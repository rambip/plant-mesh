//! Demonstrates how to define and use specialized mesh pipeline
//!
//! This example shows how to use the built-in [`SpecializedMeshPipeline`]
//! functionality with a custom [`RenderCommand`] to allow custom mesh rendering with
//! more flexibility than the material api.
//!
//! [`SpecializedMeshPipeline`] let's you customize the entire pipeline used when rendering a mesh.

use bevy::core_pipeline::oit::OrderIndependentTransparencySettingsOffset;
use bevy::core_pipeline::prepass::MotionVectorPrepass;
use bevy::utils::Parallel;
use bevy_ecs::system::lifetimeless::Read;
use bevy_pbr::{ExtractMeshesSet, MeshBindGroups, MeshCullingDataBuffer, MeshFlags, MeshInputUniform, MeshTransforms, MeshUniform, MeshViewBindGroup, NotShadowCaster, NotShadowReceiver, PreviousGlobalTransform, RenderLightmaps, RenderMeshInstanceGpuQueues, SkinIndices, TransmittedShadowReceiver};
use bevy_render::render_phase::AddRenderCommand;
use bevy_ecs::system::SystemState;
use bevy_render::render_resource::{BindGroupLayout, ShaderType};
use bevy_render::renderer::{RenderDevice, RenderQueue};
use bevy_render::texture::DefaultImageSampler;
//mod mesh_pipeline;
use bevy_render::batching::{gpu_preprocessing, GetBatchData, GetFullBatchData, NoAutomaticBatching};
//use bevy_pbr::{
//    MeshLayouts, MeshPipelineKey, MeshPipelineViewLayoutKey, MeshPipelineViewLayouts, SetMeshViewBindGroup
//};

use bevy_pbr::{MeshPipeline, MeshPipelineKey, MeshPipelineViewLayoutKey};

use bevy::{
    core_pipeline::core_3d::{Opaque3d, Opaque3dBinKey, CORE_3D_DEPTH_FORMAT},
    prelude::*,
    render::{
        extract_component::{ExtractComponent, ExtractComponentPlugin},
        mesh::{MeshVertexBufferLayoutRef, RenderMesh, RenderMeshBufferInfo},
        render_asset::RenderAssets,
        render_phase::{
            BinnedRenderPhaseType, DrawFunctions, SetItemPipeline,
            ViewBinnedRenderPhases, TrackedRenderPass, RenderCommand, RenderCommandResult,
            PhaseItem,
        },
        render_resource::{
            ColorTargetState, ColorWrites, CompareFunction, DepthStencilState, Face, FragmentState,
            FrontFace, MultisampleState, PipelineCache, PolygonMode, PrimitiveState,
            RenderPipelineDescriptor, SpecializedMeshPipeline, SpecializedMeshPipelineError,
            SpecializedMeshPipelines, TextureFormat, VertexState,
        },
        view::{self, ExtractedView, RenderVisibleEntities, ViewTarget, VisibilitySystems},
        Render, RenderApp, RenderSet,
    },
    math::Affine3,
};

use bevy_ecs::{
    system::{
        lifetimeless::SRes,
        SystemParamItem
    },
    query::ROQueryItem,
};
use bevy_render::view::{RenderVisibilityRanges, ViewUniformOffset, VisibilityRange};
use bevy_render::Extract;
use bevy_render::{
    batching::gpu_preprocessing::{IndirectParameters, IndirectParametersBuffer},
    mesh::allocator::MeshAllocator,
    sync_world::{MainEntity, MainEntityHashMap}
};

use bytemuck::{Pod, Zeroable};
use smallvec::{SmallVec, smallvec};
use nonmax::NonMaxU32;

/// Information that the render world keeps about each entity that contains a
/// mesh, when using GPU mesh instance data building.
#[derive(Default, Deref, DerefMut, Resource)]
pub struct RenderMeshInstancesCpu(MainEntityHashMap<RenderMeshInstanceCpu>);

pub struct RenderMeshQueueData {
    /// The translation of the mesh instance.
    pub mesh_asset_id: AssetId<Mesh>,
}

impl RenderMeshInstancesCpu {
    fn mesh_asset_id(&self, entity: MainEntity) -> Option<AssetId<Mesh>> {
        self.get(&entity)
            .map(|render_mesh_instance| render_mesh_instance.mesh_asset_id)
    }

    fn render_mesh_queue_data(&self, entity: MainEntity) -> Option<RenderMeshQueueData> {
        self.get(&entity)
            .map(|render_mesh_instance| RenderMeshQueueData {
                mesh_asset_id: render_mesh_instance.mesh_asset_id,
            })
    }
}

//#[derive(Component)]
//pub struct MeshTransforms {
//    pub world_from_local: Affine3,
//    pub previous_world_from_local: Affine3,
//}
//
/// CPU data that the render world keeps for each entity, when *not* using GPU
/// mesh uniform building.
pub struct RenderMeshInstanceCpu {
    /// The transform of the mesh.
    ///
    /// This will be written into the [`MeshUniform`] at the appropriate time.
    pub transforms: MeshTransforms,
    pub mesh_asset_id: AssetId<Mesh>,
}

/// Information that the render world keeps about each entity that contains a
/// mesh.
///
/// The set of information needed is different depending on whether CPU or GPU
/// [`MeshUniform`] building is in use.
#[derive(Resource, Default)]
pub struct RenderMeshInstances(
    /// Information needed when using CPU mesh instance data building.
    RenderMeshInstancesCpu,
);

pub struct DrawMesh;
impl RenderCommand<Opaque3d> for DrawMesh {
    type Param = (
        SRes<RenderAssets<RenderMesh>>,
        SRes<RenderMeshInstances>,
        SRes<IndirectParametersBuffer>,
        SRes<MeshAllocator>,
    );
    type ViewQuery = ();
    type ItemQuery = ();
    #[inline]
    fn render<'w>(
        item: &Opaque3d,
        _: ROQueryItem<Self::ViewQuery>,
        _item_query: Option<()>,
        (
            meshes,
            mesh_instances,
            indirect_parameters_buffer,
            mesh_allocator,
        ): SystemParamItem<'w, '_, Self::Param>,
        pass: &mut TrackedRenderPass<'w>,
    ) -> RenderCommandResult {
        // If we're using GPU preprocessing, then we're dependent on that
        // compute shader having been run, which of course can only happen if
        // it's compiled. Otherwise, our mesh instance data won't be present.

        let meshes = meshes.into_inner();
        let mesh_instances = mesh_instances.into_inner();
        let indirect_parameters_buffer = indirect_parameters_buffer.into_inner();
        let mesh_allocator = mesh_allocator.into_inner();

        let Some(mesh_asset_id) = 
            mesh_instances.0.mesh_asset_id(item.main_entity())
            else {
                return RenderCommandResult::Skip;
        }
            ;

    let Some(gpu_mesh) = meshes.get(mesh_asset_id) else {
        return RenderCommandResult::Skip;
    };
    let Some(vertex_buffer_slice) = mesh_allocator.mesh_vertex_slice(&mesh_asset_id) else {
        return RenderCommandResult::Skip;
    };

        // Calculate the indirect offset, and look up the buffer.
        let indirect_parameters = match item.extra_index().as_indirect_parameters_index() {
            None => None,
            Some(index) => match indirect_parameters_buffer.buffer() {
                None => {
                    warn!("Not rendering mesh because indirect parameters buffer wasn't present");
                    return RenderCommandResult::Skip;
                }
                Some(buffer) => Some((
                    index as u64 * size_of::<IndirectParameters>() as u64,
                    buffer,
                )),
            },
        };

        pass.set_vertex_buffer(0, vertex_buffer_slice.buffer.slice(..));

        let batch_range = item.batch_range();

        // Draw either directly or indirectly, as appropriate.
        match &gpu_mesh.buffer_info {
            RenderMeshBufferInfo::Indexed {
                index_format,
                count,
            } => {
                let Some(index_buffer_slice) = mesh_allocator.mesh_index_slice(&mesh_asset_id)
                else {
                    return RenderCommandResult::Skip;
                };

                pass.set_index_buffer(index_buffer_slice.buffer.slice(..), 0, *index_format);

                match indirect_parameters {
                    None => {
                        pass.draw_indexed(
                            index_buffer_slice.range.start
                                ..(index_buffer_slice.range.start + *count),
                            vertex_buffer_slice.range.start as i32,
                            batch_range.clone(),
                        );
                    }
                    Some((indirect_parameters_offset, indirect_parameters_buffer)) => pass
                        .draw_indexed_indirect(
                            indirect_parameters_buffer,
                            indirect_parameters_offset,
                        ),
                }
            }
            RenderMeshBufferInfo::NonIndexed => match indirect_parameters {
                None => {
                    pass.draw(vertex_buffer_slice.range, batch_range.clone());
                }
                Some((indirect_parameters_offset, indirect_parameters_buffer)) => {
                    pass.draw_indirect(indirect_parameters_buffer, indirect_parameters_offset);
                }
            },
        }
        RenderCommandResult::Success
    }
}

pub struct SetMeshViewBindGroup<const I: usize>;
impl<P: PhaseItem, const I: usize> RenderCommand<P> for SetMeshViewBindGroup<I> {
    type Param = ();
    type ViewQuery = (
        Read<ViewUniformOffset>,
        Read<MeshViewBindGroup>,
        Option<Read<OrderIndependentTransparencySettingsOffset>>,
    );
    type ItemQuery = ();

    #[inline]
    fn render<'w>(
        _item: &P,
        (
            view_uniform,
            mesh_view_bind_group,
            maybe_oit_layers_count_offset,
        ): ROQueryItem<'w, Self::ViewQuery>,
        _entity: Option<()>,
        _: SystemParamItem<'w, '_, Self::Param>,
        pass: &mut TrackedRenderPass<'w>,
    ) -> RenderCommandResult {
        let mut offsets: SmallVec<[u32; 8]> = smallvec![
            view_uniform.offset,
        ];
        if let Some(layers_count_offset) = maybe_oit_layers_count_offset {
            offsets.push(layers_count_offset.offset);
        }
        pass.set_bind_group(I, &mesh_view_bind_group.value, &offsets);

        RenderCommandResult::Success
    }
}

pub struct SetMeshBindGroup<const I: usize>;
impl<P: PhaseItem, const I: usize> RenderCommand<P> for SetMeshBindGroup<I> {
    type Param = (
        SRes<MeshBindGroups>,
        SRes<RenderMeshInstances>,
    );
    type ViewQuery = Has<MotionVectorPrepass>;
    type ItemQuery = ();

    #[inline]
    fn render<'w>(
        item: &P,
        has_motion_vector_prepass: bool,
        _item_query: Option<()>,
        (bind_groups, mesh_instances): SystemParamItem<
            'w,
            '_,
            Self::Param,
        >,
        pass: &mut TrackedRenderPass<'w>,
    ) -> RenderCommandResult {
        let bind_groups = bind_groups.into_inner();
        let mesh_instances = mesh_instances.into_inner();

        let entity = &item.main_entity();

        let Some(mesh_asset_id) = mesh_instances.0.mesh_asset_id(*entity) else {
            return RenderCommandResult::Success;
        };

        let Some(bind_group) = bind_groups.get(
            mesh_asset_id,
            None,
            Default::default(),
            Default::default(),
            has_motion_vector_prepass,
        ) else {
            return RenderCommandResult::Failure(
                "The MeshBindGroups resource wasn't set in the render phase. \
                It should be set by the prepare_mesh_bind_group system.\n\
                This is a bevy bug! Please open an issue.",
            );
        };

        let mut dynamic_offsets: [u32; 3] = Default::default();
        let mut offset_count = 0;
        if let Some(dynamic_offset) = item.extra_index().as_dynamic_offset() {
            dynamic_offsets[offset_count] = dynamic_offset.get();
            offset_count += 1;
        }

        pass.set_bind_group(I, bind_group, &dynamic_offsets[0..offset_count]);

        RenderCommandResult::Success
    }
}

//#[derive(ShaderType, Clone)]
//pub struct MeshUniform {
//    // Affine 4x3 matrices transposed to 3x4
//    pub world_from_local: [Vec4; 3],
//    pub previous_world_from_local: [Vec4; 3],
//    // 3x3 matrix packed in mat2x4 and f32 as:
//    //   [0].xyz, [1].x,
//    //   [1].yz, [2].xy
//    //   [2].z
//    pub local_from_world_transpose_a: [Vec4; 2],
//    pub local_from_world_transpose_b: f32,
//    pub flags: u32,
//    // Four 16-bit unsigned normalized UV values packed into a `UVec2`:
//    //
//    //                         <--- MSB                   LSB --->
//    //                         +---- min v ----+ +---- min u ----+
//    //     lightmap_uv_rect.x: vvvvvvvv vvvvvvvv uuuuuuuu uuuuuuuu,
//    //                         +---- max v ----+ +---- max u ----+
//    //     lightmap_uv_rect.y: VVVVVVVV VVVVVVVV UUUUUUUU UUUUUUUU,
//    //
//    // (MSB: most significant bit; LSB: least significant bit.)
//    pub lightmap_uv_rect: UVec2,
//    /// The index of this mesh's first vertex in the vertex buffer.
//    ///
//    /// Multiple meshes can be packed into a single vertex buffer (see
//    /// [`MeshAllocator`]). This value stores the offset of the first vertex in
//    /// this mesh in that buffer.
//    pub first_vertex_index: u32,
//    /// Padding.
//    pub pad_a: u32,
//    /// Padding.
//    pub pad_b: u32,
//    /// Padding.
//    pub pad_c: u32,
//}
//
//impl MeshUniform {
//    pub fn new(
//        mesh_transforms: &MeshTransforms,
//        first_vertex_index: u32,
//    ) -> Self {
//        let (local_from_world_transpose_a, local_from_world_transpose_b) =
//            mesh_transforms.world_from_local.inverse_transpose_3x3();
//        Self {
//            world_from_local: mesh_transforms.world_from_local.to_transpose(),
//            previous_world_from_local: mesh_transforms.previous_world_from_local.to_transpose(),
//            local_from_world_transpose_a,
//            local_from_world_transpose_b,
//            lightmap_uv_rect: UVec2::ZERO,
//            flags: Default::default(),
//            first_vertex_index,
//            pad_a: 0,
//            pad_b: 0,
//            pad_c: 0,
//        }
//    }
//}
//
///// Information that has to be transferred from CPU to GPU in order to produce
///// the full [`MeshUniform`].
/////
///// This is essentially a subset of the fields in [`MeshUniform`] above.
//#[derive(ShaderType, Pod, Zeroable, Clone, Copy)]
//#[repr(C)]
//pub struct MeshInputUniform {
//    /// Affine 4x3 matrix transposed to 3x4.
//    pub world_from_local: [Vec4; 3],
//    /// Four 16-bit unsigned normalized UV values packed into a `UVec2`:
//    ///
//    /// ```text
//    ///                         <--- MSB                   LSB --->
//    ///                         +---- min v ----+ +---- min u ----+
//    ///     lightmap_uv_rect.x: vvvvvvvv vvvvvvvv uuuuuuuu uuuuuuuu,
//    ///                         +---- max v ----+ +---- max u ----+
//    ///     lightmap_uv_rect.y: VVVVVVVV VVVVVVVV UUUUUUUU UUUUUUUU,
//    ///
//    /// (MSB: most significant bit; LSB: least significant bit.)
//    /// ```
//    pub lightmap_uv_rect: UVec2,
//    /// Various [`MeshFlags`].
//    pub flags: u32,
//    /// The index of this mesh's [`MeshInputUniform`] in the previous frame's
//    /// buffer, if applicable.
//    ///
//    /// This is used for TAA. If not present, this will be `u32::MAX`.
//    pub previous_input_index: u32,
//    /// The index of this mesh's first vertex in the vertex buffer.
//    ///
//    /// Multiple meshes can be packed into a single vertex buffer (see
//    /// [`MeshAllocator`]). This value stores the offset of the first vertex in
//    /// this mesh in that buffer.
//    pub first_vertex_index: u32,
//    /// Padding.
//    pub pad_a: u32,
//    /// Padding.
//    pub pad_b: u32,
//    /// Padding.
//    pub pad_c: u32,
//}
//
///// All data needed to construct a pipeline for rendering 3D meshes.
//#[derive(Resource, Clone)]
//pub struct MeshPipeline {
//    /// A reference to all the mesh pipeline view layouts.
//    pub view_layouts: MeshPipelineViewLayouts,
//    // This dummy white texture is to be used in place of optional StandardMaterial textures
//    pub mesh_layouts: MeshLayouts,
//}
//
//impl FromWorld for MeshPipeline {
//    fn from_world(world: &mut World) -> Self {
//        let mut system_state: SystemState<(
//            Res<RenderDevice>,
//            Res<DefaultImageSampler>,
//            Res<RenderQueue>,
//            Res<MeshPipelineViewLayouts>,
//        )> = SystemState::new(world);
//        let (render_device, _, _, view_layouts) =
//            system_state.get_mut(world);
//
//
//        MeshPipeline {
//            view_layouts: view_layouts.clone(),
//            mesh_layouts: MeshLayouts::new(&render_device),
//        }
//    }
//}
//
//impl MeshPipeline {
//    pub fn get_view_layout(&self, layout_key: MeshPipelineViewLayoutKey) -> &BindGroupLayout {
//        self.view_layouts.get_view_layout(layout_key)
//    }
//}
//
//impl GetBatchData for MeshPipeline {
//    type Param = (
//        SRes<RenderMeshInstancesCpu>,
//        SRes<RenderAssets<RenderMesh>>,
//        SRes<MeshAllocator>,
//    );
//    // The material bind group ID, the mesh ID, and the lightmap ID,
//    // respectively.
//    type CompareData = AssetId<Mesh>;
//
//    type BufferData = MeshUniform;
//
//    fn get_batch_data(
//        (mesh_instances, _, mesh_allocator): &SystemParamItem<Self::Param>,
//        (_entity, main_entity): (Entity, MainEntity),
//    ) -> Option<(Self::BufferData, Option<Self::CompareData>)> {
//        let mesh_instance = mesh_instances.get(&main_entity)?;
//        let first_vertex_index =
//            match mesh_allocator.mesh_vertex_slice(&mesh_instance.mesh_asset_id) {
//                Some(mesh_vertex_slice) => mesh_vertex_slice.range.start,
//                None => 0,
//            };
//        Some((
//            MeshUniform::new(
//                &mesh_instance.transforms,
//                first_vertex_index,
//            ),
//            Some(
//                mesh_instance.mesh_asset_id,
//            ),
//        ))
//    }
//}
//
//impl GetFullBatchData for MeshPipeline {
//    type BufferInputData = MeshInputUniform;
//
//    fn get_index_and_compare_data(
//        _: &SystemParamItem<Self::Param>,
//        _: (Entity, MainEntity),
//    ) -> Option<(NonMaxU32, Option<Self::CompareData>)> {
//            error!(
//                "`get_index_and_compare_data` should never be called in CPU mesh uniform building \
//                mode"
//            );
//            None
//    }
//
//    fn get_binned_batch_data(
//        (mesh_instances, _, mesh_allocator): &SystemParamItem<Self::Param>,
//        (_entity, main_entity): (Entity, MainEntity),
//    ) -> Option<Self::BufferData> {
//        let mesh_instance = mesh_instances.get(&main_entity)?;
//        let first_vertex_index =
//            match mesh_allocator.mesh_vertex_slice(&mesh_instance.mesh_asset_id) {
//                Some(mesh_vertex_slice) => mesh_vertex_slice.range.start,
//                None => 0,
//            };
//
//        Some(MeshUniform::new(
//            &mesh_instance.transforms,
//            first_vertex_index,
//        ))
//    }
//
//    fn get_binned_index(
//        _: &SystemParamItem<Self::Param>,
//        _: (Entity, MainEntity),
//    ) -> Option<NonMaxU32> {
//        // This should only be called during GPU building.
//            error!(
//                "`get_binned_index` should never be called in CPU mesh uniform \
//                building mode"
//            );
//        None
//    }
//
//    fn get_batch_indirect_parameters_index(
//        _: &SystemParamItem<Self::Param>,
//        _: &mut IndirectParametersBuffer,
//        _: (Entity, MainEntity),
//        _: u32,
//    ) -> Option<NonMaxU32> {
//            error!(
//                "`get_batch_indirect_parameters_index` should never be called in CPU mesh uniform \
//                building mode"
//            );
//            None
//    }
//}

/// Extracts meshes from the main world into the render world, populating the
/// [`RenderMeshInstances`].
///
/// This is the variant of the system that runs when we're *not* using GPU
/// [`MeshUniform`] building.
pub fn extract_meshes_for_cpu_building(
    mut render_mesh_instances: ResMut<RenderMeshInstances>,
    render_visibility_ranges: Res<RenderVisibilityRanges>,
    mut render_mesh_instance_queues: Local<Parallel<Vec<(Entity, RenderMeshInstanceCpu)>>>,
    meshes_query: Extract<
        Query<(
            Entity,
            &ViewVisibility,
            &GlobalTransform,
            Option<&PreviousGlobalTransform>,
            &Mesh3d,
            Has<NotShadowReceiver>,
            Has<TransmittedShadowReceiver>,
            Has<NotShadowCaster>,
            Has<NoAutomaticBatching>,
            Has<VisibilityRange>,
        )>,
    >,
) {
    meshes_query.par_iter().for_each_init(
        || render_mesh_instance_queues.borrow_local_mut(),
        |queue,
         (
            entity,
            view_visibility,
            transform,
            previous_transform,
            mesh,
            not_shadow_receiver,
            transmitted_receiver,
            not_shadow_caster,
            no_automatic_batching,
            visibility_range,
        )| {
            if !view_visibility.get() {
                return;
            }

            let mut lod_index = None;
            if visibility_range {
                lod_index = render_visibility_ranges.lod_index_for_entity(entity.into());
            }

            let mesh_flags = MeshFlags::NONE;

            let world_from_local = transform.affine();
            queue.push((
                entity,
                RenderMeshInstanceCpu {
                    transforms: MeshTransforms {
                        world_from_local: (&world_from_local).into(),
                        previous_world_from_local: (&previous_transform
                            .map(|t| t.0)
                            .unwrap_or(world_from_local))
                            .into(),
                        flags: mesh_flags.bits(),
                    },
                    mesh_asset_id: mesh.id(),
                },
            ));
        },
    );

    // Collect the render mesh instances.
    render_mesh_instances.0.clear();
    for queue in render_mesh_instance_queues.iter_mut() {
        for (entity, render_mesh_instance) in queue.drain(..) {
            render_mesh_instances.0.insert_unique_unchecked(entity.into(), render_mesh_instance);

        }
    }
}

pub const MESH_SHADER_HANDLE: Handle<Shader> =
    Handle::weak_from_u128(13828845428412094821);

const MESH_SHADER: &str = include_str!("shader.wgsl");

// When writing custom rendering code it's generally recommended to use a plugin.
// The main reason for this is that it gives you access to the finish() hook
// which is called after rendering resources are initialized.
pub struct CustomRenderedMeshPipelinePlugin;
impl Plugin for CustomRenderedMeshPipelinePlugin {
    fn build(&self, app: &mut App) {
        let mut shaders = app.world_mut().resource_mut::<Assets<Shader>>();
        shaders.insert(
            &MESH_SHADER_HANDLE,
            Shader::from_wgsl(MESH_SHADER, file!()),
        );
        app.add_plugins(ExtractComponentPlugin::<CustomRenderedEntity>::default())
            .add_systems(
                PostUpdate,
                // Make sure to tell Bevy to check our entity for visibility. Bevy won't
                // do this by default, for efficiency reasons.
                // This will do things like frustum culling and hierarchy visibility
                view::check_visibility::<WithCustomRenderedEntity>
                .in_set(VisibilitySystems::CheckVisibility),
            )

            .add_systems(
                ExtractSchedule,
                extract_meshes_for_cpu_building.in_set(ExtractMeshesSet),
            )
            ;
        // We make sure to add these to the render app, not the main app.
        let Some(render_app) = app.get_sub_app_mut(RenderApp) else {
            return;
        };
        render_app
            // This is needed to tell bevy about your custom pipeline
            .init_resource::<SpecializedMeshPipelines<CustomMeshPipeline>>()
            // We need to use a custom draw command so we need to register it
            .add_render_command::<Opaque3d, DrawSpecializedPipelineCommands>()
            .add_systems(Render, queue_custom_mesh_pipeline.in_set(RenderSet::Queue));
    }

    fn finish(&self, app: &mut App) {
        let Some(render_app) = app.get_sub_app_mut(RenderApp) else {
            return;
        };
        // Creating this pipeline needs the RenderDevice and RenderQueue
        // which are only available once rendering plugins are initialized.
        render_app
            //.insert_resource(RenderMeshInstancesCpu::default())
            .insert_resource(RenderMeshInstances::default())
            .init_resource::<CustomMeshPipeline>();
    }
}

/// A marker component that represents an entity that is to be rendered using
/// our specialized pipeline.
///
/// Note the [`ExtractComponent`] trait implementation. This is necessary to
/// tell Bevy that this object should be pulled into the render world.
#[derive(Clone, Component, ExtractComponent)]
pub struct CustomRenderedEntity;

/// The custom draw commands that Bevy executes for each entity we enqueue into
/// the render phase.
type DrawSpecializedPipelineCommands = (
    // Set the pipeline
    SetItemPipeline,
    // Set the view uniform at bind group 0
    SetMeshViewBindGroup<0>,
    // Set the mesh uniform at bind group 1
    SetMeshBindGroup<1>,
    // Draw the mesh
    DrawMesh,
);

/// A query filter that tells [`view::check_visibility`] about our custom
/// rendered entity.
type WithCustomRenderedEntity = With<CustomRenderedEntity>;

#[derive(Resource)]
struct CustomMeshPipeline(MeshPipeline);

impl FromWorld for CustomMeshPipeline {
    fn from_world(world: &mut World) -> Self {
        // Load the shader
        Self( MeshPipeline::from_world(world),)
    }
}

impl SpecializedMeshPipeline for CustomMeshPipeline {
    /// Pipeline use keys to determine how to specialize it.
    /// The key is also used by the pipeline cache to determine if
    /// it needs to create a new pipeline or not
    ///
    /// In this example we just use the base `MeshPipelineKey` defined by bevy, but this could be anything.
    /// For example, if you want to make a pipeline with a procedural shader you could add the Handle<Shader> to the key.
    type Key = MeshPipelineKey;

    fn specialize(
        &self,
        mesh_key: Self::Key,
        layout: &MeshVertexBufferLayoutRef,
    ) -> Result<RenderPipelineDescriptor, SpecializedMeshPipelineError> {
        // Define the vertex attributes based on a standard bevy [`Mesh`]
        let mut vertex_attributes = Vec::new();
        if layout.0.contains(Mesh::ATTRIBUTE_POSITION) {
            // Make sure this matches the shader location
            vertex_attributes.push(Mesh::ATTRIBUTE_POSITION.at_shader_location(0));
        }
        if layout.0.contains(Mesh::ATTRIBUTE_COLOR) {
            // Make sure this matches the shader location
            vertex_attributes.push(Mesh::ATTRIBUTE_COLOR.at_shader_location(1));
        }
        if layout.0.contains(Mesh::ATTRIBUTE_NORMAL) {
            // Make sure this matches the shader location
            vertex_attributes.push(Mesh::ATTRIBUTE_NORMAL.at_shader_location(2));
        }
        // This will automatically generate the correct `VertexBufferLayout` based on the vertex attributes
        let vertex_buffer_layout = layout.0.get_layout(&vertex_attributes)?;

        Ok(RenderPipelineDescriptor {
            label: Some("Specialized Mesh Pipeline".into()),
            layout: vec![
                // Bind group 0 is the view uniform
                self.0
                    .get_view_layout(MeshPipelineViewLayoutKey::from(mesh_key))
                    .clone(),
                // Bind group 1 is the mesh uniform
                self.0.mesh_layouts.model_only.clone(),
            ],
            push_constant_ranges: vec![],
            vertex: VertexState {
                shader: MESH_SHADER_HANDLE,
                shader_defs: vec![],
                entry_point: "vertex".into(),
                // Customize how to store the meshes' vertex attributes in the vertex buffer
                buffers: vec![vertex_buffer_layout],
            },
            fragment: Some(FragmentState {
                shader: MESH_SHADER_HANDLE,
                shader_defs: vec![],
                entry_point: "fragment".into(),
                targets: vec![Some(ColorTargetState {
                    // This isn't required, but bevy supports HDR and non-HDR rendering
                    // so it's generally recommended to specialize the pipeline for that
                    format: if mesh_key.contains(MeshPipelineKey::HDR) {
                        ViewTarget::TEXTURE_FORMAT_HDR
                    } else {
                        TextureFormat::bevy_default()
                    },
                    // For this example we only use opaque meshes,
                    // but if you wanted to use alpha blending you would need to set it here
                    blend: None,
                    write_mask: ColorWrites::ALL,
                })],
            }),
            primitive: PrimitiveState {
                topology: mesh_key.primitive_topology(),
                front_face: FrontFace::Ccw,
                cull_mode: Some(Face::Back),
                polygon_mode: PolygonMode::Fill,
                ..default()
            },
            // Note that if your view has no depth buffer this will need to be
            // changed.
            depth_stencil: Some(DepthStencilState {
                format: CORE_3D_DEPTH_FORMAT,
                depth_write_enabled: true,
                depth_compare: CompareFunction::GreaterEqual,
                stencil: default(),
                bias: default(),
            }),
            // It's generally recommended to specialize your pipeline for MSAA,
            // but it's not always possible
            multisample: MultisampleState {
                count: mesh_key.msaa_samples(),
                ..MultisampleState::default()
            },
            zero_initialize_workgroup_memory: false,
        })
    }
}


/// A render-world system that enqueues the entity with custom rendering into
/// the opaque render phases of each view.
#[allow(clippy::too_many_arguments)]
fn queue_custom_mesh_pipeline(
    pipeline_cache: Res<PipelineCache>,
    custom_mesh_pipeline: Res<CustomMeshPipeline>,
    mut opaque_render_phases: ResMut<ViewBinnedRenderPhases<Opaque3d>>,
    opaque_draw_functions: Res<DrawFunctions<Opaque3d>>,
    mut specialized_mesh_pipelines: ResMut<SpecializedMeshPipelines<CustomMeshPipeline>>,
    views: Query<(Entity, &RenderVisibleEntities, &ExtractedView, &Msaa), With<ExtractedView>>,
    render_meshes: Res<RenderAssets<RenderMesh>>,
    render_mesh_instances: Res<RenderMeshInstances>,
) {
    // Get the id for our custom draw function
    let draw_function_id = opaque_draw_functions
        .read()
        .id::<DrawSpecializedPipelineCommands>();

    // Render phases are per-view, so we need to iterate over all views so that
    // the entity appears in them. (In this example, we have only one view, but
    // it's good practice to loop over all views anyway.)
    for (view_entity, view_visible_entities, view, msaa) in views.iter() {
        let Some(opaque_phase) = opaque_render_phases.get_mut(&view_entity) else {
            continue;
        };

        // Create the key based on the view. In this case we only care about MSAA and HDR
        let view_key = MeshPipelineKey::from_msaa_samples(msaa.samples())
            | MeshPipelineKey::from_hdr(view.hdr);

        // Find all the custom rendered entities that are visible from this
        // view.
        for &(render_entity, visible_entity) in view_visible_entities
            .get::<WithCustomRenderedEntity>()
            .iter()
        {
            // Get the mesh instance
            let Some(mesh_instance) = render_mesh_instances.0.render_mesh_queue_data(visible_entity)
            else {
                continue;
            };

            // Get the mesh data
            let Some(mesh) = render_meshes.get(mesh_instance.mesh_asset_id) else {
                continue;
            };

            // Specialize the key for the current mesh entity
            // For this example we only specialize based on the mesh topology
            // but you could have more complex keys and that's where you'd need to create those keys
            let mut mesh_key = view_key;
            mesh_key |= MeshPipelineKey::from_primitive_topology(mesh.primitive_topology());

            // Finally, we can specialize the pipeline based on the key
            let pipeline_id = specialized_mesh_pipelines
                .specialize(
                    &pipeline_cache,
                    &custom_mesh_pipeline,
                    mesh_key,
                    &mesh.layout,
                )
                // This should never with this example, but if your pipeline specialization
                // can fail you need to handle the error here
                .expect("Failed to specialize mesh pipeline");

            // Add the mesh with our specialized pipeline
            opaque_phase.add(
                Opaque3dBinKey {
                    draw_function: draw_function_id,
                    pipeline: pipeline_id,
                    // The asset ID is arbitrary; we simply use [`AssetId::invalid`],
                    // but you can use anything you like. Note that the asset ID need
                    // not be the ID of a [`Mesh`].
                    asset_id: AssetId::<Mesh>::invalid().untyped(),
                    material_bind_group_id: None,
                    lightmap_image: None,
                },
                (render_entity, visible_entity),
                // This example supports batching, but if your pipeline doesn't
                // support it you can use `BinnedRenderPhaseType::UnbatchableMesh`
                BinnedRenderPhaseType::BatchableMesh,
            );
        }
    }
}
