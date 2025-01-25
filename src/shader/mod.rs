//! Demonstrates how to define and use specialized mesh pipeline
//!
//! This example shows how to use the built-in [`SpecializedMeshPipeline`]
//! functionality with a custom [`RenderCommand`] to allow custom mesh rendering with
//! more flexibility than the material api.
//!
//! [`SpecializedMeshPipeline`] let's you customize the entire pipeline used when rendering a mesh.

use bevy_render::{batching::gpu_preprocessing::{IndirectParameters, IndirectParametersBuffer}, mesh::{allocator::MeshAllocator, PrimitiveTopology, RenderMeshBufferInfo}, render_phase::{PhaseItem, RenderCommand, RenderCommandResult, TrackedRenderPass}, render_resource::MultisampleState, renderer::RenderAdapter, view::ViewUniform};
use bevy::{
    core_pipeline::core_3d::{Opaque3d, Opaque3dBinKey, CORE_3D_DEPTH_FORMAT},
    prelude::*,
    render::{
        extract_component::{ExtractComponent, ExtractComponentPlugin},
        mesh::{MeshVertexBufferLayoutRef, RenderMesh},
        render_asset::RenderAssets,
        render_phase::{
            AddRenderCommand, BinnedRenderPhaseType, DrawFunctions, SetItemPipeline,
            ViewBinnedRenderPhases,
        },
        render_resource::{
            ColorTargetState, ColorWrites, CompareFunction, DepthStencilState, Face, FragmentState,
            FrontFace, PipelineCache, PolygonMode, PrimitiveState,
            RenderPipelineDescriptor, SpecializedMeshPipeline, SpecializedMeshPipelineError,
            SpecializedMeshPipelines, TextureFormat, VertexState,
        },
        view::{self, ExtractedView, RenderVisibleEntities, VisibilitySystems},
        Render, RenderApp, RenderSet,
    },
};

use bevy_ecs::{query::ROQueryItem, system::{lifetimeless::SRes, SystemParamItem, SystemState}};
use bevy_pbr::{
    layout_entries, MeshPipeline, MeshPipelineKey, MeshPipelineViewLayoutKey, MeshUniform, RenderMeshInstances, SetMeshBindGroup, SetMeshViewBindGroup
};
use bevy_render::{render_resource::{BindGroupLayout, BindGroupLayoutEntries, GpuArrayBuffer, ShaderStages}, renderer::RenderDevice};

const SHADER_ASSET_PATH: &str = "shader.wgsl";


pub struct MyDrawMesh;
impl<P: PhaseItem> RenderCommand<P> for MyDrawMesh {
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
        item: &P,
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

        let Some(mesh_asset_id) = mesh_instances.mesh_asset_id(item.main_entity()) else {
            return RenderCommandResult::Skip;
        };
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

// When writing custom rendering code it's generally recommended to use a plugin.
// The main reason for this is that it gives you access to the finish() hook
// which is called after rendering resources are initialized.
pub struct CustomRenderedMeshPipelinePlugin;
impl Plugin for CustomRenderedMeshPipelinePlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(ExtractComponentPlugin::<CustomRenderedEntity>::default())
            .add_systems(
                PostUpdate,
                // Make sure to tell Bevy to check our entity for visibility. Bevy won't
                // do this by default, for efficiency reasons.
                // This will do things like frustum culling and hierarchy visibility
                view::check_visibility::<WithCustomRenderedEntity>
                    .in_set(VisibilitySystems::CheckVisibility),
            );

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
        render_app.init_resource::<CustomMeshPipeline>();
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
    MyDrawMesh,
);

pub const CLUSTERED_FORWARD_STORAGE_BUFFER_COUNT: u32 = 3;
pub const VISIBILITY_RANGES_STORAGE_BUFFER_COUNT: u32 = 4;

/// A query filter that tells [`view::check_visibility`] about our custom
/// rendered entity.
type WithCustomRenderedEntity = With<CustomRenderedEntity>;

// This contains the state needed to speciazlize a mesh pipeline
#[derive(Resource)]
struct CustomMeshPipeline {
    /// Stores the shader used for this pipeline directly on the pipeline.
    /// This isn't required, it's only done like this for simplicity.
    shader_handle: Handle<Shader>,
    view_layout: BindGroupLayout,
    model: BindGroupLayout,

}
impl FromWorld for CustomMeshPipeline {
    fn from_world(world: &mut World) -> Self {
        // Load the shader
        let shader_handle: Handle<Shader> = world.resource::<AssetServer>().load(SHADER_ASSET_PATH);

        let mut system_state: SystemState<( Res<RenderDevice>, Res<RenderAdapter>)> = SystemState::new(world);
        let (render_device, render_adapter) = system_state.get_mut(world);
        let model = render_device.create_bind_group_layout(
            "mesh_layout",
            &BindGroupLayoutEntries::single(
                ShaderStages::empty(),
                GpuArrayBuffer::<MeshUniform>::binding_layout(&render_device)
                .visibility(ShaderStages::VERTEX_FRAGMENT).clone()
            ),
        );


        let clustered_forward_buffer_binding_type = render_device
            .get_supported_read_only_binding_type(CLUSTERED_FORWARD_STORAGE_BUFFER_COUNT);
        let visibility_ranges_buffer_binding_type = render_device
            .get_supported_read_only_binding_type(VISIBILITY_RANGES_STORAGE_BUFFER_COUNT);

        let key = MeshPipelineKey::from_primitive_topology(PrimitiveTopology::TriangleList);
        let entries = layout_entries(
            clustered_forward_buffer_binding_type,
            visibility_ranges_buffer_binding_type,
            key.into(),
            &render_device,
            &render_adapter,
        );

        let view_layout = render_device.create_bind_group_layout(
            "mesh_view_layout",
            &entries.to_vec()
        );

        Self {
            shader_handle,
            view_layout,
            model
        }
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
                self.view_layout.clone(),
                // Bind group 1 is the mesh uniform
                self.model.clone(),
            ],
            push_constant_ranges: vec![],
            vertex: VertexState {
                shader: self.shader_handle.clone(),
                shader_defs: vec![],
                entry_point: "vertex".into(),
                // Customize how to store the meshes' vertex attributes in the vertex buffer
                buffers: vec![vertex_buffer_layout],
            },
            fragment: Some(FragmentState {
                shader: self.shader_handle.clone(),
                shader_defs: vec![],
                entry_point: "fragment".into(),
                targets: vec![Some(ColorTargetState {
                    // This isn't required, but bevy supports HDR and non-HDR rendering
                    // so it's generally recommended to specialize the pipeline for that
                    format: TextureFormat::bevy_default(),
                    // For this example we only use opaque meshes,
                    // but if you wanted to use alpha blending you would need to set it here
                    blend: None,
                    write_mask: ColorWrites::ALL,
                })],
            }),
            primitive: PrimitiveState {
                topology: PrimitiveTopology::TriangleList,
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

        // Find all the custom rendered entities that are visible from this
        // view.
        for &(render_entity, visible_entity) in view_visible_entities
            .get::<WithCustomRenderedEntity>()
            .iter()
        {
            // Get the mesh instance
            let Some(mesh_instance) = render_mesh_instances.render_mesh_queue_data(visible_entity)
            else {
                continue;
            };

            // Get the mesh data
            let Some(mesh) = render_meshes.get(mesh_instance.mesh_asset_id) else {
                continue;
            };

            let mesh_key = MeshPipelineKey::from_msaa_samples(msaa.samples())
            | MeshPipelineKey::from_hdr(view.hdr)
            | MeshPipelineKey::from_primitive_topology(mesh.primitive_topology());

            // Finally, we can specialize the pipeline based on the key
            let pipeline_id = specialized_mesh_pipelines
                .specialize(
                    &pipeline_cache,
                    &custom_mesh_pipeline,
                    mesh_key.into(),
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
