//! Demonstrates how to enqueue custom draw commands in a render phase.
//!
//! This example shows how to use the built-in
//! [`bevy_render::render_phase::BinnedRenderPhase`] functionality with a
//! custom [`RenderCommand`] to allow inserting arbitrary GPU drawing logic
//! into Bevy's pipeline. This is not the only way to add custom rendering code
//! into Bevy—render nodes are another, lower-level method—but it does allow
//! for better reuse of parts of Bevy's built-in mesh rendering logic.

use bevy_app::{App, Plugin};
use bevy_asset::{AssetId, Handle};
use bevy_core_pipeline::core_3d::{Transparent3d, CORE_3D_DEPTH_FORMAT};
use bevy_ecs::{
    prelude::*,
    query::ROQueryItem,
    system::{
        lifetimeless::{Read, SRes},
        SystemParamItem,
    },
};
use bevy_render::{
    mesh::{
        allocator::MeshAllocator, MeshVertexBufferLayouts, PrimitiveTopology, RenderMesh,
        RenderMeshBufferInfo,
    },
    prelude::{Mesh, Mesh3d, Msaa},
    render_asset::RenderAssets,
    render_phase::{
        AddRenderCommand, DrawFunctions, PhaseItem, RenderCommand, RenderCommandResult,
        SetItemPipeline, TrackedRenderPass,
    },
    render_phase::{PhaseItemExtraIndex, ViewSortedRenderPhases},
    render_resource::{
        binding_types::uniform_buffer, BindGroup, BindGroupLayout, BindGroupLayoutEntry,
        DynamicBindGroupEntries, DynamicBindGroupLayoutEntries, Face, ShaderStages,
    },
    render_resource::{
        ColorTargetState, ColorWrites, CompareFunction, DepthStencilState, FragmentState,
        MultisampleState, PipelineCache, PrimitiveState, RenderPipelineDescriptor, Shader,
        SpecializedRenderPipeline, SpecializedRenderPipelines, TextureFormat, VertexState,
    },
    renderer::RenderDevice,
    sync_world::MainEntity,
    view::{ExtractedView, RenderVisibleEntities},
    view::{ViewUniform, ViewUniforms},
    Extract, ExtractSchedule, Render, RenderApp, RenderSet,
};

use std::collections::HashMap;

#[derive(Component)]
pub struct CustomEntity;

#[derive(Clone)]
pub struct MeshInstance {
    id: AssetId<Mesh>,
}

#[derive(Clone, Default, Resource)]
pub struct MeshInstances(pub HashMap<MainEntity, MeshInstance>);

/// Holds a reference to our shader.
///
/// This is loaded at app creation time.
#[derive(Resource)]
struct CustomMeshPipeline {
    shader: Handle<Shader>,
    view_layout: BindGroupLayout,
}

/// A [`RenderCommand`] that binds the vertex and index buffers and issues the
/// draw command for our custom phase item.
struct DrawCustomMesh;

/// FIXME: prepare the view uniforms before passing them to the pass ?
/// https://docs.rs/bevy_pbr/0.15.1/src/bevy_pbr/render/mesh_view_bindings.rs.html#715
/// I don't know why SRes does not work with ViewUniforms
impl<P> RenderCommand<P> for DrawCustomMesh
where
    P: PhaseItem,
{
    type Param = (
        SRes<RenderAssets<RenderMesh>>,
        SRes<MeshInstances>,
        SRes<MeshAllocator>,
    );

    type ViewQuery = Read<ViewBindGroup>;

    type ItemQuery = ();

    fn render<'w>(
        item: &P,
        views: ROQueryItem<'w, Self::ViewQuery>,
        _: Option<ROQueryItem<'w, Self::ItemQuery>>,
        (meshes, mesh_instances, mesh_allocator): SystemParamItem<'w, '_, Self::Param>,
        pass: &mut TrackedRenderPass<'w>,
    ) -> RenderCommandResult {
        let Some(MeshInstance { id, .. }) = mesh_instances.0.get(&item.main_entity()) else {
            return RenderCommandResult::Skip;
        };

        let mesh_allocator = mesh_allocator.into_inner();

        let Some(gpu_mesh) = meshes.into_inner().get(*id) else {
            return RenderCommandResult::Skip;
        };
        let Some(vertex_buffer_slice) = mesh_allocator.mesh_vertex_slice(&id) else {
            return RenderCommandResult::Skip;
        };

        let RenderMeshBufferInfo::Indexed {
            index_format,
            count,
        } = &gpu_mesh.buffer_info
        else {
            return RenderCommandResult::Failure("mesh buffer is not indexed wtf");
        };

        let index_buffer_slice = mesh_allocator.mesh_index_slice(&id).unwrap();

        // Tell the GPU where the vertices are.
        pass.set_vertex_buffer(0, vertex_buffer_slice.buffer.slice(..));

        pass.set_index_buffer(index_buffer_slice.buffer.slice(..), 0, *index_format);

        pass.set_bind_group(0, &views.0, &[0]);

        pass.draw_indexed(
            index_buffer_slice.range.start..(index_buffer_slice.range.start + count),
            vertex_buffer_slice.range.start as i32,
            0..1,
        );

        RenderCommandResult::Success
    }
}

/// The custom draw commands that Bevy executes for each entity we enqueue into
/// the render phase.
type DrawCustomMeshCommands = (SetItemPipeline, DrawCustomMesh);

pub struct CustomMeshPipelinePlugin {
    pub shader: Handle<Shader>,
}

#[derive(Resource)]
struct CustomMeshPipelineShader(Handle<Shader>);

impl Plugin for CustomMeshPipelinePlugin {
    fn build(&self, app: &mut App) {
        let Some(render_app) = app.get_sub_app_mut(RenderApp) else {
            return;
        };
        // We make sure to add these to the render app, not the main app.
        render_app
            .insert_resource(CustomMeshPipelineShader(self.shader.clone()))
            .init_resource::<SpecializedRenderPipelines<CustomMeshPipeline>>()
            .init_resource::<MeshInstances>()
            .add_render_command::<Transparent3d, DrawCustomMeshCommands>()
            .add_systems(ExtractSchedule, extract_meshes)
            .add_systems(
                Render,
                (
                    prepare_view_bind_groups.in_set(RenderSet::Prepare),
                    queue_meshes.in_set(RenderSet::Queue),
                ),
            );
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

/// A render-world system that enqueues the entity with custom rendering into
/// the opaque render phases of each view.
fn queue_meshes(
    pipeline_cache: Res<PipelineCache>,
    custom_mesh_pipeline: Res<CustomMeshPipeline>,
    mut transparent_render_phases: ResMut<ViewSortedRenderPhases<Transparent3d>>,
    transparent_draw_functions: Res<DrawFunctions<Transparent3d>>,
    mut specialized_render_pipelines: ResMut<SpecializedRenderPipelines<CustomMeshPipeline>>,
    views: Query<(Entity, &RenderVisibleEntities, &Msaa), With<ExtractedView>>,
) {
    let draw_function = transparent_draw_functions
        .read()
        .id::<DrawCustomMeshCommands>();

    // Render phases are per-view, so we need to iterate over all views so that
    // the entity appears in them. (In this example, we have only one view, but
    // it's good practice to loop over all views anyway.)
    for (view_entity, view_visible_entities, msaa) in views.iter() {
        let Some(transparent_phase) = transparent_render_phases.get_mut(&view_entity) else {
            continue;
        };

        // Find all the custom rendered entities that are visible from this
        // view.
        for &entity in view_visible_entities.get::<With<Mesh3d>>().iter() {
            // Ordinarily, the [`SpecializedRenderPipeline::Key`] would contain
            // some per-view settings, such as whether the view is HDR, but for
            // simplicity's sake we simply hard-code the view's characteristics,
            // with the exception of number of MSAA samples.
            let pipeline = specialized_render_pipelines.specialize(
                &pipeline_cache,
                &custom_mesh_pipeline,
                *msaa,
            );

            // Add the custom render item. We use the
            // [`BinnedRenderPhaseType::NonMesh`] type to skip the special
            // handling that Bevy has for meshes (preprocessing, indirect
            // draws, etc.)
            //
            // The asset ID is arbitrary; we simply use [`AssetId::invalid`],
            // but you can use anything you like. Note that the asset ID need
            // not be the ID of a [`Mesh`].

            transparent_phase.add(Transparent3d {
                entity,
                draw_function,
                pipeline,
                distance: 0.,
                batch_range: 0..1,
                extra_index: PhaseItemExtraIndex::NONE,
            });
        }
    }
}

impl SpecializedRenderPipeline for CustomMeshPipeline {
    type Key = Msaa;

    fn specialize(&self, msaa: Self::Key) -> RenderPipelineDescriptor {
        // TODO: lazy static and check the layout
        let mut l: MeshVertexBufferLayouts = Default::default();
        let vertex_layout = Mesh::new(PrimitiveTopology::TriangleList, Default::default())
            .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, vec![[0., 0., 0.]])
            .with_inserted_attribute(Mesh::ATTRIBUTE_NORMAL, vec![[0., 0., 0.]])
            .with_inserted_attribute(Mesh::ATTRIBUTE_COLOR, vec![[0., 0., 0., 0.]])
            .get_mesh_vertex_buffer_layout(&mut l)
            .0
            .layout()
            .clone();

        RenderPipelineDescriptor {
            label: Some("custom render pipeline".into()),
            layout: vec![self.view_layout.clone()],
            push_constant_ranges: vec![],
            vertex: VertexState {
                shader: self.shader.clone(),
                shader_defs: vec![],
                entry_point: "vertex".into(),
                buffers: vec![vertex_layout],
            },
            fragment: Some(FragmentState {
                shader: self.shader.clone(),
                shader_defs: vec![],
                entry_point: "fragment".into(),
                targets: vec![Some(ColorTargetState {
                    // Ordinarily, you'd want to check whether the view has the
                    // HDR format and substitute the appropriate texture format
                    // here, but we omit that for simplicity.
                    format: TextureFormat::Rgba8UnormSrgb,
                    blend: None,
                    write_mask: ColorWrites::ALL,
                })],
            }),
            primitive: PrimitiveState {
                cull_mode: Some(Face::Back),
                ..PrimitiveState::default()
            },
            // Note that if your view has no depth buffer this will need to be
            // changed.
            depth_stencil: Some(DepthStencilState {
                format: CORE_3D_DEPTH_FORMAT,
                depth_write_enabled: true,
                depth_compare: CompareFunction::Greater,
                stencil: Default::default(),
                bias: Default::default(),
            }),
            multisample: MultisampleState {
                count: msaa.samples(),
                mask: !0,
                alpha_to_coverage_enabled: false,
            },
            zero_initialize_workgroup_memory: false,
        }
    }
}

#[derive(Component)]
struct ViewBindGroup(BindGroup);

// FIXME: init view uniforms
fn prepare_view_bind_groups(
    mut commands: Commands,
    render_device: Res<RenderDevice>,
    mesh_pipeline: Res<CustomMeshPipeline>,
    view_uniforms: Res<ViewUniforms>,
    views: Query<Entity>,
) {
    if let Some(view_binding) = view_uniforms.uniforms.binding() {
        for entity in &views {
            let entries = DynamicBindGroupEntries::new_with_indices(((0, view_binding.clone()),));
            let layout = &mesh_pipeline.view_layout;

            commands
                .entity(entity)
                .insert(ViewBindGroup(render_device.create_bind_group(
                    "view_bind_group",
                    layout,
                    &entries,
                )));
        }
    }
}

fn extract_meshes(
    mut mesh_instances: ResMut<MeshInstances>,
    meshes: Extract<Query<(Entity, &Mesh3d), With<CustomEntity>>>,
) {
    mesh_instances.0.clear();
    for (entity, m) in &meshes {
        mesh_instances
            .0
            .insert(entity.into(), MeshInstance { id: m.0.id() });
    }
}

impl FromWorld for CustomMeshPipeline {
    fn from_world(world: &mut World) -> Self {
        // Load and compile the shader in the background.
        let render_device = world.resource::<RenderDevice>();
        let shader = world.resource::<CustomMeshPipelineShader>();

        let view_layout_entries: Vec<BindGroupLayoutEntry> =
            DynamicBindGroupLayoutEntries::new_with_indices(
                ShaderStages::FRAGMENT,
                ((
                    0,
                    uniform_buffer::<ViewUniform>(true).visibility(ShaderStages::VERTEX_FRAGMENT),
                ),),
            )
            .to_vec();

        let view_layout =
            render_device.create_bind_group_layout("mesh_view_layout", &view_layout_entries);

        CustomMeshPipeline {
            view_layout,
            shader: shader.0.clone(),
        }
    }
}
