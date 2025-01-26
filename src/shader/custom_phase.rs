//! Demonstrates how to enqueue custom draw commands in a render phase.
//!
//! This example shows how to use the built-in
//! [`bevy_render::render_phase::BinnedRenderPhase`] functionality with a
//! custom [`RenderCommand`] to allow inserting arbitrary GPU drawing logic
//! into Bevy's pipeline. This is not the only way to add custom rendering code
//! into Bevy—render nodes are another, lower-level method—but it does allow
//! for better reuse of parts of Bevy's built-in mesh rendering logic.

use bevy::{
    asset::{AssetEvents, RenderAssetUsages}, core_pipeline::core_3d::{Opaque3d, Opaque3dBinKey, CORE_3D_DEPTH_FORMAT}, ecs::{
        query::ROQueryItem,
        system::{lifetimeless::SRes, SystemParamItem},
    }, math::{vec3, Vec3A}, prelude::*, render::{
        extract_component::{ExtractComponent, ExtractComponentPlugin},
        primitives::Aabb,
        render_phase::{
            AddRenderCommand, BinnedRenderPhaseType, DrawFunctions, PhaseItem, RenderCommand,
            RenderCommandResult, SetItemPipeline, TrackedRenderPass, ViewBinnedRenderPhases,
        },
        render_resource::{
            BufferUsages, ColorTargetState, ColorWrites, CompareFunction, DepthStencilState,
            FragmentState, IndexFormat, MultisampleState, PipelineCache, PrimitiveState,
            RawBufferVec, RenderPipelineDescriptor, SpecializedRenderPipeline,
            SpecializedRenderPipelines, TextureFormat, VertexAttribute, VertexBufferLayout,
            VertexFormat, VertexState, VertexStepMode,
        },
        renderer::{RenderDevice, RenderQueue},
        view::{self, ExtractedView, RenderVisibleEntities, VisibilitySystems},
        Render, RenderApp, RenderSet,
    }
};
use bevy_ecs::system::lifetimeless::Read;
use bevy_render::{mesh::{Indices, PrimitiveTopology}, render_resource::{binding_types::uniform_buffer, BindGroup, BindGroupLayout, BindGroupLayoutEntry, DynamicBindGroupEntries, DynamicBindGroupLayoutEntries, ShaderStages}, view::{ViewUniform, ViewUniforms}};
use bytemuck::{Pod, Zeroable};

use crate::Vertex;

/// A marker component that represents an entity that is to be rendered using
/// our custom phase item.
///
/// Note the [`ExtractComponent`] trait implementation. This is necessary to
/// tell Bevy that this object should be pulled into the render world.
#[derive(Clone, Component, ExtractComponent)]
struct CustomRenderedEntity;

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
struct DrawCustomPhaseItem;


/// FIXME: prepare the view uniforms before passing them to the pass ?
/// https://docs.rs/bevy_pbr/0.15.1/src/bevy_pbr/render/mesh_view_bindings.rs.html#715
/// I don't know why SRes does not work with ViewUniforms
impl<P> RenderCommand<P> for DrawCustomPhaseItem
where
    P: PhaseItem,
{
    type Param = SRes<CustomPhaseItemBuffers>;

    type ViewQuery = Read<ViewBindGroup>;

    type ItemQuery = ();

    fn render<'w>(
        _: &P,
        views: ROQueryItem<'w, Self::ViewQuery>,
        _: Option<ROQueryItem<'w, Self::ItemQuery>>,
        custom_phase_item_buffers: SystemParamItem<'w, '_, Self::Param>,
        pass: &mut TrackedRenderPass<'w>,
    ) -> RenderCommandResult {
        // Borrow check workaround.
        let custom_phase_item_buffers = custom_phase_item_buffers.into_inner();

        //let Some(mesh) = meshes else {return RenderCommandResult::Skip};


        println!("render the actual buffer");

        // Tell the GPU where the vertices are.
        pass.set_vertex_buffer(
            0,
            custom_phase_item_buffers
                .vertices
                .buffer()
                .unwrap()
                .slice(..),
        );

        // Tell the GPU where the indices are.
        pass.set_index_buffer(
            custom_phase_item_buffers
                .indices
                .buffer()
                .unwrap()
                .slice(..),
            0,
            IndexFormat::Uint32,
        );

        pass.set_bind_group(0, &views.0, &[0]);

        // Draw one triangle (3 vertices).
        pass.draw_indexed(0..6, 0, 0..1);

        RenderCommandResult::Success
    }
}

/// The GPU vertex and index buffers for our custom phase item.
///
/// As the custom phase item is a single triangle, these are uploaded once and
/// then left alone.
#[derive(Resource)]
struct CustomPhaseItemBuffers {
    /// The vertices for the single triangle.
    ///
    /// This is a [`RawBufferVec`] because that's the simplest and fastest type
    /// of GPU buffer, and [`Vertex`] objects are simple.
    vertices: RawBufferVec<Vertex>,

    /// The indices of the single triangle.
    ///
    /// As above, this is a [`RawBufferVec`] because `u32` values have trivial
    /// size and alignment.
    indices: RawBufferVec<u32>,
}

//impl Default for CustomPhaseItemBuffers {
//    fn default() -> Self {
//        CustomPhaseItemBuffers {
//            vertices: RawBufferVec::new(BufferUsages::VERTEX),
//            indices: RawBufferVec::new(BufferUsages::INDEX),
//        }
//    }
//}

/// The custom draw commands that Bevy executes for each entity we enqueue into
/// the render phase.
type DrawCustomPhaseItemCommands = (SetItemPipeline, DrawCustomPhaseItem);

/// A query filter that tells [`view::check_visibility`] about our custom
/// rendered entity.
type WithCustomRenderedEntity = With<CustomRenderedEntity>;

/// A single sraure's worth of vertices, for demonstration purposes.

/// The entry point.
pub fn test_main() {
    let mut app = App::new();
    app.add_plugins(DefaultPlugins)
        .add_plugins(CustomPhasePipelinePlugin)
        .add_systems(Startup, setup);

    app.run();
}

fn create_dummy_mesh() -> Mesh {
    Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::default())
    .with_inserted_attribute(
        Mesh::ATTRIBUTE_POSITION,
        vec![
            // top (facing towards +y)
            [-0.5, 0.5, -0.5], // vertex with index 0
            [0.5, 0.5, -0.5], // vertex with index 1
            [0.5, 0.5, 0.5], // etc. until 23
            [-0.5, 0.5, 0.5],
        ],
    )
    .with_inserted_indices(Indices::U32(vec![
        0,1,2 , 1,0,2, // triangles making up the top (+y) facing side.
    ]))
}

static VERTICES: [Vertex; 4] = [
    Vertex::new(vec3(-0.866, -0.5, 0.5), vec3(1.0, 0.0, 0.0)),
    Vertex::new(vec3(0.866, -0.5, 0.5), vec3(0.0, 1.0, 0.0)),
    Vertex::new(vec3(0.0, 1.0, 0.5), vec3(0.0, 0.0, 1.0)),
    Vertex::new(vec3(0.0, -1.0, 0.5), vec3(1.0, 1.0, 1.0)),
];

/// Spawns the objects in the scene.
fn setup(
    mut commands: Commands, 
    render_device: Res<RenderDevice>,
    render_queue: Res<RenderQueue>,
    mut meshes: ResMut<Assets<Mesh>>,
    ) {
    println!("setup scene");
    // Spawn a single entity that has custom rendering. It'll be extracted into
    // the render world via [`ExtractComponent`].
    commands.spawn((
        Visibility::default(),
        Transform::default(),
        // This `Aabb` is necessary for the visibility checks to work.
        Aabb {
            center: Vec3A::ZERO,
            half_extents: Vec3A::splat(0.5),
        },
        CustomRenderedEntity,
    ));
    //commands.spawn((
    //    Visibility::default(),
    //    Transform::default(),
    //    Mesh3d(meshes.add(create_dummy_mesh())),
    //    CustomRenderedEntity,
    //));

    // Spawn the camera.
    commands.spawn((
        Camera3d::default(),
        bevy::core_pipeline::tonemapping::Tonemapping::None,
        Transform::from_xyz(0.0, 0.0, 10.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    //commands.insert_resource(CustomPhaseItemBuffers::demo(&render_device, &render_queue));
}

pub struct CustomPhasePipelinePlugin;

impl Plugin for CustomPhasePipelinePlugin {
    fn build(&self, app: &mut App) {
        app
            .add_plugins(ExtractComponentPlugin::<CustomRenderedEntity>::default())
            .add_systems(
                PostUpdate,
                view::check_visibility::<WithCustomRenderedEntity>
                    .in_set(VisibilitySystems::CheckVisibility),
            );
        //render_app
        //    .init_resource::<SpecializedRenderPipelines<CustomMeshPipeline>>()
        //    .add_render_command::<Opaque3d, DrawCustomPhaseItemCommands>()
            //.add_systems(
            //    PostUpdate,
            //    // Make sure to tell Bevy to check our entity for visibility. Bevy won't
            //    // do this by default, for efficiency reasons.
            //    // This will do things like frustum culling and hierarchy visibility
            //    view::check_visibility::<WithCustomRenderedEntity>
            //    .in_set(VisibilitySystems::CheckVisibility),
            //)
            ;

    }
    fn finish(&self, app: &mut App) {
        // We make sure to add these to the render app, not the main app.
        let Some(render_app) = app.get_sub_app_mut(RenderApp) else {
            return;
        };

        render_app
            .init_resource::<CustomMeshPipeline>()
            .init_resource::<SpecializedRenderPipelines<CustomMeshPipeline>>()
            .add_render_command::<Opaque3d, DrawCustomPhaseItemCommands>()
            .add_systems(
                Render,
                (prepare_custom_phase_item_buffers.in_set(RenderSet::Prepare),
                prepare_view_bind_groups.in_set(RenderSet::Prepare),
                )
            )
            .add_systems(Render, queue_custom_phase_item.in_set(RenderSet::Queue));
        ;
    }
}

///// Creates the [`CustomPhaseItemBuffers`] resource.
/////
///// This must be done in a startup system because it needs the [`RenderDevice`]
///// and [`RenderQueue`] to exist, and they don't until [`App::run`] is called.
fn prepare_custom_phase_item_buffers(
    mut commands: Commands,
    //meshes: Query<&Mesh3d>,
    ) {
    commands.init_resource::<CustomPhaseItemBuffers>();
}

/// A render-world system that enqueues the entity with custom rendering into
/// the opaque render phases of each view.
fn queue_custom_phase_item(
    pipeline_cache: Res<PipelineCache>,
    custom_phase_pipeline: Res<CustomMeshPipeline>,
    mut opaque_render_phases: ResMut<ViewBinnedRenderPhases<Opaque3d>>,
    opaque_draw_functions: Res<DrawFunctions<Opaque3d>>,
    mut specialized_render_pipelines: ResMut<SpecializedRenderPipelines<CustomMeshPipeline>>,
    views: Query<(Entity, &RenderVisibleEntities, &Msaa), With<ExtractedView>>,
) {
    println!("queuing item");
    let draw_custom_phase_item = opaque_draw_functions
        .read()
        .id::<DrawCustomPhaseItemCommands>();

    // Render phases are per-view, so we need to iterate over all views so that
    // the entity appears in them. (In this example, we have only one view, but
    // it's good practice to loop over all views anyway.)
    for (view_entity, view_visible_entities, msaa) in views.iter() {
        let Some(opaque_phase) = opaque_render_phases.get_mut(&view_entity) else {
            continue;
        };

        // Find all the custom rendered entities that are visible from this
        // view.
        for &entity in view_visible_entities
            .get::<WithCustomRenderedEntity>()
            .iter()
        {
            println!("specialize pipeline for visible entity");
            // Ordinarily, the [`SpecializedRenderPipeline::Key`] would contain
            // some per-view settings, such as whether the view is HDR, but for
            // simplicity's sake we simply hard-code the view's characteristics,
            // with the exception of number of MSAA samples.
            let pipeline_id = specialized_render_pipelines.specialize(
                &pipeline_cache,
                &custom_phase_pipeline,
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
            opaque_phase.add(
                Opaque3dBinKey {
                    draw_function: draw_custom_phase_item,
                    pipeline: pipeline_id,
                    asset_id: AssetId::<Mesh>::invalid().untyped(),
                    material_bind_group_id: None,
                    lightmap_image: None,
                },
                entity,
                BinnedRenderPhaseType::NonMesh,
            );
        }
    }
}

impl SpecializedRenderPipeline for CustomMeshPipeline {
    type Key = Msaa;

    fn specialize(&self, msaa: Self::Key) -> RenderPipelineDescriptor {
        RenderPipelineDescriptor {
            label: Some("custom render pipeline".into()),
            layout: vec![
                self.view_layout.clone()
            ],
            push_constant_ranges: vec![],
            vertex: VertexState {
                shader: self.shader.clone(),
                shader_defs: vec![],
                entry_point: "vertex".into(),
                buffers: vec![VertexBufferLayout {
                    array_stride: size_of::<Vertex>() as u64,
                    step_mode: VertexStepMode::Vertex,
                    // This needs to match the layout of [`Vertex`].
                    attributes: vec![
                        VertexAttribute {
                            format: VertexFormat::Float32x3,
                            offset: 0,
                            shader_location: 0,
                        },
                        VertexAttribute {
                            format: VertexFormat::Float32x3,
                            offset: 16,
                            shader_location: 1,
                        },
                    ],
                }],
            },
            fragment: Some(FragmentState {
                shader: self.shader.clone(),
                shader_defs: vec![],
                entry_point: "fragment".into(),
                targets: vec![Some(ColorTargetState {
                    // Ordinarily, you'd want to check whether the view has the
                    // HDR format and substitute the appropriate texture format
                    // here, but we omit that for simplicity.
                    format: TextureFormat::bevy_default(),
                    blend: None,
                    write_mask: ColorWrites::ALL,
                })],
            }),
            primitive: PrimitiveState::default(),
            // Note that if your view has no depth buffer this will need to be
            // changed.
            depth_stencil: Some(DepthStencilState {
                format: CORE_3D_DEPTH_FORMAT,
                depth_write_enabled: false,
                depth_compare: CompareFunction::Always,
                stencil: default(),
                bias: default(),
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

//static VERTICES: [Vertex; 4] = [
//    Vertex::new(vec3(-0.866, -0.5, 0.5), vec3(1.0, 0.0, 0.0)),
//    Vertex::new(vec3(0.866, -0.5, 0.5), vec3(0.0, 1.0, 0.0)),
//    Vertex::new(vec3(0.0, 1.0, 0.5), vec3(0.0, 0.0, 1.0)),
//    Vertex::new(vec3(0.0, -1.0, 0.5), vec3(1.0, 1.0, 1.0)),
//];
//
fn extract_meshes_test(
    mut buffer: ResMut<CustomPhaseItemBuffers>,
    meshes_handles: Query<&Mesh3d>,
    meshes: Option<Res<Assets<Mesh>>>,
    render_device: Res<RenderDevice>,
    render_queue: Res<RenderQueue>,
){
    if meshes.is_some() {
        println!("I can access meshes");
    }
    if meshes_handles.iter().count() > 0 {
        println!("I have a mesh");
    }
}

fn extract_meshes(
    mut buffer: ResMut<CustomPhaseItemBuffers>,
    meshes_handles: Query<&Mesh3d>,
    meshes: Res<Assets<Mesh>>,
    render_device: Res<RenderDevice>,
    render_queue: Res<RenderQueue>,
){
    println!("extracting stuff ...");

    let mut vbo = RawBufferVec::new(BufferUsages::VERTEX);
    let mut ibo = RawBufferVec::new(BufferUsages::INDEX);
    for mesh in &meshes_handles {
        println!("extracting a mesh ...");
        let inner_mesh: &Mesh = match meshes.get(mesh) {
            Some(m) => m,
            _ => return
        };
        let count = inner_mesh.count_vertices();
        let new_indices = 
            inner_mesh.indices()
            .unwrap()
            .iter()
            .map(|x| (x+ count) as u32);
        ibo.extend(new_indices);
        // TODO: extract info

    }

    vbo.extend(VERTICES);

    vbo.write_buffer(&render_device, &render_queue);
    ibo.write_buffer(&render_device, &render_queue);

    *buffer = CustomPhaseItemBuffers {
        vertices: vbo,
        indices: ibo,
    }
}

impl CustomPhaseItemBuffers {
    fn demo(render_device: &RenderDevice, render_queue: &RenderQueue) -> Self {
        println!("creating custom buffer with data");
        let mut vbo = RawBufferVec::new(BufferUsages::VERTEX);
        let mut ibo = RawBufferVec::new(BufferUsages::INDEX);

        for vertex in &VERTICES {
            vbo.push(*vertex);
        }
        ibo.extend([0u32, 1, 2, 1, 0, 3]);

        // These two lines are required in order to trigger the upload to GPU.
        vbo.write_buffer(&render_device, &render_queue);
        ibo.write_buffer(&render_device, &render_queue);

        CustomPhaseItemBuffers {
            vertices: vbo,
            indices: ibo,
        }
    }
}

impl FromWorld for CustomPhaseItemBuffers {
    fn from_world(world: &mut World) -> Self {
        let render_device = world.resource::<RenderDevice>();
        let render_queue = world.resource::<RenderQueue>();

        // Create the vertex and index buffers.
        let mut vbo = RawBufferVec::new(BufferUsages::VERTEX);
        let mut ibo = RawBufferVec::new(BufferUsages::INDEX);

        for vertex in &VERTICES {
            vbo.push(*vertex);
        }
        ibo.extend([0u32, 1, 2, 1, 0, 3]);

        // These two lines are required in order to trigger the upload to GPU.
        vbo.write_buffer(render_device, render_queue);
        ibo.write_buffer(render_device, render_queue);

        CustomPhaseItemBuffers {
            vertices: vbo,
            indices: ibo,
        }
    }
}

#[derive(Component)]
struct ViewBindGroup( BindGroup );

// FIXME: init view uniforms
fn prepare_view_bind_groups(
    mut commands: Commands,
    render_device: Res<RenderDevice>,
    mesh_pipeline: Res<CustomMeshPipeline>,
    view_uniforms: Res<ViewUniforms>,
    views: Query< Entity>,
) {
    println!("give view data to shader");
    if let Some(view_binding) = view_uniforms.uniforms.binding() {
        for entity in &views {
        let entries = DynamicBindGroupEntries::new_with_indices(((0, view_binding.clone()),));
        let layout = &mesh_pipeline.view_layout;

        commands.entity(entity).insert(ViewBindGroup (
            render_device.create_bind_group("view_bind_group", layout, &entries),
        ));
        }
    }
}

impl FromWorld for CustomMeshPipeline {
    fn from_world(world: &mut World) -> Self {
        // Load and compile the shader in the background.
        let asset_server = world.resource::<AssetServer>();
        let render_device = world.resource::<RenderDevice>();

        let view_layout_entries: Vec<BindGroupLayoutEntry> = DynamicBindGroupLayoutEntries::new_with_indices(
            ShaderStages::FRAGMENT,
            (
                (0, uniform_buffer::<ViewUniform>(true).visibility(ShaderStages::VERTEX_FRAGMENT),
                    ),)
            ).to_vec()
        ;

        let view_layout = render_device.create_bind_group_layout(
            "mesh_view_layout",
            &view_layout_entries,
        );

        CustomMeshPipeline {
            shader: asset_server.load("custom_phase_item.wgsl"),
            view_layout,
        }
    }
}
