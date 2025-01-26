//! A simple 3D scene with light shining over a cube sitting on a plane.

use bevy::math::vec3;
pub(crate) use bevy::prelude::*;
use bevy_render::render_resource::{BufferUsages, RawBufferVec};
use bytemuck::{Pod, Zeroable};
use std::f32::consts::PI;

/// The CPU-side structure that describes a single vertex of the triangle.
#[derive(Clone, Copy, Pod, Zeroable)]
#[repr(C)]
pub struct Vertex {
    /// The 3D position of the triangle vertex.
    position: Vec3,
    /// Padding.
    pad0: u32,
    /// The color of the triangle vertex.
    color: Vec3,
    /// Padding.
    pad1: u32,
}

impl Vertex {
    /// Creates a new vertex structure.
    pub(crate) const fn new(position: Vec3, color: Vec3) -> Vertex {
        Vertex {
            position,
            color,
            pad0: 0,
            pad1: 0,
        }
    }
}

#[derive(Component)]
pub struct CustomMesh {
    /// The vertices for the single triangle.
    ///
    /// This is a [`RawBufferVec`] because that's the simplest and fastest type
    /// of GPU buffer, and [`Vertex`] objects are simple.
    vertices: Vec<Vertex>,

    /// The indices of the single triangle.
    ///
    /// As above, this is a [`RawBufferVec`] because `u32` values have trivial
    /// size and alignment.
    indices: Vec<u32>,
}

static VERTICES: [Vertex; 4] = [
    Vertex::new(vec3(-0.866, -0.5, 0.5), vec3(1.0, 0.0, 0.0)),
    Vertex::new(vec3(0.866, -0.5, 0.5), vec3(0.0, 1.0, 0.0)),
    Vertex::new(vec3(0.0, 1.0, 0.5), vec3(0.0, 0.0, 1.0)),
    Vertex::new(vec3(0.0, -1.0, 0.5), vec3(1.0, 1.0, 1.0)),
];

#[derive(Component)]
struct Tree {
    plant_graph: growing::PlantNode,
    // including root
    node_count: usize,
    particle_per_leaf: usize,
    cache: meshing::MeshBuilder,
}

mod meshing;
mod growing;

mod shader;
//use shader::{CustomRenderedMeshPipelinePlugin, CustomRenderedEntity};
use shader::custom_phase::test_main;

fn main(){
    test_main()
}

//fn main() {
//    App::new()
//        .add_plugins(DefaultPlugins)
//        //.add_plugins(bevy_pbr::MeshRenderPlugin {use_gpu_instance_buffer_builder: false})
//        //.add_plugins(bevy_pbr::PbrPlugin::default())
//        .add_plugins(CustomRenderedMeshPipelinePlugin)
//        .init_resource::<CameraSettings>()
//        .add_systems(Startup, setup)
//        .add_systems(Update, draw_tree)
//        .add_systems(Update, (handle_input, update_view))
//        .run();
//}

/// set up a simple 3D scene
fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    //mut materials: ResMut<Assets<StandardMaterial>>,
    camera_settings: Res<CameraSettings>,
) {
    // draw a floor
    commands.spawn((
        Mesh3d(meshes.add(Cuboid::new(100., 100., 0.1))),
        //MeshMaterial3d(materials.add(Color::srgb_u8(200, 200, 200))),
        Transform::from_translation(-Vec3::Z * 0.1)
    ));
    commands.spawn((
        Mesh3d(meshes.add(Circle::new(4.0))),
        //MeshMaterial3d(materials.add(Color::WHITE)),
    ));
    // light
    //commands.spawn((
    //    PointLight {
    //        shadows_enabled: true,
    //        ..default()
    //    },
    //    Transform::from_xyz(5.0, 0.0, 20.0),
    //));
    // camera
    commands.spawn((
        Camera3d::default(),
        bevy::core_pipeline::tonemapping::Tonemapping::None,
        camera_settings.transform(),
    ));
    commands.spawn((
        Tree::default(),
        Mesh3d::default(),
        NeedRender(true),
        //MeshMaterial3d(materials.add(Color::srgba(0.5, 1.0, 0.3, 0.8))),
        //CustomRenderedEntity,
    ));
}

#[derive(Resource, Debug)]
struct CameraSettings {
    pub orbit_distance: f32,
    pub orbit_angle: f32,
    pub sensibility: f32,
    z: f32,
}

impl Default for CameraSettings {
    fn default() -> Self {
        Self {
            orbit_distance: 20.,
            orbit_angle: 0.,
            sensibility: 1.,
            z: 5.,
        }
    }

}
impl CameraSettings {
    fn transform(&self) -> Transform {
        let mut camera = Transform::default();
        camera.rotation = Quat::from_rotation_z(self.orbit_angle)
                        * Quat::from_rotation_x(0.5*PI) // swap y and z
        ;

        // Adjust the translation to maintain the correct orientation toward the orbit target.
        // In our example it's a static target, but this could easily be customized.
        let target = Vec3::Z * self.z;
        camera.translation = target - camera.forward() * self.orbit_distance;
        camera
    }
}

fn update_view(
    camera_settings: Res<CameraSettings>,
    mut camera: Single<&mut Transform, With<Camera>>,
){
    **camera = camera_settings.transform();
}

fn handle_input(
    mut camera_settings: ResMut<CameraSettings>,
    keyboard: Res<ButtonInput<KeyCode>>,
    mut renders: Query<&mut NeedRender>,
    time: Res<Time>,
) {



    if keyboard.pressed(KeyCode::ArrowRight) {
        camera_settings.orbit_angle += camera_settings.sensibility * time.delta_secs()
    }
    if keyboard.pressed(KeyCode::ArrowLeft) {
        camera_settings.orbit_angle -= camera_settings.sensibility * time.delta_secs()
    }
    if keyboard.just_pressed(KeyCode::NumpadAdd) {
        camera_settings.orbit_distance -= 1.;
    }
    if keyboard.just_pressed(KeyCode::NumpadSubtract) {
        camera_settings.orbit_distance += 1.;
    }
    if keyboard.just_pressed(KeyCode::Space) {
        for mut r in &mut renders {
            r.0 = true
        }
    }

}

#[derive(Component)]
struct NeedRender(bool);

fn draw_tree(
    mut gizmos: Gizmos,
    mut trees: Query<(&mut Mesh3d, &mut Tree, &mut NeedRender)>,
    mut meshes: ResMut<Assets<Mesh>>,
    time: Res<Time>,
    ) {


    for (mut mesh, mut tree, mut need_render) in trees.iter_mut() {
        let t = time.elapsed_secs();
        let r: f32 = t/5. - (t / 5.).floor();
        tree.cache.set_triangle_proportion(r);

        mesh.0 = meshes.add(tree.render_mesh(need_render.0));
        need_render.0 = false;
        // only debug the tree after trying to render it
        //tree.debug(&mut gizmos);
    }
}


impl Default for Tree {
    fn default() -> Self {
        let plant_graph = growing::PlantNode::demo();
        let node_count = plant_graph.count_node();
        let cache = meshing::MeshBuilder::new(&plant_graph, node_count);
        Self {
            plant_graph,
            node_count,
            particle_per_leaf: 50,
            cache,
        }
    }
}

impl Tree {
    // TODO: don't use pbr but my own rendering pipeline.
    // https://github.com/bevyengine/bevy/blob/main/examples/shader/specialized_mesh_pipeline.rs
    pub fn render_mesh(&mut self, recompute: bool) -> Mesh {
        if recompute {
            self.cache = meshing::MeshBuilder::new(&self.plant_graph, self.node_count);
            self.cache.compute_trajectories(self.particle_per_leaf);
            self.cache.compute_each_branch();
        }

        self.cache.to_mesh()
    }


    // TODO: optimization = return an iterator

    pub fn debug(&self, gizmos: &mut Gizmos) {
        self.plant_graph.debug(gizmos);
        self.cache.debug(gizmos);

    }
}

// TODO: multiple materials for one single tree
