//! A simple 3D scene with light shining over a cube sitting on a plane.

pub(crate) use bevy::prelude::*;
use std::f32::consts::PI;

use bevy_gizmos::prelude::Gizmos;

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

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(bevy_sprite::SpritePlugin {})
        .add_plugins(bevy_gizmos::GizmoPlugin)
        .add_plugins(shader::CustomMeshPipelinePlugin)
        .init_resource::<CameraSettings>()
        .add_systems(Startup, setup)
        .add_systems(Update, draw_tree)
        .add_systems(Update, (handle_input, update_view))
        .run();
}

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
        Transform::from_translation(-Vec3::Z * 0.1)
    ));
    commands.spawn((
        Mesh3d(meshes.add(Circle::new(4.0))),
    ));
    commands.spawn((
        Camera3d::default(),
        bevy::core_pipeline::tonemapping::Tonemapping::None,
        camera_settings.transform(0.),
    ));
    commands.spawn((
        Tree::default(),
        Mesh3d::default(),
        NeedRender(true),
        Visibility::default(),
        shader::CustomEntity,
    ));
}

#[derive(Resource, Debug)]
struct CameraSettings {
    pub orbit_distance: f32,
    pub orbit_angle: f32,
    pub sensibility: f32,
    z: f32,
    debug_shapes: bool,
    animate: bool,
}

impl Default for CameraSettings {
    fn default() -> Self {
        Self {
            orbit_distance: 20.,
            orbit_angle: 0.,
            sensibility: 1.,
            z: 5.,
            debug_shapes: false,
            animate: true,
        }
    }

}
impl CameraSettings {
    fn transform(&self, time: f32) -> Transform {
        let mut camera = Transform::default();
        let angle = self.orbit_angle + if self.animate {0.5*time} else {0.};
        camera.rotation = Quat::from_rotation_z(angle)
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
    time: Res<Time>,
){
    **camera = camera_settings.transform(time.elapsed_secs());
}

fn handle_input(
    mut camera_settings: ResMut<CameraSettings>,
    keyboard: Res<ButtonInput<KeyCode>>,
    mut renders: Query<&mut NeedRender>,
    time: Res<Time>,
) {



    if keyboard.pressed(KeyCode::ArrowRight) {
        camera_settings.orbit_angle -= camera_settings.sensibility * time.delta_secs()
    }
    if keyboard.pressed(KeyCode::ArrowLeft) {
        camera_settings.orbit_angle += camera_settings.sensibility * time.delta_secs()
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
    if keyboard.just_pressed(KeyCode::KeyD) {
        camera_settings.debug_shapes ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyA) {
        camera_settings.animate ^= true;
    }

}

#[derive(Component)]
struct NeedRender(bool);

fn draw_tree(
    mut gizmos: Gizmos,
    mut trees: Query<(&mut Mesh3d, &mut Tree, &mut NeedRender)>,
    mut meshes: ResMut<Assets<Mesh>>,
    camera_settings: Res<CameraSettings>,
    ) {


    for (mut mesh, mut tree, mut need_render) in trees.iter_mut() {
        mesh.0 = meshes.add(tree.render_mesh(need_render.0));
        need_render.0 = false;
        // only debug the tree after trying to render it
        if camera_settings.debug_shapes {
            tree.debug(&mut gizmos);
        }
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
