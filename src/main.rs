//! A simple 3D scene with light shining over a cube sitting on a plane.

use bevy::prelude::*;
use std::f32::consts::PI;

#[derive(Component)]
struct Tree {
    plant_graph: growing::PlantNode,
    // including root
    node_count: usize,
    particle_per_leaf: usize,
    cache: Option<meshing::MeshBuilder>,
}

mod meshing;
mod growing;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
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
    mut materials: ResMut<Assets<StandardMaterial>>,
    camera_settings: Res<CameraSettings>,
) {
    // draw a floor
    commands.spawn((
        Mesh3d(meshes.add(Cuboid::new(100., 100., 0.1))),
        MeshMaterial3d(materials.add(Color::srgb_u8(200, 200, 200))),
        Transform::from_translation(-Vec3::Z * 0.1)
    ));
    commands.spawn((
        Mesh3d(meshes.add(Circle::new(4.0))),
        MeshMaterial3d(materials.add(Color::WHITE)),
    ));
    // light
    commands.spawn((
        PointLight {
            shadows_enabled: true,
            ..default()
        },
        Transform::from_xyz(5.0, 0.0, 20.0),
    ));
    // camera
    commands.spawn((
        Camera3d::default(),
        camera_settings.transform(),
    ));
    commands.spawn((
        Tree::default(),
        Mesh3d::default(),
        NeedRender(true),
        MeshMaterial3d(materials.add(Color::srgba(0.5, 1.0, 0.3, 0.8))),
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

}

#[derive(Component)]
struct NeedRender(bool);

fn draw_tree(
    mut gizmos: Gizmos,
    mut trees: Query<(&mut Mesh3d, &mut Tree, &mut NeedRender)>,
    mut meshes: ResMut<Assets<Mesh>>,
    ) {


    for (mut mesh, mut tree, mut need_render) in trees.iter_mut() {
        if need_render.0 {
            mesh.0 = meshes.add(tree.render_mesh());
            need_render.0 = false;
        }
        // only debug the tree after trying to render it
        tree.debug(&mut gizmos);
    }
}


impl Default for Tree {
    fn default() -> Self {
        let plant_graph = growing::PlantNode::demo();
        let node_count = plant_graph.count_node();
        Self {
            plant_graph,
            node_count,
            particle_per_leaf: 50,
            cache: None,
        }
    }
}

impl Tree {
    pub fn render_mesh(&mut self) -> Mesh {
        let mut cache = meshing::MeshBuilder::new(self);
        cache.compute_trajectories(self.particle_per_leaf);
        cache.compute_each_branch();

        // FIXME: duplicate points at nodes
        let mesh = cache.to_mesh();
        self.cache = Some(cache);

        mesh
    }


    // TODO: optimization = return an iterator

    pub fn debug(&self, gizmos: &mut Gizmos) {
        self.plant_graph.debug(gizmos);
        self.cache.as_ref().unwrap().debug(gizmos);

    }
}

// TODO: multiple materials for one single tree
