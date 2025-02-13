//! A simple 3D scene with light shining over a cube sitting on a plane.

use crate::meshing::VolumetricTree;
use bevy::input::{gestures::PinchGesture, mouse::{MouseMotion, MouseWheel}};
pub(crate) use bevy::prelude::*;
use growing::{PlantNode, TreeSkeleton};
use meshing::GeometryData;
use rand::{rngs::StdRng, Rng, SeedableRng};
use shader::CustomEntity;
use smallvec::SmallVec;
use std::f32::consts::PI;

use bevy_gizmos::prelude::Gizmos;


#[derive(Copy, Clone, Default, Debug, Resource)]
struct DebugFlags {
    triangles: bool,
    strands: bool,
    skeleton: bool,
    other: bool,
    contours: bool,
}

#[derive(Clone)]
pub struct NodeInfo {
    pub depth: usize,
    pub parent: Option<usize>,
    // id in prefix traversal order of tree
    pub id: usize,
    // TODO: bigger children first
    pub children: SmallVec<[usize; 2]>,
}

#[derive(Copy, Clone, Debug, Default)]
pub struct PlantNodeProps {
    pub position: Vec3,
    pub radius: f32,
    pub orientation: Quat,
}

trait VisualDebug {
    fn debug<R: Rng + Clone>(&self,
        gizmos: &mut Gizmos,
        rng: R,
        debug_flags: DebugFlags
        );
}


mod meshing;
mod growing;

#[derive(Component)]
struct TreeConfig {
    particle_per_leaf: usize,
}

impl TreeConfig {
    fn to_tree(&self) -> TreeSkeleton {
        PlantNode::basic_random().to_tree()
    }
}

mod shader;
mod tools;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                canvas: Some("#bevy".into()),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(bevy_sprite::SpritePlugin {})
        .add_plugins(bevy_gizmos::GizmoPlugin)
        .add_plugins(shader::CustomMeshPipelinePlugin)
        .init_resource::<CameraSettings>()
        .init_resource::<DebugFlags>()
        .add_systems(Startup, setup)
        .add_systems(Update, draw_tree)
        .add_systems(Update, (handle_input, update_view, visual_debug))
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
        Transform::from_translation(-Vec3::Z * 0.1),
    ));
    commands.spawn((Mesh3d(meshes.add(Circle::new(4.0))),));
    commands.spawn((
        Camera3d::default(),
        bevy::core_pipeline::tonemapping::Tonemapping::None,
        camera_settings.transform(0.),
    ));
    commands.spawn((
        TreeConfig::default(),
        Mesh3d::default(),
        NeedRender(true),
        shader::CustomEntity,
    ));
}

#[derive(Resource, Debug)]
struct CameraSettings {
    pub orbit_distance: f32,
    pub orbit_angle: f32,
    pub sensibility: f32,
    z: f32,
    animate: bool,
    show_mesh: bool,
}

impl Default for CameraSettings {
    fn default() -> Self {
        Self {
            orbit_distance: 20.,
            orbit_angle: 0.,
            sensibility: 1.,
            z: 5.,
            animate: true,
            show_mesh: true,
        }
    }
}
impl CameraSettings {
    fn transform(&self, time: f32) -> Transform {
        let mut camera = Transform::default();
        let angle = self.orbit_angle + if self.animate { 0.5 * time } else { 0. };
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
) {
    **camera = camera_settings.transform(time.elapsed_secs());
}

fn handle_input(
    mut camera_settings: ResMut<CameraSettings>,
    keyboard: Res<ButtonInput<KeyCode>>,
    mouse: Res<ButtonInput<MouseButton>>,
    mut evr_motion: EventReader<MouseMotion>,
    mut evr_scroll: EventReader<MouseWheel>,
    mut evr_gesture_pinch: EventReader<PinchGesture>,
    mut renders: Query<&mut NeedRender>,
    mut flags: ResMut<DebugFlags>,
    time: Res<Time>,
) {
    if keyboard.pressed(KeyCode::ArrowRight) {
        camera_settings.orbit_angle -= camera_settings.sensibility * time.delta_secs()
    }
    if keyboard.pressed(KeyCode::ArrowLeft) {
        camera_settings.orbit_angle += camera_settings.sensibility * time.delta_secs()
    }
    for ev in evr_motion.read() {
        if mouse.pressed(MouseButton::Left) {
            camera_settings.animate = false;
            camera_settings.z += 0.0005*ev.delta.y * camera_settings.orbit_distance;
            camera_settings.orbit_angle -= 0.0005*ev.delta.x * camera_settings.orbit_distance;
        }
    }
    for ev in evr_scroll.read() {
        #[cfg(target_family="wasm")]
        {
            camera_settings.orbit_distance -= 0.03*ev.y;
        }
        #[cfg(not(target_family="wasm"))]
        {
            camera_settings.orbit_distance -= ev.y;
        }
    }
    for ev in evr_gesture_pinch.read() {
        camera_settings.orbit_distance -= ev.0;
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
    if keyboard.just_pressed(KeyCode::Enter) {
        camera_settings.show_mesh ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyW) {
        flags.triangles ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyS) {
        flags.strands ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyG) {
        flags.skeleton ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyC) {
        flags.contours ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyD) {
        flags.other ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyA) {
        camera_settings.animate ^= true;
    }
}

#[derive(Component)]
struct NeedRender(bool);

fn draw_tree(
    mut commands: Commands,
    mut trees: Query<(Entity, &mut Mesh3d, &TreeConfig, &mut NeedRender)>,
    mut meshes: ResMut<Assets<Mesh>>,
    camera_settings: Res<CameraSettings>,
) {
    for (e, mut mesh, tree_config, mut need_render) in trees.iter_mut() {
        if need_render.0 {
            let seed = rand::random::<u64>();
            let rng = StdRng::seed_from_u64(seed);
            let tree = tree_config.to_tree();
            let mut mesh_builder = GeometryData::default();
            let strands = VolumetricTree::from_tree(
                tree.clone(), 
                rng.clone(),
                tree_config.particle_per_leaf
            );
            strands.compute_branches(&mut mesh_builder, rng.clone());

            mesh.0 = meshes.add(mesh_builder.to_mesh());
            need_render.0 = false;

            commands.entity(e).insert(tree);
            commands.entity(e).insert(strands);
            commands.entity(e).insert(mesh_builder);
        }


        if camera_settings.show_mesh {
            commands.entity(e).insert(CustomEntity);
        } else {
            commands.entity(e).remove::<CustomEntity>();
        }
    }
}

fn visual_debug(
    query: Query<
        (
            Option<&TreeSkeleton>,
            Option<&VolumetricTree>,
            Option<&GeometryData>,
        )>,
    flags: Res<DebugFlags>,
    mut gizmos: Gizmos,
) {
    let rng = StdRng::seed_from_u64(42);
    for (s, t, g) in &query {
        s.map(|x| x.debug(&mut gizmos, rng.clone(), *flags));
        t.map(|x| x.debug(&mut gizmos, rng.clone(), *flags));
        g.map(|x| x.debug(&mut gizmos, rng.clone(), *flags));
    }
}


impl Default for TreeConfig {
    fn default() -> Self {
        Self {
            particle_per_leaf: 30,
        }
    }
}
