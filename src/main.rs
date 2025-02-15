//! A simple 3D scene with light shining over a cube sitting on a plane.

pub(crate) use bevy::prelude::*;
use bevy::{
    asset::{AsyncReadExt, LoadContext, AssetMetaCheck},
    input::{
        gestures::PinchGesture,
        mouse::{MouseMotion, MouseWheel},
    },
};
use rand::{rngs::StdRng, SeedableRng};
use serde::{Deserialize, Serialize};
use shader::CustomEntity;
use std::f32::consts::PI;

use bevy::asset::AssetLoader;
use bevy_gizmos::prelude::Gizmos;


mod shader;

use plant_mesh::{
    GrowConfig,
    StrandsConfig,
    MeshConfig,
    DebugFlags,
    Seed,
    GeometryData,
    TreeSkeleton,
    Grow,
    PlantNode,
    TreeSkeletonDebugData,
    VisualDebug,
    VolumetricTree,
    TrajectoryBuilder
};

#[derive(Component, Serialize, Deserialize, TypePath, Asset)]
struct TreeConfig {
    grow: GrowConfig,
    strands: StrandsConfig,
    mesh: MeshConfig,
}

// TODO: easiest way to do it ?
#[derive(Component)]
struct TreeConfigHandle(Handle<TreeConfig>);

struct TreeConfigLoader;

impl AssetLoader for TreeConfigLoader {
    type Asset = TreeConfig;

    type Settings = ();

    type Error = String;

    fn load(
        &self,
        reader: &mut dyn bevy::asset::io::Reader,
        _: &Self::Settings,
        _: &mut LoadContext,
    ) -> impl bevy::utils::ConditionalSendFuture<Output = Result<Self::Asset, Self::Error>> {
        async {
            // Read the content into a String
            let mut content = String::new();
            reader.read_to_string(&mut content).await.unwrap();

            // Parse TOML content into our Config struct
            let config: TreeConfig = toml::from_str(&content).unwrap();
            Ok(config)
        }
    }
    fn extensions(&self) -> &[&str] {
        &["toml"]
    }
}

fn main() {
    App::new()
        .add_plugins(
            DefaultPlugins
                .set(WindowPlugin {
                    primary_window: Some(Window {
                        canvas: Some("#bevy".into()),
                        ..default()
                    }),
                    ..default()
                })
                .set(AssetPlugin {
                    meta_check: AssetMetaCheck::Never,
                    ..default()
                }),
        )
        .init_asset::<TreeConfig>()
        .register_asset_loader(TreeConfigLoader)
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
    server: Res<AssetServer>,
) {
    let config_handle: Handle<TreeConfig> = server.load("tree_config.toml");
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
        TreeConfigHandle(config_handle),
        NeedRender(true),
        shader::CustomEntity,
    ));
}

#[derive(Resource, Debug)]
struct CameraSettings {
    orbit_distance: f32,
    orbit_angle: f32,
    sensibility: f32,
    tilt: f32,
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
            tilt: 0.,
            animate: true,
            show_mesh: true,
        }
    }
}
impl CameraSettings {
    fn transform(&self, time: f32) -> Transform {
        let mut camera = Transform::default();
        let angle = self.orbit_angle + if self.animate { 0.5 * time } else { 0. };
        camera.rotation =
              Quat::from_rotation_z(angle)
            * Quat::from_rotation_x(0.5*PI + self.tilt) // swap y and z
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
    if keyboard.just_pressed(KeyCode::ArrowUp) {
        camera_settings.tilt += 0.1 * std::f32::consts::PI;
    }
    if keyboard.just_pressed(KeyCode::ArrowDown) {
        camera_settings.tilt -= 0.1 * std::f32::consts::PI;
    }
    for ev in evr_motion.read() {
        if mouse.pressed(MouseButton::Left) {
            camera_settings.animate = false;
            camera_settings.z += 0.0005 * ev.delta.y * camera_settings.orbit_distance;
            camera_settings.orbit_angle -= 0.0005 * ev.delta.x * camera_settings.orbit_distance;
        }
    }
    for ev in evr_scroll.read() {
        #[cfg(target_family = "wasm")]
        {
            camera_settings.orbit_distance -= 0.03 * ev.y;
        }
        #[cfg(not(target_family = "wasm"))]
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
    mut trees: Query<(Entity, &TreeConfigHandle, &mut NeedRender)>,
    mut meshes: ResMut<Assets<Mesh>>,
    camera_settings: Res<CameraSettings>,
    configs: Res<Assets<TreeConfig>>,
    flags: Res<DebugFlags>,
) {
    for (e, tree_config_handle, mut need_render) in trees.iter_mut() {
        if camera_settings.show_mesh {
            commands.entity(e).insert(CustomEntity);
        } else {
            commands.entity(e).remove::<CustomEntity>();
        }

        if !need_render.0 {
            return;
        }
        need_render.0 = false;

        let rng = StdRng::seed_from_u64(rand::random::<u64>());

        let mut plant_builder = rng.clone().into();
        let mut skeleton_builder = rng.clone().into();
        let mut particle_builder = rng.clone().into();
        let mut mesh_builder = rng.clone().into();

        let Some(tree_config) = configs.get(&tree_config_handle.0) else {
            return;
        };
        let strands = Seed
            .grow::<PlantNode>(&tree_config.grow, &mut plant_builder)
            .grow::<TreeSkeleton>(&(), &mut skeleton_builder)
            .grow::<VolumetricTree>(&tree_config.strands, &mut particle_builder);

        commands
            .entity(e)
            .insert((skeleton_builder, particle_builder));

        if camera_settings.show_mesh || flags.need_render_mesh() {
            let tree_mesh = strands.grow::<Mesh>(&tree_config.mesh, &mut mesh_builder);

            let mesh = meshes.add(tree_mesh);
            commands.entity(e).insert(Mesh3d(mesh));

            commands.entity(e).insert(mesh_builder);
        }
    }
}

// TODO: more readable
fn visual_debug(
    query: Query<(
        Option<&TreeSkeletonDebugData>,
        Option<&TrajectoryBuilder>,
        Option<&GeometryData>,
    )>,
    flags: Res<DebugFlags>,
    mut gizmos: Gizmos,
) {
    for (a, b, c) in &query {
        a.debug(&mut gizmos, *flags);
        b.debug(&mut gizmos, *flags);
        c.debug(&mut gizmos, *flags);
    }
}
