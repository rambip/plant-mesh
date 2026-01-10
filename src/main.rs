pub(crate) use bevy::prelude::*;
use bevy::{
    asset::{AssetMetaCheck, AsyncReadExt, LoadContext},
    input::{
        gestures::PinchGesture,
        mouse::{MouseMotion, MouseWheel},
    },
};
use bevy_simple_graphics::SimpleMeshPipelinePlugin;
use rand::{rngs::StdRng, SeedableRng};
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;

use bevy::asset::AssetLoader;
use bevy_gizmos::prelude::Gizmos;

use plant_mesh::{
    BevyMesh, GeometryData, Grow, MeshDebugFlags, PlantNode, Seed, TrajectoryBuilder, TreeConfig,
    TreeSkeleton, VisualDebug,
};
use plant_core::meshing::VolumetricTree;

#[derive(Copy, Clone, Default, Debug, Resource)]
pub struct DebugFlags {
    pub strands: bool,
    pub skeleton: bool,
    pub other: bool,
    pub mesh: MeshDebugFlags,
}

#[derive(Component)]
struct Tree {
    config: Handle<TreeConfig>,
    last_render_time: f32,
    need_render: bool,
    seed: u64,
}

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
    let mut app = App::new();
    app.add_plugins(
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
    .add_plugins(bevy_gizmos::GizmoPlugin)
    .add_plugins(SimpleMeshPipelinePlugin {
        shader_path: "shader.wgsl",
    })
    .insert_resource(ClearColor(Color::srgb(0.2, 0.25, 0.2)))
    .init_resource::<CameraSettings>()
    .init_resource::<DebugFlags>()
    .add_systems(Startup, setup)
    .add_systems(
        Update,
        (
            handle_input,
            update_view,
            automatic_mode,
            draw_tree,
            visual_debug,
        ),
    )
    .run();
}

/// set up a simple 3D scene
fn setup(mut commands: Commands, camera_settings: Res<CameraSettings>, server: Res<AssetServer>) {
    let config_handle: Handle<TreeConfig> = server.load("tree_config.toml");
    commands.spawn((
        Camera3d::default(),
        bevy::core_pipeline::tonemapping::Tonemapping::None,
        camera_settings.transform(0.),
    ));
    commands.spawn((
        Tree {
            config: config_handle,
            last_render_time: 0.,
            need_render: true,
            seed: 0,
        },
        Visibility::Visible,
    ));
}

#[derive(Resource, Debug)]
struct CameraSettings {
    orbit_distance: f32,
    orbit_angle: f32,
    sensibility: f32,
    tilt: f32,
    z: f32,
    automatic_mode: bool,
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
            automatic_mode: true,
            show_mesh: true,
        }
    }
}
impl CameraSettings {
    fn transform(&self, time: f32) -> Transform {
        let mut camera = Transform::default();
        let angle = if self.automatic_mode {
            0.5 * time
        } else {
            self.orbit_angle
        };
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
    mut trees: Query<&mut Tree>,
    mut flags: ResMut<DebugFlags>,
    time: Res<Time>,
) {
    if keyboard.just_pressed(KeyCode::KeyA) {
        camera_settings.automatic_mode = true;
    }
    if keyboard.pressed(KeyCode::ArrowRight) {
        camera_settings.automatic_mode = false;
        camera_settings.orbit_angle -= camera_settings.sensibility * time.delta_secs()
    }
    if keyboard.pressed(KeyCode::ArrowLeft) {
        camera_settings.automatic_mode = false;
        camera_settings.orbit_angle += camera_settings.sensibility * time.delta_secs()
    }
    if keyboard.just_pressed(KeyCode::ArrowUp) {
        camera_settings.automatic_mode = false;
        camera_settings.tilt += 0.1 * std::f32::consts::PI;
    }
    if keyboard.just_pressed(KeyCode::ArrowDown) {
        camera_settings.automatic_mode = false;
        camera_settings.tilt -= 0.1 * std::f32::consts::PI;
    }
    for ev in evr_motion.read() {
        if mouse.pressed(MouseButton::Left) {
            camera_settings.automatic_mode = false;
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
        for mut t in &mut trees {
            t.need_render = true;
            t.seed += 1;
        }
    }
    if keyboard.just_pressed(KeyCode::Enter) {
        camera_settings.show_mesh ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyW) {
        flags.mesh.triangles ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyS) {
        flags.strands ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyG) {
        flags.skeleton ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyC) {
        flags.mesh.contours ^= true;
    }
    if keyboard.just_pressed(KeyCode::KeyD) {
        flags.other ^= true;
    }
}

fn automatic_mode(
    time: Res<Time>,
    camera_settings: Res<CameraSettings>,
    mut trees: Query<&mut Tree>,
) {
    const UPDATE_TIME: f32 = 5.0;
    if camera_settings.automatic_mode {
        for mut t in &mut trees {
            if time.elapsed_secs() - t.last_render_time > UPDATE_TIME {
                t.seed += 1;
                t.need_render = true;
            }
        }
    }
}

fn draw_tree(
    mut commands: Commands,
    mut trees: Query<(Entity, &mut Tree)>,
    mut meshes: ResMut<Assets<Mesh>>,
    camera_settings: Res<CameraSettings>,
    configs: Res<Assets<TreeConfig>>,
    time: Res<Time>,
) {
    for (e, mut tree) in trees.iter_mut() {
        commands.entity(e).insert(if camera_settings.show_mesh {
            Visibility::Visible
        } else {
            Visibility::Hidden
        });

        if !tree.need_render {
            return;
        }
        tree.need_render = false;
        tree.last_render_time = time.elapsed_secs();

        let rng = StdRng::seed_from_u64(tree.seed);

        let mut plant_builder = rng.clone();
        let mut skeleton_builder = ();
        let mut particle_builder = TrajectoryBuilder::new(rng.clone());
        let mut mesh_builder = GeometryData::new(rng.clone());

        let Some(tree_config) = configs.get(&tree.config) else {
            return;
        };
        let BevyMesh(tree_mesh) = Seed
            .grow::<PlantNode>(&tree_config.grow, &mut plant_builder)
            .grow::<TreeSkeleton>(&(), &mut skeleton_builder)
            .grow::<VolumetricTree>(&tree_config.strands, &mut particle_builder)
            .grow::<BevyMesh>(&tree_config.mesh, &mut mesh_builder);

        let mesh = meshes.add(tree_mesh);
        commands.entity(e).insert(Mesh3d(mesh));

        commands.entity(e).insert((mesh_builder, particle_builder));
    }
}

fn visual_debug(
    query: Query<(
        Option<&TreeSkeleton>,
        Option<&TrajectoryBuilder>,
        Option<&GeometryData>,
    )>,
    flags: Res<DebugFlags>,
    mut gizmos: Gizmos,
) {
    for (a, b, c) in query.iter() {
        if let Some(a) = a {
            a.debug(&mut gizmos, flags.skeleton);
        }
        if let Some(b) = b {
            b.debug(&mut gizmos, flags.strands);
        }
        if let Some(c) = c {
            c.debug(&mut gizmos, flags.mesh);
        }
    }
}
