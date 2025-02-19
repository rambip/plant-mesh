use rand::Rng;
use rand::{rngs::StdRng, SeedableRng};
use bevy::prelude::*;
use bevy_gizmos::gizmos::Gizmos;
use plant_mesh::meshing::{
    StrandsConfig,
    particles::{UniformDisk, spread_points}
};

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(bevy_sprite::SpritePlugin {})
        .add_plugins(bevy_gizmos::GizmoPlugin)
        .add_systems(Startup, setup)
        .add_systems(Update, (handle_input, simulate, draw))
        .run();
}

#[derive(Resource)]
struct StartRng(StdRng);

#[derive(Resource)]
struct SimulationSetup(StrandsConfig);

#[derive(Resource)]
struct NeedCompute(bool);

#[derive(Component)]
struct SimulationResult {
    points: Vec<Vec2>,
}

fn setup(
    mut commands: Commands
) {
    commands.spawn((
        Camera3d::default(),
        bevy::core_pipeline::tonemapping::Tonemapping::None,
        Transform::from_xyz(0., 0., 3.)
    ));
    commands.spawn(
        SimulationResult {
            points: vec![]
        }
    );

    commands.insert_resource(SimulationSetup (
        StrandsConfig {
            alpha: 1.0,
            contour_attraction: 0.,
            jump: 1,
            particles_per_leaf : 1000,
            wall_repulsion : 0.01,
            repulsion : 0.1,
            dt : 0.03,
            n_steps : 1,
            max_velocity_factor: 0.05,
            interaction_radius: 0.1,
        }
    ));
    commands.insert_resource(NeedCompute(true));
    commands.insert_resource(StartRng(StdRng::seed_from_u64(42)));
}

fn handle_input(
    mut simulation_params: ResMut<SimulationSetup>,
    keyboard: Res<ButtonInput<KeyCode>>,
    mut need_compute: ResMut<NeedCompute>,
) {
    if keyboard.just_pressed(KeyCode::ArrowRight) {
        simulation_params.0.n_steps += 1;
        println!("n steps: {}", simulation_params.0.n_steps);
        need_compute.0 = true;
    }
    if keyboard.just_pressed(KeyCode::ArrowLeft) {
        if simulation_params.0.n_steps > 0 {
            simulation_params.0.n_steps -= 1;
        }
        println!("n steps: {}", simulation_params.0.n_steps);
        need_compute.0 = true;
    }
}

fn simulate(
    config: Res<SimulationSetup>,
    rng: Res<StartRng>,
    mut result: Query<&mut SimulationResult>,
    mut need_compute: ResMut<NeedCompute>
) {
    for mut r in &mut result {
        if !need_compute.0 {
            return
        }
        need_compute.0 = false;

        let cloud1 = rng.0
            .clone()
            .sample_iter(UniformDisk::new(Vec2::ZERO, 0.5))
            .take(config.0.particles_per_leaf);

        let cloud2 = rng.0
            .clone()
            .sample_iter(UniformDisk::new(Vec2::X, 0.5))
            .take(config.0.particles_per_leaf);

        let mut cloud = cloud1.chain(cloud2).collect();

        spread_points(
            &mut cloud,
            1.0,
            &config.0,
        );

        r.points = cloud;
    }
}

fn draw(
    mut gizmos: Gizmos,
    result: Query<&SimulationResult>,
) {
    for r in &result {
        let color = Color::srgb(1.0, 0.0, 1.0);
        gizmos.circle(Isometry3d::default(), 1.0, color);
        let n = r.points.len()/2;
        for &p in &r.points[0..n] {
            let color = Color::srgb(0.0, 1.0, 0.0);
            gizmos.cross(p.extend(0.), 0.05, color);
        }
        for &p in &r.points[n..] {
            let color = Color::srgb(1.0, 0.0, 0.0);
            gizmos.cross(p.extend(0.), 0.05, color);
        }
    }
}
