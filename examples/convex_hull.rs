use bevy::prelude::*;
use bevy_gizmos::gizmos::Gizmos;
use plant_mesh::meshing::{algorithms::convex_hull_graham, particles::UniformDisk};
use rand::Rng;
use rand::{rngs::StdRng, SeedableRng};

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(bevy_gizmos::GizmoPlugin)
        .add_systems(Startup, setup)
        .add_systems(Update, (handle_input, simulate, draw))
        .run();
}

#[derive(Resource)]
struct SimulationSetup {
    n_points: usize,
    min_angle: f32,
    seed: u64,
}

#[derive(Resource)]
struct NeedCompute(bool);

#[derive(Component, Default)]
struct SimulationResult {
    points: Vec<Vec2>,
    border: Vec<usize>,
}

fn setup(mut commands: Commands) {
    commands.spawn((
        Camera3d::default(),
        bevy::core_pipeline::tonemapping::Tonemapping::None,
        Transform::from_xyz(0., 0., 3.),
    ));
    commands.spawn(SimulationResult::default());

    commands.insert_resource(SimulationSetup {
        n_points: 20000,
        min_angle: 2.0,
        seed: 0,
    });
    commands.insert_resource(NeedCompute(true));
}

fn handle_input(
    mut simulation_params: ResMut<SimulationSetup>,
    keyboard: Res<ButtonInput<KeyCode>>,
    mut need_compute: ResMut<NeedCompute>,
) {
    if keyboard.just_pressed(KeyCode::Space) {
        simulation_params.seed += 1;
        println!("seed: {}", simulation_params.seed);
        need_compute.0 = true;
    }
}

fn simulate(
    config: Res<SimulationSetup>,
    mut result: Query<&mut SimulationResult>,
    mut need_compute: ResMut<NeedCompute>,
) {
    for mut r in &mut result {
        if !need_compute.0 {
            return;
        }
        need_compute.0 = false;

        let rng = StdRng::seed_from_u64(config.seed);
        let cloud1 = rng
            .clone()
            .sample_iter(UniformDisk::new(Vec2::ZERO, 0.4))
            .take(config.n_points);

        let cloud2 = rng
            .clone()
            .sample_iter(UniformDisk::new(Vec2::X, 0.6))
            .take(config.n_points);

        let cloud: Vec<Vec2> = cloud1.chain(cloud2).collect();

        let result = convex_hull_graham(&cloud, Some(config.min_angle));

        r.points = cloud.to_vec();
        r.border = result.collect();
    }
}

fn draw(mut gizmos: Gizmos, result: Query<&SimulationResult>) {
    for r in &result {
        //for &p in &r.points {
        //    let color = Color::srgba(0.0, 0.5, 0.0, 0.3);
        //    gizmos.cross(p.extend(0.), 0.05, color);
        //}
        let n = r.border.len();
        for i in 0..n {
            let p1 = r.points[r.border[i]];
            let p2 = r.points[r.border[(i + 1) % n]];
            let color = Color::srgb(0.5, 0.5, 0.5);
            gizmos.line(p1.extend(0.), p2.extend(0.), color);
        }
    }
}
