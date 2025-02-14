use bevy::{
    color::Color,
    math::{FloatPow, Vec2, Vec3},
    prelude::Component,
};
use bevy_gizmos::gizmos::Gizmos;
use rand::{prelude::Distribution, rngs::StdRng, Rng};

use crate::{growing::TreeSkeleton, VisualDebug};

use super::{
    algorithms::{extended_catmull_spline, SplineIndex},
    StrandsConfig,
};

struct UniformDisk {
    radius: f32,
}

impl UniformDisk {
    fn new(radius: f32) -> Self {
        Self { radius }
    }
}

impl Distribution<Vec2> for UniformDisk {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec2 {
        Vec2::from_angle(rng.gen_range(0f32..std::f32::consts::TAU))
            * self.radius
            * rng.gen_range(0f32..1f32).sqrt()
    }
}

fn spread_points(points: &mut Vec<Vec2>, radius: f32, config: &StrandsConfig) {
    let n = points.len();
    let mut velocities = vec![Vec2::ZERO; n];
    let wall_repulsion = config.wall_repulsion * radius;
    let repulsion = config.repulsion * radius.squared() / (n as f32).sqrt();
    let typical_distance = radius * (n as f32).sqrt();
    let max_force = config.max_force_factor * repulsion / typical_distance.squared() / config.dt;

    let max_radius = points.iter().map(|x| x.length()).reduce(f32::max).unwrap();

    let scale = 0.9 * radius / max_radius;

    for i in 0..n {
        points[i] *= scale;
    }

    for _ in 0..config.n_steps {
        for i in 0..n {
            velocities[i] = Vec2::ZERO;
            for j in 0..n {
                if i != j {
                    let relative_pos = points[i] - points[j];
                    let force = f32::min(max_force, repulsion / relative_pos.length_squared());
                    velocities[i] += force * relative_pos.normalize();
                }
            }
        }
        for i in 0..n {
            let target_position = points[i] + config.dt * velocities[i];
            if target_position.length() > radius {
                points[i] -= config.dt * wall_repulsion * target_position.normalize()
            } else {
                points[i] = target_position;
            }
        }
    }
}

#[derive(Component)]
pub struct TrajectoryBuilder {
    pub particles_per_node: Vec<Vec<usize>>,
    pub trajectories: Vec<Vec<Vec3>>,
    rng: StdRng,
}

impl TrajectoryBuilder {
    pub fn clear_for_tree(&mut self, tree: &TreeSkeleton) {
        let n = tree.node_count();
        self.particles_per_node = vec![vec![]; n];
        self.trajectories = vec![];
    }
    fn register_particles_for_leaf(&mut self, tree: &TreeSkeleton, node: usize, cloud: &[Vec2]) {
        for &point in cloud {
            let particle_id = self.trajectories.len();
            self.particles_per_node[node].push(particle_id);
            let mut empty_trajectory = vec![Vec3::ZERO; tree.depth(node)];
            empty_trajectory.push(tree.plane_to_space(node, point));

            self.trajectories.push(empty_trajectory);
        }
    }

    fn register_particles_for_node(
        &mut self,
        tree: &TreeSkeleton,
        node: usize,
        cloud: &[Vec2],
        particle_ids: &[usize],
    ) {
        assert_eq!(cloud.len(), particle_ids.len());

        let n = particle_ids.len();
        for i in 0..n {
            let d = tree.depth(node);
            self.trajectories[particle_ids[i]][d] = tree.plane_to_space(node, cloud[i]);
            self.particles_per_node[node].push(particle_ids[i]);
        }
    }

    fn project_particles(
        &mut self,
        tree: &TreeSkeleton,
        parent: usize,
        child: usize,
        offset: Vec3,
    ) -> Vec<Vec2> {
        let origin = tree.position(parent);
        let normal = tree.normal(parent);
        let orientation = tree.orientation(parent);
        let p_child = tree.position(child);

        let d = (origin - p_child).normalize();

        let projected = |particle_id: &usize| {
            let pos_particle = self.trajectories[*particle_id][tree.depth(child)];
            let u = pos_particle - origin;

            let l = normal.dot(u) / normal.dot(d.normalize());
            if l.is_nan() {
                return u.truncate();
            }
            //assert!(!l.is_nan());
            (orientation.inverse() * offset + (u - l * d.normalize())).truncate()
        };

        // TODO: no collect for opti
        self.particles_per_node[child]
            .iter()
            .map(projected)
            .collect()
    }

    pub fn compute_trajectories(
        &mut self,
        tree: &TreeSkeleton,
        root: usize,
        config: &StrandsConfig,
    ) {
        assert!(self.particles_per_node[root].len() == 0);
        let radius = tree.radius(root);

        match tree.children(root) {
            [] => {
                let radius = tree.radius(root);

                let mut cloud = self
                    .rng
                    .clone()
                    .sample_iter(UniformDisk::new(radius))
                    .take(config.particles_per_leaf)
                    .collect();

                spread_points(&mut cloud, radius, config);
                self.register_particles_for_leaf(tree, root, &cloud);
            }
            &[child] => {
                self.compute_trajectories(tree, child, config);
                let mut cloud = self.project_particles(tree, root, child, Vec3::ZERO);
                spread_points(&mut cloud, radius, config);
                let particle_ids = self.particles_per_node[child].clone();
                self.register_particles_for_node(tree, root, &cloud, &particle_ids);
            }
            &[m_child, s_child] => {
                self.compute_trajectories(tree, m_child, config);
                self.compute_trajectories(tree, s_child, config);
                let normal = tree.normal(root);
                let offset_direction = {
                    let u = tree.position(m_child) - tree.position(s_child);
                    (u - u.dot(normal) * normal).normalize()
                };
                let m_cloud = self.project_particles(
                    tree,
                    root,
                    m_child,
                    offset_direction * tree.radius(m_child),
                );
                let s_cloud = self.project_particles(
                    tree,
                    root,
                    s_child,
                    -offset_direction * tree.radius(s_child),
                );

                let mut cloud = m_cloud;
                cloud.extend(s_cloud);

                let mut particle_ids = self.particles_per_node[m_child].clone();
                particle_ids.extend(&self.particles_per_node[s_child]);

                spread_points(&mut cloud, radius, config);

                self.register_particles_for_node(tree, root, &cloud, &particle_ids);
            }
            _ => panic!("did not expect more than 2 childs for node {root}"),
        }
    }
}

impl From<StdRng> for TrajectoryBuilder {
    fn from(rng: StdRng) -> Self {
        Self {
            particles_per_node: vec![],
            trajectories: vec![],
            rng,
        }
    }
}

impl VisualDebug for TrajectoryBuilder {
    fn debug(&self, gizmos: &mut Gizmos, debug_flags: crate::DebugFlags) {
        if debug_flags.strands {
            let mut rng = self.rng.clone();
            for (i_t, traj) in self.trajectories.iter().enumerate() {
                let a: f32 = i_t as f32 / self.trajectories.len() as f32;
                let b: f32 = rng.gen();
                let color = Color::srgb(1., 0.3 + 0.5 * a, 0.3 + 0.5 * b);
                for i in 1..100 {
                    let t1 = i as f32 / 100.;
                    let t2 = (i + 1) as f32 / 100.;
                    let pos1 = extended_catmull_spline(traj, SplineIndex::Global(t1));
                    let pos2 = extended_catmull_spline(traj, SplineIndex::Global(t2));
                    gizmos.line(pos1, pos2, color);
                }
            }
        }
    }
}
