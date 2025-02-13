use bevy::math::{FloatPow, Vec2, Vec3};
use rand::{prelude::Distribution, Rng};

use crate::growing::TreeSkeleton;

struct UniformDisk {
    radius: f32
}

impl UniformDisk {
    fn new(radius: f32) -> Self {
        Self {radius}
    }
}

impl Distribution<Vec2> for UniformDisk {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec2 {
        Vec2::from_angle(rng.gen_range(0f32..std::f32::consts::TAU))
            * self.radius
            * rng.gen_range(0f32..1f32).sqrt()
    }
}

struct ParticleSimulationConfig {
    repulsion: f32,
    dt: f32,
    n_steps: usize,
    // TODO: compute depending on number of particles
    particle_size: f32,
}

const DEFAULT_SIM_CONFIG: ParticleSimulationConfig = ParticleSimulationConfig {
    repulsion: 0.3,
    dt: 0.001,
    n_steps: 10,
    particle_size: 0.01,
};

fn spread_points(
    points: &mut Vec<Vec2>,
    radius: f32,
    config: ParticleSimulationConfig,
) {
    let n = points.len();
    let mut velocities = vec![Vec2::ZERO; n];
    let repulsion = config.repulsion * radius;

    let max_radius = points
        .iter()
        .map(|x| x.length())
        .reduce(f32::max)
        .unwrap();

    let scale = 0.9*radius / max_radius;

    for i in 0..n {
        points[i] *= scale;
    }

    for _ in 0..config.n_steps {
        for i in 0..n {
            velocities[i] = Vec2::ZERO;
            for j in 0..n {
                if i != j {
                    let relative_pos = points[i] - points[j];
                    let d = f32::max(
                        relative_pos.length_squared(),
                        config.particle_size.squared(),
                    );
                    velocities[i] += repulsion * relative_pos / d;
                }
            }
        }
        for i in 0..n {
            let target_position = points[i] + config.dt * velocities[i];
            if target_position.length() > radius {
                points[i] -= config.dt * velocities[i];
            } else {
                points[i] = target_position;
            }
        }
    }
}

pub struct TrajectoryBuilder<'a> {
    particles_per_node: Vec<Vec<usize>>,
    trajectories: Vec<Vec<Vec3>>,
    tree: &'a TreeSkeleton,
}

impl<'a> TrajectoryBuilder<'a> {
    pub fn new(tree: &'a TreeSkeleton) -> Self {
        let n = tree.node_count();
        Self {
            particles_per_node: vec![vec![]; n],
            trajectories: vec![],
            tree
        }
    }

    pub fn extract(self) -> (Vec<Vec<Vec3>>, Vec<Vec<usize>>) {
        (self.trajectories, self.particles_per_node)
    }

    fn register_particles_for_leaf(&mut self, node: usize, cloud: &[Vec2]) {
        for &point in cloud {
            let particle_id = self.trajectories.len();
            self.particles_per_node[node].push(particle_id);
            let mut empty_trajectory = vec![Vec3::ZERO; self.tree.depth(node)];
            empty_trajectory.push(self.tree.plane_to_space(node, point));

            self.trajectories.push(empty_trajectory);
        }
    }

    fn register_particles_for_node(&mut self, node: usize, cloud: &[Vec2], particle_ids: &[usize]) {
        assert_eq!(cloud.len(), particle_ids.len());

        let n = particle_ids.len();
        for i in 0..n {
            let d = self.tree.depth(node);
            self.trajectories[particle_ids[i]][d] = self.tree.plane_to_space(node, cloud[i]);
            self.particles_per_node[node].push(particle_ids[i]);
        }
    }

    fn project_particles(&mut self, parent: usize, child: usize, offset: Vec3) -> Vec<Vec2> {
        let origin = self.tree.position(parent);
        let normal = self.tree.normal(parent);
        let orientation = self.tree.orientation(parent);
        let p_child = self.tree.position(child);

        let d = (origin - p_child).normalize();

        let projected = |particle_id: &usize| {
            let pos_particle = self.trajectories[*particle_id][self.tree.depth(child)];
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

    pub fn compute_trajectories(&mut self, 
        root: usize, 
        rng: &mut impl Rng,
        particle_per_leaf: usize
        ) {
        assert!(self.particles_per_node[root].len() == 0);
        let radius = self.tree.radius(root);

        match self.tree.children(root) {
            [] => {
                let radius = self.tree.radius(root);

                let mut cloud = rng.sample_iter(UniformDisk::new(radius))
                    .take(particle_per_leaf)
                    .collect();

                spread_points(&mut cloud, radius, DEFAULT_SIM_CONFIG);
                self.register_particles_for_leaf(root, &cloud);

            }
            &[child] => {
                self.compute_trajectories(child, rng, particle_per_leaf);
                let mut cloud = self.project_particles(root, child, Vec3::ZERO);
                spread_points(&mut cloud, radius, DEFAULT_SIM_CONFIG);
                let particle_ids = self.particles_per_node[child].clone();
                self.register_particles_for_node(root, &cloud, &particle_ids);
            }
            &[m_child, s_child] => {
                self.compute_trajectories(m_child, rng, particle_per_leaf);
                self.compute_trajectories(s_child, rng, particle_per_leaf);
                let normal = self.tree.normal(root);
                let offset_direction = {
                    let u = self.tree.position(m_child) - self.tree.position(s_child);
                    (u - u.dot(normal) * normal).normalize()
                };
                let m_cloud = self.project_particles(
                    root,
                    m_child,
                    offset_direction * self.tree.radius(m_child)
                );
                let s_cloud = self.project_particles(
                    root,
                    s_child,
                    -offset_direction * self.tree.radius(s_child)
                );

                let mut cloud = m_cloud;
                cloud.extend(s_cloud);

                let mut particle_ids = self.particles_per_node[m_child].clone();
                particle_ids.extend(&self.particles_per_node[s_child]);

                spread_points(&mut cloud, radius, DEFAULT_SIM_CONFIG);

                self.register_particles_for_node(root, &cloud, &particle_ids);
            }
            _ => panic!("did not expect more than 2 childs for node {root}"),
        }
    }
}
