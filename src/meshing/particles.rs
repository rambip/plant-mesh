use bevy::math::{Quat, Vec2, Vec3};

use crate::growing::{NodeInfo, PlantNodeProps};

fn sample_uniform_circle(radius: f32) -> Vec2 {
    Vec2::from_angle(rand::random::<f32>() * std::f32::consts::TAU) * radius
}

pub struct TrajectoryBuilder<'a> {
    particles_per_node: Vec<Vec<usize>>,
    trajectories: Vec<Vec<Vec3>>,
    node_props: &'a [PlantNodeProps],
    node_info: &'a [NodeInfo],
}

impl<'a> TrajectoryBuilder<'a> {
    pub fn new(node_props: &'a [PlantNodeProps], node_info: &'a [NodeInfo]) -> Self {
        let n = node_props.len();
        Self {
            particles_per_node: vec![vec![]; n],
            trajectories: vec![],
            node_props,
            node_info,
        }
    }

    pub fn extract(self) -> (Vec<Vec<Vec3>>, Vec<Vec<usize>>) {
        (self.trajectories, self.particles_per_node)
    }

    fn register_particle_trajectory(&mut self, leaf_id: usize) {
        let PlantNodeProps {
            position,
            radius,
            orientation,
        } = self.node_props[leaf_id];
        let particle_id = self.trajectories.len();
        let rotation = Quat::from_rotation_arc(Vec3::Z, orientation);
        let relative_pos = rotation * sample_uniform_circle(radius).extend(0.);
        self.particles_per_node[leaf_id].push(particle_id);
        let mut empty_trajectory = vec![Vec3::ZERO; self.node_info[leaf_id].depth];
        empty_trajectory.push(position + relative_pos);

        self.trajectories.push(empty_trajectory);
    }

    fn register_particle_position_for_node(
        &mut self,
        particle_id: usize,
        position: Vec3,
        current_node: usize,
    ) {
        let d = self.node_info[current_node].depth;
        self.trajectories[particle_id][d] = position;
        self.particles_per_node[current_node].push(particle_id);
    }

    fn project_particles(&mut self, parent: usize, child: usize, offset: Vec3) {
        let origin = self.node_props[parent].position;
        let normal = self.node_props[parent].orientation;
        let r_parent = self.node_props[parent].radius;
        let p_child = self.node_props[child].position;

        let d = (origin - p_child).normalize();

        // FIXME: no clone
        for p in self.particles_per_node[child].clone() {
            let pos_particle = self.trajectories[p][self.node_info[child].depth];
            let u = pos_particle - origin;

            let l = normal.dot(u) / normal.dot(d.normalize());
            assert!(!l.is_nan());
            let projected = origin + r_parent * (offset + (u - l * d.normalize())).normalize();

            self.register_particle_position_for_node(p, projected, parent);
        }
    }

    pub fn compute_trajectories(&mut self, root: usize, particle_per_leaf: usize) {
        assert!(self.particles_per_node[root].len() == 0);

        match &self.node_info[root].children[..] {
            [] => {
                for _ in 0..particle_per_leaf {
                    self.register_particle_trajectory(root);
                }
            }
            &[child] => {
                self.compute_trajectories(child, particle_per_leaf);
                self.project_particles(root, child, Vec3::ZERO);
            }
            &[m_child, s_child] => {
                self.compute_trajectories(m_child, particle_per_leaf);
                self.compute_trajectories(s_child, particle_per_leaf);
                let normal = self.node_props[root].orientation;
                let offset_direction = {
                    let u = self.node_props[m_child].position - self.node_props[s_child].position;
                    (u - u.dot(normal) * normal).normalize()
                };
                self.project_particles(
                    root,
                    m_child,
                    offset_direction * self.node_props[m_child].radius,
                );
                self.project_particles(
                    root,
                    s_child,
                    -offset_direction * self.node_props[s_child].radius,
                );
            }
            _ => panic!("did not expect more than 2 childs for node {root}"),
        }
    }
}
