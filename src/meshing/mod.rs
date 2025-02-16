mod algorithms;
mod mesh_builder;
pub mod particles;
pub use mesh_builder::GeometryData;

use bevy::color::Color;
use bevy::math::{Mat3, Quat, Vec2, Vec3};
use bevy::prelude::{Component, Mesh};
use rand::prelude::Distribution;
use rand::Rng;

use crate::growing::{BranchSectionPosition, TreeSkeleton};
use crate::tools::{split_slice_circular, FloatProducer};
use crate::TreePipelinePhase;

use algorithms::{convex_hull_graham, extended_catmull_spline, mesh_between_contours, SplineIndex};
pub use particles::TrajectoryBuilder;

#[derive(Component)]
pub struct VolumetricTree {
    pub particles_per_node: Vec<Vec<usize>>,
    pub trajectories: Vec<Vec<Vec3>>,
    pub tree: TreeSkeleton,
}

#[derive(Copy, Clone, serde::Serialize, serde::Deserialize)]
pub struct StrandsConfig {
    pub particles_per_leaf: usize,
    pub repulsion: f32,
    pub wall_repulsion: f32,
    pub dt: f32,
    pub n_steps: usize,
    pub max_velocity_factor: f32,
    pub jump: usize,
}

impl TreePipelinePhase for VolumetricTree {
    type Previous = TreeSkeleton;
    type Config = StrandsConfig;
    type Builder = TrajectoryBuilder;
    fn generate_from(
        prev: Self::Previous,
        config: &Self::Config,
        builder: &mut Self::Builder,
    ) -> Self {
        builder.clear_for_tree(&prev);
        builder.compute_trajectories(&prev, 0, config);
        Self {
            trajectories: builder.trajectories.clone(),
            particles_per_node: builder.particles_per_node.clone(),
            tree: prev,
        }
    }
}

#[derive(Copy, Clone, serde::Serialize, serde::Deserialize)]
pub struct MeshConfig {
    leaf_size: f32,
    leaf_angle: f32,
    interior_angle: f32,
}

impl TreePipelinePhase for Mesh {
    type Previous = VolumetricTree;
    type Config = MeshConfig;
    type Builder = GeometryData;
    fn generate_from(
        prev: Self::Previous,
        config: &Self::Config,
        builder: &mut Self::Builder,
    ) -> Self {
        let pos_root = BranchSectionPosition::new(prev.tree.root(), 0.);
        let root_section = prev.register_branch_contour(pos_root, builder, config);
        prev.compute_each_branch_recursive(pos_root, root_section, builder, config);

        builder.to_mesh()
    }
}

impl VolumetricTree {
    fn particles_on_section(&self, pos: BranchSectionPosition) -> Vec<Vec3> {
        let depth = self.tree.depth(pos.node);
        if pos.length == 0. {
            let spline_index = SplineIndex::Local(depth, 0.);
            self.particles_per_node[pos.node]
                .iter()
                .map(|&particle| extended_catmull_spline(&self.trajectories[particle], spline_index))
                .collect()
        } else if pos.length < 0. {
            let parent = self.tree.parent(pos.node).unwrap();
            let branch_len = (self.tree.position(pos.node) - self.tree.position(parent)).length();
            let t = (branch_len + pos.length) / branch_len;
            let spline_index = SplineIndex::Local(depth-1, t);
            self.particles_per_node[pos.node]
                .iter()
                .map(|&particle| extended_catmull_spline(&self.trajectories[particle], spline_index))
                .collect()
        } else {
            match self.tree.children(pos.node) {
                &[child] => {
                    let t = pos.length / self.tree.branch_length_to_parent(child);
                    let spline_index = SplineIndex::Local(depth, t);
                    self.particles_per_node[pos.node]
                        .iter()
                        .map(|&particle| extended_catmull_spline(&self.trajectories[particle], spline_index))
                        .collect()
                },
                &[child1, child2] => {
                    let t1 = pos.length / self.tree.branch_length_to_parent(child1);
                    let t2 = pos.length / self.tree.branch_length_to_parent(child2);
                    let spline_index_1 = SplineIndex::Local(depth, t1);
                    let spline_index_2 = SplineIndex::Local(depth, t2);
                    let part_1 = self.particles_per_node[child1]
                        .iter()
                        .map(|&particle| extended_catmull_spline(&self.trajectories[particle], spline_index_1));
                    let part_2 = self.particles_per_node[child2]
                        .iter()
                        .map(|&particle| extended_catmull_spline(&self.trajectories[particle], spline_index_2));

                    part_1.chain(part_2).collect()
                },
                _ => panic!("leaf does not have outgoing branch")
            }
        }

    }

    fn branch_is_spliting(&self, pos: BranchSectionPosition) -> bool {
        if self.tree.children(pos.node).len() != 2 {
            return false;
        };
        let m_child = self.tree.children(pos.node)[0];
        let s_child = self.tree.children(pos.node)[1];

        let t = pos.length / self.tree.average_branch_length_to_children(pos.node);
        let m_center = self
            .tree
            .position(pos.node)
            .lerp(self.tree.position(m_child), t);
        let s_center = self
            .tree
            .position(pos.node)
            .lerp(self.tree.position(s_child), t);

        let pos_along_dir = |a: Vec3| a.dot(s_center - m_center);

        let m_pos = BranchSectionPosition::new(
            m_child,
            f32::min(0., pos.length - self.tree.branch_length_to_parent(m_child)),
        );
        let s_pos = BranchSectionPosition::new(
            s_child,
            f32::min(0., pos.length - self.tree.branch_length_to_parent(s_child)),
        );
        let m_part = self.particles_on_section(m_pos);
        let s_part = self.particles_on_section(s_pos);

        let m_relative_pos = m_part.into_iter().map(pos_along_dir);
        let s_relative_pos = s_part.into_iter().map(pos_along_dir);

        m_relative_pos.reduce(f32::max) < s_relative_pos.reduce(f32::min)
    }

    fn register_branch_contour(
        &self,
        pos: BranchSectionPosition,
        mesh: &mut GeometryData,
        config: &MeshConfig,
    ) -> Vec<usize> {
        let contour = self.branch_contour(pos, mesh, config);
        mesh.register_points_trunk(&contour)
    }

    fn branch_contour(&self, pos: BranchSectionPosition, mesh: &mut GeometryData, config: &MeshConfig) -> Vec<Vec3> {
        let parent = pos.node;
        let points = self.particles_on_section(pos);

        if points.len() == 0 {
            println!("node {} has no particles", pos.node);
        }

        let projected_points: Vec<Vec2> = points
            .iter()
            .map(|x: &Vec3| self.tree.space_to_plane(parent, *x))
            .collect();

        let center = self.tree.branch_section_center(pos);
        mesh.add_debug(center, Color::srgb(1.0, 1.0, 1.0));
        let center = self.tree.space_to_plane(parent, center);
        let result: Vec<Vec3> = convex_hull_graham(
            Some(center),
            &projected_points,
            Some(config.interior_angle),
        )
        .into_iter()
        // TODO: project point to avoid strands with too much diff ?
        .map(|i| points[i])
        .collect();
        mesh.add_contour(&result);

        result
    }

    fn compute_branch_join(
        &self,
        pos: BranchSectionPosition,
        previous_contour: Vec<usize>,
        mesh: &mut GeometryData,
        config: &MeshConfig,
    ) -> (
        (BranchSectionPosition, Vec<usize>),
        (BranchSectionPosition, Vec<usize>),
    ) {
        assert!(self.tree.children(pos.node).len() == 2);
        let normal = self.tree.normal(pos.node);

        let m_child = self.tree.children(pos.node)[0];
        let s_child = self.tree.children(pos.node)[1];

        let mut compute_properties = |child| {
            let branch_length = self.tree.branch_length_to_parent(child);
            let pos = BranchSectionPosition::new(child, f32::min(0., pos.length - branch_length));
            let center = self.tree.branch_section_center(pos);
            let contour = self.register_branch_contour(pos, mesh, config);
            (pos, center, contour)
        };

        let (m_p, m_c, m_cont) = compute_properties(m_child);
        let (s_p, s_c, s_cont) = compute_properties(s_child);

        assert!(
            m_cont.len() >= 3 && s_cont.len() >= 3,
            "not enough particles in the branch to compute join"
        );

        let m_dist_center = |i: &usize| (mesh.point(i) - m_c).length();
        let s_dist_center = |i: &usize| (mesh.point(i) - s_c).length();

        let i_m_furthest = m_cont.iter().map(s_dist_center).arg_min().unwrap();
        let i_s_furthest = s_cont.iter().map(m_dist_center).arg_min().unwrap();

        let center = 0.5 * (
              mesh.point(&m_cont[i_m_furthest])
            + mesh.point(&s_cont[i_s_furthest])
        );
        let dist_center = |i: &usize| (mesh.point(i) - center).length();

        let side =
            |i: &usize| Mat3::from_cols(m_c - s_c, normal, mesh.point(i) - center).determinant();


        let i_a = previous_contour
            .iter()
            .map(|i| side(i) / dist_center(i))
            .arg_min()
            .unwrap() as i32;
        let i_b = previous_contour
            .iter()
            .map(|i| side(i) / dist_center(i))
            .arg_max()
            .unwrap() as i32;

        let (m_junction, m_above) =
            split_slice_circular(&m_cont, i_m_furthest as i32 - 1, i_m_furthest as i32 + 1);
        let (mut s_junction, s_above) =
            split_slice_circular(&s_cont, i_s_furthest as i32 - 1, i_s_furthest as i32 + 1);
        s_junction.reverse();
        let (m_under, s_under) = split_slice_circular(&previous_contour, i_b, i_a);

        mesh.register_triangles(&mesh_between_contours(
            &mesh.points(),
            &m_under,
            &m_above,
            false,
        ));
        mesh.register_triangles(&mesh_between_contours(
            &mesh.points(),
            &s_under,
            &s_above,
            false,
        ));

        mesh.register_triangles(&mesh_between_contours(
            &mesh.points(),
            &s_junction,
            &m_junction,
            false,
        ));

        mesh.register_triangles(&[
            previous_contour[i_a as usize],
            s_junction[0],
            m_junction[0],
            previous_contour[i_b as usize],
            m_junction[2],
            s_junction[2],
        ]);

        ((m_p, m_cont), (s_p, s_cont))
    }

    fn compute_branch_while(
        &self,
        pos: BranchSectionPosition,
        previous_contour: &mut Vec<usize>,
        dz: f32,
        condition: impl Fn(BranchSectionPosition, &Self) -> bool,
        mesh: &mut GeometryData,
        config: &MeshConfig,
    ) -> BranchSectionPosition {
        let mut p = pos;
        assert!(dz != 0.);
        while condition(p, self) {
            if p.length > 1000. * dz {
                panic!("stopping, the branch at position {p:?} is already too long")
            }
            let current_contour = self.register_branch_contour(p, mesh, config);
            let triangles =
                mesh_between_contours(&mesh.points(), &previous_contour, &current_contour, true);
            mesh.register_triangles(&triangles);
            *previous_contour = current_contour;
            p += dz;
        }
        p
    }

    fn compute_each_branch_recursive(
        &self,
        pos: BranchSectionPosition,
        mut previous_contour: Vec<usize>,
        mesh: &mut GeometryData,
        config: &MeshConfig,
    ) {
        let root = pos.node;
        let radius = self.tree.radius(root);
        let dz = 0.2 * radius;

        match self.tree.children(root) {
            [] => {
                let leaf =
                    self.tree.position(root) + self.tree.radius(root) * self.tree.normal(root);
                let i_end = mesh.register_points_trunk(&[leaf])[0];
                let n = previous_contour.len();
                for i in 0..n {
                    mesh.register_triangles(&[
                        previous_contour[i],
                        previous_contour[(i + 1) % n],
                        i_end,
                    ]);
                }
                let random_point_distrib = RandomLeafPoint {
                    length: config.leaf_size,
                    max_angle: config.leaf_angle,
                    normal: self.tree.normal(root),
                };

                let p1 = leaf + mesh.rng.sample(random_point_distrib);
                let p2 = leaf + mesh.rng.sample(random_point_distrib);
                let i_leaf = mesh.register_points_leaf(&[leaf, p1, p2])[0];
                mesh.register_triangles(&[
                    i_leaf,
                    i_leaf + 1,
                    i_leaf + 2,
                    i_leaf,
                    i_leaf + 2,
                    i_leaf + 1,
                ]);
            }
            &[child] => {
                let branch_length = self.tree.branch_length_to_parent(child);
                let new_pos = BranchSectionPosition {
                    node: child,
                    length: pos.length - branch_length
                };
                let last_pos = self.compute_branch_while(
                    new_pos,
                    &mut previous_contour,
                    dz,
                    |p, _| p.length < 0.,
                    mesh,
                    config
                );
                self.compute_each_branch_recursive(last_pos, previous_contour, mesh, config)
            }
            &[m_child, s_child] => {
                let max_branch_length = f32::max(
                    self.tree.branch_length_to_parent(m_child),
                    self.tree.branch_length_to_parent(s_child),
                );
                let pos_split = self.compute_branch_while(
                    pos,
                    &mut previous_contour,
                    dz,
                    |p, me| !me.branch_is_spliting(p) && p.length < max_branch_length,
                    mesh,
                    config
                );

                let ((m_pos, mut m_cont), (s_pos, mut s_cont)) =
                    self.compute_branch_join(pos_split, previous_contour, mesh, config);

                let last_pos = self.compute_branch_while(m_pos+dz, &mut m_cont, dz, |p, _| p.length < 0., mesh, config);
                self.compute_each_branch_recursive(last_pos, m_cont, mesh, config);

                let last_pos = self.compute_branch_while(s_pos+dz, &mut s_cont, dz, |p, _| p.length < 0., mesh, config);
                self.compute_each_branch_recursive(last_pos, s_cont, mesh, config);
            }
            _ => panic!("did not expect more than 2 childs for node {root}"),
        }
    }
}

#[derive(Copy, Clone)]
struct RandomLeafPoint {
    max_angle: f32,
    length: f32,
    normal: Vec3,
}

impl Distribution<Vec3> for RandomLeafPoint {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        let a: f32 = rng.gen_range(0f32..self.max_angle);
        let c: f32 = rng.gen_range(0f32..2. * std::f32::consts::PI);
        let rot = Quat::from_rotation_z(c) * Quat::from_rotation_x(a) * Quat::from_rotation_z(-c);
        rot * (self.length * self.normal)
    }
}
