mod algorithms;
mod mesh_builder;
pub mod particles;
pub use mesh_builder::GeometryData;

use bevy::color::Color;
use bevy::math::{FloatExt, Mat3, Quat, Vec2, Vec3};
use bevy::prelude::{Component, Mesh};
use rand::prelude::Distribution;
use rand::Rng;

use crate::growing::TreeSkeleton;
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
    spacing: f32,
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
        let spline_index = SplineIndex::Local(0, 0.);
        let points_base = 
            prev.particles_per_node[0]
                .iter()
                .map(|&particle| extended_catmull_spline(&prev.trajectories[particle], spline_index));
        let contour = builder.register_points_trunk(points_base);
        prev.compute_each_branch_recursive(0, 0., contour, builder, config);

        builder.to_mesh()
    }
}

impl VolumetricTree {
    fn branch_is_spliting(&self, 
        parent: usize,
        m_child: usize,
        s_child: usize,
        l: f32,
        ) -> bool {
        let depth = self.tree.depth(parent);
        let m_t = l / self.tree.branch_length_to_parent(m_child);
        let m_spline_index = SplineIndex::Local(depth, m_t);

        let s_t = l / self.tree.branch_length_to_parent(s_child);
        let s_spline_index = SplineIndex::Local(depth, s_t);

        let radius = self.tree.radius(parent);

        let m_part: Vec<Vec3> = self.particles_per_node[m_child]
            .iter()
            .map(|&particle| extended_catmull_spline(&self.trajectories[particle], m_spline_index))
            .collect();

        let s_part: Vec<Vec3> = self.particles_per_node[s_child]
            .iter()
            .map(|&particle| extended_catmull_spline(&self.trajectories[particle], s_spline_index))
            .collect();

        let mut min_dist = f32::INFINITY;
        for &m in &m_part {
            for &s in &s_part {
                min_dist = f32::min(min_dist, (m-s).length());
            }
        }

        min_dist >= 0.1 * radius
    }

    fn compute_branch_join(
        &self,
        parent: usize,
        m_child: usize,
        s_child: usize,
        l: f32,
        previous_contour: Vec<usize>,
        mesh: &mut GeometryData,
        config: &MeshConfig,
    ) -> ( Vec<usize>, Vec<usize>) {
        let normal = self.tree.normal(parent);
        let depth = self.tree.depth(parent);

        let mut compute_properties = |child: usize| {
            let branch_len = self.tree.branch_length_to_parent(child);
            let t = l / branch_len;
            let center = self.tree.position(parent).lerp(
                self.tree.position(child), t);
            let spline_index = SplineIndex::Local(depth, t);
            let points: Vec<Vec3> = self.particles_per_node[child]
                .iter()
                .map(|&particle| extended_catmull_spline(&self.trajectories[particle], spline_index))
                .collect();
            let projected_points: Vec<Vec2> = points
                .iter()
                .map(|&x: &Vec3| self.tree.space_to_plane(parent, x))
                .collect();

            let result: Vec<Vec3> = convex_hull_graham(
                Some(self.tree.space_to_plane(child, center)),
                &projected_points,
                Some(config.interior_angle),
            )
            .into_iter()
            .map(|i| points[i])
            .collect();
            let contour = mesh.register_points_trunk(result);
            mesh.add_contour(&contour);
            (center, contour)
        };

        let (m_c, m_cont) = compute_properties(m_child);
        let (s_c, s_cont) = compute_properties(s_child);

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

        (m_cont, s_cont)
    }
    
    fn join_contours(&self, 
        parent: usize,
        previous_contour: &mut Vec<usize>,
        point_cloud: impl IntoIterator<Item=Vec3>,
        center: Vec3,
        mesh: &mut GeometryData,
        config: &MeshConfig,
        ) {
        let points: Vec<Vec3> = point_cloud.into_iter().collect();
        mesh.add_debug(center, Color::srgb(1.0, 1.0, 1.0));
        let center = self.tree.space_to_plane(parent, center);

        let projected_points: Vec<Vec2> = points
            .iter()
            .map(|&x: &Vec3| self.tree.space_to_plane(parent, x))
            .collect();

        let result: Vec<Vec3> = convex_hull_graham(
            Some(center),
            &projected_points,
            Some(config.interior_angle),
        )
        .into_iter()
        .map(|i| points[i])
        .collect();

        let current_contour = mesh.register_points_trunk(result);
        mesh.add_contour(&current_contour);
        let triangles = mesh_between_contours(
            &mesh.points(),
            &previous_contour,
            &current_contour, 
            true
        );
        mesh.register_triangles(&triangles);
        *previous_contour = current_contour;
    }

    fn compute_each_branch_recursive(
        &self,
        root: usize,
        mut l: f32,
        mut previous_contour: Vec<usize>,
        mesh: &mut GeometryData,
        config: &MeshConfig,
    ) {
        let depth = self.tree.depth(root);
        let parent_radius = self.tree.radius(root);

        match self.tree.children(root) {
            [] => {
                let leaf =
                    self.tree.position(root) + parent_radius * self.tree.normal(root);
                let i_end = mesh.register_points_trunk([leaf])[0];
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
                let i_leaf = mesh.register_points_leaf([leaf, p1, p2])[0];
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
                let child_radius = self.tree.radius(child);
                let branch_length = self.tree.branch_length_to_parent(child);
                while l < branch_length {
                    let t = l / branch_length;
                    let dz = (config.spacing*parent_radius).lerp(config.spacing*child_radius, t);
                    let spline_index = SplineIndex::Local(depth, t);
                    let points = self.particles_per_node[child] .iter() .map(|&particle| extended_catmull_spline(&self.trajectories[particle], spline_index));
                    let center = self.tree.position(root).lerp(
                        self.tree.position(child), t);
                    self.join_contours(root, &mut previous_contour, points, center, mesh, config);
                    l += dz;
                }
                self.compute_each_branch_recursive(child, l-branch_length, previous_contour, mesh, config)
            }
            &[m_child, s_child] => {
                let m_branch_length = self.tree.branch_length_to_parent(m_child);
                let s_branch_length = self.tree.branch_length_to_parent(s_child);
                let min_branch_length = f32::min(m_branch_length, s_branch_length);
                let dz = config.spacing*parent_radius;
                let m_radius = self.tree.radius(m_child);
                let s_radius = self.tree.radius(s_child);
                while l < min_branch_length && !self.branch_is_spliting(root, m_child, s_child, l) {
                    let m_t = l / m_branch_length;
                    let s_t = l / s_branch_length;
                    let m_spline_index = SplineIndex::Local(depth, m_t);
                    let s_spline_index = SplineIndex::Local(depth, s_t);
                    let m_points = self.particles_per_node[m_child] .iter() .map(|&particle| extended_catmull_spline(&self.trajectories[particle], m_spline_index));
                    let s_points = self.particles_per_node[s_child] .iter() .map(|&particle| extended_catmull_spline(&self.trajectories[particle], s_spline_index));
                    let points = m_points.chain(s_points);
                    let m_center = self.tree.position(root).lerp(
                        self.tree.position(m_child), m_t);
                    let s_center = self.tree.position(root).lerp(
                        self.tree.position(s_child), s_t);
                    let center = 0.5*(m_center + s_center);
                    self.join_contours(root, &mut previous_contour, points, center, mesh, config);
                    l += dz;
                }
                let (mut m_cont, mut s_cont) =
                    self.compute_branch_join(root, m_child, s_child, l, previous_contour, mesh, config);

                let mut m_l = l+dz;
                while m_l < m_branch_length {
                    let m_t = m_l / m_branch_length;
                    let m_dz = (config.spacing*parent_radius).lerp(config.spacing*m_radius, (m_l - l)/(m_branch_length - l));
                    let m_spline_index = SplineIndex::Local(depth, m_t);
                    let m_points = self.particles_per_node[m_child] .iter() .map(|&particle| extended_catmull_spline(&self.trajectories[particle], m_spline_index));
                    let m_center = self.tree.position(root).lerp(
                        self.tree.position(m_child), m_t);
                    self.join_contours(root, &mut m_cont, m_points, m_center, mesh, config);
                    m_l += m_dz;
                }
                self.compute_each_branch_recursive(m_child, m_l - m_branch_length, m_cont, mesh, config);

                let mut s_l = l+dz;
                while s_l < s_branch_length {
                    let s_t = s_l / s_branch_length;
                    let s_dz = (config.spacing*parent_radius).lerp(config.spacing*s_radius, (s_l - l)/(s_branch_length - l));
                    let s_spline_index = SplineIndex::Local(depth, s_t);
                    let s_points = self.particles_per_node[s_child] .iter() .map(|&particle| extended_catmull_spline(&self.trajectories[particle], s_spline_index));
                    let s_center = self.tree.position(root).lerp(
                        self.tree.position(s_child), s_t);
                    self.join_contours(root, &mut s_cont, s_points, s_center, mesh, config);
                    s_l += s_dz;
                }
                self.compute_each_branch_recursive(s_child, s_l - s_branch_length, s_cont, mesh, config)            }
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
