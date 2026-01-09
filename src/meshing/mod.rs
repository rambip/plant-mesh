use bevy::prelude::Mesh;
use bevy::math::{FloatExt, Mat3, Quat, Vec2, Vec3};
use bevy::prelude::Component;
use rand::Rng;
use rand::prelude::Distribution;

use plant_core::growing::TreeSkeleton;
use plant_core::meshing::algorithms::{convex_hull_graham, extended_catmull_spline, mesh_between_contours, SplineIndex};
use plant_core::meshing::particles::TrajectoryBuilder;
use plant_core::meshing::mesh_builder::GeometryData;
use plant_core::TreePipelinePhase;
pub use plant_core::StrandsConfig;

use crate::tools::{split_slice_circular, FloatProducer};

pub mod algorithms {
    pub use plant_core::meshing::algorithms::*;
}
pub mod particles {
    pub use plant_core::meshing::particles::*;
}

#[derive(Component)]
pub struct VolumetricTree {
    pub particles_per_node: Vec<Vec<usize>>,
    pub trajectories: Vec<Vec<Vec3>>,
    pub tree: TreeSkeleton,
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
        builder.compute_trajectories(&prev, prev.root(), config);
        Self {
            trajectories: builder.trajectories.clone(),
            particles_per_node: builder.particles_per_node.clone(),
            tree: prev,
        }
    }
}

pub use plant_core::meshing::mesh_builder::MeshConfig;

pub struct BevyMesh(pub Mesh);

impl TreePipelinePhase for BevyMesh {
    type Previous = VolumetricTree;
    type Config = MeshConfig;
    type Builder = GeometryData;
    fn generate_from(
        prev: Self::Previous,
        config: &Self::Config,
        builder: &mut Self::Builder,
    ) -> Self {
        let spline_index = SplineIndex::Local(0, 0.);
        let points_base = prev.particles_per_node[0]
            .iter()
            .map(|&particle| extended_catmull_spline(&prev.trajectories[particle], spline_index));
        
        let mut rng = builder.rng.clone();
        let contour = builder.register_points_trunk(points_base, &mut rng);
        builder.rng = rng;
        
        builder.add_contour(&contour);
        prev.compute_each_branch_recursive(0, 0., contour, builder, config);
        
        // Convert GeometryData to Bevy Mesh
        let mut mesh = Mesh::new(
            bevy::render::mesh::PrimitiveTopology::TriangleList,
            bevy::asset::RenderAssetUsages::default(),
        )
        .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, builder.points.clone())
        .with_inserted_attribute(
            Mesh::ATTRIBUTE_COLOR,
            builder.colors.clone()
        )
        .with_inserted_indices(bevy::render::mesh::Indices::U32(
            builder.triangles.clone(),
        ));
        mesh.compute_smooth_normals();
        BevyMesh(mesh)
    }
}

impl VolumetricTree {
    fn branch_is_spliting(&self, parent: usize, m_child: usize, s_child: usize, l: f32) -> bool {
        let depth = self.tree.depth(parent);
        let m_t = l / (self.tree.position(m_child) - self.tree.position(parent)).length();
        let m_spline_index = SplineIndex::Local(depth, m_t);

        let s_t = l / (self.tree.position(s_child) - self.tree.position(parent)).length();
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
                min_dist = f32::min(min_dist, (m - s).length());
            }
        }

        min_dist >= 0.1 * radius
    }

    fn find_branch_split(
        &self,
        root: usize,
        m_child: usize,
        s_child: usize,
        dz: f32,
    ) -> Option<f32> {
        let m_len = (self.tree.position(m_child) - self.tree.position(root)).length();
        let s_len = (self.tree.position(s_child) - self.tree.position(root)).length();
        let min_branch_length = f32::min(m_len, s_len);
        
        if !self.branch_is_spliting(root, m_child, s_child, min_branch_length) {
            return None;
        }
        let mut interval = 0f32..min_branch_length;
        while (interval.end - interval.start) > dz {
            let middle = 0.5 * (interval.start + interval.end);
            if self.branch_is_spliting(root, m_child, s_child, middle) {
                interval = interval.start..middle
            } else {
                interval = middle..interval.end
            }
        }
        Some(interval.start)
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
    ) -> (Vec<usize>, Vec<usize>) {
        let normal = self.tree.normal(parent);
        let depth = self.tree.depth(parent);

        let mut compute_properties = |child: usize| {
            let branch_len = (self.tree.position(child) - self.tree.position(parent)).length();
            let t = l / branch_len;
            let center = self
                .tree
                .position(parent)
                .lerp(self.tree.position(child), t);
            let spline_index = SplineIndex::Local(depth, t);
            let points: Vec<Vec3> = self.particles_per_node[child]
                .iter()
                .map(|&particle| {
                    extended_catmull_spline(&self.trajectories[particle], spline_index)
                })
                .collect();
            let projected_points: Vec<Vec2> = points
                .iter()
                .map(|&x: &Vec3| self.tree.space_to_plane(parent, x))
                .collect();

            let result: Vec<Vec3> =
                convex_hull_graham(&projected_points, Some(config.interior_angle))
                    .into_iter()
                    .map(|i| points[i])
                    .collect();
            
            let mut rng = mesh.rng.clone();
            let contour = mesh.register_points_trunk(result, &mut rng);
            mesh.rng = rng;
            
            mesh.add_contour(&contour);
            (center, contour)
        };

        let (m_c, m_cont) = compute_properties(m_child);
        let (s_c, s_cont) = compute_properties(s_child);

        assert!(
            m_cont.len() >= 3 && s_cont.len() >= 3,
            "not enough particles in the branch to compute join"
        );

        let m_dist_center = |i: &usize| (mesh.point(*i) - m_c).length();
        let s_dist_center = |i: &usize| (mesh.point(*i) - s_c).length();

        let i_m_furthest = m_cont.iter().map(s_dist_center).arg_min().unwrap();
        let i_s_furthest = s_cont.iter().map(m_dist_center).arg_min().unwrap();

        let center = 0.5 * (mesh.point(m_cont[i_m_furthest]) + mesh.point(s_cont[i_s_furthest]));
        let dist_center = |i: &usize| (mesh.point(*i) - center).length();

        let side =
            |i: &usize| Mat3::from_cols(m_c - s_c, normal, mesh.point(*i) - center).determinant();

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
            &mesh.points,
            &m_under,
            &m_above,
            false,
        ));
        mesh.register_triangles(&mesh_between_contours(
            &mesh.points,
            &s_under,
            &s_above,
            false,
        ));

        mesh.register_triangles(&mesh_between_contours(
            &mesh.points,
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

    fn join_contours(
        &self,
        parent: usize,
        previous_contour: &mut Vec<usize>,
        point_cloud: impl IntoIterator<Item = Vec3>,
        center: Vec3,
        mesh: &mut GeometryData,
        config: &MeshConfig,
    ) {
        let points: Vec<Vec3> = point_cloud.into_iter().collect();
        mesh.debug_points.push((center, [1.0, 1.0, 1.0, 1.0]));

        let projected_points: Vec<Vec2> = points
            .iter()
            .map(|&x: &Vec3| self.tree.space_to_plane(parent, x))
            .collect();

        let result: Vec<Vec3> = convex_hull_graham(&projected_points, Some(config.interior_angle))
            .into_iter()
            .map(|i| points[i])
            .collect();

        let mut rng = mesh.rng.clone();
        let current_contour = mesh.register_points_trunk(result, &mut rng);
        mesh.rng = rng;
        
        mesh.add_contour(&current_contour);
        let triangles =
            mesh_between_contours(&mesh.points, &previous_contour, &current_contour, true);
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

        let children = self.tree.children(root);
        match children.len() {
            0 => {
                let leaf = self.tree.position(root) + parent_radius * self.tree.normal(root);
                let mut rng = mesh.rng.clone();
                let i_end = mesh.register_points_trunk([leaf], &mut rng)[0];
                mesh.rng = rng;
                
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
            1 => {
                let child = children[0];
                let branch_length = (self.tree.position(child) - self.tree.position(root)).length();
                let child_radius = self.tree.radius(child);
                while l < branch_length {
                    let t = l / branch_length;
                    let dz =
                        (config.spacing * parent_radius).lerp(config.spacing * child_radius, t);
                    let spline_index = SplineIndex::Local(depth, t);
                    let points = self.particles_per_node[child].iter().map(|&particle| {
                        extended_catmull_spline(&self.trajectories[particle], spline_index)
                    });
                    let center = self.tree.position(root).lerp(self.tree.position(child), t);
                    self.join_contours(root, &mut previous_contour, points, center, mesh, config);
                    l += dz;
                }
                self.compute_each_branch_recursive(
                    child,
                    l - branch_length,
                    previous_contour,
                    mesh,
                    config,
                )
            }
            2 => {
                let m_child = children[0];
                let s_child = children[1];
                let m_branch_length = (self.tree.position(m_child) - self.tree.position(root)).length();
                let s_branch_length = (self.tree.position(s_child) - self.tree.position(root)).length();
                let min_branch_length = f32::min(m_branch_length, s_branch_length);
                let dz = config.spacing * parent_radius;
                let m_radius = self.tree.radius(m_child);
                let s_radius = self.tree.radius(s_child);

                let l_split = self
                    .find_branch_split(root, m_child, s_child, dz)
                    .map(|x| x-3.*dz)
                    .unwrap_or(min_branch_length);

                while l < l_split {
                    let m_t = l / m_branch_length;
                    let s_t = l / s_branch_length;
                    let m_spline_index = SplineIndex::Local(depth, m_t);
                    let s_spline_index = SplineIndex::Local(depth, s_t);
                    let m_points = self.particles_per_node[m_child].iter().map(|&particle| {
                        extended_catmull_spline(&self.trajectories[particle], m_spline_index)
                    });
                    let s_points = self.particles_per_node[s_child].iter().map(|&particle| {
                        extended_catmull_spline(&self.trajectories[particle], s_spline_index)
                    });
                    let points = m_points.chain(s_points);
                    let m_center = self
                        .tree
                        .position(root)
                        .lerp(self.tree.position(m_child), m_t);
                    let s_center = self
                        .tree
                        .position(root)
                        .lerp(self.tree.position(s_child), s_t);
                    let center = 0.5 * (m_center + s_center);
                    self.join_contours(root, &mut previous_contour, points, center, mesh, config);
                    l += dz;
                }
                let (mut m_cont, mut s_cont) = self.compute_branch_join(
                    root,
                    m_child,
                    s_child,
                    l,
                    previous_contour,
                    mesh,
                    config,
                );

                let mut m_l = l + dz;
                while m_l < m_branch_length {
                    let m_t = m_l / m_branch_length;
                    let m_dz = (config.spacing * parent_radius)
                        .lerp(config.spacing * m_radius, (m_l - l) / (m_branch_length - l));
                    let m_spline_index = SplineIndex::Local(depth, m_t);
                    let m_points = self.particles_per_node[m_child].iter().map(|&particle| {
                        extended_catmull_spline(&self.trajectories[particle], m_spline_index)
                    });
                    let m_center = self
                        .tree
                        .position(root)
                        .lerp(self.tree.position(m_child), m_t);
                    self.join_contours(root, &mut m_cont, m_points, m_center, mesh, config);
                    m_l += m_dz;
                }
                self.compute_each_branch_recursive(
                    m_child,
                    m_l - m_branch_length,
                    m_cont,
                    mesh,
                    config,
                );

                let mut s_l = l + dz;
                while s_l < s_branch_length {
                    let s_t = s_l / s_branch_length;
                    let s_dz = (config.spacing * parent_radius)
                        .lerp(config.spacing * s_radius, (s_l - l) / (s_branch_length - l));
                    let s_spline_index = SplineIndex::Local(depth, s_t);
                    let s_points = self.particles_per_node[s_child].iter().map(|&particle| {
                        extended_catmull_spline(&self.trajectories[particle], s_spline_index)
                    });
                    let s_center = self
                        .tree
                        .position(root)
                        .lerp(self.tree.position(s_child), s_t);
                    self.join_contours(root, &mut s_cont, s_points, s_center, mesh, config);
                    s_l += s_dz;
                }
                self.compute_each_branch_recursive(
                    s_child,
                    s_l - s_branch_length,
                    s_cont,
                    mesh,
                    config,
                )
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
