mod algorithms;
mod mesh_builder;
mod particles;
pub use mesh_builder::GeometryData;

use bevy::color::Color;
use bevy::math::{Mat3, Vec2, Vec3};
use bevy::prelude::{Component, Mesh};

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
    pub dt: f32,
    pub n_steps: usize,
    // TODO: compute depending on number of particles
    pub particle_size: f32,
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

impl TreePipelinePhase for Mesh {
    type Previous = VolumetricTree;
    type Config = ();
    type Builder = GeometryData;
    fn generate_from(prev: Self::Previous, _: &Self::Config, builder: &mut Self::Builder) -> Self {
        let pos_root = BranchSectionPosition::new(prev.tree.root(), 0.);
        let root_section = prev.register_branch_contour(pos_root, builder);
        prev.compute_each_branch_recursive(prev.tree.root(), root_section, builder);

        builder.to_mesh()
    }
}

impl VolumetricTree {
    fn particles_on_section(&self, pos: BranchSectionPosition) -> Vec<Vec3> {
        let t = if pos.length < 0. {
            let parent = self.tree.parent(pos.node).unwrap();
            let branch_len = (self.tree.position(pos.node) - self.tree.position(parent)).length();
            (branch_len + pos.length) / branch_len
        } else {
            pos.length / self.tree.branch_length_to_main_children(pos.node)
        };
        let depth = if pos.length < 0. {
            self.tree.depth(pos.node) - 1
        } else {
            self.tree.depth(pos.node)
        };
        let spline_index = SplineIndex::Local(depth, t);

        self.particles_per_node[pos.node]
            .iter()
            .map(|&particle| extended_catmull_spline(&self.trajectories[particle], spline_index))
            .collect()
    }

    fn branch_is_spliting(&self, pos: BranchSectionPosition) -> bool {
        if self.tree.children(pos.node).len() != 2 {
            return false;
        };
        let m_child = self.tree.children(pos.node)[0];
        let s_child = self.tree.children(pos.node)[1];

        let t = pos.length / self.tree.branch_length_to_main_children(pos.node);
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
            pos.length - self.tree.branch_length_to_parent(m_child),
        );
        let s_pos = BranchSectionPosition::new(
            s_child,
            pos.length - self.tree.branch_length_to_parent(s_child),
        );
        let m_relative_pos = self
            .particles_on_section(m_pos)
            .into_iter()
            .map(pos_along_dir);
        let s_relative_pos = self
            .particles_on_section(s_pos)
            .into_iter()
            .map(pos_along_dir);

        m_relative_pos.reduce(f32::max) < s_relative_pos.reduce(f32::min)
    }

    fn register_branch_contour(
        &self,
        pos: BranchSectionPosition,
        mesh: &mut GeometryData,
    ) -> Vec<usize> {
        let contour = self.branch_contour(pos, mesh);
        mesh.register_points(&contour)
    }

    fn branch_contour(&self, pos: BranchSectionPosition, mesh: &mut GeometryData) -> Vec<Vec3> {
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
            Some(0.9 * std::f32::consts::PI),
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
            let pos = BranchSectionPosition::new(child, pos.length - branch_length);
            let center = self.tree.branch_section_center(pos);
            let contour = self.register_branch_contour(pos, mesh);
            (pos, center, contour)
        };

        let (m_p, m_c, m_cont) = compute_properties(m_child);
        let (s_p, s_c, s_cont) = compute_properties(s_child);

        let center = 0.5 * (m_c + s_c);

        let m_dist_center = |i: &usize| (mesh.point(i) - m_c).length();
        let s_dist_center = |i: &usize| (mesh.point(i) - s_c).length();
        let dist_center = |i: &usize| (mesh.point(i) - center).length();

        let side =
            |i: &usize| Mat3::from_cols(m_c - s_c, normal, mesh.point(i) - center).determinant();

        let i_m_furthest = m_cont.iter().map(s_dist_center).arg_min().unwrap() as i32;
        let i_s_furthest = s_cont.iter().map(m_dist_center).arg_min().unwrap() as i32;

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
            split_slice_circular(&m_cont, i_m_furthest - 1, i_m_furthest + 1);
        let (mut s_junction, s_above) =
            split_slice_circular(&s_cont, i_s_furthest - 1, i_s_furthest + 1);
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
    ) -> BranchSectionPosition {
        let mut p = pos;
        assert!(dz != 0.);
        p += dz;
        while condition(p, self) {
            if p.length > 1000. * dz {
                panic!("stopping, the branch at position {p:?} is already too long")
            }
            let current_contour = self.register_branch_contour(p, mesh);
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
        root: usize,
        mut previous_contour: Vec<usize>,
        mesh: &mut GeometryData,
    ) {
        let pos_root = BranchSectionPosition::new(root, 0.);

        let radius = self.tree.radius(root);
        let dz = 0.2 * radius;

        match self.tree.children(root) {
            [] => {
                let leaf =
                    self.tree.position(root) + self.tree.radius(root) * self.tree.normal(root);
                let i_end = mesh.register_points(&vec![leaf])[0];
                let n = previous_contour.len();
                for i in 0..n {
                    mesh.register_triangles(&[
                        previous_contour[i],
                        previous_contour[(i + 1) % n],
                        i_end,
                    ]);
                }
            }
            &[child] => {
                let branch_length = self.tree.branch_length_to_parent(child);
                self.compute_branch_while(
                    pos_root,
                    &mut previous_contour,
                    dz,
                    |p, _| p.length < branch_length,
                    mesh,
                );
                self.compute_each_branch_recursive(child, previous_contour, mesh)
            }
            &[m_child, s_child] => {
                let pos_split = self.compute_branch_while(
                    pos_root,
                    &mut previous_contour,
                    dz,
                    |p, me| !me.branch_is_spliting(p),
                    mesh,
                );

                let ((m_pos, mut m_cont), (s_pos, mut s_cont)) =
                    self.compute_branch_join(pos_split, previous_contour, mesh);

                self.compute_branch_while(m_pos, &mut m_cont, dz, |p, _| p.length < 0., mesh);
                self.compute_each_branch_recursive(m_child, m_cont, mesh);

                self.compute_branch_while(s_pos, &mut s_cont, dz, |p, _| p.length < 0., mesh);
                self.compute_each_branch_recursive(s_child, s_cont, mesh);
            }
            _ => panic!("did not expect more than 2 childs for node {root}"),
        }
    }
}
