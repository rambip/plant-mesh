use std::ops::{Add, AddAssign};

use bevy::math::{Quat, Vec2, Vec3};
use bevy::prelude::{Mesh, Color};
use bevy::asset::RenderAssetUsages;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy_gizmos::prelude::Gizmos;

use crate::growing::{PlantNode, PlantNodeProps, NodeInfo};
mod meshing;

use meshing::{extended_catmull_spline, SplineIndex};

#[derive(Copy, Clone, Debug)]
struct BranchSectionPosition {
    // the node being considered
    node: usize,
    // the linear parameter.
    // If `t=0.`, we consider the branch section of the node.
    // If `t=1.`, we consider the branch sections of the direct children
    // If `t=-1.`, we consider the branch section of the parent.
    t: f32,
}

impl BranchSectionPosition {
    fn decimal_part(&self) -> f32 {
        self.t - self.t.floor()
    }
}

impl Add<f32> for BranchSectionPosition {
    type Output = BranchSectionPosition;

    fn add(self, rhs: f32) -> Self::Output {
        Self {
            node: self.node,
            t: self.t+rhs
        }
    }
}

impl AddAssign<f32> for BranchSectionPosition {
    fn add_assign(&mut self, rhs: f32) {
        self.t += rhs
    }
}

impl BranchSectionPosition {
    fn new(node: usize, t: f32) -> Self {
        Self {node, t}
    }
}

pub struct MeshBuilder {
    node_info: Vec<NodeInfo>,
    // from bottom to top
    trajectories: Vec<Vec<Vec3>>,
    particles_per_node: Vec<Vec<usize>>,
    node_props: Vec<PlantNodeProps>,
    mesh_points: Vec<Vec3>,
    debug_points: Vec<(Vec3, Color)>,
    mesh_triangles: Vec<usize>,
    mesh_normals: Vec<Vec3>,
    mesh_colors: Vec<[f32; 4]>,
    tree_depth: usize,
}

impl MeshBuilder {
    pub fn new(plant_graph: &PlantNode) -> Self {
        let mut node_props = Vec::new();
        plant_graph.register_node_properties(&mut node_props);
        let mut node_info = Vec::new();
        plant_graph.register_node_info(&mut node_info, 0);

        let node_count = node_info.len();
        let tree_depth = plant_graph.compute_depth();

        Self {
            node_props,
            node_info,
            tree_depth,
            trajectories: vec![],
            particles_per_node: vec![vec![]; node_count],
            mesh_points: vec![],
            mesh_triangles: vec![],
            mesh_normals: vec![],
            mesh_colors: vec![],
            debug_points: vec![],
        }
    }

    fn max_depth(&self) -> usize {
        self.tree_depth
    }

    fn depth(&self, node_id: usize) -> usize {
        self.node_info[node_id].depth
    }

    fn position(&self, node_id: usize) -> Vec3 {
        self.node_props[node_id].position
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
        let mut empty_trajectory = vec![Vec3::ZERO; self.depth(leaf_id)];
        empty_trajectory.push(position+relative_pos);

        self.trajectories.push(empty_trajectory);
    }

    fn register_particle_position_for_node(&mut self, particle_id: usize, position: Vec3, current_node: usize) {
        self.trajectories[particle_id][self.node_info[current_node].depth] = position;
        self.particles_per_node[current_node].push(particle_id);
    }

    fn register_points_on_contour(&mut self, points: &[Vec3], pos: BranchSectionPosition) -> Vec<usize> {
        let orientation = self.node_props[pos.node].orientation;
        let i0 = self.mesh_points.len();
        let n = points.len();
        self.mesh_points.extend(points);

        let r = (self.depth(pos.node) as f32 + pos.t) / self.max_depth() as f32;
        let color = [1. - r, 0.5+0.5*r, 0.2, 1.0];
        self.mesh_colors.extend(vec![color; n]);

        for i in 0..n {
            let v1 = points[(n+i-1)%n] - points[(n+i-2)%n];
            let v2 = points[(n+i+0)%n] - points[(n+i-1)%n];
            let v3 = points[(n+i+1)%n] - points[(n+i+0)%n];
            let v4 = points[(n+i+2)%n] - points[(n+i+1)%n];

            let curv1 = v1.normalize() - v2.normalize();
            let curv2 = v2.normalize() - v3.normalize();
            let curv3 = v3.normalize() - v4.normalize();

            let mut curv = curv1+2.*curv2+curv3;
            curv = curv - curv.dot(orientation)*orientation;

            self.mesh_normals.push(curv.normalize());
        }

        (i0..i0+n)
            .into_iter()
            .collect()
    }

    fn project_particles(&mut self, parent: usize, child: usize, offset: bool) {
        let origin = self.node_props[parent].position;
        let normal = self.node_props[parent].orientation;
        let r_parent = self.node_props[parent].radius;
        let p_child = self.node_props[child].position;

        let d = (origin - p_child).normalize();

        // FIXME: no clone
        for p in self.particles_per_node[child].clone() {
            let pos_particle = self.trajectories[p][self.depth(child)];
            let u = pos_particle - origin;

            let projected_offset = p_child - origin - (p_child - origin).dot(normal)*normal;

            let l = normal.dot(u) / normal.dot(d.normalize());
            assert!(!l.is_nan());
            let o = if offset {0.5*projected_offset} else {Vec3::ZERO};
            let projected = origin + r_parent*(o + (u - l * d.normalize())).normalize();

            self.register_particle_position_for_node(p, projected, parent);
        }
    }

    fn particles_on_section(&self, pos: BranchSectionPosition) -> Vec<Vec3> {
        let floor = pos.t.floor();
        let depth = (self.depth(pos.node) as i32) + floor as i32;
        let spline_index = SplineIndex::Local(depth as usize, pos.t - floor);
        //if (self.node_info[pos.node].children.len() == 0) {
        //    return vec![]
        //}
        println!("
            node is {}, at depth {}.
            It has an offset of {}.
            The calculated spline index is {:?}

        ", pos.node, self.depth(pos.node), pos.t, spline_index);

        self.particles_per_node[pos.node]
            .iter()
            .map(|&particle| extended_catmull_spline(&self.trajectories[particle], spline_index))
            .collect()
    }

    fn split(&mut self, pos: BranchSectionPosition) -> bool {
        if self.node_info[pos.node].children.len() != 2 {return false};
        let child1 = self.node_info[pos.node].children[0];
        let child2 = self.node_info[pos.node].children[1];

        let center_1 = self.position(pos.node).lerp(self.position(child1), pos.t);
        let center_2 = self.position(pos.node).lerp(self.position(child2), pos.t);

        let pos_along_dir = |a: Vec3| a.dot(center_2 - center_1);
        let relative_pos_1 = self.particles_on_section(BranchSectionPosition::new(child1, pos.t-1.)).into_iter().map(pos_along_dir);
        let relative_pos_2 = self.particles_on_section(BranchSectionPosition::new(child2, pos.t-1.)).into_iter().map(pos_along_dir);

        relative_pos_1.reduce(f32::max) < relative_pos_2.reduce(f32::min)
    }


    pub fn compute_trajectories(&mut self, root: usize, particle_per_leaf: usize) {
        assert!(self.particles_per_node[root].len()==0);

        match &self.node_info[root].children[..] {
            [] => {
                for _ in 0..particle_per_leaf {
                    self.register_particle_trajectory(root);
                }
            },
            &[child] => {
                self.compute_trajectories(child, particle_per_leaf);
                self.project_particles(root, child, false);
            },
            // TODO: assume child1 is the main branch
            &[child1, child2] => {
                let mut big_child = child1;
                let mut small_child = child2;
                self.compute_trajectories(child1, particle_per_leaf);
                self.compute_trajectories(child2, particle_per_leaf);
                if self.particles_per_node[big_child].len() < self.particles_per_node[small_child].len() {
                    std::mem::swap(&mut big_child, &mut small_child);
                }
                self.project_particles(root, small_child, false);
                self.project_particles(root, big_child, true);
            }
            _ => panic!("did not expect more than 2 childs for node {root}")
        }
    }

    fn branch_contour(&self, pos: BranchSectionPosition) -> Vec<Vec3> {
        let parent = pos.node;
        let points = self.particles_on_section(pos);

        if points.len() == 0 {
            println!("node {} has no particles", pos.node);
        }

        let projected_points: Vec<Vec2> = points
            .iter()
            // TODO: smarter projection
            .map(|&x| Quat::from_rotation_arc(self.node_props[parent].orientation, Vec3::Z) * x)
            .map(|x: Vec3| x.truncate())
            .collect();

        meshing::convex_hull_graham(&projected_points)
            .into_iter()
            .map(|i| points[i])
            .collect()
    }

    fn register_triangles(&mut self, triangles: &[usize]) {
        self.mesh_triangles.extend(triangles)
    }

    pub fn compute_each_branch(&mut self) {
        let pos_root = BranchSectionPosition::new(0, 0.);
        let root_section = self.register_points_on_contour(&self.branch_contour(pos_root), pos_root);
        self.compute_each_branch_recursive(0, root_section)
    }

    fn compute_branch_until(&mut self, 
        pos: BranchSectionPosition, 
        mut previous_contour: Vec<usize>,
        dt:f32,
        condition: impl Fn(BranchSectionPosition, &mut Self)->bool
        ) -> Result<Vec<usize>, BranchSectionPosition> {
        let mut p = pos;
        while p.decimal_part()+dt <= 1. {
            if condition(p, self) {return Err(p)};
            let current_contour = self.register_points_on_contour(&self.branch_contour(p), p);
            let triangles = meshing::mesh_between_contours(&self.mesh_points, &previous_contour, &current_contour); 
            self.register_triangles(&triangles);
            previous_contour = current_contour;
            p += dt;
        }
        Ok(previous_contour)
    }

    fn compute_each_branch_recursive(&mut self, root: usize, start_section: Vec<usize>) {
        let previous_contour = start_section;

        let pos_root = BranchSectionPosition::new(root, 0.);

        let radius = self.node_props[root].radius;
        let dz = 2.0*radius * std::f32::consts::PI / previous_contour.len() as f32;

        match &self.node_info[root].children[..] {
            [] => {},
            &[child] => {
                let branch_length = (self.position(root) - self.position(child)).length();
                let last_contour = self.compute_branch_until(pos_root, previous_contour.clone(), dz/branch_length, |_,_| false).unwrap();
                self.compute_each_branch_recursive(child, last_contour)
            },
            // TODO: assume child1 is the main branch
            &[child1, child2] => {
                let branch_length = 0.5*((self.position(root) - self.position(child1)).length()
                    + (self.position(root) - self.position(child2)).length());
                let pos_split = self.compute_branch_until(pos_root, previous_contour.clone(), dz/branch_length, |t, me| 
                    // FIXME: don't pass self as argument
                    me.split(t)
                    ).err().unwrap();

                for c in [child1, child2] {
                    let p = BranchSectionPosition::new(c, pos_split.t - 1.);
                    let current_contour = self.register_points_on_contour(&self.branch_contour(p), p);
                    let last_contour = self.compute_branch_until(p, current_contour, dz/branch_length, |_,_| false).unwrap();
                    self.compute_each_branch_recursive(c, last_contour);
                }
            },
            _ => panic!("did not expect more than 2 childs for node {root}")
        }
    }

    pub fn to_mesh(&self) -> Mesh {
        let mut result = Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::default())
            .with_inserted_attribute(
                Mesh::ATTRIBUTE_POSITION,
                self.mesh_points.clone(),
            )
            //.with_inserted_attribute(
            //    Mesh::ATTRIBUTE_NORMAL,
            //    self.mesh_normals.clone(),
            //)
            .with_inserted_attribute(
                Mesh::ATTRIBUTE_COLOR,
                self.mesh_colors.clone(),
            )
            .with_inserted_indices(Indices::U32(
                    self.mesh_triangles.iter().map(|x| *x as u32).collect()
            ));
            result.compute_smooth_normals();
            result
    }
    
    pub fn debug(&self, gizmos: &mut Gizmos, debug_flags: crate::DebugFlags) {
        if debug_flags.mesh {
            for i in 0..self.mesh_triangles.len()/3 {
                let (ia, ib, ic) = (self.mesh_triangles[3*i], self.mesh_triangles[3*i+1], self.mesh_triangles[3*i+2]);
                let (pa, pb, pc) = (self.mesh_points[ia], self.mesh_points[ib], self.mesh_points[ic]);
                let color = Color::srgb(0., 0.4, 0.);
                gizmos.line(pa, pb,  color);
                gizmos.line(pb, pc,  color);
                gizmos.line(pc, pa,  color);
            }
        }

        if debug_flags.normals {
            for i in 0..self.mesh_normals.len() {
                let p = self.mesh_points[i];
                let normal = self.mesh_normals[i];
                let color = Color::srgb(1., 0.0, 0.);
                gizmos.line(p, p+0.3*normal, color);
            }
        }
        if debug_flags.strands {
            for traj in &self.trajectories {
                for i in 1..100 {
                    let t1 = i as f32 / 100.;
                    let t2 = (i+1) as f32 / 100.;
                    let pos1 = meshing::extended_catmull_spline(traj, SplineIndex::Global(t1));
                    let pos2 = meshing::extended_catmull_spline(traj, SplineIndex::Global(t2));
                    let color = Color::srgb(1., 0.5, 0.5);
                    gizmos.line(pos1, pos2, color);
                }

            }
        }
        if debug_flags.other {
            for &(p, c) in &self.debug_points {
                gizmos.cross(p, 0.2, c);
            }
        }
    }
}


fn sample_uniform_disk(radius: f32) -> Vec2 {
    Vec2::from_angle(rand::random::<f32>()*std::f32::consts::TAU) * radius*(rand::random::<f32>()).sqrt()
}

fn sample_uniform_circle(radius: f32) -> Vec2 {
    Vec2::from_angle(rand::random::<f32>()*std::f32::consts::TAU) * radius
}

