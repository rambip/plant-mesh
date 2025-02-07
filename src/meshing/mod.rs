use std::ops::{Add, AddAssign};
use std::sync::{Arc, Mutex};

use bevy::math::{Mat3, Quat, Vec2, Vec3};
use bevy::prelude::{Mesh, Color};
use bevy::asset::RenderAssetUsages;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy_gizmos::prelude::Gizmos;
use rand::prelude::*;

use crate::growing::{PlantNode, PlantNodeProps, NodeInfo};
use crate::tools::FloatProducer;
mod meshing;

use meshing::{extended_catmull_spline, SplineIndex};

#[derive(Copy, Clone, Debug)]
// TODO: enum ?
struct BranchSectionPosition {
    // the node being considered
    node: usize,
    // the distance we traveled from the node to the leaves
    // it can be negative, in this case we consider the parent
    length: f32,
}

impl Add<f32> for BranchSectionPosition {
    type Output = BranchSectionPosition;

    fn add(self, rhs: f32) -> Self::Output {
        Self {
            node: self.node,
            length: self.length+rhs
        }
    }
}

impl AddAssign<f32> for BranchSectionPosition {
    fn add_assign(&mut self, rhs: f32) {
        self.length += rhs
    }
}

impl BranchSectionPosition {
    fn new(node: usize, length: f32) -> Self {
        Self {node, length}
    }
}

fn split_contour(source: &[usize], start: i32, end: i32) -> (Vec<usize>, Vec<usize>) {
    let mut direct = Vec::new();
    let mut indirect = Vec::new();
    let n = source.len() as i32;

    let id = |x| ((x+n)%n) as usize;

    let mut i = start;
    while id(i) != id(end) {
        direct.push(source[id(i)]);
        i += 1;
    }
    direct.push(source[id(end)]);
    while id(i) != id(start) {
        indirect.push(source[id(i)]);
        i += 1;
    }
    indirect.push(source[id(start)]);
    (direct, indirect)
}


pub struct MeshBuilder {
    node_info: Vec<NodeInfo>,
    // from bottom to top
    trajectories: Vec<Vec<Vec3>>,
    particles_per_node: Vec<Vec<usize>>,
    node_props: Vec<PlantNodeProps>,
    mesh_points: Vec<Vec3>,
    debug_points: Arc<Mutex<Vec<(Vec3, Color)>>>,
    mesh_triangles: Vec<usize>,
    mesh_normals: Vec<Vec3>,
    mesh_colors: Vec<[f32; 4]>,
    contours: Arc<Mutex<Vec<Vec<Vec3>>>>,
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
            contours: Arc::new(vec![].into()),
            node_props,
            node_info,
            tree_depth,
            trajectories: vec![],
            particles_per_node: vec![vec![]; node_count],
            mesh_points: vec![],
            mesh_triangles: vec![],
            mesh_normals: vec![],
            mesh_colors: vec![],
            debug_points: Arc::new(vec![].into()),
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
    fn radius(&self, node_id: usize) -> f32 {
        self.node_props[node_id].radius
    }
    fn orientation(&self, node_id: usize) -> Vec3 {
        self.node_props[node_id].orientation
    }
    fn main_children(&self, node_id: usize) -> Option<usize> {
        self.node_info[node_id].children.get(0).copied()
    }
    fn parent(&self, node_id: usize) -> Option<usize> {
        self.node_info[node_id].parent
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
        let d = self.depth(current_node);
        self.trajectories[particle_id][d] = position;
        self.particles_per_node[current_node].push(particle_id);
    }

    fn branch_length_to_parent(&self, child: usize) -> f32 {
        let parent = self.parent(child).expect("there is no branch under root");
        (self.position(child) - self.position(parent)).length()
    }

    fn branch_length_to_main_children(&self, node: usize) -> f32 {
        let child = self.main_children(node).expect("a leaf does not have a length");
        (self.position(node) - self.position(child)).length()
    }

    fn register_points_at_position(&mut self, points: &[Vec3], pos: BranchSectionPosition, compute_normals: bool) -> Vec<usize> {
        let orientation = self.node_props[pos.node].orientation;
        let i0 = self.mesh_points.len();
        let n = points.len();
        self.mesh_points.extend(points);

        let r = (self.depth(pos.node) as f32) / self.max_depth() as f32;
        let color = [1. - r, 0.5+0.5*r, 0.2, 1.0];
        self.mesh_colors.extend(vec![color; n]);

        if compute_normals {

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
        }

        (i0..i0+n)
            .into_iter()
            .collect()
    }

    fn project_particles(&mut self, parent: usize, child: usize, offset: Vec3) {
        let origin = self.node_props[parent].position;
        let normal = self.node_props[parent].orientation;
        let r_parent = self.node_props[parent].radius;
        let p_child = self.node_props[child].position;

        let d = (origin - p_child).normalize();

        // FIXME: no clone
        for p in self.particles_per_node[child].clone() {
            let pos_particle = self.trajectories[p][self.depth(child)];
            let u = pos_particle - origin;

            let l = normal.dot(u) / normal.dot(d.normalize());
            assert!(!l.is_nan());
            let projected = origin + r_parent*(offset + (u - l * d.normalize())).normalize();

            self.register_particle_position_for_node(p, projected, parent);
        }
    }

    fn particles_on_section(&self, pos: BranchSectionPosition) -> Vec<Vec3> {
        let t = if pos.length < 0. {
            let parent = self.parent(pos.node).unwrap();
            let branch_len = (self.position(pos.node) - self.position(parent)).length();
            (branch_len + pos.length)/branch_len
        } else {
            pos.length / self.branch_length_to_main_children(pos.node)
        };
        let depth = if pos.length < 0. {self.depth(pos.node) - 1} else {self.depth(pos.node)};
        let spline_index = SplineIndex::Local(depth, t);

        self.particles_per_node[pos.node]
            .iter()
            .map(|&particle| extended_catmull_spline(&self.trajectories[particle], spline_index))
            .collect()
    }

    fn split(&self, pos: BranchSectionPosition) -> bool {
        if self.node_info[pos.node].children.len() != 2 {return false};
        let m_child = self.node_info[pos.node].children[0];
        let s_child = self.node_info[pos.node].children[1];

        let t = pos.length / self.branch_length_to_main_children(pos.node);
        let m_center = self.position(pos.node).lerp(self.position(m_child), t);
        let s_center = self.position(pos.node).lerp(self.position(s_child), t);

        let pos_along_dir = |a: Vec3| a.dot(s_center - m_center);

        let m_pos = BranchSectionPosition::new(m_child, pos.length-self.branch_length_to_parent(m_child));
        let s_pos = BranchSectionPosition::new(s_child, pos.length-self.branch_length_to_parent(s_child));
        let m_relative_pos = self.particles_on_section(m_pos).into_iter().map(pos_along_dir);
        let s_relative_pos = self.particles_on_section(s_pos).into_iter().map(pos_along_dir);

        m_relative_pos.reduce(f32::max) < s_relative_pos.reduce(f32::min)
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
                self.project_particles(root, child, Vec3::ZERO);
            },
            &[m_child, s_child] => {
                self.compute_trajectories(m_child, particle_per_leaf);
                self.compute_trajectories(s_child, particle_per_leaf);
                let normal = self.node_props[root].orientation;
                let offset_direction = {
                    let u = self.position(m_child) - self.position(s_child);
                    (u - u.dot(normal) * normal).normalize() 
                };
                self.project_particles(root, m_child, offset_direction * self.radius(m_child));
                self.project_particles(root, s_child, -offset_direction * self.radius(s_child));
            }
            _ => panic!("did not expect more than 2 childs for node {root}")
        }
    }

    fn branch_section_center(&self, pos: BranchSectionPosition) -> Vec3 {
        if pos.length < 0. {
            let parent = self.parent(pos.node).expect("there is no branch under root");
            let length = self.branch_length_to_parent(pos.node);
            self.node_props[pos.node].position.lerp(self.node_props[parent].position, -pos.length / length)
        }
        else {
            let length = self.branch_length_to_main_children(pos.node);
            let child = *self.node_info[pos.node].children.get(0).expect("there is no branch after this leaf");
            self.node_props[pos.node].position.lerp(self.node_props[child].position, pos.length / length)
        }
    }

    fn branch_contour(&self, pos: BranchSectionPosition) -> Vec<Vec3> {
        let parent = pos.node;
        let points = self.particles_on_section(pos);

        if points.len() == 0 {
            println!("node {} has no particles", pos.node);
        }

        // TODO: smarter projection
        let to_plane = |x: Vec3| 
            (Quat::from_rotation_arc(self.node_props[parent].orientation, Vec3::Z) * x)
            .truncate();


        let projected_points: Vec<Vec2> = points
            .iter()
            .map(|x: &Vec3| to_plane(*x))
            .collect();

        let center = self.branch_section_center(pos);
        self.debug_points.lock().unwrap().push((center, Color::srgb(1.0, 1.0, 1.0)));
        let center = to_plane(center);
        let result: Vec<Vec3> = meshing::convex_hull_graham(Some(center), &projected_points, Some(0.9*std::f32::consts::PI))
            .into_iter()
            // TODO: project point to avoid strands with too much diff ?
            .map(|i| points[i])
            .collect();
        self.contours.lock().unwrap().push(result.clone());
        result
    }

    fn register_triangles(&mut self, triangles: &[usize]) {
        self.mesh_triangles.extend(triangles)
    }

    pub fn compute_each_branch(&mut self) {
        let pos_root = BranchSectionPosition::new(0, 0.);
        let root_section = self.register_points_at_position(&self.branch_contour(pos_root), pos_root, true);
        self.compute_each_branch_recursive(0, root_section)
    }

    fn compute_branch_join(&mut self, p: BranchSectionPosition, previous_contour: Vec<usize>, dz: f32) -> ((BranchSectionPosition, Vec<usize>), (BranchSectionPosition, Vec<usize>)){
        assert!(self.node_info[p.node].children.len() == 2);
        let normal = self.node_props[p.node].orientation;

        let m_child = self.node_info[p.node].children[0];
        let s_child = self.node_info[p.node].children[1];

        let mut compute_properties = |child| {
            let branch_length = self.branch_length_to_parent(child);
            let pos = BranchSectionPosition::new(child, p.length - branch_length + dz);
            let center = self.branch_section_center(pos);
            let contour = self.register_points_at_position(&self.branch_contour(pos), pos, true);
            (pos, center, contour)
        };

        let (m_p, m_c, m_cont) = compute_properties(m_child);
        let (s_p, s_c, s_cont) = compute_properties(s_child);

        let center = 0.5*(m_c + s_c);

        let m_dist_center = |i: &usize| (self.mesh_points[*i] - m_c).length();
        let s_dist_center = |i: &usize| (self.mesh_points[*i] - s_c).length();
        let dist_center = |i: &usize| (self.mesh_points[*i] - center).length();

        let side = |i: &usize| Mat3::from_cols(m_c - s_c, normal, self.mesh_points[*i] - center).determinant();

        let i_m_furthest = m_cont.iter().map(s_dist_center).arg_min().unwrap() as i32;
        let i_s_furthest = s_cont.iter().map(m_dist_center).arg_min().unwrap() as i32;

        let i_a = previous_contour.iter().map(|i| side(i) / dist_center(i))
            .arg_min()
            .unwrap() as i32;
        let i_b = previous_contour.iter().map(|i| side(i) / dist_center(i))
            .arg_max()
            .unwrap() as i32;

        self.debug_points.lock().unwrap().push((self.mesh_points[previous_contour[i_a as usize]], Color::srgb(1.0, 1.0, 0.8)));
        self.debug_points.lock().unwrap().push((self.mesh_points[previous_contour[i_b as usize]], Color::srgb(0.8, 1.0, 1.0)));

        let (m_junction, m_above) = split_contour(&m_cont, i_m_furthest-1, i_m_furthest+1);
        let (mut s_junction, s_above) = split_contour(&s_cont, i_s_furthest-1, i_s_furthest+1);
        s_junction.reverse();
        let (m_under, s_under) = split_contour(&previous_contour, i_b, i_a);

        self.mesh_triangles.extend(
            meshing::mesh_between_contours(&self.mesh_points, 
            &m_under, &m_above,
            false
        ));
        self.mesh_triangles.extend(
            meshing::mesh_between_contours(&self.mesh_points, 
            &s_under, &s_above,
            false)
        );

        self.mesh_triangles.extend(
            meshing::mesh_between_contours(&self.mesh_points, 
            &s_junction, &m_junction,
            false)
        );

        self.mesh_triangles.extend([
            m_junction[0], s_junction[0], previous_contour[i_a as usize],
            s_junction[2], m_junction[2], previous_contour[i_b as usize],
        ]);

        ((m_p, m_cont), (s_p, s_cont))
    }

    fn compute_branch_while(&mut self, 
        pos: BranchSectionPosition, 
        previous_contour: &mut Vec<usize>,
        dz:f32,
        condition: impl Fn(BranchSectionPosition, &Self)->bool
        ) -> BranchSectionPosition {
        let mut p = pos;
        assert!(dz != 0.);
        while condition(p, self) {
            if p.length > 1000.*dz {
                panic!("stopping, the branch at position {p:?} is already too long")
            }
            let current_contour = self.register_points_at_position(&self.branch_contour(p), p, true);
            let triangles = meshing::mesh_between_contours(&self.mesh_points, &previous_contour, &current_contour, true); 
            self.register_triangles(&triangles);
            *previous_contour = current_contour;
            p += dz;
        }
        p
    }

    fn compute_each_branch_recursive(&mut self, root: usize, mut previous_contour: Vec<usize>) {
        let pos_root = BranchSectionPosition::new(root, 0.);

        let radius = self.radius(root);
        let dz = 0.2*radius;

        match &self.node_info[root].children[..] {
            [] => {
                let leaf = self.position(root) + self.radius(root)*self.orientation(root);
                let i_end = self.register_points_at_position(&vec![leaf], pos_root, false)[0];
                let n = previous_contour.len();
                for i in 0..n {
                    self.mesh_triangles.extend([
                        previous_contour[(i+1)%n], previous_contour[i], i_end
                    ]);
                }
            },
            &[child] => {
                let branch_length = self.branch_length_to_parent(child);
                dbg!(branch_length);
                self.compute_branch_while(pos_root, &mut previous_contour, dz, |p,_| p.length < branch_length);
                self.compute_each_branch_recursive(child, previous_contour)
            },
            &[m_child, s_child] => {
                let pos_split = self.compute_branch_while(pos_root, &mut previous_contour, dz, |p, me| 
                    // FIXME: don't pass self as argument
                    !me.split(p+dz)
                );

                let ((m_pos, mut m_cont), (s_pos, mut s_cont)) = self.compute_branch_join(pos_split, previous_contour, dz);

                self.compute_branch_while(m_pos, &mut m_cont, dz, |p,_| p.length < 0.);
                self.compute_each_branch_recursive(m_child, m_cont);

                self.compute_branch_while(s_pos, &mut s_cont, dz, |p,_| p.length < 0.);
                self.compute_each_branch_recursive(s_child, s_cont);
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
        if debug_flags.contours {
            let mut rng = StdRng::seed_from_u64(42);
            for c in self.contours.lock().unwrap().iter() {
                let t = rng.gen();
                let n = c.len();
                for i in 0..n {
                    let r = i as f32 / n as f32; 
                    let color = Color::srgb(t, 0.3+0.2*r, 0.3+0.2*r);
                    let pos1 = c[i];
                    let pos2 = c[(i+1)%n];
                    gizmos.line(pos1, pos2, color);
                }

            }
        }
        if debug_flags.other {
            for (p, c) in self.debug_points.lock().unwrap().iter() {
                gizmos.cross(*p, 0.2, *c);
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

