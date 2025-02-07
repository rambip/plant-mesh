use std::ops::{Add, AddAssign};
use std::sync::{Arc, Mutex};

use bevy::math::{Mat3, Quat, Vec2, Vec3};
use bevy::prelude::{Mesh, Color};
use bevy::asset::RenderAssetUsages;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy_gizmos::prelude::Gizmos;
use rand::prelude::*;

use crate::growing::{PlantNode, PlantNodeProps, NodeInfo};
use crate::tools::min_by_key;
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

    fn branch_length_parent(&self, child: usize) -> f32 {
        let parent = self.node_info[child].parent.unwrap();
        (self.position(child) - self.position(parent)).length()
    }

    fn branch_length_smallest_children(&self, node: usize) -> f32 {
        match &self.node_info[node].children[..] {
            &[] => panic!("a leaf does not have a length"),
            &[child] => (self.position(node) - self.position(child)).length(),
            &[child1, child2] => f32::min(
                (self.position(node) - self.position(child1)).length(),
                (self.position(node) - self.position(child2)).length()
            ),
            _ => panic!("node has more than 2 children")
        }
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
        let t = if pos.length < 0. {
            let parent = self.node_info[pos.node].parent.unwrap();
            let branch_len = (self.position(pos.node) - self.position(parent)).length();
            (branch_len + pos.length)/branch_len
        } else {
            pos.length / self.branch_length_smallest_children(pos.node)
        };
        let depth = if pos.length < 0. {self.depth(pos.node) - 1} else {self.depth(pos.node)};
        let spline_index = SplineIndex::Local(depth, t);

        self.particles_per_node[pos.node]
            .iter()
            .map(|&particle| extended_catmull_spline(&self.trajectories[particle], spline_index))
            .collect()
    }

    fn split(&mut self, pos: BranchSectionPosition) -> bool {
        if self.node_info[pos.node].children.len() != 2 {return false};
        let child1 = self.node_info[pos.node].children[0];
        let child2 = self.node_info[pos.node].children[1];

        let t = pos.length / self.branch_length_smallest_children(pos.node);
        let center_1 = self.position(pos.node).lerp(self.position(child1), t);
        let center_2 = self.position(pos.node).lerp(self.position(child2), t);

        let pos_along_dir = |a: Vec3| a.dot(center_2 - center_1);
        let relative_pos_1 = self.particles_on_section(BranchSectionPosition::new(child1, pos.length-self.branch_length_parent(child1))).into_iter().map(pos_along_dir);
        let relative_pos_2 = self.particles_on_section(BranchSectionPosition::new(child2, pos.length-self.branch_length_parent(child2))).into_iter().map(pos_along_dir);

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

    fn branch_section_center(&self, pos: BranchSectionPosition) -> Vec3 {
        if pos.length < 0. {
            let parent = self.node_info[pos.node].parent.unwrap();
            let length = self.branch_length_parent(pos.node);
            self.node_props[pos.node].position.lerp(self.node_props[parent].position, -pos.length / length)
        }
        else {
            let length = self.branch_length_smallest_children(pos.node);

            let pos_according_to_child = |child: &usize| {
                self.node_props[pos.node].position.lerp(self.node_props[*child].position, pos.length / length)
            };

            let children = &self.node_info[pos.node].children;
            children.iter()
                .map(pos_according_to_child)
                .sum::<Vec3>() / children.len() as f32
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
        let c1 = self.node_info[p.node].children[0];
        let c2 = self.node_info[p.node].children[1];

        let p1 = BranchSectionPosition::new(c1, p.length - self.branch_length_parent(c1)+2.*dz);
        let p2 = BranchSectionPosition::new(c2, p.length - self.branch_length_parent(c2)+2.*dz);
        let cont1 = self.register_points_at_position(&self.branch_contour(p1), p1, true);
        let n1 = cont1.len();
        let cont2 = self.register_points_at_position(&self.branch_contour(p2), p2, true);
        let n2 = cont2.len();

        let center_1 = self.branch_section_center(p1);
        let center_2 = self.branch_section_center(p2);
        let center = 0.5*(center_1 + center_2);

        let normal = self.node_props[p.node].orientation;

        let i_p = min_by_key(0..n1, |&i| (self.mesh_points[cont1[i]] - center_2 ).length()).unwrap() as i32;
        let i_q = min_by_key(0..n2, |&i| (self.mesh_points[cont2[i]] - center_1).length()).unwrap() as i32;

        let (m_junction, m_above) = split_contour(&cont1, i_p-2, i_p+2);
        let (mut s_junction, s_above) = split_contour(&cont2, i_q-2, i_q+2);
        s_junction.reverse();
        assert_eq!(m_junction.len(), 5);
        assert_eq!(s_junction.len(), 5);

        let side = |p: Vec3| Mat3::from_cols(center_1 - center_2, normal, p - center).determinant();
        let mut side_a = Vec::new();
        let mut side_b = Vec::new();
        for i in 0..previous_contour.len() {
            let point = self.mesh_points[previous_contour[i]];
            if side(point) >= 0. {
                side_b.push(i);
                //self.debug_points.lock().unwrap().push((point, Color::srgb(1.0, 1.0, 0.5)));
            }
            else {
                side_a.push(i);
                //self.debug_points.lock().unwrap().push((point, Color::srgb(0.5, 1.0, 1.0)));
            }
        }

        let i_a = min_by_key(side_a,
            |&i| (self.mesh_points[previous_contour[i]] - self.mesh_points[m_junction[0]]).length() 
               + (self.mesh_points[previous_contour[i]] - self.mesh_points[s_junction[0]]).length()
        ).unwrap() as i32;
        let i_b = min_by_key(side_b,
            |&i| (self.mesh_points[previous_contour[i]] - self.mesh_points[m_junction[4]]).length() 
               + (self.mesh_points[previous_contour[i]] - self.mesh_points[s_junction[4]]).length()
        ).unwrap() as i32;

        self.debug_points.lock().unwrap().push((self.mesh_points[previous_contour[i_a as usize]], Color::srgb(1.0, 1.0, 0.8)));
        self.debug_points.lock().unwrap().push((self.mesh_points[previous_contour[i_b as usize]], Color::srgb(0.8, 1.0, 1.0)));

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
            s_junction[4], m_junction[4], previous_contour[i_b as usize],
        ]);

        ((p1, cont1), (p2, cont2))
    }

    fn compute_branch_until(&mut self, 
        pos: BranchSectionPosition, 
        previous_contour: &mut Vec<usize>,
        dz:f32,
        condition: impl Fn(BranchSectionPosition, &mut Self)->bool
        ) -> Result<(), BranchSectionPosition> {
        let mut p = pos;
        // TODO: remove this strange condition by changing the logic in
        // `compute_each_branch_recursive` ?
        let end_length = if p.length < 0. {0.} else {self.branch_length_smallest_children(pos.node)};
        while p.length <= end_length {
            if condition(p, self) {return Err(p)};
            let current_contour = self.register_points_at_position(&self.branch_contour(p), p, true);
            let triangles = meshing::mesh_between_contours(&self.mesh_points, &previous_contour, &current_contour, true); 
            self.register_triangles(&triangles);
            *previous_contour = current_contour;
            p += dz;
        }
        Ok(())
    }

    fn compute_each_branch_recursive(&mut self, root: usize, mut previous_contour: Vec<usize>) {
        let pos_root = BranchSectionPosition::new(root, 0.);

        let radius = self.node_props[root].radius;
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
                self.compute_branch_until(pos_root, &mut previous_contour, dz, |_,_| false).unwrap();
                self.compute_each_branch_recursive(child, previous_contour)
            },
            // TODO: assume child1 is the main branch
            &[child1, child2] => {
                let pos_split = self.compute_branch_until(pos_root, &mut previous_contour, dz, |p, me| 
                    // FIXME: don't pass self as argument
                    me.split(p+2.*dz)
                    ).err().expect("2 children but branch does not split !");

                let ((p1, mut cont1), (p2, mut cont2)) = self.compute_branch_join(pos_split, previous_contour, dz);

                self.compute_branch_until(p1, &mut cont1, dz, |_,_| false).unwrap();
                self.compute_each_branch_recursive(child1, cont1);

                self.compute_branch_until(p2, &mut cont2, dz, |_,_| false).unwrap();
                self.compute_each_branch_recursive(child2, cont2);
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

