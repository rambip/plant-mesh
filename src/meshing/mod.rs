use std::ops::{Add, AddAssign};
use std::sync::{Arc, Mutex};

use bevy::math::{Quat, Vec2, Vec3};
use bevy::prelude::{Mesh, Color};
use bevy::asset::RenderAssetUsages;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy_gizmos::prelude::Gizmos;

use crate::growing::{PlantNode, PlantNodeProps, NodeInfo};
use crate::tools::{max_by_key, min_by_key};
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

    fn branch_section_center(&self, pos: BranchSectionPosition) -> Vec3 {
        match (self.node_info[pos.node].parent, &self.node_info[pos.node].children[..], pos.t) {
            (_, &[a], t) if t >= 0. => self.node_props[pos.node].position.lerp(self.node_props[a].position, t),
            // a is main branch
            (_, &[a, _], t) if t >= 0. => self.node_props[pos.node].position.lerp(self.node_props[a].position, t),
            (Some(p), _, t) if t <= 0. => self.node_props[pos.node].position.lerp(self.node_props[p].position, -t),
            (_, &[], _) => panic!("node {} has no children so it has no center", pos.node),
            _ => panic!("unable to compute branch center position {pos:?}. It has {:?} childrens and its parent is {:?}", self.node_info[pos.node].children, self.node_info[pos.node].parent)
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

        let center = to_plane(self.branch_section_center(pos));
        self.debug_points.lock().unwrap().push((self.branch_section_center(pos), Color::srgb(1.0, 1.0, 1.0)));
        meshing::convex_hull_graham(Some(center), &projected_points)
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

    fn compute_branch_join(&mut self, p: BranchSectionPosition, previous_contour: Vec<usize>) -> ((BranchSectionPosition, Vec<usize>), (BranchSectionPosition, Vec<usize>)){
        assert!(self.node_info[p.node].children.len() == 2);
        let c1 = self.node_info[p.node].children[0];
        let c2 = self.node_info[p.node].children[1];

        let p1 = BranchSectionPosition::new(c1, p.t - 1.);
        let p2 = BranchSectionPosition::new(c2, p.t - 1.);
        let cont1 = self.register_points_on_contour(&self.branch_contour(p1), p1);
        let n1 = cont1.len();
        let cont2 = self.register_points_on_contour(&self.branch_contour(p2), p2);
        let n2 = cont2.len();

        let center_1 = self.branch_section_center(p1);
        let center_2 = self.branch_section_center(p2);

        let i_p = min_by_key(0..n1, |&i| (self.mesh_points[cont1[i]] - center_2 ).length()).unwrap() as i32;
        let i_q = min_by_key(0..n2, |&i| (self.mesh_points[cont2[i]] - center_1).length()).unwrap() as i32;

        let (m_junction, m_above) = split_contour(&cont1, i_p-1, i_p+1);
        let (mut s_junction, s_above) = split_contour(&cont2, i_q-1, i_q+1);
        s_junction.reverse();
        assert_eq!(m_junction.len(), 3);
        assert_eq!(s_junction.len(), 3);

        let i_a = min_by_key(0..previous_contour.len(),
            |&i| (self.mesh_points[previous_contour[i]] - self.mesh_points[m_junction[0]]).length() 
               + (self.mesh_points[previous_contour[i]] - self.mesh_points[s_junction[0]]).length()
        ).unwrap() as i32;
        let i_b = min_by_key(0..previous_contour.len(),
            |&i| (self.mesh_points[previous_contour[i]] - self.mesh_points[m_junction[2]]).length() 
               + (self.mesh_points[previous_contour[i]] - self.mesh_points[s_junction[2]]).length()
        ).unwrap() as i32;

        let (m_under, s_under) = split_contour(&previous_contour, i_b, i_a);

        let mut debug_points = self.debug_points.lock().unwrap();
        for &i in &m_junction {
            debug_points.push((self.mesh_points[i], Color::srgb(1.0, 1.0, 0.0)));
        }
        for &i in &s_junction {
            debug_points.push((self.mesh_points[i], Color::srgb(1.0, 0.0, 1.0)));
        }

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

        ((p1, cont1), (p2, cont2))
    }

    fn compute_branch_until(&mut self, 
        pos: BranchSectionPosition, 
        previous_contour: &mut Vec<usize>,
        dt:f32,
        condition: impl Fn(BranchSectionPosition, &mut Self)->bool
        ) -> Result<(), BranchSectionPosition> {
        let mut p = pos;
        while p.decimal_part()+dt <= 1. {
            if condition(p, self) {return Err(p)};
            let current_contour = self.register_points_on_contour(&self.branch_contour(p), p);
            let triangles = meshing::mesh_between_contours(&self.mesh_points, &previous_contour, &current_contour, true); 
            self.register_triangles(&triangles);
            *previous_contour = current_contour;
            p += dt;
        }
        Ok(())
    }

    fn compute_each_branch_recursive(&mut self, root: usize, start_section: Vec<usize>) {
        let mut previous_contour = start_section;

        let pos_root = BranchSectionPosition::new(root, 0.);

        let radius = self.node_props[root].radius;
        let dz = 2.0*radius * std::f32::consts::PI / previous_contour.len() as f32;

        match &self.node_info[root].children[..] {
            [] => {},
            &[child] => {
                let branch_length = (self.position(root) - self.position(child)).length();
                self.compute_branch_until(pos_root, &mut previous_contour, dz/branch_length, |_,_| false).unwrap();
                self.compute_each_branch_recursive(child, previous_contour)
            },
            // TODO: assume child1 is the main branch
            &[child1, child2] => {
                let branch_length = 0.5*((self.position(root) - self.position(child1)).length()
                    + (self.position(root) - self.position(child2)).length());
                let pos_split = self.compute_branch_until(pos_root, &mut previous_contour, dz/branch_length, |t, me| 
                    // FIXME: don't pass self as argument
                    me.split(t)
                    ).err().unwrap();

                let ((p1, mut cont1), (p2, mut cont2)) = self.compute_branch_join(pos_split, previous_contour);

                self.compute_branch_until(p1, &mut cont1, dz/branch_length, |_,_| false).unwrap();
                self.compute_each_branch_recursive(child1, cont1);

                self.compute_branch_until(p2, &mut cont2, dz/branch_length, |_,_| false).unwrap();
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

