use bevy::math::{Quat, Vec2, Vec3};
use bevy::prelude::{Mesh, Color};
use bevy::asset::RenderAssetUsages;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy_gizmos::prelude::Gizmos;
use meshing::lerp;

use crate::growing::{PlantNode, PlantNodeProps};
mod meshing;

// TODO: no heap allocation
pub struct MeshBuilder {
    node_count: usize,
    // reversed: from bottom to top
    trajectories: Vec<Vec<Vec3>>,
    depths: Vec<usize>,
    parents: Vec<usize>,
    leaves: Vec<usize>,
    particles_per_node: Vec<Vec<usize>>,
    node_props: Vec<PlantNodeProps>,
    mesh_points: Vec<Vec3>,
    debug_points: Vec<Vec3>,
    mesh_triangles: Vec<usize>,
    mesh_normals: Vec<Vec3>,
    mesh_colors: Vec<[f32; 4]>,
    tree_depth: usize,
    children: Vec<Vec<usize>>,
    topological_sort: Vec<usize>,
}

impl MeshBuilder {
    pub fn new(plant_graph: &PlantNode, node_count: usize) -> Self {
        let mut parents = Vec::new();
        plant_graph.register_parents(&mut parents, 0);

        let mut node_props = Vec::new();
        plant_graph.register_node_properties(&mut node_props);

        let mut depths = Vec::new();
        plant_graph.register_depths(&mut depths, 0);

        let mut leaves = Vec::new();
        plant_graph.register_leaves(&mut 0, &mut leaves);

        let mut topological_sort = Vec::new();
        plant_graph.register_node_postfix(&mut 0, &mut topological_sort);

        let mut children = Vec::new();
        plant_graph.register_node_children(&mut children);

        Self {
            node_count,
            depths,
            parents,
            node_props,
            leaves,
            tree_depth: plant_graph.compute_depth(),
            trajectories: vec![],
            particles_per_node: vec![vec![]; node_count],
            mesh_points: vec![],
            mesh_triangles: vec![],
            mesh_normals: vec![],
            mesh_colors: vec![],
            debug_points: vec![],
            topological_sort,
            children,
        }
    }

    fn depth(&self) -> usize {
        self.tree_depth
    }

    fn register_particle_position_for_leaf(&mut self, leaf_id: usize) {
        let PlantNodeProps {
            position,
            radius,
            orientation,
        } = self.node_props[leaf_id];
        let particle_id = self.trajectories.len();
        let rotation = Quat::from_rotation_arc(Vec3::Z, orientation);
        let relative_pos = rotation * sample_uniform_circle(radius).extend(0.);
        self.particles_per_node[leaf_id].push(particle_id);
        let mut empty_trajectory = vec![Vec3::ZERO; self.depths[leaf_id]];
        empty_trajectory.push(position+relative_pos);

        self.trajectories.push(empty_trajectory);
    }

    fn register_particle_position_for_node(&mut self, particle_id: usize, position: Vec3, current_node: usize) {
        self.trajectories[particle_id][self.depths[current_node]] = position;
        self.particles_per_node[current_node].push(particle_id);
    }

    fn global_t(&self, bottom_node: usize, t: f32) -> f32 {
        self.depths[bottom_node] as f32 + t
    }

    fn particle_position(&self, particle_id: usize, bottom_node: usize, t: f32) -> Vec3 {
        meshing::extended_catmull_spline(&self.trajectories[particle_id], 
            self.global_t(bottom_node, t)
        )
    }
    pub fn register_points_on_contour(&mut self, points: &[Vec3], global_depth: f32, orientation: Vec3) -> Vec<usize> {
        let i0 = self.mesh_points.len();
        let n = points.len();
        self.mesh_points.extend(points);
        //self.cache.debug_points.push(points[0]);
        let r = global_depth / self.depth() as f32;
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
            let pos_particle = self.trajectories[p][self.depths[child]];
            let u = pos_particle - origin;

            let projected_offset = p_child - origin - (p_child - origin).dot(normal)*normal;

            let l = normal.dot(u) / normal.dot(d.normalize());
            assert!(!l.is_nan());
            let o = if offset {0.5*projected_offset} else {Vec3::ZERO};
            let projected = origin + r_parent*(o + (u - l * d.normalize())).normalize();

            self.register_particle_position_for_node(p, projected, parent);
            //self.trajectories[p].push(projected);
        }
    }


    pub fn compute_trajectories(&mut self, particle_per_leaf: usize) {
        // FIXME: don't clone
        for parent in self.topological_sort.clone() {
            println!("parent: {}", parent);
            assert!(self.particles_per_node[parent].len()==0);
            match &self.children[parent][..] {
                [] => {
                    for _ in 0..particle_per_leaf {
                        self.register_particle_position_for_leaf(parent);
                    }
                },
                &[child] => {
                    self.project_particles(parent, child, false);
                },
                &[child1, child2] => {
                    let mut big_child = child1;
                    let mut small_child = child2;
                    if self.particles_per_node[big_child].len() < self.particles_per_node[small_child].len() {
                        std::mem::swap(&mut big_child, &mut small_child);
                    }
                    self.project_particles(parent, small_child, false);
                    self.project_particles(parent, big_child, true);
                }
                _ => panic!("did not expect more than 2 childs for node {parent}")
            }
            dbg!(&self.particles_per_node);
        }
    }
    pub fn branch_contour(&self, top_node: usize, t: f32) -> Vec<Vec3> {
        let parent = self.parents[top_node];
        let particles = &self.particles_per_node[top_node];

        let points: Vec<Vec3> = particles
            .iter()
            .map(|&p| self.particle_position(p, parent, t))
            .collect();

        if points.len() == 0 {
            println!("node {} has no particles", top_node);
        }
        meshing::convex_hull_graham(&points)
            .into_iter()
            .map(|i| points[i])
            .collect()
    }

    fn register_triangles(&mut self, triangles: &[usize]) {
        self.mesh_triangles.extend(triangles)
    }

    pub fn compute_each_branch(&mut self) {
        for child in 1..self.node_count{

            let parent = self.parents[child];


            let new_points = self.branch_contour(child, 0.);
            let mut previous_contour_ids = self.register_points_on_contour(
                &new_points, 
                self.global_t(parent, 0.),
                self.node_props[parent].orientation
            );

            let radius = self.node_props[parent].radius;
            let dz = 2.0*radius * std::f32::consts::PI / new_points.len() as f32;
            let branch_length = (self.node_props[parent].position - self.node_props[child].position).length();
            let n_steps = (branch_length / dz) as usize;

            // FIXME: duplicate points at nodes
            for i in 1..=n_steps {
                let t = (i as f32) / n_steps as f32;
                let orientation = lerp(self.node_props[parent].orientation, self.node_props[child].orientation, t);

                let current_contour = self.branch_contour(child, t);
                let current_contour_ids = self.register_points_on_contour(
                    &current_contour,
                    self.global_t(parent, t),
                    orientation
                );

                let triangles = meshing::mesh_between_contours(&self.mesh_points, &current_contour_ids, &previous_contour_ids); 

                self.register_triangles(&triangles);

                previous_contour_ids = current_contour_ids.clone();
            }
        }
    }

    pub fn to_mesh(&self) -> Mesh {
        Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::default())
            .with_inserted_attribute(
                Mesh::ATTRIBUTE_POSITION,
                self.mesh_points.clone(),
            )
            .with_inserted_attribute(
                Mesh::ATTRIBUTE_NORMAL,
                self.mesh_normals.clone(),
            )
            .with_inserted_attribute(
                Mesh::ATTRIBUTE_COLOR,
                self.mesh_colors.clone(),
            )
            .with_inserted_indices(Indices::U32(
                    self.mesh_triangles.iter().map(|x| *x as u32).collect()
            ))
    }
    
    pub fn debug(&self, gizmos: &mut Gizmos) {
        for i in 0..self.mesh_triangles.len()/3 {
            let (ia, ib, ic) = (self.mesh_triangles[3*i], self.mesh_triangles[3*i+1], self.mesh_triangles[3*i+2]);
            let (pa, pb, pc) = (self.mesh_points[ia], self.mesh_points[ib], self.mesh_points[ic]);
            let color = Color::srgb(0., 0.4, 0.);
            gizmos.line(pa, pb,  color);
            gizmos.line(pb, pc,  color);
            gizmos.line(pc, pa,  color);
        }
        for i in 0..self.mesh_normals.len() {
            let p = self.mesh_points[i];
            let normal = self.mesh_normals[i];
            let color = Color::srgb(1., 0.0, 0.);
            //gizmos.cross(Isometry3d::from_translation(p), 0.1, color);
            gizmos.line(p, p+0.3*normal, color);
        }
    }
}


fn sample_uniform_disk(radius: f32) -> Vec2 {
    Vec2::from_angle(rand::random::<f32>()*std::f32::consts::TAU) * radius*(rand::random::<f32>()).sqrt()
}

fn sample_uniform_circle(radius: f32) -> Vec2 {
    Vec2::from_angle(rand::random::<f32>()*std::f32::consts::TAU) * radius
}

