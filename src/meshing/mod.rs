use bevy::math::{Vec3, Vec2, Isometry3d};
use bevy::prelude::{Mesh, Color};
use bevy::asset::RenderAssetUsages;
use bevy::render::mesh::{Indices, PrimitiveTopology};

use crate::growing::{PlantNode, PlantNodeProps};
mod meshing;

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
    max_triangle_proportion: f32,
    mesh_normals: Vec<Vec3>,
    mesh_colors: Vec<[f32; 4]>,
}

impl MeshBuilder {
    pub fn new(plant_graph: &PlantNode, node_count: usize) -> Self {
        let mut parents = vec![0; node_count];
        plant_graph.register_parents(&mut parents, 0);

        let mut node_props = vec![PlantNodeProps::default(); node_count];
        plant_graph.register_node_properties(&mut node_props);

        let mut depths = vec![0; node_count];
        plant_graph.register_depths(&mut depths, 0);

        let mut leaves = vec![];
        plant_graph.register_leaves(&mut leaves);

        Self {
            node_count,
            depths,
            parents,
            node_props,
            leaves,
            trajectories: vec![vec![]; node_count],
            particles_per_node: vec![vec![]; node_count],
            mesh_points: vec![],
            mesh_triangles: vec![],
            max_triangle_proportion: 1.,
            mesh_normals: vec![],
            mesh_colors: vec![],
            debug_points: vec![],
        }
    }
    fn register_particle_position(&mut self, particle_id: usize, trajectory: &mut Vec<Vec3>, current_node: usize) {
        let PlantNodeProps {
            position,
            radius,
            ..
        } = self.node_props[current_node];
        let p = position + uniform_disk(radius).extend(0.);
        trajectory.push(p);
        self.particles_per_node[current_node].push(particle_id);
    }
    fn particle_position(&self, particle_id: usize, bottom_node: usize, t: f32) -> Vec3 {
        let global_t = self.depths[bottom_node] as f32 + t;
        meshing::extended_catmull_spline(&self.trajectories[particle_id], 
            global_t
            )
    }
    pub fn register_points_on_contour(&mut self, points: &[Vec3]) -> Vec<usize> {
        let i0 = self.mesh_points.len();
        let n = points.len();
        self.mesh_points.extend(points);
        //self.cache.debug_points.push(points[0]);
        self.mesh_colors.extend(vec![[1.0, 1.0, 1.0, 1.]; n]);

        for i in 0..n {
            let mean = 0.4*(
                points[(n+i-2)%n] 
                + points[(n+i-1)%n]
                + points[(n+i+1)%n]
                + points[(n+i+2)%n]
            ).normalize();

            self.mesh_normals.push((points[i] - mean).normalize());
        }

        (i0..i0+n)
            .into_iter()
            .collect()
    }
    pub fn compute_trajectories(&mut self, particle_per_leaf: usize) {
        for l in self.leaves.clone() {
            for _ in 0..particle_per_leaf {
                let particle_id = self.trajectories.len();

                let mut new_traj = Vec::new();
                let mut node = l;

                self.register_particle_position(particle_id, &mut new_traj, node);

                while node != 0 {
                    node = self.parents[node];
                    self.register_particle_position(particle_id, &mut new_traj, node);
                }
                new_traj.reverse();
                self.trajectories.push(new_traj);
            }
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

    fn triangles(&self) -> &[usize] {
        let target_n_triangles : f32 = 
            // at least one triangle (3 points)
            f32::max(1., self.max_triangle_proportion*self.mesh_triangles.len() as f32 / 3.
        );
        &self.mesh_triangles[0..3*(target_n_triangles as usize)]
    }

    pub fn compute_each_branch(&mut self) {
        for child in 1..self.node_count{

            let parent = self.parents[child];

            let new_points = self.branch_contour(child, 0.);
            let mut previous_contour_ids = self.register_points_on_contour(&new_points);

            let radius = self.node_props[parent].radius;
            let dz = radius * std::f32::consts::PI / new_points.len() as f32;
            let branch_length = (self.node_props[parent].position - self.node_props[child].position).length();
            let n_steps = (branch_length / dz) as usize;

            // FIXME: duplicate points at nodes
            for i in 1..=n_steps {
                let t = (i as f32) / n_steps as f32;

                let current_contour = self.branch_contour(child, t);
                let current_contour_ids = self.register_points_on_contour(&current_contour);

                let triangles = meshing::mesh_between_contours(&self.mesh_points, &current_contour_ids, &previous_contour_ids); 

                self.register_triangles(&triangles);

                previous_contour_ids = current_contour_ids.clone();
            }
        }
    }

    //pub fn keep_only_propotion_of_triangles(&mut self, r: f32) {
    //    let target_length = r * self.mesh_triangles.len() as f32;
    //    self.mesh_triangles.truncate(target_length as usize);
    //}
    pub fn set_triangle_proportion(&mut self, r: f32){
        self.max_triangle_proportion = r;
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
                    self.triangles().iter().map(|x| *x as u32).collect()
            ))
    }
    
    pub fn debug(&self, gizmos: &mut bevy::gizmos::gizmos::Gizmos) {
        for i in 0..self.triangles().len()/3 {
            let (ia, ib, ic) = (self.mesh_triangles[3*i], self.mesh_triangles[3*i+1], self.mesh_triangles[3*i+2]);
            let (pa, pb, pc) = (self.mesh_points[ia], self.mesh_points[ib], self.mesh_points[ic]);
            let color = Color::srgb(0., 0.4, 0.);
            gizmos.line(pa, pb,  color);
            gizmos.line(pb, pc,  color);
            gizmos.line(pc, pa,  color);
        }
        for &p in &self.debug_points {
            let color = Color::srgb(1., 0.0, 0.);
            gizmos.cross(Isometry3d::from_translation(p), 0.1, color);
        }
    }
}


fn uniform_disk(radius: f32) -> Vec2 {
    Vec2::from_angle(rand::random::<f32>()*std::f32::consts::TAU) * radius*(rand::random::<f32>()).sqrt()
}


