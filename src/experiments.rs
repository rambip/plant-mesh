use bevy::prelude::*;
use bevy::asset::RenderAssetUsages;
use bevy::render::mesh::{Indices, PrimitiveTopology};

#[derive(Component)]
pub struct Tree {
    need_render: bool,
    plant_graph: PlantNode,
    // including root
    node_count: usize,
    particle_per_leaf: usize,
    cache: TreeComputations
}

struct TreeComputations {
    trajectories: Vec<Vec<Vec3>>
}

struct PlantNode {
    position: Vec3,
    children: Vec<PlantNode>,
    id: usize,

}

fn lerp<T>(a: T, b: T, t: f32) -> T 
where T: std::ops::Mul<f32, Output=T> + std::ops::Add<Output=T>
{
    a*(1.-t) + b*t
}

fn extended_catmull_spline(points: &[Vec3], r: f32) -> Vec3 {
    let n = points.len();
    if r==1. {return points[n-2]};

    let position_of_point = r*(n as f32 - 1.);
    // edge case, we might get an index error
    let i0 = usize::min(position_of_point as usize, n-2);

    let points_to_interpolate: [Vec3; 4] = {
        let mut p = [Vec3::ZERO; 4];
        p[0] = if i0 == 0 {2.*points[1] - points[0]} else {points[i0-1]};
        p[1] = points[i0+0];
        p[2] = points[i0+1];
        p[3] = if i0 == n-2 {2.*points[n-2] - points[n-1]} else {points[i0+2]};
        p
    };
    let mut knot_sequence = [0.; 4];
    knot_sequence[0] = 0.;
    for i in 1..4 {
        knot_sequence[i] = knot_sequence[i-1] +
            (points_to_interpolate[i] - points_to_interpolate[i-1]).length()
            .sqrt();
    }
    let t = lerp(knot_sequence[1], knot_sequence[2], position_of_point - (i0 as f32));

    let ratio = |i: usize, j: usize| (t - knot_sequence[i]) / (knot_sequence[j] - knot_sequence[i]);

    let mut a_points = [Vec3::ZERO; 3];
    for i in 0..3 {
        a_points[i] = lerp(points_to_interpolate[i], points_to_interpolate[i+1], ratio(i, i+1));
    }

    let b1 = lerp(a_points[0], a_points[1], ratio(0, 2));
    let b2 = lerp(a_points[1], a_points[2], ratio(1, 3));

    lerp(b1, b2, ratio(1, 2))
}

fn uniform_square() -> Vec2 {
    Vec2::new(
        rand::random(),
        rand::random(),
    )
}

impl PlantNode {
    fn demo() -> Self {
        Self {
            position: Vec3::ZERO,
            id: 0,
            children: vec![
                PlantNode {
                    position: Vec3::new(0., 0., 2.),
                    id: 6,
                    children: vec![
                        PlantNode {
                            position: Vec3::new(1., 0., 3.),
                            id: 2,
                            children: vec![
                                PlantNode {
                                    position: Vec3::new(1.5, 0., 4.),
                                    id: 4,
                                    children: vec![],
                                },
                                PlantNode {
                                    position: Vec3::new(1.5, -1., 5.),
                                    id: 5,
                                    children: vec![],
                                }
                            ],
                        },
                        PlantNode {
                            position: Vec3::new(-1., 0., 4.),
                            id: 1,
                            children: vec![
                                PlantNode {
                                    position: Vec3::new(-0.5, 1., 6.),
                                    id: 3,
                                    children: vec![],
                                } ],
                        } ],
                }],
        }
    }


    fn debug(&self, gizmos: &mut Gizmos) {
        gizmos.circle(Isometry3d::from_translation(self.position), 0.2, Color::srgb(0., 0.8, 0.5));
        for c in &self.children {
            c.debug(gizmos)
        }

    }

    fn register_node_positions(&self, acc: &mut Vec<Vec3>) {
        acc[self.id] = self.position;
        for c in &self.children {
            c.register_node_positions(acc)
        }
    }

    fn register_parents(&self, acc: &mut Vec<usize>, parent: usize) {
        acc[self.id] = parent;
        for c in &self.children {
            c.register_parents(acc, self.id)
        }
    }

    fn register_leaves(&self, acc: &mut Vec<usize>) {
        if self.children.len() == 0 {
            acc.push(self.id);
        }
        for c in &self.children {
            c.register_leaves(acc);
        }
    }
}

impl Default for Tree {
    fn default() -> Self {
        Self {
            need_render: true,
            plant_graph: PlantNode::demo(),
            // FIXME: compute it
            node_count: 7,
            particle_per_leaf: 10,
            cache: TreeComputations {trajectories: Vec::new()}
        }
    }
}

fn dummy_mesh() -> Mesh {
    Mesh::new(PrimitiveTopology::TriangleList, RenderAssetUsages::default())
        // Add 4 vertices, each with its own position attribute (coordinate in
        // 3D space), for each of the corners of the parallelogram.
        .with_inserted_attribute(
            Mesh::ATTRIBUTE_POSITION,
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 2.0], [2.0, 0.0, 2.0]],
        )
        // Assign a UV coordinate to each vertex.
        .with_inserted_attribute(
            Mesh::ATTRIBUTE_UV_0,
            vec![[0.0, 1.0], [0.5, 0.0], [1.0, 0.0]]
        )
        // Assign normals (everything points outwards)
        .with_inserted_attribute(
            Mesh::ATTRIBUTE_NORMAL,
            vec![[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]
        )
        // After defining all the vertices and their attributes, build each triangle using the
        // indices of the vertices that make it up in a counter-clockwise order.
        .with_inserted_indices(Indices::U32(vec![ 0, 1, 2, ]))
}

impl Tree {

    fn render_mesh(&mut self) -> Mesh {
        self.cache.trajectories = self.compute_trajectories();
        // TODO:
        // Etape 1: calcul des trajectoires
        // - calculer pour chaque noeud l'ensemble des points de la trajectoire.
        // - stocker pour chaque noeud l'ensemble des points qui y passent, dans l'ordre croissant 
        // des ID des particules -> inutile ?
        // - penser à faire les collisions et tout
        //   - PBR = ?
        //
        // Etape 2: précalculs
        // - construire le graphe des parents
        // - calculer la profondeur de chaque noeud
        // - établir une fonction (temps constant) qui retourne la position interpolée d'une particule entre le
        // noeud a et le noeud b avec une progression t
        //
        // Etape 2: mesh
        // - pour chaque noeud:
        //   - si c'est une feuille
        //      - ne rien faire
        //   - si ce n'est pas une feuille
        //      - lister les particules dans ce noeud (linéaire)
        //      - pour t entre 0 et 1
        //        - calculer la position interpolée de chaque particule
        //        - calculer l'envelope convexe des positions interpolées
        //          - extension: delaunay avec seuil
        //        - ajouter ces points à la mesh
        //        - ajouter des triangles pour relier (couplage de graphe ?)
        dummy_mesh()
    }

    fn compute_trajectories(&self) -> Vec<Vec<Vec3>> {
        let parents : Vec<usize> = {
            let mut p = vec![0; self.node_count];
            self.plant_graph.register_parents(&mut p, 0);
            p
        };
        let node_positions: Vec<Vec3> = {
            let mut po = vec![Vec3::ZERO; self.node_count];
            self.plant_graph.register_node_positions(&mut po);
            po
        };

        let leaves: Vec<usize> = {
            let mut l = Vec::new();
            self.plant_graph.register_leaves(&mut l);
            l
        };

        let mut trajectories = Vec::new();
        for l in leaves {
            for _ in 0..self.particle_per_leaf {
                let mut new_traj = Vec::new();
                let start = node_positions[l] + (uniform_square()-Vec2::new(0.5, 0.5)).extend(0.);
                new_traj.push(start);
                let mut node = l;
                while node != 0 {
                    node = parents[node];
                    let pos = node_positions[node] + (uniform_square()-Vec2::new(0.5, 0.5)).extend(0.);
                    new_traj.push(pos);
                }
                trajectories.push(new_traj);
            }
        }

        trajectories
    }

    fn debug(&self, gizmos: &mut Gizmos) {
        self.plant_graph.debug(gizmos);

        for traj in &self.cache.trajectories {
            for i in 0..50 {
                let t = (i as f32) / 50.;
                let pos = extended_catmull_spline(traj, t);
                let color = Color::srgb(t, 0.8, 0.);
                gizmos.cross(Isometry3d::from_translation(pos), 0.1, color);
            }
        }
    }
}

// TODO: multiple materials for one single tree
pub fn draw_tree(
    mut gizmos: Gizmos,
    mut trees: Query<(&mut Mesh3d, &mut Tree)>,
    mut meshes: ResMut<Assets<Mesh>>,
    ) {


    for (mut mesh, mut tree) in trees.iter_mut() {
        tree.debug(&mut gizmos);
        if tree.need_render {
            println!("update mesh");
            mesh.0 = meshes.add(tree.render_mesh());
            tree.need_render = false;
        }
    }
}
