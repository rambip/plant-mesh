use std::ops::{Add, AddAssign};

use bevy::color::Color;
use bevy::math::{Isometry3d, Quat, Vec2, Vec3};
use bevy::prelude::Component;
use bevy_gizmos::prelude::Gizmos;
use rand::Rng;
use smallvec::SmallVec;

use crate::{NodeInfo, PlantNodeProps, VisualDebug};


#[derive(Debug)]
pub struct PlantNode {
    children: Vec<PlantNode>,
    props: PlantNodeProps,
}


impl PlantNode {
    pub fn basic_random() -> Self {
        let mut rng = rand::rngs::ThreadRng::default();
        let root = PlantNodeProps {
            radius: 1.,
            orientation: Quat::default(),
            position: Vec3::ZERO,
        };
        grow_tree_basic(&mut rng, root, 0, 0.1)
    }
    pub fn demo() -> Self {
        Self {
            props: PlantNodeProps::new(Vec3::new(0., 0., -2.), 1.1, Vec3::new(0., 0., 1.)),
            children: vec![PlantNode {
                props: PlantNodeProps::new(Vec3::new(0., 0., 2.), 0.9, Vec3::new(0., 0., 1.)),
                children: vec![
                    PlantNode {
                        props: PlantNodeProps::new(
                            Vec3::new(1., 0., 3.),
                            0.5,
                            Vec3::new(0., 0., 1.),
                        ),
                        children: vec![
                            PlantNode {
                                props: PlantNodeProps::new(
                                    Vec3::new(1.5, 0., 6.),
                                    0.3,
                                    Vec3::new(0., 0., 1.),
                                ),
                                children: vec![],
                            },
                            PlantNode {
                                props: PlantNodeProps::new(
                                    Vec3::new(2., -2., 5.),
                                    0.3,
                                    Vec3::new(0., -2., 1.),
                                ),
                                children: vec![],
                            },
                        ],
                    },
                    PlantNode {
                        props: PlantNodeProps::new(
                            Vec3::new(-1., 0., 3.5),
                            0.3,
                            Vec3::new(0., 1., 1.),
                        ),
                        children: vec![PlantNode {
                            props: PlantNodeProps::new(
                                Vec3::new(-0.5, 1., 5.),
                                0.1,
                                Vec3::new(0., 0., 1.),
                            ),
                            children: vec![],
                        }],
                    },
                ],
            }],
        }
    }

    pub fn _skew() -> Self {
        Self {
            props: PlantNodeProps::new(Vec3::new(0., 0., 3.), 2.0, Vec3::new(0., 0., 1.)),
            children: vec![PlantNode {
                props: PlantNodeProps::new(Vec3::new(3., 0., 5.), 1.0, Vec3::new(0., 0., 1.)),
                children: vec![],
            }],
        }
    }


    pub fn register_node_info(&self, acc: &mut Vec<NodeInfo>, parent_id: usize) {
        let id = acc.len();
        let parent = acc.get(parent_id);
        acc.push(NodeInfo {
            depth: parent.map(|x| x.depth + 1).unwrap_or(0),
            parent: parent.map(|x| x.id),
            id,
            children: SmallVec::new(),
        });
        for c in &self.children {
            let children_id = acc.len();
            acc[id].children.push(children_id);
            c.register_node_info(acc, id)
        }
    }

    pub fn register_node_properties(&self, acc: &mut Vec<PlantNodeProps>) {
        acc.push(self.props);
        for c in &self.children {
            c.register_node_properties(acc)
        }
    }

    pub fn compute_depth(&self) -> usize {
        1 + self
            .children
            .iter()
            .map(|x| x.compute_depth())
            .max()
            .unwrap_or_default()
    }

    pub fn to_tree(&self) -> TreeSkeleton {
        let mut node_props = Vec::new();
        self.register_node_properties(&mut node_props);
        let mut node_info = Vec::new();
        self.register_node_info(&mut node_info, 0);

        TreeSkeleton {
            node_info,
            node_props,
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub struct BranchSectionPosition {
    // the node being considered
    pub node: usize,
    // the distance we traveled from the node to the leaves
    // it can be negative, in this case we consider the parent
    pub length: f32,
}

impl Add<f32> for BranchSectionPosition {
    type Output = BranchSectionPosition;

    fn add(self, rhs: f32) -> Self::Output {
        Self {
            node: self.node,
            length: self.length + rhs,
        }
    }
}

impl AddAssign<f32> for BranchSectionPosition {
    fn add_assign(&mut self, rhs: f32) {
        self.length += rhs
    }
}

impl BranchSectionPosition {
    pub fn new(node: usize, length: f32) -> Self {
        Self { node, length }
    }
}

#[derive(Component, Clone)]
pub struct TreeSkeleton {
    pub node_props: Vec<PlantNodeProps>,
    pub node_info: Vec<NodeInfo>,
}

impl VisualDebug for TreeSkeleton {
    fn debug(&self, gizmos: &mut Gizmos, debug_flags: crate::DebugFlags) {
        if debug_flags.skeleton {
            for i in 0..self.node_count() {
                let isometry = Isometry3d {
                    translation: self.position(i).into(),
                    rotation: self.orientation(i),
                };
                gizmos.circle(isometry, 1.1 * self.radius(i), Color::srgb(0., 0.8, 0.5));
                for &c in self.children(i) {
                    gizmos.line(
                        self.position(i),
                        self.position(c),
                        Color::srgb(0.1, 0.1, 0.1),
                    );
                }
            }
        }
    }
}


impl TreeSkeleton {
    pub fn root(&self) -> usize {
        0
    }
    pub fn depth(&self, node_id: usize) -> usize {
        self.node_info[node_id].depth
    }
    pub fn position(&self, node_id: usize) -> Vec3 {
        self.node_props[node_id].position
    }
    pub fn radius(&self, node_id: usize) -> f32 {
        self.node_props[node_id].radius
    }
    pub fn orientation(&self, node_id: usize) -> Quat {
        self.node_props[node_id].orientation
    }
    pub fn normal(&self, node_id: usize) -> Vec3 {
        let result = self.orientation(node_id)*Vec3::Z;
        assert!(result != Vec3::ZERO);
        result
    }
    pub fn plane_to_space(&self, node_id: usize, v: Vec2) -> Vec3 {
        self.position(node_id) + self.orientation(node_id)*v.extend(0.)
    }
    pub fn space_to_plane(&self, node_id: usize, v: Vec3) -> Vec2 {
        (self.orientation(node_id).inverse() * (v - self.position(node_id)))
            .truncate()
    }
    pub fn main_children(&self, node_id: usize) -> Option<usize> {
        self.node_info[node_id].children.get(0).copied()
    }
    pub fn parent(&self, node_id: usize) -> Option<usize> {
        self.node_info[node_id].parent
    }

    pub fn branch_length_to_parent(&self, child: usize) -> f32 {
        let parent = self.parent(child).expect("there is no branch under root");
        (self.position(child) - self.position(parent)).length()
    }

    pub fn branch_length_to_main_children(&self, node: usize) -> f32 {
        let child = self
            .main_children(node)
            .expect("a leaf does not have a length");
        (self.position(node) - self.position(child)).length()
    }

    pub(crate) fn node_count(&self) -> usize {
        self.node_props.len()
    }

    pub(crate) fn children(&self, root: usize) -> &[usize] {
        &self.node_info[root].children[..]
    }

    pub fn branch_section_center(&self, pos: BranchSectionPosition) -> Vec3 {
        if pos.length < 0. {
            let parent = self
                .parent(pos.node)
                .expect("there is no branch under root");
            let length = self.branch_length_to_parent(pos.node);
            self.node_props[pos.node]
                .position
                .lerp(self.node_props[parent].position, -pos.length / length)
        } else {
            let length = self.branch_length_to_main_children(pos.node);
            let child = *self.node_info[pos.node]
                .children
                .get(0)
                .expect("there is no branch after this leaf");
            self.node_props[pos.node]
                .position
                .lerp(self.node_props[child].position, pos.length / length)
        }
    }
}

fn sample_random_children_rotations(rng: &mut impl Rng, min_angle: f32, max_angle: f32) -> (Quat, Quat) {
    let a: f32 = rng.gen_range(min_angle..max_angle);
    let b: f32 = rng.gen_range(min_angle..max_angle);
    let c: f32 = rng.gen_range(0f32..2.*std::f32::consts::PI);
    let rot1 = Quat::from_rotation_z(c)
        * Quat::from_rotation_x(a)
        * Quat::from_rotation_z(-c);
    let rot2 = Quat::from_rotation_z(c)
        * Quat::from_rotation_x(-b)
        * Quat::from_rotation_z(-c);
    (rot1, rot2)

}

fn sample_radius(rng: &mut impl Rng, parent_radius: f32) -> f32 {
    rng.gen_range(parent_radius*0.4..parent_radius*0.6)
}

fn sample_size(rng: &mut impl Rng, parent_size: f32) -> f32 {
    parent_size * rng.gen_range(2f32..4f32)
}

fn grow_tree_basic(
    rng: &mut impl Rng,
    root: PlantNodeProps,
    depth: usize,
    min_radius: f32,
    ) -> PlantNode {

    // TODO: why 

    let location = |rotation, radius, size| {
        PlantNodeProps {
            position: root.position + root.orientation*rotation*(size*Vec3::Z),
            orientation: root.orientation * rotation,
            radius,
        }
    };

    let children = if root.radius < min_radius {
        vec![]
    }
    else if rand::random::<f32>() > 1./f32::powf(2., depth as f32) {
        let (rot1, rot2) = sample_random_children_rotations(rng, 0.5, 0.8);
        let (r1, r2) = (sample_radius(rng, root.radius), sample_radius(rng, root.radius));
        let (s1, s2) = (sample_size(rng, root.radius), sample_size(rng, root.radius));
        vec![
            grow_tree_basic(rng, location(rot1, r1, s1), depth+1, min_radius),
            grow_tree_basic(rng, location(rot2, r2, s2), depth+1, min_radius)
        ]
    }
    else {
        let (rot, _) = sample_random_children_rotations(rng, 0., 0.2);
        let r = root.radius;
        vec![
            grow_tree_basic(rng, location(rot, r, 3.*root.radius), depth+1, min_radius)
        ]
    };

    PlantNode {
        props: root,
        children,
    }

}
