use crate::TreePipelinePhase;
use std::ops::{Add, AddAssign};

use bevy::color::Color;
use bevy::math::{Isometry3d, Quat, Vec2, Vec3};
use bevy::prelude::Component;
use bevy_gizmos::prelude::Gizmos;
use smallvec::SmallVec;

use crate::VisualDebug;

#[derive(Clone)]
pub struct NodeInfo {
    pub depth: usize,
    pub parent: Option<usize>,
    // id in prefix traversal order of tree
    pub id: usize,
    // TODO: bigger children first
    pub children: SmallVec<[usize; 2]>,
}

mod generation;


#[derive(Copy, Clone, Debug, Default)]
pub struct PlantNodeProps {
    pub position: Vec3,
    pub radius: f32,
    pub orientation: Quat,
}

impl PlantNodeProps {
    fn new(position: Vec3, radius: f32, orientation: Vec3) -> Self {
        Self {
            position,
            radius,
            orientation: Quat::from_rotation_arc(Vec3::Z, orientation.normalize()),
        }
    }
}

#[derive(Debug)]
pub struct PlantNode {
    children: Vec<PlantNode>,
    props: PlantNodeProps,
}

pub struct Seed;

pub struct GrowConfig { }

#[derive(Default, Component)]
pub struct TreeSkeletonDebugData {
    copy: TreeSkeleton,
}

impl TreePipelinePhase for PlantNode {
    type Previous = Seed;
    type Config = GrowConfig;
    type DebugCache = ();
    fn generate_from(
        _: Self::Previous, 
        config: &Self::Config, 
        _: &mut Self::DebugCache, 
        mut rng: rand::prelude::StdRng) -> Self {
        let root = PlantNodeProps {
            radius: 0.5,
            orientation: Quat::default(),
            position: Vec3::ZERO,
        };
        generation::grow_tree_basic(&mut rng, root, 0, 0.05)
    }
}

impl TreePipelinePhase for TreeSkeleton {
    type Previous = PlantNode;
    type Config = ();
    type DebugCache = TreeSkeletonDebugData;
    fn generate_from(prev: Self::Previous, 
        _: &Self::Config, 
        cache: &mut Self::DebugCache, 
        _rng: rand::prelude::StdRng) -> Self {

        let mut node_props = Vec::new();
        prev.register_node_properties(&mut node_props);
        let mut node_info = Vec::new();
        prev.register_node_info(&mut node_info, 0);

        let result = TreeSkeleton {
            node_info,
            node_props,
        };
        cache.copy = result.clone();
        result
    }
}


impl PlantNode {
    pub fn _demo() -> Self {
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

#[derive(Clone, Default)]
pub struct TreeSkeleton {
    pub node_props: Vec<PlantNodeProps>,
    pub node_info: Vec<NodeInfo>,
}

impl VisualDebug for TreeSkeletonDebugData {
    fn debug<R: rand::Rng + Clone>(&self, 
        gizmos: &mut Gizmos,
        _rng: R,
        debug_flags: crate::DebugFlags
) {
        if debug_flags.skeleton {
            for i in 0..self.copy.node_count() {
                let isometry = Isometry3d {
                    translation: self.copy.position(i).into(),
                    rotation: self.copy.orientation(i),
                };
                gizmos.circle(isometry, 1.1 * self.copy.radius(i), Color::srgb(0., 0.8, 0.5));
                for &c in self.copy.children(i) {
                    gizmos.line(
                        self.copy.position(i),
                        self.copy.position(c),
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

