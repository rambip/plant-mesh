use crate::TreePipelinePhase;

use bevy::color::Color;
use bevy::math::{Isometry3d, Quat, Vec2, Vec3};
use bevy::prelude::Component;
use bevy_gizmos::prelude::Gizmos;
use rand::rngs::StdRng;
use smallvec::SmallVec;

use crate::VisualDebug;

#[derive(Copy, Clone, Debug, Default)]
pub struct PlantNodeProps {
    pub position: Vec3,
    pub radius: f32,
    pub orientation: Quat,
}

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
pub use generation::GrowConfig;

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

#[derive(Component)]
pub struct TreeSkeletonDebugData {
    copy: TreeSkeleton,
}

impl TreePipelinePhase for PlantNode {
    type Previous = Seed;
    type Config = GrowConfig;
    type Builder = StdRng;
    fn generate_from(_: Self::Previous, config: &Self::Config, rng: &mut Self::Builder) -> Self {
        let root = PlantNodeProps {
            radius: 0.5,
            orientation: Quat::default(),
            position: Vec3::ZERO,
        };
        generation::grow_tree_basic(config, rng, root, 0)
    }
}

impl TreePipelinePhase for TreeSkeleton {
    type Previous = PlantNode;
    type Config = ();
    type Builder = TreeSkeletonDebugData;
    fn generate_from(prev: Self::Previous, _: &Self::Config, cache: &mut Self::Builder) -> Self {
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


#[derive(Clone)]
pub struct TreeSkeleton {
    pub node_props: Vec<PlantNodeProps>,
    pub node_info: Vec<NodeInfo>,
}

impl From<StdRng> for TreeSkeletonDebugData {
    fn from(_: StdRng) -> Self {
        Self {
            copy: TreeSkeleton {
                node_props: vec![],
                node_info: vec![],
            },
        }
    }
}

impl VisualDebug for TreeSkeletonDebugData {
    fn debug(&self, gizmos: &mut Gizmos, debug_flags: crate::DebugFlags) {
        if debug_flags.skeleton {
            for i in 0..self.copy.node_count() {
                let isometry = Isometry3d {
                    translation: self.copy.position(i).into(),
                    rotation: self.copy.orientation(i),
                };
                gizmos.circle(
                    isometry,
                    1.1 * self.copy.radius(i),
                    Color::srgb(0., 0.8, 0.5),
                );
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
        let result = self.orientation(node_id) * Vec3::Z;
        assert!(result != Vec3::ZERO);
        result
    }
    pub fn plane_to_space(&self, node_id: usize, v: Vec2) -> Vec3 {
        self.position(node_id) + self.orientation(node_id) * v.extend(0.)
    }
    pub fn space_to_plane(&self, node_id: usize, v: Vec3) -> Vec2 {
        (self.orientation(node_id).inverse() * (v - self.position(node_id))).truncate()
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

    pub fn average_branch_length_to_children(&self, node: usize) -> f32 {
        let sum: f32 = self.children(node)
            .iter()
            .map(|&child| (self.position(node) - self.position(child)).length())
            .sum();
        sum / (self.children(node).len() as f32)
    }
    pub fn min_branch_length_to_children(&self, node: usize) -> f32 {
        self.children(node)
            .iter()
            .map(|&child| (self.position(node) - self.position(child)).length())
            .reduce(f32::min)
            .unwrap()
    }

    pub(crate) fn node_count(&self) -> usize {
        self.node_props.len()
    }

    pub(crate) fn children(&self, root: usize) -> &[usize] {
        &self.node_info[root].children[..]
    }
}
