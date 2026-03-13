#[cfg(feature = "bevy")]
use bevy::prelude::Component;

use crate::VisualDebug;
use glam::{Quat, Vec2, Vec3};
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;

#[derive(Copy, Clone, Debug, Default, Serialize, Deserialize)]
#[cfg_attr(feature = "bevy", derive(Component))]
pub struct PlantNodeProps {
    pub position: Vec3,
    pub radius: f32,
    pub orientation: Quat,
}

impl PlantNodeProps {
    pub fn new(position: Vec3, radius: f32, orientation: Vec3) -> Self {
        Self {
            position,
            radius,
            orientation: Quat::from_rotation_arc(Vec3::Z, orientation.normalize()),
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct NodeInfo {
    pub depth: usize,
    pub parent: Option<usize>,
    pub id: usize,
    pub children: SmallVec<[usize; 2]>,
}

#[derive(Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "bevy", derive(Component))]
pub struct PlantNode {
    pub children: Vec<PlantNode>,
    pub props: PlantNodeProps,
}

impl PlantNode {
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

    pub fn grow_skeleton(&self) -> TreeSkeleton {
        let mut node_props = Vec::new();
        self.register_node_properties(&mut node_props);
        let mut node_info = Vec::new();
        self.register_node_info(&mut node_info, 0);
        TreeSkeleton {
            node_info,
            node_props,
        }
    }

    pub fn grow_skeleton_debug<F>(&self, mut callback: F) -> TreeSkeleton
    where
        F: FnMut(crate::DebugGeometry),
    {
        let skeleton = self.grow_skeleton();
        let debug_data = TreeSkeletonDebugData {
            copy: skeleton.clone(),
        };
        callback(debug_data.debug_data());
        skeleton
    }

    // TODO: remove
    pub fn grow_skeleton_with_debug(&self) -> (TreeSkeleton, crate::DebugGeometry) {
        let skeleton = self.grow_skeleton();
        let debug_data = TreeSkeletonDebugData {
            copy: skeleton.clone(),
        };
        (skeleton, debug_data.debug_data())
    }
}

#[derive(Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "bevy", derive(Component))]
pub struct TreeSkeleton {
    pub node_props: Vec<PlantNodeProps>,
    pub node_info: Vec<NodeInfo>,
}

impl TreeSkeleton {
    pub fn root(&self) -> usize {
        0
    }
    pub fn node_count(&self) -> usize {
        self.node_props.len()
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
    pub fn children(&self, root: usize) -> &[usize] {
        &self.node_info[root].children[..]
    }
    pub fn parent(&self, node_id: usize) -> Option<usize> {
        self.node_info[node_id].parent
    }
    pub fn normal(&self, node_id: usize) -> Vec3 {
        self.orientation(node_id) * Vec3::Z
    }
    pub fn plane_to_space(&self, node_id: usize, v: Vec2) -> Vec3 {
        self.position(node_id) + self.orientation(node_id) * v.extend(0.)
    }
    pub fn space_to_plane(&self, node_id: usize, v: Vec3) -> Vec2 {
        (self.orientation(node_id).inverse() * (v - self.position(node_id))).truncate()
    }

    pub fn grow_strands(
        &self,
        config: &crate::StrandsConfig,
        rng: rand::rngs::StdRng,
    ) -> crate::VolumetricTree {
        let mut cache = crate::TrajectoryBuilder::new(rng);
        cache.clear_for_tree(self);
        cache.compute_trajectories(self, self.root(), config);
        super::meshing::VolumetricTree {
            trajectories: cache.trajectories.clone(),
            particles_per_node: cache.particles_per_node.clone(),
            tree: self.clone(),
        }
    }

    pub fn grow_strands_debug<F>(
        &self,
        config: &crate::StrandsConfig,
        rng: rand::rngs::StdRng,
        mut callback: F,
    ) -> crate::VolumetricTree
    where
        F: FnMut(crate::DebugGeometry),
    {
        let rng_for_debug = rng.clone();
        let volumetric = self.grow_strands(config, rng);
        let debug_data = crate::TrajectoryBuilder {
            particles_per_node: volumetric.particles_per_node.clone(),
            trajectories: volumetric.trajectories.clone(),
            rng: rng_for_debug,
        };
        callback(debug_data.debug_data());
        volumetric
    }
}

#[cfg_attr(feature = "bevy", derive(bevy::prelude::Component))]
pub struct TreeSkeletonDebugData {
    pub copy: TreeSkeleton,
}

impl Default for TreeSkeletonDebugData {
    fn default() -> Self {
        Self::new()
    }
}

impl TreeSkeletonDebugData {
    pub fn new() -> Self {
        Self {
            copy: TreeSkeleton {
                node_props: Default::default(),
                node_info: Default::default(),
            },
        }
    }
}

impl crate::VisualDebug for TreeSkeletonDebugData {
    fn debug_data(&self) -> crate::DebugGeometry {
        let mut out = crate::DebugGeometry::new();
        let sk = &self.copy;
        for i in 0..sk.node_count() {
            let circle = crate::Circle {
                position: sk.position(i),
                orientation: sk.orientation(i),
                radius: 1.1 * sk.radius(i),
            };
            out.circles
                .push((circle, crate::DebugColor::rgb(0., 0.8, 0.5)));
            for &c in sk.children(i) {
                out.lines.push((
                    sk.position(i),
                    sk.position(c),
                    crate::DebugColor::rgb(0.1, 0.1, 0.1),
                ));
            }
        }
        out
    }
}

pub mod generation;
