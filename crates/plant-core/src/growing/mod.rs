use glam::{Quat, Vec2, Vec3};
use smallvec::SmallVec;
use serde::{Serialize, Deserialize};
#[cfg(feature = "bevy")]
use bevy_gizmos::prelude::Gizmos;
#[cfg(feature = "bevy")]
use bevy_color::Color;
#[cfg(feature = "bevy")]
use bevy_math::Isometry3d;

#[derive(Copy, Clone, Debug, Default, Serialize, Deserialize)]
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
}

#[derive(Clone, Serialize, Deserialize)]
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
}

impl crate::VisualDebug for TreeSkeleton {
    type Flags = bool;
    #[cfg(feature = "bevy")]
    fn debug(&self, gizmos: &mut Gizmos, debug_flags: Self::Flags) {
        if debug_flags {
            for i in 0..self.node_count() {
                let isometry = Isometry3d {
                    translation: self.position(i).into(),
                    rotation: self.orientation(i),
                };
                gizmos.circle(
                    isometry,
                    1.1 * self.radius(i),
                    Color::srgb(0., 0.8, 0.5),
                );
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

pub mod generation;
