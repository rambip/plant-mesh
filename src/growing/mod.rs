use bevy::math::{Isometry3d, Quat, Vec3};
use bevy_gizmos::prelude::Gizmos;
use bevy::color::Color;
use smallvec::SmallVec;

#[derive(Clone)]
pub struct NodeInfo {
    pub depth: usize,
    pub parent: Option<usize>,
    // id in prefix traversal order of tree
    pub id: usize,
    pub children: SmallVec<[usize; 2]>,
}

//mod simple_generation;

pub struct PlantNode {
    children: Vec<PlantNode>,
    props: PlantNodeProps,
}

#[derive(Copy, Clone, Debug, Default)]
pub struct PlantNodeProps {
    pub position: Vec3,
    pub radius: f32,
    pub orientation: Vec3,
}

impl PlantNodeProps {
    fn new(position: Vec3, radius: f32, orientation: Vec3) -> Self {
        Self {position, radius, orientation: orientation.normalize()}
    }
}

impl PlantNode {
    pub fn demo() -> Self {
        Self {
            props: PlantNodeProps::new(Vec3::new(0., 0., -2.), 1.1, Vec3::new(0., 0., 1.)),
            children: vec![
                PlantNode {
                    props: PlantNodeProps::new(Vec3::new(0., 0., 2.), 0.9, Vec3::new(0., 0., 1.)),
                    children: vec![
                        PlantNode {
                            props: PlantNodeProps::new(Vec3::new(1., 0., 3.), 0.5, Vec3::new(0., 0., 1.)),
                            children: vec![
                                PlantNode {
                                    props: PlantNodeProps::new(Vec3::new(1.5, 0., 6.), 0.3, Vec3::new(0., 0., 1.)),
                                    children: vec![],
                                },
                                PlantNode {
                                    props: PlantNodeProps::new(Vec3::new(2., -2., 5.), 0.3, Vec3::new(0., -2., 1.)),
                                    children: vec![],
                                }
                            ],
                        },
                        PlantNode {
                            props: PlantNodeProps::new(Vec3::new(-1., 0., 4.5), 0.2, Vec3::new(0., 1., 1.)),
                            children: vec![
                                PlantNode {
                                    props: PlantNodeProps::new(Vec3::new(-0.5, 1., 6.), 0.1, Vec3::new(0., 0., 1.)),
                                    children: vec![],
                                } ],
                        } ],
                }],
        }
    }

    pub fn _skew() -> Self {
        Self {
            props: PlantNodeProps::new(Vec3::new(0., 0., 3.), 2.0, Vec3::new(0., 0., 1.)),
            children: vec![
                PlantNode {
                    props: PlantNodeProps::new(Vec3::new(3., 0., 5.), 1.0, Vec3::new(0., 0., 1.)),
                    children: vec![]
                }
            ],
        }
    }

    pub fn debug(&self, gizmos: &mut Gizmos) {
        let isometry = Isometry3d {
            translation: self.props.position.into(),
            rotation: Quat::from_rotation_arc(Vec3::Z, self.props.orientation)
        };
        gizmos.circle(isometry, 1.1*self.props.radius, Color::srgb(0., 0.8, 0.5));
        for c in &self.children {
            gizmos.line(self.props.position, c.props.position, Color::srgb(0.1, 0.1, 0.1));
            c.debug(gizmos)
        }

    }

    pub fn register_node_info(&self, acc: &mut Vec<NodeInfo>, parent_id: usize) {
        let id = acc.len();
        let parent = acc.get(parent_id);
        acc.push(NodeInfo {
            depth: parent.map(|x| x.depth+1).unwrap_or(0),
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
        1 + self.children
            .iter()
            .map(|x| x.compute_depth())
            .max()
            .unwrap_or_default()
    }
}
