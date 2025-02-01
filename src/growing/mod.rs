use bevy::math::{Isometry3d, Quat, Vec3};
use bevy_gizmos::prelude::Gizmos;
use bevy::color::Color;

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

    pub fn register_node_postfix(&self, count: &mut usize, acc: &mut Vec<usize>) {
        let my_id = *count;
        *count+=1;
        for c in &self.children {
            c.register_node_postfix(count, acc)
        }
        acc.push(my_id);
    }

    pub fn register_node_children(&self, acc: &mut Vec<Vec<usize>>) {
        let self_id = acc.len();
        acc.push(Vec::new());
        for c in &self.children {
            let children_id = acc.len();
            acc[self_id].push(children_id);
            c.register_node_children(acc)
        }
    }

    pub fn register_node_properties(&self, acc: &mut Vec<PlantNodeProps>) {
        acc.push(self.props);
        for c in &self.children {
            c.register_node_properties(acc)
        }
    }

    pub fn register_parents(&self, acc: &mut Vec<usize>, parent: usize) {
        let my_id = acc.len();
        acc.push(parent);
        for c in &self.children {
            c.register_parents(acc, my_id)
        }
    }

    pub fn register_leaves(&self, count: &mut usize, acc: &mut Vec<usize>) {
        if self.children.len() == 0 {
            acc.push(*count);
        }
        *count += 1;
        for c in &self.children {
            c.register_leaves(count, acc);
        }
    }

    pub fn compute_depth(&self) -> usize {
        1 + self.children
            .iter()
            .map(|x| x.compute_depth())
            .max()
            .unwrap_or_default()
    }

    pub fn register_depths(&self, acc: &mut Vec<usize>, start_depth: usize) {
        acc.push(start_depth);
        for c in &self.children {
            c.register_depths(acc, start_depth+1)
        }
    }

    pub fn count_node(&self) -> usize {
        let n_node_child: usize = 
            self.children.iter()
            .map(|x| x.count_node())
            .sum();

        n_node_child + 1
    }
}
