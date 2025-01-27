use bevy::math::{Vec3, Isometry3d};
use bevy::gizmos::gizmos::Gizmos;
use bevy::color::Color;


pub struct PlantNode {
    children: Vec<PlantNode>,
    props: PlantNodeProps,
    id: usize,

}

#[derive(Copy, Clone, Debug, Default)]
pub struct PlantNodeProps {
    pub position: Vec3,
    pub radius: f32,
}

impl PlantNodeProps {
    fn new(position: Vec3, radius: f32) -> Self {
        Self {position, radius }
    }
}

impl PlantNode {
    pub fn demo() -> Self {
        Self {
            id: 0,
            props: PlantNodeProps::new(Vec3::ZERO, 1.1),
            children: vec![
                PlantNode {
                    id: 1,
                    props: PlantNodeProps::new(Vec3::new(0., 0., 2.), 0.9),
                    children: vec![
                        PlantNode {
                            id: 2,
                            props: PlantNodeProps::new(Vec3::new(1., 0., 3.), 0.5),
                            children: vec![
                                PlantNode {
                                    id: 3,
                                    props: PlantNodeProps::new(Vec3::new(1.5, 0., 6.), 0.3),
                                    children: vec![],
                                },
                                PlantNode {
                                    id: 4,
                                    props: PlantNodeProps::new(Vec3::new(2., -2., 5.), 0.3),
                                    children: vec![],
                                }
                            ],
                        },
                        PlantNode {
                            id: 5,
                            props: PlantNodeProps::new(Vec3::new(-1., 0., 4.5), 0.2),
                            children: vec![
                                PlantNode {
                                    id: 6,
                                    props: PlantNodeProps::new(Vec3::new(-0.5, 1., 6.), 0.1),
                                    children: vec![],
                                } ],
                        } ],
                }],
        }
    }

    pub fn _skew() -> Self {
        Self {
            id: 0,
            props: PlantNodeProps::new(Vec3::new(0., 0., 3.), 2.0),
            children: vec![
                PlantNode {
                    id: 1,
                    props: PlantNodeProps::new(Vec3::new(3., 0., 5.), 1.0),
                    children: vec![]
                }
            ],
        }
    }

    pub fn debug(&self, gizmos: &mut Gizmos) {
        gizmos.circle(Isometry3d::from_translation(self.props.position), 1.1*self.props.radius, Color::srgb(0., 0.8, 0.5));
        for c in &self.children {
            c.debug(gizmos)
        }

    }

    pub fn register_node_properties(&self, acc: &mut Vec<PlantNodeProps>) {
        acc[self.id] = self.props;
        for c in &self.children {
            c.register_node_properties(acc)
        }
    }

    pub fn register_parents(&self, acc: &mut Vec<usize>, parent: usize) {
        acc[self.id] = parent;
        for c in &self.children {
            c.register_parents(acc, self.id)
        }
    }

    pub fn register_leaves(&self, acc: &mut Vec<usize>) {
        if self.children.len() == 0 {
            acc.push(self.id);
        }
        for c in &self.children {
            c.register_leaves(acc);
        }
    }

    pub fn register_depths(&self, acc: &mut Vec<usize>, start_depth: usize) {
        acc[self.id] = start_depth;
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
