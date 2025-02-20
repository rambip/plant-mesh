use bevy::math::{Quat, Vec3};
use rand::{prelude::Distribution, Rng};

use super::{PlantNode, PlantNodeProps};

#[derive(Debug)]
struct ChildrenBranches {
    max_turn_angle: f32,
    length: f32,
    variability: f32,
    min_distance: f32
}

impl Distribution<((Quat, f32), (Quat, f32))> for ChildrenBranches {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> ((Quat, f32), (Quat, f32)) {
        for _ in 0..100 {
            let a: f32 = rng.gen_range(0f32..self.max_turn_angle);
            let b: f32 = rng.gen_range(0f32..self.max_turn_angle);
            let c: f32 = rng.gen_range(0f32..2. * std::f32::consts::PI);
            let rot1 = Quat::from_rotation_z(c) * Quat::from_rotation_x(a) * Quat::from_rotation_z(-c);
            let rot2 = Quat::from_rotation_z(c) * Quat::from_rotation_x(-b) * Quat::from_rotation_z(-c);
            let branch_range = (self.length - self.variability) ..(self.length + self.variability);
            let l1 = rng.gen_range(branch_range.clone());
            let l2 = rng.gen_range(branch_range.clone());
            let pos1 = rot1 * l1;
            let pos2 = rot2 * l2;
            if (pos1 - pos2).length() > self.min_distance {
                return ((rot1, l1), (rot2, l2))
            }
        }
        panic!("not able to generate branch children with these parameters: {self:?}");
    }
}

struct ChildBranch {
    max_angle: f32,
    length: f32,
    variability: f32,
}

impl Distribution<(Quat, f32)> for ChildBranch {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> (Quat, f32) {
        let a: f32 = rng.gen_range(0f32..self.max_angle);
        let c: f32 = rng.gen_range(0f32..2. * std::f32::consts::PI);
        let rot = Quat::from_rotation_z(c) * Quat::from_rotation_x(a) * Quat::from_rotation_z(-c);
        let branch_range = (self.length - self.variability) ..(self.length + self.variability);
        let l = rng.gen_range(branch_range.clone());
        (rot, l)
    }
}

struct ChildrenBranchRadiuses {
    parent_radius: f32,
    main_children_factor: f32,
    secondary_children_factor: f32,
    variance: f32,
}

impl Distribution<(f32, f32)> for ChildrenBranchRadiuses {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> (f32, f32) {
        (
            self.parent_radius
                * rng.gen_range(
                    self.main_children_factor - self.variance
                        ..self.main_children_factor + self.variance,
                ),
            self.parent_radius
                * rng.gen_range(
                    self.secondary_children_factor - self.variance
                        ..self.secondary_children_factor + self.variance,
                ),
        )
    }
}

#[derive(Copy, Clone, serde::Serialize, serde::Deserialize)]
pub struct GrowConfig {
    pub min_radius: f32,
    pub max_turn_angle: f32,
    pub max_zigzag_angle: f32,
    pub radius_to_branch_ratio: f32,
    pub branch_variance: f32,
    pub main_children_radius_factor: f32,
    pub secondary_children_radius_factor: f32,
    pub radius_variance: f32,
    pub birth_coefficient: f32,
    pub birth_power: f32,
    pub up_attraction_factor: f32,
    pub base_radius: f32,
}

impl GrowConfig {
    fn branch_rotation_2_child(&self, radius: f32) -> ChildrenBranches {
        ChildrenBranches {
            max_turn_angle: self.max_turn_angle,
            length: radius*self.radius_to_branch_ratio,
            variability: radius*self.branch_variance,
            min_distance: radius
        }
    }
    fn branch_rotation_1_child(&self, radius: f32) -> ChildBranch {
        ChildBranch {
            max_angle: self.max_zigzag_angle,
            length: radius*self.radius_to_branch_ratio,
            variability: radius*self.branch_variance,
        }
    }
    fn branch_radius(&self, radius: f32) -> ChildrenBranchRadiuses {
        ChildrenBranchRadiuses {
            parent_radius: radius,
            main_children_factor: self.main_children_radius_factor,
            secondary_children_factor: self.secondary_children_radius_factor,
            variance: self.radius_variance,
        }
    }
}

pub fn grow_tree_basic(
    config: &GrowConfig,
    rng: &mut impl Rng,
    root: PlantNodeProps,
    depth: usize,
) -> PlantNode {
    let location = |rotation, radius, size| PlantNodeProps {
        position: root.position + root.orientation * rotation * (size * Vec3::Z),
        orientation: root.orientation.lerp(Quat::IDENTITY, config.up_attraction_factor) * rotation,
        radius,
    };

    let children = if root.radius < config.min_radius {
        vec![]
    } else if rng.gen_range(0f32..1f32) > 1. / config.birth_coefficient / (depth as f32).powf(config.birth_power) {
        let ((rot1, s1), (rot2, s2)) = rng.sample(config.branch_rotation_2_child(root.radius));
        let (r1, r2) = rng.sample(config.branch_radius(root.radius));
        vec![
            grow_tree_basic(config, rng, location(rot1, r1, s1), depth + 1),
            grow_tree_basic(config, rng, location(rot2, r2, s2), depth + 1),
        ]
    } else {
        let (rot, s) = rng.sample(config.branch_rotation_1_child(root.radius));
        let r = root.radius;
        vec![grow_tree_basic(config, rng, location(rot, r, s), depth + 1)]
    };

    PlantNode {
        props: root,
        children,
    }
}
