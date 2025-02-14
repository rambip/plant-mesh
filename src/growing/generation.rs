use bevy::math::{Quat, Vec3};
use rand::{distributions::Uniform, prelude::Distribution, Rng};

use super::{PlantNode, PlantNodeProps};

struct ChildrenBranchRotations {
    max_turn_angle: f32,
    min_split_angle: f32,
}

impl Distribution<(Quat, Quat)> for ChildrenBranchRotations {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> (Quat, Quat) {
        let a: f32 = rng.gen_range(0f32..0.5*self.max_turn_angle);
        let b: f32 = rng.gen_range((self.min_split_angle-a)..self.max_turn_angle);
        let c: f32 = rng.gen_range(0f32..2. * std::f32::consts::PI);
        let rot1 = Quat::from_rotation_z(c) * Quat::from_rotation_x(a) * Quat::from_rotation_z(-c);
        let rot2 = Quat::from_rotation_z(c) * Quat::from_rotation_x(-b) * Quat::from_rotation_z(-c);
        (rot1, rot2)
    }
}

struct ChildBranchRotation {
    max_angle: f32
}

impl Distribution<Quat> for ChildBranchRotation {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Quat {
        let a: f32 = rng.gen_range(0f32..self.max_angle);
        let c: f32 = rng.gen_range(0f32..2. * std::f32::consts::PI);
        Quat::from_rotation_z(c) * Quat::from_rotation_x(a) * Quat::from_rotation_z(-c)
    }
}

struct ChildrenBranchRadiuses {
    radius: f32,
    main_children_factor: f32,
    secondary_children_factor: f32,
    variance: f32,
}

impl Distribution<(f32, f32)> for ChildrenBranchRadiuses {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> (f32, f32) {
        (
            self.radius*rng.gen_range(
                self.main_children_factor-self.variance
                ..self.main_children_factor+self.variance
                ),
            self.radius*rng.gen_range(
                self.secondary_children_factor-self.variance
                ..self.secondary_children_factor+self.variance
                ),
        )
    }
}

struct ChildrenBranchSize {
    radius: f32,
}

impl Distribution<f32> for ChildrenBranchSize {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f32 {
        self.radius * rng.gen_range(3f32..5f32)
    }
}

#[derive(Copy, Clone, serde::Serialize, serde::Deserialize)]
pub struct GrowConfig {
    pub min_radius: f32,
    pub max_turn_angle: f32,
    pub min_split_angle: f32,
    pub max_zigzag_angle: f32,
    pub radius_to_branch_ratio: f32,
    pub branch_variance: f32,
    pub main_children_radius_factor: f32,
    pub secondary_children_radius_factor: f32,
    pub radius_variance: f32

}

impl GrowConfig {
    fn branch_rotation_2_child(&self) -> ChildrenBranchRotations {
        ChildrenBranchRotations {
            min_split_angle: self.min_split_angle,
            max_turn_angle: self.max_turn_angle,
        }
    }
    fn branch_rotation_1_child(&self) -> ChildBranchRotation {
        ChildBranchRotation {
            max_angle: self.max_zigzag_angle,
        }
    }
    fn branch_size(&self, radius: f32) -> Uniform<f32> {
        Uniform::from(
            radius*(self.radius_to_branch_ratio - self.branch_variance)..
            radius*(self.radius_to_branch_ratio + self.branch_variance)
        )
    }
    fn branch_radius(&self, radius: f32) -> ChildrenBranchRadiuses {
        ChildrenBranchRadiuses {
            radius,
            main_children_factor: self.main_children_radius_factor,
            secondary_children_factor: self.secondary_children_radius_factor,
            variance: self.radius_variance
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
        orientation: root.orientation * rotation,
        radius,
    };

    let children = if root.radius < config.min_radius {
        vec![]
    } else if rng.gen_range(0f32..1f32) > 1. / f32::powf(2., depth as f32) {
        let (rot1, rot2) = rng.sample(config.branch_rotation_2_child());
        let (r1, r2) = rng.sample(config.branch_radius(root.radius));
        let (s1, s2) = (
            rng.sample(config.branch_size(root.radius)),
            rng.sample(config.branch_size(root.radius))
        );
        vec![
            grow_tree_basic(config, rng, location(rot1, r1, s1), depth + 1),
            grow_tree_basic(config, rng, location(rot2, r2, s2), depth + 1),
        ]
    } else {
        let rot = rng.sample(config.branch_rotation_1_child());
        let r = root.radius;
        let s = rng.sample(ChildrenBranchSize {
            radius: root.radius,
        });
        vec![grow_tree_basic(
            config,
            rng,
            location(rot, r, s),
            depth + 1,
        )]
    };

    PlantNode {
        props: root,
        children,
    }
}
