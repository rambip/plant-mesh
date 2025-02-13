use bevy::math::{Quat, Vec3};
use rand::{prelude::Distribution, Rng};

use super::{PlantNode, PlantNodeProps};

struct ChildrenBranchRotations { }

impl Distribution<(Quat, Quat)> for ChildrenBranchRotations {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> (Quat, Quat) {
        let a: f32 = rng.gen_range(0f32..0.5);
        let b: f32 = rng.gen_range(0.3f32..1f32);
        let c: f32 = rng.gen_range(0f32..2.*std::f32::consts::PI);
        let rot1 = Quat::from_rotation_z(c)
            * Quat::from_rotation_x(a)
            * Quat::from_rotation_z(-c);
        let rot2 = Quat::from_rotation_z(c)
            * Quat::from_rotation_x(-b)
            * Quat::from_rotation_z(-c);
        (rot1, rot2)
    }
}

struct ChildBranchRotation { }

impl Distribution<Quat> for ChildBranchRotation {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Quat {
        let a: f32 = rng.gen_range(0f32..0.2);
        let c: f32 = rng.gen_range(0f32..2.*std::f32::consts::PI);
        Quat::from_rotation_z(c)
            * Quat::from_rotation_x(a)
            * Quat::from_rotation_z(-c)
    }
}

struct ChildrenBranchRadiuses {radius: f32}

impl Distribution<(f32, f32)> for ChildrenBranchRadiuses {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> (f32, f32) {
        (
            rng.gen_range(self.radius*0.6..self.radius*0.8),
            rng.gen_range(self.radius*0.4..self.radius*0.7),
        )
    }
}

struct ChildrenBranchSize {radius: f32}

impl Distribution<f32> for ChildrenBranchSize {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f32 {
        self.radius * rng.gen_range(3f32..5f32)
    }
}

pub fn grow_tree_basic(
    rng: &mut impl Rng,
    root: PlantNodeProps,
    depth: usize,
    min_radius: f32,
    ) -> PlantNode {

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
    else if rng.gen_range(0f32..1f32) > 1./f32::powf(2., depth as f32) {
        let (rot1, rot2) = rng.sample(ChildrenBranchRotations {
        });
        let (r1, r2) = rng.sample(ChildrenBranchRadiuses {radius: root.radius});
        let (s1, s2) = 
            (rng.sample(ChildrenBranchSize{radius: root.radius}),
             rng.sample(ChildrenBranchSize{radius: root.radius})
            );
        vec![
            grow_tree_basic(rng, location(rot1, r1, s1), depth+1, min_radius),
            grow_tree_basic(rng, location(rot2, r2, s2), depth+1, min_radius)
        ]
    }
    else {
        let rot = rng.sample(ChildBranchRotation {});
        let r = root.radius;
        let s = rng.sample(ChildrenBranchSize {radius: root.radius});
        vec![
            grow_tree_basic(rng, location(rot, r, s), depth+1, min_radius)
        ]
    };

    PlantNode {
        props: root,
        children,
    }

}
