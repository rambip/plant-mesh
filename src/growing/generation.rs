use bevy::math::{Quat, Vec3};
use rand::Rng;

use super::{PlantNode, PlantNodeProps};


fn sample_random_children_rotations(rng: &mut impl Rng, min_angle: f32, max_angle: f32) -> (Quat, Quat) {
    let a: f32 = rng.gen_range(min_angle..max_angle);
    let b: f32 = rng.gen_range(min_angle..max_angle);
    let c: f32 = rng.gen_range(0f32..2.*std::f32::consts::PI);
    let rot1 = Quat::from_rotation_z(c)
        * Quat::from_rotation_x(a)
        * Quat::from_rotation_z(-c);
    let rot2 = Quat::from_rotation_z(c)
        * Quat::from_rotation_x(-b)
        * Quat::from_rotation_z(-c);
    (rot1, rot2)

}

fn sample_radius(rng: &mut impl Rng, parent_radius: f32) -> f32 {
    rng.gen_range(parent_radius*0.4..parent_radius*0.7)
}

fn sample_size(rng: &mut impl Rng, parent_size: f32) -> f32 {
    parent_size * rng.gen_range(2f32..4f32)
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
    else if rand::random::<f32>() > 1./f32::powf(2., depth as f32) {
        let (rot1, rot2) = sample_random_children_rotations(rng, 0.5, 0.8);
        let (r1, r2) = (sample_radius(rng, root.radius), sample_radius(rng, root.radius));
        let (s1, s2) = (sample_size(rng, root.radius), sample_size(rng, root.radius));
        vec![
            grow_tree_basic(rng, location(rot1, r1, s1), depth+1, min_radius),
            grow_tree_basic(rng, location(rot2, r2, s2), depth+1, min_radius)
        ]
    }
    else {
        let (rot, _) = sample_random_children_rotations(rng, 0., 0.2);
        let r = root.radius;
        vec![
            grow_tree_basic(rng, location(rot, r, 3.*root.radius), depth+1, min_radius)
        ]
    };

    PlantNode {
        props: root,
        children,
    }

}
