use std::ops::Deref;

use bevy::asset::RenderAssetUsages;
use bevy::color::ColorToComponents;
use bevy::math::Vec3;
use bevy::prelude::{Color, Component, Mesh};
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy_gizmos::prelude::Gizmos;
use rand::prelude::*;

use crate::VisualDebug;

#[derive(Default, Component)]
pub struct GeometryData {
    contours: Vec<Vec<Vec3>>,
    debug_points: Vec<(Vec3, Color)>,
    triangles: Vec<usize>,
    colors: Vec<Color>,
    points: Vec<Vec3>,
}

impl GeometryData {
    pub fn to_mesh(&self) -> Mesh {
        let mut result = Mesh::new(
            PrimitiveTopology::TriangleList,
            RenderAssetUsages::default(),
        )
        .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, self.points.clone())
        .with_inserted_attribute(Mesh::ATTRIBUTE_COLOR, self.colors
            .iter()
            .map(|a| a.to_linear().to_f32_array())
            .collect::<Vec<_>>())
        .with_inserted_indices(Indices::U32(
            self.triangles.iter().map(|x| *x as u32).collect(),
        ));
        result.compute_smooth_normals();
        result
    }

    pub fn register_points(
        &mut self,
        points: &[Vec3],
    ) -> Vec<usize> {
        let i0 = self.points.len();
        let n = points.len();
        self.points.extend(points);

        for _ in 0..n {
            let a: f32 = rand::random();
            let color = Color::srgb(0.35, 0.2, 0.1+0.01*a);
            self.colors.push(color);
        }

        (i0..i0 + n).into_iter().collect()
    }

    pub fn point(&self, i: impl Deref<Target=usize>) -> Vec3 {
        self.points[*i.deref()]
    }

    pub fn register_triangles(&mut self, triangles: &[usize]) {
        self.triangles.extend(triangles)
    }

    pub fn mark_debug(&mut self, id: impl Deref<Target=usize>, color: Color) {
        self.debug_points.push(
            (self.point(id), color)
        );
    }

    pub fn add_debug(&mut self, pos: Vec3, color: Color) {
        self.debug_points.push(
            (pos, color)
        );
    }

    pub fn points(&self) -> &[Vec3] {
        &self.points
    }

    pub fn add_contour(&mut self, points: &[Vec3]) {
        self.contours.push(points.to_vec())
    }

}

impl VisualDebug for GeometryData {
    fn debug(&self, gizmos: &mut Gizmos, debug_flags: crate::DebugFlags) {
        if debug_flags.triangles {
            for i in 0..self.triangles.len() / 3 {
                let (ia, ib, ic) = (
                    self.triangles[3 * i],
                    self.triangles[3 * i + 1],
                    self.triangles[3 * i + 2],
                );
                let (pa, pb, pc) = (
                    self.points[ia],
                    self.points[ib],
                    self.points[ic],
                );
                let color = Color::srgb(0., 0.4, 0.);
                gizmos.line(pa, pb, color);
                gizmos.line(pb, pc, color);
                gizmos.line(pc, pa, color);
            }
        }

        // TODO: store normals
        //if debug_flags.normals {
        //    let mesh = self.to_mesh();
        //    let normals = mesh.attribute(Mesh::ATTRIBUTE_NORMAL).unwrap()
        //        .as_float3().unwrap();
        //    for i in 0..normals.len() {
        //        let p = self.normals[i];
        //        let normal = Vec3::from(normals[i]);
        //        let color = Color::srgb(1., 0.0, 0.);
        //        gizmos.line(p, p + 0.3 * normal, color);
        //    }
        //}
        if debug_flags.contours {
            let mut rng = StdRng::seed_from_u64(42);
            for c in self.contours.iter() {
                let t = rng.gen();
                let n = c.len();
                for i in 0..n {
                    let r = i as f32 / n as f32;
                    let color = Color::srgb(t, 0.3 + 0.2 * r, 0.3 + 0.2 * r);
                    let pos1 = c[i];
                    let pos2 = c[(i + 1) % n];
                    gizmos.line(pos1, pos2, color);
                }
            }
        }
        if debug_flags.other {
            for (p, c) in self.debug_points.iter() {
                gizmos.cross(*p, 0.2, *c);
            }
        }
    }
}

