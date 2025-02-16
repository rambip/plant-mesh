use std::ops::Deref;

use bevy::asset::RenderAssetUsages;
use bevy::color::ColorToComponents;
use bevy::math::Vec3;
use bevy::prelude::{Color, Component, Mesh};
use bevy_gizmos::prelude::Gizmos;
use bevy_render::mesh::{Indices, PrimitiveTopology};
use rand::rngs::StdRng;
use rand::Rng;

use crate::VisualDebug;

#[derive(Component)]
pub struct GeometryData {
    contours: Vec<Vec<usize>>,
    debug_points: Vec<(Vec3, Color)>,
    triangles: Vec<usize>,
    colors: Vec<Color>,
    points: Vec<Vec3>,
    pub rng: StdRng,
}

impl GeometryData {
    pub fn to_mesh(&self) -> Mesh {
        let mut result = Mesh::new(
            PrimitiveTopology::TriangleList,
            RenderAssetUsages::default(),
        )
        .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, self.points.clone())
        .with_inserted_attribute(
            Mesh::ATTRIBUTE_COLOR,
            self.colors
                .iter()
                .map(|a| a.to_linear().to_f32_array())
                .collect::<Vec<_>>(),
        )
        .with_inserted_indices(Indices::U32(
            self.triangles.iter().map(|x| *x as u32).collect(),
        ));
        result.compute_smooth_normals();
        result
    }

    pub fn register_points_trunk(&mut self, points: impl IntoIterator<Item=Vec3>) -> Vec<usize> {
        let mut result = Vec::new();

        for p in points {
            let blue: f32 = self.rng.gen_range(0.1f32..0.13);
            let color = Color::srgb(0.35, 0.2, blue);
            self.colors.push(color);
            result.push(self.points.len());
            self.points.push(p);
        }

        result
    }
    pub fn register_points_leaf(&mut self, points: impl IntoIterator<Item=Vec3>) -> Vec<usize> {
        let mut result = Vec::new();

        for p in points {
            let color = Color::srgb(0.1, 0.8, 0.3);
            self.colors.push(color);
            result.push(self.points.len());
            self.points.push(p);
        }

        result
    }

    pub fn point(&self, i: impl Deref<Target = usize>) -> Vec3 {
        self.points[*i.deref()]
    }

    pub fn register_triangles(&mut self, triangles: &[usize]) {
        self.triangles.extend(triangles)
    }

    pub fn mark_debug(&mut self, id: impl Deref<Target = usize>, color: Color) {
        self.debug_points.push((self.point(id), color));
    }

    pub fn add_debug(&mut self, pos: Vec3, color: Color) {
        self.debug_points.push((pos, color));
    }

    pub fn points(&self) -> &[Vec3] {
        &self.points
    }

    pub fn add_contour(&mut self, indices: &[usize]) {
        self.contours.push(indices.to_vec())
    }
}

impl From<StdRng> for GeometryData {
    fn from(rng: StdRng) -> Self {
        Self {
            contours: vec![],
            debug_points: vec![],
            triangles: vec![],
            colors: vec![],
            points: vec![],
            rng,
        }
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
                let (pa, pb, pc) = (self.points[ia], self.points[ib], self.points[ic]);
                let color = Color::srgb(0., 0.4, 0.);
                gizmos.line(pa, pb, color);
                gizmos.line(pb, pc, color);
                gizmos.line(pc, pa, color);
            }
        }

        if debug_flags.contours {
            let mut rng = self.rng.clone();
            for c in self.contours.iter() {
                let t = rng.gen();
                let n = c.len();
                for i in 0..n {
                    let r = i as f32 / n as f32;
                    let color = Color::srgb(t, 0.3 + 0.2 * r, 0.3 + 0.2 * r);
                    let pos1 = self.points[c[i]];
                    let pos2 = self.points[c[(i + 1) % n]];
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
