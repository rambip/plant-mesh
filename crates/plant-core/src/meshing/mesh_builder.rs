#[cfg(feature = "bevy")]
use bevy::prelude::Component;
#[cfg(feature = "bevy")]
use bevy_color::Color;
#[cfg(feature = "bevy")]
use bevy_gizmos::prelude::Gizmos;
use glam::Vec3;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
#[cfg_attr(feature = "bevy", derive(Component))]
pub struct GeometryData {
    pub points: Vec<Vec3>,
    pub normals: Vec<Vec3>,
    pub colors: Vec<[f32; 4]>,
    pub triangles: Vec<u32>,
    pub contours: Vec<Vec<usize>>,
    pub debug_points: Vec<(Vec3, [f32; 4])>,
    #[serde(skip, default = "default_rng")]
    pub rng: rand::rngs::StdRng,
}

fn default_rng() -> rand::rngs::StdRng {
    rand::SeedableRng::seed_from_u64(0)
}

impl GeometryData {
    pub fn new(rng: rand::rngs::StdRng) -> Self {
        Self {
            points: Vec::new(),
            normals: Vec::new(),
            colors: Vec::new(),
            triangles: Vec::new(),
            contours: Vec::new(),
            debug_points: Vec::new(),
            rng,
        }
    }

    /// Computes smooth normals for the mesh using angle-weighted face normals.
    ///
    /// This method iterates over all triangles, calculates the face normal,
    /// and accumulates it into the vertex normals weighted by the angle at each vertex.
    /// Finally, it normalizes all accumulated vectors.
    pub fn compute_smooth_normals(&mut self) {
        let mut normals = vec![Vec3::ZERO; self.points.len()];

        for chunk in self.triangles.chunks_exact(3) {
            let (ia, ib, ic) = (chunk[0] as usize, chunk[1] as usize, chunk[2] as usize);
            let (pa, pb, pc) = (self.points[ia], self.points[ib], self.points[ic]);

            let dir_ab = (pb - pa).normalize_or_zero();
            let dir_ac = (pc - pa).normalize_or_zero();
            let dir_ba = (pa - pb).normalize_or_zero();
            let dir_bc = (pc - pb).normalize_or_zero();
            let dir_ca = (pa - pc).normalize_or_zero();
            let dir_cb = (pb - pc).normalize_or_zero();

            let weight_a = dir_ab.dot(dir_ac).clamp(-1.0, 1.0).acos();
            let weight_b = dir_ba.dot(dir_bc).clamp(-1.0, 1.0).acos();
            let weight_c = dir_ca.dot(dir_cb).clamp(-1.0, 1.0).acos();

            let face_normal = (pb - pa).cross(pc - pa).normalize_or_zero();

            normals[ia] += face_normal * weight_a;
            normals[ib] += face_normal * weight_b;
            normals[ic] += face_normal * weight_c;
        }

        self.normals = normals.into_iter().map(|n| n.normalize_or_zero()).collect();
    }

    pub fn register_points_trunk(
        &mut self,
        points: impl IntoIterator<Item = Vec3>,
        rng: &mut impl rand::Rng,
    ) -> Vec<usize> {
        let mut result = Vec::new();
        for p in points {
            let blue: f32 = rng.gen_range(0.1f32..0.2);
            self.colors.push([0.4, 0.3, blue, 1.0]);
            result.push(self.points.len());
            self.points.push(p);
        }
        result
    }

    pub fn register_points_leaf(&mut self, points: impl IntoIterator<Item = Vec3>) -> Vec<usize> {
        let mut result = Vec::new();
        for p in points {
            self.colors.push([0.1, 0.8, 0.3, 1.0]);
            result.push(self.points.len());
            self.points.push(p);
        }
        result
    }

    pub fn register_triangles(&mut self, triangles: &[usize]) {
        self.triangles.extend(triangles.iter().map(|&i| i as u32));
    }

    pub fn add_contour(&mut self, indices: &[usize]) {
        self.contours.push(indices.to_vec());
    }

    pub fn point(&self, i: usize) -> Vec3 {
        self.points[i]
    }

    pub fn mark_debug(&mut self, id: usize, color: [f32; 4]) {
        self.debug_points.push((self.point(id), color));
    }

    pub fn add_debug(&mut self, pos: Vec3, color: [f32; 4]) {
        self.debug_points.push((pos, color));
    }
}

impl Default for GeometryData {
    fn default() -> Self {
        Self::new(default_rng())
    }
}

#[derive(Copy, Clone, Debug, Default, Serialize, Deserialize)]
pub struct MeshDebugFlags {
    pub triangles: bool,
    pub contours: bool,
    pub other: bool,
}

impl crate::VisualDebug for GeometryData {
    type Flags = MeshDebugFlags;
    #[cfg(feature = "bevy")]
    fn debug(&self, gizmos: &mut Gizmos, debug_flags: Self::Flags) {
        if debug_flags.triangles {
            for i in 0..self.triangles.len() / 3 {
                let (ia, ib, ic) = (
                    self.triangles[3 * i] as usize,
                    self.triangles[3 * i + 1] as usize,
                    self.triangles[3 * i + 2] as usize,
                );
                let (pa, pb, pc) = (self.points[ia], self.points[ib], self.points[ic]);
                let color = Color::srgb(0., 0.4, 0.);
                gizmos.line(pa, pb, color);
                gizmos.line(pb, pc, color);
                gizmos.line(pc, pa, color);
            }
        }

        if debug_flags.contours {
            for c in self.contours.iter() {
                let n = c.len();
                for i in 0..n {
                    let r = i as f32 / n as f32;
                    let color = Color::srgb(0.5, 0.3 + 0.2 * r, 0.3 + 0.2 * r);
                    let pos1 = self.points[c[i]];
                    let pos2 = self.points[c[(i + 1) % n]];
                    gizmos.line(pos1, pos2, color);
                }
            }
        }

        if debug_flags.other {
            for (p, c) in self.debug_points.iter() {
                let color = Color::linear_rgba(c[0], c[1], c[2], c[3]);
                gizmos.cross(*p, 0.2, color);
            }
        }
    }
}

#[derive(Copy, Clone, Serialize, Deserialize)]
pub struct MeshConfig {
    pub leaf_size: f32,
    pub leaf_angle: f32,
    pub interior_angle: f32,
    pub spacing: f32,
}
