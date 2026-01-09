use glam::Vec3;
use serde::{Serialize, Deserialize};
#[cfg(feature = "bevy")]
use bevy_gizmos::prelude::Gizmos;
#[cfg(feature = "bevy")]
use bevy_color::Color;

#[derive(Serialize, Deserialize)]
pub struct GeometryData {
    pub points: Vec<Vec3>,
    pub colors: Vec<[f32; 4]>,
    pub triangles: Vec<u32>,
    pub contours: Vec<Vec<usize>>,
    pub debug_points: Vec<(Vec3, [f32; 4])>,
}

impl GeometryData {
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            colors: Vec::new(),
            triangles: Vec::new(),
            contours: Vec::new(),
            debug_points: Vec::new(),
        }
    }

    pub fn to_mesh_tools(&self) -> mesh_tools::models::Mesh {
        mesh_tools::models::Mesh::default()
    }

    pub fn register_points_trunk(&mut self, points: impl IntoIterator<Item = Vec3>, rng: &mut impl rand::Rng) -> Vec<usize> {
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
    pub interaction_radius: f32,
    pub repulsion: f32,
    pub n_steps: usize,
    pub dt: f32,
    pub particles_per_leaf: usize,
}

impl From<super::particles::StrandsConfig> for MeshConfig {
    fn from(s: super::particles::StrandsConfig) -> Self {
        Self {
            interaction_radius: s.interaction_radius,
            repulsion: s.repulsion,
            n_steps: s.n_steps,
            dt: s.dt,
            particles_per_leaf: s.particles_per_leaf,
        }
    }
}

impl From<MeshConfig> for super::particles::StrandsConfig {
    fn from(m: MeshConfig) -> Self {
        Self {
            interaction_radius: m.interaction_radius,
            repulsion: m.repulsion,
            n_steps: m.n_steps,
            dt: m.dt,
            particles_per_leaf: m.particles_per_leaf,
        }
    }
}
