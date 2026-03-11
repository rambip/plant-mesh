#[cfg(feature = "bevy")]
use bevy::prelude::Component;

#[cfg(feature = "python")]
use crate::export::Expr;
#[cfg(feature = "python")]
use crate::TreeEncoder;
use crate::VisualDebug;
use glam::Vec3;
#[cfg(feature = "python")]
use numpy::IntoPyArray;
#[cfg(feature = "python")]
use numpy::PyArrayDyn;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
#[cfg_attr(feature = "bevy", derive(Component))]
#[cfg_attr(feature = "python", pyo3::pyclass)]
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

#[cfg(feature = "python")]
#[pyo3::pymethods]
impl GeometryData {
    #[getter]
    fn points<'py>(&self, py: pyo3::Python<'py>) -> pyo3::Bound<'py, PyArrayDyn<f32>> {
        let array = ndarray::Array::from_shape_fn((self.points.len(), 3), |(i, j)| match j {
            0 => self.points[i].x,
            1 => self.points[i].y,
            _ => self.points[i].z,
        })
        .into_dyn();

        array.into_pyarray(py)
    }
    #[getter]
    fn normals<'py>(&self, py: pyo3::Python<'py>) -> pyo3::Bound<'py, PyArrayDyn<f32>> {
        let array = ndarray::Array::from_shape_fn((self.normals.len(), 3), |(i, j)| match j {
            0 => self.normals[i].x,
            1 => self.normals[i].y,
            _ => self.normals[i].z,
        })
        .into_dyn();

        array.into_pyarray(py)
    }
    #[getter]
    fn colors<'py>(&self, py: pyo3::Python<'py>) -> pyo3::Bound<'py, PyArrayDyn<f32>> {
        let array =
            ndarray::Array::from_shape_fn((self.colors.len(), 4), |(i, j)| self.colors[i][j])
                .into_dyn();

        array.into_pyarray(py)
    }
    #[getter]
    fn triangles<'py>(&self, py: pyo3::Python<'py>) -> pyo3::Bound<'py, PyArrayDyn<u32>> {
        let array =
            ndarray::Array::from_shape_fn((self.triangles.len() / 3, 3), |(i, j)| match j {
                0 => self.triangles[3 * i],
                1 => self.triangles[3 * i + 1],
                _ => self.triangles[3 * i + 2],
            })
            .into_dyn();

        array.into_pyarray(py)
    }

    fn to_json(&self, include_debug: bool) -> String {
        use crate::export::Expr;
        use std::collections::HashMap;

        let mut encoder = TreeEncoder::new();

        let quant_scale = 256i32;

        let vx: Vec<i32> = self
            .points
            .iter()
            .map(|p| (p.x * quant_scale as f32) as i32)
            .collect();
        let vy: Vec<i32> = self
            .points
            .iter()
            .map(|p| (p.y * quant_scale as f32) as i32)
            .collect();
        let vz: Vec<i32> = self
            .points
            .iter()
            .map(|p| (p.z * quant_scale as f32) as i32)
            .collect();

        let nx: Vec<i32> = self.normals.iter().map(|n| (n.x * 127.0) as i32).collect();
        let ny: Vec<i32> = self.normals.iter().map(|n| (n.y * 127.0) as i32).collect();
        let nz: Vec<i32> = self.normals.iter().map(|n| (n.z * 127.0) as i32).collect();

        let cr: Vec<i32> = self.colors.iter().map(|c| (c[0] * 255.0) as i32).collect();
        let cg: Vec<i32> = self.colors.iter().map(|c| (c[1] * 255.0) as i32).collect();
        let cb: Vec<i32> = self.colors.iter().map(|c| (c[2] * 255.0) as i32).collect();
        let ca: Vec<i32> = self.colors.iter().map(|c| (c[3] * 255.0) as i32).collect();

        let n_triangles = self.triangles.len() / 3;

        let mut ix = Vec::with_capacity(n_triangles);
        let mut iy = Vec::with_capacity(n_triangles);
        let mut iz = Vec::with_capacity(n_triangles);

        for i in 0..n_triangles {
            ix.push(self.triangles[3 * i] as i32);
            iy.push(self.triangles[3 * i + 1] as i32);
            iz.push(self.triangles[3 * i + 2] as i32);
        }

        // Delta encode indices
        let mut delta_ix = Vec::with_capacity(n_triangles);
        let mut delta_iy = Vec::with_capacity(n_triangles);
        let mut delta_iz = Vec::with_capacity(n_triangles);
        let mut prev_ix = 0i32;
        let mut prev_iy = 0i32;
        let mut prev_iz = 0i32;

        for i in 0..n_triangles {
            delta_ix.push(ix[i] - prev_ix);
            delta_iy.push(iy[i] - prev_iy);
            delta_iz.push(iz[i] - prev_iz);
            prev_ix = ix[i];
            prev_iy = iy[i];
            prev_iz = iz[i];
        }

        encoder.add_delta_buffer("ix", &delta_ix, 2);
        encoder.add_delta_buffer("iy", &delta_iy, 2);
        encoder.add_delta_buffer("iz", &delta_iz, 2);

        let vertices = Expr::divp2(
            Expr::vec3(Expr::var("vx"), Expr::var("vy"), Expr::var("vz")),
            8,
        );
        let normals = Expr::divp2(
            Expr::vec3(Expr::var("nx"), Expr::var("ny"), Expr::var("nz")),
            7,
        );
        let triangles = Expr::triangle(
            Expr::cumsum(Expr::var("ix")),
            Expr::cumsum(Expr::var("iy")),
            Expr::cumsum(Expr::var("iz")),
        );

        encoder.set_vertices(vertices);
        encoder.set_normals(normals);
        encoder.set_triangles(triangles);

        if !self.colors.is_empty() {
            let colors = Expr::divp2(
                Expr::vec4(
                    Expr::var("cr"),
                    Expr::var("cg"),
                    Expr::var("cb"),
                    Expr::var("ca"),
                ),
                8,
            );
            encoder.set_colors(colors);
        }

        if include_debug {
            let debug_data = self.debug_data();
            let mut debug_layers = crate::DebugLayers {
                layers: HashMap::new(),
            };
            debug_layers.layers.insert(
                "mesh".to_string(),
                encoder.add_debug_layer("mesh", &debug_data),
            );
            encoder.set_debug(debug_layers);
        }

        encoder.to_json()
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
    fn debug_data(&self) -> crate::DebugGeometry {
        let mut out = crate::DebugGeometry::new();
        for i in 0..self.triangles.len() / 3 {
            let (ia, ib, ic) = (
                self.triangles[3 * i] as usize,
                self.triangles[3 * i + 1] as usize,
                self.triangles[3 * i + 2] as usize,
            );
            let (pa, pb, pc) = (self.points[ia], self.points[ib], self.points[ic]);
            let color = crate::DebugColor::rgb(0., 0.4, 0.);
            out.lines.push((pa, pb, color));
            out.lines.push((pb, pc, color));
            out.lines.push((pc, pa, color));
        }

        for c in self.contours.iter() {
            let n = c.len();
            for i in 0..n {
                let r = i as f32 / n as f32;
                let color = crate::DebugColor::rgb(0.5, 0.3 + 0.2 * r, 0.3 + 0.2 * r);
                let pos1 = self.points[c[i]];
                let pos2 = self.points[c[(i + 1) % n]];
                out.lines.push((pos1, pos2, color));
            }
        }

        for (p, c) in self.debug_points.iter() {
            let color = crate::DebugColor(*c);
            out.points.push((*p, color));
        }

        out
    }
}

#[derive(Copy, Clone, Serialize, Deserialize)]
pub struct MeshConfig {
    pub leaf_size: f32,
    pub leaf_angle: f32,
    pub interior_angle: f32,
    pub spacing: f32,
    pub smoothing_iters: u32,
}

impl Default for MeshConfig {
    fn default() -> Self {
        Self {
            leaf_size: 0.5,
            leaf_angle: 0.5,
            interior_angle: 0.3,
            spacing: 0.5,
            smoothing_iters: 2,
        }
    }
}
