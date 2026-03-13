use numpy::{IntoPyArray, PyArrayDyn};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rand::SeedableRng;
use std::path::{Path, PathBuf};

use crate::DebugGeometry;

const DEFAULT_RNG_SEED: u64 = 42;

#[pyclass(name = "Seed")]
pub struct PySeed;

#[pymethods]
impl PySeed {
    #[new]
    fn new() -> Self {
        PySeed
    }

    #[pyo3(signature = (*, min_radius=0.05, max_turn_angle=0.9, max_zigzag_angle=0.3, radius_to_branch_ratio=1.8, branch_variance=0.5, main_children_radius_factor=0.75, secondary_children_radius_factor=0.65, radius_variance=0.1, birth_coefficient=0.8, birth_power=1.0, up_attraction_factor=0.15, base_radius=1.0))]
    fn grow_plant(
        &self,
        min_radius: f32,
        max_turn_angle: f32,
        max_zigzag_angle: f32,
        radius_to_branch_ratio: f32,
        branch_variance: f32,
        main_children_radius_factor: f32,
        secondary_children_radius_factor: f32,
        radius_variance: f32,
        birth_coefficient: f32,
        birth_power: f32,
        up_attraction_factor: f32,
        base_radius: f32,
    ) -> PyPlantNode {
        let config = crate::GrowConfig {
            min_radius,
            max_turn_angle,
            max_zigzag_angle,
            radius_to_branch_ratio,
            branch_variance,
            main_children_radius_factor,
            secondary_children_radius_factor,
            radius_variance,
            birth_coefficient,
            birth_power,
            up_attraction_factor,
            base_radius,
        };
        let mut rng = rand::rngs::StdRng::seed_from_u64(DEFAULT_RNG_SEED);
        let node = crate::Seed::grow_plant(&config, &mut rng);
        PyPlantNode(node)
    }
}

#[pyclass(name = "PlantNode")]
pub struct PyPlantNode(crate::PlantNode);

#[pymethods]
impl PyPlantNode {
    fn grow_skeleton(&self) -> PyTreeSkeleton {
        PyTreeSkeleton(self.0.grow_skeleton())
    }
}

#[pyclass(name = "Skeleton")]
pub struct PyTreeSkeleton(crate::TreeSkeleton);

#[pymethods]
impl PyTreeSkeleton {
    #[pyo3(signature = (*, repulsion=0.5, interaction_radius=0.1, n_steps=10, dt=0.01, particles_per_leaf=5))]
    fn grow_strands(
        &self,
        repulsion: f32,
        interaction_radius: f32,
        n_steps: usize,
        dt: f32,
        particles_per_leaf: usize,
    ) -> PyVolumetricTree {
        let config = crate::StrandsConfig {
            repulsion,
            interaction_radius,
            n_steps,
            dt,
            particles_per_leaf,
        };
        let rng = rand::rngs::StdRng::seed_from_u64(DEFAULT_RNG_SEED);
        let mut debug = None;
        let volumetric = self.0.grow_strands_debug(&config, rng, |x| debug = Some(x));
        PyVolumetricTree {
            data: volumetric,
            debug: debug.unwrap(),
        }
    }

    fn debug(&self) -> PyDebugGeometry {
        let skeleton_debug = crate::TreeSkeletonDebugData {
            copy: self.0.clone(),
        };
        PyDebugGeometry(crate::VisualDebug::debug_data(&skeleton_debug))
    }
}

#[pyclass(name = "VolumetricTree")]
pub struct PyVolumetricTree {
    data: crate::VolumetricTree,
    debug: DebugGeometry,
}

#[pymethods]
impl PyVolumetricTree {
    #[pyo3(signature = (*, spacing=0.15, leaf_size=0.5, smoothing_iters=2, leaf_angle=0.4, interior_angle=2.5))]
    fn build_mesh(
        &self,
        spacing: f32,
        leaf_size: f32,
        smoothing_iters: u32,
        leaf_angle: f32,
        interior_angle: f32,
    ) -> PyGeometryData {
        let config = crate::MeshConfig {
            spacing,
            leaf_size,
            leaf_angle,
            interior_angle,
            smoothing_iters,
        };
        let rng = rand::rngs::StdRng::seed_from_u64(DEFAULT_RNG_SEED);
        let mut debug = None;
        let mut geometry = self
            .data
            .build_mesh_debug(&config, rng, |x| debug = Some(x));
        geometry.compute_smooth_normals();
        PyGeometryData {
            data: geometry,
            debug: debug.unwrap(),
        }
    }

    fn debug(&self) -> PyDebugGeometry {
        // FIXME: no clone
        PyDebugGeometry(self.debug.clone())
    }
}

#[pyclass(name = "TreeMesh")]
pub struct PyGeometryData {
    data: crate::GeometryData,
    debug: crate::DebugGeometry,
}

#[pymethods]
impl PyGeometryData {
    fn to_json(&self, include_debug: bool) -> String {
        let mut encoder = crate::TreeEncoder::new();
        let scale = 256i32;

        // Add quantized geometry components with custom prefixes for flexibility
        encoder.add_vec3_components("v", &self.data.points, scale);
        encoder.add_vec3_components("n", &self.data.normals, scale);

        let indices: Vec<i32> = self.data.triangles.iter().map(|&i| i as i32).collect();
        encoder.add_immediate_buffer("indices", &indices);

        // Set required geometry outputs.
        encoder.set_vertices(encoder.make_scaled_vec3("v", 8));
        encoder.set_normals(encoder.make_scaled_vec3("n", 8));

        // Colors are optional in the decoder path. Avoid emitting divp2(vec4(...))
        // when there are no colors, since divp2 on an empty vec4 would fail at runtime.
        if !self.data.colors.is_empty() {
            encoder.add_color_components("c", &self.data.colors, scale);
            encoder.set_colors(encoder.make_scaled_vec4("c", 8));
        }

        encoder.set_triangles(crate::export::Expr::var("indices"));

        if include_debug {
            let debug_geom = crate::VisualDebug::debug_data(&self.data);
            let mut debug_layers = crate::DebugLayers {
                layers: std::collections::HashMap::new(),
            };
            debug_layers.layers.insert(
                "mesh".to_string(),
                encoder.add_debug_layer("mesh", &debug_geom),
            );
            encoder.set_debug(debug_layers);
        }

        encoder.to_json()
    }

    #[getter]
    fn points(&self) -> Vec<[f32; 3]> {
        self.data.points.iter().map(|v| v.to_array()).collect()
    }

    #[getter]
    fn normals(&self) -> Vec<[f32; 3]> {
        self.data.normals.iter().map(|v| v.to_array()).collect()
    }

    #[getter]
    fn colors(&self) -> Vec<[f32; 4]> {
        self.data.colors.iter().copied().collect()
    }

    #[getter]
    fn triangles(&self) -> Vec<u32> {
        self.data.triangles.clone()
    }

    fn debug(&self) -> PyDebugGeometry {
        PyDebugGeometry(self.debug.clone())
    }

    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        let json = self.to_json(false);
        embed_viewer(py, &json)
    }
}

#[pyclass(name = "DebugData")]
pub struct PyDebugGeometry(crate::DebugGeometry);

#[pymethods]
impl PyDebugGeometry {
    fn to_json(&self, layer_name: &str) -> String {
        let mut encoder = crate::TreeEncoder::new();
        let debug_layer = encoder.add_debug_layer(layer_name, &self.0);
        let mut debug_layers = crate::DebugLayers {
            layers: std::collections::HashMap::new(),
        };
        debug_layers
            .layers
            .insert(layer_name.to_string(), debug_layer);
        encoder.set_debug(debug_layers);
        encoder.to_json()
    }

    #[getter]
    fn points<'py>(&self, py: pyo3::Python<'py>) -> pyo3::Bound<'py, PyArrayDyn<f32>> {
        let array = ndarray::Array::from_shape_fn((self.0.points.len(), 3), |(i, j)| match j {
            0 => self.0.points[i].0.x,
            1 => self.0.points[i].0.y,
            _ => self.0.points[i].0.z,
        })
        .into_dyn();

        array.into_pyarray(py)
    }

    #[getter]
    fn lines<'py>(&self, py: pyo3::Python<'py>) -> pyo3::Bound<'py, PyArrayDyn<f32>> {
        let array =
            ndarray::Array::from_shape_fn((self.0.lines.len(), 2, 3), |(i, side, j)| {
                match (side, j) {
                    (0, 0) => self.0.lines[i].0.x,
                    (0, 1) => self.0.lines[i].0.y,
                    (0, _) => self.0.lines[i].0.z,
                    (_, 0) => self.0.lines[i].1.x,
                    (_, 1) => self.0.lines[i].1.y,
                    (_, _) => self.0.lines[i].1.z,
                }
            })
            .into_dyn();

        array.into_pyarray(py)
    }

    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        let json = self.to_json("debug");
        embed_viewer(py, &json)
    }
}

fn bundled_render_js_path(py: Python<'_>) -> PyResult<PathBuf> {
    let tubulin = py.import("tubulin")?;
    let module_file: String = tubulin.getattr("__file__")?.extract()?;
    let module_path = Path::new(&module_file);
    let module_dir = module_path.parent().ok_or_else(|| {
        PyRuntimeError::new_err(format!(
            "Could not resolve package directory from tubulin.__file__={module_file}"
        ))
    })?;
    Ok(module_dir.join("render.js"))
}

fn load_render_js(py: Python<'_>) -> PyResult<String> {
    let path = bundled_render_js_path(py)?;
    std::fs::read_to_string(&path).map_err(|e| {
        PyRuntimeError::new_err(format!(
            "Failed to read bundled viewer JS at {}: {e}",
            path.display()
        ))
    })
}

fn embed_viewer(py: Python<'_>, json_data: &str) -> PyResult<String> {
    let render_js = load_render_js(py)?;
    Ok(format!(
        r#"<div id='main' style='width:100%;height:400px;'></div>
<script type="module">{}</script>
<script type="module">initTreeViewer({});</script>"#,
        render_js, json_data
    ))
}

#[pyclass(name = "DemoMesh")]
pub struct PyDemoMesh {}

#[pymethods]
impl PyDemoMesh {
    #[new]
    pub fn new() -> PyDemoMesh {
        PyDemoMesh {}
    }
    pub fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        let mut encoder = crate::TreeEncoder::new();

        let t_scale = 1024i32;
        let n_splines = 8;
        let n_samples = 8;
        let ncp = 4;

        let cp_z = [-0.5f32, -1.0 / 6.0, 1.0 / 6.0, 0.5];
        let cp_r = [0.5f32, 0.25, 0.25, 0.5];

        let mut spine_t = Vec::with_capacity(n_samples);
        for i in 0..n_samples {
            spine_t.push(
                ((i as f32 / (n_samples - 1) as f32) * (ncp - 1) as f32 * t_scale as f32) as i32,
            );
        }
        encoder.add_immediate_buffer("spine_t", &spine_t);

        for s in 0..n_splines {
            let angle = (s as f32 / n_splines as f32) * std::f32::consts::PI * 2.0;
            let cos_a = angle.cos();
            let sin_a = angle.sin();

            let cpx: Vec<i32> = cp_r.iter().map(|r| (cos_a * r * 256.0) as i32).collect();
            let cpy: Vec<i32> = cp_r.iter().map(|r| (sin_a * r * 256.0) as i32).collect();
            let cpz: Vec<i32> = cp_z.iter().map(|z| (z * 256.0) as i32).collect();

            let cnx: Vec<i32> = cp_r.iter().map(|_| (cos_a * 256.0) as i32).collect();
            let cny: Vec<i32> = cp_r.iter().map(|_| (sin_a * 256.0) as i32).collect();
            let cnz: Vec<i32> = vec![0; ncp];

            let ccr: Vec<i32> = cp_z.iter().map(|z| ((z + 0.5) * 256.0) as i32).collect();
            let ccg: Vec<i32> = cp_z
                .iter()
                .map(|z| ((0.3 + 0.4 * (z + 0.5)) * 256.0) as i32)
                .collect();
            let ccb: Vec<i32> = cp_z
                .iter()
                .map(|z| ((0.5 + 0.5 * (z + 0.5)) * 256.0) as i32)
                .collect();
            let cca: Vec<i32> = vec![256; ncp];

            let delta = |arr: &[i32]| -> Vec<i32> {
                arr.iter()
                    .enumerate()
                    .map(|(i, &v)| if i == 0 { v } else { v - arr[i - 1] })
                    .collect()
            };

            encoder.add_delta_buffer(&format!("sp{}_x", s), &delta(&cpx), 2);
            encoder.add_delta_buffer(&format!("sp{}_y", s), &delta(&cpy), 2);
            encoder.add_delta_buffer(&format!("sp{}_z", s), &delta(&cpz), 2);
            encoder.add_immediate_buffer(&format!("sp{}_nx", s), &cnx);
            encoder.add_immediate_buffer(&format!("sp{}_ny", s), &cny);
            encoder.add_immediate_buffer(&format!("sp{}_nz", s), &cnz);
            encoder.add_immediate_buffer(&format!("sp{}_cr", s), &ccr);
            encoder.add_immediate_buffer(&format!("sp{}_cg", s), &ccg);
            encoder.add_immediate_buffer(&format!("sp{}_cb", s), &ccb);
            encoder.add_immediate_buffer(&format!("sp{}_ca", s), &cca);
        }

        let mut vert_exprs = Vec::new();
        let mut norm_exprs = Vec::new();
        let mut color_exprs = Vec::new();

        let t_expr = crate::export::Expr::divp2(crate::export::Expr::var("spine_t"), 10);

        for s in 0..n_splines {
            let bx = format!("sp{}_x", s);
            let by = format!("sp{}_y", s);
            let bz = format!("sp{}_z", s);
            let bnx = format!("sp{}_nx", s);
            let bny = format!("sp{}_ny", s);
            let bnz = format!("sp{}_nz", s);
            let bcr = format!("sp{}_cr", s);
            let bcg = format!("sp{}_cg", s);
            let bcb = format!("sp{}_cb", s);
            let bca = format!("sp{}_ca", s);

            let cp_vec3 = |xb: &str, yb: &str, zb: &str| -> crate::export::Expr {
                crate::export::Expr::divp2(
                    crate::export::Expr::vec3(
                        crate::export::Expr::cumsum(crate::export::Expr::var(xb)),
                        crate::export::Expr::cumsum(crate::export::Expr::var(yb)),
                        crate::export::Expr::cumsum(crate::export::Expr::var(zb)),
                    ),
                    8,
                )
            };

            let cp_vec4 = |rb: &str, gb: &str, bb: &str, ab: &str| -> crate::export::Expr {
                crate::export::Expr::divp2(
                    crate::export::Expr::vec4(
                        crate::export::Expr::cumsum(crate::export::Expr::var(rb)),
                        crate::export::Expr::cumsum(crate::export::Expr::var(gb)),
                        crate::export::Expr::cumsum(crate::export::Expr::var(bb)),
                        crate::export::Expr::cumsum(crate::export::Expr::var(ab)),
                    ),
                    8,
                )
            };

            let spline_expr = crate::export::Expr::spline(cp_vec3(&bx, &by, &bz), t_expr.clone());
            let nspline_expr =
                crate::export::Expr::spline(cp_vec3(&bnx, &bny, &bnz), t_expr.clone());
            let cspline_expr =
                crate::export::Expr::spline(cp_vec4(&bcr, &bcg, &bcb, &bca), t_expr.clone());

            vert_exprs.push(spline_expr);
            norm_exprs.push(nspline_expr);
            color_exprs.push(cspline_expr);
        }

        let vertices = crate::export::Expr::interleave(vert_exprs);
        let normals = crate::export::Expr::interleave(norm_exprs);
        let colors = crate::export::Expr::interleave(color_exprs);

        encoder.set_vertices(vertices);
        encoder.set_normals(normals);
        encoder.set_colors(colors);

        let _n_verts = n_samples * n_splines;
        let mut indices_i = Vec::new();
        let mut indices_j = Vec::new();
        let mut indices_k = Vec::new();

        for i in 0..n_samples - 1 {
            for s in 0..n_splines {
                let a = i * n_splines + s;
                let b = i * n_splines + (s + 1) % n_splines;
                let c = (i + 1) * n_splines + s;
                let d = (i + 1) * n_splines + (s + 1) % n_splines;

                indices_i.push(a as i32);
                indices_j.push(c as i32);
                indices_k.push(b as i32);

                indices_i.push(b as i32);
                indices_j.push(c as i32);
                indices_k.push(d as i32);
            }
        }

        let mut delta_i = Vec::new();
        let mut delta_j = Vec::new();
        let mut delta_k = Vec::new();
        let mut prev_i = 0i32;
        let mut prev_j = 0i32;
        let mut prev_k = 0i32;

        for t in 0..indices_i.len() {
            delta_i.push(indices_i[t] - prev_i);
            delta_j.push(indices_j[t] - prev_j);
            delta_k.push(indices_k[t] - prev_k);
            prev_i = indices_i[t];
            prev_j = indices_j[t];
            prev_k = indices_k[t];
        }

        encoder.add_delta_buffer("indices_i", &delta_i, 2);
        encoder.add_delta_buffer("indices_j", &delta_j, 2);
        encoder.add_delta_buffer("indices_k", &delta_k, 2);

        let triangles = crate::export::Expr::triangle(
            crate::export::Expr::cumsum(crate::export::Expr::var("indices_i")),
            crate::export::Expr::cumsum(crate::export::Expr::var("indices_j")),
            crate::export::Expr::cumsum(crate::export::Expr::var("indices_k")),
        );
        encoder.set_triangles(triangles);

        embed_viewer(py, &encoder.to_json())
    }
}

#[pymodule]
fn _tubulin(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySeed>()?;
    m.add_class::<PyPlantNode>()?;
    m.add_class::<PyTreeSkeleton>()?;
    m.add_class::<PyVolumetricTree>()?;
    m.add_class::<PyGeometryData>()?;
    m.add_class::<PyDebugGeometry>()?;
    m.add_class::<PyDemoMesh>()?;
    Ok(())
}
