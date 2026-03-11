use pyo3::prelude::*;
use rand::SeedableRng;

#[pyclass(name = "Seed")]
pub struct PySeed;

#[pymethods]
impl PySeed {
    #[new]
    fn new() -> Self {
        PySeed
    }

    fn grow_plant(&self) -> PyPlantNode {
        let config = crate::GrowConfig::default();
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
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
    #[pyo3(signature = (*, repulsion=0.5, interaction_radius=0.1, n_steps=10, dt=0.1, particles_per_leaf=20))]
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
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let volumetric = self.0.grow_strands(&config, rng);
        PyVolumetricTree(volumetric)
    }

    fn debug_data(&self) -> PyDebugData {
        let debug = crate::TreeSkeletonDebugData::default();
        let debug_geom = crate::VisualDebug::debug_data(&debug);
        PyDebugData(debug_geom)
    }

    fn to_json(&self) -> String {
        let grow_config = crate::GrowConfig::default();
        let strands_config = crate::StrandsConfig::default();
        let mesh_config = crate::MeshConfig::default();

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let node = crate::Seed::grow_plant(&grow_config, &mut rng);
        let skeleton = node.grow_skeleton();
        let volumetric = skeleton.grow_strands(&strands_config, rng.clone());
        let geometry = volumetric.build_mesh(&mesh_config, rng);

        let mut encoder = crate::TreeEncoder::new();

        let scale = 256i32;
        let vx: Vec<i32> = geometry
            .points
            .iter()
            .map(|p| (p.x * scale as f32) as i32)
            .collect();
        let vy: Vec<i32> = geometry
            .points
            .iter()
            .map(|p| (p.y * scale as f32) as i32)
            .collect();
        let vz: Vec<i32> = geometry
            .points
            .iter()
            .map(|p| (p.z * scale as f32) as i32)
            .collect();
        let nx: Vec<i32> = geometry
            .normals
            .iter()
            .map(|n| (n.x * scale as f32) as i32)
            .collect();
        let ny: Vec<i32> = geometry
            .normals
            .iter()
            .map(|n| (n.y * scale as f32) as i32)
            .collect();
        let nz: Vec<i32> = geometry
            .normals
            .iter()
            .map(|n| (n.z * scale as f32) as i32)
            .collect();
        let cr: Vec<i32> = geometry
            .colors
            .iter()
            .map(|c| (c[0] * scale as f32) as i32)
            .collect();
        let cg: Vec<i32> = geometry
            .colors
            .iter()
            .map(|c| (c[1] * scale as f32) as i32)
            .collect();
        let cb: Vec<i32> = geometry
            .colors
            .iter()
            .map(|c| (c[2] * scale as f32) as i32)
            .collect();
        let ca: Vec<i32> = geometry
            .colors
            .iter()
            .map(|c| (c[3] * scale as f32) as i32)
            .collect();
        let indices: Vec<i32> = geometry.triangles.iter().map(|&i| i as i32).collect();

        encoder.add_immediate_buffer("vx", &vx);
        encoder.add_immediate_buffer("vy", &vy);
        encoder.add_immediate_buffer("vz", &vz);
        encoder.add_immediate_buffer("nx", &nx);
        encoder.add_immediate_buffer("ny", &ny);
        encoder.add_immediate_buffer("nz", &nz);
        encoder.add_immediate_buffer("cr", &cr);
        encoder.add_immediate_buffer("cg", &cg);
        encoder.add_immediate_buffer("cb", &cb);
        encoder.add_immediate_buffer("ca", &ca);
        encoder.add_immediate_buffer("indices", &indices);

        let vertices = crate::export::Expr::divp2(
            crate::export::Expr::vec3(
                crate::export::Expr::var("vx"),
                crate::export::Expr::var("vy"),
                crate::export::Expr::var("vz"),
            ),
            8,
        );
        let normals = crate::export::Expr::divp2(
            crate::export::Expr::vec3(
                crate::export::Expr::var("nx"),
                crate::export::Expr::var("ny"),
                crate::export::Expr::var("nz"),
            ),
            8,
        );
        let colors = crate::export::Expr::divp2(
            crate::export::Expr::vec4(
                crate::export::Expr::var("cr"),
                crate::export::Expr::var("cg"),
                crate::export::Expr::var("cb"),
                crate::export::Expr::var("ca"),
            ),
            8,
        );
        let triangles = crate::export::Expr::var("indices");

        encoder.set_vertices(vertices);
        encoder.set_normals(normals);
        encoder.set_colors(colors);
        encoder.set_triangles(triangles);

        let debug = crate::TreeSkeletonDebugData::default();
        let debug_geom = crate::VisualDebug::debug_data(&debug);
        let mut debug_layers = crate::DebugLayers {
            layers: std::collections::HashMap::new(),
        };
        debug_layers.layers.insert(
            "skeleton".to_string(),
            encoder.add_debug_layer("skeleton", &debug_geom),
        );
        encoder.set_debug(debug_layers);
        encoder.to_json()
    }

    fn _repr_html_(&self) -> String {
        "TODO: Skeleton viewer".to_string()
    }
}

#[pyclass(name = "VolumetricTree")]
pub struct PyVolumetricTree(crate::VolumetricTree);

#[pymethods]
impl PyVolumetricTree {
    #[pyo3(signature = (*, spacing=0.5, leaf_size=0.5, smoothing_iters=2))]
    fn build_mesh(&self, spacing: f32, leaf_size: f32, smoothing_iters: u32) -> PyGeometryData {
        let config = crate::MeshConfig {
            spacing,
            leaf_size,
            leaf_angle: 0.5,
            interior_angle: 0.3,
            smoothing_iters,
        };
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let geometry = self.0.build_mesh(&config, rng);
        PyGeometryData(geometry)
    }
}

#[pyclass(name = "TreeMesh")]
pub struct PyGeometryData(crate::GeometryData);

#[pymethods]
impl PyGeometryData {
    fn to_json(&self, include_debug: bool) -> String {
        let mut encoder = crate::TreeEncoder::new();

        let scale = 256i32;

        let vx: Vec<i32> = self
            .0
            .points
            .iter()
            .map(|p| (p.x * scale as f32) as i32)
            .collect();
        let vy: Vec<i32> = self
            .0
            .points
            .iter()
            .map(|p| (p.y * scale as f32) as i32)
            .collect();
        let vz: Vec<i32> = self
            .0
            .points
            .iter()
            .map(|p| (p.z * scale as f32) as i32)
            .collect();

        let nx: Vec<i32> = self
            .0
            .normals
            .iter()
            .map(|n| (n.x * scale as f32) as i32)
            .collect();
        let ny: Vec<i32> = self
            .0
            .normals
            .iter()
            .map(|n| (n.y * scale as f32) as i32)
            .collect();
        let nz: Vec<i32> = self
            .0
            .normals
            .iter()
            .map(|n| (n.z * scale as f32) as i32)
            .collect();

        let cr: Vec<i32> = self
            .0
            .colors
            .iter()
            .map(|c| (c[0] * scale as f32) as i32)
            .collect();
        let cg: Vec<i32> = self
            .0
            .colors
            .iter()
            .map(|c| (c[1] * scale as f32) as i32)
            .collect();
        let cb: Vec<i32> = self
            .0
            .colors
            .iter()
            .map(|c| (c[2] * scale as f32) as i32)
            .collect();
        let ca: Vec<i32> = self
            .0
            .colors
            .iter()
            .map(|c| (c[3] * scale as f32) as i32)
            .collect();

        let indices: Vec<i32> = self.0.triangles.iter().map(|&i| i as i32).collect();

        encoder.add_immediate_buffer("vx", &vx);
        encoder.add_immediate_buffer("vy", &vy);
        encoder.add_immediate_buffer("vz", &vz);
        encoder.add_immediate_buffer("nx", &nx);
        encoder.add_immediate_buffer("ny", &ny);
        encoder.add_immediate_buffer("nz", &nz);
        encoder.add_immediate_buffer("cr", &cr);
        encoder.add_immediate_buffer("cg", &cg);
        encoder.add_immediate_buffer("cb", &cb);
        encoder.add_immediate_buffer("ca", &ca);
        encoder.add_immediate_buffer("indices", &indices);

        let vertices = crate::export::Expr::divp2(
            crate::export::Expr::vec3(
                crate::export::Expr::var("vx"),
                crate::export::Expr::var("vy"),
                crate::export::Expr::var("vz"),
            ),
            8,
        );

        let normals = crate::export::Expr::divp2(
            crate::export::Expr::vec3(
                crate::export::Expr::var("nx"),
                crate::export::Expr::var("ny"),
                crate::export::Expr::var("nz"),
            ),
            8,
        );

        let colors = crate::export::Expr::divp2(
            crate::export::Expr::vec4(
                crate::export::Expr::var("cr"),
                crate::export::Expr::var("cg"),
                crate::export::Expr::var("cb"),
                crate::export::Expr::var("ca"),
            ),
            8,
        );

        let triangles = crate::export::Expr::var("indices");

        encoder.set_vertices(vertices);
        encoder.set_normals(normals);
        encoder.set_colors(colors);
        encoder.set_triangles(triangles);

        if include_debug {
            let debug_geom = crate::VisualDebug::debug_data(&self.0);
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
    fn points(&self) -> Vec<Vec<f32>> {
        self.0.points.iter().map(|v| vec![v.x, v.y, v.z]).collect()
    }

    #[getter]
    fn normals(&self) -> Vec<Vec<f32>> {
        self.0.normals.iter().map(|v| vec![v.x, v.y, v.z]).collect()
    }

    #[getter]
    fn colors(&self) -> Vec<Vec<f32>> {
        self.0
            .colors
            .iter()
            .map(|c| vec![c[0], c[1], c[2], c[3]])
            .collect()
    }

    #[getter]
    fn triangles(&self) -> Vec<u32> {
        self.0.triangles.clone()
    }

    fn _repr_html_(&self) -> String {
        let json = self.to_json(false);
        format!(
            r#"<iframe srcdoc="{}" width="100%" height="500" frameborder="0"></iframe>"#,
            html_escape(&json)
        )
    }
}

#[pyclass(name = "DebugData")]
pub struct PyDebugData(crate::DebugGeometry);

#[pymethods]
impl PyDebugData {
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
}

fn html_escape(s: &str) -> String {
    s.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace("\"", "&quot;")
        .replace("'", "&#39;")
}

#[pyfunction]
pub fn build_demo_tree() -> PyGeometryData {
    let geometry = crate::GeometryData::build_from_config(
        &crate::TreeConfig {
            grow: crate::GrowConfig::default(),
            strands: crate::StrandsConfig::default(),
            mesh: crate::MeshConfig::default(),
        },
        42,
    );
    PyGeometryData(geometry)
}

#[pyfunction]
pub fn demo_mesh() -> String {
    let mut encoder = crate::TreeEncoder::new();

    let t_scale = 1024i32;
    let n_splines = 8;
    let n_samples = 8;
    let ncp = 4;

    let cp_z = [-0.5f32, -1.0 / 6.0, 1.0 / 6.0, 0.5];
    let cp_r = [0.5f32, 0.25, 0.25, 0.5];

    let mut spine_t = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        spine_t
            .push(((i as f32 / (n_samples - 1) as f32) * (ncp - 1) as f32 * t_scale as f32) as i32);
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

        let ccr: Vec<i32> = cp_z.iter().map(|z| (((z + 0.5) * 256.0) as i32)).collect();
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
        let nspline_expr = crate::export::Expr::spline(cp_vec3(&bnx, &bny, &bnz), t_expr.clone());
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

    let n_verts = n_samples * n_splines;
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

    let debug_geom = crate::DebugGeometry::new();
    let mut debug_layers = crate::DebugLayers {
        layers: std::collections::HashMap::new(),
    };
    debug_layers.layers.insert(
        "skeleton".to_string(),
        encoder.add_debug_layer("skeleton", &debug_geom),
    );
    encoder.set_debug(debug_layers);

    encoder.to_json()
}

#[pyfunction]
pub fn debug_to_json(debug: &PyDebugData, layer_name: &str) -> String {
    debug.to_json(layer_name)
}

#[pymodule]
fn _tubulin(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySeed>()?;
    m.add_class::<PyPlantNode>()?;
    m.add_class::<PyTreeSkeleton>()?;
    m.add_class::<PyVolumetricTree>()?;
    m.add_class::<PyGeometryData>()?;
    m.add_class::<PyDebugData>()?;
    m.add_function(wrap_pyfunction!(build_demo_tree, m)?)?;
    m.add_function(wrap_pyfunction!(demo_mesh, m)?)?;
    m.add_function(wrap_pyfunction!(debug_to_json, m)?)?;
    Ok(())
}
