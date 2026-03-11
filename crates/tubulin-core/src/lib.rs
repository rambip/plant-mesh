pub mod export;
pub mod growing;
pub mod meshing;
pub mod utils;

pub use growing::generation::GrowConfig;
pub use growing::NodeInfo;
pub use growing::PlantNode;
pub use growing::PlantNodeProps;
pub use growing::TreeSkeleton;
pub use growing::TreeSkeletonDebugData;
pub use meshing::mesh_builder::MeshConfig;
pub use meshing::mesh_builder::MeshDebugFlags;
pub use meshing::particles::TrajectoryBuilder;
pub use meshing::GeometryData;
pub use meshing::SplineIndex;
pub use meshing::StrandsConfig;

use glam::Quat;
use glam::Vec3;
pub use meshing::VolumetricTree;
use rand::SeedableRng;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Default)]
pub struct DebugColor(pub [f32; 4]);

impl DebugColor {
    pub fn rgb(r: f32, g: f32, b: f32) -> Self {
        Self([r, g, b, 1.0])
    }
    pub fn rgba(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self([r, g, b, a])
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Circle {
    pub position: Vec3,
    pub orientation: Quat,
    pub radius: f32,
}

#[derive(Debug, Default)]
pub struct DebugGeometry {
    pub lines: Vec<(Vec3, Vec3, DebugColor)>,
    pub circles: Vec<(Circle, DebugColor)>,
    pub points: Vec<(Vec3, DebugColor)>,
}

impl DebugGeometry {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn layer_name(&self) -> &'static str {
        "unknown"
    }
}

pub use export::{DebugLayer, DebugLayers, DebugLines, DebugPoints, TreeEncoder};

pub trait VisualDebug {
    fn debug_data(&self) -> DebugGeometry;
}

#[derive(Serialize, Deserialize, Clone)]
#[cfg_attr(
    feature = "bevy",
    derive(bevy::prelude::TypePath, bevy::prelude::Asset)
)]
pub struct TreeConfig {
    pub grow: GrowConfig,
    pub strands: StrandsConfig,
    pub mesh: MeshConfig,
}

impl GeometryData {
    pub fn build_from_config(config: &TreeConfig, seed: u64) -> Self {
        let rng = rand::rngs::StdRng::seed_from_u64(seed);

        let mut plant_builder = rng.clone();
        let rng_for_strands = rng.clone();
        let rng_for_mesh = rng.clone();

        let mut geometry = Seed
            .grow::<PlantNode>(&config.grow, &mut plant_builder)
            .grow_skeleton()
            .grow_strands(&config.strands, rng_for_strands)
            .build_mesh(&config.mesh, rng_for_mesh);

        geometry.compute_smooth_normals();
        geometry
    }
}

pub trait Grow {
    fn grow<Next>(self, config: &Next::Config, builder: &mut Next::Builder) -> Next
    where
        Next: TreePipelinePhase<Previous = Self>;
}

impl<T> Grow for T {
    fn grow<Next>(
        self,
        config: &<Next as TreePipelinePhase>::Config,
        builder: &mut <Next as TreePipelinePhase>::Builder,
    ) -> Next
    where
        Next: TreePipelinePhase<Previous = T>,
    {
        Next::generate_from(self, config, builder)
    }
}

pub trait TreePipelinePhase {
    type Previous;
    type Config;
    type Builder;
    fn generate_from(
        prev: Self::Previous,
        config: &Self::Config,
        builder: &mut Self::Builder,
    ) -> Self;
}

pub struct Seed;

impl Seed {
    pub fn grow_plant(config: &GrowConfig, rng: &mut rand::rngs::StdRng) -> PlantNode {
        let root = PlantNodeProps {
            radius: config.base_radius,
            orientation: glam::Quat::IDENTITY,
            position: glam::Vec3::ZERO,
        };
        growing::generation::grow_tree_basic(config, rng, root, 0)
    }
}

impl TreePipelinePhase for PlantNode {
    type Previous = Seed;
    type Config = GrowConfig;
    type Builder = rand::rngs::StdRng;
    fn generate_from(_: Self::Previous, config: &Self::Config, rng: &mut Self::Builder) -> Self {
        Seed::grow_plant(config, rng)
    }
}

impl TreePipelinePhase for TreeSkeleton {
    type Previous = PlantNode;
    type Config = ();
    type Builder = TreeSkeletonDebugData;
    fn generate_from(prev: Self::Previous, _: &Self::Config, cache: &mut Self::Builder) -> Self {
        let mut node_props = Vec::new();
        prev.register_node_properties(&mut node_props);
        let mut node_info = Vec::new();
        prev.register_node_info(&mut node_info, 0);
        cache.copy = TreeSkeleton {
            node_info,
            node_props,
        };
        cache.copy.clone()
    }
}

impl TreePipelinePhase for TrajectoryBuilder {
    type Previous = TreeSkeleton;
    type Config = StrandsConfig;
    type Builder = rand::rngs::StdRng;
    fn generate_from(prev: Self::Previous, config: &Self::Config, rng: &mut Self::Builder) -> Self {
        let mut builder = TrajectoryBuilder::new(rng.clone());
        builder.clear_for_tree(&prev);
        builder.compute_trajectories(&prev, prev.root(), config);
        builder
    }
}

#[cfg(feature = "python")]
#[pyo3::pymodule]
mod _tubulin {
    use crate::{GeometryData, TreeConfig, TreeEncoder};
    #[pyo3::pyfunction]
    pub fn build_demo_tree(birth_power: f32) -> GeometryData {
        let config_str = include_str!("../../../assets/tree_config.toml");
        let config: TreeConfig = toml::from_str(config_str).unwrap();

        return GeometryData::build_from_config(&config, 42);
    }

    #[pyo3::pyfunction]
    pub fn demo_mesh() -> String {
        let mut encoder = TreeEncoder::new();

        // Simple cylinder: 8 splines, 8 samples each = 64 vertices
        let pos_scale = 256i32;
        let t_scale = 1024i32;
        let n_splines = 8;
        let n_samples = 8;
        let ncp = 4;

        // z positions and radii
        let cp_z = [-0.5f32, -1.0 / 6.0, 1.0 / 6.0, 0.5];
        let cp_r = [0.5f32, 0.25, 0.25, 0.5];

        // t values
        let mut spine_t = Vec::with_capacity(n_samples);
        for i in 0..n_samples {
            spine_t.push(
                ((i as f32 / (n_samples - 1) as f32) * (ncp - 1) as f32 * t_scale as f32) as i32,
            );
        }
        encoder.add_immediate_buffer("spine_t", &spine_t);

        // Generate splines
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

        // Build expressions for vertices, normals, colors
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

        // Interleave splines for ring-major layout
        let vertices = crate::export::Expr::interleave(vert_exprs);
        let normals = crate::export::Expr::interleave(norm_exprs);
        let colors = crate::export::Expr::interleave(color_exprs);

        encoder.set_vertices(vertices);
        encoder.set_normals(normals);
        encoder.set_colors(colors);

        // Triangle indices (ring-major layout)
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

        // Delta encode indices
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

        // Debug layer: skeleton lines
        let skel_start_x = vec![0i32, 100, 200, 50, 150];
        let skel_start_y = vec![0i32, 50, 50, 100, 100];
        let skel_start_z = vec![0i32, 80, 80, 120, 120];
        let skel_end_x = vec![100i32, 200, 150, 150, 200];
        let skel_end_y = vec![50i32, 50, 100, 100, 100];
        let skel_end_z = vec![80i32, 80, 120, 120, 120];
        let skel_color_r = vec![0i32, 128, 255, 128, 255];
        let skel_color_g = vec![200i32, 100, 0, 150, 50];
        let skel_color_b = vec![100i32, 200, 100, 50, 200];

        encoder.add_immediate_buffer("skel_start_x", &skel_start_x);
        encoder.add_immediate_buffer("skel_start_y", &skel_start_y);
        encoder.add_immediate_buffer("skel_start_z", &skel_start_z);
        encoder.add_immediate_buffer("skel_end_x", &skel_end_x);
        encoder.add_immediate_buffer("skel_end_y", &skel_end_y);
        encoder.add_immediate_buffer("skel_end_z", &skel_end_z);
        encoder.add_immediate_buffer("skel_color_r", &skel_color_r);
        encoder.add_immediate_buffer("skel_color_g", &skel_color_g);
        encoder.add_immediate_buffer("skel_color_b", &skel_color_b);

        let debug_geom = crate::DebugGeometry {
            lines: Vec::new(),
            circles: Vec::new(),
            points: Vec::new(),
        };
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
}
