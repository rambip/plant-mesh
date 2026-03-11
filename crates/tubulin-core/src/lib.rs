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
use meshing::VolumetricTree;
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

pub use export::{DebugLayer, DebugLayers, DebugLines, DebugPoints, TreeMeshExporter};

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
        let mut skeleton_builder = TreeSkeletonDebugData::new();
        let mut particle_builder = TrajectoryBuilder::new(rng.clone());
        let mut mesh_builder = GeometryData::new(rng.clone());

        let mut geometry = Seed
            .grow::<PlantNode>(&config.grow, &mut plant_builder)
            .grow::<TreeSkeleton>(&(), &mut skeleton_builder)
            .grow::<VolumetricTree>(&config.strands, &mut particle_builder)
            .grow::<GeometryData>(&config.mesh, &mut mesh_builder);

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

impl TreePipelinePhase for PlantNode {
    type Previous = Seed;
    type Config = GrowConfig;
    type Builder = rand::rngs::StdRng;
    fn generate_from(_: Self::Previous, config: &Self::Config, rng: &mut Self::Builder) -> Self {
        let root = PlantNodeProps {
            radius: config.base_radius,
            orientation: glam::Quat::IDENTITY,
            position: glam::Vec3::ZERO,
        };
        growing::generation::grow_tree_basic(config, rng, root, 0)
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
mod plant_core {
    use crate::{GeometryData, TreeConfig};
    #[pyo3::pyfunction]
    pub fn build_demo_tree(birth_power: f32) -> GeometryData {
        let config_str = include_str!("../../../assets/tree_config.toml");
        let mut config: TreeConfig = toml::from_str(config_str).unwrap();
        // config.grow.birth_power = birth_power;

        return GeometryData::build_from_config(&config, 42);
    }
}
