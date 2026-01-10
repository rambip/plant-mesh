pub mod growing;
pub mod meshing;
pub mod utils;

pub use growing::generation::GrowConfig;
pub use growing::NodeInfo;
pub use growing::PlantNode;
pub use growing::PlantNodeProps;
pub use growing::TreeSkeleton;
pub use meshing::mesh_builder::MeshConfig;
pub use meshing::mesh_builder::MeshDebugFlags;
pub use meshing::particles::TrajectoryBuilder;
pub use meshing::GeometryData;
pub use meshing::SplineIndex;
pub use meshing::StrandsConfig;

use meshing::VolumetricTree;
use rand::SeedableRng;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone)]
#[cfg_attr(feature = "bevy", derive(bevy::prelude::TypePath, bevy::prelude::Asset))]
pub struct TreeConfig {
    pub grow: GrowConfig,
    pub strands: StrandsConfig,
    pub mesh: MeshConfig,
}

impl GeometryData {
    pub fn build_from_config(config: &TreeConfig, seed: u64) -> Self {
        let rng = rand::rngs::StdRng::seed_from_u64(seed);

        let mut plant_builder = rng.clone();
        let mut particle_builder = TrajectoryBuilder::new(rng.clone());
        let mut mesh_builder = GeometryData::new(rng.clone());

        let mut geometry = Seed
            .grow::<PlantNode>(&config.grow, &mut plant_builder)
            .grow::<TreeSkeleton>(&(), &mut ())
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

pub trait VisualDebug {
    type Flags;
    #[cfg(feature = "bevy")]
    fn debug(&self, gizmos: &mut bevy_gizmos::prelude::Gizmos, debug_flags: Self::Flags);
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
    type Builder = ();
    fn generate_from(prev: Self::Previous, _: &Self::Config, _: &mut Self::Builder) -> Self {
        let mut node_props = Vec::new();
        prev.register_node_properties(&mut node_props);
        let mut node_info = Vec::new();
        prev.register_node_info(&mut node_info, 0);

        TreeSkeleton {
            node_info,
            node_props,
        }
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
