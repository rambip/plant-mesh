pub mod meshing;
pub mod growing;

pub use meshing::StrandsConfig;
pub use meshing::SplineIndex;
pub use meshing::GeometryData;
pub use meshing::mesh_builder::MeshDebugFlags;
pub use meshing::mesh_builder::MeshConfig;
pub use meshing::particles::TrajectoryBuilder;
pub use growing::generation::GrowConfig;
pub use growing::TreeSkeleton;
pub use growing::PlantNode;
pub use growing::PlantNodeProps;
pub use growing::NodeInfo;

pub trait VisualDebug {
    type Flags;
    #[cfg(feature = "bevy")]
    fn debug(&self, gizmos: &mut bevy_gizmos::prelude::Gizmos, debug_flags: Self::Flags);
}

pub trait TreePipelinePhase {
    type Previous;
    type Config;
    type Builder;
    fn generate_from(prev: Self::Previous, config: &Self::Config, builder: &mut Self::Builder) -> Self;
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
