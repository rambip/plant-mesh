use bevy::prelude::Resource;
use bevy_gizmos::gizmos::Gizmos;
use rand::rngs::StdRng;
use serde::{Deserialize, Serialize};

#[derive(Copy, Clone, Default, Debug, Resource)]
pub struct DebugFlags {
    pub triangles: bool,
    pub strands: bool,
    pub skeleton: bool,
    pub other: bool,
    pub contours: bool,
}

impl DebugFlags {
    pub fn need_render_mesh(&self) -> bool {
        self.contours || self.triangles || self.other
    }
}

pub trait VisualDebug {
    fn debug(&self, gizmos: &mut Gizmos, debug_flags: DebugFlags);
}

impl VisualDebug for () {
    fn debug(&self, _: &mut Gizmos, _: DebugFlags) {}
}

impl<T> VisualDebug for Option<&T>
where
    T: VisualDebug,
{
    fn debug(&self, gizmos: &mut Gizmos, debug_flags: DebugFlags) {
        self.as_ref().map(|x| x.debug(gizmos, debug_flags));
    }
}

impl VisualDebug for StdRng {
    fn debug(&self, _: &mut Gizmos, _: DebugFlags) {}
}

pub trait TreePipelinePhase {
    type Previous;
    type Config: Copy + Serialize + Deserialize<'static>;
    type Builder: VisualDebug + From<StdRng>;

    fn generate_from(
        prev: Self::Previous,
        config: &Self::Config,
        builder: &mut Self::Builder,
    ) -> Self;
}

pub trait Grow {
    fn grow<Next>(self, config: &Next::Config, cache: &mut Next::Builder) -> Next
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

mod tools;
pub mod growing;
pub mod meshing;
pub mod shader;

pub use growing::{
    GrowConfig,
    Seed,
    TreeSkeleton,
    PlantNode,
    TreeSkeletonDebugData
};
pub use meshing::{
    StrandsConfig, MeshConfig, GeometryData,  VolumetricTree, TrajectoryBuilder,
    particles::spread_points
};
pub use shader::{CustomEntity, CustomMeshPipelinePlugin};
