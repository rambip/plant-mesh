use bevy_gizmos::gizmos::Gizmos;
use rand::rngs::StdRng;
use serde::{Deserialize, Serialize};

pub trait VisualDebug {
    type Flags: Copy + std::fmt::Debug + Default;
    fn debug(&self, gizmos: &mut Gizmos, debug_flags: Self::Flags);
}

impl VisualDebug for () {
    type Flags = ();
    fn debug(&self, _: &mut Gizmos, (): ()) {}
}

impl<T> VisualDebug for Option<&T>
where
    T: VisualDebug,
{
    type Flags = T::Flags;
    fn debug(&self, gizmos: &mut Gizmos, debug_flags: Self::Flags) {
        self.as_ref().map(|x| x.debug(gizmos, debug_flags));
    }
}

impl VisualDebug for StdRng {
    type Flags = ();
    fn debug(&self, _: &mut Gizmos, (): ()) {}
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

pub mod growing;
pub mod meshing;
mod tools;

pub use growing::{GrowConfig, PlantNode, Seed, TreeSkeleton, TreeSkeletonDebugData};
pub use meshing::{
    particles::spread_points, GeometryData, MeshConfig, MeshDebugFlags, StrandsConfig,
    TrajectoryBuilder, VolumetricTree,
};
