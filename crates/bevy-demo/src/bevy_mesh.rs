use bevy::prelude::Mesh;
use tubulin_core::meshing::{GeometryData, MeshConfig, VolumetricTree};
use tubulin_core::TreePipelinePhase;

pub mod algorithms {
    pub use tubulin_core::meshing::algorithms::*;
}
pub mod particles {
    pub use tubulin_core::meshing::particles::*;
}

pub struct BevyMesh(pub Mesh);

impl TreePipelinePhase for BevyMesh {
    type Previous = VolumetricTree;
    type Config = MeshConfig;
    type Builder = GeometryData;
    fn generate_from(
        prev: Self::Previous,
        config: &Self::Config,
        builder: &mut Self::Builder,
    ) -> Self {
        *builder = GeometryData::generate_from(prev, config, builder);
        builder.compute_smooth_normals();

        // Convert GeometryData to Bevy Mesh
        let mesh = Mesh::new(
            bevy::render::mesh::PrimitiveTopology::TriangleList,
            bevy::asset::RenderAssetUsages::default(),
        )
        .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, builder.points.clone())
        .with_inserted_attribute(Mesh::ATTRIBUTE_NORMAL, builder.normals.clone())
        .with_inserted_attribute(Mesh::ATTRIBUTE_COLOR, builder.colors.clone())
        .with_inserted_indices(bevy::render::mesh::Indices::U32(builder.triangles.clone()));
        BevyMesh(mesh)
    }
}
