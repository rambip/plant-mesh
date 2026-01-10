use bevy::prelude::Mesh;
use plant_core::meshing::{VolumetricTree, GeometryData, MeshConfig};
use plant_core::TreePipelinePhase;

pub mod algorithms {
    pub use plant_core::meshing::algorithms::*;
}
pub mod particles {
    pub use plant_core::meshing::particles::*;
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
        let builder = GeometryData::generate_from(prev, config, builder);
        
        // Convert GeometryData to Bevy Mesh
        let mut mesh = Mesh::new(
            bevy::render::mesh::PrimitiveTopology::TriangleList,
            bevy::asset::RenderAssetUsages::default(),
        )
        .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, builder.points.clone())
        .with_inserted_attribute(
            Mesh::ATTRIBUTE_COLOR,
            builder.colors.clone()
        )
        .with_inserted_indices(bevy::render::mesh::Indices::U32(
            builder.triangles.clone(),
        ));
        mesh.compute_smooth_normals();
        BevyMesh(mesh)
    }
}
