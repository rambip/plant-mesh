use bevy::prelude::Mesh;

pub mod algorithms {
    pub use tubulin_core::meshing::algorithms::*;
}
pub mod particles {
    pub use tubulin_core::meshing::particles::*;
}

pub struct BevyMesh(pub Mesh);

impl BevyMesh {
    pub fn from_geometry_data(mut geometry: tubulin_core::GeometryData) -> Self {
        geometry.compute_smooth_normals();

        let mesh = Mesh::new(
            bevy::render::mesh::PrimitiveTopology::TriangleList,
            bevy::asset::RenderAssetUsages::default(),
        )
        .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, geometry.points.clone())
        .with_inserted_attribute(Mesh::ATTRIBUTE_NORMAL, geometry.normals.clone())
        .with_inserted_attribute(Mesh::ATTRIBUTE_COLOR, geometry.colors.clone())
        .with_inserted_indices(bevy::render::mesh::Indices::U32(geometry.triangles.clone()));
        BevyMesh(mesh)
    }

    pub fn from_geometry(
        volumetric: &tubulin_core::VolumetricTree,
        config: &tubulin_core::MeshConfig,
        rng: rand::rngs::StdRng,
    ) -> Self {
        let mut cache = volumetric.build_mesh(config, rng);
        cache.compute_smooth_normals();

        let mesh = Mesh::new(
            bevy::render::mesh::PrimitiveTopology::TriangleList,
            bevy::asset::RenderAssetUsages::default(),
        )
        .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, cache.points.clone())
        .with_inserted_attribute(Mesh::ATTRIBUTE_NORMAL, cache.normals.clone())
        .with_inserted_attribute(Mesh::ATTRIBUTE_COLOR, cache.colors.clone())
        .with_inserted_indices(bevy::render::mesh::Indices::U32(cache.triangles.clone()));
        BevyMesh(mesh)
    }
}
