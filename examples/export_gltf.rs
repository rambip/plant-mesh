use mesh_tools::{Gltf, GltfBuilder};
use plant_core::meshing::GeometryData;
use plant_core::TreeConfig;
use std::collections::HashMap;
use std::fs;

fn main() {
    // 1. Load config from assets/tree_config.toml
    let config_str = fs::read_to_string("assets/tree_config.toml")
        .expect("Failed to read assets/tree_config.toml");
    let config: TreeConfig = toml::from_str(&config_str).expect("Failed to parse tree_config.toml");

    // 2. Generate the tree geometry
    let geometry = GeometryData::build_from_config(&config, 42);
    println!(
        "Generated geometry: {} points, {} triangles",
        geometry.points.len(),
        geometry.triangles.len() / 3
    );

    // 3. Flatten geometry data
    let positions: Vec<f32> = geometry
        .points
        .iter()
        .flat_map(|p| [p.x, p.y, p.z])
        .collect();

    let normals: Vec<f32> = geometry
        .normals
        .iter()
        .flat_map(|n| [n.x, n.y, n.z])
        .collect();

    let colors: Vec<f32> = geometry
        .colors
        .iter()
        .flat_map(|c| [c[0], c[1], c[2], c[3]])
        .collect();

    let indices = geometry.triangles;

    // Convert to bytes
    let pos_bytes: Vec<u8> = positions.iter().flat_map(|f| f.to_le_bytes()).collect();
    let norm_bytes: Vec<u8> = normals.iter().flat_map(|f| f.to_le_bytes()).collect();
    let color_bytes: Vec<u8> = colors.iter().flat_map(|f| f.to_le_bytes()).collect();
    let idx_bytes: Vec<u8> = indices.iter().flat_map(|i| i.to_le_bytes()).collect();

    // Combine into one buffer
    let mut buffer_data = pos_bytes.clone();
    let normals_offset = buffer_data.len();
    buffer_data.extend(&norm_bytes);
    let colors_offset = buffer_data.len();
    buffer_data.extend(&color_bytes);
    let indices_offset = buffer_data.len();
    buffer_data.extend(&idx_bytes);

    // Create buffer
    let buffer = mesh_tools::Buffer {
        byte_length: buffer_data.len(),
        uri: None,
    };

    // Buffer views
    let pos_buffer_view = mesh_tools::BufferView {
        buffer: 0,
        byte_offset: 0,
        byte_length: pos_bytes.len(),
        byte_stride: None,
        target: Some(34962), // ARRAY_BUFFER
    };

    let norm_buffer_view = mesh_tools::BufferView {
        buffer: 0,
        byte_offset: normals_offset,
        byte_length: norm_bytes.len(),
        byte_stride: None,
        target: Some(34962),
    };

    let color_buffer_view = mesh_tools::BufferView {
        buffer: 0,
        byte_offset: colors_offset,
        byte_length: color_bytes.len(),
        byte_stride: None,
        target: Some(34962),
    };

    let idx_buffer_view = mesh_tools::BufferView {
        buffer: 0,
        byte_offset: indices_offset,
        byte_length: idx_bytes.len(),
        byte_stride: None,
        target: Some(34963), // ELEMENT_ARRAY_BUFFER
    };

    // Calculate min/max for positions
    let min_x = positions
        .iter()
        .step_by(3)
        .cloned()
        .fold(f32::INFINITY, f32::min);
    let min_y = positions
        .iter()
        .skip(1)
        .step_by(3)
        .cloned()
        .fold(f32::INFINITY, f32::min);
    let min_z = positions
        .iter()
        .skip(2)
        .step_by(3)
        .cloned()
        .fold(f32::INFINITY, f32::min);
    let max_x = positions
        .iter()
        .step_by(3)
        .cloned()
        .fold(f32::NEG_INFINITY, f32::max);
    let max_y = positions
        .iter()
        .skip(1)
        .step_by(3)
        .cloned()
        .fold(f32::NEG_INFINITY, f32::max);
    let max_z = positions
        .iter()
        .skip(2)
        .step_by(3)
        .cloned()
        .fold(f32::NEG_INFINITY, f32::max);

    // Accessors
    let pos_accessor = mesh_tools::Accessor {
        buffer_view: 0,
        byte_offset: Some(0),
        component_type: 5126, // FLOAT
        count: positions.len() / 3,
        type_: "VEC3".to_string(),
        min: Some(vec![min_x, min_y, min_z]),
        max: Some(vec![max_x, max_y, max_z]),
        normalized: None,
    };

    let norm_accessor = mesh_tools::Accessor {
        buffer_view: 1,
        byte_offset: Some(0),
        component_type: 5126, // FLOAT
        count: normals.len() / 3,
        type_: "VEC3".to_string(),
        min: None,
        max: None,
        normalized: None,
    };

    let color_accessor = mesh_tools::Accessor {
        buffer_view: 2,
        byte_offset: Some(0),
        component_type: 5126, // FLOAT
        count: colors.len() / 4,
        type_: "VEC4".to_string(),
        min: None,
        max: None,
        normalized: None,
    };

    let idx_accessor = mesh_tools::Accessor {
        buffer_view: 3,
        byte_offset: Some(0),
        component_type: 5125, // UNSIGNED_INT
        count: indices.len(),
        type_: "SCALAR".to_string(),
        min: None,
        max: None,
        normalized: None,
    };

    let primitive = mesh_tools::Primitive {
        attributes: {
            let mut attrs = HashMap::new();
            attrs.insert("POSITION".to_string(), 0);
            attrs.insert("NORMAL".to_string(), 1);
            attrs.insert("COLOR_0".to_string(), 2);
            attrs
        },
        indices: Some(3),
        material: None,
        mode: Some(4), // TRIANGLES
    };

    let mesh = mesh_tools::Mesh {
        primitives: vec![primitive],
        name: None,
    };

    let node = mesh_tools::Node {
        name: None,
        mesh: Some(0),
        translation: None,
        rotation: None,
        scale: None,
        matrix: None,
        children: None,
    };

    let scene = mesh_tools::Scene {
        nodes: Some(vec![0]),
        name: None,
    };

    let asset = mesh_tools::Asset {
        version: "2.0".to_string(),
        generator: Some("Custom Rust Generator".to_string()),
        copyright: None,
    };

    let gltf = Gltf {
        asset,
        scene: Some(0),
        scenes: Some(vec![scene]),
        nodes: Some(vec![node]),
        meshes: Some(vec![mesh]),
        accessors: Some(vec![
            pos_accessor,
            norm_accessor,
            color_accessor,
            idx_accessor,
        ]),
        buffer_views: Some(vec![
            pos_buffer_view,
            norm_buffer_view,
            color_buffer_view,
            idx_buffer_view,
        ]),
        buffers: Some(vec![buffer]),
        materials: None,
        textures: None,
        images: None,
        samplers: None,
        animations: None,
        extensions: None,
        extensions_used: None,
        extensions_required: None,
        extras: None,
    };

    let glb_builder = GltfBuilder {
        gltf,
        buffer_data: buffer_data,
    };

    glb_builder
        .export_glb("tree.glb")
        .expect("Failed to generate GLB file");

    println!("Successfully exported tree.glb with normals and colors");
}
