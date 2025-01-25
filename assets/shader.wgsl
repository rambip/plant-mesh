//! Very simple shader used to demonstrate how to get the world position and pass data
//! between the vertex and fragment shader. Also shows the custom vertex layout.

struct ColorGrading {
    balance: mat3x3<f32>,
    saturation: vec3<f32>,
    contrast: vec3<f32>,
    gamma: vec3<f32>,
    gain: vec3<f32>,
    lift: vec3<f32>,
    midtone_range: vec2<f32>,
    exposure: f32,
    hue: f32,
    post_saturation: f32,
}
struct View {
    clip_from_world: mat4x4<f32>,
    unjittered_clip_from_world: mat4x4<f32>,
    world_from_clip: mat4x4<f32>,
    world_from_view: mat4x4<f32>,
    view_from_world: mat4x4<f32>,
    clip_from_view: mat4x4<f32>,
    view_from_clip: mat4x4<f32>,
    world_position: vec3<f32>,
    exposure: f32,
    // viewport(x_origin, y_origin, width, height)
    viewport: vec4<f32>,
    frustum: array<vec4<f32>, 6>,
    color_grading: ColorGrading,
    mip_bias: f32,
};

struct Mesh {
    // Affine 4x3 matrices transposed to 3x4
    // Use bevy_render::maths::affine3_to_square to unpack
    world_from_local: mat3x4<f32>,
    previous_world_from_local: mat3x4<f32>,
    // 3x3 matrix packed in mat2x4 and f32 as:
    // [0].xyz, [1].x,
    // [1].yz, [2].xy
    // [2].z
    // Use bevy_pbr::mesh_functions::mat2x4_f32_to_mat3x3_unpack to unpack
    local_from_world_transpose_a: mat2x4<f32>,
    local_from_world_transpose_b: f32,
    // 'flags' is a bit field indicating various options. u32 is 32 bits so we have up to 32 options.
    flags: u32,
    lightmap_uv_rect: vec2<u32>,
    // The index of the mesh's first vertex in the vertex buffer.
    first_vertex_index: u32,
    pad_a: u32,
    pad_b: u32,
    pad_c: u32,
};

@group(0) @binding(0) var<uniform> view: View;
#ifdef PER_OBJECT_BUFFER_BATCH_SIZE
@group(1) @binding(0) var<uniform> mesh: array<Mesh, #{PER_OBJECT_BUFFER_BATCH_SIZE}u>;
#else
@group(1) @binding(0) var<storage> mesh: array<Mesh>;
#endif

fn affine3_to_square(affine: mat3x4<f32>) -> mat4x4<f32> {
    return transpose(mat4x4<f32>(
        affine[0],
        affine[1],
        affine[2],
        vec4<f32>(0.0, 0.0, 0.0, 1.0),
    ));
}

fn position_world_to_clip(world_pos: vec3<f32>) -> vec4<f32> {
    let clip_pos = view.clip_from_world * vec4(world_pos, 1.0);
    return clip_pos;
}


fn get_world_from_local(instance_index: u32) -> mat4x4<f32> {
    return affine3_to_square(mesh[instance_index].world_from_local);
}
fn mesh_position_local_to_world(world_from_local: mat4x4<f32>, vertex_position: vec4<f32>) -> vec4<f32> {
    return world_from_local * vertex_position;
}

struct Vertex {
    // This is needed if you are using batching and/or gpu preprocessing
    // It's a built in so you don't need to define it in the vertex layout
    @builtin(instance_index) instance_index: u32,
    // Like we defined for the vertex layout
    // position is at location 0
    @location(0) position: vec3<f32>,
    // and color at location 1
    @location(1) color: vec4<f32>,
    @location(2) normal: vec4<f32>,
};

// This is the output of the vertex shader and we also use it as the input for the fragment shader
struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) world_position: vec4<f32>,
    @location(1) color: vec3<f32>,
};


@vertex
fn vertex(vertex: Vertex) -> VertexOutput {
    var out: VertexOutput;
    // This is how bevy computes the world position
    // The vertex.instance_index is very important. Especially if you are using batching and gpu preprocessing
    var world_from_local = get_world_from_local(vertex.instance_index);
    out.world_position = mesh_position_local_to_world(world_from_local, vec4(vertex.position, 1.0));
    out.clip_position = position_world_to_clip(out.world_position.xyz);

    // We just use the raw vertex color
    out.color = vertex.color.rgb;

    return out;
}

@fragment
fn fragment(in: VertexOutput) -> @location(0) vec4<f32> {
    // output the color directly
    return vec4(in.color, 1.0);
}
