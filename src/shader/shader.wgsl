
#import bevy_render::view::View

@group(0) @binding(0) var<uniform> view: View;

// The GPU-side vertex structure.
struct Vertex {
    // The world-space position of the vertex.
    @location(0) position: vec3<f32>,
    // The color of the vertex.
    @location(1) normal: vec3<f32>,
    @location(2) color: vec4<f32>,
};

// Information passed from the vertex shader to the fragment shader.
struct VertexOutput {
    // The clip-space position of the vertex.
    @builtin(position) clip_position: vec4<f32>,
    // The color of the vertex.
    @location(0) color: vec4<f32>,
};

// The vertex shader entry point.
@vertex
fn vertex(vertex: Vertex) -> VertexOutput {
    // Use an orthographic projection.
    var vertex_output: VertexOutput;
    let clip = view.clip_from_world * vec4(vertex.position.xyz, 1.0);
    vertex_output.clip_position = clip;

    // hardcoded light direction
    var diffuse_light: f32 = 0.2*dot(vertex.normal, vec3(0., 0.6, 0.6));
    var intensity : f32 = 0.5+diffuse_light;
    vertex_output.color = vec4(intensity*vertex.color.rgb, vertex.color.a);
    return vertex_output;
}

// The fragment shader entry point.
@fragment
fn fragment(vertex_output: VertexOutput) -> @location(0) vec4<f32> {
    return vertex_output.color;
}
