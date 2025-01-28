
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
    //var intensity: f32 = 1.+dot(vertex.normal, vec3(0., 0.5, 0.5));
    // FIXME: always same z position
//    var intensity: f32 = 0.1*vertex_output.clip_position.z;
    //vertex_output.color = intensity*vertex.color;
    vertex_output.color = vec4(0.1*clip.z, clip.w, 0., 1.);
    return vertex_output;
}

// The fragment shader entry point.
@fragment
fn fragment(vertex_output: VertexOutput) -> @location(0) vec4<f32> {
    return vertex_output.color;
}
