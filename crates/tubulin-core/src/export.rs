use serde::{ser::SerializeStruct, Serializer};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Buffer {
    pub k: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub offset: Option<i32>,
    pub data: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub length: Option<usize>,
}

#[derive(Debug, Clone)]
pub enum Expr {
    Var(String),
    Cumsum(Box<Expr>),
    Divp2(Box<Expr>, u32),
    Vec3(Box<Expr>, Box<Expr>, Box<Expr>),
    Vec4(Box<Expr>, Box<Expr>, Box<Expr>, Box<Expr>),
    Spline(Box<Expr>, Box<Expr>),
    Concat(Vec<Expr>),
    Interleave(Vec<Expr>),
    Triangle(Box<Expr>, Box<Expr>, Box<Expr>),
}

impl Serialize for Expr {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        if let Expr::Var(name) = self {
            return s.serialize_str(name);
        }
        let mut st = s.serialize_struct("Expr", 2)?;
        match self {
            Expr::Var(_) => unreachable!(),
            Expr::Cumsum(e) => {
                st.serialize_field("op", "cumsum")?;
                st.serialize_field("args", &[e])?;
            }
            Expr::Divp2(e, n) => {
                st.serialize_field("op", "divp2")?;
                st.serialize_field("args", &(e, n))?; // tuple → array
            }
            Expr::Vec3(a, b, c) => {
                st.serialize_field("op", "vec3")?;
                st.serialize_field("args", &[a, b, c])?;
            }
            Expr::Vec4(a, b, c, d) => {
                st.serialize_field("op", "vec4")?;
                st.serialize_field("args", &[a, b, c, d])?;
            }
            Expr::Spline(a, b) => {
                st.serialize_field("op", "spline")?;
                st.serialize_field("args", &[a, b])?;
            }
            Expr::Concat(v) => {
                st.serialize_field("op", "concat")?;
                st.serialize_field("args", v)?;
            }
            Expr::Interleave(v) => {
                st.serialize_field("op", "interleave")?;
                st.serialize_field("args", v)?;
            }
            Expr::Triangle(a, b, c) => {
                st.serialize_field("op", "triangle")?;
                st.serialize_field("args", &[a, b, c])?;
            }
        }
        st.end()
    }
}

impl Expr {
    pub fn var(name: &str) -> Self {
        Expr::Var(name.to_string())
    }

    pub fn cumsum(self) -> Self {
        Expr::Cumsum(Box::new(self))
    }

    pub fn divp2(self, exp: u32) -> Self {
        Expr::Divp2(Box::new(self), exp)
    }

    pub fn vec3(x: Self, y: Self, z: Self) -> Self {
        Expr::Vec3(Box::new(x), Box::new(y), Box::new(z))
    }

    pub fn vec4(r: Self, g: Self, b: Self, a: Self) -> Self {
        Expr::Vec4(Box::new(r), Box::new(g), Box::new(b), Box::new(a))
    }

    pub fn spline(points: Self, times: Self) -> Self {
        Expr::Spline(Box::new(points), Box::new(times))
    }

    pub fn concat(args: Vec<Self>) -> Self {
        Expr::Concat(args)
    }

    pub fn interleave(args: Vec<Self>) -> Self {
        Expr::Interleave(args)
    }

    pub fn triangle(i: Self, j: Self, k: Self) -> Self {
        Expr::Triangle(Box::new(i), Box::new(j), Box::new(k))
    }
}

fn zigzag_encode(n: i32) -> u32 {
    if n >= 0 {
        (n * 2) as u32
    } else {
        ((-n) * 2 - 1) as u32
    }
}

pub fn rice_encode(values: &[i32], k: u32) -> Vec<u8> {
    let mut bits: Vec<u8> = Vec::new();

    for &signed in values {
        let value = zigzag_encode(signed);
        let quotient = value >> k;
        let remainder = value & ((1 << k) - 1);

        if quotient > 0 {
            bits.extend(std::iter::repeat_n(1, quotient as usize));
        }
        bits.push(0);
        for i in (0..k).rev() {
            bits.push(((remainder >> i) & 1) as u8);
        }
    }

    while !bits.len().is_multiple_of(8) {
        bits.push(0);
    }

    let mut bytes = vec![0u8; bits.len() / 8];
    for (i, chunk) in bits.chunks(8).enumerate() {
        let mut byte = 0u8;
        for &bit in chunk.iter() {
            byte = (byte << 1) | bit;
        }
        bytes[i] = byte;
    }

    bytes
}

pub fn encode_int32(values: &[i32]) -> String {
    let bytes: Vec<u8> = values.iter().flat_map(|&v| v.to_le_bytes()).collect();
    base64_encode(&bytes)
}

fn base64_encode(bytes: &[u8]) -> String {
    const CHARS: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut result = String::new();

    for chunk in bytes.chunks(3) {
        let b0 = chunk[0] as usize;
        let b1 = chunk.get(1).copied().unwrap_or(0) as usize;
        let b2 = chunk.get(2).copied().unwrap_or(0) as usize;

        result.push(CHARS[b0 >> 2] as char);
        result.push(CHARS[((b0 & 0x03) << 4) | (b1 >> 4)] as char);

        if chunk.len() > 1 {
            result.push(CHARS[((b1 & 0x0f) << 2) | (b2 >> 6)] as char);
        } else {
            result.push('=');
        }

        if chunk.len() > 2 {
            result.push(CHARS[b2 & 0x3f] as char);
        } else {
            result.push('=');
        }
    }

    result
}

pub struct TreeEncoder {
    pub buffers: HashMap<String, Buffer>,
    pub outputs: Outputs,
    pub debug: Option<DebugLayers>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Outputs {
    pub vertices: Expr,
    pub normals: Expr,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub colors: Option<Expr>,
    pub triangles: Expr,
}

#[derive(Debug, Clone, Serialize)]
pub struct DebugLayers {
    #[serde(flatten)]
    pub layers: HashMap<String, DebugLayer>,
}

#[derive(Debug, Clone, Serialize)]
pub struct DebugLayer {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub points: Option<DebugPoints>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lines: Option<DebugLines>,
}

#[derive(Debug, Clone, Serialize)]
pub struct DebugPoints {
    pub positions: Expr,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub colors: Option<Expr>,
}

#[derive(Debug, Clone, Serialize)]
pub struct DebugLines {
    pub starts: Expr,
    pub ends: Expr,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub colors: Option<Expr>,
}

impl TreeEncoder {
    pub fn new() -> Self {
        Self {
            buffers: HashMap::from([(
                "empty".to_string(),
                Buffer {
                    k: 0,
                    offset: Some(0),
                    length: Some(0),
                    data: String::new(),
                },
            )]),
            outputs: Outputs {
                vertices: Expr::var("empty"),
                normals: Expr::var("empty"),
                colors: None,
                triangles: Expr::var("empty"),
            },
            debug: None,
        }
    }

    pub fn add_buffer(&mut self, name: &str, k: u32, data: String, length: usize) {
        self.buffers.insert(
            name.to_string(),
            Buffer {
                k,
                offset: None,
                data,
                length: Some(length),
            },
        );
    }

    pub fn add_delta_buffer(&mut self, name: &str, values: &[i32], k: u32) {
        let encoded = rice_encode(values, k);
        let data = base64_encode(&encoded);
        self.add_buffer(name, k, data, values.len());
    }

    pub fn add_immediate_buffer(&mut self, name: &str, values: &[i32]) {
        let data = encode_int32(values);
        self.add_buffer(name, 0, data, values.len());
    }

    pub fn set_vertices(&mut self, expr: Expr) {
        self.outputs.vertices = expr;
    }

    pub fn set_normals(&mut self, expr: Expr) {
        self.outputs.normals = expr;
    }

    pub fn set_colors(&mut self, expr: Expr) {
        self.outputs.colors = Some(expr);
    }

    pub fn set_triangles(&mut self, expr: Expr) {
        self.outputs.triangles = expr;
    }

    pub fn set_debug(&mut self, debug: DebugLayers) {
        self.debug = Some(debug);
    }
    pub fn add_debug_layer(
        &mut self,
        layer_name: &str,
        geometry: &crate::DebugGeometry,
    ) -> DebugLayer {
        let mut layer = DebugLayer {
            points: None,
            lines: None,
        };
        // Convert lines to debug lines
        if !geometry.lines.is_empty() {
            let n = geometry.lines.len();

            let mut starts_x = Vec::with_capacity(n);
            let mut starts_y = Vec::with_capacity(n);
            let mut starts_z = Vec::with_capacity(n);
            let mut ends_x = Vec::with_capacity(n);
            let mut ends_y = Vec::with_capacity(n);
            let mut ends_z = Vec::with_capacity(n);
            let mut colors_r = Vec::with_capacity(n);
            let mut colors_g = Vec::with_capacity(n);
            let mut colors_b = Vec::with_capacity(n);
            let mut colors_a = Vec::with_capacity(n);
            for (start, end, color) in &geometry.lines {
                starts_x.push((start.x * 256.0) as i32);
                starts_y.push((start.y * 256.0) as i32);
                starts_z.push((start.z * 256.0) as i32);
                ends_x.push((end.x * 256.0) as i32);
                ends_y.push((end.y * 256.0) as i32);
                ends_z.push((end.z * 256.0) as i32);
                colors_r.push((color.0[0] * 255.0) as i32);
                colors_g.push((color.0[1] * 255.0) as i32);
                colors_b.push((color.0[2] * 255.0) as i32);
                colors_a.push((color.0[3] * 255.0) as i32);
            }
            let prefix = format!("debug_{}_", layer_name);
            let sx = format!("{}start_x", prefix);
            let sy = format!("{}start_y", prefix);
            let sz = format!("{}start_z", prefix);
            let ex = format!("{}end_x", prefix);
            let ey = format!("{}end_y", prefix);
            let ez = format!("{}end_z", prefix);
            let cr = format!("{}color_r", prefix);
            let cg = format!("{}color_g", prefix);
            let cb = format!("{}color_b", prefix);
            let ca = format!("{}color_a", prefix);
            self.add_immediate_buffer(&sx, &starts_x);
            self.add_immediate_buffer(&sy, &starts_y);
            self.add_immediate_buffer(&sz, &starts_z);
            self.add_immediate_buffer(&ex, &ends_x);
            self.add_immediate_buffer(&ey, &ends_y);
            self.add_immediate_buffer(&ez, &ends_z);
            self.add_immediate_buffer(&cr, &colors_r);
            self.add_immediate_buffer(&cg, &colors_g);
            self.add_immediate_buffer(&cb, &colors_b);
            self.add_immediate_buffer(&ca, &colors_a);
            let starts = Expr::divp2(
                Expr::vec3(Expr::var(&sx), Expr::var(&sy), Expr::var(&sz)),
                8,
            );
            let ends = Expr::divp2(
                Expr::vec3(Expr::var(&ex), Expr::var(&ey), Expr::var(&ez)),
                8,
            );
            let colors = Expr::divp2(
                Expr::vec4(
                    Expr::var(&cr),
                    Expr::var(&cg),
                    Expr::var(&cb),
                    Expr::var(&ca),
                ),
                8,
            );
            layer.lines = Some(DebugLines {
                starts,
                ends,
                colors: Some(colors),
            });
        }
        // Convert points to debug points
        if !geometry.points.is_empty() {
            let n = geometry.points.len();

            let mut positions_x = Vec::with_capacity(n);
            let mut positions_y = Vec::with_capacity(n);
            let mut positions_z = Vec::with_capacity(n);
            let mut colors_r = Vec::with_capacity(n);
            let mut colors_g = Vec::with_capacity(n);
            let mut colors_b = Vec::with_capacity(n);
            let mut colors_a = Vec::with_capacity(n);
            for (pos, color) in &geometry.points {
                positions_x.push((pos.x * 256.0) as i32);
                positions_y.push((pos.y * 256.0) as i32);
                positions_z.push((pos.z * 256.0) as i32);
                colors_r.push((color.0[0] * 255.0) as i32);
                colors_g.push((color.0[1] * 255.0) as i32);
                colors_b.push((color.0[2] * 255.0) as i32);
                colors_a.push((color.0[3] * 255.0) as i32);
            }
            let prefix = format!("debug_{}_pt_", layer_name);
            let px = format!("{}x", prefix);
            let py = format!("{}y", prefix);
            let pz = format!("{}z", prefix);
            let cr = format!("{}color_r", prefix);
            let cg = format!("{}color_g", prefix);
            let cb = format!("{}color_b", prefix);
            let ca = format!("{}color_a", prefix);
            self.add_immediate_buffer(&px, &positions_x);
            self.add_immediate_buffer(&py, &positions_y);
            self.add_immediate_buffer(&pz, &positions_z);
            self.add_immediate_buffer(&cr, &colors_r);
            self.add_immediate_buffer(&cg, &colors_g);
            self.add_immediate_buffer(&cb, &colors_b);
            self.add_immediate_buffer(&ca, &colors_a);
            let positions = Expr::divp2(
                Expr::vec3(Expr::var(&px), Expr::var(&py), Expr::var(&pz)),
                8,
            );
            let colors = Expr::divp2(
                Expr::vec4(
                    Expr::var(&cr),
                    Expr::var(&cg),
                    Expr::var(&cb),
                    Expr::var(&ca),
                ),
                8,
            );
            layer.points = Some(DebugPoints {
                positions,
                colors: Some(colors),
            });
        }
        layer
    }
}

impl Default for TreeEncoder {
    fn default() -> Self {
        Self::new()
    }
}

fn serialize_expr(expr: &Expr) -> serde_json::Value {
    match expr {
        Expr::Var(name) => serde_json::Value::String(name.clone()),
        Expr::Cumsum(inner) => serde_json::json!({
            "op": "cumsum",
            "args": [serialize_expr(inner)]
        }),
        Expr::Divp2(inner, exp) => serde_json::json!({
            "op": "divp2",
            "args": [serialize_expr(inner), *exp]
        }),
        Expr::Vec3(x, y, z) => serde_json::json!({
            "op": "vec3",
            "args": [serialize_expr(x), serialize_expr(y), serialize_expr(z)]
        }),
        Expr::Vec4(r, g, b, a) => serde_json::json!({
            "op": "vec4",
            "args": [serialize_expr(r), serialize_expr(g), serialize_expr(b), serialize_expr(a)]
        }),
        Expr::Spline(points, times) => serde_json::json!({
            "op": "spline",
            "args": [serialize_expr(points), serialize_expr(times)]
        }),
        Expr::Concat(args) => serde_json::json!({
            "op": "concat",
            "args": args.iter().map(serialize_expr).collect::<Vec<_>>()
        }),
        Expr::Interleave(args) => serde_json::json!({
            "op": "interleave",
            "args": args.iter().map(serialize_expr).collect::<Vec<_>>()
        }),
        Expr::Triangle(i, j, k) => serde_json::json!({
            "op": "triangle",
            "args": [serialize_expr(i), serialize_expr(j), serialize_expr(k)]
        }),
    }
}

impl TreeEncoder {
    pub fn add_vec3_components(&mut self, prefix: &str, vectors: &[glam::Vec3], scale: i32) {
        let (x, y, z) = crate::utils::quantize_vec3_components(vectors, scale);
        self.add_immediate_buffer(&format!("{}_x", prefix), &x);
        self.add_immediate_buffer(&format!("{}_y", prefix), &y);
        self.add_immediate_buffer(&format!("{}_z", prefix), &z);
    }

    pub fn add_color_components(&mut self, prefix: &str, colors: &[[f32; 4]], scale: i32) {
        let (r, g, b, a) = crate::utils::quantize_color_components(colors, scale);
        self.add_immediate_buffer(&format!("{}_r", prefix), &r);
        self.add_immediate_buffer(&format!("{}_g", prefix), &g);
        self.add_immediate_buffer(&format!("{}_b", prefix), &b);
        self.add_immediate_buffer(&format!("{}_a", prefix), &a);
    }

    /// Creates scaled Vec3 expression from quantized components.
    ///
    /// Assumes buffers named `{prefix}_x`, `{prefix}_y`, `{prefix}_z` exist.
    pub fn make_scaled_vec3(&self, prefix: &str, scale_bits: u32) -> Expr {
        Expr::divp2(
            Expr::vec3(
                Expr::var(&format!("{}_x", prefix)),
                Expr::var(&format!("{}_y", prefix)),
                Expr::var(&format!("{}_z", prefix)),
            ),
            scale_bits,
        )
    }

    /// Creates scaled Vec4 expression from quantized components.
    ///
    /// Assumes buffers named `{prefix}_r`, `{prefix}_g`, `{prefix}_b`, `{prefix}_a` exist.
    pub fn make_scaled_vec4(&self, prefix: &str, scale_bits: u32) -> Expr {
        Expr::divp2(
            Expr::vec4(
                Expr::var(&format!("{}_r", prefix)),
                Expr::var(&format!("{}_g", prefix)),
                Expr::var(&format!("{}_b", prefix)),
                Expr::var(&format!("{}_a", prefix)),
            ),
            scale_bits,
        )
    }

    /// Helper to add geometry data and return debug geometry layer expression.
    pub fn add_geometry_with_debug(
        &mut self,
        geom: &crate::GeometryData,
        prefix: &str,
        scale: i32,
    ) -> DebugLayer {
        self.add_vec3_components(&format!("{}_v", prefix), &geom.points, scale);
        self.add_vec3_components(&format!("{}_n", prefix), &geom.normals, scale);
        self.add_color_components(&format!("{}_c", prefix), &geom.colors, scale);

        let indices: Vec<i32> = geom.triangles.iter().map(|&i| i as i32).collect();
        self.add_immediate_buffer(&format!("{}_indices", prefix), &indices);

        let debug_geom = crate::VisualDebug::debug_data(geom);
        self.add_debug_layer(prefix, &debug_geom)
    }

    /// Sets standard geometry outputs from component prefixes.
    ///
    /// Assumes buffers follow the naming convention: `{v_prefix}_x/y/z`, etc.
    pub fn set_geometry_outputs(
        &mut self,
        v_prefix: &str,
        n_prefix: &str,
        c_prefix: &str,
        scale_bits: u32,
    ) {
        self.set_vertices(self.make_scaled_vec3(v_prefix, scale_bits));
        self.set_normals(self.make_scaled_vec3(n_prefix, scale_bits));
        self.set_colors(self.make_scaled_vec4(c_prefix, scale_bits));
    }

    pub fn to_json(&self) -> String {
        let mut map = serde_json::Map::new();
        map.insert(
            "treemesh".to_string(),
            serde_json::Value::String("0.1".to_string()),
        );
        map.insert(
            "spline_convention".to_string(),
            serde_json::Value::String("reflect".to_string()),
        );
        map.insert(
            "buffers".to_string(),
            serde_json::to_value(&self.buffers).unwrap(),
        );

        let mut outputs_map = serde_json::Map::new();
        outputs_map.insert(
            "vertices".to_string(),
            serialize_expr(&self.outputs.vertices),
        );
        outputs_map.insert("normals".to_string(), serialize_expr(&self.outputs.normals));
        if let Some(colors) = &self.outputs.colors {
            outputs_map.insert("colors".to_string(), serialize_expr(colors));
        }
        outputs_map.insert(
            "triangles".to_string(),
            serialize_expr(&self.outputs.triangles),
        );
        map.insert(
            "outputs".to_string(),
            serde_json::Value::Object(outputs_map),
        );

        if let Some(debug) = &self.debug {
            let debug_json = serde_json::to_value(debug).unwrap();
            map.insert("debug".to_string(), debug_json);
        }

        serde_json::to_string_pretty(&serde_json::Value::Object(map)).unwrap()
    }
}

pub fn create_cylinder_mesh() -> String {
    let mut exporter = TreeEncoder::new();

    let pos_scale = 256i32;
    let t_scale = 1024i32;
    let color_scale = 256i32;

    let n_splines = 8;
    let n_samples = 8;
    let ncp1 = 4;

    let cp_z: Vec<f32> = vec![-0.5, -1.0 / 6.0, 1.0 / 6.0, 0.5];
    let cp_r: Vec<f32> = vec![0.5, 0.25, 0.25, 0.5];

    let mut spine_t = Vec::new();
    for i in 0..n_samples {
        spine_t.push(
            ((i as f32 / (n_samples - 1) as f32) * (ncp1 - 1) as f32 * t_scale as f32).round()
                as i32,
        );
    }
    exporter.add_immediate_buffer("spine_t", &spine_t);

    let mut spine_spline_exprs = Vec::new();
    let mut spine_nspline_exprs = Vec::new();
    let mut spine_colors_exprs = Vec::new();

    for s in 0..n_splines {
        let angle = (s as f32 / n_splines as f32) * std::f32::consts::PI * 2.0;
        let cos = angle.cos();
        let sin = angle.sin();

        let cpx: Vec<i32> = cp_r
            .iter()
            .map(|r| (cos * r * pos_scale as f32).round() as i32)
            .collect();
        let cpy: Vec<i32> = cp_r
            .iter()
            .map(|r| (sin * r * pos_scale as f32).round() as i32)
            .collect();
        let cpz: Vec<i32> = cp_z
            .iter()
            .map(|z| (z * pos_scale as f32).round() as i32)
            .collect();

        let cnx: Vec<i32> = cp_r
            .iter()
            .map(|_| (cos * pos_scale as f32).round() as i32)
            .collect();
        let cny: Vec<i32> = cp_r
            .iter()
            .map(|_| (sin * pos_scale as f32).round() as i32)
            .collect();
        let cnz: Vec<i32> = vec![0; ncp1];

        let ccr: Vec<i32> = cp_z
            .iter()
            .map(|z| ((z + 0.5) * color_scale as f32).round() as i32)
            .collect();
        let ccg: Vec<i32> = cp_z
            .iter()
            .map(|z| ((0.3 + 0.4 * (z + 0.5)) * color_scale as f32).round() as i32)
            .collect();
        let ccb: Vec<i32> = cp_z
            .iter()
            .map(|z| ((0.5 + 0.5 * (z + 0.5)) * color_scale as f32).round() as i32)
            .collect();
        let cca: Vec<i32> = vec![color_scale; ncp1];

        let delta = |arr: &[i32]| -> Vec<i32> {
            arr.iter()
                .enumerate()
                .map(|(i, &v)| if i == 0 { v } else { v - arr[i - 1] })
                .collect()
        };

        let name = format!("sp{}_", s);
        let xb = format!("{}x", name);
        let yb = format!("{}y", name);
        let zb = format!("{}z", name);
        let nxb = format!("{}nx", name);
        let nyb = format!("{}ny", name);
        let nzb = format!("{}nz", name);
        let crb = format!("{}cr", name);
        let cgb = format!("{}cg", name);
        let cbb = format!("{}cb", name);
        let cab = format!("{}ca", name);

        exporter.add_delta_buffer(&xb, &delta(&cpx), 2);
        exporter.add_delta_buffer(&yb, &delta(&cpy), 2);
        exporter.add_delta_buffer(&zb, &delta(&cpz), 2);
        exporter.add_immediate_buffer(&nxb, &cnx);
        exporter.add_immediate_buffer(&nyb, &cny);
        exporter.add_immediate_buffer(&nzb, &cnz);
        exporter.add_immediate_buffer(&crb, &ccr);
        exporter.add_immediate_buffer(&cgb, &ccg);
        exporter.add_immediate_buffer(&cbb, &ccb);
        exporter.add_immediate_buffer(&cab, &cca);

        let cp_vec3 = |xb: &str, yb: &str, zb: &str| -> Expr {
            Expr::divp2(
                Expr::vec3(
                    Expr::cumsum(Expr::var(xb)),
                    Expr::cumsum(Expr::var(yb)),
                    Expr::cumsum(Expr::var(zb)),
                ),
                8,
            )
        };

        let cp_vec4 = |rb: &str, gb: &str, bb: &str, ab: &str| -> Expr {
            Expr::divp2(
                Expr::vec4(
                    Expr::cumsum(Expr::var(rb)),
                    Expr::cumsum(Expr::var(gb)),
                    Expr::cumsum(Expr::var(bb)),
                    Expr::cumsum(Expr::var(ab)),
                ),
                8,
            )
        };

        let t_expr = Expr::divp2(Expr::var("spine_t"), 10);

        spine_spline_exprs.push(Expr::spline(cp_vec3(&xb, &yb, &zb), t_expr.clone()));
        spine_nspline_exprs.push(Expr::spline(cp_vec3(&nxb, &nyb, &nzb), t_expr.clone()));
        spine_colors_exprs.push(Expr::spline(
            cp_vec4(&crb, &cgb, &cbb, &cab),
            t_expr.clone(),
        ));
    }

    let _n_skel_lines = 5;
    let skel_start_x = vec![0i32, 100, 200, 50, 150];
    let skel_start_y = vec![0i32, 50, 50, 100, 100];
    let skel_start_z = vec![0i32, 80, 80, 120, 120];
    let skel_end_x = vec![100i32, 200, 150, 150, 200];
    let skel_end_y = vec![50i32, 50, 100, 100, 100];
    let skel_end_z = vec![80i32, 80, 120, 120, 120];
    let skel_color_r = vec![0i32, 128, 255, 128, 255];
    let skel_color_g = vec![200i32, 100, 0, 150, 50];
    let skel_color_b = vec![100i32, 200, 100, 50, 200];

    exporter.add_immediate_buffer("skel_start_x", &skel_start_x);
    exporter.add_immediate_buffer("skel_start_y", &skel_start_y);
    exporter.add_immediate_buffer("skel_start_z", &skel_start_z);
    exporter.add_immediate_buffer("skel_end_x", &skel_end_x);
    exporter.add_immediate_buffer("skel_end_y", &skel_end_y);
    exporter.add_immediate_buffer("skel_end_z", &skel_end_z);
    exporter.add_immediate_buffer("skel_color_r", &skel_color_r);
    exporter.add_immediate_buffer("skel_color_g", &skel_color_g);
    exporter.add_immediate_buffer("skel_color_b", &skel_color_b);

    exporter.add_immediate_buffer("empty", &[]);

    let vertices = Expr::concat(vec![Expr::interleave(spine_spline_exprs)]);

    let normals = Expr::concat(vec![Expr::interleave(spine_nspline_exprs)]);

    let colors = Expr::concat(vec![Expr::interleave(spine_colors_exprs)]);

    exporter.set_vertices(vertices);
    exporter.set_normals(normals);
    exporter.set_colors(colors);
    exporter.set_triangles(Expr::var("empty"));

    let mut debug_layers = DebugLayers {
        layers: HashMap::new(),
    };

    let skeleton_lines = DebugLines {
        starts: Expr::divp2(
            Expr::vec3(
                Expr::var("skel_start_x"),
                Expr::var("skel_start_y"),
                Expr::var("skel_start_z"),
            ),
            8,
        ),
        ends: Expr::divp2(
            Expr::vec3(
                Expr::var("skel_end_x"),
                Expr::var("skel_end_y"),
                Expr::var("skel_end_z"),
            ),
            8,
        ),
        colors: Some(Expr::divp2(
            Expr::vec3(
                Expr::var("skel_color_r"),
                Expr::var("skel_color_g"),
                Expr::var("skel_color_b"),
            ),
            8,
        )),
    };

    debug_layers.layers.insert(
        "skeleton".to_string(),
        DebugLayer {
            points: None,
            lines: Some(skeleton_lines),
        },
    );

    exporter.set_debug(debug_layers);

    exporter.to_json()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_zigzag_encode() {
        assert_eq!(zigzag_encode(0), 0);
        assert_eq!(zigzag_encode(1), 2);
        assert_eq!(zigzag_encode(-1), 1);
        assert_eq!(zigzag_encode(2), 4);
        assert_eq!(zigzag_encode(-2), 3);
    }
    #[test]
    fn test_rice_encode_roundtrip() {
        let values = vec![0i32, 1, 2, 3, -1, -2, 100, -100];
        let encoded = rice_encode(&values, 2);
        let decoded = decode_rice(&encoded, 2, values.len());
        assert_eq!(values, decoded);
    }
    fn decode_rice(data: &[u8], k: u32, length: usize) -> Vec<i32> {
        let mut result = Vec::new();
        let mut bit_buffer = 0u8;
        let mut bit_count = 0;
        let mut byte_idx = 0;
        let mut read_bit = || -> u8 {
            if bit_count == 0 {
                if byte_idx >= data.len() {
                    return 0;
                }
                bit_buffer = data[byte_idx];
                byte_idx += 1;
                bit_count = 8;
            }
            let bit = (bit_buffer >> 7) & 1;
            bit_buffer <<= 1;
            bit_count -= 1;
            bit
        };
        while result.len() < length {
            let mut quotient = 0u32;
            while read_bit() == 1 {
                quotient += 1;
            }
            let mut remainder = 0u32;
            for _ in 0..k {
                remainder = (remainder << 1) | read_bit() as u32;
            }
            let zigzag = (quotient << k) + remainder;
            let signed = if (zigzag & 1) == 0 {
                (zigzag >> 1) as i32
            } else {
                -((zigzag >> 1) as i32) - 1
            };
            result.push(signed);
        }
        result
    }
    #[test]
    fn test_base64_encode() {
        assert_eq!(base64_encode(&[0]), "AA==");
        assert_eq!(base64_encode(&[1]), "AQ==");
        assert_eq!(base64_encode(&[0, 1]), "AAE=");
        assert_eq!(base64_encode(&[0, 1, 2]), "AAEC");
    }
    #[test]
    fn test_expr_var() {
        let expr = Expr::var("test_buffer");
        let json = serialize_expr(&expr);
        assert_eq!(json, serde_json::Value::String("test_buffer".to_string()));
    }
    #[test]
    fn test_expr_cumsum() {
        let expr = Expr::cumsum(Expr::var("buffer"));
        let json = serialize_expr(&expr);
        let expected = serde_json::json!({
            "op": "cumsum",
            "args": ["buffer"]
        });
        assert_eq!(json, expected);
    }
    #[test]
    fn test_expr_divp2() {
        let expr = Expr::var("buffer").divp2(8);
        let json = serialize_expr(&expr);
        let expected = serde_json::json!({
            "op": "divp2",
            "args": ["buffer", 8]
        });
        assert_eq!(json, expected);
    }
    #[test]
    fn test_expr_vec3() {
        let expr = Expr::vec3(Expr::var("x"), Expr::var("y"), Expr::var("z"));
        let json = serialize_expr(&expr);
        let expected = serde_json::json!({
            "op": "vec3",
            "args": ["x", "y", "z"]
        });
        assert_eq!(json, expected);
    }
    #[test]
    fn test_expr_concat() {
        let expr = Expr::concat(vec![Expr::var("a"), Expr::var("b")]);
        let json = serialize_expr(&expr);
        let expected = serde_json::json!({
            "op": "concat",
            "args": ["a", "b"]
        });
        assert_eq!(json, expected);
    }
    #[test]
    fn test_tree_mesh_exporter() {
        let mut exporter = TreeEncoder::new();
        exporter.add_immediate_buffer("test", &[1, 2, 3]);
        exporter.set_vertices(Expr::var("test"));
        exporter.set_normals(Expr::var("test"));
        exporter.set_triangles(Expr::var("test"));

        let json = exporter.to_json();
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();

        assert_eq!(parsed["treemesh"], "0.1");
        assert!(parsed["buffers"]["test"].is_object());
    }
    #[test]
    fn test_tree_mesh_exporter_with_colors() {
        let mut exporter = TreeEncoder::new();
        exporter.add_immediate_buffer("verts", &[1, 2, 3]);
        exporter.add_immediate_buffer("colors", &[255, 0, 0, 255]);

        exporter.set_vertices(Expr::var("verts"));
        exporter.set_normals(Expr::var("verts"));
        exporter.set_colors(Expr::var("colors"));
        exporter.set_triangles(Expr::var("verts"));

        let json = exporter.to_json();
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();

        assert_eq!(parsed["outputs"]["colors"], "colors");
    }
    #[test]
    fn test_debug_layers() {
        let mut exporter = TreeEncoder::new();

        let mut debug_layers = DebugLayers {
            layers: HashMap::new(),
        };

        let lines = DebugLines {
            starts: Expr::vec3(Expr::var("sx"), Expr::var("sy"), Expr::var("sz")),
            ends: Expr::vec3(Expr::var("ex"), Expr::var("ey"), Expr::var("ez")),
            colors: None,
        };

        debug_layers.layers.insert(
            "skeleton".to_string(),
            DebugLayer {
                points: None,
                lines: Some(lines),
            },
        );

        exporter.set_debug(debug_layers);

        let json = exporter.to_json();
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();

        assert!(parsed["debug"]["skeleton"]["lines"].is_object());
    }
    #[test]
    fn test_create_cylinder_mesh() {
        let json = create_cylinder_mesh();
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();

        assert_eq!(parsed["treemesh"], "0.1");
        assert!(parsed["buffers"].is_object());
        assert!(parsed["outputs"]["vertices"].is_object());
        assert!(parsed["debug"].is_object());
    }
    #[test]
    fn test_add_debug_layer() {
        use crate::{DebugColor, DebugGeometry};

        let mut exporter = TreeEncoder::new();
        let mut geometry = DebugGeometry::new();

        // Add a line
        geometry.lines.push((
            glam::Vec3::new(0.0, 0.0, 0.0),
            glam::Vec3::new(1.0, 1.0, 1.0),
            DebugColor::rgb(1.0, 0.0, 0.0),
        ));

        // Add a point
        geometry.points.push((
            glam::Vec3::new(0.5, 0.5, 0.5),
            DebugColor::rgb(0.0, 1.0, 0.0),
        ));

        let layer = exporter.add_debug_layer("test_layer", &geometry);

        // Need to set debug with the layer
        let mut debug_layers = DebugLayers {
            layers: HashMap::new(),
        };
        debug_layers.layers.insert("test_layer".to_string(), layer);
        exporter.set_debug(debug_layers);

        let json = exporter.to_json();
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();

        // Check debug section has the layer
        assert!(parsed["debug"]["test_layer"]["lines"].is_object());
        assert!(parsed["debug"]["test_layer"]["points"].is_object());
    }
}
