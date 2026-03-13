pub trait FloatProducer: ExactSizeIterator<Item = f32> + Sized {
    fn arg_min(self) -> Option<usize> {
        let (mut i_min, mut v_min) = (0, f32::NAN);
        for (i, v) in self.enumerate() {
            if v < v_min || v_min.is_nan() {
                (i_min, v_min) = (i, v)
            }
        }
        if v_min.is_nan() {
            None
        } else {
            Some(i_min)
        }
    }

    fn arg_max(self) -> Option<usize> {
        let (mut i_max, mut v_max) = (0, f32::NAN);
        for (i, v) in self.enumerate() {
            if v > v_max || v_max.is_nan() {
                (i_max, v_max) = (i, v)
            }
        }
        if v_max.is_nan() {
            None
        } else {
            Some(i_max)
        }
    }
}

impl<I> FloatProducer for I where I: ExactSizeIterator<Item = f32> {}

pub fn split_slice_circular(source: &[usize], start: i32, end: i32) -> (Vec<usize>, Vec<usize>) {
    let mut direct = Vec::new();
    let mut indirect = Vec::new();
    let n = source.len() as i32;

    let id = |x| ((x + n) % n) as usize;

    let mut i = start;
    while id(i) != id(end) {
        direct.push(source[id(i)]);
        i += 1;
    }
    direct.push(source[id(end)]);
    while id(i) != id(start) {
        indirect.push(source[id(i)]);
        i += 1;
    }
    indirect.push(source[id(start)]);
    (direct, indirect)
}

/// Quantizes Vec3 components to i32 arrays with optional component names.
///
/// Extracts x, y, z components from a slice of Vec3s and quantizes them by
/// multiplying by the scale factor and truncating to i32.
///
/// # Arguments
/// * `vectors` - Slice of Vec3 to quantize
/// * `scale` - Scale factor (e.g., 256 for 8-bit fixed point)
///
/// # Returns
/// Tuple of (x_quantized, y_quantized, z_quantized)
pub fn quantize_vec3_components(
    vectors: &[glam::Vec3],
    scale: i32,
) -> (Vec<i32>, Vec<i32>, Vec<i32>) {
    let mut x = Vec::with_capacity(vectors.len());
    let mut y = Vec::with_capacity(vectors.len());
    let mut z = Vec::with_capacity(vectors.len());

    for v in vectors {
        x.push((v.x * scale as f32) as i32);
        y.push((v.y * scale as f32) as i32);
        z.push((v.z * scale as f32) as i32);
    }

    (x, y, z)
}

/// Quantizes color components to i32 arrays.
///
/// Extracts r, g, b, a components from a slice of [f32; 4] colors and
/// quantizes them by multiplying by the scale factor and truncating to i32.
///
/// # Arguments
/// * `colors` - Slice of [f32; 4] colors to quantize
/// * `scale` - Scale factor (e.g., 256 for 8-bit fixed point)
///
/// # Returns
/// Tuple of (r_quantized, g_quantized, b_quantized, a_quantized)
pub fn quantize_color_components(
    colors: &[[f32; 4]],
    scale: i32,
) -> (Vec<i32>, Vec<i32>, Vec<i32>, Vec<i32>) {
    let mut r = Vec::with_capacity(colors.len());
    let mut g = Vec::with_capacity(colors.len());
    let mut b = Vec::with_capacity(colors.len());
    let mut a = Vec::with_capacity(colors.len());

    for c in colors {
        r.push((c[0] * scale as f32) as i32);
        g.push((c[1] * scale as f32) as i32);
        b.push((c[2] * scale as f32) as i32);
        a.push((c[3] * scale as f32) as i32);
    }

    (r, g, b, a)
}
