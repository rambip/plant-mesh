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
