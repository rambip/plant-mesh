// TODO: rename module

pub trait FloatProducer: ExactSizeIterator<Item=f32> + Sized {
    fn arg_min(self) -> Option<usize> 
    {
        let (mut i_min, mut v_min) = (0, f32::NAN);
        for (i, v) in self.enumerate() {
            if v < v_min || v_min.is_nan() {
                (i_min, v_min) = (i, v)
            }
        }
        if v_min.is_nan() {
            None
        }
        else {
            Some(i_min)
        }
    }

    fn arg_max(self) -> Option<usize> 
    {
        let (mut i_max, mut v_max) = (0, f32::NAN);
        for (i, v) in self.enumerate() {
            if v > v_max || v_max.is_nan() {
                (i_max, v_max) = (i, v)
            }
        }
        if v_max.is_nan() {
            None
        }
        else {
            Some(i_max)
        }
    }
}

impl<I> FloatProducer for I 
where I: ExactSizeIterator<Item=f32> { }

