// TODO: rename

pub fn min_by_key<I, T, F>(iter: I, mut f: F) -> Option<T>
where
    I: IntoIterator<Item = T>,
    F: FnMut(&T) -> f32,
{
    iter.into_iter().reduce(|a, b| {
        let a_key = f(&a);
        let b_key = f(&b);
        
        // Handle NaN cases - treat NaN as greater than everything
        match (a_key.is_nan(), b_key.is_nan()) {
            (true, true) => a,   // If both are NaN, keep first
            (true, false) => b,  // If a is NaN, choose b
            (false, true) => a,  // If b is NaN, choose a
            (false, false) => {
                if a_key <= b_key { a } else { b }
            }
        }
    })
}

pub fn max_by_key<I, T, F>(iter: I, mut f: F) -> Option<T>
where
    I: IntoIterator<Item = T>,
    F: FnMut(&T) -> f32,
{
    iter.into_iter().reduce(|a, b| {
        let a_key = f(&a);
        let b_key = f(&b);
        
        // Handle NaN cases - treat NaN as greater than everything
        match (a_key.is_nan(), b_key.is_nan()) {
            (true, true) => a,   // If both are NaN, keep first
            (true, false) => b,  // If a is NaN, choose b
            (false, true) => a,  // If b is NaN, choose a
            (false, false) => {
                if a_key >= b_key { a } else { b }
            }
        }
    })
}
