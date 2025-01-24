use bevy::math::Vec3;


fn lerp<T>(a: T, b: T, t: f32) -> T 
where T: std::ops::Mul<f32, Output=T> + std::ops::Add<Output=T>
{
    a*(1.-t) + b*t
}

// r is between 0 and points.len()-1
pub fn extended_catmull_spline(points: &[Vec3], r: f32) -> Vec3 {
    let n = points.len();

    // edge case, we might get an index error
    let i0 = usize::min(r as usize, n-2);

    let points_to_interpolate: [Vec3; 4] = {
        let mut p = [Vec3::ZERO; 4];
        p[0] = if i0 == 0 {2.*points[1] - points[0]} else {points[i0-1]};
        p[1] = points[i0+0];
        p[2] = points[i0+1];
        p[3] = if i0 == n-2 {2.*points[n-2] - points[n-1]} else {points[i0+2]};
        p
    };
    let mut knot_sequence = [0.; 4];
    knot_sequence[0] = 0.;
    for i in 1..4 {
        knot_sequence[i] = knot_sequence[i-1] +
            (points_to_interpolate[i] - points_to_interpolate[i-1]).length()
            .sqrt();
    }
    let t = lerp(knot_sequence[1], knot_sequence[2], r - (i0 as f32));

    let ratio = |i: usize, j: usize| (t - knot_sequence[i]) / (knot_sequence[j] - knot_sequence[i]);

    let mut a_points = [Vec3::ZERO; 3];
    for i in 0..3 {
        a_points[i] = lerp(points_to_interpolate[i], points_to_interpolate[i+1], ratio(i, i+1));
    }

    let b1 = lerp(a_points[0], a_points[1], ratio(0, 2));
    let b2 = lerp(a_points[1], a_points[2], ratio(1, 3));

    lerp(b1, b2, ratio(1, 2))
}

fn cmp(a: f32, b: f32) -> std::cmp::Ordering {
    a.partial_cmp(&b).expect("cannot compare NaN")
}

fn det(a: Vec3, b: Vec3) -> f32 {
    a.x * b.y - a.y * b.x
}

/// returns the set of points in the convex hull of `points`,
/// once projected on a 2D plane.
pub fn convex_hull_graham(points: &[Vec3]) -> Vec<usize> {
    let n = points.len();
    let p0 = (0..n)
        .into_iter()
        .min_by(|&i, &j| cmp(points[i].y, points[j].y)).unwrap();

    // probably broken
    let sorted_points : Vec<usize> = {
        let mut result: Vec<_> = (0..n).into_iter()
            .filter(|p| *p != p0)
            .collect();

        result.sort_by(|i: &usize, j: &usize| 
            cmp(0., det(points[*i] - points[p0], points[*j] - points[p0]))
            );
        result
    };

    let mut result: Vec<usize> = vec![p0, sorted_points[0]];

    for &i in &sorted_points[1..] {
        loop {
            let n_hull = result.len();
            let p_1 = result[n_hull-2];
            let p_2 = result[n_hull-1];
            if det(points[p_2] - points[p_1], points[i] - points[p_1]) >= 0. || result.len() < 3 {
                break
            }
            result.pop().unwrap();
        }
        result.push(i);
    }

    result
}

// c1 is above c2, both contours are clockwise
pub fn mesh_between_contours(points: &[Vec3], c1: &[usize], c2: &[usize]) -> Vec<usize> {
    let mut f1 = c1.iter();
    let mut f2 = c2.iter();
    let mut result: Vec<usize> = Vec::new();

    enum CountourId {
        C1,
        C2
    }

    use CountourId::*;

    let current = |f: &std::slice::Iter<_>| f.as_slice()[0];
    let next = |f: &std::slice::Iter<_>| f.as_slice()[1];
    let current_p = |f: &std::slice::Iter<_>| points[f.as_slice()[0]];
    let next_p = |f: &std::slice::Iter<_>| points[f.as_slice()[1]];

    let mut add_point = |f1: &mut std::slice::Iter<_>, f2: &mut std::slice::Iter<_>, c: CountourId| {
        match c {
            CountourId::C1 => {
                result.extend([current(f1), current(f2), next(f1)]);
                let _ = f1.next();
            },
            CountourId::C2 => {
                result.extend([current(f1), current(f2), next(f2)]);
                let _ = f2.next();
            },
        }
    };

    let i1_0 = current(&f1);
    let i2_0 = current(&f2);

    while f1.as_slice().len() > 1 || f2.as_slice().len() > 1 {
        // TODO: faire la jointure
        if f2.as_slice().len() == 1 {add_point(&mut f1, &mut f2, C1)}
        else if f1.as_slice().len() == 1 {add_point(&mut f1, &mut f2, C2)}

        else {

            //let a1 = cosine_sim(current_p(&f1) - next_p(&f1) , current_p(&f2) - next_p(&f1));
            //let a2 = cosine_sim(current_p(&f1) - next_p(&f2) , current_p(&f2) - next_p(&f2));
            let d1 = f32::max(0. , (current_p(&f2) - next_p(&f1)).length());
            let d2 = f32::max((current_p(&f1) - next_p(&f2)).length() , 0.);

            if d1 < d2 { add_point(&mut f1, &mut f2, C1) }
            else       { add_point(&mut f1, &mut f2, C2) }
        }

    }
    result.extend([current(&f1), current(&f2), i1_0]);
    result.extend([i1_0, current(&f2), i2_0]);

    result
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::prelude::*;

    #[test]
    fn convex_hull_square(){
        let points: Vec<Vec3> = [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.8, 0.8),
            (1.0, 1.0),
            (0.0, 1.0),
        ].into_iter()
            .map(|(a, b)| Vec3::new(a, b, 0.))
            .collect();

        let convex = convex_hull_graham(&points);
        assert_eq!(convex, vec![0, 1, 3, 4]);
    }

    #[test]
    fn convex_hull_polygon(){
        let points: Vec<Vec3> = [
            (2.0, 2.0),
            (1.0, 1.0),
            (0.8, 1.0),
            (-1.0, 1.0),
            (-0.9, 0.0),
            (-1.0, -1.1),
            (1.0, -1.0),
            (2.0, 0.0),
        ].into_iter()
            .map(|(a, b)| Vec3::new(a-10., b, 0.))
            .collect();

        let convex = convex_hull_graham(&points);
        assert_eq!(convex, vec![5, 6, 7, 0, 3]);
    }

    #[test]
    fn test_random() {
        for i in 0..10 {
            let mut rng = StdRng::seed_from_u64(i);
            let points: Vec<Vec3> = (0..10)
                .into_iter()
                .map(|_| Vec3::new(rng.gen_range(0..10) as f32, rng.gen_range(0..10) as f32, rng.gen_range(0..10) as f32))
                .collect();

            dbg!(&points);

            let result = convex_hull_graham(&points);

            let mut turns = (0..result.len()-2)
                .map(|i| (result[i+0], result[i+1], result[i+2]))
                .map(|(p1, p2, p3)| det(points[p2] - points[p1], points[p3] - points[p1]));

            assert!(turns.all(|x| x >= 0.));
        }
    }
}
