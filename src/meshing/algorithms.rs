use bevy::math::{FloatExt, Vec2, Vec3};

use crate::tools::FloatProducer;

#[derive(Copy, Clone, Debug)]
pub enum SplineIndex {
    /// real number between 0 and 1
    Global(f32),
    /// index of the point pair to consider and real number between 0 and 1
    Local(usize, f32),
}

// r is between 0 and points.len()-1
pub fn extended_catmull_spline(points: &[Vec3], pos: SplineIndex) -> Vec3 {
    let n = points.len();

    // edge case, we might get an index error
    let (i0, r) = match pos {
        SplineIndex::Global(1.) => (n - 2, 1.),
        SplineIndex::Global(t) => {
            let step = f32::floor((n - 1) as f32 * t);
            (step as usize, t * (n - 1) as f32 - step)
        }
        SplineIndex::Local(i, 0.) if i == n - 1 => (n - 2, 1.),
        SplineIndex::Local(i0, r) => (i0, r),
    };

    assert!(i0 <= n - 2, "{pos:?}, n={n}");

    let points_to_interpolate = [
        if i0 == 0 {
            2. * points[0] - points[1]
        } else {
            points[i0 - 1]
        },
        points[i0 + 0],
        points[i0 + 1],
        if i0 == n - 2 {
            2. * points[n - 1] - points[n - 2]
        } else {
            points[i0 + 2]
        },
    ];
    let mut knot_sequence = [0.; 4];
    knot_sequence[0] = 0.;
    for i in 1..4 {
        knot_sequence[i] = knot_sequence[i - 1]
            + (points_to_interpolate[i] - points_to_interpolate[i - 1])
                .length()
                .sqrt();
    }
    let t = knot_sequence[1].lerp(knot_sequence[2], r);

    let ratio = |i: usize, j: usize| (t - knot_sequence[i]) / (knot_sequence[j] - knot_sequence[i]);

    let mut a_points = [Vec3::ZERO; 3];
    for i in 0..3 {
        a_points[i] = points_to_interpolate[i].lerp(points_to_interpolate[i + 1], ratio(i, i + 1));
    }

    let b1 = a_points[0].lerp(a_points[1], ratio(0, 2));
    let b2 = a_points[1].lerp(a_points[2], ratio(1, 3));

    b1.lerp(b2, ratio(1, 2))
}

pub fn angle_to(a: Vec2, b: Vec2) -> f32 {
    if a.length_squared() * b.length_squared() == 0. {
        return 0.;
    }

    // TODO: more efficient
    let angle = f32::acos(0.99 * a.dot(b) / f32::sqrt(a.length_squared() * b.length_squared()));

    if a.perp_dot(b) >= 0. {
        angle
    } else {
        2. * std::f32::consts::PI - angle
    }
}

/// returns the set of points in the convex hull of `points`,
/// once projected on a 2D plane.
/// we suppose that the barycenter is in the convex hull
pub fn convex_hull_graham(
    pivot: Option<Vec2>,
    points: &[Vec2],
    min_angle: Option<f32>,
) -> Vec<usize> {
    let min_angle = min_angle.unwrap_or(std::f32::consts::PI);

    let n = points.len();
    let i_min_y = points.into_iter().map(|point| point.y).arg_min().unwrap();
    let pivot = pivot.unwrap_or((1. / n as f32) * points.iter().sum::<Vec2>());
    let u0 = points[i_min_y] - pivot;

    let compare_angle = |a: &usize, b: &usize| {
        angle_to(u0, points[*a] - pivot)
            .partial_cmp(&angle_to(u0, points[*b] - pivot))
            .unwrap()
    };

    assert!(n >= 3);
    let sorted_points: Vec<usize> = {
        let mut result: Vec<usize> = (0..n).into_iter().collect();
        result.sort_by(compare_angle);
        result
    };

    let mut result: Vec<usize> = sorted_points[..2].into();
    let is_sharp_angle = |result: &Vec<usize>, i| {
        let n_hull = result.len();
        let p_1 = result[n_hull - 2];
        let p_2 = result[n_hull - 1];
        let angle = angle_to(points[p_1] - points[p_2], points[i] - points[p_2]);
        angle < min_angle
    };

    for &i in sorted_points.iter().cycle().skip(2).take(n - 1) {
        while result.len() >= 2 && is_sharp_angle(&result, i) {
            result.pop().unwrap();
        }
        result.push(i);
    }
    result.pop();
    result
}

// c2 is above c1, both contours are clockwise
pub fn mesh_between_contours(
    points: &[Vec3],
    c1: &[usize],
    c2: &[usize],
    close_contour: bool,
) -> Vec<usize> {
    assert!(c1.len() > 0);
    assert!(c2.len() > 0);
    let mut f1 = c1.iter();
    let mut f2 = c2.iter();
    let mut result: Vec<usize> = Vec::new();

    enum CountourId {
        C1,
        C2,
    }

    use CountourId::*;

    let current = |f: &std::slice::Iter<_>| f.as_slice()[0];
    let next = |f: &std::slice::Iter<_>| f.as_slice()[1];
    let current_p = |f: &std::slice::Iter<_>| points[f.as_slice()[0]];
    let next_p = |f: &std::slice::Iter<_>| points[f.as_slice()[1]];

    let mut add_point =
        |f1: &mut std::slice::Iter<_>, f2: &mut std::slice::Iter<_>, c: CountourId| match c {
            CountourId::C1 => {
                result.extend([current(f2), current(f1), next(f1)]);
                let _ = f1.next();
            }
            CountourId::C2 => {
                result.extend([current(f2), current(f1), next(f2)]);
                let _ = f2.next();
            }
        };

    let i1_0 = current(&f1);
    let i2_0 = current(&f2);

    while f1.as_slice().len() > 1 || f2.as_slice().len() > 1 {
        // TODO: faire la jointure
        if f2.as_slice().len() == 1 {
            add_point(&mut f1, &mut f2, C1)
        } else if f1.as_slice().len() == 1 {
            add_point(&mut f1, &mut f2, C2)
        } else {
            let d1 = (current_p(&f2) - next_p(&f1)).length();
            let d2 = (current_p(&f1) - next_p(&f2)).length();

            if d1 < d2 {
                add_point(&mut f1, &mut f2, C1)
            } else {
                add_point(&mut f1, &mut f2, C2)
            }
        }
    }

    if close_contour {
        result.extend([current(&f2), current(&f1), i1_0]);
        result.extend([current(&f2), i1_0, i2_0]);
    }

    result
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::prelude::*;

    #[test]
    fn convex_hull_square() {
        let points: Vec<Vec2> = [(0.0, 0.0), (1.0, 0.0), (0.8, 0.8), (1.0, 1.0), (0.0, 1.0)]
            .into_iter()
            .map(|(a, b)| Vec2::new(a, b))
            .collect();

        let convex = convex_hull_graham(Some(points[0]), &points, None);
        assert_eq!(convex, vec![0, 1, 3, 4]);
    }

    #[test]
    fn convex_hull_polygon() {
        let points: Vec<Vec2> = [
            (2.0, 2.0),
            (1.0, 1.0),
            (0.8, 1.0),
            (-1.0, 1.0),
            (-0.9, 0.0),
            (-1.0, -1.1),
            (1.0, -1.0),
            (2.0, 0.0),
        ]
        .into_iter()
        .map(|(a, b)| Vec2::new(a - 10., b))
        .collect();

        let convex = convex_hull_graham(Some(points[5]), &points, None);
        assert_eq!(convex, vec![0, 3, 5, 6, 7]);
    }

    #[test]
    fn test_random() {
        const N: u64 = 10;
        for i in 0..N {
            let mut rng = StdRng::seed_from_u64(i);
            let points: Vec<Vec2> = (0..10)
                .into_iter()
                .map(|_| Vec2::new(rng.gen_range(0..10) as f32, rng.gen_range(0..10) as f32))
                .collect();

            dbg!(&points);

            let result = convex_hull_graham(None, &points, None);

            let mut turns = (0..result.len() - 2)
                .map(|i| (result[i + 0], result[i + 1], result[i + 2]))
                .map(|(p1, p2, p3)| (points[p2] - points[p1]).perp_dot(points[p3] - points[p1]));

            assert!(turns.all(|x| x >= 0.));
        }
    }
}
