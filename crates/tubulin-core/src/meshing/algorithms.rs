use std::collections::VecDeque;
use glam::{Vec2, Vec3};

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
    let t = knot_sequence[1] + (knot_sequence[2] - knot_sequence[1]) * r;

    let ratio = |i: usize, j: usize| (t - knot_sequence[i]) / (knot_sequence[j] - knot_sequence[i]);

    let mut a_points = [Vec3::ZERO; 3];
    for i in 0..3 {
        a_points[i] = points_to_interpolate[i].lerp(points_to_interpolate[i + 1], ratio(i, i + 1));
    }

    let b1 = a_points[0].lerp(a_points[1], ratio(0, 2));
    let b2 = a_points[1].lerp(a_points[2], ratio(1, 3));

    b1.lerp(b2, ratio(1, 2))
}

fn angle_to(a: Vec2, b: Vec2) -> f32 {
    let mut result = a_angle_to(a, b);
    if result < 0. {
        result += std::f32::consts::TAU
    }
    result
}

fn a_angle_to(a: Vec2, b: Vec2) -> f32 {
    f32::atan2(a.x * b.y - a.y * b.x, a.x * b.x + a.y * b.y)
}

/// returns the set of points in the convex hull of `points`,
/// once projected on a 2D plane.
/// we suppose that the barycenter is in the convex hull
pub fn convex_hull_graham(points: &[Vec2], min_angle: Option<f32>) -> impl ExactSizeIterator<Item=usize> {
    let min_angle = min_angle.unwrap_or(std::f32::consts::PI);

    assert!(!points.is_empty());
    let n = points.len();
    let pivot = points.iter().sum::<Vec2>() / n as f32;

    let u = Vec2::X;

    let compare_angle = |&a: &usize, &b: &usize| {
        let angle_a = angle_to(u, points[a] - pivot);
        let angle_b = angle_to(u, points[b] - pivot);
        angle_a.partial_cmp(&angle_b).unwrap()
    };

    assert!(n >= 3, "not enough particles to compute convex hull");
    let mut sorted_points = {
        let mut result: Vec<usize> = (0..n).into_iter().collect();
        result.sort_by(compare_angle);
        result.into_iter().cycle()
    };

    let is_sharp_angle = |result: &VecDeque<usize>, i| {
        let n_hull = result.len();
        let p_1 = result[n_hull - 2];
        let p_2 = result[n_hull - 1];
        let angle = angle_to(points[p_1] - points[p_2], points[i] - points[p_2]);
        if angle > 0.99 * std::f32::consts::TAU {
            // very strange: the convex hull is very pointy.
            // it is probably a numerical instability.
            return true;
        }
        angle < min_angle
    };

    let mut result = VecDeque::new();
    result.push_back(sorted_points.next().unwrap());
    // FIXME: stop when the 2 last points are the 2 first points
    for i in sorted_points {
        while result.len() >= 2 && is_sharp_angle(&result, i) {
            result.pop_back().unwrap();
        }
        if result.len() >= 3 && result.get(1) == Some(&i) {
            let previous = result.pop_front().unwrap();
            if previous == result[result.len()-1] {
                break
            }
        }
        result.push_back(i);
    }
    // we remove the last `i0` we added to the list,
    // because it was already added at the begining.
    result.into_iter()
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
    let mut f1_idx = 0;
    let mut f2_idx = 0;
    let mut result: Vec<usize> = Vec::new();

    let current_p = |c: &[usize], idx: usize| points[c[idx]];
    let next_p = |c: &[usize], idx: usize| points[c[idx + 1]];

    let i1_0 = c1[0];
    let i2_0 = c2[0];

    while f1_idx < c1.len() - 1 || f2_idx < c2.len() - 1 {
        if f2_idx == c2.len() - 1 {
            result.extend([c2[f2_idx], c1[f1_idx], c1[f1_idx + 1]]);
            f1_idx += 1;
        } else if f1_idx == c1.len() - 1 {
            result.extend([c2[f2_idx], c1[f1_idx], c2[f2_idx + 1]]);
            f2_idx += 1;
        } else {
            let d1 = (current_p(c2, f2_idx) - next_p(c1, f1_idx)).length();
            let d2 = (current_p(c1, f1_idx) - next_p(c2, f2_idx)).length();

            if d1 < d2 {
                result.extend([c2[f2_idx], c1[f1_idx], c1[f1_idx + 1]]);
                f1_idx += 1;
            } else {
                result.extend([c2[f2_idx], c1[f1_idx], c2[f2_idx + 1]]);
                f2_idx += 1;
            }
        }
    }

    if close_contour {
        result.extend([c2[f2_idx], c1[f1_idx], i1_0]);
        result.extend([c2[f2_idx], i1_0, i2_0]);
    }

    result
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn convex_hull_square() {
        let points: Vec<Vec2> = [(0.0, 0.0), (1.0, 0.0), (0.8, 0.8), (1.0, 1.0), (0.0, 1.0)]
            .into_iter()
            .map(|(a, b)| Vec2::new(a, b))
            .collect();

        let mut convex: Vec<usize> = convex_hull_graham(&points, None).collect();
        convex.sort();
        assert_eq!(convex, vec![0, 1, 3, 4]);
    }
}
