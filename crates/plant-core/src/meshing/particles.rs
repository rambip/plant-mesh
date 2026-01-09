use glam::{Vec2, Vec3};
use rand::{prelude::Distribution, Rng};
use smallvec::SmallVec;
use serde::{Serialize, Deserialize};
#[cfg(feature = "bevy")]
use bevy_gizmos::prelude::Gizmos;
#[cfg(feature = "bevy")]
use bevy_color::Color;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StrandsConfig {
    pub repulsion: f32,
    pub interaction_radius: f32,
    pub n_steps: usize,
    pub dt: f32,
    pub particles_per_leaf: usize,
}

#[derive(Debug)]
struct HGrid {
    corner_a: Vec2,
    corner_b: Vec2,
    cells: Vec<SmallVec<[usize; 32]>>,
    n_x: usize,
    n_y: usize,
}

impl std::ops::IndexMut<(usize, usize)> for HGrid {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.cells[index.1 * self.n_x + index.0]
    }
}

impl std::ops::Index<(usize, usize)> for HGrid {
    type Output = SmallVec<[usize; 32]>;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.cells[index.1 * self.n_x + index.0]
    }
}

impl HGrid {
    fn new_in_region(corner_a: Vec2, corner_b: Vec2, n_x: usize, n_y: usize) -> Self {
        let cells = vec![SmallVec::new(); n_x * n_y];
        Self {
            corner_a,
            corner_b,
            cells,
            n_x,
            n_y,
        }
    }
    fn clear(&mut self) {
        for i_x in 0..self.n_x {
            for i_y in 0..self.n_y {
                self[(i_x, i_y)].clear()
            }
        }
    }
    fn width(&self) -> f32 {
        (self.corner_b - self.corner_a).x
    }
    fn height(&self) -> f32 {
        (self.corner_b - self.corner_a).y
    }
    fn point_to_cell(&self, point: Vec2) -> (usize, usize) {
        let x = f32::min(
            (point.x - self.corner_a.x) / self.width() * self.n_x as f32,
            self.n_x as f32 - 1.,
        );
        let y = f32::min(
            (point.y - self.corner_a.y) / self.height() * self.n_y as f32,
            self.n_y as f32 - 1.,
        );
        (x as usize, y as usize)
    }
    fn insert(&mut self, point: Vec2, id: usize) {
        let i_grid = self.point_to_cell(point);
        self[i_grid].push(id)
    }
    fn query_rect(
        &self,
        mut corner_a: Vec2,
        mut corner_b: Vec2,
        condition: impl Fn(usize) -> bool,
    ) -> SmallVec<[usize; 32]> {
        let mut result = SmallVec::new();
        corner_a.x = f32::max(corner_a.x, self.corner_a.x);
        corner_a.y = f32::max(corner_a.y, self.corner_a.y);
        corner_b.x = f32::min(corner_b.x, self.corner_b.x);
        corner_b.y = f32::min(corner_b.y, self.corner_b.y);
        let cell_a = self.point_to_cell(corner_a);
        let cell_b = self.point_to_cell(corner_b);
        for i_x in cell_a.0..=cell_b.0 {
            for i_y in cell_a.1..=cell_b.1 {
                for &i in &self[(i_x, i_y)] {
                    if condition(i) {
                        result.push(i)
                    }
                }
            }
        }
        result
    }
}

struct CubicKernel {
    h: f32,
}

impl CubicKernel {
    fn f(&self, q: f32) -> f32 {
        if q < 1e0 {
            1e0 - q.powi(2) * 1.5 + q.powi(3) * 0.75
        } else if q < 2e0 {
            0.25 * (2e0 - q).powi(3)
        } else {
            0.
        }
    }

    fn c(&self) -> f32 {
        10e0 / (7e0 * std::f32::consts::PI * self.h.powi(2))
    }

    fn w(&self, l: f32) -> f32 {
        self.c() * self.f(2. * l / self.h)
    }
}

pub struct UniformDisk {
    center: Vec2,
    radius: f32,
}

impl UniformDisk {
    pub fn new(center: Vec2, radius: f32) -> Self {
        Self { center, radius }
    }
}

impl Distribution<Vec2> for UniformDisk {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec2 {
        let angle = rng.gen_range(0f32..std::f32::consts::TAU);
        self.center
            + Vec2::new(angle.cos(), angle.sin())
                * self.radius
                * rng.gen_range(0f32..1f32).sqrt()
    }
}

fn average<I: IntoIterator<Item = Vec2>>(i: I) -> Option<Vec2> {
    let mut result = Vec2::ZERO;
    let mut n = 0.;
    for x in i.into_iter() {
        result += x;
        n += 1.
    }
    if n == 0. {
        None
    } else {
        Some(result / n)
    }
}

pub fn spread_points(points: &mut Vec<Vec2>, radius: f32, config: &StrandsConfig) {
    let n = points.len();
    let mut pressures = vec![Vec2::ZERO; n];
    let repulsion = config.repulsion * radius / config.interaction_radius / (n as f32).sqrt();

    let d = radius * config.interaction_radius;
    let n_cells = (1. / config.interaction_radius) as usize + 1;
    let diag = Vec2::new(d, d);

    let mut grid = HGrid::new_in_region(
        Vec2::new(-radius, -radius),
        Vec2::new(radius, radius),
        n_cells,
        n_cells,
    );
    let mut neighbourghs = vec![SmallVec::<[usize; 32]>::new(); n];

    let kernel = CubicKernel { h: d };

    let max_radius = points.iter().map(|x| x.length()).reduce(f32::max).unwrap();
    for i in 0..n {
        points[i] *= 0.9 * radius / max_radius
    }

    for _ in 0..config.n_steps {
        grid.clear();
        for i in 0..n {
            grid.insert(points[i], i)
        }
        for i in 0..n {
            neighbourghs[i] = grid.query_rect(points[i] - diag, points[i] + diag, |j| {
                (points[i] - points[j]).length() < d
            });
        }
        for i in 0..n {
            pressures[i] = repulsion
                * average(
                    neighbourghs[i]
                        .iter()
                        .filter(|&&j| points[i] != points[j])
                        .map(|&j| (points[i] - points[j]).normalize()),
                )
                .unwrap_or(Vec2::ZERO);
        }
        for i in 0..n {
            let mut velocity = Vec2::ZERO;
            for &j in &neighbourghs[i] {
                let l = (points[i] - points[j]).length();
                velocity += kernel.w(l) * pressures[j];
            }
            points[i] += config.dt * velocity;
            if points[i].length() > radius {
                points[i] = 0.99 * radius * points[i].normalize()
            }
        }
    }
    assert!(points.iter().all(|x| !x.is_nan()));
}

pub struct TrajectoryBuilder {
    pub particles_per_node: Vec<Vec<usize>>,
    pub trajectories: Vec<Vec<Vec3>>,
    pub rng: rand::rngs::StdRng,
}

impl TrajectoryBuilder {
    pub fn new(rng: rand::rngs::StdRng) -> Self {
        Self {
            particles_per_node: Vec::new(),
            trajectories: Vec::new(),
            rng,
        }
    }

    pub fn clear_for_tree(&mut self, tree: &super::super::growing::TreeSkeleton) {
        let n = tree.node_count();
        self.particles_per_node = vec![vec![]; n];
        self.trajectories = vec![];
    }

    pub fn register_particles_for_leaf(&mut self, tree: &super::super::growing::TreeSkeleton, node: usize, cloud: &[Vec2]) {
        for &point in cloud {
            let particle_id = self.trajectories.len();
            self.particles_per_node[node].push(particle_id);
            let mut empty_trajectory = vec![Vec3::ZERO; tree.depth(node)];
            empty_trajectory.push(tree.plane_to_space(node, point));
            self.trajectories.push(empty_trajectory);
        }
    }

    pub fn register_particles_for_node(
        &mut self,
        tree: &super::super::growing::TreeSkeleton,
        node: usize,
        cloud: &[Vec2],
        particle_ids: &[usize],
    ) {
        assert_eq!(cloud.len(), particle_ids.len());
        let n = particle_ids.len();
        for i in 0..n {
            let d = tree.depth(node);
            self.trajectories[particle_ids[i]][d] = tree.plane_to_space(node, cloud[i]);
            self.particles_per_node[node].push(particle_ids[i]);
        }
    }

    pub fn project_particles(&mut self, tree: &super::super::growing::TreeSkeleton, child: usize, offset: Vec3) -> Vec<Vec2> {
        let mut result = Vec::new();
        for &particle_id in &self.particles_per_node[child] {
            let pos_particle = self.trajectories[particle_id][tree.depth(child)];
            result.push(tree.space_to_plane(child, offset + pos_particle));
        }
        result
    }

    pub fn compute_trajectories(
        &mut self,
        tree: &super::super::growing::TreeSkeleton,
        root: usize,
        config: &StrandsConfig,
    ) {
        assert!(self.particles_per_node[root].is_empty());
        let radius = tree.radius(root);

        let children = tree.children(root);
        match children.len() {
            0 => {
                let mut cloud = Vec::new();
                let disk = UniformDisk::new(Vec2::ZERO, radius);
                let mut rng = self.rng.clone();
                for _ in 0..config.particles_per_leaf {
                    cloud.push(disk.sample(&mut rng));
                }
                spread_points(&mut cloud, radius, config);
                self.register_particles_for_leaf(tree, root, &cloud);
            }
            1 => {
                let child = children[0];
                self.compute_trajectories(tree, child, config);
                let mut cloud = self.project_particles(tree, child, Vec3::ZERO);
                spread_points(&mut cloud, radius, config);
                let particle_ids = self.particles_per_node[child].clone();
                self.register_particles_for_node(tree, root, &cloud, &particle_ids);
            }
            2 => {
                let m_child = children[0];
                let s_child = children[1];
                self.compute_trajectories(tree, m_child, config);
                self.compute_trajectories(tree, s_child, config);
                let normal = tree.normal(root);
                let offset_direction = {
                    let u = tree.position(m_child) - tree.position(s_child);
                    (u - u.dot(normal) * normal).normalize()
                };
                let m_cloud = self.project_particles(tree, m_child, offset_direction * tree.radius(m_child));
                let s_cloud = self.project_particles(tree, s_child, -offset_direction * tree.radius(s_child));

                let mut cloud = m_cloud;
                cloud.extend(s_cloud);

                let mut particle_ids = self.particles_per_node[m_child].clone();
                particle_ids.extend(&self.particles_per_node[s_child]);

                spread_points(&mut cloud, radius, config);
                self.register_particles_for_node(tree, root, &cloud, &particle_ids);
            }
            _ => panic!("did not expect more than 2 childs for node {root}"),
        }
    }
}

impl crate::VisualDebug for TrajectoryBuilder {
    type Flags = bool;
    #[cfg(feature = "bevy")]
    fn debug(&self, gizmos: &mut Gizmos, debug_flags: bool) {
        if debug_flags {
            let mut rng = self.rng.clone();
            for (i_t, traj) in self.trajectories.iter().enumerate() {
                let a: f32 = i_t as f32 / self.trajectories.len() as f32;
                let b: f32 = rng.gen();
                let color = Color::srgb(1., 0.3 + 0.5 * a, 0.3 + 0.5 * b);
                let positions = (0..50).map(|i| {
                    super::algorithms::extended_catmull_spline(traj, super::algorithms::SplineIndex::Global((i as f32 / 49.).powi(2)))
                });
                gizmos.linestrip(positions, color);
            }
        }
    }
}
