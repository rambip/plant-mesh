# Project Roadmap: plant-mesh

This document outlines the development path for the `plant-mesh` generator, focusing on the transition from a Bevy-centric prototype to a high-performance, multi-platform geometry engine.

## Phase 1: Core Logic Decoupling
Extract geometric algorithms into a "headless" `plant-core` crate to enable multi-platform usage.

- **Crate Extraction**
    - Move `src/meshing/algorithms.rs` (Convex Hull, Catmull-Rom) and `src/meshing/particles.rs` (HGrid, Simulation).
    - Move `src/growing/` (Skeleton generation logic).
    - Replace `bevy::math` with `glam` for zero-overhead math (sticking to `f32`).
- **Mesh Abstraction via `mesh-tools`**
    - Adopt `mesh-tools` as the internal mesh representation to maintain engine-agnosticism.
    - Use `mesh-tools::compat` for seamless, efficient conversion to `bevy_render::Mesh`.
    - Leverage `mesh-tools` built-in processing for normal calculation and attribute management.
- **Feature-Gating & Traits**
    - Feature-gate `VisualDebug` and `Gizmos` behind a `bevy` flag.
    - Ensure `RawMesh` (via `mesh-tools`) supports vertex colors (`Vec<[f32; 4]>`).
- **Remaining Questions & Clarifications**
    - Does the current Bevy version (0.14/0.15) align with `mesh-tools` dependency requirements?

## Phase 2: Python Integration (`maturin`)
Create a native Python extension to expose the generation pipeline.

- **Core Python Binding (The "Bridge")**
    - Use `PyO3` to expose `generate_mesh(config: Dict) -> (Vertices, Indices, Normals, Colors)`.
    - Release the Global Interpreter Lock (GIL) using `Python::allow_threads` during simulation.
- **Progress Reporting**
    - Implement a callback mechanism in the Rust simulation loop to allow Python-side progress bars (e.g., `tqdm`).
- **Numpy Extension (Performance Layer)**
    - Implement an optional `numpy` feature in Rust using `rust-numpy` for zero-copy `PyArray` returns.
- **Remaining Questions & Clarifications**
    - How to efficiently map `mesh-tools` attribute buffers to NumPy arrays via PyO3?

## Phase 3: Interactive Widget
Develop a visualization tool for rapid parameter tuning.

- **Existing Ecosystem Integration**
    - Primary target: **PyVista** or `ipyvolume` for robust 3D support in Jupyter/Notebooks.
- **Custom Extension**
    - Future: Build a specialized Three.js widget using the **anywidget** framework for a tailored UI.
- **Draft Mode Logic**
    - Implement a "Low-Fidelity" simulation mode (fewer particles, larger `dt`) for real-time slider feedback.
- **Remaining Questions & Clarifications**
    - How to best synchronize the "Skeleton" view vs "Mesh" view in third-party widgets?

## Phase 4: Performance & Scaling
Optimize the simulation to handle complex forests.

- **Rayon Parallelism**
    - Parallelize the `HGrid` neighbor search and force calculation (CPU-bound).
- **Memory Management**
    - Optimize `SmallVec` sizes in the `HGrid` to minimize heap fragmentation.
- **Remaining Questions & Clarifications**
    - Is there a threshold where a GPU-based compute shader (via `wgpu`) becomes necessary?
