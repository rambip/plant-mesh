# Project Roadmap: plant-mesh

This document outlines the development path for the `plant-mesh` generator, focusing on the transition from a Bevy-centric prototype to a high-performance, multi-platform geometry engine.

## Phase 1: Core Logic Decoupling (Current Focus)
The goal is to extract the geometric algorithms into a "headless" crate that can be used by both the Bevy renderer and the Python extension.

- **Extract `plant-core` crate**
    - Move `src/meshing/algorithms.rs` (Convex Hull, Catmull-Rom).
    - Move `src/meshing/particles.rs` (HGrid, Particle Simulation).
    - Move `src/growing/` (Skeleton generation logic).
    - **Nested Details:**
        - Replace `bevy::math` with `glam` for zero-overhead math.
        - Implement a custom `NormalSolver` to replace Bevy's built-in mesh normal calculation.
        - Feature-gate `serde` for configuration serialization.
- **Refactor for Pure Functions**
    - Remove all `Query` and `Res` types from the simulation logic.
    - Ensure the simulation takes a `Config` struct and returns a `RawMesh` (Vec of vertices, indices, and normals).

## Phase 2: Python Integration (`maturin`)
Create a high-performance Python package to allow researchers and artists to generate meshes within Python environments (Jupyter, Blender, etc.).

- **Create `plant-python` wrapper**
    - Use `PyO3` to expose Rust structs to Python.
    - Map Python dictionaries to Rust `GrowConfig` and `MeshConfig`.
- **Performance Optimizations**
    - **Rayon Integration:** Parallelize the `HGrid` neighbor search in the particle simulation.
    - **GIL Management:** Release the Python Global Interpreter Lock during the simulation to keep the host application responsive.
- **Data Transfer**
    - Return mesh data as NumPy-compatible buffers to avoid expensive string parsing or list copying.

## Phase 3: Interactive Python Widget
Develop a lightweight visualization tool for rapid parameter tuning.

- **Backend:** The `plant-python` extension.
- **Frontend Options:**
    - **PyVista/IPyWidgets:** For integration into Jupyter Lab.
    - **PySide6 + WebEngine:** To reuse the existing WASM/WebGL build for a native desktop feel.
- **Nested Details:**
    - Implement a "Draft Mode" in the simulation (fewer particles/steps) for real-time slider feedback.
    - Add a "Live Export" button to save `.obj` or `.glb` files directly from the widget.

## Phase 4: Algorithmic Improvements
Enhance the quality and diversity of the generated plants.

- **Advanced Branching Models**
    - Implement biology-based growth (e.g., apical dominance, gravitropism).
    - Support non-circular branch profiles (elliptical or custom polygons).
- **Mesh Quality**
    - Improve the "Branch Fusion" logic to handle 3+ children at a single node.
    - Implement strand-swapping to reduce "mangling" in high-curvature areas.

## Phase 5: Distribution
- **PyPI:** Publish the `plant-mesh` package with pre-compiled binaries for Windows, macOS, and Linux.
- **Documentation:** Create a `docs/` site with interactive examples using the Python widget.
