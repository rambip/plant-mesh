# Agent Guidelines - plant-mesh

## Project Scope & Goal
This project is a procedural plant mesh generator developed for a computer geometry course. It uses a pipeline-based approach to grow tree skeletons and generate volumetric meshes using Bevy for rendering and geometric algorithms like convex hulls and particle-based meshing. See [ARCHITECTURE.md](ARCHITECTURE.md) for a detailed breakdown of the pipeline and internal crates.

## Build & Test
- **Build:** `cargo build`
- **Lint:** `cargo clippy`
- **Test All:** `cargo test`
- **Single Test:** `cargo test -- <test_name_substring>`
- **Run Example:** `cargo run --example <example_name>` (e.g., `convex_hull`)
- **Compile Report:** `typst compile docs/report.typ docs/report.pdf` (requires `typst`)

## Planning & Documentation
- **Roadmap:** Read and update [ROADMAP.md](ROADMAP.md) regularly to track progress and align with the long-term vision.
- **Concision:** Prefer concise, actionable updates over verbose descriptions. Keep documentation and comments high-density.

## Code Style
- **Formatting:** Standard `rustfmt`. Use `cargo fmt` before committing.
- **Naming:** `snake_case` for functions/variables, `PascalCase` for types/traits.
- **Imports:** Grouped by: `std`, external crates (e.g., `bevy`, `rand`), then local modules (`super`, `self`).
- **Types:** Use `f32` for geometric calculations. Prefer `bevy::math` types (`Vec3`, `Quat`).
- **Error Handling:** Use `panic!` for unreachable states in generation logic; otherwise, prefer `Option` or `Result`.
- **Conventions:**
    - Implement `VisualDebug` for types requiring gizmo rendering.
    - Follow the `TreePipelinePhase` pattern for generation steps.
    - Keep shaders in `assets/` and use `.wgsl`.
    - Use `serde` for configuration structs (`GrowConfig`, `MeshConfig`).
