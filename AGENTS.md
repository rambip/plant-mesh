# tubulin

Procedural plant mesh generator using volumetric invigoration. A WIP with a Rust core, Python bindings, and a JS viewer.

## Repository Structure

```
tubulin/
├── crates/
│   ├── tubulin-core/        # Pure geometry library: growing pipeline, meshing, splines
│   │                        # Feature flags: "bevy" (Bevy types), "python" (pyo3 bindings)
│   ├── bevy-demo/           # Bevy interactive viewer (binary). Depends on tubulin-core.
│   ├── bevy-gizmos/         # Vendored Bevy gizmos crate (upstream fork)
│   └── bevy-simple-graphics/ # Minimal Bevy render pipeline helper
├── python/
│   └── tubulin/             # Python package source
│       └── __init__.py      # Will expose TreeMesh class with _repr_html_()
├── js/                      # JS decoder + Three.js viewer
│   ├── src/
│   │   ├── decoder.js       # TreeMesh format decoder (Rice, spline eval, operators)
│   │   ├── generate.js      # Geometry generator (dev/test, outputs geometry.json)
│   │   └── render.js        # Three.js renderer
│   └── dist/                # Built JS bundle (gitignored, regenerate with bun)
├── Cargo.toml               # Workspace manifest only (no [package])
├── pyproject.toml           # Maturin build config — points to crates/tubulin-core
└── GEO_SPEC.md              # TreeMesh intermediate representation spec
```

## JS Decoder (WIP)

The `js/` folder contains a JavaScript decoder that consumes mesh data from the Rust core and renders it in the browser using Three.js.

> **IMPORTANT for agents:** Read [js/DEVELOP.md](js/DEVELOP.md) before making any changes to the JS decoder, renderer, or geometry generation. It documents the encoding format, common failure modes, and their root causes — skipping it will likely result in subtle bugs that are hard to diagnose.

See [js/DEVELOP.md](js/DEVELOP.md) for stack, commands, and format details.

# Agent Guidelines

## Project Scope & Goal
Procedural plant mesh generator for a computer geometry course. Pipeline: grow tree skeleton → volumetric meshing → export as TreeMesh JSON → decode in JS or Python for rendering. See [ARCHITECTURE.md](ARCHITECTURE.md) for the full pipeline breakdown.

## Build & Test
- **Build all:** `cargo build`
- **Build demo:** `cargo build -p bevy-demo`
- **Lint:** `cargo clippy`
- **Test:** `cargo test`
- **Single test:** `cargo test -- <test_name_substring>`
- **Run demo:** `cargo run -p bevy-demo`
- **Run example:** `cargo run --example <example_name>`
- **Compile report:** `typst compile docs/report.typ docs/report.pdf` (requires `typst`)
- **Python build:** `maturin develop` (from repo root, requires Python venv)

## Planning & Documentation
- **Roadmap:** Read and update [ROADMAP.md](ROADMAP.md) regularly.
- **Concision:** Prefer concise, actionable updates. Keep documentation high-density.

## Code Style
- **Formatting:** `rustfmt`. Run `cargo fmt` before committing.
- **Naming:** `snake_case` for functions/variables, `PascalCase` for types/traits, `kebab-case` for crate names.
- **Imports:** Grouped: `std`, external crates, then local (`super`, `self`).
- **Types:** `f32` for geometry. Prefer `bevy::math` types (`Vec3`, `Quat`) where available.
- **Error handling:** `panic!` for unreachable states in generation; otherwise `Option`/`Result`.
- **Conventions:**
    - Implement `VisualDebug` for types requiring gizmo rendering.
    - Follow the `TreePipelinePhase` pattern for generation steps.
    - Keep shaders in `assets/` as `.wgsl`.
    - Use `serde` for config structs (`GrowConfig`, `MeshConfig`).
