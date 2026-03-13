# Architecture - plant-mesh

## Overview

`plant-mesh` (tubulin) is a procedural tree mesh generator with three frontends:

- Rust core (`crates/tubulin-core`) for growth, volumetric strands, meshing, and TreeMesh export.
- Bevy demo (`crates/bevy-demo`) for interactive native/wasm visualization.
- Python + JS viewer for notebooks, where Rust emits TreeMesh JSON and a bundled JS renderer displays it.

The main flow is:

`Seed -> PlantNode -> TreeSkeleton -> VolumetricTree -> GeometryData -> TreeMesh JSON`

## Workspace Layout

- `crates/tubulin-core/`: geometry and serialization core.
- `crates/bevy-demo/`: interactive app, camera controls, debug toggles.
- `crates/bevy-gizmos/`, `crates/bevy-simple-graphics/`: rendering support crates used by the demo.
- `python/tubulin/`: Python package containing `_tubulin` extension and bundled `render.js`.
- `js/src/`: TreeMesh decoder and viewer source.
- `GEO_SPEC.md`: TreeMesh IR format and operator semantics.

## Pipeline Stages

### 1) Growth (`PlantNode`)

- `Seed::grow_plant` creates a recursive branching structure.
- Nodes carry `PlantNodeProps` (`position`, `radius`, `orientation`) plus children.
- Config: `GrowConfig`.

### 2) Skeleton (`TreeSkeleton`)

- `PlantNode::grow_skeleton` flattens recursion into array-based data:
  - `node_props: Vec<PlantNodeProps>`
  - `node_info: Vec<NodeInfo>` (parent/children/depth/index)
- This representation is used by downstream algorithms and debug export.

### 3) Strands (`VolumetricTree`)

- `TreeSkeleton::grow_strands` builds particle trajectories with `TrajectoryBuilder`.
- Leaves spawn particle clouds in local branch disks.
- At each parent, child clouds are projected/merged and relaxed via `spread_points`.
- A uniform grid (`HGrid`) accelerates local neighbor queries for repulsion.
- Result:
  - `particles_per_node`: particle IDs active on each skeleton node
  - `trajectories`: per-particle 3D path across tree depth

### 4) Surface Meshing (`GeometryData`)

- `VolumetricTree::build_mesh` samples trajectories with centripetal Catmull-Rom (`extended_catmull_spline`).
- Per-slice contours are computed from projected particle clouds using convex hull.
- Neighboring contours are stitched with triangle strips (`mesh_between_contours`).
- Branch joins handle split/merge regions using local contour partitioning.
- Output: positions, normals, colors, triangles.

### 5) Export (`TreeEncoder`)

- `TreeEncoder` emits TreeMesh JSON program:
  - compressed integer `buffers` (Rice + zigzag for signed deltas)
  - expression `outputs` (`vec3`, `divp2`, `cumsum`, `triangle`, `spline`, ...)
  - optional `debug` layers (`points`, `lines`, and `show` initial visibility flag)
- `Buffer.length` is emitted to avoid Rice padding decode drift.

## Debug Architecture

- Rust types implement `VisualDebug` and produce `DebugGeometry` (`lines`, `circles`, `points`).
- Python bindings expose debug objects and JSON conversion.
- TreeMesh debug layers can be toggled in viewer UI.
- `show` on each debug layer controls initial visible/hidden state.

## Frontends

### Bevy Demo

- Loads `TreeConfig` from TOML.
- Regenerates trees on demand.
- Supports camera orbit/tilt/zoom and debug toggles.

### Python Notebooks

- pyo3 module `_tubulin` exposes pipeline objects (`Seed`, `Skeleton`, `VolumetricTree`, `TreeMesh`, `DebugData`).
- `_repr_html_` embeds viewer HTML and loads bundled JS at runtime from package path (`tubulin.__file__` + `render.js`).

### JS Decoder/Renderer

- `js/src/decoder.js` evaluates TreeMesh expressions into typed arrays.
- `js/src/render.js` uploads geometry, manages camera interaction, and renders debug layers.
- The notebook bundle target is `python/tubulin/render.js`.

## Build and Deployment Notes

- Python package build: `.venv/bin/maturin develop`.
- JS bundle for Python package: `bun build js/src/render.js --outfile=python/tubulin/render.js`.
- GitHub Pages workflow builds the root Trunk app and deploys `dist/` artifact.
