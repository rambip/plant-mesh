# Tubulin Python API — Roadmap

## Goal

Expose the plant mesh pipeline as an inspectable, composable Python API for use in
educational notebooks. Each pipeline stage is a first-class object: configurable,
runnable in isolation, and able to emit debug geometry that renders interactively in
the browser.

---

## Pipeline Overview

```
Grow  →  Strands  →  Mesh
 ↓           ↓         ↓
skeleton  strand     surface
           cloud      mesh
```

## Export Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        Python Builder                           │
│  Tubulin().grow().strands().mesh().debug().run()               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼ calls Rust pipeline
┌─────────────────────────────────────────────────────────────────┐
│                     Rust Pipeline                               │
│  PlantNode → TreeSkeleton → TrajectoryBuilder → GeometryData   │
│                                           │                    │
│                    VisualDebug::debug_data()                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    TreeMeshExporter                             │
│  - add_buffer() for encoded data                               │
│  - add_debug_layer(DebugGeometry) → JSON expressions           │
│  - to_json() → TreeMesh JSON string                            │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼ passes JSON
┌─────────────────────────────────────────────────────────────────┐
│                   Python Wrapper                                │
│  GeometryData.to_json() → str                                  │
│  _repr_html_() renders JS viewer with JSON                     │
└─────────────────────────────────────────────────────────────────┘
```

Each stage takes the output of the previous one. All stages are lazy: nothing runs
until `.run()` is called on the pipeline.

---

## API Design

### Builder pattern

```python
from tubulin import Tubulin

tree = (
    Tubulin()
    .grow(base_radius=1.0, up_attraction=0.15, branch_variance=0.5)
    .strands(repulsion=0.5, n_steps=10)
    .mesh(spacing=0.5, leaf_size=0.5)
    .run()
)

tree  # renders Three.js viewer inline in the notebook
```

### Stage-level inspection

Each stage can be run independently, returning a typed intermediate result:

```python
skeleton    = Tubulin().grow(base_radius=1.0).run()
strand_cloud = Tubulin().grow(...).strands(...).run()
```

Intermediate types (`Skeleton`, `StrandCloud`, `TreeMesh`) are real Python objects
backed by Rust via pyo3. They expose data lazily as numpy arrays:

```python
skeleton.points        # np.ndarray, shape (N, 3) — node positions
skeleton.radii         # np.ndarray, shape (N,)
strand_cloud.points    # np.ndarray, shape (M, 3)
mesh.vertices          # np.ndarray, shape (V, 3)
mesh.triangles         # np.ndarray, shape (T, 3), dtype int
```

Each intermediate type has a `_repr_html_()` that renders it in the notebook with
the Three.js viewer.

### Debug layers

Call `.debug()` on any stage to include its debug geometry in the final output:

```python
tree = (
    Tubulin()
    .grow(base_radius=1.0).debug()      # emits skeleton layer
    .strands(repulsion=0.5).debug()     # emits strand cloud layer
    .mesh(spacing=0.5)
    .run()
)
```

The viewer shows a checkbox per stage: **Skeleton**, **Strands**, **Mesh**. Layers
not opted into are absent from the payload entirely.

---

## TreeMesh JSON Format — Debug Extension

The existing format gains an optional `debug` top-level key. Each stage that ran
with `.debug()` contributes a named sub-object. All arrays use the existing
buffer/expr encoding (Rice-encoded, same as `vertices`/`triangles`).

```json
{
  "treemesh": "0.1",
  "spline_convention": "reflect",
  "buffers": { "...": "..." },
  "outputs": {
    "vertices":  "<expr>",
    "triangles": "<expr>"
  },
  "debug": {
    "skeleton": {
      "lines":   { "starts": "<expr>", "ends": "<expr>", "colors": "<expr>" },
      "circles": { "origins": "<expr>", "orientations": "<expr>", "radii": "<expr>", "colors": "<expr>" }
    },
    "strands": {
      "points": { "positions": "<expr>", "colors": "<expr>" }
    },
    "mesh": {
      "lines": { "starts": "<expr>", "ends": "<expr>", "colors": "<expr>" }
    }
  }
}
```

All `debug` fields are optional. The JS viewer iterates the keys present and creates
one checkbox per key. No changes needed to the buffer encoding or Rice decoder.

---

## Rust Side

### `DebugGeometry` — pure data, no Bevy dependency

```rust
pub struct DebugFrame { pub origin: Vec3, pub orientation: Quat }

pub struct DebugGeometry {
    pub lines:   Vec<(Vec3, Vec3, [f32; 4])>,
    pub circles: Vec<(DebugFrame, f32, [f32; 4])>,
    pub points:  Vec<(Vec3, [f32; 4])>,
}

pub trait VisualDebug {
    fn debug_data(&self) -> DebugGeometry;
}
```

The pipeline collects debug output by calling `debug_data()` on each stage,
producing a `DebugGeometry` which is then passed to `add_debug_layer()` with
a layer name string. The existing Bevy gizmo rendering becomes a thin adapter.

### Intermediate types exposed via pyo3

`Skeleton`, `StrandCloud`, and `TreeMesh` are wrapped with `#[pyclass]`. Numpy
arrays are returned via `rust-numpy` (`PyArray`), computed on first access and
cached. No Bevy types cross the Python boundary.

---

## JS Viewer

- Reads `debug` key from the JSON payload if present
- Creates one `THREE.Group` per layer, toggled by checkbox
- Lines → `THREE.LineSegments`, circles → `THREE.RingGeometry` oriented by
  quaternion, points → `THREE.Points`
- Checkbox state is local (no round-trip to Python needed)

---

## Milestones

### M1 — `DebugGeometry` infrastructure ✓
- [x] Define `DebugGeometry` and `VisualDebug` trait in `tubulin-core` (no Bevy feature flag)
- [x] Implement `VisualDebug` for `TreeSkeleton`
- [x] Adapt Bevy gizmo rendering to consume `DebugGeometry` (thin adapter, no logic change)
- [x] Serialize `DebugGeometry` into the `debug` key of the TreeMesh JSON

### M2 — Python builder API
- [ ] Implement `Tubulin` builder with `.grow()`, `.strands()`, `.mesh()`, `.debug()`, `.run()`
- [ ] Wrap `GeometryData` as `#[pyclass]` with `to_json()` method
- [ ] Expose `.vertices`, `.normals`, `.colors`, `.triangles` as numpy arrays via `rust-numpy`
- [ ] `_repr_html_()` on `GeometryData` (reuse JS viewer)

### M3 — JS debug layer rendering ✓
- [x] Parse `debug` key in JS decoder
- [x] Render lines, circles, points as Three.js objects
- [x] Checkbox UI per layer (one per key present in `debug`)

### M4 — `VisualDebug` for pipeline stages ✓
- [x] `TreeSkeletonDebugData` (skeleton layer)
- [x] `TrajectoryBuilder` (strands layer)
- [x] `GeometryData` (mesh layer)

### M5 — Export integration
- [ ] Add `to_json(include_debug: bool)` method to `GeometryData`
- [ ] Wire pipeline stages → VisualDebug → TreeMeshExporter JSON output
- [ ] Test full pipeline export

### M6 — Polish
- [ ] Consistent color conventions across stages
- [ ] Notebook example demonstrating full pipeline inspection
- [ ] Docs / docstrings on all public Python API surface

---

## Known Bugs

- **JS Decode error**: `TypeError: can't access property "length", x is undefined` — occurs when the full tree mesh is rendered in the browser. The JS decoder expects buffers to be present but some path in the evaluation chain returns undefined. Need to add assertions in decoder to catch this early and identify which expression/buffer is missing.

---

## Open Questions

- Should `Skeleton.points` include only node positions, or also interpolated spline
  points? (Affects usefulness for manual numpy analysis.)
- Default colors per stage, or user-configurable in `.debug(color=...)`?

---

## Technical Pain Points

### P1: Dual-purpose Rust structs for pyo3

Problem: We want `Seed` to be a unit struct for normal Rust usage, but pyo3 requires a struct with fields to generate Python constructors.

Attempted solutions:
- `#[cfg_attr(feature = "python", pyo3::pyclass)]` on a unit struct - doesn't work, pyo3 needs fields
- Separate `PySeed` wrapper - adds duplication
- `cfg_attr` with field - complex attribute syntax, doesn't compile well

**Solution needed**: Find a pattern that avoids duplicating the struct definition while satisfying pyo3's requirements.

### P2: Cross-boundary types

Problem: To return debug data to Python, we need types like `DebugGeometry` and `GrowConfig` to implement pyo3 traits.

Currently:
- `DebugGeometry` - contains tuples which pyo3 can't serialize directly
- `GrowConfig` - needs `#[pyclass]` and proper derives

**Solution needed**: Either wrap the return type entirely in Rust (generate JSON in Rust), or implement proper pyo3 serialization for these types.

### P3: Format strings with self fields

Problem: In `format!()` macros, you can't use `self.field` - need to use positional args or intermediate variables.

```rust
// Doesn't work:
format!("{}", self.field)

// Works:
format!("{}", &self.field)
format!("{}", self.some_method())
```

**Solution needed**: Be aware of this limitation when writing pymethods.

### P4: Method resolution with cfg features

Problem: When both `#[cfg(feature = "python")]` impl blocks and regular impl blocks define methods with the same name, Rust gets confused about which to use.

**Solution needed**: Use distinct method names for Python vs Rust (e.g., `grow_plant_py`).
