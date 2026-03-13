# Tubulin Roadmap

## Direction

The next phase focuses on making Tubulin pleasant and reliable for Python-first
users (notebooks, teaching demos, quick experiments) while keeping the Rust core
as the single source of truth for geometry logic.

Priorities:

1. Better Python API ergonomics
2. More inspectable intermediate data
3. Cleaner Rust/Python boundary (`pyo3` isolation)
4. Performance improvements across generation + export
5. Strong documentation and examples

---

## Product Goals

### G1 — Python API that feels intentional

- Pipeline objects should be easy to discover and chain.
- Method names should match user intent (`grow`, `build_mesh`, `to_json`, `debug`).
- Configuration should be kwargs-first with stable defaults.
- Errors should be actionable from notebooks (no cryptic panic surfaces).

### G2 — Richer access to intermediate data

- Expose enough structure for learning/debugging without leaking internal clutter.
- Keep return types numeric and analysis-friendly (`numpy` arrays / small typed wrappers).
- Make stage outputs self-describing (counts, shapes, debug metadata).

### G3 — Robust Rust/Python boundary

- Keep all `pyo3` usage in `crates/tubulin-core/src/python.rs`.
- Core modules (`growing`, `meshing`, `export`) remain `pyo3`-free.
- Keep conversion/adaptation logic near bindings, not mixed in domain code.

### G4 — Faster end-to-end iteration

- Improve performance for typical notebook loops (generate, inspect, render).
- Reduce avoidable allocations and redundant recomputation.
- Keep JSON/viewer paths efficient for large meshes.

### G5 — Documentation that teaches and unblocks

- Clear API docs with examples at each stage.
- Up-to-date architecture and format docs.
- One canonical notebook walkthrough from seed to final mesh.

---

## Milestones

## M1 — Python API Cleanup

Target: coherent, stable stage API for notebook users.

- [ ] Audit current Python surface and remove naming inconsistencies.
- [ ] Standardize stage classes and constructor patterns:
  - `Seed`
  - `PlantNode`
  - `Skeleton`
  - `VolumetricTree`
  - `TreeMesh`
  - `DebugData`
- [ ] Normalize kwargs for stage configuration with explicit defaults.
- [ ] Add clear repr/help text for each exposed class.
- [ ] Ensure every stage has predictable `debug()` and/or `to_json()` behavior.

Deliverable: stable API shape documented with short code snippets.

## M2 — Expose More Intermediate Data

Target: users can inspect internals without patching Rust.

- [ ] `Skeleton`:
  - [ ] expose node positions/radii/orientations
  - [ ] expose parent/children topology in analysis-friendly form
- [ ] `VolumetricTree`:
  - [ ] expose trajectory counts and per-node particle mapping
  - [ ] expose sampled trajectory points for inspection
- [ ] `TreeMesh`:
  - [ ] expose vertices, normals, colors, triangles with consistent shapes
  - [ ] expose lightweight mesh stats (counts, bounds)
- [ ] Keep API read-only by default to preserve Rust invariants.

Deliverable: notebook users can inspect all major stage outputs directly.

## M3 — pyo3 Boundary Isolation

Target: no `pyo3` usage outside `python.rs`.

- [ ] Remove remaining `#[cfg(feature = "python")]` binding code from core modules.
- [ ] Move Python-only helper methods into `python.rs` wrappers.
- [ ] Keep core structs pure Rust + serde + math traits only.
- [ ] Add CI check/documented convention to prevent pyo3 leakage.

Deliverable: architecture clean split between domain logic and bindings.

## M4 — Performance Pass

Target: noticeably faster notebook workflows.

- [ ] Profile hotspots for representative cases (small/medium/large trees).
- [ ] Optimize geometry generation path:
  - [ ] reduce temporary allocations in trajectory/mesh steps
  - [ ] avoid unnecessary clones across stage transitions
- [ ] Optimize export path:
  - [ ] avoid repeated quantization work when unchanged
  - [ ] minimize JSON assembly overhead where practical
- [ ] Optimize viewer payload path for large debug layers.
- [ ] Add benchmark script and record baseline vs improved timings.

Deliverable: measured speedup and documented benchmark results.

## M5 — Docs and Examples

Target: users can succeed without reading Rust internals.

- [ ] Python API guide:
  - [ ] quickstart
  - [ ] stage-by-stage walkthrough
  - [ ] debug visualization usage
- [ ] Update `ARCHITECTURE.md` and keep it synchronized with implementation.
- [ ] Keep `GEO_SPEC.md` aligned with encoder/decoder behavior.
- [ ] Add notebook example(s):
  - [ ] minimal end-to-end
  - [ ] intermediate data inspection
  - [ ] debug-layer exploration
- [ ] Add troubleshooting section for common decode/render errors.

Deliverable: docs + notebook examples usable by new contributors/students.

---

## Cross-Cutting Engineering Tasks

- [ ] Improve error quality across Rust ↔ JS ↔ Python boundaries.
- [ ] Keep viewer bundle workflow explicit and reproducible.
- [ ] Add regression tests for known decode failure modes.
- [ ] Track API changes in changelog notes for notebook users.

---

## Definition of Done (for this roadmap)

This roadmap is complete when:

- Python API is coherent and documented.
- Intermediate data access covers skeleton, strands, and mesh inspection.
- `pyo3` usage is isolated to `python.rs`.
- Performance improvements are benchmarked and measurable.
- Docs and examples are current and sufficient for first-time users.

---

## Out of Scope (for now)

- Major TreeMesh format version bump.
- Full GUI editor for all generation parameters.
- Non-notebook frontend redesign.

---

## End-of-Roadmap Cleanup

- [ ] Remove `bevy` feature/dependencies from `tubulin-core` by moving Bevy-specific derives/wrappers into `crates/bevy-demo`, so core stays engine-agnostic.
