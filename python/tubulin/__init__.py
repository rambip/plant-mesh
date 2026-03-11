# tubulin — procedural plant mesh generator
# Python bindings are provided by the Rust extension _tubulin (built with maturin).
from ._tubulin import (
    Seed,
    PlantNode,
    Skeleton,
    VolumetricTree,
    TreeMesh,
    DebugData,
    debug_to_json,
    demo_mesh,
)
