# tubulin — procedural plant mesh generator
# Python bindings are provided by the Rust extension _tubulin (built with maturin).
from ._tubulin import build_demo_tree

import os


def _get_viewer_html(json_data: str) -> str:
    """Generate HTML to display the tree mesh."""

    viewer_dir = os.path.dirname(os.path.abspath(__file__))
    decoder_path = os.path.join(viewer_dir, "viewer", "decoder.js")
    render_path = os.path.join(viewer_dir, "viewer", "render.js")

    with open(decoder_path, "r") as f:
        decoder_js = f.read()

    with open(render_path, "r") as f:
        render_js = f.read()

    return f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <style>
    body {{ margin: 0; padding: 0; font-family: sans-serif; }}
    #main {{ width: 100%; height: 500px; }}
  </style>
</head>
<body>
  <div id="main"></div>
  <script type="module">
    {decoder_js}
  </script>
  <script type="module">
    {render_js}
    initTreeViewer({json_data}, 'main');
  </script>
</body>
</html>"""


class TreeMesh:
    """Wrapper around GeometryData with HTML representation."""

    def __init__(self, geometry_data):
        self._data = geometry_data

    def to_json(self, include_debug: bool = False) -> str:
        """Export to TreeMesh JSON format."""
        return self._data.to_json(include_debug)

    def _repr_html_(self) -> str:
        """Render in Jupyter notebook."""
        json_str = self._data.to_json(False)
        return _get_viewer_html(json_str)

    @property
    def points(self):
        """Vertex positions as numpy array."""
        return self._data.points

    @property
    def normals(self):
        """Vertex normals as numpy array."""
        return self._data.normals

    @property
    def colors(self):
        """Vertex colors as numpy array."""
        return self._data.colors

    @property
    def triangles(self):
        """Triangle indices as numpy array."""
        return self._data.triangles


def grow(**config) -> TreeMesh:
    """Generate a tree with given config (placeholder - uses demo tree for now)."""
    birth_power = config.get("birth_power", 0.5)
    data = build_demo_tree(birth_power)
    return TreeMesh(data)
