# TreeMesh Intermediate Representation — Language Specification

Version: 0.1  
Status: Draft

---

## Overview

A TreeMesh file is a JSON document that encodes a tree mesh as a small **program** operating on compressed integer buffers. A compliant decoder executes the program to produce a vertex position array and a triangle index array, suitable for direct upload to a GPU.

The format has two parts:

- **`buffers`** — named, compressed integer arrays
- **`outputs`** — a program with two required outputs: `vertices` and `triangles`

---

## Top-Level Structure

```json
{
  "treemesh": "0.1",
  "spline_convention": "reflect",
  "buffers": { ... },
  "outputs": {
    "vertices":  <expr>,
    "triangles": <expr>
  }
}
```

### Fields

| Field | Type | Required | Description |
|---|---|---|---|
| `treemesh` | string | yes | Format version. Must be `"0.1"` |
| `spline_convention` | string | yes | Border handling for spline evaluation. Must be `"reflect"` in this version |
| `buffers` | object | yes | Named compressed buffers (see Buffers) |
| `outputs.vertices` | expr | yes | Expression evaluating to `Buffer<Vec3>` |
| `outputs.triangles` | expr | yes | Expression evaluating to `Buffer<Int>`, flat `[i0,j0,k0, i1,j1,k1,...]` |

---

## Buffers

Each buffer is a named Rice-decoded integer array. Buffers are the only source of raw data in the program; all other values are derived from them.

```json
"buffers": {
  "<name>": {
    "k":      <int>,      // Golomb-Rice parameter, k >= 0
    "offset": <int>,      // Optional offset added after decoding (default: 0)
    "data":   "<base64>"  // Rice-encoded, base64-encoded byte string
  }
}
```

### Decoding

1. Base64-decode `data` to a byte array
2. Apply Golomb-Rice decoding with parameter `k` to produce a sequence of non-negative integers
3. Add `offset` to each value (if present, default is 0)

The result is a `Buffer<Int>` bound to `<name>` in the program's scope.

The `offset` field allows encoding negative values by shifting them into the non-negative range before Rice encoding. For example, if the original data ranges from -128 to +127, use `"offset": 128` to shift it to 0-255.

### Golomb-Rice Encoding Convention

- All encoded values are **non-negative** integers. Signed values must be zigzag-encoded before Rice encoding: `n >= 0 → 2n`, `n < 0 → -2n - 1`
- Bit order within bytes: **most significant bit first**
- The Rice quotient is encoded in **unary** (sequence of 1-bits terminated by a 0-bit), followed by `k` remainder bits

---

## Types

The language has three types:

| Type | Description |
|---|---|
| `Int` | A 32-bit signed integer |
| `Scalar` | A 32-bit float |
| `Vec3` | A triple of `Scalar` values `(x, y, z)` |

All types exist as **buffers** (fixed-length sequences). All operators are vectorized: they operate elementwise over their buffer arguments. Buffer lengths must match unless otherwise specified.

---

## Expressions

An expression is either:

- A **string** — a reference to a named buffer from `buffers`, evaluating to `Buffer<Int>`
- An **object** with an `"op"` field — an operator application

Arguments to operators are always written as:

- `"args": [...]` — a positional list of arguments

---

## Operators

### `cumsum` — Prefix Sum

```json
{ "op": "cumsum", "args": [<Buffer<Int>>] }
```

Computes the cumulative sum of the input buffer. Used to undo delta encoding.

```
out[0] = in[0]
out[i] = out[i-1] + in[i]
```

**Input:** `Buffer<Int>`  
**Output:** `Buffer<Int>`, same length

---

### `divp2` — Divide by Power of Two

```json
{ "op": "divp2", "args": [<Buffer<Int> | Buffer<Vec3>>, <int>] }
```

Divides every element by `2^args[1]`. `args[1]` is a **literal integer constant**, not a buffer reference.

If the input is `Buffer<Int>`, the output is `Buffer<Scalar>`.  
If the input is `Buffer<Vec3>`, the output is `Buffer<Vec3>` with each component divided.

```
out[i] = in[i] / (2 ^ args[1])
```

**Input:** `Buffer<Int>` or `Buffer<Vec3>`  
**Output:** `Buffer<Scalar>` or `Buffer<Vec3>`

---

### `vec3` — Construct Vec3 Buffer

```json
{ "op": "vec3", "args": [<Buffer<Int>>, <Buffer<Int>>, <Buffer<Int>>] }
```

Zips three integer buffers into a buffer of Vec3. The three buffers are `x`, `y`, `z` in order. Typically followed by `divp2` for dequantization.

```
out[i] = (args[0][i], args[1][i], args[2][i])
```

**Input:** three `Buffer<Int>` of equal length  
**Output:** `Buffer<Vec3>`

---

### `spline` — Evaluate Centripetal Catmull-Rom Spline

```json
{ "op": "spline", "args": [<Buffer<Vec3>>, <Buffer<Scalar>>] }
```

Evaluates a centripetal Catmull-Rom spline defined by the control points in `args[0]`, at each parameter value in `args[1]`.

- `args[0]` is a `Buffer<Vec3>` of `K` control points (`K >= 2`)
- `args[1]` is a `Buffer<Scalar>` of query parameter values, each in `[0, K-1]`
- `floor(t)` identifies the segment index `i0`; the fractional part `r = t - i0` is the local parameter in `[0, 1]`

For each query `t`, evaluation uses the 4-point window `[p[i0-1], p[i0], p[i0+1], p[i0+2]]`. The border convention `"reflect"` synthesizes out-of-range points:

```
p[-1] = 2*p[0]   - p[1]
p[K]  = 2*p[K-1] - p[K-2]
```

The knot sequence for the 4-point window is **centripetal** (knot spacing = sqrt of chord length):

```
knot[0] = 0
knot[i] = knot[i-1] + sqrt(length(p[i] - p[i-1]))   for i in 1..3
```

The local parameter is mapped into the knot interval `[knot[1], knot[2]]`:

```
t_knot = knot[1] + (knot[2] - knot[1]) * r
```

Evaluation uses de Boor's algorithm with the centripetal knots (3 levels of lerp). See `crates/plant-core/src/meshing/algorithms.rs:extended_catmull_spline` for the reference implementation.

```
out[i] = spline(args[0], args[1][i])
```

**Input:** `Buffer<Vec3>` (length K), `Buffer<Scalar>` (length N)  
**Output:** `Buffer<Vec3>` (length N)

---

### `concat` — Concatenate Buffers

```json
{ "op": "concat", "args": [<Buffer<T>>, <Buffer<T>>, ...] }
```

Appends any number of buffers of the same type sequentially. The output contains all elements of `args[0]`, followed by all elements of `args[1]`, and so on. Works for any buffer type (`Int`, `Scalar`, `Vec3`).

Used as the final step for `vertices`: collects per-spline position buffers into one contiguous vertex buffer.

```
out = [args[0][0], ..., args[0][N0-1], args[1][0], ..., args[1][N1-1], ...]
```

**Input:** one or more `Buffer<T>` of the same type  
**Output:** `Buffer<T>`

---

### `interleave` — Zip Buffers Elementwise

```json
{ "op": "interleave", "args": [<Buffer<T>>, <Buffer<T>>, ...] }
```

Interleaves elements from all input buffers elementwise. All inputs must have the same length N. The output has length `N * len(args)`, with elements from each buffer alternating.

```
out = [args[0][0], args[1][0], ..., args[0][1], args[1][1], ...]
```

**Input:** one or more `Buffer<T>` of equal length  
**Output:** `Buffer<T>` of length `N * len(args)`

---

### `triangle` — Interleave Triangle Index Buffers

```json
{ "op": "triangle", "args": [<Buffer<Int>>, <Buffer<Int>>, <Buffer<Int>>] }
```

Interleaves three index buffers `i`, `j`, `k` (each of length T) into a flat buffer of length `3T`, in the order `[i0, j0, k0, i1, j1, k1, ...]`. The result is a standard flat triangle index buffer for GPU upload.

```
out[3*n]   = args[0][n]
out[3*n+1] = args[1][n]
out[3*n+2] = args[2][n]
```

**Input:** three `Buffer<Int>` of equal length T  
**Output:** `Buffer<Int>` of length 3T

---

## Vertex Index Convention

Triangle index buffers reference vertices by their position in the output of the top-level `vertices` expression (after `interleave`). If spline 0 contributes `N0` vertices and spline 1 contributes `N1` vertices, then spline 1's vertices are indexed starting at `N0`. The encoder is responsible for writing delta-encoded triangle indices consistent with this ordering.

---

## Worked Example (3 splines)

```json
{
  "treemesh": "0.1",
  "spline_convention": "reflect",

  "buffers": {
    "spline_0_x": { "k": 3, "data": "..." },
    "spline_0_y": { "k": 3, "data": "..." },
    "spline_0_z": { "k": 3, "data": "..." },
    "spline_0_t": { "k": 2, "data": "..." },
    "spline_1_x": { "k": 3, "data": "..." },
    "spline_1_y": { "k": 3, "data": "..." },
    "spline_1_z": { "k": 3, "data": "..." },
    "spline_1_t": { "k": 2, "data": "..." },
    "spline_2_x": { "k": 3, "data": "..." },
    "spline_2_y": { "k": 3, "data": "..." },
    "spline_2_z": { "k": 3, "data": "..." },
    "spline_2_t": { "k": 2, "data": "..." },
    "triangle_i": { "k": 1, "data": "..." },
    "triangle_j": { "k": 1, "data": "..." },
    "triangle_k": { "k": 1, "data": "..." }
  },

  "outputs": {
    "vertices": {
      "op": "interleave",
      "args": [
        {
          "op": "spline",
          "args": [
            { "op": "divp2", "args": [{ "op": "vec3", "args": ["spline_0_x", "spline_0_y", "spline_0_z"] }, 4096] },
            { "op": "divp2", "args": ["spline_0_t", 1024] }
          ]
        },
        {
          "op": "spline",
          "args": [
            { "op": "divp2", "args": [{ "op": "vec3", "args": ["spline_1_x", "spline_1_y", "spline_1_z"] }, 4096] },
            { "op": "divp2", "args": ["spline_1_t", 1024] }
          ]
        },
        {
          "op": "spline",
          "args": [
            { "op": "divp2", "args": [{ "op": "vec3", "args": ["spline_2_x", "spline_2_y", "spline_2_z"] }, 4096] },
            { "op": "divp2", "args": ["spline_2_t", 1024] }
          ]
        }
      ]
    },

    "triangles": {
      "op": "triangle",
      "args": [
        { "op": "cumsum", "args": ["triangle_i"] },
        { "op": "cumsum", "args": ["triangle_j"] },
        { "op": "cumsum", "args": ["triangle_k"] }
      ]
    }
  }
}
```

---

## Decoder Algorithm (Pseudocode)

```
function decode(doc):
  // 1. Decode all buffers
  scope = {}
  for name, buf in doc.buffers:
    scope[name] = rice_decode(base64_decode(buf.data), buf.k)

  // 2. Evaluate outputs
  vertices  = eval(doc.outputs.vertices, scope)   // → Float32Array (x,y,z interleaved)
  triangles = eval(doc.outputs.triangles, scope)  // → Uint32Array  (i,j,k interleaved)
  return { vertices, triangles }

function eval(expr, scope):
  if expr is string:
    return scope[expr]                            // Buffer<Int>

  switch expr.op:
    case "cumsum":
      return prefix_sum(eval(arg(expr, 0), scope))

    case "divp2":
      return map(eval(arg(expr, 0), scope), x => x / (1 << expr.arg_1))

    case "vec3":
      x, y, z = eval(arg(expr, 0)), eval(arg(expr, 1)), eval(arg(expr, 2))
      return zip(x, y, z)

    case "spline":
      points = eval(expr.args[0], scope)          // Buffer<Vec3>
      times  = eval(expr.args[1], scope)          // Buffer<Scalar>
      return times.map(t => eval_spline(points, t))

    case "concat":
      return append(expr.args.map(a => eval(a, scope)))

    case "interleave":
      buffers = expr.args.map(a => eval(a, scope))
      return zip_elementwise(buffers)

    case "triangle":
      i, j, k = eval(expr.args[0]), eval(expr.args[1]), eval(expr.args[2])
      return interleave3(i, j, k)
```

---

## Open Questions / Future Extensions

- **Quantization precision** — `divp2` exponents are currently chosen per-file by the encoder. A future version may standardize them or derive them from bounding box metadata.
- **Multiple spline conventions** — only `"reflect"` is specified in v0.1. Future versions may add `"clamp"` and `"natural"`.
- **LOD** — no level-of-detail mechanism is specified. A future version may allow multiple `outputs` blocks at different resolutions.
- **Normals** — not part of the format. The decoder is expected to recompute normals from the triangle mesh.
