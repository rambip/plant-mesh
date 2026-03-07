# JS Viewer Development

## Overview

JavaScript viewer for TreeMesh format (see `GEO_SPEC.md` in project root).

## Files

- `src/decoder.js` - TreeMesh decoder (Rice decode, expression eval)
- `src/generate.js` - Generates geometry.json in TreeMesh format
- `src/render.js` - Three.js renderer, fetches geometry.json
- `index.html` - Entry point, shows JSON + 3D view
- `dist/` - Built files (generated)

## Commands

```bash
# Generate geometry.json (run after any change to generate.js)
bun run js/src/generate.js

# Build render.js (run after any change to render.js or decoder.js)
bun build js/src/render.js --outdir=js/dist

# Serve (from project root)
python -m http.server 8081
```

## TreeMesh Format

See `GEO_SPEC.md` in project root. Key points:

- Buffers: `{ k: <int>, offset?: <int>, data: "<base64>" }`
- k=0: raw 4-byte little-endian ints
- k>0: Golomb-Rice encoded
- offset: added after decode (for negative values)

### Operators

- `divp2` - Divide by power of 2 (dequantization)
- `vec3` - Zip 3 int buffers to Vec3
- `cumsum` - Prefix sum (undo delta encoding)
- `triangle` - Interleave 3 index buffers
- `interleave` - Concatenate Vec3 buffers

### Outputs

- `vertices` - Buffer<Vec3>, dequantized
- `normals` - Buffer<Vec3>, dequantized  
- `triangles` - Buffer<Int>, flat [i,j,k...]

## Geometry Generation Tips

1. **Order vertices for compression**: Alternate between circles (bottom→top per segment) so consecutive indices are monotonic
2. **Delta-encode indices**: Store differences, use cumsum to recover
3. **Quantize to int then Rice**: Scale floats, round to int, use Rice encoding (k=2-4 works well for small deltas)
4. **Zigzag-encode signed deltas**: Rice coding only handles non-negatives. Map signed→unsigned via zigzag before encoding: `n >= 0 ? n*2 : -n*2-1`. Decode symmetrically. **Do not use `offset` for index deltas** — offset shifts all values uniformly and cannot handle mixed-sign deltas.

## Common Failures

### `WebGL warning: drawElementsInstanced: Indexed vertex fetch requires N vertices, but attribs only supply M`
Indices reference vertices that don't exist. Root causes:
- **Signed index deltas not zigzag-encoded**: Negative deltas passed to `riceEncode` produce garbage via JS bitwise ops on negatives. Fix: zigzag-encode before Rice, zigzag-decode after Rice.
- **Stale bundle**: `js/dist/render.js` is not rebuilt after editing `src/`. Always run `bun build js/src/render.js --outdir=js/dist` and hard-refresh (Ctrl+Shift+R).
- **Buffer length not stored**: Without `length` in the JSON buffer metadata, the Rice decoder reads padding bits as extra values, inflating buffer lengths and misaligning cumsum results.

### `Decode error: Index out of range` / `Negative index`
Thrown by assertions in `decoder.js`. Check:
- Index delta buffers have the correct `length` field in the JSON
- The encoder uses zigzag encoding for all Rice-coded signed data
- `cumsum` is applied to the correct (delta) buffer, not to already-absolute indices

### `vec3 component length mismatch` / `triangle index length mismatch`
Coordinate or index buffers have different lengths. Check that all three component buffers (x/y/z or i/j/k) are generated from the same loop and stored with correct `length` metadata.

### `Buffer "X": expected N values, got M`
Rice decoder produced wrong count. Check:
- `length` field is set correctly in the JSON for Rice-coded buffers (`k > 0`)
- The encoder's Rice output matches the declared length