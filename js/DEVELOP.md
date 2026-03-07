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
4. **Use offset for negatives**: If data ranges -128 to +127, add offset=128 before Rice encoding