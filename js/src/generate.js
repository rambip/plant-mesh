/**
 * Generates geometry JSON data in TreeMesh format (GEO_SPEC.md)
 *
 * Vertex layout: perimeter order — for each height level j, all `segments`
 * vertices around the circumference. Index = j * segments + i.
 * This keeps consecutive triangle indices close together (small deltas).
 *
 * Buffers: one x/y/z triple per height ring, interleaved in the output expr.
 */

function zigzagEncode(n) {
  return n >= 0 ? n * 2 : -n * 2 - 1;
}

function riceEncode(arr, k) {
  const bits = [];

  for (const signed of arr) {
    const value = zigzagEncode(signed);
    const quotient = value >> k;
    const remainder = value & ((1 << k) - 1);

    for (let i = 0; i < quotient; i++) {
      bits.push(1);
    }
    bits.push(0);

    for (let i = k - 1; i >= 0; i--) {
      bits.push((remainder >> i) & 1);
    }
  }

  while (bits.length % 8 !== 0) {
    bits.push(0);
  }

  const bytes = new Uint8Array(bits.length / 8);
  for (let i = 0; i < bytes.length; i++) {
    let byte = 0;
    for (let j = 0; j < 8; j++) {
      byte = (byte << 1) | bits[i * 8 + j];
    }
    bytes[i] = byte;
  }

  let binary = '';
  for (let i = 0; i < bytes.length; i++) {
    binary += String.fromCharCode(bytes[i]);
  }
  return btoa(binary);
}

function encodeInt32(arr) {
  const int32 = new Int32Array(arr);
  const bytes = new Uint8Array(int32.buffer);
  let binary = '';
  for (let i = 0; i < bytes.length; i++) {
    binary += String.fromCharCode(bytes[i]);
  }
  return btoa(binary);
}

function createCylinderData() {
  const radius = 0.5;
  const height = 1.0;
  const segments = 16;   // points around circumference
  const rings = 5;       // height levels
  const scale = 256;

  // Perimeter order: for each ring j, emit all `segments` points.
  // Vertex index = j * segments + i
  const vertices = [];  // flat [x,y,z, x,y,z, ...]
  const normals = [];

  for (let j = 0; j < rings; j++) {
    const y = -height / 2 + (height * j / (rings - 1));
    for (let i = 0; i < segments; i++) {
      const angle = (i / segments) * Math.PI * 2;
      const x = Math.cos(angle) * radius;
      const z = Math.sin(angle) * radius;
      vertices.push(y, x, z);
      normals.push(0, Math.cos(angle), Math.sin(angle));
    }
  }

  // Triangles: quad strip between ring j and j+1, wrapping in i.
  // For each quad: vertices at (j,i), (j,i+1 mod seg), (j+1,i), (j+1,i+1 mod seg)
  // Two triangles per quad, wound CCW from outside.
  const indices = [];
  for (let j = 0; j < rings - 1; j++) {
    for (let i = 0; i < segments; i++) {
      const a = j * segments + i;
      const b = j * segments + (i + 1) % segments;
      const c = (j + 1) * segments + i;
      const d = (j + 1) * segments + (i + 1) % segments;
      indices.push(a, c, b);
      indices.push(b, c, d);
    }
  }

  // --- Per-ring buffers ---
  // For each ring j: one buffer triple (ring_j_x, ring_j_y, ring_j_z).
  // Values are delta-encoded within each ring (consecutive points on a circle
  // have small coordinate deltas after quantization).
  const buffers = {};
  const ringVec3Exprs = [];   // for vertices interleave
  const ringNVec3Exprs = [];  // for normals interleave

  for (let j = 0; j < rings; j++) {
    const vx = [], vy = [], vz = [];
    const nx = [], ny = [], nz = [];

    for (let i = 0; i < segments; i++) {
      const idx = (j * segments + i) * 3;
      vx.push(Math.round(vertices[idx]     * scale));
      vy.push(Math.round(vertices[idx + 1] * scale));
      vz.push(Math.round(vertices[idx + 2] * scale));
      nx.push(Math.round(normals[idx]      * scale));
      ny.push(Math.round(normals[idx + 1]  * scale));
      nz.push(Math.round(normals[idx + 2]  * scale));
    }

    // Delta-encode each ring's coordinates (consecutive values on a circle
    // differ by small amounts, giving good Rice compression).
    const deltaEncode = (arr) => arr.map((v, i) => v - (i === 0 ? 0 : arr[i - 1]));

    const dvx = deltaEncode(vx), dvy = deltaEncode(vy), dvz = deltaEncode(vz);
    const dnx = deltaEncode(nx), dny = deltaEncode(ny), dnz = deltaEncode(nz);

    buffers[`ring${j}_x`] = { k: 2, data: riceEncode(dvx, 2), length: segments };
    buffers[`ring${j}_y`] = { k: 0, data: encodeInt32(vy), length: segments };
    buffers[`ring${j}_z`] = { k: 2, data: riceEncode(dvz, 2), length: segments };
    buffers[`ring${j}_nx`] = { k: 2, data: riceEncode(dnx, 2), length: segments };
    buffers[`ring${j}_ny`] = { k: 0, data: encodeInt32(ny), length: segments };
    buffers[`ring${j}_nz`] = { k: 2, data: riceEncode(dnz, 2), length: segments };

    // x and z use cumsum to undo delta; y is stored absolute (constant per ring).
    ringVec3Exprs.push({
      op: 'vec3',
      args: [
        { op: 'cumsum', args: [`ring${j}_x`] },
        `ring${j}_y`,
        { op: 'cumsum', args: [`ring${j}_z`] },
      ],
    });
    ringNVec3Exprs.push({
      op: 'vec3',
      args: [
        { op: 'cumsum', args: [`ring${j}_nx`] },
        `ring${j}_ny`,
        { op: 'cumsum', args: [`ring${j}_nz`] },
      ],
    });
  }

  // --- Index buffers ---
  // Delta-encode i, j, k components of triangle indices independently.
  const indicesI = [], indicesJ = [], indicesK = [];
  let prevI = 0, prevJ = 0, prevK = 0;
  for (let t = 0; t < indices.length; t += 3) {
    indicesI.push(indices[t]     - prevI);
    indicesJ.push(indices[t + 1] - prevJ);
    indicesK.push(indices[t + 2] - prevK);
    prevI = indices[t];
    prevJ = indices[t + 1];
    prevK = indices[t + 2];
  }
  const nTriangles = indices.length / 3;
  buffers['indices_i'] = { k: 2, data: riceEncode(indicesI, 2), length: nTriangles };
  buffers['indices_j'] = { k: 2, data: riceEncode(indicesJ, 2), length: nTriangles };
  buffers['indices_k'] = { k: 2, data: riceEncode(indicesK, 2), length: nTriangles };

  return {
    treemesh: "0.1",
    spline_convention: "reflect",
    buffers,
    outputs: {
      vertices: {
        op: 'divp2',
        args: [{ op: 'interleave', args: ringVec3Exprs }, 8],
      },
      normals: {
        op: 'divp2',
        args: [{ op: 'interleave', args: ringNVec3Exprs }, 8],
      },
      triangles: {
        op: 'triangle',
        args: [
          { op: 'cumsum', args: ['indices_i'] },
          { op: 'cumsum', args: ['indices_j'] },
          { op: 'cumsum', args: ['indices_k'] },
        ],
      },
    },
  };
}

const fs = require('fs');
const path = require('path');

const data = createCylinderData();
const outputPath = path.join(__dirname, '..', 'dist', 'geometry.json');
fs.writeFileSync(outputPath, JSON.stringify(data, null, 2));
console.log('Wrote geometry.json to', outputPath);
