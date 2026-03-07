/**
 * Generates geometry JSON data in TreeMesh format (GEO_SPEC.md)
 *
 * Mesh 1 — longitudinal spine cylinder (z axis):
 *   8 splines, one per meridian angle, each running bottom→top along z.
 *   4 control points per spline; radius is smaller in the middle (waist).
 *   After evaluation, interleave (zip) all 8 splines so the vertex layout is
 *   ring-major: [s0[0],s1[0],...,s7[0], s0[1],s1[1],...,s7[1], ...]
 *   Triangles connect adjacent rings, wrapping around the 8 splines.
 *
 * Mesh 2 — latitudinal ring cylinder (x axis):
 *   5 rings, each a 4-control-point circular spline evaluated at 16 points.
 *   concat-ed into one vertex buffer.
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

    for (let i = 0; i < quotient; i++) bits.push(1);
    bits.push(0);
    for (let i = k - 1; i >= 0; i--) bits.push((remainder >> i) & 1);
  }

  while (bits.length % 8 !== 0) bits.push(0);

  const bytes = new Uint8Array(bits.length / 8);
  for (let i = 0; i < bytes.length; i++) {
    let byte = 0;
    for (let j = 0; j < 8; j++) byte = (byte << 1) | bits[i * 8 + j];
    bytes[i] = byte;
  }

  let binary = '';
  for (let i = 0; i < bytes.length; i++) binary += String.fromCharCode(bytes[i]);
  return btoa(binary);
}

function encodeInt32(arr) {
  const int32 = new Int32Array(arr);
  const bytes = new Uint8Array(int32.buffer);
  let binary = '';
  for (let i = 0; i < bytes.length; i++) binary += String.fromCharCode(bytes[i]);
  return btoa(binary);
}

function encodeFloat32(arr) {
  const f32 = new Float32Array(arr);
  const bytes = new Uint8Array(f32.buffer);
  let binary = '';
  for (let i = 0; i < bytes.length; i++) binary += String.fromCharCode(bytes[i]);
  return btoa(binary);
}

function createCylinderData() {
  const posScale = 256;  // quantization scale for positions
  const tScale = 1024;   // quantization scale for spline t parameter

  const buffers = {};

  // -------------------------------------------------------------------------
  // Mesh 1: longitudinal spine cylinder, z axis, centered at origin
  //
  // 8 splines at angles 0°, 45°, 90°, ..., 315°. Each spline has 4 control
  // points along z, with a smaller radius in the middle (waist):
  //   cp[0]: z=-0.5, r=0.5   (bottom)
  //   cp[1]: z=-1/6, r=0.25  (lower waist)
  //   cp[2]: z=+1/6, r=0.25  (upper waist)
  //   cp[3]: z=+0.5, r=0.5   (top)
  //
  // After evaluation at nSamples t-values, interleave (zip) all 8 splines so
  // the vertex layout is ring-major: for sample i, all 8 spline points appear
  // consecutively → vertex index = i * nSplines + s.
  // -------------------------------------------------------------------------

  const nSplines = 8;
  const nSamples = 8;   // vertices per spline (rings in the output mesh)
  const ncp1 = 4;       // control points per spline

  // z positions and radii for the 4 control points
  const cpZ  = [-0.5, -1/6,  1/6,  0.5];
  const cpR  = [ 0.5,  0.25, 0.25, 0.5];

  // t values: uniform in [0, ncp1-1], shared across all 8 splines
  const spine_t = [];
  for (let i = 0; i < nSamples; i++)
    spine_t.push(Math.round((i / (nSamples - 1)) * (ncp1 - 1) * tScale));
  buffers['spine_t'] = { k: 0, data: encodeInt32(spine_t), length: nSamples };

  const spineSplineExprs = [];
  const spineNSplineExprs = [];

  for (let s = 0; s < nSplines; s++) {
    const angle = (s / nSplines) * Math.PI * 2;
    const cos = Math.cos(angle), sin = Math.sin(angle);

    // Control points: (x, y, z) = (cos*r, sin*r, z)
    const cpx = cpR.map(r => Math.round(cos * r * posScale));
    const cpy = cpR.map(r => Math.round(sin * r * posScale));
    const cpz = cpZ.map(z => Math.round(z * posScale));

    // Normals point radially outward: (cos, sin, 0) at each cp
    const cnx = cpR.map(() => Math.round(cos * posScale));
    const cny = cpR.map(() => Math.round(sin * posScale));
    const cnz = cpR.map(() => 0);

    const delta = (arr) => arr.map((v, i) => v - (i === 0 ? 0 : arr[i - 1]));

    buffers[`sp${s}_x`] = { k: 2, data: riceEncode(delta(cpx), 2), length: ncp1 };
    buffers[`sp${s}_y`] = { k: 2, data: riceEncode(delta(cpy), 2), length: ncp1 };
    buffers[`sp${s}_z`] = { k: 2, data: riceEncode(delta(cpz), 2), length: ncp1 };
    buffers[`sp${s}_nx`] = { k: 0, data: encodeInt32(cnx), length: ncp1 };
    buffers[`sp${s}_ny`] = { k: 0, data: encodeInt32(cny), length: ncp1 };
    buffers[`sp${s}_nz`] = { k: 0, data: encodeInt32(cnz), length: ncp1 };

    const cpVec3 = (xb, yb, zb) => ({
      op: 'divp2',
      args: [{ op: 'vec3', args: [
        { op: 'cumsum', args: [xb] },
        { op: 'cumsum', args: [yb] },
        { op: 'cumsum', args: [zb] },
      ]}, 8],
    });

    const tExpr = { op: 'divp2', args: ['spine_t', 10] };

    spineSplineExprs.push({ op: 'spline', args: [cpVec3(`sp${s}_x`, `sp${s}_y`, `sp${s}_z`), tExpr] });
    spineNSplineExprs.push({ op: 'spline', args: [cpVec3(`sp${s}_nx`, `sp${s}_ny`, `sp${s}_nz`), tExpr] });
  }

  // Triangles for mesh 1.
  // After interleave, vertex index = i * nSplines + s.
  // Quads between ring i and i+1, wrapping around splines.
  const spineIndices = (offset) => {
    const idx = [];
    for (let i = 0; i < nSamples - 1; i++) {
      for (let s = 0; s < nSplines; s++) {
        const a = offset + i * nSplines + s;
        const b = offset + i * nSplines + (s + 1) % nSplines;
        const c = offset + (i + 1) * nSplines + s;
        const d = offset + (i + 1) * nSplines + (s + 1) % nSplines;
        idx.push(a, c, b);
        idx.push(b, c, d);
      }
    }
    return idx;
  };

  const vertsPerSpineCyl = nSamples * nSplines;  // 8*8 = 64

  // -------------------------------------------------------------------------
  // Mesh 2: latitudinal ring cylinder, x axis, offset to x=[0.6, 1.6]
  // 5 rings × 4 control points, evaluated at 16 samples per ring.
  // -------------------------------------------------------------------------

  const rings = 5;
  const segments = 16;
  const ncp2 = 5;  // 4 control points + wrap
  const radii = [0.5, 0.35, 0.25, 0.35, 0.5];

  const tValues = [];
  for (let i = 0; i < segments; i++)
    tValues.push(Math.round((i / segments) * 4 * tScale));
  buffers['ring_t'] = { k: 0, data: encodeInt32(tValues), length: segments };

  const ring2SplineExprs = [];
  const ring2NSplineExprs = [];

  for (let j = 0; j < rings; j++) {
    const x = 0.6 + j / (rings - 1);
    const radius = radii[j];

    const cpx = [], cpy = [], cpz = [];
    const cnx = [], cny = [], cnz = [];
    for (let i = 0; i < ncp2; i++) {
      const angle = (i / 4) * Math.PI * 2;
      cpx.push(Math.round(x                        * posScale));
      cpy.push(Math.round(Math.cos(angle) * radius * posScale));
      cpz.push(Math.round(Math.sin(angle) * radius * posScale));
      cnx.push(0);
      cny.push(Math.round(Math.cos(angle) * posScale));
      cnz.push(Math.round(Math.sin(angle) * posScale));
    }

    const delta = (arr) => arr.map((v, i) => v - (i === 0 ? 0 : arr[i - 1]));

    buffers[`ring${j}_cp_x`] = { k: 0, data: encodeInt32(cpx), length: ncp2 };
    buffers[`ring${j}_cp_y`] = { k: 0, data: encodeInt32(cpy), length: ncp2 };
    buffers[`ring${j}_cp_z`] = { k: 2, data: riceEncode(delta(cpz), 2), length: ncp2 };
    buffers[`ring${j}_cn_x`] = { k: 0, data: encodeInt32(cnx), length: ncp2 };
    buffers[`ring${j}_cn_y`] = { k: 0, data: encodeInt32(cny), length: ncp2 };
    buffers[`ring${j}_cn_z`] = { k: 2, data: riceEncode(delta(cnz), 2), length: ncp2 };

    const cpVec3r = (xbuf, ybuf, zbuf) => ({
      op: 'divp2',
      args: [{ op: 'vec3', args: [xbuf, ybuf, { op: 'cumsum', args: [zbuf] }]}, 8],
    });

    const tExpr2 = { op: 'divp2', args: ['ring_t', 10] };

    ring2SplineExprs.push({ op: 'spline', args: [cpVec3r(`ring${j}_cp_x`, `ring${j}_cp_y`, `ring${j}_cp_z`), tExpr2] });
    ring2NSplineExprs.push({ op: 'spline', args: [cpVec3r(`ring${j}_cn_x`, `ring${j}_cn_y`, `ring${j}_cn_z`), tExpr2] });
  }

  const ringIndices = (offset) => {
    const idx = [];
    for (let j = 0; j < rings - 1; j++) {
      for (let i = 0; i < segments; i++) {
        const a = offset + j * segments + i;
        const b = offset + j * segments + (i + 1) % segments;
        const c = offset + (j + 1) * segments + i;
        const d = offset + (j + 1) * segments + (i + 1) % segments;
        idx.push(a, c, b);
        idx.push(b, c, d);
      }
    }
    return idx;
  };

  const vertsPerRingCyl = rings * segments;  // 5*16 = 80

  const indices = [
    ...spineIndices(0),
    ...ringIndices(vertsPerSpineCyl),
  ];

  // Delta-encode all triangle indices together
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
        op: 'concat',
        args: [
          // Mesh 1: interleave (zip) all 8 splines so layout is ring-major
          { op: 'interleave', args: spineSplineExprs },
          // Mesh 2: concat rings sequentially
          { op: 'concat', args: ring2SplineExprs },
        ],
      },
      normals: {
        op: 'concat',
        args: [
          { op: 'interleave', args: spineNSplineExprs },
          { op: 'concat', args: ring2NSplineExprs },
        ],
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
