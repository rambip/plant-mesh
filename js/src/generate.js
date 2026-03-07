/**
 * Generates geometry JSON data in TreeMesh format (GEO_SPEC.md)
 */

function riceEncode(arr, k) {
  const bits = [];
  
  for (const value of arr) {
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

function createCubeData() {
  const radius = 0.5;
  const height = 1.0;
  const segments = 16;
  const scale = 256;

  const vertices = [];
  const normals = [];
  const indices = [];

  for (let i = 0; i <= segments; i++) {
    const angle = (i / segments) * Math.PI * 2;
    const x = Math.cos(angle) * radius;
    const z = Math.sin(angle) * radius;
    
    vertices.push(-height / 2, x, z);
    normals.push(0, Math.cos(angle), Math.sin(angle));
    
    vertices.push(height / 2, x, z);
    normals.push(0, Math.cos(angle), Math.sin(angle));
  }

  for (let i = 0; i < segments; i++) {
    const base = i * 2;
    indices.push(base, base + 1, base + 2);
    indices.push(base + 2, base + 1, base + 3);
  }

  const verticesX = [];
  const verticesY = [];
  const verticesZ = [];
  for (let i = 0; i < vertices.length; i += 3) {
    verticesX.push(Math.round(vertices[i] * scale));
    verticesY.push(Math.round(vertices[i + 1] * scale));
    verticesZ.push(Math.round(vertices[i + 2] * scale));
  }

  const normalsX = [];
  const normalsY = [];
  const normalsZ = [];
  for (let i = 0; i < normals.length; i += 3) {
    normalsX.push(Math.round(normals[i] * scale));
    normalsY.push(Math.round(normals[i + 1] * scale));
    normalsZ.push(Math.round(normals[i + 2] * scale));
  }

  const indicesI = [];
  const indicesJ = [];
  const indicesK = [];
  let prevI = 0, prevJ = 0, prevK = 0;
  for (let i = 0; i < indices.length; i += 3) {
    indicesI.push(indices[i] - prevI);
    indicesJ.push(indices[i + 1] - prevJ);
    indicesK.push(indices[i + 2] - prevK);
    prevI = indices[i];
    prevJ = indices[i + 1];
    prevK = indices[i + 2];
  }

  return {
    treemesh: "0.1",
    spline_convention: "reflect",
    buffers: {
      vertices_x: { k: 0, data: encodeInt32(verticesX) },
      vertices_y: { k: 0, data: encodeInt32(verticesY) },
      vertices_z: { k: 0, data: encodeInt32(verticesZ) },
      normals_x: { k: 0, data: encodeInt32(normalsX) },
      normals_y: { k: 0, data: encodeInt32(normalsY) },
      normals_z: { k: 0, data: encodeInt32(normalsZ) },
      indices_i: { k: 2, data: riceEncode(indicesI, 2) },
      indices_j: { k: 2, data: riceEncode(indicesJ, 2) },
      indices_k: { k: 2, data: riceEncode(indicesK, 2) },
    },
    outputs: {
      vertices: {
        op: "divp2",
        args: [
          {
            op: "vec3",
            args: ["vertices_x", "vertices_y", "vertices_z"]
          },
          8
        ]
      },
      normals: {
        op: "divp2",
        args: [
          {
            op: "vec3",
            args: ["normals_x", "normals_y", "normals_z"]
          },
          8
        ]
      },
      triangles: {
        op: "triangle",
        args: [
          { op: "cumsum", args: ["indices_i"] },
          { op: "cumsum", args: ["indices_j"] },
          { op: "cumsum", args: ["indices_k"] }
        ]
      }
    }
  };
}

const fs = require('fs');
const path = require('path');

const data = createCubeData();
const outputPath = path.join(__dirname, '..', 'dist', 'geometry.json');
fs.writeFileSync(outputPath, JSON.stringify(data, null, 2));
console.log('Wrote geometry.json to', outputPath);