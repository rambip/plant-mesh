/**
 * @module decoder
 */

/**
 * @typedef {Object} GeometryOptions
 * @property {Float32Array} [vertexBuffer]
 * @property {Float32Array} [normalBuffer]
 * @property {Uint32Array} [indexBuffer]
 */

/**
 * Represents a 3D geometry with vertex, normal, and index buffers.
 * @class
 */
class Geometry {
  /**
   * @type {Float32Array}
   */
  vertexBuffer;

  /**
   * @type {Float32Array}
   */
  normalBuffer;

  /**
   * @type {Uint32Array}
   */
  indexBuffer;

  /**
   * @param {GeometryOptions} [options]
   */
  constructor(options = {}) {
    this.vertexBuffer = options.vertexBuffer || new Float32Array(0);
    this.normalBuffer = options.normalBuffer || new Float32Array(0);
    this.indexBuffer = options.indexBuffer || new Uint32Array(0);
  }
}

/**
 * @typedef {Object} TreeMeshInput
 * @property {string} treemesh
 * @property {string} spline_convention
 * @property {Object} buffers
 * @property {Object} outputs
 */

/**
 * @class
 */
class TreeMeshDecoder {
  /**
   * @param {TreeMeshInput} jsonInput
   * @returns {Geometry}
   */
  decode(jsonInput) {
    const { buffers, outputs } = jsonInput;

    const scope = {};
    for (const [name, buf] of Object.entries(buffers)) {
      const offset = buf.offset !== undefined ? buf.offset : 0;
      const decoded = this._decodeBuffer(buf.data, buf.k, offset, buf.length);
      if (buf.length !== undefined && decoded.length !== buf.length) {
        throw new Error(`Buffer "${name}": expected ${buf.length} values, got ${decoded.length}. ` +
          `Check Rice k parameter and zigzag encoding in the encoder.`);
      }
      scope[name] = decoded;
    }

    try {
      const vertices = this._eval(outputs.vertices, scope);
      const normals = this._eval(outputs.normals, scope);
      const triangles = this._eval(outputs.triangles, scope);

      if (vertices.length !== normals.length) {
        throw new Error(`Vertex/normal count mismatch: ${vertices.length} vertices vs ${normals.length} normals.`);
      }

      const nVertices = vertices.length;
      const maxIndex = Math.max(...triangles);
      if (maxIndex >= nVertices) {
        throw new Error(
          `Index out of range: max index ${maxIndex} >= vertex count ${nVertices}. ` +
          `Likely cause: signed index deltas encoded without zigzag, or cumsum applied to wrong buffer.`
        );
      }
      const minIndex = Math.min(...triangles);
      if (minIndex < 0) {
        throw new Error(
          `Negative index ${minIndex} in triangle buffer. ` +
          `Likely cause: signed deltas not zigzag-encoded, causing cumsum to produce negative values.`
        );
      }
      if (triangles.length % 3 !== 0) {
        throw new Error(`Triangle buffer length ${triangles.length} is not a multiple of 3.`);
      }

      const vertexBuffer = this._vec3ToInterleaved(vertices);
      const normalBuffer = this._vec3ToInterleaved(normals);

      return new Geometry({
        vertexBuffer,
        normalBuffer,
        indexBuffer: triangles,
      });
    } catch (e) {
      console.error('Decode error:', e);
      throw e;
    }
  }

  _decodeBuffer(base64, k, offset = 0, length = undefined) {
    const binary = atob(base64);
    const byteArray = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i++) {
      byteArray[i] = binary.charCodeAt(i);
    }

    let result;
    if (k === 0) {
      const count = length !== undefined ? length : byteArray.length / 4;
      result = new Int32Array(byteArray.buffer, byteArray.byteOffset, count);
    } else {
      result = this._riceDecode(byteArray, k, length);
    }

    if (offset !== 0) {
      for (let i = 0; i < result.length; i++) {
        result[i] += offset;
      }
    }

    return result;
  }

  _riceDecode(byteArray, k, length = undefined) {
    const result = [];
    let bitBuffer = 0;
    let bitCount = 0;
    let byteIndex = 0;

    const readBit = () => {
      if (bitCount === 0) {
        if (byteIndex >= byteArray.length) return 0;
        bitBuffer = byteArray[byteIndex++];
        bitCount = 8;
      }
      const bit = (bitBuffer >> 7) & 1;
      bitBuffer <<= 1;
      bitCount--;
      return bit;
    };

    const readBits = (n) => {
      let value = 0;
      for (let i = 0; i < n; i++) {
        value = (value << 1) | readBit();
      }
      return value;
    };

    const limit = length !== undefined ? length : Infinity;
    while ((byteIndex < byteArray.length || bitCount > 0) && result.length < limit) {
      let quotient = 0;
      while (readBit() === 1) {
        quotient++;
      }
      const remainder = readBits(k);
      const zigzag = (quotient << k) + remainder;
      // zigzag decode: even -> positive, odd -> negative
      const signed = (zigzag & 1) === 0 ? zigzag >> 1 : -(zigzag >> 1) - 1;
      result.push(signed);
    }

    return new Int32Array(result);
  }

  _eval(expr, scope) {
    if (typeof expr === 'string') {
      return scope[expr];
    }

    if (typeof expr === 'number') {
      return expr;
    }

    const op = expr.op;
    
    if (expr.args) {
      const args = expr.args.map(a => this._eval(a, scope));
      
      switch (op) {
        case 'cumsum':
          return this._cumsum(args[0]);
        case 'divp2':
          return this._divp2(args[0], args[1]);
        case 'vec3':
          return this._vec3(args[0], args[1], args[2]);
        case 'triangle':
          return this._triangle(args[0], args[1], args[2]);
        case 'concat':
          return this._concat(args);
        case 'interleave':
          return this._interleave(args);
        case 'spline':
          return this._spline(args[0], args[1]);
        default:
          throw new Error(`Unknown operator: ${op}`);
      }
    }

    throw new Error(`Missing args for operator: ${op}`);
  }

  _cumsum(buffer) {
    const result = new Int32Array(buffer.length);
    let sum = 0;
    for (let i = 0; i < buffer.length; i++) {
      sum += buffer[i];
      result[i] = sum;
    }
    return result;
  }

  _divp2(buffer, exponent) {
    const divisor = Math.pow(2, exponent);
    if (buffer instanceof Int32Array) {
      const result = new Float32Array(buffer.length);
      for (let i = 0; i < buffer.length; i++) {
        result[i] = buffer[i] / divisor;
      }
      return result;
    } else if (Array.isArray(buffer) && buffer.length > 0 && typeof buffer[0] === 'object' && 'x' in buffer[0]) {
      return buffer.map(v => ({
        x: v.x / divisor,
        y: v.y / divisor,
        z: v.z / divisor,
      }));
    }
    throw new Error('Unsupported buffer type for divp2');
  }

  _vec3(x, y, z) {
    if (x.length !== y.length || x.length !== z.length) {
      throw new Error(`vec3 component length mismatch: x=${x.length}, y=${y.length}, z=${z.length}. ` +
        `All three coordinate buffers must have the same length.`);
    }
    const result = [];
    for (let i = 0; i < x.length; i++) {
      result.push({ x: x[i], y: y[i], z: z[i] });
    }
    return result;
  }

  _triangle(i, j, k) {
    if (i.length !== j.length || i.length !== k.length) {
      throw new Error(`triangle index length mismatch: i=${i.length}, j=${j.length}, k=${k.length}. ` +
        `All three index buffers must have the same length (one entry per triangle).`);
    }
    const result = new Uint32Array(i.length * 3);
    for (let n = 0; n < i.length; n++) {
      result[n * 3] = i[n];
      result[n * 3 + 1] = j[n];
      result[n * 3 + 2] = k[n];
    }
    return result;
  }

  // Centripetal Catmull-Rom spline. Matches extended_catmull_spline() in
  // crates/plant-core/src/meshing/algorithms.rs.
  // points: Array of {x,y,z}, length K >= 2
  // times:  Float32Array, values in [0, K-1]
  _spline(points, times) {
    const K = points.length;
    if (K < 2) throw new Error(`spline: need at least 2 control points, got ${K}`);

    // Reflect border: virtual points at index -1 and K
    const p = (i) => {
      if (i === -1)  return { x: 2*points[0].x - points[1].x,
                              y: 2*points[0].y - points[1].y,
                              z: 2*points[0].z - points[1].z };
      if (i === K)   return { x: 2*points[K-1].x - points[K-2].x,
                              y: 2*points[K-1].y - points[K-2].y,
                              z: 2*points[K-1].z - points[K-2].z };
      return points[i];
    };

    const dist = (a, b) => {
      const dx = b.x-a.x, dy = b.y-a.y, dz = b.z-a.z;
      return Math.sqrt(Math.sqrt(dx*dx + dy*dy + dz*dz)); // sqrt of length = centripetal
    };

    const lerp3 = (a, b, t) => ({
      x: a.x + (b.x - a.x) * t,
      y: a.y + (b.y - a.y) * t,
      z: a.z + (b.z - a.z) * t,
    });

    const evalAt = (t) => {
      // Clamp t=K-1 to last segment
      let i0, r;
      if (t >= K - 1) { i0 = K - 2; r = 1.0; }
      else            { i0 = Math.floor(t); r = t - i0; }

      const pts = [p(i0 - 1), p(i0), p(i0 + 1), p(i0 + 2)];

      // Centripetal knot sequence
      const knots = [0, 0, 0, 0];
      for (let i = 1; i < 4; i++) {
        knots[i] = knots[i-1] + dist(pts[i-1], pts[i]);
      }

      const tk = knots[1] + (knots[2] - knots[1]) * r;

      const ratio = (i, j) => {
        const d = knots[j] - knots[i];
        return d === 0 ? 0 : (tk - knots[i]) / d;
      };

      // Level 1: 3 points
      const a = [
        lerp3(pts[0], pts[1], ratio(0, 1)),
        lerp3(pts[1], pts[2], ratio(1, 2)),
        lerp3(pts[2], pts[3], ratio(2, 3)),
      ];
      // Level 2: 2 points
      const b1 = lerp3(a[0], a[1], ratio(0, 2));
      const b2 = lerp3(a[1], a[2], ratio(1, 3));
      // Level 3: result
      return lerp3(b1, b2, ratio(1, 2));
    };

    return Array.from(times).map(evalAt);
  }

  _concat(buffers) {
    // Sequentially appends buffers of the same type (Int32Array, Float32Array, or Vec3 array).
    if (buffers[0] instanceof Int32Array || buffers[0] instanceof Float32Array) {
      const total = buffers.reduce((s, b) => s + b.length, 0);
      const result = new buffers[0].constructor(total);
      let offset = 0;
      for (const b of buffers) { result.set(b, offset); offset += b.length; }
      return result;
    }
    // Vec3 array
    const result = [];
    for (const buf of buffers) for (const v of buf) result.push(v);
    return result;
  }

  _interleave(buffers) {
    // Elementwise zip: all buffers must have equal length N.
    // Output length = N * len(buffers), elements alternating.
    const N = buffers[0].length;
    for (const buf of buffers) {
      if (buf.length !== N) {
        throw new Error(`interleave: all buffers must have equal length, got ${N} and ${buf.length}`);
      }
    }
    if (buffers[0] instanceof Int32Array || buffers[0] instanceof Float32Array) {
      const result = new buffers[0].constructor(N * buffers.length);
      for (let i = 0; i < N; i++)
        for (let b = 0; b < buffers.length; b++)
          result[i * buffers.length + b] = buffers[b][i];
      return result;
    }
    // Vec3 array
    const result = [];
    for (let i = 0; i < N; i++)
      for (const buf of buffers) result.push(buf[i]);
    return result;
  }

  _vec3ToInterleaved(vec3Array) {
    const result = new Float32Array(vec3Array.length * 3);
    for (let i = 0; i < vec3Array.length; i++) {
      result[i * 3] = vec3Array[i].x;
      result[i * 3 + 1] = vec3Array[i].y;
      result[i * 3 + 2] = vec3Array[i].z;
    }
    return result;
  }
}

module.exports = { Geometry, TreeMeshDecoder };