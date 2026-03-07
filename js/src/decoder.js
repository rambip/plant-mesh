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
      scope[name] = this._decodeBuffer(buf.data, buf.k, offset);
    }

    try {
      const vertices = this._eval(outputs.vertices, scope);
      const normals = this._eval(outputs.normals, scope);
      const triangles = this._eval(outputs.triangles, scope);

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

  _decodeBuffer(base64, k, offset = 0) {
    const binary = atob(base64);
    const byteArray = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i++) {
      byteArray[i] = binary.charCodeAt(i);
    }

    let result;
    if (k === 0) {
      result = new Int32Array(byteArray.buffer, byteArray.byteOffset, byteArray.length / 4);
    } else {
      result = this._riceDecode(byteArray, k);
    }

    if (offset !== 0) {
      for (let i = 0; i < result.length; i++) {
        result[i] += offset;
      }
    }

    return result;
  }

  _riceDecode(byteArray, k) {
    const result = [];
    let bitBuffer = 0;
    let bitCount = 0;
    let byteIndex = 0;

    const readBit = () => {
      if (bitCount === 0) {
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

    while (byteIndex < byteArray.length || bitCount > 0) {
      let quotient = 0;
      while (readBit() === 1) {
        quotient++;
      }
      const remainder = readBits(k);
      result.push((quotient << k) + remainder);
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
        case 'interleave':
          return this._interleave(args);
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
    const result = [];
    for (let i = 0; i < x.length; i++) {
      result.push({ x: x[i], y: y[i], z: z[i] });
    }
    return result;
  }

  _triangle(i, j, k) {
    const result = new Uint32Array(i.length * 3);
    for (let n = 0; n < i.length; n++) {
      result[n * 3] = i[n];
      result[n * 3 + 1] = j[n];
      result[n * 3 + 2] = k[n];
    }
    return result;
  }

  _interleave(buffers) {
    const totalLength = buffers.reduce((sum, b) => sum + b.length, 0);
    const result = [];
    for (const buffer of buffers) {
      for (const v of buffer) {
        result.push(v);
      }
    }
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