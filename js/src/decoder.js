/**
 * @module decoder
 */

/**
 * @typedef {Object} GeometryOptions
 * @property {Float32Array} [vertexBuffer]
 * @property {Float32Array} [normalBuffer]
 * @property {Uint16Array} [indexBuffer]
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
   * @type {Uint16Array}
   */
  indexBuffer;

  /**
   * @param {GeometryOptions} [options]
   */
  constructor(options = {}) {
    this.vertexBuffer = options.vertexBuffer || new Float32Array(0);
    this.normalBuffer = options.normalBuffer || new Float32Array(0);
    this.indexBuffer = options.indexBuffer || new Uint16Array(0);
  }
}

/**
 * @typedef {Object} GeoVM3dInput
 * @property {number} n_vertices
 * @property {number} n_triangles
 * @property {string} vertices
 * @property {string} normals
 * @property {string} indices
 */

/**
 * @class
 */
class GeoVM3d {
  /**
   * @param {GeoVM3dInput} jsonInput
   * @returns {Geometry}
   */
  decode(jsonInput) {
    const { n_vertices, n_triangles, vertices, normals, indices } = jsonInput;

    const vertexData = this._base64ToFloat32Array(vertices, n_vertices * 3);
    const normalData = this._base64ToFloat32Array(normals, n_vertices * 3);
    const indexData = this._base64ToUint16Array(indices, n_triangles * 3);

    return new Geometry({
      vertexBuffer: vertexData,
      normalBuffer: normalData,
      indexBuffer: indexData,
    });
  }

  /**
   * @param {string} base64
   * @param {number} length
   * @returns {Float32Array}
   */
  _base64ToFloat32Array(base64, length) {
    const binary = atob(base64);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i++) {
      bytes[i] = binary.charCodeAt(i);
    }
    return new Float32Array(bytes.buffer, bytes.byteOffset, length);
  }

  /**
   * @param {string} base64
   * @param {number} length
   * @returns {Uint16Array}
   */
  _base64ToUint16Array(base64, length) {
    const binary = atob(base64);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i++) {
      bytes[i] = binary.charCodeAt(i);
    }
    return new Uint16Array(bytes.buffer, bytes.byteOffset, length);
  }
}

module.exports = { Geometry, GeoVM3d };
