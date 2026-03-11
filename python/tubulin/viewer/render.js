var __create = Object.create;
var __getProtoOf = Object.getPrototypeOf;
var __defProp = Object.defineProperty;
var __getOwnPropNames = Object.getOwnPropertyNames;
var __hasOwnProp = Object.prototype.hasOwnProperty;
function __accessProp(key) {
  return this[key];
}
var __toESMCache_node;
var __toESMCache_esm;
var __toESM = (mod, isNodeMode, target) => {
  var canCache = mod != null && typeof mod === "object";
  if (canCache) {
    var cache = isNodeMode ? __toESMCache_node ??= new WeakMap : __toESMCache_esm ??= new WeakMap;
    var cached = cache.get(mod);
    if (cached)
      return cached;
  }
  target = mod != null ? __create(__getProtoOf(mod)) : {};
  const to = isNodeMode || !mod || !mod.__esModule ? __defProp(target, "default", { value: mod, enumerable: true }) : target;
  for (let key of __getOwnPropNames(mod))
    if (!__hasOwnProp.call(to, key))
      __defProp(to, key, {
        get: __accessProp.bind(mod, key),
        enumerable: true
      });
  if (canCache)
    cache.set(mod, to);
  return to;
};
var __commonJS = (cb, mod) => () => (mod || cb((mod = { exports: {} }).exports, mod), mod.exports);

// js/src/decoder.js
var require_decoder = __commonJS((exports, module) => {
  class Geometry {
    vertexBuffer;
    normalBuffer;
    colorBuffer;
    indexBuffer;
    debugLayers;
    constructor(options = {}) {
      this.vertexBuffer = options.vertexBuffer || new Float32Array(0);
      this.normalBuffer = options.normalBuffer || new Float32Array(0);
      this.colorBuffer = options.colorBuffer || new Float32Array(0);
      this.indexBuffer = options.indexBuffer || new Uint32Array(0);
      this.debugLayers = options.debugLayers || {};
    }
  }

  class TreeMeshDecoder {
    decode(jsonInput) {
      const { buffers, outputs, debug } = jsonInput;
      const scope = {};
      for (const [name, buf] of Object.entries(buffers)) {
        const offset = buf.offset !== undefined ? buf.offset : 0;
        const decoded = this._decodeBuffer(buf.data, buf.k, offset, buf.length);
        if (buf.length !== undefined && decoded.length !== buf.length) {
          throw new Error(`Buffer "${name}": expected ${buf.length} values, got ${decoded.length}. ` + `Check Rice k parameter and zigzag encoding in the encoder.`);
        }
        scope[name] = decoded;
      }
      try {
        const vertices = this._eval(outputs.vertices, scope);
        const normals = this._eval(outputs.normals, scope);
        const colors = outputs.colors ? this._eval(outputs.colors, scope) : null;
        const triangles = this._eval(outputs.triangles, scope);
        if (vertices.length !== normals.length) {
          throw new Error(`Vertex/normal count mismatch: ${vertices.length} vertices vs ${normals.length} normals.`);
        }
        const nVertices = vertices.length;
        const maxIndex = Math.max(...triangles);
        if (maxIndex >= nVertices) {
          throw new Error(`Index out of range: max index ${maxIndex} >= vertex count ${nVertices}. ` + `Likely cause: signed index deltas encoded without zigzag, or cumsum applied to wrong buffer.`);
        }
        const minIndex = Math.min(...triangles);
        if (minIndex < 0) {
          throw new Error(`Negative index ${minIndex} in triangle buffer. ` + `Likely cause: signed deltas not zigzag-encoded, causing cumsum to produce negative values.`);
        }
        if (triangles.length % 3 !== 0) {
          throw new Error(`Triangle buffer length ${triangles.length} is not a multiple of 3.`);
        }
        const vertexBuffer = this._vec3ToInterleaved(vertices);
        const normalBuffer = this._vec3ToInterleaved(normals);
        const colorBuffer = colors ? this._rgbaToInterleaved(colors) : new Float32Array(0);
        const debugLayers = this._decodeDebug(debug, scope);
        return new Geometry({
          vertexBuffer,
          normalBuffer,
          colorBuffer,
          indexBuffer: triangles,
          debugLayers
        });
      } catch (e) {
        console.error("Decode error:", e);
        throw e;
      }
    }
    _decodeBuffer(base64, k, offset = 0, length = undefined) {
      const binary = atob(base64);
      const byteArray = new Uint8Array(binary.length);
      for (let i = 0;i < binary.length; i++) {
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
        for (let i = 0;i < result.length; i++) {
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
          if (byteIndex >= byteArray.length)
            return 0;
          bitBuffer = byteArray[byteIndex++];
          bitCount = 8;
        }
        const bit = bitBuffer >> 7 & 1;
        bitBuffer <<= 1;
        bitCount--;
        return bit;
      };
      const readBits = (n) => {
        let value = 0;
        for (let i = 0;i < n; i++) {
          value = value << 1 | readBit();
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
        const signed = (zigzag & 1) === 0 ? zigzag >> 1 : -(zigzag >> 1) - 1;
        result.push(signed);
      }
      return new Int32Array(result);
    }
    _eval(expr, scope) {
      if (typeof expr === "string") {
        return scope[expr];
      }
      if (typeof expr === "number") {
        return expr;
      }
      if (Array.isArray(expr)) {
        throw new Error(`Expected expression object, got array: ${JSON.stringify(expr).slice(0, 100)}`);
      }
      if (typeof expr !== "object" || expr === null) {
        throw new Error(`Invalid expression: ${expr}, expected object with 'op' field`);
      }
      const op = expr.op;
      if (!op) {
        throw new Error(`Expression missing 'op' field: ${JSON.stringify(expr).slice(0, 200)}`);
      }
      if (!expr.args) {
        throw new Error(`Operator '${op}' missing 'args' field`);
      }
      if (!Array.isArray(expr.args)) {
        throw new Error(`Operator '${op}' 'args' must be an array, got: ${typeof expr.args}`);
      }
      const args = expr.args.map((a) => this._eval(a, scope));
      switch (op) {
        case "cumsum":
          return this._cumsum(args[0]);
        case "divp2":
          return this._divp2(args[0], args[1]);
        case "vec3":
          return this._vec3(args[0], args[1], args[2]);
        case "vec4":
          return this._vec4(args[0], args[1], args[2], args[3]);
        case "triangle":
          return this._triangle(args[0], args[1], args[2]);
        case "concat":
          return this._concat(args);
        case "interleave":
          return this._interleave(args);
        case "spline":
          return this._spline(args[0], args[1]);
        default:
          throw new Error(`Unknown operator: ${op}`);
      }
    }
    _cumsum(buffer) {
      const result = new Int32Array(buffer.length);
      let sum = 0;
      for (let i = 0;i < buffer.length; i++) {
        sum += buffer[i];
        result[i] = sum;
      }
      return result;
    }
    _divp2(buffer, exponent) {
      const divisor = Math.pow(2, exponent);
      if (buffer instanceof Int32Array) {
        const result = new Float32Array(buffer.length);
        for (let i = 0;i < buffer.length; i++) {
          result[i] = buffer[i] / divisor;
        }
        return result;
      } else if (Array.isArray(buffer) && buffer.length > 0 && typeof buffer[0] === "object" && "x" in buffer[0]) {
        return buffer.map((v) => ({
          x: v.x / divisor,
          y: v.y / divisor,
          z: v.z / divisor,
          w: v.w !== undefined ? v.w / divisor : undefined
        }));
      }
      throw new Error("Unsupported buffer type for divp2");
    }
    _vec3(x, y, z) {
      if (x.length !== y.length || x.length !== z.length) {
        throw new Error(`vec3 component length mismatch: x=${x.length}, y=${y.length}, z=${z.length}. ` + `All three coordinate buffers must have the same length.`);
      }
      const result = [];
      for (let i = 0;i < x.length; i++) {
        result.push({ x: x[i], y: y[i], z: z[i] });
      }
      return result;
    }
    _vec4(x, y, z, w) {
      if (x.length !== y.length || x.length !== z.length || x.length !== w.length) {
        throw new Error(`vec4 component length mismatch: x=${x.length}, y=${y.length}, z=${z.length}, w=${w.length}. ` + `All four component buffers must have the same length.`);
      }
      const result = [];
      for (let i = 0;i < x.length; i++) {
        result.push({ x: x[i], y: y[i], z: z[i], w: w[i] });
      }
      return result;
    }
    _triangle(i, j, k) {
      if (i.length !== j.length || i.length !== k.length) {
        throw new Error(`triangle index length mismatch: i=${i.length}, j=${j.length}, k=${k.length}. ` + `All three index buffers must have the same length (one entry per triangle).`);
      }
      const result = new Uint32Array(i.length * 3);
      for (let n = 0;n < i.length; n++) {
        result[n * 3] = i[n];
        result[n * 3 + 1] = j[n];
        result[n * 3 + 2] = k[n];
      }
      return result;
    }
    _spline(points, times) {
      const K = points.length;
      if (K < 2)
        throw new Error(`spline: need at least 2 control points, got ${K}`);
      const p = (i) => {
        if (i === -1)
          return {
            x: 2 * points[0].x - points[1].x,
            y: 2 * points[0].y - points[1].y,
            z: 2 * points[0].z - points[1].z
          };
        if (i === K)
          return {
            x: 2 * points[K - 1].x - points[K - 2].x,
            y: 2 * points[K - 1].y - points[K - 2].y,
            z: 2 * points[K - 1].z - points[K - 2].z
          };
        return points[i];
      };
      const dist = (a, b) => {
        const dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
        return Math.sqrt(Math.sqrt(dx * dx + dy * dy + dz * dz));
      };
      const lerp3 = (a, b, t) => ({
        x: a.x + (b.x - a.x) * t,
        y: a.y + (b.y - a.y) * t,
        z: a.z + (b.z - a.z) * t
      });
      const evalAt = (t) => {
        let i0, r;
        if (t >= K - 1) {
          i0 = K - 2;
          r = 1;
        } else {
          i0 = Math.floor(t);
          r = t - i0;
        }
        const pts = [p(i0 - 1), p(i0), p(i0 + 1), p(i0 + 2)];
        const knots = [0, 0, 0, 0];
        for (let i = 1;i < 4; i++) {
          knots[i] = knots[i - 1] + dist(pts[i - 1], pts[i]);
        }
        const tk = knots[1] + (knots[2] - knots[1]) * r;
        const ratio = (i, j) => {
          const d = knots[j] - knots[i];
          return d === 0 ? 0 : (tk - knots[i]) / d;
        };
        const a = [
          lerp3(pts[0], pts[1], ratio(0, 1)),
          lerp3(pts[1], pts[2], ratio(1, 2)),
          lerp3(pts[2], pts[3], ratio(2, 3))
        ];
        const b1 = lerp3(a[0], a[1], ratio(0, 2));
        const b2 = lerp3(a[1], a[2], ratio(1, 3));
        return lerp3(b1, b2, ratio(1, 2));
      };
      return Array.from(times).map(evalAt);
    }
    _concat(buffers) {
      if (buffers[0] instanceof Int32Array || buffers[0] instanceof Float32Array) {
        const total = buffers.reduce((s, b) => s + b.length, 0);
        const result2 = new buffers[0].constructor(total);
        let offset = 0;
        for (const b of buffers) {
          result2.set(b, offset);
          offset += b.length;
        }
        return result2;
      }
      const result = [];
      for (const buf of buffers)
        for (const v of buf)
          result.push(v);
      return result;
    }
    _interleave(buffers) {
      const N = buffers[0].length;
      for (const buf of buffers) {
        if (buf.length !== N) {
          throw new Error(`interleave: all buffers must have equal length, got ${N} and ${buf.length}`);
        }
      }
      if (buffers[0] instanceof Int32Array || buffers[0] instanceof Float32Array) {
        const result2 = new buffers[0].constructor(N * buffers.length);
        for (let i = 0;i < N; i++)
          for (let b = 0;b < buffers.length; b++)
            result2[i * buffers.length + b] = buffers[b][i];
        return result2;
      }
      const result = [];
      for (let i = 0;i < N; i++)
        for (const buf of buffers)
          result.push(buf[i]);
      return result;
    }
    _vec3ToInterleaved(vec3Array) {
      const result = new Float32Array(vec3Array.length * 3);
      for (let i = 0;i < vec3Array.length; i++) {
        result[i * 3] = vec3Array[i].x;
        result[i * 3 + 1] = vec3Array[i].y;
        result[i * 3 + 2] = vec3Array[i].z;
      }
      return result;
    }
    _rgbaToInterleaved(rgbaArray) {
      const result = new Float32Array(rgbaArray.length * 4);
      for (let i = 0;i < rgbaArray.length; i++) {
        const c = rgbaArray[i];
        result[i * 4] = c.x;
        result[i * 4 + 1] = c.y;
        result[i * 4 + 2] = c.z;
        result[i * 4 + 3] = c.w !== undefined ? c.w : 1;
      }
      return result;
    }
    _decodeDebug(debug, scope) {
      const layers = {};
      if (!debug)
        return layers;
      for (const [layerName, layerData] of Object.entries(debug)) {
        const layer = {};
        if (layerData.points) {
          const positions = this._eval(layerData.points.positions, scope);
          const colors = layerData.points.colors ? this._eval(layerData.points.colors, scope) : null;
          layer.points = {
            positions: this._vec3ToInterleaved(positions),
            colors: colors ? this._rgbaToInterleaved(colors) : null
          };
        }
        if (layerData.lines) {
          const starts = this._eval(layerData.lines.starts, scope);
          const ends = this._eval(layerData.lines.ends, scope);
          const colors = layerData.lines.colors ? this._eval(layerData.lines.colors, scope) : null;
          layer.lines = {
            starts: this._vec3ToInterleaved(starts),
            ends: this._vec3ToInterleaved(ends),
            colors: colors ? this._rgbaToInterleaved(colors) : null
          };
        }
        if (Object.keys(layer).length > 0) {
          layers[layerName] = layer;
        }
      }
      return layers;
    }
  }
  module.exports = { Geometry, TreeMeshDecoder };
});

// js/src/render.js
var import_decoder = __toESM(require_decoder(), 1);
import * as THREE from "https://unpkg.com/three@0.160.0/build/three.module.js";
function init(geometryData) {
  const container = document.getElementById("main");
  if (!container) {
    console.error('Element with id "main" not found');
    return;
  }
  const scene = new THREE.Scene;
  scene.background = new THREE.Color(2236962);
  const camera = new THREE.PerspectiveCamera(75, container.clientWidth / container.clientHeight, 0.1, 1000);
  camera.position.z = 2;
  const renderer = new THREE.WebGLRenderer;
  renderer.setSize(container.clientWidth, container.clientHeight);
  container.appendChild(renderer.domElement);
  const ambientLight = new THREE.AmbientLight(4210752);
  scene.add(ambientLight);
  const directionalLight = new THREE.DirectionalLight(16777215, 1);
  directionalLight.position.set(1, 1, 1);
  scene.add(directionalLight);
  const treeMeshDecoder = new import_decoder.TreeMeshDecoder;
  const geometry = treeMeshDecoder.decode(geometryData);
  const threeGeometry = new THREE.BufferGeometry;
  threeGeometry.setAttribute("position", new THREE.BufferAttribute(geometry.vertexBuffer, 3));
  threeGeometry.setAttribute("normal", new THREE.BufferAttribute(geometry.normalBuffer, 3));
  if (geometry.colorBuffer && geometry.colorBuffer.length > 0) {
    threeGeometry.setAttribute("color", new THREE.BufferAttribute(geometry.colorBuffer, 4));
  }
  threeGeometry.setIndex(new THREE.BufferAttribute(geometry.indexBuffer, 1));
  const material = new THREE.MeshBasicMaterial({
    color: 16777215,
    vertexColors: true,
    wireframe: true
  });
  const mesh = new THREE.Mesh(threeGeometry, material);
  scene.add(mesh);
  const debugGroups = {};
  for (const [layerName, layer] of Object.entries(geometry.debugLayers)) {
    const group = new THREE.Group;
    group.visible = false;
    if (layer.points && layer.points.positions) {
      const positions = layer.points.positions;
      const colors = layer.points.colors;
      const nPoints = positions.length / 3;
      const pointsGeo = new THREE.BufferGeometry;
      pointsGeo.setAttribute("position", new THREE.BufferAttribute(positions, 3));
      if (colors && colors.length > 0) {
        pointsGeo.setAttribute("color", new THREE.BufferAttribute(colors, 4));
      }
      const pointsMat = new THREE.PointsMaterial({
        size: 0.05,
        vertexColors: colors && colors.length > 0
      });
      const pointsObj = new THREE.Points(pointsGeo, pointsMat);
      group.add(pointsObj);
    }
    if (layer.lines && layer.lines.starts && layer.lines.ends) {
      const starts = layer.lines.starts;
      const ends = layer.lines.ends;
      const colors = layer.lines.colors;
      const linePositions = new Float32Array(starts.length + ends.length);
      linePositions.set(starts, 0);
      linePositions.set(ends, starts.length);
      const lineGeo = new THREE.BufferGeometry;
      lineGeo.setAttribute("position", new THREE.BufferAttribute(linePositions, 3));
      const lineColors = colors ? new Float32Array(starts.length / 3 * 2 * 4) : null;
      if (colors) {
        for (let i = 0;i < starts.length / 3; i++) {
          lineColors[i * 8] = colors[i * 4];
          lineColors[i * 8 + 1] = colors[i * 4 + 1];
          lineColors[i * 8 + 2] = colors[i * 4 + 2];
          lineColors[i * 8 + 3] = colors[i * 4 + 3];
          lineColors[i * 8 + 4] = colors[i * 4];
          lineColors[i * 8 + 5] = colors[i * 4 + 1];
          lineColors[i * 8 + 6] = colors[i * 4 + 2];
          lineColors[i * 8 + 7] = colors[i * 4 + 3];
        }
        lineGeo.setAttribute("color", new THREE.BufferAttribute(lineColors, 4));
      }
      const lineMat = new THREE.LineBasicMaterial({
        vertexColors: colors !== null
      });
      const linesObj = new THREE.LineSegments(lineGeo, lineMat);
      group.add(linesObj);
    }
    scene.add(group);
    debugGroups[layerName] = group;
  }
  if (Object.keys(debugGroups).length > 0) {
    const controlsDiv = document.createElement("div");
    controlsDiv.style.cssText = "position:absolute;top:10px;right:10px;background:rgba(0,0,0,0.7);padding:10px;color:white;font-family:sans-serif;";
    for (const layerName of Object.keys(debugGroups)) {
      const label = document.createElement("label");
      label.style.cssText = "display:block;margin:5px 0;";
      const checkbox = document.createElement("input");
      checkbox.type = "checkbox";
      checkbox.onchange = (e) => {
        debugGroups[layerName].visible = e.target.checked;
      };
      label.appendChild(checkbox);
      label.appendChild(document.createTextNode(" " + layerName));
      controlsDiv.appendChild(label);
    }
    document.body.appendChild(controlsDiv);
  }
  function animate() {
    requestAnimationFrame(animate);
    mesh.rotation.y += 0.01;
    renderer.render(scene, camera);
  }
  animate();
  window.addEventListener("resize", () => {
    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(container.clientWidth, container.clientHeight);
  });
}
async function main() {
  const response = await fetch("./dist/geometry.json");
  const geometryData = await response.json();
  document.getElementById("json").textContent = JSON.stringify(geometryData, null, 2);
  init(geometryData);
}
main();
