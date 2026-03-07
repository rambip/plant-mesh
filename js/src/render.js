import * as THREE from 'https://unpkg.com/three@0.160.0/build/three.module.js';
import { GeoVM3d } from './decoder.js';

/**
 * @typedef {Object} CubeData
 * @property {number} n_vertices
 * @property {number} n_triangles
 * @property {string} vertices
 * @property {string} normals
 * @property {string} indices
 */

/**
 * Encodes a number array to base64 little-endian
 * @param {number[]} arr
 * @returns {string}
 */
function encodeFloat32(arr) {
  const float32 = new Float32Array(arr);
  const bytes = new Uint8Array(float32.buffer);
  let binary = '';
  for (let i = 0; i < bytes.length; i++) {
    binary += String.fromCharCode(bytes[i]);
  }
  return btoa(binary);
}

/**
 * Encodes a number array to base64 little-endian
 * @param {number[]} arr
 * @returns {string}
 */
function encodeUint16(arr) {
  const uint16 = new Uint16Array(arr);
  const bytes = new Uint8Array(uint16.buffer);
  let binary = '';
  for (let i = 0; i < bytes.length; i++) {
    binary += String.fromCharCode(bytes[i]);
  }
  return btoa(binary);
}

function createCubeJson() {
  const s = 0.5;

  const vertices = [
    -s, -s, -s,  s, -s, -s,  s,  s, -s, -s,  s, -s,
    -s, -s,  s,  s, -s,  s,  s,  s,  s, -s,  s,  s,
    -s, -s, -s, -s,  s, -s, -s,  s,  s, -s, -s,  s,
     s, -s, -s,  s,  s, -s,  s,  s,  s,  s, -s,  s,
    -s, -s, -s, -s, -s,  s,  s, -s,  s,  s, -s, -s,
    -s,  s, -s, -s,  s,  s,  s,  s,  s,  s,  s, -s,
  ];

  const normals = [
     0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,
     0,  0,  1,  0,  0,  1,  0,  0,  1,  0,  0,  1,
    -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,  0,
     1,  0,  0,  1,  0,  0,  1,  0,  0,  1,  0,  0,
     0, -1,  0,  0, -1,  0,  0, -1,  0,  0, -1,  0,
     0,  1,  0,  0,  1,  0,  0,  1,  0,  0,  1,  0,
  ];

  const indices = [
    0, 1, 2, 0, 2, 3,
    4, 5, 6, 4, 6, 7,
    8, 9, 10, 8, 10, 11,
    12, 13, 14, 12, 14, 15,
    16, 17, 18, 16, 18, 19,
    20, 21, 22, 20, 22, 23,
  ];

  /** @type {CubeData} */
  return {
    n_vertices: 24,
    n_triangles: 12,
    vertices: encodeFloat32(vertices),
    normals: encodeFloat32(normals),
    indices: encodeUint16(indices),
  };
}

function init() {
  const container = document.getElementById('main');
  if (!container) {
    console.error('Element with id "main" not found');
    return;
  }

  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x222222);

  const camera = new THREE.PerspectiveCamera(
    75,
    container.clientWidth / container.clientHeight,
    0.1,
    1000
  );
  camera.position.z = 2;

  const renderer = new THREE.WebGLRenderer();
  renderer.setSize(container.clientWidth, container.clientHeight);
  container.appendChild(renderer.domElement);

  const ambientLight = new THREE.AmbientLight(0x404040);
  scene.add(ambientLight);

  const directionalLight = new THREE.DirectionalLight(0xffffff, 1);
  directionalLight.position.set(1, 1, 1);
  scene.add(directionalLight);

  const cubeJson = createCubeJson();
  const geoVM3d = new GeoVM3d();
  const geometry = geoVM3d.decode(cubeJson);

  const threeGeometry = new THREE.BufferGeometry();
  threeGeometry.setAttribute(
    'position',
    new THREE.BufferAttribute(geometry.vertexBuffer, 3)
  );
  threeGeometry.setAttribute(
    'normal',
    new THREE.BufferAttribute(geometry.normalBuffer, 3)
  );
  threeGeometry.setIndex(new THREE.BufferAttribute(geometry.indexBuffer, 1));

  const material = new THREE.MeshPhongMaterial({
    color: 0x00ff88,
    side: THREE.DoubleSide,
  });

  const mesh = new THREE.Mesh(threeGeometry, material);
  scene.add(mesh);

  function animate() {
    requestAnimationFrame(animate);
    mesh.rotation.x += 0.01;
    mesh.rotation.y += 0.01;
    renderer.render(scene, camera);
  }

  animate();

  window.addEventListener('resize', () => {
    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(container.clientWidth, container.clientHeight);
  });
}

init();
