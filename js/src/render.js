import * as THREE from 'https://unpkg.com/three@0.160.0/build/three.module.js';
import { TreeMeshDecoder } from './decoder.js';

/**
 * @typedef {Object} CubeData
 * @property {number} n_vertices
 * @property {number} n_triangles
 * @property {string} vertices
 * @property {string} normals
 * @property {string} indices
 */

function init(geometryData) {
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

  const treeMeshDecoder = new TreeMeshDecoder();
  const geometry = treeMeshDecoder.decode(geometryData);

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

  const material = new THREE.MeshBasicMaterial({
    color: 0x00ff88,
    wireframe: true,
  });

  const mesh = new THREE.Mesh(threeGeometry, material);
  scene.add(mesh);

  function animate() {
    requestAnimationFrame(animate);
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

async function main() {
  const response = await fetch('./dist/geometry.json');
  const geometryData = await response.json();
  
  document.getElementById('json').textContent = JSON.stringify(geometryData, null, 2);
  
  init(geometryData);
}

main();