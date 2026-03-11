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
  if (geometry.colorBuffer && geometry.colorBuffer.length > 0) {
    threeGeometry.setAttribute(
      'color',
      new THREE.BufferAttribute(geometry.colorBuffer, 4)
    );
  }
  threeGeometry.setIndex(new THREE.BufferAttribute(geometry.indexBuffer, 1));

  const material = new THREE.MeshBasicMaterial({
    color: 0xffffff,
    vertexColors: true,
    wireframe: true,
  });

  const mesh = new THREE.Mesh(threeGeometry, material);
  scene.add(mesh);

  const debugGroups = {};
  for (const [layerName, layer] of Object.entries(geometry.debugLayers)) {
    const group = new THREE.Group();
    group.visible = false;
    
    if (layer.points && layer.points.positions) {
      const positions = layer.points.positions;
      const colors = layer.points.colors;
      const nPoints = positions.length / 3;
      
      const pointsGeo = new THREE.BufferGeometry();
      pointsGeo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
      if (colors && colors.length > 0) {
        pointsGeo.setAttribute('color', new THREE.BufferAttribute(colors, 4));
      }
      
      const pointsMat = new THREE.PointsMaterial({
        size: 0.05,
        vertexColors: colors && colors.length > 0,
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
      
      const lineGeo = new THREE.BufferGeometry();
      lineGeo.setAttribute('position', new THREE.BufferAttribute(linePositions, 3));
      
      const lineColors = colors ? new Float32Array((starts.length / 3) * 2 * 4) : null;
      if (colors) {
        for (let i = 0; i < starts.length / 3; i++) {
          lineColors[i * 8] = colors[i * 4];
          lineColors[i * 8 + 1] = colors[i * 4 + 1];
          lineColors[i * 8 + 2] = colors[i * 4 + 2];
          lineColors[i * 8 + 3] = colors[i * 4 + 3];
          lineColors[i * 8 + 4] = colors[i * 4];
          lineColors[i * 8 + 5] = colors[i * 4 + 1];
          lineColors[i * 8 + 6] = colors[i * 4 + 2];
          lineColors[i * 8 + 7] = colors[i * 4 + 3];
        }
        lineGeo.setAttribute('color', new THREE.BufferAttribute(lineColors, 4));
      }
      
      const lineMat = new THREE.LineBasicMaterial({
        vertexColors: colors !== null,
      });
      const linesObj = new THREE.LineSegments(lineGeo, lineMat);
      group.add(linesObj);
    }
    
    scene.add(group);
    debugGroups[layerName] = group;
  }
  
  if (Object.keys(debugGroups).length > 0) {
    const controlsDiv = document.createElement('div');
    controlsDiv.style.cssText = 'position:absolute;top:10px;right:10px;background:rgba(0,0,0,0.7);padding:10px;color:white;font-family:sans-serif;';
    
    for (const layerName of Object.keys(debugGroups)) {
      const label = document.createElement('label');
      label.style.cssText = 'display:block;margin:5px 0;';
      const checkbox = document.createElement('input');
      checkbox.type = 'checkbox';
      checkbox.onchange = (e) => {
        debugGroups[layerName].visible = e.target.checked;
      };
      label.appendChild(checkbox);
      label.appendChild(document.createTextNode(' ' + layerName));
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

  window.addEventListener('resize', () => {
    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(container.clientWidth, container.clientHeight);
  });
}

export { init as initTreeViewer };