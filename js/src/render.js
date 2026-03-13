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

  const cameraSettings = {
    orbitDistance: 2,
    orbitAngle: 0,
    sensibility: 1,
    tilt: 0,
    z: 0,
    automaticMode: true,
  };

  const xAxis = new THREE.Vector3(1, 0, 0);
  const zAxis = new THREE.Vector3(0, 0, 1);
  const target = new THREE.Vector3();
  const rotationZ = new THREE.Quaternion();
  const rotationX = new THREE.Quaternion();
  const cameraRotation = new THREE.Quaternion();
  const cameraForward = new THREE.Vector3();

  function updateCamera(timeSeconds) {
    const angle = cameraSettings.automaticMode
      ? 0.5 * timeSeconds
      : cameraSettings.orbitAngle;

    rotationZ.setFromAxisAngle(zAxis, angle);
    rotationX.setFromAxisAngle(xAxis, 0.5 * Math.PI + cameraSettings.tilt);
    cameraRotation.copy(rotationZ).multiply(rotationX);

    camera.quaternion.copy(cameraRotation);

    target.set(0, 0, cameraSettings.z);
    cameraForward.set(0, 0, -1).applyQuaternion(cameraRotation);
    camera.position.copy(target).addScaledVector(cameraForward, -cameraSettings.orbitDistance);
  }

  const renderer = new THREE.WebGLRenderer();
  renderer.setSize(container.clientWidth, container.clientHeight);
  container.appendChild(renderer.domElement);

  let pointerDragging = false;
  let lastPointerX = 0;
  let lastPointerY = 0;

  renderer.domElement.addEventListener('pointerdown', (event) => {
    if (event.button !== 0) {
      return;
    }
    pointerDragging = true;
    lastPointerX = event.clientX;
    lastPointerY = event.clientY;
    cameraSettings.automaticMode = false;
    renderer.domElement.setPointerCapture(event.pointerId);
  });

  renderer.domElement.addEventListener('pointermove', (event) => {
    if (!pointerDragging) {
      return;
    }

    const dx = event.clientX - lastPointerX;
    const dy = event.clientY - lastPointerY;
    lastPointerX = event.clientX;
    lastPointerY = event.clientY;

    cameraSettings.z += 0.0005 * dy * cameraSettings.orbitDistance;
    cameraSettings.orbitAngle -=
      0.0005 * dx * cameraSettings.orbitDistance * cameraSettings.sensibility;
  });

  renderer.domElement.addEventListener('pointerup', (event) => {
    if (event.button !== 0) {
      return;
    }
    pointerDragging = false;
    renderer.domElement.releasePointerCapture(event.pointerId);
  });

  renderer.domElement.addEventListener('pointercancel', (event) => {
    pointerDragging = false;
    renderer.domElement.releasePointerCapture(event.pointerId);
  });

  renderer.domElement.addEventListener('wheel', (event) => {
    event.preventDefault();
    cameraSettings.orbitDistance += 0.001 * event.deltaY;
    cameraSettings.orbitDistance = Math.max(0.2, cameraSettings.orbitDistance);
  }, { passive: false });

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
      const nLines = starts.length / 3;
      
      const linePositions = new Float32Array((nLines * 2) * 3);
      for (let i = 0; i < nLines; i++) {
        const src = i * 3;
        const dst = i * 6;
        linePositions[dst] = starts[src];
        linePositions[dst + 1] = starts[src + 1];
        linePositions[dst + 2] = starts[src + 2];
        linePositions[dst + 3] = ends[src];
        linePositions[dst + 4] = ends[src + 1];
        linePositions[dst + 5] = ends[src + 2];
      }
      
      const lineGeo = new THREE.BufferGeometry();
      lineGeo.setAttribute('position', new THREE.BufferAttribute(linePositions, 3));
      
      const hasLineColors = colors && colors.length > 0;
      const lineColors = hasLineColors ? new Float32Array(nLines * 2 * 4) : null;
      if (colors) {
        for (let i = 0; i < nLines; i++) {
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
        vertexColors: hasLineColors,
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
    updateCamera(performance.now() / 1000);
    renderer.render(scene, camera);
  }

  animate();

  window.addEventListener('resize', () => {
    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(container.clientWidth, container.clientHeight);
  });
}

window.initTreeViewer = init;
