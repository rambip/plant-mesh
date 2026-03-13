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

  function showDecodeError(error, inputData) {
    const previous = document.getElementById('treemesh-decode-error');
    if (previous) {
      previous.remove();
    }

    const panel = document.createElement('div');
    panel.id = 'treemesh-decode-error';
    panel.style.cssText = [
      'position:fixed',
      'left:12px',
      'right:12px',
      'top:12px',
      'z-index:99999',
      'max-height:45vh',
      'overflow:auto',
      'background:#2b0f12',
      'border:1px solid #ff5f67',
      'color:#ffd7d9',
      'padding:12px',
      'font:12px/1.45 monospace',
      'white-space:pre-wrap',
      'border-radius:8px',
      'box-shadow:0 6px 24px rgba(0,0,0,0.35)',
    ].join(';');

    const outputKeys = inputData && inputData.outputs
      ? Object.keys(inputData.outputs)
      : [];
    const bufferCount = inputData && inputData.buffers
      ? Object.keys(inputData.buffers).length
      : 0;

    const errorName = error && error.name ? error.name : 'Error';
    const errorMessage = error && error.message ? error.message : String(error);

    const details = [
      '[TreeMesh decode failed]',
      '',
      `${errorName}: ${errorMessage}`,
      '',
      String(error && error.stack ? error.stack : error),
      '',
      `treemesh=${inputData && inputData.treemesh ? inputData.treemesh : 'unknown'}`,
      `buffers=${bufferCount}`,
      `outputs=${outputKeys.join(', ')}`,
    ];

    panel.textContent = details.join('\n');
    document.body.appendChild(panel);
  }

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
    targetX: 0,
    targetY: 0,
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

    target.set(cameraSettings.targetX, cameraSettings.targetY, cameraSettings.z);
    cameraForward.set(0, 0, -1).applyQuaternion(cameraRotation);
    camera.position.copy(target).addScaledVector(cameraForward, -cameraSettings.orbitDistance);
  }

  function expandBoundsFromFlat(bounds, flatPositions) {
    if (!flatPositions || flatPositions.length < 3) {
      return;
    }
    for (let i = 0; i + 2 < flatPositions.length; i += 3) {
      bounds.expandByPoint(
        new THREE.Vector3(flatPositions[i], flatPositions[i + 1], flatPositions[i + 2])
      );
    }
  }

  function fitInitialCamera(bounds) {
    if (bounds.isEmpty()) {
      return;
    }

    const size = new THREE.Vector3();
    const center = new THREE.Vector3();
    bounds.getSize(size);
    bounds.getCenter(center);

    const radius = Math.max(size.length() * 0.5, 0.001);
    const vFov = THREE.MathUtils.degToRad(camera.fov);
    const hFov = 2 * Math.atan(Math.tan(vFov * 0.5) * camera.aspect);
    const minFov = Math.max(Math.min(vFov, hFov), 0.01);
    const distance = (radius / Math.sin(minFov * 0.5)) * 1.15;

    cameraSettings.orbitDistance = Math.max(distance, 0.2);
    cameraSettings.targetX = center.x;
    cameraSettings.targetY = center.y;
    cameraSettings.z = center.z;
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
    const distance = Math.max(cameraSettings.orbitDistance, 0.2);
    cameraSettings.orbitDistance += 0.0011 * event.deltaY * Math.sqrt(distance);
    cameraSettings.orbitDistance = Math.max(0.2, cameraSettings.orbitDistance);
  }, { passive: false });

  const ambientLight = new THREE.AmbientLight(0x404040);
  scene.add(ambientLight);

  const directionalLight = new THREE.DirectionalLight(0xffffff, 1);
  directionalLight.position.set(1, 1, 1);
  scene.add(directionalLight);

  const treeMeshDecoder = new TreeMeshDecoder();
  let geometry;
  try {
    geometry = treeMeshDecoder.decode(geometryData);
  } catch (error) {
    console.error('TreeMesh decode failed', error, geometryData);
    showDecodeError(error, geometryData);
    return;
  }

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

  const initialBounds = new THREE.Box3();
  const meshPositions = threeGeometry.getAttribute('position');
  if (meshPositions && meshPositions.count > 0) {
    initialBounds.setFromBufferAttribute(meshPositions);
  }

  const debugGroups = {};
  for (const [layerName, layer] of Object.entries(geometry.debugLayers)) {
    const group = new THREE.Group();
    group.visible = layer.show === true;
    
    if (layer.points && layer.points.positions) {
      const positions = layer.points.positions;
      const colors = layer.points.colors;
      const nPoints = positions.length / 3;
      expandBoundsFromFlat(initialBounds, positions);
      
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
      expandBoundsFromFlat(initialBounds, starts);
      expandBoundsFromFlat(initialBounds, ends);
      
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

  fitInitialCamera(initialBounds);
  
  if (Object.keys(debugGroups).length > 0) {
    const controlsDiv = document.createElement('div');
    controlsDiv.style.cssText = 'position:absolute;top:10px;right:10px;background:rgba(0,0,0,0.7);padding:10px;color:white;font-family:sans-serif;';
    
    for (const layerName of Object.keys(debugGroups)) {
      const label = document.createElement('label');
      label.style.cssText = 'display:block;margin:5px 0;';
      const checkbox = document.createElement('input');
      checkbox.type = 'checkbox';
      checkbox.checked = debugGroups[layerName].visible;
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
