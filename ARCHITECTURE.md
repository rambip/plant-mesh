# Architecture - plant-mesh

## Overview
This project implements a procedural plant generation pipeline based on the paper *"Interactive Invigoration: Volumetric Modeling of Trees with Strands"*. It generates 3D meshes from scratch through a multi-step process.

## Generation Pipeline
The pipeline follows a bottom-up/top-down hybrid approach:

1.  **Skeleton Simulation (Growing):** Generates a binary tree structure (`PlantNode`) where each node contains position, radius, and orientation (stored as a `Quat`).
2.  **Particle Flow & Interpolation:** 
    - Particles are generated at the leaves and projected down towards the root.
    - A particle-based simulation (with repulsion forces) ensures uniform distribution within the branch volume.
    - **Strands:** Catmull-Rom splines interpolate particle positions across nodes to create continuous "fibers".
3.  **Contours Calculation:**
    - Detects branch splits by measuring distances between particle clouds.
    - Computes 2D boundaries of particle sections using a custom $O(n \log n)$ convex hull algorithm (Graham scan variant).
4.  **Meshing & Branch Fusion:**
    - Joins consecutive contours using a greedy triangulation algorithm that minimizes triangle perimeters.
    - **Branch Fusion:** Handles the complex geometry where two child branches merge into a single parent trunk by finding nearest points and creating a seamless transition.

## Core Libraries & Frameworks

### Bevy
The project uses the [Bevy](https://bevyengine.org/) engine for rendering, asset management, and windowing.

### bevy_simple_graphics
A lightweight, custom rendering plugin located in `bevy_simple_graphics/`. It provides a simplified 3D rendering pipeline for Bevy without the full overhead of `bevy_pbr`. It handles:
- Custom mesh vertex layouts (Position, Normal, Color).
- Transparent 3d render phase integration.
- Efficient mesh instance extraction and queuing.

### bevy_gizmos
A specialized fork of Bevy's gizmo system located in `bevy_gizmos/`. It has been modified to operate without the standard Bevy ECS overhead where possible, providing an immediate-mode drawing API for:
- **Visual Debugging:** Drawing skeletons, particle clouds, and contours in real-time.
- **Line Rendering:** Efficiently rendering large numbers of lines and joints (miter, round, bevel) using specialized WGSL shaders (`lines.wgsl`, `line_joints.wgsl`).

## Data Structures
- **PlantNode:** Recursive pointer-based structure for growth.
- **Array-based Skeleton:** Used in later pipeline stages for performance and easier traversal.
- **Quadtree:** Used in particle simulations to optimize $O(n^2)$ neighbor searches.
