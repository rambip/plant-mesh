#align(center, [
#set text(size: 20pt)

])

= Related work


3 mains steps in a lot of plant mesh generation tasks:
- creating the skeleton (macro-level)
- generating details (micro-level)
- generating the geometry

== Skeleton generation

- L-system
- biology


- P Prusinkiewicz. 1986. Graphical Applications of L-systems

- T. Hädrich, J. Scheffczyk, W. Palubicki, S. Pirk, and D. L. Michels. 2020. Interactive Wood Fracture. In Eurographics/ ACM SIGGRAPH Symposium on Computer Animation Posters. The Eurographics Association.

== Micro level

- particles
X. Chen, B. Neubert, Y.-Q. Xu, O. Deussen, and S. B. Kang. 2008. Sketch-Based Tree Modeling Using Markov Random Field. ACM TOG 27, 5, Article 109 (Dec. 2008).

idea of strand: TODO

== Geometry

Established techniques:
- delaunay
- splines
- graph traversal




== particles


== meshing

- idea of strands

= Implementation choices

- rust
    - bevy


= Algorithmic choices

- overview of the pipeline
    1) defining the plant graph attributes
    2) simulating the growth of the tree
    3) particle flow
    4) interpolation
    5) convex hull and points selection
    6) meshing

== plant-graph simulation

== branch meshing

we go from bottom to top(s)



- delaunay is not ideal
    - complex and costly
    - need to recompute the convex hull after
    - what we need is an algorithm for convex hull with added constraints

- convex hull is better
    - problem 1: the points are not on a plane
        - solution: project the points on a plane. PCA is possible but inefficient
          to minimize the vertical distance, we can use the pseudo-inverse. works if points not almost vertical.
          At the end, the choice was to 

    - problem 2: we don't want too long vertices
        - solution: modified graham scan to ignore vertices when length too long

    - problem 3: multiple connected components
        

- triangulation
    - max vs sum

- non-uniform sampling: too much triangles in some parts
    - adapt dt according to radius and length of branch

- normals: compute using neighbourgs. strategy: average of 4 neighbourgs

- particle positioning
    - 1) projection
    - 2) offsets
    - 3) particle based dynamics


- recurring problem: strands "demangling"
    - the strands are based on particle positions, and they can easily mix together
        - particles inside move outside
        - particles a and b are swaped

    - problem in the meshing algorithm: vertices can appear and disappear.
        - rectangles when 2 particles are present before and after
        - triangle when two particles merge
        - triangle when one particle merge into two
        - case when one particle is replaced or 2 particles swap: we hope it will not appear ?


- parametrization: if 2 branches are not the same size, it's very hard to merge them.
We must find the right step size for both
