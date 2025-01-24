#align(center, [
#set text(size: 20pt)

])

= Related work

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

- delaunay is not ideal
    - complex and costly
    - need to recompute the convex hull after
    - what we need is an algorithm for convex hull with added constraints

- convex hull is better
    - problem 1: the points are not on a plane
        - solution: project the points on a plane. PCA is possible but inefficient
          to minimize the vertical distance, we can use the pseudo-inverse. works if points not almost vertical.
          We can just drop the z coordinates

    - problem 2: we don't want too long vertices
        - solution: modified graham scan to ignore vertices when length too long


- triangulation
    - max vs sum

- non-uniform sampling: too much triangles in some parts
    - adapt dt according to radius and length of branch

- normals: compute laplacian using neighbourgs. strategy: average of 4 neighbourgs
