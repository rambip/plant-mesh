#align(center, [
#set text(size: 20pt)

])

= Related work

= Implementation choices

- rust
    - bevy


= Algorithmic choices

== plant-graph simulation

== branch meshing

- delaunay is not ideal
    - complex and costly
    - need to recompute the convex hull after
    - what we need is an algorithm for convex hull with added constraints

- convex hull is better
    - problem 1: the points are not on a plane
        - solution: project the points on a plane. PCA is possible but inefficient
          to minimize the vertical distance, we can use the pseudo-inverse. works if points not almost vertical

    - problem 2: we don't want too long vertices
        - solution: modified graham scan to ignore vertices when length too long
