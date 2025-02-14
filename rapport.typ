#import "@preview/fletcher:0.5.4" as fletcher: diagram, node, edge
#import "@preview/note-me:0.3.0": *


#align(center, [
#set text(size: 20pt)
    = Final semester project IGR
    == Interactive Invigoration: Volumetric Modeling of Trees with Strands

    PERONNET Antonin

])

In this project, I tried to reproduce and extend the work done in:

Bosheng Li, Nikolas A. Schwarz, Wojtek Pałubicki, Sören Pirk, and Bedrich
Benes. 2018. Interactive Invigoration: Volumetric Modeling of Trees with
Strands. In . ACM, New York, NY, USA, 13 pages

= Related work

The goal of this project is to generate the entire plant geometry (the mesh including colors) from no data.

In most plant mesh generation pipelines, there are 3 main steps:
- creating the skeleton (macro-level representation)
- computing the details (micro-level)
- generating the geometry

== Skeleton generation

The first approach that was used to generate skeleton of trees was fractal modeling. The basic idea is that a branch contains other sub-branches, which look a lot like the main branch, but scaled down. The idea of L-systems #footnote[P Prusinkiewicz. 1986. Graphical Applications of L-systems] generalizes this idea. It consists of rules (similar to formal grammars) to generate segments recursively.

Nowadays, biological models of plants are used to simulate the growth. A very sophisticated one has been made to include wood fracture for example. #footnote[- T. Hädrich, J. Scheffczyk, W. Palubicki, S. Pirk, and D. L. Michels. 2020. Interactive Wood Fracture. In Eurographics/ ACM SIGGRAPH Symposium on Computer Animation Posters. The Eurographics Association.]

The paper I will use depends on the work from an older paper, that simulate the growth of the roots of a tree#footnote[Li, J. Klein, D. L. Michels, B. Benes, S. Pirk, and W. Pałubicki. 2023. Rhizomorph: The Coordinated Function of Shoots and Roots. ACM Trans. Graph. 42, 4, Article 59 (jul 2023), 16 pages.]

== Micro level

Genereating the micro-level gemoetry of a tree is very resource-intensive. To have a perfect simulation, one would need to simulate every cell, which can interact in complicated ways (such as the patterns in brak and dendrites).

As for many complex system, an succesful approach is particule-based#footnote[X. Chen, B. Neubert, Y.-Q. Xu, O. Deussen, and S. B. Kang. 2008. Sketch-Based Tree Modeling Using Markov Random Field. ACM TOG 27, 5, Article 109 (Dec. 2008).]. Particles represent unitary volumes in the tree, in he same way as they can represent unit volumes of water in fluid simulations.

The main innovation of the paper is to use strands, sort of threads that interpolate between particles.


== Geometry

Very often, the geometry part uses very well established techniques in the field of computer graphics. The most used techniques are:
- delaunay triangulation
- convex hull algorithms
- cubic splines
- mesh techniques such as smoothed normal calculation, mesh smoothing.


= Implementation

From a very high level, the pipeline is the following:



#align(center, diagram(
    node((0, 0), name:<start>),
    node((0, 1), [1) Plant Graph simulation (growing)], name: <p1>, stroke: 1pt),
    node((0, 2), [2) Particle flow and interpolation], name: <p2>, stroke: 1pt),
    node((0, 3), [3) Contours calculation], name: <p3>, stroke: 1pt),
    node((0, 4), [4) Meshing and branch fusion], name: <p4>, stroke: 1pt),
    node((0, 5), name: <end>),
    edge(<start>, <p1>, [tree caracteristics], "->", label-side: left),
    edge(<p1>, <p2>, [skeleton], "->", label-side: left),
    edge(<p2>, <p3>, [strands], "->", label-side: left),
    edge(<p3>, <p4>, [contours (branch sections)], "->", label-side: left),
    edge(<p4>, <end>, [3d mesh], "->", label-side: left),
))


I focused on the particle and meshing parts (2, 3, 4) and spend very little time on the skeleton calculation (1). Thus, I will present the steps 2), 3), 4) and 1) at the very end.


For this project, I used #link("https://bevyengine")[bevy], a rust cross-platfrom game design library. I decided to use it because of the following features:
- a `gizmos` tool that allow to draw 3d points and lines, vey useful for debuging
- a very good keyboard and mouse support
- a dedicated way to store assets, resources and properties
- the ability to compile to `wasm` and `webgl`


My work is available #link("https://github.com/rambip/plant-mesh")[on github] under the MIT license and can be tested on a web interface here: #link("https://rambip.github.io/plant-mesh/")[https://rambip.github.io/plant-mesh/]


#pagebreak()

= Step 0: Defining the skeleton

The skeleton is a simple tree structure with some properties stored in each node:

```rust
pub struct PlantNode {
    children: Vec<PlantNode>,
    props: PlantNodeProps,
}

pub struct PlantNodeProps {
    // a 3d vector
    pub position: Vec3,
    pub radius: f32,
    // a normalized 3d vector
    pub orientation: Vec3,
}
```


#note[
An optimization I did later to make the code easier to understand is to represent the orientation each node by a quaternion. This way, the transformation to go from an absolute plane to a branch and back are just multiplication by the quaternion.
]

I decided to restrict the tree to a binary tree (each node has 0, 1 or 2 children) like in the paper. It makes the particle projection and the branch fusions a lot easier.

I took the convention that the first child is the *main* branch, and the other is the *secondary* branch (if it exists).

This pointer-based representation of the tree is ideal for the skeleton growing phase, because it is recursive in nature.
For the latter stages of the pipeline, I decided to use an array based representation.


#{
set text(size: 10pt)
set align(center)

grid(columns: (1fr),
inset: 10pt,
stroke: black,
[=== Pointer-based],
diagram(
  spacing: (15pt, 15pt),
  node-stroke: 1pt,
  edge-stroke: 1pt+blue,
  mark-scale: 60%,
  
  // Root node and its properties
  node((0,0), $A$, shape:circle, name: <a>),
  node((-1,0), $v_a$, shape:rect, name: <va>),
  edge((0,0), (-1,0), "..|>"),
  
  // Left child (B)
  node((-1,1), $B$, shape:circle, name: <b>),
  edge(<a>, <b>, "-|>"),
  node((-2,1), $v_b$, shape:rect, name: <vb>),
  edge(<b>, <vb>, "..|>"),
  
  // Right child (C)
  node((1,1), $C$, shape:circle, name: <c>),
  edge(<a>, <c>, "-|>"),
  node((0,1), $v_c$, shape:rect, name: <vc>),
  edge(<c>, <vc>, "..|>"),
),

[],
[=== Array based],
diagram(
    spacing: (12pt, 10pt),
    node-stroke: 1pt,
    edge-stroke: 1pt+blue,
    mark-scale: 60%,

    // Node A [index 0]
    node((0,0), "parent", shape: rect, name: <pa>),
    node((1,0), "children", shape: rect, name: <ca>),
    node((2,0), $v_a$, shape: rect, name: <va>),
    node((1, -1), "A", enclose: ((1, -1), <pa>, <ca>, <va>), name: <a>),


    node((4,1), "parent", shape: rect, name: <pb>),
    node((5,1), "children", shape: rect, name: <cb>),
    node((6,1), $v_b$, shape: rect, name: <vb>),
    node((5, 0), "B", enclose: ((5, 0), <pb>, <cb>, <vb>), name: <b>),


    node((7,2), "parent", shape: rect, name: <pc>),
    node((8,2), "children", shape: rect, name: <cc>),
    node((9,2), $v_c$, shape: rect, name: <vc>),
    node((8, 1), "C", enclose: ((8, 1), <pc>, <cc>, <vc>), name: <c>),

    edge(<ca>, "d,r,r", "->"),
    edge(<ca>, "d,d,r,r,r,r,r", "->"),

    edge(<pb>, "u,u,l", "->"),
    edge(<pc>, "u,u,u,l,l,l,l,l", "->"),
))
}


= Step 1

As I said, I will detail the plant growing strategy later.

For all the tests and illustrations you will see, I used a manually crafted tree with 7 nodes:


#align(center, image("images/manually_crafted_tree.png", width: 30%))


= Step 2: Particle flow and interpolation

*code*: `src/meshing/particles.rs`

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

= Step 3: Contours calculation

*code*: `src/meshing/mod.rs`

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


- parametrization: if 2 branches are not the same size, it's very hard to merge them.
We must find the right step size for both

#image("images/tree_branch_non_uniform.svg")
#image("images/mesh_non_uniform_triangulation.png")
#image("images/mesh_uniform_triangulation.png")


but it is still an issue for branch fusion. To have a really regular mesh, the solution would be to have a good reparametrization of the spline. I decided not to go this way but to approximate it branch section by branch section. It works until there is a large ratio between the lengths of the 2 branches.


setting the right parameter for the convex hull / delaunay treshold is hard.

= Step 4: Meshing and branch fusion


- non-uniform sampling: too much triangles in some parts
    - adapt dt according to radius and length of branch
        


- normals: compute using neighbourgs. strategy: average of 4 neighbourgs


- a lot of edge cases:
    - when not enough points to merge
    - when two children branches do not split
    - when the branches split just before or just after a node in the tree


- branch fusion


#let branches() = {
        let n = 20
        let points_m = range(2, n)
        .map(i => {
                let θ = i*360deg/n
                let x = calc.cos(θ)-1.5
                let y = calc.sin(θ)
                (x, y)
                })

        let points_s = points_m.map(pos => {
            (- pos.at(0), pos.at(1))
        })

        let points_aam = range(2, n)
        .map(i => {
                let θ = 45deg + i*275deg/n
                let x = 3*calc.cos(θ)-1.5
                let y = 3*calc.sin(θ)
                (x, y)
                })

        let points_aas = points_aam
        .map(pos => (- pos.at(0), pos.at(1)))

        let points_am = range(n - 2)
            .map(i => (
                0.5*(points_m.at(i).at(0) + points_aam.at(i).at(0)),
                0.5*(points_m.at(i).at(1) + points_aam.at(i).at(1)),
            ))

        let points_as = range(n - 2)
            .map(i => (
                0.5*(points_s.at(i).at(0) + points_aas.at(i).at(0)),
                0.5*(points_s.at(i).at(1) + points_aas.at(i).at(1)),
            ))

        let n = 19
        for i in range(n - 2) {
        let pos = points_m.at(i)
            node(pos, [$times$], name: "s"+str(i))
        }

        for i in range(n - 2) {
        let pos = points_s.at(i)
            node(pos, [$times$], name: "m"+str(i))
        }
        for i in range(n - 2) {
        let pos = points_am.at(i)
            node(pos, [$times$], name: "am"+str(i))
        }
        for i in range(n - 2) {
        let pos = points_as.at(i)
            node(pos, [$times$], name: "as"+str(i))
        }

        for i in range(n - 3) {
        // refer to nodes by label, e.g., <1>
            edge(label("s"+str(i)), label("s"+str(i+1)), "-")
            edge(label("m"+str(i)), label("m"+str(i+1)), "-")
            edge(label("am"+str(i)), label("am"+str(i+1)), "-")
            edge(label("as"+str(i)), label("as"+str(i+1)), "-")
        }
}


#diagram(branches())

#diagram(branches(),
    node((0, 0), box({set text(stroke: blue); $times$}))
)
