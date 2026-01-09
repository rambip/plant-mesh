#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge
#import "@preview/note-me:0.5.0": *

#show link: underline
#show image: it => {
    set align(center)
    it
}
#show heading: it => {
    v(15pt)
    it
}

#show raw: it => {
    set text(8pt)
    it
}

#let scale_center(x) = {
    set align(center)
    scale(70%, reflow: true, x)
}


#align(center, [
#set text(size: 30pt)
    = IGR \ Final semester project
    == Procedural geometry generation \ for plants and trees

    #v(50pt)

#set text(size: 20pt)
    Based on the paper:

#link("https://dl.acm.org/doi/10.1145/3658206")[_Bosheng Li, Nikolas A. Schwarz, Wojtek Pałubicki, Sören Pirk, and Bedrich Benes. 2018. Interactive Invigoration: Volumetric Modeling of Trees with Strands. In . ACM, New York, NY, USA, 13 pages_]

    #v(100pt)
    PERONNET Antonin

])

#pagebreak()

#set heading(numbering:"1.1")

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


#pagebreak()

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

#v(30pt)


For this project, I used #link("https://bevyengine")[bevy], a rust cross-platfrom game design library. I decided to use it because of the following features:
- a `gizmos` tool that allow to draw 3d points and lines, vey useful for debuging
- a very good keyboard and mouse support
- a dedicated way to store assets, resources and properties
- the ability to compile to `wasm` and `webgl`


My work is available #link("https://github.com/rambip/plant-mesh")[on github] under the MIT license and can be tested on a web interface here: https://rambip.github.io/plant-mesh/

In the online demo, a tree is generated automatically every 5s. The user can visualize the different steps (see @gallery). I did not implement the ability to change the parameters, but it can easily be added#footnote[In the desktop version, the parameters are loaded from a config file.]


#pagebreak()

== Step 0: Defining the skeleton

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
  spacing: (15pt, 10pt),
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

grid.cell(fill: black, inset: -2pt,[]),
[=== Array based],
diagram(
    spacing: (12pt, 8pt),
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

== Step 1

As I said, I will detail the plant growing strategy later.

For all the tests and illustrations you will see, I used a manually crafted tree with 7 nodes:


#align(center, image("images/manually_crafted_tree.png", width: 30%))


== Step 2: Particle flow and interpolation

*code*: `src/meshing/particles.rs`


The pseudo-code for this step is the following:

```
compute_paticles(node):
    if node is a leaf:
        create a set of particles for this leaf
    else:
        for each child:
            compute_particles(child)
        project particles onto parent
        do particle based collision detection
        return the set of particles
```

As you see, they are 2 main parts in this algorithm:
- projecting the particles
- doing the particle simulation

#note[
The paper proposed to use arbitrary shaped for the section, I only implemented circles. It may improve slightly the diversity of possible trees, but most trees have circular branches.
]

#pagebreak()

=== Particle projection

When there is only one child, the projection operation is pretty easy. For 2 children, it is more complicated.
The paper did not describe how they defined the projection entirely, so I had to interpret what they did.

We want to compute the positions of the 2 clouds of points (orange and red) on the root parent plane:

#image("images/projection_1.svg", height: 200pt)

The first step is to use the position of the parent and the 2 children to compute the important directions:
- from child1 to parent
- from child2 to parent
- from child1 to child2, in the parent plane

#image("images/projection_2.svg", height: 200pt)


Then, the a point on the cloud of child1 will be projected on the plane parallel to $D_1$ and moved with an offset in the direction of -$D_{12}$. Similarly, a point in the cloud of child2 will be projected parralel to $D_2$ and offset in the direction $D_12$.

The offsets are calculated such that the two set of points are adjacent in the parent plane.


#pagebreak()

After projection and scaling, the particle trajectories look like:
#image("images/strands_without_collision.png", height: 150pt)

In order to give the realistic result of trunks, we need to make the particles interact with each other.

=== Particle simulation

Again, the paper was not very precise regarding the type of particle simulation they used.

I tried to use an electrostatic-like repulsion force with an euler integratino scheme, with collisions on the border.

The force is defined as:

$
arrow(F_(a b)) = -k/(A B)^2 arrow(u_(a b))
$

My simulation has normalized parameters like:
- the repulsion force between particles
- the repulsion force with the wall
- the time step

For numeric stability, I had to update the parameters depending on the radius and the number of particles. For example, the repulsion is proportionnal to the radius squared.

Here is the result on an example:

#let images = range(1,8).map(i => image("images/particles_visu_a"+str(i)+".png"))
#grid(columns: 8, ..images)

It works, but as you can see there is a problem: there are too much particles on the contours and not enough inside. This is due to the fact that the initial speed of the particles is too high.

To solve this problem, I introduced a maximum velocity. I also multiplied the repulsion constant by $1/sqrt(N)$ where $N$ is the number of particles. This way, the average speed is the same no matter the number of particles. After these adjustments, I got a more uniform particle density:

#let images = range(1,8).map(i => image("images/particles_visu_b"+str(i)+".png"))
#grid(columns: 8, ..images)

#pagebreak()

It made a huge difference for the phases after that:

#grid(columns: (1fr, 1fr),
figure(image("images/base_strands_non_uniform.png", height: 140pt),
caption: "particles without max veolcity"),

figure(image("images/base_strands_uniform.png", height: 140pt),
caption: "particles with max velocity")
)

#v(50pt)


The particle simulation step is also the most resource-intensive part, because of the $O(n^2)$ complexity to compute the interactions between particles. (see @performance). This is why I tried to optimize it, by constructing the neighbourgs less frequently.

Building the neighbourgs depending on some interaction radius every 5 steps worked quite well. To improve the performance further, I used a quadtree implementation for efficient retrieval. See @performance for more details.


#{
set align(horizon)
figure(image("images/numerical_instability.png", height: 200pt),
caption: "if dt is not small enough, it lead to numerical instability")
}



#pagebreak()

=== Strands

Once we know the particle positions on each node, we can interpolate them using splines.

The paper proposed *B-splines*, but as this type of spline does not interploate the points, \ I chose the #link("https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline")[Catmull-Rom spline].

A very frequent problem in this project is strands "mangling". After the projection operations and the simulation, strands can:
- go inside or outside the branch
- turn one around another


#grid(columns: (1fr, 1fr),
figure(image("images/strands_zigzag.png", height: 120pt), caption: "strands can zigzag"),
figure(image("images/strands_mangling.png", height: 120pt), caption: "strands can mangle")
)

This can cause a number of problems for the next steps, like contours calculation and triangulation. The paper proposed some techniques like swapping the points positions when it is minimizes the distance, but I did not implement it.



#pagebreak()
== Step 3: Contours calculation

*code*: `src/meshing/mod.rs`


The pseudo-code for this part is the following:

```
compute_contours(node):
    if node is a leaf:
        stop
    else if one children:
        for each position from parent to children:
            compute the positions of the point cloud
            select the boundary of the point cloud
        compute_contours(children)

    else if 2 children:
        compute the position when the branch splits
        for each position before the split:
            compute the positions of the point cloud
            select the boundary of the point cloud

        for each position between the split and first children:
            compute the positions of the point cloud
            select the boundary of the point cloud
        compute_contours(first children)

        for each position between the split and second children:
            compute the positions of the point cloud
            select the boundary of the point cloud
        compute_contours(second children)
```

Note that this algorithm is bottom-up, contrary to the previous one which is top-down.


=== Detecting if branch is spliting

To detect if the branch is splitting, I compute the minimum distance between the points in the first cloud and the points in the second cloud. If this distance is bigger than a treshold, that means the branch is splitting.#footnote[I realized I can do it by dichotomy]

Before the detected branch split, I use the union of the 2 points clouds to compute the contour. After the split, I compute the contours separately.

#image("images/contours_with_branch_split.png", height: 100pt)

#warning[
In some cases, the branch never splits. It is an edge case I did not imagine at first, but it can happen with lot's of branches. In this case, we stop at the shortest branch.
]



#pagebreak()

=== Points selection

This was probably the hardest part of the entire project.

The paper proposed delaynay triangulation to compute the boundary of the points, but early on I thought this was not the best idea. Indeed:
- it is very costly, $O(n^2)$ in the worst case even more complex in 3d.
- it builds the entire triangulation of the section, but uses only the contour
- it does not give the contour directly, one must do a graph traversal afterwards.
- the paper removed the edges from the triangulation when they are longer than a treshold, but don't indicate how to set this treshold.

#note[
    To convert from 3d points (splines in space) to 2d, I considered using a minimization algorithm, to find the 2d plane that miminimize the error from the reality. After, I figured out I can simply project the points using the orientation vector of the nearest parent.

]

So my first intuition was to use a convex hull algorithm like #link("https://en.wikipedia.org/wiki/Graham_scan")[Graham scan]. It worked, but the issue is that it does not follow the contour when 2 branch merge (see below). I had to adapt the algorithm.


After a lot of work and debuging, I ended up with this algorithm:

$
 "select_contours" & ((P_i)): \
    & O <- 1/n sum P_i \
    & i_0 <- "argmin"_i (e_y dot P_i) \
    & L <- [i_0] \
    & sigma <- "sort points depending on" P_i |-> P_(i_0) O P_i  \
    & "for" i in sigma : \
    & | | "add i to" L \
    & | | "remove the points in L depending on angle" \
    & "return L"
$

#v(20pt)

To illustrate:
#grid(columns: (1fr, 1fr),
figure(image("images/convex_hull_1.svg", height: 150pt), caption: "First, sort the points around the center"),
figure(image("images/convex_hull_2.svg", height: 150pt), caption: [Each time one point is added, \ remove the previous one depending on angle])
)


This algorithm is $O(n log n)$ in the worst case, because of the sort.

By varying the angle treshold, we can select different boundaries:

#grid(columns: (1fr, 1fr, 1fr),
image("images/convex_hull_no_treshold.png"),
image("images/convex_hull_treshold_1.png"),
image("images/convex_hull_treshold_2.png"),
)

But setting the right treshold for the convex hull is not trivial.
- very high treshold: the contour does not follow the shape of the strands, it is especially problematic for meshing (see @meshing)
- very low treshold: too much noise and not enough regularity in the contour that is seleted.

#warning[
When computing the mesh of a large number of points, I noticed some strange behaviours:

#image("images/branch_join_convex_bug.png", height: 100pt)

I was able to reproduce it and to understand where the issue was coming from:

#figure(image("images/convex_hull_instability.png", height: 130pt), caption: "some points are added, but they are not on the boundary.")

After checking my code again and again, I discovered that the issue comes from the way the angle is calculated in the math library:
```rust
pub fn angle_to(self, rhs: Self) -> f32 {
    let angle = math::acos_approx(
        self.dot(rhs) / math::sqrt(self.length_squared() * rhs.length_squared())
    );
    angle * math::signum(self.perp_dot(rhs))
}
```

The floating points error caused the sign to flip when the points are almost colinear, and thus to be treated as almost 360deg.
]

#pagebreak()

=== Parametrization

The paper did not mention how the step between 2 consecutive sections is calculated.

At the start, I used a fixed number of steps for each branch. This resulted in a non-uniform mesh:

#image("images/mesh_non_uniform_triangulation.png", height: 180pt)

Other issue: When the two children branches are not the same length, it creates a gap between the 2 sides:

#image("images/tree_branch_non_uniform.svg", height: 150pt)

TO solve these issues, I changed the approach and used a step $d z$ computed according to the radius of the branch. To get the spline parameter, I use $t = k (d z) / L$ with $L$ the total branch length.

It was a lot better:

#image("images/mesh_uniform_triangulation.png", height: 180pt)

Unfrtunately, it did not solve all the problems. In particular, there was an irregularity between the last contour of the branch and the first of the next branch.

#figure(image("images/parametrization_gap.png", height: 200pt), caption: "the distance between contours is too small")

The solution was to pass the offset from the end of the branch in the recursive call to the function.

After all of this, I got a good parametrization:

#image("images/good_parametrization.png", height: 300pt)


#note[
    I could have used a parametrization of the spline directly, but it would have been unpractical to compute the branch joins.
]




== Step 4: Meshing and branch fusion <meshing>

The final step is to comupte the mesh (vertices, colors, normals and triangles) of the entire tree

The basic idea is to take 2 contours, and join them together using triangles on the lateral surface. Again, it was not entirely clear from the paper so I had to adapt.

=== Joining contours

I looked at the problem as a minimazation problem.

Given two lists of points $P_i$ and $Q_i$, you want to interlace them in the way that minimizes the sum of the perimeters of the triangles.

The perimeter of the triangle at one moment of the algotihm is one of 2 possibilities:
- $P_i P_(i+1) + P_(i+1) Q_j + Q_j P_i$
- $P_i Q_(j) + Q_j Q_(j+1) + Q_(j+1) P_i$

#figure(image("images/contour_algorithm_1.svg", width: 400pt), caption: "the 2 possible choices of triangles")

Since the terms $P_i P_(i+1)$ and $Q_j Q_(j+1)$ will be added either way, we just have to minimize the sum of $P_i Q_j$.

This leads to the following greedy algorithm:

```
join_contours(P, Q):
    (i, j) = (0, 0)
    if d(P[i], Q[j+1]) < d(P[i+1], Q[j]):
        add the triangle P[i],Q[j],Q[j+1]
        i += 1
    else:
        add the triangle P[i],P[i+1],Q[j]
        j += 1
```

#note[
    At the end, $P_n$ and $P_0$, $Q_n$ and $Q_0$ must be linked together to close the lateral surface.
]

#pagebreak()

=== Branch fusion

There is a very particular case of meshing when 2 branch split: there is one contour above and 2 contours over.

#let n = 20

#let branches() = {
        let points_m = range(1, n - 1)
        .map(i => {
                let θ = i*360deg/n + 180deg/n
                let x = calc.cos(θ)-1.5
                let y = calc.sin(θ)
                (x, y)
                })

        let points_s = points_m.map(pos => {
            (- pos.at(0), pos.at(1))
        })

        let points_aam = range(2, n)
        .map(i => {
                let θ = 45deg + i*260deg/n
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

        for i in range(n - 2) {
        let pos = points_m.at(i)
            node(pos, [$times$], name: "m"+str(i))
        }

        for i in range(n - 2) {
        let pos = points_s.at(i)
            node(pos, [$times$], name: "s"+str(i))
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


The situation can be represented by:
#scale_center(diagram(
    branches(),
    node((0, 1.5), [x]),
    node((0, -1.5), [x]),
    node((-0.5, 0.0), [x]),
    node((0.5, 0.0), [x]),
))

In this situation, we want to mesh everything such that the result has no hole. In order to create the mesh, I start by finding the 2 nearest points from the 2 branches:

#scale_center(diagram(
    branches(),
    node((-0.5, 0.0), [P], fill: red),
    node((0.5, 0.0), [Q], fill: red),
    node((0, 1.5), [x]),
    node((0, -1.5), [x]),
))

And then I find the points on the above contour that are the nearest from them:


#scale_center(diagram(
    branches(),
    node((-0.5, 0.0), [P], fill: red),
    node((0.5, 0.0), [Q], fill: red),
    node((0, 1.5), [A], fill: blue),
    node((0, -1.5), [B], fill: blue),
))

#note[
    To be sure the 2 selected points $A$ and $B$ are not on the same side, I compute the determinant $det(P Q, Q A)$ and det $(P Q, Q B)$ . They must be opposite signs.
]

Once I found all these special points, I can use the previous algorithm to create 3 strips. There are 2 special triangles, represented in gray on the figure.

#diagram(
    mark-scale: 10%,
    branches(),
    node((-0.5, 0.0), [P], fill: red, name: <p>),
    node((0.5, 0.0), [Q], fill: red, name: <q>),
    node((0, 1.5), [A], fill: blue, name: <a>),
    node((0, -1.5), [B], fill: blue, name: <b>),
    edge(<s1>, <as0>, "-", stroke: yellow),
    edge(<m1>, <am0>, "-", stroke: green),
    for i in range(1, n - 3) {
    // refer to nodes by label, e.g., <1>
        edge(label("m"+str(i)), label("am"+str(i)), "-", stroke: green)
        edge(label("m"+str(i)), label("am"+str(i+1)), "-", stroke: green)
        edge(label("s"+str(i)), label("as"+str(i)), "-", stroke: yellow)
        edge(label("s"+str(i)), label("as"+str(i+1)), "-", stroke: yellow)
    },
    edge(<a>, <m1>, stroke: gray),
    edge(<a>, <s1>, stroke: gray),
    edge(<m1>, <s1>, stroke: gray),

    edge(<b>, <m16>, stroke: gray),
    edge(<b>, <s16>, stroke: gray),
    edge(<m16>, <s16>, stroke: gray),

    edge(<m17>, <s17>, stroke: purple),
    edge(<m0>, <s0>, stroke: purple),
    edge(<p>, <q>, stroke: purple),

    edge(<m17>, <s16>, stroke: purple),
    edge(<s0>, <m1>, stroke: purple),
    edge(<p>, <s17>, stroke: purple),
    edge(<m0>, <q>, stroke: purple),
)

=== Normals

The first version of my algorithm computed the normals section by section, computing the bisector for each segment and interpolating. The result was very unstable and did not work when the branch changed direction abruptly. To solve this issue I computed the normals when my mesh was completely generated, using a pre-existing algorithm.

=== Adding leaves

For the visual aspect, I added leaves by creating random triangles at the end of each trunk.

#image("images/green_leaves.png", height: 150pt)


#pagebreak()

== Skeleton generation

For the skeleton generation, I used a simpler version than what has done the paper. I adapted the ideas from another paper that used L-systems to generate the skeleton.

I create the tree from the root, and add 1 or 2 children with random rotations and random branch lengths. Some details:
- the probablity of having 2 childrens increase with depth
- I make sure that the angle between 2 children is not too small, because in reality 2 children branches interact one with another
- the radiuses of the 2 children don't have the same probability distribution.
- I stop when the branch radius is smaller than a treshold.


As I said, I chose to focus on the mesh generation and performance, and not the skeleton generation part.

#pagebreak()
= Analysis

== Performance <performance>

With the quadtree implementation, my algorithm took around 15s to generate 200_000 particles and 16_000 strands. This is slower than the paper. For reference:#footnote[PC is time for particle simulation in seconds, M is time for mesh generation in second]

#image("images/stats_from_paper.png")

I think my particle simulation code is slower, but my mesh generation code is faster. Also, I don't know what machine they used.

Iw would be very valuable to benchmark each phase to know what is the bottleneck, but I did not have time to do it.

== Limitations and improvements

My implementation suffer from some edge cases, and some of them can make the app panic.
    - when there are not enough points to merge
    - when two children branches do not split (solved)
    - when the branches split just before or just after a node in the tree


The most useful improvements in the current state of the project are:
- using a biology-based model to generate the skeleton
- being able to specify custom branch profiles, not just circles.
- benchmarking each phase to know which phase should be optimized
- testing other data-structures to improve the performance of the simulation.
- using paralelization






#pagebreak()
= Gallery <gallery>

I spent some time adding visualization code for each phase of the generation. They can be enabled or disabled in the web interface also.

#for i in range(1, 5) {
    grid(columns: (1fr, 1fr),
    image("images/tree_"+str(i)+"_skeleton.png"),
    image("images/tree_"+str(i)+"_strands.png"),
    image("images/tree_"+str(i)+"_contours.png"),
    image("images/tree_"+str(i)+"_render.png"),
    )
}
