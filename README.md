# BATS: Basic Applied Topology Subprograms

This repository contains header files for
* Representing topological spaces on a computer
* Computing algebraic invariants of topological spaces, and diagrams of topological spaces.

The contents of this repository were originally developed to implement the computational framework in the paper
[Persistent and Zigzag Homology: A Matrix Factorization Viewpoint](https://arxiv.org/abs/1911.10693) by Gunnar Carlsson, Anjan Dwaraknath, and Bradley J. Nelson.  Portions of the code base is still under development - please contact the authors if you would like to have certain functionality to support an application!

## Features

The key features that distinguish this package are

* Computation of general induced maps on homology
* Computing persistent and zigzag homology using induced maps

We hope to add more attractions in future work. The code base is written with the following principles in mind:

* Header-only - simply include in your project to get started.
* Modularity - topological spaces, chain complexes, quiver representations etc. are all separate classes - if you want to use an optimized version of something, you simply need to provide the right interfaces
* Parallelism - We use OpenMP to take advantage of multiple cores.
* Templates - (almost) everything is templated.  Want to compute with different field coefficients?  Simply change a template parameter!

# Getting Started

To include the BATS headers in your project, simply make sure the project include folder is on your path, then
```cpp
#include <bats.h>
```
When you compile, use `-std=c++17 -fopenmp`.  You can use `demo/makefile` to get started, and modify for your needs.

## Demo

Some examples are in the `demo` folder.  To make a demo, simply go to the folder and type something like
```bash
make zigzag_cover.out
```
where `zigzag_cover.out` pattern matches `zigzag_cover.cpp`.

We'll now look a bit more closely at how thing work

## Build a diagram of spaces

We first build a diagram of spaces.  This could potentially be done in many ways (depending on the application).  In `zigzag_cover.cpp`, we use a cover of the cube, and consider inclusion maps from sets to unions with the adjacent sets
```cpp
// this just creates data
auto x = sample_cube<double>(d, 1000);

/// project onto first coordinate
auto p = coordinate_projection(x, d, 0);

// create open cover using projection
auto cover = uniform_interval_cover(p, 20);

// Diagram of Sets and inclusions
auto SetDgm = linear_cover_union_diagram(cover);

// Diagram of Spaces and maps
auto TopDgm = Rips(SetDgm, x, d, rmax, 2);
```
The diagram of spaces is now held in `TopDgm`.

## Compute a diagram of vector spaces via homology
We will now get the quiver representation of the induced maps on homology.  This goes through two steps - first, we create a diagram of chain complexes.  Then we compute induced maps on homology (in dimesnion 1 in this example).

```cpp
// diagram in Chain
auto ChainDgm = Chain<MT>(TopDgm);

// diagram in Homology in dimension 1
auto HkDgm = Hom(ChainDgm, 1);
```

The matrices used in the chain complex and reduction are determined by the template parameter `MT`, which at the beginning of the program was defined via
```cpp
#define FT ModP<int, 3>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>
```
So homology will be computed with coefficients in `F3`.





# Contributing
Code should use C++17 standard

Much of the code should be written in a templated manner in a header file.
