# BATS: Basic Applied Topology Subprograms

[![Build Status](https://travis-ci.com/bnels/BATS.svg?branch=master)](https://travis-ci.com/bnels/BATS)

This repository contains header files for
* Representing topological spaces on a computer
* Computing algebraic invariants of topological spaces, and diagrams of topological spaces.

The contents of this repository were originally developed to implement the computational framework in the paper
[Persistent and Zigzag Homology: A Matrix Factorization Viewpoint](https://arxiv.org/abs/1911.10693) by Gunnar Carlsson, Anjan Dwaraknath, and Bradley J. Nelson.  Portions of the code base is still under development - please contact the authors (or submit an issue) if you would like to have certain functionality to support an application!

Python bindings for a subset of BATS can be found in a separate repository: [BATS.py](https://github.com/bnels/BATS.py).  

Documentation:

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
auto X = sample_cube<double>(d, 1000);

/// project onto first coordinate
auto p = coordinate_projection(X, 0);

// create open cover using projection
auto cover = uniform_interval_cover(p, 20);

// Diagram of Sets and inclusions
auto SetDgm = linear_cover_union_diagram(cover);

// Diagram of Spaces and maps
auto TopDgm = Rips(SetDgm, X, Euclidean(), rmax, 2);
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


# Fields

BATS linear algebra routines are templated over the choice of field, so you can choose between provided field types, or implement your own.

## Field mod a prime

BATS contains a type for mod-p arithmetic: `ModP<IntT, P>`.  `IntT` is the type of the underlying data, and `P` should be a prime number (we check with a `static_assert` that it is indeed prime at compile time). For example,
```cpp
ModP<int, 3>
```
will perform arithmetic with the integers mod 3.

## Rationals

BATS supports rational arithmetic with an provided `Rational<IntT>` template.

When computing with the in-built rational type, it is possible that the factorizations used in the quiver will create larger and larger numerators and denominators, eventually leading to overflow.  To avoid this, we recommend `Rational<int64_t>` types (instead of 32 bit-types).  You may also consider using big-int type packages which will hold arbitrarily large integers.

If you want to use the in-built rational type, you can pass in a compiler flag `-DWARN_RATIONAL_OVERFLOW` which will check whether `Rational<int64_t>` and `Rational<int32_t>` types are in danger of overflow (other integer types will not perform checks).

We should note that on "real" problems, induced maps on homology don't seem to run into these overflow issues.  These issues have been observed when testing the quiver algorithms on random matrices.  If you find an example in the wild, we're interested to know!

# Parallelism

This package supports OpenMP by default.  If you want to specify the number of threads (instead of using your environment default), try setting the environment variable `OMP_NUM_THREADS`, e.g.
```
time OMP_NUM_THREADS=2 ./zigzag_cylinder.out
```

# IO

We support reading/writing various objects in the following formats

Our convention is to save object `X` using a `save` method.
```cpp
std::string fname = "...";
X.save(fname);
```

all file names will begin with the type of the object.

## SimplicialComplexes

You can write a simplicial complex to a text file using
```cpp
X.save(fname); // std::string fname
```

The format is a comma separated value file with a simplex on each line.  E.g. the simplicial complex with simplices `{0}, {1}, {0,1}` will be written to a file
```
SimplicialComplex
0,
1,
0,1,
```
by default, `.scpx` is used as the file extension

## Sparse Matrices

Sparse matrices can be stored in index-value format `i:v`, where each column of the matrix gets its own row.  The size of the matrix is encoded on the first line as `m, n`. For example, the matrix
```
0  1
1  0
-1  0
```
would be stored as
```
Sparse Matrix
3,2
1:1,2:-1,
0:1
```

By default, `.smtx` is used as the file extension.

## CellularMap
A CellularMap is saved as several sparse matrices in a single file.

If the 0-map is n columns, then after n rows of the file, an entry for the 1-map will begin.

## Diagrams
If the nodes and edges of a diagram have a `save` method on the objects, then the diagram can be saved.

Diagrams are saved in a `*.dgm` folder.  Every node will be saved as `*.dgm/node<i>`, where `<i>` is replaced by the node index in the diagram.  Every edge will be saved as `*.dgm/edge<j>` where `<j>` is the edge index in the diagram.  Metadata for the diagram can be found in `*.dgm/metadata`, which is a text file where the first line is the number of nodes, the second line is the number of edges, and then a line for each edge with contains the source and target of each edge as a comma-separated pair.

For example, a diagram of two simplicial complexes and cellular maps saved in `test.dgm` will have the following directory structure
```
test.dgm/
	metadata
	node0
	node1
	edge1
```
where `node0` and `node1` each contain a saved `SimplicialComplex`, and `edge1` is the saved `CellularMap`.  If the edge is from node 0 to node 1, then the metadata file will contain
```
2
1
0,1
```

# Filtrations

BATS provides a filtration wrapper for complexes, chain complexes, and reduced chain complexes.


```cpp
std::vector<size_t> spx;
Filtration<double, SimplicialComplex> F;
spx = {0}; F.add(0.0, spx);
spx = {1}; F.add(0.0, spx);
spx = {2}; F.add(0.0, spx);
spx = {0,1}; F.add(1.0, spx);
spx = {0,2}; F.add(1.0, spx);
spx = {1,2}; F.add(1.0, spx);
```

Extracting homology
```cpp
auto FC = FilteredChainComplex<double, MT>(F);

auto RFC = ReducedFilteredChainComplex(FC);

// persistence pairs for H1
auto ps = RFC.persistence_pairs(1);

for (auto p : ps) {
	std::cout << p.str() << std::endl;
}
```


# Contributing
Code should use C++17 standard

Much of the code should be written in a templated manner in a header file.
