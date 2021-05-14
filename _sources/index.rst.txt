.. Basic Applied Topology Subprograms documentation master file, created by
   sphinx-quickstart on Wed May 12 21:00:14 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Basic Applied Topology Subprograms
==================================

The Basic Applied Topology Subprograms (BATS) are a C++ header-only library for
applied algebraic topology.  It includes functionality for

* Creation of simplicial, cubical, and cell complexes, as well as cellular maps
* Implementation of chain and homology functors
* Filtered complexes and persistent homology
* Diagrams of spaces and maps, and computation of zigzag homology
* Topological constructions such as Vietoris-Rips complexes, Witness complexes, Nerves of covers
* Sparse linear algebra over (finite) fields
* And more!

BATS attempts to balance performance with extensibility.  The library uses template
metaprogramming to make it easy to swap in and out different constructions
while maintaining a consistent core functionality.

Why BATS?
---------
There are many very high-performance libraries for computing things like persistent
homology that have been developed over the past decade.  Unlike many of these libraries
BATS is focused on functorality, and provides functionality to handle maps between
topological spaces, chain maps, and induced maps on homology.  The goal is to make
it easier for researchers and practitioners to implement and explore the vast back
catalog of algebraic topology while also providing applied functionality.



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   quickstart
   tutorials/tut_menu
   examples/ex_menu
   api/library_root
   about



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
