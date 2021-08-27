.. quick start guide:

Quick Start Guide
=================

This guide will cover several basic use cases for BATS.
More specialized functionality is covered in the tutorials and examples.
Ultimately, you can use the API reference.

Simplicial Complexes and Homology
---------------------------------

BATS offers two implementations of simplicial complexes: SimplicialComplex
and LightSimplicialComplex.  While the internal representations differ, they
both have the same interface which can be used.  When dealing with stand-alone
simplices, BATS uses :code:`std::vector<size_t>` to represent simplices, from which
the vertex set and dimension of the simplex can be extracted.

.. warning::
    Simplices in bats should generally be assumed to be *ordered*, meaning that
    :code:`{0,1,2}` is not the same as :code:`{1,2,0}`.  If you want to use
    *unordered* simplices, you can either add vertices in sorted order, or use
    a sorting algorithm before adding simplices to complexes.

The two main methods for adding simplices to simplicial complexes are :code:`add`,
which assumes you have already added the boundary of a simplex, and :code:`add_recursive`
which will add any faces that are not already present.

.. code-block:: cpp


   bats::SimplcialComplex X();
   X.add_recursive({0,1,2});
   X.add_recursive({2,3});
   X.add({1,3});

   X.print_summary();

The above code will create a SimplicialComplex with a single connected component
and a single hole.  The call to :code:`X.print_summary()` will produce

.. code-block:: none

   SimplicialComplex, maxdim = 2
      dim 0 : 4 cells
      dim 1 : 5 cells
      dim 2 : 1 cells
   10 total

Let's now compute homology.

.. code-block:: cpp


   using F2 = ModP<int, 2>;
   auto R = bats::Reduce(X, F2());

   R.print_summary();

The output of :code:`bats::Reduce` will be a ReducedChainComplex with F2 coefficients,
which holds information used to compute homology.
The output of :code:`R.print_summary()` will be

.. code-block:: none

   ReducedChainComplex, maxdim =  2
      dim 0: 4, betti_0: 1
      dim 1: 5, betti_1: 1
      dim 2: 1, betti_2: 0



Persistent Homology
-------------------

A filtration in BATS is a class which is templated over the
type of the filtration values as well as the type of the underlying complex.

.. code-block:: cpp


   bats::Filtration<double, bats::SimplicialComplex> F;
   std::vector<size_t> spx;
   spx = {0,1,2}; F.add_recursive(0.0, spx);
   spx = {2,3};   F.add_recursive(1.0, spx);
   spx = {1,3};   F.add(2.0, spx);

   F.complex().print_summary();

The underlying SimplicialComplex is the same as in the previous example:

.. code-block:: none

   SimplicialComplex, maxdim = 2
      dim 0 : 4 cells
      dim 1 : 5 cells
      dim 2 : 1 cells
   10 total

Again, we can use the :code:`Reduce` function to compute homology.  Because we
are using a filtration as input, persistent homology will be computed, returning
a ReducedFilteredChainComplex.

.. code-block:: cpp


   using F2 = ModP<int, 2>;
   auto R = bats::Reduce(F, F2());

   for (size_t k = 0; k < R.maxdim(); ++k) {
      std::cout << "\n" << k << "-dimensional barcode:\n";
      for ( auto p : R.persistence_pairs(k)) {
         std::cout << p.str() << std::endl;
      }
   }

The output will show one persistent 0-dimensional homology class
as well as one persistent 1-dimensional homology class

.. code-block:: none

   0-dimensional barcode:
   0 : (0,inf) <0,-1>
   0 : (0,0) <1,0>
   0 : (0,0) <2,1>
   0 : (1,1) <3,3>

   1-dimensional barcode:
   1 : (0,0) <2,0>
   1 : (2,inf) <4,-1>


The output of :code:`R.persistence_pairs(k)` is a vector of PersistencePairs
for k-dimensional persistent homology.

A PersistencePair includes 5 pieces of information:
* The dimension of the homology class.
* The birth and death parameters of the homology class.
* The simplex indices responsible for birth and death.

Maps
----

BATS makes dealing with maps between topological spaces and associated chain maps
and induced maps on homology easy.  The relevant class is a CellularMap which
keeps track of what cells in one complex map to what cells in another.

We'll just look at a wrapper for CellularMap, called SimplcialMap which can be used
to extend a map on the vertex set of a SimplicialComplex to a map of simplices.

First, we'll build identical simplicial complexes X and Y which are both cycle graphs
on four vertices.

.. code-block:: cpp


   bats::SimplicialComplex X;
   X.add_recursive({0,1});
   X.add_recursive({1,2});
   X.add_recursive({2,3});
   X.add_recursive({0,3});
   bats::SimplicialComplex Y = X; // copy

We then build a simplicial map from X to Y which is extended
from a reflection of the vertices.

.. code-block:: cpp


   std::vector<size_t> f0{2,1,0,3}; // reflection map
   auto F = bats::SimplicialMap(X, Y, f0);

The map is extended by sending vertex :code:`i` in X to vertex
:code:`f0[i]` in Y.  Next, we can apply the chain functor.  We'll use
F3 coefficients.

.. code-block:: cpp


   // apply the chain functor
   using F3 = ModP<int, 3>;
   auto CX = bats::Chain(X, F3());
   auto CY = bats::Chain(Y, F3());
   auto CF = bats::Chain(F, F3());

Finally, we can compute homology of the chain complexes and the induced maps.

.. code-block:: cpp


   auto RX = bats::ReducedChainComplex(CX);
   auto RY = bats::ReducedChainComplex(CY);
   RX.print_summary();
   RY.print_summary();

   for (size_t k = 0; k < 2; k++) {
      std::cout << "\nInduced map in dimension " << k << std::endl;
      auto HF = bats::induced_map(CF, RX, RY, k);
      HF.print();
   }

The following output will be produced:

.. code-block:: none

   ReducedChainComplex, maxdim =  1
      dim 0: 4, betti_0: 1
      dim 1: 4, betti_1: 1
   ReducedChainComplex, maxdim =  1
      dim 0: 4, betti_0: 1
      dim 1: 4, betti_1: 1

   Induced map in dimension 0
   [0x7fff6f336460] : 1 x 1 ColumnMatrix
      1

   Induced map in dimension 1
   [0x7fff6f336460] : 1 x 1 ColumnMatrix
      2

As expected, we see that X and Y both have 1-dimensional homology in
dimensions 0 and 1.  The induced map in dimension 0 is the identity, and the
induced map in dimension 1 is multiplication by -1 (in mod-3 coefficients).


Zigzag Homology
---------------

We'll now compute a simple zigzag barcode, using the above example.  We'll consider
a diagram with two (identical) spaces, connected by a single edge which applies
the reflection map in the above example.

.. code-block:: cpp


   bats::Diagram<bats::SimplicialComplex, bats::CellularMap> D(2,1);

   bats::SimplicialComplex X;
   X.add_recursive({0,1});
   X.add_recursive({1,2});
   X.add_recursive({2,3});
   X.add_recursive({0,3});

   std::vector<size_t> f0{2,1,0,3}; // reflection map
   auto F = bats::SimplicialMap(X, X, f0);

   D.set_node(0, X);
   D.set_node(1, X);
   D.set_edge(0, 0, 1, F); // edge 0: (0,1)

We can then apply the Chain and Hom functors, to obtain a diagram
of homology spaces and maps between them

.. code-block:: cpp


   using F3 = ModP<int, 3>;
   auto CD = bats::Chain(D, F3());

   auto HD = bats::Hom(CD, (size_t) 1); // homology in dimension 1

Finally, we extract the barcode from the homology diagram

.. code-block:: cpp


   auto ps = barcode(HD, 1);
   for (auto p : ps) {
      std::cout << p.str() << std::endl;
   }

The output should look like: `1 : (0,1) <0,0>`.
This indicates there is a 1-dimensional homology bar, which is born in the space
with index 0 and survives until the space with index 1.  The `<0,0>` indicates which
generators are associated with the homology class in the diagram.


Source Code
-----------

.. literalinclude:: ../../demo/quickstart.cpp
   :language: c++
   :linenos:
