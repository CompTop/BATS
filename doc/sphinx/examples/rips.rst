Rips Complexes
==============

A Rips complex is constructed on a finite metric space with
a fixed distance parameter.  All pairs of points which have a smaller
distance than the specified parameter are connected by an edge, and the
maximal simplicial complex (up to a specified dimension) is constructed
on this edge set (i.e. a flag complex or clique complex).

First, let's generate a dataset:

.. code-block: cpp

   // in header
   using namespace bats;

   size_t d = 2; // dimension of Euclidean Space
   size_t n = 350; // size of point cloud
   DataSet<double> x = = sample_cube<double>(d, n, 0); // seed with 0

Standard Rips Complexes
-----------------------

now we can construct a Rips complex.  It can either be a SimplicialComplex
or LightSimplicialComplex.

.. code-block:: cpp

   // using CpxT = SimplicialComplex;
   using CpxT = LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;

   size_t maxdim = 3; // maximum simplex dimension
   double rmax = 0.2; // connectivity radius

   auto X = RipsComplex<CpxT>(x, LInfDist(), rmax, maxdim);

Rips Filtrations
----------------

A Rips filtration can be constructed with an almost identical call.  The difference
is that now you keep track of the parameter at which each simplex appears in the filtration
and can compute persistent homology.

.. code-block:: cpp

   auto F = bats::RipsFiltration<CpxT>(X, dist, rmax, maxdim);
