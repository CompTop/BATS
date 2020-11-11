#pragma once

// TV is type of values stored in matrix

/*
TODO:
specify functions that any matrix class should implement
*/


/*
AbstractMatrix

Provides a CRTP base class that specifies operations that all matrix
implementations should support
*/
//template <class Derived>
class AbstractMatrix
{
public:

  // one constructor should take in a csc matrix

  // should support getindex

  // should be able to dump to dense array

  // matrix-vector multiply
  // template <class VecT>
  // VecT matvec(VecT &v) {
  //   return static_cast<Derived *>(this)->matvec(v);
  // }


};

// identity matrix
// template <class TM>
// TM identity(size_t n);
