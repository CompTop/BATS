#pragma once
/*
  Abstract cell complex class that all derived complexes should inherit from
*/

// #include <vector>
#include <cstddef>
#include <string>
#include <iostream>
#include <sstream>

// indicates dimension and index of a cell
struct cell_ind {
public:
    size_t dim;
    size_t ind;

    cell_ind() {}
    cell_ind(size_t dim, size_t ind) : dim(dim), ind(ind) {}

    std::string str() {
        std::ostringstream oss;
        oss << '(' << dim << ',' << ind << ')';
        return oss.str();
    }
};

// //#include "csc_matrix.h"
// //#include "sparse_vector.h"
//
// class AbstractComplex
// {
// public:
//   // adds cell with index s to complex in dimension dim
//   // add_cell(size_t s, size_t dim);
//   // will likely need to extend methods
//   // examples: if boundary and dimension can be inferred, can leave out arguments
//   //           if bdr should take coefficents
//
//   // return number of cells in dimension dim
//   size_t ncells(size_t dim);
//
//   // return maximal dimension of a cell in complex
//   size_t maxdim();
//
//   // return boundary of i'th cell in dimension dim in SparseVector<int> format
//   //SparseVector<size_t, int> boundary_int(size_t dim, size_t i);
//
//   // return boundary of i'th cell in dimension dim in VecT format
//   //template <class VecT>
//   //VecT boundary(size_t dim, size_t i) = VecT(boundary_int(dim, i));
//
//   // return dimension dim boundary in MatT format
//   //CSCMatrix<size_t, int> boundary_csc(size_t dim);
//
//   // return dimension dim boundary in MatT format.
//   // default to converting csc integer matrix,
//   // but specialized implementation should be provided if appropriate
//   //template <class MatT>
//   //MatT boundary(size_t dim) = MatT(boundary_csc(dim));
//
// };
