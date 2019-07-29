#pragma once

#include "abstract_complex.h"
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <correspondence/correspondence.h>
#include <util/sorted.h>
#include <util/trie.h>
#include <vector>
#include <cstddef>
#include <iostream>
#include <map>
#include <limits>
#include <algorithm>
#include <util/simplex.h> // for hash function


// no index
#define NO_IND std::numeric_limits<size_t>::max();

/*
simplicial complex implementation
*/
class SimplicialComplex : public AbstractComplex
{
private:
  // maybe make template parameter for container?
  // using spx_map = std::unordered_map<std::vector<size_t>, size_t, SimplexHasher>;
  // using spx_map = std::map<std::vector<size_t>, size_t>;
  using spx_map = SparseTrie<size_t, size_t>;

  // container for simplices of a fixed dimension
  using spx_container = std::vector<std::vector<size_t>>;
  // using spx_container = SimplexContainer;
  // hold number of 0-cells
  size_t ncells0;
  // spx_list[k][i] holds i'th simplex in dimension k+1
  std::vector<spx_container> spx_list;
  // keeps track of how many cells are in given dimension
  std::vector<size_t> _ncells;
  // map to find simplices when forming boundary
  spx_map spx_to_idx;
  // TODO: if forming full clique complex
  // can give all simplices unique location for lookup

  // returns index of simplex
  size_t find_idx(std::vector<size_t> &s) {
    size_t dim = s.size() - 1;
    if (dim == 0) {
      return s[0] < ncells0 ? s[0] : NO_IND;
    } else {
      if (spx_to_idx.count(s)) {
        return spx_to_idx[s];
      }
    }
    return NO_IND;
  }

public:

  // default constructor
  SimplicialComplex() : ncells0(0) {}

  // constructor that initializes to set dimension
  SimplicialComplex(size_t maxdim) : ncells0(0), _ncells(maxdim+1, 0) {
    spx_list.reserve(maxdim);
    for (size_t i = 0; i < maxdim; i++) {
      // spx_list.emplace_back(spx_container(i+1));
      spx_list.emplace_back(spx_container());
    }
  }

   bool set_ncells0(size_t n) {
     ncells0 = n;
     return true;
   }


  // adds a simplex to the complex without doing any checks
  bool add_unsafe(std::vector<size_t> &s) {
    size_t dim = s.size() - 1;
    if (dim == 0) {
      ncells0 = std::max(s[0] + 1, ncells0);
      spx_to_idx.emplace(s, _ncells[dim]);
    } else {
      // add simplex to appropriate dimension
      spx_list[dim-1].emplace_back(s);
      // set reverse map
      spx_to_idx.emplace(s, _ncells[dim]);
    }
    _ncells[dim]++;
    return true;
  }

  // add simplex to complex with appropriate checks
  bool add(std::vector<size_t> &s) {
    size_t dim = s.size() - 1;
    if (dim == 0){
      ncells0 = std::max(s[0] + 1, ncells0);
      return true;
    }

    // ensure simplex is sorted
    std::sort(s.begin(), s.end());

    // add dimensions if necessary
    while (spx_list.size() < dim) {
      // spx_list.push_back(spx_container(spx_list.size() + 1));
      spx_list.emplace_back(spx_container());
    }

    // check if simplex is already in complex
    if (spx_to_idx.count(s) > 0) {
      // simplex is already in complex
      return false;
    }

    return add_unsafe(s);
  }

  void add_dimension_recursive_flag_unsafe(
          std::vector<std::vector<size_t>> &nbrs,
          size_t d, size_t maxd,
          std::vector<size_t> &iter_idxs,
          std::vector<size_t> &spx_idxs) {

    if (d == maxd) {
      // maximum dimension
      for (size_t j : iter_idxs) {
        // add each simplex
        spx_idxs.emplace_back(j);
        std::vector<size_t> s(spx_idxs);
        std::sort(s.begin(), s.end()); // TODO: is this sorted?
        add_unsafe(s);
        spx_idxs.pop_back();
        // no recursion
      }
    } else {
      std::vector<size_t> iter_idxs2; // allocate once
      iter_idxs2.reserve(iter_idxs.size());
      for (size_t j : iter_idxs) {
        spx_idxs.emplace_back(j);
        std::vector<size_t> s(spx_idxs);
        std::sort(s.begin(), s.end());
        add_unsafe(s);

        // recurse
        // get list of indices to recurse on
        intersect_sorted_lt(iter_idxs, nbrs[j], j, iter_idxs2);
        if (!(iter_idxs2.empty())) {
          // only recurse if there is something to do
          add_dimension_recursive_flag_unsafe(nbrs, d+1, maxd, iter_idxs2, s);
        }

        // ready for next simplex
        spx_idxs.pop_back();
      }
    }

  }

  // return 0-skeleton of cell
  // std::vector<size_t> skeleton0(size_t dim, size_t i) {
  //   if (dim == 0) {
  //     return {i};
  //   }
  //   return spx_list[dim-1][i];
  // }


  const size_t ncells() const {
    size_t ct = 0;
    for (size_t k = 0; k < maxdim() + 1; k++) {
      ct += ncells(k);
    }
    return ct;
  }

  const size_t ncells(size_t dim) const {
    return dim == 0 ? ncells0 : _ncells[dim];
  }

  const size_t maxdim() const {
    return spx_list.size();
  }

  void print_dims() {
    std::cout << "maxdim = " << maxdim() << std::endl;
    for (size_t dim = 0; dim < maxdim() + 1; dim ++) {
      std::cout << "dim " << dim << " : " << ncells(dim) << std::endl;
    }
  }

  // void print_cell(size_t dim, size_t i) {
  //   if (dim == 0) {
  //     std::cout << i << std::endl;
  //     return;
  //   }
  //   for (size_t j = 0; j < dim; j++) {
  //     std::cout << spx_list[dim - 1][i][j] << ",";
  //   }
  //   std::cout << spx_list[dim - 1][i][dim] << std::endl;
  // }

  // void print(size_t dim) {
  //   if (dim == 0) {
  //     for (size_t i = 0; i < ncells0; i++) {
  //       print_cell(0, i);
  //     }
  //   } else {
  //     for (size_t i = 0; i < ncells(dim); i++) {
  //       print_cell(dim, i);
  //     }
  //   }
  // }

  // void print() {
  //   std::cout << "SimplicialComplex of dimension " << maxdim() << std::endl;
  //   for (size_t dim = 0; dim < maxdim() + 1; dim ++) {
  //     std::cout << "dim " << dim << " : " << ncells(dim) << " cells" << std::endl;
  //     print(dim);
  //   }
  // }


  // get boundary of simplex i in dimension dim
  // template <typename TV>
  // SparseVector<size_t, TV> boundary(size_t dim, size_t i) {
  //   //assert(dim > 0);
  //   std::vector<size_t> bdr_ind;
  //   std::vector<int> bdr_val;
  //   std::vector<size_t> face;
  //   int coeff = -1;
  //   // loop over faces in lexicographical order
  //   for (size_t k = 0; k < dim+1; k++) {
  //     size_t k2 = dim-k; // index to skip
  //     face.clear();
  //     for (size_t j = 0; j < k2; j++) {
  //       face.push_back(spx_list[dim-1][i][j]);
  //     }
  //     for (size_t j = k2+1; j < dim+1; j++) {
  //       face.push_back(spx_list[dim-1][i][j]);
  //     }
  //
  //     bdr_ind.push_back(find_idx(face));
  //     bdr_val.push_back(coeff);
  //     coeff = -coeff;
  //   }
  //   return SparseVector<size_t, TV>(bdr_ind, bdr_val);
  // }

  // template over column type
  // template <class TVec>
  // ColumnMatrix<TVec> boundary(size_t dim) {
  //   //assert(dim > 0);
  //   std::vector<TVec> col;
  //   std::vector<size_t> bdr_ind;
  //   std::vector<int> bdr_val;
  //   std::vector<size_t> face;
  //   for (size_t i = 0; i < ncells(dim); i++) {
  //     bdr_ind.clear();
  //     bdr_val.clear();
  //     int coeff = -1;
  //     // loop over faces in lexicographical order
  //     for (size_t k = 0; k < dim+1; k++) {
  //       size_t k2 = dim-k; // index to skip
  //       face.clear();
  //       for (size_t j = 0; j < k2; j++) {
  //         face.push_back(spx_list[dim-1][i][j]);
  //       }
  //       for (size_t j = k2+1; j < dim+1; j++) {
  //         face.push_back(spx_list[dim-1][i][j]);
  //       }
  //
  //       bdr_ind.push_back(find_idx(face));
  //       bdr_val.push_back(coeff);
  //       coeff = -coeff;
  //     }
  //     col.push_back(TVec(bdr_ind, bdr_val));
  //   }
  //   return ColumnMatrix<TVec>(col);
  // }

  // return indices of simplices in dimension dim whose vertex set is all in vtx_list
  std::vector<size_t> sub_complex(std::vector<size_t> &vtx_list, size_t dim);
  // step 1: sort vtx_list
  // step 2: run through all dim simplices and see if simplex is a subset of vtx_list

  // construct a quotient complex by a relation
  // x_1 ~ x_2 if f(x_1) = f(x_2)
  SimplicialComplex quotient(Function &R);

};

// return a Flag Complex using an adjacency list of vertices.
// complex is constructed up to maxdim
// ASSUME:
// if (i,j) in edge set and j < i, then j in nbrs[i], and i not in nbrs[j]
SimplicialComplex FlagComplex(std::vector<std::vector<size_t>> &nbrs, size_t maxdim) {
  SimplicialComplex X = SimplicialComplex(maxdim);
  // sets 0-cells
  X.set_ncells0(nbrs.size());

  std::vector<size_t> spx_idxs(2);
  std::vector<size_t> iter_idxs;
  iter_idxs.reserve(nbrs.size()); // maximum size

  for (size_t i = 0; i < nbrs.size(); i++) {
    std::cout << "i = " << i << std::endl;
    // make sure neighbor list is sorted
    std::sort(nbrs[i].begin(), nbrs[i].end());
    for (size_t j : nbrs[i]) {
      if (j > i) {break;}
      std::cout << "  j = " << j << std::endl;
      spx_idxs[0] = j;
      spx_idxs[1] = i;
      X.add_unsafe(spx_idxs);

      // get common neighbors for higher dim simplices
      intersect_sorted_lt(nbrs[i], nbrs[j], i, iter_idxs);
      std::cout << "iter idxs:" << std::endl;
      for (auto val : iter_idxs) {
        std::cout << val << std::endl;
      }
      X.add_dimension_recursive_flag_unsafe(nbrs, 2, maxdim, iter_idxs, spx_idxs);
    }
  }
  return X;
}

// Flag complex using list of edges
// (edges[2*k], edges[2*k+1]) = (i, j) is an edge
// n - number of vertices
// maxdim - maximum dimension of simplices
SimplicialComplex FlagComplex(std::vector<size_t> &edges, size_t n, size_t maxdim) {
  SimplicialComplex X = SimplicialComplex(maxdim);
  // sets 0-cells
  X.set_ncells0(n);

  std::vector<std::vector<size_t>> nbrs(n);

  std::vector<size_t> spx_idxs(2);
  std::vector<size_t> iter_idxs;
  iter_idxs.reserve(n); // maximum size

  size_t m = edges.size() / 2;
  for (size_t k = 0; k < m; k++) {
    size_t i = edges[2*k];
    size_t j = edges[2*k + 1];
    spx_idxs[0] = i;
    spx_idxs[1] = j;
    X.add_unsafe(spx_idxs);

    intersect_sorted(nbrs[i], nbrs[j], iter_idxs);

    nbrs[i].emplace_back(j);
    std::sort(nbrs[i].begin(), nbrs[i].end());
    nbrs[j].emplace_back(i);
    std::sort(nbrs[j].begin(), nbrs[j].end());

    X.add_dimension_recursive_flag_unsafe(nbrs, 2, maxdim, iter_idxs, spx_idxs);
  }

  return X;
}
//
// // construct a Clique Complex on n vertices up to dimension k
// SimplicialComplex CliqueComplex(size_t n, size_t k) {
//
// }
