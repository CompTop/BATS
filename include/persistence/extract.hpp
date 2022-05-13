#pragma once
// extract barcode pairs from reduced matrix

#include <map>
#include <vector>
#include <linalg/sparse_vector.hpp>
#include <linalg/col_matrix.hpp>
#include <util/barcode.hpp>
#include "barcode.hpp"

namespace bats {
/**
Extract persistence pairs

@param R reduced matrix in dimension k
@param p2c pivot-to-column map for reduced dimension k+1
@param valsk filtration values in dimension k
@param valsk1 fitration values in dimension k+1
@param k (optional) dimension
*/
template <typename T, typename MT>
std::vector<PersistencePair<T>> extract_pairs(
  const MT &R,
  const std::vector<size_t>& p2c,
  const std::vector<T>& valsk,
  const std::vector<T>& valsk1,
  size_t k=0
)
{
  std::vector<PersistencePair<T>> pairs;
  pairs.reserve(R.ncol());
  for (size_t i =0; i < R.ncol(); i++) {
    if (R[i].nnz() == 0) {
      // homology generated
      if (p2c[i] == bats::NO_IND)  {
        // infinite bar
        pairs.emplace_back(
          k, i, bats::NO_IND,
          valsk[i], std::numeric_limits<T>::infinity()
        );
      } else {
        size_t j = p2c[i];
        // finite bar
        pairs.emplace_back(
          k, i, j,
          valsk[i], valsk1[j]
        );
      }
    }
  }
  pairs.shrink_to_fit();
  return pairs;
}

} // namespace bats
