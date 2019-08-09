#pragma once
// extract barcode pairs from reduced matrix

#include <map>
#include <vector>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <util/barcode.h>

/*
returns vector of pairs
two entries = two pairs
*/
std::vector<size_t> extract_pairs(
  ColumnMatrix<TVec> &M,
  std::map<size_t, size_t> pivot_to_col)
{
    std::vector<size_t> p;
}
