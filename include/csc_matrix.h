#pragma once

#include <vector>
#include <cstddef>
#include "abstract_matrix.h"

// template over index type and value type
template <typename TI, typename TV>
class CSCMatrix : public AbstractMatrix
{
public:
  std::vector<TI> colptr;
  std::vector<TI> rowind;
  std::vector<TV> val;

};
