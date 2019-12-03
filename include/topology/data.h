#pragma once

#include <linalg/naive_dense.h>

// data matrix - d x n
// where d is dimension, n is number of points
template<typename T>
using Matrix = A<Dense<T,ColMaj>>;
