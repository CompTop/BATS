#include <complex/simplicial_complex.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <filtration.h>
#include <iostream>

int main() {

  Filtration<SimplicialComplex, float> F;

  F.add({0}, 0.0);
  F.add({1}, 1.0);
  F.add({0,1}, 2.0);

  F.sort();
  F.print();

  return 0;
}
