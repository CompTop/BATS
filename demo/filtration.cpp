#include <complex/simplicial_complex.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <filtration.h>
#include <iostream>
#include <homology_reduction.h>

#include <linalg/field.h>
#define F2 ModP<int, 2>

int main() {

  Filtration<SimplicialComplex, float> F;

  F.add({0}, 0.0);
  F.add({1}, 1.0);
  F.add({2}, 2.0);
  F.add({0,1}, 2.0);
  F.add({0,2}, 3.0);
  F.add({1,2}, 4.0);

  F.sort();
  F.print();

  std::cout << "\n1 boundary:" << std::endl;
  auto B = F.boundary<SparseVector<size_t, F2>>(1);
  B.print();

  // auto v = B[0];
  // v.print();

  std::cout << "\nreduced\n" << std::endl;
  auto p2c = reduce_matrix(B);
  B.print();





  return 0;
}
