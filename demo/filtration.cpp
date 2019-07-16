#include <complex/simplicial_complex.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <filtration.h>
#include <iostream>
#include <homology_reduction.h>

#include <linalg/field.h>
#define F2 ModP<int, 2>
#define VecT SparseVector<size_t, F2>
#define MatT ColumnMatrix<VecT>

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
  auto B = F.boundary<VecT>(1);
  B.print();

  // auto v = B[0];
  // v.print();

  std::cout << "\nreduced:" << std::endl;
  auto p2c = reduce_matrix(B);
  B.print();

  std::cout << "\nidentity:" << std::endl;
  MatT U = identity<VecT>(3);
  U.print();

  B = F.boundary<VecT>(1);
  p2c = reduce_matrix(B, U);
  U.print();

  auto y = gemv(U, VecT(2));
  y.print_row();



  return 0;
}
