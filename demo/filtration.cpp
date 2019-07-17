#include <complex/simplicial_complex.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
#include <filtration.h>
#include <iostream>
#include <homology_reduction.h>

#include <linalg/field.h>
#define FF ModP<int, 3>
#define VecT SparseVector<size_t, FF>
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
  std::cout << "\nU:" << std::endl;
  U.print();

  std::vector<size_t> inds({0,2});
  std::vector<FF> vals({FF(1), FF(1)});
  VecT x(inds, vals);
  x.print_row();
  auto y = gemv(U, x);
  y.print_row();

  auto z = ut_solve(U, x);
  z.print_row();

  MatT I = identity<VecT>(3);
  auto Uinv = ut_solve(U, I);
  std::cout << "\nUinv:" << std::endl;
  Uinv.print();

  MatT A = MatT(4, 4);
  A.print();

  return 0;
}
