#include <simplicial_complex.h>
#include <sparse_vector.h>
#include <col_matrix.h>
#include <filtration.h>
#include <iostream>

int main() {

  SimplicialComplex X = SimplicialComplex();
  X.add({0});
  X.add({1});
  X.add({0, 1});
  X.add({2});
  X.add({0,2});
  X.add({2,1});
  X.add({0,1,2});
  X.print();

  // std::cout << "{1} loc : " << X.find_idx({1}) << std::endl;
  // std::cout << "{3} loc : " << X.find_idx({3}) << std::endl;
  // std::cout << "{0,2} loc : " << X.find_idx({0,2}) << std::endl;
  // std::cout << "{0,3} loc : " << X.find_idx({0,3}) << std::endl;

  std::cout << "boundary of edge 0" << std::endl;
  auto v = X.boundary<int>(1, 0);
  v.print();

  std::cout << "boundary of edge 1" << std::endl;
  v = X.boundary<int>(1, 1);
  v.print();

  std::cout << "boundary of triangle" << std::endl;
  v = X.boundary<int>(2, 0);
  v.print();

  auto B = X.boundary<SparseVector<size_t, int>>(1);

  auto F = LowerStarFiltration<SimplicialComplex, float>(X, {1., 2., 3.});

  F.print();

  return 0;

}
