#include <complex/simplicial_complex.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>
//#include <filtration/filtration.h>
#include <iostream>

int main() {

  std::vector<std::vector<size_t>> nbrs(4);

  std::vector<size_t> edges = {0,1,\
                               0,2,\
                               0,3,\
                               1,2,\
                               1,3,\
                               2,3};
  SimplicialComplex Y = FlagComplex(edges, 4, 3);
  Y.print();

  // X.add({0});
  // X.add({1});
  // auto ret0 = X.add({0, 1});
  // std::cout << ret0 << std::endl;
  // auto ret = X.add({0, 1});
  // std::cout << ret << std::endl;
  // X.add({2});
  // X.add({0,2});
  // X.add({2,1});
  // X.add({0,1,2});
  // X.print();
  //
  // SimplicialComplex Y = SimplicialComplex(3);
  // Y.print();
  //
  // // std::cout << "{1} loc : " << X.find_idx({1}) << std::endl;
  // // std::cout << "{3} loc : " << X.find_idx({3}) << std::endl;
  // // std::cout << "{0,2} loc : " << X.find_idx({0,2}) << std::endl;
  // // std::cout << "{0,3} loc : " << X.find_idx({0,3}) << std::endl;
  //
  // std::cout << "boundary of edge 0" << std::endl;
  // auto v = X.boundary<int>(1, 0);
  // v.print();
  //
  // std::cout << "boundary of edge 1" << std::endl;
  // v = X.boundary<int>(1, 1);
  // v.print();
  //
  // std::cout << "boundary of triangle" << std::endl;
  // v = X.boundary<int>(2, 0);
  // v.print();
  //
  // auto B = X.boundary<SparseVector<size_t, int>>(1);
  // B.print();
  //
  // auto F = LowerStarFiltration<SimplicialComplex, float>(X, {2., 1., 3.});
  //
  // F.print();
  //
  // B = F.boundary<SparseVector<size_t, int>>(1);
  // B.print();
  //
  // B = F.boundary<SparseVector<size_t, int>>(2);
  // B.print();

  return 0;

}
