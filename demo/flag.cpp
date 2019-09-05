#include <filtration/flag.h>
#include <iostream>

#include <linalg/field.h>
#include <linalg/sparse_vector.h>
#include <linalg/col_matrix.h>

#define FF ModP<int, 3>
#define VecT SparseVector<size_t, FF>
#define MatT ColumnMatrix<VecT>

int main() {

  std::vector<size_t> edges = {0, 1,\
                               1, 2,\
                               0, 2,\
                               1, 3,\
                               2, 3,\
                               0, 3};
  std::vector<float> t = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

  auto F = FlagFiltration(edges, t, 0.0f, 4, 3);

  F.sort();
  F.print();

  std::cout << F.ncells() << std::endl;
  std::cout << F.maxdim() << std::endl;

  std::cout << "\n1 boundary:" << std::endl;
  auto B1 = F.boundary<VecT>(1);
  B1.print();

  std::cout << "\n2 boundary:" << std::endl;
  auto B2 = F.boundary<VecT>(2);
  B2.print();

  std::cout << "\n3 boundary:" << std::endl;
  auto B3 = F.boundary<VecT>(3);
  B3.print();


  return 0;
}
