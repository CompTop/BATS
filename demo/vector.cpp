#include <linalg/field.h>
#include <linalg/sparse_vector.h>
#include <vector>

#define F2 ModP<int, 2>
#define VecT SparseVector<size_t, F2>

int main() {

  std::vector<size_t> inds({1,3});
  std::vector<F2> vals({F2(1), F2(1)});

  // for (size_t j = 0; j < 2; j++) {
  //   std::cout << inds[j] << " : " << vals[j] << std::endl;
  //   auto p = std::make_pair(inds[j], vals[j]);
  //   std::cout << p.first << " : " << p.second << std::endl;
  // }

  VecT x(inds, vals);
  x.print_row();

  VecT y({0, 1}, {F2(1), F2(1)});
  y.print_row();

  y.axpy(F2(1), x);
  y.print_row();

  VecT z(3);
  z.print_row();


  return 0;
}
