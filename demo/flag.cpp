#include <filtration/flag.h>
#include <iostream>

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

  return 0;
}
