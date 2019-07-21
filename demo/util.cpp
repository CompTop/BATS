#include <util/sorted.h>
#include <vector>
#include <iostream>

int main() {
  std::vector<int> a = {1,2,3};
  std::vector<int> b = {1, 3, 5};
  std::vector<int> c;
  intersect_sorted(a, b, c);
  std::cout << "c is size " << c.size() << std::endl;
  for (auto i : c) {
    std::cout << i << std::endl;
  }

  intersect_sorted_lt(a, b, 3, c);
  std::cout << "c is size " << c.size() << std::endl;
  for (auto i : c) {
    std::cout << i << std::endl;
  }
  return 0;
}
