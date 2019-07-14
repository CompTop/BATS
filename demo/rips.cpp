#include <iostream>
#include <filtration/rips.h>

#define TE tedge<float, size_t>

int main() {

  auto a = make_edge(0, 1, 1.0);
  std::cout << a << std::endl;

  return 0;
}
