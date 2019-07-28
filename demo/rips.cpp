#include <iostream>
#include <filtration/flag.h>
#include <filtration/rips.h>
#include <vector>
#include <random>
#include <chrono>

#define TE tedge<float, size_t>

int main() {

  std::default_random_engine generator;
  std::uniform_real_distribution<float> distribution(0.0,1.0);

  size_t n = 100;
  std::vector<float> v(n);
  for (size_t i = 0; i < n; i++) {
      v[i] = distribution(generator);
  }

  std::vector<size_t> edges;
  std::vector<float> t;
  auto start = std::chrono::high_resolution_clock::now();
  rips_edges(v, edges, t);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "microseconds: " << duration.count() << std::endl;

  // std::cout << "entering Flag construction" << std::endl;
  start = std::chrono::high_resolution_clock::now();
  auto F1 = FlagFiltration(edges, t, 0.0f, v.size(), 2);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "microseconds: " << duration.count() << std::endl;

  // std::cout << "exiting Flag construction" << std::endl;

  std::cout << F1.ncells() << std::endl;
  std::cout << F1.maxdim() << std::endl;

  return 0;
}
