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
  size_t maxdim = 2;
  std::vector<float> v(n);
  for (size_t i = 0; i < n; i++) {
      v[i] = distribution(generator);
  }

  std::vector<size_t> edges;
  std::vector<float> t;
  auto start = std::chrono::high_resolution_clock::now();
  rips_edges(v, edges, t);
  sort_edges(edges, t);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "microseconds: " << duration.count() << std::endl;
  std::cout << "n edges =" << t.size() << std::endl;

  // std::cout << "entering Flag construction" << std::endl;
  start = std::chrono::high_resolution_clock::now();
  // SimplicialComplex X;
  // Filtration<float, SimplicialComplex> F;
  auto [X, F] = FlagFiltration(edges, t, v.size(), maxdim, 0.0f);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "microseconds: " << duration.count() << std::endl;

  // std::cout << "exiting Flag construction" << std::endl;

  for (size_t k = 0; k < X.maxdim() + 1; k++) {
      std::cout << "cells in dim " << k << " = " << X.ncells(k) << std::endl;
  }
  std::cout << F.maxdim() << std::endl;

  return 0;
}
