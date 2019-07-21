#pragma once
/*
Some operations on sorted containers
*/
#include <vector>
#include <iostream>

/*
sets c = intersection(a, b)
over-writes c
TODO: template over container type as well
*/
template <typename T>
void intersect_sorted(const std::vector<T> &a, const std::vector<T> &b, std::vector<T> &c) {
  c.clear();
  auto ia = a.cbegin();
  auto ib = b.cbegin();
  while (ia < a.cend() && ib < b.cend()) {
    if (*ia < *ib) {
      // std::cout << "a small" << std::endl;
      ++ia;
    } else if (*ib < *ia) {
      // std::cout << "b small" << std::endl;
      ++ib;
    } else {
      // std::cout << "agree" << std::endl;
      // std::cout << *ia << std::endl;
      // *ia == *ib
      c.emplace_back(*ia);
      ++ia;
      ++ib;
    }
  }
  // for (auto i : c) {
  //   std::cout << i << std::endl;
  // }
}

/*
sets c = intersection(a, b, (-inf, maxval))
over-writes c
TODO: template over container type as well
*/
template <typename T>
void intersect_sorted_lt(const std::vector<T> &a, const std::vector<T> &b, T maxval, std::vector<T> &c) {
  c.clear();
  auto ia = a.cbegin();
  auto ib = b.cbegin();
  while (ia < a.cend() && ib < b.cend()) {
    if (*ia < *ib) {
      ++ia;
      if (!(*ia < maxval)) {break;}
    } else if (*ib < *ia) {
      ++ib;
      if (!(*ib < maxval)) {break;}
    } else {
      // *ia == *ib
      c.emplace_back(*ia);
      ++ia;
      ++ib;
      if (!(*ia < maxval) || !(*ib < maxval)) {break;}
    }
  }
}
