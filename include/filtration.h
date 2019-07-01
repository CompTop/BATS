#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include "abstract_complex.h"

// template over
//  TC - complex type
//  TF - filtration type
template <class TC, typename TF>
class Filtration
{
private:
  TC cpx;
  std::vector<std::vector<TF>> val; // filtration value for each cell
  std::vector<std::vector<size_t>> filt_perm; // filtration perumation for each dimension
  std::vector<std::vector<size_t>> inv_perm; // inverse permutation for each dimension
public:
  // initialize with only complex
  Filtration(TC cpx) : cpx(cpx) {}

  Filtration(TC cpx,
    std::vector<std::vector<TF>> val) : cpx(cpx), val(val) {}

  // assume filtration value has been set for each cell.  Find sort permutaion in each dimension
  void sort() {
    filt_perm.resize(cpx.maxdim() + 1);
    inv_perm.resize(cpx.maxdim() + 1);
    // loop over dimension
    for (size_t dim = 0; dim < cpx.maxdim() + 1; dim++) {
      // std::cout << "sorting dim = " << dim << std::endl;
      size_t ncells_dim = cpx.ncells(dim);
      // std::cout << "ncells dim = " << ncells_dim << std::endl;
      filt_perm[dim] = std::vector<size_t>(ncells_dim);
      // std::cout << "created vector" << std::endl;
      std::iota(filt_perm[dim].begin(), filt_perm[dim].end(), 0);
      // std::cout << "filled iota" << std::endl;
      std::sort(filt_perm[dim].begin(), filt_perm[dim].end(),
        [&](const int& a, const int& b) {
          return (val[dim][a] < val[dim][b]);
        }
      );
      // std::cout << "sorted" << std::endl;
      // now filt_perm[dim] holds the sort permutation on cells for the filtration
      inv_perm[dim] = std::vector<size_t>(ncells_dim);
      for (size_t i = 0; i < ncells_dim; i++) {
        inv_perm[dim][filt_perm[dim][i]] = i;
      }
      // std::cout << "inv sorted" << std::endl;
      // now inv_perm[dim] holds the inverse permutation
    }
    return;
  }

  void print(size_t dim) {
    for (auto i : filt_perm[dim]) {
      std::cout << val[dim][i] << " | ";
      cpx.print_cell(dim, i);
    }
  }

  void print() {
    for (size_t dim = 0; dim < cpx.maxdim() + 1; dim ++) {
      std::cout << "dim " << dim << " : " << cpx.ncells(dim) << " cells" << std::endl;
      print(dim);
    }
  }

  // create boundary in dimension dim

};


// create a lower-star filtraiton on a cell complex of type TC
template <class TC, typename TF>
Filtration<TC,TF> LowerStarFiltration(TC cpx, std::vector<TF> f) {

  // assert f.size() == cpx.ncells(0);

  // initialize filtration
  std::vector<std::vector<TF>> val(cpx.maxdim()+1);
  val[0] = f;
  // std::cout << "extending filtration..." << std::endl;
  for (size_t dim = 1; dim < cpx.maxdim() + 1; dim++) {
    size_t ncells_dim = cpx.ncells(dim);
    val[dim] = std::vector<TF>(ncells_dim);
    for (size_t i = 0; i < ncells_dim; i++) {
      std::vector<size_t> skel0 = cpx.skeleton0(dim, i);
      size_t imax = *std::max_element(skel0.begin(),skel0.end(),
        [&f](int i1, int i2){return f[i1]<f[i2];}
      );
      val[dim][i] = f[imax];
    }
  }

  // std::cout << "creating object" << std::endl;
  auto F = Filtration<TC, TF>(cpx, val);
  // std::cout << "sorting" << std::endl;
  F.sort();
  return F;
}
