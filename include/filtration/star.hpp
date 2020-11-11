
namespace bats {


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

} // namespace bats
