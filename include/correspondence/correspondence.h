#pragma once
#include <map>
#include <vector>

// subclass for correspondences (multiset)
// subclass for function (single point)
// subclass for inclusion (injective map)
// note that inclusion can be defined just on 0-cells, then extended

class AbstractCorrespondence {
  /*
  correspondence between points
  */
public:
  AbstractCorrespondence() {}

  // add y to correspondence set of x
  void add(const size_t x, const size_t y);
};

class Correspondence : public AbstractCorrespondence {
private:
  size_t nx;
  size_t ny;
  std::vector<std::vector<size_t>> vtx_map;
public:
  //Correspondence() : {}

  Correspondence(size_t nx, size_t ny) : nx(nx), ny(ny) {
    vtx_map.resize(nx);
  }

  // operator []

  // add y to correspondence set of x
  void add(const size_t x, const size_t y) {
    vtx_map[x].push_back(y);
  }

  // get(x)
  std::vector<size_t> get(const size_t x) {
    return vtx_map[x];
  }

  // transpose
  Correspondence transpose();

  bool is_surjective() const;
  bool is_function() const;
  bool is_inclusion() const;
};

class Function : public AbstractCorrespondence {
private:
  size_t nx;
  size_t ny;
  std::vector<size_t> vtx_map;
public:
  Function(size_t nx, size_t ny) : nx(nx), ny(ny) {
    vtx_map.resize(nx);
  }

  // add y to correspondence set of x
  void add(const size_t x, const size_t y) {
    vtx_map[x] = y;
  }

  // get(x)
  size_t get(const size_t x) {
    return vtx_map[x];
  }

  Correspondence transpose();

  bool is_surjective() const;
  inline bool is_function() const {return true;}
  bool is_inclusion() const ;

};

class Inclusion : private AbstractCorrespondence {
private:
  size_t nx;
  size_t ny;
  std::vector<size_t> vtx_map;

public:
  Inclusion(size_t nx, size_t ny) : nx(nx), ny(ny) {
    vtx_map.resize(nx);
  }

  // add y to correspondence set of x
  void add(const size_t x, const size_t y) {
    vtx_map[x] = y;
  }

  // get(x)
  size_t get(const size_t x) {
    return vtx_map[x];
  }

  bool is_surjective() const;
  inline bool is_function() const {return true;}
  inline bool is_inclusion() const {return true;}


};
