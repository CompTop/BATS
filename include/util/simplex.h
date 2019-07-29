// utilities for simplices
#include <functional> // for hash
#include <vector>



// implement has function for vectors

struct SimplexHasher
{
  std::size_t operator()(const std::vector<size_t>& k) const
  {
    using std::size_t;
    using std::hash;

    // Compute individual hash values for first,
    // second and third and combine them using XOR
    // and bit shifting:
    size_t ret = 0;
    for (auto i : k) {
      // ret = ret ^ (i);
      ret *= i;
      // ret += ((i+1) << 6);
    }
    return ret;
  }
};


class SimplexContainer {
private:
  std::vector<size_t> data;
  size_t k; // number of vertices in this simplex size
public:

  // empty constructor - not available
  // constructor with dimension specified
  SimplexContainer(size_t d) : k(d+1) {}

  // emplace_back with vector
  void emplace_back(std::vector<size_t> &s) {
    if (data.size() != k) {throw "unexpected simplex size!";}
    for (auto i : s) {
      data.emplace_back(i);
    }
  }

  // size - number of simplices
  size_t size() const {
    return data.size() / k;
  }
  // operator[] - return vector
  // reserve - reserve space for a fixed number of simplices
  // dim - dimension
  inline size_t dim() const {
    return k-1;
  }
  // function to return iterator over vertices of simplex i
};
