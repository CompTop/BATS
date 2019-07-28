// utilities for simplices
#include <functional> // for hash



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
