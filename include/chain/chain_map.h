// Chain map class
#include <vector>
#include <cstddef>


template <typename R, class TM>
class ChainMap<R> {
/**
Chain map over a ring R

stores matrices to encode maps
**/
private:
  size_t maxdim;
  std::vector<TM> map;
public:

  ChainMap(size_t maxdim) : maxdim(maxdim) {
    map.reserve(maxdim + 1);
  }

  // set map in dimension k
  void set_map(size_t k, TM &F) {
    while (map.size() < k+1) {
      map.push_back(TM());
    }
    map[k] = F;
  }

};
