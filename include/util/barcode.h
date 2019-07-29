// barcode utilities
#include <limits>

// template over filtration type
template <typename T>
class BarcodePair {
public:
  T birth;
  T death;

  BarcodePair()    : birth(T(0)), death(T(0)) {}
  BarcodePair(T b) : birth(b), death(std::numeric_limits<T>::infinity()) {}
};
