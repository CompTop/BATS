// barcode utilities
#include <limits>
#include <iostream>

namespace bats {
  
// template over filtration type
template <typename T>
class BarcodePair {
public:
  T birth;
  T death;

  BarcodePair()    : birth(T(0)), death(T(0)) {}
  BarcodePair(T b) : birth(b), death(std::numeric_limits<T>::infinity()) {}
  BarcodePair(T b, T d) : birth(b), death(d) {}

  void print() {
    std::cout << "( " << birth << " , " << death << " )" << std::endl;
  }
};

} // namespace bats
