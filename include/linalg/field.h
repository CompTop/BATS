/*
Implementation of fields

TODO: this could be swapped out for some library
*/

#include <iostream>

template<typename IntT, IntT P>
class ModP {
private:
  IntT val;

public:

  ModP(IntT val) : val((val + P) % P) {}

  ModP operator+( const ModP &b) const {
    return ModP(val + b.val);
  }

  ModP operator-( const ModP &b) const  {
    return ModP(val - b.val + P);
  }

  ModP operator-() const {
    return ModP(P - val);
  }

  ModP operator*( const ModP &b) const {
    return ModP(val * b.val);
  }

  inline bool operator==( const ModP &b) const {
    return val == b.val;
  }

  inline bool operator<( const ModP &b) const {
    return false;
  }

  ModP inv() const {
    IntT b = 1;
    IntT c = val;
    while (c % P != 1) {
      b++;
      c += val;
    }
    return ModP(b);
  }

  ModP operator/( const ModP &b) const {
    ModP binv = b.inv();
    return (val * binv.val);
  }

  bool iszero() const {
    // std::cout << val << " == 0:" << (val == 0) << std::endl;
    return val == 0;
  }

  friend std::ostream& operator<<( std::ostream& os, const ModP &x) {
    os << x.val << " mod " << P;
    return os;
  }

  void print() {
    std::cout << val << " mod " << P << std::endl;
  }
};
