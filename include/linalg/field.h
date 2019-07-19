/*
Implementation of fields

TODO: this could be swapped out for some library
*/

#include <iostream>


template<class Derived>
class AbstractField {
public:
  // addition
  inline Derived operator+(const Derived &b) {
    return static_cast<Derived *>(this)->operator+(b);
  }

  // subtraction
  inline Derived operator-(const Derived &b) {
    return static_cast<Derived *>(this)->operator-(b);
  }

  // additive inverse
  inline Derived operator-() {
    return static_cast<Derived *>(this)->operator-();
  }

  // multiplication
  inline Derived operator*(const Derived &b) {
    return static_cast<Derived *>(this)->operator*(b);
  }

  // division
  inline Derived operator/(const Derived &b) {
    return static_cast<Derived *>(this)->operator/(b);
  }

  // multiplicative inverse
  inline Derived inv() {
    return static_cast<Derived *>(this)->inv();
  }

  // equalivalence
  inline Derived operator==(const Derived &b) {
    return static_cast<Derived *>(this)->operator==(b);
  }

  friend std::ostream& operator<<( std::ostream& os, AbstractField &x) {
    os << *static_cast<Derived*>(&x);
    return os;
  }

  void print() {
    std::cout << *this << std::endl;
  }


};

template<typename IntT, IntT P>
class ModP : public AbstractField<ModP<IntT, P>> {
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

};

// specialized Mod2 implementation
template<typename IntT>
class Mod2 : public AbstractField<Mod2<IntT>> {
private:
  IntT val;

public:
  Mod2(IntT val) : val(val) {}

  Mod2 operator+( const Mod2 &b) const {
    return Mod2(val ^ b.val); // xor
  }

  Mod2 operator-( const Mod2 &b) const {
    return Mod2(val ^ b.val); // xor
  }

  Mod2 operator-() const {
    return Mod2(val); // xor
  }

  Mod2 operator*(const Mod2 &b) const {
    return Mod2(val | b.val); // or
  }

  Mod2 operator/(const Mod2 &b) const {
    if (b.val & 0x1 == 0) {throw "Division by zero!";}
    return Mod2(val); // or
  }

  Mod2 inv() const {
    if (val & 0x1 == 0) {throw "Inversion of zero!";}
    return Mod2(0x1);
  }

  bool operator==(const Mod2 &b) const {
    return (val & 0x1) == (b.val == 0x1);
  }

  friend std::ostream& operator<<( std::ostream& os, const Mod2 &x) {
    os << (x.val & 0x1) << " mod " << 2;
    return os;
  }

};
