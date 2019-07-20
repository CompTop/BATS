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

template<typename IntT, unsigned P>
class ModP : public AbstractField<ModP<IntT, P>> {
private:
  IntT val;
  // Check whether P is prime before compiling
  // only checks for small numbers. doesn't cover all primes
  static_assert(!(P == 1), "not prime!");
  static_assert(!((P > 2) && (P % 2 == 0)), "not prime!");
  static_assert(!((P > 3) && (P % 3 == 0)), "not prime!");
  static_assert(!((P > 2) && (P % 2 == 0)), "not prime!");
  static_assert(!((P > 5) && (P % 5 == 0)), "not prime!");
  static_assert(!((P > 7) && (P % 7 == 0)), "not prime!");
  static_assert(!((P > 11) && (P % 11 == 0)), "not prime!");

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

  ModP inv() const;

  ModP operator/( const ModP &b) const {
    if (b.val == 0) {throw "Division by zero!";}
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

// helper function for inversion
template<typename IntT, unsigned P>
IntT ff_inv(IntT val) {
  if (val == 0) {throw "Inversion of zero!";}
  IntT b = 1;
  IntT c = val;
  while (c % P != 1) {
    b++;
    c += val;
  }
  return b;
}


// implementation of field inversion
template<typename IntT, unsigned P>
ModP<IntT, P> ModP<IntT,P>::inv() const {
  return ModP<IntT, P>(ff_inv<IntT, P>(val));
}

// specialized Mod2 implementation
template<typename IntT>
class ModP<IntT, 2> : public AbstractField<ModP<IntT, 2>> {
private:
  IntT val;

public:
  ModP(IntT val) : val(val) {}

  ModP operator+( const ModP &b) const {
    return ModP(val ^ b.val); // xor
  }

  ModP operator-( const ModP &b) const {
    return ModP(val ^ b.val); // xor
  }

  ModP operator-() const {
    return ModP(val); // xor
  }

  ModP operator*(const ModP &b) const {
    return ModP(val | b.val); // or
  }

  ModP operator/(const ModP &b) const {
    if (b.val & 0x1 == 0) {throw "Division by zero!";}
    return ModP(val); // or
  }

  ModP inv() const {
    if (val & 0x1 == 0) {throw "Inversion of zero!";}
    return ModP(0x1);
  }

  bool operator==(const ModP &b) const {
    return (val & 0x1) == (b.val == 0x1);
  }

  friend std::ostream& operator<<( std::ostream& os, const ModP &x) {
    os << (x.val & 0x1) << " mod " << 2;
    return os;
  }

};
