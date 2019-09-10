#pragma once
/*
Implementation of fields

TODO: this could be swapped out for some library
*/

#include <iostream>
#include <numeric>


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

  // integer equivalence
  inline Derived operator==(const int b) {
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

constexpr bool isprime(unsigned p) {
     if (p < 2) { return false; }
     for (unsigned i = 2; i < (p>>1) + 1; i++) {
         if (p % i == 0) { return false; }
     }
     return true;
}

template<typename IntT, unsigned P>
class ModP : public AbstractField<ModP<IntT, P>> {
private:
  IntT val;
  // Check whether P is prime before compiling
  static_assert(isprime(P), "not prime!");

public:


  ModP() : val(0) {}
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

  inline bool operator==( const int b) const {
    return val == b;
  }

  inline bool operator<( const ModP &b) const {
    return true;
  }

  // ModP operator=(const int &a) {
  //   return ModP(a);
  // }

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
IntT ff_inv(const IntT val) {
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

  ModP() : val(0) {}
  ModP(IntT val) : val(val) {}

  ModP operator+( const ModP &b) const {
    return ModP(val ^ b.val); // xor
  }

  ModP operator-( const ModP &b) const {
    return ModP(val ^ b.val); // xor
  }

  inline ModP& operator-() {
    return *this; // no-op
  }

  ModP operator*(const ModP &b) const {
    return ModP(val | b.val); // or
  }

  ModP operator/(const ModP &b) const {
    if ((b.val & 0x1) == 0) {throw "Division by zero!";}
    return ModP(val); // or
  }

  ModP inv() const {
    if ((val & 0x1) == 0) {throw "Inversion of zero!";}
    return ModP(0x1);
  }

  bool operator==(const ModP &b) const {
    return (val & 0x1) == (b & 0x1);
  }

  inline bool operator==( const int b) const {
    return (val & 0x1) == b;
  }

  ModP operator=(const int &a) {
    return ModP(a);
  }

  bool iszero() const {
    // std::cout << val << " == 0:" << (val == 0) << std::endl;
    return (val & 0x1) == 0;
  }

  inline bool operator<( const ModP &b) const {
    return true;
  }

  friend std::ostream& operator<<( std::ostream& os, const ModP &x) {
    os << (x.val & 0x1) << " mod " << 2;
    return os;
  }

};

// sign function
template <typename T>
inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename IntT>
class Rational : public AbstractField<Rational<IntT>> {
private:
  IntT n; // numerator
  IntT d; // denominator

public:

  Rational(IntT n0, IntT d0) {
    IntT gcd = std::gcd(n0, d0);
    int s = sgn(d0);
    n = s * (n0 / gcd);
    d = s * (d0 / gcd);
  }

  Rational(IntT n) : n(n), d(1) {}

  Rational operator+( const Rational &b) const {
    return Rational(n * b.d + b.n * d, d * b.d);
  }

  Rational operator-( const Rational &b) const  {
    return Rational(n * b.d - b.n * d, d * b.d);
  }

  Rational operator-() const {
    return Rational(-n, d);
  }

  Rational operator*( const Rational &b) const {
    return Rational(n * b.n, d * b.d);
  }

  inline bool operator==( const Rational &b) const {
    return (n == b.n) && (d == b.n);
  }

  inline bool operator==( const int b) const {
    return (n == b) && (d == 1);
  }

  inline bool operator<( const Rational &b) const {
    return (n * b.d) < (b.n * d);
  }

  Rational inv() const {
    if (n == 0) {throw "Inversion of zero!";}
    return Rational(d, n);
  }

  Rational operator/( const Rational &b) const {
    if (b.n == 0) {throw "Division by zero!";}
    return Rational(n * b.d, d * b.n);
  }

  bool iszero() const {
    return n == 0;
  }

  friend std::ostream& operator<<( std::ostream& os, const Rational &x) {
    os << x.n << "//" << x.d;
    return os;
  }

};
