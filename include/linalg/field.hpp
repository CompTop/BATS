#pragma once
/*
Implementation of fields

TODO: this could be swapped out for some library
*/

#include <iostream>
#include <numeric>
#include <string>
#include <stdexcept>


template<class Derived>
class AbstractField {
public:
  // addition
  inline Derived operator+(const Derived &b) {
    return static_cast<Derived *>(this)->operator+(b);
  }

  // +=, -=, *=, /=

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

  std::string str() {
    std::ostringstream oss;
    oss << *this;
    return oss.str();
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

  inline IntT to_int() const { return val; }

  ModP operator+( const ModP &b) const {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    return ModP(val + b.val);
  }

  ModP& operator+=( const ModP &b) {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    val = (val + b.val) % P;
    return *this;
  }

  ModP operator-( const ModP &b) const  {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    return ModP(val - b.val + P);
  }

  ModP& operator-=( const ModP &b) {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    val = (val - b.val + P) % P;
    return *this;
  }

  ModP operator-() const {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    return ModP(P - val);
  }

  ModP operator*( const ModP &b) const {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    return ModP(val * b.val);
  }

  ModP& operator*=( const ModP &b) {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    val = (val * b.val) % P;
    return *this;
  }

  inline bool operator==( const ModP &b) const {
    return val == b.val;
  }

  inline bool operator!=( const ModP &b) const {
    return val != b.val;
  }

  inline bool operator==( const int b) const {
    return ((b - val + P) % P) == 0;
  }

  inline bool operator!=( const int b) const {
    return ((val - b + P) % P) != 0;
  }

  inline bool operator<( const ModP &b) const {
    return (val % P) < (b.val % P);
  }

  // ModP operator=(const int &a) {
  //   return ModP(a);
  // }

  ModP inv() const;

  ModP operator/( const ModP &b) const {
    if (b.val == 0) {throw std::invalid_argument("division by 0");}
    ModP binv = b.inv();
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    return (val * binv.val);
  }

  ModP& operator/=( const ModP &b) {
    if (b.val == 0) {throw std::invalid_argument("division by 0");}
    ModP binv = b.inv();
    val = (val * binv.val) % P;
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    return *this;
  }

  bool iszero() const {
    // std::cout << val << " == 0:" << (val == 0) << std::endl;
    return val == 0;
  }

  friend std::ostream& operator<<( std::ostream& os, const ModP &x) {
    //os << x.val << " mod " << P;
    os << x.val;
    return os;
  }

};

// helper function for inversion
template<typename IntT, unsigned P>
IntT ff_inv(const IntT val) {
  if (val == 0) {throw std::invalid_argument("inversion of 0");}
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

  inline IntT to_int() const { return val & 0x1; }

  ModP operator+( const ModP &b) const {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    return ModP(val ^ b.val); // xor
  }

  ModP& operator+=( const ModP &b ) {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    val = val ^ b.val; // xor
    return *this;
  }

  ModP operator-( const ModP &b) const {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    return ModP(val ^ b.val); // xor
  }

  ModP& operator-=( const ModP &b ) {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    val = val ^ b.val; // xor
    return *this;
  }

  inline ModP& operator-() {
    return *this; // no-op
  }

  ModP operator*(const ModP &b) const {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    return ModP(val & b.val); // or
  }

  ModP& operator*=(const ModP &b) {
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    val = val & b.val; // or
    return *this;
  }

  ModP operator/(const ModP &b) const {
    if ((b.val & 0x1) == 0) {throw std::runtime_error("Division by zero!");}
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    return ModP(val); // or
  }

  ModP& operator/=(const ModP &b) {
    if ((b.val & 0x1) == 0) {throw std::runtime_error("Division by zero!");}
    #ifdef BATS_OPCOUNT
    bats::global_ops++;
    #endif
    return *this; // no-op
  }

  ModP inv() const {
    if ((val & 0x1) == 0) {throw std::runtime_error("Inversion of zero!");}
    return ModP(0x1);
  }

  inline bool operator==( const ModP &b ) const {
    return (val & 0x1) == (b.val & 0x1);
  }

  inline bool operator!=( const ModP &b ) const {
    return (val & 0x1) != (b.val & 0x1);
  }


  inline bool operator==( const int b) const {
    return (val & 0x1) == (b & 0x1);
  }

  inline bool operator!=( const int b) const {
    return (val & 0x1) != b;
  }

  ModP operator=(const int &a) {
    return ModP(a);
  }

  bool iszero() const {
    // std::cout << val << " == 0:" << (val == 0) << std::endl;
    return (val & 0x1) == 0;
  }

  inline bool operator<( const ModP &b) const {
    return (val & 0x1) < (b.val & 0x1);
  }

  friend std::ostream& operator<<( std::ostream& os, const ModP &x) {
    // os << (x.val & 0x1) << " mod " << 2;
    os << (x.val & 0x1);
    return os;
  }

};

// sign function
template <typename T>
inline T sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


template<typename IntT>
class Rational : public AbstractField<Rational<IntT>> {
private:
  IntT n; // numerator
  IntT d; // denominator

  // reduction of numerator and denominator
  void reduce() {
    IntT gcd = std::gcd(n, d);
    IntT s = sgn(d);
    n = s * (n / gcd);
    d = s * (d / gcd);
    #ifdef WARN_RATIONAL_OVERFLOW
    if constexpr(std::is_same<IntT,int32_t>::value) {
      const int32_t mask = 0x7F00; // ignore sign bit
      auto n2 = n < 0 ? -n : n;  // make numerator positive
      if ( ((n2 |d) & mask) > 0) {
        std::cerr<< "Rational overflow possibility! " << std::endl;
      }
    } else if constexpr(std::is_same<IntT, int64_t>::value) {
      const int64_t mask = 0x7FFF0000; // ignore sign bit
      auto n2 = n < 0 ? -n : n;
      if ( ((n2|d) & mask) != 0) {
        std::cerr<< "Rational overflow possibility! " << std::endl;
      }
    }
    #endif
  }

public:

  Rational() : n(0), d(1) {}

  Rational(IntT n0, IntT d0) : n(n0), d(d0) {
    reduce();
  }

  Rational(IntT n) : n(n), d(1) {}

  inline IntT to_int() const { return n / d; }

  Rational operator+( const Rational &b) const {
    return Rational(n * b.d + b.n * d, d * b.d);
  }

  Rational& operator+=(const Rational &b) {
    n = n * b.d + b.n * d;
    d *= b.d;
    reduce();
    return *this;
  }

  Rational operator-( const Rational &b) const  {
    return Rational(n * b.d - b.n * d, d * b.d);
  }

  Rational& operator-=( const Rational &b) {
    n = n * b.d - b.n * d;
    d *= b.d;
    reduce();
    return *this;
  }

  Rational operator-() const {
    return Rational(-n, d);
  }

  Rational operator*( const Rational &b) const {
    return Rational(n * b.n, d * b.d);
  }

  Rational& operator*=( const Rational &b) {
    n *= b.n;
    d *= b.d;
    reduce();
    return *this;
  }

  inline bool operator==( const Rational &b) const {
    return (n == b.n) && (d == b.d);
  }

  inline bool operator!=( const Rational &b) const {
    return (n != b.n) || (d != b.d);
  }

  inline bool operator==( const int b) const {
    return (n == b) && (d == 1);
  }

  inline bool operator!=( const int b) const {
    return (n != b) || (d != 1);
  }

  inline bool operator<( const Rational &b) const {
    return (n * b.d) < (b.n * d);
  }

  Rational inv() const {
    if (n == 0) {throw std::runtime_error("Inversion of zero!");}
    return Rational(d, n);
  }

  Rational operator/( const Rational &b) const {
    if (b.n == 0) {throw std::runtime_error("Division by zero!");}
    return Rational(n * b.d, d * b.n);
  }

  Rational& operator/=( const Rational &b) {
    if (b.n == 0) {throw std::runtime_error("Division by zero!");}
    n *= b.d;
    d *= b.n;
    reduce();
    return *this;
  }

  bool iszero() const {
    return n == 0;
  }

  friend std::ostream& operator<<( std::ostream& os, const Rational &x) {
    os << x.n << "/" << x.d;
    return os;
  }

};
