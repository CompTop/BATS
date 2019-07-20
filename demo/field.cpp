#include <linalg/field.h>
#include <iostream>

#define F2 ModP<bool, 2>
#define F3 ModP<int, 3>
// should throw compile-time error because not prime
//#define FF ModP<int, 10>
#define Q Rational<int>

int main() {

  //FF t(1);
  Q t(1,2);
  t.print();
  Q s(1,-3);
  s.print();
  Q r = t + s;
  r.print();
  r = t - s;
  r.print();
  r = t * s;
  r.print();
  r = t/s;
  r.print();

  F2 x(1);
  F2 y(0);
  F2 z(2);
  x.print();
  y.print();
  z.print();

  F3 a(0);
  F3 b(2);
  F3 c(4);
  a.print();
  b.print();
  c.print();

  auto d = c + a;
  d.print();

  d = c - b;
  d.print();

  d = -d;
  d.print();

  auto e = -d.inv();
  e.print();

  e = e * e.inv();
  e.print();

  std::cout << e << std::endl;

  return 0;
}
