#pragma once
#include <vector>
#include <cmath>
using namespace std;

inline mpz_class gcd_iter(const mpz_class& a, const mpz_class& b)
{
  mpz_class t, u, v;
  u = a;
  v = b;
  auto zero = mpz_class(0);

  while (v != zero)
  {
    t = u;
    u = v;
    v = t % v;
  }

  return abs(u);
}
