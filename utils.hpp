#pragma once
#include <vector>
#include <cmath>
#include <time.h>
using namespace std;

namespace
{
  // random int from min to max inclusive
  mpz_class rand_range(const mpz_class& min, const mpz_class& max)
  {
    static gmp_randclass r (gmp_randinit_default);
    r.seed(time(NULL));
    return min + (r.get_z_range(max - min + 1));
  }
}

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
