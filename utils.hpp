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

// return true if it PROBABLY is prime
inline bool miller_rabin_prime(const mpz_class & n, unsigned int k)
{
  if (mpz_odd_p(n.get_mpz_t()) == 0)
  {
    return false;
  }

  mpz_class m, d;
  d = m = n - 1;
  mpz_class zero = mpz_class(0);
  mpz_class one = mpz_class(1);
  unsigned int s = 0;

  while (mpz_even_p(d.get_mpz_t()))
  {
    s++;
    d /= 2;
  }

  for (unsigned int i = 0; i < k; i++)
  {
    mpz_class a = rand_range(2, n - 2);
    mpz_class x;
    // x = a^d mod n
    mpz_powm(x.get_mpz_t(), a.get_mpz_t(), d.get_mpz_t(), n.get_mpz_t());

    if (x == one || x == m)
    {
      continue;
    }

    bool maybe_prime = false;
    for (unsigned int j = 0; j < s - 1; j++)
    {
      mpz_powm_ui(x.get_mpz_t(), x.get_mpz_t(), 2, n.get_mpz_t());
      if (x == one)
      {
        return false;
      }
      if (x == m)
      {
        maybe_prime = true;
        break;
      }
    }

    if (maybe_prime)
    {
      continue;
    }
    return false;
  }

  return true;
}
