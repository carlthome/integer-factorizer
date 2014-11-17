#include "legendre_symbol.hpp"

namespace
{
  inline mpz_class pow_mpz(unsigned long int base, unsigned long int e)
  {
    mpz_class ret;
    mpz_ui_pow_ui(ret.get_mpz_t(), base, e);
    return ret;
  }
}

// solves the congruence x^2 = n (mod p)
inline pair<unsigned long int, unsigned long int> tonelli_shanks(unsigned long int a, unsigned long int p)
{
  if (p == 2)
  {
    return {a, a};
  }
  if (p % 4 == 3)
  {
    auto x = powm(a, (p + 1) / 4, p);
    return {x, a - x};
  }

  unsigned long int s = p - 1;
  unsigned long int e = 0;
  while (s % 2 == 0)
  {
    s /= 2;
    e+= 1;
  }

  unsigned long int n = 2;
  while (legendre_symbol(n, p) != -1)
  {
    n++;
  }

  unsigned long int x = powm(a, (s + 1) / 2, p);
  unsigned long int b = powm(a, s, p);
  unsigned long int g = powm(n, s, p);
  unsigned long int r = e;

  while (true)
  {
    unsigned long int t = b;
    unsigned long int m = 0;
    for (m = 0; m < r; m++)
    {
      if (t == 1)
      {
        break;
      }
      t = powm<unsigned long int>(t, 2, p);
    }

    if (m == 0)
    {
      return {x, a - x};
    }

    unsigned long int gs = powm<unsigned long int>(g, pow(2, r - m - 1), p);
    g = (gs * gs) % p;
    x = (x * gs) % p;
    b = (b * g) % p;
    r = m;
  }
}
