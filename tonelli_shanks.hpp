#pragma once
#include <utility>
#include "legendre_symbol.hpp"
using namespace std;

// solves the congruence x^2 = n (mod p)
inline pair<unsigned int, unsigned int> tonelli_shanks(unsigned int n_, unsigned int p_)
{
  if (p_ == 2)
  {
    return {n_, n_};
  }
  unsigned long int p = p_;
  unsigned long int n = n_;

  unsigned long int S = 0;
  unsigned long int Q = p - 1;

  while (Q % 2 == 0)
  {
    Q /= 2;
    S++;
  }

  unsigned long int z = 2;

  while (legendre_symbol(z, p) != -1)
  {
    z++;
  }

  unsigned long int c = powm(z, Q, p);
  unsigned long int R = powm(n, (Q + 1) / 2, p);
  unsigned long int t = powm(n, Q, p);
  unsigned long int M = S;

  while (t % p != 1)
  {
    unsigned long int i = 1;
    while (powm<unsigned long int>(t, pow(2, i), p) != 1)
    {
      i++;
    }

    unsigned long int b = powm<unsigned long int>(c, pow(2, M - i - 1), p);
    R = R * b % p;
    t = t * b * b % p;
    c = b * b % p;
    M = i;
  }

  return {R, p - R};
}
