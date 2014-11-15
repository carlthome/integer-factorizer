#pragma once
#include <utility>
#include "legendre_symbol.hpp"
using namespace std;

// solves the congruence x^2 = n (mod p)
inline pair<unsigned int, unsigned int> tonelli_shanks(unsigned int n, unsigned int p)
{
  if (p == 2)
  {
    return {n, n};
  }

  unsigned int S = 0;
  unsigned int Q = p - 1;

  while (Q % 2 == 0)
  {
    Q /= 2;
    S++;
  }

  unsigned int z = 2;

  while (legendre_symbol(z, p) != -1)
  {
    z++;
  }

  unsigned int c = powm(z, Q, p);
  unsigned int R = powm(n, (Q + 1) / 2, p);
  unsigned int t = powm(n, Q, p);
  unsigned int M = S;

  while (t % p != 1)
  {
    unsigned int i = 1;
    while (powm(t, (unsigned int)pow(2, i), p) != 1)
    {
      i++;
    }

    unsigned int b = powm(c, (unsigned int)pow(2, M - i - 1), p);
    R = R * b % p;
    t = t * b * b % p;
    c = b * b % p;
    M = i;
  }

  return {R, p - R};
}
