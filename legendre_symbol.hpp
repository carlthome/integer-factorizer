#pragma once
#include "powm.hpp"
using namespace std;

template <class Num>
inline int legendre_symbol(const Num& a, const Num& p)
{
  auto t = powm<Num>(a, (p - 1) / 2, p);

  if (t > 1)
  {
  	return -1;
  }
  return t;
}

// TODO: remove this once the other stuff works
template <>
inline int legendre_symbol(const mpz_class& a, const mpz_class& p)
{
	return mpz_legendre(a.get_mpz_t(), p.get_mpz_t());
}
