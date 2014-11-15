#pragma once
using namespace std;

template <class Num>
inline Num powm(Num a, Num b, const Num& m)
{
  Num result = 1;
  while (b > 0)
  {
    if (b % 2 == 1)
    {
      result = result * a % m;
    }
    // divide by two
    b >>= 1;

    a = a * a % m;
  }

  return result;
}

// TODO: remove this once the other stuff works
template <>
inline mpz_class powm(mpz_class a, mpz_class b, const mpz_class& m)
{
  mpz_class result;
  mpz_powm(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), m.get_mpz_t());
  return result;
}
