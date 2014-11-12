#pragma once
#include <vector>
#include "utils.hpp"
#include <iostream>
using namespace std;

namespace
{
  unsigned int optimal_B(const mpz_class& n)
  {
    // reference: http://www.andrew-g-west.com/docs/thesis_slides.pdf
    auto digits = mpz_sizeinbase(n.get_mpz_t(), 10)
    if (digits <= 35)
      return 17023;
    else if (digits <= 40)
      return 29723;
    else if (digits <= 45)
      return 70575;
    else if (digits <= 50)
      return 137134;
    else if (digits <= 55)
      return 263382;
    else if (digits <= 60)
      return 568835;
    else if (digits <= 65)
      return 807450;
    else if (digits <= 70)
      return 1582110;
    else
      return 60000000; // definitely black magic by now, tweak this
  }
}

inline vector<mpz_class> quadratic_sieve(const mpz_class& n)
{
  vector<mpz_class> factors;
  mpz_class a, tmp;
  mpz_sqrt(a.get_mpz_t(), n.get_mpz_t());
  // round up instead of down
  a++;
  // ~magic number~, needs to be different for different n
  B = optimal_B(n);
  


  return factors;
}
