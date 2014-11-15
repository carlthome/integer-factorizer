#pragma once
#include <vector>
#include <utility>
#include "utils.hpp"
#include "tonelli_shanks.hpp"
#include <iostream>
using namespace std;

namespace
{
  const unsigned int CHUNK_SIZE = 65536;

  unsigned int optimal_B(const mpz_class& n)
  {
    // reference: http://www.andrew-g-west.com/docs/thesis_slides.pdf
    auto digits = mpz_sizeinbase(n.get_mpz_t(), 10);
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
  auto factors = vector<mpz_class>();
  mpz_class a, tmp;
  mpz_sqrt(a.get_mpz_t(), n.get_mpz_t());
  // round up instead of down
  a++;
  // ~magic number~, needs to be different for different n
  auto B = optimal_B(n);
  // sieving range
  auto M = B * B * B;

  auto primes = sieve_of_eratosthenes(B);

  // TODO: try this when it works
  //primes.insert(primes.begin(), -1);

  auto factor_base = vector<unsigned int>();
  for (unsigned int i = 0; i < primes.size(); i++)
  {
    mpz_class p = primes[i];
    if (mpz_legendre(n.get_mpz_t(), p.get_mpz_t()) == 1)
    {
      factor_base.push_back(primes[i]);
    }
  }

  auto factor_base_indices = vector<vector<unsigned int>>(2, vector<unsigned int>(factor_base.size()));
  auto sqrt_n = sqrt(n);

  // find sieving indices for each factor base number
  for (unsigned int p = 0; p < factor_base.size(); p++)
  {
    mpz_class tmp = n % mpz_class(factor_base[p]);
    unsigned int i0, i1;
    tie(i0, i1) = tonelli_shanks(tmp.get_ui(), factor_base[p]);

    tmp = i0 - sqrt_n;
    tmp = ((tmp % factor_base[p]) + factor_base[p]) % factor_base[p];
    factor_base_indices[0][p] = tmp.get_ui();

    tmp = i1 - sqrt_n;
    tmp = ((tmp % factor_base[p]) + factor_base[p]) % factor_base[p];
    factor_base_indices[1][p] = tmp.get_ui();
  }

  auto Y = vector<float>(M * 2);

  auto X = vector<unsigned int>();
  auto smooth = vector<vector<unsigned int>>();
  unsigned int start = 0;
  unsigned int end = CHUNK_SIZE;
  // we want at least 20 more smooth numbers than "needed",
  // to make sure that we get a dependancy that factorizes n
  while (smooth.size() < (factor_base.size() + 20))
  {
    for (unsigned int t = 1; t < CHUNK_SIZE; t++)
    {
      // y = (sqrt(n) + x)^2 - n
      mpz_class y = (sqrt_n + t + start) * (sqrt_n + t + start) - n;

      // we estimate the 2 logarithm by counting the number of bits
      Y[t] = mpz_sizeinbase(y.get_mpz_t(), 2);
    }

    // time to sieve!
    for (unsigned int p = 0; p < factor_base.size(); p++)
    {
      float lg = log(factor_base[p]) / log(2);

      for (unsigned int t = 0; t < 2; t++)
      {
        while (factor_base_indices[t][p] < end)
        {
          Y[factor_base_indices[t][p] - start] -= lg;
          factor_base_indices[t][p] += factor_base[p];
        }

        // if p = 2 we only have one modular root
        if (factor_base[p] == 2)
        {
          break;
        }
      }
    }

    float threshold = log(factor_base.back()) / log(2);

    for (unsigned int i = 0; i < CHUNK_SIZE; i++)
    {
      if (abs(Y[i]) < threshold)
      {
        mpz_class y = (sqrt_n + i + start) * (sqrt_n + i + start) - n;
        smooth.emplace_back();

        for (unsigned int p = 0; p < factor_base.size(); p++)
        {
          while (mpz_divisible_ui_p(y.get_mpz_t(), factor_base[p]))
          {
            mpz_divexact_ui(y.get_mpz_t(), y.get_mpz_t(), factor_base[p]);
            smooth.back().push_back(p);
          }
        }

        if (y == 1)
        {
          // smooth! woop woop

          X.push_back(i + start);
          if (smooth.size() >= (factor_base.size() + 20))
          {
            break;
          }
        }
        else
        {
          smooth.pop_back();
        }
      }
    }

    start += CHUNK_SIZE;
    end += CHUNK_SIZE;
  }

  return factors;
}
