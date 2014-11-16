#pragma once
#include <vector>
#include <utility>
#include <algorithm>
#include "utils.hpp"
#include "tonelli_shanks.hpp"
#include "gf2.hpp"
#include "trial_division.hpp"
#include <iostream>
using namespace std;

namespace
{
  const unsigned int CHUNK_SIZE = 655360;

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

inline vector<mpz_class> quadratic_sieve(const mpz_class& n_)
{
  auto n = n_;
  auto factors = sieve_of_eratosthenes_factorization(n);
  n = factors.back();
  factors.pop_back();
  unsigned int factorized_i = factors.size();
  cout << "n is now " << n << endl;

  // ~magic number~, needs to be different for different n
  auto B = optimal_B(n);
  // sieving range
  auto M = B * B * B;

  auto primes = sieve_of_eratosthenes(B);

  // TODO: try this when it works
  //primes.insert(primes.begin(), -1);

  cout << "calculating factor base from " << primes.size() << " primes" << endl;
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

  cout << "finding sieving indices" << endl;
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

  auto Y = vector<float>(CHUNK_SIZE);

  auto X = vector<unsigned int>();
  auto smooth = vector<vector<unsigned int>>();
  unsigned int start = 0;
  unsigned int end = CHUNK_SIZE;
  // we want at least 20 more smooth numbers than "needed",
  // to make sure that we get a dependancy that factorizes n
  while (smooth.size() < (factor_base.size() + 20))
  {
    cout << "sieving, chunk " << start / CHUNK_SIZE << ", ";
    cout << smooth.size() << " of " << factor_base.size() + 20 << endl;
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
  cout << "transferring to matrix" << endl;

  auto matrix = gf2(smooth.size(), factor_base.size());
  // for each smooth number
  for (unsigned int s = 0; s < smooth.size(); s++)
  {
    // and for each factor in that smooth number
    for (unsigned int p = 0; p < smooth[s].size(); p++)
    {
      matrix.add_bit(s, smooth[s][p]);
    }
  }

  cout << "doing gauss elimination" << endl;

  auto dependencies = matrix.fast_gauss();

  cout << "finding squares" << endl;
  for (unsigned int d = 0; d < dependencies.size(); d++)
  {
    mpz_class R = 1;
    for (unsigned int i = 0; i < dependencies[d].size(); i++)
    {
      R *= X[dependencies[d][i]];
    }
    R = sqrt(R);

    mpz_class b = 1;
    for (unsigned int i = 0; i < smooth[dependencies[d][0]].size(); i++)
    {
      b *= smooth[dependencies[d][0]][i];
    }
    auto gcd = gcd_iter(R - b, n);
    if (gcd != 1 && gcd != n)
    {
      factors.push_back(gcd);
      //break;
    }
  }

  for (; factorized_i < factors.size(); factorized_i++)
  {
    n /= factors[factorized_i];
  }
  if (n != 1)
  {
    factors.push_back(n);
  }

  return factors;
}
