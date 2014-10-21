#pragma once
#include <vector>
#include <climits>
#include <iostream>
using namespace std;

inline vector<mpz_class> sieve_of_eratosthenes(const mpz_class& n)
{
  mpz_class remaining = n;
  vector<mpz_class> factors;

  unsigned int limit = pow(2, 26);
  mpz_class tmp;
  mpz_root(tmp.get_mpz_t(), n.get_mpz_t(), 2);
  tmp++;
  if (tmp.fits_uint_p() && tmp < limit)
  {
    limit = tmp.get_ui();
    cout << "limit now: " << limit << endl;
  }

  vector<bool> is_prime (limit, true);
  
  is_prime[0] = false;
  is_prime[1] = false;
  for (int i = 2; i < limit; i++)
  {
    if (is_prime[i])
    {
      while (remaining % i == 0)
      {
        remaining /= i;
        factors.push_back(i);
      }

      for (int j = i*i; j < limit; j += i)
      {
        is_prime[j] = false;
      }
    }
  }

  // Don't forget the last factor if one exists.
  if (remaining > 1) factors.push_back(remaining);

  return factors;
}
