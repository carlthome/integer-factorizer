#pragma once
#include <vector>
#include <climits>
using namespace std;

inline vector<mpz_class> sieve_of_eratosthenes(const mpz_class& n)
{
  auto remaining = n;
  vector<mpz_class> factors;
  vector<bool> is_prime (n.get_ui(), true);
  
  is_prime[0] = false;
  is_prime[1] = false;
  for (int i = 2; i <= sqrt(n); i++)
  {
    if (is_prime[i])
    {
      while (remaining % i == 0)
      {
        remaining /= i;
        factors.push_back(i);
      }

      for (int j = i*i; j <= n; j += i)
      {
        is_prime[j] = false;
      }
    }
  }

  // Don't forget the last factor if one exists.
  if (remaining > 1) factors.push_back(remaining);

  return factors;
}