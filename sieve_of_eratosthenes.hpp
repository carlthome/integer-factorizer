#include <vector>
#include <cmath>

#ifndef INTEGERFACTORIZER_SIEVEOFERATOSTHENES
#define INTEGERFACTORIZER_SIEVEOFERATOSTHENES

template<class Num>
inline vector<Num> sieve_of_eratosthenes(const Num& n)
{
  auto remaining = n;
  vector<Num> factors;
  vector<bool> is_prime (n, true);
  
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


#endif