#pragma once
#include <vector>
#include <climits>
using namespace std;

inline vector<mpz_class> sieve_of_atkin(const mpz_class& n)
{
  mpz_class x, y, i;
  
  auto remaining = n;
  vector<mpz_class> factors;
  auto factorize = [&remaining, &factors] (mpz_class factor) {
    while (remaining % factor == 0)
    {
      remaining /= factor;
      factors.push_back(factor);
    }
  };  
  
  vector<bool> is_prime(n.get_ui(), false);

  mpz_class seq[] = {2, 4};
  mpz_class k1 = 0, k = 0;
 
  x = 1;
  y = 0;
  while (x < sqrt(n / 4) + 1)
  {
    i = 0;
    k1 = 4 * x * x;
    y = 1;
    if (x % 3 == 0)
    {
      while (true)
      {
        k = k1 + y * y;
        if (k >= n) break;
        is_prime[k.get_ui()] = !is_prime[k.get_ui()];
        mpz_class tmp = ++i & 1;
        y += seq[tmp.get_ui()];
      }
    }
    else
    {
      while (true)
      {
        k = k1 + y * y;
        if (k >= n) break;
        is_prime[k.get_ui()] = !is_prime[k.get_ui()];
        y += 2;
      }
    }
    x++;
  }
 
  x = 1;
  y = 0;
  while (x < sqrt(n / 3) + 1)
  {
    i = 1;
    k1 = 3 * x * x;
    y = 2;
    while (true)
    {
      k = k1 + y * y;
      if (k >= n) break;
      is_prime[k.get_ui()] = !is_prime[k.get_ui()];
        mpz_class tmp = ++i & 1;
        y += seq[tmp.get_ui()];
    }
    x += 2;
  }
 
  x = 1;
  y = 0;
  while (x < sqrt(n))
  {
    k1 = 3 * x * x;
    if ((x & 1) == 0)
    {
      y = 1;
      i = 0;
    } else
    {
      y = 2;
      i = 1;
    }
    while (y < x)
    {
      k = k1 - y * y;
      if (k < n) is_prime[k.get_ui()] = !is_prime[k.get_ui()];
        mpz_class tmp = ++i & 1;
        y += seq[tmp.get_ui()];
    }
    x++;
  }
 
  is_prime[2] = true;
  factorize(2);
  is_prime[3] = true;
  factorize(3);
  for (int j = 5; j <= sqrt(n) + 1; ++j)
  {
    if (is_prime[j])
    {
      factorize(j);
      mpz_class n2 = n * n;
      for (k = n2; k < n; k += n2)
      {
        is_prime[k.get_ui()] = false;
      }
    }
  }

   // Don't forget the last factor if one exists.
  if (remaining > 1) factors.push_back(remaining);

  return factors;
}
