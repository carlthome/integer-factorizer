#include "miller_rabin.hpp"

void go(const mpz_class& n, vector<mpz_class>& factors)
{
  if (n == one) return;
  if (miller_rabin_prime(n, 64))
  {
    factors.push_back(n);
  }
  else
  {
    static auto g = [=] (const mpz_class& x) -> mpz_class
    {
      return (x * x + 1) % n;
    };

    auto x = two;
    auto h = one;
    auto x_fixed = two;
    auto cycle_size = two;
    while (h == one)
    {
      auto count = one;

      while (count <= cycle_size && h == one)
      {
        x = g(x);
        count += 1;
        h = gcd_iter(x - x_fixed, n);
      }

      if (h != one)
      {
        break;
      }

      cycle_size *= 2;
      x_fixed = x;
    }

    go(h, factors);
    go(n / h, factors);
  }
}

// Repeatedly apply Pollard's Rho algorithm to factor the input number n.
vector<mpz_class> pollard_rho(const mpz_class& n)
{
  vector<mpz_class> factors;
  go(n, factors);
  return factors;
}
