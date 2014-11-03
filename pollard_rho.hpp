#pragma once
#include <vector>
#include "utils.hpp"
using namespace std;

inline vector<mpz_class> pollard_rho(const mpz_class& n)
{
  vector<mpz_class> factors;

  auto g = [=] (const mpz_class& x) -> mpz_class { return (x * x + 1) % n; };
  
  auto one = mpz_class(1);
  auto two = mpz_class(2);
  auto x_fixed = two;
  auto cycle_size = two;
  auto x = two;
  auto h = one;

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

  factors.push_back(h);
  if (n != h)
  {
    factors.push_back(n / h);
  }
  return factors;
}
