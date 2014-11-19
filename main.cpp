#include <iostream>
#include <mpirxx.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <set>
#include <stack>
#include <utility>
#include <functional>
#include <cassert>
#include "gf2.hpp"
using namespace std;
typedef mpz_class num;
typedef function<num(num, num)> polynomial;
const long SIEVE_WINDOW = 1000000;
const long TRIAL_BOUND = 100000000;
vector<num> primes;

// Number of bits in a is approximately log_2(a). 
inline float log(const num& a) { return mpz_sizeinbase(a.get_mpz_t(), 2); }

// Solve the congruence x^2 = n (mod p).
inline pair<num, num> tonelli_shanks(const num& n, const num& p)
{
  if (p == 2) return{n, n};

  num S = 0, Q = p - 1;

  while (Q % 2 == 0)
  {
    Q /= 2;
    ++S;
  }

  num z = 2;

  while (mpz_legendre(z.get_mpz_t(), p.get_mpz_t()) != -1) ++z;

  static auto pow = [](const num& base, const num& exponent)
  {
    assert(mpz_fits_ui_p(exponent.get_mpz_t()));
    num r;
    mpz_pow_ui(r.get_mpz_t(), base.get_mpz_t(), exponent.get_ui());
    return r;
  };

  static auto pow_mod = [](const num& base, const num& exponent, const num& modulo)
  {
    num r;
    mpz_powm(r.get_mpz_t(), base.get_mpz_t(), exponent.get_mpz_t(), modulo.get_mpz_t());
    return r;
  };

  num c = pow_mod(z, Q, p);
  num R = pow_mod(n, (Q + 1) / 2, p);
  num t = pow_mod(n, Q, p);
  num M = S;

  while (t % p != 1)
  {
    num i = 1;
    while (pow_mod(t, pow(2, i), p) != 1) ++i;
    num b = pow_mod(c, pow(2, M - i - 1), p);
    R = R * b % p;
    t = t * b * b % p;
    c = b * b % p;
    M = i;
  }

  return{R, p - R};
}

// Factorize n.
inline num quadratic_sieve(const num& n)
{
  // Set smoothness bound with Torbjörn Granlund's magic formula.
  const num B = ceil(3 * exp(0.5 * sqrt(log(n) * log(log(n)))));
  cerr << "    Smoothness bound: " << B << endl;

  // Generate the factor base.
  vector<num> factor_base;
  for (const auto& prime : primes)
  {
    if (B <= prime) break;
    if (mpz_legendre(n.get_mpz_t(), prime.get_mpz_t()) == 1)
    {
      factor_base.push_back(prime);
    }
  }
  cerr << "    Factor base: " << factor_base.size() << endl;

  // Calculate sieve index (where to start the sieve) for each factor base number.
  vector<pair<size_t, size_t>> indexes(0);
  for (unsigned long p = 0; p < factor_base.size(); ++p)
  {
    // Solve the congruence x^2 = n (mod p) to find out where to start sieving.
    auto r = tonelli_shanks(n % factor_base[p], factor_base[p]);
    num idx1 = (((r.first - sqrt(n)) % factor_base[p]) + factor_base[p]) % factor_base[p];
    num idx2 = (((r.second - sqrt(n)) % factor_base[p]) + factor_base[p]) % factor_base[p];
    assert(mpz_fits_ui_p(idx1.get_mpz_t()));
    assert(mpz_fits_ui_p(idx2.get_mpz_t()));
    indexes.push_back({idx1.get_ui(), idx2.get_ui()});
  }
  cerr << "    Determined sieve indices per factor base number." << endl;

  // Sieve new chunks until we have enough smooth numbers.
  cerr << "    Sieving:" << endl;
  float last_estimate = 0;
  size_t next_estimate = 1;
  size_t min_x = 0, max_x = SIEVE_WINDOW;
  vector<size_t> X;
  vector<float> Y(SIEVE_WINDOW, 0);
  vector<vector<long>> smooth; //TODO Replace with parity bits directly.
  while (smooth.size() < factor_base.size() + 20)
  {
    static auto square = [](const num& n) { return n * n; };
    static vector<polynomial> polynomials = {
      [](const num& n, const num& x) -> num { return square((sqrt(1 * n) + x)) - n; },
      [](const num& n, const num& x) -> num { return square((sqrt(2 * n) + x)) - n; },
      [](const num& n, const num& x) -> num { return square((sqrt(3 * n) + x)) - n; },
      [](const num& n, const num& x) -> num { return square((sqrt(4 * n) + x)) - n; },
      [](const num& n, const num& x) -> num { return square((sqrt(5 * n) + x)) - n; },
      [](const num& n, const num& x) -> num { return square((sqrt(6 * n) + x)) - n; },
      [](const num& n, const num& x) -> num { return square((sqrt(7 * n) + x)) - n; },
      [](const num& n, const num& x) -> num { return square((sqrt(8 * n) + x)) - n; },
      [](const num& n, const num& x) -> num { return square((sqrt(9 * n) + x)) - n; },
    };
    for (const auto& Q : polynomials)
    {
      // Calculate Y vector of log approximations.
      for (unsigned long i = 1; i < SIEVE_WINDOW; ++i)
      {
        auto x = i + min_x;

        // Only calculate log estimates if necessary.
        if (next_estimate <= x)
        {
          auto y = Q(n, x);
          last_estimate = log(y);
          next_estimate *= 2;
        }

        Y[i] = last_estimate;
      }

      // Sieve
      auto i = 0;
      for (const auto& p : factor_base)
      {
        auto& idx1 = indexes[i].first;
        auto& idx2 = indexes[i].second;
        auto sieve = [&](size_t& idx)
        {
          while (idx < max_x)
          {
            Y[idx - min_x] -= log(p);
            assert(mpz_fits_ui_p(p.get_mpz_t()));
            idx += p.get_ui();
          }
        };
        sieve(idx1);
        sieve(idx2);
        ++i;
      }

      // Factor all values whose logarithms were reduced to approximately zero using trial division.
      for (unsigned long i = 0; i < SIEVE_WINDOW; ++i)
      {
        auto x = i + min_x;
        if (abs(Y[i]) < 0.1*log(n)) //TODO Set proper tolerance.
        {
          //TODO Avoid duplicates from other polynomials.
          auto y = Q(n, x);

          vector<long> factorization;
          for (unsigned long p = 0; p < factor_base.size(); ++p)
          {
            while (y % factor_base[p] == 0)
            {
              y /= factor_base[p];
              factorization.push_back(p);
            }
          }

          // Keep factorization if smooth.
          if (y == 1)
          {
            X.push_back(i + min_x);
            smooth.push_back(factorization);
          }
        }
      }
    }

    min_x += SIEVE_WINDOW;
    max_x += SIEVE_WINDOW;
    cerr << "      Smooth numbers found = " << smooth.size() << " out of " << factor_base.size() + 20 << "" << endl;
  }
  cerr << "    Sieving completed" << endl;

  //TODO Replace below with fast gauss.

  // Create parity matrix of exponent vectors by going through each factor in each smooth number.
  auto matrix = gf2(smooth.size(), factor_base.size());

  for (unsigned long s = 0; s < smooth.size(); ++s)
    for (unsigned long p = 0; p < smooth[s].size(); ++p)
      matrix.add_bit(s, smooth[s][p]);
  cerr << "    Created matrix" << endl;

  auto dependencies = matrix.fast_gauss();

  cerr << "    Performed Gauss elimination" << endl;

  auto factors = set<num>();
  for (unsigned long d = 0; d < dependencies.size(); d++)
  {
    num R = 1;
    for (unsigned long i = 0; i < dependencies[d].size(); i++)
    {
      R *= X[dependencies[d][i]];
    }
    R = sqrt(R);

    // TODO: I think the code might be wrong from this point onwards
    // it doesn't always find dependencies even though there are easy
    // factors available
    num b = 1;
    for (unsigned int i = 0; i < smooth[dependencies[d][0]].size(); i++)
    {
      b *= smooth[dependencies[d][0]][i];
    }

    num gcd = R - b;
    mpz_gcd(gcd.get_mpz_t(), gcd.get_mpz_t(), n.get_mpz_t());
    if (gcd != 1 && gcd != n)
    {
     factors.insert(gcd);
    }
  }
  if (factors.size() == 0)
    cerr << "FAILED on " << n << endl;

  cerr << "  Found factor " << (*factors.rbegin()) << " with quadratic sieve." << endl;

  // TODO: return all factors
  // right now it returns the biggest (last in the ordered set) factor
  return factors.size() > 0 ? (*factors.rbegin()) : 0;
}

int main()
{
  // Sieve out some small primes.
  bool *is_prime = new bool[TRIAL_BOUND];
  memset(is_prime, true, TRIAL_BOUND);
  is_prime[0] = is_prime[1] = false;
  for (int i = 2; i <= sqrt(TRIAL_BOUND + 1); i++)
  {
    if (!is_prime[i]) continue;
    primes.push_back(i);
    for (int j = i + i; j <= TRIAL_BOUND; j += i) is_prime[j] = false;
  }

  // Factor numbers.
  num number;
  while (cin >> number)
  {
    num n = number, product(1);
    stack<num> factors;
    factors.push(n);

factorize:
    while (!factors.empty())
    {
      num factor = factors.top(); factors.pop();
      cerr << endl << "Factoring " << factor << ":" << endl;

      // Print probable primes.
      if (mpz_probab_prime_p(factor.get_mpz_t(), 64))
      {
        cerr << "  Number was prime." << endl;
        cout << factor << endl;
        product *= factor;
        goto factorize;
      }
      else
      {
        //TODO Test Pollard's rho algorithm for a limited time.

        // Trial divide small primes.
        for (const auto& prime : primes)
        {
          if (mpz_divisible_p(factor.get_mpz_t(), prime.get_mpz_t()))
          {
            cerr << "  Trial division found prime " << prime << "." << endl;
            factors.push(prime);
            factors.push(factor / prime);
            goto factorize;
          }
        }

        // Avoid perfect powers by dividing and continuing per base factor.
        if (mpz_perfect_power_p(factor.get_mpz_t()))
        {
          cerr << "  Avoiding perfect powers." << endl;
          for (unsigned long n = 2; n < sqrt(factor); ++n)
          {
            num root, rem;
            mpz_rootrem(root.get_mpz_t(), rem.get_mpz_t(), factor.get_mpz_t(), n);
            if (rem == 0) for (unsigned long i = 0; i < n; ++i) factors.push(root);
          }
        }
        else
        {
          cerr << "  Using quadratic sieve:" << endl;
          num f = quadratic_sieve(factor);
          factors.push(f);
          factors.push(factor / f);
        }
      }
    }
    assert(number == product);
    cout << endl;
  }

  return 0;
}
