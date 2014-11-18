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
using namespace std;
typedef long long ll;
typedef mpz_class num;
typedef function<num(num, num)> polynomial;
const ll SIEVE_WINDOW = 100000;
const ll TRIAL_BOUND = 100000000;
vector<ll> primes;

inline num square(const num& n)
{
  return n * n;
}

inline double log(const num& n)
{
  //TODO * log(2)
  return mpz_sizeinbase(n.get_mpz_t(), 2); // Number of bits is approximately log_2.
}

inline num pow_mod(const num& base, const num& exponent, const num& modulo)
{
  num r;
  mpz_powm(r.get_mpz_t(), base.get_mpz_t(), exponent.get_mpz_t(), modulo.get_mpz_t());
  return r;
}

inline num legendre_symbol(const num& a, const num& p)
{
  num exp = (p - 1) / 2;
  num r = pow_mod(a, exp, p);
  return (r > 1) ? num(-1) : r;
}

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

  while (legendre_symbol(z, p) != -1) ++z;

  num c = pow_mod(z, Q, p);
  num R = pow_mod(n, (Q + 1) / 2, p);
  num t = pow_mod(n, Q, p);
  num M = S;

  while (t % p != 1)
  {
    num i = 1;
    while (pow_mod(t, num(pow(2, i.get_ui())), p) != 1) ++i; //TODO Can i become too large?
    num exp = M - i - 1;
    num b = pow_mod(c, num(pow(2, exp.get_ui())), p);
    R = R * b % p;
    t = t * b * b % p;
    c = b * b % p;
    M = i;
  }

  return{R, p - R};
}

// Get the i:th bit in row.
inline ll get_bit(ll i, ll *row)
{
  return (row[i / sizeof(ll)] & (1 << (i % sizeof(ll)))) != 0;
}

// Set the i:th bit in row.
inline void set_bit(ll i, ll *row)
{
  row[i / sizeof(ll)] |= (1 << (i % sizeof(ll)));
}

// Unset the i:th bit in row.
inline void unset_bit(ll i, ll *row)
{
  row[i / sizeof(ll)] &= ~(1 << (i % sizeof(ll)));
}

// Flip the i'th bit in row.
inline void toggle_bit(ll i, ll *row)
{
  row[i / sizeof(ll)] ^= (1 << (i % sizeof(ll)));
}

// A quadratic sieve implementation for integers up to 100 bits. n must be composite.
inline num quadratic_sieve(const num& n)
{
  // Set smoothness bound with Torbjörn Granlund's magic formula.
  const num B = ceil(3 * exp(0.5 * sqrt(log(n) * log(log(n)))));
  cerr << "    Smoothness bound: " << B << endl;

  // Generate the factor base.
  vector<num> factor_base;
  for (const auto& prime : primes)
  {
    if (prime >= B) break;
    if (mpz_legendre(n.get_mpz_t(), num(prime).get_mpz_t()) == 1)
      factor_base.push_back(num(prime));
  }

  cerr << "    Factor base: " << factor_base.size() << endl;

  // Calculate sieve index (where to start the sieve) for each factor base number.
  vector<pair<size_t, size_t>> indexes(0);
  for (auto p = 0; p < factor_base.size(); ++p)
  {
    // Solve the congruence x^2 = n (mod p) to find out where to start sieving.
    auto r = tonelli_shanks(n % factor_base[p], factor_base[p]);
    num idx1 = (((r.first - sqrt(n)) % factor_base[p]) + factor_base[p]) % factor_base[p];
    num idx2 = (((r.second - sqrt(n)) % factor_base[p]) + factor_base[p]) % factor_base[p];
    indexes.push_back({idx1.get_ui(), idx2.get_ui()});
  }
  cerr << "    Determined sieve indices per factor base number." << endl;

  // Sieve new chunks until we have enough smooth numbers.
  cerr << "    Sieving:" << endl;
  float last_estimate = 0;
  ll next_estimate = 1;
  size_t min_x = 0, max_x = SIEVE_WINDOW;
  vector<ll> X;
  vector<double> Y(SIEVE_WINDOW, 0);
  vector<vector<ll>> smooth; //TODO Replace with parity bits directly.
  while (smooth.size() < factor_base.size() + 20)
  {
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
      for (auto i = 1; i < SIEVE_WINDOW; ++i)
      {
        auto x = i + min_x;

        // Only calculate log estimates if necessary.
        if (next_estimate <= x)
        {
          auto y = Q(n, x);
          last_estimate = log(y);
          next_estimate *= 1.8 + 1; //TODO Document properly. The higher t gets, the less the logarithm of Y[t] changes.
        }

        Y[i] = last_estimate;
      }

      // Sieve
      for (auto p = 0; p < factor_base.size(); ++p)
      {
        auto& idx1 = indexes[p].first;
        auto& idx2 = indexes[p].second;
        const auto val = log(factor_base[p].get_ui()) / log(2);
        while (idx1 < max_x) Y[idx1 - min_x] -= val, idx1 += factor_base[p].get_ui();
        while (idx2 < max_x) Y[idx2 - min_x] -= val, idx2 += factor_base[p].get_ui();
        //TODO if (factor_base[p] == 2) break; // p has only one modular root
      }

      // Factor all values whose logarithms were reduced to approximately zero using trial division.
      for (auto i = 0; i < SIEVE_WINDOW; ++i)
      {
        auto x = i + min_x;
        if (abs(Y[i]) < 0.1*log(n)) //TODO Set proper tolerance.
        {
          //TODO Avoid duplicates from other polynomials.
          auto y = Q(n, x);

          vector<ll> factorization;
          for (auto p = 0; p < factor_base.size(); ++p)
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

  // Create parity matrix of exponent vectors by going through each factor in each smooth number.
  auto **matrix = new ll*[factor_base.size()];
  auto row_words = (smooth.size() + sizeof(ll)) / sizeof(ll);
  for (auto i = 0; i < factor_base.size(); ++i)
  {
    matrix[i] = new ll[row_words];
    memset(matrix[i], 0, row_words * sizeof(ll));
  }
  for (auto s = 0; s < smooth.size(); ++s)
    for (auto p = 0; p < smooth[s].size(); ++p)
      toggle_bit(s, matrix[smooth[s][p]]);
  cerr << "    Created matrix" << endl;

  // Gauss elimination.
  ll i = 0, j = 0;
  while (i < factor_base.size() && j < (smooth.size() + 1))
  {
    ll maxi = i;

    // Find pivot element.
    for (auto k = i + 1; k < factor_base.size(); ++k)
    {
      if (get_bit(j, matrix[k]) == 1)
      {
        maxi = k;
        break;
      }
    }
    if (get_bit(j, matrix[maxi]) == 1)
    {
      swap(matrix[i], matrix[maxi]);

      for (auto u = i + 1; u < factor_base.size(); ++u)
      {
        if (get_bit(j, matrix[u]) == 1)
        {
          for (auto w = 0; w < row_words; ++w)
            matrix[u][w] ^= matrix[i][w];
        }
      }
      ++i;
    }
    ++j;
  }
  cerr << "    Performed Gauss elimination" << endl;

  // A copy of matrix that we'll perform back-substitution on.
  num a, b;
  ll **back_matrix = new ll*[factor_base.size()];
  for (auto i = 0; i < factor_base.size(); ++i)
    back_matrix[i] = new ll[row_words];

  ll *x = new ll[smooth.size()];
  ll *combination = new ll[factor_base.size()];

  // Loop until a != +/- b (mod n) to find a non-trivial factor.
  while (a % n == b % n || a % n == (-b) % n + n)
  {
    // Copy the gauss eliminated matrix.
    for (auto i = 0; i < factor_base.size(); ++i) memcpy(back_matrix[i], matrix[i], row_words * sizeof(ll));

    // Clear the x vector.
    memset(x, 0, smooth.size() * sizeof(ll));

    // Perform back-substitution.
    ll i = factor_base.size() - 1;
    while (i >= 0)
    {
      // Count non-zero elements in current row.
      ll count = 0, current = -1;
      for (auto c = 0; c < smooth.size(); ++c)
      {
        count += get_bit(c, back_matrix[i]);
        current = get_bit(c, back_matrix[i]) ? c : current;
      }

      // Empty row, advance to next.
      if (count == 0)
      {
        --i;
        continue;
      }

      // Underdermined, pick x[current] randomly.
      x[current] = count > 1 ? rand() % 2 : get_bit(smooth.size(), back_matrix[i]);

      for (auto u = 0; u <= i; ++u)
      {
        if (get_bit(current, back_matrix[u]) == 1)
        {
          if (x[current] == 1) toggle_bit(smooth.size(), back_matrix[u]);
          unset_bit(current, back_matrix[u]);
        }
      }

      if (count == 1) --i;
    }

    // Combine factor base to find square.
    a = 1, b = 1;
    memset(combination, 0, sizeof(ll) * factor_base.size());
    for (auto i = 0; i < smooth.size(); ++i)
    {
      if (x[i] == 1)
      {
        for (auto p = 0; p < smooth[i].size(); ++p) ++combination[smooth[i][p]];
        b *= (X[i] + sqrt(n));
      }
    }

    for (auto p = 0; p < factor_base.size(); ++p)
    {
      for (auto i = 0; i < (combination[p] / 2); ++i)
        a *= factor_base[p];
    }
  }
  b -= a;

  num factor;
  mpz_gcd(factor.get_mpz_t(), b.get_mpz_t(), n.get_mpz_t());

  cerr << "  Found factor " << factor << " with quadratic sieve." << endl;

  return factor;
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
          if (mpz_divisible_p(factor.get_mpz_t(), num(prime).get_mpz_t()))
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
          for (auto n = 2; n < sqrt(factor); ++n)
          {
            num root, rem;
            mpz_rootrem(root.get_mpz_t(), rem.get_mpz_t(), factor.get_mpz_t(), n);
            if (rem == 0) for (auto i = 0; i < n; ++i) factors.push(root);
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