using namespace std;
#include <iostream>
#include <mpirxx.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <set>
#include <stack>
#include <cassert>
typedef long long ll;
typedef mpz_class num;

const ll SIEVE_WINDOW = 100000;
const ll TRIAL_BOUND = 100000000;
vector<ll> primes;

// Fast modular exponentiation.
inline ll pow_mod(ll a, ll b, ll m)
{
  auto r = 1;
  while (b > 0)
  {
    if (b & 1)
      r = r * a % m;
    b >>= 1;
    a = a * a % m;
  }

  return r;
}

inline ll legendre_symbol(ll a, ll p)
{
  auto t = pow_mod(a, (p - 1) / 2, p);
  return t > 1 ? -1 : t;
}

// Solve the congruence x^2 = n (mod p).
inline void tonelli_shanks(ll n, ll p, ll *result)
{
  if (p == 2)
  {
    result[0] = result[1] = n;
    return;
  }

  ll S = 0, Q = p - 1;

  while (Q % 2 == 0)
  {
    Q /= 2;
    ++S;
  }

  ll z = 2;

  while (legendre_symbol(z, p) != -1)
    ++z;

  ll c = pow_mod(z, Q, p);
  ll R = pow_mod(n, (Q + 1) / 2, p);
  ll t = pow_mod(n, Q, p);
  ll M = S;

  while (t % p != 1)
  {
    ll i = 1;
    while (pow_mod(t, pow(2, i), p) != 1)
      ++i;

    ll b = pow_mod(c, pow(2, M - i - 1), p);
    R = R * b % p;
    t = t * b * b % p;
    c = b * b % p;
    M = i;
  }

  result[0] = R;
  result[1] = p - R;
}

// Get the i'th bit in row.
inline ll get_bit(ll i, ll *row)
{
  return (row[i / sizeof(ll)] & (1 << (i % sizeof(ll)))) != 0;
}

// Set the i'th bit in row to 1.
inline void set_bit(ll i, ll *row)
{
  row[i / sizeof(ll)] |= (1 << (i % sizeof(ll)));
}

// Set the i'th bit in row to 0.
inline void unset_bit(ll i, ll *row)
{
  row[i / sizeof(ll)] &= ~(1 << (i % sizeof(ll)));
}

// Toggle the i'th bit in row.
inline void toggle_bit(ll i, ll *row)
{
  row[i / sizeof(ll)] ^= (1 << (i % sizeof(ll)));
}

// A quadratic sieve implementation for integers up to 100 bits. n must be composite.
inline num quadratic_sieve(num& n)
{
  // Polynomial
  auto Q = [] (const num& x, const num& n) -> num
  {
    return (sqrt(n) + x) * (sqrt(n) + x) - n;
  };

  // Set smoothness bound with Torbjörn Granlund's magic formula.
  const float log_n = mpz_sizeinbase(n.get_mpz_t(), 2) * log(2);
  const ll B = ceil(100 + 3 * exp(0.5 * sqrt(log_n * log(log_n))));
  cerr << "    Set smoothness bound to " << B << "." << endl;

  // Generate the factor base.
  vector<ll> factor_base;
  for (const auto& prime : primes)
  {
    if (prime >= B) break;
    if (mpz_legendre(n.get_mpz_t(), num(prime).get_mpz_t()) == 1)
      factor_base.push_back(prime);
  }
  cerr << "    Generated factor base of size " << factor_base.size() << "." << endl;

  // Calculate sieve index (where to start the sieve) for each factor base number.
  ll **fb_indexes = new ll*[2];
  fb_indexes[0] = new ll[factor_base.size()];
  fb_indexes[1] = new ll[factor_base.size()];
  for (auto p = 0; p < factor_base.size(); ++p)
  {
    // At what indexes do we start this sieve? Solve the congruence x^2 = n (mod p) to find out.
    // Results in two solutions, so we do two sieve iterations for each prime in the factor base.
    ll idxs[2];
    num temp = n % num(factor_base[p]);
    tonelli_shanks(temp.get_ui(), factor_base[p], idxs);

    temp = idxs[0] - sqrt(n);
    temp = ((temp % factor_base[p]) + factor_base[p]) % factor_base[p];
    fb_indexes[0][p] = temp.get_ui();

    temp = idxs[1] - sqrt(n);
    temp = ((temp % factor_base[p]) + factor_base[p]) % factor_base[p];
    fb_indexes[1][p] = temp.get_ui();
  }
  cerr << "    Determined sieve indices per factor base number." << endl;

  // Sieve new chunks until we have enough smooth numbers.
  cerr << "    Sieving:" << endl;
  float last_estimate = 0;
  ll next_estimate = 1;
  ll min_x = 0, max_x = SIEVE_WINDOW;
  vector<ll> X;
  vector<float> Y(SIEVE_WINDOW, 0);
  vector<vector<ll>> smooth;
  while (smooth.size() < factor_base.size() + 20)
  {
    // Calculate Y vector of log approximations.
    for (auto i = 1; i < SIEVE_WINDOW; ++i)
    {
      auto x = i + min_x;

      // Only calculate log estimates if neccessary.
      if (next_estimate <= x)
      {
        auto y = Q(x, n);
        last_estimate = mpz_sizeinbase(y.get_mpz_t(), 2); // Number of bits is approximately log_2.
        next_estimate *= 1.8 + 1; //TODO Document properly. The higher t gets, the less the logarithm of Y[t] changes.
      }

      Y[i] = last_estimate;
    }

    // Sieve
    for (auto p = 0; p < factor_base.size(); ++p)
    {
      for (auto t = 0; t < 2; ++t)
      {
        while (fb_indexes[t][p] < max_x)
        {
          Y[fb_indexes[t][p] - min_x] -= log(factor_base[p]) / log(2);
          fb_indexes[t][p] += factor_base[p];
        }
        if (factor_base[p] == 2) break; // p = 2 only has one modular root
      }
    }

    // Factor all values whose logarithms were reduced to approximately zero using trial division.
    float threshold = log(factor_base.back()) / log(2);
    for (auto i = 0; i < SIEVE_WINDOW; ++i)
    {
      auto x = i + min_x;

      if (fabs(Y[i]) < threshold)
      {
        auto y = Q(x, n);
        smooth.push_back(vector<ll>());

        for (auto p = 0; p < factor_base.size(); ++p)
        {
          while (mpz_divisible_ui_p(y.get_mpz_t(), factor_base[p]))
          {
            mpz_divexact_ui(y.get_mpz_t(), y.get_mpz_t(), factor_base[p]);
            smooth.back().push_back(p);
          }
        }

        // Keep smooth number.
        if (y == 1) X.push_back(i + min_x);
        else smooth.pop_back(); // Not smooth. Remove!
      }
    }

    min_x += SIEVE_WINDOW;
    max_x += SIEVE_WINDOW;
    cerr << "      Sieved (smooth numbers found = " << smooth.size() << " out of " << factor_base.size() + 20 << "." << endl;
  }
  cerr << "    Sieving completed." << endl;

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
  cerr << "    Created matrix." << endl;

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
  cerr << "    Performed Gauss elimination." << endl;

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
            cerr << "  Trial dividing prime " << prime << "." << endl;
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
          cerr << "  Using quadratic sieve." << endl;
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