inline vector<mpz_class> sieve_of_eratosthenes_factorization(const mpz_class& n)
{
  mpz_class remaining = n;
  vector<mpz_class> factors;

  unsigned int limit = pow(2, 26);
  mpz_class tmp;
  mpz_sqrt(tmp.get_mpz_t(), n.get_mpz_t());
  tmp++;
  if (tmp.fits_uint_p() && tmp < limit)
  {
    limit = tmp.get_ui();
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

inline vector<int> sieve_of_eratosthenes(unsigned int B)
{
  vector<int> primes;
  vector<bool> is_prime (B, true);

  is_prime[0] = false;
  is_prime[1] = false;
  for (int i = 2; i < B; i++)
  {
    if (is_prime[i])
    {
      primes.push_back(i);

      for (int j = i*i; j < B; j += i)
      {
        is_prime[j] = false;
      }
    }
  }

  return primes;
}
