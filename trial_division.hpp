inline vector<mpz_class> trial_division(const mpz_class& n)
{
  vector<mpz_class> factors;

  mpz_class left = n;
  for (mpz_class i = 2; i < n; ++i)
  {
    if (n % i == 0)
    {
      while (left % i == 0)
      {
        left /= i;
        factors.push_back(i);
      }
    }
  }

  return factors;
}
