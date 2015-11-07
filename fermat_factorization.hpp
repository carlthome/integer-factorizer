inline vector<mpz_class> fermat_factorization(const mpz_class& n)
{
  vector<mpz_class> factors;
  if (mpz_odd_p(n.get_mpz_t()) == 0)
  {
    cerr << "n must be odd to apply fermat's method" << endl << "n: " << n << endl;
    return factors;
  }

  mpz_class a, tmp;
  mpz_sqrt(a.get_mpz_t(), n.get_mpz_t());
  a++;

  mpz_class b2 = a * a - n;
  while (mpz_root(tmp.get_mpz_t(), b2.get_mpz_t(), 2) == 0)
  {
    a++;
    b2 = a * a - n;
  }

  mpz_sqrt(b2.get_mpz_t(), b2.get_mpz_t());
  factors.push_back(a - b2);
  factors.push_back(a + b2);

  return factors;
}
