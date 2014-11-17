const auto zero = mpz_class(0);
const auto one = mpz_class(1);
const auto two = mpz_class(2);
const auto three = mpz_class(3);

inline mpz_class gcd_iter(const mpz_class& a, const mpz_class& b)
{
  mpz_class t, u, v;
  u = a;
  v = b;

  while (v != zero)
  {
    t = u;
    u = v;
    v = t % v;
  }

  return abs(u);
}
