#include <vector>
#include <cmath>

#ifndef INTEGERFACTORIZER_TRIALDIVISION
#define INTEGERFACTORIZER_TRIALDIVISION

template<class Num>
inline vector<Num> trial_division(const Num& n)
{
  vector<Num> factors;

  Num left = n;
  for (Num i = 2; i < n; ++i)
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


#endif