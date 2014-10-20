#pragma once
#include <vector>
#include <cmath>
using namespace std;

template <typename T>
inline vector<T> trial_division(const T& n)
{
  vector<T> factors;

  T left = n;
  for (T i = 2; i < n; ++i)
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
