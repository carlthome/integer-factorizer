#pragma once
#include <vector>
#include <climits>
using namespace std;

inline unsigned int powm(unsigned int a, unsigned int b, unsigned int m)
{
  unsigned int result = 1;
  while (b > 0)
  {
  	// fast way to check b mod 2 == 1
  	if (b & 1)
  	{
  		result = result * a % m;
  	}
  	// divide by two
  	b >>= 1;

  	a = a * a % m;
  }

  return result;
}
