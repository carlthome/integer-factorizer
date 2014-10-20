#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>
using namespace std;

#include "sieve_of_eratosthenes.hpp"

typedef unsigned long long factor; //TODO Use GMP.

vector<factor> trial_division(const factor& n)
{
  vector<factor> factors;

  factor left = n;
  for (factor i = 2; i < n; ++i)
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

clock_t start = clock();
int main()
{ 
  factor n = 100;
  auto factors = sieve_of_eratosthenes<factor>(n);
  cout << "Composite number: " << n << endl;
  cout << "No. of factors found: " << factors.size() << endl;
  cout << "Runtime: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
  cout << "Factors:" << endl;
  for (auto f : factors) cout << f << endl;
}
