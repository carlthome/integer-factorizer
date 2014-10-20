#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>
using namespace std;

typedef unsigned long long factor; //TODO Use GMP.

vector<factor> sieve_of_eratosthenes(const factor& n)
{
  auto remaining = n;
  vector<factor> factors;
  vector<bool> is_prime (n, true);
  
  is_prime[0] = false;
  is_prime[1] = false;
  for (int i = 2; i <= sqrt(n); i++)
  {
    if (is_prime[i])
    {
      while (remaining % i == 0)
      {
        remaining /= i;
        factors.push_back(i);
      }

      for (int j = i*i; j <= n; j += i)
      {
        is_prime[j] = false;
      }
    }
  }

  // Don't forget the last factor if one exists.
  if (remaining > 1) factors.push_back(remaining);

  return factors;
}

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
  auto factors = sieve_of_eratosthenes(n);
  cout << "Composite number: " << n << endl;
  cout << "No. of factors found: " << factors.size() << endl;
  cout << "Runtime: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
  cout << "Factors:" << endl;
  for (auto f : factors) cout << f << endl;
}
