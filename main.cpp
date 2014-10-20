#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>
using namespace std;

#include "sieve_of_eratosthenes.hpp"
#include "trial_division.hpp"

typedef unsigned long long factor; //TODO Use GMP.

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
