#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>
#include "sieve_of_eratosthenes.hpp"
#include "trial_division.hpp"
using namespace std;

typedef unsigned long long factor; //TODO Use GMP.

clock_t start = clock();
int main(int argc, char* argv[])
{
  if (argc != 2) cerr << "Missing argument." << endl;
  string arg1 = argv[1];
  factor n = strtol(arg1.c_str(), NULL, 10); // TODO Use GMP.
  auto factors = sieve_of_eratosthenes<factor>(n);
  cout << "Composite number: " << n << endl;
  cout << "No. of factors found: " << factors.size() << endl;
  cout << "Runtime: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
  cout << "Factors:" << endl;
  for (auto f : factors) cout << f << endl;
}
