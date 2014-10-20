#include <iostream>
#include <vector>
#include <functional>
#include <utility>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "trial_division.hpp"
#include "sieve_of_eratosthenes.hpp"
#include "sieve_of_atkin.hpp"
using namespace std;

typedef unsigned long long factor; //TODO Use GMP.

clock_t start = clock();
int main(int argc, char* argv[])
{
  if (argc != 2) cerr << "Missing argument." << endl;
  string arg1 = argv[1];
  factor n = strtol(arg1.c_str(), NULL, 10); // TODO Use GMP.
  
  vector<pair<string, function<vector<factor>(factor)>>> factorization_algorithms = {
    {"Trial Division", trial_division<factor>},
    {"Sieve of Eratosthenes", sieve_of_eratosthenes<factor>},
    {"Sieve of Atkin", sieve_of_atkin<factor>},
  };

  for (auto p : factorization_algorithms)
  {
    const auto& name = p.first;
    const auto& factorization_algorithm = p.second;
    cout << name << ':' << endl;

    vector<factor> factors = factorization_algorithm(n);

    cout << "  - Composite number: " << n << endl;
    cout << "  - No. of factors found: " << factors.size() << endl;
    cout << "  - Runtime: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
    cout << "  - Factors:";
    for (auto f : factors) cout << ' ' << f;
    cout << endl << endl;
  }
  
}
