#include <iostream>
#include <vector>
#include <functional>
#include <utility>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <gmpxx.h>
#include "trial_division.hpp"
#include "sieve_of_eratosthenes.hpp"
#include "sieve_of_atkin.hpp"
using namespace std;

clock_t start = clock();
int main(int argc, char* argv[])
{
  if (argc != 2) cerr << "Missing argument." << endl;
  string arg1 = argv[1];
  mpz_class n(arg1.c_str()); // TODO Use GMP.
  
  vector<pair<string, function<vector<mpz_class>(mpz_class)>>> factorization_algorithms = {
    {"Trial Division", trial_division},
    {"Sieve of Eratosthenes", sieve_of_eratosthenes},
    {"Sieve of Atkin", sieve_of_atkin<},
  };

  for (auto p : factorization_algorithms)
  {
    const auto& name = p.first;
    const auto& factorization_algorithm = p.second;
    cout << name << ':' << endl;

    vector<mpz_class> factors = factorization_algorithm(n);

    cout << "  - Composite number: " << n << endl;
    cout << "  - No. of factors found: " << factors.size() << endl;
    cout << "  - Runtime: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
    cout << "  - Factors:";
    for (auto f : factors) cout << ' ' << f;
    cout << endl << endl;
  }
  
}
