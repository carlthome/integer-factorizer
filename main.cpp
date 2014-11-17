using namespace std;
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <functional>
#include <utility>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <mpirxx.h>
#include "utils.hpp"
#include "trial_division.hpp"
#include "sieve_of_eratosthenes.hpp"
#include "fermat_factorization.hpp"
#include "pollard_rho.hpp"
#include "quadratic_sieve.hpp"

typedef function<vector<mpz_class>(mpz_class)> factorization_algorithm;
map<string, factorization_algorithm> factorization_algorithms = {
  { "Trial Division", trial_division },
  { "Sieve of Eratosthenes", sieve_of_eratosthenes_factorization },
  { "Fermat Factorization", fermat_factorization },
  { "Pollard Rho", pollard_rho },
  { "Quadratic Sieve", quadratic_sieve },
};

int main(int argc, char* argv[])
{
bool interactive_mode = argc < 2;

start:
  string arg1;
  
  if (interactive_mode)
  {
    cout << "Enter number: ";
    cin >> arg1;
    cout << endl << "Entered: " << arg1 << endl;
  }
  else
  {
    arg1 = argv[1];
  }

  const auto n = mpz_class(arg1, 10);

  auto go = [] (mpz_class n, string name, factorization_algorithm f)
  {
    auto start = clock();
    auto factors = f(n);
    auto runtime = (clock() - start) / (double)CLOCKS_PER_SEC;
    cout << name << ':' << endl;
    cout << "  - Composite number: " << n << endl;
    cout << "  - No. of factors found: " << factors.size() << endl;
    cout << "  - Runtime: " << runtime << " seconds" << endl;
    cout << "  - Factors: ";
    for (auto f : factors) cout << f << ", ";
    cout << endl << endl;
  };

  if (argc > 2)
  {
    for (int i = 2; i < argc; ++i)
    {
      string s = argv[i];
      auto f = factorization_algorithms[s];
      go(n, s, f);
    }
  }
  else
  {
    //TODO Don't factor the number with each algorithm. Instead combine the algorithms to factor the number once. Start by trying pollard rho for a while, then if not factored in time: use quadratic sieve.
    for (auto p : factorization_algorithms)
    {
      go(n, p.first, p.second);
    }
  }

  if (interactive_mode)
  {
    string answer;
    cout << "Done?" << endl;
    cin >> answer;
    if (answer == "y") return 0;
    else goto start;
  }
}
