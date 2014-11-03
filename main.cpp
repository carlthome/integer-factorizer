#include <iostream>
#include <vector>
#include <map>
#include <functional>
#include <utility>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <mpirxx.h>
#include "trial_division.hpp"
#include "sieve_of_eratosthenes.hpp"
#include "fermat_factorization.hpp"
#include "pollard_rho.hpp"
#include "utils.hpp"
using namespace std;
typedef mpz_class factor;
map<string, function<vector<factor>(factor)>> factorization_algorithms = {
  {"Trial Division", trial_division},
  {"Sieve of Eratosthenes", sieve_of_eratosthenes},
  {"Fermat Factorization", fermat_factorization},
  {"Pollard Rho", pollard_rho},
};

clock_t start = clock();
int main(int argc, char* argv[])
{
start:
  string arg1;
  if (argc < 2)
  {
    cout << "Enter number: ";
    cin >> arg1;
    cout << endl << "Entered " << arg1 << endl;
  }
  else
  {
    arg1 = argv[1];
  }

  const factor n = mpz_class(arg1, 10);

  auto go = [](factor n, string name, function<vector<factor>(factor)> factorization_algorithm)
  {
    vector<factor> factors = factorization_algorithm(n);
    cout << name << ':' << endl;
    cout << "  - Composite number: " << n << endl;
    cout << "  - No. of factors found: " << factors.size() << endl;
    cout << "  - Runtime: " << (double)(clock() - start) / CLOCKS_PER_SEC << " seconds" << endl;
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
    for (auto p : factorization_algorithms)
    {
      go(n, p.first, p.second);
    }
  }

  string answer;
  cout << "Done?" << endl;
  cin >> answer;
  if (answer == "y") return 0;
  else goto start;
}
