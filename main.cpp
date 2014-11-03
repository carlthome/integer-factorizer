#include <iostream>
#include <vector>
#include <queue>
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

clock_t start = clock();
int main(int argc, char* argv[])
{
  queue<mpz_class> todo;
  cout << "Factoring number: " << argv[1] << endl;
  todo.push(mpz_class(argv[1]));

  while (!todo.empty())
  {
    auto& next = todo.front();
    auto factors = pollard_rho(next);

    for (auto& fac : factors)
    {
      if (miller_rabin_prime(fac, 50))
      {
        cout << fac << endl;
      }
      else
      {
        todo.push(fac);
      }
    }

    todo.pop();
  }

  cout << "Done!" << endl;
}
