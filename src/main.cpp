#include <cstdlib>
#include "GeneratorOmegaPiPi.hpp"

int main(int argc, char* argv[]) {
  int n = std::atoi(argv[1]);
  double energy = std::atof(argv[2]);
  GeneratorOmegaPiPi q(argv[3], 100000, 4.);  
  q.generate(CFourVector(0., 0., -0.5 * energy, 0.5 * energy),
	     CFourVector(0., 0., 0.5 * energy, 0.5 * energy));
  return 0;
}
