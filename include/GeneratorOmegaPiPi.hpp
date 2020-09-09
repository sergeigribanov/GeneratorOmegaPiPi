#ifndef __GENERATOR_OMEGAPIPI_HPP__
#define __GENERATOR_OMEGAPIPI_HPP__
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <MEUtils/CFourVector.hpp>
#include <MEUtils/MElParticle.hpp>
#include <RandomStar/RandomStar.hpp>
#include <EventOmegaPiPi/EventOmegaPiPi.hpp>

class GeneratorOmegaPiPi {
public:
  GeneratorOmegaPiPi(const std::string&, int, double = 2);
  virtual ~GeneratorOmegaPiPi();
  std::string getPath() const;
  int getN() const;
  void setPath(const std::string&);
  void setN(int);
  void generate(const CFourVector&, const CFourVector&);
private:
  static const int _nMajorant;
  static const int _nInterpMajorant;
  double getMajorant(double) const;
  double calcMajorant(double);
  void interpolateMajorant();
  void genOmegaPiPi(const CFourVector&);
  std::string _path;
  int _n;
  double _maxCMEnergy;
  std::vector<MElParticle*> _particles;
  RandomStar _star;
  gsl_spline* _spline;
  gsl_interp_accel* _acc;
  double* _energies;
  double* _majorants;
  double _m;
};

#endif
