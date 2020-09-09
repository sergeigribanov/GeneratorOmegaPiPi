#include <iostream>
#include <limits>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <MElCalc/MElCalc.hpp>
#include "GeneratorOmegaPiPi.hpp"

const int GeneratorOmegaPiPi::_nMajorant = 10000;
const int GeneratorOmegaPiPi::_nInterpMajorant = 200;

GeneratorOmegaPiPi::GeneratorOmegaPiPi(const std::string& path, int n, double maxCMEnergy) :
  _path(path), _n(n), _maxCMEnergy(maxCMEnergy),
  _particles({
    new MElParticle(-211),
    new MElParticle(-211),
    new MElParticle(211),
    new MElParticle(211),
    new MElParticle(111)}),
  _star(_particles),
  _spline(gsl_spline_alloc(gsl_interp_cspline, _nInterpMajorant)),
  _acc(gsl_interp_accel_alloc()),
  _energies(new double[_nInterpMajorant]),
  _majorants(new double[_nInterpMajorant]) {
  _m = 0;
  for (const auto& el : _particles) {
    _m += el->getMass();
  }
  interpolateMajorant();
}

GeneratorOmegaPiPi::~GeneratorOmegaPiPi() {
  gsl_spline_free(_spline);
  gsl_interp_accel_free(_acc);
  delete [] _energies;
  delete [] _majorants;
}

std::string GeneratorOmegaPiPi::getPath() const {
  return _path;
}

int GeneratorOmegaPiPi::getN() const {
  return _n;
}

void GeneratorOmegaPiPi::setPath(const std::string& path) {
  _path = path;
}

void GeneratorOmegaPiPi::setN(int n) {
  _n = n;
}

void GeneratorOmegaPiPi::generate(const CFourVector& eMiP, const CFourVector& ePlP) {
  CFourVector p = eMiP + ePlP;
  // TODO : make boost of p to CM framework
  // TODO : compare CM energy with the maximum CM energy (exception)
  auto file = TFile::Open(_path.c_str(), "recreate");
  file->cd();
  auto event = new EventOmegaPiPi();
  auto tree = new TTree("tree", "");
  tree->Branch("ev", &event, 16000, 99);
  TLorentzVector piMi0;
  TLorentzVector piMi1;
  TLorentzVector piPl0;
  TLorentzVector piPl1;
  TLorentzVector pi0;
  int l = _n / 100;
  for (int i = 0; i < _n; ++i) {
    if (i % l == 0) {
      std::cout << "event : " << i << std::endl;
    }
    genOmegaPiPi(p);
    for (int j = 0; j < 4; ++j) {
      piMi0[j] = _star.getMomenta()[0][j].real();
      piMi1[j] = _star.getMomenta()[1][j].real();
      piPl0[j] = _star.getMomenta()[2][j].real();
      piPl1[j] = _star.getMomenta()[3][j].real();
      pi0[j] = _star.getMomenta()[4][j].real();
    }
    event->setPiMi0(piMi0);
    event->setPiMi1(piMi1);
    event->setPiPl0(piPl0);
    event->setPiPl1(piPl1);
    event->setPi0(pi0);
    tree->Fill();
  }
  file->Write();
  file->Close();
  delete event;
}

double GeneratorOmegaPiPi::calcMajorant(double cmEnergy) {
  if (cmEnergy <= _m) {
    return 0;
  }
  CFourVector p(0., 0., 0., cmEnergy);  
  _star.setInitialMomentum(p);
  double result = 0;
  double tmp;
  for (int i = 0; i < _nMajorant; ++i) {
    _star.generate();
    tmp = MElCalc::getOmega2PiMEl2
      (std::make_pair(_star.getMomenta()[0], _star.getMomenta()[1]),
       std::make_pair(_star.getMomenta()[2], _star.getMomenta()[3]),
       _star.getMomenta()[4]) * 
      _star.getPhaseSpace();
    if (tmp > result) {
      result = tmp;
    }
  }
  return 1.3 * result;
}

void GeneratorOmegaPiPi::interpolateMajorant() {
  double h = (_maxCMEnergy - _m) / (_nInterpMajorant - 1);
  double energy = _m;
  for (int i = 0; i < _nInterpMajorant; ++i) {
    std::cout << "interpolate : " << i << " / " << _nInterpMajorant << std::endl;
    energy += h;
    _energies[i] = energy;
    _majorants[i] = calcMajorant(energy);
  }
  gsl_spline_init(_spline, _energies, _majorants, _nInterpMajorant);
}

double GeneratorOmegaPiPi::getMajorant(double cmEnergy) const {
  return gsl_spline_eval(_spline, cmEnergy, _acc);
}

void GeneratorOmegaPiPi::genOmegaPiPi(const CFourVector& p) {
  double val = std::numeric_limits<double>::infinity();
  double tmp = 0;
  _star.setInitialMomentum(p);
  while (tmp < val) {
    val = _star.getRndGen().Rndm() * getMajorant(p.getE().real());
    _star.generate();
    tmp = MElCalc::getOmega2PiMEl2
      (std::make_pair(_star.getMomenta()[0], _star.getMomenta()[1]),
       std::make_pair(_star.getMomenta()[2], _star.getMomenta()[3]),
       _star.getMomenta()[4]) * _star.getPhaseSpace();
  }
}
