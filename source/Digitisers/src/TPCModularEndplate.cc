
#include "TPCModularEndplate.h"

#include "marlin/VerbosityLevels.h"
#include <DD4hep/DD4hepUnits.h>

#include <math.h>

TPCModularEndplate::TPCModularEndplate(const dd4hep::rec::FixedPadSizeTPCData* tpc) : _tpc(tpc) {}

void TPCModularEndplate::addModuleRing(unsigned nModules, double phi0) {
  if (isInitialized) {
    throw std::runtime_error("TPCModularEndplate: addModuleRing() called after initialize() ");
  }

  _moduleRings.push_back({nModules, phi0, 0, 0, 0});
}

void TPCModularEndplate::initialize() {
  unsigned nRing = _moduleRings.size();

  double deltaR = (_tpc->rMaxReadout / dd4hep::mm - _tpc->rMinReadout / dd4hep::mm) / nRing;

  double rMin = _tpc->rMinReadout / dd4hep::mm;

  double rMax = rMin + deltaR;

  for (auto& modRing : _moduleRings) {
    modRing.rMin = rMin;
    modRing.rMax = rMax;

    modRing.deltaPhi = 2. * M_PI / modRing.nModules;

    rMin += deltaR;
    rMax += deltaR;
  }

  isInitialized = true;
}

double TPCModularEndplate::computeDistanceRPhi(const dd4hep::rec::Vector3D& hit) {
  if (!isInitialized) {
    throw std::runtime_error("TPCModularEndplate: computeDistanceRPhi() called before initialize() ");
  }

  unsigned indexR = std::floor(_moduleRings.size() * (hit.rho() - _tpc->rMinReadout / dd4hep::mm) /
                               (_tpc->rMaxReadout / dd4hep::mm - _tpc->rMinReadout / dd4hep::mm));

  if (indexR > _moduleRings.size() - 1) {
    streamlog_out(WARNING) << " wrong index  : " << indexR << " for point " << hit << std::endl;

    return 1e6;
  }

  auto modRing = _moduleRings.at(indexR);

  double phi = hit.phi() - modRing.phi0;

  double deltaPhi = std::fabs(std::fmod(phi, modRing.deltaPhi));

  if (deltaPhi > modRing.deltaPhi / 2.) {
    deltaPhi = modRing.deltaPhi - deltaPhi;
  }

  return hit.rho() * deltaPhi;
}
