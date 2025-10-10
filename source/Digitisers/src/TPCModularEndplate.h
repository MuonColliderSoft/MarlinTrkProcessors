#ifndef TPCModularEndplate_h
#define TPCModularEndplate_h

#include <DDRec/DetectorData.h>
#include <DDRec/Vector3D.h>

#include <vector>

/** Helper class for defining a modular TPC endplate with regular modules.
 *  Computes the distance from a module boundary for hits.
 *
 * @author F.Gaede, DESY
 * @date June, 2017
 */

class TPCModularEndplate {
public:
  /// internal helper struct
  struct ModuleRing {
    unsigned nModules;
    double phi0;
    double rMin;
    double rMax;
    double deltaPhi;
  };

  /// no default c'tor
  TPCModularEndplate() = delete;
  TPCModularEndplate(const TPCModularEndplate&) = delete;
  TPCModularEndplate& operator=(const TPCModularEndplate&) = delete;

  /// intitialize for the gice TPC data
  TPCModularEndplate(const dd4hep::rec::FixedPadSizeTPCData* tpc);

  /// add a module ring with the given number of modules and phi0 - modules are added inside out.
  void addModuleRing(unsigned nModules, double phi0);

  /// inititalize with the modules added so far - after this no more modules can be added.
  void initialize();

  /// compute the distance in the rphi-plane from a module boundary in mm.
  double computeDistanceRPhi(const dd4hep::rec::Vector3D& hit);

protected:
  std::vector<ModuleRing> _moduleRings{};

  const dd4hep::rec::FixedPadSizeTPCData* _tpc{};

  bool isInitialized{false};
};

#endif
