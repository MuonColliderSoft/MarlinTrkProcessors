#ifndef SetTrackerHitExtensions_h
#define SetTrackerHitExtensions_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <string>

#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitField64.h>

namespace EVENT {
class MCParticle;
class Track;
} // namespace EVENT

namespace IMPL {
class TrackImpl;
}

namespace UTIL {
class LCRelationNavigator;
}

namespace MarlinTrk {
class IMarlinTrkSystem;
}

/**  Set the LCIO Extensions to relate SimTrackerHits to TrackerHits via a pointer.
 *

 *  <h4>Input - Prerequisites</h4>
 *  Needs a collections of LCIO TrackerHits.
 *
 * @param TrackerHitsInputCollections Name of the tracker hit input collections <br>
 * (default value: FTDTrackerHits SITTrackerHits TPCTrackerHits VXDTrackerHits )
 *
 * @param TrackerHitsRelInputCollections Name of the lcrelation collections, that link the TrackerHits to their
 SimTrackerHits.
 * Have to be in same order as TrackerHitsInputCollections!!! <br>
 * (default value: FTDTrackerHitRelations SITTrackerHitRelations TPCTrackerHitRelations VXDTrackerHitRelations )
 *
 *
 *
 * @author S. J. Aplin, DESY
 *
 */
class SetTrackerHitExtensions : public marlin::Processor {
public:
  virtual marlin::Processor* newProcessor() { return new SetTrackerHitExtensions; }

  SetTrackerHitExtensions();
  SetTrackerHitExtensions(const SetTrackerHitExtensions&) = delete;
  SetTrackerHitExtensions& operator=(const SetTrackerHitExtensions&) = delete;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

protected:
  const LCObjectVec* getSimHits(TrackerHit* trkhit, const FloatVec* weights = NULL);

  UTIL::BitField64* _encoder{nullptr};
  int getDetectorID(TrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::subdet()];
  }
  int getSideID(TrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::side()];
  };
  int getLayerID(TrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::layer()];
  };
  int getModuleID(TrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::module()];
  };
  int getSensorID(TrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::sensor()];
  };

  /** helper function to get collection using try catch block */
  LCCollection* GetCollection(LCEvent* evt, std::string colName);

  /** helper function to get relations using try catch block */
  LCRelationNavigator* GetRelations(LCEvent* evt, std::string RelName);

  /** sets up the different collections */
  void SetupInputCollections(LCEvent* evt);

  /** input TrackerHit collections
   */
  std::vector<std::string> _colNamesTrackerHits{};

  /** input relation collections
   */
  std::vector<std::string> _colNamesTrackerHitRelations{};

  //   int _nEventPrintout ;
  int _n_run{};
  int _n_evt{};
  int _current_evt_number{};

  std::vector<LCCollection*> _colTrackerHits{};
  std::vector<LCRelationNavigator*> _navTrackerHitRel{};
};

#endif
