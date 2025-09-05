#ifndef ExtrToTracker_h
#define ExtrToTracker_h 1

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <marlin/Processor.h>
#include <math.h>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitField64.h>
#include <UTIL/BitSet32.h>
#include <UTIL/LCRelationNavigator.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCRelation.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <IMPL/TrackerHitImpl.h>

#include "KiTrack/SubsetHopfieldNN.h"
#include "KiTrack/SubsetSimple.h"

#include "marlin/Global.h"

#include "DDRec/DetectorSurfaces.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceHelper.h"
#include "DDRec/SurfaceManager.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"

using namespace KiTrack;
namespace MarlinTrk {
class IMarlinTrkSystem;
}

class ExtrToTracker : public marlin::Processor {
public:
  virtual marlin::Processor* newProcessor() { return new ExtrToTracker; }

  ExtrToTracker();
  ExtrToTracker(const ExtrToTracker&) = delete;
  ExtrToTracker& operator=(const ExtrToTracker&) = delete;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(lcio::LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(lcio::LCEvent* evt);

  virtual void check(lcio::LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

  int FitInit2(Track* track, MarlinTrk::IMarlinTrack* _marlinTrk);

  struct compare_r {
    bool operator()(EVENT::TrackerHit* a, EVENT::TrackerHit* b) const {
      double r_a_sqd = a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1];
      double r_b_sqd = b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1];
      return (r_a_sqd < r_b_sqd);
    }
  };

  // TrackerHitPlane* getSiHit(LCCollection*& sitHitsCol, int fitElID, MarlinTrk::IMarlinTrack*& marlin_trk, int&
  // nHitsOnDetEl); TrackerHit* getSiHit(LCCollection*& sitHitsCol, int fitElID, MarlinTrk::IMarlinTrack*& marlin_trk,
  // int& nHitsOnDetEl); bool getSiHit(LCCollection*& sitHitsCol, int fitElID, MarlinTrk::IMarlinTrack*& marlin_trk,
  // TrackerHit*& selectedHit);

  TrackerHitPlane* getSiHit(std::vector<TrackerHitPlane*>& hitsOnDetEl, MarlinTrk::IMarlinTrack*& marlin_trk);

  TrackerHitPlane* getSiHit(std::vector<dd4hep::CellID>& vecElID,
                            std::map<int, std::vector<TrackerHitPlane*>>& mapElHits,
                            MarlinTrk::IMarlinTrack*& marlin_trk);

  /* void getNeighbours(int elID, std::vector<int >& vecIDs, std::string cellIDEcoding, std::map<int , int >
   * mapLayerNModules); */

  void fillMapElHits(std::vector<LCCollection*>& vecHitCol,
                     std::vector<std::map<int, std::vector<TrackerHitPlane*>>>& vecMaps);

  /* void addHitOnNextElID(int elementID, MarlinTrk::IMarlinTrack*& marlin_trk, EVENT::TrackerHitVec& trkHits,
   * LCCollection*& sitHitsCol, LCCollection*& otHitsCol, int& iL, int& nSITR, int& TotalSITHits, int& SITHitsPerTrk,
   * int& SITHitsFitted, int& SITHitsNonFitted); */

  void fillVecSubdet(lcio::LCEvent*& evt);

  void getGeoInfo();

  void FindAndAddHit(size_t& idet, int& elID, MarlinTrk::IMarlinTrack*& mtrk, EVENT::TrackerHitVec& trkHits,
                     int& SITHitsPerTrk, int& layer);

  void getCellID0AndPositionInfo(LCCollection*& col);

protected:
  /* helper function to get collection using try catch block */
  lcio::LCCollection* GetCollection(lcio::LCEvent* evt, std::string colName);

  /* /\* helper function to get relations using try catch block *\/ */
  /* lcio::LCRelationNavigator* GetRelations(lcio::LCEvent * evt, std::string RelName ) ; */

  /** Input track collection name for refitting.
   */
  std::string _input_track_col_name{};

  /** output collection name for the not used hits.
   */
  std::string _output_not_used_col_name{};

  /** output track collection name.
   */
  std::string _output_track_col_name{};

  /** Output track relations name for refitting.
   */
  std::string _output_track_rel_name{};

  /** pointer to the IMarlinTrkSystem instance
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem{nullptr};

  /* std::string _mcParticleCollectionName ; */

  bool _MSOn{};
  bool _ElossOn{};
  bool _SmoothOn{};
  double _Max_Chi2_Incr{};
  double _searchSigma{};

  int _n_run{};
  int _n_evt{};
  int SITHitsFitted{};
  int SITHitsNonFitted{};
  int TotalSITHits{};
  int _nHitsChi2{};

  float _bField{};

  bool _performFinalRefit{};

  bool _extrapolateForward{};

  const dd4hep::rec::SurfaceMap* _map{nullptr};

  // processor parameters

  StringVec _vecDigiHits{};
  StringVec _vecSubdetName{};
  std::vector<bool> _vecSubdetIsBarrel{};
  std::vector<int> _vecSubdetNLayers{};
  std::vector<int> _vecSubdetID{};
  std::vector<LCCollection*> _vecDigiHitsCol{};
  std::vector<std::map<dd4hep::CellID, std::vector<dd4hep::CellID>>*> _vecMapNeighbours{};

  std::vector<std::map<int, std::vector<TrackerHitPlane*>>> _vecMapsElHits{};

  std::vector<std::vector<TrackerHitPlane*>> _vecvecHitsInCol{};
};

#endif
