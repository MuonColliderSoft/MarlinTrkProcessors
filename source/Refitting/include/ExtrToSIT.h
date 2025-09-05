#ifndef ExtrToSIT_h
#define ExtrToSIT_h 1

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

class ExtrToSIT : public marlin::Processor {
public:
  virtual marlin::Processor* newProcessor() { return new ExtrToSIT; }

  ExtrToSIT();
  ExtrToSIT(const ExtrToSIT&) = delete;
  ExtrToSIT& operator=(const ExtrToSIT&) = delete;

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

  void SelectBestCandidate(EVENT::TrackerHitVec& HitsInLayer, const float*& pivot, EVENT::TrackerHit*& BestHit,
                           bool& BestHitFound, int& pointer);

  void SelectBestCandidateLimited(EVENT::TrackerHitVec& HitsInLayer, const float*& pivot, EVENT::TrackerHit*& BestHit,
                                  const FloatVec& covLCIO, double& radius, bool& BestHitFound, double& sigma,
                                  int& pointer, int& PossibleHits, float& dU, float& dV, double& DimDist,
                                  TrackerHitVec& usedSiHits);

  int FitInit(std::vector<TrackerHit*> trackerHits, MarlinTrk::IMarlinTrack* _marlinTrk);
  int FitInit2(Track* track, MarlinTrk::IMarlinTrack* _marlinTrk);

  struct compare_r {
    bool operator()(EVENT::TrackerHit* a, EVENT::TrackerHit* b) const {
      double r_a_sqd = a->getPosition()[0] * a->getPosition()[0] + a->getPosition()[1] * a->getPosition()[1];
      double r_b_sqd = b->getPosition()[0] * b->getPosition()[0] + b->getPosition()[1] * b->getPosition()[1];
      return (r_a_sqd < r_b_sqd);
    }
  };

protected:
  /* helper function to get collection using try catch block */
  lcio::LCCollection* GetCollection(lcio::LCEvent* evt, std::string colName);

  /* helper function to get relations using try catch block */
  lcio::LCRelationNavigator* GetRelations(lcio::LCEvent* evt, std::string RelName);

  /** Input track collection name for refitting.
   */
  std::string _input_track_col_name{};

  /** Input track relations name for refitting.
   */
  std::string _input_track_rel_name{};

  /** Input SIT tracker summer hit collection.
   */
  std::string _sitColName{};

  /** Input VXD tracker summer hit collection.
   */
  std::string _vxdColName{};

  /** refitted track collection name.
   */
  std::string _output_track_col_name{};

  /** Output track relations name for refitting.
   */
  std::string _output_track_rel_name{};

  /** Output silicon track collection.
   */
  std::string _siTrkColName{};

  /** pointer to the IMarlinTrkSystem instance
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem{nullptr};
  std::string _trkSystemName{};

  std::string _mcParticleCollectionName{};

  bool _MSOn{};
  bool _ElossOn{};
  bool _SmoothOn{};
  bool _dirInsideOut{};

  double _Max_Chi2_Incr{};
  int _tpcHitsCut{};
  float _chi2NDoFCut{};
  float _DoCut{};
  float _ZoCut{};
  double _searchSigma{};
  bool _isSpacePoints{};
  int _propToLayer{};

  int _n_run{};
  int _n_evt{};
  int SITHitsFitted{};
  int SITHitsNonFitted{};
  int TotalSITHits{};
  int _nHitsChi2{};

  float _maxChi2PerHit{};
  float _bField{};

  unsigned int _nSITLayers{};

  StringVec _colNamesTrackerHitRelations{};
};

#endif
