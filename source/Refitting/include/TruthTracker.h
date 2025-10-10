#ifndef TruthTracker_h
#define TruthTracker_h 1

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
class LCEvent;
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

/**  Track creation based on MC truth.
 *

 *  <h4>Input - Prerequisites</h4>
 *  Needs a collections of LCIO TrackerHits.
 *
 *  <h4>Output</h4>
 *  LCIO Track Collection
 *
 * @param TrackerHitsInputCollections Name of the tracker hit input collections <br>
 * (default value: FTDTrackerHits SITTrackerHits TPCTrackerHits VXDTrackerHits )
 *
 * @param TrackerHitsRelInputCollections Name of the lcrelation collections, that link the TrackerHits to their
 SimTrackerHits.
 * Have to be in same order as TrackerHitsInputCollections!!! <br>
 * (default value: FTDTrackerHitRelations SITTrackerHitRelations TPCTrackerHitRelations VXDTrackerHitRelations )
 *
 * @param OutputTrackCollectionName Name of the output Track collection <br>
 * (default value: TruthTracks )
 *
 * @param OutputTrackRelCollection Name of the MCParticle-Track Relations collection for output tracks <br>
 * (default value: TruthTracksMCP )
 *
 * @param MCpThreshold Transverse Momentum Threshold MC particles which will produce tracks GeV <br>
 * (default value: 0.1 )
 *
 * @param FitTracksWithMarlinTrk Fit the Tracks with MarlinTrk, otherwise take track parameters from MCParticle <br>
 * (default value: true )
 *
 * @param MultipleScatteringOn Use MultipleScattering in Fit <br>
 * (default value: true )
 *
 * @param EnergyLossOn Use Energy Loss in Fit <br>
 * (default value: true )
 *
 * @param SmoothOn Smooth All Mesurement Sites in Fit <br>
 * (default value: false )
 *
 *
 * @author S. J. Aplin, DESY ; R. Glattauer, HEPHY
 *
 */
class TruthTracker : public marlin::Processor {
public:
  virtual marlin::Processor* newProcessor() { return new TruthTracker; }

  TruthTracker();
  TruthTracker(const TruthTracker&) = delete;
  TruthTracker& operator=(const TruthTracker&) = delete;

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

  struct SimTrackerHitSortPredicate {
    bool operator()(std::pair<SimTrackerHit*, TrackerHit*> p1, std::pair<SimTrackerHit*, TrackerHit*> p2) {
      SimTrackerHit* simHit1 = p1.first;
      SimTrackerHit* simHit2 = p2.first;

      if (simHit1->getMCParticle() == simHit2->getMCParticle()) {
        return simHit1->getTime() < simHit2->getTime();
      } else {
        return simHit1->getMCParticle() < simHit2->getMCParticle();
      }
    }
  };

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

  void createTrack(MCParticle* mcp, UTIL::BitField64& cellID_encoder,
                   std::vector<std::pair<SimTrackerHit*, TrackerHit*>>& hit_list);

  void createTrack_old(MCParticle* mcp, UTIL::BitField64& cellID_encoder, std::vector<TrackerHit*>& hit_list);

  void createTrack_iterative(MCParticle* mcp, UTIL::BitField64& cellID_encoder,
                             std::vector<std::pair<SimTrackerHit*, TrackerHit*>>& hit_list);

  void drawEvent();

  /** input MCParticle collection
   */
  std::string _colNameMCParticles{};

  /** input TrackerHit collections
   */
  std::vector<std::string> _colNamesTrackerHits{};

  /** input relation collections
   */
  std::vector<std::string> _colNamesTrackerHitRelations{};

  std::vector<LCCollection*> _colTrackerHits{};
  std::vector<LCRelationNavigator*> _navTrackerHitRel{};

  /** output track collection
   */
  std::string _output_track_col_name{};
  LCCollectionVec* _trackVec{nullptr};

  /** Output track relations
   */
  std::string _output_track_rel_name{};
  LCCollectionVec* _trackRelVec{nullptr};

  /** output track segments collection, used for tracks which cannot be formed from a single fit
   */
  std::string _output_track_segments_col_name{};
  LCCollectionVec* _trackSegmentsVec{nullptr};

  /** Output track segments relations, used for tracks which cannot be formed from a single fit
   */
  std::string _output_track_segment_rel_name{};
  LCCollectionVec* _trackSegmentsRelVec{nullptr};

  int _nMCP{};

  int _n_run{};
  int _n_evt{};

  float _MCpThreshold{};

  bool _useMCParticleParametersFotInitOfFit{};

  /** pointer to the IMarlinTrkSystem instance
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem{nullptr};
  bool _runMarlinTrkDiagnostics{};
  std::string _MarlinTrkDiagnosticsName{};

  bool _FitTracksWithMarlinTrk{};
  bool _create_prefit_using_MarlinTrk{};

  bool _MSOn{};
  bool _ElossOn{};
  bool _SmoothOn{};

  float _initialTrackError_d0{};
  float _initialTrackError_phi0{};
  float _initialTrackError_omega{};
  float _initialTrackError_z0{};
  float _initialTrackError_tanL{};

  bool _UseIterativeFitting{};
  bool _UseEventDisplay{};

  double _maxChi2PerHit{};

  double _Bz{};

  unsigned _nCreatedTracks{};

  EVENT::LCEvent* _current_event{nullptr};

  int _detector_model_for_drawing{};
  std::vector<int> _colours{};
  float _helix_max_r{};

  std::string _trkSystemName{};

  int _fitDirection{};
};

#endif
