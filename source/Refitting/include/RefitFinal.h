#ifndef RefitFinal_h
#define RefitFinal_h 1

#include <marlin/Processor.h>

#include <UTIL/BitField64.h>

#include <EVENT/Track.h>
#include <IMPL/TrackImpl.h>
#include <cfloat>

namespace MarlinTrk {
class IMarlinTrkSystem;
class IMarlinTrack;
} // namespace MarlinTrk

class RefitFinal : public marlin::Processor {
  /// Nhits cut struct
  struct NHitsCut {
    std::vector<int> detIDs{}; // list of detectors IDs from which hits will be summed
    int nHits_min{};           // minimum number of hits
  };

public:
  virtual marlin::Processor* newProcessor() { return new RefitFinal; }

  RefitFinal();
  RefitFinal(const RefitFinal&) = delete;
  RefitFinal& operator=(const RefitFinal&) = delete;

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

protected:
  int FitInit2(Track* track, MarlinTrk::IMarlinTrack* _marlinTrk);

  /* helper function to get collection using try catch block */
  lcio::LCCollection* GetCollection(lcio::LCEvent* evt, std::string colName);

  /** Input track collection name for refitting.
   */
  std::string _input_track_col_name = "TruthTracks";

  /** output track collection name.
   */
  std::string _output_track_col_name = "RefittedTracks";

  /** Input track relations name.
   */
  std::string _input_track_rel_name = "SiTrackRelations";

  /** Output track relations name for refitting.
   */
  std::string _output_track_rel_name = "RefittedRelation";

  /** pointer to the IMarlinTrkSystem instance
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem = nullptr;

  int _n_run = -1;
  int _n_evt = -1;

  bool _MSOn = true;
  bool _ElossOn = true;

  StringVec _NHitsCutsStr{};
  std::vector<NHitsCut> _NHitsCuts{};
  double _ReducedChi2Cut = -1.0;

  bool _SmoothOn = false;
  double _Max_Chi2_Incr = DBL_MAX;
  int _refPoint = -1;

  float _bField = 0.0;

  bool _extrapolateForward = true;
  int _minClustersOnTrackAfterFit = 0;

  std::shared_ptr<UTIL::BitField64> _encoder{};
};

#endif
