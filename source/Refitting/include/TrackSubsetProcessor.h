#ifndef TrackSubsetProcessor_h
#define TrackSubsetProcessor_h 1

#include "EVENT/Track.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "lcio.h"
#include "marlin/Processor.h"

#include "Math/ProbFunc.h"

using namespace lcio;
using namespace marlin;

/**  Processor that takes tracks from multiple sources and outputs them (or modified versions, or a subset of them)
 * as one track collection.
 *
 *  <h4>Input - Prerequisites</h4>
 *  Track collections
 *
 *  <h4>Output</h4>
 *  A single track collection
 *
 * @param TrackInputCollections A vector of the input track collections <br>
 * (default value ForwardTracks SiTracks )
 *
 * @param TrackOutputCollection Name of the output track collection <br>
 * (default value SubsetTracks )
 *
 * @param MultipleScatteringOn Whether to take multiple scattering into account when fitting the tracks<br>
 * (default value true )
 *
 * @param EnergyLossOn Whether to take energyloss into account when fitting the tracks<br>
 * (default value true )
 *
 * @param SmoothOn Whether to smooth all measurement sites in fit<br>
 * (default value false )
 *
 * @param Omega The parameter omega for the HNN. Controls the influence of the quality indicator. Between 0 and 1:
 * 1 means high influence of quality indicator, 0 means no influence.
 *
 * @author Robin Glattauer, HEPHY
 *
 */

class TrackSubsetProcessor : public Processor {
public:
  virtual Processor* newProcessor() { return new TrackSubsetProcessor; }

  TrackSubsetProcessor();
  TrackSubsetProcessor(const TrackSubsetProcessor&) = delete;
  TrackSubsetProcessor& operator=(const TrackSubsetProcessor&) = delete;

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
  /** helper method that removes short tracks from the list that have the same hits as another, longer track in the list
   *  by setting the corresponding pointer to NULL.
   */
  void removeShortTracks(std::vector<EVENT::Track*>& tracks);

  /** Input collection names */
  std::vector<std::string> _trackInputColNames{};

  /** Output collection name */
  std::string _trackOutputColName{};

  MarlinTrk::IMarlinTrkSystem* _trkSystem{nullptr};
  std::string _trkSystemName{};

  bool _MSOn{};
  bool _ElossOn{};
  bool _SmoothOn{};
  bool _removeShortTracks{};

  float _initialTrackError_d0{};
  float _initialTrackError_phi0{};
  float _initialTrackError_omega{};
  float _initialTrackError_z0{};
  float _initialTrackError_tanL{};

  double _maxChi2PerHit{};

  float _bField{};

  int _nRun{};
  int _nEvt{};

  double _omega{};
};

/** A functor to return whether two tracks are compatible: The criterion is if the share a TrackerHit or more */
class TrackCompatibility {
public:
  inline bool operator()(Track* trackA, Track* trackB) {
    std::vector<TrackerHit*> hitsA = trackA->getTrackerHits();
    std::vector<TrackerHit*> hitsB = trackB->getTrackerHits();

    for (unsigned i = 0; i < hitsA.size(); i++) {
      for (unsigned j = 0; j < hitsB.size(); j++) {
        if (hitsA[i] == hitsB[j])
          return false; // a hit is shared -> incompatible
      }
    }

    return true;
  }
};

/** A functor to return the quality of a track, which is currently the chi2 probability. */
class TrackQI {
public:
  /** @param trkSystem a pointer to an IMarlinTrkSystem, needed for fitting of tracks */
  TrackQI(MarlinTrk::IMarlinTrkSystem* trkSystem) : _trkSystem(trkSystem) {}

  inline double operator()(Track* track) { return ROOT::Math::chisquared_cdf_c(track->getChi2(), track->getNdf()); }

protected:
  MarlinTrk::IMarlinTrkSystem* _trkSystem{nullptr};
};

#endif
