#ifndef FilterConeHits_h
#define FilterConeHits_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <string>
#include <vector>

#include <TH1F.h>

using namespace lcio;
using namespace marlin;

/** Utility processor that selects and saves the tracker hits that are included in
 *  a DeltaR cone around the MC particle direction along with the corresponding
 *  sim hits and the reco-sim relations.
 *
 *  @parameter MCParticleCollection name of the MCParticle collection
 *  @parameter TrackerHitInputCollections name of the tracker hit input collections
 *  @parameter TrackerSimHitInputCollections name of the tracker simhit input collections
 *  @parameter TrackerHitInputRelations name of the tracker hit relation input collections
 *  @parameter TrackerHitOutputCollections name of the tracker hit output collections
 *  @parameter TrackerSimHitOutputCollections name of the tracker simhit output collections
 *  @parameter TrackerHitOutputRelations name of the tracker hit relation output collections
 *  @parameter DeltaRCut maximum angular distance between the hits and the particle direction
 *  @parameter FillHistograms flag to fill the diagnostic histograms
 *  @parameter ConeAroundStatus list of Pythia generator statuses that we can choose to cone around
 *
 * @author M. Casarsa, INFN Trieste
 * @date  22 January 2021
 * @version $Id: $
 */

class FilterConeHits : public Processor {
public:
  virtual Processor* newProcessor() { return new FilterConeHits; }

  FilterConeHits();
  FilterConeHits(const FilterConeHits&) = delete;
  FilterConeHits& operator=(const FilterConeHits&) = delete;

  virtual void init();

  virtual void processRunHeader(LCRunHeader* run);

  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  virtual void end();

protected:
  // --- Input/output collection names:
  std::string m_inputMCParticlesCollName{};
  std::vector<std::string> m_inputTrackerHitsCollNames{};
  std::vector<std::string> m_inputTrackerSimHitsCollNames{};
  std::vector<std::string> m_inputTrackerHitRelNames{};
  std::vector<std::string> m_outputTrackerHitsCollNames{};
  std::vector<std::string> m_outputTrackerSimHitsCollNames{};
  std::vector<std::string> m_outputTrackerHitRelNames{};

  // --- Processor parameters:
  bool m_fillHistos{};
  double m_deltaRCut{};
  double m_dist3DCut{};

  // Coning around MCParticles whose generator status is in this list, given in steering file
  std::vector<int> m_coneAroundStatus{};

  // --- Diagnostic histograms:
  TH1F* m_distXY = nullptr;
  TH1F* m_distZ = nullptr;
  TH1F* m_dist3D = nullptr;
  TH1F* m_angle = nullptr;
  TH1F* m_pathLength = nullptr;
  TH1F* m_time = nullptr;

  // --- Magneti field value
  double m_magneticField{};
  const double trackerOuterRadius = 1500.; // mm

  // --- Run and event counters:
  int _nRun{};
  int _nEvt{};
};

#endif
