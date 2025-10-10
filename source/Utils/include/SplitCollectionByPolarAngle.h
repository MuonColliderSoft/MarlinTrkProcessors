#ifndef SplitCollectionByPolarAngle_h
#define SplitCollectionByPolarAngle_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <string>
#include <vector>

#include <TH1F.h>

using namespace lcio;
using namespace marlin;

/** Utility processor that selects and saves the tracker hits with a polar angle
 *  within a given range (defined by the lower and upper limit), along with
 *  the corresponding sim hits and the reco-sim relations.
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
 *
 * @author F. Meloni, DESY; R. Simoniello, CERN; A. Montella, Università degli studi di Trieste/INFN Trieste
 * @date  18 March 2021
 * @version $Id: SplitCollectionByPolarAngle.h,v 0.2 2021-03-18 17:00:00 amontella Exp $
 */

class SplitCollectionByPolarAngle : public Processor {
public:
  virtual Processor* newProcessor() { return new SplitCollectionByPolarAngle; }

  SplitCollectionByPolarAngle();
  SplitCollectionByPolarAngle(const SplitCollectionByPolarAngle&) = delete;
  SplitCollectionByPolarAngle& operator=(const SplitCollectionByPolarAngle&) = delete;

  virtual void init();

  virtual void processRunHeader(LCRunHeader* run);

  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  virtual void end();

protected:
  // --- Input/output collection names:
  std::vector<std::string> m_inputTrackerHitsCollNames{};
  std::vector<std::string> m_inputTrackerSimHitsCollNames{};
  std::vector<std::string> m_inputTrackerHitRelNames{};
  std::vector<std::string> m_outputTrackerHitsCollNames{};
  std::vector<std::string> m_outputTrackerSimHitsCollNames{};
  std::vector<std::string> m_outputTrackerHitRelNames{};

  // --- Processor parameters:
  bool m_fillHistos{};
  double m_theta_min{};
  double m_theta_max{};

  // --- Diagnostic histograms:
  TH1F* m_theta = nullptr;

  // --- Run and event counters:
  int _nRun{};
  int _nEvt{};
};

#endif
