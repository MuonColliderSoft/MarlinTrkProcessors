#ifndef FilterTimeHits_h
#define FilterTimeHits_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <string>
#include <vector>

#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include <TH1F.h>
#include <TMath.h>

using namespace lcio;
using namespace marlin;

/** Utility processor that selects and saves the tracker hits with a polar angle
 *  within a given range (defined by the lower and upper limit), along with
 *  the corresponding sim hits and the reco-sim relations.
 *
 *  @parameter TrackerHitInputCollections name of the tracker hit input collections
 *  @parameter TrackerSimHitInputCollections name of the tracker simhit input collections
 *  @parameter TrackerHitInputRelations name of the tracker hit relation input collections
 *  @parameter TrackerHitOutputCollections name of the tracker hit output collections
 *  @parameter TrackerSimHitOutputCollections name of the tracker simhit output collections
 *  @parameter TrackerHitOutputRelations name of the tracker hit relation output collections
 *  @parameter TargetBeta target beta=v/c for hit time of flight correction
 *  @parameter TimeLowerLimit lower limit on the corrected hit time in ns
 *  @parameter TimeUpperLimit upper limit on the corrected hit time in ns
 *  @parameter FillHistograms flag to fill the diagnostic histograms
 *
 * @author F. Meloni, DESY
 * @date  8 June 2021
 * @version $Id: FilterTimeHits.h,v 0.1 2021-06-08 08:58:00 fmeloni Exp $
 */

class FilterTimeHits : public Processor {
public:
  virtual Processor* newProcessor() { return new FilterTimeHits; }

  FilterTimeHits();
  FilterTimeHits(const FilterTimeHits&) = delete;
  FilterTimeHits& operator=(const FilterTimeHits&) = delete;

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
  double m_beta{};
  double m_time_min{};
  double m_time_max{};

  // --- Diagnostic histograms:
  TH1F* m_corrected_time = nullptr;

  // --- Run and event counters:
  int _nRun{};
  int _nEvt{};
};

#endif
