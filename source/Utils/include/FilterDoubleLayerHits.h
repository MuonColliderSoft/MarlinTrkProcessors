#ifndef FilterDoubleLayerHits_h
#define FilterDoubleLayerHits_h 1

#include "marlin/EventModifier.h"
#include "marlin/Processor.h"

#include "DDRec/SurfaceManager.h"
#include <EVENT/TrackerHitPlane.h>

#include "lcio.h"
#include <string>

#include <TH1.h>

using namespace lcio;
using namespace marlin;

/** Utility processor that removes tracker hits in double layers if they don't have
 *  a corresponding close-by hit in the other sublayer.
 *  Pairs of considered layers are configurable and extracted from the cellID word.
 *  Works for all four lcio hit classes.
 *
 *  @parameter InputCollection name of the hit collection with (Sim)TrackerHits/(Sim)CalorimeterHits
 *  @parameter OutputCollections ( ColName  StartLayer EndLayer )
 *
 * @author N. Bartosik, INFN Torino
 * @date  17 June 2020
 * @version $Id: $
 */

class FilterDoubleLayerHits : public Processor {
protected:
  static const size_t NHITS_MAX = 10000000;

  struct SensorPosition {
    unsigned int layer;
    unsigned int side;
    unsigned int ladder;
    unsigned int module;

    bool operator<(const SensorPosition& rhs) const {
      return std::tie(layer, side, ladder, module) < std::tie(rhs.layer, rhs.side, rhs.ladder, rhs.module);
    }
  };

  /// Double layer cut struct
  struct DoubleLayerCut {
    unsigned int layer0;
    unsigned int layer1;
    double dPhi_max;
    double dTheta_max;
  };

public:
  virtual Processor* newProcessor() { return new FilterDoubleLayerHits; }

  FilterDoubleLayerHits();
  FilterDoubleLayerHits(const FilterDoubleLayerHits&) = delete;
  FilterDoubleLayerHits& operator=(const FilterDoubleLayerHits&) = delete;

  virtual const std::string& name() const { return Processor::name(); }

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
  dd4hep::rec::Vector2D globalToLocal(long int cellID, const dd4hep::rec::Vector3D& posGlobal,
                                      dd4hep::rec::ISurface** surf);

  ////Input collection name.
  std::string _inColName{};

  ////Output collection name.
  std::string _outColName{};

  ////Maximum time difference between hits in a doublet
  double _dtMax{};

  ////Double layer cuts configuration
  StringVec _dlCutConfigs{};

  ////Whether to fill diagnostic histograms
  bool _fillHistos{false};

  ////Subdetector name (needed to get the sensor surface manager)
  std::string _subDetName{};

  std::vector<DoubleLayerCut> _dlCuts{};

  ////Surface map for getting local hit positions at sensor surface
  const dd4hep::rec::SurfaceMap* _map{nullptr};

  ////Array of flags for hits to be accepted
  bool _hitAccepted[NHITS_MAX];

  ////Map of vectors of hits grouped by position in the detector
  std::map<SensorPosition, std::vector<size_t>> _hitsGrouped{};

  ////Monitoring histograms
  std::map<std::string, TH1*> _histos{};

  int _nRun{};
  int _nEvt{};
};

#endif
