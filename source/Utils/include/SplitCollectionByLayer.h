#ifndef SplitCollectionByLayer_h
#define SplitCollectionByLayer_h 1

#include "marlin/EventModifier.h"
#include "marlin/Processor.h"

#include "lcio.h"
#include <string>

using namespace lcio;
using namespace marlin;

/** Utility processor that allows to split a collection of Hits into
 *  several collections based on the layer information in the cellID word.
 *  Works for all four lcio hit classes.
 *
 *  @parameter InputCollection name of the hit collection with (Sim)TrackerHits/(Sim)CalorimeterHits
 *  @parameter OutputCollections ( ColName  StartLayer EndLayer )
 *
 * @author F. Gaede, CERN/DESY
 * @date  30 Oct 2014
 * @version $Id: $
 */

class SplitCollectionByLayer : public Processor {
protected:
  /// helper struct
  struct OutColInfo {
    std::string name{};
    unsigned layer0{};
    unsigned layer1{};
    LCCollection* collection{nullptr};
  };

  /// Enum used for hit types
  enum HitType { SimTrackerHitType = 1, TrackerHitType, SimCalorimeterHitType, CalorimeterHitType, UnkownType };

public:
  virtual Processor* newProcessor() { return new SplitCollectionByLayer; }

  SplitCollectionByLayer();
  SplitCollectionByLayer(const SplitCollectionByLayer&) = delete;
  SplitCollectionByLayer& operator=(const SplitCollectionByLayer&) = delete;

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
  ////Input collection name.
  std::string _colName{};

  /// Output collections and layers:
  StringVec _outColAndLayers{};

  std::vector<OutColInfo> _outCols{};

  HitType _type{};

  int _nRun{};
  int _nEvt{};
};

#endif
