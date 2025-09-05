#include "SplitCollectionByLayer.h"
#include <cstdlib>
#include <iostream>

#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>

#include <UTIL/BitField64.h>
// #include "UTIL/LCTrackerConf.h"
// #include <UTIL/LCTOOLS.h>

using namespace lcio;
using namespace marlin;

SplitCollectionByLayer aSplitCollectionByLayer;

SplitCollectionByLayer::SplitCollectionByLayer() : Processor("SplitCollectionByLayer"), _nRun(0), _nEvt(0) {
  // modify processor description
  _description = "split a hit collection based on the layer number of the hits ";

  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("InputCollection", "Name of the input collection with hits", _colName,
                             std::string("FTDCollection"));

  StringVec outColEx;
  outColEx.push_back("FTD_PIXELCollection");
  outColEx.push_back("0");
  outColEx.push_back("1");
  outColEx.push_back("FTD_STRIPCollection");
  outColEx.push_back("2");
  outColEx.push_back("6");

  registerProcessorParameter("OutputCollections", "Name of the output collection with start and end layer number",
                             _outColAndLayers, outColEx);
}

template <class T>
long cellIDFromHit(const LCObject* o) {
  long id = -1;
  const T* h = dynamic_cast<const T*>(o);
  if (h)
    id = (h->getCellID0() & 0xffffffff) | ((long(h->getCellID1()) << 32) & 0xffffffff00000000);
  return id;
}

void SplitCollectionByLayer::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;

  // usually a good idea to
  printParameters();

  _outCols.resize(_outColAndLayers.size() / 3);

  unsigned i = 0, index = 0;
  while (i < _outColAndLayers.size()) {
    _outCols[index].name = _outColAndLayers[i++];
    _outCols[index].layer0 = std::atoi(_outColAndLayers[i++].c_str());
    _outCols[index].layer1 = std::atoi(_outColAndLayers[i++].c_str());
    _outCols[index].collection = 0;

    ++index;
  }

  _nRun = 0;
  _nEvt = 0;
}

void SplitCollectionByLayer::processRunHeader(LCRunHeader*) { _nRun++; }

void SplitCollectionByLayer::processEvent(LCEvent* evt) {
  LCCollection* col = 0;

  try {
    col = evt->getCollection(_colName);

  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(DEBUG5) << " input collection not in event : " << _colName << "   - nothing to do    !!! "
                          << std::endl;
    return;
  }

  // remember the type of the hit collection
  if (col->getTypeName() == lcio::LCIO::SIMTRACKERHIT)
    _type = SimTrackerHitType;
  else if (col->getTypeName() == lcio::LCIO::TRACKERHIT)
    _type = TrackerHitType;
  else if (col->getTypeName() == lcio::LCIO::SIMCALORIMETERHIT)
    _type = SimCalorimeterHitType;
  else if (col->getTypeName() == lcio::LCIO::CALORIMETERHIT)
    _type = CalorimeterHitType;
  else
    _type = UnkownType;

  std::string encoderString = col->getParameters().getStringVal("CellIDEncoding");

  UTIL::BitField64 encoder(encoderString);

  unsigned layerIndex = encoder.index("layer");

  //---- create output collections
  for (unsigned i = 0, N = _outCols.size(); i < N; ++i) {
    LCCollectionVec* newCol = new LCCollectionVec(col->getTypeName());

    newCol->setSubset(true);

    newCol->parameters().setValue("CellIDEncoding", encoderString);

    outCol.second.collection =  newCol ;

    _outCols[i].collection = newCol;

    streamlog_out(DEBUG5) << " create new output collection " << _outCols[i].name << " of type " << col->getTypeName()
                          << std::endl;
    streamlog_out(DEBUG5) << " create new output collection " << _outCols[i].name << " of type " << col->getTypeName()
                          << " to the event  " << std::endl;
  }

  //---- loop over hits

  int nHit = col->getNumberOfElements();

  for (int iHit = 0; iHit < nHit; iHit++) {
    lcio::LCObject* h = col->getElementAt(iHit);

    long id = -1;

    switch (_type) {
    case SimTrackerHitType:
      id = cellIDFromHit<SimTrackerHit>(h);
      break;
    case TrackerHitType:
      id = cellIDFromHit<TrackerHit>(h);
      break;
    case SimCalorimeterHitType:
      id = cellIDFromHit<SimCalorimeterHit>(h);
      break;
    case CalorimeterHitType:
      id = cellIDFromHit<CalorimeterHit>(h);
      break;
    case UnkownType:
      continue;
    }

    encoder.setValue(id);

    unsigned int layerID = encoder[layerIndex];

    // check if we have an output collection for this layer
    for (int i = 0, N = _outCols.size(); i < N; ++i) {
      if ((_outCols[i].layer0 <= layerID) && (layerID <= _outCols[i].layer1)) {
        _outCols[i].collection->addElement(h);

        streamlog_out(DEBUG0) << " adding hit for layer " << layerID << " to collection : " << _outCols[i].name
                              << std::endl;
      }
    }
  }

  // add non empty collections to the event
  for (unsigned i = 0, N = _outCols.size(); i < N; ++i) {
    LCCollection* newCol = _outCols[i].collection;

    if (newCol->getNumberOfElements() > 0) {
      evt->addCollection(newCol, _outCols[i].name);

      streamlog_out(DEBUG5) << " output collection " << _outCols[i].name << " of type " << col->getTypeName()
                            << " added to the event  " << std::endl;
    }
  }

  streamlog_out(DEBUG3) << "   processing event: " << evt->getEventNumber() << "   in run:  " << evt->getRunNumber()
                        << std::endl;

  _nEvt++;
}

void SplitCollectionByLayer::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void SplitCollectionByLayer::end() {
  streamlog_out(MESSAGE) << "SplitCollectionByLayer::end()  " << name() << " processed " << _nEvt << " events in "
                         << _nRun << " runs " << std::endl;
}
