#include "SplitCollectionByLayer.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>

#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
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

  registerProcessorParameter("KeepEmptyCollections",
                             "Whether collections should be added to the event even if they are empty",
                             _addEmptyCollections, false);
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

  unsigned i = 0;
  while (i < _outColAndLayers.size()) {
    std::string name = _outColAndLayers[i++];
    size_t layerStart = std::atoi(_outColAndLayers[i++].c_str());
    size_t layerEnd = std::atoi(_outColAndLayers[i++].c_str());

    // Ensuring that a collection setup is created for this name
    if (_outCols.find(name) == _outCols.end()) {
      _outCols[name] = OutColInfo();
    }

    // Adding layers for the collection setup
    for (size_t layer = layerStart; layer <= layerEnd; ++layer) {
      _outCols[name].layers.push_back(layer);
    }
  }

  streamlog_out(MESSAGE) << "Will split collection " << _colName
                         << " by layers into the following collections:" << std::endl;
  for (const auto& outCol : _outCols) {
    streamlog_out(MESSAGE) << "  " << outCol.first << ":  ";
    for (const auto& layer : outCol.second.layers) {
      streamlog_out(MESSAGE) << layer << " ";
    }
    streamlog_out(MESSAGE) << std::endl;
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
  else if (col->getTypeName() == lcio::LCIO::TRACKERHITPLANE)
    _type = TrackerHitPlaneType;
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
  for (auto& outCol : _outCols) {
    LCCollectionVec* newCol = new LCCollectionVec(col->getTypeName());

    newCol->setSubset(true);

    newCol->parameters().setValue("CellIDEncoding", encoderString);

    outCol.second.collection = newCol;

    streamlog_out(DEBUG5) << " created new output collection " << outCol.first << " of type " << col->getTypeName()
                          << std::endl;
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
    case TrackerHitPlaneType:
      id = cellIDFromHit<TrackerHitPlane>(h);
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

    size_t layerID = encoder[layerIndex];

    // check if we have an output collection for this layer
    for (auto& outCol : _outCols) {
      if (std::find(outCol.second.layers.begin(), outCol.second.layers.end(), layerID) == outCol.second.layers.end())
        continue;
      outCol.second.collection->addElement(h);
      streamlog_out(DEBUG0) << " adding hit for layer " << layerID << " to collection : " << outCol.first << std::endl;
    }
  }

  // add non empty collections to the event
  for (const auto& outCol : _outCols) {
    LCCollection* newCol = outCol.second.collection;

    if (_addEmptyCollections || newCol->getNumberOfElements() > 0) {
      evt->addCollection(newCol, outCol.first);
      streamlog_out(DEBUG5) << " Output collection " << outCol.first << " of type " << col->getTypeName()
                            << " added to the event" << std::endl;
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
