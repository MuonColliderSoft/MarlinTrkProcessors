#include "FilterDoubleLayerHits.h"
#include <cstdlib>
#include <iostream>
#include <string.h>

#include <EVENT/TrackerHitPlane.h>
#include <IMPL/LCCollectionVec.h>

#include "marlin/AIDAProcessor.h"

#include <UTIL/CellIDDecoder.h>

#include "DD4hep/DD4hepUnits.h"

#include <TH1F.h>
#include <TH2F.h>

using namespace lcio;
using namespace marlin;

FilterDoubleLayerHits aFilterDoubleLayerHits;

FilterDoubleLayerHits::FilterDoubleLayerHits() : Processor("FilterDoubleLayerHits"), _dtMax(-1.0), _nRun(0), _nEvt(0) {
  // modify processor description
  _description = "remove hits in double layers if they don't have a corresponding close-by hit in the other sublayer";

  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("InputCollection", "Name of the input collection with hits", _inColName,
                             std::string("VXDTrackerHitPlanes"));

  registerProcessorParameter("OutputCollection", "Name of the output collection with filtered hits", _outColName,
                             std::string("VXDTrackerHitPlanes_DLFiltered"));

  registerProcessorParameter("SubDetectorName", "Name of dub detector", _subDetName, std::string("Vertex"));

  registerProcessorParameter("FillHistograms", "Whether to fill diagnostic histograms", _fillHistos, false);

  registerProcessorParameter("DeltaTimeMax", "Maximum time difference between hits in a doublet [ns]", _dtMax, -1.0);

  StringVec dlCutConfigsEx;
  dlCutConfigsEx.push_back("0");
  dlCutConfigsEx.push_back("1");
  dlCutConfigsEx.push_back("0.5");
  dlCutConfigsEx.push_back("0.05");

  registerProcessorParameter("DoubleLayerCuts",
                             "Layer IDs and angular cuts [mrad] to be applied: <layer 0> <layer 1> <dPhi> <dTheta>",
                             _dlCutConfigs, dlCutConfigsEx);
}

dd4hep::rec::Vector2D FilterDoubleLayerHits::globalToLocal(long int cellID, const dd4hep::rec::Vector3D& posGlobal,
                                                           dd4hep::rec::ISurface** surfptr = nullptr) {
  dd4hep::rec::ISurface* surf;
  // Using directly the provided surface object if available
  if (surfptr && *surfptr) {
    surf = *surfptr;
  } else {
    // Finding the surface corresponding to the cellID
    dd4hep::rec::SurfaceMap::const_iterator surfIt = _map->find(cellID);
    if (surfIt == _map->end()) {
      std::stringstream err;
      err << " FilterDoubleLayerHits::processEvent(): no surface found for cellID: " << cellID;
      throw Exception(err.str());
    }
    surf = surfIt->second;
    // Saving the surface object outside the function to be reused for the same cellID
    if (surfptr)
      *surfptr = surf;
  }
  // Converting global position to local in [cm]
  dd4hep::rec::Vector2D posLocal = surf->globalToLocal(dd4hep::mm * posGlobal);

  return dd4hep::rec::Vector2D(posLocal.u() / dd4hep::mm, posLocal.v() / dd4hep::mm);
}

void FilterDoubleLayerHits::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;
  printParameters();

  // Extracting double-layer cut configurations
  _dlCuts.resize(_dlCutConfigs.size() / 4);

  unsigned i = 0, index = 0;
  while (i < _dlCutConfigs.size()) {
    _dlCuts[index].layer0 = std::atoi(_dlCutConfigs[i++].c_str());
    _dlCuts[index].layer1 = std::atoi(_dlCutConfigs[i++].c_str());
    _dlCuts[index].dPhi_max = std::atof(_dlCutConfigs[i++].c_str()) / 1e3;   // converting mrad -> rad
    _dlCuts[index].dTheta_max = std::atof(_dlCutConfigs[i++].c_str()) / 1e3; // converting mrad -> rad
    ++index;
  }

  // Get the surface map from the SurfaceManager
  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>();
  dd4hep::DetElement det = theDetector.detector(_subDetName);

  _map = surfMan.map(det.name());

  if (!_map) {
    std::stringstream err;
    err << " Could not find surface map for detector: " << _subDetName << " in SurfaceManager";
    throw Exception(err.str());
  }

  // Booking diagnostic histograms for each configured cut
  AIDAProcessor::histogramFactory(this);
  char hname[100];
  for (size_t iCut = 0, nCuts = _dlCuts.size(); iCut < nCuts; ++iCut) {
    const DoubleLayerCut& cut = _dlCuts.at(iCut);
    if (_fillHistos) {
      // Properties of the closest hits across 2 sublayers
      sprintf(hname, "h_dU_layers_%d_%d", cut.layer0, cut.layer1);
      _histos[std::string(hname)] = new TH1F(hname, ";#DeltaU [mm]; Hit pairs", 2000, -5, 5);
      sprintf(hname, "h2_dU_dPhi_layers_%d_%d", cut.layer0, cut.layer1);
      _histos[std::string(hname)] = new TH2F(hname, ";#DeltaU [mm]; #Delta#phi [mrad]", 500, -5, 5, 500, 0, 50);
      sprintf(hname, "h_dTheta_layers_%d_%d", cut.layer0, cut.layer1);
      _histos[std::string(hname)] = new TH1F(hname, ";#Delta#Theta [mrad]; Hit pairs", 3000, -30, 30);
      sprintf(hname, "h_dPhi_layers_%d_%d", cut.layer0, cut.layer1);
      _histos[std::string(hname)] = new TH1F(hname, ";|#Delta#phi| [mrad]; Hit pairs", 1000, 0, 50);
      sprintf(hname, "h_dt_layers_%d_%d", cut.layer0, cut.layer1);
      _histos[std::string(hname)] = new TH1F(hname, ";|#Deltat| [ns]; Hit pairs", 2000, -10, 10);
      // Hit properties for each individual layer
      std::vector<unsigned int> layers{cut.layer0, cut.layer1};
      for (auto layer : layers) {
        sprintf(hname, "h2_posUV_rejected_layer_%d", layer);
        _histos[std::string(hname)] = new TH2F(hname, ";U [mm]; V [mm]", 500, -100, 100, 1000, -200, 200);
      }
    }
    // Printing the configured cut
    streamlog_out(DEBUG5) << iCut << ". layers: " << cut.layer0 << " >> " << cut.layer1 << ";  dPhi: " << cut.dPhi_max
                          << " rad;  dTheta: " << cut.dTheta_max << " rad" << std::endl;
  }

  _nRun = 0;
  _nEvt = 0;
}

void FilterDoubleLayerHits::processRunHeader(LCRunHeader*) { _nRun++; }

void FilterDoubleLayerHits::processEvent(LCEvent* evt) {
  LCCollection* col(0);

  try {
    col = evt->getCollection(_inColName);
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(WARNING) << " input collection not in event : " << _inColName << "   - nothing to do    !!! "
                           << std::endl;
    return;
  }

  std::string encoderString = col->getParameters().getStringVal("CellIDEncoding");
  UTIL::CellIDDecoder<TrackerHitPlane> decoder(encoderString);

  //---- create the output collection
  LCCollectionVec* colOut = new LCCollectionVec(col->getTypeName());
  colOut->setSubset(true);
  colOut->parameters().setValue("CellIDEncoding", encoderString);

  // Set acceptance flags for all hits to FALSE
  const size_t nHit = col->getNumberOfElements();
  memset(&_hitAccepted, false, nHit);

  // Splitting hits by sensor ids for faster association
  _hitsGrouped.clear();
  for (size_t iHit = 0; iHit < nHit; iHit++) {
    TrackerHitPlane* h = (TrackerHitPlane*)col->getElementAt(iHit);

    unsigned int layerID = decoder(h)["layer"];
    unsigned int sideID = decoder(h)["side"];
    unsigned int ladderID = decoder(h)["module"];
    unsigned int moduleID = decoder(h)["sensor"];

    SensorPosition sensPos = {layerID, sideID, ladderID, moduleID};
    if (_hitsGrouped.find(sensPos) == _hitsGrouped.end()) {
      _hitsGrouped[sensPos] = std::vector<size_t>();
      _hitsGrouped[sensPos].reserve(nHit);
    }
    _hitsGrouped[sensPos].push_back(iHit);
  }

  //---- loop over hits
  for (size_t iHit = 0; iHit < nHit; iHit++) {
    // Skipping hits that are already accepted
    if (_hitAccepted[iHit])
      continue;

    TrackerHitPlane* h = (TrackerHitPlane*)col->getElementAt(iHit);

    unsigned int layerID = decoder(h)["layer"];
    unsigned int sideID = decoder(h)["side"];
    unsigned int ladderID = decoder(h)["module"];
    unsigned int moduleID = decoder(h)["sensor"];
    streamlog_out(DEBUG5) << " Checking 1st hit " << iHit << " / " << nHit << " at layer: " << layerID
                          << "  ladder: " << ladderID << "  module: " << moduleID << std::endl;

    const SensorPosition sensPos = {layerID, sideID, ladderID, moduleID};

    // Checking if the hit is at the inner double layer to be filtered
    const DoubleLayerCut* dlCut(0);
    for (int iCut = 0, nCuts = _dlCuts.size(); iCut < nCuts; ++iCut) {
      const DoubleLayerCut& cut = _dlCuts.at(iCut);
      if ((layerID != cut.layer0) && (layerID != cut.layer1))
        continue;
      dlCut = &cut;
      break;
    }

    // Accepting hit immediately if it's not affected by any double-layer cut
    if (!dlCut) {
      _hitAccepted[iHit] = true;
      continue;
    }

    // Skipping the hit from the first pass if it belongs to the outer sublayer of the cut
    if (layerID == dlCut->layer1)
      continue;

    // Getting local and global hit positions
    dd4hep::rec::Vector3D posGlobal(h->getPosition()[0], h->getPosition()[1], h->getPosition()[2]);
    dd4hep::rec::Vector2D posLocal = globalToLocal(decoder(h).getValue(), posGlobal);

    // Setting the values for closest hits
    double dR_min(999.0);
    double dU_closest(0.0);
    double dTheta_closest(0.0);
    double dPhi_closest(0.0);
    double dt_closest(0.0);

    // Looking for the compliment hits in the 2nd sublayer
    size_t nCompatibleHits(0);
    SensorPosition sensPos2 = sensPos;
    sensPos2.layer = dlCut->layer1;
    dd4hep::rec::ISurface* surf = nullptr;
    // Checking if there are any hits in the corresponding sensor at the other sublayer
    if (_hitsGrouped.find(sensPos2) == _hitsGrouped.end())
      continue;
    for (size_t iHit2 : _hitsGrouped.at(sensPos2)) {
      TrackerHitPlane* h2 = (TrackerHitPlane*)col->getElementAt(iHit2);
      unsigned int layerID2 = decoder(h2)["layer"];

      // Checking whether hit is in the time acceptance window
      double dt = h2->getTime() - h->getTime();
      if (_dtMax >= 0.0 && std::fabs(dt) > _dtMax)
        continue;

      // Getting the local and global hit positions
      dd4hep::rec::Vector3D posGlobal2(h2->getPosition()[0], h2->getPosition()[1], h2->getPosition()[2]);
      dd4hep::rec::Vector2D posLocal2 = globalToLocal(decoder(h2).getValue(), posGlobal2, &surf);

      // Checking whether hit is close enough to the 1st one
      double dU = posLocal2.u() - posLocal.u();
      double dTheta = posGlobal2.theta() - posGlobal.theta();
      double dPhi = std::fabs(posGlobal2.phi() - posGlobal.phi());
      if (dPhi > dd4hep::pi)
        dPhi = dd4hep::twopi - dPhi;
      double dR = sqrt(dPhi * dPhi + dTheta * dTheta);
      streamlog_out(DEBUG0) << " Checking 2nd hit at layer: " << layerID2 << ";  dPhi: " << dPhi
                            << ";  dTheta: " << dTheta << std::endl;

      // Updating the minimal values
      if (dR < dR_min) {
        dR_min = dR;
        dU_closest = dU;
        dTheta_closest = dTheta;
        dPhi_closest = dPhi;
        dt_closest = dt;
      }

      // Skipping if the hit is outside the cut window
      if (std::fabs(dPhi) > dlCut->dPhi_max)
        continue;
      if (std::fabs(dTheta) > dlCut->dTheta_max)
        continue;

      nCompatibleHits++;
      _hitAccepted[iHit2] = true;
      streamlog_out(DEBUG5) << " Accepted 2nd hit at layer: " << layerID2 << ";  dPhi: " << dPhi
                            << ";  dTheta: " << dTheta << std::endl;
    }
    // Filling diagnostic histograms
    if (_fillHistos && dR_min < 998) {
      char hname[100];
      sprintf(hname, "h_dU_layers_%d_%d", dlCut->layer0, dlCut->layer1);
      _histos[hname]->Fill(dU_closest);
      sprintf(hname, "h2_dU_dPhi_layers_%d_%d", dlCut->layer0, dlCut->layer1);
      _histos[hname]->Fill(dU_closest, dPhi_closest * 1e3);
      sprintf(hname, "h_dTheta_layers_%d_%d", dlCut->layer0, dlCut->layer1);
      _histos[hname]->Fill(dTheta_closest * 1e3);
      sprintf(hname, "h_dPhi_layers_%d_%d", dlCut->layer0, dlCut->layer1);
      _histos[hname]->Fill(dPhi_closest * 1e3);
      sprintf(hname, "h_dt_layers_%d_%d", dlCut->layer0, dlCut->layer1);
      _histos[hname]->Fill(dt_closest);
    }

    // Accepting the first hit if it has at least one compatible pair
    if (nCompatibleHits > 0) {
      _hitAccepted[iHit] = true;
      streamlog_out(DEBUG5) << " Accepted 1st hit at layer: " << layerID << std::endl;
    }
  }

  // Adding accepted hits to the output collection
  size_t nHitsAccepted(0);
  for (size_t iHit = 0; iHit < nHit; iHit++) {
    if (!_hitAccepted[iHit]) {
      // Filling the positions of rejected hits
      if (_fillHistos) {
        TrackerHitPlane* h = (TrackerHitPlane*)col->getElementAt(iHit);
        unsigned int layerID = decoder(h)["layer"];

        // Getting local hit position
        dd4hep::rec::Vector3D posGlobal(h->getPosition()[0], h->getPosition()[1], h->getPosition()[2]);
        dd4hep::rec::Vector2D posLocal = globalToLocal(decoder(h).getValue(), posGlobal);

        char hname[100];
        sprintf(hname, "h2_posUV_rejected_layer_%d", layerID);
        _histos[hname]->Fill(posLocal.u(), posLocal.v());
      }
      continue;
    }
    colOut->addElement(col->getElementAt(iHit));
    nHitsAccepted++;
  }
  streamlog_out(MESSAGE) << " " << nHitsAccepted << " hits added to collection: " << _outColName << std::endl;

  evt->addCollection(colOut, _outColName);
  streamlog_out(DEBUG5) << " output collection " << _outColName << " of type " << col->getTypeName()
                        << " added to the event  " << std::endl;

  _nEvt++;
}

void FilterDoubleLayerHits::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void FilterDoubleLayerHits::end() {
  streamlog_out(MESSAGE) << "FilterDoubleLayerHits::end()  " << name() << " processed " << _nEvt << " events in "
                         << _nRun << " runs " << std::endl;
}
