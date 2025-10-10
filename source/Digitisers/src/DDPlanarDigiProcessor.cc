/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "DDPlanarDigiProcessor.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitSet32.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
#include <UTIL/LCRelationNavigator.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"

#include <TMath.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "AIDA/AIDA.h"

#include "marlin/AIDAProcessor.h"
#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"

#include "CLHEP/Vector/TwoVector.h"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace lcio;
using namespace marlin;
using namespace std;

DDPlanarDigiProcessor aDDPlanarDigiProcessor;

DDPlanarDigiProcessor::DDPlanarDigiProcessor() : Processor("DDPlanarDigiProcessor") {
  // modify processor description
  _description = "DDPlanarDigiProcessor creates TrackerHits from SimTrackerHits, smearing their position and time "
                 "according to the input parameters."
                 "The geoemtry of the surface is taken from the DDRec::Surface associated to the hit via the cellID";

  // register steering parameters: name, description, class-variable, default value

  FloatVec resUEx;
  resUEx.push_back(0.0040);

  registerProcessorParameter(
      "ResolutionU", "resolution in direction of u - either one per layer or one for all layers ", _resU, resUEx);

  FloatVec resVEx;
  resVEx.push_back(0.0040);

  registerProcessorParameter(
      "ResolutionV", "resolution in direction of v - either one per layer or one for all layers ", _resV, resVEx);

  FloatVec resTEx;
  resTEx.push_back(-1);

  registerProcessorParameter("ResolutionT",
                             "resolution of time - either one per layer or one for all layers. if the single entry is "
                             "negative, disable time smearing. ",
                             _resT, resTEx);

  registerProcessorParameter("IsStrip", "whether hits are 1D strip hits", _isStrip, bool(false));

  registerProcessorParameter("SubDetectorName", "Name of dub detector", _subDetName, std::string("VXD"));

  // Input collections
  registerInputCollection(LCIO::SIMTRACKERHIT, "SimTrackHitCollectionName",
                          "Name of the Input SimTrackerHit collection", _inColName, std::string("VXDCollection"));

  // Output collections
  registerOutputCollection(LCIO::TRACKERHITPLANE, "TrackerHitCollectionName",
                           "Name of the TrackerHit output collection", _outColName, std::string("VTXTrackerHits"));

  registerOutputCollection(LCIO::LCRELATION, "SimTrkHitRelCollection",
                           "Name of TrackerHit SimTrackHit relation collection", _outRelColName,
                           std::string("VTXTrackerHitRelations"));

  registerProcessorParameter("ForceHitsOntoSurface",
                             "Project hits onto the surface in case they are not yet on the surface (default: false)",
                             _forceHitsOntoSurface, bool(false));

  registerProcessorParameter("MinimumEnergyPerHit", "Minimum Energy (in GeV!) to accept hits, other hits are ignored",
                             _minEnergy, double(0.0));

  registerProcessorParameter("CorrectTimesForPropagation",
                             "Correct hit time for the propagation: radial distance/c (default: false)",
                             _correctTimesForPropagation, bool(false));

  registerProcessorParameter(
      "UseTimeWindow", "Only accept hits with time (after smearing) within the specified time window (default: false)",
      _useTimeWindow, bool(false));

  FloatVec timeWindow_min;
  timeWindow_min.push_back(-1e9);
  registerProcessorParameter(
      "TimeWindowMin",
      "Minimum time a hit must have after smearing to be accepted [ns] - either one per layer or one for all layers",
      _timeWindow_min, timeWindow_min);

  FloatVec timeWindow_max;
  timeWindow_max.push_back(1e9);
  registerProcessorParameter(
      "TimeWindowMax",
      "Maximum time a hit must have after smearing to be accepted [ns] - either one per layer or one for all layers",
      _timeWindow_max, timeWindow_max);

  // setup the list of supported detectors
}

enum { hu = 0, hv, hT, hitE, hitsAccepted, diffu, diffv, diffT, hSize };

void DDPlanarDigiProcessor::init() {
  // usually a good idea to
  printParameters();

  _nRun = 0;
  _nEvt = 0;

  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);

  _h.resize(hSize);

  Global::EVENTSEEDER->registerProcessor(this);

  if (_resU.size() != _resV.size()) {
    std::stringstream ss;
    ss << name() << "::init() - Inconsistent number of resolutions given for U and V coordinate: "
       << "ResolutionU  :" << _resU.size() << " != ResolutionV : " << _resV.size();

    throw EVENT::Exception(ss.str());
  }

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  //===========  get the surface map from the SurfaceManager ================

  dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>();

  dd4hep::DetElement det = theDetector.detector(_subDetName);

  _map = surfMan.map(det.name());

  if (!_map) {
    std::stringstream err;
    err << " Could not find surface map for detector: " << _subDetName << " in SurfaceManager ";
    throw Exception(err.str());
  }

  streamlog_out(DEBUG3) << " DDPlanarDigiProcessor::init(): found " << _map->size()
                        << " surfaces for detector:" << _subDetName << std::endl;

  streamlog_out(MESSAGE) << " *** DDPlanarDigiProcessor::init(): creating histograms" << std::endl;

  AIDAProcessor::histogramFactory(this); //->createHistogram1D( "hMCPEnergy", "energy of the MCParticles", 100 ) ;

  _h[hu] = new TH1F("hu", "smearing u", 50, -5., +5.);
  _h[hv] = new TH1F("hv", "smearing v", 50, -5., +5.);
  _h[hT] = new TH1F("hT", "smearing time", 50, -5., +5.);

  _h[diffu] = new TH1F("diffu", "diff u", 1000, -5., +5.);
  _h[diffv] = new TH1F("diffv", "diff v", 1000, -5., +5.);
  _h[diffT] = new TH1F("diffT", "diff time", 1000, -5., +5.);

  _h[hitE] = new TH1F("hitE", "hitEnergy in keV", 1000, 0, 200);
  _h[hitsAccepted] = new TH1F("hitsAccepted", "Fraction of accepted hits [%]", 201, 0, 100.5);
}

void DDPlanarDigiProcessor::processRunHeader(LCRunHeader*) { ++_nRun; }

void DDPlanarDigiProcessor::processEvent(LCEvent* evt) {
  gsl_rng_set(_rng, Global::EVENTSEEDER->getSeed(this));
  streamlog_out(DEBUG4) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;

  LCCollection* STHcol = 0;
  try {
    STHcol = evt->getCollection(_inColName);
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << "Collection " << _inColName.c_str() << " is unavailable in event " << _nEvt << std::endl;
  }

  if (STHcol != 0) {
    unsigned nCreatedHits = 0;
    unsigned nDismissedHits = 0;

    LCCollectionVec* trkhitVec = new LCCollectionVec(LCIO::TRACKERHITPLANE);

    CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder(lcio::LCTrackerCellID::encoding_string(), trkhitVec);

    // Relation collection TrackerHit, SimTrackerHit
    LCCollection* thsthcol = 0;
    UTIL::LCRelationNavigator thitNav = UTIL::LCRelationNavigator(LCIO::TRACKERHITPLANE, LCIO::SIMTRACKERHIT);

    CellIDDecoder<SimTrackerHit> cellid_decoder(STHcol);

    int nSimHits = STHcol->getNumberOfElements();

    streamlog_out(DEBUG4) << " processing collection " << _inColName << " with " << nSimHits << " hits ... "
                          << std::endl;

    for (int i = 0; i < nSimHits; ++i) {
      SimTrackerHit* simTHit = dynamic_cast<SimTrackerHit*>(STHcol->getElementAt(i));

      _h[hitE]->Fill(simTHit->getEDep() * (dd4hep::GeV / dd4hep::keV));

      if (simTHit->getEDep() < _minEnergy) {
        streamlog_out(DEBUG) << "Hit with insufficient energy " << simTHit->getEDep() * (dd4hep::GeV / dd4hep::keV)
                             << " keV" << std::endl;
        continue;
      }

      const int cellID0 = simTHit->getCellID0();

      //***********************************************************
      // get the measurement surface for this hit using the CellID
      //***********************************************************

      dd4hep::rec::SurfaceMap::const_iterator sI = _map->find(cellID0);

      if (sI == _map->end()) {
        std::cout << " DDPlanarDigiProcessor::processEvent(): no surface found for cellID : "
                  << cellid_decoder(simTHit).valueString() << std::endl;

        std::stringstream err;
        err << " DDPlanarDigiProcessor::processEvent(): no surface found for cellID : "
            << cellid_decoder(simTHit).valueString();
        throw Exception(err.str());
      }

      const dd4hep::rec::ISurface* surf = sI->second;

      int layer = cellid_decoder(simTHit)["layer"];

      dd4hep::rec::Vector3D oldPos(simTHit->getPosition()[0], simTHit->getPosition()[1], simTHit->getPosition()[2]);

      dd4hep::rec::Vector3D newPos;

      //************************************************************
      // Check if Hit is inside sensitive
      //************************************************************

      if (!surf->insideBounds(dd4hep::mm * oldPos)) {
        streamlog_out(DEBUG3) << "  hit at " << oldPos << " " << cellid_decoder(simTHit).valueString()
                              << " is not on surface " << *surf << " distance: " << surf->distance(dd4hep::mm * oldPos)
                              << std::endl;

        if (_forceHitsOntoSurface) {
          dd4hep::rec::Vector2D lv = surf->globalToLocal(dd4hep::mm * oldPos);

          dd4hep::rec::Vector3D oldPosOnSurf = (1. / dd4hep::mm) * surf->localToGlobal(lv);

          streamlog_out(DEBUG3) << " moved to " << oldPosOnSurf << " distance " << (oldPosOnSurf - oldPos).r()
                                << std::endl;

          oldPos = oldPosOnSurf;

        } else {
          ++nDismissedHits;

          continue;
        }
      }

      //***************************************************************
      // Smear time of the hit and apply the time window cut if needed
      //***************************************************************

      double hitT = simTHit->getTime();

      // Smearing time of the hit
      if (_resT.size() and _resT[0] > 0.0) {
        float resT = _resT.size() > 1 ? _resT.at(layer) : _resT.at(0);
        double tSmear = resT > 0.0 ? gsl_ran_gaussian(_rng, resT) : 0.0;
        _h[hT]->Fill(resT > 0.0 ? tSmear / resT : 0.0);
        _h[diffT]->Fill(tSmear);

        hitT += tSmear;
        streamlog_out(DEBUG3) << "smeared hit at T: " << simTHit->getTime() << " ns to T: " << hitT
                              << " ns according to resolution: " << resT << " ns" << std::endl;
      }

      // Correcting for the propagation time
      if (_correctTimesForPropagation) {
        double dt = oldPos.r() / (TMath::C() / 1e6);
        hitT -= dt;
        streamlog_out(DEBUG3) << "corrected hit at R: " << oldPos.r() << " mm by propagation time: " << dt
                              << " ns to T: " << hitT << " ns" << std::endl;
      }

      // Skipping the hit if its time is outside the acceptance time window
      if (_useTimeWindow) {
        float timeWindow_min = _timeWindow_min.size() > 1 ? _timeWindow_min.at(layer) : _timeWindow_min.at(0);
        float timeWindow_max = _timeWindow_max.size() > 1 ? _timeWindow_max.at(layer) : _timeWindow_max.at(0);
        if (hitT < timeWindow_min || hitT > timeWindow_max) {
          streamlog_out(DEBUG4) << "hit at T: " << simTHit->getTime() << " smeared to: " << hitT
                                << " is outside the time window: hit dropped" << std::endl;
          ++nDismissedHits;
          continue;
        }
      }

      //*********************************************************************************
      // Try to smear the hit position but ensure the hit is inside the sensitive region
      //*********************************************************************************

      dd4hep::rec::Vector3D u = surf->u();
      dd4hep::rec::Vector3D v = surf->v();

      // get local coordinates on surface
      dd4hep::rec::Vector2D lv = surf->globalToLocal(dd4hep::mm * oldPos);
      double uL = lv[0] / dd4hep::mm;
      double vL = lv[1] / dd4hep::mm;

      bool accept_hit = false;
      unsigned tries = 0;
      static const unsigned MaxTries = 10;

      float resU = (_resU.size() > 1 ? _resU.at(layer) : _resU.at(0));
      float resV = (_resV.size() > 1 ? _resV.at(layer) : _resV.at(0));

      while (tries < MaxTries) {
        if (tries > 0)
          streamlog_out(DEBUG0) << "retry smearing for " << cellid_decoder(simTHit).valueString() << " : retries "
                                << tries << std::endl;

        double uSmear = gsl_ran_gaussian(_rng, resU);
        double vSmear = gsl_ran_gaussian(_rng, resV);

        // dd4hep::rec::Vector3D newPosTmp = oldPos +  uSmear * u ;
        // if( ! _isStrip )  newPosTmp = newPosTmp +  vSmear * v ;

        dd4hep::rec::Vector3D newPosTmp;
        if (_isStrip) {
          if (_subDetName == "SET") {
            double xStripPos, yStripPos, zStripPos;
            // Find intersection of the strip with the z=centerOfSensor plane to set it as the center of the SET strip
            dd4hep::rec::Vector3D simHitPosSmeared =
                (1. / dd4hep::mm) * (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, 0.)));
            zStripPos = surf->origin()[2] / dd4hep::mm;
            double lineParam = (zStripPos - simHitPosSmeared[2]) / v[2];
            xStripPos = simHitPosSmeared[0] + lineParam * v[0];
            yStripPos = simHitPosSmeared[1] + lineParam * v[1];
            newPosTmp = dd4hep::rec::Vector3D(xStripPos, yStripPos, zStripPos);
          } else {
            newPosTmp =
                (1. / dd4hep::mm) * (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, 0.)));
          }
        } else {
          newPosTmp =
              (1. / dd4hep::mm) *
              (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, (vL + vSmear) * dd4hep::mm)));
        }

        streamlog_out(DEBUG1) << " hit at    : " << oldPos << " smeared to: " << newPosTmp << " uL: " << uL
                              << " vL: " << vL << " uSmear: " << uSmear << " vSmear: " << vSmear << std::endl;

        if (surf->insideBounds(dd4hep::mm * newPosTmp)) {
          accept_hit = true;
          newPos = newPosTmp;

          _h[hu]->Fill(uSmear / resU);
          _h[hv]->Fill(vSmear / resV);

          _h[diffu]->Fill(uSmear);
          _h[diffv]->Fill(vSmear);

          break;

        } else {
          streamlog_out(DEBUG1) << "  hit at " << newPosTmp << " " << cellid_decoder(simTHit).valueString()
                                << " is not on surface "
                                << " distance: " << surf->distance(dd4hep::mm * newPosTmp) << std::endl;
        }

        ++tries;
      }

      if (accept_hit == false) {
        streamlog_out(DEBUG4) << "hit could not be smeared within ladder after " << MaxTries << "  tries: hit dropped"
                              << std::endl;
        ++nDismissedHits;
        continue;
      }

      //**************************************************************************
      // Store hit variables to TrackerHitPlaneImpl
      //**************************************************************************

      TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl;

      const int cellID1 = simTHit->getCellID1();
      trkHit->setCellID0(cellID0);
      trkHit->setCellID1(cellID1);

      trkHit->setPosition(newPos.const_array());
      trkHit->setTime(hitT);
      trkHit->setEDep(simTHit->getEDep());

      float u_direction[2];
      u_direction[0] = u.theta();
      u_direction[1] = u.phi();

      float v_direction[2];
      v_direction[0] = v.theta();
      v_direction[1] = v.phi();

      streamlog_out(DEBUG0) << " U[0] = " << u_direction[0] << " U[1] = " << u_direction[1]
                            << " V[0] = " << v_direction[0] << " V[1] = " << v_direction[1] << std::endl;

      trkHit->setU(u_direction);
      trkHit->setV(v_direction);

      trkHit->setdU(resU);

      if (_isStrip) {
        // store the resolution from the length of the wafer - in case a fitter might want to treat this as 2d hit ....
        double stripRes = (surf->length_along_v() / dd4hep::mm) / std::sqrt(12.);
        trkHit->setdV(stripRes);

      } else {
        trkHit->setdV(resV);
      }

      if (_isStrip) {
        trkHit->setType(UTIL::set_bit(trkHit->getType(), UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL));
      }

      //**************************************************************************
      // Set Relation to SimTrackerHit
      //**************************************************************************

      // Set relation with LCRelationNavigator
      thitNav.addRelation(trkHit, simTHit);

      //**************************************************************************
      // Add hit to collection
      //**************************************************************************

      trkhitVec->addElement(trkHit);

      ++nCreatedHits;

      streamlog_out(DEBUG3) << "-------------------------------------------------------" << std::endl;
    }

    // Create relation collection
    thsthcol = thitNav.createLCCollection();

    // Filling the fraction of accepted hits in the event
    float accFraction = nSimHits > 0 ? float(nCreatedHits) / float(nSimHits) * 100.0 : 0.0;
    _h[hitsAccepted]->Fill(accFraction);

    //**************************************************************************
    // Add collection to event
    //**************************************************************************

    evt->addCollection(trkhitVec, _outColName);
    evt->addCollection(thsthcol, _outRelColName);

    streamlog_out(DEBUG4) << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits  dismissed\n";
  }
  _nEvt++;
}

void DDPlanarDigiProcessor::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void DDPlanarDigiProcessor::end() {
  gsl_rng_free(_rng);

  streamlog_out(MESSAGE) << " end()  " << name() << " processed " << _nEvt << " events in " << _nRun << " runs "
                         << std::endl;
}
