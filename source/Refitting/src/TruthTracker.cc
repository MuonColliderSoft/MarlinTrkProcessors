#include "TruthTracker.h"
#include <iostream>

#include <algorithm>
#include <vector>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>

#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <IMPL/TrackImpl.h>

#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/IMarlinTrack.h"

#include "MarlinTrk/MarlinTrkDiagnostics.h"
#ifdef MARLINTRK_DIAGNOSTICS_ON
#include "MarlinTrk/DiagnosticsController.h"
#endif

#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/LCIOTrackPropagators.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitField64.h>
#include <UTIL/BitSet32.h>
#include <UTIL/ILDConf.h>

#include "MarlinCED.h"
#include <UTIL/Operators.h>

using namespace lcio;
using namespace marlin;

using namespace MarlinTrk;

TruthTracker aTruthTracker;

TruthTracker::TruthTracker() : Processor("TruthTracker") {
  // modify processor description
  _description = "Creates Track Collection from MC Truth. Can handle composite spacepoints as long as they consist of "
                 "two TrackerHits";

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection("MCParticle", "MCParticleCollectionName", "Name of the MCParticle input collection",
                          _colNameMCParticles, std::string("MCParticle"));

  StringVec trackerHitsRelInputColNamesDefault;
  trackerHitsRelInputColNamesDefault.push_back("VXDTrackerHitRelations");
  trackerHitsRelInputColNamesDefault.push_back("SITTrackerHitRelations");
  trackerHitsRelInputColNamesDefault.push_back("FTDPixelTrackerHitRelations");
  trackerHitsRelInputColNamesDefault.push_back("FTDSpacePointRelations");
  trackerHitsRelInputColNamesDefault.push_back("TPCTrackerHitRelations");
  trackerHitsRelInputColNamesDefault.push_back("SETTrackerHitRelations");

  registerInputCollections("LCRelation", "TrackerHitsRelInputCollections",
                           "Name of the lcrelation collections, that link the TrackerHits to their SimTrackerHits. "
                           "Have to be in same order as TrackerHitsInputCollections!!!",
                           _colNamesTrackerHitRelations, trackerHitsRelInputColNamesDefault);

  StringVec trackerHitsInputColNamesDefault;

  trackerHitsInputColNamesDefault.push_back("VXDTrackerHits");
  trackerHitsInputColNamesDefault.push_back("SITTrackerHits");
  trackerHitsInputColNamesDefault.push_back("FTDPixelTrackerHits");
  trackerHitsInputColNamesDefault.push_back("FTDSpacePointRelations");
  trackerHitsInputColNamesDefault.push_back("TPCTrackerHits");
  trackerHitsInputColNamesDefault.push_back("SETTrackerHits");

  registerInputCollections("TrackerHit", "TrackerHitsInputCollections", "Name of the tracker hit input collections",
                           _colNamesTrackerHits, trackerHitsInputColNamesDefault);

  registerOutputCollection(LCIO::TRACK, "OutputTrackCollectionName", "Name of the output track collection",
                           _output_track_col_name, std::string("TruthTracks"));

  registerOutputCollection(LCIO::LCRELATION, "OutputTrackRelCollection",
                           "Name of the MCParticle-Track Relations collection for output tracks",
                           _output_track_rel_name, std::string("TruthTracksMCP"));

  registerOutputCollection(LCIO::TRACK, "OutputTrackSegmentCollectionName",
                           "Name of the output track segment collection", _output_track_segments_col_name,
                           std::string("TruthTrackSegments"));

  registerOutputCollection(LCIO::LCRELATION, "OutputTrackSegmentRelCollection",
                           "Name of the MCParticle-Track Relations collection for output track segments",
                           _output_track_segment_rel_name, std::string("TruthTrackSegmentsMCP"));

  //   registerProcessorParameter( "nEventPrintout",
  //                              "Print out progress every N events ",
  //                              _nEventPrintout,
  //                              int(1000));

  registerProcessorParameter("MCpThreshold", "Transverse Momentum Threshold MC particles which will produce tracks GeV",
                             _MCpThreshold, float(0.1));

  registerProcessorParameter("FitTracksWithMarlinTrk",
                             "Fit the Tracks with MarlinTrk, otherwise take track parameters from MCParticle",
                             _FitTracksWithMarlinTrk, bool(true));

  registerProcessorParameter("MultipleScatteringOn", "Use MultipleScattering in Fit", _MSOn, bool(true));

  registerProcessorParameter("EnergyLossOn", "Use Energy Loss in Fit", _ElossOn, bool(true));

  registerProcessorParameter("SmoothOn", "Smooth All Mesurement Sites in Fit", _SmoothOn, bool(false));

  registerProcessorParameter(
      "UseMCParticleParametersFotInitOfFit",
      "When fitting take track parameters from MCParticle for the initialisation of the Track Fit",
      _useMCParticleParametersFotInitOfFit, bool(false));

  registerProcessorParameter(
      "CreatePrefitUsingMarlinTrk",
      "When fitting take track parameters a full pre-fit for the initialisation of the Track Fit",
      _create_prefit_using_MarlinTrk, bool(false));

  registerProcessorParameter("InitialTrackErrorD0", "Value used for the initial d0 variance of the trackfit",
                             _initialTrackError_d0, float(1.e6));

  registerProcessorParameter("InitialTrackErrorPhi0", "Value used for the initial phi0 variance of the trackfit",
                             _initialTrackError_phi0, float(1.e2));

  registerProcessorParameter("InitialTrackErrorOmega", "Value used for the initial omega variance of the trackfit",
                             _initialTrackError_omega, float(1.e-4));

  registerProcessorParameter("InitialTrackErrorZ0", "Value used for the initial z0 variance of the trackfit",
                             _initialTrackError_z0, float(1.e6));

  registerProcessorParameter("InitialTrackErrorTanL", "Value used for the initial tanL variance of the trackfit",
                             _initialTrackError_tanL, float(1.e2));

  registerProcessorParameter("MaxChi2PerHit", "Maximum Chi-squared value allowed when assigning a hit to a track",
                             _maxChi2PerHit, double(1.e2));

  registerProcessorParameter(
      "UseIterativeFitting",
      "Fit the Tracks with MarlinTrk using iterative approach. If 3 consecutive hits fail to be included then the "
      "current fit is written out and a new fit started. Use instead of FitTracksWithMarlinTrk.",
      _UseIterativeFitting, bool(false));

  registerProcessorParameter("UseEventDisplay",
                             "When using UseIterativeFitting show status of each track fit using CED event display.",
                             _UseEventDisplay, bool(false));

  registerProcessorParameter("DetectorTypeForDraw", "Detector type sent to MarlinCED for drawing",
                             _detector_model_for_drawing, int(0));

  registerProcessorParameter("HelixMaxR", "Max R (mm) Extent for drawing Helix if UseTPCForLimitsOfHelix false",
                             _helix_max_r, float(2000.0));

  registerProcessorParameter("TrackSystemName",
                             "Name of the track fitting system to be used (KalTest, DDKalTest, aidaTT, ... )",
                             _trkSystemName, std::string("KalTest"));

  registerProcessorParameter("FitDirection", "Fit direction: -1: backward [default], +1: forward", _fitDirection,
                             int(-1));

#ifdef MARLINTRK_DIAGNOSTICS_ON

  registerOptionalParameter(
      "RunMarlinTrkDiagnostics",
      "Run MarlinTrk Diagnostics. MarlinTrk must be compiled with MARLINTRK_DIAGNOSTICS_ON defined",
      _runMarlinTrkDiagnostics, bool(false));

  registerOptionalParameter("DiagnosticsName", "Name of the root file and root tree if running Diagnostics",
                            _MarlinTrkDiagnosticsName, std::string("TruthTrackerDiagnostics"));

#endif

  _n_run = 0;
  _n_evt = 0;

  _current_event = 0;
}

void TruthTracker::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;

  _encoder = new UTIL::BitField64(lcio::LCTrackerCellID::encoding_string());

  // usually a good idea to
  printParameters();

  if (_UseEventDisplay) {
    MarlinCED::init(this);
  }

  _colours.push_back(0xff00ff);
  _colours.push_back(0xffff00);
  _colours.push_back(0x0000ff);
  _colours.push_back(0xff00ff);
  _colours.push_back(0x00ffff);
  _colours.push_back(0xffffff);

  _colours.push_back(0xff88ff);
  _colours.push_back(0xffff88);
  _colours.push_back(0x8888ff);
  _colours.push_back(0xff88ff);
  _colours.push_back(0x88ffff);
  _colours.push_back(0xffffff);

  // set up the trk system
  _trksystem = MarlinTrk::Factory::createMarlinTrkSystem(_trkSystemName, 0, "");

  if (_trksystem == 0) {
    throw EVENT::Exception(std::string("  Cannot initialize MarlinTrkSystem of Type: ") + _trkSystemName);
  }

  _trksystem->setOption(IMarlinTrkSystem::CFG::useQMS, _MSOn);
  _trksystem->setOption(IMarlinTrkSystem::CFG::usedEdx, _ElossOn);
  _trksystem->setOption(IMarlinTrkSystem::CFG::useSmoothing, _SmoothOn);
  _trksystem->init();

#ifdef MARLINTRK_DIAGNOSTICS_ON

  void* dcv = _trksystem->getDiagnositicsPointer();
  DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
  dc->init(_MarlinTrkDiagnosticsName, _MarlinTrkDiagnosticsName, _runMarlinTrkDiagnostics);

#endif

  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
  double bFieldVec[3];
  lcdd.field().magneticField({0, 0, 0}, bFieldVec); // get the magnetic field vector from DD4hep
  _Bz = bFieldVec[2] / dd4hep::tesla;               // z component at (0,0,0)
}

void TruthTracker::processRunHeader(LCRunHeader*) { ++_n_run; }

void TruthTracker::processEvent(LCEvent* evt) {
  // set the correct configuration for the tracking system for this event
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(_trksystem, _MSOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(_trksystem, _ElossOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);

  _current_event = evt;

  streamlog_out(DEBUG3) << "   processing event: " << _n_evt << std::endl;

  _nMCP = 0;

  _colTrackerHits.clear();
  _navTrackerHitRel.clear();
  _nCreatedTracks = 0;

  /**********************************************************************************************/
  /*                Prepare the collections                                                     */
  /**********************************************************************************************/

  // get the input collections and fill the vectors
  this->SetupInputCollections(evt);

  // establish the track collection that will be created
  _trackVec = new LCCollectionVec(LCIO::TRACK);

  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  _trackVec->setFlag(trkFlag.getFlag());

  // establish the track relations collection that will be created
  _trackRelVec = new LCCollectionVec(LCIO::LCRELATION);

  // track segments

  _trackSegmentsVec = new LCCollectionVec(LCIO::TRACK);

  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkSegFlag(0);
  trkSegFlag.setBit(LCIO::TRBIT_HITS);
  _trackSegmentsVec->setFlag(trkSegFlag.getFlag());

  // establish the track relations collection that will be created
  _trackSegmentsRelVec = new LCCollectionVec(LCIO::LCRELATION);

  // create the encoder to decode cellID0
  UTIL::BitField64 cellID_encoder(LCTrackerCellID::encoding_string());

  /**********************************************************************************************/
  /*                Store the TrackerHits and SimTrackerHits as pairs                           */
  /**********************************************************************************************/

  std::vector<std::pair<SimTrackerHit*, TrackerHit*>> simHitTrkHit;

  for (unsigned iCol = 0; iCol < _colTrackerHits.size(); iCol++) {
    LCCollection* trackerHitCol = _colTrackerHits[iCol];
    int nHits = trackerHitCol->getNumberOfElements();

    LCRelationNavigator* nav = _navTrackerHitRel[iCol];

    for (int j = 0; j < nHits; j++) {
      if (trackerHitCol->getElementAt(j) == 0) {
        streamlog_out(DEBUG0) << "Pointer to TrackerHit" << j << " is NULL " << std::endl;
      }

      TrackerHit* trkhit = dynamic_cast<TrackerHit*>(trackerHitCol->getElementAt(j));

      if (trkhit == 0) {
        std::stringstream errorMsg;
        errorMsg << "dynamic_cast to TrackerHit for hit " << j
                 << " failed. Pointer = " << trackerHitCol->getElementAt(j) << std::endl;

        throw lcio::Exception(errorMsg.str());
      }

      const LCObjectVec& to = nav->getRelatedToObjects(trkhit);

      if (BitSet32(trkhit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) { // it is a composite spacepoint

        if (to.size() == 2) {
          SimTrackerHit* simhitA = dynamic_cast<SimTrackerHit*>(to.at(0));
          SimTrackerHit* simhitB = dynamic_cast<SimTrackerHit*>(to.at(1));

          MCParticle* mcpA = simhitA->getMCParticle();
          MCParticle* mcpB = simhitB->getMCParticle();

          if (mcpA == 0 || mcpB == 0) {
            streamlog_out(DEBUG5) << "spacepoint discarded, because MCParticle from simHit is NULL: mcpA = " << mcpA
                                  << " mcpB = " << mcpB << "\n";
          }

          // Check if the simHits are from the same particel in order to avoid ghost hits
          if (mcpA == mcpB)
            simHitTrkHit.push_back(std::make_pair(simhitA, trkhit));
          else
            streamlog_out(DEBUG0) << "spacepoint discarded, because simHits are not equal " << mcpA << " != " << mcpB
                                  << "\n";

          // #ifdef MARLINTRK_DIAGNOSTICS_ON
          //
          //           // set the pointer to the simhit via lcio extention MCTruth4HitExt
          //
          //           //Split it up and add both hits to the MarlinTrk
          //           const LCObjectVec rawObjects = trkhit->getRawHits();
          //
          //           for( unsigned k=0; k< rawObjects.size(); k++ ){
          //
          //             TrackerHit* rawHit = dynamic_cast< TrackerHit* >( rawObjects[k] );
          //             if( rawHit ){
          //
          //               if( rawHit->getCellID0() == simhitA->getCellID0() ) {
          //                 streamlog_out( DEBUG4 ) << "link simhit = " << simhitA << " Cell ID = " <<
          //                 simhitA->getCellID0() << " with trkhit = " << rawHit << " Cell ID = " <<
          //                 rawHit->getCellID0() << std::endl; rawHit->ext<MarlinTrk::MCTruth4HitExt>() = new
          //                 MarlinTrk::MCTruth4HitExtStruct; rawHit->ext<MarlinTrk::MCTruth4HitExt>()->simhit =
          //                 simhitA;
          //               }
          //               if( rawHit->getCellID0() == simhitB->getCellID0() ) {
          //                 streamlog_out( DEBUG4 ) << "link simhit = " << simhitB << " Cell ID = " <<
          //                 simhitB->getCellID0() << " with trkhit = " << rawHit << " Cell ID = " <<
          //                 rawHit->getCellID0() << std::endl; rawHit->ext<MarlinTrk::MCTruth4HitExt>() = new
          //                 MarlinTrk::MCTruth4HitExtStruct; rawHit->ext<MarlinTrk::MCTruth4HitExt>()->simhit =
          //                 simhitB;
          //               }
          //
          //             }
          //           }
          // #endif

        } else {
          streamlog_out(DEBUG0) << "spacepoint discarded, because it is related to " << to.size()
                                << "SimTrackerHits. It should be 2!\n";
        }

      } else { // no composite spacepoint

        if (to.size() == 1) { // only take trackerHits, that have only one related SimHit

          SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(to.at(0));
          simHitTrkHit.push_back(std::make_pair(simhit, trkhit));

          // #ifdef MARLINTRK_DIAGNOSTICS_ON
          //           trkhit->ext<MarlinTrk::MCTruth4HitExt>() = new MarlinTrk::MCTruth4HitExtStruct;
          //           trkhit->ext<MarlinTrk::MCTruth4HitExt>()->simhit = simhit;
          // #endif
        } else {
          streamlog_out(DEBUG0) << "TrackerHit discarded, because it is related to " << to.size()
                                << "SimTrackerHits. It should be 1!\n";
        }
      }
    }
  }
  streamlog_out(DEBUG4) << "Number of Tracker hits = " << simHitTrkHit.size() << std::endl;

  // now order the hits by MCParticle and then in time
  std::sort(simHitTrkHit.begin(), simHitTrkHit.end(), SimTrackerHitSortPredicate());

  // some information output
  for (unsigned int i = 0; i < simHitTrkHit.size(); ++i) {
    SimTrackerHit* simHit = simHitTrkHit[i].first;
    TrackerHit* trackerHit = simHitTrkHit[i].second;

    streamlog_out(DEBUG2) << "Tracker hit: [" << i << "] = " << trackerHit << "  mcp = " << simHit->getMCParticle()
                          << " time = " << simHit->getTime() << " cellid (trkHit) = " << trackerHit->getCellID0()
                          << " (de" << getDetectorID(trackerHit) << ",si" << getSideID(trackerHit) << ",la"
                          << getLayerID(trackerHit) << ",mo" << getModuleID(trackerHit) << ",se"
                          << getSensorID(trackerHit) << ")" << std::endl;
  }

  /**********************************************************************************************/
  /*                Take hits with the same MCP and create tracks from them                         */
  /**********************************************************************************************/

  streamlog_out(DEBUG4) << "Add Tracker hits to Tracks" << std::endl;

  //  std::vector<TrackerHit*> hit_list;
  std::vector<std::pair<SimTrackerHit*, TrackerHit*>> hit_list;

  if (simHitTrkHit.size() > 0) {
    MCParticle* mcplast = NULL;

    for (unsigned int i = 0; i < simHitTrkHit.size(); ++i) {
      SimTrackerHit* simhit = simHitTrkHit[i].first;

      MCParticle* mcp = simhit->getMCParticle();

      if (mcp == 0) {
        streamlog_out(DEBUG5) << "hit discarded, because MCParticle from simHit is NULL: mcp = " << mcp << "\n";
        continue;
      }

      double const* p = mcp->getMomentum();
      float const pmag2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
      //      float  const pt2   = p[0]*p[0] + p[1]*p[1] ;

      if (i == 0) {
        mcplast = mcp;
      }

      if (mcp != mcplast) {
        // new MCParticle

        if (hit_list.size() >= 3) {
          // create track from vector of hits
          streamlog_out(DEBUG2) << "Create New Track for MCParticle " << mcplast << std::endl;

          if (_UseIterativeFitting) {
            this->createTrack_iterative(mcplast, cellID_encoder, hit_list);
          } else {
            this->createTrack(mcplast, cellID_encoder, hit_list);
          }
        }

        hit_list.clear(); // clear the list for the new mcp
      }

      if (pmag2 > (_MCpThreshold * _MCpThreshold)) { // if momentum is greater than cut add hit to list

        streamlog_out(DEBUG2) << "Add hit from det " << simhit->getCellID0() << " to track from MCParticle " << mcp
                              << " : current number of hits = " << hit_list.size() << std::endl;
        std::pair<SimTrackerHit*, TrackerHit*> pair(simHitTrkHit[i]);
        hit_list.push_back(pair);
      }

      // set last mcparticle
      mcplast = mcp;

    } // end of loop over hits

    // check if there is still a track to be created
    if (hit_list.size() >= 3) {
      // then create a new track
      streamlog_out(DEBUG3) << "Create New Track for Last MCParticle " << mcplast << std::endl;
      if (_UseIterativeFitting) {
        this->createTrack_iterative(mcplast, cellID_encoder, hit_list);
      } else {
        this->createTrack(mcplast, cellID_encoder, hit_list);
      }

      hit_list.clear();
    }
  }

  evt->addCollection(_trackVec, _output_track_col_name);
  evt->addCollection(_trackRelVec, _output_track_rel_name);
  evt->addCollection(_trackSegmentsVec, _output_track_segments_col_name);
  evt->addCollection(_trackSegmentsRelVec, _output_track_segment_rel_name);

  streamlog_out(DEBUG4) << "Created " << _nCreatedTracks << " truth tracks\n";

  // clear the navigator relations
  for (unsigned i = 0; i < _navTrackerHitRel.size(); ++i) {
    delete _navTrackerHitRel[i];
  }

  ++_n_evt;
}

void TruthTracker::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void TruthTracker::end() {
  streamlog_out(DEBUG4) << "TruthTracker::end()  " << name() << " processed " << _n_evt << " events in " << _n_run
                        << " runs " << std::endl;

  delete _encoder;

  //  delete _trksystem ;
}

void TruthTracker::SetupInputCollections(LCEvent* evt) {
  // Check if there are as many tracker hit input collections as relation collections
  if (_colNamesTrackerHits.size() != _colNamesTrackerHitRelations.size()) {
    streamlog_out(ERROR)
        << "There must be as many input collections of tracker Hits as of relations. At the moment, there are "
        << _colNamesTrackerHits.size() << " tracker hit collections and " << _colNamesTrackerHitRelations.size()
        << " relation collections passed as steering paremeters!\n";

    exit(1);
  }

  for (unsigned i = 0; i < _colNamesTrackerHits.size(); i++) {
    // the tracker hits
    LCCollection* colTrkHits = GetCollection(evt, _colNamesTrackerHits[i]);
    if (colTrkHits == NULL)
      continue;

    // the relations of them
    LCRelationNavigator* nav = GetRelations(evt, _colNamesTrackerHitRelations[i]);
    if (nav == NULL)
      continue;

    _colTrackerHits.push_back(colTrkHits);
    _navTrackerHitRel.push_back(nav);
  }
}

void TruthTracker::createTrack(MCParticle* mcp, UTIL::BitField64& cellID_encoder,
                               std::vector<std::pair<SimTrackerHit*, TrackerHit*>>& hit_list) {
  ///////////////////////////////////////////////////////
  // check inputs
  ///////////////////////////////////////////////////////

  if (mcp == 0) {
    throw EVENT::Exception(std::string("TruthTracker::createTrack: MCParticle == NULL "));
  }

  if (hit_list.empty())
    return;

  HelixTrack hel(mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), _Bz);
  streamlog_out(DEBUG3) << "\n MCParticle paramters: "
                        << " d0 " << hel.getD0() << " phi0 " << hel.getPhi0() << " omega " << hel.getOmega() << " z0 "
                        << hel.getZ0() << " tanl " << hel.getTanLambda()
                        << " total number of hits = " << hit_list.size() << "\n"
                        << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////
  // Only fit the hits before looping over SJA:FIXME: FTD may need special treatment
  //////////////////////////////////////////////////////////////////////////////////////

  std::vector<TrackerHit*> hit_list_inner_r;

  hit_list_inner_r.reserve(300);

  float delta_phi = 0.0;

  HelixTrack hel_copy(mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), _Bz);

  for (unsigned ihit = 0; ihit < hit_list.size(); ++ihit) {
    float x = hit_list[ihit].second->getPosition()[0];
    float y = hit_list[ihit].second->getPosition()[1];
    float z = hit_list[ihit].second->getPosition()[2];

    delta_phi += fabs(hel_copy.moveRefPoint(x, y, z));

    if (delta_phi < M_PI) {
      hit_list_inner_r.push_back(hit_list[ihit].second);
    }
  }

  if (hit_list_inner_r.size() < 3) {
    streamlog_out(DEBUG2) << " Reject Track as the number of hits before the turn " << hit_list_inner_r.size()
                          << " is less than 3 hits. Total delta phi before last hit = " << delta_phi
                          << " . Total number of hits = " << hit_list.size() << std::endl;
    return;
  }

  TrackImpl* Track = new TrackImpl;

  if (_FitTracksWithMarlinTrk) {
    // setup initial dummy covariance matrix
    EVENT::FloatVec covMatrix;
    covMatrix.resize(15);

    for (unsigned icov = 0; icov < covMatrix.size(); ++icov) {
      covMatrix[icov] = 0;
    }

    covMatrix[0] = (_initialTrackError_d0);    // sigma_d0^2
    covMatrix[2] = (_initialTrackError_phi0);  // sigma_phi0^2
    covMatrix[5] = (_initialTrackError_omega); // sigma_omega^2
    covMatrix[9] = (_initialTrackError_z0);    // sigma_z0^2
    covMatrix[14] = (_initialTrackError_tanL); // sigma_tanl^2

    streamlog_out(DEBUG3) << "Create track with " << hit_list_inner_r.size() << " hits" << std::endl;

    std::vector<std::pair<float, EVENT::TrackerHit*>> r2_values;
    r2_values.reserve(hit_list_inner_r.size());

    for (TrackerHitVec::iterator it = hit_list_inner_r.begin(); it != hit_list_inner_r.end(); ++it) {
      EVENT::TrackerHit* h = *it;
      float r2 = h->getPosition()[0] * h->getPosition()[0] + h->getPosition()[1] * h->getPosition()[1];
      r2_values.push_back(std::make_pair(r2, *it));
    }

    sort(r2_values.begin(), r2_values.end());

    hit_list_inner_r.clear();
    hit_list_inner_r.reserve(r2_values.size());

    for (std::vector<std::pair<float, EVENT::TrackerHit*>>::iterator it = r2_values.begin(); it != r2_values.end();
         ++it) {
      hit_list_inner_r.push_back(it->second);
    }

    TrackStateImpl* prefit_trackState = 0;

    // bool fit_backwards = IMarlinTrack::backward;

    bool fit_direction = ((_fitDirection < 0) ? IMarlinTrack::backward : IMarlinTrack::forward);

    streamlog_out(DEBUG1) << "TruthTracker::createTrack: fit direction used for fit (-1:backward,+1forward) : "
                          << _fitDirection << std::endl;

    MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();

    try {
      int error = 0;

      if (_useMCParticleParametersFotInitOfFit) {
        hel.moveRefPoint(hit_list_inner_r.front()->getPosition()[0], hit_list_inner_r.front()->getPosition()[1],
                         hit_list_inner_r.front()->getPosition()[2]);

        const float referencePoint[3] = {float(hel.getRefPointX()), float(hel.getRefPointY()),
                                         float(hel.getRefPointZ())};

        prefit_trackState = new TrackStateImpl(lcio::TrackState::AtIP, hel.getD0(), hel.getPhi0(), hel.getOmega(),
                                               hel.getZ0(), hel.getTanLambda(), covMatrix, referencePoint);

        error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, hit_list_inner_r, Track, fit_direction,
                                                    prefit_trackState, _Bz, _maxChi2PerHit);

      } else {
        error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, hit_list_inner_r, Track, fit_direction, covMatrix, _Bz,
                                                    _maxChi2PerHit);
      }

      if (error != IMarlinTrack::success || Track->getNdf() < 0) {
        streamlog_out(DEBUG2) << "TruthTracker::createTrack: EVENT: << " << _current_event->getEventNumber()
                              << " >> Track fit returns error code " << error << " NDF = " << Track->getNdf()
                              << ". Number of hits = " << hit_list_inner_r.size() << std::endl;
        return;
      }

#ifdef MARLINTRK_DIAGNOSTICS_ON
      if (error != IMarlinTrack::success) {
        void* dcv = _trksystem->getDiagnositicsPointer();
        DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
        dc->skip_current_track();
      }
#endif

    } catch (...) {
      streamlog_out(ERROR) << "TruthTracker::createTrack: EVENT: << " << _current_event->getEventNumber()
                           << " >> exception caught and rethown. MCParticle = " << mcp << std::endl;

      //      delete Track;
      //      delete marlinTrk;

      throw;
    }

    std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
    std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;
    std::vector<TrackerHit*> all_hits;
    all_hits.reserve(300);

    marlinTrk->getHitsInFit(hits_in_fit);

    for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      all_hits.push_back(hits_in_fit[ihit].first);
    }

    MarlinTrk::addHitNumbersToTrack(Track, all_hits, true, cellID_encoder);

    marlinTrk->getOutliers(outliers);

    for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      all_hits.push_back(outliers[ihit].first);
    }

    MarlinTrk::addHitNumbersToTrack(Track, all_hits, false, cellID_encoder);

    delete marlinTrk;

  } else {
    // use mcp pos and mom to set the track parameters
    HelixTrack hel_tmp(mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), _Bz);

    hel_tmp.moveRefPoint(0.0, 0.0, 0.0);

    TrackStateImpl* ts = new TrackStateImpl();

    ts->setD0(hel_tmp.getD0());
    ts->setPhi(hel_tmp.getPhi0());
    ts->setOmega(hel_tmp.getOmega());
    ts->setZ0(hel_tmp.getZ0());
    ts->setTanLambda(hel_tmp.getTanLambda());

    float ref[3];

    ref[0] = hel_tmp.getRefPointX();
    ref[1] = hel_tmp.getRefPointY();
    ref[2] = hel_tmp.getRefPointZ();

    ts->setReferencePoint(ref);

    ts->setLocation(lcio::TrackState::AtIP);

    Track->addTrackState(ts);

    std::vector<EVENT::TrackerHit*> added_hits;

    for (unsigned i = 0; i < hit_list.size(); ++i) {
      added_hits.push_back(hit_list[i].second);
    }

    MarlinTrk::addHitNumbersToTrack(Track, added_hits, true, cellID_encoder);
    MarlinTrk::addHitNumbersToTrack(Track, added_hits, false, cellID_encoder);

    for (unsigned i = 0; i < hit_list.size(); i++) {
      Track->addHit(hit_list[i].second);
    }
  }

  streamlog_out(DEBUG3) << "Add Track " << Track << " to collection related to mcp -> " << mcp << std::endl;

  _trackVec->addElement(Track);

  LCRelationImpl* rel = new LCRelationImpl;
  rel->setFrom(Track);
  rel->setTo(mcp);
  rel->setWeight(1.0);
  _trackRelVec->addElement(rel);

  _nCreatedTracks++;
}

void TruthTracker::createTrack_iterative(MCParticle* mcp, UTIL::BitField64& cellID_encoder,
                                         std::vector<std::pair<SimTrackerHit*, TrackerHit*>>& hit_list) {
  ///////////////////////////////////////////////////////
  // check inputs
  ///////////////////////////////////////////////////////

  if (mcp == 0) {
    throw EVENT::Exception(std::string("TruthTracker::createTrack: MCParticle == NULL "));
  }

  if (hit_list.empty()) {
    throw EVENT::Exception(std::string("TruthTracker::createTrack: Hitlist is empty"));
  }

  int layer = 9;
  int size = 3;
  int marker = 1;
  // int ml     = 0 ;
  float helix_max_r = 0;
  float helix_max_z = 0;
  int color = 0;

  if (_UseEventDisplay) {
    MarlinCED::newEvent(this, _detector_model_for_drawing);

    //    CEDPickingHandler &pHandler=CEDPickingHandler::getInstance();
    //
    //    pHandler.update(_current_event);

    helix_max_z = fabs(mcp->getEndpoint()[2]);

    streamlog_out(MESSAGE) << "Draw MCParticle : " << *mcp << std::endl;

    MarlinCED::add_layer_description("MCParticle_For_Fit", layer);

    MarlinCED::drawHelix(_Bz, mcp->getCharge(), mcp->getVertex()[0], mcp->getVertex()[1], mcp->getVertex()[2],
                         mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2], layer, size, 0x7af774,
                         0.0, _helix_max_r, helix_max_z, mcp->id());

    const std::string colName = "Hits_For_Fit";

    size = 10;
    layer = 11;
    //    ml = marker | ( layer << CED_LAYER_SHIFT ) ;

    // ced_describe_layer( colName.c_str() ,layer);
    MarlinCED::add_layer_description(colName, layer);

    color = 0xFFFFFF;

    for (std::vector<std::pair<SimTrackerHit*, TrackerHit*>>::const_iterator it = hit_list.begin();
         it != hit_list.end(); it++) {
      TrackerHit* trkhit = (*it).second;
      ced_hit_ID(trkhit->getPosition()[0], trkhit->getPosition()[1], trkhit->getPosition()[2], marker, layer, size,
                 color, trkhit->id());

    } // hits
  }

  std::vector<EVENT::Track*> track_segments;
  std::vector<IMPL::LCRelationImpl*> track_segments_rels;

  track_segments.reserve(10);

  TrackStateImpl* prefit_trackState = 0;

  // setup initial dummy covariance matrix
  EVENT::FloatVec initial_cov_matrix;
  initial_cov_matrix.resize(15);

  for (unsigned icov = 0; icov < initial_cov_matrix.size(); ++icov) {
    initial_cov_matrix[icov] = 0;
  }

  initial_cov_matrix[0] = (_initialTrackError_d0);    // sigma_d0^2
  initial_cov_matrix[2] = (_initialTrackError_phi0);  // sigma_phi0^2
  initial_cov_matrix[5] = (_initialTrackError_omega); // sigma_omega^2
  initial_cov_matrix[9] = (_initialTrackError_z0);    // sigma_z0^2
  initial_cov_matrix[14] = (_initialTrackError_tanL); // sigma_tanl^2

  // the hits are already ordered in time so there is no need to re-order them

  MarlinTrk::IMarlinTrack* marlinTrk = 0;

  TrackerHitVec added_hits;
  added_hits.reserve(500);
  TrackerHitVec outlier_hits;
  outlier_hits.reserve(20);

  // used to keep track of whether hits are accepted by the fitter before fitting
  bool isSuccessful = false;

  // fit the last hit

  // to initialise we need to add at least one hit.
  // given that we need to also work with COMPOSITE_SPACEPOINT hits we will simply add both strip hits in this case.
  // fit will then be called on the hits added and then addandFit will be used for the remaining hits to provide more
  // feedback and control

  bool fit_running = false;
  unsigned index_of_last_added_hit = 0;
  unsigned running_number_of_rejected_hits = 0;

  // start by trying to add the last hit and break at the first accepted hit
  for (int j = hit_list.size() - 1; j >= 0; --j) {
    TrackerHit* trkHit = hit_list[j].second;

    if (fit_running == false) { // try to start fit

      marlinTrk = _trksystem->createTrack();

      running_number_of_rejected_hits = 0;

      // Check for spacepoints
      if (BitSet32(trkHit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) { // it is a composite spacepoint

        // Split it up and add both hits to the MarlinTrk
        const LCObjectVec rawObjects = trkHit->getRawHits();

        for (unsigned k = 0; k < rawObjects.size(); k++) {
          TrackerHit* rawHit = dynamic_cast<TrackerHit*>(rawObjects[k]);
          if (marlinTrk->addHit(rawHit) == IMarlinTrack::success) {
            isSuccessful = true; // if at least one hit from the spacepoint gets added
          }
        }
      } else { // normal non-composite tracker hit

        isSuccessful = marlinTrk->addHit(trkHit) == IMarlinTrack::success;
      }

      if (isSuccessful) {
        added_hits.push_back(trkHit);
        // at this point one hit added successfully so break the for loop
        fit_running = true;

        // use the simhit to get the track parameters

        SimTrackerHit* simhit_last = hit_list[j].first;

        double mcp_p_last[3];

        mcp_p_last[0] = simhit_last->getMomentum()[0];
        mcp_p_last[1] = simhit_last->getMomentum()[1];
        mcp_p_last[2] = simhit_last->getMomentum()[2];

        HelixTrack hel_at_end(simhit_last->getPosition(), mcp_p_last, mcp->getCharge(), _Bz);

        // set up the initial track and fit parameters

        const float referencePoint[3] = {float(hel_at_end.getRefPointX()), float(hel_at_end.getRefPointY()),
                                         float(hel_at_end.getRefPointZ())};

        prefit_trackState = new TrackStateImpl(lcio::TrackState::AtLastHit, hel_at_end.getD0(), hel_at_end.getPhi0(),
                                               hel_at_end.getOmega(), hel_at_end.getZ0(), hel_at_end.getTanLambda(),
                                               initial_cov_matrix, referencePoint);

        // We now have either one 2-D hit or two 1-D hits added.
        // So proceed with initialisation.

        streamlog_out(DEBUG3) << "\n Helix for prefit: "
                              << " d0 =  " << prefit_trackState->getD0() << " phi0 =  " << prefit_trackState->getPhi()
                              << " omega =  " << prefit_trackState->getOmega()
                              << " z0 =  " << prefit_trackState->getZ0()
                              << " tanl =  " << prefit_trackState->getTanLambda()
                              << " ref =  " << prefit_trackState->getReferencePoint()[0] << " "
                              << prefit_trackState->getReferencePoint()[1] << " "
                              << prefit_trackState->getReferencePoint()[2] << "\n"
                              << std::endl;

        // set the initial track state for the track
        marlinTrk->initialise(*prefit_trackState, _Bz, IMarlinTrack::backward);

        if (_UseEventDisplay) {
          double pt = _Bz * 3e-4 / std::abs(prefit_trackState->getOmega());
          double charge = (prefit_trackState->getOmega() > 0. ? 1. : -1.);

          double px = pt * std::cos(prefit_trackState->getPhi());
          double py = pt * std::sin(prefit_trackState->getPhi());
          double pz = pt * prefit_trackState->getTanLambda();

          // start point for drawing ( PCA to reference point )

          double xs =
              prefit_trackState->getReferencePoint()[0] - prefit_trackState->getD0() * sin(prefit_trackState->getPhi());
          double ys =
              prefit_trackState->getReferencePoint()[1] + prefit_trackState->getD0() * cos(prefit_trackState->getPhi());
          double zs = prefit_trackState->getReferencePoint()[2] + prefit_trackState->getZ0();

          helix_max_z = 2500.0;

          streamlog_out(MESSAGE) << " Draw TrackState for Prefit : " << *prefit_trackState << std::endl;

          layer = 7;
          size = 2;

          MarlinCED::add_layer_description("prefit_trackState", layer);

          MarlinCED::drawHelix(_Bz, charge, xs, ys, zs, px, py, pz, layer, size, 0xFFFFFF, 0.0, _helix_max_r,
                               helix_max_z, prefit_trackState->id());
        }

        delete prefit_trackState;

        // filter the first 1 or 2 hits
        int fit_status = marlinTrk->fit(_maxChi2PerHit);

        streamlog_out(DEBUG4) << "fit_status = " << fit_status << std::endl;

        // check that first hit is accepted by the fit, if this fails we bail here as there is little chance of
        // recovering.

        if (fit_status != IMarlinTrack::success) { // no chance to initialise print warning and return
          streamlog_out(WARNING) << "TruthTracker::createTrack Initial Hit not accepted by the fit, track droped."
                                 << std::endl;
          delete marlinTrk;
          return;
        }

        else {
          index_of_last_added_hit = j;
        }
      }
    }
    // continue to loop over the hits adding them to the fit
    else {
      isSuccessful = false;

      // Check for spacepoints
      if (BitSet32(trkHit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) { // it is a composite spacepoint

        // Split it up and add both hits to the MarlinTrk
        const LCObjectVec rawObjects = trkHit->getRawHits();

        for (int k = rawObjects.size() - 1; k >= 0; --k) {
          TrackerHit* rawHit = dynamic_cast<TrackerHit*>(rawObjects[k]);

          double chi2inc = 0;

          if (marlinTrk->addAndFit(rawHit, chi2inc, _maxChi2PerHit) == IMarlinTrack::success) {
            isSuccessful = true; // if at least one hit from the spacepoint gets added
          }
        }
      }
      // normal non composite tracker hit
      else {
        double chi2inc = 0;

        isSuccessful = marlinTrk->addAndFit(trkHit, chi2inc, _maxChi2PerHit) == IMarlinTrack::success;
      }

      if (isSuccessful) {
        added_hits.push_back(trkHit);
        index_of_last_added_hit = j;
        running_number_of_rejected_hits = 0;
      } else {
        ++running_number_of_rejected_hits;

        // if the running number of rejected hits reaches 3 try to save current track and start a new one
        if (running_number_of_rejected_hits > 3) {
          if (added_hits.size() > 3) {
            if (_SmoothOn)
              marlinTrk->smooth();

            IMPL::TrackImpl* Track = new IMPL::TrackImpl();
            int error = MarlinTrk::finaliseLCIOTrack(marlinTrk, Track, added_hits, IMarlinTrack::backward);

            std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
            std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;
            std::vector<TrackerHit*> all_hits;
            all_hits.reserve(300);

            marlinTrk->getHitsInFit(hits_in_fit);

            for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
              all_hits.push_back(hits_in_fit[ihit].first);
            }

            MarlinTrk::addHitNumbersToTrack(Track, all_hits, true, cellID_encoder);

            marlinTrk->getOutliers(outliers);

            for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
              all_hits.push_back(outliers[ihit].first);
            }

            MarlinTrk::addHitNumbersToTrack(Track, all_hits, false, cellID_encoder);

            delete marlinTrk;
            marlinTrk = 0;

            fit_running = false;

            if (error == IMarlinTrack::success) {
              streamlog_out(DEBUG3) << "Add Track " << Track << " to collection related to mcp -> " << mcp << std::endl;

              track_segments.push_back(Track);
              //              _trackVec->addElement(Track);

              LCRelationImpl* rel = new LCRelationImpl;
              rel->setFrom(Track);
              rel->setTo(mcp);
              rel->setWeight(1.0);

              track_segments_rels.push_back(rel);
              //              _trackRelVec->addElement(rel);

              _nCreatedTracks++;

              if (_UseEventDisplay) {
                const EVENT::TrackState* trkStateIP = Track->getTrackState(EVENT::TrackState::AtIP);

                double pt = _Bz * 3e-4 / std::abs(trkStateIP->getOmega());
                double charge = (trkStateIP->getOmega() > 0. ? 1. : -1.);

                double px = pt * std::cos(trkStateIP->getPhi());
                double py = pt * std::sin(trkStateIP->getPhi());
                double pz = pt * trkStateIP->getTanLambda();

                // start point for drawing ( PCA to reference point )

                double xs = trkStateIP->getReferencePoint()[0] - trkStateIP->getD0() * sin(trkStateIP->getPhi());
                double ys = trkStateIP->getReferencePoint()[1] + trkStateIP->getD0() * cos(trkStateIP->getPhi());
                double zs = trkStateIP->getReferencePoint()[2] + trkStateIP->getZ0();

                helix_max_z = fabs(hit_list.back().second->getPosition()[2]);

                streamlog_out(MESSAGE) << "Draw Partial TrackState : " << *trkStateIP << std::endl;

                layer = 8;
                size = 4;

                //                MarlinCED::add_layer_description("TrackState_From_Fit", layer);

                color = _colours[j % _colours.size()];

                MarlinCED::drawHelix(_Bz, charge, xs, ys, zs, px, py, pz, layer, size, color, 0.0, helix_max_r,
                                     helix_max_z, trkStateIP->id());

                //                ml = marker | ( layer+10 << CED_LAYER_SHIFT ) ;

                for (std::vector<TrackerHit*>::const_iterator it = added_hits.begin(); it != added_hits.end(); it++) {
                  ced_hit_ID((*it)->getPosition()[0], (*it)->getPosition()[1], (*it)->getPosition()[2], marker,
                             layer + 10, 4, color, (*it)->id());

                } // hits
              }
            }
          }
          // reject track
          else {
            streamlog_out(DEBUG3) << "Reject Track : number of hits = " << added_hits.size() << " for MCParticle -> "
                                  << mcp << std::endl;

#ifdef MARLINTRK_DIAGNOSTICS_ON

            void* dcv = _trksystem->getDiagnositicsPointer();
            DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
            dc->skip_current_track();

#endif
            delete marlinTrk;
            marlinTrk = 0;
          }

          // set the counter back to the one after the last hit added
          j = index_of_last_added_hit;
          running_number_of_rejected_hits = 0;
          added_hits.clear();
          fit_running = false;
        }
      }
    }

    if (j == 0) {
      if (added_hits.size() > 3) {
        if (_SmoothOn)
          marlinTrk->smooth();

        IMPL::TrackImpl* Track = new IMPL::TrackImpl();
        int error = MarlinTrk::finaliseLCIOTrack(marlinTrk, Track, added_hits, IMarlinTrack::backward);

        std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
        std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;
        std::vector<TrackerHit*> all_hits;
        all_hits.reserve(300);

        marlinTrk->getHitsInFit(hits_in_fit);

        for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
          all_hits.push_back(hits_in_fit[ihit].first);
        }

        MarlinTrk::addHitNumbersToTrack(Track, all_hits, true, cellID_encoder);

        marlinTrk->getOutliers(outliers);

        for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
          all_hits.push_back(outliers[ihit].first);
        }

        MarlinTrk::addHitNumbersToTrack(Track, all_hits, false, cellID_encoder);

        delete marlinTrk;
        marlinTrk = 0;

        if (error == IMarlinTrack::success) {
          streamlog_out(DEBUG3) << "Add Track " << Track << " to collection related to mcp -> " << mcp << std::endl;

          track_segments.push_back(Track);
          //          _trackVec->addElement(Track);

          LCRelationImpl* rel = new LCRelationImpl;
          rel->setFrom(Track);
          rel->setTo(mcp);
          rel->setWeight(1.0);

          //          _trackRelVec->addElement(rel);
          track_segments_rels.push_back(rel);

          _nCreatedTracks++;

          if (_UseEventDisplay) {
            const EVENT::TrackState* trkStateIP = Track->getTrackState(EVENT::TrackState::AtIP);

            double pt = _Bz * 3e-4 / std::abs(trkStateIP->getOmega());
            double charge = (trkStateIP->getOmega() > 0. ? 1. : -1.);

            double px = pt * std::cos(trkStateIP->getPhi());
            double py = pt * std::sin(trkStateIP->getPhi());
            double pz = pt * trkStateIP->getTanLambda();

            // start point for drawing ( PCA to reference point )

            double xs = trkStateIP->getReferencePoint()[0] - trkStateIP->getD0() * sin(trkStateIP->getPhi());
            double ys = trkStateIP->getReferencePoint()[1] + trkStateIP->getD0() * cos(trkStateIP->getPhi());
            double zs = trkStateIP->getReferencePoint()[2] + trkStateIP->getZ0();

            helix_max_z = fabs(hit_list.back().second->getPosition()[2]);

            streamlog_out(MESSAGE) << " Draw Single TrackState : " << *trkStateIP << std::endl;

            layer = 8;
            size = 2;

            color = _colours[j % _colours.size()];

            MarlinCED::add_layer_description("TrackState_From_Fit", layer);

            MarlinCED::drawHelix(_Bz, charge, xs, ys, zs, px, py, pz, layer, size, color, 0.0, _helix_max_r,
                                 helix_max_z, trkStateIP->id());

            // ml = marker | ( layer+10 << CED_LAYER_SHIFT ) ;

            for (std::vector<TrackerHit*>::const_iterator it = added_hits.begin(); it != added_hits.end(); it++) {
              ced_hit_ID((*it)->getPosition()[0], (*it)->getPosition()[1], (*it)->getPosition()[2], marker, layer + 10,
                         4, color, (*it)->id());

            } // hits
          }
        }
      }
      // reject track
      else {
        streamlog_out(DEBUG3) << "Reject Track : number of hits = " << added_hits.size() << " for MCParticle -> " << mcp
                              << std::endl;

#ifdef MARLINTRK_DIAGNOSTICS_ON

        void* dcv = _trksystem->getDiagnositicsPointer();
        DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
        dc->skip_current_track();

#endif
        delete marlinTrk;
        marlinTrk = 0;
      }
    }

  } // end reverse loop over trakerhits

  // now check how many track segments were created
  if (track_segments.size() != track_segments_rels.size()) {
    throw EVENT::Exception(std::string(
        "TruthTracker::createTrack_iterative: missmatch in size of track_segments and track_segments_rels"));
  }

  if (track_segments.size() == 1) { // only 1 track element so just add it to the collection
    _trackVec->addElement(track_segments[0]);
    _trackRelVec->addElement(track_segments_rels[0]);
  } else if (track_segments.size() != 0) {
    // as the hits are looped over in reverse then the first track segment will be the one which should
    // be used for the fit at the calo face, the last will be the one used for the IP

    // create a track which will be composed from the segements
    IMPL::TrackImpl* Track = new IMPL::TrackImpl();

    for (unsigned itrkseg = 0; itrkseg < track_segments.size(); ++itrkseg) {
      EVENT::Track* seg = track_segments[itrkseg];
      EVENT::LCRelation* rel = track_segments_rels[itrkseg];
      _trackSegmentsVec->addElement(seg);
      _trackSegmentsRelVec->addElement(rel);

      Track->addTrack(seg);

      EVENT::TrackerHitVec hits = seg->getTrackerHits();

      for (unsigned ihit = 0; ihit < hits.size(); ++ihit) {
        Track->addHit(hits[ihit]);
      }

      EVENT::IntVec hitNums = seg->getSubdetectorHitNumbers();
      Track->subdetectorHitNumbers().resize(hitNums.size());

      for (unsigned index = 0; index < hitNums.size(); ++index) {
        Track->subdetectorHitNumbers()[index] += hitNums[index];
      }
    }

    EVENT::Track* innerMostSegment = track_segments.back();
    EVENT::Track* outerMostSegment = track_segments.front();

    Track->setNdf(innerMostSegment->getNdf());
    Track->setChi2(innerMostSegment->getChi2());

    IMPL::TrackStateImpl* atIP = new IMPL::TrackStateImpl(*innerMostSegment->getTrackState(EVENT::TrackState::AtIP));
    Track->addTrackState(atIP);

    IMPL::TrackStateImpl* atFirstHit =
        new IMPL::TrackStateImpl(*innerMostSegment->getTrackState(EVENT::TrackState::AtFirstHit));
    Track->addTrackState(atFirstHit);

    Track->setRadiusOfInnermostHit(innerMostSegment->getRadiusOfInnermostHit());

    IMPL::TrackStateImpl* atLastHit =
        new IMPL::TrackStateImpl(*outerMostSegment->getTrackState(EVENT::TrackState::AtLastHit));

    Track->addTrackState(atLastHit);

    IMPL::TrackStateImpl* atCalo =
        new IMPL::TrackStateImpl(*outerMostSegment->getTrackState(EVENT::TrackState::AtCalorimeter));

    Track->addTrackState(atCalo);

    _trackVec->addElement(Track);

    LCRelationImpl* rel = new LCRelationImpl;
    rel->setFrom(Track);
    rel->setTo(mcp);
    rel->setWeight(1.0);

    _trackRelVec->addElement(rel);
  }

  if (_UseEventDisplay) {
    this->drawEvent();
  }

  delete marlinTrk; // just to be sure
}

void TruthTracker::createTrack_old(MCParticle* mcp, UTIL::BitField64& cellID_encoder,
                                   std::vector<TrackerHit*>& hit_list) {
  TrackImpl* Track = new TrackImpl;

  streamlog_out(DEBUG3) << "Create track with " << hit_list.size() << " hits" << std::endl;

  // check that there are actually hits to fit
  if (hit_list.empty())
    return;

  HelixTrack hel(mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), _Bz);
  streamlog_out(DEBUG3) << "\n MCParticle paramters: "
                        << " d0 " << hel.getD0() << " phi0 " << hel.getPhi0() << " omega " << hel.getOmega() << " z0 "
                        << hel.getZ0() << " tanl " << hel.getTanLambda() << "\n"
                        << std::endl;

  double d0;
  double phi0;
  double omega;
  double z0;
  double tanL;

  double chi2 = 0;
  int ndf = 0;

  float ref[3];

  TrackerHitVec added_hits;

  if (_FitTracksWithMarlinTrk) {
    TrackStateImpl* prefit_trackState = 0;

    // setup initial dummy covariance matrix
    EVENT::FloatVec covMatrix;
    covMatrix.resize(15);

    for (unsigned icov = 0; icov < covMatrix.size(); ++icov) {
      covMatrix[icov] = 0;
    }

    covMatrix[0] = (_initialTrackError_d0);    // sigma_d0^2
    covMatrix[2] = (_initialTrackError_phi0);  // sigma_phi0^2
    covMatrix[5] = (_initialTrackError_omega); // sigma_omega^2
    covMatrix[9] = (_initialTrackError_z0);    // sigma_z0^2
    covMatrix[14] = (_initialTrackError_tanL); // sigma_tanl^2

    // First sort the hits in R, so here we are assuming that the track came from the IP and that we want to fit out to
    // in.
    std::vector<std::pair<float, EVENT::TrackerHit*>> r2_values;

    for (TrackerHitVec::iterator it = hit_list.begin(); it != hit_list.end(); ++it) {
      EVENT::TrackerHit* h = *it;
      float r2 = h->getPosition()[0] * h->getPosition()[0] + h->getPosition()[1] * h->getPosition()[1];
      r2_values.push_back(std::make_pair(r2, *it));
    }

    sort(r2_values.begin(), r2_values.end());

    hit_list.clear();
    hit_list.resize(r2_values.size());

    for (std::vector<std::pair<float, EVENT::TrackerHit*>>::iterator it = r2_values.begin(); it != r2_values.end();
         ++it) {
      hit_list.push_back(it->second);
    }

    // loop over all the hits and create a list consisting only 2D hits

    TrackerHitVec twoD_hits;

    for (unsigned ihit = 0; ihit < hit_list.size(); ++ihit) {
      // check if this a space point or 2D hit
      if (BitSet32(hit_list[ihit]->getType())[UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL] == false) {
        // then add to the list
        twoD_hits.push_back(hit_list[ihit]);
      }
    }

    // now either take the MCParticle parameters directly, use a full prefit or create a helix from 3 2-D hits

    if (_useMCParticleParametersFotInitOfFit) {
      hel.moveRefPoint(hit_list.front()->getPosition()[0], hit_list.front()->getPosition()[1],
                       hit_list.front()->getPosition()[2]);

      const float referencePoint[3] = {float(hel.getRefPointX()), float(hel.getRefPointY()), float(hel.getRefPointZ())};

      prefit_trackState = new TrackStateImpl(lcio::TrackState::AtIP, hel.getD0(), hel.getPhi0(), hel.getOmega(),
                                             hel.getZ0(), hel.getTanLambda(), covMatrix, referencePoint);

    } else { // use 2-D hits to start fit either simple helix construction or full prefit

      // check that there are enough 2-D hits to create a helix
      if (twoD_hits.size() < 3) { // no chance to initialise print warning and return
        streamlog_out(WARNING) << "TruthTracker::createTrack Cannot create helix from less than 3 2-D hits"
                               << std::endl;
        delete Track;
        return;
      }

      // make a helix from 3 hits to get a trackstate
      // SJA:FIXME: this may not be the optimal 3 hits to take in certain cases where the 3 hits are not well spread
      // over the track length
      const double* x1 = twoD_hits[0]->getPosition();
      const double* x2 = twoD_hits[twoD_hits.size() / 2]->getPosition();
      const double* x3 = twoD_hits.back()->getPosition();

      HelixTrack helixTrack(x1, x2, x3, _Bz, HelixTrack::forwards);

      helixTrack.moveRefPoint(hit_list.back()->getPosition()[0], hit_list.back()->getPosition()[1],
                              hit_list.back()->getPosition()[2]);

      const float referencePoint[3] = {float(helixTrack.getRefPointX()), float(helixTrack.getRefPointY()),
                                       float(helixTrack.getRefPointZ())};

      prefit_trackState =
          new TrackStateImpl(lcio::TrackState::AtIP, helixTrack.getD0(), helixTrack.getPhi0(), helixTrack.getOmega(),
                             helixTrack.getZ0(), helixTrack.getTanLambda(), covMatrix, referencePoint);

      // if performing full pre fit use prefit_trackState and improve with full fit

      if (_create_prefit_using_MarlinTrk) {
        MarlinTrk::IMarlinTrack* prefit_trk = _trksystem->createTrack();

#ifdef MARLINTRK_DIAGNOSTICS_ON

        void* dcv = _trksystem->getDiagnositicsPointer();
        DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
        dc->skip_current_track();

#endif

        EVENT::TrackerHitVec::iterator it = hit_list.begin();

        streamlog_out(DEBUG2) << "Start Pre Fit: AddHits: number of hits to fit " << hit_list.size() << std::endl;

        int pre_fit_ndof_added = 0;
        TrackerHitVec pre_fit_added_hits;

        for (it = hit_list.begin(); it != hit_list.end(); ++it) {
          TrackerHit* trkHit = *it;
          bool isSuccessful = false;

          if (BitSet32(
                  trkHit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) { // it is a composite spacepoint

            // Split it up and add both hits to the MarlinTrk
            const LCObjectVec rawObjects = trkHit->getRawHits();

            for (unsigned k = 0; k < rawObjects.size(); k++) {
              TrackerHit* rawHit = dynamic_cast<TrackerHit*>(rawObjects[k]);

              if (prefit_trk->addHit(rawHit) == IMarlinTrack::success) {
                isSuccessful = true; // if at least one hit from the spacepoint gets added
                ++pre_fit_ndof_added;
              }
            }
          } else { // normal non composite hit

            if (prefit_trk->addHit(trkHit) == IMarlinTrack::success) {
              isSuccessful = true;
              pre_fit_ndof_added += 2;
            }
          }

          if (isSuccessful) {
            pre_fit_added_hits.push_back(trkHit);
          }

          else {
            streamlog_out(DEBUG4) << "Hit " << it - hit_list.begin() << " Dropped " << std::endl;
          }
        }

        if (pre_fit_ndof_added < 8) {
          streamlog_out(DEBUG3)
              << "TruthTracker::Prefit failed: Cannot fit less with less than 8 degrees of freedom. Number of hits =  "
              << pre_fit_added_hits.size() << " ndof = " << pre_fit_ndof_added << std::endl;
          delete prefit_trk;
          delete prefit_trackState;
          return;
        }

        prefit_trk->initialise(*prefit_trackState, _Bz, IMarlinTrack::backward);

        delete prefit_trackState;

        int fit_status = prefit_trk->fit();

        if (fit_status != 0) {
          streamlog_out(DEBUG3) << "TruthTracker::Prefit failed: fit_status = " << fit_status << std::endl;
          delete prefit_trk;
          return;
        } else {
          prefit_trackState = new TrackStateImpl();

          const Vector3D point(hit_list.back()->getPosition()[0], hit_list.back()->getPosition()[1],
                               hit_list.back()->getPosition()[2]); // position of outermost hit

          int return_code = prefit_trk->propagate(point, *prefit_trackState, chi2, ndf);

          // make sure that the track state can be propagated to the IP and that the NDF is not less than 0
          if (return_code != 0 || ndf < 0) {
            streamlog_out(WARNING)
                << "TruthTracker::createTrack Prefit: Track droped return_code for prefit propagation = " << return_code
                << " NDF = " << ndf << std::endl;
            delete prefit_trackState;
            delete prefit_trk;
            return;
          }

          prefit_trackState->setLocation(lcio::TrackState::AtLastHit);

          // blow up covariance matrix

          prefit_trackState->setCovMatrix(covMatrix);
        }

        delete prefit_trk;
      }
    }

    MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();

    // used to keep track of whether hits are accepted by the fitter before fitting
    bool isSuccessful = false;

    // needed to keep track of which hit should used for the first addandfit call
    int index_of_next_hit_to_fit = hit_list.size() - 1;

    // to initialise we need to add at least one hit.
    // given that we need to also work with COMPOSITE_SPACEPOINT hits we will simply added both strip hits in this case.
    // fit will then be called on the hits added and then addandFit will be used for the remaining hits to provide more
    // feedback and control

    // start by trying to add the last hit and break at the first accepted hit
    for (int j = index_of_next_hit_to_fit; j >= 0; --j) {
      TrackerHit* trkHit = hit_list[j];

      --index_of_next_hit_to_fit; // first decrement the index of the next hit to fit

      // Check for spacepoints
      if (BitSet32(trkHit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) { // it is a composite spacepoint

        // Split it up and add both hits to the MarlinTrk
        const LCObjectVec rawObjects = trkHit->getRawHits();

        for (unsigned k = 0; k < rawObjects.size(); k++) {
          TrackerHit* rawHit = dynamic_cast<TrackerHit*>(rawObjects[k]);
          if (marlin_trk->addHit(rawHit) == IMarlinTrack::success) {
            isSuccessful = true; // if at least one hit from the spacepoint gets added
          }
        }
      } else { // normal non-composite tracker hit

        isSuccessful = marlin_trk->addHit(trkHit) == IMarlinTrack::success;
      }

      if (isSuccessful) {
        added_hits.push_back(trkHit);
        // at this point one hit added successfully so break the for loop
        break;
      }

    } // end of reverse loop over hits

    // We now have either one 2-D hit or two 1-D hits added.
    // So proceed with initialisation.

    streamlog_out(DEBUG3) << "\n Helix for prefit: "
                          << " d0 =  " << prefit_trackState->getD0() << " phi0 =  " << prefit_trackState->getPhi()
                          << " omega =  " << prefit_trackState->getOmega() << " z0 =  " << prefit_trackState->getZ0()
                          << " tanl =  " << prefit_trackState->getTanLambda()
                          << " ref =  " << prefit_trackState->getReferencePoint()[0] << " "
                          << prefit_trackState->getReferencePoint()[1] << " "
                          << prefit_trackState->getReferencePoint()[2] << "\n"
                          << std::endl;

    // set the initial track state for the track
    marlin_trk->initialise(*prefit_trackState, _Bz, IMarlinTrack::backward);

    delete prefit_trackState;

    // filter the first 1 or 2 hits
    int fit_status = marlin_trk->fit();

    streamlog_out(DEBUG4) << "fit_status = " << fit_status << std::endl;

    // check that first hit is accepted by the fit, if this fails we bail here as there is little chance of recovering.

    if (fit_status != IMarlinTrack::success) { // no chance to initialise print warning and return
      streamlog_out(WARNING) << "TruthTracker::createTrack Initial Hit not accepted by the fit, track droped."
                             << std::endl;
      delete Track;
      delete marlin_trk;
      return;
    }

    // now used addAndFit to add the remaining hits
    for (int j = index_of_next_hit_to_fit; j >= 0; --j) {
      TrackerHit* trkHit = hit_list[j];

      isSuccessful = false;

      // Check for spacepoints
      if (BitSet32(trkHit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) { // it is a composite spacepoint

        // Split it up and add both hits to the MarlinTrk
        const LCObjectVec rawObjects = trkHit->getRawHits();

        for (int k = rawObjects.size() - 1; k >= 0; --k) {
          TrackerHit* rawHit = dynamic_cast<TrackerHit*>(rawObjects[k]);

          double chi_incr = DBL_MAX;

          if (marlin_trk->addAndFit(rawHit, chi_incr) == IMarlinTrack::success) {
            isSuccessful = true; // if at least one hit from the spacepoint gets added
          }
        }
      } else { // normal non composite tracker hit

        double chi_incr = DBL_MAX;

        isSuccessful = marlin_trk->addAndFit(trkHit, chi_incr) == IMarlinTrack::success;
      }
      if (isSuccessful)
        added_hits.push_back(trkHit);
    }

    const Vector3D point(0., 0., 0.); // nominal IP

    TrackStateImpl* trkStateIP = new TrackStateImpl();

    int return_code = marlin_trk->propagate(point, *trkStateIP, chi2, ndf);

    // make sure that the track state can be propagated to the IP and that the NDF is not less than 0
    if (return_code != 0 || ndf < 0) {
      streamlog_out(WARNING) << "TruthTracker::createTrack: Track droped return_code for propagation = " << return_code
                             << " NDF = " << ndf << std::endl;
      delete trkStateIP;
      delete marlin_trk;
      return;
    }

    trkStateIP->setLocation(lcio::TrackState::AtIP);
    Track->trackStates().push_back(trkStateIP);
    Track->setChi2(chi2);
    Track->setNdf(ndf);

    std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;

    marlin_trk->getHitsInFit(hits_in_fit);

    EVENT::TrackerHit* first_trkhit = hits_in_fit.front().first;
    EVENT::TrackerHit* last_trkhit = hits_in_fit.back().first;

    // now loop over the hits and get the track State there

    for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      EVENT::TrackerHit* trkhit = hits_in_fit[ihit].first;

      if (trkhit == 0) {
        throw EVENT::Exception(std::string("TruthTracker::createTrack: TrackerHit pointer == NULL "));
      }

      TrackStateImpl* ts = new TrackStateImpl;
      return_code = marlin_trk->getTrackState(trkhit, *ts, chi2, ndf);

      if (return_code != MarlinTrk::IMarlinTrack::success) {
        streamlog_out(DEBUG5) << "  >>>>>>>>>>> TruthTracker::createTrack:  could not get TrackState at Hit index "
                              << ihit << std::endl;
      }

      // SJA:FIXME: if the hit is 1D then the pivot will be set to 0,0,0 : what should we do? Probably propagate to the
      // Space Point if there is one .... set the location for the first hit
      if (trkhit == first_trkhit) {
        ts->setLocation(lcio::TrackState::AtFirstHit);
      }
      // set the location for the last hit
      else if (trkhit == last_trkhit) {
        ts->setLocation(lcio::TrackState::AtLastHit);
      }
      // set the location for the ihit'th  hit with an offset of lcio::TrackState::AtFirstHit + 1000
      else {
        ts->setLocation(lcio::TrackState::AtFirstHit + 1000 + ihit);
      }

      Track->trackStates().push_back(ts);
    }

    // now get the state at the calo face

    TrackStateImpl* trkStateCalo = new TrackStateImpl;

    UTIL::BitField64 encoder(lcio::LCTrackerCellID::encoding_string());
    encoder.reset(); // reset to 0

    encoder[lcio::LCTrackerCellID::subdet()] = lcio::ILDDetID::ECAL;
    encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::barrel;
    encoder[lcio::LCTrackerCellID::layer()] = 0;

    int detElementID = 0;
    return_code = marlin_trk->propagateToLayer(encoder.lowWord(), last_trkhit, *trkStateCalo, chi2, ndf, detElementID,
                                               IMarlinTrack::modeForward);

    if (return_code == MarlinTrk::IMarlinTrack::no_intersection) { // try forward or backward
      if (trkStateIP->getTanLambda() > 0) {
        encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::fwd;
      } else {
        encoder[lcio::LCTrackerCellID::side()] = lcio::ILDDetID::bwd;
      }
      return_code = marlin_trk->propagateToLayer(encoder.lowWord(), last_trkhit, *trkStateCalo, chi2, ndf, detElementID,
                                                 IMarlinTrack::modeForward);
    }

    if (return_code != MarlinTrk::IMarlinTrack::success) {
      streamlog_out(DEBUG5) << "  >>>>>>>>>>> FinalRefit :  could not get TrackState at Calo Face: return_code = "
                            << return_code << std::endl;
    }

    trkStateCalo->setLocation(lcio::TrackState::AtCalorimeter);
    Track->trackStates().push_back(trkStateCalo);

    delete marlin_trk;

  } else {
    // use mcp pos and mom to set the track parameters
    HelixTrack hel_tmp(mcp->getVertex(), mcp->getMomentum(), mcp->getCharge(), _Bz);
    d0 = hel_tmp.getD0();
    phi0 = hel_tmp.getPhi0();
    omega = hel_tmp.getOmega();
    z0 = hel_tmp.getZ0();
    tanL = hel_tmp.getTanLambda();

    ref[0] = hel_tmp.getRefPointX();
    ref[1] = hel_tmp.getRefPointY();
    ref[2] = hel_tmp.getRefPointZ();

    Track->setD0(d0);
    Track->setPhi(phi0);
    Track->setOmega(omega);
    Track->setZ0(z0);
    Track->setTanLambda(tanL);

    Track->setReferencePoint(ref);

    added_hits = hit_list;
  }

  std::map<int, int> hitNumbers;

  hitNumbers[lcio::ILDDetID::VXD] = 0;
  hitNumbers[lcio::ILDDetID::SIT] = 0;
  hitNumbers[lcio::ILDDetID::FTD] = 0;
  hitNumbers[lcio::ILDDetID::TPC] = 0;
  hitNumbers[lcio::ILDDetID::SET] = 0;
  hitNumbers[lcio::ILDDetID::ETD] = 0;

  for (unsigned int j = 0; j < added_hits.size(); ++j) {
    Track->addHit(added_hits.at(j));

    cellID_encoder.setValue(added_hits.at(j)->getCellID0());
    int detID = cellID_encoder[LCTrackerCellID::subdet()];
    ++hitNumbers[detID];
    //    streamlog_out( DEBUG1 ) << "Hit from Detector " << detID << std::endl;
  }

  // SJA:FIXME no distiction made for hits in fit or not
  Track->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 2] = hitNumbers[lcio::ILDDetID::VXD];
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::FTD - 2] = hitNumbers[lcio::ILDDetID::FTD];
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::SIT - 2] = hitNumbers[lcio::ILDDetID::SIT];
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::TPC - 2] = hitNumbers[lcio::ILDDetID::TPC];
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::SET - 2] = hitNumbers[lcio::ILDDetID::SET];
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::ETD - 2] = hitNumbers[lcio::ILDDetID::ETD];
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 1] = hitNumbers[lcio::ILDDetID::VXD];
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::FTD - 1] = hitNumbers[lcio::ILDDetID::FTD];
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::SIT - 1] = hitNumbers[lcio::ILDDetID::SIT];
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::TPC - 1] = hitNumbers[lcio::ILDDetID::TPC];
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::SET - 1] = hitNumbers[lcio::ILDDetID::SET];
  Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::ETD - 1] = hitNumbers[lcio::ILDDetID::ETD];

  _trackVec->addElement(Track);

  LCRelationImpl* rel = new LCRelationImpl;
  rel->setFrom(Track);
  rel->setTo(mcp);
  rel->setWeight(1.0);
  _trackRelVec->addElement(rel);

  _nCreatedTracks++;
}

LCCollection* TruthTracker::GetCollection(LCEvent* evt, std::string colName) {
  LCCollection* col = NULL;

  try {
    col = evt->getCollection(colName.c_str());
    streamlog_out(DEBUG4) << " --> " << colName.c_str()
                          << " collection found, number of elements = " << col->getNumberOfElements() << std::endl;
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << " --> " << colName.c_str() << " collection absent" << std::endl;
  }

  return col;
}

LCRelationNavigator* TruthTracker::GetRelations(LCEvent* evt, std::string RelName) {
  LCRelationNavigator* nav = NULL;
  LCCollection* col = NULL;

  try {
    col = evt->getCollection(RelName.c_str());
    nav = new LCRelationNavigator(col);
    streamlog_out(DEBUG4) << " --> " << RelName
                          << " track relation collection found, number of elements = " << col->getNumberOfElements()
                          << std::endl;
  } catch (DataNotAvailableException& e) {
    streamlog_out(ERROR) << " --> " << RelName.c_str() << " track relation collection absent" << std::endl;
  }

  return nav;
}

void TruthTracker::drawEvent() {
  LCCollection* colMCP = GetCollection(_current_event, _colNameMCParticles);

  if (colMCP) {
    for (int iMCP = 0; iMCP < colMCP->getNumberOfElements(); ++iMCP) {
      MCParticle* mcp = dynamic_cast<MCParticle*>(colMCP->getElementAt(iMCP));

      float pmag2 = mcp->getMomentum()[0] * mcp->getMomentum()[0] + mcp->getMomentum()[1] * mcp->getMomentum()[1] +
                    mcp->getMomentum()[2] * mcp->getMomentum()[2];

      if (fabs(mcp->getCharge()) > 0.01 && pmag2 > _MCpThreshold * _MCpThreshold) {
        int layer = 0;
        int size = 1;

        float helix_max_r =
            sqrt(mcp->getEndpoint()[0] * mcp->getEndpoint()[0] + mcp->getEndpoint()[1] * mcp->getEndpoint()[1]);

        helix_max_r = _helix_max_r;

        float helix_max_z = fabs(mcp->getEndpoint()[2]);

        MarlinCED::add_layer_description(_colNameMCParticles, layer);

        MarlinCED::drawHelix(_Bz, mcp->getCharge(), mcp->getVertex()[0], mcp->getVertex()[1], mcp->getVertex()[2],
                             mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2], layer, size,
                             _colours[iMCP % _colours.size()], 0.0, helix_max_r, helix_max_z, mcp->id());
      }
    }
  }

  for (unsigned iCol = 0; iCol < _colTrackerHits.size(); iCol++) {
    LCCollection* trackerHitCol = _colTrackerHits[iCol];

    int color = 0xee0044;

    int layer = iCol + 1;
    int marker = 1;
    int size = 10;
    const std::string& colName = _colNamesTrackerHits[iCol];

    // ced_describe_layer( colName.c_str() ,layer);
    MarlinCED::add_layer_description(colName, layer);

    // draw a marker at hit position
    LCTypedVector<TrackerHit> v(trackerHitCol);
    MarlinCED::drawObjectsWithPosition(v.begin(), v.end(), marker, size, color, layer);
  }

  int wait_for_keyboard = 1;
  //++++++++++++++++++++++++++++++++++++
  MarlinCED::draw(this, wait_for_keyboard);
  //++++++++++++++++++++++++++++++++++++
}
