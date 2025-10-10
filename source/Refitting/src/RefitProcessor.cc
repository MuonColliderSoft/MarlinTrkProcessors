#include "RefitProcessor.h"

#include <algorithm>

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>

#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/Operators.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"

using namespace lcio;
using namespace marlin;

using namespace MarlinTrk;

RefitProcessor aRefitProcessor;

RefitProcessor::RefitProcessor() : Processor("RefitProcessor") {
  // modify processor description
  _description = "RefitProcessor refits an input track collection, producing a new collection of tracks.";

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection(LCIO::TRACK, "InputTrackCollectionName", "Name of the input track collection",
                          _input_track_col_name, std::string("TruthTracks"));

  registerInputCollection(LCIO::LCRELATION, "InputTrackRelCollection",
                          "Name of the MCParticle-Track Relations collection for input tracks", _input_track_rel_name,
                          std::string("TruthTracksMCP"));

  registerOutputCollection(LCIO::TRACK, "OutputTrackCollectionName", "Name of the output track collection",
                           _output_track_col_name, std::string("RefittedTracks"));

  registerOutputCollection(LCIO::LCRELATION, "OutputTrackRelCollection",
                           "Name of the MCParticle-Track Relations collection for output tracks",
                           _output_track_rel_name, std::string("RefittedTracksMCP"));

  registerProcessorParameter("MultipleScatteringOn", "Use MultipleScattering in Fit", _MSOn, bool(true));

  registerProcessorParameter("EnergyLossOn", "Use Energy Loss in Fit", _ElossOn, bool(true));

  registerProcessorParameter("SmoothOn", "Smooth All Mesurement Sites in Fit", _SmoothOn, bool(false));

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
                             _maxChi2PerHit, float(1.e2));

  registerProcessorParameter("TrackSystemName",
                             "Name of the track fitting system to be used (KalTest, DDKalTest, aidaTT, ... )",
                             _trkSystemName, std::string("KalTest"));

  registerProcessorParameter("InitialTrackState",
                             "TrackState to use for initialization of the fit: -1: refit from hits [default], 1: AtIP, "
                             "2: AtFirstHit, 3: AtLastHit, 4:AtCalo",
                             _initialTrackState, int(-1));

  registerProcessorParameter("FitDirection", "Fit direction: -1: backward [default], +1: forward", _fitDirection,
                             int(-1));

  registerProcessorParameter("ParticleMass",
                             "particle mass that is used in the fit - default is the pion mass: 0.13957018 )", _mass,
                             double(0.13957018));
}

void RefitProcessor::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;

  // usually a good idea to
  printParameters();

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  double bFieldVec[3];
  theDetector.field().magneticField({0, 0, 0}, bFieldVec); // get the magnetic field vector from DD4hep
  _bField = bFieldVec[2] / dd4hep::tesla;                  // z component at (0,0,0)

  //----
  // set up the geometery needed for tracking

  _trksystem = MarlinTrk::Factory::createMarlinTrkSystem(_trkSystemName, 0, "");

  if (_trksystem == 0) {
    throw EVENT::Exception(std::string("  Cannot initialize MarlinTrkSystem of Type: ") + _trkSystemName);
  }

  _trksystem->setOption(IMarlinTrkSystem::CFG::useQMS, _MSOn);
  _trksystem->setOption(IMarlinTrkSystem::CFG::usedEdx, _ElossOn);
  _trksystem->setOption(IMarlinTrkSystem::CFG::useSmoothing, _SmoothOn);
  _trksystem->init();

  _n_run = 0;
  _n_evt = 0;
}

void RefitProcessor::processRunHeader(LCRunHeader*) { ++_n_run; }

void RefitProcessor::processEvent(LCEvent* evt) {
  // set the correct configuration for the tracking system for this event
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(_trksystem, _MSOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(_trksystem, _ElossOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);

  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG4) << "   processing event: " << _n_evt << std::endl;

  // get input collection and relations
  LCCollection* input_track_col = this->GetCollection(evt, _input_track_col_name);

  auto input_track_rels = this->GetRelations(evt, _input_track_rel_name);

  if (input_track_col != 0) {
    // establish the track collection that will be created
    LCCollectionVec* trackVec = new LCCollectionVec(LCIO::TRACK);

    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0);
    trkFlag.setBit(LCIO::TRBIT_HITS);
    trackVec->setFlag(trkFlag.getFlag());

    // establish the track relations collection that will be created
    auto trackRelNav = UTIL::LCRelationNavigator(LCIO::TRACK, LCIO::MCPARTICLE);

    int nTracks = input_track_col->getNumberOfElements();
    streamlog_out(DEBUG4) << "Processing input collection " << _input_track_col_name << " with " << nTracks
                          << " tracks\n";

    // loop over the input tacks and refit using KalTest
    for (int i = 0; i < nTracks; ++i) {
      Track* track_to_refit = dynamic_cast<Track*>(input_track_col->getElementAt(i));

      EVENT::TrackerHitVec trkHits = track_to_refit->getTrackerHits();

      std::vector<std::pair<float, EVENT::TrackerHit*>> r2_values;
      r2_values.reserve(trkHits.size());

      for (TrackerHitVec::iterator it = trkHits.begin(); it != trkHits.end(); ++it) {
        EVENT::TrackerHit* h = *it;
        float r2 = h->getPosition()[0] * h->getPosition()[0] + h->getPosition()[1] * h->getPosition()[1];
        r2_values.push_back(std::make_pair(r2, *it));
      }

      sort(r2_values.begin(), r2_values.end());

      trkHits.clear();
      trkHits.reserve(r2_values.size());

      UTIL::BitField64 cellID_encoder(lcio::LCTrackerCellID::encoding_string());

      for (std::vector<std::pair<float, EVENT::TrackerHit*>>::iterator it = r2_values.begin(); it != r2_values.end();
           ++it) {
        streamlog_out(DEBUG0) << " -- added tracker hit : " << *it->second << std::endl;

        trkHits.push_back(it->second);
      }

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

      bool fit_direction = ((_fitDirection < 0) ? IMarlinTrack::backward : IMarlinTrack::forward);

      streamlog_out(DEBUG1) << "TruthTracker::createTrack: fit direction used for fit (-1:backward,+1forward) : "
                            << _fitDirection << std::endl;

      MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();

      marlinTrk->setMass(_mass);

      TrackImpl* refittedTrack = new TrackImpl;

      try {
        int error = 0;

        if (_initialTrackState < 0) { // initialize the track from three hits

          // error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, refittedTrack, fit_direction, covMatrix,
          // _bField, _maxChi2PerHit);

          // call with empty pre_fit  -> should use default initialisation of the implementation, e.g.
          // use an internal pre fit in aidaTT
          error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, refittedTrack, fit_direction, 0, _bField,
                                                      _maxChi2PerHit);

        } else { // use the specified track state

          EVENT::TrackState* ts = const_cast<EVENT::TrackState*>(track_to_refit->getTrackState(_initialTrackState));

          if (!ts) {
            std::stringstream ess;
            ess << "  Could not get track state at " << _initialTrackState << " from track to refit ";
            throw EVENT::Exception(ess.str());
          }

          IMPL::TrackStateImpl pre_fit(*ts);
          pre_fit.setCovMatrix(covMatrix);

          error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, refittedTrack, fit_direction, &pre_fit,
                                                      _bField, _maxChi2PerHit);
        }

        if (error != IMarlinTrack::success || refittedTrack->getNdf() < 0) {
          streamlog_out(DEBUG6) << "in event << " << evt->getEventNumber() << " >> track refit returns error code "
                                << error << "; NDF = " << refittedTrack->getNdf()
                                << ". Number of hits = " << trkHits.size() << std::endl;
        }

      } catch (...) {
        streamlog_out(ERROR) << "RefitProcessor::processEvent: EVENT: << " << evt->getEventNumber()
                             << " >> exception caught and rethown. Track = " << track_to_refit << std::endl;

        delete refittedTrack;

        throw;
      }

      // fitting finished get hit in the fit

      std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
      std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;

      // remember the hits are ordered in the order in which they were fitted
      // here we are fitting inwards to the first is the last and vice verse

      marlinTrk->getHitsInFit(hits_in_fit);

      if (hits_in_fit.size() < 3) {
        streamlog_out(DEBUG6) << "RefitProcessor: Less than 3 hits in fit: number of hits =  " << trkHits.size()
                              << std::endl;
      }

      std::vector<TrackerHit*> all_hits;
      all_hits.reserve(300);

      for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
        all_hits.push_back(hits_in_fit[ihit].first);
      }

      //      UTIL::BitField64 cellID_encoder( lcio::LCTrackerCellID::encoding_string() ) ;

      MarlinTrk::addHitNumbersToTrack(refittedTrack, all_hits, true, cellID_encoder);

      marlinTrk->getOutliers(outliers);

      for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
        all_hits.push_back(outliers[ihit].first);
      }

      MarlinTrk::addHitNumbersToTrack(refittedTrack, all_hits, false, cellID_encoder);

      delete marlinTrk;

      int nhits_in_vxd = refittedTrack->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 1];
      int nhits_in_ftd = refittedTrack->subdetectorHitNumbers()[2 * lcio::ILDDetID::FTD - 1];
      int nhits_in_sit = refittedTrack->subdetectorHitNumbers()[2 * lcio::ILDDetID::SIT - 1];
      int nhits_in_tpc = refittedTrack->subdetectorHitNumbers()[2 * lcio::ILDDetID::TPC - 1];
      int nhits_in_set = refittedTrack->subdetectorHitNumbers()[2 * lcio::ILDDetID::SET - 1];

      streamlog_out(DEBUG3) << " Hit numbers for Track " << refittedTrack->id() << ": "
                            << " vxd hits = " << nhits_in_vxd << " ftd hits = " << nhits_in_ftd
                            << " sit hits = " << nhits_in_sit << " tpc hits = " << nhits_in_tpc
                            << " set hits = " << nhits_in_set << std::endl;

      if (nhits_in_vxd > 0)
        refittedTrack->setTypeBit(lcio::ILDDetID::VXD);
      if (nhits_in_ftd > 0)
        refittedTrack->setTypeBit(lcio::ILDDetID::FTD);
      if (nhits_in_sit > 0)
        refittedTrack->setTypeBit(lcio::ILDDetID::SIT);
      if (nhits_in_tpc > 0)
        refittedTrack->setTypeBit(lcio::ILDDetID::TPC);
      if (nhits_in_set > 0)
        refittedTrack->setTypeBit(lcio::ILDDetID::SET);

      trackVec->addElement(refittedTrack);

      // assign the relations previously assigned to the input tracks
      if (input_track_rels) {
        const auto& objVec = input_track_rels->getRelatedToObjects(track_to_refit);
        const auto& weights = input_track_rels->getRelatedToWeights(track_to_refit);

        for (unsigned int irel = 0; irel < objVec.size(); ++irel) {
          trackRelNav.addRelation(refittedTrack, objVec[irel], weights[irel]);
        }
      }
    }

    evt->addCollection(trackVec, _output_track_col_name);
    auto trackRelVec = trackRelNav.createLCCollection();
    evt->addCollection(trackRelVec, _output_track_rel_name);
  }
  ++_n_evt;
}

void RefitProcessor::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void RefitProcessor::end() {
  streamlog_out(DEBUG) << "RefitProcessor::end()  " << name() << " processed " << _n_evt << " events in " << _n_run
                       << " runs " << std::endl;
}

LCCollection* RefitProcessor::GetCollection(LCEvent* evt, std::string colName) {
  LCCollection* col = NULL;

  try {
    col = evt->getCollection(colName.c_str());
    streamlog_out(DEBUG4) << " --> " << colName.c_str() << " track collection found in event = " << col
                          << " number of elements " << col->getNumberOfElements() << std::endl;
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << " --> " << colName.c_str() << " collection absent in event" << std::endl;
  }

  return col;
}

std::unique_ptr<LCRelationNavigator> RefitProcessor::GetRelations(LCEvent* evt, std::string RelName) {
  std::unique_ptr<UTIL::LCRelationNavigator> nav = nullptr;

  try {
    nav.reset(new LCRelationNavigator(evt->getCollection(RelName.c_str())));
    streamlog_out(DEBUG2) << "RefitProcessor --> " << RelName << " track relation collection in event = " << nav.get()
                          << std::endl;
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG2) << "RefitProcessor --> " << RelName.c_str() << " track relation collection absent in event"
                          << std::endl;
  }

  return nav;
}
