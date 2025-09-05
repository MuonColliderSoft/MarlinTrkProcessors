#include "RefitFinal.h"

#include <marlin/Exceptions.h>
#include <marlin/Global.h>
#include <marlin/VerbosityLevels.h>

#include <MarlinTrk/Factory.h>
#include <MarlinTrk/IMarlinTrack.h>
#include <MarlinTrk/MarlinTrkUtils.h>

#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <UTIL/BitField64.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/Operators.h>

#include <algorithm>

using namespace lcio;
using namespace marlin;

RefitFinal aRefitFinal;

RefitFinal::RefitFinal() : Processor("RefitFinal") {
  // modify processor description
  _description = "Refit processor that calls finaliseLCIOTrack after taking "
                 "the trackstate from the existing track. No re-sorting of "
                 "hits is done";

  // register steering parameters: name, description, class-variable, default
  // value

  registerInputCollection(LCIO::TRACK, "InputTrackCollectionName", "Name of the input track collection",
                          _input_track_col_name, _input_track_col_name);

  registerInputCollection(LCIO::LCRELATION, "InputRelationCollectionName",
                          "Name of the input track to MCParticle relation collection", _input_track_rel_name,
                          _input_track_rel_name);

  registerOutputCollection(LCIO::TRACK, "OutputTrackCollectionName", "Name of the output track collection",
                           _output_track_col_name, _output_track_col_name);

  registerOutputCollection(LCIO::LCRELATION, "OutputRelationCollectionName",
                           "Refit Track to MCParticle relation collection Name", _output_track_rel_name,
                           _output_track_rel_name);

  registerProcessorParameter("MultipleScatteringOn", "Use MultipleScattering in Fit", _MSOn, bool(true));

  registerProcessorParameter("EnergyLossOn", "Use Energy Loss in Fit", _ElossOn, bool(true));

<<<<<<< HEAD

  registerProcessorParameter("SmoothOn", "Smooth All Mesurement Sites in Fit",
                             _SmoothOn, bool(false));

  registerProcessorParameter("Max_Chi2_Incr",
                             "maximum allowable chi2 increment when moving from one site to another",
=======
  registerProcessorParameter("SmoothOn", "Smooth All Mesurement Sites in Fit", _SmoothOn, bool(false));

  registerProcessorParameter("Max_Chi2_Incr", "maximum allowable chi2 increment when moving from one site to another",
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
                             _Max_Chi2_Incr, _Max_Chi2_Incr);

  registerProcessorParameter("ReferencePoint",
                             "Identifier of the reference point to use for the "
                             "fit initialisation, -1 means at 0 0 0",
                             _refPoint, _refPoint);

  registerProcessorParameter("extrapolateForward",
                             "if true extrapolation in the forward direction "
                             "(in-out), otherwise backward (out-in)",
                             _extrapolateForward, _extrapolateForward);

  registerProcessorParameter("MinClustersOnTrackAfterFit", "Final minimum number of track clusters",
                             _minClustersOnTrackAfterFit, int(4));
<<<<<<< HEAD

  registerProcessorParameter("ReducedChi2Cut", "Cut on the reduced chi square", _ReducedChi2Cut,
                             double(-1.));

  registerProcessorParameter("NHitsCuts", "Cuts on Nhits: <detID>,<detID>,... <N hits min> ", _NHitsCutsStr,
                             StringVec(0));

=======
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
}

void RefitFinal::init() {
  // usually a good idea to
  printParameters();

<<<<<<< HEAD
  // Extracting nHits cuts
  _NHitsCuts.resize( _NHitsCutsStr.size() / 2  ) ;

  unsigned i=0, index=0 ;
  while( i < _NHitsCutsStr.size() ) {
    std::vector<int> detIds;
    // Extracting the comma-separated list of detector IDs
    std::string split;
    std::istringstream ss(_NHitsCutsStr[ i++ ].c_str());
    while ( std::getline(ss, split, ',') )
      _NHitsCuts[index].detIDs.push_back(std::atoi(split.c_str()));
    // Storing the nHits value
    _NHitsCuts[index].nHits_min = std::atoi( _NHitsCutsStr[ i++ ].c_str() );

    streamlog_out(DEBUG9) << "nHits cut " << index << ":  " << _NHitsCuts[index].nHits_min << " from detIds: ";
    for (int detId : _NHitsCuts[index].detIDs) streamlog_out(DEBUG9) << " " << detId;
    streamlog_out(DEBUG9) << std::endl;

    ++index ;
  }

  _trksystem =
      MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
=======
  _trksystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90

  ///////////////////////////////

  _encoder = std::make_shared<UTIL::BitField64>(lcio::LCTrackerCellID::encoding_string());

  if (not _trksystem) {
    throw EVENT::Exception("Cannot initialize MarlinTrkSystem of Type: DDKalTest");
  }

  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useQMS, _MSOn);
  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, _ElossOn);
  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, _SmoothOn);
  _trksystem->init();

  _n_run = 0;
  _n_evt = 0;
}

void RefitFinal::processRunHeader(LCRunHeader*) { ++_n_run; }

<<<<<<< HEAD
void RefitFinal::processEvent(LCEvent *evt) {

  // set the correct configuration for the tracking system for this event
  MarlinTrk::TrkSysConfig< MarlinTrk::IMarlinTrkSystem::CFG::useQMS>       mson( _trksystem,  _MSOn ) ;
  MarlinTrk::TrkSysConfig< MarlinTrk::IMarlinTrkSystem::CFG::usedEdx>      elosson( _trksystem,_ElossOn) ;
  MarlinTrk::TrkSysConfig< MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon( _trksystem,_SmoothOn) ;
=======
void RefitFinal::processEvent(LCEvent* evt) {
  // set the correct configuration for the tracking system for this event
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(_trksystem, _MSOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(_trksystem, _ElossOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90

  ++_n_evt;

  // get input collection and relations
  LCCollection* input_track_col = this->GetCollection(evt, _input_track_col_name);
  if (not input_track_col) {
    return;
  }

  // establish the track collection that will be created
<<<<<<< HEAD
  LCCollectionVec *trackVec = new LCCollectionVec(LCIO::TRACK);
  UTIL::BitField64 encoder(lcio::LCTrackerCellID::encoding_string());
  encoder.reset(); // reset to 0
=======
  LCCollectionVec* trackVec = new LCCollectionVec(LCIO::TRACK);
  _encoder->reset();
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  trackVec->setFlag(trkFlag.getFlag());

  LCCollection* input_rel_col = this->GetCollection(evt, _input_track_rel_name);
  LCCollectionVec* trackRelationCollection = nullptr;
  std::shared_ptr<LCRelationNavigator> relation;

  if (not input_rel_col) {
    streamlog_out(DEBUG9) << "No input relation collection, not creating one either" << std::endl;
  } else {
    trackRelationCollection = new LCCollectionVec(LCIO::LCRELATION);
    relation = std::make_shared<LCRelationNavigator>(input_rel_col);
  }

  const int nTracks = input_track_col->getNumberOfElements();

  streamlog_out(MESSAGE4) << " Number of Tracks " << nTracks << std::endl;

  int counter = 0;
  std::map<int, int> hitInSubDet;
  // loop over the input tracks and refit
  for (int iTrack = 0; iTrack < nTracks; ++iTrack) {
    Track* track = static_cast<Track*>(input_track_col->getElementAt(iTrack));

    auto marlin_trk = std::unique_ptr<MarlinTrk::IMarlinTrack>(_trksystem->createTrack());
    EVENT::TrackerHitVec trkHits = track->getTrackerHits();

    streamlog_out(DEBUG5) << "---- track n = " << iTrack << "  n hits = " << trkHits.size() << std::endl;

    const int nHitsTrack = trkHits.size();

    hitInSubDet.clear();

    for (int iHit = 0; iHit < nHitsTrack && iHit < nHitsTrack; ++iHit) {
      marlin_trk->addHit(trkHits[iHit]);
      encoder.setValue(trkHits[iHit]->getCellID0()) ;
      ++hitInSubDet[encoder[UTIL::LCTrackerCellID::subdet()]];      
    }

    int init_status = FitInit2(track, marlin_trk.get());

    if (init_status != 0) {
      continue;
    }

<<<<<<< HEAD
    // Checking numbers of hits in subdetectors
    if (_NHitsCuts.size() > 0) {
      bool skipTrack(false);
      for (const NHitsCut& cut : _NHitsCuts) {
        // Summing the number of hits from each detector
        int nHits = 0;
        for (int detId : cut.detIDs) {
          nHits += hitInSubDet[detId];
        }
        // Checking if the cut is satisfied _before_ the fit
        if (nHits < cut.nHits_min) {
          streamlog_out(DEBUG5) << "Skip track " << track->id() << ": "
                                << " not enough detector hits: " << nHits << "/" << cut.nHits_min << std::endl;
        skipTrack = true;
        counter++;
        break;
        }
      }
      if (skipTrack) continue;
    } 

    streamlog_out(DEBUG4) << "Refit: Trackstate after initialisation\n"
                          << marlin_trk->toString() << std::endl;
=======
    streamlog_out(DEBUG4) << "Refit: Trackstate after initialisation\n" << marlin_trk->toString() << std::endl;
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90

    streamlog_out(DEBUG5) << "track initialised " << std::endl;

    int fit_status = marlin_trk->fit();

    streamlog_out(DEBUG4) << "RefitHit: Trackstate after fit()\n" << marlin_trk->toString() << std::endl;

    if (fit_status != 0) {
      continue;
    }

    auto lcio_trk = std::unique_ptr<IMPL::TrackImpl>(new IMPL::TrackImpl());

    const bool fit_direction = MarlinTrk::IMarlinTrack::forward;
    int return_code = finaliseLCIOTrack(marlin_trk.get(), lcio_trk.get(), trkHits, fit_direction);

    if (return_code != MarlinTrk::IMarlinTrack::success) {
      streamlog_out(DEBUG3) << "finaliseLCIOTrack failed" << std::endl;
      continue;
    }

    streamlog_out(DEBUG5) << " *** created finalized LCIO track - return code " << return_code << std::endl
                          << *lcio_trk << std::endl;

    // fit finished - get hits in the fit
    std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
    std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;

    // remember the hits are ordered in the order in which they were fitted

    marlin_trk->getHitsInFit(hits_in_fit);

    if (int(hits_in_fit.size()) < _minClustersOnTrackAfterFit) {
<<<<<<< HEAD
      streamlog_out(DEBUG5) << "Less than " << _minClustersOnTrackAfterFit << " hits in fit: Track "
=======
      streamlog_out(DEBUG3) << "Less than " << _minClustersOnTrackAfterFit
                            << " hits in fit: Track "
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
                               "Discarded. Number of hits =  "
                            << trkHits.size() << std::endl;
      continue;
    }

    marlin_trk->getOutliers(outliers);

    std::vector<TrackerHit*> all_hits;
    all_hits.reserve(hits_in_fit.size() + outliers.size());

    for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      all_hits.push_back(hits_in_fit[ihit].first);
    }

    for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      all_hits.push_back(outliers[ihit].first);
    }

    UTIL::BitField64 encoder2(lcio::LCTrackerCellID::encoding_string());
    encoder2.reset(); // reset to 0
    MarlinTrk::addHitNumbersToTrack(lcio_trk.get(), all_hits, false, encoder2);
    MarlinTrk::addHitNumbersToTrack(lcio_trk.get(), hits_in_fit, true, encoder2);

    streamlog_out(DEBUG4) << "processEvent: Hit numbers for track " << lcio_trk->id() << ":  " << std::endl;
    int detID = 0;
    for (size_t ip = 0; ip < lcio_trk->subdetectorHitNumbers().size(); ip = ip + 2) {
      detID++;
      streamlog_out(DEBUG4) << "  det id " << detID << " , nhits in track = " << lcio_trk->subdetectorHitNumbers()[ip]
                            << " , nhits in fit = " << lcio_trk->subdetectorHitNumbers()[ip + 1] << std::endl;
      if (lcio_trk->subdetectorHitNumbers()[ip] > 0)
        lcio_trk->setTypeBit(detID);
    }

    auto lcioTrkPtr = lcio_trk.release();

    // if required apply the ReducedChi2 cut
    if (_ReducedChi2Cut > 0. && lcioTrkPtr->getChi2()/lcioTrkPtr->getNdf() > _ReducedChi2Cut) {
      streamlog_out(DEBUG5) << "Skip track " << lcioTrkPtr->id() << ": "
                            << "Chi2/ndof " << lcioTrkPtr->getChi2()/lcioTrkPtr->getNdf() << std::endl;
      counter++; 
      continue;
    }

    trackVec->addElement(lcioTrkPtr);

    if (input_rel_col) {
      auto mcParticleVec = relation->getRelatedToObjects(track);
      auto weightVec = relation->getRelatedToWeights(track);
      for (size_t i = 0; i < mcParticleVec.size(); ++i) {
        LCRelationImpl* relationTrack = new LCRelationImpl(lcioTrkPtr, mcParticleVec[i], weightVec[i]);
        trackRelationCollection->addElement(relationTrack);
      }
    }

  } // for loop to the tracks
  streamlog_out(MESSAGE4) << " Final number of tracks " << trackVec->getNumberOfElements() 
                          << " Skipped tracks: " << counter << std::endl;


  evt->addCollection(trackVec, _output_track_col_name);
  if (input_rel_col) {
    evt->addCollection(trackRelationCollection, _output_track_rel_name);
  }
}

void RefitFinal::check(LCEvent*) {}

void RefitFinal::end() {}

LCCollection* RefitFinal::GetCollection(LCEvent* evt, std::string colName) {
  LCCollection* col = nullptr;

  try {
    col = evt->getCollection(colName.c_str());
    streamlog_out(DEBUG3) << " --> " << colName.c_str() << " track collection found in event = " << col
                          << " number of elements " << col->getNumberOfElements() << std::endl;
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG3) << " --> " << colName.c_str() << " collection absent in event" << std::endl;
  }

  return col;
}

int RefitFinal::FitInit2(Track* track, MarlinTrk::IMarlinTrack* marlinTrk) {
  TrackStateImpl trackState;

  if (_refPoint == -1) {
    trackState = TrackStateImpl(TrackState::AtOther, track->getD0(), track->getPhi(), track->getOmega(), track->getZ0(),
                                track->getTanLambda(), track->getCovMatrix(), track->getReferencePoint());
  } else {
    const TrackState* trackAtHit = track->getTrackState(_refPoint);
    if (not trackAtHit) {
      streamlog_out(ERROR) << "Cannot find trackstate for " << _refPoint << std::endl;
      return MarlinTrk::IMarlinTrack::error;
    }
    trackState = TrackStateImpl(*trackAtHit);
  }

  const bool direction = _extrapolateForward ? MarlinTrk::IMarlinTrack::forward : MarlinTrk::IMarlinTrack::backward;
  marlinTrk->initialise(trackState, _bField, direction);

  return MarlinTrk::IMarlinTrack::success;
}

