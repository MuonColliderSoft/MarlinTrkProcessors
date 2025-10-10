#include "ClonesAndSplitTracksFinder.h"
#include "HitsSorterAndDebugger.h"

#include <marlin/Exceptions.h>
#include <marlin/Global.h>
#include <marlin/VerbosityLevels.h>

#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/MarlinTrkDiagnostics.h"
#include <MarlinTrk/Factory.h>
#include <MarlinTrk/IMarlinTrack.h>
#include <MarlinTrk/MarlinTrkUtils.h>

#include <marlinutil/GeometryUtil.h>

#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <UTIL/BitField64.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/Operators.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/SurfaceManager.h"

#include <algorithm>

using namespace lcio;
using namespace marlin;
using namespace std;

ClonesAndSplitTracksFinder aClonesAndSplitTracksFinder;

ClonesAndSplitTracksFinder::ClonesAndSplitTracksFinder() : Processor("ClonesAndSplitTracksFinder") {
  // modify processor description
  _description =
      "ClonesAndSplitTracksFinder takes the track collection, checks for doubles and merges them in an output track "
      "collection";

  // register steering parameters: name, description, class-variable, default
  // value

  registerInputCollection(LCIO::TRACK, "InputTrackCollectionName", "Name of the input track collection",
                          _input_track_col_name, std::string("SiTracks"));

  registerOutputCollection(LCIO::TRACK, "OutputTrackCollectionName", "Name of the output track collection",
                           _output_track_col_name, std::string("SiTracksMerged"));

  registerProcessorParameter("MultipleScatteringOn", "Use MultipleScattering in Fit", _MSOn, bool(true));

  registerProcessorParameter("EnergyLossOn", "Use Energy Loss in Fit", _ElossOn, bool(true));

  registerProcessorParameter("SmoothOn", "Smooth All Measurement Sites in Fit", _SmoothOn, bool(false));

  registerProcessorParameter("extrapolateForward",
                             "if true extrapolation in the forward direction "
                             "(in-out), otherwise backward (out-in)",
                             _extrapolateForward, bool(true));

  registerProcessorParameter("minTrackPt", "minimum track pt for merging (in GeV/c)", _minPt, double(1.0));

  registerProcessorParameter("maxSignificanceTheta", "maximum significance separation in tanLambda",
                             _maxSignificanceTheta, double(3.0));

  registerProcessorParameter("maxSignificancePhi", "maximum significance separation in phi", _maxSignificancePhi,
                             double(3.0));

  registerProcessorParameter("maxSignificancePt", "maximum significance separation in pt", _maxSignificancePt,
                             double(2.0));

  registerProcessorParameter("mergeSplitTracks", "if true, the merging of split tracks is performed", _mergeSplitTracks,
                             bool(false));
}

void ClonesAndSplitTracksFinder::init() {
  // usually a good idea to
  printParameters();

  _trksystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");

  _magneticField = MarlinUtil::getBzAtOrigin();
  ///////////////////////////////

  _encoder = std::make_shared<UTIL::BitField64>(lcio::LCTrackerCellID::encoding_string());

  if (not _trksystem) {
    throw EVENT::Exception("Cannot initialize MarlinTrkSystem of Type: DDKalTest");
  }

  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useQMS, _MSOn);
  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, _ElossOn);
  _trksystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, _SmoothOn);
  _trksystem->init();

  // Put default values for track fitting
  _initialTrackError_d0 = 1.e6;
  _initialTrackError_phi0 = 1.e2;
  _initialTrackError_omega = 1.e-4;
  _initialTrackError_z0 = 1.e6;
  _initialTrackError_tanL = 1.e2;
  _maxChi2perHit = 1.e2;

  _n_run = 0;
  _n_evt = 0;
}

void ClonesAndSplitTracksFinder::processRunHeader(LCRunHeader*) { ++_n_run; }

void ClonesAndSplitTracksFinder::processEvent(LCEvent* evt) {
  // set the correct configuration for the tracking system for this event
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(_trksystem, _MSOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(_trksystem, _ElossOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);

  ++_n_evt;

  // get input collection and relations
  LCCollection* input_track_col = this->GetCollection(evt, _input_track_col_name);
  if (not input_track_col) {
    return;
  }
  const int nTracks = input_track_col->getNumberOfElements();
  streamlog_out(DEBUG5) << " >> ClonesAndSplitTracksFinder starts with " << nTracks << " tracks." << std::endl;

  // establish the track collection that will be created
  auto trackVec = std::unique_ptr<LCCollectionVec>(new LCCollectionVec(LCIO::TRACK));
  _encoder->reset();
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  trackVec->setFlag(trkFlag.getFlag());

  //------------
  // FIRST STEP: REMOVE CLONES
  //------------

  EVENT::TrackVec tracksWithoutClones;
  removeClones(tracksWithoutClones, input_track_col);
  const int ntracksWithoutClones = tracksWithoutClones.size();
  streamlog_out(DEBUG5) << " >> ClonesAndSplitTracksFinder found " << ntracksWithoutClones << " tracks without clones."
                        << std::endl;

  if (_mergeSplitTracks && ntracksWithoutClones > 1) {
    streamlog_out(DEBUG5) << " Try to merge tracks ..." << std::endl;

    //------------
    // SECOND STEP: MERGE TRACKS
    //------------

    mergeSplitTracks(trackVec, input_track_col, tracksWithoutClones);
  } else {
    streamlog_out(DEBUG5) << " Not even try to merge tracks ..." << std::endl;
    for (UInt_t iTrk = 0; iTrk < tracksWithoutClones.size(); iTrk++) {
      TrackImpl* trackFinal = new TrackImpl;
      fromTrackToTrackImpl(tracksWithoutClones.at(iTrk), trackFinal);
      trackVec->addElement(trackFinal);
    }
  }

  evt->addCollection(trackVec.release(), _output_track_col_name);
}

void ClonesAndSplitTracksFinder::check(LCEvent*) {}

void ClonesAndSplitTracksFinder::end() {}

LCCollection* ClonesAndSplitTracksFinder::GetCollection(LCEvent* evt, std::string colName) {
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

// Function to check if two KDtracks contain several hits in common
int ClonesAndSplitTracksFinder::overlappingHits(const Track* track1, const Track* track2) {
  int nHitsInCommon = 0;
  const EVENT::TrackerHitVec trackVec1 = track1->getTrackerHits();
  const EVENT::TrackerHitVec trackVec2 = track2->getTrackerHits();
  for (size_t hit = 0; hit < trackVec1.size(); hit++) {
    if (std::find(trackVec2.begin(), trackVec2.end(), trackVec1.at(hit)) != trackVec2.end())
      nHitsInCommon++;
  }
  return nHitsInCommon;
}

void ClonesAndSplitTracksFinder::fromTrackToTrackImpl(const Track* track, TrackImpl*& trackFinal) {
  const TrackState* ts_atOther = 0;
  ts_atOther = track->getTrackState(TrackState::AtOther);
  if (ts_atOther)
    trackFinal->addTrackState(new TrackStateImpl(*ts_atOther));
  const TrackState* ts_atIP = 0;
  ts_atIP = track->getTrackState(TrackState::AtIP);
  if (ts_atIP)
    trackFinal->addTrackState(new TrackStateImpl(*ts_atIP));
  const TrackState* ts_atFirstHit = 0;
  ts_atFirstHit = track->getTrackState(TrackState::AtFirstHit);
  if (ts_atFirstHit)
    trackFinal->addTrackState(new TrackStateImpl(*ts_atFirstHit));
  const TrackState* ts_atLastHit = 0;
  ts_atLastHit = track->getTrackState(TrackState::AtLastHit);
  if (ts_atLastHit)
    trackFinal->addTrackState(new TrackStateImpl(*ts_atLastHit));
  const TrackState* ts_atCalorimeter = 0;
  ts_atCalorimeter = track->getTrackState(TrackState::AtCalorimeter);
  if (ts_atCalorimeter)
    trackFinal->addTrackState(new TrackStateImpl(*ts_atCalorimeter));
  const TrackState* ts_atVertex = 0;
  ts_atVertex = track->getTrackState(TrackState::AtVertex);
  if (ts_atVertex)
    trackFinal->addTrackState(new TrackStateImpl(*ts_atVertex));
  const TrackState* ts_atLastLocation = 0;
  ts_atLastLocation = track->getTrackState(TrackState::LastLocation);
  if (ts_atLastLocation)
    trackFinal->addTrackState(new TrackStateImpl(*ts_atLastLocation));

  for (UInt_t i = 0; i < (track->getTrackerHits()).size(); i++) {
    trackFinal->addHit(track->getTrackerHits().at(i));
  }
  trackFinal->setRadiusOfInnermostHit(track->getRadiusOfInnermostHit());
  trackFinal->setChi2(track->getChi2());
  trackFinal->setNdf(track->getNdf());
  trackFinal->setdEdx(track->getdEdx());
  trackFinal->setdEdxError(track->getdEdxError());
  trackFinal->subdetectorHitNumbers() = track->getSubdetectorHitNumbers();
}

void ClonesAndSplitTracksFinder::removeClones(EVENT::TrackVec& tracksWithoutClones, LCCollection*& input_track_col) {
  streamlog_out(DEBUG8) << "ClonesAndSplitTracksFinder::removeClones " << std::endl;
  const int nTracks = input_track_col->getNumberOfElements();

  // loop over the input tracks

  std::multimap<int, std::pair<int, Track*>> candidateClones;

  for (int iTrack = 0; iTrack < nTracks; ++iTrack) { // first loop over tracks
    int countClones = 0;
    Track* track_i = static_cast<Track*>(input_track_col->getElementAt(iTrack));

    for (int jTrack = 0; jTrack < nTracks; ++jTrack) { // second loop over tracks

      Track* track_j = static_cast<Track*>(input_track_col->getElementAt(jTrack));
      if (track_i != track_j) { // track1 != track2

        const unsigned int nOverlappingHits = overlappingHits(track_i, track_j);
        if (nOverlappingHits >= 2) { // clones
          countClones++;
          Track* bestTrack;
          bestInClones(track_i, track_j, nOverlappingHits, bestTrack);
          candidateClones.insert(make_pair(iTrack, make_pair(jTrack, bestTrack)));
        } else {
          continue;
        }
      }

    } // end second track loop

    if (countClones == 0) {
      tracksWithoutClones.push_back(track_i);
    }

  } // end first track loop

  filterClonesAndMergedTracks(candidateClones, input_track_col, tracksWithoutClones, true);
}

void ClonesAndSplitTracksFinder::mergeSplitTracks(std::unique_ptr<LCCollectionVec>& trackVec,
                                                  LCCollection*& input_track_col,
                                                  EVENT::TrackVec& tracksWithoutClones) {
  streamlog_out(DEBUG8) << "ClonesAndSplitTracksFinder::mergeSplitTracks " << std::endl;

  std::multimap<int, std::pair<int, Track*>> mergingCandidates;
  std::set<int> iter_duplicates;

  for (UInt_t iTrack = 0; iTrack < tracksWithoutClones.size(); ++iTrack) {
    int countMergingPartners = 0;
    bool toBeMerged = false;
    Track* track_i = static_cast<Track*>(tracksWithoutClones.at(iTrack));

    double pt_i = 0.3 * _magneticField / (fabs(track_i->getOmega()) * 1000.);
    double theta_i = (M_PI / 2 - atan(track_i->getTanLambda())) * 180. / M_PI;
    double phi_i = track_i->getPhi() * 180. / M_PI;

    // Merge only tracks with min pt
    // Try to avoid merging loopers for now
    if (pt_i < _minPt) {
      streamlog_out(DEBUG5) << " Track #" << iTrack << ": pt = " << pt_i << ", theta = " << theta_i
                            << ", phi = " << phi_i << std::endl;
      streamlog_out(DEBUG5) << " Track #" << iTrack << " does not fulfil min pt requirement." << std::endl;
      streamlog_out(DEBUG5) << " TRACK STORED" << std::endl;

      TrackImpl* trackFinal = new TrackImpl;
      fromTrackToTrackImpl(track_i, trackFinal);
      trackVec->addElement(trackFinal);
      continue;
    }

    for (UInt_t jTrack = iTrack + 1; jTrack < tracksWithoutClones.size(); ++jTrack) {
      Track* track_j = static_cast<Track*>(tracksWithoutClones.at(jTrack));
      bool isCloseInTheta = false, isCloseInPhi = false, isCloseInPt = false;

      if (track_j != track_i) {
        double pt_j = 0.3 * _magneticField / (fabs(track_j->getOmega() * 1000.));
        double theta_j = (M_PI / 2 - atan(track_j->getTanLambda())) * 180. / M_PI;
        double phi_j = track_j->getPhi() * 180. / M_PI;
        streamlog_out(DEBUG5) << " Track #" << iTrack << ": pt = " << pt_i << ", theta = " << theta_i
                              << ", phi = " << phi_i << std::endl;
        streamlog_out(DEBUG5) << " Track #" << jTrack << ": pt = " << pt_j << ", theta = " << theta_j
                              << ", phi = " << phi_j << std::endl;

        if (pt_j < _minPt) {
          streamlog_out(DEBUG5) << " Track #" << jTrack << " does not fulfil min pt requirement. Skip. " << std::endl;
          continue;
        }

        double significanceTanLambda = calculateSignificanceTanLambda(track_i, track_j);
        double significancePhi = calculateSignificancePhi(track_i, track_j);
        double significancePt = calculateSignificancePt(track_i, track_j);

        streamlog_out(DEBUG5) << " -> tanLambda significance = " << significanceTanLambda << " with cut at "
                              << _maxSignificanceTheta << std::endl;
        if (significanceTanLambda < _maxSignificanceTheta) {
          isCloseInTheta = true;
          streamlog_out(DEBUG5) << " Tracks are close in theta " << std::endl;
        }

        streamlog_out(DEBUG5) << " -> phi significance = " << significancePhi << " with cut at " << _maxSignificancePhi
                              << std::endl;
        if (significancePhi < _maxSignificancePhi) {
          isCloseInPhi = true;
          streamlog_out(DEBUG5) << " Tracks are close in phi " << std::endl;
        }

        streamlog_out(DEBUG5) << " -> pt significance = " << significancePt << " with cut at " << _maxSignificancePt
                              << std::endl;
        if (significancePt < _maxSignificancePt) {
          isCloseInPt = true;
          streamlog_out(DEBUG5) << " Tracks are close in pt  " << std::endl;
        }

        if (streamlog::out.write<streamlog::DEBUG5>()) {
          streamlog_out(DEBUG5) << " Track #" << iTrack << ": " << std::endl;
          printHits(track_i);
          streamlog_out(DEBUG5) << " Track #" << jTrack << ": " << std::endl;
          printHits(track_j);
        }

        toBeMerged = isCloseInTheta && isCloseInPhi && isCloseInPt;

        if (toBeMerged) { // merging, refitting, storing in a container of mergingCandidates (multimap <*track1,
                          // pair<*track2,*trackMerged>>)
          EVENT::Track* lcioTrkPtr = nullptr;
          mergeAndFit(track_i, track_j, lcioTrkPtr);
          if (not lcioTrkPtr) {
            continue;
          }
          mergingCandidates.insert(make_pair(iTrack, make_pair(jTrack, lcioTrkPtr)));
          countMergingPartners++;
          iter_duplicates.insert(iTrack);
          iter_duplicates.insert(jTrack);
        } else { // no merging conditions met
          continue;
        }
      }

    } // end loop on jTracks

    // Track was already found as duplicate
    const bool is_in = iter_duplicates.find(iTrack) != iter_duplicates.end();
    if (countMergingPartners == 0 && !is_in) { // if track_i has no merging partner, store it in the output vec
      streamlog_out(DEBUG5) << " Track #" << iTrack << " has no merging partners, so TRACK STORED." << std::endl;

      TrackImpl* trackFinal = new TrackImpl;
      fromTrackToTrackImpl(track_i, trackFinal);
      trackVec->addElement(trackFinal);

    } else {
      streamlog_out(DEBUG5) << " TRACK NOT STORED" << std::endl;
    }
    if (countMergingPartners != 0)
      streamlog_out(DEBUG5) << " possible merging partners for track #" << iTrack << " are = " << countMergingPartners
                            << std::endl;

  } // end loop on iTracks

  EVENT::TrackVec finalTracks;
  filterClonesAndMergedTracks(mergingCandidates, input_track_col, finalTracks, false);

  for (UInt_t iTrk = 0; iTrk < finalTracks.size(); iTrk++) {
    streamlog_out(DEBUG5) << " TRACK STORED" << std::endl;
    TrackImpl* trackFinal = new TrackImpl;
    fromTrackToTrackImpl(finalTracks.at(iTrk), trackFinal);
    trackVec->addElement(trackFinal);
  }
}

double ClonesAndSplitTracksFinder::calculateSignificancePt(const Track* first, const Track* second) {
  float omegaFirst = first->getOmega();
  float omegaSecond = second->getOmega();

  double ptFirst = 0.3 * _magneticField / (fabs(first->getOmega() * 1000.));
  double ptSecond = 0.3 * _magneticField / (fabs(second->getOmega() * 1000.));

  const float sigmaPOverPFirst = sqrt(first->getCovMatrix()[5]) / fabs(omegaFirst);
  const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5]) / fabs(omegaSecond);
  const float sigmaPtFirst = ptFirst * sigmaPOverPFirst;
  const float sigmaPtSecond = ptSecond * sigmaPOverPSecond;

  const double significance = calculateSignificance(ptFirst, ptSecond, sigmaPtFirst, sigmaPtSecond);

  return significance;
}

double ClonesAndSplitTracksFinder::calculateSignificancePhi(const Track* first, const Track* second) {
  float phiFirst = first->getPhi();
  float phiSecond = second->getPhi();
  float deltaPhi = (M_PI - std::abs(std::abs(phiFirst - phiSecond) - M_PI));

  const float sigmaPhiFirst = sqrt(first->getCovMatrix()[2]);
  const float sigmaPhiSecond = sqrt(second->getCovMatrix()[2]);

  const double significance = calculateSignificance(deltaPhi, 0.0, sigmaPhiFirst, sigmaPhiSecond);

  return significance;
}

double ClonesAndSplitTracksFinder::calculateSignificanceTanLambda(const Track* first, const Track* second) {
  float tanLambdaFirst = first->getTanLambda();
  float tanLambdaSecond = second->getTanLambda();

  const float sigmaTanLambdaFirst = sqrt(first->getCovMatrix()[14]);
  const float sigmaTanLambdaSecond = sqrt(second->getCovMatrix()[14]);

  const double significance =
      calculateSignificance(tanLambdaFirst, tanLambdaSecond, sigmaTanLambdaFirst, sigmaTanLambdaSecond);

  return significance;
}

double ClonesAndSplitTracksFinder::calculateSignificance(const double firstPar, const double secondPar,
                                                         const double firstPar_sigma, const double secondPar_sigma) {
  const float delta = fabs(firstPar - secondPar);
  const float sigmaDelta = sqrt(firstPar_sigma * firstPar_sigma + secondPar_sigma * secondPar_sigma);

  return delta / sigmaDelta;
}

void ClonesAndSplitTracksFinder::filterClonesAndMergedTracks(std::multimap<int, std::pair<int, Track*>>& candidates,
                                                             LCCollection*& inputTracks, TrackVec& trackVecFinal,
                                                             bool clones) {
  std::vector<TrackerHitVec> savedHitVec;

  for (const auto& iter : candidates) {
    int track_a_id = iter.first;
    int track_b_id = iter.second.first;
    Track* track_final = iter.second.second;
    int countConnections = candidates.count(track_a_id);
    bool multiConnection = (countConnections > 1);

    if (!multiConnection) { // if only 1 connection

      if (clones) { // clones: compare the track pointers
        auto it_trk = find(trackVecFinal.begin(), trackVecFinal.end(), track_final);
        if (it_trk != trackVecFinal.end()) { // if the track is already there, do nothing
          continue;
        }
        trackVecFinal.push_back(track_final);
      } else { // mergeable tracks: compare the sets of tracker hits

        TrackerHitVec track_final_hits = track_final->getTrackerHits();
        bool toBeSaved = true;

        for (const auto& hitsVec : savedHitVec) {
          if (equal(hitsVec.begin(), hitsVec.end(), track_final_hits.begin())) {
            toBeSaved = false;
            break;
          }
        }

        if (toBeSaved) {
          savedHitVec.push_back(track_final_hits);
          trackVecFinal.push_back(track_final);
        } else {
          delete track_final;
        }
      }

    } else { // if more than 1 connection, clones and mergeable tracks have to be treated a little different

      if (clones) { // clones

        // look at the elements with equal range. If their bestTrack is the same, store it (if not already in). If their
        // bestTrack is different, don't store it
        auto ret =
            candidates.equal_range(track_a_id); // a std::pair of iterators on the multimap [
                                                // std::pair<std::multimap<Track*,std::pair<Track*,Track*>>::iterator,
                                                // std::multimap<Track*,std::pair<Track*,Track*>>::iterator> ]
        TrackVec bestTracksMultiConnections;
        for (std::multimap<int, std::pair<int, Track*>>::iterator it = ret.first; it != ret.second; ++it) {
          Track* track_best = it->second.second;
          bestTracksMultiConnections.push_back(track_best);
        }
        if (std::adjacent_find(bestTracksMultiConnections.begin(), bestTracksMultiConnections.end(),
                               std::not_equal_to<Track*>()) ==
            bestTracksMultiConnections.end()) { // one best track with the same track key
          auto it_trk = find(trackVecFinal.begin(), trackVecFinal.end(), bestTracksMultiConnections.at(0));
          if (it_trk != trackVecFinal.end()) { // if the track is already there, do nothing
            continue;
          }
          trackVecFinal.push_back(bestTracksMultiConnections.at(0));

        } else { // multiple best tracks with the same track key
          continue;
        }

      } // end of clones

      else {                // mergeable tracks -- at the moment they are all stored (very rare anyways)
        delete track_final; // not using the mergedTracks, so delete it

        Track* track_a = static_cast<Track*>(inputTracks->getElementAt(track_a_id));
        Track* track_b = static_cast<Track*>(inputTracks->getElementAt(track_b_id));

        auto trk1 = find(trackVecFinal.begin(), trackVecFinal.end(), track_a);

        if (trk1 != trackVecFinal.end()) { // if the track1 is already there
          continue;
        }
        // otherwise store the two tracks
        trackVecFinal.push_back(track_a);
        trackVecFinal.push_back(track_b);

      } // end of mergeable tracks
    }
  }
}

void ClonesAndSplitTracksFinder::mergeAndFit(Track* track_i, Track* track_j, Track*& lcioTrkPtr) {
  streamlog_out(DEBUG8) << "ClonesAndSplitTracksFinder::mergeAndFit " << std::endl;
  EVENT::TrackerHitVec trkHits_i = track_i->getTrackerHits();
  EVENT::TrackerHitVec trkHits_j = track_j->getTrackerHits();

  EVENT::TrackerHitVec trkHits;
  for (UInt_t iHits = 0; iHits < trkHits_i.size(); iHits++) {
    trkHits.push_back(trkHits_i.at(iHits));
  }
  // Remove common hits while filling for the second track
  for (UInt_t jHits = 0; jHits < trkHits_j.size(); jHits++) {
    if (std::find(trkHits.begin(), trkHits.end(), trkHits_j.at(jHits)) != trkHits.end()) {
      streamlog_out(DEBUG8) << " This hit is already in the track" << std::endl;
      continue;
    } else {
      trkHits.push_back(trkHits_j.at(jHits));
    }
  }
  std::sort(trkHits.begin(), trkHits.end(), sort_by_radius);
  if (streamlog::out.write<streamlog::DEBUG5>()) {
    streamlog_out(DEBUG8) << " Hits in track to be merged: " << std::endl;
    printHits(trkHits);
  }

  auto mergedTrack = std::unique_ptr<TrackImpl>(new TrackImpl);

  auto marlin_trk = std::unique_ptr<MarlinTrk::IMarlinTrack>(_trksystem->createTrack());

  // Make an initial covariance matrix with very broad default values
  EVENT::FloatVec covMatrix(15, 0);          // Size 15, filled with 0s
  covMatrix[0] = (_initialTrackError_d0);    // sigma_d0^2
  covMatrix[2] = (_initialTrackError_phi0);  // sigma_phi0^2
  covMatrix[5] = (_initialTrackError_omega); // sigma_omega^2
  covMatrix[9] = (_initialTrackError_z0);    // sigma_z0^2
  covMatrix[14] = (_initialTrackError_tanL); // sigma_tanl^2

  const bool direction = _extrapolateForward ? MarlinTrk::IMarlinTrack::forward : MarlinTrk::IMarlinTrack::backward;

  int fit_status = MarlinTrk::createFinalisedLCIOTrack(marlin_trk.get(), trkHits, mergedTrack.get(), direction,
                                                       covMatrix, _magneticField, _maxChi2perHit);

  if (fit_status != 0) {
    streamlog_out(DEBUG4) << "Fit failed with error status " << fit_status << std::endl;
    return;
  }
  streamlog_out(DEBUG8) << " >> Fit not failed ! " << std::endl;

  // fit finished - get hits in the fit
  std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
  std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;

  // remember the hits are ordered in the order in which they were fitted

  marlin_trk->getHitsInFit(hits_in_fit);
  if (hits_in_fit.size() < 3) {
    streamlog_out(DEBUG4) << "Less than 3 hits in fit: Track discarded. Number of hits =  " << trkHits.size()
                          << std::endl;
    return;
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
  MarlinTrk::addHitNumbersToTrack(mergedTrack.get(), all_hits, false, encoder2);
  MarlinTrk::addHitNumbersToTrack(mergedTrack.get(), hits_in_fit, true, encoder2);

  if (streamlog::out.write<streamlog::DEBUG5>()) {
    streamlog_out(DEBUG5) << " Merged track : " << std::endl;
    printHits(&*mergedTrack);
  }

  lcioTrkPtr = mergedTrack.release();
}

void ClonesAndSplitTracksFinder::bestInClones(Track* track_a, Track* track_b, int nOverlappingHits, Track*& bestTrack) {
  // This function compares two tracks which have a certain number of overlapping hits and returns the best track
  // The best track is chosen based on length (in terms of number of hits) and chi2/ndf requirements
  // In general, the longest track is preferred. When clones have same length, the one with best chi2/ndf is chosen

  TrackerHitVec trackerHit_a = track_a->getTrackerHits();
  TrackerHitVec trackerHit_b = track_b->getTrackerHits();

  int trackerHit_a_size = trackerHit_a.size();
  int trackerHit_b_size = trackerHit_b.size();

  double b_chi2 = track_b->getChi2() / track_b->getNdf();
  double a_chi2 = track_a->getChi2() / track_a->getNdf();

  if (nOverlappingHits == trackerHit_a_size) { // if the second track is the first track + segment
    bestTrack = track_b;
  } else if (nOverlappingHits == trackerHit_b_size) { // if the second track is a subtrack of the first track
    bestTrack = track_a;
  } else if (trackerHit_b_size == trackerHit_a_size) { // if the two tracks have the same length
    if (b_chi2 <= a_chi2) {
      bestTrack = track_b;
    } else {
      bestTrack = track_a;
    }
  } else if (trackerHit_b_size > trackerHit_a_size) { // if the second track is longer
    bestTrack = track_b;
  } else if (trackerHit_b_size < trackerHit_a_size) { // if the second track is shorter
    bestTrack = track_a;
  }
}
