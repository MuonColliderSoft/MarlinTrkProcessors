#include "ExtrToTracker.h"

#include <IMPL/TrackerHitPlaneImpl.h>

#include "MarlinTrk/MarlinTrkUtils.h"

#include "UTIL/Operators.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include "DDRec/DetectorData.h"

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitField64.h>
#include <UTIL/BitSet32.h>
#include <UTIL/ILDConf.h>

#include <algorithm>

// CxxUtils/
#include "fpcompare.h"

using namespace lcio;
using namespace marlin;
using namespace MarlinTrk;

using namespace dd4hep::rec;

//------------------------------------------------------------------------------------------

struct PtSort { // sort tracks wtr to pt - largest first
  inline bool operator()(const lcio::LCObject* l, const lcio::LCObject* r) {
    return CxxUtils::fpcompare::less(std::abs(((const lcio::Track*)l)->getOmega()),
                                     std::abs(((const lcio::Track*)r)->getOmega()));
  }
};
//------------------------------------------------------------------------------------------

struct InversePtSort { // sort tracks wtr to pt - smallest first
  inline bool operator()(const lcio::LCObject* l, const lcio::LCObject* r) {
    return CxxUtils::fpcompare::greater(std::abs(((const lcio::Track*)l)->getOmega()),
                                        std::abs(((const lcio::Track*)r)->getOmega()));
  }
};

//------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------

struct ZSort { // sort track segments wtr to Z - smallest first
  inline bool operator()(const lcio::LCObject* l, const lcio::LCObject* r) {
    return CxxUtils::fpcompare::less(std::abs(((const lcio::Track*)l)->getZ0()),
                                     std::abs(((const lcio::Track*)r)->getZ0()));
  }
};

//------------------------------------------------------------------------------------------

ExtrToTracker aExtrToTracker;

ExtrToTracker::ExtrToTracker() : Processor("ExtrToTracker") {
  // modify processor description
  _description = "ExtrToTracker refits an input VXD track collection and uses IMarlinTrk tools to propagate it to the "
                 "main tracker";

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection(LCIO::TRACK, "InputTrackCollectionName", "Name of the input track collection",
                          _input_track_col_name, std::string("TruthTracks"));

  /////////////////////////

  StringVec vecDigiHitsDefault;
  vecDigiHitsDefault.push_back("ITrackerHits");
  vecDigiHitsDefault.push_back("OTrackerHits");
  vecDigiHitsDefault.push_back("ITrackerEndcapHits");
  vecDigiHitsDefault.push_back("OTrackerEndcapHits");

  registerInputCollections(LCIO::TRACKERHITPLANE, "vecDigiHits",
                           "vector of name of the digi hits collection - need to be syncro with vecSubdetName!",
                           _vecDigiHits, vecDigiHitsDefault);

  StringVec vecSubdetNameDefault;
  vecSubdetNameDefault.push_back("InnerTrackerBarrel");
  vecSubdetNameDefault.push_back("OuterTrackerBarrel");
  vecSubdetNameDefault.push_back("InnerTrackerEndcap");
  vecSubdetNameDefault.push_back("OuterTrackerEndcap");

  registerProcessorParameter("vecSubdetName", "vector of names of all subdetector to exrapolate to", _vecSubdetName,
                             vecSubdetNameDefault);

  /////////////////////////

  registerOutputCollection(LCIO::TRACKERHITPLANE, "OutputNotUsedHitCollectionName",
                           "Name of the output collection with the not used hits", _output_not_used_col_name,
                           std::string("NotUsedHits"));

  registerOutputCollection(LCIO::TRACK, "OutputTrackCollectionName", "Name of the output track collection",
                           _output_track_col_name, std::string("ExtrTracks"));

  registerProcessorParameter("MultipleScatteringOn", "Use MultipleScattering in Fit", _MSOn, bool(true));

  registerProcessorParameter("EnergyLossOn", "Use Energy Loss in Fit", _ElossOn, bool(true));

  registerProcessorParameter("SmoothOn", "Smooth All Mesurement Sites in Fit", _SmoothOn, bool(false));

  registerProcessorParameter("Max_Chi2_Incr", "maximum allowable chi2 increment when moving from one site to another",
                             _Max_Chi2_Incr, double(1000));

  registerProcessorParameter("SearchSigma", "times d0(Z0) acceptable from track extrapolation point", _searchSigma,
                             double(3));

  registerProcessorParameter("PerformFinalRefit", "perform a final refit of the extrapolated track", _performFinalRefit,
                             bool(false));

  registerProcessorParameter("extrapolateForward",
                             "if true extrapolation in the forward direction (in-out), otherwise backward (out-in)",
                             _extrapolateForward, bool(true));
}

void ExtrToTracker::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;

  // usually a good idea to
  printParameters();

  getGeoInfo();

  _trksystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", 0, "");

  ///////////////////////////////

  if (_trksystem == 0) {
    throw EVENT::Exception(std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("DDKalTest"));
  }

  _trksystem->setOption(IMarlinTrkSystem::CFG::useQMS, _MSOn);
  _trksystem->setOption(IMarlinTrkSystem::CFG::usedEdx, _ElossOn);
  _trksystem->setOption(IMarlinTrkSystem::CFG::useSmoothing, _SmoothOn);
  _trksystem->init();

  _n_run = 0;
  _n_evt = 0;
  SITHitsFitted = 0;
  SITHitsNonFitted = 0;
  TotalSITHits = 0;
}

void ExtrToTracker::processRunHeader(LCRunHeader*) { ++_n_run; }

void ExtrToTracker::processEvent(LCEvent* evt) {
  // set the correct configuration for the tracking system for this event
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(_trksystem, _MSOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(_trksystem, _ElossOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);

  // printParameters();

  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG1) << "   processing event: " << _n_evt << std::endl;

  // get input collection and relations
  LCCollection* input_track_col = this->GetCollection(evt, _input_track_col_name);

  if (input_track_col != 0) {
    ////////////////////////

    fillVecSubdet(evt);
    fillMapElHits(_vecDigiHitsCol, _vecMapsElHits);

    ////////////////////////

    // establish the track collection that will be created
    LCCollectionVec* trackVec = new LCCollectionVec(LCIO::TRACK);

    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0);
    trkFlag.setBit(LCIO::TRBIT_HITS);
    trackVec->setFlag(trkFlag.getFlag());

    int nTracks = input_track_col->getNumberOfElements();

    streamlog_out(DEBUG4) << " ######### NO OF TRACKS $$$$$$$$$$ " << nTracks << std::endl;

    LCCollectionVec* inputTrackVec = new LCCollectionVec(LCIO::TRACK);
    inputTrackVec->setSubset(true);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int n = 0; n < nTracks; ++n) {
      Track* testTrack = dynamic_cast<Track*>(input_track_col->getElementAt(n));
      inputTrackVec->addElement(testTrack);
    }

    std::sort(inputTrackVec->begin(), inputTrackVec->end(), PtSort());
    // std::sort( inputTrackVec->begin() , inputTrackVec->end() ,  InversePtSort()  ) ;

    // loop over the input tracks and refit using KalTest
    for (int i = 0; i < nTracks; ++i) {
      int SITHitsPerTrk = 0;

      Track* track = dynamic_cast<Track*>(inputTrackVec->getElementAt(i));

      MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();

      EVENT::TrackerHitVec trkHits = track->getTrackerHits();

      streamlog_out(DEBUG2) << "---- tracks n = " << i << "   n VTX hits = " << trkHits.size() << std::endl;

      // sort the hits in R, so here we are assuming that the track came from the IP

      sort(trkHits.begin(), trkHits.end(), ExtrToTracker::compare_r());

      double z_last_hit = (trkHits.back())->getPosition()[2]; // mm

      EVENT::TrackerHitVec::iterator it = trkHits.begin();

      for (it = trkHits.begin(); it != trkHits.end(); ++it) {
        marlin_trk->addHit(*it);
      }

      int init_status = FitInit2(track, marlin_trk);

      if (init_status == 0) {
        streamlog_out(DEBUG4) << "track initialised " << std::endl;

        int fit_status = marlin_trk->fit();

        if (fit_status == 0) {
          double chi2 = 0;
          int ndf = 0;
          TrackStateImpl trkState0;

          UTIL::BitField64 encoder(lcio::LCTrackerCellID::encoding_string());

          encoder.reset(); // reset to 0

          int layerID = encoder.lowWord();
          int elementID = 0;

          //________________________________________________________________________________________________________
          //
          // starting loop on subdetectors and loop on each subdetector layer
          //________________________________________________________________________________________________________

          for (size_t idet = 0; idet < _vecSubdetName.size(); idet++) {
            streamlog_out(DEBUG4) << "LOOP - detID = " << _vecSubdetID.at(idet) << " begins " << std::endl;

            if (_vecDigiHitsCol.at(idet) != 0) {
              for (int iL = 0; iL < _vecSubdetNLayers.at(idet); iL++) {
                streamlog_out(DEBUG4) << "LOOP - layer = " << iL << " begins " << std::endl;

                encoder[UTIL::LCTrackerCellID::subdet()] = _vecSubdetID.at(idet);

                int side = 0;
                if (_vecSubdetIsBarrel.at(idet))
                  side = 0; // lcio::ILDDetID::barrel;
                else {
                  // for disks: side +1 corresposds to positive z and -1 to negative z
                  // decide to go to either the left or right side according to the z position of the last hit in the
                  // track
                  if (z_last_hit > 0.)
                    side = 1;
                  else
                    side = -1;
                }
                encoder[UTIL::LCTrackerCellID::side()] = side;

                int layer = 0;
                if (_extrapolateForward)
                  layer = iL;
                else
                  layer = _vecSubdetNLayers.at(idet) - (iL + 1);
                encoder[UTIL::LCTrackerCellID::layer()] = layer;

                layerID = encoder.lowWord();
                streamlog_out(DEBUG4) << "layerID = " << layerID << std::endl;

                ///////////////////////////////////////////////////////////

                if (marlin_trk->propagateToLayer(layerID, trkState0, chi2, ndf, elementID, IMarlinTrack::modeClosest) ==
                    MarlinTrk::IMarlinTrack::success) {
                  const FloatVec& covLCIO = trkState0.getCovMatrix();
                  const float* pivot = trkState0.getReferencePoint();
                  double r = sqrt(pivot[0] * pivot[0] + pivot[1] * pivot[1]);

                  streamlog_out(DEBUG4) << " kaltest track parameters: "
                                        << " chi2/ndf " << chi2 / ndf << " chi2 " << chi2 << std::endl

                                        << "\t D0 " << trkState0.getD0() << "[+/-" << sqrt(covLCIO[0]) << "] "
                                        << "\t Phi :" << trkState0.getPhi() << "[+/-" << sqrt(covLCIO[2]) << "] "
                                        << "\t Omega " << trkState0.getOmega() << "[+/-" << sqrt(covLCIO[5]) << "] "
                                        << "\t Z0 " << trkState0.getZ0() << "[+/-" << sqrt(covLCIO[9]) << "] "
                                        << "\t tan(Lambda) " << trkState0.getTanLambda() << "[+/-" << sqrt(covLCIO[14])
                                        << "] "

                                        << "\t pivot : [" << pivot[0] << ", " << pivot[1] << ", " << pivot[2]
                                        << " - r: " << r << "]" << std::endl;

                  streamlog_out(DEBUG4) << " layer " << iL << " max search distances Z : Rphi "
                                        << _searchSigma * sqrt(covLCIO[9]) << " : " << _searchSigma * sqrt(covLCIO[0])
                                        << std::endl;

                  FindAndAddHit(idet, elementID, marlin_trk, trkHits, SITHitsPerTrk, iL);

                } // successful propagation to layer

              } // loop to all subdetector layers

            } // end colDigi not empty

          } // end loop on subdetectors
          streamlog_out(DEBUG4) << " no of hits in the track (after adding SIT hits) " << trkHits.size()
                                << " SIT hits added " << SITHitsPerTrk << " event " << _n_evt << std::endl;

          //==============================================================================================================

          IMPL::TrackImpl* lcio_trk = new IMPL::TrackImpl();

          IMarlinTrack* marlinTrk = 0;

          if (!_performFinalRefit) {
            // fg: ------ here we just create a final LCIO track from the extrapolation :

            marlinTrk = marlin_trk;

            bool fit_direction = IMarlinTrack::forward;
            int return_code = finaliseLCIOTrack(marlin_trk, lcio_trk, trkHits, fit_direction);

            streamlog_out(DEBUG) << " *** created finalized LCIO track - return code " << return_code << std::endl
                                 << *lcio_trk << std::endl;

          } else { // fg: ------- perform a final refit - does not work right now ...

            // refitted track collection creation

            TrackStateImpl* trkState = new TrackStateImpl();
            double chi2_fin = 0.;
            int ndf_fin = 0;

            marlin_trk->getTrackState(*trkState, chi2_fin, ndf_fin);
            // const FloatVec& covMatrix = trkState->getCovMatrix();

            //////////////////////////////////////////////////////////////////////////////////

            sort(trkHits.begin(), trkHits.end(), ExtrToTracker::compare_r());

            bool fit_backwards = IMarlinTrack::backward;
            // bool fit_forwards = IMarlinTrack::forward;

            marlinTrk = _trksystem->createTrack();

            std::vector<EVENT::TrackerHit*> vec_hits;
            vec_hits.clear();
            for (unsigned int ih = 0; ih < trkHits.size(); ih++) {
              vec_hits.push_back(trkHits.at(ih));
            } // end loop on hits
            streamlog_out(DEBUG1) << " --- vec_hits.size() = " << vec_hits.size() << std::endl;

            // int ndf_test_0;
            // int return_error_0 = marlinTrk->getNDF(ndf_test_0);
            // streamlog_out(DEBUG1) << "++++ 0 - getNDF returns " << return_error_0 << std::endl;
            // streamlog_out(DEBUG1) << "++++ 0 - getNDF returns ndf = " << ndf_test_0 << std::endl;

            // Kalman filter smoothing - fit track from out to in
            int error_fit = createFit(vec_hits, marlinTrk, trkState, _bField, fit_backwards, _Max_Chi2_Incr);
            streamlog_out(DEBUG) << "---- createFit - error_fit = " << error_fit << std::endl;

            bool fit_direction = fit_backwards;

            if (error_fit == 0) {
              int error = finaliseLCIOTrack(marlinTrk, lcio_trk, vec_hits, fit_direction);
              streamlog_out(DEBUG) << "---- finalisedLCIOTrack - error = " << error << std::endl;

              // int ndf_test;
              // int return_error = marlinTrk->getNDF(ndf_test);
              // streamlog_out(DEBUG3) << "++++ getNDF returns " << return_error << std::endl;
              // streamlog_out(DEBUG3) << "++++ getNDF returns ndf = " << ndf_test << std::endl;

              if (error != 0) {
                streamlog_out(DEBUG3) << "Error from finaliseLCIOTrack non zero! deleting tracks. error=" << error
                                      << " noHits: " << trkHits.size() << " marlinTrk: " << marlinTrk
                                      << " lcio_trk: " << lcio_trk << std::endl;

                delete lcio_trk;
                delete marlinTrk;
                continue;
              }
            } else {
              streamlog_out(DEBUG3) << "Error from createFit non zero! deleting tracks. error_fit=" << error_fit
                                    << std::endl;

              delete lcio_trk;
              delete marlinTrk;
              continue;
            }

            delete trkState;

          } // !_perfomFinalRefit

          // fit finished - get hits in the fit

          std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
          std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;

          // remember the hits are ordered in the order in which they were fitted

          marlinTrk->getHitsInFit(hits_in_fit);

          if (hits_in_fit.size() < 3) {
            streamlog_out(DEBUG3) << "RefitProcessor: Less than 3 hits in fit: Track Discarded. Number of hits =  "
                                  << trkHits.size() << std::endl;
            delete marlinTrk;
            delete lcio_trk;
            continue;
          }

          marlinTrk->getOutliers(outliers);

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
          MarlinTrk::addHitNumbersToTrack(lcio_trk, all_hits, false, encoder2);
          MarlinTrk::addHitNumbersToTrack(lcio_trk, hits_in_fit, true, encoder2);

          streamlog_out(DEBUG4) << "ExtrToTracker::processEvent - Hit numbers for track " << lcio_trk->id() << ":  "
                                << std::endl;
          int detID = 0;
          for (size_t ip = 0; ip < lcio_trk->subdetectorHitNumbers().size(); ip = ip + 2) {
            detID++;
            streamlog_out(DEBUG4) << "  det id " << detID
                                  << " , nhits in track = " << lcio_trk->subdetectorHitNumbers()[ip]
                                  << " , nhits in fit = " << lcio_trk->subdetectorHitNumbers()[ip + 1] << std::endl;
            if (lcio_trk->subdetectorHitNumbers()[ip] > 0)
              lcio_trk->setTypeBit(detID);
          }

          trackVec->addElement(lcio_trk);

          if (_performFinalRefit)
            delete marlinTrk;

        } // good fit status
      } // good initialisation status

      delete marlin_trk;

    } // for loop to the tracks

    //-------------------------------------------------------------------------------------------------------

    evt->addCollection(trackVec, _output_track_col_name);

    delete inputTrackVec;

    /////////////////////////////////////////////////////////////////
    // Save not used hits in a collection for possible further use //
    /////////////////////////////////////////////////////////////////

    LCCollectionVec* notUsedHitsVec = new LCCollectionVec(LCIO::TRACKERHITPLANE);
    CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder(
        lcio::LCTrackerCellID::encoding_string(),
        notUsedHitsVec); // do not change it, code will not work with a different encoder
    notUsedHitsVec->setSubset(true);

    for (size_t iDet = 0; iDet < _vecMapsElHits.size(); iDet++) {
      for (std::map<int, std::vector<TrackerHitPlane*>>::iterator iterator = _vecMapsElHits.at(iDet).begin();
           iterator != _vecMapsElHits.at(iDet).end(); iterator++) {
        std::vector<TrackerHitPlane*> vecHits = iterator->second;
        for (size_t iHitOnEl = 0; iHitOnEl < vecHits.size(); iHitOnEl++) {
          notUsedHitsVec->addElement(vecHits.at(iHitOnEl));
        } // end loop on hits on each det element
      } // end loop on map detEl <--> vector of hits on the detEl
    } // end loops on vector of maps - one for each subdetector

    evt->addCollection(notUsedHitsVec, _output_not_used_col_name);

    // delete notUsedHitsVec;

    /////////////////////////////////////////////////////////////////

  } // track collection no empty

  ++_n_evt;

  // cout << " event " << _n_evt << std::endl ;
}

void ExtrToTracker::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void ExtrToTracker::end() {
  streamlog_out(DEBUG) << "ExtrToTracker::end()  " << name() << " processed " << _n_evt << " events in " << _n_run
                       << " runs " << std::endl;

  streamlog_out(DEBUG4) << " SIT hits considered for track-hit association " << TotalSITHits
                        << " how many of them were matched and fitted successfully ? " << SITHitsFitted
                        << " for how many the fit failed ? " << SITHitsNonFitted << std::endl;
}

LCCollection* ExtrToTracker::GetCollection(LCEvent* evt, std::string colName) {
  LCCollection* col = NULL;

  try {
    col = evt->getCollection(colName.c_str());
    streamlog_out(DEBUG3) << " --> " << colName.c_str() << " track collection found in event = " << col
                          << " number of elements " << col->getNumberOfElements() << std::endl;
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG3) << " --> " << colName.c_str() << " collection absent in event" << std::endl;
  }

  return col;
}

int ExtrToTracker::FitInit2(Track* track, MarlinTrk::IMarlinTrack* _marlinTrk) {
  // EVENT::FloatVec covMatrix = track->getCovMatrix();

  TrackStateImpl trackState(TrackState::AtOther, track->getD0(), track->getPhi(), track->getOmega(), track->getZ0(),
                            track->getTanLambda(), track->getCovMatrix(), track->getReferencePoint());

  bool direction;
  if (_extrapolateForward)
    direction = IMarlinTrack::forward;
  else
    direction = IMarlinTrack::backward;
  _marlinTrk->initialise(trackState, _bField, direction);

  return IMarlinTrack::success;
}

void ExtrToTracker::fillVecSubdet(lcio::LCEvent*& evt) {
  _vecDigiHitsCol.clear();

  for (size_t idet = 0; idet < _vecSubdetName.size(); idet++) {
    streamlog_out(DEBUG4) << "idet = " << idet << std::endl;

    if (_vecDigiHits.at(idet) != "") {
      LCCollection* colDigi = 0;
      try {
        colDigi = evt->getCollection(_vecDigiHits.at(idet));
        getCellID0AndPositionInfo(colDigi);

      } catch (DataNotAvailableException& e) {
        streamlog_out(DEBUG4) << "Collection " << _vecDigiHits.at(idet).c_str() << " is unavailable " << std::endl;
      }

      _vecDigiHitsCol.push_back(colDigi);

    } // diginame !=0

    else
      _vecDigiHitsCol.push_back(0);

  } // end loop on subdetectors

} // end

void ExtrToTracker::getCellID0AndPositionInfo(LCCollection*& col) {
  std::string cellIDEcoding0 = col->getParameters().getStringVal("CellIDEncoding");
  UTIL::BitField64 cellid_decoder(cellIDEcoding0);

  for (int i = 0; i < col->getNumberOfElements(); i++) {
    TrackerHitPlane* trackerHit = dynamic_cast<TrackerHitPlane*>(col->getElementAt(i));

    const dd4hep::CellID id = trackerHit->getCellID0();
    cellid_decoder.setValue(id);

    int subdet = cellid_decoder[UTIL::LCTrackerCellID::subdet()];
    int side = cellid_decoder[UTIL::LCTrackerCellID::side()];
    int layer = cellid_decoder[UTIL::LCTrackerCellID::layer()];
    int module = cellid_decoder[UTIL::LCTrackerCellID::module()];
    int sensor = cellid_decoder[UTIL::LCTrackerCellID::sensor()];

    streamlog_out(DEBUG2) << " hit" << i << " ( subdetector: " << subdet << ", side: " << side << ", layer: " << layer
                          << ", module: " << module << ", sensor: " << sensor << " ) ( r: "
                          << sqrt(pow(trackerHit->getPosition()[0], 2) + pow(trackerHit->getPosition()[1], 2))
                          << " , phi: " << atan(trackerHit->getPosition()[1] / trackerHit->getPosition()[0])
                          << " , z: " << trackerHit->getPosition()[2] << " ) \n";

    // SurfaceMap::const_iterator si = _surfMap.find(id);

    // if( si == _surfMap.end()) {
    //   streamlog_out(MESSAGE2) << "----- IT IS A SURFACES PROBLEM "  <<std::endl;
    // } else {
    //   streamlog_out(MESSAGE2) << "----- OK I HAVE SOMETHING "  <<std::endl;
    // }

    // ISurface* surf = (si != _surfMap.end() ?  si->second  : 0);

    // streamlog_out(MESSAGE2) << " surf = " << surf <<std::endl;

    // const double* hit_pos = trackerHit->getPosition();
    // DDSurfaces::Vector3D hit_global(hit_pos[0],hit_pos[1],hit_pos[2]);
    // DDSurfaces::Vector2D hit_local = surf->globalToLocal( dd4hep::mm * hit_global );

    // streamlog_out(MESSAGE2) << "----- hit global x, y, r = " << hit_global[0] <<"   "<< hit_global[1] <<"   "<<
    // sqrt(pow(hit_global[0],2)+pow(hit_global[1],2)) <<std::endl; streamlog_out(MESSAGE2) << "----- hit local U, V = "
    // << hit_local[0] <<"   "<< hit_local[1] <<std::endl;
  }

  return;
}

TrackerHitPlane* ExtrToTracker::getSiHit(std::vector<TrackerHitPlane*>& hitsOnDetEl,
                                         MarlinTrk::IMarlinTrack*& marlin_trk) {
  double min = 9999999.;
  double testChi2 = 0.;
  size_t nHitsOnDetEl = hitsOnDetEl.size();
  int index = -1; // index of the selected hits

  streamlog_out(DEBUG2) << "-- number of hits on the same detector element: " << nHitsOnDetEl << std::endl;

  if (nHitsOnDetEl == 0)
    return 0;

  for (size_t i = 0; i < nHitsOnDetEl; i++) {
    // if ( marlin_trk->testChi2Increment(hitsOnDetEl.at(i), testChi2) == MarlinTrk::IMarlinTrack::success ) {..} // do
    // not do this the testChi2 internally call the addandfir setting as max acceptable chi2 a very small number to be
    // sure to never add the the hit to the fit (so it always fails)
    marlin_trk->testChi2Increment(hitsOnDetEl.at(i), testChi2);
    streamlog_out(DEBUG2) << "-- testChi2: " << testChi2 << std::endl;
    if (min > testChi2 && testChi2 > 0) {
      min = testChi2;
      index = i;
    } // end min chi2
  } // end loop on hits on same detector element

  streamlog_out(DEBUG2) << "-- index of the selected hit: " << index << std::endl;

  if (index == -1)
    return 0;
  else {
    TrackerHitPlane* selectedHit = dynamic_cast<TrackerHitPlane*>(hitsOnDetEl.at(index));

    if (nHitsOnDetEl > 1)
      std::iter_swap(hitsOnDetEl.begin() + index, hitsOnDetEl.begin() + nHitsOnDetEl - 1);
    hitsOnDetEl.pop_back();
    // hitsOnDetEl.erase(hitsOnDetEl.begin()+index);

    return selectedHit;
  }

} // end getSiHit

TrackerHitPlane* ExtrToTracker::getSiHit(std::vector<dd4hep::CellID>& vecElID,
                                         std::map<int, std::vector<TrackerHitPlane*>>& mapElHits,
                                         MarlinTrk::IMarlinTrack*& marlin_trk) {
  double min = 9999999.;
  double testChi2 = 0.;
  size_t nElID = vecElID.size();
  int indexel = -1; // element index (in vecElID) of the selected hits
  int index = -1;   // index of the selected hits

  for (size_t ie = 0; ie < nElID; ie++) {
    int elID = vecElID.at(ie);
    size_t nHitsOnDetEl = 0;
    if (mapElHits.count(elID) > 0) {
      nHitsOnDetEl = mapElHits[elID].size();
    }

    // streamlog_out(MESSAGE2) << "-- elID at index = "<< ie <<" / " << nElID << " : " << vecElID.at(ie) << std::endl ;
    streamlog_out(DEBUG3) << "-- number of hits on the same detector element: " << nHitsOnDetEl << std::endl;

    for (size_t i = 0; i < nHitsOnDetEl; i++) {
      marlin_trk->testChi2Increment(mapElHits[elID].at(i), testChi2);
      // streamlog_out(MESSAGE2) << "-- trackerhit at index = " << i << " / "<< nHitsOnDetEl << std::endl ;
      // streamlog_out(MESSAGE2) << "-- testChi2: " << testChi2 << std::endl ;
      if (min > testChi2 && testChi2 > 0) {
        min = testChi2;
        indexel = elID;
        index = i;
      }
    } // end loop on hits on the same elID

  } // end loop on elIDs

  if (index == -1 || indexel == -1)
    return 0;
  else {
    streamlog_out(DEBUG3) << "-- hit added " << std::endl;
    // if ( vecElID.at(0) != indexel) streamlog_out(MESSAGE2) << "-- but not from the first elementID " << std::endl ;

    TrackerHitPlane* selectedHit = dynamic_cast<TrackerHitPlane*>(mapElHits[indexel].at(index));

    int nHitsOnSelectedEl = mapElHits[indexel].size();
    if (nHitsOnSelectedEl > 1)
      std::iter_swap(mapElHits[indexel].begin() + index, mapElHits[indexel].begin() + nHitsOnSelectedEl - 1);
    mapElHits[indexel].pop_back();
    ////(mapElHits[indexel].erase((mapElHits[indexel].begin()+index);

    return selectedHit;
  }

  return 0;

} // end getSiHit

void ExtrToTracker::fillMapElHits(std::vector<LCCollection*>& vecHitCol,
                                  std::vector<std::map<int, std::vector<TrackerHitPlane*>>>& vecMaps) {
  // fill map (el - vector of hits) for each subdtector

  vecMaps.clear();
  _vecvecHitsInCol.clear();

  for (size_t icol = 0; icol < vecHitCol.size(); icol++) {
    std::map<int, std::vector<TrackerHitPlane*>> map_el_hits;

    std::vector<TrackerHitPlane*> vecHelper;
    vecHelper.clear();

    if (vecHitCol.at(icol) != NULL) {
      int nhits = vecHitCol.at(icol)->getNumberOfElements();
      streamlog_out(DEBUG2) << "  nhits = " << nhits << std::endl;

      for (int ihit = 0; ihit < nhits; ihit++) {
        TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(vecHitCol.at(icol)->getElementAt(ihit));

        vecHelper.push_back(hit);

        int cellID0 = hit->getCellID0();
        UTIL::BitField64 encoder0(
            lcio::LCTrackerCellID::encoding_string()); // do not change it, code will not work with a different encoder
        encoder0.reset();                              // reset to 0
        encoder0.setValue(cellID0);
        int hitElID = encoder0.lowWord();

        map_el_hits[hitElID].push_back(hit);

      } // end loop on hits

    } // collection not empty

    vecMaps.push_back(map_el_hits);

    _vecvecHitsInCol.push_back(vecHelper);

  } // end loop on collections of hits

} // end fillMapElHits

void ExtrToTracker::getGeoInfo() {
  _vecSubdetID.clear();
  _vecSubdetNLayers.clear();
  _vecSubdetIsBarrel.clear();
  _vecMapNeighbours.clear();

  dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();

  // alternative way
  //  const std::vector< dd4hep::DetElement>& barrelDets = dd4hep::DetectorSelector(theDetector).detectors(  (
  //  DD4hep::DetType::TRACKER | DD4hep::DetType::BARREL )) ;

  // streamlog_out( MESSAGE2 ) << " --- flag = " << (DD4hep::DetType::TRACKER | DD4hep::DetType::BARREL) <<std::endl;
  // streamlog_out( MESSAGE2 ) << " --- number of barrel detectors = " << barrelDets.size() <<std::endl;

  const double pos[3] = {0, 0, 0};
  double bFieldVec[3] = {0, 0, 0};
  mainDetector.field().magneticField(pos, bFieldVec); // get the magnetic field vector from DD4hep
  _bField = bFieldVec[2] / dd4hep::tesla;             // z component at (0,0,0)

  streamlog_out(DEBUG2) << " - _bField = " << _bField << std::endl;

  for (size_t i = 0; i < _vecSubdetName.size(); i++) {
    int detID = 0;
    int nlayers = 0;
    bool isBarrel = false;

    try {
      const dd4hep::DetElement& theDetector = mainDetector.detector(_vecSubdetName.at(i));
      auto detType = dd4hep::DetType(theDetector.typeFlag());
      isBarrel = detType.is(dd4hep::DetType::BARREL);

      // streamlog_out( DEBUG1 ) << " is an endcap " << std::boolalpha << detType.is( DD4hep::DetType::ENDCAP ) <<
      // std::endl;
      streamlog_out(DEBUG1) << " is a barrel " << std::boolalpha << detType.is(dd4hep::DetType::BARREL) << std::endl;

      detID = theDetector.id();
      streamlog_out(DEBUG2) << " --- subdet: " << _vecSubdetName.at(i) << " - id = " << detID << std::endl;

      if (isBarrel) {
        dd4hep::rec::ZPlanarData* theExtension = 0;
        theExtension = theDetector.extension<dd4hep::rec::ZPlanarData>();

        nlayers = theExtension->layers.size();

        streamlog_out(DEBUG2) << " - n layers = " << nlayers << std::endl;

      } // end barrel type
      else {
        dd4hep::rec::ZDiskPetalsData* theExtension = 0;
        theExtension = theDetector.extension<dd4hep::rec::ZDiskPetalsData>();

        nlayers = theExtension->layers.size();

        streamlog_out(DEBUG2) << " - n layers = " << nlayers << std::endl;

      } // end endcap type

      dd4hep::rec::NeighbourSurfacesData* neighbourSurfaces = 0;
      neighbourSurfaces = theDetector.extension<dd4hep::rec::NeighbourSurfacesData>();
      _vecMapNeighbours.push_back(&(neighbourSurfaces->sameLayer));

    } catch (std::runtime_error& exception) {
      streamlog_out(WARNING)
          << "ExtrToTracker::getGeoInfo - exception in retriving number of modules per layer for subdetector : "
          << _vecSubdetName.at(i) << " : " << exception.what() << std::endl;
    }

    _vecSubdetID.push_back(detID);
    _vecSubdetNLayers.push_back(nlayers);
    _vecSubdetIsBarrel.push_back(isBarrel);

  } // end loop on subdetector names

} // end getGeoInfo

void ExtrToTracker::FindAndAddHit(size_t& idet, int& elID, MarlinTrk::IMarlinTrack*& mtrk,
                                  EVENT::TrackerHitVec& trkHits, int& SITHitsPerTrk, int& layer) {
  //_______________________________________________________________________________________
  //

  streamlog_out(DEBUG2) << " element ID " << elID << std::endl;

  if (elID != 0) {
    bool isSuccessfulFit = false;

    std::vector<dd4hep::CellID> vecIDs;
    vecIDs = _vecMapNeighbours.at(idet)->find(elID)->second;

    vecIDs.insert(std::begin(vecIDs), elID);

    TrackerHitPlane* BestHit;
    BestHit = getSiHit(vecIDs, _vecMapsElHits.at(idet), mtrk);

    // BestHit = getSiHit(_vecvecHitsInCol.at(idet), mtrk);

    if (BestHit != 0) {
      streamlog_out(DEBUG2) << " --- Best hit found " << std::endl;

      double chi2_increment = 0.;
      // isSuccessfulFit = mtrk->addAndFit( BestHit, chi2_increment, _Max_Chi2_Incr*(layer+1) ) == IMarlinTrack::success
      // ;
      isSuccessfulFit = mtrk->addAndFit(BestHit, chi2_increment, _Max_Chi2_Incr) == IMarlinTrack::success;

      streamlog_out(DEBUG4) << " --- layer+1 = " << layer + 1 << std::endl;
      streamlog_out(DEBUG4) << " --- _Max_Chi2_Incr = " << _Max_Chi2_Incr << std::endl;
      streamlog_out(DEBUG4) << " --- increment in the chi2 = " << chi2_increment << "  , max chi2 to accept the hit "
                            << _Max_Chi2_Incr * (layer + 1) << std::endl;
      streamlog_out(DEBUG4) << " --- isSuccessfulFit = " << isSuccessfulFit << std::endl;

      // TotalSITHits++;

      if (isSuccessfulFit) {
        trkHits.push_back(BestHit);

        streamlog_out(DEBUG4) << " +++ hit added " << BestHit << std::endl;

        SITHitsPerTrk++;
        // SITHitsFitted++;

        streamlog_out(DEBUG4) << " --- SITHitsPerTrk " << SITHitsPerTrk << std::endl;

      } // end successful fit
      else {
        // SITHitsNonFitted++;
        streamlog_out(DEBUG4) << " +++ HIT NOT ADDED " << std::endl;
      }

    } // besthit found

    // streamlog_out(DEBUG4) << " -- SITHitsFitted = " << SITHitsFitted << std::endl;
    // streamlog_out(DEBUG4) << " -- SITHitsNonFitted = " << SITHitsNonFitted << std::endl;

  } // elID !=0

  return;
}
