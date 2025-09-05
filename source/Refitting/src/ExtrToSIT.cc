#include "ExtrToSIT.h"

#include <IMPL/TrackerHitPlaneImpl.h>
#include <UTIL/ILDConf.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"

using namespace lcio;
using namespace marlin;
using namespace MarlinTrk;

//----------------------------------------------------------------
struct Distance3D2 {
  Vector3D _pos;
  Distance3D2(const Vector3D& pos) : _pos(pos) {}
  template <class T>
  double operator()(const T* t) {
    Vector3D p(t->getPosition());
    return (p - _pos).r2();
  }
};
//----------------------------------------------------------------

//------------------------------------------------------------------------------------------

struct PtSort { // sort tracks wtr to pt - largest first
  inline bool operator()(const lcio::LCObject* l, const lcio::LCObject* r) {
    return (std::abs(((const lcio::Track*)l)->getOmega()) <
            std::abs(((const lcio::Track*)r)->getOmega())); // pt ~ 1./omega
  }
};

//------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------

struct ZSort { // sort TPC track segments wtr to Z - smallest first
  inline bool operator()(const lcio::LCObject* l, const lcio::LCObject* r) {
    return (std::abs(((const lcio::Track*)l)->getZ0()) < std::abs(((const lcio::Track*)r)->getZ0())); // pt ~ 1./omega
  }
};

//------------------------------------------------------------------------------------------

ExtrToSIT aExtrToSIT;

ExtrToSIT::ExtrToSIT() : Processor("ExtrToSIT") {
  // modify processor description
  _description =
      "ExtrToSIT refits an input track collection (TPC or VXD), and used IMarlinTrk tools to propagate it to SIT";

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection(LCIO::TRACK, "InputTrackCollectionName", "Name of the input track collection",
                          _input_track_col_name, std::string("TruthTracks"));

  registerInputCollection(LCIO::LCRELATION, "InputTrackRelCollection",
                          "Name of the MCParticle-Track Relations collection for input tracks", _input_track_rel_name,
                          std::string("TruthTracksMCP"));

  registerInputCollection(LCIO::TRACKERHITPLANE, "digitisedVXDHits", "Name of the VTXTrackerHit collection",
                          _vxdColName, std::string("VTXTrackerHits"));

  registerInputCollection(LCIO::TRACKERHITPLANE, "digitisedSITHits", "Name of the SITTrackerHit collection",
                          _sitColName, std::string("SITTrackerHits"));

  registerOutputCollection(LCIO::TRACK, "OutputTrackCollectionName", "Name of the output track collection",
                           _output_track_col_name, std::string("RefittedTracks"));

  registerOutputCollection(LCIO::TRACK, "SiliconCollectionName", "Name of the output silicon track collection",
                           _siTrkColName, std::string("SiliconTracks"));

  registerOutputCollection(LCIO::LCRELATION, "OutputTrackRelCollection",
                           "Name of the MCParticle-Track Relations collection for output tracks",
                           _output_track_rel_name, std::string("RefittedTracksMCP"));

  registerProcessorParameter("MultipleScatteringOn", "Use MultipleScattering in Fit", _MSOn, bool(true));

  registerProcessorParameter("EnergyLossOn", "Use Energy Loss in Fit", _ElossOn, bool(true));

  registerProcessorParameter("SmoothOn", "Smooth All Mesurement Sites in Fit", _SmoothOn, bool(false));

  registerProcessorParameter("TrackSystemName",
                             "Name of the track fitting system to be used (KalTest, DDKalTest, aidaTT, ... )",
                             _trkSystemName, std::string("DDKalTest"));

  registerProcessorParameter(
      "DirInsideOut", "direction for the extrapolation. if true it means we extrapolate from VXD, otherwise from TPC",
      _dirInsideOut, bool(true));

  registerProcessorParameter("Max_Chi2_Incr", "maximum allowable chi2 increment when moving from one site to another",
                             _Max_Chi2_Incr, double(1000));

  registerProcessorParameter("PropagateToLayer", "Which layer should the seed be propagated to", _propToLayer, int(4));

  registerProcessorParameter("TPCHitsCut", "minimum acceptable no of hits for the TPC tracks", _tpcHitsCut, int(6));

  registerProcessorParameter("Chi2NDoFCut", "maximum acceptable chi2/ndof for the TPC tracks", _chi2NDoFCut,
                             float(100));

  registerProcessorParameter("DoCut", "maximum acceptable D0 at IP", _DoCut, float(2));

  registerProcessorParameter("ZoCut", "maximum acceptable Z0 at IP", _ZoCut, float(5));

  registerProcessorParameter("SearchSigma", "times d0(Z0) acceptable from track extrapolation point", _searchSigma,
                             double(3));

  registerProcessorParameter("IsSpacePoints", "If we use space points rather than hits (SIT)", _isSpacePoints,
                             bool(false));

  registerProcessorParameter(
      "NHitsChi2", "Maximal number of hits for which a track with n hits is better than one with n-1hits. (defaut 5)",
      _nHitsChi2, int(5));

  StringVec trackerHitsRelInputColNamesDefault;
  trackerHitsRelInputColNamesDefault.push_back("VXDTrackerHitRelations");
  trackerHitsRelInputColNamesDefault.push_back("SITTrackerHitRelations");

  registerInputCollections("LCRelation", "TrackerHitsRelInputCollections",
                           "Name of the lcrelation collections, that link the TrackerHits to their SimTrackerHits.",
                           _colNamesTrackerHitRelations, trackerHitsRelInputColNamesDefault);

  registerInputCollection(LCIO::MCPARTICLE, "MCParticleCollection", "Name of the MCParticle input collection",
                          _mcParticleCollectionName, std::string("MCParticle"));
}

void ExtrToSIT::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;

  // usually a good idea to
  printParameters();

  // set up the geometery needed by DDKalTest
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
  SITHitsFitted = 0;
  SITHitsNonFitted = 0;
  TotalSITHits = 0;

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  double bFieldVec[3];
  theDetector.field().magneticField({0, 0, 0}, bFieldVec); // get the magnetic field vector from DD4hep
  _bField = bFieldVec[2] / dd4hep::tesla;                  // z component at (0,0,0)

  dd4hep::DetElement sitDE = theDetector.detector("SIT");
  dd4hep::rec::ZPlanarData* sit = sitDE.extension<dd4hep::rec::ZPlanarData>();
  _nSITLayers = sit->layers.size();

  _maxChi2PerHit = 100;
}

void ExtrToSIT::processRunHeader(LCRunHeader*) { ++_n_run; }

void ExtrToSIT::processEvent(LCEvent* evt) {
  // set the correct configuration for the tracking system for this event
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(_trksystem, _MSOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(_trksystem, _ElossOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);

  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG1) << "   processing event: " << _n_evt << std::endl;

  // get input collection and relations
  LCCollection* input_track_col = this->GetCollection(evt, _input_track_col_name);
  LCCollection* sitHitsCol = this->GetCollection(evt, _sitColName);

  if (input_track_col != 0) {
    // establish the track collection that will be created
    LCCollectionVec* trackVec = new LCCollectionVec(LCIO::TRACK);
    // LCCollectionVec* SiTrkCol  = new LCCollectionVec(LCIO::TRACK);

    // relation collections

    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0);
    trkFlag.setBit(LCIO::TRBIT_HITS);
    trackVec->setFlag(trkFlag.getFlag());

    int nTracks = input_track_col->getNumberOfElements();

    streamlog_out(DEBUG4) << " ######### NO OF TRACKS $$$$$$$$$$ " << nTracks << std::endl;

    LCCollectionVec* inputTrackVec = new LCCollectionVec(LCIO::TRACK);

    for (int n = 0; n < nTracks; ++n) {
      Track* testTrack = dynamic_cast<Track*>(input_track_col->getElementAt(n));
      inputTrackVec->addElement(testTrack);
    }

    std::sort(inputTrackVec->begin(), inputTrackVec->end(), PtSort());

    int sitHits = 0;
    if (sitHitsCol != 0) {
      sitHits = sitHitsCol->getNumberOfElements();
    }

    EVENT::TrackerHitVec usedSiHits;
    std::vector<IMPL::TrackImpl*> trackCandidates;

    // loop over the input tracks and refit using KalTest
    for (int iTrack = 0; iTrack < nTracks; ++iTrack) {
      int SITHitsPerTrk = 0;

      // Track* track = dynamic_cast<Track*>( input_track_col->getElementAt( i ) ) ;
      Track* track = static_cast<Track*>(inputTrackVec->getElementAt(iTrack));

      MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();

      EVENT::TrackerHitVec trkHits = track->getTrackerHits();

      // apply a chi2/ndf track sample selection
      float chisquare = track->getChi2();
      int ndof = track->getNdf();

      // d0, Z0 values of TPC tracks

      float chi2_ndf = chisquare / (1.0 * ndof);

      // accept only tracks with > minimum hits at the TPC & chi2/ndf < maximum cut & Pt cut?
      if (chi2_ndf < _chi2NDoFCut) {
        streamlog_out(DEBUG4) << "%%%% Checking track " << track->id() << " chi2/ndf " << chi2_ndf << " cut "
                              << _chi2NDoFCut << " D0 " << track->getD0() << " Z0 " << track->getZ0() << std::endl;

        // sort the hits in R, so here we are assuming that the track came from the IP and that we want to fit out to
        // in.
        sort(trkHits.begin(), trkHits.end(), ExtrToSIT::compare_r());

        for (EVENT::TrackerHitVec::iterator it = trkHits.begin(); it != trkHits.end(); ++it) {
          marlin_trk->addHit(*it);
        }

        // int init_status = marlin_trk->initialise( IMarlinTrack::backward ) ;
        // int init_status = FitInit(trkHits, marlin_trk) ;                   // alternative way to initialise the
        // marlin_trk object
        int init_status = FitInit2(track, marlin_trk);

        if (init_status == 0) {
          streamlog_out(DEBUG4) << "track initialised " << std::endl;

          int fit_status = marlin_trk->fit();

          if (fit_status == 0) {
            int testFlag = 0;

            // some testing
            TrackStateImpl testTS;
            double chi2_test = 0;
            int ndf_test = 0;
            int testTrackState = marlin_trk->getTrackState(testTS, chi2_test, ndf_test);
            if (testTrackState == 0) {
              const float* testpivot = testTS.getReferencePoint();
              streamlog_out(DEBUG4) << "test pivot    " << testpivot[0] << ", " << testpivot[1] << ", " << testpivot[2]
                                    << " - r: " << sqrt(testpivot[0] * testpivot[0] + testpivot[1] * testpivot[1])
                                    << std::endl;
            }

            streamlog_out(DEBUG3) << "###########$$$$$$$$$$##############" << std::endl;

            // now we have a track constisting only of TPC hits

            // 1) extrapolate to outer most SIT layer

            // 2) select best hit candidate

            // 3) marlin_trk->addAndFit(hit,max_chi2_increment);

            // 3a) select best hit candidate on n+1 layer

            // 3b) marlin_trk->addAndFit(hit,max_chi2_increment);

            // 4) loop and repeat for other layers, plus VXD

            // 5) all hits have been added, propogate to IP for instance

            //_________________________________________________________

            // marlin_trk->smooth(trkHits.back());

            double chi2 = 0;
            int ndf = 0;

            Vector3D xing_point;

            UTIL::BitField64 encoder(lcio::LCTrackerCellID::encoding_string());

            encoder.reset(); // reset to 0

            int layerID = encoder.lowWord();
            int elementID = 0;

            //________________________________________________________________________________________________________
            //
            // starting loop to SIT layers
            //________________________________________________________________________________________________________

            // for loop to all SIT layers

            for (unsigned int iL = 0; iL < _nSITLayers; iL++) {
              if (sitHitsCol != 0) {
                streamlog_out(DEBUG4) << " Do I come into the loop " << std::endl;

                streamlog_out(DEBUG4) << "LOOP" << iL << " begins " << std::endl;

                encoder[lcio::LCTrackerCellID::subdet()] = ILDDetID::SIT;

                if (_dirInsideOut)
                  encoder[lcio::LCTrackerCellID::layer()] = iL; // in case we propagate outwards from VXD
                else
                  encoder[lcio::LCTrackerCellID::layer()] = 3 - iL; //  in case we propagate inwards from TPC

                layerID = encoder.lowWord();

                TrackStateImpl trkState;

                if (marlin_trk->propagateToLayer(layerID, trkState, chi2, ndf, elementID, IMarlinTrack::modeClosest) ==
                    MarlinTrk::IMarlinTrack::success) {
                  const FloatVec& covLCIO = trkState.getCovMatrix();
                  const float* pivot = trkState.getReferencePoint();

                  streamlog_out(DEBUG4) << " kaltest track parameters: "
                                        << " chi2/ndf " << chi2 / ndf << " chi2 " << chi2 << std::endl

                                        << "\t D0 " << trkState.getD0() << "[+/-" << sqrt(covLCIO[0]) << "] "
                                        << "\t Phi :" << trkState.getPhi() << "[+/-" << sqrt(covLCIO[2]) << "] "
                                        << "\t Omega " << trkState.getOmega() << "[+/-" << sqrt(covLCIO[5]) << "] "
                                        << "\t Z0 " << trkState.getZ0() << "[+/-" << sqrt(covLCIO[9]) << "] "
                                        << "\t tan(Lambda) " << trkState.getTanLambda() << "[+/-" << sqrt(covLCIO[14])
                                        << "] "

                                        << "\t pivot : [" << pivot[0] << ", " << pivot[1] << ", " << pivot[2]
                                        << " - r: " << sqrt(pivot[0] * pivot[0] + pivot[1] * pivot[1]) << "]"
                                        << std::endl;

                  double r = sqrt(pivot[0] * pivot[0] + pivot[1] * pivot[1]);

                  streamlog_out(DEBUG4) << " layer " << 3 - iL << " max search distances Z : Rphi "
                                        << _searchSigma * sqrt(covLCIO[9]) << " : " << _searchSigma * sqrt(covLCIO[0])
                                        << std::endl;

                  //_______________________________________________________________________________________
                  //
                  // track - hit association (just for fun)

                  streamlog_out(DEBUG2) << " element ID " << elementID << std::endl;

                  if (elementID != 0) {
                    testFlag = 1;

                    EVENT::TrackerHitVec HitsInLayer;

                    float dU_spres = 0;
                    float dV_spres = 0;

                    double chi2_increment = 0;

                    for (int iSit = 0; iSit < sitHits; iSit++) {
                      TrackerHit* hit = dynamic_cast<TrackerHit*>(sitHitsCol->getElementAt(iSit));

                      streamlog_out(DEBUG2) << " type = " << hit->getType() << std::endl;

                      const int celId = hit->getCellID0();

                      streamlog_out(DEBUG1) << " hit cellid0 = " << celId << std::endl;

                      int layerNumber = 0;
                      int ladderNumber = 0;
                      int sideTest = 0;
                      int sensorNumber = 0;

                      // UTIL::BitField64 encoder( lcio::LCTrackerCellID::encoding_string() ) ;
                      encoder.setValue(celId);
                      layerNumber = encoder[lcio::LCTrackerCellID::layer()];
                      ladderNumber = encoder[lcio::LCTrackerCellID::module()];
                      sideTest = encoder[lcio::LCTrackerCellID::side()];
                      sensorNumber = encoder[lcio::LCTrackerCellID::sensor()];
                      encoder.reset();

                      // Just to check the element matching between the hit (coming from the digitiser) and the track
                      // extrapolation element (coming from Mokka)

                      // int mokkaLayerNumber = 0 ;
                      // int mokkaLadderNumber = 0 ;
                      // int mokkaSideTest = 0 ;
                      // int mokkaSensorNumber = 0 ;

                      encoder.setValue(elementID);
                      // mokkaLayerNumber  = encoder[lcio::LCTrackerCellID::layer()] ;
                      // mokkaLadderNumber = encoder[lcio::LCTrackerCellID::module()] ;
                      // mokkaSideTest = encoder[lcio::LCTrackerCellID::side()] ;
                      // mokkaSensorNumber = encoder[lcio::LCTrackerCellID::sensor()] ;
                      encoder.reset();

                      streamlog_out(DEBUG2)
                          << " checking hit : type = " << hit->getType() << " cell ID = " << celId
                          << " side = " << sideTest << " layer = " << layerNumber << " ladder = " << ladderNumber
                          << " sensor = " << sensorNumber << std::endl;

                      streamlog_out(DEBUG2)
                          << " the element id where the hit is found " << celId
                          << " the element id where the track is extrapolated " << elementID << std::endl;

                      if (celId == elementID) { // cause elementID does not give sensor-side infos...
                        // if (mokkaLayerNumber==layerNumber && mokkaLadderNumber==ladderNumber){
                        streamlog_out(DEBUG2)
                            << " We found a hit at the right element with type : " << hit->getType()
                            << " side = " << sideTest << " layer = " << layerNumber << " ladder = " << ladderNumber
                            << " sensor = " << sensorNumber << std::endl;
                        streamlog_out(DEBUG3)
                            << " We found a hit at the right element with type : " << hit->getType() << std::endl;
                        HitsInLayer.push_back(hit);
                      }
                    } // end loop on hits
                    /*
                    // ------------------ in order to obtain sensors single point resolution ---------------------
                    if (_pixelSIT){
                      for (int jj=0;jj<sitHits;jj++){

                        TrackerHitPlane* hitplane = dynamic_cast<TrackerHitPlane*>( sitHitsCol->getElementAt( jj ) ) ;
                        if (hitplane){
                          const int celId_spres = hitplane->getCellID0() ;

                          if (celId_spres==elementID){
                            dU_spres = hitplane->getdU();
                            dV_spres = hitplane->getdV();
                            //streamlog_out(DEBUG4) << " U resolution = " << dU_spres << " V resolution = " << dV_spres
                    << std::endl ;
                          }
                        }
                      }
                    }
                    */
                    // if (!_pixelSIT){
                    //  Alternative way to assign the sensor s.p. resolution until I understand what's going on with
                    //  that
                    dU_spres = 0.007;
                    dV_spres = 0.05;
                    //}
                    // ---------------------------------------------------------------------------------------------

                    if (HitsInLayer.size() != 0) {
                      TrackerHit* BestHit;
                      bool BestHitFound = false;

                      int pointer = 0;
                      int PossibleHits = 0;
                      double DimDist = 0;

                      streamlog_out(DEBUG4)
                          << " calling selectbestcandidatelimited with value for possible hits = " << PossibleHits
                          << std::endl;

                      SelectBestCandidateLimited(HitsInLayer, pivot, BestHit, covLCIO, r, BestHitFound, _searchSigma,
                                                 pointer, PossibleHits, dU_spres, dV_spres, DimDist, usedSiHits);

                      if (BestHitFound) { // when selection is restricted to an area maybe we will not find an
                                          // appropriate hit

                        bool isSuccessful2 = false;

                        streamlog_out(DEBUG4) << " Best hit found: call add and fit " << std::endl;

                        TrackerHitPlaneImpl* TestHitPlane = new TrackerHitPlaneImpl;

                        TestHitPlane->setCellID0(BestHit->getCellID0());
                        TestHitPlane->setPosition(BestHit->getPosition());
                        TestHitPlane->setdU(dU_spres);
                        TestHitPlane->setdV(dV_spres);

                        isSuccessful2 = marlin_trk->addAndFit(TestHitPlane, chi2_increment, _Max_Chi2_Incr) ==
                                        IMarlinTrack::success;
                        // isSuccessful2 = marlin_trk->addAndFit( BestHit, chi2_increment, _Max_Chi2_Incr ) ==
                        // IMarlinTrack::success;

                        TotalSITHits++;

                        if (isSuccessful2) {
                          streamlog_out(DEBUG4) << " successful fit " << std::endl;

                          trkHits.push_back(BestHit);

                          // usedSiHits.push_back(BestHit) ;

                          streamlog_out(DEBUG4) << " hit added " << BestHit << std::endl;

                          SITHitsPerTrk++;
                          SITHitsFitted++;

                          HitsInLayer.erase(HitsInLayer.begin() + pointer);

                        }

                        else {
                          SITHitsNonFitted++;
                        }

                        delete TestHitPlane;

                      } // condition for found hit ends
                      HitsInLayer.clear();
                    } // end of loop on layer i

                    streamlog_out(DEBUG4) << "LOOP" << iL << " ends " << std::endl;
                    streamlog_out(DEBUG4) << "###########$$$$$$$$$$##############" << std::endl;
                  }
                  // track - hit association fini

                } // successful propagation

              } // loop to all SIT layers

            } // condition for SIT digitised hits collection

            streamlog_out(DEBUG4) << " no of hits in the track (after adding SIT hits) " << trkHits.size()
                                  << " SIT hits added " << SITHitsPerTrk << " event " << _n_evt << std::endl;

            // refitted track collection creation

            if (testFlag == 1) {
              TrackStateImpl* trkState = new TrackStateImpl();
              double chi2_fin = 0;
              int ndf_fin = 0;

              marlin_trk->getTrackState(*trkState, chi2_fin, ndf_fin);

              const FloatVec& covMatrix = trkState->getCovMatrix();

              IMPL::TrackImpl* refittedTrack = new IMPL::TrackImpl();

              //-----------------------------------------------------------------------------------------------------------------------
              // final refit of the track
              std::vector<std::pair<float, EVENT::TrackerHit*>> r2_values;
              r2_values.reserve(trkHits.size());

              streamlog_out(DEBUG4) << " size of hits vector thats gonna be added to the output track "
                                    << trkHits.size() << std::endl;

              for (TrackerHitVec::iterator it = trkHits.begin(); it != trkHits.end(); ++it) {
                EVENT::TrackerHit* h = *it;
                float r2 = h->getPosition()[0] * h->getPosition()[0] + h->getPosition()[1] * h->getPosition()[1];
                r2_values.push_back(std::make_pair(r2, *it));
              }

              sort(r2_values.begin(), r2_values.end());

              trkHits.clear();
              trkHits.reserve(r2_values.size());

              for (std::vector<std::pair<float, EVENT::TrackerHit*>>::iterator it = r2_values.begin();
                   it != r2_values.end(); ++it) {
                trkHits.push_back(it->second);
              }

              // bool fit_backwards = IMarlinTrack::backward;
              bool fit_forwards = IMarlinTrack::forward;
              MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();
              // int error = 0;

              try {
                // int error =
                MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, refittedTrack, fit_forwards, covMatrix, _bField,
                                                    _maxChi2PerHit);

              } catch (...) {
                delete refittedTrack;
                delete marlinTrk;

                throw;
              }

              // fitting finished get hit in the fit

              std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
              std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;

              // remember the hits are ordered in the order in which they were fitted
              // here we are fitting inwards to the first is the last and vice verse

              marlinTrk->getHitsInFit(hits_in_fit);

              if (hits_in_fit.size() < 3) {
                streamlog_out(DEBUG3) << "RefitProcessor: Less than 3 hits in fit: Track Discarded. Number of hits =  "
                                      << trkHits.size() << std::endl;
                delete marlinTrk;
                delete refittedTrack;
                continue;
              }

              std::vector<TrackerHit*> all_hits;
              all_hits.reserve(300);

              for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
                all_hits.push_back(hits_in_fit[ihit].first);
              }

              UTIL::BitField64 cellID_encoder(lcio::LCTrackerCellID::encoding_string());

              MarlinTrk::addHitNumbersToTrack(refittedTrack, all_hits, true, cellID_encoder);

              marlinTrk->getOutliers(outliers);

              for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
                all_hits.push_back(outliers[ihit].first);
              }

              MarlinTrk::addHitNumbersToTrack(refittedTrack, all_hits, false, cellID_encoder);

              delete marlinTrk;

              int nhits_in_vxd = refittedTrack->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 2];
              int nhits_in_ftd = refittedTrack->subdetectorHitNumbers()[2 * lcio::ILDDetID::FTD - 2];
              int nhits_in_sit = refittedTrack->subdetectorHitNumbers()[2 * lcio::ILDDetID::SIT - 2];
              int nhits_in_tpc = refittedTrack->subdetectorHitNumbers()[2 * lcio::ILDDetID::TPC - 2];
              int nhits_in_set = refittedTrack->subdetectorHitNumbers()[2 * lcio::ILDDetID::SET - 2];

              streamlog_out(DEBUG4) << " Hit numbers for Track " << refittedTrack->id() << ": "
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

              trackCandidates.push_back(
                  refittedTrack); // trackCandidates vector stores all the candidate tracks of the event

            } // end of the creation of the refitted track collection
          } // good fit status
        } // good initialisation status
      } // minimum acceptable TPC hits

      delete marlin_trk;

    } // for loop to the tracks

    //-------------------------------------------------------------------------------------------------------

    for (unsigned int i = 0; i < trackCandidates.size(); i++) {
      TrackImpl* SiliconRefTrack = dynamic_cast<TrackImpl*>(trackCandidates[i]);

      trackVec->addElement(SiliconRefTrack);
    }

    // for debugging reasons
    /*
    for (int ii=0; ii<usedSiHits.size(); ii++){
      streamlog_out(DEBUG0) << " ii " << ii << " hit added " << usedSiHits.at(ii) << std::endl;
    }
    */
    // PtTest->clear();
    // evt->addCollection( SiTrkCol , _siTrkColName ) ;
    evt->addCollection(trackVec, _output_track_col_name);
    // evt->addCollection( trackRelVec , _output_track_rel_name) ;
  } // track collection no empty

  ++_n_evt;

  // cout << " event " << _n_evt << std::endl ;
}

void ExtrToSIT::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void ExtrToSIT::end() {
  streamlog_out(DEBUG) << "ExtrToSIT::end()  " << name() << " processed " << _n_evt << " events in " << _n_run
                       << " runs " << std::endl;

  streamlog_out(DEBUG4) << " SIT hits considered for track-hit association " << TotalSITHits
                        << " how many of them were matched and fitted successfully ? " << SITHitsFitted
                        << " for how many the fit failed ? " << SITHitsNonFitted << std::endl;
}

void ExtrToSIT::SelectBestCandidateLimited(EVENT::TrackerHitVec& HitsInLayer, const float*& pivot,
                                           EVENT::TrackerHit*& BestHit, const FloatVec& covLCIO, double& radius,
                                           bool& BestHitFound, double& sigma, int& pointer, int& PossibleHits,
                                           float& dU, float& dV, double& DimDist, TrackerHitVec& /*usedSiHits*/) {
  BestHitFound = false;

  int NoOfHits = HitsInLayer.size();

  double MaxDistZ = sigma * (sqrt(covLCIO[9] + dV * dV));
  double MaxDistRphi = sigma * (sqrt(covLCIO[0] + dU * dU));
  // double 3Ddist = 0;

  streamlog_out(DEBUG4) << " sensor single point resolution Z : Rphi " << dV << " : " << dU << " at radius " << radius
                        << std::endl;

  for (int i = 0; i < NoOfHits; i++) {
    EVENT::TrackerHit* CandidateHit = HitsInLayer.at(i);

    streamlog_out(DEBUG1) << " Checking candidate hit " << CandidateHit << std::endl;

    double distZ = 0;
    double distRphi = 0;

    double posX = CandidateHit->getPosition()[0];
    double posY = CandidateHit->getPosition()[1];
    double posZ = CandidateHit->getPosition()[2];

    distZ = fabs(posZ - pivot[2]);
    distRphi = sqrt((posX - pivot[0]) * (posX - pivot[0]) + (posY - pivot[1]) * (posY - pivot[1]));

    streamlog_out(DEBUG2) << " maximum acceptable distances for track-hit matching Z : Rphi " << MaxDistZ << " : "
                          << MaxDistRphi << " at radius " << radius << std::endl;

    if (distZ < MaxDistZ && distRphi < MaxDistRphi) {
      PossibleHits++;
      BestHit = CandidateHit;
      BestHitFound = true;
      streamlog_out(DEBUG4) << " found appropriate hit at distance in Z : Rphi " << distZ << " : " << distRphi
                            << std::endl;
      MaxDistZ = distZ;
      MaxDistRphi = distRphi;
      pointer = i;
      DimDist = sqrt((distZ * distZ) + (distRphi * distRphi));
    }
  }

  streamlog_out(DEBUG4) << " No of candidate hits to be associated to the track = " << PossibleHits << std::endl;
}

LCCollection* ExtrToSIT::GetCollection(LCEvent* evt, std::string colName) {
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

LCRelationNavigator* ExtrToSIT::GetRelations(LCEvent* evt, std::string RelName) {
  LCRelationNavigator* nav = NULL;

  try {
    nav = new LCRelationNavigator(evt->getCollection(RelName.c_str()));
    streamlog_out(DEBUG2) << "ExtrToSIT --> " << RelName << " track relation collection in event = " << nav
                          << std::endl;
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG2) << "ExtrToSIT --> " << RelName.c_str() << " track relation collection absent in event"
                          << std::endl;
  }

  return nav;
}

int ExtrToSIT::FitInit(std::vector<TrackerHit*> trackerHits, MarlinTrk::IMarlinTrack* _marlinTrk) {
  if (trackerHits.size() < 3) {
    streamlog_out(ERROR) << "<<<<<< FitInit: Shortage of Hits! nhits = " << trackerHits.size() << " >>>>>>>"
                         << std::endl;
    return IMarlinTrack::error;
  }

  // initialise with space-points not strips
  // make a helix from 3 hits to get a trackstate
  const double* x1 = trackerHits[0]->getPosition();
  const double* x2 = trackerHits[trackerHits.size() / 2]->getPosition();
  const double* x3 = trackerHits.back()->getPosition();

  double r1 = sqrt(x1[0] * x1[0] + x1[1] * x1[1]);
  double r2 = sqrt(x2[0] * x2[0] + x2[1] * x2[1]);
  double r3 = sqrt(x3[0] * x3[0] + x3[1] * x3[1]);
  streamlog_out(DEBUG4) << " Radii of hits used for initialisation: " << r1 << ", " << r2 << " and " << r3 << std::endl;

  HelixTrack helixTrack(x1, x2, x3, _bField, HelixTrack::forwards);

  helixTrack.moveRefPoint(0.0, 0.0, 0.0);

  const float referencePoint[3] = {float(helixTrack.getRefPointX()), float(helixTrack.getRefPointY()),
                                   float(helixTrack.getRefPointZ())};

  EVENT::FloatVec covMatrix;

  covMatrix.resize(15);

  for (unsigned icov = 0; icov < covMatrix.size(); ++icov) {
    covMatrix[icov] = 0;
  }

  covMatrix[0] = (1.e6);  // sigma_d0^2
  covMatrix[2] = (1.e2);  // sigma_phi0^2
  covMatrix[5] = (1.e-4); // sigma_omega^2
  covMatrix[9] = (1.e6);  // sigma_z0^2
  covMatrix[14] = (1.e2); // sigma_tanl^2

  TrackStateImpl trackState(TrackState::AtOther, helixTrack.getD0(), helixTrack.getPhi0(), helixTrack.getOmega(),
                            helixTrack.getZ0(), helixTrack.getTanLambda(), covMatrix, referencePoint);

  _marlinTrk->initialise(trackState, _bField, IMarlinTrack::backward);

  return IMarlinTrack::success;
}

int ExtrToSIT::FitInit2(Track* track, MarlinTrk::IMarlinTrack* _marlinTrk) {
  // EVENT::FloatVec covMatrix = track->getCovMatrix();

  TrackStateImpl trackState(TrackState::AtOther, track->getD0(), track->getPhi(), track->getOmega(), track->getZ0(),
                            track->getTanLambda(), track->getCovMatrix(), track->getReferencePoint());

  _marlinTrk->initialise(trackState, _bField, IMarlinTrack::forward);

  return IMarlinTrack::success;
}
