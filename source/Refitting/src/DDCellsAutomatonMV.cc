#include "DDCellsAutomatonMV.h"

// using namespace MarlinTrk ;

#include <algorithm>
#include <climits>
#include <cmath>
#include <iostream>

#include <UTIL/BitField64.h>
// #include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>

// #include "MarlinTrk/Factory.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

using namespace lcio;
using namespace marlin;

using namespace MarlinTrk;

// Used to fedine the quality of the track output collection
const int DDCellsAutomatonMV::_output_track_col_quality_GOOD = 1;
const int DDCellsAutomatonMV::_output_track_col_quality_FAIR = 2;
const int DDCellsAutomatonMV::_output_track_col_quality_POOR = 3;

const double DDCellsAutomatonMV::TWOPI = 2 * M_PI;

DDCellsAutomatonMV aDDCellsAutomatonMV;

DDCellsAutomatonMV::DDCellsAutomatonMV() : Processor("DDCellsAutomatonMV") {
  registerProcessorParameter("NDivisionsInPhi", "Number of divisions in Phi", _nDivisionsInPhi, int(80));

  registerProcessorParameter("NDivisionsInTheta", "Number of divisions in Theta", _nDivisionsInTheta, int(80));

  FloatVec resUEx;
  resUEx.push_back(0.0040);

  registerProcessorParameter("PixRes", "VXD pixels resolution", _resU, resUEx);

  registerProcessorParameter("MultipleScatteringOn", "Use MultipleScattering in Fit", _MSOn, bool(true));

  registerProcessorParameter("EnergyLossOn", "Use Energy Loss in Fit", _ElossOn, bool(true));

  registerProcessorParameter("SmoothOn", "Smooth All Mesurement Sites in Fit", _SmoothOn, bool(true));
  /*
  registerProcessorParameter("UseMiddleLayer",
                             "Use middle layer for the cellular automaton",
                             _middleLayer,
                             bool(true));
  */
  registerProcessorParameter("UseSIT", "Use SIT", _useSIT, int(1));

  registerProcessorParameter("HitsPerTrackMin", "The minimum number of hits to create a track", _hitsPerTrackMin,
                             int(4));

  registerProcessorParameter("StepMax", "The maximum step between layers", _layerStepMax, int(2));

  registerProcessorParameter("LastLayerToIP", "The maximum step between layers", _lastLayerToIP, int(3));

  registerProcessorParameter("HelixFitMax", "The maximum chi2/Ndf that is allowed as result of a helix fit",
                             _helixFitMax, double(500));

  registerProcessorParameter("kalFitMax", "The maximum chi2/Ndf that is allowed as result of a Kalman fit",
                             _chi2OverNdfCut, double(500));

  registerProcessorParameter("Chi2ProbCut", "Tracks with a chi2 probability below this will get sorted out",
                             _chi2ProbCut, double(0.005));

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

  registerProcessorParameter("MaxHitsPerSector",
                             "Maximal number of hits allowed on a sector. More will cause drop of hits in sector",
                             _maxHitsPerSector, int(1000));

  registerProcessorParameter("BestSubsetFinder",
                             "The method used to find the best non overlapping subset of tracks. Available are: "
                             "SubsetHopfieldNN, SubsetSimple and None",
                             _bestSubsetFinder, std::string("SubsetHopfieldNN"));

  // Security checks to prevent combinatorial disasters

  registerProcessorParameter(
      "MaxConnectionsAutomaton",
      "If the automaton has more connections than this it will be redone with the next set of cut off parameters",
      _maxConnectionsAutomaton, int(100000));

  registerProcessorParameter("MaxDistance",
                             "The maximum distance between two hits in adjacent layers in order to form a minivector",
                             _maxDist, double(20));

  registerProcessorParameter(
      "MVHitsThetaDifference_Adjacent",
      "The difference in polar angle (in degrees)  between two hits in adjacent layers in order to form a minivector",
      _hitPairThDiff, double(0.5));

  registerProcessorParameter(
      "MVHitsThetaDifference",
      "The difference in polar angle (in degrees)  between two hits in the INNER layer in order to form a minivector",
      _hitPairThDiffInner, double(0.1));

  registerProcessorParameter(
      "NHitsChi2", "Maximal number of hits for which a track with n hits is better than one with n-1hits. (defaut 5)",
      _nHitsChi2, int(5));

  registerProcessorParameter("VXDName", "Name of the vertex detector element", _detElVXDName,
                             std::string("VertexBarrel"));

  registerProcessorParameter("InnerTrackerName", "Name of the inner tracker detector element", _detElITName,
                             std::string("InnerTrackerBarrel"));

  registerProcessorParameter("OuterTrackerName", "Name of the outer tracker detector element", _detElOTName,
                             std::string("OuterTrackerBarrel"));

  // Input Collections

  registerInputCollection(LCIO::TRACKERHITPLANE, "VTXHitCollectionName", "VTX Hit Collection Name", _VTXHitCollection,
                          std::string("VTXTrackerHits"));

  registerInputCollection(LCIO::TRACKERHIT, "SITHitCollectionName", "SIT Hit Collection Name", _SITHitCollection,
                          std::string("SITSpacePoints"));

  // Output Collections

  registerOutputCollection(LCIO::TRACK, "CATrackCollection",
                           "Name of the Cellular Automaton Tracking output collection", _CATrackCollection,
                           std::string("CATracks"));

  // The Criteria for the Cellular Automaton:

  std::vector<std::string> allCriteria = Criteria::getAllCriteriaNamesVec();

  registerProcessorParameter(
      "Criteria",
      "A vector of the criteria that are going to be used. For every criterion a min and max needs to be set!!!",
      _criteriaNames, allCriteria);

  // Now set min and max values for all the criteria
  for (unsigned i = 0; i < _criteriaNames.size(); i++) {
    std::vector<float> emptyVec;

    std::string critMinString = _criteriaNames[i] + "_min";

    registerProcessorParameter(critMinString, "The minimum of " + _criteriaNames[i], _critMinima[_criteriaNames[i]],
                               emptyVec);

    std::string critMaxString = _criteriaNames[i] + "_max";

    registerProcessorParameter(critMaxString, "The maximum of " + _criteriaNames[i], _critMaxima[_criteriaNames[i]],
                               emptyVec);
  }
}

void DDCellsAutomatonMV::init() {
  this->setupGeom();

  /**********************************************************************************************/
  /*       Make a SectorSystemVXD                                                             */
  /**********************************************************************************************/

  // The SectorSystemVXD is the object translating the sectors of the hits into layers, modules etc. and vice versa

  MiniVectors_sectors = 0;
  MiniVectors_CutSelection = 0;

  int nDivisionsInPhi = _nDivisionsInPhi;
  int nDivisionsInTheta = _nDivisionsInTheta;
  if (_useSIT == 1) {
    _nLayers = _nLayersVTX + _nLayersSIT + 1; // + 1 adding a virtual layer for the IP
  } else
    _nLayers = _nLayersVTX + 1; // + 1 adding a virtual layer for the IP
  _sectorSystemVXD = new SectorSystemVXD(_nLayers, nDivisionsInPhi, nDivisionsInTheta);

  _dPhi = TWOPI / _nDivisionsInPhi;
  _dTheta = 2.0 / _nDivisionsInTheta;

  /**********************************************************************************************/
  /*       Initialise the MarlinTrkSystem, needed by the tracks for fitting                     */
  /**********************************************************************************************/

  // set upt the geometry

  _trkSystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", 0, "");

  if (_trkSystem == 0)
    throw EVENT::Exception(std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("KalTest"));

  // set the options
  _trkSystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useQMS, _MSOn);           // multiple scattering
  _trkSystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, _ElossOn);       // energy loss
  _trkSystem->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, _SmoothOn); // smoothing

  // initialise the tracking system
  _trkSystem->init();
}

void DDCellsAutomatonMV::processRunHeader(LCRunHeader*) {}

void DDCellsAutomatonMV::processEvent(LCEvent* evt) {
  // set the correct configuration for the tracking system for this event
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(_trkSystem, _MSOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(_trkSystem, _ElossOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trkSystem, _SmoothOn);

  TrackerHitVec HitsTemp;     // Hits to be deleted at the end
  std::vector<IHit*> hitsTBD; // Hits to be deleted at the end
  // std::vector< IHit* > MiniVectorsTemp;

  // reset the hit map
  _map_sector_spacepoints.clear();
  _map_sector_hits.clear();

  InitialiseVTX(evt, HitsTemp);

  unsigned round = 0; // the round we are in
  std::vector<RawTrack> rawTracks;

  /**********************************************************************************************/
  /*                Check if no sector is overflowing with hits                                 */
  /**********************************************************************************************/

  std::map<int, EVENT::TrackerHitVec>::iterator it;

  for (it = _map_sector_spacepoints.begin(); it != _map_sector_spacepoints.end(); it++) {
    int nHits = it->second.size();
    streamlog_out(DEBUG2) << "Number of hits in sector " << it->first << " = " << nHits << "\n";

    if (nHits > _maxHitsPerSector) {
      it->second.clear(); // delete the hits in this sector, it will be dropped

      streamlog_out(ERROR)
          << " ### EVENT " << evt->getEventNumber() << " :: RUN " << evt->getRunNumber()
          << " \n ### Number of Hits in VXD Sector " << it->first << ": " << nHits << " > " << _maxHitsPerSector
          << " (MaxHitsPerSector)\n : This sector will be dropped from track search, and QualityCode set to \"Poor\" "
          << std::endl;

      _output_track_col_quality =
          _output_track_col_quality_POOR; // We had to drop hits, so the quality of the result is decreased
    }
  }

  /**********************************************************************************************/
  /*                Add the IP as a virtual hit                                                 */
  /**********************************************************************************************/

  IHit* virtualIPHitForward = createVirtualIPHit(_sectorSystemVXD);
  // HitsTemp.push_back( virtualIPHitForward );
  _map_sector_hits[0].push_back(virtualIPHitForward);
  streamlog_out(DEBUG4) << " sector of the IP hit = " << virtualIPHitForward->getSector() << std::endl;

  //*********************************************************************************************************************

  /**********************************************************************************************/
  /*                        Create Mini-Vectors                                                 */
  /**********************************************************************************************/

  for (std::map<int, EVENT::TrackerHitVec>::iterator itSecHit = _map_sector_spacepoints.begin();
       itSecHit != _map_sector_spacepoints.end(); itSecHit++) { // over all sectors

    int sector = itSecHit->first;
    CreateMiniVectors(sector); // Process one VXD sector
  }

  while (setCriteria(round)) {
    streamlog_out(DEBUG4) << " DO I ENTER IN THE GAME " << std::endl;

    round++; // count up the round we are in

    SegmentBuilder segBuilder(_map_sector_hits);

    segBuilder.addCriteria(_crit2Vec); // Add the criteria on when to connect two hits. The vector has been filled by
                                       // the method setCriteria

    unsigned layerStepMax = _layerStepMax; // how many layers to go at max
    unsigned lastLayerToIP = _lastLayerToIP;
    // bool MiddleLayer = _middleLayer ;
    // VXDSectorConnector secCon( _sectorSystemVXD , layerStepMax, lastLayerToIP, MiddleLayer );
    VXDSectorConnector secCon(_sectorSystemVXD, layerStepMax, lastLayerToIP);

    streamlog_out(DEBUG4) << " _layerStepMax =  " << _layerStepMax << std::endl;
    streamlog_out(DEBUG4) << " _lastLayerToIP = " << _lastLayerToIP << std::endl;

    segBuilder.addSectorConnector(&secCon); // Add the sector connector (so the SegmentBuilder knows what hits from
                                            // different sectors it is allowed to look for connections)

    // And get out the Cellular Automaton with the 1-segments
    Automaton automaton = segBuilder.get1SegAutomaton();

    streamlog_out(DEBUG4) << " automaton.getNumberOfConnections() = " << automaton.getNumberOfConnections()
                          << std::endl;

    // Check if there are not too many connections
    if (automaton.getNumberOfConnections() > unsigned(_maxConnectionsAutomaton)) {
      streamlog_out(DEBUG4) << "Redo the Automaton with different parameters, because there are too many connections:\n"
                            << "\tconnections( " << automaton.getNumberOfConnections()
                            << " ) > MaxConnectionsAutomaton( " << _maxConnectionsAutomaton << " )\n";
      continue;
    }

    /**********************************************************************************************/
    /*                Automaton                                                                   */
    /**********************************************************************************************/

    if (automaton.getNumberOfConnections() > 0) {
      streamlog_out(DEBUG4) << "\t\t---Automaton---\n";

      // if( _useCED ) KiTrackMarlin::drawAutomatonSegments( automaton ); // draws the 1-segments (i.e. hits)

      /*******************************/
      /*      2-hit segments         */
      /*******************************/

      streamlog_out(DEBUG4) << "\t\t--2-hit-Segments--\n";

      //          streamlog_out(DEBUG4) << "Automaton has " << automaton.getTracks( 3 ).size() << " track candidates\n";
      //          //should be commented out, because it takes time

      automaton.clearCriteria();
      automaton.addCriteria(_crit3Vec); // Add the criteria for 3 hits (i.e. 2 2-hit segments )

      // Let the automaton lengthen its 1-hit-segments to 2-hit-segments
      automaton.lengthenSegments();

      // So now we have 2-hit-segments and are ready to perform the Cellular Automaton.

      // Perform the automaton
      automaton.doAutomaton();

      // Clean segments with bad states
      automaton.cleanBadStates();

      // Reset the states of all segments
      automaton.resetStates();

      //          streamlog_out(DEBUG4) << "Automaton has " << automaton.getTracks( 3 ).size() << " track candidates\n";
      //          //should be commented out, because it takes time

      // Check if there are not too many connections
      if (automaton.getNumberOfConnections() > unsigned(_maxConnectionsAutomaton)) {
        streamlog_out(DEBUG4)
            << "Redo the Automaton with different parameters, because there are too many connections:\n"
            << "\tconnections( " << automaton.getNumberOfConnections() << " ) > MaxConnectionsAutomaton( "
            << _maxConnectionsAutomaton << " )\n";
        continue;
      }

      streamlog_out(DEBUG4) << " _hitsPerTrackMin = " << _hitsPerTrackMin << std::endl;

      // get the raw tracks (raw track = just a vector of hits, the most rudimentary form of a track)
      rawTracks = automaton.getTracks(_hitsPerTrackMin);

      break; // if we reached this place all went well and we don't need another round --> exit the loop
    }
  }

  streamlog_out(DEBUG4) << "Automaton returned " << rawTracks.size() << " raw tracks \n";

  //*************************************************************************************************************
  // Track fitting similar to forward tracking
  //*************************************************************************************************************

  std::vector<ITrack*> trackCandidates;

  // for all raw tracks we got from the automaton
  for (unsigned i = 0; i < rawTracks.size(); i++) {
    RawTrack rawTrack = rawTracks[i];

    if (rawTrack.size() < unsigned(_hitsPerTrackMin)) {
      streamlog_out(DEBUG4) << "Trackversion discarded, too few hits: only " << rawTrack.size() << " < "
                            << _hitsPerTrackMin << "(hitsPerTrackMin)\n";
      continue;
    }

    VXDTrack* trackCand = new VXDTrack(_trkSystem);

    // add the hits to the track
    for (unsigned k = 0; k < rawTrack.size(); k++) {
      IMiniVector* mvHit = dynamic_cast<IMiniVector*>(rawTrack[k]); // cast to IMiniVectors, as needed for a VXDTrack
      if (mvHit != NULL)
        trackCand->addHit(mvHit);
      else
        streamlog_out(DEBUG4) << "Hit " << rawTrack[k] << " could not be casted to IMiniVector\n";
    }

    /*-----------------------------------------------*/
    /*                Helix Fit                      */
    /*-----------------------------------------------*/

    streamlog_out(DEBUG2) << "Fitting with Helix Fit\n";
    try {
      VXDHelixFitter helixFitter(trackCand->getLcioTrack());

      TrackerHitVec testVec = trackCand->getLcioTrack()->getTrackerHits();

      streamlog_out(DEBUG2) << " $$$$ fitting track with " << testVec.size() << " hits " << std::endl;

      float chi2OverNdf = helixFitter.getChi2() / float(helixFitter.getNdf());
      streamlog_out(DEBUG2) << "chi2OverNdf = " << chi2OverNdf << "\n";

      if (chi2OverNdf > _helixFitMax) {
        streamlog_out(DEBUG2) << "Discarding track because of bad helix fit: chi2/ndf = " << chi2OverNdf << "\n";

        // debug
        // streamlog_out(DEBUG4) << " pre-fitting: deleting track " << trackCand << std::endl ;
        delete trackCand;
        continue;

      } else {
        streamlog_out(DEBUG2) << "Keeping track because of good helix fit: chi2/ndf = " << chi2OverNdf << "\n";
      }
    } catch (VXDHelixFitterException& e) {
      streamlog_out(DEBUG2) << "Track rejected, because fit failed: " << e.what() << "\n";
      delete trackCand;
      continue;
    }

    /*-----------------------------------------------*/
    /*                Kalman Fit                      */
    /*-----------------------------------------------*/

    streamlog_out(DEBUG3) << "Fitting with Kalman Filter\n";
    try {
      trackCand->fit();

      streamlog_out(DEBUG3) << " Track " << trackCand << " chi2Prob = " << trackCand->getChi2Prob()
                            << "( chi2=" << trackCand->getChi2() << ", Ndf=" << trackCand->getNdf() << " )\n";

      double test = trackCand->getChi2() / (1.0 * trackCand->getNdf()); // FIXME: give a less dull name to this var. YV
      float NoOfHitsTimes2 = 2.0 * (trackCand->getHits().size());

      // if ( trackCand->getChi2Prob() >= _chi2ProbCut ){
      if (test < _chi2OverNdfCut) {
        streamlog_out(DEBUG2) << "Track accepted (chi2prob " << trackCand->getChi2Prob() << " >= " << _chi2ProbCut
                              << " chi2 over ndf " << test << " No of hits (x4) " << NoOfHitsTimes2
                              << " Overall sorting variable " << NoOfHitsTimes2 / test << "\n";

      } else {
        streamlog_out(DEBUG2) << "Track rejected (chi2prob " << trackCand->getChi2Prob() << " < " << _chi2ProbCut
                              << " chi2 over ndf " << test << "\n";

        // debug
        // streamlog_out(DEBUG4) << " Kalman fitting: deleting track " << trackCand << std::endl ;

        delete trackCand;

        continue;
      }

    } catch (FitterException& e) {
      streamlog_out(DEBUG4) << "Track rejected, because fit failed: " << e.what() << "\n";
      delete trackCand;
      continue;
    }

    // Kalman fitting over
    //____________________________________________________________________________________________________________

    streamlog_out(DEBUG1) << "------------ trackCand = " << trackCand << "\n";

    trackCandidates.push_back(trackCand);
  }

  // FTD like track fitting over.
  //_________________________________________________________________________________________________

  // establish the track collection that will be created
  LCCollectionVec* trackVec = new LCCollectionVec(LCIO::TRACK);

  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  trackVec->setFlag(trkFlag.getFlag());

  // Selection of best consistent track subsample
  //*************************************************************************************************************

  std::vector<ITrack*> GoodTracks;
  std::vector<ITrack*> RejectedTracks;

  TrackCompatibilityShare1_MV comp;

  // Various ways to define the quality of a track. They are defined in  the header file
  // TrackQIChi2Prob_MV whatever;
  // TrackQISpecial_MV JustDoIt ;
  // TrackQI trackQI;
  // MaxHits MaxLength;
  Test test;

  streamlog_out(DEBUG4) << " best subset finder = " << _bestSubsetFinder << " no of tracks fed to the nnets "
                        << trackCandidates.size() << std::endl;

  // FIX ME: implement properly HopfieldNN also
  /*
  if( _bestSubsetFinder == "SubsetHopfieldNN" ){

    streamlog_out(DEBUG4) << "Use SubsetHopfieldNN for getting the best subset\n" << std::endl ;

    SubsetHopfieldNN< ITrack* > subset_tracks ;
    subset_tracks.add( trackCandidates );
    subset_tracks.calculateBestSet( comp, trackQI );
    GoodTracks = subset_tracks.getAccepted();
    RejectedTracks = subset_tracks.getRejected();

  }
  */

  if (_bestSubsetFinder == "SubsetSimple") {
    SubsetSimple<ITrack*> subset_tracks;
    subset_tracks.add(trackCandidates);
    // subset_tracks.calculateBestSet( comp, trackQI );
    // subset_tracks.calculateBestSet( comp, MaxLength );
    subset_tracks.calculateBestSet(comp, test);
    GoodTracks = subset_tracks.getAccepted();
    RejectedTracks = subset_tracks.getRejected();
  }

  if (_bestSubsetFinder == "NoSelection") { // in any other case take all tracks

    streamlog_out(DEBUG3) << "Input for subset = \"" << _bestSubsetFinder << "\". All tracks are kept\n";
    GoodTracks = trackCandidates;
  }

  streamlog_out(DEBUG4) << "End of Sorting, Good tracks number: " << GoodTracks.size() << std::endl;

  //******************************************************************************************************************
  // best consistent track subsample selection ends here

  // Finalise the tracks

  for (unsigned int i = 0; i < GoodTracks.size(); i++) {
    VXDTrack* myTrack = dynamic_cast<VXDTrack*>(GoodTracks[i]);

    if (myTrack != NULL) {
      TrackImpl* trackImpl = new TrackImpl(*(myTrack->getLcioTrack())); // possible leak

      try {
        finaliseTrack(trackImpl);

        // Applying a x2/ndf cut on final tracks

        if (((1.0 * trackImpl->getChi2()) / (1.0 * trackImpl->getNdf())) < 10.0) {
          trackVec->addElement(trackImpl);

          streamlog_out(DEBUG0) << "DDCellsAutomatonMV: trackImpl added to trackVec\n";
        }
      }

      catch (FitterException& e) {
        streamlog_out(DEBUG4) << "DDCellsAutomatonMV: track couldn't be finalized due to fitter error: " << e.what()
                              << "\n";
        delete trackImpl;
      }
    }
  }
  // Finalisation ends

  streamlog_out(DEBUG4) << "DDCellsAutomatonMV: _CATrackCollection = " << _CATrackCollection
                        << "     trackVec->getNumberOfElements() = " << trackVec->getNumberOfElements() << "\n";

  evt->addCollection(trackVec, _CATrackCollection);

  /**********************************************************************************************/
  /*                Clean up                                                                    */
  /**********************************************************************************************/

  // delete all the created IHits
  for (unsigned i = 0; i < HitsTemp.size(); i++)
    delete HitsTemp[i];

  if (MiniVectorsTemp.size() > 0) {
    for (unsigned i = 0; i < MiniVectorsTemp.size(); i++) {
      delete MiniVectorsTemp[i];
      MiniVectorsTemp[i] = NULL;
    }
  }

  if (TestMiniVectorsTemp.size() > 0) {
    for (unsigned i = 0; i < TestMiniVectorsTemp.size(); i++) {
      delete TestMiniVectorsTemp[i];
      TestMiniVectorsTemp[i] = NULL;
    }
  }

  // cleanup of tracks
  // if ( _bestSubsetFinder != "NoSelection") for (unsigned int i=0; i < GoodTracks.size(); i++){ delete GoodTracks[i];
  // } for ( unsigned i=0; i<RejectedTracks.size(); i++){ delete RejectedTracks[i]; } for ( unsigned i=0;
  // i<trackCandidates.size(); i++){ delete trackCandidates[i]; }
}

void DDCellsAutomatonMV::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void DDCellsAutomatonMV::end() {
  for (unsigned i = 0; i < _crit2Vec.size(); i++)
    delete _crit2Vec[i];
  for (unsigned i = 0; i < _crit3Vec.size(); i++)
    delete _crit3Vec[i];
  for (unsigned i = 0; i < _crit4Vec.size(); i++)
    delete _crit4Vec[i];
  _crit2Vec.clear();
  _crit3Vec.clear();
  _crit4Vec.clear();

  delete _sectorSystemVXD;
  _sectorSystemVXD = NULL;

  // streamlog_out(DEBUG4) << " no of candidate minivectors created with sectors apprach " << MiniVectors_sectors << "
  // no of candidate minivectors created applying a cut selection " <<  MiniVectors_CutSelection  << std::endl ;
}

void DDCellsAutomatonMV::InitialiseVTX(LCEvent* evt, EVENT::TrackerHitVec HitsTemp) {
  // Reading out VTX Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    LCCollection* hitCollection = evt->getCollection(_VTXHitCollection.c_str());

    int nelem = hitCollection->getNumberOfElements();

    streamlog_out(DEBUG4) << "Number of VTX hits = " << nelem << std::endl;

    for (int ielem = 0; ielem < nelem; ++ielem) {
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem));

      dd4hep::rec::Vector3D U(1.0, hit->getU()[1], hit->getU()[0], dd4hep::rec::Vector3D::spherical);
      dd4hep::rec::Vector3D V(1.0, hit->getV()[1], hit->getV()[0], dd4hep::rec::Vector3D::spherical);
      dd4hep::rec::Vector3D Z(0.0, 0.0, 1.0);

      const float eps = 1.0e-07;
      // V must be the global z axis
      if (fabs(1.0 - V.dot(Z)) > eps) {
        streamlog_out(ERROR) << "SiTrkHopfieldNN: VXD Hit measurment vectors V is not equal to the global Z axis. \n\n "
                                " exit(1) called from file "
                             << __FILE__ << " and line " << __LINE__ << std::endl;
        exit(1);
      }

      if (fabs(U.dot(Z)) > eps) {
        streamlog_out(ERROR) << "SiTrkHopfieldNN: VXD Hit measurment vectors U is not in the global X-Y plane. \n\n "
                                "exit(1) called from file "
                             << __FILE__ << " and line " << __LINE__ << std::endl;
        exit(1);
      }

      double pos[3];
      double radius = 0;

      for (int i = 0; i < 3; ++i) {
        pos[i] = hit->getPosition()[i];
        radius += pos[i] * pos[i];
      }

      radius = sqrt(radius);

      double cosTheta = pos[2] / radius;
      double Phi = atan2(pos[1], pos[0]);

      if (Phi < 0.)
        Phi = Phi + TWOPI;

      // get the layer number

      int celId = hit->getCellID0();

      UTIL::BitField64 encoder(lcio::LCTrackerCellID::encoding_string());
      encoder.setValue(celId);
      int layer = encoder[lcio::LCTrackerCellID::layer()] + 1; // +1 cause IP is considered layer 0

      int iPhi = int(Phi / _dPhi);
      int iTheta = int((cosTheta + double(1.0)) / _dTheta);
      int iCode = layer + _nLayers * iPhi + _nLayers * _nDivisionsInPhi * iTheta;

      streamlog_out(DEBUG4) << " CA: making a VXD hit at layer " << layer << " phi sector " << iPhi << " theta sector "
                            << iTheta << " theta angle " << acos(cosTheta) * (180.0 / M_PI) << " sector code " << iCode
                            << " total layers " << _nLayers << " no of phi sectors " << _nDivisionsInPhi
                            << " no of theta sectors " << _nDivisionsInTheta << std::endl;

      // Make an VXDHit01 from the TrackerHit
      // VXDHit01* vxdHit = new VXDHit01 ( hit , _sectorSystemVXD );   // Don't need to create VXDHits, we stick to mini
      // - vectors
      HitsTemp.push_back(hit); // so we can easily delete every created hit afterwards

      _map_sector_spacepoints[iCode].push_back(hit);
    }

  }

  catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << " collection not found : " << _VTXHitCollection.c_str() << std::endl;
  }

  //________________________________________________________________________________________________________________
  // readout of SIT collection                                                                                     |
  //_______________________________________________________________________________________________________________|

  if (_useSIT == 1) {
    try {
      LCCollection* hitCollection = evt->getCollection(_SITHitCollection.c_str());

      int nelem = hitCollection->getNumberOfElements();

      streamlog_out(DEBUG4) << "Number of SIT hits = " << nelem << std::endl;

      TrackerHit* trkhit = 0;
      TrackerHitPlane* trkhit_P = 0;
      TrackerHitZCylinder* trkhit_C = 0;

      // double drphi(NAN);
      // double dz(NAN);

      for (int ielem = 0; ielem < nelem; ++ielem) {
        // hit could be of the following type
        // 1) TrackerHit, either ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT or just standard TrackerHit
        // 2) TrackerHitPlane, either 1D or 2D
        // 3) TrackerHitZCylinder, if coming from a simple cylinder design as in the LOI

        // Establish which of these it is in the following order of likelyhood
        //    i)   ILDTrkHitTypeBit::ONE_DIMENSIONAL (TrackerHitPlane) Should Never Happen: SpacePoints Must be Used
        //    Instead ii)  ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT (TrackerHit) iii) TrackerHitPlane (Two dimentional)
        //    iv)  TrackerHitZCylinder
        //    v)   Must be standard TrackerHit

        trkhit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));

        int celId_SIT = trkhit->getCellID0();

        UTIL::BitField64 encoder(lcio::LCTrackerCellID::encoding_string());
        encoder.setValue(celId_SIT);
        int layer = encoder[lcio::LCTrackerCellID::layer()] + 1; // + 1 cause IP is considered layer 0

        // VXD and SIT are treated as one system so SIT layers start from _nLayersVTX
        layer = layer + _nLayersVTX;

        if (layer < 0 || layer >= _nLayers) {
          streamlog_out(ERROR) << "SiliconTracking_MarlinTrk => fatal error in SIT : layer is outside allowed range : "
                               << layer << std::endl;
          exit(1);
        }

        // first check that we have not been given 1D hits by mistake, as they won't work here
        if (BitSet32(trkhit->getType())[UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL]) {
          streamlog_out(ERROR)
              << "SiliconTracking_MarlinTrk: SIT Hit cannot be of type UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL "
                 "COMPOSITE SPACEPOINTS must be use instead. \n\n  exit(1) called from file "
              << __FILE__ << " and line " << __LINE__ << std::endl;
          exit(1);

        }
        // most likely case: COMPOSITE_SPACEPOINT hits formed from stereo strip hits
        else if (BitSet32(trkhit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) {
          streamlog_out(DEBUG1) << " We deal with composite spacepoints " << std::endl;

          // drphi =  2 * sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]);
          // dz    =      sqrt(trkhit->getCovMatrix()[5]);

        }
        // or a PIXEL based SIT, using 2D TrackerHitPlane like the VXD above
        else if ((trkhit_P = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem)))) {
          // first we need to check if the measurement vectors are aligned with the global coordinates
          dd4hep::rec::Vector3D U(1.0, trkhit_P->getU()[1], trkhit_P->getU()[0], dd4hep::rec::Vector3D::spherical);
          dd4hep::rec::Vector3D V(1.0, trkhit_P->getV()[1], trkhit_P->getV()[0], dd4hep::rec::Vector3D::spherical);
          dd4hep::rec::Vector3D Z(0.0, 0.0, 1.0);

          const float eps = 1.0e-07;
          // V must be the global z axis
          if (fabs(1.0 - V.dot(Z)) > eps) {
            streamlog_out(ERROR) << "SiliconTracking_MarlinTrk: PIXEL SIT Hit measurment vectors V is not equal to the "
                                    "global Z axis. \n\n  exit(1) called from file "
                                 << __FILE__ << " and line " << __LINE__ << std::endl;
            exit(1);
          }

          // U must be normal to the global z axis
          if (fabs(U.dot(Z)) > eps) {
            streamlog_out(ERROR) << "SiliconTracking_MarlinTrk: PIXEL SIT Hit measurment vectors U is not in the "
                                    "global X-Y plane. \n\n exit(1) called from file "
                                 << __FILE__ << " and line " << __LINE__ << std::endl;
            exit(1);
          }

          // drphi = trkhit_P->getdU();
          // dz    = trkhit_P->getdV();

        }
        // or a simple cylindrical design, as used in the LOI
        else if ((trkhit_C = dynamic_cast<TrackerHitZCylinder*>(hitCollection->getElementAt(ielem)))) {
          // drphi = trkhit_C->getdRPhi();
          // dz    = trkhit_C->getdZ();

        }
        // this would be very unlikely, but who knows ... just an ordinary TrackerHit, which is not a
        // COMPOSITE_SPACEPOINT
        else {
          // drphi =  2 * sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]);
          // dz =     sqrt(trkhit->getCovMatrix()[5]);
        }

        double pos[3];
        double radius = 0;

        for (int i = 0; i < 3; ++i) {
          pos[i] = trkhit->getPosition()[i];
          radius += pos[i] * pos[i];
        }

        radius = sqrt(radius);

        double cosTheta = pos[2] / radius;
        double Phi = atan2(pos[1], pos[0]);

        if (Phi < 0.)
          Phi = Phi + TWOPI;

        int iPhi = int(Phi / _dPhi);
        int iTheta = int((cosTheta + double(1.0)) / _dTheta);
        int iCode = layer + _nLayers * iPhi + _nLayers * _nDivisionsInPhi * iTheta;

        streamlog_out(DEBUG2) << " CA: making an SIT hit at layer " << layer << " phi sector " << iPhi
                              << " theta sector " << iTheta << " sector code " << iCode << " total layers " << _nLayers
                              << " no of phi sectors " << _nDivisionsInPhi << std::endl;

        // Make an VXDHit01 from the TrackerHit
        // VXDHit01* vxdHit = new VXDHit01 ( trkhit , _sectorSystemVXD );   // Don't need to create VXDHits, we stick to
        // mini - vectors
        HitsTemp.push_back(trkhit); // so we can easily delete every created hit afterwards

        _map_sector_spacepoints[iCode].push_back(trkhit);
      }

    } catch (DataNotAvailableException& e) {
      streamlog_out(DEBUG4) << " collection not found : " << _SITHitCollection.c_str() << std::endl;
    }
  }
}

//________________________________________________________________________________________________________
//                                                                                                       |
// Create all possible mini-vectors of the sector                                                        |
//                                                                                                       |
//--------------------------------------------------------------------------------------------------------

// Taking the hits on the outer side of each VXD double layer we try to connect them with hits on the inner side.
// The connections considered valid if it satisfies a condition (for example the dtheta difference of the hits or their
// distance didvided with  the distance of the two sides of the layer. The search is confined in  a number of target
// sectors of the inner side of the VXD layer.

void DDCellsAutomatonMV::CreateMiniVectors(int sector) {
  int thisTheta = sector / (_nLayers * _nDivisionsInPhi);
  int thisPhi = ((sector - (thisTheta * _nLayers * _nDivisionsInPhi)) / _nLayers);
  int layer = sector - (thisTheta * _nLayers * _nDivisionsInPhi) - (thisPhi * _nLayers);

  streamlog_out(DEBUG4) << " Taking sector " << sector << " of layer " << layer << " Phi sector " << thisPhi
                        << " Theta sector " << thisTheta << std::endl;

  int iPhi_Up = thisPhi + 2;
  int iPhi_Low = thisPhi - 2;
  int iTheta_Up = thisTheta + 2;
  int iTheta_Low = thisTheta - 2;
  if (iTheta_Low < 0)
    iTheta_Low = 0;
  if (iTheta_Up >= _nDivisionsInTheta)
    iTheta_Up = _nDivisionsInTheta - 1;

  TrackerHitVec VXDHits = _map_sector_spacepoints[sector];

  for (TrackerHitVec::iterator iter = VXDHits.begin(); iter != VXDHits.end(); ++iter) {
    TrackerHit* fromHit = *iter; // Starting hit

    streamlog_out(DEBUG2) << " hit to initiate a MV: " << fromHit << std::endl;

    int celID = fromHit->getCellID0();
    int detID = 0;
    // int layerID = 0 ;

    UTIL::BitField64 encoder(lcio::LCTrackerCellID::encoding_string());
    encoder.reset(); // reset to 0
    encoder.setValue(celID);

    detID = encoder[lcio::LCTrackerCellID::subdet()];

    layer = encoder[lcio::LCTrackerCellID::layer()] + 1; // + 1 if we consider the IP hit

    if (detID == lcio::ILDDetID::VXD) {
      // if (layer==5 || layer==3 || layer==1) {
      if (layer == 6 || layer == 4 || layer == 2) { // in case of considering an IP hit

        for (int iPhi = iPhi_Low; iPhi < iPhi_Up; iPhi++) {
          // construct mini-vectors applying a delta theta cut
          //***************************************************************************

          int iTheta_Up_mod = thisTheta + 2;
          int iTheta_Low_mod = thisTheta - 2;
          if (iTheta_Low_mod < 0)
            iTheta_Low_mod = 0;
          if (iTheta_Up_mod >= _nDivisionsInTheta)
            iTheta_Up_mod = _nDivisionsInTheta - 1;

          for (int iTheta = iTheta_Low_mod; iTheta < iTheta_Up_mod; iTheta++) {
            int target_sector = (layer - 1) + _nLayers * iPhi + _nLayers * _nDivisionsInPhi * iTheta;
            TrackerHitVec targetHitsMod = _map_sector_spacepoints[target_sector];

            streamlog_out(DEBUG2) << " Checking with TARGET sector " << target_sector << " of layer " << layer - 1
                                  << " Phi sector " << iPhi << " Theta sector " << iTheta << std::endl;

            streamlog_out(DEBUG2) << " No of hits in the sector " << targetHitsMod.size() << std::endl;

            for (TrackerHitVec::iterator iterMod = targetHitsMod.begin(); iterMod != targetHitsMod.end(); ++iterMod) {
              TrackerHit* toHitMod = *iterMod; // Candidate hit to form a mini - vector with the starting hit

              if (thetaAgreementImproved(toHitMod, fromHit, layer) == true) {
                MiniVector* testMiniVectorMod = new MiniVector(fromHit, toHitMod);
                MiniVectorHit01* miniVectorHitMod = new MiniVectorHit01(testMiniVectorMod, _sectorSystemVXD);
                MiniVectors_CutSelection++;
                _map_sector_hits[sector].push_back(miniVectorHitMod);

                streamlog_out(DEBUG2) << " making the VXD mini-vector hit " << miniVectorHitMod << " at sector "
                                      << sector << std::endl;

                MiniVectorsTemp.push_back(miniVectorHitMod);
                TestMiniVectorsTemp.push_back(testMiniVectorMod);
              }
            }

            targetHitsMod.clear();
          }
        }
      }
    } // end of VXD hits - creation of mini-vectors

    // **** SIT ****

    // We don't use SIT for time being. FIX ME

    if (_useSIT == 1) {
      if (detID == lcio::ILDDetID::SIT) {
        layer = layer + _nLayersVTX;

        // if (layer==8) {
        if (layer == 7 || layer == 9) {
          // double from_x = fromHit->getPosition()[0] ;
          // double from_y = fromHit->getPosition()[1] ;

          for (int iPhi = iPhi_Low; iPhi < iPhi_Up; iPhi++) {
            int iTheta_Up_mod = thisTheta + 10;
            int iTheta_Low_mod = thisTheta - 10;
            if (iTheta_Low_mod < 0)
              iTheta_Low_mod = 0;
            if (iTheta_Up_mod >= _nDivisionsInTheta)
              iTheta_Up_mod = _nDivisionsInTheta - 1;

            for (int iTheta = iTheta_Low; iTheta < iTheta_Up; iTheta++) {
              int target_sector = (layer - 1) + _nLayers * iPhi + _nLayers * _nDivisionsInPhi * iTheta;

              // double ThetaAngle = acos(((2*iTheta)/(_nDivisionsInTheta*1.0)) - 1);

              TrackerHitVec targetHits = _map_sector_spacepoints[target_sector];

              streamlog_out(DEBUG3) << " How many hits on the target sector " << targetHits.size() << std::endl;

              for (TrackerHitVec::iterator iter2 = targetHits.begin(); iter2 != targetHits.end(); ++iter2) {
                TrackerHit* toHit = *iter2;

                if (Dist(fromHit, toHit) < _maxDist) {
                  MiniVector* sitTestMiniVector = new MiniVector(fromHit, toHit);
                  MiniVectorHit01* sitMiniVectorHit = new MiniVectorHit01(sitTestMiniVector, _sectorSystemVXD);

                  streamlog_out(DEBUG4) << " making a SIT mini-vector hit at sector " << sector << std::endl;

                  // MiniVectors_sectors++;

                  _map_sector_hits[sector].push_back(sitMiniVectorHit);
                }
              }

              targetHits.clear();
            }
          }
        }
      }
    }

  } // end of looping on VXD -  SIT trackerhits

  VXDHits.clear();
}

void DDCellsAutomatonMV::setupGeom() {
  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  const double pos[3] = {0, 0, 0};
  double bFieldVec[3] = {0, 0, 0};
  theDetector.field().magneticField(pos, bFieldVec); // get the magnetic field vector from DD4hep
  _bField = bFieldVec[2] / dd4hep::tesla;            // z component at (0,0,0)

  _nLayersVTX = 0;
  dd4hep::DetElement vtxDE = theDetector.detector(_detElVXDName);
  dd4hep::rec::ZPlanarData* vtx = vtxDE.extension<dd4hep::rec::ZPlanarData>();
  _nLayersVTX = vtx->layers.size();

  _nLayersSIT = 0;

  int nIT = 0;
  dd4hep::DetElement itDE = theDetector.detector(_detElITName);
  dd4hep::rec::ZPlanarData* it = itDE.extension<dd4hep::rec::ZPlanarData>();
  nIT = it->layers.size();

  int nOT = 0;
  dd4hep::DetElement otDE = theDetector.detector(_detElOTName);
  dd4hep::rec::ZPlanarData* ot = otDE.extension<dd4hep::rec::ZPlanarData>();
  nOT = ot->layers.size();

  _nLayersSIT = nIT + nOT;

  streamlog_out(DEBUG0) << "###setupGeom: _bField = " << _bField << std::endl;
  streamlog_out(DEBUG0) << "###setupGeom: _nLayersVTX = " << _nLayersVTX << std::endl;
  streamlog_out(DEBUG0) << "###setupGeom: _nLayersSilicon (Inner + Outer) = " << _nLayersSIT << std::endl;
}

bool DDCellsAutomatonMV::setCriteria(unsigned round) {
  // delete the old ones
  for (unsigned i = 0; i < _crit2Vec.size(); i++)
    delete _crit2Vec[i];
  for (unsigned i = 0; i < _crit3Vec.size(); i++)
    delete _crit3Vec[i];
  for (unsigned i = 0; i < _crit4Vec.size(); i++)
    delete _crit4Vec[i];
  _crit2Vec.clear();
  _crit3Vec.clear();
  _crit4Vec.clear();

  streamlog_out(DEBUG4) << " Calling setCriteria function - criteria vector size =  " << _criteriaNames.size()
                        << std::endl;

  bool newValuesGotUsed = false; // if new values are used

  for (unsigned i = 0; i < _criteriaNames.size(); i++) {
    std::string critName = _criteriaNames[i];

    float min = _critMinima[critName].back();
    float max = _critMaxima[critName].back();

    streamlog_out(DEBUG4) << " iterator " << i << " criterio " << critName << " min. value " << min << " max. value "
                          << max << std::endl;

    // use the value corresponding to the round, if there are no new ones for this criterion, just do nothing (the
    // previous value stays in place)
    if (round + 1 <= _critMinima[critName].size()) {
      min = _critMinima[critName][round];
      newValuesGotUsed = true;
    }

    if (round + 1 <= _critMaxima[critName].size()) {
      max = _critMaxima[critName][round];
      newValuesGotUsed = true;
    }

    ICriterion* crit = Criteria::createCriterion(critName, min, max);

    // Some debug output about the created criterion
    std::string type = crit->getType();

    streamlog_out(DEBUG4) << "Added: Criterion " << critName << " (type =  " << type << " ). Min = " << min
                          << ", Max = " << max << ", round " << round << "\n";

    // Add the new criterion to the corresponding vector
    if (type == "2Hit") {
      _crit2Vec.push_back(crit);

    } else if (type == "3Hit") {
      _crit3Vec.push_back(crit);

    } else if (type == "4Hit") {
      _crit4Vec.push_back(crit);

    } else
      delete crit;
  }

  return newValuesGotUsed;
}

void DDCellsAutomatonMV::finaliseTrack(TrackImpl* trackImpl) {
  // Fitter fitter( trackImpl , _trkSystem ); //it gives problem at 90deg: sometimes the hits are taken in the inverse
  // order resulting in a fitted track in the opposite quadrant of the hits
  Fitter fitter(trackImpl, _trkSystem, 1); // it forces the hits ordering according to the radius (problem for very bent
                                           // tracks that are coming back - TO STUDY)

  trackImpl->trackStates().clear();

  TrackStateImpl* trkStateIP = new TrackStateImpl(*fitter.getTrackState(lcio::TrackState::AtIP));
  trkStateIP->setLocation(TrackState::AtIP);
  trackImpl->addTrackState(trkStateIP);

  TrackStateImpl* trkStateFirstHit = new TrackStateImpl(*fitter.getTrackState(TrackState::AtFirstHit));
  trkStateFirstHit->setLocation(TrackState::AtFirstHit);
  trackImpl->addTrackState(trkStateFirstHit);

  TrackStateImpl* trkStateLastHit = new TrackStateImpl(*fitter.getTrackState(TrackState::AtLastHit));
  trkStateLastHit->setLocation(TrackState::AtLastHit);
  trackImpl->addTrackState(trkStateLastHit);

  // TrackStateImpl* trkStateAtCalo = new TrackStateImpl( *fitter.getTrackState( TrackState::AtCalorimeter ) ) ;
  // trkStateAtCalo->setLocation( TrackState::AtCalorimeter );
  // trackImpl->addTrackState( trkStateAtCalo );

  trackImpl->setChi2(fitter.getChi2(TrackState::AtIP));
  trackImpl->setNdf(fitter.getNdf(TrackState::AtIP));

  const float* p = trkStateFirstHit->getReferencePoint();
  trackImpl->setRadiusOfInnermostHit(sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]));

  std::map<int, int> hitNumbers;

  hitNumbers[lcio::ILDDetID::VXD] = 0;
  hitNumbers[lcio::ILDDetID::SIT] = 0;
  hitNumbers[lcio::ILDDetID::FTD] = 0;
  hitNumbers[lcio::ILDDetID::TPC] = 0;
  hitNumbers[lcio::ILDDetID::SET] = 0;
  hitNumbers[lcio::ILDDetID::ETD] = 0;

  std::vector<TrackerHit*> trackerHits = trackImpl->getTrackerHits();

  streamlog_out(DEBUG0) << "DDCellsAutomatonMV: finaliseTrack - trackerHits.size() = " << trackerHits.size() << "\n";

  for (unsigned j = 0; j < trackerHits.size(); j++) {
    UTIL::BitField64 encoder(LCTrackerCellID::encoding_string());
    encoder.setValue(trackerHits[j]->getCellID0());
    int subdet = encoder[lcio::LCTrackerCellID::subdet()];

    streamlog_out(DEBUG0) << "DDCellsAutomatonMV: finaliseTrack - subdet = " << subdet << "\n";

    ++hitNumbers[subdet];
  }

  trackImpl->subdetectorHitNumbers().resize(2 * lcio::ILDDetID::ETD);
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 2] = hitNumbers[lcio::ILDDetID::VXD];
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::FTD - 2] = hitNumbers[lcio::ILDDetID::FTD];
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::SIT - 2] = hitNumbers[lcio::ILDDetID::SIT];
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::TPC - 2] = hitNumbers[lcio::ILDDetID::TPC];
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::SET - 2] = hitNumbers[lcio::ILDDetID::SET];
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::ETD - 2] = hitNumbers[lcio::ILDDetID::ETD];
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 1] = hitNumbers[lcio::ILDDetID::VXD];
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::FTD - 1] = hitNumbers[lcio::ILDDetID::FTD];
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::SIT - 1] = hitNumbers[lcio::ILDDetID::SIT];
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::TPC - 1] = hitNumbers[lcio::ILDDetID::TPC];
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::SET - 1] = hitNumbers[lcio::ILDDetID::SET];
  trackImpl->subdetectorHitNumbers()[2 * lcio::ILDDetID::ETD - 1] = hitNumbers[lcio::ILDDetID::ETD];

  return;
}

double DDCellsAutomatonMV::Dist(EVENT::TrackerHit* toHit, EVENT::TrackerHit* fromHit) {
  const double* posOuter = fromHit->getPosition();
  const double* posInner = toHit->getPosition();

  double dist = sqrt((posOuter[0] - posInner[0]) * (posOuter[0] - posInner[0]) +
                     (posOuter[1] - posInner[1]) * (posOuter[1] - posInner[1]) +
                     (posOuter[2] - posInner[2]) * (posOuter[2] - posInner[2]));

  // streamlog_out(DEBUG4) << " distance between minivectors = " << dist << std::endl ;

  return dist;
}

bool DDCellsAutomatonMV::thetaAgreementImproved(EVENT::TrackerHit* toHit, EVENT::TrackerHit* fromHit, int layer) {
  // Improved error estimation of polar angle theta, taking into account the uncertainty of the radius as well...
  // also the single point resolution is not hard-coded any more

  bool agreement = false;

  // double resolution = 0.004 ; // just temporarily here
  // double MPS_factor = 0.0 ;  // just for test, and only for vertical muons of 4GeV

  double pos_outer[3];
  double pos_inner[3];

  pos_outer[0] = fromHit->getPosition()[0];
  pos_outer[1] = fromHit->getPosition()[1];
  pos_outer[2] = fromHit->getPosition()[2];
  double radout = sqrt(pos_outer[0] * pos_outer[0] + pos_outer[1] * pos_outer[1]);
  //  double theta_out = (180.0 * atan(radout/pos_outer[2])) / M_PI ;
  // BEFORE BAD RANGE - at 90deg possible flip of sign -- fixed track ineff, still bad track
  double theta_out = atan2(radout, pos_outer[2]);
  theta_out = 2 * M_PI - theta_out;
  // theta_out = theta_out % (2*M_PI);
  // theta_out = fmod(theta_out,2*M_PI);
  while (theta_out >= 2 * M_PI) {
    theta_out = theta_out - 2 * M_PI;
  }
  theta_out = theta_out * 180. / M_PI;

  pos_inner[0] = toHit->getPosition()[0];
  pos_inner[1] = toHit->getPosition()[1];
  pos_inner[2] = toHit->getPosition()[2];
  double radinn = sqrt(pos_inner[0] * pos_inner[0] + pos_inner[1] * pos_inner[1]);
  // double theta_inn = (180.0 *atan(radinn/pos_inner[2])) / M_PI;
  // BEFORE BAD RANGE - at 90deg possible flip of sign -- fixed track ineff, still bad track
  double theta_inn = atan2(radinn, pos_inner[2]);
  theta_inn = 2 * M_PI - theta_inn;
  while (theta_inn >= 2 * M_PI) {
    theta_inn = theta_inn - 2 * M_PI;
  }
  theta_inn = theta_inn * 180. / M_PI;

  double diff_theta = fabs(theta_out - theta_inn);

  streamlog_out(DEBUG4) << "Theta of outer hit " << theta_out << " Theta of inner hit " << theta_inn
                        << " Theta diff between the adjacent hits " << diff_theta << " maximum value " << _hitPairThDiff
                        << std::endl;

  /*
  int celId_inner = toHit->getCellID0() ;
  UTIL::BitField64 encoder( lcio::LCTrackerCellID::encoding_string() ) ;
  encoder.setValue(celId_inner) ;
  int inner_layer  = encoder[lcio::LCTrackerCellID::layer()] ;

  double resolution = _resU[inner_layer] ;

  //double deltaTheta = 3*resolution*sqrt( (radinn/((radinn*radinn) +
  (pos_inner[2]*pos_inner[2])))*(radinn/((radinn*radinn) + (pos_inner[2]*pos_inner[2])))  +
  ((pos_inner[0]*pos_inner[2])/(((radinn*radinn) +
  (pos_inner[2]*pos_inner[2]))*radinn))*((pos_inner[0]*pos_inner[2])/(((radinn*radinn) +
  (pos_inner[2]*pos_inner[2]))*radinn))  +  ((pos_inner[1]*pos_inner[2])/(((radinn*radinn) +
  (pos_inner[2]*pos_inner[2]))*radinn))*((pos_inner[1]*pos_inner[2])/(((radinn*radinn) +
  (pos_inner[2]*pos_inner[2]))*radinn)) ) ;

  double deltazent = radinn/((radinn*radinn) + (pos_inner[2]*pos_inner[2])) ;
  double deltaxi = (pos_inner[0]*pos_inner[2])/(((radinn*radinn) + (pos_inner[2]*pos_inner[2]))*radinn) ;
  double deltapsi = (pos_inner[1]*pos_inner[2])/(((radinn*radinn) + (pos_inner[2]*pos_inner[2]))*radinn) ;

  double deltaTheta = 5*resolution*sqrt((deltaxi*deltaxi) + (deltapsi*deltapsi) + (deltazent*deltazent)) ;

  if ( theta_inn < (theta_out + (deltaTheta + MPS_factor)) &&  theta_inn > (theta_out - (deltaTheta + MPS_factor)) ){
    agreement = true ;
  }

  streamlog_out(DEBUG3) << " outer radius " << radout << " inner radius " << radinn << " resolution " <<
  _resU[inner_layer] << " outer theta " << theta_out*(180.0/M_PI) << " inner theta " << theta_inn*(180.0/M_PI) << "
  deltaTheta "  << deltaTheta*(180.0/M_PI)  << " Agreement " << agreement << " deltaxi = " << deltaxi  << " deltapsi = "
  << deltapsi  << " deltazent = " << deltazent << std::endl ;
  */

  if (layer == 6 || layer == 4) {
    if (diff_theta < _hitPairThDiff) {
      agreement = true;
    }
  }

  if (layer == 2) {
    if (diff_theta < _hitPairThDiffInner) {
      agreement = true;
    }
  }

  return agreement;
}
