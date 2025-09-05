#include "FullLDCTracking_MarlinTrk.h"

#include <EVENT/LCCollection.h>
#include <EVENT/LCObject.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHitPlane.h>
#include <EVENT/TrackerHitZCylinder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackImpl.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCTOOLS.h>
#include <marlin/Global.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>

#include "ClusterShapes.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "MarlinTrk/Factory.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/LCIOTrackPropagators.h"
#include "MarlinTrk/MarlinTrkDiagnostics.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#ifdef MARLINTRK_DIAGNOSTICS_ON
#include "MarlinTrk/DiagnosticsController.h"
#endif

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include <climits>
#include <cmath>

#include "UTIL/LCTrackerConf.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

using namespace lcio;
using namespace marlin;
using namespace MarlinTrk;

FullLDCTracking_MarlinTrk aFullLDCTracking_MarlinTrk;

/** debug printout helper method */
std::string toString(int iTrk, Track* tpcTrack, float bField = 3.5) {
  int nHits = int(tpcTrack->getTrackerHits().size());
  float d0TPC = tpcTrack->getD0();
  float z0TPC = tpcTrack->getZ0();
  float omegaTPC = tpcTrack->getOmega();
  float phi0TPC = tpcTrack->getPhi();
  float tanLTPC = tpcTrack->getTanLambda();
  float Chi2TPC = tpcTrack->getChi2() / float(tpcTrack->getNdf());
  int ndfTPC = tpcTrack->getNdf();

  int nlinkedTracks = tpcTrack->getTracks().size();

  HelixClass helixTPC;

  helixTPC.Initialize_Canonical(phi0TPC, d0TPC, z0TPC, omegaTPC, tanLTPC, bField);

  char strg[200];

  float pxTPC = helixTPC.getMomentum()[0];
  float pyTPC = helixTPC.getMomentum()[1];
  float pzTPC = helixTPC.getMomentum()[2];
  const float ptot = sqrt(pxTPC * pxTPC + pyTPC * pyTPC + pzTPC * pzTPC);

  sprintf(strg, "%3i   %5i %9.3f  %9.3f  %9.3f  %7.2f  %7.2f  %7.2f %4i %4i %8.3f %8i", iTrk, tpcTrack->id(), ptot,
          d0TPC, z0TPC, pxTPC, pyTPC, pzTPC, nHits, ndfTPC, Chi2TPC, nlinkedTracks);

  return std::string(strg);
}

FullLDCTracking_MarlinTrk::FullLDCTracking_MarlinTrk() : Processor("FullLDCTracking_MarlinTrk") {
  _description = "Performs full tracking in ILD detector";

  // Input tracker hit collections

  registerInputCollection(LCIO::TRACKERHITPLANE, "FTDPixelHitCollectionName", "FTD Pixel Hit Collection Name",
                          _FTDPixelHitCollection, std::string("FTDPixelTrackerHits"));

  registerInputCollection(LCIO::TRACKERHIT, "FTDSpacePointCollectionName", "FTD FTDSpacePoint Collection Name",
                          _FTDSpacePointCollection, std::string("FTDSpacePoints"));

  registerInputCollection(LCIO::TRACKERHITPLANE, "VTXHitCollection", "VTX Hit Collection Name",
                          _VTXTrackerHitCollection, std::string("VTXTrackerHits"));

  registerInputCollection(LCIO::TRACKERHIT, "SITHitCollection", "SIT Hit Collection Name", _SITTrackerHitCollection,
                          std::string("SITTrackerHits"));

  registerInputCollection(LCIO::TRACKERHIT, "SETHitCollection", "SET Hit Collection Name", _SETTrackerHitCollection,
                          std::string("SETTrackerHits"));

  registerInputCollection(LCIO::TRACKERHIT, "ETDHitCollection", "ETD Hit Collection Name", _ETDTrackerHitCollection,
                          std::string("ETDTrackerHits"));

  registerInputCollection(LCIO::TRACKERHIT, "TPCHitCollection", "TPC Hit Collection Name", _TPCTrackerHitCollection,
                          std::string("TPCTrackerHits"));

  // Input track collections

  registerInputCollection(LCIO::TRACK, "TPCTracks", "TPC Track Collection", _TPCTrackCollection,
                          std::string("TPCTracks"));

  registerInputCollection(LCIO::TRACK, "SiTracks", "Si Track Collection", _SiTrackCollection, std::string("SiTracks"));

  // Input relation collections

  registerInputCollection(LCIO::LCRELATION, "TPCTracksMCPRelColl", "TPC Track to MCP Relation Collection Name",
                          _TPCTrackMCPCollName, std::string("TPCTracksMCP"));

  registerInputCollection(LCIO::LCRELATION, "SiTracksMCPRelColl", "Si Track to Collection", _SiTrackMCPCollName,
                          std::string("SiTracksMCP"));

  // Output track collection
  registerOutputCollection(LCIO::TRACK, "LDCTrackCollection", "LDC track collection name", _LDCTrackCollection,
                           std::string("LDCTracks"));

  // steering parameters
  registerProcessorParameter("VertexDetectorName", "Name of the vertex (barrel) tracking detector in the geometry",
                             m_vtxDetName, std::string("VXD"));

  registerProcessorParameter("ForwardTrackerDetectorName", "Name of the forward tracking detector in the geometry",
                             m_forwardTrackerDetName, std::string("FTD"));

  registerProcessorParameter("InnerTrackerDetectorName", "Name of the inner (barrel) tracking detector in the geometry",
                             m_innerBarrelTrackerDetName, std::string("SIT"));

  registerProcessorParameter("VertexEndcapDetectorName", "Name of the vertex endcap tracking detector in the geometry",
                             m_vertexEndcapDetName, std::string("VertexEndcap"));

  registerProcessorParameter("SiliconExternalDetectorName",
                             "Name of the silicon external tracking detector in the geometry",
                             m_siliconExternalTrackerDetName, std::string("SET"));

  registerProcessorParameter("D0CutForMerging", "Cut on D0 difference for merging of Si and TPC segments",
                             _d0CutForMerging, float(500.0));

  registerProcessorParameter("Z0CutForMerging", "Cut on Z0 difference for merging of Si and TPC segments",
                             _z0CutForMerging, float(1000.0));

  registerProcessorParameter("OmegaCutForMerging", "Cut on Omega difference for merging Si and TPC segments",
                             _dOmegaForMerging, float(0.25));

  registerProcessorParameter("AngleCutForMerging", "Cut on Opening Angle for merging Si and TPC segments",
                             _angleForMerging, float(0.10));

  registerProcessorParameter("Chi2FitCut", "Cut on fit Chi2", _chi2FitCut, float(100.0));

  registerProcessorParameter("Debug", "Activate debugging?", _debug, int(0));

  registerProcessorParameter("ForceSiTPCMerging", "Force merging of Si and TPC segments?", _forceMerging, int(0));

  registerProcessorParameter("D0CutForForcedMerging", "Cut on D0 difference for forced merging of Si and TPC segments",
                             _d0CutForForcedMerging, float(50.));

  registerProcessorParameter("Z0CutForForcedMerging", "Cut on Z0 difference for forced merging of Si and TPC segments",
                             _z0CutForForcedMerging, float(200.));

  registerProcessorParameter("OmegaCutForForcedMerging",
                             "Cut on Omega difference for forced merging of Si and TPC segments",
                             _dOmegaForForcedMerging, float(0.15));

  registerProcessorParameter("AngleCutForForcedMerging",
                             "Cut on Opening Angle for forced merging of Si and TPC segments", _angleForForcedMerging,
                             float(0.05));

  registerProcessorParameter("ForceTPCSegmentsMerging", "Force merging of TPC Segments?", _mergeTPCSegments, int(0));

  registerProcessorParameter("D0CutToMergeTPCSegments", "Cut on D0 difference for merging TPC segments",
                             _d0CutToMergeTPC, float(100.));

  registerProcessorParameter("Z0CutToMergeTPCSegments", "Cut on Z0 difference for merging TPC segments",
                             _z0CutToMergeTPC, float(5000.0));

  registerProcessorParameter("DeltaPCutToMergeTPCSegments", "Cut on dP/P difference for merging TPC segments",
                             _dPCutToMergeTPC, float(0.1));

  registerProcessorParameter("PtCutToMergeTPCSegments", "Cut on Pt of tracks for merging TPC segments",
                             _PtCutToMergeTPC, float(1.2));

  registerProcessorParameter("cosThetaCutHighPtMerge",
                             "Cut on cos theta between the two momentum vectors when considering "
                             "merger of high Pt tracks",
                             _cosThetaCutHighPtMerge, float(0.99));

  registerProcessorParameter("cosThetaCutSoftHighPtMerge",
                             "cut on cos theta between the two momentum vectors when considering "
                             "merger of high Pt tracks for softer dp/p cut",
                             _cosThetaCutSoftHighPtMerge, float(0.998));

  registerProcessorParameter("momDiffCutHighPtMerge", "cut on dp/p when considering merger of high Pt tracks",
                             _momDiffCutHighPtMerge, float(0.01));

  registerProcessorParameter("momDiffCutSoftHighPtMerge",
                             "softer cut on dp/p when considering merger of high Pt tracks", _momDiffCutSoftHighPtMerge,
                             float(0.25));

  registerProcessorParameter("hitDistanceCutHighPtMerge",
                             "cut on 3D distance between hit and helix extrapolation when considering "
                             "merger of high Pt tracks",
                             _hitDistanceCutHighPtMerge, float(25.0));

  registerProcessorParameter("maxHitDistanceCutHighPtMerge",
                             "cut for max 3D distance between any hit and helix extrapolation when "
                             "considering merger of high Pt tracks",
                             _maxHitDistanceCutHighPtMerge, float(50.0));

  registerProcessorParameter("maxFractionOfOutliersCutHighPtMerge",
                             "cut on maximum fraction of outliers when "
                             "considering merger of high Pt tracks",
                             _maxFractionOfOutliersCutHighPtMerge, float(0.95));

  registerProcessorParameter("CutOnTPCHits", "Cut on the number of the TPC hits for tracks with no Si hits",
                             _cutOnTPCHits, int(35));

  registerProcessorParameter("CutOnSiHits", "Cut on the number of the Si hits for tracks with no TPC hits",
                             _cutOnSiHits, int(4));

  registerProcessorParameter("AssignVTXHits", "Assign left over VTX hits", _assignVTXHits, int(1));

  registerProcessorParameter("AssignFTDHits", "Assign left over FTD hits", _assignFTDHits, int(1));

  registerProcessorParameter("AssignSITHits", "Assign left over SIT hits", _assignSITHits, int(1));

  registerProcessorParameter("AssignTPCHits", "Assign left over TPC hits", _assignTPCHits, int(1));

  registerProcessorParameter("AssignSETHits", "Assign SET Hits", _assignSETHits, int(1));

  registerProcessorParameter("AssignETDHits", "Assign ETD Hits", _assignETDHits, int(1));

  registerProcessorParameter("NHitsExtrapolation", "number of hits for outward extrapolation", _nHitsExtrapolation,
                             int(35));

  registerProcessorParameter("VTXHitToTrackDistance", "Cut on distance between track and VTX hits", _distCutForVTXHits,
                             float(1.5));

  registerProcessorParameter("FTDHitToTrackDistance", "Cut on distance between track and FTD hits", _distCutForFTDHits,
                             float(2.0));

  registerProcessorParameter("SITHitToTrackDistance", "Cut on distance between track and SIT hits", _distCutForSITHits,
                             float(2.0));

  registerProcessorParameter("SETHitToTrackDistance", "Cut on distance between track and SET hits", _distCutForSETHits,
                             float(2.0));

  registerProcessorParameter("ETDHitToTrackDistance", "Cut on distance between track and ETD hits", _distCutForETDHits,
                             float(10.0));

  registerProcessorParameter("TPCHitToTrackDistance", "Cut on distance between track and TPC hits", _distCutForTPCHits,
                             float(15.0));

  registerProcessorParameter("CutOnTrackD0", "Cut on the track parameter D0", _d0TrkCut, float(500.));

  registerProcessorParameter("CutOnTrackZ0", "Cut on the track parameter Z0", _z0TrkCut, float(500.));

  registerProcessorParameter("ForbidOverlapInZTPC", "Forbid overlap in Z for the merged TPC segments",
                             _forbidOverlapInZTPC, int(0));

  registerProcessorParameter("ForbidOverlapInZComb",
                             "Forbid overlap in Z for combining TPC segments "
                             "with tracks having Si hits",
                             _forbidOverlapInZComb, int(0));

  registerProcessorParameter("MultipleScatteringOn", "Use MultipleScattering in Fit", _MSOn, bool(true));

  registerProcessorParameter("EnergyLossOn", "Use Energy Loss in Fit", _ElossOn, bool(true));

  registerProcessorParameter("SmoothOn", "Smooth All Mesurement Sites in Fit", _SmoothOn, bool(true));

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

  registerProcessorParameter("MinChi2ProbForSiliconTracks", "Minimum Chi-squared P value allowed for Silicon Tracks.",
                             _minChi2ProbForSiliconTracks, double(1.e-03));

  registerProcessorParameter("VetoMergeMomentumCut", "Minimum momentum for which Veto is applicable",
                             _vetoMergeMomentumCut, float(2.5));

  registerProcessorParameter("MaxAllowedPercentageOfOutliersForTrackCombination",
                             "Maximum number of outliers allowed before track combination is vetoed",
                             _maxAllowedPercentageOfOutliersForTrackCombination, float(0.3));

  registerProcessorParameter("MaxAllowedSiHitRejectionsForTrackCombination",
                             "Maximum number of outliers allowed before track combination is vetoed",
                             _maxAllowedSiHitRejectionsForTrackCombination, int(2));

  registerProcessorParameter("EnergyLossErrorTPC", "Fractional error of dEdx in the TPC", _energyLossErrorTPC,
                             float(0.05));

  registerProcessorParameter("TrackSystemName",
                             "Name of the track fitting system to be used "
                             "(KalTest, DDKalTest, aidaTT, ... )",
                             _trkSystemName, std::string("KalTest"));

#ifdef MARLINTRK_DIAGNOSTICS_ON

  registerOptionalParameter("RunMarlinTrkDiagnostics",
                            "Run MarlinTrk Diagnostics. MarlinTrk must be "
                            "compiled with MARLINTRK_DIAGNOSTICS_ON defined",
                            _runMarlinTrkDiagnostics, bool(false));

  registerOptionalParameter("DiagnosticsName", "Name of the root file and root tree if running Diagnostics",
                            _MarlinTrkDiagnosticsName, std::string("FullLDCTrackingDiagnostics"));

#endif
}

void FullLDCTracking_MarlinTrk::init() {
  printParameters();
  _nRun = -1;
  _nEvt = 0;
  PI = acos(-1.);
  PIOVER2 = 0.5 * PI;
  TWOPI = 2 * PI;

  _encoder = new UTIL::BitField64(lcio::LCTrackerCellID::encoding_string());

  // set upt the geometry
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

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  this->setupGeom(theDetector);
}

void FullLDCTracking_MarlinTrk::processRunHeader(LCRunHeader* run) {
  _nRun++;
  _nEvt = 0;
  streamlog_out(DEBUG5) << std::endl;
  streamlog_out(DEBUG5) << "FullLDCTracking_MarlinTrk ---> new run : run number = " << run->getRunNumber() << std::endl;
}

void FullLDCTracking_MarlinTrk::processEvent(LCEvent* evt) {
  // set the correct configuration for the tracking system for this event
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(_trksystem, _MSOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(_trksystem, _ElossOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);

  _evt = evt;

  streamlog_out(DEBUG5) << std::endl;
  streamlog_out(DEBUG5) << "FullLDCTracking_MarlinTrk -> run = " << evt->getRunNumber()
                        << "  event = " << evt->getEventNumber() << std::endl;
  streamlog_out(DEBUG5) << std::endl;

  prepareVectors(evt);
  streamlog_out(DEBUG5) << "************************************PrepareVectors done..." << std::endl;

  streamlog_out(DEBUG5) << "************************************Merge TPC/Si ..." << std::endl;

  MergeTPCandSiTracks();
  streamlog_out(DEBUG5) << "************************************Merging done ..." << std::endl;

  MergeTPCandSiTracksII();
  streamlog_out(DEBUG5) << "************************************Merging II done ..." << std::endl;

  Sorting(_allCombinedTracks);
  streamlog_out(DEBUG5) << "************************************Sorting by Chi2/NDF done ..." << std::endl;

  streamlog_out(DEBUG5) << "************************************Selection of "
                           "all 2 track combininations ..."
                        << std::endl;
  SelectCombinedTracks();
  streamlog_out(DEBUG5) << "************************************Selection of "
                           "all 2 track combininations done ..."
                        << std::endl;

  streamlog_out(DEBUG5) << "************************************Trying non combined tracks ..." << std::endl;
  AddNotCombinedTracks();
  streamlog_out(DEBUG5) << "************************************Non combined tracks added ..." << std::endl;
  // CheckTracks( );

  streamlog_out(DEBUG5) << "************************************Add Non assigned hits ..." << std::endl;
  AddNotAssignedHits();
  streamlog_out(DEBUG5) << "************************************Non assigned hits added ..." << std::endl;

  AddTrackColToEvt(evt, _trkImplVec, _LDCTrackCollection);
  streamlog_out(DEBUG5) << "Collections added to event ..." << std::endl;
  CleanUp();
  streamlog_out(DEBUG5) << "Cleanup is done." << std::endl;
  _nEvt++;
  //  getchar();
  streamlog_out(DEBUG5) << std::endl;
  streamlog_out(DEBUG5) << std::endl;
}

void FullLDCTracking_MarlinTrk::AddTrackColToEvt(LCEvent* evt, TrackExtendedVec& trkVec, std::string TrkColName) {
  LCCollectionVec* colTRK = new LCCollectionVec(LCIO::TRACK);
  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  colTRK->setFlag(trkFlag.getFlag());

  streamlog_out(DEBUG5) << "AddTrackColToEvt: Collection " << TrkColName << " is being added to event " << std::endl;

  //  LCCollectionVec * colRel = NULL;

  int nTrkCand = int(trkVec.size());

  int nTotTracks = 0;
  float eTot = 0.0;
  float pxTot = 0.0;
  float pyTot = 0.0;
  float pzTot = 0.0;

  // SJA:FIXME: So here we are going to do one final refit. This can certainly
  // be optimised, but rather than worry about the mememory management right now
  // lets make it work, and optimise it later ...

  for (int iTRK = 0; iTRK < nTrkCand; ++iTRK) {
    TrackExtended* trkCand = trkVec[iTRK];
    TrackerHitExtendedVec& hitVec = trkCand->getTrackerHitExtendedVec();
    EVENT::TrackerHitVec trkHits;

    streamlog_out(DEBUG2) << " Trying to add track " << trkCand << " to final lcio collection " << std::endl;

    int nHits = int(hitVec.size());

    streamlog_out(DEBUG2) << " Trying to add track " << trkCand << " to final lcio collection nHits = " << nHits
                          << std::endl;

    for (int ihit = 0; ihit < nHits; ++ihit) {
      if (hitVec[ihit]->getUsedInFit() == false) {
        streamlog_out(DEBUG2) << "rejecting hit for track " << trkCand << " at zhit  "
                              << hitVec[ihit]->getTrackerHit()->getPosition()[2] << std::endl;
        continue;
      }

      EVENT::TrackerHit* trkHit = hitVec[ihit]->getTrackerHit();

      if (trkHit) {
        trkHits.push_back(trkHit);
      } else {
        throw EVENT::Exception(std::string("FullLDCTracking_MarlinTrk::AddTrackColToEvt: "
                                           "TrackerHit pointer == NULL "));
      }
    }

    if (trkHits.size() < 3) {
      streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AddTrackColToEvt: "
                               "Cannot fit less than 3 hits. Number of hits =  "
                            << trkHits.size() << std::endl;
      continue;
    }

    TrackImpl* Track = new TrackImpl;

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

    // get the track state at the last hit at the outer most hit

    GroupTracks* group = trkCand->getGroupTracks();

    TrackStateImpl ts_initial;

    bool prefit_set = false;

    streamlog_out(DEBUG2) << "Track Group = " << group << std::endl;

    if (group)
      streamlog_out(DEBUG2) << "Track Group size = " << group->getTrackExtendedVec().size() << std::endl;

    if (group != NULL && group->getTrackExtendedVec().size() > 0) {
      // get the second track as this must be the one furthest from the IP

      TrackExtended* te = 0;

      if (group->getTrackExtendedVec().size() == 1) {
        te = group->getTrackExtendedVec()[0];
      } else {
        te = group->getTrackExtendedVec()[1];
      }

      if (te->getTrack()->getTrackState(lcio::TrackState::AtLastHit)) {
        streamlog_out(DEBUG2) << "Initialise Fit with trackstate from last hit" << group << std::endl;

        ts_initial = *(te->getTrack()->getTrackState(lcio::TrackState::AtLastHit));

        prefit_set = true;
      }
    }

    if (!prefit_set) { // use parameters at IP

      streamlog_out(DEBUG2) << "Initialise Fit with trackstate from IP " << group << std::endl;

      ts_initial.setD0(trkCand->getD0());
      ts_initial.setPhi(trkCand->getPhi());
      ts_initial.setZ0(trkCand->getZ0());
      ts_initial.setOmega(trkCand->getOmega());
      ts_initial.setTanLambda(trkCand->getTanLambda());

      float ref[3];
      ref[0] = ref[1] = ref[2] = 0.0;

      ts_initial.setReferencePoint(ref);

      ts_initial.setLocation(lcio::TrackStateImpl::AtIP);
    }

    ts_initial.setCovMatrix(covMatrix);

    // sort hits in R
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

    for (std::vector<std::pair<float, EVENT::TrackerHit*>>::iterator it = r2_values.begin(); it != r2_values.end();
         ++it) {
      trkHits.push_back(it->second);
    }

    bool fit_backwards = IMarlinTrack::backward;

    MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();

    int error = 0;

    try {
      error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, Track, fit_backwards, &ts_initial, _bField,
                                                  _maxChi2PerHit);
    } catch (...) {
      //      delete Track;
      //      delete marlinTrk;

      throw;
    }

#ifdef MARLINTRK_DIAGNOSTICS_ON
    if (error != IMarlinTrack::success && _runMarlinTrkDiagnostics) {
      void* dcv = _trksystem->getDiagnositicsPointer();
      DiagnosticsController* dc = static_cast<DiagnosticsController*>(dcv);
      dc->skip_current_track();
    }
#endif

    std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
    std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;
    std::vector<TrackerHit*> all_hits;
    all_hits.reserve(300);

    marlinTrk->getHitsInFit(hits_in_fit);

    for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
      all_hits.push_back(hits_in_fit[ihit].first);
    }

    UTIL::BitField64 cellID_encoder(lcio::LCTrackerCellID::encoding_string());

    MarlinTrk::addHitNumbersToTrack(Track, all_hits, true, cellID_encoder);

    marlinTrk->getOutliers(outliers);

    for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      all_hits.push_back(outliers[ihit].first);
    }

    MarlinTrk::addHitNumbersToTrack(Track, all_hits, false, cellID_encoder);

    delete marlinTrk;

    if (error != IMarlinTrack::success) {
      streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AddTrackColToEvt: Track fit failed with error code " << error
                            << " track dropped. Number of hits = " << trkHits.size() << std::endl;

      delete Track;
      continue;
    }

    if (Track->getNdf() < 0) {
      streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AddTrackColToEvt: Track fit returns " << Track->getNdf()
                            << " degress of freedom track dropped. Number of hits = " << trkHits.size() << std::endl;

      delete Track;
      continue;
    }

    const TrackState* trkStateIP = Track->getTrackState(lcio::TrackState::AtIP);

    if (trkStateIP == 0) {
      streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AddTrackColToEvt: Track fit returns " << Track->getNdf()
                            << " degress of freedom track dropped. Number of hits = " << trkHits.size() << std::endl;
      throw EVENT::Exception(std::string("FullLDCTracking_MarlinTrk::AddTrackColToEvt: trkStateIP pointer == NULL "));
    }

    if (group != NULL) {
      TrackExtendedVec trkVecGrp = group->getTrackExtendedVec();
      int nGrTRK = int(trkVecGrp.size());
      for (int iGr = 0; iGr < nGrTRK; ++iGr) {
        TrackExtended* subTrack = trkVecGrp[iGr];
        Track->addTrack(subTrack->getTrack());

        // check if it is a tpc looper ...
        if (BitSet32(subTrack->getTrack()->getType())[lcio::ILDDetID::TPC]) {
          const TrackVec segments = subTrack->getTrack()->getTracks();

          if (segments.empty() == false) {
            for (unsigned iSeg = 0; iSeg < segments.size(); ++iSeg) {
              Track->addTrack(segments[iSeg]);
            }
          }
        }
      }
    }

    float d0TrkCand = trkCand->getD0();
    float z0TrkCand = trkCand->getZ0();
    //    float phi0TrkCand = trkCand->getPhi();

    int nhits_in_vxd = Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 2];
    int nhits_in_ftd = Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::FTD - 2];
    int nhits_in_sit = Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::SIT - 2];
    int nhits_in_tpc = Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::TPC - 2];
    int nhits_in_set = Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::SET - 2];

    int nHitsSi = nhits_in_vxd + nhits_in_ftd + nhits_in_sit;

    streamlog_out(DEBUG3) << " Hit numbers for Track " << Track->id() << ": "
                          << " vxd hits = " << nhits_in_vxd << " ftd hits = " << nhits_in_ftd
                          << " sit hits = " << nhits_in_sit << " tpc hits = " << nhits_in_tpc
                          << " set hits = " << nhits_in_set << std::endl;

    if (nhits_in_vxd > 0)
      Track->setTypeBit(lcio::ILDDetID::VXD);
    if (nhits_in_ftd > 0)
      Track->setTypeBit(lcio::ILDDetID::FTD);
    if (nhits_in_sit > 0)
      Track->setTypeBit(lcio::ILDDetID::SIT);
    if (nhits_in_tpc > 0)
      Track->setTypeBit(lcio::ILDDetID::TPC);
    if (nhits_in_set > 0)
      Track->setTypeBit(lcio::ILDDetID::SET);

    bool rejectTrack_onTPCHits = (nhits_in_tpc < _cutOnTPCHits) && (nHitsSi <= 0);

    bool rejectTrackonSiliconHits = ((nhits_in_tpc <= 0) && (nHitsSi < _cutOnSiHits));
    bool rejectTrackonImpactParameters = (fabs(d0TrkCand) > _d0TrkCut) || (fabs(z0TrkCand) > _z0TrkCut);

    if (rejectTrack_onTPCHits || rejectTrackonSiliconHits || rejectTrackonImpactParameters) {
      streamlog_out(DEBUG3) << " Track " << trkCand << " rejected : rejectTrack_onTPCHits = " << rejectTrack_onTPCHits
                            << " rejectTrackonSiliconHits " << rejectTrackonSiliconHits
                            << " rejectTrackonImpactParameters " << rejectTrackonImpactParameters << std::endl;

      delete Track;
    } else {
      float omega = trkStateIP->getOmega();
      float tanLambda = trkStateIP->getTanLambda();
      float phi0 = trkStateIP->getPhi();
      float d0 = trkStateIP->getD0();
      float z0 = trkStateIP->getZ0();

      HelixClass helix;
      helix.Initialize_Canonical(phi0, d0, z0, omega, tanLambda, _bField);

      float trkPx = helix.getMomentum()[0];
      float trkPy = helix.getMomentum()[1];
      float trkPz = helix.getMomentum()[2];
      float trkP = sqrt(trkPx * trkPx + trkPy * trkPy + trkPz * trkPz);

      eTot += trkP;
      pxTot += trkPx;
      pyTot += trkPy;
      pzTot += trkPz;
      nTotTracks++;

      streamlog_out(DEBUG3) << " Add Track to final Collection: ID = " << Track->id() << " for trkCand " << trkCand
                            << std::endl;
      colTRK->addElement(Track);
    }
  }

  streamlog_out(DEBUG5) << std::endl;
  streamlog_out(DEBUG5) << "Number of accepted " << TrkColName << " = " << nTotTracks << std::endl;
  streamlog_out(DEBUG5) << "Total 4-momentum of " << TrkColName << " : E = " << eTot << " Px = " << pxTot
                        << " Py = " << pyTot << " Pz = " << pzTot << std::endl;
  streamlog_out(DEBUG5) << std::endl;

  evt->addCollection(colTRK, TrkColName.c_str());
}

void FullLDCTracking_MarlinTrk::prepareVectors(LCEvent* event) {
  _allTPCHits.clear();
  _allVTXHits.clear();
  _allFTDHits.clear();
  _allSITHits.clear();
  _allSETHits.clear();
  _allETDHits.clear();
  _allTPCTracks.clear();
  _allSiTracks.clear();
  _allCombinedTracks.clear();
  _allNonCombinedTPCTracks.clear();
  _allNonCombinedSiTracks.clear();
  _trkImplVec.clear();
  _candidateCombinedTracks.clear();

  std::map<TrackerHit*, TrackerHitExtended*> mapTrackerHits;

  // Reading TPC hits
  try {
    LCCollection* col = event->getCollection(_TPCTrackerHitCollection.c_str());

    int nelem = col->getNumberOfElements();

    streamlog_out(DEBUG5) << "Number of TPC hits = " << nelem << std::endl;

    for (int ielem = 0; ielem < nelem; ++ielem) {
      TrackerHit* hit = dynamic_cast<TrackerHit*>(col->getElementAt(ielem));
      TrackerHitExtended* hitExt = new TrackerHitExtended(hit);

      // Covariance Matrix in LCIO is defined in XYZ convert to R-Phi-Z
      // For no error in r

      double tpcRPhiRes = sqrt(hit->getCovMatrix()[0] + hit->getCovMatrix()[2]);
      double tpcZRes = sqrt(hit->getCovMatrix()[5]);

      hitExt->setResolutionRPhi(float(tpcRPhiRes));
      hitExt->setResolutionZ(float(tpcZRes));

      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      hitExt->setDet(int(INT_MAX));
      _allTPCHits.push_back(hitExt);
      mapTrackerHits[hit] = hitExt;

      double pos[3];

      for (int i = 0; i < 3; ++i) {
        pos[i] = hit->getPosition()[i];
      }

      unsigned int layer = static_cast<unsigned int>(getLayerID(hit));

      streamlog_out(DEBUG1) << " TPC Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi "
                            << tpcRPhiRes << " dz " << tpcZRes << "  layer = " << layer << std::endl;
    }
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << _TPCTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  };

  // Reading in FTD Pixel Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    LCCollection* hitCollection = event->getCollection(_FTDPixelHitCollection.c_str());

    int nelem = hitCollection->getNumberOfElements();

    streamlog_out(DEBUG5) << "Number of FTD Pixel hits = " << nelem << std::endl;

    for (int ielem = 0; ielem < nelem; ++ielem) {
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem));

      TrackerHitExtended* hitExt = new TrackerHitExtended(hit);

      double point_res_rphi = sqrt(hit->getdU() * hit->getdU() + hit->getdV() * hit->getdV());
      hitExt->setResolutionRPhi(point_res_rphi);

      // SJA:FIXME why is this needed?
      hitExt->setResolutionZ(0.1);

      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      hitExt->setDet(int(INT_MAX));

      _allFTDHits.push_back(hitExt);
      mapTrackerHits[hit] = hitExt;

      // get the layer number
      unsigned int layer = static_cast<unsigned int>(getLayerID(hit));
      unsigned int petalIndex = static_cast<unsigned int>(getModuleID(hit));

      if (_petalBasedFTDWithOverlaps == true) {
        // as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
        // the layers will follow the even odd numbering of the petals
        if (petalIndex % 2 == 0) {
          layer = 2 * layer;
        } else {
          layer = 2 * layer + 1;
        }
      }

      if (layer >= _nLayersFTD) {
        streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk => fatal error in FTD : layer is outside allowed range : "
                             << layer << " number of layers = " << _nLayersFTD << std::endl;
        exit(1);
      }

      double pos[3];

      for (int i = 0; i < 3; ++i) {
        pos[i] = hit->getPosition()[i];
      }

      streamlog_out(DEBUG1) << " FTD Pixel Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi "
                            << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ()
                            << "  layer = " << layer << std::endl;
    }
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << _FTDPixelHitCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading in FTD SpacePoint Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    LCCollection* hitCollection = event->getCollection(_FTDSpacePointCollection.c_str());

    int nelem = hitCollection->getNumberOfElements();

    streamlog_out(DEBUG5) << "Number of FTD SpacePoints hits = " << nelem << std::endl;

    for (int ielem = 0; ielem < nelem; ++ielem) {
      TrackerHit* hit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));

      TrackerHitExtended* hitExt = new TrackerHitExtended(hit);

      // SJA:FIXME: fudge for now by a factor of two and ignore covariance
      double point_res_rphi = 2 * sqrt(hit->getCovMatrix()[0] + hit->getCovMatrix()[2]);

      hitExt->setResolutionRPhi(point_res_rphi);

      // SJA:FIXME why is this needed?
      hitExt->setResolutionZ(0.1);

      // type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and
      // catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));

      _allFTDHits.push_back(hitExt);
      mapTrackerHits[hit] = hitExt;

      // get the layer number
      unsigned int layer = static_cast<unsigned int>(getLayerID(hit));
      unsigned int petalIndex = static_cast<unsigned int>(getModuleID(hit));

      if (_petalBasedFTDWithOverlaps == true) {
        // as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
        // the layers will follow the even odd numbering of the petals
        if (petalIndex % 2 == 0) {
          layer = 2 * layer;
        } else {
          layer = 2 * layer + 1;
        }
      }

      if (layer >= _nLayersFTD) {
        streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk => fatal error in FTD : layer is outside allowed range : "
                             << layer << " number of layers = " << _nLayersFTD << std::endl;
        exit(1);
      }

      double pos[3];

      for (int i = 0; i < 3; ++i) {
        pos[i] = hit->getPosition()[i];
      }

      streamlog_out(DEBUG1) << " FTD SpacePoint Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi "
                            << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ()
                            << "  layer = " << layer << std::endl;
    }
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << _FTDSpacePointCollection.c_str() << " collection is unavailable" << std::endl;
  }

  try {
    LCCollection* hitCollection = event->getCollection(_SITTrackerHitCollection.c_str());

    int nelem = hitCollection->getNumberOfElements();

    streamlog_out(DEBUG5) << "Number of SIT hits = " << nelem << std::endl;

    TrackerHit* trkhit = 0;
    TrackerHitPlane* trkhit_P = 0;
    TrackerHitZCylinder* trkhit_C = 0;

    double drphi(NAN);
    double dz(NAN);

    for (int ielem = 0; ielem < nelem; ++ielem) {
      // hit could be of the following type
      // 1) TrackerHit, either ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT or just standard TrackerHit
      // 2) TrackerHitPlane, either 1D or 2D
      // 3) TrackerHitZCylinder, if coming from a simple cylinder design as in the LOI

      // Establish which of these it is in the following order of likelyhood
      //    i)   ILDTrkHitTypeBit::ONE_DIMENSIONAL (TrackerHitPlane) Should Never Happen: SpacePoints Must be Used
      //    Instead ii)  ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT (TrackerHit) iii) TrackerHitPlane (Two dimentional) iv)
      //    TrackerHitZCylinder v)   Must be standard TrackerHit

      trkhit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));

      int layer = getLayerID(trkhit);

      if (layer < 0 || (unsigned)layer >= _nLayersSIT) {
        streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk => fatal error in SIT : hit layer number " << layer
                             << " is outside of the total number of layers in the given geometry: " << _nLayersSIT
                             << std::endl;
        exit(1);
      }

      // first check that we have not been given 1D hits by mistake, as they won't work here
      if (BitSet32(trkhit->getType())[UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL]) {
        streamlog_out(ERROR)
            << "FullLDCTracking_MarlinTrk: SIT Hit cannot be of type UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL COMPOSITE "
               "SPACEPOINTS must be use instead. \n\n  exit(1) called from file "
            << __FILE__ << " and line " << __LINE__ << std::endl;
        exit(1);
      }
      // most likely case: COMPOSITE_SPACEPOINT hits formed from stereo strip hits
      else if (BitSet32(trkhit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) {
        // SJA:FIXME: fudge for now by a factor of two and ignore covariance
        drphi = 2 * sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]);
        dz = sqrt(trkhit->getCovMatrix()[5]);
      }
      // or a PIXEL based SIT, using 2D TrackerHitPlane like the VXD above
      else if ((trkhit_P = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem)))) {
        // first we need to check if the measurement vectors are aligned with the global coordinates
        Vector3D U(1.0, trkhit_P->getU()[1], trkhit_P->getU()[0], Vector3D::spherical);
        Vector3D V(1.0, trkhit_P->getV()[1], trkhit_P->getV()[0], Vector3D::spherical);
        Vector3D Z(0.0, 0.0, 1.0);

        const float eps = 1.0e-07;
        // V must be the global z axis
        if (!(fabs(1.0 - V.dot(Z)) > eps || fabs(-1.0 - V.dot(Z)) > eps)) {
          streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk: PIXEL SIT Hit measurement vectors V is not equal to the "
                                  "global Z axis. \n\n  exit(1) called from file "
                               << __FILE__ << " and line " << __LINE__ << std::endl;
          exit(1);
        }

        // U must be normal to the global z axis
        if (fabs(U.dot(Z)) > eps) {
          streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk: PIXEL SIT Hit measurment vectors U is not in the global "
                                  "X-Y plane. \n\n exit(1) called from file "
                               << __FILE__ << " and line " << __LINE__ << std::endl;
          exit(1);
        }

        drphi = trkhit_P->getdU();
        dz = trkhit_P->getdV();
      }
      // or a simple cylindrical design, as used in the LOI
      else if ((trkhit_C = dynamic_cast<TrackerHitZCylinder*>(hitCollection->getElementAt(ielem)))) {
        drphi = trkhit_C->getdRPhi();
        dz = trkhit_C->getdZ();
      }
      // this would be very unlikely, but who knows ... just an ordinary TrackerHit, which is not a COMPOSITE_SPACEPOINT
      else {
        // SJA:FIXME: fudge for now by a factor of two and ignore covariance
        drphi = 2 * sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]);
        dz = sqrt(trkhit->getCovMatrix()[5]);
      }

      // now that the hit type has been established carry on and create a

      TrackerHitExtended* hitExt = new TrackerHitExtended(trkhit);

      // SJA:FIXME: just use planar res for now
      hitExt->setResolutionRPhi(drphi);
      hitExt->setResolutionZ(dz);

      // set type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try
      // and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));

      _allSITHits.push_back(hitExt);
      mapTrackerHits[trkhit] = hitExt;

      double pos[3];

      for (int i = 0; i < 3; ++i) {
        pos[i] = trkhit->getPosition()[i];
      }

      streamlog_out(DEBUG1) << " SIT Hit " << trkhit->id() << " added : @ " << pos[0] << " " << pos[1] << " " << pos[2]
                            << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ()
                            << "  layer = " << layer << std::endl;
    }
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << " collection not found : " << _SITTrackerHitCollection.c_str() << std::endl;
  }

  try {
    LCCollection* hitCollection = event->getCollection(_SETTrackerHitCollection.c_str());

    int nelem = hitCollection->getNumberOfElements();

    streamlog_out(DEBUG5) << "Number of SET hits = " << nelem << std::endl;

    TrackerHit* trkhit = 0;
    TrackerHitPlane* trkhit_P = 0;
    TrackerHitZCylinder* trkhit_C = 0;

    double drphi(NAN);
    double dz(NAN);

    for (int ielem = 0; ielem < nelem; ++ielem) {
      // hit could be of the following type
      // 1) TrackerHit, either ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT or just standard TrackerHit
      // 2) TrackerHitPlane, either 1D or 2D
      // 3) TrackerHitZCylinder, if coming from a simple cylinder design as in the LOI

      // Establish which of these it is in the following order of likelyhood
      //    i)   ILDTrkHitTypeBit::ONE_DIMENSIONAL (TrackerHitPlane) Should Never Happen: SpacePoints Must be Used
      //    Instead ii)  ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT (TrackerHit) iii) TrackerHitPlane (Two dimentional) iv)
      //    TrackerHitZCylinder v)   Must be standard TrackerHit

      trkhit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));

      int layer = getLayerID(trkhit);

      if (layer < 0 || (unsigned)layer >= _nLayersSET) {
        streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk => fatal error in SET : layer is outside allowed range : "
                             << layer << std::endl;
        exit(1);
      }

      // first check that we have not been given 1D hits by mistake, as they won't work here
      if (BitSet32(trkhit->getType())[UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL]) {
        streamlog_out(ERROR)
            << "FullLDCTracking_MarlinTrk: SET Hit cannot be of type UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL COMPOSITE "
               "SPACEPOINTS must be use instead. \n\n  exit(1) called from file "
            << __FILE__ << " and line " << __LINE__ << std::endl;
        exit(1);
      }
      // most likely case: COMPOSITE_SPACEPOINT hits formed from stereo strip hits
      else if (BitSet32(trkhit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) {
        // SJA:FIXME: fudge for now by a factor of two and ignore covariance
        drphi = 2 * sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]);
        dz = sqrt(trkhit->getCovMatrix()[5]);
      }
      // or a PIXEL based SET, using 2D TrackerHitPlane like the VXD above
      else if ((trkhit_P = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem)))) {
        // first we need to check if the measurement vectors are aligned with the global coordinates
        Vector3D U(1.0, trkhit_P->getU()[1], trkhit_P->getU()[0], Vector3D::spherical);
        Vector3D V(1.0, trkhit_P->getV()[1], trkhit_P->getV()[0], Vector3D::spherical);
        Vector3D Z(0.0, 0.0, 1.0);

        const float eps = 1.0e-07;
        // V must be the global z axis
        if (fabs(1.0 - V.dot(Z)) > eps) {
          streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk: PIXEL SET Hit measurment vectors V is not equal to the "
                                  "global Z axis. \n\n  exit(1) called from file "
                               << __FILE__ << " and line " << __LINE__ << std::endl;
          exit(1);
        }

        // U must be normal to the global z axis
        if (fabs(U.dot(Z)) > eps) {
          streamlog_out(ERROR) << "FullLDCTracking_MarlinTrk: PIXEL SET Hit measurment vectors U is not in the global "
                                  "X-Y plane. \n\n exit(1) called from file "
                               << __FILE__ << " and line " << __LINE__ << std::endl;
          exit(1);
        }

        drphi = trkhit_P->getdU();
        dz = trkhit_P->getdV();
      }
      // or a simple cylindrical design, as used in the LOI
      else if ((trkhit_C = dynamic_cast<TrackerHitZCylinder*>(hitCollection->getElementAt(ielem)))) {
        drphi = trkhit_C->getdRPhi();
        dz = trkhit_C->getdZ();
      }
      // this would be very unlikely, but who knows ... just an ordinary TrackerHit, which is not a COMPOSITE_SPACEPOINT
      else {
        // SJA:FIXME: fudge for now by a factor of two and ignore covariance
        drphi = 2 * sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]);
        dz = sqrt(trkhit->getCovMatrix()[5]);
      }

      // now that the hit type has been established carry on and create a

      TrackerHitExtended* hitExt = new TrackerHitExtended(trkhit);

      // SJA:FIXME: just use planar res for now
      hitExt->setResolutionRPhi(drphi);
      hitExt->setResolutionZ(dz);

      // set type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try
      // and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));

      _allSETHits.push_back(hitExt);
      mapTrackerHits[trkhit] = hitExt;

      double pos[3];

      for (int i = 0; i < 3; ++i) {
        pos[i] = trkhit->getPosition()[i];
      }

      streamlog_out(DEBUG1) << " SET Hit " << trkhit->id() << " added : @ " << pos[0] << " " << pos[1] << " " << pos[2]
                            << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ()
                            << "  layer = " << layer << std::endl;
    }
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << " collection not found : " << _SETTrackerHitCollection.c_str() << std::endl;
  }

  // Reading VTX Hits
  try {
    LCCollection* col = event->getCollection(_VTXTrackerHitCollection.c_str());

    int nelem = col->getNumberOfElements();

    streamlog_out(DEBUG5) << "Number of VXD hits = " << nelem << std::endl;

    for (int ielem = 0; ielem < nelem; ++ielem) {
      TrackerHitPlane* trkhit = dynamic_cast<TrackerHitPlane*>(col->getElementAt(ielem));
      TrackerHitExtended* hitExt = new TrackerHitExtended(trkhit);

      // SJA:FIXME: just use planar res for now
      hitExt->setResolutionRPhi(trkhit->getdU());
      hitExt->setResolutionZ(trkhit->getdV());

      // type and det are no longer used, set to INT_MAX to try and catch any missuse
      hitExt->setType(int(INT_MAX));
      hitExt->setDet(int(INT_MAX));
      _allVTXHits.push_back(hitExt);
      mapTrackerHits[trkhit] = hitExt;

      double pos[3];

      for (int i = 0; i < 3; ++i) {
        pos[i] = trkhit->getPosition()[i];
      }

      int layer = getLayerID(trkhit);

      streamlog_out(DEBUG1) << " VXD Hit " << trkhit->id() << " added : @ " << pos[0] << " " << pos[1] << " " << pos[2]
                            << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ()
                            << "  layer = " << layer << std::endl;
    }
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << _VTXTrackerHitCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading TPC Tracks
  try {
    LCCollection* col = event->getCollection(_TPCTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    streamlog_out(DEBUG5) << std::endl;
    streamlog_out(DEBUG5) << "Number of TPC Tracks = " << nelem << " in " << _TPCTrackCollection.c_str() << std::endl;
    streamlog_out(DEBUG5)
        << " Trk    ID        p          D0         Z0       Px       Py       Pz   ntpc ndf  Chi2/ndf nlinkedTracks"
        << std::endl;
    //                       "  0             1.111   0.059      0.022    -0.54     0.61    -0.45    0.185

    for (int iTrk = 0; iTrk < nelem; ++iTrk) {
      Track* tpcTrack = dynamic_cast<Track*>(col->getElementAt(iTrk));

      TrackerHitVec hitVec = tpcTrack->getTrackerHits();
      int nHits = int(hitVec.size());

      streamlog_out(DEBUG5) << toString(iTrk, tpcTrack, _bField) << std::endl;

      TrackExtended* trackExt = new TrackExtended(tpcTrack);

      trackExt->setOmega(tpcTrack->getOmega());
      trackExt->setTanLambda(tpcTrack->getTanLambda());
      trackExt->setPhi(tpcTrack->getPhi());
      trackExt->setD0(tpcTrack->getD0());
      trackExt->setZ0(tpcTrack->getZ0());
      float cov[15];
      //      float reso[4];
      // float param[5];
      // param[0] = tpcTrack->getOmega();
      // param[1] = tpcTrack->getTanLambda();
      // param[2] = tpcTrack->getPhi();
      // param[3] = tpcTrack->getD0();
      // param[4] = tpcTrack->getZ0();

      const FloatVec Cov = tpcTrack->getCovMatrix();
      int NC = int(Cov.size());
      for (int ic = 0; ic < NC; ic++) {
        cov[ic] = Cov[ic];
      }

      trackExt->setCovMatrix(cov);
      trackExt->setNDF(tpcTrack->getNdf());
      trackExt->setChi2(tpcTrack->getChi2());

      for (int iHit = 0; iHit < nHits; ++iHit) {
        TrackerHit* hit = hitVec[iHit];
        TrackerHitExtended* hitExt = mapTrackerHits[hit];
        hitExt->setTrackExtended(trackExt);
        trackExt->addTrackerHitExtended(hitExt);
      }

      _allTPCTracks.push_back(trackExt);
    }
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG5) << _TPCTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }

  // Reading Si Tracks
  try {
    LCCollection* col = event->getCollection(_SiTrackCollection.c_str());
    int nelem = col->getNumberOfElements();
    streamlog_out(DEBUG5) << std::endl;
    streamlog_out(DEBUG5) << "Number of Si Tracks = " << nelem << std::endl;
    streamlog_out(DEBUG5)
        << " Trk    ID        p          D0         Z0       Px       Py       Pz   hitsSi ndf Chi2/ndf" << std::endl;

    for (int iTrk = 0; iTrk < nelem; ++iTrk) {
      Track* siTrack = dynamic_cast<Track*>(col->getElementAt(iTrk));

      double prob = (siTrack->getNdf() > 0 ? gsl_cdf_chisq_Q(siTrack->getChi2(), (double)siTrack->getNdf()) : 0.);

      if (prob < _minChi2ProbForSiliconTracks) {
        streamlog_out(DEBUG5) << "Si Tracks " << siTrack << " id : " << siTrack->id() << " rejected with prob " << prob
                              << " < " << _minChi2ProbForSiliconTracks << std::endl;
        continue;
      }

      TrackExtended* trackExt = new TrackExtended(siTrack);
      TrackerHitVec hitVec = siTrack->getTrackerHits();
      int nHits = int(hitVec.size());
      trackExt->setOmega(siTrack->getOmega());
      trackExt->setTanLambda(siTrack->getTanLambda());
      trackExt->setPhi(siTrack->getPhi());
      trackExt->setD0(siTrack->getD0());
      trackExt->setZ0(siTrack->getZ0());
      float cov[15];
      // float param[5];
      // param[0] = siTrack->getOmega();
      // param[1] = siTrack->getTanLambda();
      // param[2] = siTrack->getPhi();
      // param[3] = siTrack->getD0();
      // param[4] = siTrack->getZ0();

      const FloatVec Cov = siTrack->getCovMatrix();
      int NC = int(Cov.size());
      for (int ic = 0; ic < NC; ic++) {
        cov[ic] = Cov[ic];
      }
      //      }
      trackExt->setCovMatrix(cov);
      trackExt->setNDF(siTrack->getNdf());
      trackExt->setChi2(siTrack->getChi2());
      char strg[200];
      HelixClass helixSi;
      for (int iHit = 0; iHit < nHits; ++iHit) {
        TrackerHit* hit = hitVec[iHit];
        TrackerHitExtended* hitExt = mapTrackerHits[hit];
        hitExt->setTrackExtended(trackExt);

        trackExt->addTrackerHitExtended(hitExt);
      }

      float d0Si = trackExt->getD0();
      float z0Si = trackExt->getZ0();
      float omegaSi = trackExt->getOmega();
      float phi0Si = trackExt->getPhi();
      float tanLSi = trackExt->getTanLambda();
      helixSi.Initialize_Canonical(phi0Si, d0Si, z0Si, omegaSi, tanLSi, _bField);
      float pxSi = helixSi.getMomentum()[0];
      float pySi = helixSi.getMomentum()[1];
      float pzSi = helixSi.getMomentum()[2];
      const float pTot = sqrt(pxSi * pxSi + pySi * pySi + pzSi * pzSi);
      const int ndfSi = trackExt->getNDF();
      float Chi2Si = trackExt->getChi2() / float(trackExt->getNDF());
      sprintf(strg, "%3i   %5i %9.3f  %9.3f  %9.3f  %7.2f  %7.2f  %7.2f %4i %4i %8.3f", iTrk, siTrack->id(), pTot, d0Si,
              z0Si, pxSi, pySi, pzSi, nHits, ndfSi, Chi2Si);
      streamlog_out(DEBUG5) << strg << std::endl;

      if (nHits > 0) {
        _allSiTracks.push_back(trackExt);
      } else {
        delete trackExt;
      }
    }

    streamlog_out(DEBUG5) << std::endl;
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG5) << _SiTrackCollection.c_str() << " collection is unavailable" << std::endl;
  }
}

void FullLDCTracking_MarlinTrk::CleanUp() {
  int nNonCombTpc = int(_allNonCombinedTPCTracks.size());
  for (int i = 0; i < nNonCombTpc; ++i) {
    TrackExtended* trkExt = _allNonCombinedTPCTracks[i];
    GroupTracks* group = trkExt->getGroupTracks();
    delete group;
  }
  _allNonCombinedTPCTracks.clear();

  int nNonCombSi = int(_allNonCombinedSiTracks.size());
  for (int i = 0; i < nNonCombSi; ++i) {
    TrackExtended* trkExt = _allNonCombinedSiTracks[i];
    GroupTracks* group = trkExt->getGroupTracks();
    delete group;
  }
  _allNonCombinedSiTracks.clear();

  int nSITHits = int(_allSITHits.size());
  for (int i = 0; i < nSITHits; ++i) {
    TrackerHitExtended* hitExt = _allSITHits[i];
    delete hitExt;
  }
  _allSITHits.clear();

  int nSETHits = int(_allSETHits.size());
  for (int i = 0; i < nSETHits; ++i) {
    TrackerHitExtended* hitExt = _allSETHits[i];
    delete hitExt;
  }
  _allSETHits.clear();

  int nTPCHits = int(_allTPCHits.size());
  for (int i = 0; i < nTPCHits; ++i) {
    TrackerHitExtended* hitExt = _allTPCHits[i];
    delete hitExt;
  }
  _allTPCHits.clear();

  int nFTDHits = int(_allFTDHits.size());
  for (int i = 0; i < nFTDHits; ++i) {
    TrackerHitExtended* hitExt = _allFTDHits[i];
    delete hitExt;
  }
  _allFTDHits.clear();

  int nETDHits = int(_allETDHits.size());
  for (int i = 0; i < nETDHits; ++i) {
    TrackerHitExtended* hitExt = _allETDHits[i];
    delete hitExt;
  }
  _allETDHits.clear();

  int nVTXHits = int(_allVTXHits.size());
  for (int i = 0; i < nVTXHits; ++i) {
    TrackerHitExtended* hitExt = _allVTXHits[i];
    delete hitExt;
  }
  _allVTXHits.clear();

  int nSiTrk = int(_allSiTracks.size());
  for (int i = 0; i < nSiTrk; ++i) {
    TrackExtended* trkExt = _allSiTracks[i];
    delete trkExt;
  }
  _allSiTracks.clear();

  int nTPCTrk = int(_allTPCTracks.size());
  for (int i = 0; i < nTPCTrk; ++i) {
    TrackExtended* trkExt = _allTPCTracks[i];
    delete trkExt;
  }
  _allTPCTracks.clear();

  int nCombTrk = int(_allCombinedTracks.size());
  for (int i = 0; i < nCombTrk; ++i) {
    TrackExtended* trkExt = _allCombinedTracks[i];
    GroupTracks* group = trkExt->getGroupTracks();
    delete group;
    delete trkExt;
  }
  _allCombinedTracks.clear();

  //   int nImplTrk = int(_trkImplVec.size());
  //   for (int i=0;i<nImplTrk;++i) {
  //     TrackExtended * trkImpl = _trkImplVec[i];
  //     delete trkImpl;
  //   }
  _trkImplVec.clear();

  // AS: Dont delete the individual entries, some of them are cleared elsewhere, I think
  _candidateCombinedTracks.clear();
}

/*

 Compare all TPC tracks with all Silicon tracks
 and compare using CompareTrkII which really only compares delta d0 and delta z0
 then calculate dOmega and the angle between the tracks, which are checked against cuts here.

 if CompareTrkII and the cuts on dOmega and angle succecede try and full fit of the hits using CombineTracks
 assuming that the merger was not vetoed VetoMerge
 i.e. if the momentum of either track is less than 2.5 GeV
 or if following a full fit the NDF+10 of the combined tracks is less than the NDF_first + NDF_second

 if this is successful the track is added to _allCombinedTracks and the TPC and Si Segements are added to
 _candidateCombinedTracks

 */

void FullLDCTracking_MarlinTrk::MergeTPCandSiTracks() {
  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks = int(_allSiTracks.size());

  streamlog_out(DEBUG3) << " MergeTPCandSiTracks called nTPC tracks " << nTPCTracks << " - nSiTracks " << nSiTracks
                        << std::endl;

  for (int iTPC = 0; iTPC < nTPCTracks; ++iTPC) {
    TrackExtended* tpcTrackExt = _allTPCTracks[iTPC];
    for (int iSi = 0; iSi < nSiTracks; ++iSi) {
      TrackExtended* siTrackExt = _allSiTracks[iSi];
      int iComp = 0;
      float angle = 0;

      streamlog_out(DEBUG2) << " compare tpc trk " << toString(iTPC, tpcTrackExt->getTrack(), _bField) << std::endl;
      streamlog_out(DEBUG2) << "    to si trk    " << toString(iSi, siTrackExt->getTrack(), _bField) << std::endl;

      float dOmega = CompareTrkII(siTrackExt, tpcTrackExt, _d0CutForMerging, _z0CutForMerging, iComp, angle);

      if ((dOmega < _dOmegaForMerging) && (angle < _angleForMerging) && !VetoMerge(tpcTrackExt, siTrackExt)) {
        streamlog_out(DEBUG2) << " call CombineTracks for tpc trk " << tpcTrackExt << " si trk " << siTrackExt
                              << std::endl;

        TrackExtended* combinedTrack =
            CombineTracks(tpcTrackExt, siTrackExt, _maxAllowedPercentageOfOutliersForTrackCombination, false);

        streamlog_out(DEBUG2) << " combinedTrack returns " << combinedTrack << std::endl;

        if (combinedTrack != NULL) {
          _allCombinedTracks.push_back(combinedTrack);
          _candidateCombinedTracks.insert(tpcTrackExt);
          _candidateCombinedTracks.insert(siTrackExt);

          //          streamlog_out(DEBUG3) << " combinedTrack successfully added to _allCombinedTracks : " << toString(
          //          0,combinedTrack->getTrack(), _bField  )  << std::endl;
          streamlog_out(DEBUG3) << " *** combinedTrack successfully added to _allCombinedTracks : tpc " << iTPC
                                << " si " << iSi << std::endl;

          if (_debug >= 3) {
            int iopt = 1;
            PrintOutMerging(tpcTrackExt, siTrackExt, iopt);
          }
        } else {
          if (_debug >= 3) {
            int iopt = 6;
            PrintOutMerging(tpcTrackExt, siTrackExt, iopt);
          }
        }
      } else {
        if (_debug >= 3) {
          int iopt = 6;
          PrintOutMerging(tpcTrackExt, siTrackExt, iopt);
        }
      }
    }
  }
}

/*

 All Si and TPC tracks which have been merged using MergeTPCandSiTracks are excluded from this search.

 The remaining TPC tracks and Silicon tracks will be tested using CompareTrkIII

 CompareTrkIII does the following

 i)   significance d0    < d0Cut
 ii)  significance z0    < z0Cut
 iii) pdot > 0.999

 if above three cuts are passed return:

 the significance dAngle (return by reference)
 and
 the significance dOmega

 if CompareTrkIII and the cuts on the significance of dOmega and angle succecede try and full fit of the hits using
 CombineTracks assuming that the merger was not vetoed VetoMerge i.e. if the momentum of either track is less than 2.5
 GeV and there is no overlapping of the segments

 if this is successful the track is added to _allCombinedTracks and the TPC and Si Segements are added to
 _candidateCombinedTracks

 */

void FullLDCTracking_MarlinTrk::MergeTPCandSiTracksII() {
  int nTPCTracks = int(_allTPCTracks.size());
  int nSiTracks = int(_allSiTracks.size());

  streamlog_out(DEBUG3) << " MergeTPCandSiTracksII called nTPC tracks " << nTPCTracks << " - nSiTracks " << nSiTracks
                        << std::endl;

  for (int iTPC = 0; iTPC < nTPCTracks; ++iTPC) {
    // check if the tpc track has already been merged with CompareTrkII
    TrackExtended* tpcTrackExt = _allTPCTracks[iTPC];
    if (_candidateCombinedTracks.find(tpcTrackExt) != _candidateCombinedTracks.end())
      continue;

    for (int iSi = 0; iSi < nSiTracks; ++iSi) {
      // check if the tpc track has already been merged with CompareTrkII
      TrackExtended* siTrackExt = _allSiTracks[iSi];
      if (_candidateCombinedTracks.find(siTrackExt) != _candidateCombinedTracks.end())
        continue;

      int iComp = 0;
      float angleSignificance = 0;

      float significance =
          CompareTrkIII(siTrackExt, tpcTrackExt, _d0CutForMerging, _z0CutForMerging, iComp, angleSignificance);

      streamlog_out(DEBUG2) << " MergeTPCandSiTracksII - tpctrk " << iTPC << " - " << iSi << " - significance "
                            << significance << " angleSignificance " << angleSignificance << std::endl;

      if ((significance < 10) && (angleSignificance < 5) && !VetoMerge(tpcTrackExt, siTrackExt)) {
        TrackExtended* combinedTrack =
            CombineTracks(tpcTrackExt, siTrackExt, _maxAllowedPercentageOfOutliersForTrackCombination, false);

        streamlog_out(DEBUG2) << " combinedTrack returns " << combinedTrack << std::endl;

        if (combinedTrack != NULL) {
          _allCombinedTracks.push_back(combinedTrack);
          streamlog_out(DEBUG3) << " *** combinedTrack successfully added to _allCombinedTracks : tpc " << iTPC
                                << " si " << iSi << std::endl;

          if (_debug >= 3) {
            int iopt = 1;
            PrintOutMerging(tpcTrackExt, siTrackExt, iopt);
          }
        } else {
          if (_debug >= 3) {
            int iopt = 6;
            PrintOutMerging(tpcTrackExt, siTrackExt, iopt);
          }
        }
      } else {
        if (_debug >= 3) {
          int iopt = 6;
          PrintOutMerging(tpcTrackExt, siTrackExt, iopt);
        }
      }
    }
  }
}

// if testCombinationOnly is true then hits will not be assigned to the tracks
TrackExtended* FullLDCTracking_MarlinTrk::CombineTracks(TrackExtended* tpcTrack, TrackExtended* siTrack,
                                                        float maxAllowedOutliers, bool testCombinationOnly) {
  TrackExtended* OutputTrack = NULL;

  TrackerHitExtendedVec siHitVec = siTrack->getTrackerHitExtendedVec();
  TrackerHitExtendedVec tpcHitVec = tpcTrack->getTrackerHitExtendedVec();

  int nSiHits = int(siHitVec.size());
  int nTPCHits = int(tpcHitVec.size());
  int nHits = nTPCHits + nSiHits;

  // std::cout << "FullLDCTracking_MarlinTrk::CombineTracks nSiHits = " << nSiHits << std::endl;
  // std::cout << "FullLDCTracking_MarlinTrk::CombineTracks nTPCHits = " << nTPCHits << std::endl;

  EVENT::TrackerHitVec trkHits;
  trkHits.reserve(nHits);

  for (int ih = 0; ih < nSiHits; ++ih) {
    TrackerHit* trkHit = siHitVec[ih]->getTrackerHit();
    if (trkHit) {
      trkHits.push_back(trkHit);
    } else {
      throw EVENT::Exception(std::string("FullLDCTracking_MarlinTrk::CombineTracks: TrackerHit pointer == NULL "));
    }
  }

  for (int ih = 0; ih < nTPCHits; ++ih) {
    TrackerHit* trkHit = tpcHitVec[ih]->getTrackerHit();
    if (trkHit) {
      trkHits.push_back(trkHit);
    } else {
      throw EVENT::Exception(std::string("FullLDCTracking_MarlinTrk::CombineTracks: TrackerHit pointer == NULL "));
    }
  }

  double chi2_D;
  int ndf;

  if (trkHits.size() < 3) {
    return 0;
  }

  streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: Sorting Hits " << trkHits.size() << std::endl;

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

  for (std::vector<std::pair<float, EVENT::TrackerHit*>>::iterator it = r2_values.begin(); it != r2_values.end();
       ++it) {
    trkHits.push_back(it->second);
  }

  streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: Start Fitting: AddHits: number of hits to fit "
                        << trkHits.size() << std::endl;

  std::unique_ptr<MarlinTrk::IMarlinTrack> marlin_trk_autop(_trksystem->createTrack());
  MarlinTrk::IMarlinTrack& marlin_trk = *marlin_trk_autop.get();

  IMPL::TrackStateImpl pre_fit;

  int error = IMarlinTrack::success;

  pre_fit = *(tpcTrack->getTrack()->getTrackState(EVENT::TrackState::AtLastHit));

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

  pre_fit.setCovMatrix(covMatrix);

  error = MarlinTrk::createFit(trkHits, &marlin_trk, &pre_fit, _bField, IMarlinTrack::backward, _maxChi2PerHit);

  if (error != IMarlinTrack::success) {
    streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: creation of fit fails with error " << error
                          << std::endl;
    return 0;
  }

  const Vector3D point(0., 0., 0.); // nominal IP

  TrackStateImpl trkState;
  error = marlin_trk.propagate(point, trkState, chi2_D, ndf);

  if (error != IMarlinTrack::success) {
    streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::CombineTracks: propagate to IP fails with error " << error
                          << std::endl;
    return 0;
  }

  if (ndf < 0) {
    streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: Fit failed NDF is less that zero  " << ndf
                          << std::endl;
    return 0;
  }

  float chi2Fit = chi2_D / float(ndf);

  if (chi2Fit > _chi2FitCut) {
    streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: track fail Chi2 cut of " << _chi2FitCut
                          << " chi2 of track = " << chi2Fit << std::endl;
    return 0;
  }

  streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: Check for outliers " << std::endl;

  std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;
  marlin_trk.getOutliers(outliers);

  float outlier_pct = outliers.size() / float(trkHits.size());

  if (outlier_pct > maxAllowedOutliers) {
    streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: percentage of outliers " << outlier_pct
                          << " is greater than cut maximum: " << maxAllowedOutliers << std::endl;
    return 0;
  }

  // sort the hits into outliers from TPC and Silicon as we will reject the combination if more that 2 Si hits get
  // rejected ...

  std::vector<TrackerHitExtended*> siHitInFit;
  std::vector<TrackerHitExtended*> siOutliers;
  std::vector<TrackerHitExtended*> tpcHitInFit;
  std::vector<TrackerHitExtended*> tpcOutliers;

  for (int i = 0; i < nSiHits; ++i) {
    bool hit_is_outlier = false;

    // we need to make sure that in the case of a composite hit we reject this as well
    LCObjectVec hits;

    // all hits, both the 2D tracker hit, as well as any raw hits which belong to it
    hits.push_back(siHitVec[i]->getTrackerHit());

    // add the raw hits ...
    const LCObjectVec rawHits = siHitVec[i]->getTrackerHit()->getRawHits();

    if (rawHits.empty() == false) {
      for (unsigned ihit = 0; ihit < rawHits.size(); ++ihit) {
        hits.push_back(rawHits[ihit]);
      }
    }

    // now double loop over the outliers and the hits assosiated with this TrackerHitExtended and compare
    for (unsigned ohit = 0; ohit < outliers.size(); ++ohit) {
      for (unsigned ihit = 0; ihit < hits.size(); ++ihit) {
        // compare outlier pointer to TrackerHit pointer
        if (outliers[ohit].first == hits[ihit]) {
          // silicon outlier found so add the TrackerHitExtended to the list of outliers
          hit_is_outlier = true;
          siOutliers.push_back(siHitVec[i]);
          break;
        }
        if (hit_is_outlier) {
          break;
        }
      }
    }

    if (hit_is_outlier == false) {
      // add the TrackerHitExtended to the list of silicon hits used in the fit
      siHitInFit.push_back(siHitVec[i]);
    }
  }

  // more simple as TPC hits are never composite
  for (int i = 0; i < nTPCHits; ++i) {
    bool hit_is_outlier = false;

    for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
      // compare outlier pointer to TrackerHit pointer
      if (outliers[ihit].first == tpcHitVec[i]->getTrackerHit()) {
        // tpc outlier found so add the TrackerHitExtended to the list of outliers
        hit_is_outlier = true;
        tpcOutliers.push_back(tpcHitVec[i]);
        break;
      }
    }

    if (hit_is_outlier == false) {
      // add the TrackerHitExtended to the list of tpc hits used in the fit
      tpcHitInFit.push_back(tpcHitVec[i]);
    }
  }

  streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: Check for Silicon Hit rejections ... "
                        << std::endl;

  if ((int)siOutliers.size() > _maxAllowedSiHitRejectionsForTrackCombination) {
    streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: Fit rejects " << siOutliers.size()
                          << " silicon hits : max allowed rejections = "
                          << _maxAllowedSiHitRejectionsForTrackCombination << " : Combination rejected " << std::endl;
    return 0;
  }

  float omega = trkState.getOmega();
  float tanlambda = trkState.getTanLambda();
  float phi0 = trkState.getPhi();
  float d0 = trkState.getD0();
  float z0 = trkState.getZ0();

  OutputTrack = new TrackExtended();

  GroupTracks* group = new GroupTracks();
  OutputTrack->setGroupTracks(group);

  group->addTrackExtended(siTrack);
  group->addTrackExtended(tpcTrack);

  // note OutputTrack which is of type TrackExtended, only takes fits set for ref point = 0,0,0
  OutputTrack->setOmega(omega);
  OutputTrack->setTanLambda(tanlambda);
  OutputTrack->setPhi(phi0);
  OutputTrack->setZ0(z0);
  OutputTrack->setD0(d0);
  OutputTrack->setChi2(chi2_D);
  OutputTrack->setNDF(ndf);

  float cov[15];

  for (int i = 0; i < 15; ++i) {
    cov[i] = trkState.getCovMatrix().operator[](i);
  }

  OutputTrack->setCovMatrix(cov);

  // if this is not just a test of the combination add the hits to the combined track
  if (testCombinationOnly == false) {
    for (unsigned ihit = 0; ihit < siHitInFit.size(); ++ihit) {
      TrackerHitExtended* hitExt = siHitInFit[ihit];
      OutputTrack->addTrackerHitExtended(hitExt);
      hitExt->setUsedInFit(true);
    }

    for (unsigned ihit = 0; ihit < siOutliers.size(); ++ihit) {
      TrackerHitExtended* hitExt = siOutliers[ihit];
      OutputTrack->addTrackerHitExtended(hitExt);
      hitExt->setUsedInFit(false);
    }

    for (unsigned ihit = 0; ihit < tpcHitInFit.size(); ++ihit) {
      TrackerHitExtended* hitExt = tpcHitInFit[ihit];
      OutputTrack->addTrackerHitExtended(hitExt);
      hitExt->setUsedInFit(true);
    }

    for (unsigned ihit = 0; ihit < tpcOutliers.size(); ++ihit) {
      TrackerHitExtended* hitExt = tpcOutliers[ihit];
      OutputTrack->addTrackerHitExtended(hitExt);
      hitExt->setUsedInFit(false);
    }
  }

  streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::CombineTracks: merged track created  " << OutputTrack << " with "
                        << OutputTrack->getTrackerHitExtendedVec().size() << " hits, nhits tpc " << nTPCHits
                        << " nSiHits " << nSiHits << ", testCombinationOnly = " << testCombinationOnly << std::endl;

  return OutputTrack;
}

void FullLDCTracking_MarlinTrk::SortingTrackHitPairs(TrackHitPairVec& trackHitPairVec) {
  int sizeOfVector = int(trackHitPairVec.size());
  TrackHitPair *one, *two, *Temp;

  for (int i = 0; i < sizeOfVector - 1; i++)
    for (int j = 0; j < sizeOfVector - i - 1; j++) {
      one = trackHitPairVec[j];
      two = trackHitPairVec[j + 1];
      float oneQ = one->getDistance();
      float twoQ = two->getDistance();
      if (oneQ > twoQ) {
        Temp = trackHitPairVec[j];
        trackHitPairVec[j] = trackHitPairVec[j + 1];
        trackHitPairVec[j + 1] = Temp;
      }
    }
}

/*

 Sorts all tracks in the vector by Chi2/NDF

 */

void FullLDCTracking_MarlinTrk::Sorting(TrackExtendedVec& trackVec) {
  int sizeOfVector = int(trackVec.size());
  TrackExtended *one, *two, *Temp;

  for (int i = 0; i < sizeOfVector - 1; i++)
    for (int j = 0; j < sizeOfVector - i - 1; j++) {
      one = trackVec[j];
      two = trackVec[j + 1];
      float oneQ = one->getChi2() / float(one->getNDF());
      float twoQ = two->getChi2() / float(two->getNDF());
      if (oneQ > twoQ) {
        Temp = trackVec[j];
        trackVec[j] = trackVec[j + 1];
        trackVec[j + 1] = Temp;
      }
    }
}

/*

 this function is used to select only the combined tracks which have been formed from excatly 2 sub tracks which
 themselves are not formed from any other tracks.

 */

void FullLDCTracking_MarlinTrk::SelectCombinedTracks() {
  int nCombTrk = int(_allCombinedTracks.size());

  streamlog_out(DEBUG3) << " **SelectCombinedTracks - check " << nCombTrk << " comb. tracks " << std::endl;

  // loop over all combined tracks ...
  for (int i = 0; i < nCombTrk; ++i) {
    TrackExtended* trkExt = _allCombinedTracks[i];

    streamlog_out(DEBUG1) << " **SelectCombinedTracks - Check Track " << trkExt << std::endl;

    // get the sub tracks which have been combined to form this track
    GroupTracks* group = trkExt->getGroupTracks();
    TrackExtendedVec tracks = group->getTrackExtendedVec();

    // check that there are only 2 sub tracks, as we are after Si <-> TPC mergers only
    int nTracks = int(tracks.size());

    streamlog_out(DEBUG1) << " **SelectCombinedTracks - nTracks = " << nTracks << std::endl;

    if (nTracks == 2) {
      TrackExtended* firstTrack = tracks[0];
      TrackExtended* secondTrack = tracks[1];

      // check that the two sub tracks in question are not themselves a merger of tracks
      if ((firstTrack->getGroupTracks() == NULL) && (secondTrack->getGroupTracks() == NULL)) {
        streamlog_out(DEBUG1) << " **SelectCombinedTracks - firstTrack->getGroupTracks() == NULL ... " << std::endl;

        // associate the current group to the two sub tracks
        firstTrack->setGroupTracks(group);
        secondTrack->setGroupTracks(group);

        // get the tracker hits ...
        TrackerHitExtendedVec firstVec = firstTrack->getTrackerHitExtendedVec();
        TrackerHitExtendedVec secondVec = secondTrack->getTrackerHitExtendedVec();

        // number of hits for first and second sub track
        int nFirst = int(firstVec.size());
        int nSecond = int(secondVec.size());

        // use these to store the min and max z positions for the combination.
        float edges[2];
        edges[0] = 1.0e+20;
        edges[1] = -1.0e+20;

        // get min and max z for the first sub track
        for (int iF = 0; iF < nFirst; ++iF) {
          TrackerHitExtended* trkHitExt = firstVec[iF];
          TrackerHit* trkHit = trkHitExt->getTrackerHit();
          float zpos = float(trkHit->getPosition()[2]);
          if (zpos > edges[1])
            edges[1] = zpos;
          if (zpos < edges[0])
            edges[0] = zpos;
        }

        // get min and max z for the second sub track
        for (int iS = 0; iS < nSecond; ++iS) {
          TrackerHitExtended* trkHitExt = secondVec[iS];
          TrackerHit* trkHit = trkHitExt->getTrackerHit();
          float zpos = float(trkHit->getPosition()[2]);
          if (zpos > edges[1])
            edges[1] = zpos;
          if (zpos < edges[0])
            edges[0] = zpos;
        }

        // record the z extreams to the group ...
        group->setEdges(edges);
        // ... and add the combined track to the list.
        _trkImplVec.push_back(trkExt);
        streamlog_out(DEBUG2) << " add track " << trkExt << " to combined final list " << std::endl;

        if (_debug >= 2) {
          int iopt = 1;

          // here it is assumed that the tpc tracks is the secondTrack ...
          PrintOutMerging(secondTrack, firstTrack, iopt);
        }
      }
    } else { // if(nTracks>2)

      streamlog_out(DEBUG2) << " *****************  SelectCombinedTracks: MORE THAN TWO TRACKS " << nCombTrk
                            << std::endl;
    }
  }
}

void FullLDCTracking_MarlinTrk::AddNotCombinedTracks() {
  int nTPCTrk = int(_allTPCTracks.size());
  int nSiTrk = int(_allSiTracks.size());

  // we need some buffer vector
  TrackExtendedVec allMergedTracks;
  allMergedTracks.clear();

  // forcing merging of Si and TPC track segments
  if (_forceMerging == 1) {
    // loop over all TPC tracks
    for (int i = 0; i < nTPCTrk; ++i) {
      // get the tracks assigned to this TPC track
      TrackExtended* trkExtTPC = _allTPCTracks[i];
      GroupTracks* groupTPC = trkExtTPC->getGroupTracks();

      // if no tracks have been grouped with this TPC track
      if (groupTPC == NULL) {
        float diffMin = 1.0e+20;

        TrackExtended* siTrkToAttach = NULL;

        // loop over all Silicon Tracks
        for (int j = 0; j < nSiTrk; ++j) {
          TrackExtended* trkExtSi = _allSiTracks[j];
          GroupTracks* groupSi = trkExtSi->getGroupTracks();

          // only consider ungrouped Silicon Tracks
          if (groupSi == NULL) {
            int iComp = 0;
            //      float deltaP = CompareTrk(trkExtSi,trkExtTPC,_d0CutForForcedMerging,_z0CutForForcedMerging,iComp);
            float angle(0.);
            float angleSignificance(0.);

            // try to merge tracks using looser cuts
            float dOmega =
                CompareTrkII(trkExtSi, trkExtTPC, _d0CutForForcedMerging, _z0CutForForcedMerging, iComp, angle);

            float significance = CompareTrkIII(trkExtSi, trkExtTPC, _d0CutForForcedMerging, _z0CutForForcedMerging,
                                               iComp, angleSignificance);

            //      if (deltaP < _dPCutForForcedMerging) {

            if (((dOmega < _dOmegaForForcedMerging) && (angle < _angleForForcedMerging)) ||
                ((significance < 5) && (angleSignificance < 5))) {
              float chi2O = dOmega / _dOmegaForForcedMerging;
              float chi2A = angle / _angleForForcedMerging;
              float deltaP = chi2O * chi2O + chi2A * chi2A;

              // if this is the best match (diffMin) set the possible merger
              if (deltaP < diffMin) {
                diffMin = deltaP;
                siTrkToAttach = trkExtSi;
              }
            } else {
              if (_debug >= 3) {
                int iopt = 7;
                streamlog_out(DEBUG2) << significance << " " << angleSignificance << std::endl;
                PrintOutMerging(trkExtTPC, trkExtSi, iopt);
              }
            }
          }
        }

        if (siTrkToAttach != NULL) {
          TrackExtended* trkExtSi = siTrkToAttach;
          TrackExtended* OutputTrack = new TrackExtended();
          GroupTracks* group = new GroupTracks();

          group->addTrackExtended(trkExtSi);
          group->addTrackExtended(trkExtTPC);

          OutputTrack->setGroupTracks(group);
          //        trkExtSi->setGroupTracks(group);
          //        trkExtTPC->setGroupTracks(group);

          OutputTrack->setOmega(trkExtTPC->getOmega());
          OutputTrack->setTanLambda(trkExtSi->getTanLambda());
          OutputTrack->setPhi(trkExtSi->getPhi());
          OutputTrack->setZ0(trkExtSi->getZ0());
          OutputTrack->setD0(trkExtSi->getD0());

          float covMatTPC[15];
          float covMatSi[15];
          float covMat[15];
          for (int iCov = 0; iCov < 15; ++iCov) {
            covMatTPC[iCov] = trkExtTPC->getCovMatrix()[iCov];
            covMatSi[iCov] = trkExtSi->getCovMatrix()[iCov];
            covMat[iCov] = covMatSi[iCov];
          }

          float scaling = sqrt(covMatTPC[5] / covMatSi[5]);
          covMat[5] = covMatTPC[5];
          covMat[3] = scaling * covMatSi[3];
          covMat[4] = scaling * covMatSi[4];
          covMat[8] = scaling * covMatSi[8];
          covMat[12] = scaling * covMatSi[12];

          OutputTrack->setCovMatrix(covMat);
          TrackerHitExtendedVec tpcHitVec = trkExtTPC->getTrackerHitExtendedVec();
          TrackerHitExtendedVec siHitVec = trkExtSi->getTrackerHitExtendedVec();

          int nTPCHits = int(tpcHitVec.size());
          int nSiHits = int(siHitVec.size());

          float edges[2];
          edges[0] = 1.0e+20;
          edges[1] = -1.0e+20;

          // find the max and min z extents from hits
          for (int iH = 0; iH < nSiHits; ++iH) {
            TrackerHitExtended* hitExt = siHitVec[iH];
            OutputTrack->addTrackerHitExtended(hitExt);
            hitExt->setUsedInFit(true);
            TrackerHit* hit = hitExt->getTrackerHit();
            float zpos = float(hit->getPosition()[2]);
            if (zpos < edges[0])
              edges[0] = zpos;
            if (zpos > edges[1])
              edges[1] = zpos;
          }
          for (int iH = 0; iH < nTPCHits; ++iH) {
            TrackerHitExtended* hitExt = tpcHitVec[iH];
            OutputTrack->addTrackerHitExtended(hitExt);
            hitExt->setUsedInFit(true);
            TrackerHit* hit = hitExt->getTrackerHit();
            float zpos = float(hit->getPosition()[2]);
            if (zpos < edges[0])
              edges[0] = zpos;
            if (zpos > edges[1])
              edges[1] = zpos;
          }

          group->setEdges(edges);
          OutputTrack->setChi2(diffMin); // will be replaced if necessary
          OutputTrack->setNDF(int(1));   // will be replaced if necessary

          _allCombinedTracks.push_back(OutputTrack);
          allMergedTracks.push_back(OutputTrack);
        }
      }

    } // end of loop over TPC tracks

    // check that there are some merged tracks to process
    int nMerged = int(allMergedTracks.size());

    if (nMerged > 0) {
      // sort all merged tracks by Chi2/NDF
      // although due to the fact that NDF for these tracks is set to the value 1 above,
      // it is really only a sort on Chi2 which was really the weighted difference in angle and omega

      Sorting(allMergedTracks);

      // loop over all merged tracks
      for (int iM = 0; iM < nMerged; ++iM) {
        TrackExtended* mergedTrack = allMergedTracks[iM];
        GroupTracks* grpTrk = mergedTrack->getGroupTracks();

        TrackExtendedVec trkVec = grpTrk->getTrackExtendedVec();
        TrackExtended* trkTPC = NULL;
        TrackExtended* trkSi = NULL;

        int nT = int(trkVec.size());

        // only consider tracks which have been composed of excactly 2 sub tracks
        if (nT == 2) {
          trkTPC = trkVec[0];
          trkSi = trkVec[1];

          GroupTracks* groupTPC = trkTPC->getGroupTracks();
          GroupTracks* groupSi = trkSi->getGroupTracks();

          // check that both the TPC and SI track have not already been combined with other tracks ...
          if (groupTPC == NULL && groupSi == NULL) {
            // set the grouping, meaning that these tracks will not be considered further
            trkTPC->setGroupTracks(grpTrk);
            trkSi->setGroupTracks(grpTrk);

            TrackerHitExtendedVec hitVec = mergedTrack->getTrackerHitExtendedVec();

            int nhits = int(hitVec.size());

            int totNdf = 2 * nhits - 5;
            float totChi2 = trkTPC->getChi2() + trkSi->getChi2();

            mergedTrack->setNDF(totNdf);
            mergedTrack->setChi2(totChi2);

            if (_debug >= 2) {
              int iopt = 2;
              PrintOutMerging(trkTPC, trkSi, iopt);
            }
            _trkImplVec.push_back(mergedTrack);
          }
        }
      }
    }
  } // end of _forceMerging

  // clear buffer vector
  allMergedTracks.clear();

  // merging split up TPC segments
  if (_mergeTPCSegments) {
    std::vector<GroupTracks*> TPCSegments;
    TPCSegments.clear();

    int nNonAssignedTPCSeg = 0;

    // loop over all TPC Tracks
    for (int i = 0; i < nTPCTrk; ++i) {
      TrackExtended* trkExt = _allTPCTracks[i];
      GroupTracks* group = trkExt->getGroupTracks();

      streamlog_out(DEBUG2) << " *****************  AddNotCombinedTracks: Check track " << trkExt
                            << " id = " << trkExt->getTrack()->id() << std::endl;

      // only consider those tracks which have not yet been combined
      if (group == NULL) {
        // find the min and max z extents using the hits
        TrackerHitExtendedVec currentVec = trkExt->getTrackerHitExtendedVec();
        int nCur = int(currentVec.size());
        float zmin = 1e+20;
        float zmax = -1e+20;

        for (int iCur = 0; iCur < nCur; ++iCur) {
          TrackerHitExtended* curTrkHitExt = currentVec[iCur];
          TrackerHit* curTrkHit = curTrkHitExt->getTrackerHit();
          float zpos = float(curTrkHit->getPosition()[2]);
          if (zpos < zmin)
            zmin = zpos;
          if (zpos > zmax)
            zmax = zpos;
        }

        nNonAssignedTPCSeg++;

        // current number of TPC segment groupings
        int nGroups = int(TPCSegments.size());

        float dPtMin = 1.0e+10;
        GroupTracks* groupToAttach = NULL;
        TrackExtended* trkToAttach = NULL;

        // loop over the current TPC segment groupings
        for (int iG = 0; iG < nGroups; ++iG) {
          GroupTracks* segments = TPCSegments[iG];
          TrackExtendedVec segVec = segments->getTrackExtendedVec();

          // number of segments with the candidate group
          int nTrk = int(segVec.size());
          bool consider = true;

          if (_forbidOverlapInZTPC == 1) { // if overlap in Z of the two segments is forbidden

            // loop over all tracks in the current grouping
            for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
              TrackExtended* trkInGroup = segVec[iTrk];

              // get the number of hits from the track
              TrackerHitExtendedVec hitInGroupVec = trkInGroup->getTrackerHitExtendedVec();
              int nHitsInGrp = int(hitInGroupVec.size());

              // loop over the hits and make sure that there is no overlap in z
              for (int iHitInGrp = 0; iHitInGrp < nHitsInGrp; iHitInGrp++) {
                TrackerHitExtended* xTrkExt = hitInGroupVec[iHitInGrp];
                TrackerHit* xTrk = xTrkExt->getTrackerHit();

                float xZ = float(xTrk->getPosition()[2]);
                if (xZ > zmin && xZ < zmax) {
                  consider = false;
                  break;
                }
              }
              if (!consider)
                break; // if the candiate track's min and max z are within that of the group
            }
          }

          if (consider) {
            // again loop over the tracks in the current group
            for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
              TrackExtended* trkInGroup = segVec[iTrk];
              int iComp = 1;

              // compare the tracks ...
              float dPt = CompareTrk(trkExt, trkInGroup, _d0CutToMergeTPC, _z0CutToMergeTPC, iComp);

              // check that this tracks has the lowest delta pt  and vetomerge (i.e. fullfill momentum cut and makes
              // sure a large number of hits have not been lost in the merger)
              if (dPt < dPtMin && !VetoMerge(trkExt, trkInGroup)) {
                dPtMin = dPt;
                groupToAttach = segments;
                trkToAttach = trkInGroup;
                if (_debug >= 3) {
                  int iopt = 5;
                  PrintOutMerging(trkExt, trkInGroup, iopt);
                }
              } else {
                if (_debug >= 3) {
                  int iopt = 9;
                  PrintOutMerging(trkExt, trkInGroup, iopt);
                }
              }
            }
          } else {
            if (_debug >= 3) {
              int iopt = 9;
              for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
                TrackExtended* trkInGroup = segVec[iTrk];
                int iComp = 1;
                float dPt = CompareTrk(trkExt, trkInGroup, _d0CutToMergeTPC, _z0CutToMergeTPC, iComp);
                if (dPt >= dPtMin) {
                  PrintOutMerging(trkExt, trkInGroup, iopt);
                }
              }
            }
          }
        }

        // check the pt cut for merging and that a group has been found to match ..
        if (dPtMin < _dPCutToMergeTPC && groupToAttach != NULL) {
          // add the track to the group
          groupToAttach->addTrackExtended(trkExt);
          trkExt->setGroupTracks(groupToAttach);

          // set the min and max z extents
          float zminGroup = groupToAttach->getEdges()[0];
          float zmaxGroup = groupToAttach->getEdges()[1];

          float edges[2];
          edges[0] = zmin;
          if (zminGroup < zmin)
            edges[0] = zminGroup;
          edges[1] = zmax;
          if (zmaxGroup > zmax)
            edges[1] = zmaxGroup;
          groupToAttach->setEdges(edges);
          if (_debug >= 3) {
            int iopt = 3;
            PrintOutMerging(trkExt, trkToAttach, iopt);
          }
        } else {
          // create a new group of segments
          GroupTracks* newSegment = new GroupTracks(trkExt);
          trkExt->setGroupTracks(newSegment);
          streamlog_out(DEBUG2) << " *****************  AddNotCombinedTracks: Create new TPC Segment Group for track "
                                << trkExt << " id = " << trkExt->getTrack()->id() << std::endl;

          TPCSegments.push_back(newSegment);
          float edges[2];
          edges[0] = zmin;
          edges[1] = zmax;
          newSegment->setEdges(edges);
        }
      }
    }

    // At this stage all tpc segements will have been grouped.

    // Now try to combine the groups of TPC segments with the
    // reconstructed tracks which have already been combined into full tracks
    // containing both Si and TPC segments

    int nCombTrk = int(_trkImplVec.size());
    int nSegments = int(TPCSegments.size());

    //    std::cout << "Combined tracks = " << nCombTrk << std::endl;
    //    std::cout << "nSegments = " << nSegments << std::endl;

    // loop over all the TPC segment collections
    for (int iS = 0; iS < nSegments; ++iS) {
      GroupTracks* segments = TPCSegments[iS];
      TrackExtendedVec segVec = segments->getTrackExtendedVec();

      float zminTPCSeg = segments->getEdges()[0];
      float zmaxTPCSeg = segments->getEdges()[1];

      int nTrk = int(segVec.size());
      TrackExtended* CombTrkToAttach = NULL;
      TrackExtended* keyTrack = NULL;

      float deltaPtMin = _dPCutToMergeTPC;

      // search over the combined (good) tracks
      for (int iCTrk = 0; iCTrk < nCombTrk; ++iCTrk) {
        TrackExtended* combTrk = _trkImplVec[iCTrk];
        GroupTracks* groupComb = combTrk->getGroupTracks();

        bool consider = true;

        if (_forbidOverlapInZComb == 1) { // if overlap in Z of the two segments is forbidden
          float zminComb = groupComb->getEdges()[0];
          float zmaxComb = groupComb->getEdges()[1];
          consider = (zminTPCSeg > zmaxComb) || (zmaxTPCSeg < zminComb);
        }

        // if there are not overlaps in z, if _forbidOverlapInZComb is set above
        if (consider) {
          // loop over the TPC segments in the group
          for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
            TrackExtended* trk = segVec[iTrk];
            int iopt = 0;

            // test for compatibility
            float dPt = CompareTrk(trk, combTrk, _d0CutToMergeTPC, _z0CutToMergeTPC, iopt);
            float angleSignificance(0.);
            float significance =
                CompareTrkIII(trk, combTrk, _d0CutToMergeTPC, _z0CutToMergeTPC, iopt, angleSignificance);

            // check if this is a better match than any before
            if ((dPt < deltaPtMin || significance < 5)) {
              if (VetoMerge(trk, combTrk) == false) {
                // asign the track to be attached
                CombTrkToAttach = combTrk;
                keyTrack = trk;
                deltaPtMin = dPt;
              }
            } else {
              if (_debug >= 3) {
                GroupTracks* groupCur = combTrk->getGroupTracks();
                TrackExtended* dummySi = groupCur->getTrackExtendedVec()[0];
                int iopt_temp = 8;
                PrintOutMerging(trk, dummySi, iopt_temp);
              }
            }
          }
        } else {
          if (_debug >= 3) {
            for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
              TrackExtended* trk = segVec[iTrk];
              int iopt = 0;
              float dPt = CompareTrk(trk, combTrk, _d0CutToMergeTPC, _z0CutToMergeTPC, iopt);
              if (dPt > deltaPtMin) {
                GroupTracks* groupCur = combTrk->getGroupTracks();
                TrackExtended* dummySi = groupCur->getTrackExtendedVec()[0];
                int iopt_temp = 8;
                PrintOutMerging(trk, dummySi, iopt_temp);
              }
            }
          }
        }
      }

      if (CombTrkToAttach != NULL) { // attach TPC segment to existing Comb Track
        GroupTracks* groupToAttach = CombTrkToAttach->getGroupTracks();
        TrackExtended* SiCombTrk = groupToAttach->getTrackExtendedVec()[0];
        TrackExtended* TpcCombTrk = groupToAttach->getTrackExtendedVec()[1];

        if (_debug >= 3) {
          int iopt = 4;
          PrintOutMerging(keyTrack, SiCombTrk, iopt);
          iopt = 5;
          PrintOutMerging(keyTrack, TpcCombTrk, iopt);
        }

        for (int iTrk = 0; iTrk < nTrk; iTrk++) {
          TrackExtended* segmentTrack = segVec[iTrk];
          groupToAttach->addTrackExtended(segmentTrack);
          segmentTrack->setGroupTracks(groupToAttach);

          TrackerHitExtendedVec hitVec = segmentTrack->getTrackerHitExtendedVec();

          int nHitSeg = int(hitVec.size());

          for (int iHS = 0; iHS < nHitSeg; ++iHS) {
            // take the hit from the segment and attach it to CombTrkToAttach
            // flagging it as not to be used in the fit
            TrackerHitExtended* hitExt = hitVec[iHS];
            hitExt->setUsedInFit(false);
            CombTrkToAttach->addTrackerHitExtended(hitExt);
          }
        }
      } else {
        if (nTrk == 1) { // create a new group
          GroupTracks* newGrp = new GroupTracks();
          segVec[0]->setGroupTracks(newGrp);
          newGrp->addTrackExtended(segVec[0]);
          TrackerHitExtendedVec TpcHitVec = segVec[0]->getTrackerHitExtendedVec();
          int nTpcH = int(TpcHitVec.size());
          for (int iTpcH = 0; iTpcH < nTpcH; ++iTpcH) {
            TpcHitVec[iTpcH]->setUsedInFit(true);
          }
          _trkImplVec.push_back(segVec[0]);
          _allNonCombinedTPCTracks.push_back(segVec[0]);
        } else { // several segments

          float zMin = 1.0e+20;
          TrackExtended* chosenTrack = NULL;

          // loop over the segments
          for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
            TrackExtended* segment = segVec[iTrk];

            // get the lcio track which is behind this segemnt
            Track* track = segment->getTrack();
            TrackerHitVec hitVec = track->getTrackerHits();

            streamlog_out(DEBUG1) << "Group of orphaned TPC tracks: trying track " << track->id() << std::endl;

            int nHits = int(hitVec.size());

            // loop over it's hits
            for (int iH = 0; iH < nHits; ++iH) {
              float zPosi = fabs(hitVec[iH]->getPosition()[2]);

              // if this segment has the hit closest to the IP so far
              if (zPosi < zMin) {
                // take this as the chosen track and break
                chosenTrack = segment;
                zMin = zPosi;
                break;
              }
            }
          }

          if (chosenTrack != NULL) { // can't really ever be null.

            streamlog_out(DEBUG2) << "Group of orphaned TPC tracks: chosen track taken as "
                                  << chosenTrack->getTrack()->id() << std::endl;

            // create a new group of tracks
            GroupTracks* newGroup = new GroupTracks();

            // first add the chosen track
            chosenTrack->setGroupTracks(newGroup);
            newGroup->addTrackExtended(chosenTrack);

            // loop over the segments ...
            for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
              TrackExtended* segment = segVec[iTrk];

              TrackerHitExtendedVec hitVecS = segment->getTrackerHitExtendedVec();
              int nHitS = int(hitVecS.size());

              // loop over the hits for the current segment
              for (int iH = 0; iH < nHitS; ++iH) {
                TrackerHitExtended* trkHitExt = hitVecS[iH];

                if (segment != chosenTrack) { // ... then don't add the hits to the fit
                  // set the relation between group and track
                  segment->setGroupTracks(newGroup);
                  newGroup->addTrackExtended(segment);
                  trkHitExt->setUsedInFit(false);
                  chosenTrack->addTrackerHitExtended(trkHitExt);
                } else {
                  trkHitExt->setUsedInFit(true);
                }
              }
            }
            _allNonCombinedTPCTracks.push_back(chosenTrack);
            _trkImplVec.push_back(chosenTrack);
          }
        }
      }
    }
    for (int iS = 0; iS < nSegments; ++iS) {
      GroupTracks* segments = TPCSegments[iS];
      delete segments;
    }
    TPCSegments.clear();
  } else { // adding all TPC segments to the list of tracks (track splitting is allowed)
    for (int i = 0; i < nTPCTrk; ++i) {
      TrackExtended* trkExt = _allTPCTracks[i];
      Track* track = trkExt->getTrack();
      GroupTracks* group = trkExt->getGroupTracks();

      if (group == NULL) {
        streamlog_out(DEBUG2) << " *****************  AddNotCombinedTracks: _mergeTPCSegments = " << _mergeTPCSegments
                              << " : Add non combined TPC track " << trkExt << " id = " << track->id() << std::endl;

        TrackerHitExtendedVec hitVec = trkExt->getTrackerHitExtendedVec();
        int nHTPC = int(hitVec.size());

        for (int iHTPC = 0; iHTPC < nHTPC; ++iHTPC) {
          hitVec[iHTPC]->setUsedInFit(true);
        }
        _trkImplVec.push_back(trkExt);
        _allNonCombinedTPCTracks.push_back(trkExt);

        GroupTracks* newGrp = new GroupTracks();
        newGrp->addTrackExtended(trkExt);
        trkExt->setGroupTracks(newGrp);
      }
    }
  }

  for (int i = 0; i < nSiTrk; ++i) { // adding left-over Si segments to the list of tracks
    TrackExtended* trkExt = _allSiTracks[i];
    Track* track = trkExt->getTrack();
    GroupTracks* group = trkExt->getGroupTracks();

    if (group == NULL) {
      streamlog_out(DEBUG2) << " *****************  AddNotCombinedTracks: Add non combined Silicon Track : " << trkExt
                            << " id = " << track->id() << std::endl;

      TrackerHitExtendedVec hitVec = trkExt->getTrackerHitExtendedVec();
      int nHSi = int(hitVec.size());

      for (int iHSi = 0; iHSi < nHSi; ++iHSi) {
        hitVec[iHSi]->setUsedInFit(true);
      }

      _trkImplVec.push_back(trkExt);
      _allNonCombinedSiTracks.push_back(trkExt);

      GroupTracks* newGrp = new GroupTracks();
      newGrp->addTrackExtended(trkExt);
      trkExt->setGroupTracks(newGrp);
    }
  }
}

void FullLDCTracking_MarlinTrk::CheckTracks() {
  for (unsigned int i = 0; i < _trkImplVec.size(); i++) {
    TrackExtended* first = _trkImplVec[i];
    if (first == NULL)
      continue;
    float d0First = first->getD0();
    float z0First = first->getZ0();
    float omegaFirst = first->getOmega();
    float tanLFirst = first->getTanLambda();
    float phiFirst = first->getPhi();
    HelixClass helixFirst;
    helixFirst.Initialize_Canonical(phiFirst, d0First, z0First, omegaFirst, tanLFirst, _bField);
    float momFirst[3];
    momFirst[0] = helixFirst.getMomentum()[0];
    momFirst[1] = helixFirst.getMomentum()[1];
    momFirst[2] = helixFirst.getMomentum()[2];
    float pFirst = sqrt(momFirst[0] * momFirst[0] + momFirst[1] * momFirst[1] + momFirst[2] * momFirst[2]);
    if (std::isnan(pFirst))
      continue;
    TrackerHitExtendedVec firstHitVec = first->getTrackerHitExtendedVec();
    if (firstHitVec.size() < 1)
      continue;

    for (unsigned int j = i + 1; j < _trkImplVec.size(); j++) {
      TrackExtended* second = _trkImplVec[j];
      if (second == NULL)
        continue;
      float d0Second = second->getD0();
      float z0Second = second->getZ0();
      float omegaSecond = second->getOmega();
      float tanLSecond = second->getTanLambda();
      float phiSecond = second->getPhi();
      HelixClass helixSecond;
      helixSecond.Initialize_Canonical(phiSecond, d0Second, z0Second, omegaSecond, tanLSecond, _bField);
      float momSecond[3];
      momSecond[0] = helixSecond.getMomentum()[0];
      momSecond[1] = helixSecond.getMomentum()[1];
      momSecond[2] = helixSecond.getMomentum()[2];
      float pSecond = sqrt(momSecond[0] * momSecond[0] + momSecond[1] * momSecond[1] + momSecond[2] * momSecond[2]);
      if (std::isnan(pSecond))
        continue;
      TrackerHitExtendedVec secondHitVec = second->getTrackerHitExtendedVec();
      if (secondHitVec.size() < 1)
        continue;
      if (firstHitVec.size() + secondHitVec.size() < 10)
        continue;

      float pdot =
          (momFirst[0] * momSecond[0] + momFirst[1] * momSecond[1] + momFirst[2] * momSecond[2]) / pFirst / pSecond;
      if (pdot < 0.999)
        continue;
      // const float sigmaPOverPFirst  = sqrt(first->getCovMatrix()[5])/fabs(omegaFirst);
      // const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5])/fabs(omegaSecond);
      // const float deltaP = fabs(pFirst-pSecond);
      // const float sigmaDeltaP =
      // sqrt(pFirst*sigmaPOverPFirst*pFirst*sigmaPOverPFirst+pSecond*sigmaPOverPSecond*pSecond*sigmaPOverPSecond);
      //      const float significance = deltaP/sigmaDeltaP;

      TrackExtended* combinedTrack =
          CombineTracks(first, second, _maxAllowedPercentageOfOutliersForTrackCombination, true);
      if (combinedTrack != NULL) {
        const int minHits = std::min(firstHitVec.size(), secondHitVec.size());
        const int maxHits = std::max(firstHitVec.size(), secondHitVec.size());

        if (combinedTrack->getNDF() <= 2 * maxHits + minHits - 5) {
          delete combinedTrack->getGroupTracks();
          delete combinedTrack;
          continue;
        }

        float d0 = combinedTrack->getD0();
        float z0 = combinedTrack->getZ0();
        float omega = combinedTrack->getOmega();
        float tanL = combinedTrack->getTanLambda();
        float phi = combinedTrack->getPhi();

        HelixClass helix;
        helix.Initialize_Canonical(phi, d0, z0, omega, tanL, _bField);
        // float mom[3];
        // mom[0]  = helix.getMomentum()[0];
        // mom[1]  = helix.getMomentum()[1];
        // mom[2]  = helix.getMomentum()[2];
        // float p = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
        // float chi2Sig =  (combinedTrack->getChi2() - combinedTrack->getNDF());
        // chi2Sig = chi2Sig/sqrt(combinedTrack->getNDF()*2);

        int nTpcFirst(0);
        int nUsedFirst(0);
        for (unsigned int ihit = 0; ihit < firstHitVec.size(); ihit++) {
          if (getDetectorID(firstHitVec[ihit]->getTrackerHit()) == lcio::ILDDetID::TPC)
            nTpcFirst++;

          if (firstHitVec[ihit]->getUsedInFit() == true)
            nUsedFirst++;
        }

        int nTpcSecond(0);
        int nUsedSecond(0);
        for (unsigned int ihit = 0; ihit < secondHitVec.size(); ihit++) {
          if (getDetectorID(secondHitVec[ihit]->getTrackerHit()) == lcio::ILDDetID::TPC)
            ++nTpcSecond;
          if (secondHitVec[ihit]->getUsedInFit() == true)
            ++nUsedSecond;
        }
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      }
    }
  }
}

/*

 compare the following:

 i)   delta omega < 2 * _dOmegaForMerging
 ii)  delta d0    < d0Cut
 iii) delta z0    < z0Cut

 if above three cuts are passed return:

 the angle between the two tracks (return by reference)
   and
 the difference in omega

 */

float FullLDCTracking_MarlinTrk::CompareTrkII(TrackExtended* first, TrackExtended* second, float d0Cut, float z0Cut,
                                              int /*iopt*/, float& Angle) {
  float result = 1.0e+20;
  Angle = 1.0e+20;
  float omegaFirst = first->getOmega();
  float omegaSecond = second->getOmega();
  float deltaOmega = fabs((omegaFirst - omegaSecond) / omegaSecond);
  if (deltaOmega > 2 * _dOmegaForMerging)
    return result;

  float d0First = first->getD0();
  float z0First = first->getZ0();
  float d0Second = second->getD0();
  float z0Second = second->getZ0();

  bool isCloseInIP = (fabs(d0First - d0Second) < d0Cut);
  isCloseInIP = isCloseInIP && (fabs(z0Second - z0First) < z0Cut);

  if (!isCloseInIP)
    return result;

  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  float qFirst = PIOVER2 - atan(tanLFirst);
  float qSecond = PIOVER2 - atan(tanLSecond);

  Angle = (cos(phiFirst) * cos(phiSecond) + sin(phiFirst) * sin(phiSecond)) * sin(qFirst) * sin(qSecond) +
          cos(qFirst) * cos(qSecond);
  Angle = acos(Angle);

  result = deltaOmega;

  return result;
}

/*

 compare the following:

 i)   significance d0    < d0Cut
 ii)  significance z0    < z0Cut
 iii) pdot > 0.999

 if above three cuts are passed return:

 the significance dAngle (return by reference)
 and
 the significance dOmega

 */

float FullLDCTracking_MarlinTrk::CompareTrkIII(TrackExtended* first, TrackExtended* second, float d0Cut, float z0Cut,
                                               int /*iopt*/, float& AngleSignificance) {
  float result = 1.0e+20;

  float d0First = first->getD0();
  float z0First = first->getZ0();
  float omegaFirst = first->getOmega();
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();
  float qFirst = PIOVER2 - atan(tanLFirst);

  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  float omegaSecond = second->getOmega();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();
  float qSecond = PIOVER2 - atan(tanLSecond);

  // MB 2010 03
  float d0ErrFirst = sqrt(first->getCovMatrix()[0]);
  float z0ErrFirst = sqrt(first->getCovMatrix()[9]);
  //  float omegaErrFirst = sqrt(first->getCovMatrix()[5]);
  float phiErrFirst = sqrt(first->getCovMatrix()[2]);
  float qErrFirst = sqrt(cos(qFirst) * cos(qFirst) * first->getCovMatrix()[14]);
  // MB END
  // MB 2010 03
  float d0ErrSecond = sqrt(second->getCovMatrix()[0]);
  float z0ErrSecond = sqrt(second->getCovMatrix()[9]);
  //  float omegaErrSecond = sqrt(second->getCovMatrix()[5]);
  float phiErrSecond = sqrt(second->getCovMatrix()[2]);
  float qErrSecond = sqrt(cos(qSecond) * cos(qSecond) * second->getCovMatrix()[14]);
  // MB END

  //  bool isCloseInIP = (fabs(d0First-d0Second)<d0Cut);
  // isCloseInIP = isCloseInIP && (fabs(z0Second-z0First)<z0Cut);

  // MB 2010 03
  bool isCloseInIP = (fabs(d0First - d0Second) / sqrt(d0ErrFirst * d0ErrFirst + d0ErrSecond * d0ErrSecond) < d0Cut);
  isCloseInIP =
      isCloseInIP && (fabs(z0Second - z0First) / sqrt(z0ErrFirst * z0ErrFirst + z0ErrSecond * z0ErrSecond) < z0Cut);

  if (!isCloseInIP)
    return result;

  float Angle = (cos(phiFirst) * cos(phiSecond) + sin(phiFirst) * sin(phiSecond)) * sin(qFirst) * sin(qSecond) +
                cos(qFirst) * cos(qSecond);

  HelixClass helixFirst;
  helixFirst.Initialize_Canonical(phiFirst, d0First, z0First, omegaFirst, tanLFirst, _bField);
  HelixClass helixSecond;
  helixSecond.Initialize_Canonical(phiSecond, d0Second, z0Second, omegaSecond, tanLSecond, _bField);

  float pFirst[3];
  float pSecond[3];
  float momFirst = 0;
  float momSecond = 0;

  for (int iC = 0; iC < 3; ++iC) {
    pFirst[iC] = helixFirst.getMomentum()[iC];
    pSecond[iC] = helixSecond.getMomentum()[iC];
    momFirst += pFirst[iC] * pFirst[iC];
    momSecond += pSecond[iC] * pSecond[iC];
  }
  momFirst = sqrt(momFirst);
  momSecond = sqrt(momSecond);

  float pdot = (pFirst[0] * pSecond[0] + pFirst[1] * pSecond[1] + pFirst[2] * pSecond[2]) / momFirst / momSecond;
  if (pdot < 0.999)
    return result;

  const float sigmaPOverPFirst = sqrt(first->getCovMatrix()[5]) / fabs(omegaFirst);
  const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5]) / fabs(omegaSecond);
  const float deltaP = fabs(momFirst - momSecond);
  const float sigmaPFirst = momFirst * sigmaPOverPFirst;
  const float sigmaPSecond = momSecond * sigmaPOverPSecond;
  const float sigmaDeltaP = sqrt(sigmaPFirst * sigmaPFirst + sigmaPSecond * sigmaPSecond);
  const float significance = deltaP / sigmaDeltaP;

  // MB 2010 03
  float errorAngle =
      sin(phiFirst) * sin(phiFirst) * phiErrFirst * phiErrFirst * cos(phiSecond) * cos(phiSecond) +
      sin(phiSecond) * sin(phiSecond) * phiErrSecond * phiErrSecond * cos(phiFirst) * cos(phiFirst) +
      sin(qFirst) * sin(qFirst) * qErrFirst * qErrFirst * cos(qSecond) * cos(qSecond) +
      sin(qSecond) * sin(qSecond) * qErrSecond * qErrSecond * cos(qFirst) * cos(qFirst) +
      cos(phiFirst) * cos(phiFirst) * phiErrFirst * phiErrFirst * (sin(phiSecond) * sin(qFirst) * sin(qSecond)) *
          (sin(phiSecond) * sin(qFirst) * sin(qSecond)) +
      cos(phiSecond) * cos(phiSecond) * phiErrSecond * phiErrSecond * (sin(phiFirst) * sin(qFirst) * sin(qSecond)) *
          (sin(phiFirst) * sin(qFirst) * sin(qSecond)) +
      cos(qFirst) * cos(qFirst) * qErrFirst * qErrFirst * (sin(phiFirst) * sin(phiSecond) * sin(qSecond)) *
          (sin(phiFirst) * sin(phiSecond) * sin(qSecond)) +
      cos(qSecond) * cos(qSecond) * qErrSecond * qErrSecond * (sin(phiFirst) * sin(phiSecond) * sin(qFirst)) *
          (sin(phiFirst) * sin(phiSecond) * sin(qFirst));

  if (Angle < 1.) {
    errorAngle = sqrt(1. / (1. - Angle * Angle) * errorAngle);
  } else {
    errorAngle = sqrt(errorAngle);
  }

  if (errorAngle < 1.e-6)
    errorAngle = 1.e-6;

  AngleSignificance = fabs(acos(Angle) / errorAngle);

  return significance;
}

/*

 compare the following:

 i)   delta d0    < d0Cut  and optionally (d0_1 + d0_2) < d0cut
 ii)  delta z0    < z0Cut

 if above two cuts are passed then:

 if ( (ptFirst<_PtCutToMergeTPC) && (ptSecond<_PtCutToMergeTPC) ) {

     then check the difference in momentum

 else

     check for cases where PatRec splits non-looping TPC tracks
     look for two tracks where total tpc hits are not more than total number
     of pad rows and that the hits on one track are close to the helix of the
     other track.

     Check that the angular and momentum difference meets the cuts for either hight or low  pt
     also check if the angle is very small between the tracks and that the significance of the difference in pt is less
 than 10


 ....

 note currently this is only used in AddNotCombinedTracks

 */

float FullLDCTracking_MarlinTrk::CompareTrk(TrackExtended* first, TrackExtended* second, float d0Cut, float z0Cut,
                                            int iopt) {
  float result = 1.0e+20;

  float d0First = first->getD0();
  float z0First = first->getZ0();
  float omegaFirst = first->getOmega();
  float tanLFirst = first->getTanLambda();
  float phiFirst = first->getPhi();

  float d0Second = second->getD0();
  float z0Second = second->getZ0();
  float omegaSecond = second->getOmega();
  float tanLSecond = second->getTanLambda();
  float phiSecond = second->getPhi();

  bool isCloseInIP = (fabs(d0First - d0Second) < d0Cut);

  if (iopt > 0)
    isCloseInIP = isCloseInIP || (fabs(d0First + d0Second) < d0Cut);

  isCloseInIP = isCloseInIP && (fabs(z0Second - z0First) < z0Cut);

  HelixClass helixFirst;
  helixFirst.Initialize_Canonical(phiFirst, d0First, z0First, omegaFirst, tanLFirst, _bField);
  HelixClass helixSecond;
  helixSecond.Initialize_Canonical(phiSecond, d0Second, z0Second, omegaSecond, tanLSecond, _bField);

  float pFirst[3];
  float pSecond[3];
  float dPminus[3];
  float dPplus[3];
  float momFirst = 0;
  float momSecond = 0;
  float momMinus = 0;
  float momPlus = 0;

  if (isCloseInIP) {
    for (int iC = 0; iC < 3; ++iC) {
      pFirst[iC] = helixFirst.getMomentum()[iC];
      pSecond[iC] = helixSecond.getMomentum()[iC];

      momFirst += pFirst[iC] * pFirst[iC];
      momSecond += pSecond[iC] * pSecond[iC];

      dPminus[iC] = pFirst[iC] - pSecond[iC];
      dPplus[iC] = pFirst[iC] + pSecond[iC];

      momMinus += dPminus[iC] * dPminus[iC];
      momPlus += dPplus[iC] * dPplus[iC];
    }
    momFirst = sqrt(momFirst);
    momSecond = sqrt(momSecond);

    float ptFirst = sqrt(pFirst[0] * pFirst[0] + pFirst[1] * pFirst[1]);
    float ptSecond = sqrt(pSecond[0] * pSecond[0] + pSecond[1] * pSecond[1]);

    // if both track's pt are lower than _PtCutToMergeTPC
    if ((ptFirst < _PtCutToMergeTPC) && (ptSecond < _PtCutToMergeTPC)) {
      momMinus = sqrt(momMinus);
      momPlus = sqrt(momPlus);

      float nom = momMinus;

      // get the smaller difference for the nominator
      if (momPlus < nom && iopt > 0)
        nom = momPlus;

      float den = momFirst;

      // get the smallest momentum for the denominator
      if (momSecond < momFirst)
        den = momSecond;

      result = nom / den;
    }

    else {
      // check for cases where PatRec splits non-looping TPC tracks
      // look for two tracks where total tpc hits are not more than total number
      // of pad rows and that the hits on one track are close to the helix of the
      // other track

      float dpOverP = 2.0 * fabs(momFirst - momSecond) / (momFirst + momSecond);
      const float pdot =
          (pFirst[0] * pSecond[0] + pFirst[1] * pSecond[1] + pFirst[2] * pSecond[2]) / momFirst / momSecond;
      const float sigmaPOverPFirst = sqrt(first->getCovMatrix()[5]) / fabs(omegaFirst);
      const float sigmaPOverPSecond = sqrt(second->getCovMatrix()[5]) / fabs(omegaSecond);
      const float deltaP = fabs(momFirst - momSecond);
      const float sigmaPFirst = momFirst * sigmaPOverPFirst;
      const float sigmaPSecond = momSecond * sigmaPOverPSecond;
      const float sigmaDeltaP = sqrt(sigmaPFirst * sigmaPFirst + sigmaPSecond * sigmaPSecond);
      const float significance = deltaP / sigmaDeltaP;

      // compare angle between the two vectors (cos theta) and their momentum
      //     if( ( pdot>0.99 && dpOverP<0.01 ) || ( pdot>0.998 && dpOverP<0.25 ) ){
      if ((pdot > _cosThetaCutHighPtMerge && dpOverP < _momDiffCutHighPtMerge) ||
          (pdot > _cosThetaCutSoftHighPtMerge && dpOverP < _momDiffCutSoftHighPtMerge) ||
          (pdot > 0.9999 && significance < 10)) {
        int nTrkGrpFirst = 0;
        int nTrkGrpSecond = 0;
        TrackerHitVec hitvecFirst;
        TrackerHitVec hitvecSecond;
        GroupTracks* groupFirst = first->getGroupTracks();
        GroupTracks* groupSecond = second->getGroupTracks();

        // does the first track belong to a group ...
        // if it does then get all the hits from the group
        if (groupFirst != NULL) {
          TrackExtendedVec tracksInGroupFirst = groupFirst->getTrackExtendedVec();
          nTrkGrpFirst = int(tracksInGroupFirst.size());

          for (int iTrkGrp = 0; iTrkGrp < nTrkGrpFirst; ++iTrkGrp) {
            TrackExtended* trkGrp = tracksInGroupFirst[iTrkGrp];
            TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();

            for (unsigned int i = 0; i < hitVec.size(); ++i) {
              hitvecFirst.push_back(hitVec[i]->getTrackerHit());
            }
          }
        }

        // does the second track belong to a group ...
        // if it does then get all the hits from the group
        if (groupSecond != NULL) {
          TrackExtendedVec tracksInGroupSecond = groupSecond->getTrackExtendedVec();
          nTrkGrpSecond = int(tracksInGroupSecond.size());

          for (int iTrkGrp = 0; iTrkGrp < nTrkGrpSecond; ++iTrkGrp) {
            TrackExtended* trkGrp = tracksInGroupSecond[iTrkGrp];
            TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();

            for (unsigned int i = 0; i < hitVec.size(); ++i) {
              hitvecSecond.push_back(hitVec[i]->getTrackerHit());
            }
          }
        }

        // for non-looping tracks
        int nhitsFirst = (int)hitvecFirst.size();
        int nhitsSecond = (int)hitvecSecond.size();
        int ntpcFirst = 0;
        int ntpcSecond = 0;
        float hitxyz[3];
        float dist[3];
        float maxdistFirst = 0;
        float maxdistSecond = 0;
        int ncloseFirst = 0;
        int ncloseSecond = 0;
        float zminFirst = 99999;
        float zminSecond = 99999;
        float zmaxFirst = -99999;
        float zmaxSecond = -99999;

        for (int ih = 0; ih < nhitsFirst; ++ih) {
          float x = (float)hitvecFirst[ih]->getPosition()[0];
          float y = (float)hitvecFirst[ih]->getPosition()[1];
          float z = (float)hitvecFirst[ih]->getPosition()[2];

          if (fabs(z) < zminFirst)
            zminFirst = fabs(z);
          if (fabs(z) > zmaxFirst)
            zmaxFirst = fabs(z);

          float r = sqrt(x * x + y * y);

          // count the number of hits in the TPC for the first Track
          if (r > _tpc_inner_r)
            ntpcFirst++;

          hitxyz[0] = x;
          hitxyz[1] = y;
          hitxyz[2] = z;
          helixSecond.getDistanceToPoint(hitxyz, dist);

          // compare 3D distance between hit and extrapolation
          if (dist[2] > maxdistFirst)
            maxdistFirst = dist[2];

          // count the number of hits from the first track which are closer than _hitDistanceCutHighPtMerge to the
          // second track
          if (dist[2] < _hitDistanceCutHighPtMerge)
            ncloseFirst++;
        }

        for (int ih = 0; ih < nhitsSecond; ++ih) {
          float x = (float)hitvecSecond[ih]->getPosition()[0];
          float y = (float)hitvecSecond[ih]->getPosition()[1];
          float z = (float)hitvecSecond[ih]->getPosition()[2];

          if (fabs(z) < zminSecond)
            zminSecond = fabs(z);
          if (fabs(z) > zmaxSecond)
            zmaxSecond = fabs(z);

          float r = sqrt(x * x + y * y);

          // count the number of hits in the TPC for the second Track
          if (r > _tpc_inner_r)
            ntpcSecond++;

          hitxyz[0] = x;
          hitxyz[1] = y;
          hitxyz[2] = z;
          helixFirst.getDistanceToPoint(hitxyz, dist);

          // compare 3D distance between hit and extrapolation
          if (dist[2] > maxdistSecond)
            maxdistSecond = dist[2];

          // count the number of hits from the second track which are closer than _hitDistanceCutHighPtMerge to the
          // first track
          if (dist[2] < _hitDistanceCutHighPtMerge)
            ncloseSecond++;
        }

        // calculate the fraction of hits which are closer than _hitDistanceCutHighPtMerge
        float fcloseFirst = (float)ncloseFirst / (float)nhitsFirst;
        float fcloseSecond = (float)ncloseSecond / (float)nhitsSecond;

        bool split = false;
        // std::cout << "Momenta = " << momFirst << " " << momSecond << std::endl;
        // std::cout << "MaxDist = " << maxdistSecond << " " << maxdistFirst << " " << _maxHitDistanceCutHighPtMerge <<
        // std::endl; std::cout << "close   = " << fcloseSecond << " " << fcloseFirst << " " <<
        // _maxFractionOfOutliersCutHighPtMerge << std::endl; std::cout << "ntpc    = " << ntpcFirst << " " <<
        // ntpcSecond << " " << _tpc_pad_height+10 << std::endl; std::cout << "pdot    = " << pdot << " significance "
        // << significance << std::endl;

        // SJA:FIXME: try to fit the two tracks, without checking the number of hits which are close !!!!!!!!!!!
        TrackExtended* combinedTrack =
            CombineTracks(first, second, _maxAllowedPercentageOfOutliersForTrackCombination, true);

        // check that no more than 5 hits have been discared in the fit of the combined track
        if (combinedTrack != NULL) {
          // std::cout << "CombinedTrack " << combinedTrack->getNDF() << " c.f. " << first->getNDF()+second->getNDF()+5
          // << std::endl;
          if (combinedTrack->getNDF() + 10 > first->getNDF() + second->getNDF() + 5) {
            split = true;
            dpOverP = 0;
            // std::cout << " Forcing MERGE " << std::endl;
          }
          delete combinedTrack->getGroupTracks();
          delete combinedTrack;
        } else {
          // std::cout << "Could not combine track " << std::endl;
          //  if close in pt and the fraction of hits matching the helix extraolations is greater than
          //  _maxFractionOfOutliersCutHighPtMerge
          if (significance < 5 && fcloseFirst > _maxFractionOfOutliersCutHighPtMerge) {
            split = true;
            dpOverP = 0;
            //      int overlap = SegmentRadialOverlap(first,second);
            // std::cout << " Forcing MERGE " << overlap << std::endl;
          }
        }

        // criteria for split track
        // old criterion
        // check the maximum deviation of the hits from the helix extrapolations, and fraction of hits matching the
        // extrapolation based on a distance cut
        if (maxdistSecond < _maxHitDistanceCutHighPtMerge && maxdistFirst < _maxHitDistanceCutHighPtMerge &&
            (fcloseSecond > _maxFractionOfOutliersCutHighPtMerge ||
             fcloseFirst > _maxFractionOfOutliersCutHighPtMerge)) {
          //           &&
          //           ntpcFirst+ntpcSecond < _tpc_pad_height+10.) { // SJA:FIXME: this must be a bug, first it is a
          //           double, and second this value is never initialised ...

          split = true;
        }

        if (split) {
          result = dpOverP;
        }
      }
    }
  }

  return result;
}

void FullLDCTracking_MarlinTrk::AddNotAssignedHits() {
  // currently any non-combined Silicon or TPC tracks are added to the list of tracks meaning their hits will not be
  // available to be picked up. it might be preferable to drop these tracks, at least for track in Silicon with very bad
  // Chi2 probabilities, and allow their hits to be re-used here ...

  if (_assignSETHits > 0) { // Assignment of SET Hits

    streamlog_out(DEBUG4) << "Assign SET hits *********************************" << std::endl;

    // Creating helix extrapolations
    CreateExtrapolations();

    int nSETHits = _allSETHits.size();
    std::vector<TrackerHitExtendedVec> SETHits;

    SETHits.resize(_nLayersSET);

    for (int iSET = 0; iSET < nSETHits; ++iSET) {
      TrackerHitExtended* trkHitExt = _allSETHits[iSET];
      TrackerHit* trkHit = trkHitExt->getTrackerHit();
      int layer = getLayerID(trkHit);
      if (layer >= 0 && (unsigned)layer < _nLayersSET)
        SETHits[layer].push_back(trkHitExt);
    }
    for (unsigned iL = 0; iL < _nLayersSET; ++iL) { // loop over SET layers
      TrackerHitExtendedVec hitVec = SETHits[iL];
      int refit = 1;
      if (hitVec.empty() == false)
        AssignOuterHitsToTracks(hitVec, _distCutForSETHits, refit);
    }
  }

  CleanUpExtrapolations();

  if (_assignSITHits > 0) { // Treatment of left-over SIT hits

    streamlog_out(DEBUG4) << "Assign SIT hits *********************************" << std::endl;

    std::vector<TrackerHitExtendedVec> nonAssignedSITHits;
    nonAssignedSITHits.resize(_nLayersSIT);

    int nSITHits = int(_allSITHits.size());

    // loop over all SIT hits ...
    for (int iH = 0; iH < nSITHits; ++iH) {
      TrackerHitExtended* trkHitExt = _allSITHits[iH];
      TrackExtended* trkExt = trkHitExt->getTrackExtended();

      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
        TrackerHit* trkHit = trkHitExt->getTrackerHit();

        int layer = getLayerID(trkHit);

        if (layer >= 0 && (unsigned)layer < _nLayersSIT) {
          nonAssignedSITHits[layer].push_back(trkHitExt);
        }
      }
    }

    for (int iL = _nLayersSIT - 1; iL >= 0; --iL) { // reverse loop over layers in Si

      TrackerHitExtendedVec hitVec = nonAssignedSITHits[iL];

      if (hitVec.empty() == false) {
        streamlog_out(DEBUG3) << "AddNotAssignedHits : Try to assign hits from layer " << iL
                              << " : Number of hits = " << hitVec.size() << std::endl;
        AssignSiHitsToTracks(hitVec, _distCutForSITHits);
      }
    }
  }

  if (_assignFTDHits > 0) { // Treatment of left-over FTD hits

    streamlog_out(DEBUG4) << "Assign FTD hits *********************************" << std::endl;

    std::vector<TrackerHitExtendedVec> nonAssignedFTDHits;
    nonAssignedFTDHits.resize(_nLayersFTD);

    int nFTDHits = int(_allFTDHits.size());

    // loop over all FTD hits ...
    for (int iH = 0; iH < nFTDHits; ++iH) {
      TrackerHitExtended* trkHitExt = _allFTDHits[iH];
      TrackExtended* trkExt = trkHitExt->getTrackExtended();

      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
        TrackerHit* trkHit = trkHitExt->getTrackerHit();

        // get the layer number
        int layer = getLayerID(trkHit);
        int petalIndex = getModuleID(trkHit);

        if (_petalBasedFTDWithOverlaps == true) {
          // as we are dealing with staggered petals we will use 2*nlayers in each directions +/- z
          // the layers will follow the even odd numbering of the petals
          if (petalIndex % 2 == 0) {
            layer = 2 * layer;
          } else {
            layer = 2 * layer + 1;
          }
        }

        if (layer >= 0 && layer < (int)_nLayersFTD)
          nonAssignedFTDHits[layer].push_back(trkHitExt);
      }
    }
    for (int iL = _nLayersFTD - 1; iL >= 0; --iL) {
      if (nonAssignedFTDHits[iL].empty() == false) {
        TrackerHitExtendedVec hitVec = nonAssignedFTDHits[iL];
        streamlog_out(DEBUG3) << "AddNotAssignedHits : Try to assign hits from layer " << iL
                              << " : Number of hits = " << hitVec.size() << std::endl;
        AssignSiHitsToTracks(hitVec, _distCutForFTDHits);
      }
    }
  }

  if (_assignVTXHits > 0) { // Treatment of left-over VTX hits

    streamlog_out(DEBUG4) << "Assign VXD hits *********************************" << std::endl;

    std::vector<TrackerHitExtendedVec> nonAssignedVTXHits;
    nonAssignedVTXHits.resize(_nLayersVTX);

    int nVTXHits = int(_allVTXHits.size());

    // loop over all VXD hits ...
    for (int iH = 0; iH < nVTXHits; ++iH) {
      TrackerHitExtended* trkHitExt = _allVTXHits[iH];
      TrackExtended* trkExt = trkHitExt->getTrackExtended();

      // check if this hit has not already been assigned to a track
      if (trkExt == NULL) {
        TrackerHit* trkHit = trkHitExt->getTrackerHit();

        int layer = getLayerID(trkHit);

        if (layer >= 0 && layer < (int)_nLayersVTX)
          nonAssignedVTXHits[layer].push_back(trkHitExt);
      }
    }
    for (int iL = _nLayersVTX - 1; iL >= 0; --iL) {
      if (nonAssignedVTXHits[iL].empty() == false) {
        TrackerHitExtendedVec hitVec = nonAssignedVTXHits[iL];
        streamlog_out(DEBUG3) << "AddNotAssignedHits : Try to assign hits from layer " << iL
                              << " : Number of hits = " << hitVec.size() << std::endl;
        AssignSiHitsToTracks(hitVec, _distCutForVTXHits);
      }
    }
  }

  streamlog_out(DEBUG4) << "Assign TPC hits *********************************" << std::endl;

  if (_assignTPCHits) { // Treatment of left-over TPC hits
    TrackerHitExtendedVec nonAssignedTPCHits;
    int nTPCHits = int(_allTPCHits.size());
    for (int iH = 0; iH < nTPCHits; ++iH) {
      TrackerHitExtended* trkHitExt = _allTPCHits[iH];
      TrackExtended* trkExt = trkHitExt->getTrackExtended();
      if (trkExt == NULL) {
        nonAssignedTPCHits.push_back(trkHitExt);
      }
    }
    streamlog_out(DEBUG4) << "AddNotAssignedHits : Number of Non Assigned TPC hits = " << nonAssignedTPCHits.size()
                          << std::endl;
    AssignTPCHitsToTracks(nonAssignedTPCHits, _distCutForTPCHits);
  }
}

void FullLDCTracking_MarlinTrk::CreateExtrapolations() {
  _trackExtrapolatedHelix.clear();

  int nTrk = int(_trkImplVec.size());
  for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
    TrackExtended* trk = _trkImplVec[iTrk];
    HelixClass* helix = GetExtrapolationHelix(trk);
    _trackExtrapolatedHelix[trk] = helix;
  }
}

void FullLDCTracking_MarlinTrk::CleanUpExtrapolations() {
  int nTrk = int(_trkImplVec.size());
  for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
    TrackExtended* trk = _trkImplVec[iTrk];
    HelixClass* helix = _trackExtrapolatedHelix[trk];
    delete helix;
  }
  _trackExtrapolatedHelix.clear();
}

void FullLDCTracking_MarlinTrk::AssignOuterHitsToTracks(TrackerHitExtendedVec hitVec, float dcut, int refit) {
  streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AssignOuterHitsToTracks dcut = " << dcut << std::endl;

  // get the number of hits to try, and the number of final tracks to which the tracks will be attached
  int nHits = int(hitVec.size());
  int nTrk = int(_trkImplVec.size());

  // create maps to record which tracks and tracker hits are flagged for assignment
  std::map<TrackExtended*, bool> flagTrack;
  std::map<TrackerHitExtended*, bool> flagHit;

  // vector to hold the matchups and the distance of closest approach.
  TrackHitPairVec pairs;

  flagTrack.clear();
  flagHit.clear();
  pairs.clear();

  // loop over all hits under consideration ...
  for (int iH = 0; iH < nHits; ++iH) {
    float pos[3];

    TrackerHitExtended* trkHitExt = hitVec[iH];
    TrackerHit* hit = trkHitExt->getTrackerHit();

    for (int ip = 0; ip < 3; ++ip) {
      pos[ip] = float(hit->getPosition()[ip]);
    }

    // loop over all tracks in same z-half ...
    for (int iT = 0; iT < nTrk; ++iT) {
      TrackExtended* trkExt = _trkImplVec[iT];
      float tanLambda = trkExt->getTanLambda();
      float product = pos[2] * tanLambda;
      // check that the hit and track are in the same z-half, which won't work for the rare cases of something going
      // backwards ...

      if (product > 0) {
        // use the previously created trackextrapolations for the
        HelixClass* helix = _trackExtrapolatedHelix[trkExt];

        // skip if the extrapolations failed
        if (helix == 0) {
          streamlog_out(DEBUG3) << "helix extrapolation failed for trkExt" << std::endl;
          continue;
        }

        float distance = helix->getDistanceToPoint(pos, dcut);

        streamlog_out(DEBUG1) << "for helix extrapolation " << helix << " distance = " << distance << std::endl;

        // check the distance is less than the steerable cut ...
        if (distance < dcut) {
          streamlog_out(DEBUG3) << "for helix extrapolation " << helix << " distance = " << distance << std::endl;

          // ... if so create the association and flag the hit and track
          TrackHitPair* trkHitPair = new TrackHitPair(trkExt, trkHitExt, distance);
          pairs.push_back(trkHitPair);
          flagTrack[trkExt] = true;
          flagHit[trkHitExt] = true;
        }
      }
    }
  }

  int nPairs = int(pairs.size());

  streamlog_out(DEBUG3) << "AssignOuterHitsToTracks : Number of track hit pairs to try =  " << nPairs << std::endl;

  if (nPairs > 0) {
    // sort the pairs on distance
    SortingTrackHitPairs(pairs);

    for (int iP = 0; iP < nPairs; ++iP) {
      TrackHitPair* trkHitPair = pairs[iP];
      TrackExtended* trkExt = trkHitPair->getTrackExtended();
      TrackerHitExtended* trkHitExt =

          trkHitPair->getTrackerHitExtended();

      // check if the track or hit is still free to be combined
      if (flagTrack[trkExt] && flagHit[trkHitExt]) {
        if (refit == 0) { // just set the association
          trkExt->addTrackerHitExtended(trkHitExt);
          trkHitExt->setUsedInFit(false);
          trkHitExt->setTrackExtended(trkExt);
        }

        else {
          // get all the hits already included in the track
          TrackerHitExtendedVec hitsInTrack = trkExt->getTrackerHitExtendedVec();

          int nTotH = int(hitsInTrack.size());
          int nHitsInFit = 0;

          for (int iTH = 0; iTH < nTotH; ++iTH) {
            // count the number of hits used in the fit
            TrackerHitExtended* hitInTrack = hitsInTrack[iTH];
            if (hitInTrack->getUsedInFit()) {
              nHitsInFit++;
            }
          }

          int iHitInFit = 0;

          // add the previously used hits from the track to the vectors

          EVENT::TrackerHitVec trkHits;

          for (int iHit = 0; iHit < nTotH; ++iHit) {
            TrackerHitExtended* hitInTrack = hitsInTrack[iHit];
            if (hitInTrack->getUsedInFit()) {
              TrackerHit* hit = hitInTrack->getTrackerHit();
              iHitInFit++;
              if (hit) {
                trkHits.push_back(hit);
              } else {
                throw EVENT::Exception(
                    std::string("FullLDCTracking_MarlinTrk::AssignOuterHitsToTracks: TrackerHit pointer == NULL "));
              }
            }
          }

          // add the hit to be attached to the vectors
          TrackerHit* remainHit = trkHitExt->getTrackerHit();
          iHitInFit++;
          trkHits.push_back(remainHit);

          double chi2_D;
          int ndf;

          if (trkHits.size() < 3)
            return;

          // sort the hits in R
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

          for (std::vector<std::pair<float, EVENT::TrackerHit*>>::iterator it = r2_values.begin();
               it != r2_values.end(); ++it) {
            trkHits.push_back(it->second);
          }

          streamlog_out(DEBUG3) << "AssignOuterHitsToTracks: Start Fitting: AddHits: number of hits to fit "
                                << trkHits.size() << std::endl;

          MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();

          IMPL::TrackStateImpl pre_fit;

          pre_fit.setD0(trkExt->getD0());
          pre_fit.setPhi(trkExt->getPhi());
          pre_fit.setZ0(trkExt->getZ0());
          pre_fit.setOmega(trkExt->getOmega());
          pre_fit.setTanLambda(trkExt->getTanLambda());

          float ref[3];
          ref[0] = ref[1] = ref[2] = 0.0;

          pre_fit.setReferencePoint(ref);

          pre_fit.setLocation(lcio::TrackStateImpl::AtIP);
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

          pre_fit.setCovMatrix(covMatrix);

          int error =
              MarlinTrk::createFit(trkHits, marlin_trk, &pre_fit, _bField, IMarlinTrack::backward, _maxChi2PerHit);

          if (error != IMarlinTrack::success) {
            streamlog_out(DEBUG3)
                << "FullLDCTracking_MarlinTrk::AssignOuterHitsToTracks: creation of fit fails with error " << error
                << std::endl;

            delete marlin_trk;
            continue;
          }

          std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;
          marlin_trk->getOutliers(outliers);

          float outlier_pct = outliers.size() / float(trkHits.size());

          streamlog_out(DEBUG1) << "FullLDCTracking_MarlinTrk::AssignOuterHitsToTracks: percentage of outliers "
                                << outlier_pct << std::endl;

          if (outlier_pct > _maxAllowedPercentageOfOutliersForTrackCombination) {
            streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::AssignOuterHitsToTracks: percentage of outliers "
                                  << outlier_pct << " is greater than cut maximum: "
                                  << _maxAllowedPercentageOfOutliersForTrackCombination << std::endl;
            delete marlin_trk;
            continue;
          }

          const Vector3D point(0., 0., 0.); // nominal IP
          TrackStateImpl trkState;
          // int return_code =
          marlin_trk->propagate(point, trkState, chi2_D, ndf);

          delete marlin_trk;

          if (error != IMarlinTrack::success) {
            streamlog_out(DEBUG3)
                << "FullLDCTracking_MarlinTrk::AssignOuterHitsToTracks: propagate to IP fails with error " << error
                << std::endl;

            continue;
          }

          if (ndf < 0) {
            streamlog_out(DEBUG2)
                << "FullLDCTracking_MarlinTrk::AssignOuterHitsToTracks: Fit failed : NDF is less that zero  " << ndf
                << std::endl;

            continue;
          }

          float chi2Fit = chi2_D / float(ndf);

          if (chi2Fit > _chi2FitCut) {
            streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::AssignOuterHitsToTracks: track fail Chi2 cut of "
                                  << _chi2FitCut << " chi2 of track = " << chi2Fit << std::endl;

            return;
          }

          // note trackAR which is of type TrackExtended, only takes fits with ref point = 0,0,0
          trkExt->setOmega(trkState.getOmega());
          trkExt->setTanLambda(trkState.getTanLambda());
          trkExt->setPhi(trkState.getPhi());
          trkExt->setD0(trkState.getD0());
          trkExt->setZ0(trkState.getZ0());

          float cov[15];

          for (int i = 0; i < 15; ++i) {
            cov[i] = trkState.getCovMatrix().operator[](i);
          }

          trkExt->setCovMatrix(cov);
          trkExt->setChi2(chi2_D);
          trkExt->setNDF(ndf);

          trkExt->addTrackerHitExtended(trkHitExt);
          trkHitExt->setTrackExtended(trkExt);
          trkHitExt->setUsedInFit(true);
          flagTrack[trkExt] = false;
          flagHit[trkHitExt] = false;

          streamlog_out(DEBUG2) << "AssignOuterHitsToTracks: Hit " << trkHitExt << " successfully assigned to track "
                                << trkExt << std::endl;
        }
      }
    }

    for (int iP = 0; iP < nPairs; ++iP) {
      TrackHitPair* trkHitPair = pairs[iP];
      delete trkHitPair;
    }

    pairs.clear();
  }
}

HelixClass* FullLDCTracking_MarlinTrk::GetExtrapolationHelix(TrackExtended* track) {
  streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::GetExtrapolationHelix called for track " << track << std::endl;

  HelixClass* helix = 0;

  // get the track state at the last hit of the track with ref point of greatest abs z

  GroupTracks* group = track->getGroupTracks();

  TrackExtendedVec trk_vec = group->getTrackExtendedVec();

  const EVENT::TrackState* ts_at_calo = 0;
  float z_ref_max = 0;

  streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::GetExtrapolationHelix number of tracks associated = "
                        << trk_vec.size() << std::endl;

  for (unsigned itrk = 0; itrk < trk_vec.size(); ++itrk) {
    EVENT::Track* trk_lcio = trk_vec[itrk]->getTrack();

    if (trk_lcio) {
      // use the tracks state at the calorimeter because that will have accounted for energy loss already
      if (trk_lcio->getTrackState(lcio::TrackState::AtCalorimeter)) {
        const EVENT::TrackState* ts_at_last_hit = trk_lcio->getTrackState(lcio::TrackState::AtLastHit);
        float z_ref = ts_at_last_hit->getReferencePoint()[2];

        // make sure we use the one closest to the calo face
        if (fabs(z_ref) > z_ref_max) {
          z_ref_max = fabs(z_ref);
          ts_at_calo = trk_lcio->getTrackState(lcio::TrackState::AtCalorimeter);

          streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::GetExtrapolationHelix set ts_at_calo with ref_z = "
                                << z_ref << std::endl;
        }
      }
    }
  }

  if (ts_at_calo) {
    IMPL::TrackStateImpl ts_at_calo_forIP(*ts_at_calo);

    LCIOTrackPropagators::PropagateLCIOToNewRef(ts_at_calo_forIP, 0.0, 0.0, 0.0);

    ts_at_calo_forIP.setLocation(lcio::TrackStateImpl::AtIP);

    helix = new HelixClass();

    helix->Initialize_Canonical(ts_at_calo_forIP.getPhi(), ts_at_calo_forIP.getD0(), ts_at_calo_forIP.getZ0(),
                                ts_at_calo_forIP.getOmega(), ts_at_calo_forIP.getTanLambda(), _bField);

    streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::GetExtrapolationHelix helix created at IP" << std::endl;
  }

  return helix;
}

void FullLDCTracking_MarlinTrk::AssignTPCHitsToTracks(TrackerHitExtendedVec hitVec, float dcut) {
  int nHits = int(hitVec.size());
  int nTrk = int(_trkImplVec.size());

  for (int iT = 0; iT < nTrk; ++iT) { // loop over all tracks
    TrackExtended* foundTrack = _trkImplVec[iT];
    GroupTracks* group = foundTrack->getGroupTracks();
    TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
    int nTrkGrp = int(tracksInGroup.size());
    for (int iTrkGrp = 0; iTrkGrp < nTrkGrp; ++iTrkGrp) {
      TrackExtended* trkGrp = tracksInGroup[iTrkGrp];
      TrackerHitExtendedVec hitVecGrp = trkGrp->getTrackerHitExtendedVec();
      int nHits_Grp = int(hitVecGrp.size());
      float zMin = 1.0e+20;
      float zMax = -1.0e+20;
      float startPoint[3] = {0., 0., 0.};
      float endPoint[3] = {0., 0., 0.};
      for (int iH = 0; iH < nHits_Grp; ++iH) {
        TrackerHitExtended* trkHitExt = hitVecGrp[iH];
        float pos[3] = {0., 0., 0.};
        for (int iC = 0; iC < 3; ++iC)
          pos[iC] = float(trkHitExt->getTrackerHit()->getPosition()[iC]);
        if (pos[2] > zMax) {
          zMax = pos[2];
          for (int iC = 0; iC < 3; ++iC)
            endPoint[iC] = pos[iC];
        }
        if (pos[2] < zMin) {
          zMin = pos[2];
          for (int iC = 0; iC < 3; ++iC)
            startPoint[iC] = pos[iC];
        }
      }
      trkGrp->setStart(startPoint);
      trkGrp->setEnd(endPoint);
    }
  }

  // replace previous version with faster loop ordering

  std::vector<float> minDistances(nHits, dcut);
  std::vector<TrackExtended*> tracksToAttach(nHits);
  std::vector<std::vector<float>> HitPositions(nHits);
  std::vector<int> HitSign(nHits);     // Positive or Negative side
  for (int iH = 0; iH < nHits; ++iH) { // loop over leftover TPC hits
    tracksToAttach[iH] = NULL;
    // Get all TrackerHit positions, so we only have to get them once
    TrackerHit* temphit = hitVec[iH]->getTrackerHit();
    const double* temppos = temphit->getPosition();
    HitPositions[iH].push_back(float(temppos[0]));
    HitPositions[iH].push_back(float(temppos[1]));
    HitPositions[iH].push_back(float(temppos[2]));
    HitSign[iH] = std::signbit(temppos[2]);
  }

  streamlog_out(DEBUG3) << "AssignTPCHitsToTracks: Starting loop " << nTrk << " tracks   and  " << nHits << " hits"
                        << std::endl;

  for (int iT = 0; iT < nTrk; ++iT) { // loop over all tracks
    TrackExtended* foundTrack = _trkImplVec[iT];
    int tanlambdaSign = std::signbit(foundTrack->getTanLambda()); // we only care about positive or negative
    GroupTracks* group = foundTrack->getGroupTracks();
    TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
    int nTrkGrp = int(tracksInGroup.size());

    for (int iTrkGrp = 0; iTrkGrp < nTrkGrp; ++iTrkGrp) {
      TrackExtended* trkGrp = tracksInGroup[iTrkGrp];
      float tanLambda = trkGrp->getTanLambda();
      float omega = trkGrp->getOmega();
      float d0 = trkGrp->getD0();
      float z0 = trkGrp->getZ0();
      float phi0 = trkGrp->getPhi();
      float startPointZ = trkGrp->getStart()[2];
      float endPointZ = trkGrp->getEnd()[2];
      HelixClass helix;
      helix.Initialize_Canonical(phi0, d0, z0, omega, tanLambda, _bField);
      float OnePFivehalfPeriodZ = 1.5 * fabs(acos(-1.) * tanLambda / omega);

      for (int iH = 0; iH < nHits; ++iH) { // loop over leftover TPC hits

        // check if the hit and the track or on the same side
        // xor return 1, if hits are different
        if (tanlambdaSign ^ HitSign[iH])
          continue;

        float DeltaStart = fabs(HitPositions[iH][2] - startPointZ);
        float DeltaEnd = fabs(HitPositions[iH][2] - endPointZ);
        bool consider = DeltaStart <= OnePFivehalfPeriodZ;
        consider = consider || (DeltaEnd <= OnePFivehalfPeriodZ);
        consider = consider || ((HitPositions[iH][2] >= startPointZ) && (HitPositions[iH][2] <= endPointZ));

        if (consider) {
          float distance = helix.getDistanceToPoint(HitPositions[iH], minDistances[iH]);
          if (distance < minDistances[iH]) {
            minDistances[iH] = distance;
            tracksToAttach[iH] = foundTrack;
          }
        }
      } // loop over leftover TPC hits
    } // groups in tracks
  } // loop over all tracks

  for (int iH = 0; iH < nHits; ++iH) {
    TrackerHitExtended* trkHitExt = hitVec[iH];
    if (tracksToAttach[iH] != NULL) {
      tracksToAttach[iH]->addTrackerHitExtended(trkHitExt);
      trkHitExt->setTrackExtended(tracksToAttach[iH]);
      trkHitExt->setUsedInFit(false);
    }
  }

  streamlog_out(DEBUG3) << " Fast loop done " << std::endl;

  //     for (int iH=0;iH<nHits;iH++) { // loop over leftover TPC hits
  //    TrackerHitExtended * hitExt = hitVec[iH];
  //    float pos[3];
  //    for (int ip=0;ip<3;++ip)
  //        pos[ip] = float(hitExt->getTrackerHit()->getPosition()[ip]);
  //    float minDist = 1.0e+20;
  //    TrackExtended * trackToAttach = NULL;
  //    for (int iT=0;iT<nTrk;++iT) { // loop over all tracks
  //        TrackExtended * foundTrack = _trkImplVec[iT];
  //        float tanLambdaFound = foundTrack->getTanLambda();
  //        float product = tanLambdaFound*pos[2];
  //        if (product>0) {
  //          GroupTracks * group = foundTrack->getGroupTracks();
  //          TrackExtendedVec tracksInGroup = group->getTrackExtendedVec();
  //          int nTrkGrp = int(tracksInGroup.size());
  //          for (int iTrkGrp=0;iTrkGrp<nTrkGrp;++iTrkGrp) {
  //            TrackExtended * trkGrp = tracksInGroup[iTrkGrp];
  //            float tanLambda = trkGrp->getTanLambda();
  //            float omega = trkGrp->getOmega();
  //            float d0 = trkGrp->getD0();
  //            float z0 = trkGrp->getZ0();
  //            float phi0 = trkGrp->getPhi();
  //            float dist[3];
  //            float startPoint[3];
  //            float endPoint[3];
  //            for (int iC=0;iC<3;++iC) {
  //              startPoint[iC] = trkGrp->getStart()[iC];
  //              endPoint[iC] = trkGrp->getEnd()[iC];
  //            }
  //            HelixClass helix;
  //            helix.Initialize_Canonical(phi0,d0,z0,omega,tanLambda,_bField);
  //            float halfPeriodZ = fabs(acos(-1.)*tanLambda/omega);
  //            helix.getDistanceToPoint(pos,dist);
  //            float DeltaStart = fabs(pos[2]-startPoint[2]);
  //            float DeltaEnd = fabs(pos[2]-endPoint[2]);
  //            bool consider = DeltaStart <= 1.5*halfPeriodZ;
  //            consider = consider || (DeltaEnd <= 1.5*halfPeriodZ);
  //            consider = consider || ( (pos[2]>=startPoint[2]) && (pos[2]<=endPoint[2]) );
  // //                 float ZMin = DeltaStart;
  // //                 if (DeltaEnd<ZMin)
  // //                   ZMin = DeltaEnd;
  //            if (dist[2]<dcut && consider && dist[2]<minDist) {
  //              minDist = dist[2];
  //              trackToAttach = foundTrack;
  //            }
  //          }
  //        }
  //    }

  //    if (trackToAttach!=NULL) {
  //   trackToAttach->addTrackerHitExtended(hitExt);
  //  hitExt->setTrackExtended( trackToAttach );
  //  hitExt->setUsedInFit( false );
  //  if(trackToAttach!=tracksToAttach[iH])std::cout << " Check Failed" << trackToAttach << "  " << tracksToAttach[iH]
  //  << std::endl;
  //
  //}
  // else {
  ///     std::cout << iH << " hit is not assigned : distance to closest track = " << minDist << std::endl;
  ///}
  //}
  // std::cout << " Slow loop done " << std::endl;
}

void FullLDCTracking_MarlinTrk::AssignSiHitsToTracks(TrackerHitExtendedVec hitVec, float dcut) {
  int nHits = int(hitVec.size());
  int nTrk = int(_allNonCombinedTPCTracks.size());

  streamlog_out(DEBUG3) << "AssignSiHitsToTracks : Number of hits to assign " << hitVec.size()
                        << " : Number of available tracks = " << nTrk << std::endl;

  std::map<TrackExtended*, bool> flagTrack;
  std::map<TrackerHitExtended*, bool> flagHit;
  TrackHitPairVec pairs;
  flagTrack.clear();
  flagHit.clear();
  pairs.clear();

  for (int iH = 0; iH < nHits; ++iH) {
    float pos[3];
    TrackerHitExtended* trkHitExt = hitVec[iH];
    TrackerHit* hit = trkHitExt->getTrackerHit();

    for (int ip = 0; ip < 3; ++ip) {
      pos[ip] = float(hit->getPosition()[ip]);
    }

    for (int iT = 0; iT < nTrk; ++iT) {
      TrackExtended* trkExt = _allNonCombinedTPCTracks[iT];

      float tanLambda = trkExt->getTanLambda();
      float product = pos[2] * tanLambda;

      streamlog_out(DEBUG0) << "AssignSiHitsToTracks : product =  " << product << " z hit = " << pos[2] << std::endl;

      if (product > 0) {
        float d0 = trkExt->getD0();
        float z0 = trkExt->getZ0();
        float phi0 = trkExt->getPhi();
        float omega = trkExt->getOmega();
        tanLambda = trkExt->getTanLambda();

        HelixClass helix;
        helix.Initialize_Canonical(phi0, d0, z0, omega, tanLambda, _bField);
        float distance = helix.getDistanceToPoint(pos, dcut);

        streamlog_out(DEBUG0) << "AssignSiHitsToTracks : distance =  " << distance << " cut = " << dcut << std::endl;

        if (distance < dcut) {
          TrackHitPair* trkHitPair = new TrackHitPair(trkExt, trkHitExt, distance);
          pairs.push_back(trkHitPair);
          flagTrack[trkExt] = true;
          flagHit[trkHitExt] = true;
        }
      }
    }
  }

  int nPairs = int(pairs.size());

  streamlog_out(DEBUG3) << "AssignSiHitsToTracks : Number of track hit pairs to try =  " << nPairs << std::endl;

  if (nPairs > 0) {
    SortingTrackHitPairs(pairs);

    for (int iP = 0; iP < nPairs; ++iP) {
      TrackHitPair* trkHitPair = pairs[iP];
      TrackExtended* trkExt = trkHitPair->getTrackExtended();
      TrackerHitExtended* trkHitExt =

          trkHitPair->getTrackerHitExtended();

      if (flagTrack[trkExt] && flagHit[trkHitExt]) {
        // get the hits already assigned to the track
        TrackerHitExtendedVec hitsInTrack = trkExt->getTrackerHitExtendedVec();

        int nTotH = int(hitsInTrack.size());
        int nHitsInFit = 0;

        for (int iTH = 0; iTH < nTotH; ++iTH) {
          TrackerHitExtended* hitInTrack = hitsInTrack[iTH];

          // count the number of hits used in the fit
          if (hitInTrack->getUsedInFit()) {
            nHitsInFit++;
          }
        }

        int iHitInFit = 0;

        // add the previously used hits from the track to the vectors

        EVENT::TrackerHitVec trkHits;

        for (int iHit = 0; iHit < nTotH; ++iHit) {
          TrackerHitExtended* hitInTrack = hitsInTrack[iHit];
          if (hitInTrack->getUsedInFit()) {
            TrackerHit* hit = hitInTrack->getTrackerHit();
            iHitInFit++;
            if (hit) {
              trkHits.push_back(hit);
            } else {
              throw EVENT::Exception(
                  std::string("FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: TrackerHit pointer == NULL "));
            }
          }
        }

        // add the hit to be attached to the vectors
        TrackerHit* remainHit = trkHitExt->getTrackerHit();
        iHitInFit++;
        trkHits.push_back(remainHit);

        double chi2_D;
        int ndf;

        if (trkHits.size() < 3)
          return;

        // sort the hits in R
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

        for (std::vector<std::pair<float, EVENT::TrackerHit*>>::iterator it = r2_values.begin(); it != r2_values.end();
             ++it) {
          trkHits.push_back(it->second);
        }

        streamlog_out(DEBUG2) << "AssignSiHitsToTracks: Start Fitting: AddHits: number of hits to fit "
                              << trkHits.size() << std::endl;

        MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();

        IMPL::TrackStateImpl pre_fit;

        pre_fit.setD0(trkExt->getD0());
        pre_fit.setPhi(trkExt->getPhi());
        pre_fit.setZ0(trkExt->getZ0());
        pre_fit.setOmega(trkExt->getOmega());
        pre_fit.setTanLambda(trkExt->getTanLambda());

        float ref[3];
        ref[0] = ref[1] = ref[2] = 0.0;

        pre_fit.setReferencePoint(ref);

        pre_fit.setLocation(lcio::TrackStateImpl::AtIP);

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

        pre_fit.setCovMatrix(covMatrix);

        int error =
            MarlinTrk::createFit(trkHits, marlin_trk, &pre_fit, _bField, IMarlinTrack::backward, _maxChi2PerHit);

        if (error != IMarlinTrack::success) {
          streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: creation of fit fails with error "
                                << error << std::endl;

          delete marlin_trk;
          continue;
        }

        std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;
        marlin_trk->getOutliers(outliers);

        float outlier_pct = outliers.size() / float(trkHits.size());

        streamlog_out(DEBUG1) << "FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: percentage of outliers "
                              << outlier_pct << std::endl;

        if (outlier_pct > _maxAllowedPercentageOfOutliersForTrackCombination) {
          streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: percentage of outliers "
                                << outlier_pct << " is greater than cut maximum: "
                                << _maxAllowedPercentageOfOutliersForTrackCombination << std::endl;
          delete marlin_trk;
          continue;
        }

        const Vector3D point(0., 0., 0.); // nominal IP

        TrackStateImpl trkState;
        // int return_code =
        marlin_trk->propagate(point, trkState, chi2_D, ndf);

        delete marlin_trk;

        if (error != IMarlinTrack::success) {
          streamlog_out(DEBUG3) << "FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: propagate to IP fails with error "
                                << error << std::endl;

          continue;
        }

        if (ndf < 0) {
          streamlog_out(DEBUG2)
              << "FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: Fit failed : NDF is less that zero  " << ndf
              << std::endl;

          continue;
        }

        float chi2Fit = chi2_D / float(ndf);

        if (chi2Fit > _chi2FitCut) {
          streamlog_out(DEBUG2) << "FullLDCTracking_MarlinTrk::AssignSiHitsToTracks: track fail Chi2 cut of "
                                << _chi2FitCut << " chi2 of track = " << chi2Fit << std::endl;

          continue;
        }

        // note trackAR which is of type TrackExtended, only takes fits with ref point = 0,0,0
        trkExt->setOmega(trkState.getOmega());
        trkExt->setTanLambda(trkState.getTanLambda());
        trkExt->setPhi(trkState.getPhi());
        trkExt->setD0(trkState.getD0());
        trkExt->setZ0(trkState.getZ0());

        float cov[15];

        for (int i = 0; i < 15; ++i) {
          cov[i] = trkState.getCovMatrix().operator[](i);
        }

        trkExt->setCovMatrix(cov);
        trkExt->setChi2(chi2_D);
        trkExt->setNDF(ndf);

        trkExt->addTrackerHitExtended(trkHitExt);
        trkHitExt->setTrackExtended(trkExt);
        trkHitExt->setUsedInFit(true);
        flagTrack[trkExt] = false;
        flagHit[trkHitExt] = false;

        streamlog_out(DEBUG2) << "AssignSiHitsToTracks: Hit " << trkHitExt << " successfully assigned to track "
                              << trkExt << std::endl;
      }
    }

    for (int iP = 0; iP < nPairs; ++iP) {
      TrackHitPair* trkHitPair = pairs[iP];
      delete trkHitPair;
    }

    pairs.clear();
  }
}

void FullLDCTracking_MarlinTrk::PrintOutMerging(TrackExtended* firstTrackExt, TrackExtended* secondTrackExt, int iopt) {
  // iopt = 1 merged Si and TPC merging
  // iopt = 2 merged Si and TPC using forced merging
  // iopt = 3 merged TPC segments merging
  // iopt = 4 merged Comb Si and TPC merging
  // iopt = 5 merged Comb TPC and TPC merging
  // iopt = 6 unmerged TPC and Si segments ( when using soft merging)
  // iopt = 7 unmerged TPC and Si segments ( when using forced merging)
  // iopt = 8 unmerged Comb and TPC
  // iopt = 9 unmerged TPC segments

  char strg[200];

  try {
    Track* firstTrack = firstTrackExt->getTrack();
    Track* secondTrack = secondTrackExt->getTrack();

    std::string firstColName = _TPCTrackMCPCollName;
    std::string secondColName = _TPCTrackMCPCollName;

    if (iopt == 1) {
      secondColName = _SiTrackMCPCollName;
    } else if (iopt == 2) {
      secondColName = _SiTrackMCPCollName;
    } else if (iopt == 3) {
      secondColName = _TPCTrackMCPCollName;
    } else if (iopt == 4) {
      secondColName = _SiTrackMCPCollName;
    } else if (iopt == 5) {
      secondColName = _TPCTrackMCPCollName;
    } else if (iopt == 6) {
      secondColName = _SiTrackMCPCollName;
    } else if (iopt == 7) {
      secondColName = _SiTrackMCPCollName;
    } else if (iopt == 8) {
      secondColName = _SiTrackMCPCollName;
    } else {
      secondColName = _TPCTrackMCPCollName;
    }

    // get the relation collections, this should be done only once for each event ...
    LCCollection* firstCol = _evt->getCollection(firstColName.c_str());
    LCCollection* secondCol = _evt->getCollection(secondColName.c_str());

    // get navigators
    LCRelationNavigator firstNav(firstCol);
    LCRelationNavigator secondNav(secondCol);

    // get the MCParticles with the greatest weight for the first track
    LCObjectVec firstVec = firstNav.getRelatedToObjects(firstTrack);
    FloatVec firstWeights = firstNav.getRelatedToWeights(firstTrack);

    LCObject* firstMCP = NULL;

    float firstWght = 0;
    int nObj = firstVec.size();
    for (int iObj = 0; iObj < nObj; ++iObj) {
      if (firstWeights[iObj] > firstWght) {
        firstWght = firstWeights[iObj];
        firstMCP = firstVec[iObj];
      }
    }

    // get the MCParticles with the greatest weight for the second track
    LCObjectVec secondVec = secondNav.getRelatedToObjects(secondTrack);
    FloatVec secondWeights = secondNav.getRelatedToWeights(secondTrack);

    LCObject* secondMCP = NULL;

    float secondWght = 0;
    nObj = secondVec.size();
    for (int iObj = 0; iObj < nObj; ++iObj) {
      if (secondWeights[iObj] > secondWght) {
        secondWght = secondWeights[iObj];
        secondMCP = secondVec[iObj];
      }
    }

    // get the track parameters for both tracks and get the 3-momentum using the HelixClass
    float d0First = firstTrackExt->getD0();
    float z0First = firstTrackExt->getZ0();
    float omegaFirst = firstTrackExt->getOmega();
    float tanLFirst = firstTrackExt->getTanLambda();
    float phi0First = firstTrackExt->getPhi();

    float d0Second = secondTrackExt->getD0();
    float z0Second = secondTrackExt->getZ0();
    float omegaSecond = secondTrackExt->getOmega();
    float tanLSecond = secondTrackExt->getTanLambda();
    float phi0Second = secondTrackExt->getPhi();

    HelixClass firstHelix;
    firstHelix.Initialize_Canonical(phi0First, d0First, z0First, omegaFirst, tanLFirst, _bField);
    float pxFirst = firstHelix.getMomentum()[0];
    float pyFirst = firstHelix.getMomentum()[1];
    float pzFirst = firstHelix.getMomentum()[2];

    HelixClass secondHelix;
    secondHelix.Initialize_Canonical(phi0Second, d0Second, z0Second, omegaSecond, tanLSecond, _bField);
    float pxSecond = secondHelix.getMomentum()[0];
    float pySecond = secondHelix.getMomentum()[1];
    float pzSecond = secondHelix.getMomentum()[2];

    // get the momentum differences
    float dPx = pxFirst + pxSecond;
    float dPy = pyFirst + pySecond;
    float dPz = pzFirst + pzSecond;

    float dPplus = sqrt(dPx * dPx + dPy * dPy + dPz * dPz);

    dPx = pxFirst - pxSecond;
    dPy = pyFirst - pySecond;
    dPz = pzFirst - pzSecond;

    float dPminus = sqrt(dPx * dPx + dPy * dPy + dPz * dPz);

    // get momentum for each track
    float pFirst = sqrt(pxFirst * pxFirst + pyFirst * pyFirst + pzFirst * pzFirst);
    float pSecond = sqrt(pxSecond * pxSecond + pySecond * pySecond + pzSecond * pzSecond);

    float ptFirst = sqrt(pxFirst * pxFirst + pyFirst * pyFirst);
    float ptSecond = sqrt(pxSecond * pxSecond + pySecond * pySecond);

    const float sigmaPOverPFirst = sqrt(firstTrackExt->getCovMatrix()[5]) / fabs(omegaFirst);
    const float sigmaPOverPSecond = sqrt(secondTrackExt->getCovMatrix()[5]) / fabs(omegaSecond);

    const float sigmaPFirst = pFirst * sigmaPOverPFirst;
    const float sigmaPSecond = pSecond * sigmaPOverPSecond;

    float den = pFirst;

    if (pSecond < pFirst)
      den = pSecond;

    dPplus = dPplus / den;
    dPminus = dPminus / den;

    // now check if this was a Erroneous merger ...
    if (firstMCP != secondMCP && iopt < 6) {
      if (iopt == 1) {
        streamlog_out(DEBUG4) << "Erroneous merging of Si and TPC segments (iopt=1) mcp first = " << firstMCP
                              << " mcp second = " << secondMCP << " --->" << std::endl;
      } else if (iopt == 2) {
        streamlog_out(DEBUG4) << "Erroneous merging of Si and TPC segments (iopt=2) mcp first = " << firstMCP
                              << " mcp second = " << secondMCP << " --->" << std::endl;
      } else if (iopt == 3) {
        streamlog_out(DEBUG4) << "Erroneous merging of TPC segments (iopt=3) mcp first = " << firstMCP
                              << " mcp second = " << secondMCP << " ---> " << std::endl;
      } else if (iopt == 4) {
        streamlog_out(DEBUG4) << "Erroneous merging of combSi segment with uncombTPC segment (iopt=4) mcp first = "
                              << firstMCP << " mcp second = " << secondMCP << " ---> " << std::endl;
      } else {
        streamlog_out(DEBUG4) << "Erroneous merging of combTPC segment with uncombTPC segment (iopt=5) mcp first = "
                              << firstMCP << " mcp second = " << secondMCP << " --->" << std::endl;
      }

      streamlog_out(DEBUG4) << "    p         error      pt       D0      Z0     Px      Py      Pz      wieght"
                            << std::endl;

      sprintf(strg, "%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ", pFirst, sigmaPFirst, ptFirst, d0First,
              z0First, pxFirst, pyFirst, pzFirst);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg, "  %5.3f\n", firstWght);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg, "%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ", pSecond, sigmaPSecond, ptSecond,
              d0Second, z0Second, pxSecond, pySecond, pzSecond);
      streamlog_out(DEBUG4) << strg;
      sprintf(strg, "  %5.3f\n", secondWght);
      streamlog_out(DEBUG4) << strg;

      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;

      TrackExtended* combinedTrack =
          CombineTracks(firstTrackExt, secondTrackExt, _maxAllowedPercentageOfOutliersForTrackCombination, true);

      if (combinedTrack != NULL) {
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. "
                              << firstTrackExt->getNDF() + secondTrackExt->getNDF() + 5 << std::endl;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      } else {
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt)
                            << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;

      streamlog_out(DEBUG4) << std::endl;
    }

    // ... or if it was an incorrect TPC to TPC rejection ...
    else if (firstMCP == secondMCP && ((iopt == 8) || (iopt == 9))) {
      float deltaOmega = _dOmegaForMerging;
      float deltaAngle = _angleForMerging;
      if (iopt == 8) {
        streamlog_out(DEBUG4) << "Unmerged TPC and Comb segments (iopt=8) --->" << std::endl;
      } else {
        streamlog_out(DEBUG4) << "Unmerged TPC segments (iopt=9) --->" << std::endl;
      }

      float qFirst = PIOVER2 - atan(tanLFirst);
      float qSecond = PIOVER2 - atan(tanLSecond);
      float dOmega = fabs((omegaFirst - omegaSecond) / omegaSecond);
      float angle = (cos(phi0First) * cos(phi0Second) + sin(phi0First) * sin(phi0Second)) * sin(qFirst) * sin(qSecond) +
                    cos(qFirst) * cos(qSecond);
      angle = acos(angle);

      streamlog_out(DEBUG4) << " dOmegaCut = " << deltaOmega << " AngleCut = " << deltaAngle << " dOmega = " << dOmega
                            << " angle = " << angle << std::endl;

      streamlog_out(DEBUG4) << "    p         error      pt       D0      Z0     Px      Py      Pz      wieght"
                            << std::endl;

      sprintf(strg, "%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ", pFirst, sigmaPFirst, ptFirst, d0First,
              z0First, pxFirst, pyFirst, pzFirst);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg, "  %5.3f\n", firstWght);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg, "%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ", pSecond, sigmaPSecond, ptSecond,
              d0Second, z0Second, pxSecond, pySecond, pzSecond);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg, "  %5.3f\n", secondWght);
      streamlog_out(DEBUG4) << strg;

      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;

      TrackExtended* combinedTrack =
          CombineTracks(firstTrackExt, secondTrackExt, _maxAllowedPercentageOfOutliersForTrackCombination, true);
      if (combinedTrack != NULL) {
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. "
                              << firstTrackExt->getNDF() + secondTrackExt->getNDF() + 5 << std::endl;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      } else {
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt)
                            << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;

      streamlog_out(DEBUG4) << std::endl;
    }

    // ... or if it was an incorrect TPC to Si rejection ...
    else if (firstMCP == secondMCP && ((iopt == 6) || (iopt == 7))) {
      float deltaOmega = _dOmegaForMerging;
      float deltaAngle = _angleForMerging;

      if (iopt == 6) {
        streamlog_out(DEBUG4) << "Unmerged TPC and Si segments (iopt=6) mcp first = " << firstMCP
                              << " mcp second = " << secondMCP << " --->" << std::endl;
      } else {
        streamlog_out(DEBUG4) << "Unmerged TPC and Si segments (iopt=7) mcp first = " << firstMCP
                              << " mcp second = " << secondMCP << " --->" << std::endl;
        deltaOmega = _dOmegaForForcedMerging;
        deltaAngle = _angleForForcedMerging;
      }

      float qFirst = PIOVER2 - atan(tanLFirst);
      float qSecond = PIOVER2 - atan(tanLSecond);

      float dOmega = fabs((omegaFirst - omegaSecond) / omegaSecond);
      float angle = (cos(phi0First) * cos(phi0Second) + sin(phi0First) * sin(phi0Second)) * sin(qFirst) * sin(qSecond) +
                    cos(qFirst) * cos(qSecond);
      angle = acos(angle);

      streamlog_out(DEBUG4) << " dOmegaCut = " << deltaOmega << " AngleCut = " << deltaAngle << " dOmega = " << dOmega
                            << " angle = " << angle << std::endl;

      streamlog_out(DEBUG4) << "    p         error      pt       D0      Z0     Px      Py      Pz      wieght"
                            << std::endl;

      sprintf(strg, "%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ", pFirst, sigmaPFirst, ptFirst, d0First,
              z0First, pxFirst, pyFirst, pzFirst);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg, "  %5.3f\n", firstWght);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg, "%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ", pSecond, sigmaPSecond, ptSecond,
              d0Second, z0Second, pxSecond, pySecond, pzSecond);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg, "  %5.3f\n", secondWght);
      streamlog_out(DEBUG4) << strg;

      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;

      TrackExtended* combinedTrack =
          CombineTracks(firstTrackExt, secondTrackExt, _maxAllowedPercentageOfOutliersForTrackCombination, true);

      if (combinedTrack != NULL) {
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. "
                              << firstTrackExt->getNDF() + secondTrackExt->getNDF() + 5 << std::endl;
        delete combinedTrack->getGroupTracks();
        delete combinedTrack;
      } else {
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt)
                            << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;
      streamlog_out(DEBUG4) << std::endl;
    }
    // ... or if it was an correct merger ...
    else if (firstMCP == secondMCP && iopt < 6 && _debug > 3) {
      //      return;
      if (iopt == 1) {
        streamlog_out(DEBUG4) << "Correctly combining Si and TPC segments (iopt=1) mcp first = " << firstMCP
                              << " mcp second = " << secondMCP << " --->" << std::endl;
      } else if (iopt == 2) {
        streamlog_out(DEBUG4) << "Correctly merging of Si and TPC segments (iopt=2) mcp first = " << firstMCP
                              << " mcp second = " << secondMCP << " --->" << std::endl;
      } else if (iopt == 3) {
        streamlog_out(DEBUG4) << "Correctly merging of TPC segments (iopt=3) mcp first = " << firstMCP
                              << " mcp second = " << secondMCP << " ---> " << std::endl;
      } else if (iopt == 4) {
        streamlog_out(DEBUG4) << "Correctly merging of combSi segment with uncombTPC segment (iopt=4) mcp first = "
                              << firstMCP << " mcp second = " << secondMCP << " ---> " << std::endl;
      } else {
        streamlog_out(DEBUG4) << "Correctly merging of combTPC segment with uncombTPC segment (iopt=5) mcp first = "
                              << firstMCP << " mcp second = " << secondMCP << " --->" << std::endl;
      }

      streamlog_out(DEBUG4) << "    p         error      pt       D0      Z0     Px      Py      Pz      wieght"
                            << std::endl;

      sprintf(strg, "%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ", pFirst, sigmaPFirst, ptFirst, d0First,
              z0First, pxFirst, pyFirst, pzFirst);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg, "  %5.3f\n", firstWght);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg, "%7.2f +- %7.2f   %7.2f  %7.1f %7.1f %7.2f %7.2f %7.2f  ", pSecond, sigmaPSecond, ptSecond,
              d0Second, z0Second, pxSecond, pySecond, pzSecond);
      streamlog_out(DEBUG4) << strg;

      sprintf(strg, "  %5.3f\n", secondWght);
      streamlog_out(DEBUG4) << strg;

      streamlog_out(DEBUG4) << "Difference in +P = " << dPplus << "  -P = " << dPminus << std::endl;

      TrackExtended* combinedTrack =
          CombineTracks(firstTrackExt, secondTrackExt, _maxAllowedPercentageOfOutliersForTrackCombination, true);

      if (combinedTrack != NULL) {
        streamlog_out(DEBUG4) << "CombinedTrack " << combinedTrack->getNDF() << " c.f. "
                              << firstTrackExt->getNDF() + secondTrackExt->getNDF() + 5 << std::endl;
      } else {
        streamlog_out(DEBUG4) << "Could not combine track " << std::endl;
      }
      delete combinedTrack->getGroupTracks();
      delete combinedTrack;
      streamlog_out(DEBUG4) << " Overlap = " << SegmentRadialOverlap(firstTrackExt, secondTrackExt)
                            << " veto = " << VetoMerge(firstTrackExt, secondTrackExt) << std::endl;

      streamlog_out(DEBUG4) << std::endl;
    }
  }

  catch (DataNotAvailableException& e) {
  };
}

void FullLDCTracking_MarlinTrk::GeneralSorting(int* index, float* val, int direct, int nVal) {
  /**
   Sorting of index vector in ascending (0) /descending (!=0) order of val
   */

  float valOne, valTwo, valTemp;
  int indTemp;
  for (int i = 0; i < nVal; ++i) {
    index[i] = i;
  }

  for (int i = 0; i < nVal - 1; i++) {
    for (int j = 0; j < nVal - i - 1; j++) {
      valOne = val[j];
      valTwo = val[j + 1];
      bool order = valOne > valTwo;
      if (direct > 0)
        order = valOne <= valTwo;
      if (order) {
        valTemp = val[j];
        val[j] = val[j + 1];
        val[j + 1] = valTemp;
        indTemp = index[j];
        index[j] = index[j + 1];
        index[j + 1] = indTemp;
      }
    }
  }
}

int FullLDCTracking_MarlinTrk::SegmentRadialOverlap(TrackExtended* first, TrackExtended* second) {
  int nTrkGrpFirst = 0;
  int nTrkGrpSecond = 0;
  TrackerHitVec hitvecFirst;
  TrackerHitVec hitvecSecond;
  GroupTracks* groupFirst = first->getGroupTracks();
  GroupTracks* groupSecond = second->getGroupTracks();

  if (groupFirst != NULL) {
    TrackExtendedVec tracksInGroupFirst = groupFirst->getTrackExtendedVec();
    nTrkGrpFirst = int(tracksInGroupFirst.size());

    for (int iTrkGrp = 0; iTrkGrp < nTrkGrpFirst; ++iTrkGrp) {
      TrackExtended* trkGrp = tracksInGroupFirst[iTrkGrp];
      TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();

      for (unsigned int i = 0; i < hitVec.size(); ++i) {
        hitvecFirst.push_back(hitVec[i]->getTrackerHit());
      }
    }
  }
  if (groupSecond != NULL) {
    TrackExtendedVec tracksInGroupSecond = groupSecond->getTrackExtendedVec();
    nTrkGrpSecond = int(tracksInGroupSecond.size());

    for (int iTrkGrp = 0; iTrkGrp < nTrkGrpSecond; ++iTrkGrp) {
      TrackExtended* trkGrp = tracksInGroupSecond[iTrkGrp];
      TrackerHitExtendedVec hitVec = trkGrp->getTrackerHitExtendedVec();

      for (unsigned int i = 0; i < hitVec.size(); ++i) {
        hitvecSecond.push_back(hitVec[i]->getTrackerHit());
      }
    }
  }

  int nhitsFirst = (int)hitvecFirst.size();
  int nhitsSecond = (int)hitvecSecond.size();
  int count = 0;
  for (int i = 0; i < nhitsFirst; ++i) {
    float xi = (float)hitvecFirst[i]->getPosition()[0];
    float yi = (float)hitvecFirst[i]->getPosition()[1];
    float ri = sqrt(xi * xi + yi * yi);
    if (ri < _tpc_inner_r || ri > _tpc_pad_height)
      continue;
    for (int j = 0; j < nhitsSecond; ++j) {
      float xj = (float)hitvecSecond[j]->getPosition()[0];
      float yj = (float)hitvecSecond[j]->getPosition()[1];
      float rj = sqrt(xj * xj + yj * yj);
      if (fabs(ri - rj) < _tpc_pad_height / 2.0)
        count++;
    }
  }
  return count;
}

/*

 veto merger if the momentum of either track is less than 2.5 GeV
 or if following a full fit the NDF+10 of the combined tracks is less than the NDF_first + NDF_second

 NOTE: This function will try a full fit using CombineTracks if the momentum of both tracks is greater than
 VetoMergeMomentumCut

 */

bool FullLDCTracking_MarlinTrk::VetoMerge(TrackExtended* firstTrackExt, TrackExtended* secondTrackExt) {
  streamlog_out(DEBUG1) << "FullLDCTracking_MarlinTrk::VetoMerge called for " << firstTrackExt << " and "
                        << secondTrackExt << std::endl;

  const float d0First = firstTrackExt->getD0();
  const float z0First = firstTrackExt->getZ0();
  const float omegaFirst = firstTrackExt->getOmega();
  const float tanLFirst = firstTrackExt->getTanLambda();
  const float phi0First = firstTrackExt->getPhi();

  const float d0Second = secondTrackExt->getD0();
  const float z0Second = secondTrackExt->getZ0();
  const float omegaSecond = secondTrackExt->getOmega();
  const float tanLSecond = secondTrackExt->getTanLambda();
  const float phi0Second = secondTrackExt->getPhi();

  HelixClass firstHelix;
  firstHelix.Initialize_Canonical(phi0First, d0First, z0First, omegaFirst, tanLFirst, _bField);
  const float pxFirst = firstHelix.getMomentum()[0];
  const float pyFirst = firstHelix.getMomentum()[1];
  const float pzFirst = firstHelix.getMomentum()[2];

  HelixClass secondHelix;
  secondHelix.Initialize_Canonical(phi0Second, d0Second, z0Second, omegaSecond, tanLSecond, _bField);
  const float pxSecond = secondHelix.getMomentum()[0];
  const float pySecond = secondHelix.getMomentum()[1];
  const float pzSecond = secondHelix.getMomentum()[2];
  const float pFirst = sqrt(pxFirst * pxFirst + pyFirst * pyFirst + pzFirst * pzFirst);
  const float pSecond = sqrt(pxSecond * pxSecond + pySecond * pySecond + pzSecond * pzSecond);

  if (pFirst < _vetoMergeMomentumCut || pSecond < _vetoMergeMomentumCut) {
    streamlog_out(DEBUG1) << "FullLDCTracking_MarlinTrk::VetoMerge do not veto as below momentum cut of 2.5 : pFirst = "
                          << pFirst << " pSecond = " << pSecond << std::endl;
    return false;
  }

  bool veto = false;

  bool testCombinationOnly = true;
  TrackExtended* combinedTrack = CombineTracks(firstTrackExt, secondTrackExt,
                                               _maxAllowedPercentageOfOutliersForTrackCombination, testCombinationOnly);

  if (combinedTrack != NULL) {
    // SJA:FIXME hardcoded cut: here the check is that no more than 7 hits have been rejected in the combined fit.
    if (combinedTrack->getNDF() + 15 < firstTrackExt->getNDF() + secondTrackExt->getNDF() + 5) {
      streamlog_out(DEBUG1) << "FullLDCTracking_MarlinTrk::VetoMerge fails NDF cut " << std::endl;
      veto = true;
    }

    delete combinedTrack->getGroupTracks();
    delete combinedTrack;
  } else {
    streamlog_out(DEBUG1)
        << "FullLDCTracking_MarlinTrk::VetoMerge fails CombineTracks(firstTrackExt,secondTrackExt,true) test"
        << std::endl;
    veto = true;
  }

  if (SegmentRadialOverlap(firstTrackExt, secondTrackExt) > 10) {
    streamlog_out(DEBUG1) << "FullLDCTracking_MarlinTrk::VetoMerge fails SegmentRadialOverlap test " << std::endl;
    veto = true;
  }

  return veto;
}

void FullLDCTracking_MarlinTrk::check(LCEvent*) {}

void FullLDCTracking_MarlinTrk::end() { delete _encoder; }

void FullLDCTracking_MarlinTrk::setupGeom(const dd4hep::Detector& theDetector) {
  double bFieldVec[3];
  theDetector.field().magneticField({0, 0, 0}, bFieldVec); // get the magnetic field vector from DD4hep
  _bField = bFieldVec[2] / dd4hep::tesla;                  // z component at (0,0,0)

  //-- VXD Parameters--
  _nLayersVTX = 0;

  // const std::string nameVertexBarrelDetector = "VertexBarrel";

  try {
    streamlog_out(DEBUG9) << " filling " << m_vtxDetName << " parameters " << '\n';

    dd4hep::DetElement vtxDE = theDetector.detector(m_vtxDetName);
    dd4hep::rec::ZPlanarData* vtx = vtxDE.extension<dd4hep::rec::ZPlanarData>();
    _nLayersVTX = vtx->layers.size();
  } catch (std::runtime_error& e) {
    streamlog_out(ERROR) << " " << m_vtxDetName << " detector Not Present in Compact File" << '\n';
  }

  //-- SIT Parameters--
  _nLayersSIT = 0;

  try {
    streamlog_out(DEBUG9) << " filling " << m_innerBarrelTrackerDetName << " parameters  " << '\n';

    dd4hep::DetElement sitDE = theDetector.detector(m_innerBarrelTrackerDetName);
    dd4hep::rec::ZPlanarData* sit = sitDE.extension<dd4hep::rec::ZPlanarData>();
    _nLayersSIT = sit->layers.size();
  } catch (std::runtime_error& e) {
    streamlog_out(ERROR) << " " << m_innerBarrelTrackerDetName << " detector Not Present in Compact File " << '\n';
  }

  //-- SET Parameters--
  _nLayersSET = 0;

  try {
    streamlog_out(DEBUG9) << " filling " << m_siliconExternalTrackerDetName << " parameters " << '\n';

    dd4hep::DetElement setDE = theDetector.detector(m_siliconExternalTrackerDetName);
    dd4hep::rec::ZPlanarData* set = setDE.extension<dd4hep::rec::ZPlanarData>();
    _nLayersSET = set->layers.size();
  } catch (std::runtime_error& e) {
    streamlog_out(ERROR) << " " << m_siliconExternalTrackerDetName << " detector Not Present in Compact File " << '\n';
  }

  //-- FTD Parameters--
  _petalBasedFTDWithOverlaps = false;
  _nLayersFTD = 0;

  try {
    streamlog_out(DEBUG9) << " filling " << m_forwardTrackerDetName << " parameters " << '\n';

    dd4hep::DetElement ftdDE = theDetector.detector(m_forwardTrackerDetName);
    dd4hep::rec::ZDiskPetalsData* ftd = ftdDE.extension<dd4hep::rec::ZDiskPetalsData>();

    _nLayersFTD = ftd->layers.size();

    for (unsigned int disk = 0; disk < _nLayersFTD; ++disk) {
      _zLayerFTD.push_back(ftd->layers[disk].zPosition + ftd->layers[disk].zOffsetSensitive);
      _zLayerFTD.push_back(ftd->layers[disk].zPosition - ftd->layers[disk].zOffsetSensitive);
      _petalBasedFTDWithOverlaps = true;
    }

    // SJA: Here we increase the size of _nLayersFTD as we are treating the
    _nLayersFTD = _zLayerFTD.size();

  } catch (std::runtime_error& e) {
    streamlog_out(ERROR) << " " << m_forwardTrackerDetName << " detector Not Present in Compact File" << '\n';
  }
}
