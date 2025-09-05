#include "FPCCDSiliconTracking_MarlinTrk.h"

#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <EVENT/TrackerHitZCylinder.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCTOOLS.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>

#include <marlin/Exceptions.h>
#include <marlin/Global.h>

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"

#include "marlin/AIDAProcessor.h"

//---- ROOT -----
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TStopwatch.h"

using namespace lcio;
using namespace marlin;
using namespace MarlinTrk;

using std::abs;
using std::max;
using std::min;

const int FPCCDSiliconTracking_MarlinTrk::_output_track_col_quality_GOOD = 1;
const int FPCCDSiliconTracking_MarlinTrk::_output_track_col_quality_FAIR = 2;
const int FPCCDSiliconTracking_MarlinTrk::_output_track_col_quality_POOR = 3;

const double FPCCDSiliconTracking_MarlinTrk::TWOPI = 2 * M_PI;

FPCCDSiliconTracking_MarlinTrk aFPCCDSiliconTracking_MarlinTrk;

FPCCDSiliconTracking_MarlinTrk::FPCCDSiliconTracking_MarlinTrk() : Processor("FPCCDSiliconTracking_MarlinTrk") {
  _description = "Pattern recognition in silicon trackers";

  _petalBasedFTDWithOverlaps = false;

  // zero triplet counters
  _ntriplets = _ntriplets_good = _ntriplets_2MCP = _ntriplets_3MCP = _ntriplets_1MCP_Bad = _ntriplets_bad = 0;

  std::vector<int> comb;

  comb.push_back(8);
  comb.push_back(6);
  comb.push_back(5);
  comb.push_back(8);
  comb.push_back(6);
  comb.push_back(4);
  comb.push_back(8);
  comb.push_back(5);
  comb.push_back(4);
  comb.push_back(6);
  comb.push_back(5);
  comb.push_back(4);

  registerProcessorParameter("LayerCombinations", "Combinations of Hits in Layers", _Combinations, comb);

  std::vector<int> combFTD;

  combFTD.push_back(13);
  combFTD.push_back(11);
  combFTD.push_back(9);
  combFTD.push_back(13);
  combFTD.push_back(11);
  combFTD.push_back(8);
  combFTD.push_back(13);
  combFTD.push_back(10);
  combFTD.push_back(9);
  combFTD.push_back(13);
  combFTD.push_back(10);
  combFTD.push_back(8);
  combFTD.push_back(12);
  combFTD.push_back(11);
  combFTD.push_back(9);
  combFTD.push_back(12);
  combFTD.push_back(11);
  combFTD.push_back(8);
  combFTD.push_back(12);
  combFTD.push_back(10);
  combFTD.push_back(9);
  combFTD.push_back(12);
  combFTD.push_back(10);
  combFTD.push_back(8);
  combFTD.push_back(11);
  combFTD.push_back(9);
  combFTD.push_back(7);
  combFTD.push_back(11);
  combFTD.push_back(9);
  combFTD.push_back(6);
  combFTD.push_back(11);
  combFTD.push_back(8);
  combFTD.push_back(7);
  combFTD.push_back(11);
  combFTD.push_back(8);
  combFTD.push_back(6);
  combFTD.push_back(10);
  combFTD.push_back(9);
  combFTD.push_back(7);
  combFTD.push_back(10);
  combFTD.push_back(9);
  combFTD.push_back(6);
  combFTD.push_back(10);
  combFTD.push_back(8);
  combFTD.push_back(7);
  combFTD.push_back(9);
  combFTD.push_back(7);
  combFTD.push_back(5);
  combFTD.push_back(9);
  combFTD.push_back(7);
  combFTD.push_back(4);
  combFTD.push_back(9);
  combFTD.push_back(6);
  combFTD.push_back(5);
  combFTD.push_back(9);
  combFTD.push_back(6);
  combFTD.push_back(4);
  combFTD.push_back(8);
  combFTD.push_back(7);
  combFTD.push_back(5);
  combFTD.push_back(8);
  combFTD.push_back(7);
  combFTD.push_back(4);
  combFTD.push_back(8);
  combFTD.push_back(6);
  combFTD.push_back(5);
  combFTD.push_back(8);
  combFTD.push_back(6);
  combFTD.push_back(4);
  combFTD.push_back(7);
  combFTD.push_back(5);
  combFTD.push_back(3);
  combFTD.push_back(7);
  combFTD.push_back(5);
  combFTD.push_back(2);
  combFTD.push_back(7);
  combFTD.push_back(4);
  combFTD.push_back(3);
  combFTD.push_back(7);
  combFTD.push_back(4);
  combFTD.push_back(2);
  combFTD.push_back(6);
  combFTD.push_back(5);
  combFTD.push_back(3);
  combFTD.push_back(6);
  combFTD.push_back(5);
  combFTD.push_back(2);
  combFTD.push_back(6);
  combFTD.push_back(4);
  combFTD.push_back(3);
  combFTD.push_back(6);
  combFTD.push_back(4);
  combFTD.push_back(2);
  combFTD.push_back(5);
  combFTD.push_back(3);
  combFTD.push_back(1);
  combFTD.push_back(5);
  combFTD.push_back(3);
  combFTD.push_back(0);
  combFTD.push_back(5);
  combFTD.push_back(2);
  combFTD.push_back(1);
  combFTD.push_back(5);
  combFTD.push_back(2);
  combFTD.push_back(0);
  combFTD.push_back(4);
  combFTD.push_back(3);
  combFTD.push_back(1);
  combFTD.push_back(4);
  combFTD.push_back(3);
  combFTD.push_back(0);
  combFTD.push_back(4);
  combFTD.push_back(2);
  combFTD.push_back(1);
  combFTD.push_back(4);
  combFTD.push_back(2);
  combFTD.push_back(0);

  registerProcessorParameter("LayerCombinationsFTD", "Combinations of Hits in FTD", _CombinationsFTD, combFTD);

  registerProcessorParameter("NDivisionsInPhiFTD", "Number of divisions in Phi for FTD", _nPhiFTD, int(30));

  registerProcessorParameter("NDivisionsInPhi", "Number of divisions in Phi", _nDivisionsInPhi,
                             int(1440)); // 1 sector is equivalent to 0.25 deg.

  registerProcessorParameter("NDivisionsInTheta", "Number of divisions in Theta", _nDivisionsInTheta,
                             int(640)); // 1 sector, delta(cosTheta), is equivalent to 2/640.

  registerProcessorParameter("SearchWindowForTripletInTheta", "num of theta-sectors used in making triplets", _sw_theta,
                             int(16));

  // Input Collections for debug
  // ^^^^^^^^^^^^^^^^^
  registerInputCollection(LCIO::SIMTRACKERHIT, "SimTrackerHit of VXD (VXDCollection)",
                          "This is used to calculate picking efficiency of a track.", _colNameVXDSimHit,
                          std::string("VXDCollection"));

  registerInputCollection(LCIO::SIMTRACKERHIT, "SimTrackerHit of SIT (SITCollection)",
                          "This is used to calculate picking efficiency of a track.", _colNameSITSimHit,
                          std::string("SITCollection"));

  registerInputCollection(LCIO::SIMTRACKERHIT, "SimTrackerHit of FTD_STRIP (FTD_STRIPCollection)",
                          "This is used to calculate picking efficiency of a track.", _colNameFTDspSimHit,
                          std::string("FTD_STRIPCollection"));

  registerInputCollection(LCIO::SIMTRACKERHIT, "SimTrackerHit of FTD_PIXEL (FTD_PIXELCollection)",
                          "This is used to calculate picking efficiency of a track.", _colNameFTDpixSimHit,
                          std::string("FTD_PIXELCollection"));

  // Input Collections
  registerInputCollection(LCIO::TRACKERHITPLANE, "VTXHitCollectionName", "VTX Hit Collection Name", _VTXHitCollection,
                          std::string("VXDTrackerHits"));

  registerInputCollection(LCIO::TRACKERHITPLANE, "FTDPixelHitCollectionName", "FTD Pixel Hit Collection Name",
                          _FTDPixelHitCollection, std::string("FTDPixelTrackerHits"));

  registerInputCollection(LCIO::TRACKERHIT, "FTDSpacePointCollectionName", "FTD FTDSpacePoint Collection Name",
                          _FTDSpacePointCollection, std::string("FTDSpacePoints"));

  registerInputCollection(LCIO::TRACKERHIT, "SITHitCollectionName", "SIT Hit Collection Name", _SITHitCollection,
                          std::string("SITSpacePoints"));

  // Output Collections
  // ^^^^^^^^^^^^^^^^^^
  registerOutputCollection(LCIO::TRACK, "SiTrackCollectionName", "Silicon track Collection Name", _siTrkCollection,
                           std::string("SiTracks"));

  // Steering parameters
  // ^^^^^^^^^^^^^^^^^^^
  registerProcessorParameter("Chi2WRphiTriplet", "Chi2WRphiTriplet", _chi2WRPhiTriplet, float(1.));

  registerProcessorParameter("Chi2WZTriplet", "Chi2WZTriplet", _chi2WZTriplet, float(0.5));

  registerProcessorParameter("Chi2FitCut", "Chi2 Fit Cut", _chi2FitCut, float(120.0));

  registerProcessorParameter("Chi2FitCut_kalman", "Chi2 Fit Cut for BuildTrack_KalFit", _chi2FitCut_kalman,
                             float(120.0));

  registerProcessorParameter("AngleCutForMerging", "Angle Cut For Merging", _angleCutForMerging, float(0.1));

  registerProcessorParameter("MinDistCutAttachForFTD", "MinDistCutAttachForFTD", _minDistCutAttachForFTD, float(2.5));

  registerProcessorParameter("MinDistCutAttachForVXD", "MinDistCutAttachForVXD", _minDistCutAttachForVXD,
                             float(0.5)); // unit: mm

  registerProcessorParameter("MinMissAddition", "Allowed max number of MinMissAdditions in BuildTrack",
                             _minMissAddition, int(2));

  registerProcessorParameter("FudgeFactorForSIT_rphi_SpatialResolution",
                             "SIT's spatial resolution is not calculated precisely. (default : 2.0)",
                             _fudgeFactorForSITsr_rphi, float(2.0));

  registerProcessorParameter("FudgeFactorForSIT_z_SpatialResolution",
                             "SIT's spatial resolution is not calculated precisely. (default : 1.0)",
                             _fudgeFactorForSITsr_z, float(1.0));

  registerProcessorParameter("MinLayerToAttach", "MinLayerToAttach", _minimalLayerToAttach, int(-1));

  registerProcessorParameter("CutOnD0", "cut on D0 for tracks", _cutOnD0, float(60.0));

  registerProcessorParameter("CutOnZ0", "cut on Z0 for tracks", _cutOnZ0, float(100.0));

  registerProcessorParameter("CutOnPt_For_VXD+SIT_section", "cut on Pt for TestTriplet in VXD+SIT section", _cutOnPtVXD,
                             double(0.18));

  registerProcessorParameter("CutOnPt_For_FTD_section ", "cut on Pt for TestTriplet in FTD section", _cutOnPtFTD,
                             double(0.05));

  registerProcessorParameter("MinimalHits", "minimal hits (default 3)", _minimalHits, int(3));

  registerProcessorParameter(
      "NHitsChi2", "Maximal number of hits for which a track with n hits is better than one with n-1hits. (defaut 8)",
      _nHitsChi2, int(8));

  registerProcessorParameter(
      "MaxHitsPerSector",
      "Maximal number of hits allowed in one theta-phi sector in FTD (In VXD and SIT, this is not implemented)",
      _max_hits_per_sector, int(100));

  registerProcessorParameter("AttachRemainingHitsForVXD", "1:Fast, 2:Slow, 3:VeryFast, 0:skip this process", _attachVXD,
                             int(2)); // Very Fast is underconstruction.

  registerProcessorParameter("AttachRemainingHitsForFTD", "1:Fast, 2:Slow, 0:skip this process", _attachFTD, int(2));

  registerProcessorParameter("UseSIT", "Use SIT", _useSIT, int(1));

  registerProcessorParameter("UseFTD", "Use FTD", _useFTD, int(1));

  registerProcessorParameter("NDivisions_phiRangeForBuildTrackForHighPt",
                             "set width of phi(deg) value [0,360].ex: 360deg/80 = 4.5 (This value emulates original "
                             "SiliconTracking_MarlinTrk by default) ",
                             _phiRangeForBuildTrackForHighPt, double(4.5));

  registerProcessorParameter("NDivisions_cosThetaRangeForBuildTrackForHighPt",
                             "set width of cosTheta value [0,2].ex: 2/80 = 0.025 (This value emulates original "
                             "SiliconTracking_MarlinTrk by default)",
                             _cosThetaRangeForBuildTrackForHighPt, double(0.025));

  registerProcessorParameter("UseClusterRejection",
                             "Use Cluster Rejection (for FPCCD VXD in dense pair-BG env. option)", _useClusterRejection,
                             bool(false));

  registerProcessorParameter("minDotOf2Clusters",
                             "minimum dot of 2Clusters (for FPCCD VXD in dense pair-BG env. option)",
                             _minDotOf2Clusters, float(0.4));

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

  registerProcessorParameter("MaxChi2PerHit",
                             "Maximum Chi-squared value allowed when assigning a hit to a track (for KalTest)",
                             _maxChi2PerHit, double(1.e2));

  registerProcessorParameter(
      "MaxChi2PerHit2nd",
      "Maximum Chi-squared value allowed when assigning a hit to a track (especially for KalFit in this code)",
      _maxChi2PerHit2nd, double(1.e2));

  registerProcessorParameter("CheckForDelta", "Check for Delta rays hits in hit-to-track assignment", _checkForDelta,
                             int(1));

  registerProcessorParameter("MinDistToDelta", "Minimal distance of track hit to the delta electron hit",
                             _minDistToDelta, float(0.25));

  registerProcessorParameter("MultipleScatteringOn", "Use MultipleScattering in Fit", _MSOn, bool(true));

  registerProcessorParameter("nSigma_for_Build_phiTrack",
                             "In extrapolation, search area in phi is determined by using sigma of d0, so you can "
                             "enlarge this area by this value. ",
                             _nSigmaBuild_phi, double(5.0));

  registerProcessorParameter("nSigma_for_Build_thetaTrack",
                             "In extrapolation, search area in phi is determined by using sigma of d0, so you can "
                             "enlarge this area by this value. ",
                             _nSigmaBuild_theta, double(5.0));

  registerProcessorParameter("EnergyLossOn", "Use Energy Loss in Fit", _ElossOn, bool(true));

  registerProcessorParameter("SmoothOn", "Smooth All Mesurement Sites in Fit", _SmoothOn, bool(false));

  registerProcessorParameter("SafetyPhi-Range_ratio_ForTripletSearchAndExtrapolation", "Range <= (1 + this)*Range",
                             _safetyPhiRange_ratio,
                             float(0.10)); // old 0.3

  registerProcessorParameter("SafetyPhi-Range_fix_ForTripletSearchAndExtrapolation", "Range += this",
                             _safetyPhiRange_fix,
                             int(1)); // old 0.1

  registerInputCollection("MCParticle", "MCParticleCollectionName", "Name of the MCParticle input collection",
                          _colNameMCParticles, std::string("MCParticle"));

  registerInputCollection("LCRelation", "VXDTrackerHitRelations", "Name of the LCRelation input collection",
                          _colNameVXDTrackerHitRelations, std::string("VXDTrackerHitRelations"));

  registerInputCollection("LCRelation", "SITSpacePointRelations", "Name of the LCRelation input collection",
                          _colNameSITSpacePointRelations, std::string("SITSpacePointRelations"));

  registerInputCollection("LCRelation", "FTDSpacePointRelations", "Name of the LCRelation input collection",
                          _colNameFTDSpacePointRelations, std::string("FTDSpacePointRelations"));

  registerInputCollection("LCRelation", "FTDPixelTrackerHitRelations", "Name of the LCRelation input collection",
                          _colNameFTDPixelTrackerHitRelations, std::string("FTDPixelTrackerHitRelations"));

  FloatVec PixelSizeVec; // for FPCCD VXD
  PixelSizeVec.push_back(0.005);
  PixelSizeVec.push_back(0.005);
  PixelSizeVec.push_back(0.010);
  PixelSizeVec.push_back(0.010);
  PixelSizeVec.push_back(0.010);
  PixelSizeVec.push_back(0.010);

  registerProcessorParameter("FPCCD: Each_FPCCD_pixelSize(mm)",
                             "Each ladder's FPCCD Pixel size(unit:mm) (for cluster rejection for FPCCD)", _pixelSizeVec,
                             PixelSizeVec);

  registerProcessorParameter("FPCCD: PixelHeight", "Pixel Height(mm) (for cluster rejection for FPCCD)", _pixelheight,
                             float(0.015));

  registerProcessorParameter("fudgePhiRange",
                             "fudgePhiRange for the iterative determination of search range in the extrapolation",
                             _fudgePhiRange, int(2));

  registerProcessorParameter("fudgeThetaRange",
                             "fudgeThetaRange for the iterative determination of search range in the extrapolation",
                             _fudgeThetaRange, int(2));

  registerProcessorParameter("mydebug", "mydebug (code debuger for Mori)", _mydebug, bool(false));

  registerProcessorParameter("mydebugKalFit", "mydebugKalFit (code debuger for Mori)", _mydebugKalFit, bool(false));

  registerProcessorParameter("mydebugIntersection", "mydebugIntersection (code debuger for Mori)", _mydebugIntersection,
                             bool(false));

  registerProcessorParameter("stopwatch", "stopwatch (code debuger for Mori)", _stopwatch, bool(false));

  registerProcessorParameter("mydebugstopwatch2", "mydebugstopwatch2 (code debuger for Mori)", _mydebugstopwatch2,
                             bool(false));

  registerProcessorParameter("keepCandidate", "Option for AttachRemainingVTXHitsVeryFast (please fix false for now)",
                             _keepCandidate, bool(false));

  _output_track_col_quality = _output_track_col_quality_GOOD;
}

void FPCCDSiliconTracking_MarlinTrk::init() {
  _nRun = -1;
  _nEvt = 0;

  _encoder = new UTIL::BitField64(lcio::LCTrackerCellID::encoding_string());

  _fastfitter = new MarlinTrk::HelixFit();
  _moriUtil = new moriUTIL();
  _purityUtil = new GetPurityUtil();

  printParameters();

  // this creates a directory for this processor ....
  _trksystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", 0, "");

  if (_trksystem == 0) {
    throw EVENT::Exception(std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("KalTest"));
  }

  _trksystem->setOption(IMarlinTrkSystem::CFG::useQMS, _MSOn);
  _trksystem->setOption(IMarlinTrkSystem::CFG::usedEdx, _ElossOn);
  _trksystem->setOption(IMarlinTrkSystem::CFG::useSmoothing, _SmoothOn);
  _trksystem->init();

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  this->setupGeom(theDetector);

  if (_useSIT == 0)
    _nLayers = _nLayersVTX;
  else
    _nLayers = _nLayersVTX + _nLayersSIT;

  // initialise the container to have separate vectors for up to _nHitsChi2 hits.
  _tracksWithNHitsContainer.resize(_nHitsChi2);

  _dPhi = TWOPI / _nDivisionsInPhi;
  _dTheta = 2.0 / _nDivisionsInTheta;
  _dPhiFTD = TWOPI / _nPhiFTD;
  double cutOnR_VXD = _cutOnPtVXD / (0.299792458 * _bField);
  cutOnR_VXD = 1000. * cutOnR_VXD;
  _cutOnOmegaVXD = 1 / cutOnR_VXD;
  double cutOnR_FTD = _cutOnPtFTD / (0.299792458 * _bField);
  cutOnR_FTD = 1000. * cutOnR_FTD;
  _cutOnOmegaFTD = 1 / cutOnR_FTD;

  InitVXDGeometry(theDetector);
  if (_useSIT == 1)
    InitSITGeometry(theDetector);

  _output_track_col_quality = 0;

// Range determination for TestTriplet
#if 0
  //new version (under construction)
  int ncomb = _Combinations.size();
  for(int i = 0; i < ncomb; i += 3 ){
    std::vector<int> layers(3);
    layers[0] = _Combinations[i];
    layers[1] = _Combinations[i+1];
    layers[2] = _Combinations[i+2];
    std::vector<double> phiDiff; //for middle and inner layer
    getNeededPhiSectorsVer2(_cutOnPtVXD, layers, phiDiff);
    std::vector<int> phiRange(2);
    phiRange[0] = int(std::abs(phiDiff[0]/_dPhi)*(1.0 + _safetyPhiRange_ratio)) + _safetyPhiRange_fix;
    phiRange[1] = int(std::abs(phiDiff[1]/_dPhi)*(1.0 + _safetyPhiRange_ratio)) + _safetyPhiRange_fix;
    _phiRangeForTripletVer2.insert(std::make_pair(layers,phiRange));
  }
#else
  // old version (for now, please use this)
  int ncomb = _Combinations.size();
  for (int i = 0; i < ncomb; i += 3) {
    std::pair<int, int> lyp(_Combinations[i], _Combinations[i + 2]);
    double phiDiff = getNeededPhiSectors(_cutOnPtVXD, _Combinations[i], _Combinations[i + 2]);
    int phiRange = int(std::abs(phiDiff / _dPhi) * (1.0 + _safetyPhiRange_ratio)) + _safetyPhiRange_fix;
    // check
    if (_mydebug)
      std::cout << "phiRange for Triplet (layer " << _Combinations[i] << "  " << _Combinations[i + 2]
                << " ) : " << phiRange << std::endl;
    _phiRangeForTriplet.insert(std::make_pair(lyp, phiRange));
  }

#endif
}

void FPCCDSiliconTracking_MarlinTrk::processRunHeader(LCRunHeader*) {
  _nRun++;
  _nEvt = 0;

  streamlog_out(MESSAGE) << "FPCCDSiliconTracking_MarlinTrk ---> new run : run number = " << _nRun << std::endl;
}

void FPCCDSiliconTracking_MarlinTrk::processEvent(LCEvent* evt) {
  if (_stopwatch) {
    _timer.reset();
    _timer.inAnEvt.Start();
  }

  // set the correct configuration for the tracking system for this event
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(_trksystem, _MSOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(_trksystem, _ElossOn);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(_trksystem, _SmoothOn);

  // check
  std::cout << "=========== evt num : " << _nEvt << "  =================" << std::endl;

  _current_event = evt;

  _output_track_col_quality = _output_track_col_quality_GOOD;

  // zero triplet counters
  _ntriplets = _ntriplets_good = _ntriplets_2MCP = _ntriplets_3MCP = _ntriplets_1MCP_Bad = _ntriplets_bad = 0;

  // Clearing the working containers from the previous event
  // FIXME: partly done at the end of the event, in CleanUp. Make it consistent.
  _tracksWithNHitsContainer.clear();
  _trackImplVec.clear();

  _colTrackerHits.clear();
  _colNamesTrackerHits.clear();

  streamlog_out(DEBUG4) << "FPCCDSiliconTracking_MarlinTrk -> run = " << _nRun << "  event = " << _nEvt << std::endl;

  if (_stopwatch)
    _timer.InitialiseVTX.Start();
  int successVTX = InitialiseVTX(evt);
  if (_stopwatch)
    _timer.InitialiseVTX.Stop();

  if (_stopwatch)
    _timer.InitialiseFTD.Start();
  int successFTD = InitialiseFTD(evt);
  if (_stopwatch)
    _timer.InitialiseFTD.Stop();

  if (_stopwatch)
    _timer.ProcessOneSector.Start();
  if (successVTX == 1) {
    streamlog_out(DEBUG1) << "      phi          theta        layer      nh o :   m :   i  :: o*m*i " << std::endl;
    for (int iPhi = 0; iPhi < _nDivisionsInPhi; ++iPhi) {
      for (int iTheta = 0; iTheta < _nDivisionsInTheta; ++iTheta) {
        ProcessOneSector(iPhi, iTheta); // Process one VXD sector
      }
    }
    streamlog_out(DEBUG4) << "End of Processing VXD and SIT sectors" << std::endl;
    if (_mydebugstopwatch2) {
      printf("_timer2Triplet :  RT=%.3f s, CPU=%.3f s, count=%d \n", _timer2Triplet.RealTime(),
             _timer2Triplet.CpuTime(), _timer2Triplet.Counter());
      printf("_timer2Build :  RT=%.3f s, CPU=%.3f s, count=%d \n", _timer2Build.RealTime(), _timer2Build.CpuTime(),
             _timer2Build.Counter());
    }
  }
  if (_stopwatch)
    _timer.ProcessOneSector.Stop();

  if (_stopwatch)
    _timer.TrackingInFTD.Start();
  if (successFTD == 1) {
    streamlog_out(DEBUG1) << "      phi          side        layer      nh o :   m :   i  :: o*m*i " << std::endl;
    TrackingInFTD(); // Perform tracking in the FTD
    streamlog_out(DEBUG4) << "End of Processing FTD sectors" << std::endl;
  }
  if (_stopwatch)
    _timer.TrackingInFTD.Stop();

  if (successVTX == 1 || successFTD == 1) {
    if (_stopwatch)
      _timer.Sorting.Start();
    for (int nHits = _nHitsChi2; nHits >= _minimalHits; nHits--) {
      Sorting(_tracksWithNHitsContainer.getTracksWithNHitsVec(nHits));
    }
    streamlog_out(DEBUG4) << "End of Sorting " << std::endl;
    if (_stopwatch)
      _timer.Sorting.Stop();

    if (_stopwatch)
      _timer.CreateTrack.Start();
    for (int nHits = _nHitsChi2; nHits >= _minimalHits; nHits--) {
      TrackExtendedVec& tracksWithNHits = _tracksWithNHitsContainer.getTracksWithNHitsVec(nHits);
      for (TrackExtendedVec::iterator trackIter = tracksWithNHits.begin(); trackIter < tracksWithNHits.end();
           trackIter++) {
        CreateTrack(*trackIter);
      }
      streamlog_out(DEBUG4) << "End of creating " << nHits << " hits tracks " << std::endl;
    }
    if (_stopwatch)
      _timer.CreateTrack.Stop();

    if (_stopwatch)
      _timer.AttachRemainingVXD.Start();
    if (_attachVXD == 1)
      AttachRemainingVTXHitsFast();
    else if (_attachVXD == 2)
      AttachRemainingVTXHitsSlow();
    else if (_attachVXD == 3)
      AttachRemainingVTXHitsVeryFast();
    if (_stopwatch)
      _timer.AttachRemainingVXD.Stop();

    if (_stopwatch)
      _timer.AttachRemainingFTD.Start();
    if (_attachFTD == 1)
      AttachRemainingFTDHitsFast();
    else if (_attachFTD == 2)
      AttachRemainingFTDHitsSlow();
    if (_stopwatch)
      _timer.AttachRemainingFTD.Stop();

    streamlog_out(DEBUG4) << "End of picking up remaining hits " << std::endl;

    LCCollectionVec* trkCol = new LCCollectionVec(LCIO::TRACK);
    // if we want to point back to the hits we need to set the flag
    LCFlagImpl trkFlag(0);
    trkFlag.setBit(LCIO::TRBIT_HITS);
    trkCol->setFlag(trkFlag.getFlag());

    LCCollectionVec* relCol = NULL;

    if (_stopwatch)
      _timer.FinalRefit.Start();
    FinalRefit(trkCol, relCol);
    if (_stopwatch)
      _timer.FinalRefit.Stop();

    // set the quality of the output collection
    switch (_output_track_col_quality) {
    case _output_track_col_quality_FAIR:
      trkCol->parameters().setValue("QualityCode", "Fair");
      break;

    case _output_track_col_quality_POOR:
      trkCol->parameters().setValue("QualityCode", "Poor");
      break;

    default:
      trkCol->parameters().setValue("QualityCode", "Good");
      break;
    }

    evt->addCollection(trkCol, _siTrkCollection.c_str());
  }

  CleanUp();
  streamlog_out(DEBUG4) << "Event is done " << std::endl;

  if (_stopwatch) {
    _timer.inAnEvt.Stop();
    _timer.cout();
  }

  _nEvt++;
}

void FPCCDSiliconTracking_MarlinTrk::CleanUp() {
  _tracksWithNHitsContainer.clear();

  for (int il = 0; il < _nLayers; ++il) {
    for (int ip = 0; ip < _nDivisionsInPhi; ++ip) {
      for (int it = 0; it < _nDivisionsInTheta; ++it) {
        unsigned int iCode = il + _nLayers * ip + _nLayers * _nDivisionsInPhi * it;

        if (iCode >= _sectors.size()) {
          std::cout << "iCode index out of range: iCode =   " << iCode << " _sectors.size() = " << _sectors.size()
                    << " exit(1) called from file " << __FILE__ << " line " << __LINE__ << std::endl;
          exit(1);
        }

        TrackerHitExtendedVec& hitVec = _sectors[iCode];
        int nH = int(hitVec.size());
        for (int iH = 0; iH < nH; ++iH) {
          TrackerHitExtended* hit = hitVec[iH];
          delete hit;
        }
      }
    }
  }

  for (int iS = 0; iS < 2; ++iS) {
    for (unsigned int layer = 0; layer < _nlayersFTD; ++layer) {
      for (int ip = 0; ip < _nPhiFTD; ++ip) {
        unsigned int iCode = iS + 2 * layer + 2 * _nlayersFTD * ip;

        if (iCode >= _sectorsFTD.size()) {
          std::cout << "iCode index out of range: iCode =   " << iCode << " _sectorsFTD.size() = " << _sectorsFTD.size()
                    << " exit(1) called from file " << __FILE__ << " line " << __LINE__ << std::endl;
          exit(1);
        }

        TrackerHitExtendedVec& hitVec = _sectorsFTD[iCode];
        int nH = int(hitVec.size());
        for (int iH = 0; iH < nH; ++iH) {
          TrackerHitExtended* hit = hitVec[iH];
          delete hit;
        }
      }
    }
  }
}

int FPCCDSiliconTracking_MarlinTrk::InitialiseFTD(LCEvent* evt) {
  int success = 1;

  _nTotalFTDHits = 0;
  _sectorsFTD.clear();
  _sectorsFTD.resize(2 * _nlayersFTD * _nPhiFTD);

  // Reading in FTD Pixel Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    LCCollection* hitCollection = evt->getCollection(_FTDPixelHitCollection.c_str());

    _colNamesTrackerHits[hitCollection] = _FTDPixelHitCollection;
    _colTrackerHits.push_back(hitCollection);

    int nelem = hitCollection->getNumberOfElements();

    streamlog_out(DEBUG4) << "Number of FTD Pixel Hits = " << nelem << std::endl;
    _nTotalFTDHits = nelem;

    for (int ielem = 0; ielem < nelem; ++ielem) {
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem));

      TrackerHitExtended* hitExt = new TrackerHitExtended(hit);

      Vector3D U(1.0, hit->getU()[1], hit->getU()[0], Vector3D::spherical);
      Vector3D V(1.0, hit->getV()[1], hit->getV()[0], Vector3D::spherical);
      Vector3D Z(0.0, 0.0, 1.0);

      const float eps = 1.0e-07;
      // V must be the global z axis
      if (fabs(V.dot(Z)) > eps) {
        streamlog_out(ERROR) << "FPCCDSiliconTracking_MarlinTrk: VXD Hit measurment vectors V is not in the global X-Y "
                                "plane. \n\n  exit(1) called from file "
                             << __FILE__ << " and line " << __LINE__ << std::endl;
        exit(1);
      }

      if (fabs(U.dot(Z)) > eps) {
        streamlog_out(ERROR) << "FPCCDSiliconTracking_MarlinTrk: VXD Hit measurment vectors U is not in the global X-Y "
                                "plane. \n\n exit(1) called from file "
                             << __FILE__ << " and line " << __LINE__ << std::endl;
        exit(1);
      }

      // SJA:FIXME Here dU and dV are almost certainly dX and dY ... should test ...
      double point_res_rphi = sqrt(hit->getdU() * hit->getdU() + hit->getdV() * hit->getdV());

      hitExt->setResolutionRPhi(point_res_rphi);

      // SJA:FIXME why is this needed?
      hitExt->setResolutionZ(0.1);

      // type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try and
      // catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));

      double pos[3];

      for (int i = 0; i < 3; ++i) {
        pos[i] = hit->getPosition()[i];
      }

      double Phi = atan2(pos[1], pos[0]);
      if (Phi < 0.)
        Phi = Phi + TWOPI;

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

      if (layer >= _nlayersFTD) {
        streamlog_out(ERROR)
            << "FPCCDSiliconTracking_MarlinTrk => fatal error in FTD : layer is outside allowed range : " << layer
            << " number of layers = " << _nlayersFTD << std::endl;
        exit(1);
      }

      int iPhi = int(Phi / _dPhiFTD);

      int side = getSideID(hit);
      int iSemiSphere = 0;

      if (side > 0)
        iSemiSphere = 1;

      int iCode = iSemiSphere + 2 * layer + 2 * _nlayersFTD * iPhi;
      _sectorsFTD[iCode].push_back(hitExt);

      streamlog_out(DEBUG1) << " FTD Pixel Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi "
                            << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  iPhi = " << iPhi
                            << " iSemiSphere " << iSemiSphere << " iCode = " << iCode << "  layer = " << layer
                            << std::endl;
    }
  } catch (DataNotAvailableException& e) {
    success = 0;
  }

  // Reading out FTD SpacePoint Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    LCCollection* hitCollection = evt->getCollection(_FTDSpacePointCollection.c_str());

    _colNamesTrackerHits[hitCollection] = _FTDSpacePointCollection;
    _colTrackerHits.push_back(hitCollection);

    int nelem = hitCollection->getNumberOfElements();

    streamlog_out(DEBUG4) << "Number of FTD SpacePoints = " << nelem << std::endl;
    _nTotalFTDHits += nelem;

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

      double pos[3];

      for (int i = 0; i < 3; ++i) {
        pos[i] = hit->getPosition()[i];
      }

      double Phi = atan2(pos[1], pos[0]);
      if (Phi < 0.)
        Phi = Phi + TWOPI;

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

      if (layer >= _nlayersFTD) {
        streamlog_out(ERROR)
            << "FPCCDSiliconTracking_MarlinTrk => fatal error in FTD : layer is outside allowed range : " << layer
            << " number of layers = " << _nlayersFTD << std::endl;
        exit(1);
      }

      int iPhi = int(Phi / _dPhiFTD);

      int side = getSideID(hit);
      int iSemiSphere = 0;

      if (side > 0)
        iSemiSphere = 1;

      int iCode = iSemiSphere + 2 * layer + 2 * _nlayersFTD * iPhi;
      _sectorsFTD[iCode].push_back(hitExt);

      streamlog_out(DEBUG1) << " FTD SpacePoint Hit added : @ " << pos[0] << " " << pos[1] << " " << pos[2] << " drphi "
                            << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ() << "  iPhi = " << iPhi
                            << " iSemiSphere " << iSemiSphere << " iCode = " << iCode << "  layer = " << layer
                            << std::endl;
    }
  } catch (DataNotAvailableException& e) {
    success = 0;
  }

  for (unsigned i = 0; i < _sectorsFTD.size(); ++i) {
    int nhits = _sectorsFTD[i].size();
    if (nhits != 0)
      streamlog_out(DEBUG1) << " Number of Hits in FTD Sector " << i << " = " << _sectorsFTD[i].size() << std::endl;
    if (nhits > _max_hits_per_sector) {
      for (unsigned ihit = 0; ihit < _sectorsFTD[i].size(); ++ihit) {
        delete _sectorsFTD[i][ihit];
      }
      _sectorsFTD[i].clear();
      if (nhits != 0)
        streamlog_out(ERROR) << " ### EVENT " << evt->getEventNumber() << " :: RUN " << evt->getRunNumber()
                             << " \n ### Number of Hits in FTD Sector " << i << " = " << nhits << " : Limit is set to "
                             << _max_hits_per_sector
                             << " : This sector will be dropped from track search, and QualityCode set to \"Poor\" "
                             << std::endl;

      _output_track_col_quality = _output_track_col_quality_POOR;
    }
  }

  return success;
}

int FPCCDSiliconTracking_MarlinTrk::InitialiseVTX(LCEvent* evt) {
  _nTotalVTXHits = 0;
  _nTotalSITHits = 0;
  _sectors.clear();
  _sectors.resize(_nLayers * _nDivisionsInPhi * _nDivisionsInTheta);

  // Reading out VTX Hits Collection
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  try {
    LCCollection* hitCollection = evt->getCollection(_VTXHitCollection.c_str());

    _colNamesTrackerHits[hitCollection] = _VTXHitCollection;
    _colTrackerHits.push_back(hitCollection);

    int nelem = hitCollection->getNumberOfElements();

    streamlog_out(DEBUG4) << "Number of VTX hits = " << nelem << std::endl;
    _nTotalVTXHits = nelem;

    for (int ielem = 0; ielem < nelem; ++ielem) {
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem));

      Vector3D U(1.0, hit->getU()[1], hit->getU()[0], Vector3D::spherical);
      Vector3D V(1.0, hit->getV()[1], hit->getV()[0], Vector3D::spherical);
      Vector3D Z(0.0, 0.0, 1.0);

      const float eps = 1.0e-07;
      // V must be the global z axis
      if (fabs(1.0 - V.dot(Z)) > eps) {
        streamlog_out(ERROR) << "FPCCDSiliconTracking_MarlinTrk: VXD Hit measurment vectors V is not equal to the "
                                "global Z axis. \n\n  exit(1) called from file "
                             << __FILE__ << " and line " << __LINE__ << std::endl;
        exit(1);
      }

      if (fabs(U.dot(Z)) > eps) {
        streamlog_out(ERROR) << "FPCCDSiliconTracking_MarlinTrk: VXD Hit measurment vectors U is not in the global X-Y "
                                "plane. \n\n exit(1) called from file "
                             << __FILE__ << " and line " << __LINE__ << std::endl;
        exit(1);
      }

      TrackerHitExtended* hitExt = new TrackerHitExtended(hit);

      // SJA:FIXME: just use planar res for now
      hitExt->setResolutionRPhi(hit->getdU());
      hitExt->setResolutionZ(hit->getdV());

      // set type is now only used in one place where it is set to 0 to reject hits from a fit, set to INT_MAX to try
      // and catch any missuse
      hitExt->setType(int(INT_MAX));
      // det is no longer used set to INT_MAX to try and catch any missuse
      hitExt->setDet(int(INT_MAX));

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
      int layer = getLayerID(hit);

      if (layer < 0 || layer >= _nLayers) {
        streamlog_out(ERROR)
            << "FPCCDSiliconTracking_MarlinTrk => fatal error in VTX : layer is outside allowed range : " << layer
            << std::endl;
        exit(1);
      }

      int iPhi = int(Phi / _dPhi);
      int iTheta = int((cosTheta + double(1.0)) / _dTheta);
      int iCode = layer + _nLayers * iPhi + _nLayers * _nDivisionsInPhi * iTheta;
      _sectors[iCode].push_back(hitExt);

      streamlog_out(DEBUG1) << " VXD Hit " << hit->id() << " added : @ " << pos[0] << " " << pos[1] << " " << pos[2]
                            << " drphi " << hitExt->getResolutionRPhi() << " dz " << hitExt->getResolutionZ()
                            << "  iPhi = " << iPhi << " iTheta " << iTheta << " iCode = " << iCode
                            << "  layer = " << layer << std::endl;
    }
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << " collection not found : " << _VTXHitCollection.c_str() << std::endl;
  }

  if (_useSIT > 0) {
    try {
      LCCollection* hitCollection = evt->getCollection(_SITHitCollection.c_str());

      _colNamesTrackerHits[hitCollection] = _SITHitCollection;
      _colTrackerHits.push_back(hitCollection);

      int nelem = hitCollection->getNumberOfElements();

      streamlog_out(DEBUG4) << "Number of SIT hits = " << nelem << std::endl;
      _nTotalSITHits = nelem;

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
        //    Instead ii)  ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT (TrackerHit) iii) TrackerHitPlane (Two dimentional)
        //    iv)  TrackerHitZCylinder
        //    v)   Must be standard TrackerHit

        trkhit = dynamic_cast<TrackerHit*>(hitCollection->getElementAt(ielem));

        int layer = getLayerID(trkhit);

        // VXD and SIT are treated as one system so SIT layers start from _nLayersVTX
        layer = layer + _nLayersVTX;

        if (layer < 0 || layer >= _nLayers) {
          streamlog_out(ERROR)
              << "FPCCDSiliconTracking_MarlinTrk => fatal error in SIT : layer is outside allowed range : " << layer
              << std::endl;
          exit(1);
        }

        // first check that we have not been given 1D hits by mistake, as they won't work here
        if (BitSet32(trkhit->getType())[UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL]) {
          streamlog_out(ERROR)
              << "FPCCDSiliconTracking_MarlinTrk: SIT Hit cannot be of type UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL "
                 "COMPOSITE SPACEPOINTS must be use instead. \n\n  exit(1) called from file "
              << __FILE__ << " and line " << __LINE__ << std::endl;
          exit(1);

        }
        // most likely case: COMPOSITE_SPACEPOINT hits formed from stereo strip hits
        else if (BitSet32(trkhit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) {
          // SJA:FIXME: fudge for now by a factor of two and ignore covariance
          // MORI:FIXME: fudge factor is now changeable in steering file.
          drphi = _fudgeFactorForSITsr_rphi * sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]);
          dz = _fudgeFactorForSITsr_z * sqrt(trkhit->getCovMatrix()[5]);

        }
        // or a PIXEL based SIT, using 2D TrackerHitPlane like the VXD above
        else if ((trkhit_P = dynamic_cast<TrackerHitPlane*>(hitCollection->getElementAt(ielem)))) {
          // first we need to check if the measurement vectors are aligned with the global coordinates
          Vector3D U(1.0, trkhit_P->getU()[1], trkhit_P->getU()[0], Vector3D::spherical);
          Vector3D V(1.0, trkhit_P->getV()[1], trkhit_P->getV()[0], Vector3D::spherical);
          Vector3D Z(0.0, 0.0, 1.0);

          const float eps = 1.0e-07;
          // V must be the global z axis
          if (fabs(1.0 - V.dot(Z)) > eps) {
            streamlog_out(ERROR) << "FPCCDSiliconTracking_MarlinTrk: PIXEL SIT Hit measurment vectors V is not equal "
                                    "to the global Z axis. \n\n  exit(1) called from file "
                                 << __FILE__ << " and line " << __LINE__ << std::endl;
            exit(1);
          }

          // U must be normal to the global z axis
          if (fabs(U.dot(Z)) > eps) {
            streamlog_out(ERROR) << "FPCCDSiliconTracking_MarlinTrk: PIXEL SIT Hit measurment vectors U is not in the "
                                    "global X-Y plane. \n\n exit(1) called from file "
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
        // this would be very unlikely, but who knows ... just an ordinary TrackerHit, which is not a
        // COMPOSITE_SPACEPOINT
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
        _sectors[iCode].push_back(hitExt);

        streamlog_out(DEBUG1) << " SIT Hit " << trkhit->id() << " added : @ " << pos[0] << " " << pos[1] << " "
                              << pos[2] << " drphi " << hitExt->getResolutionRPhi() << " dz "
                              << hitExt->getResolutionZ() << "  iPhi = " << iPhi << " iTheta " << iTheta
                              << " iCode = " << iCode << "  layer = " << layer << std::endl;
      }

    } catch (DataNotAvailableException& e) {
      streamlog_out(DEBUG4) << " collection not found : " << _SITHitCollection.c_str() << std::endl;
    }
  }

  return 1; // success
}

void FPCCDSiliconTracking_MarlinTrk::check(LCEvent*) {}

void FPCCDSiliconTracking_MarlinTrk::end() {
  delete _fastfitter;
  _fastfitter = 0;
  delete _encoder;
  _encoder = 0;
  // delete _trksystem ; _trksystem = 0;

  delete _moriUtil;
  delete _purityUtil;
}

void FPCCDSiliconTracking_MarlinTrk::ProcessOneSector(int iPhi, int iTheta) {
  int counter = 0;

  int nComb = int(_Combinations.size() / 3); // number of triplet combinations
  //  std::cout << iPhi << " " << iTheta << " " << _nEvt << std::endl;
  int iNC = 0;

  for (int iComb = 0; iComb < nComb; ++iComb) { // loop over triplets

    int nLR[3];

    for (int iS = 0; iS < 3; ++iS) {
      nLR[iS] = _Combinations[iNC];
      iNC++;
    }

    RangeMap::iterator rmit = _phiRangeForTriplet.find(std::make_pair(nLR[0], nLR[2]));
    if (rmit == _phiRangeForTriplet.end()) {
      std::cout << "Combinations of 3 layers for triplet may be wrong. exit." << std::endl;
      exit(1);
    }
    int iPhi_Up = iPhi + rmit->second;
    int iPhi_Low = iPhi - rmit->second;
    int iTheta_Up = iTheta + _sw_theta;
    int iTheta_Low = iTheta - _sw_theta;

    if (iTheta_Low < 0)
      iTheta_Low = 0;
    if (iTheta_Up >= _nDivisionsInTheta)
      iTheta_Up = _nDivisionsInTheta - 1;

    // index of theta-phi bin of outer most layer
    int iCode = nLR[0] + _nLayers * iPhi + _nLayers * _nDivisionsInPhi * iTheta;

    // get the all the hits in the outer most theta-phi bin

    TrackerHitExtendedVec& hitVecOuter = _sectors[iCode];

    int nHitsOuter = int(hitVecOuter.size());
    if (nHitsOuter > 0) {
      for (int ipMiddle = iPhi_Low; ipMiddle <= iPhi_Up; ipMiddle++) {       // loop over phi in the Middle
        for (int itMiddle = iTheta_Low; itMiddle <= iTheta_Up; itMiddle++) { // loop over theta in the Middle

          int iPhiMiddle = ipMiddle;

          // catch wrap-around
          if (ipMiddle < 0)
            iPhiMiddle = _nDivisionsInPhi + ipMiddle;
          if (ipMiddle >= _nDivisionsInPhi)
            iPhiMiddle = ipMiddle - _nDivisionsInPhi;

          // index of current theta-phi bin of middle layer
          iCode = nLR[1] + _nLayers * iPhiMiddle + _nLayers * _nDivisionsInPhi * itMiddle;

          // get the all the hits in the current middle theta-phi bin
          TrackerHitExtendedVec& hitVecMiddle = _sectors[iCode];

          int nHitsMiddle = int(hitVecMiddle.size());

          // determine which inner theta-phi bins to look in

          if (nHitsMiddle > 0) { // look into inner bins

            int iPhiLowInner = iPhi_Low;
            int iPhiUpInner = iPhi_Up;
            int iThetaLowInner = iTheta_Low;
            int iThetaUpInner = iTheta_Up;

            int difP = ipMiddle - iPhi;
            if (difP > 0)
              iPhiLowInner = iPhi_Low + difP;
            else if (difP < 0)
              iPhiUpInner = iPhi_Up + difP;

            int difT = itMiddle - iTheta;
            if (difT > 0)
              iThetaLowInner = iTheta_Low + difT;
            else if (difT < 0)
              iThetaUpInner = iTheta_Up + difT;

            for (int ipInner = iPhiLowInner; ipInner <= iPhiUpInner; ipInner++) {       // loop over phi in the Inner
              for (int itInner = iThetaLowInner; itInner <= iThetaUpInner; itInner++) { // loop over theta in the Inner
                int iPhiInner = ipInner;

                // catch wrap-around
                if (ipInner < 0)
                  iPhiInner = ipInner + _nDivisionsInPhi;
                if (ipInner >= _nDivisionsInPhi)
                  iPhiInner = ipInner - _nDivisionsInPhi;

                iCode = nLR[2] + _nLayers * iPhiInner + _nLayers * _nDivisionsInPhi * itInner;

                // get hit for inner bin
                TrackerHitExtendedVec& hitVecInner = _sectors[iCode];

                int nHitsInner = int(hitVecInner.size());

                if (nHitsInner > 0) {
                  streamlog_out(DEBUG1) << " " << std::setw(3) << iPhi << " " << std::setw(3) << ipMiddle << " "
                                        << std::setw(3) << ipInner << "   " << std::setw(3) << iTheta << " "
                                        << std::setw(3) << itMiddle << " " << std::setw(3) << itInner << "  "
                                        << std::setw(3) << nLR[0] << " " << std::setw(3) << nLR[1] << " "
                                        << std::setw(3) << nLR[2] << "     " << std::setw(3) << nHitsOuter << " : "
                                        << std::setw(3) << nHitsMiddle << " : " << std::setw(3) << nHitsInner
                                        << "  :: " << std::setw(3) << nHitsOuter * nHitsMiddle * nHitsInner
                                        << std::endl;
                  // test all triplets

                  for (int iOuter = 0; iOuter < nHitsOuter; ++iOuter) { // loop over hits in the outer sector
                    TrackerHitExtended* outerHit = hitVecOuter[iOuter];
                    for (int iMiddle = 0; iMiddle < nHitsMiddle; iMiddle++) { // loop over hits in the middle sector
                      TrackerHitExtended* middleHit = hitVecMiddle[iMiddle];
                      for (int iInner = 0; iInner < nHitsInner; iInner++) { // loop over hits in the inner sector
                        TrackerHitExtended* innerHit = hitVecInner[iInner];
                        HelixClass_double helix;

                        if (_mydebugstopwatch2)
                          _timer2Triplet.Start(kFALSE);
                        TrackExtended* trackAR = TestTriplet(outerHit, middleHit, innerHit, helix, 0);
                        if (_mydebugstopwatch2)
                          _timer2Triplet.Stop();

                        if (trackAR != NULL) {
                          if (_mydebugstopwatch2)
                            _timer2Build.Start(kFALSE);
                          int nHits = BuildTrack_KalFit(outerHit, middleHit, innerHit, helix, nLR[2], trackAR);
                          if (_mydebugstopwatch2)
                            _timer2Build.Stop();

                          _tracksWithNHitsContainer.getTracksWithNHitsVec(nHits).push_back(trackAR);
                          counter++;
                        }
                      } // endloop over hits in the inner sector
                    } // endloop over hits in the middle sector
                  } // endloop over hits in the outer sector
                } // endif nHitsInner > 0
              } // endloop over theta in the Inner
            } // endloop over phi in the Inner
          } // endif nHitsMiddle > 0
        } // endloop over theta in the Middle
      } // endloop over phi in the Middle
    } // endif nHitsOuter > 0
  } // endloop over triplets

  //  streamlog_out( DEBUG2 ) << " process one sectector theta,phi " << iTheta << ", " << iPhi <<
  //  "  number of loops : " << counter << std::endl  ;
}

TrackExtended* FPCCDSiliconTracking_MarlinTrk::TestTriplet(TrackerHitExtended* outerHit, TrackerHitExtended* middleHit,
                                                           TrackerHitExtended* innerHit, HelixClass_double& helix,
                                                           int omegamode = 0) {
  /*
     Methods checks if the triplet of hits satisfies helix hypothesis
   */

  // get the tracks already associated with the triplet
  TrackExtendedVec& trackOuterVec = outerHit->getTrackExtendedVec();
  TrackExtendedVec& trackMiddleVec = middleHit->getTrackExtendedVec();
  TrackExtendedVec& trackInnerVec = innerHit->getTrackExtendedVec();

  // check if all the hits are already assigned to a track
  if ((!trackOuterVec.empty()) && (!trackMiddleVec.empty()) && (!trackInnerVec.empty())) {
    TrackExtendedVec::const_iterator middleEndIter = trackMiddleVec.end();
    TrackExtendedVec::const_iterator outerEndIter = trackOuterVec.end();
    TrackExtendedVec::const_iterator innerEndIter = trackInnerVec.end();
    TrackExtendedVec::const_iterator outerBeginIter = trackOuterVec.begin();
    TrackExtendedVec::const_iterator innerBeginIter = trackInnerVec.begin();

    // loop over the tracks from the middle hit
    for (TrackExtendedVec::const_iterator middleIter = trackMiddleVec.begin(); middleIter < middleEndIter;
         ++middleIter) {
      // loop over the track from the outer hit
      for (TrackExtendedVec::const_iterator outerIter = outerBeginIter; outerIter < outerEndIter; ++outerIter) {
        // if track from the outer and middle are not the same progress
        if (*outerIter != *middleIter)
          continue;

        // loop over the tracks from the inner hit
        for (TrackExtendedVec::const_iterator innerIter = innerBeginIter; innerIter < innerEndIter; ++innerIter) {
          // no need to check against middle, it is idendical to outer here
          if (*outerIter == *innerIter) {
            // an existing track already contains all three hits
            // return a null pointer
            streamlog_out(DEBUG2) << " TestTriplet: track " << *outerIter
                                  << " already contains all three hits: Do not create new track from these hits "
                                  << std::endl;
          }

        } // for inner
      } // for outer
    } // for middle
  } // if all vectors are not empty

  // increase triplet count
  ++_ntriplets;

  // get the hit coordinates and errors
  double xh[3];
  double yh[3];
  float zh[3];
  double wrh[3];
  float wzh[3];
  float rh[3];
  float ph[3];

  float par[5];
  float epar[15];

  // first hit
  xh[0] = outerHit->getTrackerHit()->getPosition()[0];
  yh[0] = outerHit->getTrackerHit()->getPosition()[1];
  zh[0] = float(outerHit->getTrackerHit()->getPosition()[2]);
  wrh[0] = double(1.0 / (outerHit->getResolutionRPhi() * outerHit->getResolutionRPhi()));
  wzh[0] = 1.0 / (outerHit->getResolutionZ() * outerHit->getResolutionZ());

  // second hit
  xh[1] = middleHit->getTrackerHit()->getPosition()[0];
  yh[1] = middleHit->getTrackerHit()->getPosition()[1];
  zh[1] = float(middleHit->getTrackerHit()->getPosition()[2]);
  wrh[1] = double(1.0 / (middleHit->getResolutionRPhi() * middleHit->getResolutionRPhi()));
  wzh[1] = 1.0 / (middleHit->getResolutionZ() * middleHit->getResolutionZ());

  // third hit
  xh[2] = innerHit->getTrackerHit()->getPosition()[0];
  yh[2] = innerHit->getTrackerHit()->getPosition()[1];
  zh[2] = float(innerHit->getTrackerHit()->getPosition()[2]);
  wrh[2] = double(1.0 / (innerHit->getResolutionRPhi() * innerHit->getResolutionRPhi()));
  wzh[2] = 1.0 / (innerHit->getResolutionZ() * innerHit->getResolutionZ());

  // calculate r and phi for all hits
  for (int ih = 0; ih < 3; ih++) {
    rh[ih] = float(sqrt(xh[ih] * xh[ih] + yh[ih] * yh[ih]));
    ph[ih] = atan2(yh[ih], xh[ih]);
    if (ph[ih] < 0.)
      ph[ih] = TWOPI + ph[ih];
  }

  int NPT = 3;
  int iopt = 2;
  float chi2RPhi;
  float chi2Z;

  streamlog_out(DEBUG2) << " TestTriplet: Use fastHelixFit " << std::endl;

  _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh, iopt, par, epar, chi2RPhi, chi2Z);
  par[3] = par[3] * par[0] / fabs(par[0]);

  // get helix parameters
  float omega = par[0];
  float tanlambda = par[1];
  float phi0 = par[2];
  float d0 = par[3];
  float z0 = par[4];

  // chi2 is weighted here by a factor for both rphi and z
  float Chi2 = chi2RPhi * _chi2WRPhiTriplet + chi2Z * _chi2WZTriplet;
  int ndf = 2 * NPT - 5;

  // check the truth information for the triplet

  // define these outside of the ifdef so that we don't need to keep repeating it.
  std::vector<TrackerHit*> hit_list;
  std::vector<MCParticle*> mcps_imo;
  std::vector<MCParticle*> mcp_s;
  //  int nmcps   = 0;
  //  int nbadHits = 0;
  int triplet_code = 0;

  // Check if track satisfies all conditions

  bool failed = false;

  int quality_code = triplet_code * 10;

  if (std::isnormal(d0) == false) {
    if (std::isinf(d0) == true) {
      streamlog_out(DEBUG1) << "d0 is inf" << std::endl;
    } else if (std::isnan(d0) == true) {
      streamlog_out(DEBUG1) << "d0 is nan" << std::endl;
    } else {
      streamlog_out(DEBUG1) << "Something is wrong with d0" << std::endl;
    }
    failed = true;
    quality_code += 10;
  } else if (std::isnormal(z0) == false) {
    if (std::isinf(z0)) {
      streamlog_out(DEBUG1) << "z0 is inf" << std::endl;
    } else if (std::isnan(z0)) {
      streamlog_out(DEBUG1) << "z0 is nan" << std::endl;
    } else {
      streamlog_out(DEBUG1) << "Something is wrong with z0" << std::endl;
    }
    failed = true;
    quality_code += 100;
  } else if (std::isnormal(omega) == false) {
    if (std::isinf(omega)) {
      streamlog_out(DEBUG1) << "omega is inf" << std::endl;
    } else if (std::isnan(omega)) {
      streamlog_out(DEBUG1) << "omega is nan" << std::endl;
    } else {
      streamlog_out(DEBUG1) << "Something is wrong with omega" << std::endl;
    }
    failed = true;
    quality_code += 1000;
  } else if (std::isnormal(phi0) == false) {
    if (std::isinf(phi0)) {
      streamlog_out(DEBUG1) << "phi0 is inf" << std::endl;
    } else if (std::isnan(phi0)) {
      streamlog_out(DEBUG1) << "phi0 is nan" << std::endl;
    } else {
      streamlog_out(DEBUG1) << "Something is wrong with phi0" << std::endl;
    }
    failed = true;
    quality_code += 10000;
  } else if (std::isnormal(tanlambda) == false) {
    if (std::isinf(tanlambda)) {
      streamlog_out(DEBUG1) << "tanlambda is inf" << std::endl;
    } else if (std::isnan(tanlambda)) {
      streamlog_out(DEBUG1) << "tanlambda is nan" << std::endl;
    } else {
      streamlog_out(DEBUG1) << "Something is wrong with tanlambda" << std::endl;
    }
    failed = true;
    quality_code += 100000;
  }

  if (failed == false) {
    if (Chi2 / float(ndf) > _chi2FitCut) {
      streamlog_out(DEBUG1) << "Chi2/ndf = " << Chi2 / float(ndf) << " , cut = " << _chi2FitCut << std::endl;
      failed = true;
      quality_code += 1;
    } else if (Chi2 / float(ndf) < 0) {
      streamlog_out(DEBUG1) << "Chi2/ndf = " << Chi2 / float(ndf) << " , cut = " << _chi2FitCut << std::endl;
      failed = true;
      quality_code += 5;
    } else if (fabs(d0) > _cutOnD0) {
      streamlog_out(DEBUG1) << "d0 = " << d0 << " , cut = " << _cutOnD0 << std::endl;
      failed = true;
      quality_code += 2;
    } else if (fabs(z0) > _cutOnZ0) {
      streamlog_out(DEBUG1) << "z0 = " << z0 << " , cut = " << _cutOnZ0 << std::endl;
      failed = true;
      quality_code += 3;
    } else if (omegamode == 0 && fabs(omega) > _cutOnOmegaVXD) {
      streamlog_out(DEBUG1) << "omega = " << omega << " , cut = " << _cutOnOmegaVXD << std::endl;
      failed = true;
      quality_code += 4;
    } else if (omegamode == 1 && fabs(omega) > _cutOnOmegaFTD) {
      streamlog_out(DEBUG1) << "omega = " << omega << " , cut = " << _cutOnOmegaFTD << std::endl;
      failed = true;
      quality_code += 4;
    } else {
      streamlog_out(DEBUG1) << "Success !!!!!!!" << std::endl;
    }
  }

  if (failed) {
    // return a null pointer
    return 0;
  }

  helix.Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _bField);

  TrackExtended* trackAR = new TrackExtended();
  trackAR->addTrackerHitExtended(outerHit);
  trackAR->addTrackerHitExtended(middleHit);
  trackAR->addTrackerHitExtended(innerHit);
  outerHit->addTrackExtended(trackAR);
  middleHit->addTrackExtended(trackAR);
  innerHit->addTrackExtended(trackAR);
  trackAR->setD0(d0);
  trackAR->setZ0(z0);
  trackAR->setPhi(phi0);
  trackAR->setTanLambda(tanlambda);
  trackAR->setOmega(omega);
  trackAR->setChi2(Chi2);
  trackAR->setNDF(ndf);
  trackAR->setCovMatrix(epar);

  return trackAR;
}

int FPCCDSiliconTracking_MarlinTrk::BuildTrack_KalFit(TrackerHitExtended* /*outerHit*/,
                                                      TrackerHitExtended* /*middleHit*/,
                                                      TrackerHitExtended* /*innerHit*/, HelixClass_double& helix,
                                                      int innerLayer, TrackExtended* trackAR) {
  int nMisAssign = 0;
  for (int layer = innerLayer - 1; layer >= 0; layer--) { // loop over remaining layers
    if (nMisAssign > _minMissAddition) {
      return int(trackAR->getTrackerHitExtendedVec().size());
    }

    int Boundaries[4];
    int pterr = getPhiThetaRegion(trackAR, layer, Boundaries);
    // When the number different from 0 comes, then it is the number of hits as the
    // last result of buildTrack
    if (pterr != 0) {
      return pterr;
    }

    float distMin = 1.0e+20;
    TrackerHitExtended* assignedhit = NULL;

    if (Boundaries[2] < 0)
      Boundaries[2] = 0;
    if (Boundaries[3] >= _nDivisionsInTheta)
      Boundaries[3] = _nDivisionsInTheta - 1;

    for (int ipInner = Boundaries[0]; ipInner <= Boundaries[1]; ipInner++) {
      for (int itInner = Boundaries[2]; itInner <= Boundaries[3]; itInner++) {
        int iPhiInner = ipInner;
        while (iPhiInner < 0)
          iPhiInner -= _nDivisionsInPhi;
        while (iPhiInner >= _nDivisionsInPhi)
          iPhiInner -= _nDivisionsInPhi;
        int iCode = layer + _nLayers * iPhiInner + _nLayers * _nDivisionsInPhi * itInner;
        TrackerHitExtendedVec& hitVecInner = _sectors[iCode];
        int nHitsInner = int(hitVecInner.size());
        for (int iInner = 0; iInner < nHitsInner; iInner++) {
          TrackerHitExtended* currentHit = hitVecInner[iInner];
          double pos[3];
          double distance[3];
          for (int i = 0; i < 3; ++i) {
            pos[i] = double(currentHit->getTrackerHit()->getPosition()[i]);
          }

          int tiltStatus = 0;
          double dot = 0;

          if (_useClusterRejection == true) {
            TrackerHit* curInMostHit = trackAR->getTrackerHitExtendedVec().back()->getTrackerHit();
            tiltStatus = CheckTiltOf2Clusters(curInMostHit, currentHit->getTrackerHit(), 1);
            // In BuildTrack, in this case the first and second arguments are
            // always VXD, so needless to check which it is VXD or SIT hit.
            dot = DotOf2Clusters(curInMostHit, currentHit->getTrackerHit());
          }
          bool goodDot = (dot > _minDotOf2Clusters) ? true : false;
          if (_useClusterRejection == false || (tiltStatus >= 0 && goodDot == true)) {
            double time = helix.getDistanceToPoint(pos, distance);
            if (time < 1.0e+10) {
              if (distance[2] < distMin) { // distance[2] = sqrt( d0*d0 + z0*z0 )
                distMin = distance[2];
                assignedhit = currentHit;
              }
            }
          } // end of if statement of tiltStatus //MORI
        } // endloop over hits in the Inner sector
      } // endloop over theta in the Inner region
    } // endloop over phi in the Inner region

    if (distMin > _minDistCutAttachForVXD) {
      nMisAssign++;
      continue;
    } else {
      float par[5];
      float epar[15];
      int ndf;
      float Chi2;
      TrackerHitExtendedVec& hvec = trackAR->getTrackerHitExtendedVec();
      bool validCombination = 0;

      TrackerHitVec trkHits;
      trkHits.reserve(hvec.size() + 1);
      trkHits.push_back(assignedhit->getTrackerHit());
      for (int ij = int(hvec.size()) - 1; ij >= 0; ij--) {
        trkHits.push_back(hvec[ij]->getTrackerHit());
      }
      TrackerHitVec hits_in_fit;
      TrackerHitVec outliers;
      int error_KalFit = KalFit(ndf, Chi2, trkHits, hits_in_fit, outliers, par, epar, helix);

      if (error_KalFit != 0) {
        validCombination = 0;
        streamlog_out(DEBUG2) << "ERROR! error_KalFit = " << error_KalFit << std::endl;
      } else if (outliers.size() == 0) {
        validCombination = Chi2 / float(ndf) < _chi2FitCut_kalman;
      } else {
        validCombination = Chi2 / float(ndf) < _chi2FitCut_kalman;
        nMisAssign += int(outliers.size());
        if (int(outliers.size()) > _minMissAddition) {
          return int(trackAR->getTrackerHitExtendedVec().size());
        }
        TrackerHitExtendedVec oldhvec = hvec;
        trackAR->ClearTrackerHitExtendedVec();
        for (int hi = 0; hi < int(oldhvec.size()); hi++) {
          TrackerHitVec::iterator iter = std::find(outliers.begin(), outliers.end(), oldhvec[hi]->getTrackerHit());
          if (iter == outliers.end()) {
            trackAR->addTrackerHitExtended(oldhvec[hi]);
          } else {
            outliers.erase(iter);
          } //"found" = "oldhvec[hi] is outlier"
        }
        if (outliers.size() != 0) {
          TrackerHitVec::iterator iter = std::find(outliers.begin(), outliers.end(), assignedhit->getTrackerHit());
          if (iter != outliers.end()) {
            outliers.erase(iter);
            validCombination = 0;
            nMisAssign -= 1; // Unless this is here, nMisAssign will be double-counted.
          }
        }
      }

      if (validCombination == 0)
        nMisAssign++;
      else {
        // assign hit to track and track to hit, update the track parameters
        trackAR->addTrackerHitExtended(assignedhit);
        assignedhit->addTrackExtended(trackAR);
        // Mori：two conventions of fastFit and KalFit are different.
        // I made KalFit's convention match that of lcio's covMatrix.
        float parD0 = par[0];
        float parPhi0 = par[1];
        float parOmega = par[2];
        float parZ0 = par[3];
        float parTanlambda = par[4];
        helix.Initialize_Canonical(parPhi0, parD0, parZ0, parOmega, parTanlambda, _bField);
        trackAR->setD0(parD0);
        trackAR->setZ0(parZ0);
        trackAR->setPhi(parPhi0);
        trackAR->setTanLambda(parTanlambda);
        trackAR->setOmega(parOmega);
        trackAR->setChi2(Chi2);
        trackAR->setCovMatrix(epar);
        trackAR->setNDF(ndf);
      }
    }
  } // endloop over remaining layers

  TrackerHitExtendedVec& hvec = trackAR->getTrackerHitExtendedVec();
  int nTotalHits = int(hvec.size());

  return nTotalHits;
}

int FPCCDSiliconTracking_MarlinTrk::getIntersectionEasy(HelixClass_double& helix, TrackerHit* curInmos, int layer,
                                                        double* isec, double* ref) {
  // int FPCCDSiliconTracking_MarlinTrk::getIntersectionEasy(HelixClass_double& helix, TrackerHit* curInmos , int layer,
  // double* isec, double* ref){

  if (layer > 5) {
    std::cout << "getIntersectionEasy used non-VXD layer by mistake. Check source code." << std::endl;
    exit(1);
  }

  isec[0] = 0;
  isec[1] = 0;
  isec[2] = 0;
  double hlwidth = _vxd.geodata[layer].sximax;
  double Rmin = _vxd.geodata[layer].rmes;
  double Rmax = sqrt(Rmin * Rmin + hlwidth * hlwidth);
  double Rmean = (Rmin + Rmax) / 2.0;
  double point[6];
  double time = helix.getPointOnCircle(Rmean, ref, point);
  if (time < 0) {
    return -1;
  }
  double phiX = atan2(curInmos->getPosition()[1], curInmos->getPosition()[0]);
  double diffA = atan2(point[1], point[0]) - phiX;
  double diffB = atan2(point[4], point[3]) - phiX;

  if (std::abs(diffA) < std::abs(diffB)) {
    isec[0] = point[0];
    isec[1] = point[1];
    isec[2] = point[2];
  } else {
    isec[0] = point[3];
    isec[1] = point[4];
    isec[2] = point[5];
  }

  return 0;
}

int FPCCDSiliconTracking_MarlinTrk::getIntersectionEasyTest(HelixClass_double& helix, TrackerHit* basis, int layer,
                                                            std::vector<double>& isec) {
  if (layer > 8) {
    std::cout << "getIntersectionEasyTest uses only VXD or SIT layer. Check source code." << std::endl;
    exit(1);
  }
  long double d0 = helix.getD0();
  long double z0 = helix.getZ0();
  long double phi0 = helix.getPhi0();
  long double tanL = helix.getTanLambda();
  long double omega = helix.getOmega();

  if (std::isnormal(d0) == false || std::isnormal(z0) == false || std::isnormal(omega) == false ||
      std::isnormal(phi0) == false || std::isnormal(tanL) == false) {
    return -3;
  }

  isec.clear();
  isec.resize(3);

  long double Rmean;
  if (layer >= int(_nLayersVTX)) {
    Rmean = (_sit.geodata[layer - 6].rmin + _sit.geodata[layer - 5].rmin) / 2.0;
  } else {
    long double hlwidth = _vxd.geodata[layer].sximax;
    long double Rmin = _vxd.geodata[layer].rmes;
    long double Rmax = sqrt(Rmin * Rmin + hlwidth * hlwidth);
    Rmean = (Rmin + Rmax) / 2.0;
  }

  long double A = 1.0 - 0.5 * (Rmean * Rmean - d0 * d0) / (1. / omega / omega - d0 / omega);
  if (_mydebugIntersection) {
    std::cout << " == getIntersectionEasyTest check ===" << std::endl;
    std::cout << "A : " << A << std::endl;
  }
  if (-1.0 > A || A > 1.0) {
    return -1;
  } else if (std::isnormal(A) == false) {
    return -2;
  }

  long double phi1 = phi0 - acos(A);
  long double phi2 = phi0 + acos(A);

  long double* point = new long double[6];
  point[0] = -d0 * sin(phi0) + 1.0 / omega * (sin(phi0) - sin(phi1));
  point[1] = d0 * cos(phi0) - 1.0 / omega * (cos(phi0) - cos(phi1));
  point[2] = z0 + 1.0 / omega * (phi0 - phi1) * tanL;
  point[3] = -d0 * sin(phi0) + 1.0 / omega * (sin(phi0) - sin(phi2));
  point[4] = d0 * cos(phi0) - 1.0 / omega * (cos(phi0) - cos(phi2));
  point[5] = z0 + 1.0 / omega * (phi0 - phi2) * tanL;
  long double phiX = atan2(basis->getPosition()[1], basis->getPosition()[0]);
  if (phiX < 0)
    phiX += 2.0 * M_PI;
  long double p1 = atan2(point[1], point[0]);
  if (p1 < 0)
    p1 += 2.0 * M_PI;
  long double p2 = atan2(point[4], point[3]);
  if (p2 < 0)
    p2 += 2.0 * M_PI;
  long double diff1 = std::abs(p1 - phiX);
  if (diff1 > M_PI)
    diff1 = 2.0 * M_PI - diff1;
  long double diff2 = std::abs(p2 - phiX);
  if (diff2 > M_PI)
    diff2 = 2.0 * M_PI - diff2;
  isec[0] = (diff1 < diff2) ? point[0] : point[3];
  isec[1] = (diff1 < diff2) ? point[1] : point[4];
  long double tmpZ = (diff1 < diff2) ? point[2] : point[5];
  long double curZ = basis->getPosition()[2];
  long double dist = curZ - tmpZ;
  long double distOne = std::abs(1.0 / omega * tanL * 2.0 * M_PI);
  long double nCircle = dist / distOne;
  tmpZ += distOne * static_cast<long double>(static_cast<int>(nCircle));
  long double decimal = std::abs(nCircle - int(nCircle));
  if (decimal > 0.5)
    tmpZ += distOne * nCircle / std::abs(nCircle);
  isec[2] = tmpZ;

#if 0
  std::cout << "  point[0] : " << point[0] << std::endl;
  std::cout << "  point[1] : " << point[1] << std::endl;
  std::cout << "  point[2] : " << point[2] << std::endl;
  std::cout << "  point[3] : " << point[3] << std::endl;
  std::cout << "  point[4] : " << point[4] << std::endl;
  std::cout << "  point[5] : " << point[5] << std::endl;
  std::cout << "  point[6] : " << point[6] << std::endl;
  std::cout << "  point[7] : " << point[7] << std::endl;

  std::cout << "  phi1     : " << phi1    << std::endl;
  std::cout << "  phi2     : " << phi2    << std::endl;
  std::cout << "  phiX     : " << phiX    << std::endl;
  std::cout << "  diff1     : " << diff1    << std::endl;
  std::cout << "  diff2     : " << diff2    << std::endl;
  std::cout << "  isec [0] : " << isec[0] << std::endl;
  std::cout << "  isec [1] : " << isec[1] << std::endl;
  std::cout << "  distOne     : " << distOne    << std::endl;
  std::cout << "  nCircle     : " << nCircle    << std::endl;
  std::cout << "  decimal     : " << decimal    << std::endl;
  std::cout << "  isec [2] : " << isec[2] << std::endl;
#endif

  delete[] point;

  return 0;
}

double FPCCDSiliconTracking_MarlinTrk::getNeededPhiSectors(double Pt, int outly, int inly) {
  double R = Pt / (0.299792458 * _bField);
  R = 1000.0 * R;
  double omega = 1.0 / R;
  HelixClass_double helix;
  // helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
  helix.Initialize_Canonical(0, 0, 0, omega, 0.001, _bField);
  double Ro;
  double Ri;
  if (outly >= 6) {
    Ro = (_sit.geodata[outly - 6].rmin + _sit.geodata[outly - 5].rmin) / 2.0;
  } else {
    double hlwidth = _vxd.geodata[outly].sximax;
    double Rmin = _vxd.geodata[outly].rmes;
    double Rmax = sqrt(Rmin * Rmin + hlwidth * hlwidth);
    Ro = (Rmin + Rmax) / 2.0;
  }
  // For now, I assume that inly always corresponds to VXD layer due to the
  // current geometry.
  double hlwidth = _vxd.geodata[inly].sximax;
  double Rmin = _vxd.geodata[inly].rmes;
  double Rmax = sqrt(Rmin * Rmin + hlwidth * hlwidth);
  Ri = (Rmin + Rmax) / 2.0;

  double ref[3] = {0, 0, 0};
  double point[6];
  double time = helix.getPointOnCircle(Ro, ref, point);
  if (time < 0) {
    std::cout << "_cutOnPtVXD is too small to calculate range of phi-sectors needed to find triplet with Pt more than "
                 "_cutOnPtVXD."
              << std::endl;
    std::cout << "Please set larger value." << std::endl;
    exit(1);
  }

  double phiO = atan2(point[1], point[0]);
  std::cout << "==Ro==" << std::endl;
  std::cout << "  time : " << time << std::endl;
  std::cout << "  phi : " << phiO << std::endl;
  printf("point : %f,%f,%f  %f,%f,%f \n ", point[0], point[1], point[2], point[3], point[4], point[5]);

  time = helix.getPointOnCircle(Ri, ref, point);
  if (time < 0) {
    std::cout << "_cutOnPtVXD is too small to calculate range of phi-sectors needed to find triplet with Pt more than "
                 "_cutOnPtVXD."
              << std::endl;
    std::cout << "Please set larger value." << std::endl;
    exit(1);
  }

  double phiI = atan2(point[1], point[0]);
  std::cout << "==Ri==" << std::endl;
  std::cout << "  time : " << time << std::endl;
  std::cout << "  phi : " << phiI << std::endl;
  printf("point : %f,%f,%f  %f,%f,%f \n ", point[0], point[1], point[2], point[3], point[4], point[5]);
  printf("phiRange : %f \n ", phiI - phiO);
  return phiI - phiO;
}

// This is under construction
void FPCCDSiliconTracking_MarlinTrk::getNeededPhiSectorsVer2(double Pt, std::vector<int> layers,
                                                             std::vector<double>& phiRange) {
  double R = Pt / (0.299792458 * _bField);
  R = 1000.0 * R;
  double omega = 1.0 / R;
  HelixClass_double helix;
  // helix.Initialize_Canonical(phi0,d0,z0,omega,tanlambda,_bField);
  helix.Initialize_Canonical(0, 0, 0, omega, 0.001, _bField);
  double Ro, Rm, Ri;
  if (layers[0] >= 6) { // outer layer
    Ro = (_sit.geodata[layers[0] - 6].rmin + _sit.geodata[layers[0] - 5].rmin) / 2.0;
  } else {
    double hlwidth = _vxd.geodata[layers[0]].sximax;
    double Rmin = _vxd.geodata[layers[0]].rmes;
    double Rmax = sqrt(Rmin * Rmin + hlwidth * hlwidth);
    Ro = (Rmin + Rmax) / 2.0;
  }
  if (layers[1] >= 6) { // middle layer
    Rm = (_sit.geodata[layers[1] - 6].rmin + _sit.geodata[layers[1] - 5].rmin) / 2.0;
  } else {
    double hlwidth = _vxd.geodata[layers[1]].sximax;
    double Rmin = _vxd.geodata[layers[1]].rmes;
    double Rmax = sqrt(Rmin * Rmin + hlwidth * hlwidth);
    Rm = (Rmin + Rmax) / 2.0;
  }

  double hlwidth = _vxd.geodata[layers[2]].sximax;
  double Rmin = _vxd.geodata[layers[2]].rmes;
  double Rmax = sqrt(Rmin * Rmin + hlwidth * hlwidth);
  Ri = (Rmin + Rmax) / 2.0;

  double ref[3] = {0, 0, 0};
  double point[6];
  // outer phi
  double time = helix.getPointOnCircle(Ro, ref, point);
  if (time < 0) {
    std::cout << "_cutOnPtVXD is too small to calculate range of phi-sectors needed to find triplet with Pt more than "
                 "_cutOnPtVXD."
              << std::endl;
    exit(1);
  }
  double phiO = atan2(point[1], point[0]);
  std::cout << "Ro phi :" << phiO << std::endl;

  // middle phi
  time = helix.getPointOnCircle(Rm, ref, point);
  if (time < 0) {
    std::cout << "_cutOnPtVXD is too small to calculate range of phi-sectors needed to find triplet with Pt more than "
                 "_cutOnPtVXD."
              << std::endl;
    exit(1);
  }
  double phiM = atan2(point[1], point[0]);
  std::cout << "Rm phi :" << phiM << std::endl;

  // inner phi
  time = helix.getPointOnCircle(Ri, ref, point);
  if (time < 0) {
    std::cout << "_cutOnPtVXD is too small to calculate range of phi-sectors needed to find triplet with Pt more than "
                 "_cutOnPtVXD."
              << std::endl;
    exit(1);
  }
  double phiI = atan2(point[1], point[0]);
  std::cout << "Ri phi :" << phiI << std::endl;

  // Check
  printf("phiRange[0] (for middle layer) : %f \n ", phiM - phiO);
  printf("phiRange[1] (for inner layer) : %f \n ", phiI - phiO);

  // storing
  phiRange.clear();
  phiRange.push_back(std::abs(phiM - phiO));
  phiRange.push_back(std::abs(phiI - phiO));

  return;
}

void FPCCDSiliconTracking_MarlinTrk::Sorting(TrackExtendedVec& trackVec) {
  /**
    Sorting of Track Vector in ascending order of chi2/ndf
   */

  std::sort(trackVec.begin(), trackVec.end(), compare_TrackExtended());

  // also clean up? what does this do here?
  for (size_t i = 0, sizeOfVector = trackVec.size(); i < sizeOfVector; ++i) {
    TrackerHitExtendedVec& hitVec = trackVec[i]->getTrackerHitExtendedVec();
    int nHits = int(hitVec.size());

    for (int ih = 0; ih < nHits; ih++) {
      hitVec[ih]->clearTrackVec();
    }
  }
}

void FPCCDSiliconTracking_MarlinTrk::CreateTrack(TrackExtended* trackAR) {
  /**
    Method which creates Track out of TrackExtended objects. Checks for possible
    track splitting (separate track segments in VXD and FTD).
   */

  TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();
  int nHits = int(hitVec.size());

  for (int i = 0; i < nHits; ++i) {
    TrackExtendedVec& trackVec = hitVec[i]->getTrackExtendedVec();
    if (trackVec.size() != 0)
      return;
  }

  // First check if the current track is piece of the split one
  // look for matching track segment

  int found = 0;

  int nTrk = int(_trackImplVec.size());

  for (int itrk = 0; itrk < nTrk; ++itrk) {
    TrackExtended* trackOld = _trackImplVec[itrk];
    TrackerHitExtendedVec& hitVecOld = trackOld->getTrackerHitExtendedVec();

    float phiNew = trackAR->getPhi();
    float phiOld = trackOld->getPhi();
    float thetaNew = M_PI_2 - atan(trackAR->getTanLambda());
    float thetaOld = M_PI_2 - atan(trackOld->getTanLambda());

    float angle = (cos(phiNew) * cos(phiOld) + sin(phiNew) * sin(phiOld)) * sin(thetaNew) * sin(thetaOld) +
                  cos(thetaNew) * cos(thetaOld);
    angle = acos(angle);

    if (angle < _angleCutForMerging) {
      int nHitsOld = int(hitVecOld.size());
      int nTotHits = nHits + nHitsOld;
      double* xh = new double[nTotHits];
      double* yh = new double[nTotHits];
      float* zh = new float[nTotHits];
      double* wrh = new double[nTotHits];
      float* wzh = new float[nTotHits];
      float* rh = new float[nTotHits];
      float* ph = new float[nTotHits];
      float par[5];
      float epar[15];
      // float refPoint[3] = {0.,0.,0.};
      for (int ih = 0; ih < nHits; ++ih) {
        TrackerHit* trkHit = hitVec[ih]->getTrackerHit();
        float rR = hitVec[ih]->getResolutionRPhi();
        float rZ = hitVec[ih]->getResolutionZ();
        if (int(hitVec[ih]->getTrackExtendedVec().size()) != 0)
          streamlog_out(DEBUG2) << "WARNING : HIT POINTS TO TRACK " << std::endl;
        xh[ih] = trkHit->getPosition()[0];
        yh[ih] = trkHit->getPosition()[1];
        zh[ih] = float(trkHit->getPosition()[2]);
        wrh[ih] = double(1.0 / (rR * rR));
        wzh[ih] = 1.0 / (rZ * rZ);
        rh[ih] = float(sqrt(xh[ih] * xh[ih] + yh[ih] * yh[ih]));
        ph[ih] = float(atan2(yh[ih], xh[ih]));
      }
      for (int ih = 0; ih < nHitsOld; ++ih) {
        TrackerHit* trkHit = hitVecOld[ih]->getTrackerHit();
        xh[ih + nHits] = trkHit->getPosition()[0];
        yh[ih + nHits] = trkHit->getPosition()[1];
        zh[ih + nHits] = float(trkHit->getPosition()[2]);
        float rR = hitVecOld[ih]->getResolutionRPhi();
        float rZ = hitVecOld[ih]->getResolutionZ();
        wrh[ih + nHits] = double(1.0 / (rR * rR));
        wzh[ih + nHits] = 1.0 / (rZ * rZ);
        rh[ih + nHits] = float(sqrt(xh[ih + nHits] * xh[ih + nHits] + yh[ih + nHits] * yh[ih + nHits]));
        ph[ih + nHits] = float(atan2(yh[ih + nHits], xh[ih + nHits]));
      }
      int NPT = nTotHits;
      int iopt = 2;
      float chi2RPhi;
      float chi2Z;
      int ndf = 2 * NPT - 5;

      _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh, iopt, par, epar, chi2RPhi, chi2Z);
      par[3] = par[3] * par[0] / fabs(par[0]);

      float omega = par[0];
      float tanlambda = par[1];
      float phi0 = par[2];
      float d0 = par[3];
      float z0 = par[4];

      float eparmin[15];
      for (int iparam = 0; iparam < 15; ++iparam)
        eparmin[iparam] = epar[iparam];

      // float refPointMin[3];
      // for (int ipp=0;ipp<3;++ipp)
      //   refPointMin[ipp] = refPoint[ipp];

      float chi2Min = chi2RPhi * _chi2WRPhiSeptet + chi2Z * _chi2WZSeptet;
      chi2Min = chi2Min / float(ndf);

      // float chi2MinRPhi = chi2RPhi;
      // float chi2MinZ = chi2Z;

      int iBad = -1;
      if (chi2Min < _chi2FitCut && std::isnormal(chi2Min) == true && chi2Min > 0) {
        // Mori added new requirement that chi2Min be normal value on September 1, 2013.
        found = 1;
      } else { // SJA:FIXME: UH What is going on here? setting weights to 0 and refitting?
        float* wzhOld = new float[nTotHits];
        double* wrhOld = new double[nTotHits];
        for (int i = 0; i < nTotHits; ++i) {
          wzhOld[i] = wzh[i];
          wrhOld[i] = wrh[i];
        }
        for (int i = 0; i < nTotHits; ++i) {
          for (int j = 0; j < nTotHits; ++j) {
            if (i == j) {
              wrh[j] = 0.0;
              wzh[j] = 0.0;
            } else {
              wrh[j] = wrhOld[j];
              wzh[j] = wzhOld[j];
            }
          }

          _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh, iopt, par, epar, chi2RPhi, chi2Z);
          par[3] = par[3] * par[0] / fabs(par[0]);

          float chi2Cur = chi2RPhi * _chi2WRPhiSeptet + chi2Z * _chi2WZSeptet;
          chi2Cur = chi2Cur / float(ndf);

          if (chi2Cur < chi2Min && std::isnormal(chi2Min) == 1 && chi2Min > 0) {
            // Mori added new requirement that chi2Min be normal value on September 1, 2013.
            chi2Min = chi2Cur;
            // chi2MinRPhi = chi2RPhi;
            // chi2MinZ = chi2Z;
            omega = par[0];
            tanlambda = par[1];
            phi0 = par[2];
            d0 = par[3];
            z0 = par[4];
            for (int iparam = 0; iparam < 15; ++iparam)
              eparmin[iparam] = epar[iparam];
            // for (int ipp=0;ipp<3;++ipp)
            //   refPointMin[ipp] = refPoint[ipp];
            iBad = i;
          }
        }

        if (chi2Min < _chi2FitCut && std::isnormal(chi2Min) == 1 && chi2Min > 0) {
          // Mori added new requirement that chi2Min be normal value on September 1, 2013.
          found = 1;
        }
        delete[] wzhOld;
        delete[] wrhOld;
      }

      // Split track is found.
      // Attach hits belonging to the current track segment to
      // the track already created
      if (found == 1) {
        trackOld->ClearTrackerHitExtendedVec();
        for (int i = 0; i < nHits; ++i) {
          TrackerHitExtended* trkHit = hitVec[i];
          trkHit->clearTrackVec();
          if (i == iBad) {
          } else {
            trackOld->addTrackerHitExtended(trkHit);
            trkHit->addTrackExtended(trackOld);
          }
        }
        for (int i = 0; i < nHitsOld; ++i) {
          int icur = i + nHits;
          TrackerHitExtended* trkHit = hitVecOld[i];
          trkHit->clearTrackVec();
          if (icur == iBad) {
          } else {
            trackOld->addTrackerHitExtended(trkHit);
            trkHit->addTrackExtended(trackOld);
          }
        }
        trackOld->setOmega(omega);
        trackOld->setTanLambda(tanlambda);
        trackOld->setPhi(phi0);
        trackOld->setD0(d0);
        trackOld->setZ0(z0);

        //      std::cout << "Split track found " << d0 << " " << z0 << std::endl;

        // killeb:  In the original SiliconTracking this was in the NOT simple helix branch.
        // The rest of the code uses the simple helix branch, where ndf_D is never set.
        // In fact it has never been initialised or used anywhere. I think this line should not be executed.
        // ndf = ndf_D;

        trackOld->setChi2(chi2Min * float(ndf));
        trackOld->setNDF(ndf);
        trackOld->setCovMatrix(eparmin);
        //      trackOld->setReferencePoint(refPointMin);
      }

      delete[] xh;
      delete[] yh;
      delete[] zh;
      delete[] wrh;
      delete[] wzh;
      delete[] rh;
      delete[] ph;
    }
    if (found == 1)
      break;
  }

  // Candidate is a unique track
  // No other segments are found
  if (found == 0) {
    _trackImplVec.push_back(trackAR);
    for (int i = 0; i < nHits; ++i) {
      TrackerHitExtended* hit = hitVec[i];
      hit->addTrackExtended(trackAR);
    }
  }

  return;
}

void FPCCDSiliconTracking_MarlinTrk::AttachRemainingVTXHitsVeryFast() {
  int nTracks = int(_trackImplVec.size());
  for (int iTrk = 0; iTrk < nTracks; ++iTrk) {
    TrackExtended* trackAR = _trackImplVec[iTrk];
    for (int nthLayer = 0; nthLayer < _nLayers; nthLayer++) {
      if (nthLayer == 7 || nthLayer == 9)
        continue; //_nLayersSIT is always 4. 7 & 9 should be avoided.
      if (nthLayer < _minimalLayerToAttach)
        continue;
      int Boundaries[4];
      int pterr = getPhiThetaRegion(trackAR, nthLayer, Boundaries);
      if (pterr != 0)
        continue;
      TrackerHitExtended* trkhitToAttach = NULL;
      float minDist = 1.0e+6;
      for (int ip = Boundaries[0]; ip <= Boundaries[1]; ++ip) {
        for (int it = Boundaries[2]; it <= Boundaries[3]; ++it) {
          int iCode = nthLayer + _nLayers * ip + _nLayers * _nDivisionsInPhi * it;
          TrackerHitExtendedVec& hitVec = _sectors[iCode];
          for (int iH = 0; iH < int(hitVec.size()); ++iH) {
            TrackerHitExtended* hit = hitVec[iH];
            unsigned int maxTrackSize = 0;
            if (_keepCandidate == false) {
              TrackExtendedVec& trackVec = hit->getTrackExtendedVec();
              for (unsigned int itrack = 0; itrack < trackVec.size(); itrack++) {
                TrackerHitExtendedVec hitVec_tmp = trackVec[itrack]->getTrackerHitExtendedVec();
                if (hitVec_tmp.size() > maxTrackSize)
                  maxTrackSize = hitVec_tmp.size();
                if (maxTrackSize > 3)
                  break;
              }
            }
            if (maxTrackSize <= 3) {
              bool consider = true;
              if (_checkForDelta) { // I use this without any modification.
                TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
                int NHITS = int(hitVector.size());
                for (int IHIT = 0; IHIT < NHITS; ++IHIT) {
                  // Here we are trying to find if a hits are too close i.e. closer than _minDistToDelta
                  TrackerHit* trkhit1 = hit->getTrackerHit();
                  TrackerHit* trkhit2 = hitVector[IHIT]->getTrackerHit();
                  if (trkhit1->getCellID0() == trkhit2->getCellID0()) { // i.e. they are in the same sensor
                    float distance = 0.;
                    for (int iC = 0; iC < 3; ++iC) {
                      float posFirst = float(hit->getTrackerHit()->getPosition()[iC]);
                      float posSecond = float(hitVector[IHIT]->getTrackerHit()->getPosition()[iC]);
                      float deltaPos = posFirst - posSecond;
                      distance += deltaPos * deltaPos;
                    }
                    distance = sqrt(distance);
                    if (distance < _minDistToDelta) {
                      consider = false;
                      break;
                    }
                  }
                }
              }
              if (consider) {
                double phi0 = trackAR->getPhi();
                double d0 = trackAR->getD0();
                double z0 = trackAR->getZ0();
                double omega = trackAR->getOmega();
                double tanlambda = trackAR->getTanLambda();
                HelixClass_double helix;
                helix.Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _bField);
                double pos[3];
                double distance[3];
                for (int i = 0; i < 3; ++i)
                  pos[i] = hit->getTrackerHit()->getPosition()[i];
                float time = helix.getDistanceToPoint(pos, distance);
                if (time < 1.0e+10 && distance[2] < minDist) {
                  minDist = distance[2];
                  trkhitToAttach = hit;
                }
              }
            }
          }
        }
      }
      if (minDist < _minDistCutAttachForVXD && trkhitToAttach != NULL) {
        // int attachment_status = AttachHitToTrack_KalFit(trackAR,trkhitToAttach);//1 means attachment is done.
        AttachHitToTrack_KalFit(trackAR, trkhitToAttach); // 1 means attachment is done.
      }
    }
  } // end of nTrack mori

  return;
}

void FPCCDSiliconTracking_MarlinTrk::AttachRemainingVTXHitsFast() {
  std::vector<TrackerHitExtendedVec> nonAttachedHits;
  nonAttachedHits.resize(_nDivisionsInPhi * _nDivisionsInTheta);
  std::vector<TrackExtendedVec> trackVector;
  trackVector.resize(_nDivisionsInPhi * _nDivisionsInTheta);
  int nTracks = int(_trackImplVec.size());

  for (int iTrk = 0; iTrk < nTracks; ++iTrk) {
    TrackExtended* track = _trackImplVec[iTrk];
    double Phi = double(track->getPhi());
    if (Phi < 0)
      Phi = Phi + TWOPI;
    float tanlambda = track->getTanLambda();
    double cosTheta = double(tanlambda / sqrt(1 + tanlambda * tanlambda));
    int iPhi = int(Phi / _dPhi); // This way does hardly attach hits to low Pt tracks.
    int iTheta = int((cosTheta + double(1.0)) / _dTheta);
    int iCode = iPhi + _nDivisionsInPhi * iTheta;
    trackVector[iCode].push_back(track);
  }

  for (int il = 0; il < _nLayers; ++il) {
    for (int ip = 0; ip < _nDivisionsInPhi; ++ip) {
      for (int it = 0; it < _nDivisionsInTheta; ++it) {
        int iCode = il + _nLayers * ip + _nLayers * _nDivisionsInPhi * it;
        TrackerHitExtendedVec& hitVec = _sectors[iCode];
        int nH = int(hitVec.size());
        for (int iH = 0; iH < nH; ++iH) {
          TrackerHitExtended* hitExt = hitVec[iH];
          TrackExtendedVec& trackVec = hitExt->getTrackExtendedVec();
          if (trackVec.size() == 0) {
            TrackerHit* hit = hitExt->getTrackerHit();
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
            int iPhi = int(Phi / _dPhi);
            int iTheta = int((cosTheta + double(1.0)) / _dTheta);
            iCode = iPhi + _nDivisionsInPhi * iTheta;
            nonAttachedHits[iCode].push_back(hitExt);
          }
        }
      }
    }
  }

  for (int iT = 0; iT < _nDivisionsInTheta; ++iT) {
    for (int iP = 0; iP < _nDivisionsInPhi; ++iP) {
      int iCode = iP + _nDivisionsInPhi * iT;
      int nHits = int(nonAttachedHits[iCode].size());
      int iT1 = iT - 1;
      int iT2 = iT + 1;
      if (iT == 0) {
        iT1 = iT;
        iT2 = iT1 + 1;
      }
      if (iT == _nDivisionsInTheta - 1) {
        iT2 = iT;
        iT1 = iT2 - 1;
      }
      int iPHI[3];
      iPHI[0] = iP - 1;
      iPHI[1] = iP;
      iPHI[2] = iP + 1;
      if (iP == 0)
        iPHI[0] = _nDivisionsInPhi - 1;
      if (iP == _nDivisionsInPhi - 1)
        iPHI[2] = 0;

      for (int ihit = 0; ihit < nHits; ++ihit) {
        TrackerHitExtended* hit = nonAttachedHits[iCode][ihit];
        TrackExtended* trackToAttach = NULL;
        float minDist = 1.0e+6;

        for (int iTheta = iT1; iTheta < iT2 + 1; ++iTheta) {
          for (int indexP = 0; indexP < 3; ++indexP) {
            int iPhi = iPHI[indexP];
            int iCodeForTrack = iPhi + _nDivisionsInPhi * iTheta;
            int nTrk = int(trackVector[iCodeForTrack].size());
            for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
              TrackExtended* trackAR = trackVector[iCodeForTrack][iTrk];
              bool consider = true;
              if (_checkForDelta) {
                TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
                int NHITS = int(hitVector.size());
                for (int IHIT = 0; IHIT < NHITS; ++IHIT) {
                  // Here we are trying to find if a hits are too close i.e. closer than _minDistToDelta
                  TrackerHit* trkhit1 = hit->getTrackerHit();
                  TrackerHit* trkhit2 = hitVector[IHIT]->getTrackerHit();

                  if (trkhit1->getCellID0() == trkhit2->getCellID0()) { // i.e. they are in the same sensor
                    float distance = 0.;
                    for (int iC = 0; iC < 3; ++iC) {
                      float posFirst = float(hit->getTrackerHit()->getPosition()[iC]);
                      float posSecond = float(hitVector[IHIT]->getTrackerHit()->getPosition()[iC]);
                      float deltaPos = posFirst - posSecond;
                      distance += deltaPos * deltaPos;
                    }
                    distance = sqrt(distance);
                    if (distance < _minDistToDelta) {
                      consider = false;
                      break;
                    }
                  }
                }
              }
              if (consider) {
                double phi0 = trackAR->getPhi();
                double d0 = trackAR->getD0();
                double z0 = trackAR->getZ0();
                double omega = trackAR->getOmega();
                double tanlambda = trackAR->getTanLambda();
                HelixClass_double helix;
                helix.Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _bField);
                int layer = getLayerID(hit->getTrackerHit());
                if (layer > _minimalLayerToAttach) {
                  double pos[3];
                  for (int i = 0; i < 3; ++i)
                    pos[i] = hit->getTrackerHit()->getPosition()[i];
                  double distance[3];
                  double time = helix.getDistanceToPoint(pos, distance);
                  if (time < 1.0e+10) {
                    if (distance[2] < minDist) {
                      minDist = distance[2];
                      trackToAttach = trackAR;
                    }
                  }
                }
              }
            }
          }
        }
        if (minDist < _minDistCutAttachForVXD && trackToAttach != NULL) {
          int iopt = 2;
          AttachHitToTrack(trackToAttach, hit, iopt);
        }
      }
    }
  }
}

void FPCCDSiliconTracking_MarlinTrk::AttachRemainingVTXHitsSlow() {
  TrackerHitExtendedVec nonAttachedHits;
  nonAttachedHits.clear();

  for (int il = 0; il < _nLayers; ++il) {
    for (int ip = 0; ip < _nDivisionsInPhi; ++ip) {
      for (int it = 0; it < _nDivisionsInTheta; ++it) {
        int iCode = il + _nLayers * ip + _nLayers * _nDivisionsInPhi * it;
        TrackerHitExtendedVec& hitVec = _sectors[iCode];
        int nH = int(hitVec.size());
        for (int iH = 0; iH < nH; ++iH) {
          TrackerHitExtended* hit = hitVec[iH];
          TrackExtendedVec& trackVec = hit->getTrackExtendedVec();
          // if (trackVec.size()==0)
          // nonAttachedHits.push_back( hit );
          //-- allow hits that are only used in triplets to be re-attached
          unsigned int maxTrackSize = 0;
          for (unsigned int itrack = 0; itrack < trackVec.size(); itrack++) {
            TrackerHitExtendedVec hitVec_tmp = trackVec[itrack]->getTrackerHitExtendedVec();
            unsigned int isize = hitVec_tmp.size();
            if (isize > maxTrackSize)
              maxTrackSize = isize;
          }
          if (maxTrackSize <= 3) {
            streamlog_out(DEBUG1) << " Add non attached hit to list: id = " << hit->getTrackerHit()->id() << std::endl;
            nonAttachedHits.push_back(hit);
          }
        }
      }
    }
  }

  int nNotAttached = int(nonAttachedHits.size());

  int nTrk = int(_trackImplVec.size());
  for (int iHit = 0; iHit < nNotAttached; ++iHit) {
    TrackerHitExtended* hit = nonAttachedHits[iHit];
    streamlog_out(DEBUG1) << " Try hit: id = " << hit->getTrackerHit()->id() << std::endl;
    int layer = getLayerID(hit->getTrackerHit());
    if (layer > _minimalLayerToAttach) {
      double pos[3];
      for (int i = 0; i < 3; ++i)
        pos[i] = hit->getTrackerHit()->getPosition()[i];
      double minDist = 1e+10;
      TrackExtended* trackToAttach = NULL;
      for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
        TrackExtended* trackAR = _trackImplVec[iTrk];
        bool consider = true;
        if (_checkForDelta) {
          TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
          int NHITS = int(hitVector.size());
          for (int IHIT = 0; IHIT < NHITS; ++IHIT) {
            // Here we are trying to find if a hits are too close i.e. closer than _minDistToDelta
            TrackerHit* trkhit1 = hit->getTrackerHit();
            TrackerHit* trkhit2 = hitVector[IHIT]->getTrackerHit();

            if (trkhit1->getCellID0() == trkhit2->getCellID0()) { // i.e. they are in the same sensor

              double distance = 0.;
              for (int iC = 0; iC < 3; ++iC) {
                double posFirst = double(hit->getTrackerHit()->getPosition()[iC]);
                double posSecond = double(hitVector[IHIT]->getTrackerHit()->getPosition()[iC]);
                double deltaPos = posFirst - posSecond;
                distance += deltaPos * deltaPos;
              }
              distance = sqrt(distance);
              if (distance < _minDistToDelta) {
                consider = false;
                streamlog_out(DEBUG1) << " hit: id = " << hit->getTrackerHit()->id()
                                      << " condsidered delta together with hit " << trkhit2->id() << std::endl;
                break;
              }
            }
          }
        }
        if (consider) {
          HelixClass_double helix;
          double phi0 = trackAR->getPhi();
          double d0 = trackAR->getD0();
          double z0 = trackAR->getZ0();
          double omega = trackAR->getOmega();
          double tanlambda = trackAR->getTanLambda();
          helix.Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _bField);
          double distance[3];
          double time = helix.getDistanceToPoint(pos, distance);
          if (time < 1.0e+10) {
            if (distance[2] < minDist) {
              minDist = distance[2];
              trackToAttach = trackAR;
            }
          }
        }
      }
      if (minDist < _minDistCutAttachForVXD && trackToAttach != NULL) {
        int iopt = 2;
        streamlog_out(DEBUG1) << " Hit: id = " << hit->getTrackerHit()->id() << " : try attachement" << std::endl;
        AttachHitToTrack(trackToAttach, hit, iopt);
      } else {
        streamlog_out(DEBUG1) << " Hit: id = " << hit->getTrackerHit()->id() << " rejected due to distance cut of "
                              << _minDistCutAttachForVXD << " min distance = " << minDist << std::endl;
      }
    }
  }
}

void FPCCDSiliconTracking_MarlinTrk::AttachRemainingFTDHitsSlow() {
  TrackerHitExtendedVec nonAttachedHits;
  nonAttachedHits.clear();

  for (int iS = 0; iS < 2; ++iS) {
    for (unsigned int layer = 0; layer < _nlayersFTD; ++layer) {
      for (int ip = 0; ip < _nPhiFTD; ++ip) {
        int iCode = iS + 2 * layer + 2 * _nlayersFTD * ip;
        TrackerHitExtendedVec& hitVec = _sectorsFTD[iCode];
        int nH = int(hitVec.size());
        for (int iH = 0; iH < nH; ++iH) {
          TrackerHitExtended* hit = hitVec[iH];
          TrackExtendedVec& trackVec = hit->getTrackExtendedVec();
          if (trackVec.size() == 0)
            nonAttachedHits.push_back(hit);
        }
      }
    }
  }

  int nNotAttached = int(nonAttachedHits.size());

  int nTrk = int(_trackImplVec.size());
  for (int iHit = 0; iHit < nNotAttached; ++iHit) {
    TrackerHitExtended* hit = nonAttachedHits[iHit];
    double pos[3];
    for (int i = 0; i < 3; ++i)
      pos[i] = hit->getTrackerHit()->getPosition()[i];
    double minDist = 1e+10;
    TrackExtended* trackToAttach = NULL;
    for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
      TrackExtended* trackAR = _trackImplVec[iTrk];
      bool consider = true;
      TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
      int NHITS = int(hitVector.size());

      for (int IHIT = 0; IHIT < NHITS; ++IHIT) {
        // SJA:FIXME: check to see if allowing no hits in the same sensor vs no hits in the same layer works
        //        if (hit->getTrackerHit()->getType() == hitVector[IHIT]->getTrackerHit()->getType())
        if (hit->getTrackerHit()->getCellID0() == hitVector[IHIT]->getTrackerHit()->getCellID0()) {
          consider = false;
          break;
        }
      }

      if (consider) {
        HelixClass_double helix;
        double phi0 = trackAR->getPhi();
        double d0 = trackAR->getD0();
        double z0 = trackAR->getZ0();
        double omega = trackAR->getOmega();
        double tanlambda = trackAR->getTanLambda();
        if (tanlambda * double(getSideID(hit->getTrackerHit())) > 0) {
          helix.Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _bField);
          double distance[3];
          double time = helix.getDistanceToPoint(pos, distance);
          if (time < 1.0e+10) {
            if (distance[2] < minDist) {
              minDist = distance[2];
              trackToAttach = trackAR;
            }
          }
        }
      }
    }
    if (minDist < _minDistCutAttachForFTD && trackToAttach != NULL) {
      int iopt = 2;
      AttachHitToTrack(trackToAttach, hit, iopt);
    }
  }
}

void FPCCDSiliconTracking_MarlinTrk::AttachRemainingFTDHitsFast() {
  int nTrk = _trackImplVec.size();

  for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
    TrackExtended* trackAR = _trackImplVec[iTrk];
    HelixClass_double helix;
    double phi0 = trackAR->getPhi();
    double d0 = trackAR->getD0();
    double z0 = trackAR->getZ0();
    double omega = trackAR->getOmega();
    double tanlambda = trackAR->getTanLambda();
    helix.Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _bField);
    int iSemiSphere = 0;
    if (tanlambda > 0)
      iSemiSphere = 1;
    double ref[3];
    for (int i = 0; i < 3; ++i)
      ref[i] = helix.getReferencePoint()[i];
    // Start loop over FTD layers
    for (unsigned int layer = 0; layer < _nlayersFTD; layer++) {
      double ZL = _zLayerFTD[layer];
      if (iSemiSphere == 0)
        ZL = -ZL;
      double point[3];
      helix.getPointInZ(ZL, ref, point);
      double Phi = atan2(point[1], point[0]);
      if (Phi < 0)
        Phi = Phi + TWOPI;
      int iPhi = int(Phi / _dPhiFTD);
      double distMin = 1e+10;
      TrackerHitExtended* attachedHit = NULL;
      for (int iP = iPhi - 1; iP <= iPhi + 1; ++iP) {
        int iPP = iP;
        if (iP < 0)
          iPP = iP + _nPhiFTD;
        if (iP >= _nPhiFTD)
          iPP = iP - _nPhiFTD;
        int iCode = iSemiSphere + 2 * layer + 2 * _nlayersFTD * iPP;
        int nHits = int(_sectorsFTD[iCode].size());
        for (int iHit = 0; iHit < nHits; ++iHit) {
          TrackerHitExtended* hit = _sectorsFTD[iCode][iHit];
          bool consider = true;
          TrackerHitExtendedVec& hitVector = trackAR->getTrackerHitExtendedVec();
          int NHITS = int(hitVector.size());

          // SJA:FIXME: check to see if allowing no hits in the same sensor vs no hits in the same layer works
          for (int IHIT = 0; IHIT < NHITS; ++IHIT) {
            //            if (hit->getTrackerHit()->getType() == hitVector[IHIT]->getTrackerHit()->getType())
            if (hit->getTrackerHit()->getCellID0() == hitVector[IHIT]->getTrackerHit()->getCellID0()) {
              consider = false;
              break;
            }
          }

          if (consider) {
            double pos[3];
            for (int i = 0; i < 3; ++i) {
              pos[i] = hit->getTrackerHit()->getPosition()[i];
            }
            double distance[3];
            double time = helix.getDistanceToPoint(pos, distance);
            if (time < 1.0e+10) {
              if (distance[2] < distMin) {
                distMin = distance[2];
                attachedHit = hit;
              }
            }
          }
        }
      }
      if (distMin < _minDistCutAttachForFTD && attachedHit != NULL) {
        int iopt = 2;
        AttachHitToTrack(trackAR, attachedHit, iopt);
      }
    }
  }
}

void FPCCDSiliconTracking_MarlinTrk::TrackingInFTD() {
  int nComb = int(_CombinationsFTD.size()) / 3;

  for (int iComb = 0; iComb < nComb; ++iComb) {
    int nLS[3];
    nLS[0] = _CombinationsFTD[3 * iComb];
    nLS[1] = _CombinationsFTD[3 * iComb + 1];
    nLS[2] = _CombinationsFTD[3 * iComb + 2];

    for (int iS = 0; iS < 2; ++iS) { // loop over +z and -z

      //      std::cout << "Combinations : " << iS << " " << nLS[0] << " " << nLS[1] << " " << nLS[2] << std::endl;
      //      int iC = iS + 2*nLS[0];
      //      TrackerHitExtendedVec& hitVec = _sectorsFTD[iC];
      //      int nO = int(hitVec.size());
      //      iC = iS + 2*nLS[1];
      //      hitVec = _sectorsFTD[iC];
      //      int nM = int(hitVec.size());
      //      iC = iS + 2*nLS[2];
      //      hitVec = _sectorsFTD[iC];
      //      int nI = int(hitVec.size());
      //      std::cout << nO << " " << nM << " " << nI << std::endl;

      for (int ipOuter = 0; ipOuter < _nPhiFTD; ++ipOuter) {
        int ipMiddleLow = ipOuter - 1;
        int ipMiddleUp = ipOuter + 1;

        unsigned int iCodeOuter = iS + 2 * nLS[0] + 2 * _nlayersFTD * ipOuter;

        if (iCodeOuter >= _sectorsFTD.size()) {
          streamlog_out(ERROR) << "iCodeOuter index out of range: iCodeOuter =   " << iCodeOuter
                               << " _sectorsFTD.size() = " << _sectorsFTD.size() << " exit(1) called from file "
                               << __FILE__ << " line " << __LINE__ << std::endl;
          exit(1);
        }

        TrackerHitExtendedVec& hitVecOuter = _sectorsFTD[iCodeOuter];

        int nOuter = int(hitVecOuter.size());

        for (int iOuter = 0; iOuter < nOuter; ++iOuter) {
          TrackerHitExtended* hitOuter = hitVecOuter[iOuter];

          for (int ipMiddle = ipMiddleLow; ipMiddle <= ipMiddleUp; ++ipMiddle) {
            // for(int ipMiddle=0;ipMiddle<_nPhiFTD;++ipMiddle)
            int ipM = ipMiddle;
            if (ipM < 0)
              ipM = ipMiddle + _nPhiFTD;
            if (ipM >= _nPhiFTD)
              ipM = ipMiddle - _nPhiFTD;
            int iCodeMiddle = iS + 2 * nLS[1] + 2 * _nlayersFTD * ipM;

            TrackerHitExtendedVec& hitVecMiddle = _sectorsFTD[iCodeMiddle];
            int ipInnerLow, ipInnerUp;
            ipInnerLow = ipMiddle - 1;
            ipInnerUp = ipMiddle + 1;

            int nMiddle = int(hitVecMiddle.size());

            for (int iMiddle = 0; iMiddle < nMiddle; ++iMiddle) {
              TrackerHitExtended* hitMiddle = hitVecMiddle[iMiddle];
              for (int ipInner = ipInnerLow; ipInner <= ipInnerUp; ++ipInner) {
                // for (int ipInner=0;ipInner<_nPhiFTD;++ipInner)
                int ipI = ipInner;
                if (ipI < 0)
                  ipI = ipInner + _nPhiFTD;
                if (ipI >= _nPhiFTD)
                  ipI = ipInner - _nPhiFTD;
                int iCodeInner = iS + 2 * nLS[2] + 2 * _nlayersFTD * ipI;
                TrackerHitExtendedVec& hitVecInner = _sectorsFTD[iCodeInner];

                int nInner = int(hitVecInner.size());

                for (int iInner = 0; iInner < nInner; ++iInner) {
                  TrackerHitExtended* hitInner = hitVecInner[iInner];
                  HelixClass_double helix;
                  //                  std::cout << std::endl;
                  //                  std::cout << "Outer Hit Type " << hitOuter->getTrackerHit()->getType() << " z = "
                  //                  << hitOuter->getTrackerHit()->getPosition()[2]
                  //                  << "\nMiddle Hit Type "<< hitMiddle->getTrackerHit()->getType() << " z = " <<
                  //                  hitMiddle->getTrackerHit()->getPosition()[2]
                  //                  << "\nInner Hit Type "<< hitInner->getTrackerHit()->getType() << " z = " <<
                  //                  hitInner->getTrackerHit()->getPosition()[2]  << std::endl;

                  streamlog_out(DEBUG1) << " " << std::setw(3) << ipOuter << " " << std::setw(3) << ipMiddle << " "
                                        << std::setw(3) << ipInner << "       " << std::setw(3) << iS << "     "
                                        << std::setw(3) << nLS[0] << " " << std::setw(3) << nLS[1] << " "
                                        << std::setw(3) << nLS[2] << "     " << std::setw(3) << nOuter << " : "
                                        << std::setw(3) << nMiddle << " : " << std::setw(3) << nInner
                                        << "  :: " << std::setw(3) << nOuter * nMiddle * nInner << std::endl;

                  TrackExtended* trackAR = TestTriplet(hitOuter, hitMiddle, hitInner, helix, 1);
                  if (trackAR != NULL) {
                    //                    std::cout << "FTD triplet found" << std::endl;
                    int nHits = BuildTrackFTD(trackAR, nLS, iS);

                    _tracksWithNHitsContainer.getTracksWithNHitsVec(nHits).push_back(trackAR);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

int FPCCDSiliconTracking_MarlinTrk::BuildTrackFTD(TrackExtended* trackAR, int* nLR, int iS) {
  //  std::cout << "BuildTrackFTD: Layers = " << nLR[0] << " " << nLR[1] << " " << nLR[2] << std::endl;

  // initialise a helix from the track
  HelixClass_double helix;
  const double d0 = trackAR->getD0();
  const double z0 = trackAR->getZ0();
  const double phi0 = trackAR->getPhi();
  const double tanlambda = trackAR->getTanLambda();
  const double omega = trackAR->getOmega();
  helix.Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _bField);
  double ref[3] = {helix.getReferencePoint()[0], helix.getReferencePoint()[1], helix.getReferencePoint()[2]};

  for (int iL = 0; iL < static_cast<int>(_nlayersFTD); ++iL) {
    if (iL != nLR[0] && iL != nLR[1] && iL != nLR[2]) {
      double point[3];
      double ZL = _zLayerFTD[iL];
      if (iS == 0)
        ZL = -ZL;
      helix.getPointInZ(ZL, ref, point);
      //      double Phi = atan2(point[1],point[0]);
      //      int iPhi = int(Phi/_dPhiFTD);
      double distMin = 1e+6;
      TrackerHitExtended* attachedHit = NULL;
      for (int ip = 0; ip <= _nPhiFTD; ++ip) {
        int iP = ip;
        if (iP < 0)
          iP = ip + _nPhiFTD;
        if (iP >= _nPhiFTD)
          iP = ip - _nPhiFTD;
        int iCode = iS + 2 * iL + 2 * _nlayersFTD * iP;
        TrackerHitExtendedVec& hitVec = _sectorsFTD[iCode];
        int nH = int(hitVec.size());
        for (int iH = 0; iH < nH; ++iH) {
          TrackerHitExtended* hit = hitVec[iH];
          TrackerHit* trkHit = hit->getTrackerHit();
          double pos[3];
          for (int i = 0; i < 3; ++i)
            pos[i] = double(trkHit->getPosition()[i]);
          double distance[3];
          double time = helix.getDistanceToPoint(pos, distance);
          if (time < 1.0e+10) {
            if (distance[2] < distMin) {
              distMin = distance[2];
              attachedHit = hit;
            }
          }
        }
      }
      //      std::cout << "Layer = " << iL << "  distMin = " << distMin << std::endl;
      if (distMin < _minDistCutAttachForFTD && attachedHit != NULL) {
        int iopt = 2;
        AttachHitToTrack(trackAR, attachedHit, iopt);
      }
    }
  }
  TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();
  int nH = int(hitVec.size());
  return nH;
}

#if 1
// Here will be used in AttachRemainingVTXVeryFast which is under construction.

int FPCCDSiliconTracking_MarlinTrk::AttachHitToTrack_KalFit(TrackExtended* trackAR, TrackerHitExtended* hit) {
  TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();

  std::vector<std::pair<double, TrackerHit*>> r2_values;
  r2_values.reserve(hitVec.size() + 1);

  for (TrackerHitExtendedVec::iterator it = hitVec.begin(); it != hitVec.end(); ++it) {
    EVENT::TrackerHit* h = (*it)->getTrackerHit();
    double r2 = pow(h->getPosition()[0], 2) + pow(h->getPosition()[1], 2);
    r2_values.push_back(std::make_pair(r2, h));
  }

  std::sort(r2_values.begin(), r2_values.end());

  TrackerHitVec trkHits;
  trkHits.reserve(r2_values.size());
  for (std::vector<std::pair<double, EVENT::TrackerHit*>>::iterator it = r2_values.begin(); it != r2_values.end();
       ++it) {
    trkHits.push_back(it->second);
  }

  int attached = 0;
  float par[5];
  float epar[15];
  int ndf;
  float Chi2;
  TrackerHitVec hits_in_fit;
  TrackerHitVec outliers;
  HelixClass_double helix;
  int error = KalFit(ndf, Chi2, trkHits, hits_in_fit, outliers, par, epar, helix); // 0 means the fit succeed

  float d0 = par[0];
  float phi0 = par[1];
  float omega = par[2];
  float z0 = par[3];
  float tanL = par[4];

  if (error == 0 && Chi2 / float(ndf) < _chi2FitCut_kalman) {
    trackAR->addTrackerHitExtended(hit);
    if (_keepCandidate != true)
      hit->addTrackExtended(trackAR);
    /*  After AttachRemainingVTXHits, we don't use this information.
     *  If keepCandidate != true, the above hit is not remaining hit any more.
     */
    trackAR->setChi2(Chi2);
    trackAR->setOmega(omega);
    trackAR->setTanLambda(tanL);
    trackAR->setD0(d0);
    trackAR->setZ0(z0);
    trackAR->setPhi(phi0);
    trackAR->setNDF(ndf);
    trackAR->setCovMatrix(epar);
    attached = 1;
  }

  return attached;
}

#endif

int FPCCDSiliconTracking_MarlinTrk::AttachHitToTrack(TrackExtended* trackAR, TrackerHitExtended* hit, int iopt) {
  int attached = 0;
  TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();
  int nHits = int(hitVec.size());

  double* xh = new double[nHits + 1];
  double* yh = new double[nHits + 1];
  float* zh = new float[nHits + 1];
  double* wrh = new double[nHits + 1];
  float* wzh = new float[nHits + 1];
  float* rh = new float[nHits + 1];
  float* ph = new float[nHits + 1];
  float par[5];
  float epar[15];

  for (int i = 0; i < nHits; ++i) {
    TrackerHit* trkHit = hitVec[i]->getTrackerHit();
    xh[i] = double(trkHit->getPosition()[0]);
    yh[i] = double(trkHit->getPosition()[1]);
    zh[i] = float(trkHit->getPosition()[2]);
    ph[i] = float(atan2(yh[i], xh[i]));
    rh[i] = float(sqrt(xh[i] * xh[i] + yh[i] * yh[i]));
    float rR = hitVec[i]->getResolutionRPhi();
    float rZ = hitVec[i]->getResolutionZ();
    wrh[i] = double(1.0 / (rR * rR));
    wzh[i] = 1.0 / (rZ * rZ);
  }

  TrackerHit* trkHit = hit->getTrackerHit();
  xh[nHits] = double(trkHit->getPosition()[0]);
  yh[nHits] = double(trkHit->getPosition()[1]);
  zh[nHits] = float(trkHit->getPosition()[2]);
  ph[nHits] = float(atan2(yh[nHits], xh[nHits]));
  rh[nHits] = float(sqrt(xh[nHits] * xh[nHits] + yh[nHits] * yh[nHits]));

  float rR = hit->getResolutionRPhi();
  float rZ = hit->getResolutionZ();
  wrh[nHits] = double(1.0 / (rR * rR));
  wzh[nHits] = 1.0 / (rZ * rZ);

  int NPT = nHits + 1;

  // SJA:FIXME the newtonian part is giving crazy results for FTD so just use iopt 2 for simply attaching hits
  // using SIT and VXD doesn't seem to give any problems, so make it a function parameter and let the caller decide
  //  int iopt = 3;

  float chi2RPhi = 0;
  float chi2Z = 0;

  int error = _fastfitter->fastHelixFit(NPT, xh, yh, rh, ph, wrh, zh, wzh, iopt, par, epar, chi2RPhi, chi2Z);
  par[3] = par[3] * par[0] / fabs(par[0]);

  float omega = par[0];
  float tanlambda = par[1];
  float phi0 = par[2];
  float d0 = par[3];
  float z0 = par[4];
  float chi2 = FLT_MAX;
  int ndf = INT_MAX;

  if (NPT == 3) {
    chi2 = chi2RPhi * _chi2WRPhiTriplet + chi2Z * _chi2WZTriplet;
  }
  if (NPT == 4) {
    chi2 = chi2RPhi * _chi2WRPhiQuartet + chi2Z * _chi2WZQuartet;
  }
  if (NPT > 4) {
    chi2 = chi2RPhi * _chi2WRPhiSeptet + chi2Z * _chi2WZSeptet;
  }
  ndf = 2 * NPT - 5;

  bool normalcase = false;
  if (std::isnormal(chi2) && std::isnormal(ndf) && chi2 > 0 && ndf > 0) {
    normalcase = true;
  } // Mori added new requirement that chi2Min be normal value on September 1, 2013.
  if (error == 0 && chi2 / float(ndf) < _chi2FitCut && normalcase == true) {
    trackAR->addTrackerHitExtended(hit);
    hit->addTrackExtended(trackAR);
    trackAR->setChi2(chi2);
    trackAR->setOmega(omega);
    trackAR->setTanLambda(tanlambda);
    trackAR->setD0(d0);
    trackAR->setZ0(z0);
    trackAR->setPhi(phi0);
    trackAR->setNDF(ndf);
    trackAR->setCovMatrix(epar);
    attached = 1;
    streamlog_out(DEBUG1) << "Attachement succeeded chi2/float(ndf) = " << chi2 / float(ndf)
                          << "  cut = " << _chi2FitCut << " chi2RPhi = " << chi2RPhi << " chi2Z = " << chi2Z
                          << " error = " << error << std::endl;
  } else {
    streamlog_out(DEBUG1) << "Attachement failed chi2/float(ndf) = " << chi2 / float(ndf) << "  cut = " << _chi2FitCut
                          << " chi2RPhi = " << chi2RPhi << " chi2Z = " << chi2Z << " error = " << error << std::endl;
  }

  delete[] xh;
  delete[] yh;
  delete[] zh;
  delete[] wrh;
  delete[] wzh;
  delete[] rh;
  delete[] ph;

  return attached;
}

void FPCCDSiliconTracking_MarlinTrk::FinalRefit(LCCollectionVec* trk_col, LCCollectionVec* /*rel_col*/) {
  int nTracks = int(_trackImplVec.size());

  int nSiSegments = 0;
  double eTot = 0.;
  double pxTot = 0.;
  double pyTot = 0.;
  double pzTot = 0.;

  for (int iTrk = 0; iTrk < nTracks; ++iTrk) {
    TrackExtended* trackAR = _trackImplVec[iTrk];
    TrackerHitExtendedVec& hitVec = trackAR->getTrackerHitExtendedVec();

    int nHits = int(hitVec.size());

    if (nHits >= _minimalHits) {
      //    int * lh = new int[nHits];
      std::vector<int> lh;
      lh.resize(nHits);

      for (int i = 0; i < nHits; ++i) {
        lh[i] = 0;
      }

      double d0 = trackAR->getD0();
      double z0 = trackAR->getZ0();
      double omega = trackAR->getOmega();
      double tanlambda = trackAR->getTanLambda();
      double phi0 = trackAR->getPhi();

      HelixClass_double* helix = new HelixClass_double();
      helix->Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _bField);

      // get the point of closest approach to the reference point
      // here it is implicitly assumed that the reference point is the origin
      double Pos[3];
      Pos[0] = -d0 * sin(phi0);
      Pos[1] = d0 * cos(phi0);
      Pos[2] = z0;

      // at this point is is possible to have hits from the same layer ...
      // so a check is made to ensure that the hit with the smallest distance to the
      // current helix hypothosis is used, the other hit has lh set to 0

      // start loop over the hits to
      for (int ihit = 0; ihit < nHits; ++ihit) {
        lh[ihit] = 1; // only hits which have lh=1 will be used for the fit

        // get the pointer to the lcio trackerhit for this hit
        TrackerHit* trkHit = hitVec[ihit]->getTrackerHit();

        int det = getDetectorID(trkHit);

        if (det == lcio::ILDDetID::VXD || det == lcio::ILDDetID::FTD ||
            det == lcio::ILDDetID::SIT) { // only accept VXD, FTD or SIT

          //        int layer = getLayerID(trkHit);
          //        int moduleIndex = getModuleID(trkHit);

          // start a double loop over the hits which have already been checked
          for (int lhit = 0; lhit < ihit; ++lhit) {
            // get the pointer to the lcio trackerhit for the previously checked hit
            TrackerHit* trkHitS = hitVec[lhit]->getTrackerHit();

            //          int layerS = getLayerID(trkHitS);
            //          int moduleIndexS = getModuleID(trkHitS);

            // SJA:FIXME: check to see if allowing no hits in the same sensor vs no hits in the same layer works
            // if they are on the same layer and the previously checked hits has been declared good for fitting
            //          if ((trkHitS->getType() == trkHit->getType()) && (lh[lhit] == 1))
            // check if the hits have the same layer and petal number
            //          hitVec[ihit]->
            //          if ((layer == layerS) && (moduleIndex==moduleIndexS) && (lh[lhit] == 1))
            if ((trkHit->getCellID0() == trkHitS->getCellID0()) && (lh[lhit] == 1)) {
              // get the position of the hits
              double xP[3];
              double xPS[3];
              for (int iC = 0; iC < 3; ++iC) {
                xP[iC] = double(trkHit->getPosition()[iC]);
                xPS[iC] = double(trkHitS->getPosition()[iC]);
              }

              // get the intersection of the helix with the either the cylinder or plane containing the hit
              double Point[6];
              double PointS[6];

              if (det == lcio::ILDDetID::FTD) {
                // double time =
                helix->getPointInZ(xP[2], Pos, Point);
                // double time =
                helix->getPointInZ(xPS[2], Pos, PointS);

              } else {
                double RAD = sqrt(xP[0] * xP[0] + xP[1] * xP[1]);
                double RADS = sqrt(xPS[0] * xPS[0] + xPS[1] * xPS[1]);
                // double time =
                helix->getPointOnCircle(RAD, Pos, Point);
                // double time =
                helix->getPointOnCircle(RADS, Pos, PointS);
              }

              double DIST = 0;
              double DISTS = 0;

              // get the euclidean distance between the hit and the point of intersection
              for (int iC = 0; iC < 3; ++iC) {
                DIST += (Point[iC] - xP[iC]) * (Point[iC] - xP[iC]);
                DISTS += (PointS[iC] - xPS[iC]) * (PointS[iC] - xPS[iC]);
              }
              if (DIST < DISTS) {
                lh[lhit] = 0;
              } else {
                lh[ihit] = 0;
              }
              break;
            }
          }
        }
      }

      delete helix;

      EVENT::TrackerHitVec trkHits;
      EVENT::TrackerHitVec trkHits_used_inFit;

      int nFit = 0;
      for (int i = 0; i < nHits; ++i) {
        // check if the hit has been rejected as being on the same layer and further from the helix lh==0
        if (lh[i] == 1) {
          TrackerHit* trkHit = hitVec[i]->getTrackerHit();
          nFit++;
          if (trkHit) {
            trkHits.push_back(trkHit);
          } else {
            throw EVENT::Exception(
                std::string("FPCCDSiliconTracking_MarlinTrk::FinalRefit: TrackerHit pointer == NULL "));
          }
        } else { // reject hit
          // SJA:FIXME missuse of type find a better way to signal rejected hits
          hitVec[i]->setType(int(0));
        }
      }

      if (trkHits.size() < 3) {
        streamlog_out(DEBUG3)
            << "FPCCDSiliconTracking_MarlinTrk::FinalRefit: Cannot fit less than 3 hits. Number of hits =  "
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

      //      for (unsigned ihit_indx=0 ; ihit_indx < trkHits.size(); ++ihit_indx) {
      //        std::cout << " trk hit " << *trkHits[ihit_indx] << std::endl;
      //      }

      bool fit_backwards = IMarlinTrack::backward;

      MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();

      int error = 0;

      try {
        error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, Track, fit_backwards, covMatrix, _bField,
                                                    _maxChi2PerHit);

      } catch (...) {
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

      UTIL::BitField64 cellID_encoder(lcio::LCTrackerCellID::encoding_string());

      MarlinTrk::addHitNumbersToTrack(Track, all_hits, true, cellID_encoder);

      marlinTrk->getOutliers(outliers);

      for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
        all_hits.push_back(outliers[ihit].first);
      }

      MarlinTrk::addHitNumbersToTrack(Track, all_hits, false, cellID_encoder);

      delete marlinTrk;

      int nhits_in_vxd = Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::VXD - 2];
      int nhits_in_ftd = Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::FTD - 2];
      int nhits_in_sit = Track->subdetectorHitNumbers()[2 * lcio::ILDDetID::SIT - 2];

      streamlog_out(DEBUG2) << " Hit numbers for Track " << Track->id() << ": "
                            << " vxd hits = " << nhits_in_vxd << " ftd hits = " << nhits_in_ftd
                            << " sit hits = " << nhits_in_sit << std::endl;

      if (nhits_in_vxd > 0)
        Track->setTypeBit(lcio::ILDDetID::VXD);
      if (nhits_in_ftd > 0)
        Track->setTypeBit(lcio::ILDDetID::FTD);
      if (nhits_in_sit > 0)
        Track->setTypeBit(lcio::ILDDetID::SIT);

      if (error != IMarlinTrack::success) {
        delete Track;
        streamlog_out(DEBUG3) << "FPCCDSiliconTracking_MarlinTrk::FinalRefit: Track fit failed with error code "
                              << error << " track dropped. Number of hits = " << trkHits.size() << std::endl;
        continue;
      }

      if (Track->getNdf() <= 0) {
        delete Track;
        streamlog_out(DEBUG3) << "FPCCDSiliconTracking_MarlinTrk::FinalRefit: Track fit returns " << Track->getNdf()
                              << " degress of freedom track dropped. Number of hits = " << trkHits.size() << std::endl;
        continue;
      }

      const TrackState* trkStateIP = Track->getTrackState(lcio::TrackState::AtIP);

      if (trkStateIP == 0) {
        streamlog_out(DEBUG3) << "FPCCDSiliconTracking_MarlinTrk::FinalRefit: Track fit returns " << Track->getNdf()
                              << " degress of freedom track dropped. Number of hits = " << trkHits.size() << std::endl;
        throw EVENT::Exception(std::string("FPCCDSiliconTracking_MarlinTrk::FinalRefit: trkStateIP pointer == NULL "));
        delete Track;
        continue;
      }

      trk_col->addElement(Track);

      // note trackAR which is of type TrackExtended, only takes fits set for ref point = 0,0,0
      trackAR->setOmega(trkStateIP->getOmega());
      trackAR->setTanLambda(trkStateIP->getTanLambda());
      trackAR->setPhi(trkStateIP->getPhi());
      trackAR->setD0(trkStateIP->getD0());
      trackAR->setZ0(trkStateIP->getZ0());

      float cov[15];

      for (int i = 0; i < 15; ++i) {
        cov[i] = trkStateIP->getCovMatrix().operator[](i);
      }

      trackAR->setCovMatrix(cov);
      trackAR->setChi2(Track->getChi2());
      trackAR->setNDF(Track->getNdf());

      nSiSegments++;

      HelixClass_double helix_final;

      helix_final.Initialize_Canonical(trkStateIP->getPhi(), trkStateIP->getD0(), trkStateIP->getZ0(),
                                       trkStateIP->getOmega(), trkStateIP->getTanLambda(), _bField);

      double trkPx = helix_final.getMomentum()[0];
      double trkPy = helix_final.getMomentum()[1];
      double trkPz = helix_final.getMomentum()[2];
      double trkP = sqrt(trkPx * trkPx + trkPy * trkPy + trkPz * trkPz);
      eTot += trkP;
      pxTot += trkPx;
      pyTot += trkPy;
      pzTot += trkPz;
    }
  }

  streamlog_out(DEBUG4) << "FPCCDSiliconTracking_MarlinTrk -> run " << _nRun << " event " << _nEvt << std::endl;
  streamlog_out(DEBUG4) << "Number of Si tracks = " << nSiSegments << std::endl;
  streamlog_out(DEBUG4) << "Total 4-momentum of Si tracks : E = " << eTot << " Px = " << pxTot << " Py = " << pyTot
                        << " Pz = " << pzTot << std::endl;
}

void FPCCDSiliconTracking_MarlinTrk::setupGeom(const dd4hep::Detector& theDetector) {
  double bFieldVec[3];
  theDetector.field().magneticField({0, 0, 0}, bFieldVec); // get the magnetic field vector from DD4hep
  _bField = bFieldVec[2] / dd4hep::tesla;                  // z component at (0,0,0)

  //-- VXD Parameters--
  _nLayersVTX = 0;

  try {
    streamlog_out(DEBUG9) << " filling VXD parameters  " << std::endl;

    dd4hep::DetElement vtxDE = theDetector.detector("VXD");
    dd4hep::rec::ZPlanarData* vtx = vtxDE.extension<dd4hep::rec::ZPlanarData>();
    _nLayersVTX = vtx->layers.size();

  } catch (std::runtime_error& e) {
    streamlog_out(DEBUG9) << " ### VXD detector Not Present in Detector" << std::endl;
  }

  //-- SIT Parameters--
  _nLayersSIT = 0;

  try {
    streamlog_out(DEBUG9) << " filling SIT parameters  " << std::endl;

    dd4hep::DetElement sitDE = theDetector.detector("SIT");
    dd4hep::rec::ZPlanarData* sit = sitDE.extension<dd4hep::rec::ZPlanarData>();
    _nLayersSIT = sit->layers.size();
  } catch (std::runtime_error& e) {
    streamlog_out(DEBUG9) << " ###  SIT detector Not Present in Detector " << std::endl;
  }

  //-- FTD Parameters--
  _petalBasedFTDWithOverlaps = false;
  _nlayersFTD = 0;

  try {
    streamlog_out(DEBUG9) << " filling FTD parameters  " << std::endl;

    dd4hep::DetElement ftdDE = theDetector.detector("FTD");
    dd4hep::rec::ZDiskPetalsData* ftd = ftdDE.extension<dd4hep::rec::ZDiskPetalsData>();

    _nlayersFTD = ftd->layers.size();

    for (unsigned int disk = 0; disk < _nlayersFTD; ++disk) {
      _zLayerFTD.push_back(ftd->layers[disk].zPosition + ftd->layers[disk].zOffsetSensitive);
      _zLayerFTD.push_back(ftd->layers[disk].zPosition - ftd->layers[disk].zOffsetSensitive);
      _petalBasedFTDWithOverlaps = true;
    }

    // SJA: Here we increase the size of _nlayersFTD as we are treating the
    _nlayersFTD = _zLayerFTD.size();

  } catch (std::runtime_error& e) {
    streamlog_out(DEBUG9) << " ### FTD detector Not Present in Detector" << std::endl;
  }
}

void FPCCDSiliconTracking_MarlinTrk::TracksWithNHitsContainer::clear() {
  for (std::vector<TrackExtendedVec>::iterator trackVecIter = _tracksNHits.begin(); trackVecIter < _tracksNHits.end();
       trackVecIter++) {
    for (TrackExtendedVec::iterator trackIter = trackVecIter->begin(); trackIter < trackVecIter->end(); trackIter++) {
      delete *trackIter;
    }

    trackVecIter->clear();
  }
}

LCCollection* FPCCDSiliconTracking_MarlinTrk::GetCollection(LCEvent* evt, std::string colName) {
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

LCRelationNavigator* FPCCDSiliconTracking_MarlinTrk::GetRelations(LCEvent* evt, std::string RelName) {
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

void FPCCDSiliconTracking_MarlinTrk::InitVXDGeometry(const dd4hep::Detector& theDetector) {
  // Save frequently used parameters.

  dd4hep::DetElement vtxDE = theDetector.detector("VXD");
  dd4hep::rec::ZPlanarData* vtx = vtxDE.extension<dd4hep::rec::ZPlanarData>();

  _vxd.nLayer = vtx->layers.size();

  _vxd.geodata.resize(_vxd.nLayer);
  _vxd.maxLadder = 0;

  for (int ly = 0; ly < _vxd.nLayer; ly++) {
    const dd4hep::rec::ZPlanarData::LayerLayout layerVXD = vtx->layers[ly];

    _vxd.geodata[ly].nladder = layerVXD.ladderNumber; // Number of ladders in this layer
    if (_vxd.maxLadder < _vxd.geodata[ly].nladder) {
      _vxd.maxLadder = _vxd.geodata[ly].nladder;
    }
    _vxd.geodata[ly].rmin = layerVXD.distanceSensitive; // Distance of sensitive area from IP
    _vxd.geodata[ly].dphi = (2 * M_PI) / (double)_vxd.geodata[ly].nladder;
    _vxd.geodata[ly].phi0 = layerVXD.phi0; // phi offset.
    _vxd.geodata[ly].sthick = layerVXD.thicknessSensitive;
    _vxd.geodata[ly].sximin = -layerVXD.offsetSensitive - layerVXD.widthSensitive / 2.0;
    _vxd.geodata[ly].sximax = -layerVXD.offsetSensitive + layerVXD.widthSensitive / 2.0;
    // CAUTION: plus direction of xi-axis in ld = 0, ly = 0 is equal to minus direction of global x-axis and longer than
    // minus direction of xi-axis.
    _vxd.geodata[ly].hlength = layerVXD.zHalfSensitive;
    _vxd.geodata[ly].cosphi.resize(_vxd.geodata[ly].nladder);
    _vxd.geodata[ly].sinphi.resize(_vxd.geodata[ly].nladder);
    _vxd.geodata[ly].phi.resize(_vxd.geodata[ly].nladder);
    _vxd.geodata[ly].phiAtXiMin.resize(_vxd.geodata[ly].nladder);
    _vxd.geodata[ly].phiAtXiMax.resize(_vxd.geodata[ly].nladder);
    _vxd.geodata[ly].ladder_incline.resize(_vxd.geodata[ly].nladder);
    _vxd.geodata[ly].num_xi_pixel = (int)(layerVXD.widthSensitive / _pixelSizeVec[ly]);
    _vxd.geodata[ly].num_zeta_pixel = (int)(2 * layerVXD.zHalfSensitive / _pixelSizeVec[ly]);

    if (ly % 2 == 0) {
      _vxd.geodata[ly].rmes = _vxd.geodata[ly].rmin + _vxd.geodata[ly].sthick * (15.0 / 50.0) * 0.5;
    } else {
      _vxd.geodata[ly].rmes = _vxd.geodata[ly].rmin + _vxd.geodata[ly].sthick * (1.0 - (15.0 / 50.0) * 0.5);
    }

    for (int ld = 0; ld < _vxd.geodata[ly].nladder; ld++) {
      double phi = _vxd.geodata[ly].phi0 + _vxd.geodata[ly].dphi * ld;
      _vxd.geodata[ly].cosphi[ld] = cos(phi);
      _vxd.geodata[ly].sinphi[ld] = sin(phi);
      _vxd.geodata[ly].phi[ld] = phi;

      double incline = phi - (M_PI / 2.0);
      while (incline > 1.0 * M_PI || incline < -1.0 * M_PI) {
        incline > 1.0 * M_PI ? incline -= 2.0 * M_PI : incline += 2.0 * M_PI;
      }
      _vxd.geodata[ly].ladder_incline[ld] = incline; // this expresses the tilt of the ladder by phi.
      // I adjusted the range between [-180deg, 180deg]

      _vxd.geodata[ly].phiAtXiMin[ld] =
          _vxd.geodata[ly].phi[ld] + atan2(-_vxd.geodata[ly].sximin, _vxd.geodata[ly].rmes);
      _vxd.geodata[ly].phiAtXiMax[ld] =
          _vxd.geodata[ly].phi[ld] - atan2(_vxd.geodata[ly].sximax, _vxd.geodata[ly].rmes);
    }
  }
}

void FPCCDSiliconTracking_MarlinTrk::InitSITGeometry(const dd4hep::Detector& theDetector) {
  // Save frequently used parameters.

  dd4hep::DetElement sitDE = theDetector.detector("SIT");
  dd4hep::rec::ZPlanarData* sit = sitDE.extension<dd4hep::rec::ZPlanarData>();

  _sit.nLayer = sit->layers.size();

  _sit.geodata.resize(_sit.nLayer);
  _sit.maxLadder = 0;

  for (int ly = 0; ly < _sit.nLayer; ly++) {
    const dd4hep::rec::ZPlanarData::LayerLayout layerSIT = sit->layers[ly];

    _sit.geodata[ly].nladder = layerSIT.ladderNumber; // Number of ladders in this layer
    if (_sit.maxLadder < _sit.geodata[ly].nladder) {
      _sit.maxLadder = _sit.geodata[ly].nladder;
    }
    _sit.geodata[ly].rmin = layerSIT.distanceSensitive; // Distance of sensitive area from IP
    _sit.geodata[ly].phi0 = layerSIT.phi0;              // phi offset.
    _sit.geodata[ly].hlength = layerSIT.zHalfSensitive;
    _sit.geodata[ly].cosphi.resize(_sit.geodata[ly].nladder);
    _sit.geodata[ly].sinphi.resize(_sit.geodata[ly].nladder);
    _sit.geodata[ly].phi.resize(_sit.geodata[ly].nladder);
    _sit.geodata[ly].phiAtXiMin.resize(_sit.geodata[ly].nladder);
    _sit.geodata[ly].phiAtXiMax.resize(_sit.geodata[ly].nladder);
    _sit.geodata[ly].ladder_incline.resize(_sit.geodata[ly].nladder);

    for (int ld = 0; ld < _sit.geodata[ly].nladder; ld++) {
      double phi = _sit.geodata[ly].phi0 + _sit.geodata[ly].dphi * ld;
      _sit.geodata[ly].cosphi[ld] = cos(phi);
      _sit.geodata[ly].sinphi[ld] = sin(phi);
      _sit.geodata[ly].phi[ld] = phi;

      double incline = phi - (M_PI / 2.0);
      while (incline > 1.0 * M_PI || incline < -1.0 * M_PI) {
        incline += (incline > 1.0 * M_PI) ? -2.0 * M_PI : 2.0 * M_PI;
      }
      _sit.geodata[ly].ladder_incline[ld] = incline;
    }
  }
}

int FPCCDSiliconTracking_MarlinTrk::KalFit(int& ndf, float& Chi2, TrackerHitVec trkHits, TrackerHitVec& Hits_in_fit,
                                           TrackerHitVec& Outliers, float* par, float* epar, HelixClass_double& helix) {
  TrackImpl* Track = new TrackImpl;
  EVENT::FloatVec covMatrix;
  covMatrix.resize(15);
  covMatrix[0] = (_initialTrackError_d0);    // sigma_d0^2
  covMatrix[2] = (_initialTrackError_phi0);  // sigma_phi0^2
  covMatrix[5] = (_initialTrackError_omega); // sigma_omega^2
  covMatrix[9] = (_initialTrackError_z0);    // sigma_z0^2
  covMatrix[14] = (_initialTrackError_tanL); // sigma_tanl^2
  bool fit_backwards = IMarlinTrack::backward;
  MarlinTrk::IMarlinTrack* marlinTrk = _trksystem->createTrack();
  int error = 0;

  error = MarlinTrk::createFinalisedLCIOTrack(marlinTrk, trkHits, Track, fit_backwards, covMatrix, _bField,
                                              _maxChi2PerHit2nd);
  if (error != 0) {
    if (_mydebugKalFit)
      std::cout << "KalFit error code : " << error << std::endl;
    delete Track;
    delete marlinTrk;
    return error;
    /* For Reference,
       const int IMarlinTrack::success  = 0 ;  // no error
       const int IMarlinTrack::error = 1 ; // ndf is not enough
       const int IMarlinTrack::bad_intputs = 3 ;
       const int IMarlinTrack::no_intersection = 4 ; // no intersection found
       const int IMarlinTrack::site_discarded = 5 ;  // measurement discarded by the fitter
       const int IMarlinTrack::site_fails_chi2_cut = 6 ;  // measurement discarded by the fitter due to chi2 cut
       const int IMarlinTrack::all_sites_fail_fit = 7 ;   // no single measurement added to the fit
     */
  }

  ndf = Track->getNdf();
  if (std::isnormal(ndf) == false || ndf <= 0) {
    if (_mydebugKalFit)
      std::cout << "ERROR = 8!! ndf is " << ndf << " at KalFit" << std::endl;
    delete Track;
    delete marlinTrk;
    return 8;
  }

  Chi2 = Track->getChi2();
  if (std::isnormal(Chi2) == false || Chi2 < 0.0) {
    if (_mydebugKalFit)
      std::cout << "ERROR = 9!! Chi2 is " << Chi2 << " at KalFit" << std::endl;
    delete Track;
    delete marlinTrk;
    return 9;
  }

  const TrackState* trkStateIP = Track->getTrackState(lcio::TrackState::AtIP);
  if (trkStateIP == 0) {
    if (_mydebugKalFit)
      std::cout << "ERROR = 2!! trkStateIP is NULL! : KalFit" << std::endl;
    delete Track;
    delete marlinTrk;
    return 2;
  }
  par[0] = trkStateIP->getD0();
  par[1] = trkStateIP->getPhi();
  par[2] = trkStateIP->getOmega();
  par[3] = trkStateIP->getZ0();
  par[4] = trkStateIP->getTanLambda();

  if (std::isnormal(par[0]) == false || std::isnormal(par[1]) == false || std::isnormal(par[2]) == false ||
      std::isnormal(par[3]) == false || std::isnormal(par[4]) == false) {
    if (_mydebugKalFit) {
      std::cout << "ERROR = 10!! Some of track parameters output from Kalman Filter are nan or inf. " << std::endl;
      std::cout << "d0,phi0,omega,z0,tanlambda : " << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " "
                << par[4] << std::endl;
    }
    delete Track;
    delete marlinTrk;
    return 10;
  }

  FloatVec covM = trkStateIP->getCovMatrix();

  for (int i = 0; i < 15; i++)
    epar[i] = covM[i];
  helix.Initialize_Canonical(par[1], par[0], par[3], par[2], par[4], _bField);

  std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
  std::vector<std::pair<EVENT::TrackerHit*, double>> outliers;

  marlinTrk->getHitsInFit(hits_in_fit);
  Hits_in_fit.clear();
  Hits_in_fit.reserve(hits_in_fit.size());
  for (unsigned ihit = 0; ihit < hits_in_fit.size(); ++ihit) {
    Hits_in_fit.push_back(hits_in_fit[ihit].first);
  }

  marlinTrk->getOutliers(outliers);
  if (_mydebugKalFit) {
    std::cout << "outliers.size() : " << outliers.size() << std::endl;
    for (int i = 0; i < int(outliers.size()); i++) {
      std::cout << outliers[i].first << std::endl;
    }
  }
  Outliers.clear();
  Outliers.reserve(outliers.size());
  for (unsigned ihit = 0; ihit < outliers.size(); ++ihit) {
    Outliers.push_back(outliers[ihit].first);
  }

  delete Track;
  delete marlinTrk;

  return error; // 0
}

float FPCCDSiliconTracking_MarlinTrk::DotOf2Clusters(TrackerHit* A, TrackerHit* B) {
  ClusterStatus ca(A);
  ClusterStatus cb(B);

#if 0
  std::cout << "DotOf2Clusters!" << std::endl;
  std::cout << "--A--" << std::endl;
  std::cout <<  "layer : "  << ca.layer << std::endl;
  std::cout <<  "ladder : "  << ca.ladder << std::endl;
  std::cout <<  "xiwidth : "  << ca.xiwidth << std::endl;
  std::cout <<  "zetawidth : "  << ca.zetawidth << std::endl;
  std::cout <<  "nPix : "  << ca.nPix << std::endl;
  std::cout <<  "tilt : "  << ca.tilt << std::endl;
  std::cout << "--B--" << std::endl;
  std::cout <<  "layer : "  << cb.layer << std::endl;
  std::cout <<  "ladder : "  << cb.ladder << std::endl;
  std::cout <<  "xiwidth : "  << cb.xiwidth << std::endl;
  std::cout <<  "zetawidth : "  << cb.zetawidth << std::endl;
  std::cout <<  "nPix : "  << cb.nPix << std::endl;
  std::cout <<  "tilt : "  << cb.tilt << std::endl;
#endif

  double signA = 1;
  if (ca.tilt > 1)
    signA = -1;
  double signB = 1;
  if (cb.tilt > 1)
    signB = -1;
  TVector3 ldirA(_pixelSizeVec[ca.layer] * (ca.xiwidth - 1), signA * _pixelSizeVec[ca.layer] * (ca.zetawidth - 1),
                 _pixelheight);
  TVector3 revldirA(-ldirA.X(), -ldirA.Y(), ldirA.Z());
  TVector3 ldirB(_pixelSizeVec[cb.layer] * (cb.xiwidth - 1), signB * _pixelSizeVec[cb.layer] * (cb.zetawidth - 1),
                 _pixelheight);
  TVector3 gdirA = LocalToGlobal(ldirA, ca.layer, ca.ladder);
  TVector3 revgdirA = LocalToGlobal(revldirA, ca.layer, ca.ladder);
  TVector3 gdirB = LocalToGlobal(ldirB, cb.layer, cb.ladder);
  float dot1 = (revgdirA.Unit().Dot(gdirB.Unit()) > gdirA.Unit().Dot(gdirB.Unit())) ? revgdirA.Unit().Dot(gdirB.Unit())
                                                                                    : gdirA.Unit().Dot(gdirB.Unit());

  return dot1;
}

TVector3 FPCCDSiliconTracking_MarlinTrk::LocalToGlobal(TVector3 local, int layer, int ladder) {
#if 0
  std::cout << "========CHECK LocalToGlobal===========" << std::endl;
  std::cout << "check. x " << local.X() << std::endl;
  std::cout << "check. y " << local.Y() << std::endl;
  std::cout << "check. z " << local.Z() << std::endl;
  std::cout << "sximin = " << _vxd.geodata[layer].sximin << std::endl;
  std::cout << "sinphi = " << _vxd.geodata[layer].sinphi[ladder] << std::endl;
  std::cout << "cosphi = " << _vxd.geodata[layer].cosphi[ladder] << std::endl;
  std::cout << "hlength = " << _vxd.geodata[layer].hlength << std::endl << std::endl;
#endif

  TVector3 gdir(local.X() * _vxd.geodata[layer].sinphi[ladder] + local.Z() * _vxd.geodata[layer].cosphi[ladder],
                -local.X() * _vxd.geodata[layer].cosphi[ladder] + local.Z() * _vxd.geodata[layer].sinphi[ladder],
                local.Y());
  return gdir;
}

void FPCCDSiliconTracking_MarlinTrk::calcTrackParameterOfMCP(MCParticle* pmcp, double* par) {
  HelixTrack helixMC(pmcp->getVertex(), pmcp->getMomentum(), pmcp->getCharge(), _bField);
  double oldphi0 = double(helixMC.getPhi0());
  double omega = double(helixMC.getOmega());
  double tanL = double(helixMC.getTanLambda());
  double oldsinphi0 = std::sin(oldphi0);
  double oldcosphi0 = std::cos(oldphi0);
  double centx = pmcp->getVertex()[0] + 1 / omega * oldcosphi0;
  double centy = pmcp->getVertex()[1] + 1 / omega * oldsinphi0;
  double newphi0;
  if (omega > 0) {
    newphi0 = std::atan2(centy, centx);
  } else {
    newphi0 = std::atan2(-centy, -centx);
  }
  double sinphi0 = std::sin(newphi0);
  double cosphi0 = std::cos(newphi0);
  double d0 = centx * cosphi0 + centy * sinphi0 - 1.0 / omega;
  double z0 = pmcp->getVertex()[2] - 1.0 / omega * (newphi0 - oldphi0) * tanL;

#if 0
  std::cout << "Test Pivot Transformation " << std::endl;
  std::cout << "d0 : " << d0 << std::endl;
  std::cout << "z0 : " << z0 << std::endl;
  std::cout << "omega : " << omega << std::endl;
  std::cout << "phi0 : " << newphi0 << std::endl;
  std::cout << "tanL : " << tanL << std::endl;
#endif
  par[0] = d0;
  par[1] = z0;
  par[2] = omega;
  par[3] = newphi0;
  par[4] = tanL;
}

int FPCCDSiliconTracking_MarlinTrk::CheckTiltOf2Clusters(TrackerHit* A, TrackerHit* B, int /*level*/) {
  // level is set for the future where users choose the level of requirement hardness.
  // For now, this level is not alive.

  ClusterStatus ca(A);
  ClusterStatus cb(B);

#if 0
  std::cout << "--A--" << std::endl;
  std::cout <<  "xiwidth : "  << ca.xiwidth << std::endl;
  std::cout <<  "zetawidth : "  << ca.zetawidth << std::endl;
  std::cout <<  "nPix : "  << ca.nPix << std::endl;
  std::cout <<  "tilt : "  << ca.tilt << std::endl;
  std::cout << "--B--" << std::endl;
  std::cout <<  "xiwidth : "  << cb.xiwidth << std::endl;
  std::cout <<  "zetawidth : "  << cb.zetawidth << std::endl;
  std::cout <<  "nPix : "  << cb.nPix << std::endl;
  std::cout <<  "tilt : "  << cb.tilt << std::endl;
#endif

  // tilt check.
  // tilt : 0 --> straight line shape
  // tilt : 1 --> zig-zag line going from left bottom to right top
  // tilt : 2 --> zig-zag line going from right bottom to left top
  // tilt : 3 --> the others shape

  int status = 0;
  if ((ca.tilt == 1 && cb.tilt == 2) || (ca.tilt == 2 && cb.tilt == 1)) {
    status = -1;
  } else if (ca.tilt == 1 && cb.tilt == 1) {
    status = 1;
  } else if (ca.tilt == 2 && cb.tilt == 2) {
    status = 2;
  } else if (ca.tilt == 0 && cb.tilt == 0) {
    if (ca.nPix < 3 && cb.nPix < 3)
      status = 0;
    else if (ca.xiwidth == cb.xiwidth || ca.zetawidth == cb.zetawidth) {
      status = 4;
    } else
      status = -2;
  }

  return status;

  // EXPLANATION
  // status : 0 --> anything passing through all if statement.
  //        : 1 --> matched in tilt == 1
  //        : 2 --> matched in tilt == 2
  //        : 4 --> matched in tilt == 0
  //        : -1 --> mismatched because 2 tilts are different.
  //        : -2 --> mismatched because 2 tilts are different in the meaning of tilt == 0.
}

IntVec FPCCDSiliconTracking_MarlinTrk::getNHitsInSubDet(SimTrackerHitVec simvec) {
  IntVec ivec(3);
  for (int i = 0; i < int(simvec.size()); i++) {
    int detid = getDetectorID(simvec[i]);
    if (detid == lcio::ILDDetID::VXD) {
      ivec[0]++;
    } else if (detid == lcio::ILDDetID::SIT) {
      ivec[1]++;
    } else if (detid == lcio::ILDDetID::FTD) {
      ivec[2]++;
    }
  }
  return ivec;
}

// This is a debug tool for debug by mori.
MCPMap FPCCDSiliconTracking_MarlinTrk::LoadMCPMap() {
  _mcpVXD.clear();
  _mcpSIT.clear();
  _mcpVXDSIT.clear();
  _mcpFTD.clear();
  _mcpFTDSIT.clear();
  _mcpVXDFTD.clear();
  _mcpVXDFTDSIT.clear();

  //_colMCP    = GetCollection(_current_event, _colNameMCParticles);
  _simVXD = GetCollection(_current_event, _colNameVXDSimHit);
  _simSIT = GetCollection(_current_event, _colNameSITSimHit);
  _simFTDpix = GetCollection(_current_event, _colNameFTDpixSimHit);
  _simFTDsp = GetCollection(_current_event, _colNameFTDspSimHit);

  SimTrackerHitVec simVec;
  int nvxd = 0, nsit = 0, nftdpix = 0, nftdsp = 0;
  if (_simVXD != NULL) {
    nvxd = _simVXD->getNumberOfElements();
    for (int i = 0; i < nvxd; i++)
      simVec.push_back(dynamic_cast<SimTrackerHit*>(_simVXD->getElementAt(i)));
  }
  if (_simSIT != NULL) {
    nsit = _simSIT->getNumberOfElements();
    for (int i = 0; i < nsit; i++)
      simVec.push_back(dynamic_cast<SimTrackerHit*>(_simSIT->getElementAt(i)));
  }
  if (_simFTDpix != NULL) {
    nftdpix = _simFTDpix->getNumberOfElements();
    for (int i = 0; i < nftdpix; i++)
      simVec.push_back(dynamic_cast<SimTrackerHit*>(_simFTDpix->getElementAt(i)));
  }
  if (_simFTDsp != NULL) {
    nftdsp = _simFTDsp->getNumberOfElements();
    for (int i = 0; i < nftdsp; i++)
      simVec.push_back(dynamic_cast<SimTrackerHit*>(_simFTDsp->getElementAt(i)));
  }

  MCPMap mymap = _moriUtil->MakeMCPMap(simVec);
  for (MCPMap::iterator it = mymap.begin(); it != mymap.end(); it++) {
    if (it->first == NULL)
      continue;
    IntVec ivec = getNHitsInSubDet(it->second); // ivec[0] : VXD, ivec[1] : SIT, ivec[2] : FTD
    if (ivec[2] == 0) {
      if (ivec[1] == 0) {
        _mcpVXD.insert(*it);
      } else {
        if (ivec[0] == 0) {
          _mcpSIT.insert(*it);
        } else {
          _mcpVXDSIT.insert(*it);
        }
      }
    } else {
      if (ivec[0] == 0) {
        if (ivec[1] == 0) {
          _mcpFTD.insert(*it);
        } else {
          _mcpFTDSIT.insert(*it);
        }
      } else {
        if (ivec[1] == 0) {
          _mcpVXDFTD.insert(*it);
        } else {
          _mcpVXDFTDSIT.insert(*it);
        }
      }
    }
  }

  if (1) {
    std::cout << "nvxd : " << nvxd << std::endl;
    std::cout << "nsit : " << nsit << std::endl;
    std::cout << "nftdpix : " << nftdpix << std::endl;
    std::cout << "nftdsp : " << nftdsp << std::endl;
    std::cout << "nmcp contributing tracker : " << mymap.size() << std::endl;
    std::cout << "nmcp contributing VXD: " << _mcpVXD.size() << std::endl;
    std::cout << "nmcp contributing FTD: " << _mcpFTD.size() << std::endl;
    std::cout << "nmcp contributing SIT: " << _mcpSIT.size() << std::endl;
    std::cout << "nmcp contributing VXD + FTD: " << _mcpVXDFTD.size() << std::endl;
    std::cout << "nmcp contributing VXD + SIT: " << _mcpVXDSIT.size() << std::endl;
    std::cout << "nmcp contributing FTD + SIT: " << _mcpFTDSIT.size() << std::endl;
    std::cout << "nmcp contributing VXD + FTD + SIT: " << _mcpVXDFTDSIT.size() << std::endl;
  }

  return mymap;
}

int FPCCDSiliconTracking_MarlinTrk::getPhiThetaRegion(TrackExtended* trackAR, int layer, int* Boundaries) {
  TrackerHit* currentInnermostHit = trackAR->getTrackerHitExtendedVec().back()->getTrackerHit();
  double d0error = (trackAR->getTrackerHitExtendedVec().size() == 3) ? 1.0 / sqrt(trackAR->getCovMatrix()[9])
                                                                     : sqrt(trackAR->getCovMatrix()[0]);
  double omegaerror = (trackAR->getTrackerHitExtendedVec().size() == 3) ? 1.0 / sqrt(trackAR->getCovMatrix()[0])
                                                                        : sqrt(trackAR->getCovMatrix()[5]);
  double z0error = (trackAR->getTrackerHitExtendedVec().size() == 3) ? 1.0 / sqrt(trackAR->getCovMatrix()[14])
                                                                     : sqrt(trackAR->getCovMatrix()[9]);

  double phi0 = trackAR->getPhi();
  double d0 = trackAR->getD0();
  double z0 = trackAR->getZ0();
  double omega = trackAR->getOmega();
  double tanL = trackAR->getTanLambda();

  // check
  if (_mydebugIntersection) {
    std::cout << "========layer : " << layer << " ====================" << std::endl;
    printf("d0error,omegaerror,z0error : %f, %f, %f \n", d0error, omegaerror, z0error);
    printf("d0,omega,z0,phi0,tanlambda : %f, %f, %f, %f, %f \n", d0, omega, z0, phi0, tanL);
  }

  if (std::isnormal(d0error) == false || (d0error < 1e-6)) {
    if (_mydebugIntersection)
      std::cout << "FATAL ERROR of d0error: nan check MORI " << std::endl;
    if (trackAR->getTrackerHitExtendedVec().size() == 3) {
      if (_mydebugIntersection)
        printf("d0error cov[9] : %f \n", trackAR->getCovMatrix()[9]);
      d0error = -1.0;
    } else {
      if (_mydebugIntersection)
        printf("d0error cov[0] : %f \n", trackAR->getCovMatrix()[0]);
      d0error = -1.0;
    }
  }
  if (std::isnormal(omegaerror) == false || (omegaerror < 1e-6)) {
    if (_mydebugIntersection)
      std::cout << "FATAL ERROR of omegaerror: nan check MORI " << std::endl;
    if (trackAR->getTrackerHitExtendedVec().size() == 3) {
      if (_mydebugIntersection)
        printf("omegaerror cov[0] : %f \n", trackAR->getCovMatrix()[0]);
      omegaerror = -1.0;
    } else {
      if (_mydebugIntersection)
        printf("omegaerror cov[5] : %f \n", trackAR->getCovMatrix()[5]);
      omegaerror = -1.0;
    }
  }
  if (std::isnormal(z0error) == false || (z0error < 1e-6)) {
    if (_mydebugIntersection)
      std::cout << "FATAL ERROR of z0error: nan check MORI " << std::endl;
    if (trackAR->getTrackerHitExtendedVec().size() == 3) {
      if (_mydebugIntersection)
        printf("z0error cov[14] : %f \n", trackAR->getCovMatrix()[14]);
      z0error = -1.0;
    } else {
      if (_mydebugIntersection)
        printf("z0error cov[9] : %f \n", trackAR->getCovMatrix()[9]);
      z0error = -1.0;
    }
  }

  DoubleVec phiSectsTest(5);
  DoubleVec thetaSectsTest(3);

  for (int i = 0; i < 5; i++) {
    std::vector<double> iSec;
    HelixClass_double tmpHelix;
    if (i == 1 && (d0error < 0)) {
      phiSectsTest[1] = phiSectsTest[0];
      phiSectsTest[2] = phiSectsTest[0];
      i = 2;
      continue;
    }
    if (i >= 3 && (omegaerror < 0)) {
      phiSectsTest[3] = phiSectsTest[0];
      phiSectsTest[4] = phiSectsTest[0];
      break;
    }

    if (i == 0)
      tmpHelix.Initialize_Canonical(phi0, d0, z0, omega, tanL, _bField);
    else if (i == 1)
      tmpHelix.Initialize_Canonical(phi0, d0 + d0error * _nSigmaBuild_phi, z0, omega, tanL, _bField);
    else if (i == 2)
      tmpHelix.Initialize_Canonical(phi0, d0 - d0error * _nSigmaBuild_phi, z0, omega, tanL, _bField);
    else if (i == 3)
      tmpHelix.Initialize_Canonical(phi0, d0, z0, omega + omegaerror * _nSigmaBuild_phi, tanL, _bField);
    else if (i == 4)
      tmpHelix.Initialize_Canonical(phi0, d0, z0, omega - omegaerror * _nSigmaBuild_phi, tanL, _bField);
    int error = getIntersectionEasyTest(tmpHelix, currentInnermostHit, layer, iSec);

    if (error == -1) {
      // std::cout << "getIntersectionEasyTest for phi-sector couldn't find intersection. loop : " << i << std::endl;
      return int(trackAR->getTrackerHitExtendedVec().size());
    }
    if (layer < int(_nLayersVTX)) {
      if (std::abs(iSec[2]) > _vxd.geodata[layer].hlength) {
        if (_mydebugIntersection)
          std::cout << "Z length is too long. VXD at 1" << std::endl;
        return int(trackAR->getTrackerHitExtendedVec().size());
      }
    } else if (layer >= int(_nLayersVTX)) {
      if (std::abs(iSec[2]) > _sit.geodata[layer - _nLayersVTX].hlength) {
        if (_mydebugIntersection)
          std::cout << "Z length is too long. SIT at 1" << std::endl;
        return int(trackAR->getTrackerHitExtendedVec().size());
      }
    }

    double tmpphi = atan2(iSec[1], iSec[0]);
    if (tmpphi < 0.0)
      tmpphi += 2.0 * M_PI;
    phiSectsTest[i] = tmpphi / (2.0 * M_PI / _nDivisionsInPhi);
  }

  for (int i = 0; i < 3; i++) {
    std::vector<double> iSec;
    HelixClass_double tmpHelix;
    if (i == 1 && (z0error < 0)) {
      thetaSectsTest[1] = thetaSectsTest[0];
      thetaSectsTest[2] = thetaSectsTest[0];
      break;
    }
    if (i == 0)
      tmpHelix.Initialize_Canonical(phi0, d0, z0, omega, tanL, _bField);
    else if (i == 1) {
      tmpHelix.Initialize_Canonical(phi0, d0, z0 + z0error * _nSigmaBuild_theta, omega, tanL, _bField);
    } else if (i == 2) {
      tmpHelix.Initialize_Canonical(phi0, d0, z0 - z0error * _nSigmaBuild_theta, omega, tanL, _bField);
    }

    int error = getIntersectionEasyTest(tmpHelix, currentInnermostHit, layer, iSec);
    if (error == -1) {
      // std::cout << "getIntersectionEasyTest for phi-theta couldn't find intersection. loop : " << i << std::endl;
      return int(trackAR->getTrackerHitExtendedVec().size());
    }
    if (layer < int(_nLayersVTX)) {
      if (std::abs(iSec[2]) > _vxd.geodata[layer].hlength) {
        if (_mydebugIntersection)
          std::cout << "Z length is too long. VXD at 2" << std::endl;
        return int(trackAR->getTrackerHitExtendedVec().size());
      }
    } else if (layer >= int(_nLayersVTX)) {
      if (std::abs(iSec[2]) > _sit.geodata[layer - _nLayersVTX].hlength) {
        if (_mydebugIntersection)
          std::cout << "Z length is too long. SIT at 2" << std::endl;
        return int(trackAR->getTrackerHitExtendedVec().size());
      }
    }
    double tmpcostheta = iSec[2] / sqrt(iSec[0] * iSec[0] + iSec[1] * iSec[1] + iSec[2] * iSec[2]);
    thetaSectsTest[i] = (tmpcostheta + 1.0) / _dTheta;
  }
  // check
  /*
     std::cout << "phiSectsTest[0] : " << phiSectsTest[0] << std::endl;
     std::cout << "phiSectsTest[1] : " << phiSectsTest[1] << std::endl;
     std::cout << "phiSectsTest[2] : " << phiSectsTest[2] << std::endl;
     std::cout << "phiSectsTest[3] : " << phiSectsTest[3] << std::endl;
     std::cout << "phiSectsTest[4] : " << phiSectsTest[4] << std::endl;
     std::cout << "thetaSectsTest[0] : " << thetaSectsTest[0] << std::endl;
     std::cout << "thetaSectsTest[1] : " << thetaSectsTest[1] << std::endl;
     std::cout << "thetaSectsTest[2] : " << thetaSectsTest[2] << std::endl;
   */

  int iPhiUpTest =
      int(std::max(phiSectsTest[0],
                   std::max(phiSectsTest[1], std::max(phiSectsTest[2], std::max(phiSectsTest[3], phiSectsTest[4]))))) +
      1 + _fudgePhiRange;
  int iPhiLowTest =
      int(std::min(phiSectsTest[0],
                   std::min(phiSectsTest[1], std::min(phiSectsTest[2], std::min(phiSectsTest[3], phiSectsTest[4]))))) -
      1 - _fudgePhiRange;
  int iThetaUpTest =
      int(std::max(thetaSectsTest[0], std::max(thetaSectsTest[1], thetaSectsTest[2]))) + 1 + _fudgeThetaRange;
  int iThetaLowTest =
      int(std::min(thetaSectsTest[0], std::min(thetaSectsTest[1], thetaSectsTest[2]))) - 1 - _fudgeThetaRange;

  if (iPhiLowTest + _nDivisionsInPhi - iPhiUpTest < iPhiUpTest - iPhiLowTest) {
    int swap1 = iPhiUpTest;
    int swap2 = iPhiLowTest + _nDivisionsInPhi;
    iPhiUpTest = std::max(swap1, swap2);
    iPhiLowTest = std::min(swap1, swap2);
  }

  int bugCheck = 0;
  if (iThetaLowTest < 0) {
    iThetaLowTest = 0;
    bugCheck++;
  }
  if (iThetaUpTest >= _nDivisionsInTheta) {
    iThetaUpTest = _nDivisionsInTheta - 1;
    bugCheck++;
  }
  if (bugCheck == 2) {
    std::cout << "FATAL ERROR. Check code!!; however, process isn't stopped. iThetaLow,Up = 0,0 is returned.) "
              << std::endl;
    iThetaLowTest = 0;
    iThetaUpTest = 0;
  }

  // check
  if (_mydebugIntersection) {
    std::cout << "iPhiLowTest, UpTest = " << iPhiLowTest << " , " << iPhiUpTest << std::endl;
    std::cout << "iThetaLowTest, UpTest = " << iThetaLowTest << " , " << iThetaUpTest << std::endl;
  }

  Boundaries[0] = iPhiLowTest;
  Boundaries[1] = iPhiUpTest;
  Boundaries[2] = iThetaLowTest;
  Boundaries[3] = iThetaUpTest;

  return 0;
}

FPCCDSiliconTracking_MarlinTrk::MCPContributions FPCCDSiliconTracking_MarlinTrk::getMCPContribution(IntVec nsub) {
  if (nsub[0] == 0) {
    if (nsub[1] == 0) {
      if (nsub[2] == 0) {
        std::cout << "In MCPContributionChecker: No Tracker Hit!! Fatal Error. exit" << std::endl;
        exit(1);
      } else {
        return contFTD;
      }
    } else {
      if (nsub[2] == 0) {
        return contSIT;
      } else {
        return contFTDSIT;
      }
    }
  } else {
    if (nsub[1] == 0) {
      if (nsub[2] == 0) {
        return contVXD;
      } else {
        return contVXDFTD;
      }
    } else {
      if (nsub[2] == 0) {
        return contVXDSIT;
      } else {
        return contVXDFTDSIT;
      }
    }
  }
}
