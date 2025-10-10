#ifndef FPCCDFULLLDCTRACKING_H
#define FPCCDFULLLDCTRACKING_H 1

#include "ClusterExtended.h"
#include "ClusterShapes.h"
#include "GroupTracks.h"
#include "HelixClass.h"
#include "HelixClass_double.h"
#include "TrackExtended.h"
#include "TrackHitPair.h"
#include "TrackerHitExtended.h"
#include "lcio.h"
#include "marlin/Processor.h"
#include <map>
#include <marlin/Global.h>
#include <set>
#include <string>

#include "MarlinTrk/IMarlinTrack.h"

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitField64.h>

#include "../src/GetPurity.h"
#include "./MCPMap.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TTree.h"
#include <EVENT/LCRelation.h>
#include <EVENT/SimTrackerHit.h>
#include <UTIL/LCRelationNavigator.h>

using namespace lcio;
using namespace marlin;

namespace MarlinTrk {
class HelixFit;
class IMarlinTrkSystem;
} // namespace MarlinTrk

namespace dd4hep {
class Detector;
}

/** === FPCCDFullLDCTracking_MarlinTrk Processor === <br>
 *
 *
 * This processor is based on FullLDCTracking_MarlinTrk Processor (author: Steve Apline). <br>
 * Please use this processor with FPCCDSiliconTracking_MarlinTrk, then tracking performance is improved.  <br>
 * As I said in FPCCDSiliconTracking_MarlinTrk.h, users can use CMOS VXD Simulator with these tracking code in spite of
 the name FPCCD~.  <br>
 *
 * There are only two differences between FPCCDFullLDCTracking_MarlinTrk and FullLDCTracking_MarlinTrk. <br>
 * The one is the requirement of registering silicon tracks in this processor. <br>
 * If useMaxChi2RequirementForSiTrk : false, then the requirement requires <br>
 * a threshold of probability of track (old requirement). <br>
 * If true, it does chi2/ndf (new requirement).  <br>
 * The reason that I added new one is that good reconstructed low Pt tracks  <br>
 * tend to have relatively higher probability than high Pt tracks.  <br>
 * For saving those tracks, it is more effective to use the threshold of chi2/ndf.  <br>
 * The other is the strategy for reducing pair-BG tracks.  <br>
 * If FinalTrackCut_strategyA is true, then at last stage of this processor,  <br>
 * tracks seemed to be pair-BG tracks are discarded by requiring  <br>
 * SIT htis >= 1 || TPC hits >= 1 || abs(costheta) > 0.9  <br>
 * About this track requirement, please see   <br>
 * https://agenda.linearcollider.org/getFile.py/access?contribId=5&resId=0&materialId=slides&confId=6294 <br>
 *
 * @author Tatsuya Mori (Tohoku University)<br>
 *
 *
 * The following sentence is the copy of FullLDCTracking_MarlinTrk.
 *
 *
 *
 * Processor performing track finding procedure in
 * the entire LDC detector by linking track segments
 * found by the SiliconTracking module in the silicon detectors
 * and by the LEPTracking module in TPC.
 * <h4>Input collections and prerequisites</h4>
 * Processor requires collection of digitized vertex, sit, ftd, set, etd & tpc tracker hits
 * and also the collections of tracks found in the silicon detectors
 * and in TPC.
 * <h4>Output</h4>
 * Processor produces an LCIO collection of the Tracks. Each track is characterised by
 * five parameters : Omega (signed curvuture), Tan(lambda) where
 * lambda is the dip angle, Phi (azimuthal angle @ point of closest approach), D0 (signed impact parameter),
 * Z0 (displacement along z axis at the point of closest approach to IP).
 * Covariance matrix for these parameters is also provided.
 * Only lower left corner of the covariance matrix is stored. The sequence of the covariance matrix elements
 * assigned to track is the following: <br>
 * (D0,D0) <br>
 * (Phi,D0), (Phi,Phi) <br>
 * (Omega,D0), (Omega,Phi), (Omega,Omega) <br>
 * (Z0,D0), (Z0,Phi), (Z0,Omega), (Z0,Z0) <br>
 * (TanL,D0), (TanL,Phi), (TanL,Omega), (TanL,Z0), (TanL,TanL) <br>
 * The number of hits in the different subdetectors associated
 * with each track can be accessed via method Track::getSubdetectorHitNumbers().
 * This method returns vector of integers : <br>
 * number of VTX hits used in the track fit is the 1st element in this vector
 * (Track::getSubdetectorHitNumbers()[0]) <br>
 * number of FTD hits used in the track fit is the 2nd element in this vector
 * (Track::getSubdetectorHitNumbers()[1]) <br>
 * number of SIT hits used in the track fit is the 3d element in this vector
 * (Track::getSubdetectorHitNumbers()[2]) <br>
 * number of TPC hits used in the track fit is the 4th element in this vector
 * (Track::getSubdetectorHitNumbers()[3]) <br>
 * number of SET hits used in the track fit is the 5th element in this vector
 * (Track::getSubdetectorHitNumbers()[4]) <br>
 * number of ETD hits used in the track fit is the 6th element in this vector
 * (Track::getSubdetectorHitNumbers()[5]) <br>
 * total number of VTX hits in track is the 7th element in this vector
 * (Track::getSubdetectorHitNumbers()[6]) <br>
 * total number of FTD hits in track is the 8th element in this vector
 * (Track::getSubdetectorHitNumbers()[7]) <br>
 * total number of SIT hits in track is the 9th element in this vector
 * (Track::getSubdetectorHitNumbers()[8]) <br>
 * total number of TPC hits in track is the 10th element in this vector
 * (Track::getSubdetectorHitNumbers()[9]) <br>
 * total number of SET hits in track is the 11th element in this vector
 * (Track::getSubdetectorHitNumbers()[10]) <br>
 * total number of ETD hits in track is the 12th element in this vector
 * (Track::getSubdetectorHitNumbers()[11]) <br>
 * Output track collection has by default a name "LDCTracks".
 * @param VTXHitCollection name of input VTX TrackerHit collection <br>
 * (default parameter value : "VTXTrackerHits") <br>
 * @param FTDPixelHitCollectionName name of input FTD Pixel TrackerHit collection <br>
 * (default parameter value : "FTDPixelTrackerHits") <br>
 * @param FTDSpacePointCollectionName name of input FTD Space Point TrackerHit collection <br>
 * (default parameter value : "FTDSpacePoints") <br>
 * @param SITHitCollection name of input SIT TrackerHit collection <br>
 * (default parameter value : "SITTrackerHits") <br>
 * @param TPCHitCollection name of input TPC TrackerHit collection <br>
 * (default parameter value : "TPCTrackerHits") <br>
 * @param SETHitCollection name of input SET TrackerHit collection <br>
 * (default parameter value : "SETTrackerHits") <br>
 * @param ETDHitCollection name of input ETD TrackerHit collection <br>
 * (default parameter value : "ETDTrackerHits") <br>
 * @param TPCTracks collection name of TPC tracks <br>
 * (default parameter value : "TPCTracks") <br>
 * @param TPCTracksMCPRelColl Name of input TPC track to MC particle relation collection <br>
 * (default parameter value : "TPCTracksMCP") <br>
 * @param SiTracks collection name of Si tracks <br>
 * (default parameter value : "SiTracks") <br>
 * @param SiTracksMCPRelColl Name of input Si track to MC particle relation collection <br>
 * (default parameter value : "SiTracksMCP") <br>
 * @param LDCTrackCollection name of the output LDC track collection <br>
 * (default parameter value : "LDCTracks") <br>
 * @param Chi2FitCut cut on the Chi2/Ndf of the track fit <br>
 * (default parameter value : 100.0) <br>
 * @param Chi2PrefitCut cut on the prefit Chi2 of the track candidate,
 * prefit is done with the simple helix hypothesis <br>
 * (default parameter value : 1e+5) <br>
  * @param AngleCutForMerging  cut on opening angle between
 * particle momentum reconstructed with TPC and momentum reconstructed
 * with the Silicon detectors. If the opening angle is smaller that this cut
 * the track segment in Silicon trackers and in TPC are tested for their
 * compatibility <br>
 * (default parameter value : 0.10) <br>
 * @param OmegaCutForMerging  cut on the relative difference in the track Omega
 * parameter reconstructed with TPC and with Si detectors. If the relative difference is smaller
 * than this cut, the track segments in TPC and Si are tested for their compatibility <br>
 * (default parameter value : 0.25) <br>
 * @param D0CutForMerging Upper cutoff on the difference in D0 [mm] to allow for merging
 * of the Si and TPC segments <br>
 * (default parameter value : 500) <br>
 * @param Z0CutForMerging Upper cutoff on the difference in Z0 [mm] to allow for merging
 * of the Si and TPC segments <br>
 * (default parameter value : 1000) <br>
 * @param Debug flag to allow for printout of debug information,
 * if set to 1 debugging printout is activated
 * (default parameter value : 1) <br>
 * @param ForceSiTPCMerging This flag steers merging of Si and TPC track segments. If ForceMerging=1
 * Si and TPC track segments are forced to be merged if the opening angle between Si track
 * momentum and TPC track momentum
 * is less than AngleCutForForcedMerging (see below) and difference in tracks
 * parameters Omega is less than OmegaCutForForcedMerging (see below) <br>
 * (default parameter value : 0)
 * @param AngleCutForForcedMerging cut on opening angle between Si track momentum and
 * TPC track momentum. Used to steer forced merging of Si and TPC track segments. <br>
 * (default parameter value : 0.05)
 * @param OmegaCutForForcedMerging cut on the difference between Si and TPC tracks parameter
 * Omega. Used to steer forced merging of Si and TPC track segments. Relative
 * errors are compared. <br>
 * (default parameter value : 0.15) <br>
 * @param D0CutForForcedMerging Upper cutoff on the difference in D0 to allow for forced
 * merging of the Si and TPC segments <br>
 * (default parameter value : 50) <br>
 * @param Z0CutForForcedMerging Upper cutoff on the difference in Z0 to allow for forced
 * merging of the Si and TPC segments <br>
 * (default parameter value : 200) <br>
 * @param ForceTPCSegmentsMerging If this flag is set to 1, the code attempts to
 * merge TPC segments from the low pt splitted loopers <br>
 * (default parameter value : 1) <br>
 * @param D0CutToMergeTPCSegments cut on the difference in the track parameter
 * d0 [mm] to allow for merging TPC segments <br>
 * (default parameter value : 100) <br>
 * @param Z0CutToMergeTPCSegments cut on the difference in the track parameter
 * z0 [mm] to allow for merging TPC segments <br>
 * (default parameter value : 5000) <br>
 * @param DeltaPCutToMergeTPCSegments cut on the magnitude [GeV/c] of the vectorial difference
 * of the momentum vectors, associated with TPC segments, for the TPC segment's merging procedure <br>
 * (default parameter value : 0.1) <br>
 * @param PtCutToMergeTPCSegments lower cutoff on Pt of the TPC segments of the looping track for
 * the merging procedure.
 * If transverse momentum of the segments is less than cutoff the segments are allowed to be merged. <br>
 * (default parameter value : 1.2) <br>
 * @param AssignTPCHits If this flag is set to 1, the code attempts to assign left-over
 * TPC hits to the accepted track candidates. No track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param AssignETDHits If this flag is set to 1, the code attempts to assign
 * ETD hits to the accepted track candidates. No track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param AssignVTXHits If this flag is set to 1, the code attempts to assign left-over
 * VTX hits to the accepted track candidates. Track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param AssignFTDHits If this flag is set to 1, the code attempts to assign left-over
 * FTD hits to the accepted track candidates. Track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param AssignSITHits If this flag is set to 1, the code attempts to assign left-over
 * SIT hits to the accepted track candidates. Track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param AssignSETHits If this flag is set to 1, the code attempts to assign
 * SET hits to the accepted track candidates. Track refit is done in case when hit is assigned
 * to the existing track <br>
 * (default parameter value : 1) <br>
 * @param TPCHitToTrackDistance Cut on the distance between left-over TPC hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 15.0) <br>
 * @param VTXHitToTrackDistance Cut on the distance between left-over VTX hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 1.5) <br>
 * @param FTDHitToTrackDistance Cut on the distance between left-over FTD hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 2.0) <br>
 * @param SITHitToTrackDistance Cut on the distance between left-over SIT hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 2.0) <br>
 * @param SETHitToTrackDistance Cut on the distance between SET hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 2.0) <br>
 * @param ETDHitToTrackDistance Cut on the distance between ETD hit and the track helix
 * to allow for assignment of the hit with a given track <br>
 * (default parameter value : 10.0) <br>
 * @param NHitsExtrapolation Number of the last track hits for extrapolating helix
 * to the outer tracking detectors (SET, ETD) <br>
 * (default parameter value : 35) <br>
 * @param CutOnTPCHits minimal number of TPC hits, used in the track fit, which is
 * required for tracks which have no hits from the Si detectors <br>
 * (default parameter value : 35) <br>
 * @param CutOnTrackD0 cut on the d0 parameter of the track. If the d0 parameter is greater that
 * this cut, track is rejected <br>
 * (default parameter value : 500) <br>
 * @param CutOnTrackZ0 cut on the z0 parameter of the track. If the z0 parameter is greater that
 * this cut, track is rejected <br>
 * (default parameter value : 500) <br>
 * @param ForbidOverlapInZTPC If this flag is set to 1 then merging of the TPC semiloops is
 * forbiden for segment overlapping in z <br>
 * (default parameter value : 0) <br>
 * @param ForbidOverlapInZComb If this flag is set to 1 then merging of left-over TPC semiloop and
 * combined Si-TPC track is their segments overlap in z <br>
 * (default parameter value : 0) <br>
 * @param cosThetaCutHighPtMerge cut on cos theta between the two momentum vectors
 * when considering merger of high Pt tracks <br>
 * (default is 0.99) <br>
 * @param cosThetaCutSoftHighPtMerge cut on the cos theta between the two momentum vectors
 * when considering merger of high Pt tracks for softer dp/p cut <br>
 * (default is 0.998) <br>
 * @param momDiffCutHighPtMerge cut on dp/p
 * when considering merger of high Pt tracks <br>
 * (default is 0.01 [1/GeV]) <br>
 * @param momDiffCutSoftHighPtMerge softer cut on dp/p
 * when considering merger of high Pt tracks <br>
 * (default is 0.25 [1/GeV]) <br>
 * @param hitDistanceCutHighPtMerge cut on 3D distance between hit
 * and helix extrapolation when considering merger of high Pt tracks <br>
 * (default is 25.0 [mm]) <br>
 * @param maxHitDistanceCutHighPtMerge cut for max 3D distance between any hit
 * and helix extrapolation when considering merger of high Pt tracks <br>
 * (default is 50.0 [mm]) <br>
 * @param maxFractionOfOutliersCutHighPtMerge cut on maximum fraction of outliers
 * when considering merger of high Pt tracks <br>
 * (default is 0.95 ) <br>


 * @author A. Raspereza (MPI Munich)<br>
 */

class FPCCDFullLDCTracking_MarlinTrk : public Processor {
public:
  virtual Processor* newProcessor() { return new FPCCDFullLDCTracking_MarlinTrk; }
  FPCCDFullLDCTracking_MarlinTrk();
  FPCCDFullLDCTracking_MarlinTrk(const FPCCDFullLDCTracking_MarlinTrk&) = delete;
  FPCCDFullLDCTracking_MarlinTrk& operator=(const FPCCDFullLDCTracking_MarlinTrk&) = delete;
  virtual void init();
  virtual void processRunHeader(LCRunHeader* run);
  virtual void processEvent(LCEvent* evt);
  virtual void check(LCEvent* evt);
  virtual void end();

protected:
  void prepareVectors(LCEvent* evt);
  void CleanUp();
  void MergeTPCandSiTracks();
  void MergeTPCandSiTracksII();

  TrackExtended* CombineTracks(TrackExtended* tpcTrk, TrackExtended* siTrk, float maxAllowedOutliers,
                               bool testCombinationOnly);

  //  TrackExtended * TrialCombineTracks(TrackExtended * tpcTrk, TrackExtended * siTrk);

  void Sorting(TrackExtendedVec& trackVec);
  void SelectCombinedTracks();
  void AddNotCombinedTracks();
  void CheckTracks();
  void AddNotAssignedHits();
  void RemoveSplitTracks();
  void AddTrackColToEvt(LCEvent* evt, TrackExtendedVec& trkVec, std::string TrkColName);
  float CompareTrk(TrackExtended* first, TrackExtended* second, float d0Cut, float z0Cut, int iopt);

  float CompareTrkII(TrackExtended* first, TrackExtended* second, float d0Cut, float z0Cut, int iopt, float& Angle);
  float CompareTrkIII(TrackExtended* first, TrackExtended* second, float d0Cut, float z0Cut, int iopt, float& Angle);

  void SortingTrackHitPairs(TrackHitPairVec& trackHitPairVec);

  void AssignSiHitsToTracks(TrackerHitExtendedVec hitVec, float dcut);

  void AssignTPCHitsToTracks(TrackerHitExtendedVec hitVec, float dcut);

  void AssignOuterHitsToTracks(TrackerHitExtendedVec hitVec, float dcut, int refit);

  void CreateExtrapolations();

  void CleanUpExtrapolations();

  HelixClass* GetExtrapolationHelix(TrackExtended* track);

  void PrintOutMerging(TrackExtended* firstTrackExt, TrackExtended* secondTrackExt, int iopt);

  void GeneralSorting(int* index, float* val, int direct, int nVal);

  int SegmentRadialOverlap(TrackExtended* pTracki, TrackExtended* pTrackj);
  bool VetoMerge(TrackExtended* firstTrackExt, TrackExtended* secondTrackExt);

  int _nRun{};
  int _nEvt{};

  MarlinTrk::HelixFit* _fastfitter{nullptr};

  /** pointer to the IMarlinTrkSystem instance
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem{nullptr};

  bool _MSOn{}, _ElossOn{}, _SmoothOn{};

  std::string _TPCTrackCollection{};
  std::string _SiTrackCollection{};
  std::string _TPCTrackMCPCollName{};
  std::string _SiTrackMCPCollName{};

  std::string _VTXTrackerHitCollection{};
  std::string _SITTrackerHitCollection{};
  std::string _SETTrackerHitCollection{};
  std::string _FTDPixelHitCollection{};
  std::string _FTDSpacePointCollection{};
  std::string _TPCTrackerHitCollection{};
  std::string _ETDTrackerHitCollection{};

  std::string _LDCTrackCollection{};

  TrackExtendedVec _allSiTracks{};
  TrackExtendedVec _allTPCTracks{};
  TrackExtendedVec _allCombinedTracks{};
  TrackExtendedVec _allNonCombinedTPCTracks{};
  TrackExtendedVec _allNonCombinedSiTracks{};
  TrackExtendedVec _trkImplVec{};
  TrackerHitExtendedVec _allTPCHits{};
  TrackerHitExtendedVec _allVTXHits{};
  TrackerHitExtendedVec _allFTDHits{};
  TrackerHitExtendedVec _allSITHits{};
  TrackerHitExtendedVec _allSETHits{};
  TrackerHitExtendedVec _allETDHits{};

  float PI{}, PIOVER2{}, TWOPI{};

  float _bField{};
  float _chi2PrefitCut{};
  float _chi2FitCut{};

  int _debug{};

  float _dPCutForMerging{};
  float _d0CutForMerging{};
  float _z0CutForMerging{};
  float _dOmegaForMerging{};
  float _angleForMerging{};

  int _forceMerging{};
  float _dPCutForForcedMerging{};
  float _d0CutForForcedMerging{};
  float _z0CutForForcedMerging{};
  float _dOmegaForForcedMerging{};
  float _angleForForcedMerging{};

  int _mergeTPCSegments{};
  float _dPCutToMergeTPC{};
  float _PtCutToMergeTPC{};
  float _d0CutToMergeTPC{};
  float _z0CutToMergeTPC{};

  float _cosThetaCutHighPtMerge{};
  float _cosThetaCutSoftHighPtMerge{};
  float _momDiffCutHighPtMerge{};
  float _momDiffCutSoftHighPtMerge{};
  float _hitDistanceCutHighPtMerge{};
  float _maxHitDistanceCutHighPtMerge{};
  float _maxFractionOfOutliersCutHighPtMerge{};

  float _vetoMergeMomentumCut{};

  float _initialTrackError_d0{};
  float _initialTrackError_phi0{};
  float _initialTrackError_omega{};
  float _initialTrackError_z0{};
  float _initialTrackError_tanL{};

  double _maxChi2PerHit{};
  double _minChi2ProbForSiliconTracks{};
  double _maxChi2ForSiliconTracks{};
  bool _useMaxChi2ReqForSiTrk{};
  float _maxAllowedPercentageOfOutliersForTrackCombination{};
  int _maxAllowedSiHitRejectionsForTrackCombination{};

  bool _runMarlinTrkDiagnostics{};
  std::string _MarlinTrkDiagnosticsName{};

  int _nHitsExtrapolation{};

  int _cutOnTPCHits{};
  int _cutOnSiHits{};

  int _assignVTXHits{}, _assignFTDHits{}, _assignSITHits{}, _assignTPCHits{};

  int _assignSETHits{}, _assignETDHits{};

  float _distCutForVTXHits{}, _distCutForFTDHits{}, _distCutForSITHits{}, _distCutForTPCHits{};

  float _distCutForSETHits{}, _distCutForETDHits{};

  float _d0TrkCut{}, _z0TrkCut{};

  int _forbidOverlapInZTPC{}, _forbidOverlapInZComb{};

  LCEvent* _evt{nullptr};

  std::map<TrackExtended*, HelixClass*> _trackExtrapolatedHelix{};
  std::set<TrackExtended*> _candidateCombinedTracks{};

  UTIL::BitField64* _encoder{nullptr};
  int getDetectorID(TrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::subdet()];
  }
  int getSideID(TrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::side()];
  };
  int getLayerID(TrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::layer()];
  };
  int getModuleID(TrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::module()];
  };
  int getSensorID(TrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::sensor()];
  };

  int getDetectorID(SimTrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::subdet()];
  }
  int getSideID(SimTrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::side()];
  };
  int getLayerID(SimTrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::layer()];
  };
  int getModuleID(SimTrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::module()];
  };
  int getSensorID(SimTrackerHit* hit) {
    _encoder->setValue(hit->getCellID0());
    return (*_encoder)[lcio::LCTrackerCellID::sensor()];
  };

  void setupGeom(const dd4hep::Detector& theDetector);

  double _tpc_inner_r{};
  double _tpc_outer_r{};
  double _tpc_pad_height{};
  int _tpc_nrows{};

  //  struct VXD_Layer {
  //    int nLadders;
  //    double phi0;
  //    double dphi;
  //    double senRMin;
  //    double supRMin;
  //    double length;
  //    double width;
  //    double offset;
  //    double senThickness;
  //    double supThickness;
  //  };
  //  std::vector<VXD_Layer> _VXDgeo;

  unsigned int _nLayersVTX{};

  //  struct SIT_Layer {
  //    int nLadders;
  //    double phi0;
  //    double dphi;
  //    double senRMin;
  //    double supRMin;
  //    double length;
  //    double width;
  //    double offset;
  //    double senThickness;
  //    double supThickness;
  //  };
  //  std::vector<SIT_Layer> _SITgeo;

  unsigned int _nLayersSIT{};

  unsigned int _nLayersSET{};

  //  struct FTD_Disk {
  //    int nPetals;
  //    double phi0;
  //    double dphi;
  //
  //    double alpha;
  //    double rInner;
  //    double height;
  //    double innerBaseLength;
  //    double outerBaseLength;
  //    double senThickness;
  //    double supThickness;
  //
  //    double senZPos_even_petal1;
  //    double senZPos_even_petal2;
  //    double senZPos_even_petal3;
  //    double senZPos_even_petal4;
  //
  //    double supZPos_even_petal1;
  //    double supZPos_even_petal2;
  //    double supZPos_even_petal3;
  //    double supZPos_even_petal4;
  //
  //    double senZPos_odd_petal1;
  //    double senZPos_odd_petal2;
  //    double senZPos_odd_petal3;
  //    double senZPos_odd_petal4;
  //
  //    double supZPos_odd_petal1;
  //    double supZPos_odd_petal2;
  //    double supZPos_odd_petal3;
  //    double supZPos_odd_petal4;
  //
  //
  //
  //  };
  //
  //  std::vector<FTD_Disk> _FTDgeo;
  std::vector<float> _zLayerFTD{};

  unsigned int _nLayersFTD{};
  int _nPhiFTD{};
  bool _petalBasedFTDWithOverlaps{};

  /////////////Addition by Mori//////////////////

  GetPurityUtil* _purityUtil{nullptr};
  moriUTIL* _moriUtil{nullptr};
  MCPMap _mcpMapSi{};
  MCPMap _mcpMapFull{};
  IntVec getNHitsInSubDet(SimTrackerHitVec simvec);
  std::vector<LCRelationNavigator*> _naviVecSi{};
  std::vector<LCRelationNavigator*> _naviVecFull{};
  std::string _colNameVXDTrackerHitRelations{};
  std::string _colNameSITSpacePointRelations{};
  std::string _colNameFTDSpacePointRelations{};
  std::string _colNameFTDPixelTrackerHitRelations{};
  std::string _colNameTPCTrackerHitRelations{};
  std::string _colNameSETSpacePointRelations{};
  LCRelationNavigator* _navVXD{nullptr};
  LCRelationNavigator* _navSIT{nullptr};
  LCRelationNavigator* _navFTDsp{nullptr};
  LCRelationNavigator* _navFTDpix{nullptr};
  LCRelationNavigator* _navTPC{nullptr};
  LCRelationNavigator* _navSET{nullptr};
  std::string _colNameVXDSimHit{};
  std::string _colNameSITSimHit{};
  std::string _colNameFTDspSimHit{};
  std::string _colNameFTDpixSimHit{};
  std::string _colNameTPCSimHit{};
  std::string _colNameSETSimHit{};
  LCCollection* _simVXD{nullptr};
  LCCollection* _simSIT{nullptr};
  LCCollection* _simFTDsp{nullptr};
  LCCollection* _simFTDpix{nullptr};
  LCCollection* _simTPC{nullptr};
  LCCollection* _simSET{nullptr};

  bool _mydebug{};
  bool _mydebugPrintMCP{};

  bool _FinalTrackCut_strategyA{};

  /** helper function to get collection using try catch block */
  LCCollection* GetCollection(LCEvent* evt, std::string colName);

  /** helper function to get relations using try catch block */
  LCRelationNavigator* GetRelations(LCEvent* evt, std::string RelName);

  MCPMap LoadMCPMap(int mode);
  std::map<MCParticle*, SimTrackerHitVec> _mcpVXD{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpVXDFTD{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpVXDSIT{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpVXDFTDSIT{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpFTD{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpFTDSIT{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpSIT{};

  enum MCPContributions {
    contVXD,
    contFTD,
    contSIT,
    contVXDFTD,
    contVXDSIT,
    contFTDSIT,
    contVXDFTDSIT,
    contBG,
    contSize
  };
};

#endif
