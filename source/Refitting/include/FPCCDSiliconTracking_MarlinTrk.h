#ifndef FPCCDSILICONTRACKING_MarlinTrk_H
#define FPCCDSILICONTRACKING_MarlinTrk_H 1

#include "ClusterExtended.h"
#include "HelixClass_double.h"
#include "TrackExtended.h"
#include "TrackerHitExtended.h"
#include "lcio.h"
#include "marlin/Processor.h"
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackImpl.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <marlin/Global.h>
#include <string>
#include <vector>

#include "MarlinTrk/IMarlinTrack.h"

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitField64.h>

#include "../src/GetPurity.h" //by Mori
#include "./MCPMap.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TVector3.h"

namespace DiagnosticsHistograms {
class Histograms;
}

namespace DiagnosticsHistograms2D {
class Histograms;
}

using namespace lcio;
using namespace marlin;

namespace MarlinTrk {
class HelixFit;
class IMarlinTrkSystem;
} // namespace MarlinTrk

namespace UTIL {
class LCRelationNavigator;
}

namespace dd4hep {
class Detector;
}

/** === FPCCDSiliconTracking_MarlinTrk Processor === <br>
 * This processor is based on SiliconTracking_MarlinTrk Processor (author: Steve Apline). <br>
 * The major differences between FPCCDSiliconTracking_MarlinTrk and SiliconTracking_MarlinTrk are <br>
 * TestTriplet method and BuildTrack method, which is to say, making track seed process and extrapolation process. <br>
 * In more detail, please see my talk slide in LCWS 13. <br>
 * https://agenda.linearcollider.org/getFile.py/access?contribId=313&sessionId=37&resId=0&materialId=slides&confId=6000
 * <br> This slide explains the major differences and shows improvement in tracking efficiency by using
 * FPCCDSiliconTracking_MarlinTrk and <br> FPCCDFullLDCTracking_MarlinTrk. <br> In more detail than this slide, now I am
 * preparing to submit proceedings of LCWS 13. In the future, I will write the address here. <br>
 *
 * This processor is named FPCCD~, but if you don't use FPCCD VXD Simulator but CMOS with VXDPlanarDigiProcessor or
 * something like that, <br> abailable, and improves tracking efficiency and flavor tagging and reduces CPU time down to
 * 1/10 if you include pair-BG into your analysis. <br> The performance evaluation study of this processor in the case
 * including pair-BG is shown in <br>
 * https://agenda.linearcollider.org/getFile.py/access?contribId=5&resId=0&materialId=slides&confId=6294 <br>
 *
 * @author Tatsuya Mori (Tohoku University)<br>
 *
 *
 * The following sentence is the copy of SiliconTracking_MarlinTrk.
 *
 *
 * Processor performing stand-alone pattern recognition
 * in the vertex detector (VXD), forward tracking disks and SIT. <br>
 * The procedure consists of three steps : <br>
 * 1) Tracking in VXD and SIT ; <br>
 * 2) Tracking in FTD ; <br>
 * 3) Merging compatible track segments reconstructed in VXD and FTD <br>
 * STEP 1 : TRACKING IN VXD and SIT <br>
 * Algorithm starts with finding of hit triplets satisfying helix hypothesis <br>
 * in three different layers. Two layers of SIT are effectively considered as outermost <br>
 * layers of the vertex detector. To accelerate procedure, the 4-pi solid angle
 * is divided in NDivisionsInTheta and NDivisionsInPhi sectors in cosQ and Phi,
 * respectively. Triplets are looked for in 2x2 window of adjacent sectors.
 * Once triplet is found, attempt is made to associate additional hits to
 * track. Combinatin of hits is accepted for further analysis if the Chi2
 * of the fit is less than certain predefined threshold. All accepted
 * combinations are sorted in ascending order of their Chi2. First track candidate
 * in the sorted array is automatically accepted. The hits belonging to this track are
 * marked as used, and track candidates sharing these hits are discarded.
 * The procedure proceeds with increasing index of track candidate in the sorted
 * array until all track candidate have been output or discarded. <br>
 * STEP 2 : TRACKING IN FTD <br>
 * In the next step tracking in FTD is performed. The strategy of tracking in the FTD
 * is the same as used for tracking in the VXD+SIT. <br>
 * STEP 3 : MERGING TRACK SEGMENTS FOUND IN FTD AND VXD+SIT <br>
 * In the last step, track segments reconstructed in the FTD and VXD+SIT, belonging to the
 * same track  are identified and merged into one track. All possible
 * pairings are tested for their compatibility.
 * The number of pairings considered is Ntrk_VXD_SIT*Ntrk_FTD, where Ntrk_VXD_SIT is the number of
 * track segments reconstructed in the first step in VXD+SIT (segments containing solely VXD and SIT hits) and
 * Ntrk_FTD is the number of track segments reconstructed in the second step
 * (segments containing solely FTD hits).
 * Pair of segments is accepted for further examination if the angle between track segments and
 * than certain specified threshold.
 * Pairing satisfying this condition is subjected for
 * addtitional test. The fit is performed on unified array of hits belonging to both segments.
 * If the chi2 of the fit does not exceed predefined cut value two segments are unified into
 * one track.
 * <h4>Input collections and prerequisites</h4>
 * Processor requires collection of digitized vertex, sit and ftd tracker hits. <br>
 * If such a collections with the user specified names do not exist
 * processor takes no action. <br>
 * <h4>Output</h4>
 * Processor produces an LCIO collection of the Tracks. Each track is characterised by
 * five parameters : Omega (signed curvuture), Tan(lambda) where
 * lambda is the dip angle, Phi (azimuthal angle @ point of closest approach), D0 (signed impact parameter),
 * Z0 (displacement along z axis at the point of closest approach to IP). Covariance matrix for these parameters is also
 * provided. Only lower left corner of the covariance matrix is stored. The sequence of the covariance matrix elements
 * assigned to track is the following: <br>
 * (Omega,Omega) <br>
 * (Omega,TanLambda), (TanLambda,TanLambda) <br>
 * (Omega,Phi), (TanLamda,Phi), (Phi,Phi) <br>
 * (Omega,D0), (TanLambda,D0), (Phi,D0), (D0,D0) <br>
 * (Omega,Z0), (TanLambda,Z0), (Phi,Z0), (D0,Z0), (Z0,Z0) <br>
 * The number of hits in the different subdetectors associated
 * with each track can be accessed via method Track::getSubdetectorHitNumbers().
 * This method returns vector of integers : <br>
 * number of VXD hits in track is the first element in this vector
 * (Track::getSubdetectorHitNumbers()[0]) <br>
 * number of FTD hits in track is the second element in this vector
 * (Track::getSubdetectorHitNumbers()[1]) <br>
 * number of SIT hits in track is the third element in this vector
 * (Track::getSubdetectorHitNumbers()[2]) <br>
 * Output track collection has a name "SiTracks". <br>
 * @param VXDHitCollectionName name of input VXD TrackerHit collection <br>
 * (default parameter value : "VXDTrackerHits") <br>
 * @param FTDHitCollectionName name of input FTD TrackerHit collection <br>
 * (default parameter value : "FTDTrackerHits") <br>
 * @param SITHitCollectionName name of input SIT TrackerHit collection <br>
 * (default parameter value : "SITTrackerHits") <br>
 * @param SiTrackCollectionName name of the output Silicon track collection <br>
 * (default parameter value : "SiTracks") <br>
 * @param LayerCombinations combinations of layers used to search for hit triplets in VXD+SIT <br>
 * (default parameters : 6 4 3  6 4 2  6 3 2  5 4 3  5 4 2  5 3 2  4 3 2  4 3 1  4 2 1  3 2 1) <br>
 * Note that in the VXD+SIT system the first and the second layers of SIT have indicies nLayerVXD and nLayerVXD+1.
 * Combination given above means that triplets are looked first in layers 6 4 3, and then
 * in 6 4 2;  5 4 3;  6 3 2 etc. NOTE THAT LAYER INDEXING STARTS FROM 0.
 * LAYER 0 is the innermost layer  <br>
 * @param LayerCombinationsFTD combinations of layers used to search for hit triplets in FTD <br>
 * (default parameters 6 5 4  5 4 3  5 4 2  5 4 1  5 3 2  5 3 1  5 2 1  4 3 2  4 3 1
 *  4 3 0  4 2 1  4 2 0  4 1 0  3 2 1  3 2 0  3 1 0  2 1 0).
 * NOTE THAT TRACKS IN FTD ARE SEARCHED ONLY IN ONE HEMISPHERE. TRACK IS NOT
 * ALLOWED TO HAVE HITS BOTH IN BACKWARD AND FORWARD PARTS OF FTD SIMULTANEOUSLY.
 * @param NDivisionsInPhi Number of divisions in Phi for tracking in VXD+SIT <br>
 * (default value is 40) <br>
 * @param NDivisionsInTheta Number of divisions in cosQ for tracking in VXD+SIT <br>
 * (default value is 40) <br>
 * @param NDivisionsInPhiFTD Number of divisions in Phi for tracking in FTD <br>
 * (default value is 3) <br>
 * @param Chi2WRphiTriplet weight on chi2 in R-Phi plane for track with 3 hits <br>
 * (default value is 1) <br>
 * @param Chi2WZTriplet weight on chi2 in S-Z plane for track with 3 hits <br>
 * (default value is 0.5) <br>
 * @param Chi2WRphiQuartet weight on chi2 in R-Phi plane to accept track with 4 hits <br>
 * (default value is 1) <br>
 * @param Chi2WZQuartet weight on chi2 in S-Z plane for track with 4 hits <br>
 * (default value is 0.5) <br>
 * @param Chi2WRphiSeptet weight on chi2 in R-Phi plane for track with 5 and more hits <br>
 * (default value is 1) <br>
 * @param Chi2WZSeptet Cut on chi2 in S-Z plane for track with 5 and more hits <br>
 * (default value is 0.5) <br>
 * @param Chi2FitCut Cut on chi2/ndf to accept track candidate <br>
 * (default value is 100.) <br>
 * @param AngleCutForMerging cut on the angle between two track segments.
 * If the angle is greater than this cut, segments are not allowed to be merged. <br>
 * (default value is 0.1) <br>
 * @param MinDistCutAttach cut on the distance (in mm) from hit to the helix. This parameter is used
 * to decide whether hit can be attached to the track. If the distance is less than
 * cut value. The track is refitted with a given hit being added to the list of hits already
 * assigned for the track. Additional hit is assigned if chi2 of the new fit has good chi2. <br>
 * (default value is 2 ) <br>
 * @param MinLayerToAttach the minimal layer index to attach VXD hits to the found hit triplets <br>
 * (default value is -1) <br>
 * @param CutOnZ0 cut on Z0 parameter of track (in mm). If abs(Z0) is greater than the cut value, track is
 * discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 100) <br>
 * @param CutOnD0 cut on D0 parameter of track (in mm). If abs(D0) is greater than the cut value, track is
 * discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 100) <br>
 * @param CutOnPt cut on Pt (GeV/c). If Pt is less than this cut, track is discarded (used to suppress fake
 * track rate in the presence of beam induced background hits) <br>
 * (default value is 0.1) <br>
 * @param MinimalHits minimal number of hits in track required <br>
 * (default value is 3) <br>
 * @param NHitsChi2 Maximal number of hits for which a track with n hits is aways better than one with n-1 hits.
 * For tracks with equal or more than NHitsChi2 the track  with the lower \f$\chi^2\f$ is better.
 * (default value is 5) <br>
 * @param FastAttachment if this flag is set to 1, less accurate but fast procedure to merge additional hits to tracks
 * is used <br> if set to 0, a more accurate, but slower procedure is invoked <br> (default value is 0) <br>
 * @param UseSIT When this flag is set to 1, SIT is included in pattern recognition. When this flag is set
 * to 0, SIT is excluded from the procedure of pattern recognition <br>
 * (default value is 1) <br>
 * <br>
 * @author A. Raspereza (MPI Munich)<br>
 */
class FPCCDSiliconTracking_MarlinTrk : public Processor {
public:
<<<<<<< HEAD
  
  virtual Processor*  newProcessor() { return new FPCCDSiliconTracking_MarlinTrk ; }
  
  
  FPCCDSiliconTracking_MarlinTrk() ;
  FPCCDSiliconTracking_MarlinTrk(const FPCCDSiliconTracking_MarlinTrk&) = delete ;
  FPCCDSiliconTracking_MarlinTrk& operator=(const FPCCDSiliconTracking_MarlinTrk&) = delete ;
  
  /**  
=======
  virtual Processor* newProcessor() { return new FPCCDSiliconTracking_MarlinTrk; }

  FPCCDSiliconTracking_MarlinTrk();
  FPCCDSiliconTracking_MarlinTrk(const FPCCDSiliconTracking_MarlinTrk&) = delete;
  FPCCDSiliconTracking_MarlinTrk& operator=(const FPCCDSiliconTracking_MarlinTrk&) = delete;

  /**
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
   * Initialization
   */
  virtual void init();

  /** Run header processor.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Event processor.
   */
  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  /** Called after data processing for clean up.
   */
<<<<<<< HEAD
  virtual void end() ;
  
  
protected:
  
  int _nRun {};
  int _nEvt {};
  EVENT::LCEvent* _current_event{nullptr};
  
  int _nDivisionsInPhi{};
  int _nDivisionsInTheta{};
  int _nLayers{};
  
  MarlinTrk::HelixFit* _fastfitter{nullptr};
  
  /** pointer to the IMarlinTrkSystem instance 
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem {nullptr};
  bool _runMarlinTrkDiagnostics{};
  std::string _MarlinTrkDiagnosticsName{};
  
  bool _MSOn{}, _ElossOn{}, _SmoothOn{};
  
=======
  virtual void end();

protected:
  int _nRun{};
  int _nEvt{};
  EVENT::LCEvent* _current_event{nullptr};

  int _nDivisionsInPhi{};
  int _nDivisionsInTheta{};
  int _nLayers{};

  MarlinTrk::HelixFit* _fastfitter{nullptr};

  /** pointer to the IMarlinTrkSystem instance
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem{nullptr};
  bool _runMarlinTrkDiagnostics{};
  std::string _MarlinTrkDiagnosticsName{};

  bool _MSOn{}, _ElossOn{}, _SmoothOn{};

>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  float _initialTrackError_d0{};
  float _initialTrackError_phi0{};
  float _initialTrackError_omega{};
  float _initialTrackError_z0{};
  float _initialTrackError_tanL{};
<<<<<<< HEAD
  
  double _maxChi2PerHit{};
  double _maxChi2PerHit2nd{};
  
  bool  _UseEventDisplay{};
  std::vector<int> _colours{};
  
  void drawEvent();
  
  
  // histogram member variables
  
  bool  _createDiagnosticsHistograms{};
  DiagnosticsHistograms::Histograms* _histos {nullptr};

  
  int _ntriplets{}, _ntriplets_good{}, _ntriplets_2MCP{}, _ntriplets_3MCP{}, _ntriplets_1MCP_Bad{}, _ntriplets_bad{};
  
  
=======

  double _maxChi2PerHit{};
  double _maxChi2PerHit2nd{};

  bool _UseEventDisplay{};
  std::vector<int> _colours{};

  void drawEvent();

  // histogram member variables

  bool _createDiagnosticsHistograms{};
  DiagnosticsHistograms::Histograms* _histos{nullptr};

  int _ntriplets{}, _ntriplets_good{}, _ntriplets_2MCP{}, _ntriplets_3MCP{}, _ntriplets_1MCP_Bad{}, _ntriplets_bad{};

>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  /** helper function to get collection using try catch block */
  LCCollection* GetCollection(LCEvent* evt, std::string colName);

  /** helper function to get relations using try catch block */
  LCRelationNavigator* GetRelations(LCEvent* evt, std::string RelName);

  /** input MCParticle collection and threshold used for Drawing */
<<<<<<< HEAD
  std::string  _colNameMCParticles{};
  
  
=======
  std::string _colNameMCParticles{};

>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  /// Compare tracks according to their chi2/ndf
  struct compare_TrackExtended {
    // n.b.: a and b should be TrackExtended const *, but the getters are not const :-(
    bool operator()(TrackExtended* a, TrackExtended* b) const {
      if (a == b)
        return false;
      return (a->getChi2() / a->getNDF() < b->getChi2() / b->getNDF());
    }
  };
<<<<<<< HEAD
  
  
=======

>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  std::string _VTXHitCollection{};
  std::string _FTDPixelHitCollection{};
  std::string _FTDSpacePointCollection{};
  std::string _SITHitCollection{};
  std::string _siTrkCollection{};
<<<<<<< HEAD
  
  std::vector< LCCollection* > _colTrackerHits{};
  std::map< LCCollection*, std::string > _colNamesTrackerHits{};
  
  std::vector<TrackerHitExtendedVec> _sectors{};
  std::vector<TrackerHitExtendedVec> _sectorsFTD{};
  
=======

  std::vector<LCCollection*> _colTrackerHits{};
  std::map<LCCollection*, std::string> _colNamesTrackerHits{};

  std::vector<TrackerHitExtendedVec> _sectors{};
  std::vector<TrackerHitExtendedVec> _sectorsFTD{};

>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  /**
   * A helper class to allow good code readability by accessing tracks with N hits.
   * As the smalest valid track contains three hits, but the first index in a vector is 0,
   * this class hides the index-3 calculation. As the vector access is inline there should be
   * no performance penalty.
   */
  class TracksWithNHitsContainer {
  public:
    /// Empty all the vectors and delete the tracks contained in it.
    void clear();

    /// Set the size to allow a maximum of maxHit hits.
    inline void resize(size_t maxHits) {
      _tracksNHits.resize(maxHits - 2);
      _maxIndex = (maxHits - 3);
    }

    // Sort all track vectors according to chi2/nDof
    //      void sort();

    /// Returns the  TrackExtendedVec for track with n hits.
    /// In case n is larger than the maximal number the vector with the largest n ist returned.
    /// \attention The smallest valid number is three! For
    /// performance reasons there is no safety check!
    inline TrackExtendedVec& getTracksWithNHitsVec(size_t nHits) {
      // return _tracksNHits[ std::min(nHits-3, _maxIndex) ];
      //  for debugging: with boundary check
      return _tracksNHits.at(std::min(nHits - 3, _maxIndex));
    }

  protected:
<<<<<<< HEAD
    std::vector< TrackExtendedVec > _tracksNHits{};
    size_t _maxIndex{}; /// local cache variable to avoid calculation overhead
  };
  
  TracksWithNHitsContainer _tracksWithNHitsContainer{};
  
  int InitialiseVTX(LCEvent * evt);
  int InitialiseFTD(LCEvent * evt);
=======
    std::vector<TrackExtendedVec> _tracksNHits{};
    size_t _maxIndex{}; /// local cache variable to avoid calculation overhead
  };

  TracksWithNHitsContainer _tracksWithNHitsContainer{};

  int InitialiseVTX(LCEvent* evt);
  int InitialiseFTD(LCEvent* evt);
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  void ProcessOneSector(int iSectorPhi, int iSectorTheta);
  void ProcessOneSectorVer2(int iSectorPhi, int iSectorTheta);
  void CleanUp();
  TrackExtended* TestTriplet(TrackerHitExtended* outerHit, TrackerHitExtended* middleHit, TrackerHitExtended* innerHit,
                             HelixClass_double& helix, int omegamode);

<<<<<<< HEAD
  int _useBuildTrackForHighPt{};
  double _cosThetaRangeForBuildTrackForHighPt{};
  double _phiRangeForBuildTrackForHighPt{};
  void getPhiThetaRegionForHighPt(int* boundaries,TrackExtended* trackAR);
=======
  int BuildTrack_KalFit(TrackerHitExtended* outerHit, TrackerHitExtended* middleHit, TrackerHitExtended* innerHit,
                        HelixClass_double& helix, int innerlayer, TrackExtended* trackAR);
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90

  int _useBuildTrackForHighPt{};
  double _cosThetaRangeForBuildTrackForHighPt{};
  double _phiRangeForBuildTrackForHighPt{};
  void getPhiThetaRegionForHighPt(int* boundaries, TrackExtended* trackAR);

  void Sorting(TrackExtendedVec& trackVec);
  void CreateTrack(TrackExtended* trackAR);
  void AttachRemainingVTXHitsSlow();
  void AttachRemainingFTDHitsSlow();
  void AttachRemainingVTXHitsFast();
  void AttachRemainingFTDHitsFast();
  void AttachRemainingVTXHitsVeryFast();
  void TrackingInFTD();
  int BuildTrackFTD(TrackExtended* trackAR, int* nLR, int iS);
  int AttachHitToTrack(TrackExtended* trackAR, TrackerHitExtended* hit, int iopt);
  int AttachHitToTrack_KalFit(TrackExtended* trackAR, TrackerHitExtended* hit);

  void FinalRefit(LCCollectionVec* trk_col, LCCollectionVec* rel_col);
<<<<<<< HEAD
  
=======

>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  float _bField{};
  float _chi2WRPhiTriplet{};
  float _chi2WRPhiQuartet{};
  float _chi2WRPhiSeptet{};
  float _chi2WZTriplet{};
  float _chi2WZQuartet{};
  float _chi2WZSeptet{};
  float _minDistCutAttachForFTD{};
  float _minDistCutAttachForVXD{};
  int _minimalLayerToAttach{};
  int _minMissAddition{};
<<<<<<< HEAD
  
  // two pi is not a constant in cmath. Calculate it, once!
  static const double TWOPI;
  
  double _dPhi{};
  double _dTheta{};
  double _dPhiFTD{};
  

  
  std::vector<int> _Combinations{};
  std::vector<int> _CombinationsFTD{};
  
  float _resolutionRPhiVTX{};
  float _resolutionZVTX{};
  
  float _resolutionRPhiFTD{};
  float _resolutionZFTD{};
  
  float _resolutionRPhiSIT{};
  float _resolutionZSIT{};
  
  float _phiCutForMerging{};
  float _tanlambdaCutForMerging{};
  float _angleCutForMerging{};
  //float _angleCutForMerging_highPt;
  //float _angleCutForMerging_lowPt;
  
  int _print{};
  int _checkForDelta{};
  float _minDistToDelta{};
  
=======

  // two pi is not a constant in cmath. Calculate it, once!
  static const double TWOPI;

  double _dPhi{};
  double _dTheta{};
  double _dPhiFTD{};

  std::vector<int> _Combinations{};
  std::vector<int> _CombinationsFTD{};

  float _resolutionRPhiVTX{};
  float _resolutionZVTX{};

  float _resolutionRPhiFTD{};
  float _resolutionZFTD{};

  float _resolutionRPhiSIT{};
  float _resolutionZSIT{};

  float _phiCutForMerging{};
  float _tanlambdaCutForMerging{};
  float _angleCutForMerging{};
  // float _angleCutForMerging_highPt;
  // float _angleCutForMerging_lowPt;

  int _print{};
  int _checkForDelta{};
  float _minDistToDelta{};

>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  float _distRPhi{};
  float _distZ{};
  float _chi2FitCut{};
  float _chi2FitCut_lowPt{};
<<<<<<< HEAD
  
  
  TrackExtendedVec _trackImplVec{};
  
  
  float _cutOnD0{}, _cutOnZ0{}, _cutOnOmegaVXD{}, _cutOnOmegaFTD{};
  double _cutOnPtVXD{},_cutOnPtFTD{};
  
=======

  TrackExtendedVec _trackImplVec{};

  float _cutOnD0{}, _cutOnZ0{}, _cutOnOmegaVXD{}, _cutOnOmegaFTD{};
  double _cutOnPtVXD{}, _cutOnPtFTD{};

>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  int _minimalHits{};
  int _nHitsChi2{};
  int _attachVXD{};
  int _attachFTD{};
<<<<<<< HEAD
  
  int _max_hits_per_sector{};
  
  int _nTotalVTXHits{},_nTotalFTDHits{},_nTotalSITHits{};
  int _useSIT{};
  int _useFTD{};
  
  
  //  int _createMap;
  
  UTIL::BitField64* _encoder{nullptr};
  int getDetectorID(TrackerHit* hit) { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::subdet()]; }
  int getSideID(TrackerHit* hit)     { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::side()]; };
  int getLayerID(TrackerHit* hit)    { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::layer()]; };
  int getModuleID(TrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::module()]; };
  int getSensorID(TrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::sensor()]; };

  int getDetectorID(SimTrackerHit* hit) { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::subdet()]; }
  int getSideID(SimTrackerHit* hit)     { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::side()]; };
  int getLayerID(SimTrackerHit* hit)    { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::layer()]; };
  int getModuleID(SimTrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::module()]; };
  int getSensorID(SimTrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::sensor()]; };
  
  void setupGeom(const dd4hep::Detector& theDetector) ;
  
  
  unsigned int _nLayersVTX{};
  
  unsigned int _nLayersSIT{};
  
  
  std::vector<float> _zLayerFTD{};
  
  unsigned int _nlayersFTD{};
  bool _petalBasedFTDWithOverlaps{};
  int _nPhiFTD{};

=======

  int _max_hits_per_sector{};

  int _nTotalVTXHits{}, _nTotalFTDHits{}, _nTotalSITHits{};
  int _useSIT{};
  int _useFTD{};

  //  int _createMap;

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

  unsigned int _nLayersVTX{};

  unsigned int _nLayersSIT{};

  std::vector<float> _zLayerFTD{};

  unsigned int _nlayersFTD{};
  bool _petalBasedFTDWithOverlaps{};
  int _nPhiFTD{};

>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  int _output_track_col_quality{};
  static const int _output_track_col_quality_GOOD;
  static const int _output_track_col_quality_FAIR;
  static const int _output_track_col_quality_POOR;

<<<<<<< HEAD
  int _sw_theta{};//search window theta
=======
  int _sw_theta{}; // search window theta
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  float _chi2FitCut_kalman{};
  bool _useClusterRejection{};
  float _minDotOf2Clusters{};

  int getIntersectionEasy(HelixClass_double& helix, TrackerHit* curInmos, int layer, double* isec, double* ref);
  int getIntersectionEasyTest(HelixClass_double& helix, TrackerHit* basisTrkhit, int layer, std::vector<double>& vec);
  int CheckTiltOf2Clusters(TrackerHit* A, TrackerHit* B, int level);
  float DotOf2Clusters(TrackerHit* A, TrackerHit* B);
  int KalFit(int& ndf, float& Chi2, TrackerHitVec trkHits, TrackerHitVec& hits_in_fit, TrackerHitVec& outliers,
             float* par, float* epar, HelixClass_double& helix);
  int getPhiThetaRegion(TrackExtended* trackAR, int layer, int* Boundaries);

  struct GeoData_t {
    int nladder{};
<<<<<<< HEAD
    double rmin{};  // distance of inner surface of sensitive region from IP
    double dphi{};  // azimuthal angle step of each ladder
    double phi0{};  // aximuthal angle offset
    std::vector<double> cosphi{};  // cos[phi_ladder], cos_phi of each ladder
    std::vector<double> sinphi{};  // sin[phi_ladder], sin_phi of each ladder
    std::vector<double> phi{};  // phi of each ladder
    std::vector<double> phiAtXiMin{};  // phiAtXiMin of each ladder
    std::vector<double> phiAtXiMax{};  // phiAtXiMax of each ladder
    std::vector<double> ladder_incline{};//the tilt of the line of the ladder expressed by phi
    double sthick{};  // sensitive region thickness
    double sximin{};  // minimum xi of sensitive region.
    double sximax{};  // maximum xi of sensitive region
    double hlength{}; // ladder's half length in z
    int num_xi_pixel{};      // Number of xi pixel in this ladder
    int num_zeta_pixel{};    // Number of zeta pixel in this ladder
    double rmes{}; //distance in R of measurement surface
=======
    double rmin{};                        // distance of inner surface of sensitive region from IP
    double dphi{};                        // azimuthal angle step of each ladder
    double phi0{};                        // aximuthal angle offset
    std::vector<double> cosphi{};         // cos[phi_ladder], cos_phi of each ladder
    std::vector<double> sinphi{};         // sin[phi_ladder], sin_phi of each ladder
    std::vector<double> phi{};            // phi of each ladder
    std::vector<double> phiAtXiMin{};     // phiAtXiMin of each ladder
    std::vector<double> phiAtXiMax{};     // phiAtXiMax of each ladder
    std::vector<double> ladder_incline{}; // the tilt of the line of the ladder expressed by phi
    double sthick{};                      // sensitive region thickness
    double sximin{};                      // minimum xi of sensitive region.
    double sximax{};                      // maximum xi of sensitive region
    double hlength{};                     // ladder's half length in z
    int num_xi_pixel{};                   // Number of xi pixel in this ladder
    int num_zeta_pixel{};                 // Number of zeta pixel in this ladder
    double rmes{};                        // distance in R of measurement surface
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  };

  struct vxdGeoData_t {
<<<<<<< HEAD
  int nLayer{};
  int maxLadder{};
  std::vector<GeoData_t> geodata{};
  }_vxd{},_sit{};
=======
    int nLayer{};
    int maxLadder{};
    std::vector<GeoData_t> geodata{};
  } _vxd{}, _sit{};
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90

  void InitVXDGeometry(const dd4hep::Detector& theDetector);
  void InitSITGeometry(const dd4hep::Detector& theDetector);
  FloatVec _pixelSizeVec{};
  float _pixelheight{};
<<<<<<< HEAD
  TVector3 LocalToGlobal(TVector3 local,int layer,int ladder);
  //purityMCP CheckOriginOf2Clusters(TrackerHit* A, TrackerHit* B);
  void calcTrackParameterOfMCP(MCParticle* pmcp, double* par);
  
  class ClusterStatus{
    public :
    int cellid0 {};
    int cellid1{};
    unsigned int layer {};
    unsigned int ladder{};
    unsigned int   xiwidth {};
    unsigned int zetawidth {};
    unsigned int      nPix {};
    unsigned int      tilt {};
    unsigned int   quality {};
=======
  TVector3 LocalToGlobal(TVector3 local, int layer, int ladder);
  // purityMCP CheckOriginOf2Clusters(TrackerHit* A, TrackerHit* B);
  void calcTrackParameterOfMCP(MCParticle* pmcp, double* par);
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90

  class ClusterStatus {
  public:
    int cellid0{};
    int cellid1{};
    unsigned int layer{};
    unsigned int ladder{};
    unsigned int xiwidth{};
    unsigned int zetawidth{};
    unsigned int nPix{};
    unsigned int tilt{};
    unsigned int quality{};

    ClusterStatus(TrackerHit* hit) {
      cellid0 = hit->getCellID0();
      cellid1 = hit->getCellID1();
      // layer  = ( cellid0 >> 5 ) & 0x000001ff;//maybe bug.
      // ladder = (cellid0 >> 15 ) & 0x000000ff;
      xiwidth = cellid1 & 0x000001ff;
      zetawidth = (cellid1 >> 9) & 0x000001ff;
      nPix = (cellid1 >> 18) & 0x000003ff;
      tilt = (cellid1 >> 28) & 0x00000003;
      quality = (cellid1 >> 30) & 0x00000003;
      const int cellId_trk = hit->getCellID0();
      UTIL::BitField64 encoder_trk(lcio::LCTrackerCellID::encoding_string());
      encoder_trk.setValue(cellId_trk);
      layer = encoder_trk[lcio::LCTrackerCellID::layer()];
      ladder = encoder_trk[lcio::LCTrackerCellID::module()];
    }
    ClusterStatus() { ; }
  };

<<<<<<< HEAD
  //Old Ver//////////////////
  typedef std::map< std::pair< int, int >, int > RangeMap;
=======
  // Old Ver//////////////////
  typedef std::map<std::pair<int, int>, int> RangeMap;
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  RangeMap _phiRangeForTriplet{};

  double getNeededPhiSectors(double Pt, int outly, int inly); // Difference of radian is returned
  ///////////////////////////

<<<<<<< HEAD
  //New Ver////////////////////==under construction====
  typedef std::map< std::vector<int> , std::vector<int> > RangeMapVer2;
=======
  // New Ver////////////////////==under construction====
  typedef std::map<std::vector<int>, std::vector<int>> RangeMapVer2;
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  RangeMapVer2 _phiRangeForTripletVer2{};
  void getNeededPhiSectorsVer2(double Pt, std::vector<int> layers, std::vector<double>& phiDiff);
  /////////////////////////////========================

<<<<<<< HEAD
  float _safetyPhiRange_ratio{};/**
  Extra range in addition to main range used for triplet construction process and needed to find triplet 
  is calculated by getNeededPhiSectors, but safety range is not considered.
  Value of _safetyRange is used such as Range = Range*(1 + _safetyRange);
  */
=======
  float _safetyPhiRange_ratio{}; /**
   Extra range in addition to main range used for triplet construction process and needed to find triplet
   is calculated by getNeededPhiSectors, but safety range is not considered.
   Value of _safetyRange is used such as Range = Range*(1 + _safetyRange);
   */
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  int _safetyPhiRange_fix{};
  float _fudgeFactorForSITsr_rphi{};
  float _fudgeFactorForSITsr_z{};
  int _fudgePhiRange{};
  int _fudgeThetaRange{};

  std::string _colNameVXDTrackerHitRelations{};
  std::string _colNameSITSpacePointRelations{};
  std::string _colNameFTDSpacePointRelations{};
  std::string _colNameFTDPixelTrackerHitRelations{};
  LCRelationNavigator* _navVXD{nullptr};
  LCRelationNavigator* _navSIT{nullptr};
  LCRelationNavigator* _navFTDsp{nullptr};
  LCRelationNavigator* _navFTDpix{nullptr};
  std::string _colNameVXDSimHit{};
  std::string _colNameSITSimHit{};
  std::string _colNameFTDspSimHit{};
  std::string _colNameFTDpixSimHit{};
  LCCollection* _simVXD{nullptr};
  LCCollection* _simSIT{nullptr};
  LCCollection* _simFTDsp{nullptr};
  LCCollection* _simFTDpix{nullptr};

  double _nSigmaBuild_phi{};
  double _nSigmaBuild_theta{};

  bool _keepCandidate{}; // used in AttachRemainingVTXHitsVeryFast <-- under construction for now

<<<<<<< HEAD
  double _nSigmaBuild_phi{};
  double _nSigmaBuild_theta{};

  bool _keepCandidate{};//used in AttachRemainingVTXHitsVeryFast <-- under construction for now

  moriUTIL* _moriUtil{nullptr};
  GetPurityUtil* _purityUtil{nullptr};

////////////////////////////////////////////////////////////////
////from here, a lot of debug tools and variables for mori /////
////////////////////////////////////////////////////////////////
=======
  moriUTIL* _moriUtil{nullptr};
  GetPurityUtil* _purityUtil{nullptr};

  ////////////////////////////////////////////////////////////////
  ////from here, a lot of debug tools and variables for mori /////
  ////////////////////////////////////////////////////////////////
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  bool _mydebug{};
  bool _mydebugKalFit{};
  bool _mydebugIntersection{};
  bool _mydebugRoot{};
  bool _mydebugRootVXD{};
  bool _mydebugRootFTD{};
  bool _mydebugRootVXDFTD{};
  bool _mydebugVXDHits{};
  bool _mydebugSITHits{};
  bool _mydebugTriplet{};
<<<<<<< HEAD
  int  _mydebugTripletMode{};
=======
  int _mydebugTripletMode{};
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  bool _mydebugTripletVXD{};
  bool _mydebugTripletFTD{};
  bool _mydebugTripletVXDFTD{};
  bool _mydebugBuildTrack{};
<<<<<<< HEAD
  int  _mydebugBuildTrackMode{};
  bool _mydebugAttachVXD{};
  bool _mydebugPrintMCP{};
  bool _stopwatch{};
  class Timer{
   public:
=======
  int _mydebugBuildTrackMode{};
  bool _mydebugAttachVXD{};
  bool _mydebugPrintMCP{};
  bool _stopwatch{};
  class Timer {
  public:
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
    TStopwatch inAnEvt{};
    TStopwatch InitialiseVTX{};
    TStopwatch InitialiseFTD{};
    TStopwatch ProcessOneSector{};
    TStopwatch TrackingInFTD{};
    TStopwatch Sorting{};
    TStopwatch CreateTrack{};
    TStopwatch AttachRemainingVXD{};
    TStopwatch AttachRemainingFTD{};
    TStopwatch FinalRefit{};
<<<<<<< HEAD
    void reset(){
=======
    void reset() {
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
      inAnEvt.Reset();
      InitialiseVTX.Reset();
      InitialiseFTD.Reset();
      ProcessOneSector.Reset();
      TrackingInFTD.Reset();
      Sorting.Reset();
      CreateTrack.Reset();
      AttachRemainingVXD.Reset();
      AttachRemainingFTD.Reset();
      FinalRefit.Reset();
    }
    void cout() {
      printf("inAnEvt           : RT=%.3f s, CPU=%.3f s \n", inAnEvt.RealTime(), inAnEvt.CpuTime());
      printf("InitialiseVTX     : RT=%.3f s, CPU=%.3f s \n", InitialiseVTX.RealTime(), InitialiseVTX.CpuTime());
      printf("InitialiseFTD     : RT=%.3f s, CPU=%.3f s \n", InitialiseFTD.RealTime(), InitialiseFTD.CpuTime());
      printf("ProcessOneSector  : RT=%.3f s, CPU=%.3f s \n", ProcessOneSector.RealTime(), ProcessOneSector.CpuTime());
      printf("TrackingInFTD     : RT=%.3f s, CPU=%.3f s \n", TrackingInFTD.RealTime(), TrackingInFTD.CpuTime());
      printf("Sorting           : RT=%.3f s, CPU=%.3f s \n", Sorting.RealTime(), Sorting.CpuTime());
      printf("CreateTrack       : RT=%.3f s, CPU=%.3f s \n", CreateTrack.RealTime(), CreateTrack.CpuTime());
      printf("AttachRemainingVXD: RT=%.3f s, CPU=%.3f s \n", AttachRemainingVXD.RealTime(),
             AttachRemainingVXD.CpuTime());
      printf("AttachRemainingFTD: RT=%.3f s, CPU=%.3f s \n", AttachRemainingFTD.RealTime(),
             AttachRemainingFTD.CpuTime());
      printf("FinalRefit        : RT=%.3f s, CPU=%.3f s \n", FinalRefit.RealTime(), FinalRefit.CpuTime());
    }
<<<<<<< HEAD
    Timer(){reset();}
  }_timer{};
  TStopwatch _timer2Triplet{};
  TStopwatch _timer2Build{};
  bool  _mydebugstopwatch2{};

  float _currentPurity{};
  MCPMap LoadMCPMap();
  std::map< MCParticle*, SimTrackerHitVec > _mcpVXD{};
  std::map< MCParticle*, SimTrackerHitVec > _mcpVXDFTD{};
  std::map< MCParticle*, SimTrackerHitVec > _mcpVXDSIT{};
  std::map< MCParticle*, SimTrackerHitVec > _mcpVXDFTDSIT{};
  std::map< MCParticle*, SimTrackerHitVec > _mcpFTD{};
  std::map< MCParticle*, SimTrackerHitVec > _mcpFTDSIT{};
  std::map< MCParticle*, SimTrackerHitVec > _mcpSIT{};
=======
    Timer() { reset(); }
  } _timer{};
  TStopwatch _timer2Triplet{};
  TStopwatch _timer2Build{};
  bool _mydebugstopwatch2{};

  float _currentPurity{};
  MCPMap LoadMCPMap();
  std::map<MCParticle*, SimTrackerHitVec> _mcpVXD{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpVXDFTD{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpVXDSIT{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpVXDFTDSIT{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpFTD{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpFTDSIT{};
  std::map<MCParticle*, SimTrackerHitVec> _mcpSIT{};
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90

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

  enum MCPContributions getMCPContribution(IntVec nsub);

<<<<<<< HEAD
  class Triplet : public IMPL::TrackImpl{
   public :
     static int numOfProcesses;
     int id{};
     int quality_code{};
     int detID[3]{};//内側に近いところから詰める。
     int layer[3]{};
     float probability{};
     float purity{};
     SimTrackerHitVec simVecFromTrk{};//TrackerHitVecに対応するもの
     std::pair< MCParticle* , SimTrackerHitVec > mcTruth {};//Dominant MCPについての情報

     //new (2013_08_10)
     TrackerHitVec truthTrkHits {};//Dominant MCPのsimthitsに対応するTrackerHitVec.
     enum MCPContributions mcpcont{};
     IntVec nsub{};//vxd:0,sit:1,ftd:2 <--for Dominant MCP
     

     float mcp_d0{};
     float mcp_z0{};
     float mcp_omega{};
     float mcp_phi0{};
     float mcp_tanL{};
=======
  class Triplet : public IMPL::TrackImpl {
  public:
    static int numOfProcesses;
    int id{};
    int quality_code{};
    int detID[3]{}; // 内側に近いところから詰める。
    int layer[3]{};
    float probability{};
    float purity{};
    SimTrackerHitVec simVecFromTrk{};                   // TrackerHitVecに対応するもの
    std::pair<MCParticle*, SimTrackerHitVec> mcTruth{}; // Dominant MCPについての情報

    // new (2013_08_10)
    TrackerHitVec truthTrkHits{}; // Dominant MCPのsimthitsに対応するTrackerHitVec.
    enum MCPContributions mcpcont {};
    IntVec nsub{}; // vxd:0,sit:1,ftd:2 <--for Dominant MCP

    float mcp_d0{};
    float mcp_z0{};
    float mcp_omega{};
    float mcp_phi0{};
    float mcp_tanL{};
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90

    // bool availableInBuildedTrack;

<<<<<<< HEAD
     Triplet(bool incrementID = true){
       id = numOfProcesses; 
       if(incrementID){
         numOfProcesses++;
       }
       quality_code = -1;
       detID[0] = detID[1] = detID[2] = -1;
       layer[0] = layer[1] = layer[2] = -1;
       probability = 1e20;
       purity = 1e20;
       simVecFromTrk.clear();
       mcTruth.first = 0; mcTruth.second.clear();
       truthTrkHits.clear();
       mcpcont = contSize;
       nsub.clear();
       mcp_d0 = mcp_z0 = mcp_omega = mcp_phi0 = mcp_tanL = 1e20;
     }
     Triplet(const Triplet&) = default;
  };
  Triplet* _curtriplet{nullptr};
  

  enum BuildTrackResult{
    NormalEnd,
    ManyMisAssignments,
    ManyOutliers,
    PhiThetaError,
    BuildTrackResultSize
  };
  class BuildedTrack : public IMPL::TrackImpl{
   public :
     static int numOfProcesses;
     int id{};

     BuildTrackResult result{};
     int nMisAssign{};
     Triplet triplet{};
     const MCParticle* truthmcp{nullptr};
     SimTrackerHitVec simvec{};
     float probability{};
     float purity{};
     SimTrackerHitVec simVecFromTrk{};//TrackerHitVecに対応するもの
     BuildedTrack():triplet(false){
        id = numOfProcesses; 
        numOfProcesses++;
        result = BuildTrackResultSize;
        nMisAssign = -1;
        truthmcp = 0;
        simvec.clear();
        probability = 1e20;
        purity = 1e20;
        simVecFromTrk.clear();
     }
     BuildedTrack(const BuildedTrack&) = default;
     BuildedTrack& operator=(const BuildedTrack&) = default;
  };
  bool _availableInBuildedTrack{};

  typedef std::vector< BuildedTrack > BuildedTrackVec;
=======
    Triplet(bool incrementID = true) {
      id = numOfProcesses;
      if (incrementID) {
        numOfProcesses++;
      }
      quality_code = -1;
      detID[0] = detID[1] = detID[2] = -1;
      layer[0] = layer[1] = layer[2] = -1;
      probability = 1e20;
      purity = 1e20;
      simVecFromTrk.clear();
      mcTruth.first = 0;
      mcTruth.second.clear();
      truthTrkHits.clear();
      mcpcont = contSize;
      nsub.clear();
      mcp_d0 = mcp_z0 = mcp_omega = mcp_phi0 = mcp_tanL = 1e20;
    }
    Triplet(const Triplet&) = default;
  };
  Triplet* _curtriplet{nullptr};

  enum BuildTrackResult { NormalEnd, ManyMisAssignments, ManyOutliers, PhiThetaError, BuildTrackResultSize };
  class BuildedTrack : public IMPL::TrackImpl {
  public:
    static int numOfProcesses;
    int id{};

    BuildTrackResult result{};
    int nMisAssign{};
    Triplet triplet{};
    const MCParticle* truthmcp{nullptr};
    SimTrackerHitVec simvec{};
    float probability{};
    float purity{};
    SimTrackerHitVec simVecFromTrk{}; // TrackerHitVecに対応するもの
    BuildedTrack() : triplet(false) {
      id = numOfProcesses;
      numOfProcesses++;
      result = BuildTrackResultSize;
      nMisAssign = -1;
      truthmcp = 0;
      simvec.clear();
      probability = 1e20;
      purity = 1e20;
      simVecFromTrk.clear();
    }
    BuildedTrack(const BuildedTrack&) = default;
    BuildedTrack& operator=(const BuildedTrack&) = default;
  };
  bool _availableInBuildedTrack{};

  typedef std::vector<BuildedTrack> BuildedTrackVec;
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  BuildedTrackVec _buildedTrackContainer{};
  void SetBuildedTrack(TrackExtended* trackAR, BuildedTrackVec& btrackvec);
  void BuildedTrackDebuger1(BuildedTrackVec::iterator begin, BuildedTrackVec::iterator end);
  void BuildedTrackDebuger2(std::vector<const BuildedTrack*>::iterator begin,
                            std::vector<const BuildedTrack*>::iterator end);

<<<<<<< HEAD
  //時間がアレばclass templateにしてみたい。
  class MCP_BuildedTrack{
   public :
    std::pair<MCParticle*,SimTrackerHitVec> pair{};
    IntVec nsub{};//nhits in sub detectors. 0 VXD, 1 FTD, 2 SIT
=======
  // 時間がアレばclass templateにしてみたい。
  class MCP_BuildedTrack {
  public:
    std::pair<MCParticle*, SimTrackerHitVec> pair{};
    IntVec nsub{}; // nhits in sub detectors. 0 VXD, 1 FTD, 2 SIT
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
    std::vector<const BuildedTrack*> BuildedTrackVec{};
    std::vector<const BuildedTrack*> BuildedTrack2ndVec{};
    std::vector<const BuildedTrack*> BuildedTrack3rdVec{};
  };
<<<<<<< HEAD
  typedef std::map<MCParticle*,MCP_BuildedTrack> BTMap; 
  BTMap _mcpBuildedTrackContainer{};

  std::vector< std::map<MCParticle*, MCP_BuildedTrack> > _mcpRemainingBTContainerA{};
  //purity 100%のトリプレットが生成されなかった生成されなかったmcp
  std::vector< std::map<MCParticle*, MCP_BuildedTrack> > _mcpRemainingBTContainerB{};
  //purity 100%のトリプレットが生成されたものの、quality_code != 0により結局再構成されなかったmcp
  
=======
  typedef std::map<MCParticle*, MCP_BuildedTrack> BTMap;
  BTMap _mcpBuildedTrackContainer{};

  std::vector<std::map<MCParticle*, MCP_BuildedTrack>> _mcpRemainingBTContainerA{};
  // purity 100%のトリプレットが生成されなかった生成されなかったmcp
  std::vector<std::map<MCParticle*, MCP_BuildedTrack>> _mcpRemainingBTContainerB{};
  // purity 100%のトリプレットが生成されたものの、quality_code != 0により結局再構成されなかったmcp
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90

  std::vector<Triplet> _tripletContainer{};

<<<<<<< HEAD

  std::vector< Triplet > _tripletContainer{};

  class MCP_Triplet{
   public :
    std::pair<MCParticle*,SimTrackerHitVec> pair{};
    IntVec nsub{};//nhits in sub detectors. 0 VXD, 1 FTD, 2 SIT
=======
  class MCP_Triplet {
  public:
    std::pair<MCParticle*, SimTrackerHitVec> pair{};
    IntVec nsub{}; // nhits in sub detectors. 0 VXD, 1 FTD, 2 SIT
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
    std::vector<const Triplet*> tripvec{};
    std::vector<const Triplet*> trip2ndvec{};
    std::vector<const Triplet*> trip3rdvec{};
  };
<<<<<<< HEAD
  typedef std::map<MCParticle*,MCP_Triplet> TripMap; 
=======
  typedef std::map<MCParticle*, MCP_Triplet> TripMap;
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90
  TripMap _mcpTripletContainer{};

  std::vector<std::map<MCParticle*, MCP_Triplet>> _mcpRemainingTRContainerA{};
  // purity 100%のトリプレットが生成されなかった生成されなかったmcp
  std::vector<std::map<MCParticle*, MCP_Triplet>> _mcpRemainingTRContainerB{};
  // purity 100%のトリプレットが生成されたものの、quality_code != 0により結局再構成されなかったmcp

<<<<<<< HEAD
  std::vector< std::map<MCParticle*, MCP_Triplet> > _mcpRemainingTRContainerA{};
  //purity 100%のトリプレットが生成されなかった生成されなかったmcp
  std::vector< std::map<MCParticle*, MCP_Triplet> > _mcpRemainingTRContainerB{};
  //purity 100%のトリプレットが生成されたものの、quality_code != 0により結局再構成されなかったmcp
  
  void TripletDebugerWithMCPRemain(int nth, std::vector< std::map<MCParticle*, MCP_Triplet> > A);
=======
  void TripletDebugerWithMCPRemain(int nth, std::vector<std::map<MCParticle*, MCP_Triplet>> A);
>>>>>>> 665b5fed8309e84cc4197b80012c02c680e32b90

  MCPMap _mcpMap{};
  std::vector<LCRelationNavigator*> _naviVec{};
  IntVec getNHitsInSubDet(SimTrackerHitVec simvec);
  void TripletDebuger1(std::vector<Triplet>::iterator begin, std::vector<Triplet>::iterator end);
  void TripletDebuger2(std::vector<const Triplet*>::iterator begin, std::vector<const Triplet*>::iterator end);
};

int FPCCDSiliconTracking_MarlinTrk::Triplet::numOfProcesses = 0;
int FPCCDSiliconTracking_MarlinTrk::BuildedTrack::numOfProcesses = 0;

#endif
