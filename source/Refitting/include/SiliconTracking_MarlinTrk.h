#ifndef SILICONTRACKING_MarlinTrk_H
#define SILICONTRACKING_MarlinTrk_H 1

#include "marlin/Processor.h"
#include <marlin/Global.h>
#include "lcio.h"
#include <string>
#include <vector>
#include <cmath>
#include <IMPL/TrackImpl.h>
#include "ClusterExtended.h"
#include "TrackExtended.h"
#include "TrackerHitExtended.h"
#include "HelixClass.h"

#include "MarlinTrk/IMarlinTrack.h"

#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"

namespace DiagnosticsHistograms {
  class Histograms ;
}

namespace DiagnosticsHistograms2D {
  class Histograms ;
}

using namespace lcio ;
using namespace marlin ;


namespace MarlinTrk {
  class HelixFit;
  class IMarlinTrkSystem ;
}

namespace dd4hep{
    class Detector ;
}

namespace UTIL{
  class LCRelationNavigator ;
}


/** === Silicon Tracking Processor === <br>
 * Processor performing stand-alone pattern recognition
 * in the vertex detector (VTX), forward tracking disks and SIT. <br>
 * The procedure consists of three steps : <br> 
 * 1) Tracking in VTX and SIT ; <br>
 * 2) Tracking in FTD ; <br>
 * 3) Merging compatible track segments reconstructed in VTX and FTD <br>
 * STEP 1 : TRACKING IN VTX and SIT <br>
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
 * is the same as used for tracking in the VTX+SIT. <br>
 * STEP 3 : MERGING TRACK SEGMENTS FOUND IN FTD AND VTX+SIT <br>
 * In the last step, track segments reconstructed in the FTD and VTX+SIT, belonging to the
 * same track  are identified and merged into one track. All possible 
 * pairings are tested for their compatibility.
 * The number of pairings considered is Ntrk_VTX_SIT*Ntrk_FTD, where Ntrk_VTX_SIT is the number of 
 * track segments reconstructed in the first step in VTX+SIT (segments containing solely VTX and SIT hits) and
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
 * Z0 (displacement along z axis at the point of closest approach to IP). Covariance matrix for these parameters is also provided.
 * Only lower left corner of the covariance matrix is stored. The sequence of the covariance matrix elements 
 * assigned to track is the following: <br>
 * (Omega,Omega) <br>
 * (Omega,TanLambda), (TanLambda,TanLambda) <br>
 * (Omega,Phi), (TanLamda,Phi), (Phi,Phi) <br>
 * (Omega,D0), (TanLambda,D0), (Phi,D0), (D0,D0) <br>
 * (Omega,Z0), (TanLambda,Z0), (Phi,Z0), (D0,Z0), (Z0,Z0) <br>
 * The number of hits in the different subdetectors associated
 * with each track can be accessed via method Track::getSubdetectorHitNumbers().
 * This method returns vector of integers : <br>
 * number of VTX hits in track is the first element in this vector  
 * (Track::getSubdetectorHitNumbers()[0]) <br>
 * number of FTD hits in track is the second element in this vector  
 * (Track::getSubdetectorHitNumbers()[1]) <br>
 * number of SIT hits in track is the third element in this vector  
 * (Track::getSubdetectorHitNumbers()[2]) <br>
 * Output track collection has a name "SiTracks". <br>
 * @param VTXHitCollectionName name of input VTX TrackerHit collection <br>
 * (default parameter value : "VTXTrackerHits") <br>
 * @param FTDHitCollectionName name of input FTD TrackerHit collection <br>
 * (default parameter value : "FTDTrackerHits") <br>
 * @param SITHitCollectionName name of input SIT TrackerHit collection <br>
 * (default parameter value : "SITTrackerHits") <br>
 * @param SiTrackCollectionName name of the output Silicon track collection <br>
 * (default parameter value : "SiTracks") <br>
 * @param LayerCombinations combinations of layers used to search for hit triplets in VTX+SIT <br>
 * (default parameters : 6 4 3  6 4 2  6 3 2  5 4 3  5 4 2  5 3 2  4 3 2  4 3 1  4 2 1  3 2 1) <br> 
 * Note that in the VTX+SIT system the first and the second layers of SIT have indicies nLayerVTX and nLayerVTX+1. 
 * Combination given above means that triplets are looked first in layers 6 4 3, and then 
 * in 6 4 2;  5 4 3;  6 3 2 etc. NOTE THAT LAYER INDEXING STARTS FROM 0.
 * LAYER 0 is the innermost layer  <br>
 * @param LayerCombinationsFTD combinations of layers used to search for hit triplets in FTD <br>
 * (default parameters 6 5 4  5 4 3  5 4 2  5 4 1  5 3 2  5 3 1  5 2 1  4 3 2  4 3 1  
 *  4 3 0  4 2 1  4 2 0  4 1 0  3 2 1  3 2 0  3 1 0  2 1 0). 
 * NOTE THAT TRACKS IN FTD ARE SEARCHED ONLY IN ONE HEMISPHERE. TRACK IS NOT 
 * ALLOWED TO HAVE HITS BOTH IN BACKWARD AND FORWARD PARTS OF FTD SIMULTANEOUSLY. 
 * @param NDivisionsInPhi Number of divisions in Phi for tracking in VTX+SIT <br>
 * (default value is 40) <br>
 * @param NDivisionsInTheta Number of divisions in cosQ for tracking in VTX+SIT <br>
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
 * @param MinLayerToAttach the minimal layer index to attach VTX hits to the found hit triplets <br>
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
 * @param FastAttachment if this flag is set to 1, less accurate but fast procedure to merge additional hits to tracks is used <br> 
 * if set to 0, a more accurate, but slower procedure is invoked <br>
 * (default value is 0) <br>
 * @param UseSIT When this flag is set to 1, SIT is included in pattern recognition. When this flag is set
 * to 0, SIT is excluded from the procedure of pattern recognition <br>
 * (default value is 1) <br>
 * <br>
 * @author A. Raspereza (MPI Munich)<br>
 */
class SiliconTracking_MarlinTrk : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new SiliconTracking_MarlinTrk ; }
  
  
  SiliconTracking_MarlinTrk() ;
  
  /**  
   * Initialization
   */
  virtual void init() ;
  
  /** Run header processor.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Event processor.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
protected:
  
  int _nRun ;
  int _nEvt ;
  EVENT::LCEvent* _current_event;
  
  int _nDivisionsInPhi;
  int _nDivisionsInTheta;
  int _nLayers;
  
  MarlinTrk::HelixFit* _fastfitter;
  
  /** pointer to the IMarlinTrkSystem instance 
   */
  MarlinTrk::IMarlinTrkSystem* _trksystem ;
  bool _runMarlinTrkDiagnostics;
  std::string _MarlinTrkDiagnosticsName;
  
  bool _MSOn, _ElossOn, _SmoothOn ;

  /** switches: False for backward compatible (default).
                True to apply new methods.
   */
  bool _useSimpleUpdatedCoreBin;
  bool _useSimpleAttachHitToTrack;

  float _initialTrackError_d0;
  float _initialTrackError_phi0;
  float _initialTrackError_omega;
  float _initialTrackError_z0;
  float _initialTrackError_tanL;
  
  double _maxChi2PerHit;  
  
  bool  _UseEventDisplay;
  int _detector_model_for_drawing;
  std::vector<int> _colours;  
  float     _helix_max_r;
  
  void drawEvent();
  
  
  // histogram member variables
  
  bool  _createDiagnosticsHistograms;
  DiagnosticsHistograms::Histograms* _histos{nullptr};

  
  int _ntriplets, _ntriplets_good, _ntriplets_2MCP, _ntriplets_3MCP, _ntriplets_1MCP_Bad, _ntriplets_bad;
  
  
  /** helper function to get collection using try catch block */
  LCCollection* GetCollection(  LCEvent * evt, std::string colName ) ;
  
  /** helper function to get relations using try catch block */
  LCRelationNavigator* GetRelations( LCEvent * evt, std::string RelName ) ;
  
  /** input MCParticle collection and threshold used for Drawing
   */
  std::string  _colNameMCParticles;
  float _MCpThreshold ;
  
  
  /// Compare tracks according to their chi2/ndf
  struct compare_TrackExtended{
    // n.b.: a and b should be TrackExtended const *, but the getters are not const :-(
    bool operator()(TrackExtended *a, TrackExtended *b) const {
      if ( a == b ) return false;
      return (a->getChi2()/a->getNDF() < b->getChi2()/b->getNDF() );
    }
  };
  
  
  std::string _VTXHitCollection;
  std::string _FTDPixelHitCollection;
  std::string _FTDSpacePointCollection;
  std::string _SITHitCollection;
  std::string _siTrkCollection;
  
  std::vector< LCCollection* > _colTrackerHits;
  std::map< LCCollection*, std::string > _colNamesTrackerHits;
  
  std::vector<TrackerHitExtendedVec> _sectors;
  std::vector<TrackerHitExtendedVec> _sectorsFTD;
  
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
      _tracksNHits.resize(maxHits-2);
      _maxIndex=(maxHits-3);
    }
    
    // Sort all track vectors according to chi2/nDof
    //      void sort();
    
    /// Returns the  TrackExtendedVec for track with n hits. 
    /// In case n is larger than the maximal number the vector with the largest n ist returned.
    /// \attention The smallest valid number is three! For
    /// performance reasons there is no safety check!
    inline TrackExtendedVec & getTracksWithNHitsVec( size_t nHits ) {
      //return _tracksNHits[ std::min(nHits-3, _maxIndex) ];
      // for debugging: with boundary check
      return _tracksNHits.at(std::min(nHits-3, _maxIndex));
    }
    
  protected:
    std::vector< TrackExtendedVec > _tracksNHits;
    size_t _maxIndex; /// local cache variable to avoid calculation overhead
  };
  
  TracksWithNHitsContainer _tracksWithNHitsContainer;
  
  int InitialiseVTX(LCEvent * evt);
  int InitialiseFTD(LCEvent * evt);
  void ProcessOneSector(int iSectorPhi, int iSectorTheta);
  void CleanUp();
  TrackExtended * TestTriplet(TrackerHitExtended * outerHit, 
                              TrackerHitExtended * middleHit,
                              TrackerHitExtended * innerHit,
                              HelixClass & helix);
  
  int BuildTrack(TrackerHitExtended * outerHit, 
                 TrackerHitExtended * middleHit,
                 TrackerHitExtended * innerHit,
                 HelixClass & helix, 
                 int innerlayer,
                 int iPhiLow, int iPhiUp,
                 int iTheta, int iThetaUp,
                 TrackExtended * trackAR);
  
  void Sorting( TrackExtendedVec & trackVec);
  void CreateTrack(TrackExtended * trackAR );
  void AttachRemainingVTXHitsSlow();
  void AttachRemainingFTDHitsSlow();
  void AttachRemainingVTXHitsFast();
  void AttachRemainingFTDHitsFast();
  void TrackingInFTD();
  int BuildTrackFTD(TrackExtended* trackAR, int* nLR, int iS);
  int AttachHitToTrack(TrackExtended * trackAR, TrackerHitExtended * hit, int iopt);
  
  void FinalRefit(LCCollectionVec* trk_col, LCCollectionVec* rel_col);
  
  float _bField;
  float _chi2WRPhiTriplet;
  float _chi2WRPhiQuartet;
  float _chi2WRPhiSeptet;
  float _chi2WZTriplet;
  float _chi2WZQuartet;
  float _chi2WZSeptet;
  float _minDistCutAttach;
  int _minimalLayerToAttach;
  
  // two pi is not a constant in cmath. Calculate it, once!
  static const double TWOPI;
  
  double _dPhi;
  double _dTheta;
  double _dPhiFTD;
  

  
  std::vector<int> _Combinations;
  std::vector<int> _CombinationsFTD;
  
  float _resolutionRPhiVTX;
  float _resolutionZVTX;
  
  float _resolutionRPhiFTD;
  float _resolutionZFTD;
  
  float _resolutionRPhiSIT;
  float _resolutionZSIT;
  
  float _phiCutForMerging;
  float _tanlambdaCutForMerging;
  float _angleCutForMerging;
  
  int _print;
  int _checkForDelta;
  float _minDistToDelta;
  
  float _distRPhi;
  float _distZ;
  float _chi2FitCut;
  
  
  TrackExtendedVec _trackImplVec;
  
  
  float _cutOnD0, _cutOnZ0, _cutOnOmega, _cutOnPt;
  
  int _minimalHits;
  int _nHitsChi2;
  int _attachFast;
  
  int _max_hits_per_sector;
  
  int _nTotalVTXHits,_nTotalFTDHits,_nTotalSITHits;
  int _useSIT;

  std::string _trkSystemName ;
  
  //  int _createMap;
  
  UTIL::BitField64* _encoder;
  int getDetectorID(TrackerHit* hit) { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::subdet()]; }
  int getSideID(TrackerHit* hit)     { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::side()]; };
  int getLayerID(TrackerHit* hit)    { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::layer()]; };
  int getModuleID(TrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::module()]; };
  int getSensorID(TrackerHit* hit)   { _encoder->setValue(hit->getCellID0()); return (*_encoder)[lcio::LCTrackerCellID::sensor()]; };
  
  void setupGeom(const dd4hep::Detector& theDetector);
  
  
  unsigned int _nLayersVTX;
  
  unsigned int _nLayersSIT;
  
  
  std::vector<float> _zLayerFTD;
  
  unsigned int _nlayersFTD;
  bool _petalBasedFTDWithOverlaps;
  int _nPhiFTD; 

  int _output_track_col_quality;
  static const int _output_track_col_quality_GOOD;
  static const int _output_track_col_quality_FAIR;
  static const int _output_track_col_quality_POOR;

  
} ;

#endif



