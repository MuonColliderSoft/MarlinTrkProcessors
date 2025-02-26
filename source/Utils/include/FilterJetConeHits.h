#ifndef FilterJetConeHits_h
#define FilterJetConeHits_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

#include <TH1F.h>

#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "IMPL/ReconstructedParticleImpl.h"

using namespace lcio ;
using namespace marlin ;


/** Utility processor that selects and saves the tracker hits that are included in
 *  a DeltaR cone around the jet direction direction along with the corresponding
 *  sim hits and the reco-sim relations.
 *
 *  @parameter MCParticleCollection name of the MCParticle collection
 *  @parameter TrackerHitInputCollections name of the tracker hit input collections
 *  @parameter TrackerSimHitInputCollections name of the tracker simhit input collections
 *  @parameter TrackerHitInputRelations name of the tracker hit relation input collections
 *  @parameter TrackerHitOutputCollections name of the tracker hit output collections
 *  @parameter TrackerSimHitOutputCollections name of the tracker simhit output collections
 *  @parameter TrackerHitOutputRelations name of the tracker hit relation output collections
 *  @parameter DeltaRCut maximum angular distance between the hits and the particle direction
 *  @parameter FillHistograms flag to fill the diagnostic histograms
 *
 * @author L. Sestini, INFN Padova
 * @date  18 March 2021
 * @version $Id: $
 */

class FilterJetConeHits : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new FilterJetConeHits ; }
  
  
  FilterJetConeHits() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;

  bool filterJetBib(ReconstructedParticle* jet);

  void saveJet( ReconstructedParticle* jet, LCCollectionVec* jetsColl );
  
  void directionCorrection( const double* p, double* pcorr ) ;
  
 protected:

  // --- Input/output collection names:
  std::string m_inputJetCaloCollName{} ;
  std::string m_filteredJetCaloCollName{} ;
  std::vector<std::string> m_inputTrackerHitsCollNames{} ;
  std::vector<std::string> m_inputTrackerSimHitsCollNames{} ;
  std::vector<std::string> m_inputTrackerHitRelNames{} ;
  std::vector<std::string> m_outputTrackerHitsCollNames{} ;
  std::vector<std::string> m_outputTrackerSimHitsCollNames{} ;
  std::vector<std::string> m_outputTrackerHitRelNames{} ;


  // --- Processor parameters:
  bool m_fillHistos{} ;
  double m_deltaRCut{} ;

  // Jet filter parameters with BIB:
  double m_minDaughterMaxPt{} ;
  int m_minNTracks{} ;
  bool m_createFilteredJets{} ;

  // Jet direction correction params
  bool m_makeDirCorrection{} ;
  double m_corrConst{};
  double m_corrLin{};
  double m_corrQuad{};
  double m_corrCub{};

  // --- Diagnostic histograms:
  TH1F* m_dist = nullptr ;
  
  // --- Run and event counters:
  int _nRun{} ;
  int _nEvt{} ;

} ;

#endif



