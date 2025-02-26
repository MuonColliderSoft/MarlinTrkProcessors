#include "FilterJetConeHits.h"
#include <iostream>
#include <cmath>
#include <set>


#include <math.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>

#include <marlin/AIDAProcessor.h>
#include <marlinutil/GeometryUtil.h>

#include "HelixClass_double.h"

using namespace lcio ;
using namespace marlin ;


FilterJetConeHits aFilterJetConeHits ;


FilterJetConeHits::FilterJetConeHits() : Processor("FilterJetConeHits") {

  // --- Processor description:

  _description = "FilterJetConeHits selects tracker hits in a cone opened around a jet cone direction";

    
  // --- Processor parameters:
  
  registerProcessorParameter("JetCaloCollection",
			     "Name of the JetCalo collection",
			     m_inputJetCaloCollName,
			     std::string("JetCalo") );

  registerProcessorParameter("TrackerHitInputCollections",
			     "Name of the tracker hit input collections",
			     m_inputTrackerHitsCollNames,
			     {} );

  registerProcessorParameter("TrackerSimHitInputCollections",
			     "Name of the tracker simhit input collections",
			     m_inputTrackerSimHitsCollNames,
			     {} );

  registerProcessorParameter("TrackerHitInputRelations",
			     "Name of the tracker hit relation collections",
			     m_inputTrackerHitRelNames,
			     {} );

  registerProcessorParameter("TrackerHitOutputCollections",
			     "Name of the tracker hit output collections",
			     m_outputTrackerHitsCollNames,
			     {} );

  registerProcessorParameter("TrackerSimHitOutputCollections",
			     "Name of the tracker simhit output collections",
			     m_outputTrackerSimHitsCollNames,
			     {} );

  registerProcessorParameter("TrackerHitOutputRelations",
			     "Name of the tracker hit relation collections",
			     m_outputTrackerHitRelNames,
			     {} );

  registerProcessorParameter( "DeltaRCut" ,
			      "Maximum angular distance between the hits and the particle direction" ,
			      m_deltaRCut,
			      double(1.) );

  registerProcessorParameter( "FillHistograms",
			      "Flag to fill the diagnostic histograms",
			      m_fillHistos,
			      false );

  registerProcessorParameter( "MinJetTracks",
			      "Min number of tracks in jet to use it as filter",
			      m_minNTracks,
			      int(1) );  
  
  registerProcessorParameter( "MinDaughterMaxPt",
			      "Min pT of the highest-pT track in jet to use it as filter",
			      m_minDaughterMaxPt,
			      double(2.) );   

  registerProcessorParameter( "FilteredJetsCollectionName",
            "Name of the (optional) output filtered jets collection",
            m_filteredJetCaloCollName,
            std::string("FilteredJets") ); 

  registerProcessorParameter( "MakeFilteredJetsCollection",
            "Flag for producing collection of filtered jets",
            m_createFilteredJets,
            bool(false) );

  registerProcessorParameter( "MakeDirectionCorrection",
            "Flag for correcting direction of filtered jets",
            m_makeDirCorrection,
             bool(false) );
 
  registerProcessorParameter( "CorrectionConstant",
            "Direction correction constant term",
            m_corrConst,
            double(0.) );

  registerProcessorParameter( "CorrectionLinear",
            "Direction correction linear term",
            m_corrLin,
            double(1.) );

  registerProcessorParameter( "CorrectionQuadratic",
            "Direction correction quadratic term",
            m_corrQuad,
            double(0.) );

  registerProcessorParameter( "CorrectionCubic",
            "Direction correction cubic term",
            m_corrCub,
            double(0.) );
}



void FilterJetConeHits::init() { 

  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  
  // --- Print the processor parameters:

  printParameters() ;

  
  // --- Initialize the run and event counters:

  _nRun = 0 ;
  _nEvt = 0 ;


  // --- Initialize the AIDAProcessor and book the diagnostic histograms: 

  AIDAProcessor::histogramFactory(this);

  m_dist = new TH1F("m_dist", "deltaR distance between hit and jet axis", 1000, 0., 6.);

}


void FilterJetConeHits::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;

} 



void FilterJetConeHits::processEvent( LCEvent * evt ) { 

  // --- Check whether the number of input and output collections match
  bool noSimColl = false, noRelColl = false;
  if ( m_inputTrackerSimHitsCollNames.size() == 0 ) noSimColl = true;
  if ( m_inputTrackerHitRelNames.size() == 0 ) noRelColl = true;
  if ( noSimColl != noRelColl ) {
    std::stringstream err_msg;
    err_msg << "Error: RelColl are necessary to save SimColl!"
	    << std::endl ;

    throw EVENT::Exception( err_msg.str() ) ;
  }

  if ( m_inputTrackerHitsCollNames.size() != m_outputTrackerHitsCollNames.size() ||
      ( !noSimColl && m_inputTrackerSimHitsCollNames.size() != m_outputTrackerSimHitsCollNames.size() ) ||
      ( !noRelColl && m_inputTrackerHitRelNames.size() != m_outputTrackerHitRelNames.size() )  ) {

    std::stringstream err_msg;
    err_msg << "Mismatch between input and output collections"
	    << std::endl ;

    throw EVENT::Exception( err_msg.str() ) ;
  }

  if ( !noSimColl && ( m_inputTrackerHitsCollNames.size() != m_inputTrackerSimHitsCollNames.size() ||
                      m_inputTrackerHitsCollNames.size() != m_inputTrackerHitRelNames.size() ) ){

    std::stringstream err_msg;
    err_msg << "Mismatch between the reco and sim hits input collections"
	    << std::endl ;

    throw EVENT::Exception( err_msg.str() ) ;

  }

  if ( !noSimColl && ( m_outputTrackerHitsCollNames.size() != m_outputTrackerSimHitsCollNames.size() ||
                      m_outputTrackerHitsCollNames.size() != m_outputTrackerHitRelNames.size() ) ){

    std::stringstream err_msg;
    err_msg << "Mismatch between the reco and sim hits output collections"
	    << std::endl ;

    throw EVENT::Exception( err_msg.str() ) ;

  }
  

  // --- Get the JetCalo collection:

  LCCollection* m_inputJetCalo = NULL;
  try {
    m_inputJetCalo = evt->getCollection( m_inputJetCaloCollName );
  }
  catch( lcio::DataNotAvailableException& e ) {
    streamlog_out(WARNING) << m_inputJetCaloCollName << " collection not available" << std::endl;
    return;
  }


  // --- Get the input hit collections and create the corresponding output collections:

  const unsigned int nTrackerHitCol = m_inputTrackerHitsCollNames.size();
  std::vector<LCCollection*> inputHitColls(nTrackerHitCol);
  std::vector<LCCollection*> inputSimHitColls(nTrackerHitCol);
  std::vector<LCCollection*> inputHitRels(nTrackerHitCol);

  std::vector<LCCollectionVec*> outputTrackerHitColls(nTrackerHitCol);
  std::vector<LCCollectionVec*> outputTrackerSimHitColls(nTrackerHitCol);
  std::vector<LCCollectionVec*> outputTrackerHitRels(nTrackerHitCol);

  for (unsigned int icol=0; icol<nTrackerHitCol ; ++icol) {

    // get the reco hits
    try {
      inputHitColls[icol] = evt->getCollection(m_inputTrackerHitsCollNames[icol]);
    }
    catch( lcio::DataNotAvailableException& e ) {
      streamlog_out(WARNING) << m_inputTrackerHitsCollNames[icol]
			     << " collection not available" << std::endl;
      continue;
    }

    // get the sim hits
    if(!noSimColl){

      try {
        inputSimHitColls[icol] = evt->getCollection(m_inputTrackerSimHitsCollNames[icol]);
      }
      catch( lcio::DataNotAvailableException& e ) {
        streamlog_out(WARNING) << m_inputTrackerSimHitsCollNames[icol]
			     << " collection not available" << std::endl;
        continue;
      }

      try {
        inputHitRels[icol] = evt->getCollection(m_inputTrackerHitRelNames[icol]);
      }
      catch( lcio::DataNotAvailableException& e ) {
        streamlog_out(WARNING) << m_inputTrackerHitRelNames[icol]
			     << " collection not available" << std::endl;
        continue;
      }
    } 

    // reco hit output collections
    std::string encoderString = inputHitColls[icol]->getParameters().getStringVal( "CellIDEncoding" );
    outputTrackerHitColls[icol] = new LCCollectionVec( inputHitColls[icol]->getTypeName() );
    outputTrackerHitColls[icol]->parameters().setValue( "CellIDEncoding", encoderString );
    LCFlagImpl lcFlag(inputHitColls[icol]->getFlag());
    outputTrackerHitColls[icol]->setFlag(lcFlag.getFlag());
    
    // sim hit output collections
    if(!noSimColl){
      outputTrackerSimHitColls[icol] = new LCCollectionVec( inputSimHitColls[icol]->getTypeName() );
      outputTrackerSimHitColls[icol]->parameters().setValue( "CellIDEncoding", encoderString );
      LCFlagImpl lcFlag_sim(inputSimHitColls[icol]->getFlag());
      outputTrackerSimHitColls[icol]->setFlag(lcFlag_sim.getFlag());
    
      outputTrackerHitRels[icol] = new LCCollectionVec( inputHitRels[icol]->getTypeName() );
      LCFlagImpl lcFlag_rel(inputHitRels[icol]->getFlag()) ;
      outputTrackerHitRels[icol]->setFlag( lcFlag_rel.getFlag() ) ;
    }
    
  }


  // --- Loop over the JetCalo:
  
  std::vector<std::set<int> > hits_to_save(nTrackerHitCol);

  LCCollectionVec* filterJetsColl = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
  LCFlagImpl flagJets( m_inputJetCalo->getFlag() );
  filterJetsColl->setFlag( flagJets.getFlag() );

  unsigned int nGoodJets = 0;
  
  for (int ipart=0; ipart<m_inputJetCalo->getNumberOfElements(); ++ipart){

    ReconstructedParticle* part = dynamic_cast<ReconstructedParticle*>( m_inputJetCalo->getElementAt(ipart) );

    if( !FilterJetConeHits::filterJetBib(part) ) continue; //quality selection to reject fake jets
    nGoodJets++;

    if( m_createFilteredJets ) saveJet( part , filterJetsColl); //save filtered jet

    double mom[3];
    if( !m_makeDirCorrection ){
      for(int n = 0; n < 3; n++) mom[n] = part->getMomentum()[n]; //> to avoid const cast
    }
    else{
      directionCorrection( part->getMomentum(), mom );
    }

    // --- Loop over the tracker hits and select hits inside a cone around the jet axis:

    for (unsigned int icol=0; icol<inputHitColls.size(); ++icol){

      LCCollection* hit_col  =  inputHitColls[icol];
      if( !hit_col ) continue ;
      
      for (int ihit=0; ihit<hit_col->getNumberOfElements(); ++ihit){

	      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(hit_col->getElementAt(ihit));

	      // --- Skip hits that are in the opposite hemisphere w.r.t. the jet axis:
        if ( ( hit->getPosition()[0]*mom[0] +
              hit->getPosition()[1]*mom[1] +
              hit->getPosition()[2]*mom[2] ) < 0. ) continue;


        // --- Get the distance between the hit and the jet axis
        
        double jet_p = sqrt( pow(mom[0],2) + pow(mom[1],2) + pow(mom[2],2)  );
        double jet_theta = acos(mom[2]/jet_p);
        double jet_eta = -std::log(tan(jet_theta/2.));
        
        double hit_d = sqrt( pow(hit->getPosition()[0],2) + pow(hit->getPosition()[1],2) + pow(hit->getPosition()[2],2)  );
        double hit_theta = acos(hit->getPosition()[2]/hit_d);
        double hit_eta = -std::log(tan(hit_theta/2.));

        double jet_pxy = sqrt( pow(mom[0],2) + pow(mom[1],2) );
        double hit_dxy = sqrt( pow(hit->getPosition()[0],2) + pow(hit->getPosition()[1],2) );
        double deltaPhi = acos( (mom[0]*hit->getPosition()[0] + mom[1]*hit->getPosition()[1] )/jet_pxy/hit_dxy);

        double deltaR = sqrt( pow(deltaPhi,2) + pow(jet_eta-hit_eta,2)  );

        if ( m_fillHistos ){

          m_dist->Fill(deltaR);

        }
	
        if ( deltaR < m_deltaRCut )
          hits_to_save[icol].insert(ihit);
	
      } // ihit loop
       
    } // icol loop

  } // ipart loop


  // --- Add the filtered hits to the output collections:
  
  for (unsigned int icol=0; icol<inputHitColls.size(); ++icol){

    for ( auto& ihit: hits_to_save[icol] ){

      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(inputHitColls[icol]->getElementAt(ihit));
      TrackerHitPlaneImpl* hit_new = new TrackerHitPlaneImpl();

      hit_new->setCellID0(hit->getCellID0());
      hit_new->setCellID1(hit->getCellID1());
      hit_new->setType(hit->getType());
      hit_new->setPosition(hit->getPosition());
      hit_new->setU(hit->getU());         
      hit_new->setV(hit->getV());     
      hit_new->setdU(hit->getdU());
      hit_new->setdV(hit->getdV());
      hit_new->setEDep(hit->getEDep());      
      hit_new->setEDepError(hit->getEDepError()); 
      hit_new->setTime(hit->getTime());
      hit_new->setQuality(hit->getQuality());   

      outputTrackerHitColls[icol]->addElement( hit_new );

      if( noSimColl ) continue;

      LCRelation* rel = dynamic_cast<LCRelation*>(inputHitRels[icol]->getElementAt(ihit));

      
      SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(rel->getTo());
      SimTrackerHitImpl* simhit_new = new SimTrackerHitImpl();

      simhit_new->setCellID0(simhit->getCellID0());
      simhit_new->setCellID1(simhit->getCellID1());
      simhit_new->setPosition(simhit->getPosition());
      simhit_new->setEDep(simhit->getEDep());
      simhit_new->setTime(simhit->getTime());
      simhit_new->setMCParticle(simhit->getMCParticle());
      simhit_new->setMomentum(simhit->getMomentum());
      simhit_new->setPathLength(simhit->getPathLength());
      simhit_new->setQuality(simhit->getQuality());
      simhit_new->setOverlay(simhit->isOverlay());
      simhit_new->setProducedBySecondary(simhit->isProducedBySecondary());

      outputTrackerSimHitColls[icol]->addElement( simhit_new );

      
      LCRelationImpl* rel_new = new LCRelationImpl();
      
      rel_new->setFrom(hit_new);
      rel_new->setTo(simhit_new);
      rel_new->setWeight(rel->getWeight());

      outputTrackerHitRels[icol]->addElement( rel_new );

    } // ihit loop

    streamlog_out( MESSAGE ) << " " << hits_to_save[icol].size() << " hits added to the collections: "
			     << m_outputTrackerHitsCollNames[icol] << std::endl;

    evt->addCollection( outputTrackerHitColls[icol], m_outputTrackerHitsCollNames[icol] ) ;
    if(!noSimColl){
      evt->addCollection( outputTrackerSimHitColls[icol], m_outputTrackerSimHitsCollNames[icol] ) ;
      evt->addCollection( outputTrackerHitRels[icol], m_outputTrackerHitRelNames[icol] ) ;
    }

    streamlog_out( DEBUG5 ) << " output collection " << m_outputTrackerHitsCollNames[icol] << " of type "
			    << outputTrackerHitColls[icol]->getTypeName() << " added to the event" << std::endl;
		if(!noSimColl) streamlog_out( DEBUG5 )  << " output collection " << m_outputTrackerSimHitsCollNames[icol] << " of type "
			    << outputTrackerSimHitColls[icol]->getTypeName() << " added to the event \n"
			    << " output collection " << m_outputTrackerHitRelNames[icol] << " of type "
			    << outputTrackerHitRels[icol]->getTypeName() << " added to the event  "
			    << std::endl ;

  } // icol loop

  if( m_createFilteredJets ){
    evt->addCollection(filterJetsColl, m_filteredJetCaloCollName);
    streamlog_out(MESSAGE) << "Saved Filtered Jets collection, with " << filterJetsColl->getNumberOfElements() << " jets." << std::endl;
  }
  else{
    delete filterJetsColl;
  }

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() << std::endl ;
  
  _nEvt ++ ;

}



void FilterJetConeHits::check( LCEvent * evt ) { 
}



void FilterJetConeHits::end(){ 

  std::cout << "FilterJetConeHits::end()  " << name() 
   	    << " processed " << _nEvt << " events in " << _nRun << " runs "
     	    << std::endl ;

}

// jet quality selection
bool FilterJetConeHits::filterJetBib(ReconstructedParticle* jet){

  auto daughters = jet->getParticles();

  double maxPt = 0., pt = 0.;
  int ntracks = 0;

  for(unsigned int i = 0; i < daughters.size(); i++){
    
    if(daughters[i]->getCharge() == 0) continue;  //only care about tracks
    ntracks++;
    pt = sqrt( pow( daughters[i]->getMomentum()[0] ,2) + pow( daughters[i]->getMomentum()[1] ,2) );
    if(pt > maxPt) maxPt = pt;
  }

  if(ntracks < m_minNTracks || maxPt < m_minDaughterMaxPt) return false; //check quality of jet, if bad reject

  streamlog_out(DEBUG) << "Accepted jet with ntracks = " << ntracks << " and max daughter pT = " << maxPt << std::endl;
  return true;
}

// save jet in collection
void FilterJetConeHits::saveJet( ReconstructedParticle* jet, LCCollectionVec* jetsColl ){

  // this way only a pointer to the original jet is saved
  // if necessary implement here the phys copy save
  //jetsColl->addElement(jet);
  ReconstructedParticleImpl* j = new ReconstructedParticleImpl();
  j->setType( jet->getType() );

  double pcorr[3];
  if( !m_makeDirCorrection ) j->setMomentum( jet->getMomentum() );
  else{
    directionCorrection(jet->getMomentum(), pcorr);
    j->setMomentum( pcorr );
  }

  j->setEnergy( jet->getEnergy() );
  j->setMass( jet->getMass() );
  j->setCharge( jet->getCharge() );
  j->setReferencePoint( jet->getReferencePoint() );

  const EVENT::ReconstructedParticleVec rpvec = jet->getParticles();
  for(unsigned int i = 0; i < rpvec.size(); i++) j->addParticle(rpvec[i]);

  const EVENT::ClusterVec clvec = jet->getClusters();
  for(unsigned int i = 0; i < clvec.size(); i++) j->addCluster(clvec[i]);

  const EVENT::TrackVec trvec = jet->getTracks();
  for(unsigned int i = 0; i < trvec.size(); i++) j->addTrack(trvec[i]);

  jetsColl->addElement(j);

  streamlog_out( MESSAGE ) << "Saving Jet: p = ( " << j->getMomentum()[0] << " , "
                                                    << j->getMomentum()[1] << " , "
                                                    << j->getMomentum()[2] << " ) " << std::endl;

  if( m_makeDirCorrection ){
    const double* pold = jet->getMomentum();
    double oldtheta = acos( pold[2] / sqrt(pold[0]*pold[0] + pold[1]*pold[1] + pold[2]*pold[2]) ) * 180. / 3.14159265;
    const double* pnew = j->getMomentum();
    double newtheta = acos( pnew[2] / sqrt(pnew[0]*pnew[0] + pnew[1]*pnew[1] + pnew[2]*pnew[2]) ) * 180. / 3.14159265;
    streamlog_out( MESSAGE ) << "Corrected: old THETA = " << oldtheta << " => new THETA = " << newtheta << std::endl;
  }

  return;
}

void FilterJetConeHits::directionCorrection( const double* p, double* pcorr ){

  double ptot = sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
  double theta = acos( p[2] / ptot ) * 180. / 3.14159265;
  double theta_pre = theta;

  //> only in transition region, hardcoded for now...
  if(theta > 30. && theta < 60.){
    theta = m_corrConst + m_corrLin * theta + m_corrQuad * pow(theta, 2.) + m_corrCub * pow(theta, 3.);
  }
  else if(theta > 120. && theta < 150.){
    theta = 180. - theta;
    theta = m_corrConst + m_corrLin * theta + m_corrQuad * pow(theta, 2.) + m_corrCub * pow(theta, 3.);
    theta = 180. - theta;
  }

  pcorr[0] = p[0] * sin( theta * 3.14159265 / 180. ) / sin( theta_pre * 3.14159265 / 180. );
  pcorr[1] = p[1] * sin( theta * 3.14159265 / 180. ) / sin( theta_pre * 3.14159265 / 180. );
  pcorr[2] = p[2] * cos( theta * 3.14159265 / 180. ) / cos( theta_pre * 3.14159265 / 180. );

  return;
}
  