# include "FilterTracks.h"

#include <math.h>

#include <DD4hep/Detector.h>

#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTrackerConf.h>

FilterTracks aFilterTracks ;

FilterTracks::FilterTracks()
  : Processor("FilterTracks")
{
  // modify processor description
  _description = "FilterTracks processor filters a collection of tracks based on NHits and MinPt and outputs a filtered collection";

  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("BarrelOnly",
		  	     "If true, just keep tracks with only barrel hits",
			     _BarrelOnly,
			     _BarrelOnly
			      );  
  
  registerProcessorParameter("NHitsTotal",
		  	     "Minimum number of hits on track",
			     _NHitsTotal,
			     _NHitsTotal
			      );
  
  registerProcessorParameter("NHitsVertex",
		  	     "Minimum number of hits on vertex detector",
			     _NHitsVertex,
			     _NHitsVertex
			      );

  registerProcessorParameter("NHitsInner",
		  	     "Minimum number of hits on inner tracker",
			     _NHitsInner,
			     _NHitsInner
			      );

  registerProcessorParameter("NHitsOuter",
		  	     "Minimum number of hits on outer tracker",
			     _NHitsOuter,
			     _NHitsOuter
			      );

  registerProcessorParameter("MinPt",
		  	     "Minimum transverse momentum",
			     _MinPt,
			     _MinPt
		 	      );

  registerProcessorParameter("Chi2Spatial",
		  	     "Spatial chi squared",
			     _Chi2Spatial,
			     _Chi2Spatial
		 	      );

  registerInputCollection( LCIO::TRACK,
		  	   "InputTrackCollectionName" ,
			   "Name of the input collection",
			   _InputTrackCollection,
		     _InputTrackCollection
		 	    );

  registerOutputCollection( LCIO::TRACK,
		  	   "OutputTrackCollectionName" ,
			   "Name of output collection",
			   _OutputTrackCollection,
			   std::string("FilteredTracks")
			    );

}

void FilterTracks::init()
{
  // Print the initial parameters
  printParameters() ;
  buildBfield() ;
}

void FilterTracks::processRunHeader( LCRunHeader* /*run*/)
{ }

void FilterTracks::buildBfield() 
{
  // Get the magnetic field
  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
  const double position[3] = {
      0, 0,
      0};  // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3] = {
      0, 0, 0};  // initialise object to hold magnetic field
  lcdd.field().magneticField(
      position,
      magneticFieldVector);  // get the magnetic field vector from DD4hep
  _Bz = magneticFieldVector[2]/dd4hep::tesla;
}

void FilterTracks::processEvent( LCEvent * evt )
{
  // Make the output track collection
  LCCollectionVec *OutputTrackCollection = new LCCollectionVec(LCIO::TRACK);
  OutputTrackCollection->setSubset(true);

  // Get input collection
  LCCollection* InputTrackCollection  =evt->getCollection(_InputTrackCollection);

  if( InputTrackCollection->getTypeName() != lcio::LCIO::TRACK )
    { throw EVENT::Exception( "Invalid collection type: " + InputTrackCollection->getTypeName() ) ; }

  // Filter
  std::string encoderString = lcio::LCTrackerCellID::encoding_string();
  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(encoderString);

  for(int i=0; i<InputTrackCollection->getNumberOfElements(); i++) {
    EVENT::Track *trk=dynamic_cast<EVENT::Track*>(InputTrackCollection->getElementAt(i));

    int nhittotal  = trk->getTrackerHits().size();

    const EVENT::IntVec& subdetectorHits = trk->getSubdetectorHitNumbers();
    int nhitvertex = subdetectorHits[1]+subdetectorHits[2];
    int nhitinner = subdetectorHits[3]+subdetectorHits[4];
    int nhitouter = subdetectorHits[5]+subdetectorHits[6];

    float pt = fabs(0.3*_Bz/trk->getOmega()/1000);

    float chi2spatial = trk->getChi2();

    if(_BarrelOnly == true) {
      bool endcaphits = false;
      for(int j=0; j<nhittotal; ++j) {
	//Find what subdetector the hit is on 
	uint32_t systemID = decoder(trk->getTrackerHits()[j])["system"];
	if(systemID == 2 || systemID == 4 || systemID == 6) {
	  endcaphits = true;
	  break;
	}
      }
      if(endcaphits == false) { OutputTrackCollection->addElement(trk); }
    } else { // track property cuts
      if(nhittotal    > _NHitsTotal  &&
	 nhitvertex   > _NHitsVertex &&
	 nhitinner    > _NHitsInner  &&
	 nhitouter    > _NHitsOuter  &&
	 pt           > _MinPt       &&
	 chi2spatial  > _Chi2Spatial)
	{ OutputTrackCollection->addElement(trk); }
    }
  }

  // Save output track collection
  evt->addCollection(OutputTrackCollection, _OutputTrackCollection);  
}

void FilterTracks::end()
{ }
