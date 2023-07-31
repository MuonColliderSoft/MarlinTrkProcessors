# include "TrackPerf/FilterClusters.hxx"

#include <math.h>

#include <DD4hep/Detector.h>

#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/LCTrackerConf.h>
#include <IMPL/LCRelationImpl.h>

#include <marlin/AIDAProcessor.h>

#include <AIDA/ITree.h>

#include "marlin/VerbosityLevels.h"


FilterClusters aFilterClusters ;

FilterClusters::FilterClusters()
  : Processor("FilterClusters")
{
  // modify processor description
  _description = "FilterClusters processor filters a collection of tracker hits based on cluster size and outputs a filtered collection";

  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("ThetaRanges",
         "Divide theta into bins for different cluster size cuts",
         _ThetaRanges, 
         {}
          );
  registerProcessorParameter("ClusterSize",
		  	 "Maximum cluster size for each theta range",
			   _ClusterSize,
			   {}
			    );  
  registerProcessorParameter("Layers",
		  	 "Layers to be filtered",
			   _Layers,
			   {}
			    );  
  registerInputCollection( LCIO::TRACKERHIT,
		  	 "InTrackerHitCollection" ,
			   "Name of the input tracker hit collection",
			   _InTrackerHitCollection,
		     _InTrackerHitCollection
		 	    );

  registerInputCollection( LCIO::LCRELATION,
		     "InRelationCollection" ,
			   "Name of the input relation collection",
			   _InRelationCollection,
		     _InRelationCollection
		 	    );

  registerOutputCollection( LCIO::TRACKERHIT,
		  	 "OutTrackerHitCollection" ,
			   "Name of output tracker hit collection",
			   _OutTrackerHitCollection,
			   std::string("FilteredVBTrackerHits")
			    );

    registerOutputCollection( LCIO::LCRELATION,
		  	 "OutRelationCollection" ,
			   "Name of output relation collection",
			   _OutRelationCollection,
			   std::string("FilteredVBTrackerHitsRelations")
			    );

}

void FilterClusters::init()
{
  // Print the initial parameters
  printParameters() ;
}

void FilterClusters::processRunHeader( LCRunHeader* /*run*/)
{ }

void FilterClusters::processEvent( LCEvent * evt )
{
  // Make the output track collection
  LCCollectionVec *OutTrackerHitCollection = new LCCollectionVec(LCIO::TRACKERHIT);
  OutTrackerHitCollection->setSubset(true);
  LCCollectionVec *OutRelationCollection   = new LCCollectionVec(LCIO::LCRELATION);
  OutRelationCollection->setSubset(true);

  // Get input collection
  LCCollection* InTrackerHitCollection  = evt->getCollection(_InTrackerHitCollection);
  LCCollection* InRelationCollection    = evt->getCollection(_InRelationCollection);

  if( InTrackerHitCollection->getTypeName() != lcio::LCIO::TRACKERHITPLANE )
    { throw EVENT::Exception( "Invalid collection type: " + InTrackerHitCollection->getTypeName() ) ; }
    streamlog_out(DEBUG0) << "Wrong collection type for TrackerHitCollection. \n";
  if( InRelationCollection->getTypeName() != lcio::LCIO::LCRELATION )
    { throw EVENT::Exception( "Invalid collection type: " + InRelationCollection->getTypeName() ) ; }
    streamlog_out(DEBUG0) << "Wrong collection type for InRelationCollection. \n";


  streamlog_out(DEBUG0) << "Number of Elements in VB Tracker Hits Collection: " << InTrackerHitCollection->getNumberOfElements() <<std::endl;
  // Filter
  for(int i=0; i<InTrackerHitCollection->getNumberOfElements(); ++i) //loop through all hits
    {
      streamlog_out(DEBUG0) << "Loop over hits opened. \n";
      EVENT::TrackerHit *trkhit=static_cast<EVENT::TrackerHit*>(InTrackerHitCollection->getElementAt(i)); //define trkhit var, pointer to i'th element of tracker hits
      EVENT::LCRelation *rel=static_cast<EVENT::LCRelation*>(InRelationCollection->getElementAt(i));
      
      //Calculating theta 
      float x = trkhit->getPosition()[0];
      float y = trkhit->getPosition()[1];
      float z = trkhit->getPosition()[2];
      float r = sqrt(pow(x,2)+pow(y,2));
      float incidentTheta = std::atan(r/z); 

      if(incidentTheta<0)
        incidentTheta += M_PI;

      //Calculating cluster size
      const lcio::LCObjectVec &rawHits = trkhit->getRawHits();
      float max = -1000000;
      float min = 1000000;
      for (size_t j=0; j<rawHits.size(); ++j) {
        lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit*>( rawHits[j] );
        const double *localPos = hitConstituent->getPosition();
        x = localPos[0];
        y = localPos[1];
        if (y < min){
          min = y;
        }
        else if (y > max){
          max = y;          
        }
      }
      float cluster_size = (max - min)+1;

      //Get hit subdetector/layer 
      std::string _encoderString = lcio::LCTrackerCellID::encoding_string();
      UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
      uint32_t systemID = decoder(trkhit)["system"];
      uint32_t layerID = decoder(trkhit)["layer"];
      bool filter_layer = false;
      for (int i=0; i<_Layers.size(); ++i){
        if (layerID == std::stof(_Layers[i])) {
          filter_layer = true;
        }
      }
      streamlog_out(DEBUG0) << "Filter layer: " << filter_layer << std::endl;
      
      for (int i=0; i<_ThetaRanges.size()-1; ++i) {
        streamlog_out( DEBUG0 ) << "theta: " << incidentTheta << std::endl;
        float min = std::stof(_ThetaRanges[i]);
        float max = std::stof(_ThetaRanges[i+1]);
        streamlog_out( DEBUG0 ) << "theta range: " << min << ", " << max << std::endl;

        if(incidentTheta > min and incidentTheta <= max and not filter_layer){
          streamlog_out( DEBUG0 ) << "theta in range" << std::endl;
          streamlog_out( DEBUG0 ) << "cluster size cut off: " << _ClusterSize[i] << std::endl;
          streamlog_out( DEBUG0 ) << "cluster size: " << cluster_size << std::endl;
          if(cluster_size < std::stof(_ClusterSize[i])) {
            streamlog_out( DEBUG0 ) << "cluster added" << std::endl;
            OutTrackerHitCollection->addElement(trkhit); 
            OutRelationCollection->addElement(rel); }
          else {
            streamlog_out( DEBUG0 ) << "cluster rejected" << std::endl;
          }
        }
        else{
          streamlog_out( DEBUG0 ) << "theta out of range or layer filtered" << std::endl;
        }
      }
    }

  // Save output track collection
  evt->addCollection(OutTrackerHitCollection, _OutTrackerHitCollection); 
  evt->addCollection(OutRelationCollection, _OutRelationCollection); 
}

void FilterClusters::end()
{ }