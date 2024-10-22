#include "FilterClusters.h"

#include <math.h>

#include <DD4hep/Detector.h>

#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/CellIDDecoder.h>
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
  registerProcessorParameter("ThetaBins",
         "Number of bins in theta",
         _ThetaBins,
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


  streamlog_out(DEBUG0) << "Number of Elements in Tracker Hits Collection: " << InTrackerHitCollection->getNumberOfElements() <<std::endl;
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
      float ymax = -1000000;
      float xmax = -1000000;
      float ymin =  1000000;
      float xmin =  1000000;
      for (size_t j=0; j<rawHits.size(); ++j) {
        lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit*>( rawHits[j] );
        const double *localPos = hitConstituent->getPosition();
	float x_local = localPos[0];
	float y_local = localPos[1];

	if (y_local < ymin){
	  ymin = y_local;
	}
	if (y_local > ymax){
	  ymax = y_local;
	}
	
	if (x_local < xmin){
	  xmin = x_local;
	}
	if (x_local > xmax){
	  xmax = x_local;
	}
      }
      float cluster_size_y = (ymax - ymin)+1;
      float cluster_size_x = (xmax - xmin)+1;
      float cluster_size = rawHits.size();

      //Get hit subdetector/layer 
      std::string _encoderString = lcio::LCTrackerCellID::encoding_string();
      UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
      uint32_t systemID = decoder(trkhit)["system"];
      uint32_t layerID = decoder(trkhit)["layer"];
      bool filter_layer = false;

      int rows = _Layers.size(), cols = std::stoi(_ThetaBins);
      if((rows*cols != _ThetaRanges.size()) || (rows*(cols-1) != _ClusterSize.size())){
	std::cout<<"Either theta cuts or cluster cuts not provided for each layer. Please change the config, exiting now..."<<std::endl;
	return;
      }
      
      std::vector<std::vector<float>> _thetaCuts_byLayer;
      std::vector<std::vector<float>> _clusterSizeCuts_byLayer;

      for (int i = 0; i < rows; ++i) {
	std::vector<float> row;
	for (int j = 0; j < cols; ++j) {
	  row.push_back(std::stof(_ThetaRanges[j]));
	}
	_thetaCuts_byLayer.push_back(row);
      }

      for (int i = 0; i < rows; ++i) {
        std::vector<float> row;
        for (int j = 0; j < cols; ++j) {
          row.push_back(std::stof(_ClusterSize[j]));
	}
        _clusterSizeCuts_byLayer.push_back(row);
      }

	
      for (size_t j=0; j<_Layers.size(); ++j){
        if (layerID == std::stof(_Layers[j])) {
          filter_layer = true;
	  break;
        }
      }
      streamlog_out(DEBUG0) << "Filter layer: " << filter_layer << std::endl;

      for (size_t j=0; j<_thetaCuts_byLayer[layerID].size()-1; ++j) {
        streamlog_out( DEBUG0 ) << "theta: " << incidentTheta << std::endl;
	float min_theta = _thetaCuts_byLayer[layerID][j];
	float max_theta = _thetaCuts_byLayer[layerID][j+1];
        streamlog_out( DEBUG0 ) << "theta range: " << min_theta << ", " << max_theta << std::endl;
	
        if(incidentTheta >= min_theta and incidentTheta <= max_theta and filter_layer){
          streamlog_out( DEBUG0 ) << "theta in range" << std::endl;
	  streamlog_out( DEBUG0 ) << "cluster size cut off: " << _clusterSizeCuts_byLayer[layerID][j] << std::endl;
          streamlog_out( DEBUG0 ) << "cluster size: " << cluster_size << std::endl;
          if(cluster_size < _clusterSizeCuts_byLayer[layerID][j]) {
            streamlog_out( DEBUG0 ) << "cluster added" << std::endl;
            OutTrackerHitCollection->addElement(trkhit); 
            OutRelationCollection->addElement(rel); }
          else {
            streamlog_out( DEBUG0 ) << "cluster rejected" << std::endl;
          }
        }
        else{
          streamlog_out( DEBUG0 ) << "theta out of range or filtering not enabled for this layer" << std::endl;
        }
      }
    }

  // Save output track collection
  evt->addCollection(OutTrackerHitCollection, _OutTrackerHitCollection); 
  evt->addCollection(OutRelationCollection, _OutRelationCollection); 
}

void FilterClusters::end()
{ }
