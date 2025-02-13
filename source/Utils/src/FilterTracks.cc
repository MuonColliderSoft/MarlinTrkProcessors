# include "FilterTracks.h"

#include <math.h>

#include <DD4hep/Detector.h>

#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/TrackImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTrackerConf.h>

#include "TMVA/Reader.h"

#include <iostream> 

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

  registerProcessorParameter("MaxPt",
             "Max transverse momentum",
           _MaxPt,
           _MaxPt
            );

  registerProcessorParameter("MinTheta",
             "Minimum theta",
           _MinTheta,
           _MinTheta
            );

  registerProcessorParameter("MaxTheta",
             "Max theta",
           _MaxTheta,
           _MaxTheta
            );

  registerProcessorParameter("Chi2Spatial",
		  	     "iMinimum value for Spatial chi squared",
			     _Chi2Spatial,
			     _Chi2Spatial
		 	      );

  registerProcessorParameter("MinNdf",
             "Minimum value for ndf",
           _MinNdf,
           _MinNdf
            );

  registerProcessorParameter("MaxHoles",
             "Max number of holes",
           _MaxHoles,
           _MaxHoles
            );

  registerProcessorParameter("MaxOutliers",
             "Max number of outliers",
           _MaxOutl,
           _MaxOutl
          );

  registerProcessorParameter("MaxOutliersOverHits",
            "Max ratio of outliers/hits",
            _MaxOutlOverHits,
            _MaxOutlOverHits
         );

  registerProcessorParameter("Bz",
             "Magnetic field in Tesla (default: 5)",
           _Bz,
           _Bz
            );

  registerProcessorParameter("NNmethod",
           "Name of the NN method, if empty uses standard cuts",
           _NNmethod,
           std::string("")
            );

  registerProcessorParameter("NNweights",
             "File xml with the weights of the NN",
           _NNweights,
           std::string("")
            );

  registerProcessorParameter("NNvars",
             "Sorted list with the names of the variables used in the training"
              "Supported variables are: trtvhn, trtihn, trtohn, trthn, trtnh, trch2, trndf",
           _NNvars,
           _NNvars
            );

  registerProcessorParameter("NNthr",
             "NN threshold (Takes tracks with prediction > NNthr)",
           _NNthr,
           _NNthr
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
  // it requires that geometry has been already initialized
  // buildBfield() ;
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
  // do not store pointers but real tracks
  OutputTrackCollection->setSubset(false);

  // if we want to point back to the hits we need to set the flag
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  OutputTrackCollection->setFlag(trkFlag.getFlag());

  TMVA::Reader* reader = new TMVA::Reader();

  // Get input collection
  LCCollection* InputTrackCollection = evt->getCollection(_InputTrackCollection);
  
  if( InputTrackCollection->getTypeName() != lcio::LCIO::TRACK )
    { throw EVENT::Exception( "Invalid collection type: " + InputTrackCollection->getTypeName() ) ; }

  // Filter
  std::string encoderString = lcio::LCTrackerCellID::encoding_string();
  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(encoderString);

  std::unordered_map<std::string, float> vars = {
    {"trtvhn", 0},
    {"trtihn", 0},
    {"trtohn", 0},
    {"trthn", 0},
    {"trtnh", 0},
    {"trch2", 0},
    {"trndf", 0},
  };

  if ( not _NNmethod.empty() ) {

    for (unsigned i=0,N=_NNvars.size(); i<N; ++i) {
      std::string name = _NNvars[i];
      if (vars.find(name) == vars.end())
        throw EVENT::Exception( "NN variable not supported: " + name ) ;  
      reader->AddVariable( name, &vars[name] );
    }
    reader->BookMVA(_NNmethod, _NNweights);
  }

  for(int i=0; i<InputTrackCollection->getNumberOfElements(); i++) {
    EVENT::Track *trk=dynamic_cast<EVENT::Track*>(InputTrackCollection->getElementAt(i));

    vars["trtnh"]  = trk->getTrackerHits().size();

    const EVENT::IntVec& subdetectorHits = trk->getSubdetectorHitNumbers();
    vars["trtvhn"] = subdetectorHits[1]+subdetectorHits[2];
    vars["trtihn"] = subdetectorHits[3]+subdetectorHits[4];
    vars["trtohn"] = subdetectorHits[5]+subdetectorHits[6];

    vars["trthn"] = trk->getNholes();
    float pt = fabs(0.3*_Bz/trk->getOmega()/1000);
    float theta = M_PI_2-atan(trk->getTanLambda());
    vars["trch2"] = trk->getChi2();

    vars["trndf"] = trk->getNdf();

    if(_BarrelOnly == true) {
      bool endcaphits = false;
      for(int j=0; j<vars["trtnh"]; ++j) {
	      //Find what subdetector the hit is on 
	      uint32_t systemID = decoder(trk->getTrackerHits()[j])["system"];
	      if(systemID == 2 || systemID == 4 || systemID == 6) {
	        endcaphits = true;
	        break;
	      }
      }
      if (endcaphits == false) { 
        auto itrk = dynamic_cast<IMPL::TrackImpl*>(trk);
        OutputTrackCollection->addElement(new IMPL::TrackImpl(*itrk)); 
      }
    } else { // track property cuts
      if ( not _NNmethod.empty() ) { // NN cuts
        auto prediction = reader->EvaluateMVA(_NNmethod);
        if ( prediction > _NNthr ) {
          auto itrk = dynamic_cast<IMPL::TrackImpl*>(trk);
          OutputTrackCollection->addElement(new IMPL::TrackImpl(*itrk)); 
        }
      } else { // user cuts
        if (vars["trtnh"]  > _NHitsTotal  &&
	          vars["trtvhn"] > _NHitsVertex &&
	          vars["trtihn"] > _NHitsInner  &&
	          vars["trtohn"] > _NHitsOuter  &&
	          pt             > _MinPt       &&
            pt             < _MaxPt       &&
            theta          > _MinTheta    &&
            theta          < _MaxTheta    &&
	          vars["trch2"]  > _Chi2Spatial &&
            vars["trndf"]  > _MinNdf      &&
            vars["trtnh"]-vars["trndf"]/2 < _MaxOutl &&
            (vars["trtnh"]-vars["trndf"]/2) / vars["trtnh"] < _MaxOutlOverHits &&
            vars["trthn"]  < _MaxHoles)
	          { 
              auto itrk = dynamic_cast<IMPL::TrackImpl*>(trk);
              OutputTrackCollection->addElement(new IMPL::TrackImpl(*itrk)); 
            }
      }
    }
  }

  // Save output track collection
  evt->addCollection(OutputTrackCollection, _OutputTrackCollection);  
  delete reader;
}

void FilterTracks::end()
{  }
