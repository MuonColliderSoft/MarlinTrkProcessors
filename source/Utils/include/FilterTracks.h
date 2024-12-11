#pragma once

#include <marlin/Processor.h>

namespace TrackPerf
{
}

/** Utility processor that removes tracks based on different track properties.
 *
 *  @parameter InputTrackCollectionName Name of the input collection
 *  @parameter OutputTrackCollectionName Name of output collection
 *
 *  @parameter BarrelOnly Keep tracks with only barrel hits.
 *  @parameter NHitsTotal Minimum number of hits on track
 *  @parameter NHitsVertex Minimum number of hits on vertex detector
 *  @parameter NHitsInner Minimum number of hits on inner tracker
 *  @parameter NHitsOuter Minimum number of hits on outer tracker
 *  @parameter MaxOutliers Maximum number of outliers hits on track
 *  @parameter MaxHoles Maximum number of holes on track
 *  @parameter MinNdf Minimum value for ndf
 *  @parameter MinPt Minimum transverse momentum
 *  @parameter MaxPt Max transverse momentum
 *  @parameter MinTheta Minimum theta
 *  @parameter MaxTheta Max theta
 *  @parameter Chi2Spatial Spatial chi squared
 *
 * @author N. Bruhwiler
 * @author K. Krizka
 * @date  7 December 2022
 * @version $Id: $
 */
class FilterTracks : public marlin::Processor
{
public:
   virtual Processor* newProcessor() { return new FilterTracks ; }

   FilterTracks(const FilterTracks &) = delete ;
   FilterTracks& operator =(const FilterTracks &) = delete ;
   FilterTracks() ;

   /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
   virtual void init() ;

   /** Called for every run.
    */
   virtual void processRunHeader( LCRunHeader* run ) ;

   virtual void buildBfield() ;

   /** Called for every event - the working horse.
    */
   virtual void processEvent(LCEvent* evt) ;


   /** Called after data processing for clean up.
    */
   virtual void end() ;

private:
   //! Input track collection
   std::string _InputTrackCollection {};

   //! Output track collection
   std::string _OutputTrackCollection {};

   bool _BarrelOnly = false;

   //! Cut off for total number of hits
   int _NHitsTotal = 7;
   //! Cut off for number of hits in vertex detector (barrel and endcap combined)
   int _NHitsVertex = 3;
   //! Cut off for number of hits in inner tracker (barrel and endcap combined)
   int _NHitsInner = 2;
   //! Cut off for number of hits in outer tracker (barrel and endcap combined)
   int _NHitsOuter = 1;
   //! Cut off for number of holes
   int _MaxHoles = 10;

   //! Cut off for momentum (GeV)
   float _MinPt = 0.5;    //units GeV
   float _MaxPt = 1000.0; //units GeV

   //! Cut off for theta (rad)
   float _MinTheta = 0;    
   float _MaxTheta = 3.14; 

   //! Cut off for the value ndf
   int _MinNdf = 1;

   //! Cut off for outliers number   
   int _MaxOutl = 10;  

   //! Cut off for spatial and temporal chi squared values
   float _Chi2Spatial = 0;

   //! NN parameters
   std::string _NNmethod = ""; // if defined apply the NN
   std::string _NNweights = ""; // xml file with weights
   std::vector<std::string> _NNvars; // sorted list of variables used by NN
   float _NNthr = 0; // NN threshold   

   //! Default magnetic field value (Tesla)
   float _Bz = 5.0;
};

