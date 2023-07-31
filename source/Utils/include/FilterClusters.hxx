#pragma once

#include <marlin/Processor.h>

//#include <ACTSTracking/GeometryIdMappingTool.hxx>
#include "/opt/ilcsoft/muonc/ACTSTracking/v1.1.0/ACTSTracking/GeometryIdMappingTool.hxx"

namespace TrackPerf
{
}

class FilterClusters : public marlin::Processor
{
public:
   virtual Processor* newProcessor() { return new FilterClusters ; }

   FilterClusters(const FilterClusters &) = delete ;
   FilterClusters& operator =(const FilterClusters &) = delete ;
   FilterClusters() ;

   /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
   virtual void init() ;

   /** Called for every run.
    */
   virtual void processRunHeader( LCRunHeader* run ) ;

   /** Called for every event - the working horse.
    */
   virtual void processEvent(LCEvent* evt) ;


   /** Called after data processing for clean up.
    */
   virtual void end() ;

private:
   //! Input track collection
   std::string _InTrackerHitCollection {};
   std::string _InRelationCollection {};

   //! Output track collection
   std::string _OutTrackerHitCollection {};
   std::string _OutRelationCollection {};

   //! Ranges for theta
   std::vector<std::string> _ThetaRanges;

   //! Cut-offs for cluster size in various theta ranges
   std::vector<std::string> _ClusterSize;

   //! Layers to be filtered
   std::vector<std::string> _Layers;

};