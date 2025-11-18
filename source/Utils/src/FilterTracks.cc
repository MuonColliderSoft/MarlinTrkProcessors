#include "FilterTracks.h"

#include <math.h>

#include <DD4hep/Detector.h>

#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTrackerConf.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCRelationImpl.h>

#include <UTIL/LCRelationNavigator.h>

FilterTracks aFilterTracks;

FilterTracks::FilterTracks() : Processor("FilterTracks") {
  // modify processor description
  _description = "FilterTracks processor filters a collection of tracks based on NHits and MinPt and outputs a "
                 "filtered collection";

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection(LCIO::TRACK, "InputTrackCollectionName", "Name of the input collection",
                          _InputTrackCollection, _InputTrackCollection);

  registerOutputCollection(LCIO::TRACK, "OutputTrackCollectionName", "Name of output collection",
                           _OutputTrackCollection, std::string("FilteredTracks"));

  registerInputCollection(LCIO::LCRELATION, "InputRelationCollectionName",
                          "Name of the input track to MCParticle relation collection", _InputTrackRelationCollection,
                          _InputTrackRelationCollection);

  registerOutputCollection(LCIO::LCRELATION, "OutputRelationCollectionName",
                           "Refit Track to MCParticle relation collection Name", _OutputTrackRelationCollection,
                           _OutputTrackRelationCollection);

  registerProcessorParameter("BarrelOnly", "If true, just keep tracks with only barrel hits", _BarrelOnly, _BarrelOnly);

  registerProcessorParameter("HasCaloState",
                             "If true, just keep tracks that have a TrackState at the Calorimeter surface",
                             _HasCaloState, _HasCaloState);

  registerProcessorParameter("NHitsTotal", "Minimum number of hits on track", _NHitsTotal, _NHitsTotal);

  registerProcessorParameter("NHitsVertex", "Minimum number of hits on vertex detector", _NHitsVertex, _NHitsVertex);

  registerProcessorParameter("NHitsInner", "Minimum number of hits on inner tracker", _NHitsInner, _NHitsInner);

  registerProcessorParameter("NHitsOuter", "Minimum number of hits on outer tracker", _NHitsOuter, _NHitsOuter);

  registerProcessorParameter("MaxHoles", "Maximum number of holes on track", _MaxHoles, _MaxHoles);

  registerProcessorParameter("MinPt", "Minimum transverse momentum", _MinPt, _MinPt);

  registerProcessorParameter("MaxD0", "Max d0", _MaxD0, _MaxD0);

  registerProcessorParameter("MaxZ0", "Max z0", _MaxZ0, _MaxZ0);

  registerProcessorParameter("Chi2Spatial", "Spatial chi squared", _Chi2Spatial, _Chi2Spatial);
}

void FilterTracks::init() {
  // Print the initial parameters
  printParameters();
  buildBfield();
}

void FilterTracks::processRunHeader(LCRunHeader* /*run*/) {}

void FilterTracks::buildBfield() {
  // Get the magnetic field
  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
  const double position[3] = {0, 0, 0};      // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3] = {0, 0, 0}; // initialise object to hold magnetic field
  lcdd.field().magneticField(position,
                             magneticFieldVector); // get the magnetic field vector from DD4hep
  _Bz = magneticFieldVector[2] / dd4hep::tesla;
}

void FilterTracks::processEvent(LCEvent* evt) {
  // Make the output track collection
  LCCollectionVec* OutputTrackCollection = new LCCollectionVec(LCIO::TRACK);
  OutputTrackCollection->setSubset(true);

  // track-mcpart relation output collections
  UTIL::LCRelationNavigator myNav = UTIL::LCRelationNavigator(LCIO::TRACK, LCIO::MCPARTICLE); 
  
  // Get input collection
  LCCollection* InputTrackCollection = evt->getCollection(_InputTrackCollection);

  if (InputTrackCollection->getTypeName() != lcio::LCIO::TRACK) {
    throw EVENT::Exception("Invalid collection type: " + InputTrackCollection->getTypeName());
  }

  // Get input relations
  LCCollection* InputTrackRelationCollection = nullptr;
  std::shared_ptr<LCRelationNavigator> relation;
  bool has_input_relations = true;
  try {
    InputTrackRelationCollection = evt->getCollection(_InputTrackRelationCollection);
    relation = std::make_shared<LCRelationNavigator>(InputTrackRelationCollection);
  } catch (lcio::DataNotAvailableException& e) {
    streamlog_out(WARNING) << "Relation collection: " << _InputTrackRelationCollection << " not available" << std::endl;
    has_input_relations = false;
  }

  // Filter
  std::string encoderString = lcio::LCTrackerCellID::encoding_string();
  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(encoderString);

  for (int i = 0; i < InputTrackCollection->getNumberOfElements(); i++) {
    EVENT::Track* trk = dynamic_cast<EVENT::Track*>(InputTrackCollection->getElementAt(i));

    int nhittotal = trk->getTrackerHits().size();

    const EVENT::IntVec& subdetectorHits = trk->getSubdetectorHitNumbers();
    int nhitvertex = subdetectorHits[1] + subdetectorHits[2];
    int nhitinner = subdetectorHits[3] + subdetectorHits[4];
    int nhitouter = subdetectorHits[5] + subdetectorHits[6];

    float pt = fabs(0.3 * _Bz / trk->getOmega() / 1000);

    float chi2spatial = trk->getChi2();

    int nholes = trk->getNholes();

    // Check if a TrackState at the calo surface exists
    const std::vector<EVENT::TrackState*>& trackStates = trk->getTrackStates();
    const auto foundCaloState = std::find_if(trackStates.begin(), trackStates.end(), [](const auto ts) {
                                  return ts->getLocation() == EVENT::TrackState::AtCalorimeter;
                                }) != trackStates.end();
    if (_HasCaloState && !foundCaloState) {
      streamlog_out(DEBUG) << "No calo state, skipping track!" << std::endl;
      continue;
    }

    if (_BarrelOnly == true) {
      bool endcaphits = false;
      for (int j = 0; j < nhittotal; ++j) {
        // Find what subdetector the hit is on
        uint32_t systemID = decoder(trk->getTrackerHits()[j])["system"];
        if (systemID == 2 || systemID == 4 || systemID == 6) {
          endcaphits = true;
          break;
        }
      }
      if (endcaphits == false) {
        OutputTrackCollection->addElement(trk);
      }
    } else { // track property cuts
      if (nhittotal < _NHitsTotal){
        streamlog_out(DEBUG) << "Only " << nhittotal << " hits, skipping track!" << std::endl;
        continue;
      }
      if (nhitvertex < _NHitsVertex){
        streamlog_out(DEBUG) << "Only " << nhitvertex << " vertex hits, skipping track!" << std::endl;
        continue;
      }
      if (nhitinner < _NHitsInner) {
        streamlog_out(DEBUG) << "Only " << nhitinner << " inner tracker hits, skipping track!" << std::endl;
        continue;
      }
      if (nhitouter < _NHitsOuter) {
        streamlog_out(DEBUG) << "Only " << nhitouter << " outer tracker hits, skipping track!" << std::endl;
        continue;
      }
      if (pt < _MinPt){
        streamlog_out(DEBUG) << "Pt = " << pt << " GeV, skipping track!" << std::endl;
        continue;
      }
      if (chi2spatial > _Chi2Spatial){
        streamlog_out(DEBUG) << "Chi2 spatial = " << chi2spatial << ", skipping track!" << std::endl;
        continue;
      }
      if (nholes > _MaxHoles){
        streamlog_out(DEBUG) << "Nholes = " << nholes << ", skipping track!" << std::endl;
        continue;
      }
      if (fabs(trk->getD0()) > _MaxD0){
        streamlog_out(DEBUG) << "D0 = " << trk->getD0() << ", skipping track!" << std::endl;
        continue;
      }
      if (fabs(trk->getZ0()) > _MaxZ0){
        streamlog_out(DEBUG) << "Z0 = " << trk->getZ0() << ", skipping track!" << std::endl;
        continue;
      }
      
      OutputTrackCollection->addElement(trk);
      
      // Add track-mcpart relations
      if (has_input_relations) {
        auto mcParticleVec = relation->getRelatedToObjects(trk);
        auto weightVec = relation->getRelatedToWeights(trk);
        for (size_t irel = 0; irel < mcParticleVec.size(); ++irel) {
          myNav.addRelation(trk, mcParticleVec[irel], weightVec[irel]);
        }
      }
      
    }
  }

  // Save output track collection
  evt->addCollection(OutputTrackCollection, _OutputTrackCollection);
  if (has_input_relations){
    LCCollection* outputTrackHitRel = myNav.createLCCollection();
    evt->addCollection(outputTrackHitRel, _OutputTrackRelationCollection);
  }
}

void FilterTracks::end() {}
