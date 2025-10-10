/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "TruthTrackFinder.h"
#include "HitsSorterAndDebugger.h"

#include <marlinutil/HelixClass.h>

#include "MarlinTrk/Factory.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include "MarlinTrk/MarlinTrkDiagnostics.h"
#include "MarlinTrk/MarlinTrkUtils.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitSet32.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCRelationNavigator.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/SurfaceManager.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"

#include "CLHEP/Vector/TwoVector.h"

#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IHistogramFactory.h>

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace lcio;
using namespace marlin;
using namespace std;
using namespace dd4hep;
using namespace AIDA;

TruthTrackFinder aTruthTrackFinder;

/*

 This code performs a true pattern recognition by looping over all MC particles and adding all hits
 associated to them onto a track. This is then fitted and output.

*/

TruthTrackFinder::TruthTrackFinder() : Processor("TruthTrackFinder") {
  // modify processor description
  _description = "TruthTrackFinder builds tracks out of all hits associated to an MC particle";

  // Input collections - mc particles, tracker hits and the relationships between them

  std::vector<std::string> inputTrackerHitCollections;
  inputTrackerHitCollections.push_back(std::string("VXDTrackerHits"));

  registerInputCollections(LCIO::TRACKERHITPLANE, "TrackerHitCollectionNames",
                           "Name of the TrackerHit input collections", m_inputTrackerHitCollections,
                           inputTrackerHitCollections);

  std::vector<std::string> inputTrackerHitRelationCollections;
  inputTrackerHitRelationCollections.push_back(std::string("VXDTrackerHitRelations"));

  registerInputCollections(LCIO::LCRELATION, "SimTrackerHitRelCollectionNames",
                           "Name of TrackerHit SimTrackHit relation collections", m_inputTrackerHitRelationCollections,
                           inputTrackerHitRelationCollections);

  registerInputCollection(LCIO::MCPARTICLE, "MCParticleCollectionName", "Name of the MCParticle input collection",
                          m_inputParticleCollection, std::string("MCParticle"));

  // Output collections - tracks and relations
  registerOutputCollection(LCIO::TRACK, "SiTrackCollectionName", "Silicon track Collection Name",
                           m_outputTrackCollection, std::string("SiTracks"));

  registerOutputCollection(LCIO::LCRELATION, "SiTrackRelationCollectionName",
                           "Silicon track particle relation Collection Name", m_outputTrackRelationCollection,
                           std::string("SiTrackRelations"));

  // Flag parameters

  registerProcessorParameter(
      "UseTruthInPrefit",
      "If true use the truth information to initialise the helical prefit, otherwise use prefit by fitting 3 hits",
      m_useTruthInPrefit, bool(false));

  registerProcessorParameter("FitForward",
                             "If true fit 'forward' (go forward + smooth back adding last two hits with Kalman FIlter "
                             "steps), otherwise fit backward ",
                             m_fitForward, bool(false));
}

void TruthTrackFinder::init() {
  // Print the initial parameters
  printParameters();

  // Reset counters
  m_runNumber = 0;
  m_eventNumber = 0;
  m_fitFails = 0;

  // Set up the track fit factory
  trackFactory = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", nullptr, "");
  trackFactory->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useQMS, true);
  trackFactory->setOption(MarlinTrk::IMarlinTrkSystem::CFG::usedEdx, true);
  trackFactory->setOption(MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing, false);
  trackFactory->init();

  // Put default values for track fitting
  m_initialTrackError_d0 = 1.e6;
  m_initialTrackError_phi0 = 1.e2;
  m_initialTrackError_omega = 1.e-4;
  m_initialTrackError_z0 = 1.e6;
  m_initialTrackError_tanL = 1.e2;
  m_maxChi2perHit = 1.e2;

  // Get the magnetic field
  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
  const double position[3] = {0, 0, 0};      // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3] = {0, 0, 0}; // initialise object to hold magnetic field
  lcdd.field().magneticField(position, magneticFieldVector); // get the magnetic field vector from DD4hep
  m_magneticField = magneticFieldVector[2] / dd4hep::tesla;  // z component at (0,0,0)

  // Register this process
  Global::EVENTSEEDER->registerProcessor(this);

  // Initialize CellID encoder
  m_encoder = new UTIL::BitField64(lcio::LCTrackerCellID::encoding_string());
}

void TruthTrackFinder::processRunHeader(LCRunHeader*) { m_runNumber++; }

void TruthTrackFinder::processEvent(LCEvent* evt) {
  // set the correct configuration for the tracking system for this event
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useQMS> mson(trackFactory, true);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::usedEdx> elosson(trackFactory, true);
  MarlinTrk::TrkSysConfig<MarlinTrk::IMarlinTrkSystem::CFG::useSmoothing> smoothon(trackFactory, false);

  // Get the collection of MC particles
  LCCollection* particleCollection = 0;
  getCollection(particleCollection, m_inputParticleCollection, evt);
  if (particleCollection == 0)
    return;

  // Make objects to hold all of the tracker hit, simulated hit and relation collections
  std::vector<LCCollection*> trackerHitCollections;
  std::vector<LCCollection*> trackerHitRelationCollections;
  std::vector<LCRelationNavigator*> relations;

  // Loop over each input collection and get the data
  for (unsigned int collection = 0; collection < m_inputTrackerHitCollections.size(); collection++) {
    // Get the collection of tracker hits
    LCCollection* trackerHitCollection = 0;
    getCollection(trackerHitCollection, m_inputTrackerHitCollections[collection], evt);
    if (trackerHitCollection == 0)
      continue;
    trackerHitCollections.push_back(trackerHitCollection);

    // Get the collection of tracker hit relations
    LCCollection* trackerHitRelationCollection = 0;
    getCollection(trackerHitRelationCollection, m_inputTrackerHitRelationCollections[collection], evt);
    trackerHitRelationCollections.push_back(trackerHitRelationCollection);

    // Create the relations navigator
    LCRelationNavigator* relation = new LCRelationNavigator(trackerHitRelationCollection);
    relations.push_back(relation);
  }

  // Make the output track collection
  LCCollectionVec* trackCollection = new LCCollectionVec(LCIO::TRACK);

  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  trackCollection->setFlag(trkFlag.getFlag());

  // Make the output particle to track relation collection
  LCCollectionVec* trackRelationCollection = new LCCollectionVec(LCIO::LCRELATION);

  /*
Now for each MC particle we want the list of hits belonging to it. The most
efficient way is to loop over all hits once, and store the pointers in a
map, with the key a pointer to the MC particle. We can then loop over each
MC particle at the end and get all of the hits, before making a track.
*/

  // Make the container
  std::map<MCParticle*, std::vector<TrackerHit*>> particleHits;

  // Loop over all input collections
  for (unsigned int collection = 0; collection < trackerHitCollections.size(); collection++) {
    // Loop over tracker hits
    int nHits = trackerHitCollections[collection]->getNumberOfElements();
    for (int itHit = 0; itHit < nHits; itHit++) {
      // Get the hit
      TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(trackerHitCollections[collection]->getElementAt(itHit));

      // Get the related simulated hit(s)
      const LCObjectVec& simHitVector = relations[collection]->getRelatedToObjects(hit);

      // Take the first hit only (this should be changed? Yes - loop over all related simHits and add an entry for each
      // mcparticle so that this hit is in each fit)
      SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));

      // If the hit was produced by a secondary which was not saved to the MCParticle collection
      if (simHit->isProducedBySecondary())
        continue;

      // Get the particle belonging to that hit
      MCParticle* particle = simHit->getMCParticle();

      // Push back the element into the container
      particleHits[particle].push_back(hit);
    }
  }

  // Now loop over all particles and get the list of hits
  int nParticles = particleCollection->getNumberOfElements();
  for (int itP = 0; itP < nParticles; itP++) {
    // Get the particle
    MCParticle* mcParticle = static_cast<MCParticle*>(particleCollection->getElementAt(itP));

    // Get the vector of hits from the container
    if (particleHits.count(mcParticle) == 0)
      continue;
    std::vector<TrackerHit*> trackHits = particleHits[mcParticle];

    // Only make tracks with 3 or more hits
    if (trackHits.size() < 3)
      continue;

    // Sort the hits from smaller to larger radius
    std::sort(trackHits.begin(), trackHits.end(), sort_by_radius);

    // Remove the hits on the same layers (removing those with higher R)
    EVENT::TrackerHitVec trackFilteredByRHits;
    removeHitsSameLayer(trackHits, trackFilteredByRHits);
    if (trackFilteredByRHits.size() < 3)
      continue;

    /*
     Fit - this gets complicated.
     Need to pass a series of objects, including some initial states,
     covariance matrix etc. Set these up, then call the fit.
    */

    // Make the track object and relations object
    LCRelationImpl* relationTrack = new LCRelationImpl;
    TrackImpl* track = new TrackImpl;

    // IMarlinTrk used to fit track - IMarlinTrk interface to separete pattern recogition from fit implementation
    MarlinTrk::IMarlinTrack* marlinTrack = trackFactory->createTrack();
    MarlinTrk::IMarlinTrack* marlinTrackZSort = trackFactory->createTrack();

    // Save a vector of the hits to be used (why is this not attached to the track directly?? MarlinTrkUtils to be
    // updated?)
    EVENT::TrackerHitVec trackfitHits;
    for (unsigned int itTrackHit = 0; itTrackHit < trackFilteredByRHits.size(); itTrackHit++)
      trackfitHits.push_back(trackFilteredByRHits[itTrackHit]);

    // Make an initial covariance matrix with very broad default values
    EVENT::FloatVec covMatrix(15, 0);           // Size 15, filled with 0s
    covMatrix[0] = (m_initialTrackError_d0);    // sigma_d0^2
    covMatrix[2] = (m_initialTrackError_phi0);  // sigma_phi0^2
    covMatrix[5] = (m_initialTrackError_omega); // sigma_omega^2
    covMatrix[9] = (m_initialTrackError_z0);    // sigma_z0^2
    covMatrix[14] = (m_initialTrackError_tanL); // sigma_tanl^2

    bool direction = ((m_fitForward) ? MarlinTrk::IMarlinTrack::forward : MarlinTrk::IMarlinTrack::backward);

    int fitError = 0;

    if (m_useTruthInPrefit) {
      HelixTrack helix(mcParticle->getVertex(), mcParticle->getMomentum(), mcParticle->getCharge(), m_magneticField);

      double trueD0 = helix.getD0();
      double truePhi = helix.getPhi0();
      double trueOmega = helix.getOmega();
      double trueZ0 = helix.getZ0();
      double trueTanLambda = helix.getTanLambda();

      // float ref_point[3] = { 0., 0., 0. };
      helix.moveRefPoint(trackfitHits.at(0)->getPosition()[0], trackfitHits.at(0)->getPosition()[1],
                         trackfitHits.at(0)->getPosition()[2]);
      float ref_point[3] = {float(helix.getRefPointX()), float(helix.getRefPointY()), float(helix.getRefPointZ())};
      TrackStateImpl* trkState =
          new TrackStateImpl(TrackState::AtIP, trueD0, truePhi, trueOmega, trueZ0, trueTanLambda, covMatrix, ref_point);

      // int prefitError =  createFit(trackfitHits, marlinTrack, trkState, m_magneticField,  direction,
      // m_maxChi2perHit); streamlog_out(DEBUG2) << "---- createFit - error_fit = " << error_fit << std::endl;

      // if (prefitError == 0) {
      //   fitError = finaliseLCIOTrack(marlinTrack, track, trackfitHits,  direction );
      //   streamlog_out(DEBUG2) << "---- finalisedLCIOTrack - error = " << error << std::endl;
      // }

      fitError = MarlinTrk::createFinalisedLCIOTrack(marlinTrack, trackfitHits, track, direction, trkState,
                                                     m_magneticField, m_maxChi2perHit);

      // If first fit attempt fails, try a new fit with hits ordered by z

      if (fitError != 0) {
        // Sort the hits from smaller to larger z
        std::sort(trackfitHits.begin(), trackfitHits.end(), sort_by_z);

        // Removing the hits on the same layers (remove those with higher z)
        EVENT::TrackerHitVec trackFilteredByZHits;
        removeHitsSameLayer(trackfitHits, trackFilteredByZHits);
        if (trackFilteredByZHits.size() < 3)
          continue;

        // If fit with hits ordered by radius has failed, the track is probably a 'spiral' track.
        // Fitting 'backward' is very difficult for spiral track, so the default direction here is set as 'forward'
        fitError = MarlinTrk::createFinalisedLCIOTrack(marlinTrackZSort, trackFilteredByZHits, track,
                                                       MarlinTrk::IMarlinTrack::forward, trkState, m_magneticField,
                                                       m_maxChi2perHit);
      }

      delete trkState;

      ////////////////////////

    } // end: helical prefit initialised with info from truth and then track fitted and saved as a lcio track

    else {
      // DEFAULT procedure: Try to fit
      fitError = MarlinTrk::createFinalisedLCIOTrack(marlinTrack, trackfitHits, track, direction, covMatrix,
                                                     m_magneticField, m_maxChi2perHit);

      // If first fit attempt fails, try a new fit with hits ordered by z

      if (fitError != 0) {
        // we need to clean the track object
        delete track;
        track = new TrackImpl;

        // Sort the hits from smaller to larger z
        std::sort(trackfitHits.begin(), trackfitHits.end(), sort_by_z);

        // Removing the hits on the same layers (remove those with higher z)
        EVENT::TrackerHitVec trackFilteredByZHits;
        removeHitsSameLayer(trackfitHits, trackFilteredByZHits);
        if (trackFilteredByZHits.size() < 3)
          continue;

        // If fit with hits ordered by radius has failed, the track is probably a 'spiral' track.
        // Fitting 'backward' is very difficult for spiral track, so the default directiin here is set as 'forward'
        fitError = MarlinTrk::createFinalisedLCIOTrack(marlinTrackZSort, trackFilteredByZHits, track,
                                                       MarlinTrk::IMarlinTrack::forward, covMatrix, m_magneticField,
                                                       m_maxChi2perHit);
      }

    } // end: track fitted (prefit from fit from 3 hits) and saved as a lcio track

    /////////////////////////////////////////////

    streamlog_out(DEBUG2) << "TruthTrackFinder: fitError " << fitError << std::endl;

    // Check track quality - if fit fails chi2 will be 0
    if (fitError != 0) {
      delete track;
      delete relationTrack;
      delete marlinTrack;
      delete marlinTrackZSort;
      m_fitFails++;
      continue;
    }
    if (track->getChi2() <= 0.) {
      delete track;
      delete relationTrack;
      delete marlinTrack;
      delete marlinTrackZSort;
      m_fitFails++;
      continue;
    }
    if (track->getNdf() <= 0.) {
      delete track;
      delete relationTrack;
      delete marlinTrack;
      delete marlinTrackZSort;
      m_fitFails++;
      continue;
    }

    std::vector<std::pair<EVENT::TrackerHit*, double>> hits_in_fit;
    hits_in_fit.reserve(trackHits.size()); // Reserve at most the total number of hits
    marlinTrack->getHitsInFit(hits_in_fit);

    /// Fill hits associated to the track by pattern recognition and hits in fit
    m_encoder->reset(); // reset to 0
    MarlinTrk::addHitNumbersToTrack(track, trackHits, false, *m_encoder);
    MarlinTrk::addHitNumbersToTrack(track, hits_in_fit, true, *m_encoder);

    streamlog_out(DEBUG5) << "TruthTrackFinder: trackHits.size(): " << trackHits.size()
                          << " trackfitHits.size(): " << trackfitHits.size()
                          << " hits_in_fit.size(): " << hits_in_fit.size() << std::endl;

    // Push back to the output container
    trackCollection->addElement(track);

    // Make the particle to track link
    relationTrack->setFrom(track);
    relationTrack->setTo(mcParticle);
    relationTrack->setWeight(1.0);
    trackRelationCollection->addElement(relationTrack);

    delete marlinTrack;
    delete marlinTrackZSort;
  }

  // Save the output track collection
  evt->addCollection(trackCollection, m_outputTrackCollection);
  // Save the output particle to track relation collection
  evt->addCollection(trackRelationCollection, m_outputTrackRelationCollection);

  // Increment the event number
  m_eventNumber++;
  for (unsigned int collection = 0; collection < trackerHitCollections.size(); collection++)
    delete relations[collection];
}

void TruthTrackFinder::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void TruthTrackFinder::end() {
  streamlog_out(MESSAGE) << " end()  " << name() << " processed " << m_eventNumber << " events in " << m_runNumber
                         << " runs " << std::endl;

  delete m_encoder;
}

void TruthTrackFinder::getCollection(LCCollection*& collection, std::string collectionName, LCEvent* evt) {
  try {
    collection = evt->getCollection(collectionName);
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG5) << "- cannot get collection. Collection " << collectionName.c_str() << " is unavailable"
                          << std::endl;
    return;
  }
  return;
}

void TruthTrackFinder::removeHitsSameLayer(const std::vector<TrackerHit*>& trackHits,
                                           std::vector<TrackerHit*>& trackFilteredHits) {
  trackFilteredHits.push_back(*(trackHits.begin()));

  for (std::vector<TrackerHit*>::const_iterator it = trackHits.begin() + 1; it != trackHits.end(); ++it) {
    int subdet = getSubdetector(*it);
    int layer = getLayer(*it);
    if (subdet != getSubdetector(*(it - 1))) {
      trackFilteredHits.push_back(*it);
    } else if (layer != getLayer(*(it - 1))) {
      trackFilteredHits.push_back(*it);
    }
  }
}
