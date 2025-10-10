#include "FilterTimeHits.h"
#include <cmath>
#include <iostream>
#include <set>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <UTIL/LCRelationNavigator.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"

#include <marlin/AIDAProcessor.h>

using namespace lcio;
using namespace marlin;

FilterTimeHits aFilterTimeHits;

FilterTimeHits::FilterTimeHits() : Processor("FilterTimeHits") {
  // --- Processor description:

  _description = "FilterTimeHits selects tracker hits based on their time corrected for the time of flight";

  // --- Processor parameters:

  registerProcessorParameter("TrackerHitInputCollections", "Name of the tracker hit input collections",
                             m_inputTrackerHitsCollNames, {});

  registerProcessorParameter("TrackerSimHitInputCollections", "Name of the tracker simhit input collections",
                             m_inputTrackerSimHitsCollNames, {});

  registerProcessorParameter("TrackerHitInputRelations", "Name of the tracker hit relation collections",
                             m_inputTrackerHitRelNames, {});

  registerProcessorParameter("TrackerHitOutputCollections", "Name of the tracker hit output collections",
                             m_outputTrackerHitsCollNames, {});

  registerProcessorParameter("TrackerSimHitOutputCollections", "Name of the tracker simhit output collections",
                             m_outputTrackerSimHitsCollNames, {});

  registerProcessorParameter("TrackerHitOutputRelations", "Name of the tracker hit relation collections",
                             m_outputTrackerHitRelNames, {});

  registerProcessorParameter("TargetBeta", "Target beta=v/c for hit time of flight correction", m_beta, double(1.0));

  registerProcessorParameter("TimeLowerLimit", "Lower limit on the corrected hit time in ns", m_time_min,
                             double(-90.0));

  registerProcessorParameter("TimeUpperLimit", "Upper limit on the corrected hit time in ns", m_time_max, double(90.0));

  registerProcessorParameter("FillHistograms", "Flag to fill the diagnostic histograms", m_fillHistos, false);
}

void FilterTimeHits::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;

  // --- Print the processor parameters:

  printParameters();

  // --- Initialize the run and event counters:

  _nRun = 0;
  _nEvt = 0;

  // --- Initialize the AIDAProcessor and book the diagnostic histograms:

  AIDAProcessor::histogramFactory(this);

  m_corrected_time = new TH1F("m_corrected_time", "Corrected time of the hit [ps]", 1000, -250., 250.);
}

void FilterTimeHits::processRunHeader(LCRunHeader*) { _nRun++; }

void FilterTimeHits::processEvent(LCEvent* evt) {
  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() << "   in run:  " << evt->getRunNumber()
                       << std::endl;

  // --- Check whether the number of input and output collections match

  if (m_inputTrackerHitsCollNames.size() != m_inputTrackerSimHitsCollNames.size() ||
      m_inputTrackerHitsCollNames.size() != m_inputTrackerHitRelNames.size()) {
    std::stringstream err_msg;
    err_msg << "Mismatch between the reco and sim hits input collections" << std::endl;

    throw EVENT::Exception(err_msg.str());
  }

  if (m_outputTrackerHitsCollNames.size() != m_outputTrackerSimHitsCollNames.size() ||
      m_outputTrackerHitsCollNames.size() != m_outputTrackerHitRelNames.size()) {
    std::stringstream err_msg;
    err_msg << "Mismatch between the reco and sim hits output collections" << std::endl;

    throw EVENT::Exception(err_msg.str());
  }

  streamlog_out(DEBUG) << "   passed container size checks" << std::endl;

  // --- Get the input hit collections and create the corresponding output collections:

  const unsigned int nTrackerHitCol = m_inputTrackerHitsCollNames.size();
  std::vector<LCCollection*> inputHitColls(nTrackerHitCol);
  std::vector<LCCollection*> inputSimHitColls(nTrackerHitCol);
  std::vector<LCCollection*> inputHitRels(nTrackerHitCol);

  LCCollectionVec* outputTrackerHitCol = 0;
  LCCollectionVec* outputTrackerSimHitCol = 0;
  LCCollection* outputTrackerHitRel = 0;

  for (unsigned int icol = 0; icol < nTrackerHitCol; ++icol) {
    // get the reco hits
    try {
      inputHitColls[icol] = evt->getCollection(m_inputTrackerHitsCollNames[icol]);
    } catch (lcio::DataNotAvailableException& e) {
      streamlog_out(WARNING) << m_inputTrackerHitsCollNames[icol] << " collection not available" << std::endl;
      continue;
    }

    // get the sim hits
    try {
      inputSimHitColls[icol] = evt->getCollection(m_inputTrackerSimHitsCollNames[icol]);
    } catch (lcio::DataNotAvailableException& e) {
      streamlog_out(WARNING) << m_inputTrackerSimHitsCollNames[icol] << " collection not available" << std::endl;
      continue;
    }

    // get the reco-sim relations
    try {
      inputHitRels[icol] = evt->getCollection(m_inputTrackerHitRelNames[icol]);
    } catch (lcio::DataNotAvailableException& e) {
      streamlog_out(WARNING) << m_inputTrackerHitRelNames[icol] << " collection not available" << std::endl;
      continue;
    }
  }

  // --- Loop over the tracker hits and select hits inside the chosen time window:
  // std::vector<std::set<int>> hits_to_save(nTrackerHitCol);

  for (unsigned int icol = 0; icol < inputHitColls.size(); ++icol) {
    LCCollection* hit_col = inputHitColls[icol];
    if (!hit_col)
      continue;

    // reco hit output collections
    std::string encoderString = inputHitColls[icol]->getParameters().getStringVal("CellIDEncoding");
    outputTrackerHitCol = new LCCollectionVec(inputHitColls[icol]->getTypeName());
    outputTrackerHitCol->setSubset(true);
    outputTrackerHitCol->parameters().setValue("CellIDEncoding", encoderString);

    // sim hit output collections
    outputTrackerSimHitCol = new LCCollectionVec(inputSimHitColls[icol]->getTypeName());
    outputTrackerSimHitCol->parameters().setValue("CellIDEncoding", encoderString);
    outputTrackerSimHitCol->setSubset(true);

    // reco-sim relation output collections
    UTIL::LCRelationNavigator thitNav = UTIL::LCRelationNavigator(LCIO::TRACKERHITPLANE, LCIO::SIMTRACKERHIT);

    for (int ihit = 0; ihit < hit_col->getNumberOfElements(); ++ihit) {
      if (TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(hit_col->getElementAt(ihit))) {
        // Skipping the hit if its time is outside the acceptance time window
        double hitT = hit->getTime();

        dd4hep::rec::Vector3D pos = hit->getPosition();
        double hitR = pos.r();

        // Correcting for the propagation time
        double dt = hitR / (TMath::C() * m_beta / 1e6);
        hitT -= dt;
        streamlog_out(DEBUG3) << "corrected hit at R: " << hitR << " mm by propagation time: " << dt
                              << " ns to T: " << hitT << " ns" << std::endl;

        // Apply time window selection
        if (hitT < m_time_min || hitT > m_time_max) {
          streamlog_out(DEBUG4) << "hit at T: " << hitT << " ns is rejected by timing cuts" << std::endl;
          continue;
        }

        // hits_to_save[icol].insert(ihit);
        outputTrackerHitCol->addElement(hit);

        LCRelation* rel = static_cast<LCRelation*>(inputHitRels[icol]->getElementAt(ihit));
        SimTrackerHit* simhit = static_cast<SimTrackerHit*>(rel->getTo());
        outputTrackerSimHitCol->addElement(simhit);

        thitNav.addRelation(hit, simhit);

        if (m_fillHistos)
          m_corrected_time->Fill(hitT);
      }

    } // ihit loop

    streamlog_out(MESSAGE) << " " << outputTrackerHitCol->getNumberOfElements()
                           << " hits added to the collections: " << m_outputTrackerHitsCollNames[icol] << ", "
                           << m_outputTrackerSimHitsCollNames[icol] << ", " << m_outputTrackerHitRelNames[icol]
                           << std::endl;

    evt->addCollection(outputTrackerHitCol, m_outputTrackerHitsCollNames[icol]);
    evt->addCollection(outputTrackerSimHitCol, m_outputTrackerSimHitsCollNames[icol]);
    outputTrackerHitRel = thitNav.createLCCollection();
    evt->addCollection(outputTrackerHitRel, m_outputTrackerHitRelNames[icol]);

    streamlog_out(DEBUG) << " output collection " << m_outputTrackerHitsCollNames[icol] << " of type "
                         << outputTrackerHitCol->getTypeName() << " added to the event \n"
                         << " output collection " << m_outputTrackerSimHitsCollNames[icol] << " of type "
                         << outputTrackerSimHitCol->getTypeName() << " added to the event \n"
                         << " output collection " << m_outputTrackerHitRelNames[icol] << " of type "
                         << outputTrackerHitRel->getTypeName() << " added to the event  " << std::endl;

  } // icol loop

  streamlog_out(DEBUG) << "   Event processed " << std::endl;

  _nEvt++;
}

void FilterTimeHits::check(LCEvent*) {}

void FilterTimeHits::end() {
  std::cout << "FilterTimeHits::end()  " << name() << " processed " << _nEvt << " events in " << _nRun << " runs "
            << std::endl;
}
