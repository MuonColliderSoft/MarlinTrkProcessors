#include "SetTrackerHitExtensions.h"
#include <iostream>

#include <cstdlib>
#include <vector>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>

#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>

#include <IMPL/LCRelationImpl.h>
#include <UTIL/ILDConf.h>
#include <UTIL/LCRelationNavigator.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "MarlinTrk/MarlinTrkDiagnostics.h"
#ifdef MARLINTRK_DIAGNOSTICS_ON
#include "MarlinTrk/DiagnosticsController.h"
#endif

using namespace lcio;
using namespace marlin;

using namespace MarlinTrk;

SetTrackerHitExtensions aSetTrackerHitExtensions;

SetTrackerHitExtensions::SetTrackerHitExtensions() : Processor("SetTrackerHitExtensions") {
  // modify processor description
  _description = "Creates Track Collection from MC Truth. Can handle composite spacepoints as long as they consist of "
                 "two TrackerHits";

  // register steering parameters: name, description, class-variable, default value

  StringVec trackerHitsRelInputColNamesDefault;
  trackerHitsRelInputColNamesDefault.push_back("VXDTrackerHitRelations");
  trackerHitsRelInputColNamesDefault.push_back("SITTrackerHitRelations");
  trackerHitsRelInputColNamesDefault.push_back("FTDPixelTrackerHitRelations");
  trackerHitsRelInputColNamesDefault.push_back("FTDSpacePointRelations");
  trackerHitsRelInputColNamesDefault.push_back("TPCTrackerHitRelations");
  trackerHitsRelInputColNamesDefault.push_back("SETTrackerHitRelations");

  registerInputCollections("LCRelation", "TrackerHitsRelInputCollections",
                           "Name of the lcrelation collections, that link the TrackerHits to their SimTrackerHits. "
                           "Have to be in same order as TrackerHitsInputCollections!!!",
                           _colNamesTrackerHitRelations, trackerHitsRelInputColNamesDefault);

  StringVec trackerHitsInputColNamesDefault;

  trackerHitsInputColNamesDefault.push_back("VXDTrackerHits");
  trackerHitsInputColNamesDefault.push_back("SITTrackerHits");
  trackerHitsInputColNamesDefault.push_back("FTDPixelTrackerHits");
  trackerHitsInputColNamesDefault.push_back("FTDSpacePointRelations");
  trackerHitsInputColNamesDefault.push_back("TPCTrackerHits");
  trackerHitsInputColNamesDefault.push_back("SETTrackerHits");

  registerInputCollections("TrackerHit", "TrackerHitsInputCollections", "Name of the tracker hit input collections",
                           _colNamesTrackerHits, trackerHitsInputColNamesDefault);

  _n_run = 0;
  _n_evt = 0;
}

void SetTrackerHitExtensions::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;
  _encoder = new UTIL::BitField64(lcio::LCTrackerCellID::encoding_string());

  // usually a good idea to
  printParameters();
}

void SetTrackerHitExtensions::processRunHeader(LCRunHeader*) { ++_n_run; }

void SetTrackerHitExtensions::processEvent(LCEvent* evt) {
  streamlog_out(DEBUG3) << "   processing event: " << _n_evt << std::endl;

  _current_evt_number = evt->getEventNumber();

  _colTrackerHits.clear();
  _navTrackerHitRel.clear();

  /**********************************************************************************************/
  /*                Prepare the collections                                                     */
  /**********************************************************************************************/

  // get the input collections and fill the vectors
  this->SetupInputCollections(evt);

  // create the encoder to decode cellID0
  UTIL::BitField64 cellID_encoder(LCTrackerCellID::encoding_string());

  for (unsigned iCol = 0; iCol < _colTrackerHits.size(); iCol++) {
    LCCollection* trackerHitCol = _colTrackerHits[iCol];
    int nHits = trackerHitCol->getNumberOfElements();

    LCRelationNavigator* nav = _navTrackerHitRel[iCol];

    for (int j = 0; j < nHits; j++) {
      if (trackerHitCol->getElementAt(j) == 0) {
        streamlog_out(DEBUG0) << "Pointer to TrackerHit" << j << " is NULL " << std::endl;
      }

      TrackerHit* trkhit = dynamic_cast<TrackerHit*>(trackerHitCol->getElementAt(j));

      if (trkhit == 0) {
        std::stringstream errorMsg;
        errorMsg << "dynamic_cast to TrackerHit for hit " << j
                 << " failed. Pointer = " << trackerHitCol->getElementAt(j) << std::endl;

        throw lcio::Exception(errorMsg.str());
      }

      const LCObjectVec& to = nav->getRelatedToObjects(trkhit);

      if (BitSet32(trkhit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) { // it is a composite spacepoint

        if (to.size() == 2) {
#ifdef MARLINTRK_DIAGNOSTICS_ON

          SimTrackerHit* simhitA = dynamic_cast<SimTrackerHit*>(to.at(0));
          SimTrackerHit* simhitB = dynamic_cast<SimTrackerHit*>(to.at(1));

          // set the pointer to the simhit via lcio extention MCTruth4HitExt

          const LCObjectVec rawObjects = trkhit->getRawHits();

          for (unsigned k = 0; k < rawObjects.size(); k++) {
            TrackerHit* rawHit = dynamic_cast<TrackerHit*>(rawObjects[k]);
            if (rawHit) {
              if (rawHit->getCellID0() == simhitA->getCellID0()) {
                streamlog_out(DEBUG4) << "link simhit = " << simhitA << " Cell ID = " << simhitA->getCellID0()
                                      << " with trkhit = " << rawHit << " Cell ID = " << rawHit->getCellID0()
                                      << std::endl;
                rawHit->ext<MarlinTrk::MCTruth4HitExt>() = new MarlinTrk::MCTruth4HitExtStruct;
                rawHit->ext<MarlinTrk::MCTruth4HitExt>()->simhit = simhitA;
              }
              if (rawHit->getCellID0() == simhitB->getCellID0()) {
                streamlog_out(DEBUG4) << "link simhit = " << simhitB << " Cell ID = " << simhitB->getCellID0()
                                      << " with trkhit = " << rawHit << " Cell ID = " << rawHit->getCellID0()
                                      << std::endl;
                rawHit->ext<MarlinTrk::MCTruth4HitExt>() = new MarlinTrk::MCTruth4HitExtStruct;
                rawHit->ext<MarlinTrk::MCTruth4HitExt>()->simhit = simhitB;
              }
            }
          }
#endif

        } else {
          streamlog_out(DEBUG1) << "spacepoint discarded, because it is related to " << to.size()
                                << "SimTrackerHits. It should be 2!\n";
        }

      } else { // no composite spacepoint

        if (to.size() == 1) { // only take trackerHits, that have only one related SimHit

#ifdef MARLINTRK_DIAGNOSTICS_ON

          SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(to.at(0));

          trkhit->ext<MarlinTrk::MCTruth4HitExt>() = new MarlinTrk::MCTruth4HitExtStruct;
          trkhit->ext<MarlinTrk::MCTruth4HitExt>()->simhit = simhit;

          streamlog_out(DEBUG4) << "link simhit = " << simhit << " Cell ID = " << simhit->getCellID0()
                                << " with trkhit = " << trkhit << " Cell ID = " << trkhit->getCellID0() << std::endl;

#endif
        } else {
          streamlog_out(DEBUG1) << "TrackerHit discarded, because it is related to " << to.size()
                                << "SimTrackerHits. It should be 1!\n";
        }
      }
    }
  }

  ++_n_evt;
}

void SetTrackerHitExtensions::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void SetTrackerHitExtensions::end() {
  streamlog_out(DEBUG4) << "SetTrackerHitExtensions::end()  " << name() << " processed " << _n_evt << " events in "
                        << _n_run << " runs " << std::endl;

  delete _encoder;
}

void SetTrackerHitExtensions::SetupInputCollections(LCEvent* evt) {
  // Check if there are as many tracker hit input collections as relation collections
  if (_colNamesTrackerHits.size() != _colNamesTrackerHitRelations.size()) {
    streamlog_out(ERROR)
        << "There must be as many input collections of tracker Hits as of relations. At the moment, there are "
        << _colNamesTrackerHits.size() << " tracker hit collections and " << _colNamesTrackerHitRelations.size()
        << " relation collections passed as steering paremeters!\n";

    exit(1);
  }

  for (unsigned i = 0; i < _colNamesTrackerHits.size(); i++) {
    // the tracker hits
    LCCollection* colTrkHits = GetCollection(evt, _colNamesTrackerHits[i]);
    if (colTrkHits == NULL)
      continue;

    // the relations of them
    LCRelationNavigator* nav = GetRelations(evt, _colNamesTrackerHitRelations[i]);
    if (nav == NULL)
      continue;

    _colTrackerHits.push_back(colTrkHits);
    _navTrackerHitRel.push_back(nav);
  }
}

LCCollection* SetTrackerHitExtensions::GetCollection(LCEvent* evt, std::string colName) {
  LCCollection* col = NULL;

  try {
    col = evt->getCollection(colName.c_str());
    streamlog_out(DEBUG4) << " --> " << colName.c_str()
                          << " collection found, number of elements = " << col->getNumberOfElements() << std::endl;
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG4) << " --> " << colName.c_str() << " collection absent" << std::endl;
  }

  return col;
}

LCRelationNavigator* SetTrackerHitExtensions::GetRelations(LCEvent* evt, std::string RelName) {
  LCRelationNavigator* nav = NULL;
  LCCollection* col = NULL;

  try {
    col = evt->getCollection(RelName.c_str());
    nav = new LCRelationNavigator(col);
    streamlog_out(DEBUG4) << " --> " << RelName
                          << " track relation collection found, number of elements = " << col->getNumberOfElements()
                          << std::endl;
  } catch (DataNotAvailableException& e) {
    streamlog_out(ERROR) << " --> " << RelName.c_str() << " track relation collection absent" << std::endl;
  }

  return nav;
}
