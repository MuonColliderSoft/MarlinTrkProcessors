#include "CalcTrackerHitResiduals.h"
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

#include "MarlinTrk/Factory.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"

/// root headers
#include <TFile.h>
#include <TH1F.h>

using namespace lcio;
using namespace marlin;

using namespace MarlinTrk;
using namespace dd4hep::rec;

CalcTrackerHitResiduals aCalcTrackerHitResiduals;

CalcTrackerHitResiduals::CalcTrackerHitResiduals() : Processor("CalcTrackerHitResiduals") {
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

  _root_file = 0;

  _histo_map.clear();
  _histo_map_it = _histo_map.begin();
}

void CalcTrackerHitResiduals::init() {
  streamlog_out(DEBUG) << "   init called  " << std::endl;

  _encoder = new UTIL::BitField64(lcio::LCTrackerCellID::encoding_string());

  // usually a good idea to
  printParameters();

  // FIXME:SJA: if we want the surface store to be filled we need to create an instance of MarlinTrk implemented with
  // KalTest/KalDet
  MarlinTrk::IMarlinTrkSystem* trksystem = MarlinTrk::Factory::createMarlinTrkSystem("DDKalTest", 0, "");

  if (trksystem == 0) {
    throw EVENT::Exception(std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("DDKalTest"));
  }

  trksystem->init();

  /// Write histogram to file
  _root_file = new TFile("CalcTrackerHitResiduals.root", "RECREATE");
  this->createHistogramBuffers();

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

  streamlog_out(DEBUG9) << " get the surface manager from theDetector ... " << std::endl;

  const dd4hep::rec::SurfaceManager* surfMan = theDetector.extension<dd4hep::rec::SurfaceManager>();

  _surfMap = *surfMan->map("world");
}

void CalcTrackerHitResiduals::processRunHeader(LCRunHeader*) { ++_n_run; }

void CalcTrackerHitResiduals::processEvent(LCEvent* evt) {
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

  /**********************************************************************************************/
  /*                Loop over the tracker hit collections                                       */
  /**********************************************************************************************/

  std::vector<std::pair<SimTrackerHit*, TrackerHit*>> simHitTrkHit;

  for (unsigned iCol = 0; iCol < _colTrackerHits.size(); iCol++) {
    LCCollection* trackerHitCol = _colTrackerHits[iCol];
    int nHits = trackerHitCol->getNumberOfElements();

    streamlog_out(DEBUG1) << "Process " << _colNamesTrackerHits[iCol] << " collection with " << nHits << " hits \n";

    double rec_pos[3];
    double sim_pos[3];

    LCRelationNavigator* nav = _navTrackerHitRel[iCol];

    /**********************************************************************************************/
    /*                Loop over the tracker hits in this collection                               */
    /**********************************************************************************************/

    for (int j = 0; j < nHits; j++) {
      if (trackerHitCol->getElementAt(j) == 0) {
        streamlog_out(DEBUG1) << "Pointer to TrackerHit" << j << " is NULL " << std::endl;
      }

      TrackerHit* trkhit = dynamic_cast<TrackerHit*>(trackerHitCol->getElementAt(j));

      if (trkhit == 0) {
        std::stringstream errorMsg;
        errorMsg << "dynamic_cast to TrackerHit for hit " << j
                 << " failed. Pointer = " << trackerHitCol->getElementAt(j) << std::endl;

        throw lcio::Exception(errorMsg.str());
      }

      const LCObjectVec& to = nav->getRelatedToObjects(trkhit);

      const int celId = trkhit->getCellID0();

      UTIL::BitField64 encoder(lcio::LCTrackerCellID::encoding_string());

      encoder.setValue(celId);
      int side = encoder[lcio::LCTrackerCellID::side()];
      int layer = encoder[lcio::LCTrackerCellID::layer()];
      int module = encoder[lcio::LCTrackerCellID::module()];
      int sensor = encoder[lcio::LCTrackerCellID::sensor()];

      streamlog_out(DEBUG3) << "Hit = " << j << " has celId " << celId << std::endl;
      streamlog_out(DEBUG3) << "side = " << side << std::endl;
      streamlog_out(DEBUG3) << "layerNumber = " << layer << std::endl;
      streamlog_out(DEBUG3) << "moduleNumber = " << module << std::endl;
      streamlog_out(DEBUG3) << "sensorNumber = " << sensor << std::endl;

      rec_pos[0] = trkhit->getPosition()[0];
      rec_pos[1] = trkhit->getPosition()[1];
      rec_pos[2] = trkhit->getPosition()[2];

      dd4hep::rec::SurfaceMap::const_iterator si = _surfMap.find(celId);
      ISurface* ms = (si != _surfMap.end() ? si->second : 0);

      if (ms == NULL) {
        streamlog_out(DEBUG3) << " no surface found for hit  id: " << celId << std::endl;
        continue;
      }

      Vector3D globalPointRec(rec_pos[0], rec_pos[1], rec_pos[2]);
      Vector2D localPointRec = ms->globalToLocal(globalPointRec);

      if (UTIL::BitSet32(trkhit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) {
        /**********************************************************************************************/
        /*                Treat Composite Spacepoints                                                 */
        /**********************************************************************************************/

        // Check that the space point is only created from 2 strip hits
        if (to.size() == 2) {
          SimTrackerHit* simhitA = dynamic_cast<SimTrackerHit*>(to.at(0));
          SimTrackerHit* simhitB = dynamic_cast<SimTrackerHit*>(to.at(1));

          // Check if the simHits are from the same particle in order to avoid problems
          if (simhitA->getMCParticle() == simhitB->getMCParticle()) {
            sim_pos[0] = simhitA->getPosition()[0];
            sim_pos[1] = simhitA->getPosition()[1];
            sim_pos[2] = simhitA->getPosition()[2];

            Vector3D globalPointSim(sim_pos[0], sim_pos[1], sim_pos[2]);

            double dx = globalPointRec.x() - globalPointSim.x();
            double dy = globalPointRec.y() - globalPointSim.y();
            double dz = globalPointRec.z() - globalPointSim.z();

            double dr = globalPointRec.trans() - globalPointSim.trans();
            double drphi =
                globalPointRec.phi() * globalPointRec.trans() - globalPointSim.phi() * globalPointSim.trans();

            streamlog_out(DEBUG1) << "Spacepoint Residuals dx " << dx << " dy " << dy << " dz " << dz << "\n";

            fill_histo(_colNamesTrackerHits[iCol] + "_res_x", dx);
            fill_histo(_colNamesTrackerHits[iCol] + "_res_y", dy);
            fill_histo(_colNamesTrackerHits[iCol] + "_res_z", dz);

            fill_histo(_colNamesTrackerHits[iCol] + "_res_r", dr);
            fill_histo(_colNamesTrackerHits[iCol] + "_res_rphi", drphi);

          } else {
            streamlog_out(DEBUG1) << "spacepoint discarded, because simHits are not equal " << simhitA->getMCParticle()
                                  << " != " << simhitB->getMCParticle() << "\n";
          }

        } else {
          streamlog_out(DEBUG1) << "spacepoint discarded, because it is related to " << to.size()
                                << "SimTrackerHits. It should be 2!\n";
        }

      } else {
        /**********************************************************************************************/
        /*                Treat Normal Hits                                                           */
        /**********************************************************************************************/

        if (to.size() == 1) { // only take trackerHits, that have only one related SimHit

          SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(to.at(0));

          sim_pos[0] = simhit->getPosition()[0];
          sim_pos[1] = simhit->getPosition()[1];
          sim_pos[2] = simhit->getPosition()[2];

          Vector3D globalPointSim(sim_pos[0], sim_pos[1], sim_pos[2]);
          Vector2D localPointSim = ms->globalToLocal(globalPointSim);

          double dx = globalPointRec.x() - globalPointSim.x();
          double dy = globalPointRec.y() - globalPointSim.y();
          double dz = globalPointRec.z() - globalPointSim.z();

          double du = localPointRec.u() - localPointSim.u();
          double dv = localPointRec.v() - localPointSim.v();
          double dw = 0; // by construction

          streamlog_out(DEBUG1) << "TrackerHit Residuals dx " << dx << " dy " << dy << " dz " << dz << " du " << du
                                << " dv " << dv << " dw " << dw << "\n";

          fill_histo(_colNamesTrackerHits[iCol] + "_res_x", dx);
          fill_histo(_colNamesTrackerHits[iCol] + "_res_y", dy);
          fill_histo(_colNamesTrackerHits[iCol] + "_res_z", dz);
          fill_histo(_colNamesTrackerHits[iCol] + "_res_du", du);
          fill_histo(_colNamesTrackerHits[iCol] + "_res_dv", dv);
          fill_histo(_colNamesTrackerHits[iCol] + "_res_dw", dw);

        }

        else {
          streamlog_out(DEBUG1) << "TrackerHit discarded, because it is related to " << to.size()
                                << "SimTrackerHits. It should be 1!\n";
        }
      }

      // calculate residual
    }
  }
  streamlog_out(DEBUG4) << "Number of SimTracker hits = " << simHitTrkHit.size() << std::endl;

  // some information output
  for (unsigned int i = 0; i < simHitTrkHit.size(); ++i) {
    SimTrackerHit* simHit = simHitTrkHit[i].first;
    TrackerHit* trackerHit = simHitTrkHit[i].second;

    streamlog_out(DEBUG2) << "Tracker hit: [" << i << "] = " << trackerHit << "  mcp = " << simHit->getMCParticle()
                          << " time = " << simHit->getTime() << " cellid (trkHit) = " << trackerHit->getCellID0()
                          << " (de" << getDetectorID(trackerHit) << ",si" << getSideID(trackerHit) << ",la"
                          << getLayerID(trackerHit) << ",mo" << getModuleID(trackerHit) << ",se"
                          << getSensorID(trackerHit) << ")" << std::endl;
  }

  ++_n_evt;
}

void CalcTrackerHitResiduals::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void CalcTrackerHitResiduals::end() {
  streamlog_out(DEBUG4) << "CalcTrackerHitResiduals::end()  " << name() << " processed " << _n_evt << " events in "
                        << _n_run << " runs " << std::endl;

  // convert any remaining buffers to histograms

  std::map<std::string, std::list<float>*>::iterator it = _histo_buffer_map.begin();

  streamlog_out(DEBUG4) << "Number of remaining buffers = " << _histo_buffer_map.size() << std::endl;

  for (it = _histo_buffer_map.begin(); it != _histo_buffer_map.end(); ++it) {
    streamlog_out(DEBUG4) << "Write remaining buffer " << it->first << std::endl;

    this->write_buffer_to_histo(it);

    streamlog_out(DEBUG4) << "Write remaining buffer: done " << it->first << std::endl;
  }

  _root_file->Write();
  _root_file->Close();
  delete _root_file;
  _root_file = 0;

  delete _encoder;
}

void CalcTrackerHitResiduals::SetupInputCollections(LCEvent* evt) {
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

LCCollection* CalcTrackerHitResiduals::GetCollection(LCEvent* evt, std::string colName) {
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

LCRelationNavigator* CalcTrackerHitResiduals::GetRelations(LCEvent* evt, std::string RelName) {
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

// void CalcTrackerHitResiduals::bookHistograms(){
//
//   _root_file->cd();
//
//   for( unsigned iCol=0; iCol<_colNamesTrackerHits.size(); iCol++){
//
//     int nbinsx = 100;
//     float xlow = 0.0;
//     float xup  = 10.0;
//
//     TH1F* h = 0;
//     std::string name;
//
//     name = _colNamesTrackerHits[iCol] + "_res_x";
//     h = new TH1F(name.c_str(), name.c_str(), nbinsx, xlow, xup);
//     _histo_map[name] = h;
//
//     name = _colNamesTrackerHits[iCol] + "_res_y";
//     h = new TH1F(name.c_str(), name.c_str(), nbinsx, xlow, xup);
//     _histo_map[name] = h;
//
//     name = _colNamesTrackerHits[iCol] + "_res_z";
//     h = new TH1F(name.c_str(), name.c_str(), nbinsx, xlow, xup);
//     _histo_map[name] = h;
//
//
//
//
//   }
//
//
// }

void CalcTrackerHitResiduals::createHistogramBuffers() {
  for (unsigned iCol = 0; iCol < _colNamesTrackerHits.size(); iCol++) {
    std::list<float>* buffer = 0;

    std::string name;

    name = _colNamesTrackerHits[iCol] + "_res_x";
    buffer = new std::list<float>;
    _histo_buffer_map[name] = buffer;

    name = _colNamesTrackerHits[iCol] + "_res_y";
    buffer = new std::list<float>;
    _histo_buffer_map[name] = buffer;

    name = _colNamesTrackerHits[iCol] + "_res_z";
    buffer = new std::list<float>;
    _histo_buffer_map[name] = buffer;

    name = _colNamesTrackerHits[iCol] + "_res_r";
    buffer = new std::list<float>;
    _histo_buffer_map[name] = buffer;

    name = _colNamesTrackerHits[iCol] + "_res_rphi";
    buffer = new std::list<float>;
    _histo_buffer_map[name] = buffer;

    name = _colNamesTrackerHits[iCol] + "_res_du";
    buffer = new std::list<float>;
    _histo_buffer_map[name] = buffer;

    name = _colNamesTrackerHits[iCol] + "_res_dv";
    buffer = new std::list<float>;
    _histo_buffer_map[name] = buffer;

    name = _colNamesTrackerHits[iCol] + "_res_dw";
    buffer = new std::list<float>;
    _histo_buffer_map[name] = buffer;
  }
}

void CalcTrackerHitResiduals::write_buffer_to_histo(std::map<std::string, std::list<float>*>::iterator it_buffer) {
  // to avoid outliers we will discard the top and bottom 10 percent of values and calculate the rms from that

  std::string name = it_buffer->first;
  std::list<float>* buffer = it_buffer->second;

  if (buffer->empty()) {
    return;
  }

  buffer->sort();

  std::list<float>::iterator it;
  std::list<float>::iterator it_first_10pc = buffer->begin();
  std::list<float>::iterator it_last_10pc = buffer->end();

  int ten_percent = buffer->size() * 0.1;

  for (int i = 0; i < ten_percent; ++i) {
    it_first_10pc++;
    it_last_10pc--;
  }

  float sum2 = 0;
  float sum = 0;
  int i = 0;

  for (it = it_first_10pc; it != it_last_10pc; ++it) {
    sum += (*it);
    sum2 += (*it) * (*it);
    streamlog_out(MESSAGE) << "i = " << i << " value " << *it << " sum = " << sum << " sum2 = " << sum2 << std::endl;
    i++;
  }

  unsigned N = buffer->size() - 2 * ten_percent;
  float mean = sum / N;
  float rms = sqrt(sum2 / N - mean * mean);

  streamlog_out(MESSAGE) << "write_buffer_to_histo: Histo Name = " << name << " sum2 = " << sum2 << " sum = " << sum
                         << " rms = " << rms << " mean = " << mean << " N = " << N << std::endl;

  int nbinsx = 100;

  float xlow = mean - 5 * rms;
  float xup = mean + 5 * rms;

  _root_file->cd();

  TH1F* histo = new TH1F(name.c_str(), name.c_str(), nbinsx, xlow, xup);
  _histo_map[name] = histo;

  for (it = buffer->begin(); it != buffer->end(); ++it) {
    histo->Fill(*it);
  }
}

void CalcTrackerHitResiduals::fill_histo(const std::string& name, float value) {
  streamlog_out(DEBUG1) << "fill_histo: Histo Name = " << name << std::endl;

  // first check if the buffer is still being used

  std::map<std::string, std::list<float>*>::iterator it = _histo_buffer_map.begin();

  it = _histo_buffer_map.find(name);

  if (it != _histo_buffer_map.end()) {
    streamlog_out(DEBUG1) << "fill_histo: Use buffer" << std::endl;

    std::list<float>* buffer = it->second;

    buffer->push_back(value);

    // now check if the buffer is full

    if (buffer->size() == MAXBUFFERSIZE) {
      write_buffer_to_histo(it);

      delete buffer;
      _histo_buffer_map.erase(it);
    }

    // else just fill the historgram directly

  } else {
    streamlog_out(DEBUG1) << "fill_histo: Use Histo" << std::endl;

    _histo_map_it = _histo_map.find(name);

    if (_histo_map_it == _histo_map.end()) {
      for (_histo_map_it = _histo_map.begin(); _histo_map_it != _histo_map.end(); ++_histo_map_it) {
        streamlog_out(ERROR) << " Histo Name = " << _histo_map_it->first << " Pointer = " << _histo_map_it->second
                             << std::endl;
      }

      throw EVENT::Exception(std::string(" Histogram ") + name + std::string(" does not exist in the _histo_map"));
    }

    TH1F* histo = _histo_map_it->second;

    histo->Fill(value);
  }
}
