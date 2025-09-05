#include "FilterClusters.h"

#include <math.h>

#include <DD4hep/Detector.h>

#include <EVENT/SimTrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTrackerConf.h>

#include <UTIL/LCRelationNavigator.h>

#include <marlin/AIDAProcessor.h>

#include <AIDA/ITree.h>

#include "marlin/VerbosityLevels.h"

FilterClusters aFilterClusters;

FilterClusters::FilterClusters() : Processor("FilterClusters") {
  // modify processor description
  _description = "FilterClusters processor filters a collection of tracker hits based on cluster size and outputs a "
                 "filtered collection";

  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("ThetaRanges", "Divide theta into bins for different cluster size cuts", _ThetaRanges, {});
  registerProcessorParameter("ThetaBins", "Number of bins in theta", _ThetaBins, {});
  registerProcessorParameter("ClusterSize", "Maximum cluster size for each theta range", _ClusterSize, {});
  registerProcessorParameter("Layers", "Layers to be filtered", _Layers, {});
  registerInputCollection(LCIO::SIMTRACKERHIT, "InSimTrackerHitCollection",
                          "Name of the input sim tracker hit collection", _InSimTrackerHitCollection,
                          _InSimTrackerHitCollection);

  registerInputCollection(LCIO::TRACKERHITPLANE, "InTrackerHitCollection", "Name of the input tracker hit collection",
                          _InTrackerHitCollection, _InTrackerHitCollection);

  registerInputCollection(LCIO::LCRELATION, "InRelationCollection", "Name of the input relation collection",
                          _InRelationCollection, _InRelationCollection);

  registerOutputCollection(LCIO::SIMTRACKERHIT, "OutSimTrackerHitCollection",
                           "Name of output sim tracker hit collection", _OutSimTrackerHitCollection,
                           _OutSimTrackerHitCollection);

  registerOutputCollection(LCIO::TRACKERHITPLANE, "OutTrackerHitCollection", "Name of output tracker hit collection",
                           _OutTrackerHitCollection, _OutTrackerHitCollection);

  registerOutputCollection(LCIO::LCRELATION, "OutRelationCollection", "Name of output relation collection",
                           _OutRelationCollection, _OutRelationCollection);

  registerProcessorParameter("FillHistograms", "Flag to fill the diagnostic histograms", m_fillHistos, false);
}

void FilterClusters::init() {
  // Print the initial parameters
  printParameters();

  marlin::AIDAProcessor::histogramFactory(this);

  m_clusterTheta_beforeCut = new TH1F("m_ClusterTheta_before", "Cluster Theta [radian]", 32, 0., 3.2);
  m_clusterTheta_afterCut = new TH1F("m_ClusterTheta_after", "Cluster Theta [radian]", 32, 0., 3.2);
  m_clusterLayer_beforeCut = new TH1F("m_ClusterLayer_before", "Cluster Layer Index", 10, 0, 10);
  m_clusterLayer_afterCut = new TH1F("m_ClusterLayer_after", "Cluster Layer Index", 10, 0, 10);
  m_clusterSize_beforeCut = new TH1F("m_ClusterSize_before", "Cluster hit multiplicity", 20, 0, 20);
  m_clusterSize_afterCut = new TH1F("m_ClusterSize_after", "Cluster hit multiplicity", 20, 0, 20);
}

void FilterClusters::processRunHeader(LCRunHeader* /*run*/) {}

void FilterClusters::processEvent(LCEvent* evt) {
  // Get input collection
  LCCollection* InTrackerHitCollection = evt->getCollection(_InTrackerHitCollection);
  LCCollection* InRelationCollection = evt->getCollection(_InRelationCollection);
  LCCollection* InSimTrackerHitCollection = evt->getCollection(_InSimTrackerHitCollection);

  streamlog_out(DEBUG0) << "Starting processing event." << std::endl;

  if (InTrackerHitCollection->getTypeName() != lcio::LCIO::TRACKERHITPLANE) {
    throw EVENT::Exception("Invalid collection type: " + InTrackerHitCollection->getTypeName());
    streamlog_out(DEBUG0) << "Wrong collection type for TrackerHitCollection. \n";
  }
  if (InRelationCollection->getTypeName() != lcio::LCIO::LCRELATION) {
    throw EVENT::Exception("Invalid collection type: " + InRelationCollection->getTypeName());
    streamlog_out(DEBUG0) << "Wrong collection type for InRelationCollection. \n";
  }
  if (InSimTrackerHitCollection->getTypeName() != lcio::LCIO::SIMTRACKERHIT) {
    throw EVENT::Exception("Invalid collection type: " + InSimTrackerHitCollection->getTypeName());
    streamlog_out(DEBUG0) << "Wrong collection type for SimTrackerHitCollection. \n";
  }

  streamlog_out(DEBUG1) << "Number of Elements in Tracker Hits Collection: "
                        << InTrackerHitCollection->getNumberOfElements() << std::endl;

  // Make the output collections: reco hits, sim hits, reco-sim relationship
  std::string encoderString = InTrackerHitCollection->getParameters().getStringVal("CellIDEncoding");
  LCCollectionVec* OutTrackerHitCollection = new LCCollectionVec(InTrackerHitCollection->getTypeName());
  OutTrackerHitCollection->setSubset(true);
  OutTrackerHitCollection->parameters().setValue("CellIDEncoding", encoderString);

  // sim hit output collections
  LCCollectionVec* OutSimTrackerHitCollection = new LCCollectionVec(InSimTrackerHitCollection->getTypeName());
  OutSimTrackerHitCollection->parameters().setValue("CellIDEncoding", encoderString);
  OutSimTrackerHitCollection->setSubset(true);

  // reco-sim relation output collections
  UTIL::LCRelationNavigator thitNav = UTIL::LCRelationNavigator(LCIO::TRACKERHITPLANE, LCIO::SIMTRACKERHIT);

  // Filter
  for (int i = 0; i < InTrackerHitCollection->getNumberOfElements(); ++i) // loop through all hits
  {
    streamlog_out(DEBUG2) << "Loop over hits:" << i << std::endl;
    TrackerHitPlane* trkhit = static_cast<TrackerHitPlane*>(
        InTrackerHitCollection->getElementAt(i)); // define trkhit var, pointer to i'th element of tracker hits

    if (!trkhit) {
      streamlog_out(WARNING) << "Cannot retrieve valid point to cluster. Skipping it" << std::endl;
    }

    // Calculating theta
    float x = trkhit->getPosition()[0];
    float y = trkhit->getPosition()[1];
    float z = trkhit->getPosition()[2];
    float r = sqrt(pow(x, 2) + pow(y, 2));
    float incidentTheta = std::atan(r / z);

    if (incidentTheta < 0)
      incidentTheta += M_PI;

    // Calculating cluster size
    const lcio::LCObjectVec& rawHits = trkhit->getRawHits();
    float ymax = -1000000;
    float xmax = -1000000;
    float ymin = 1000000;
    float xmin = 1000000;
    streamlog_out(DEBUG1) << "Looping over hits constituents." << std::endl;
    for (size_t j = 0; j < rawHits.size(); ++j) {
      lcio::SimTrackerHit* hitConstituent = dynamic_cast<lcio::SimTrackerHit*>(rawHits[j]);
      const double* localPos = hitConstituent->getPosition();
      float x_local = localPos[0];
      float y_local = localPos[1];

      if (y_local < ymin) {
        ymin = y_local;
      }
      if (y_local > ymax) {
        ymax = y_local;
      }

      if (x_local < xmin) {
        xmin = x_local;
      }
      if (x_local > xmax) {
        xmax = x_local;
      }
    }
    // float cluster_size_y = (ymax - ymin) + 1;
    // float cluster_size_x = (xmax - xmin) + 1;
    float cluster_size = rawHits.size();

    streamlog_out(DEBUG2) << "Cluster size:" << cluster_size << std::endl;

    // Get hit subdetector/layer
    std::string _encoderString = lcio::LCTrackerCellID::encoding_string();
    UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
    uint32_t systemID = decoder(trkhit)["system"];
    uint32_t layerID = decoder(trkhit)["layer"];
    bool filter_layer = false;

    streamlog_out(DEBUG3) << "Decoded system/layer:" << systemID << "/" << layerID << std::endl;

    size_t rows = _Layers.size(), cols = std::stoi(_ThetaBins);
    if ((rows * (cols + 1) != _ThetaRanges.size()) || (rows * cols != _ClusterSize.size())) {
      std::cout
          << "Either theta cuts or cluster cuts not provided for each layer. Please change the config, exiting now..."
          << std::endl;
      return;
    }

    std::vector<std::vector<float>> _thetaBins_byLayer;
    std::vector<std::vector<float>> _clusterSizeCuts_byLayer;

    for (size_t k = 0; k < rows; ++k) {
      std::vector<float> row;
      for (size_t j = 0; j <= cols; ++j) {
        row.push_back(std::stof(_ThetaRanges[j + k * (cols + 1)]));
      }
      _thetaBins_byLayer.push_back(row);
    }

    for (size_t k = 0; k < rows; ++k) {
      std::vector<float> row;
      for (size_t j = 0; j < cols; ++j) {
        row.push_back(std::stof(_ClusterSize[j + k * cols]));
      }
      _clusterSizeCuts_byLayer.push_back(row);
    }

    int layerInd = -1;
    for (size_t j = 0; j < _Layers.size(); ++j) {
      if (layerID == std::stof(_Layers[j])) {
        filter_layer = true;
        layerInd = j;
        break;
      }
    }
    streamlog_out(DEBUG1) << "Filter layer: " << filter_layer << std::endl;

    if (m_fillHistos) {
      m_clusterTheta_beforeCut->Fill(incidentTheta);
      m_clusterLayer_beforeCut->Fill(layerID);
      m_clusterSize_beforeCut->Fill(cluster_size);
    }

    bool store_hit = true;
    if (filter_layer) {
      store_hit = false;
      for (size_t j = 0; j < _thetaBins_byLayer[layerInd].size() - 1; ++j) {
        streamlog_out(DEBUG0) << "theta: " << incidentTheta << std::endl;
        float min_theta = _thetaBins_byLayer[layerInd][j];
        float max_theta = _thetaBins_byLayer[layerInd][j + 1];
        streamlog_out(DEBUG0) << "theta range: " << min_theta << ", " << max_theta << std::endl;

        if (incidentTheta >= min_theta && incidentTheta <= max_theta && filter_layer) {
          streamlog_out(DEBUG0) << "theta in range" << std::endl;
          streamlog_out(DEBUG0) << "cluster size cut off: " << _clusterSizeCuts_byLayer[layerInd][j] << std::endl;
          streamlog_out(DEBUG0) << "cluster size: " << cluster_size << std::endl;
          if (cluster_size < _clusterSizeCuts_byLayer[layerInd][j]) {
            store_hit = true;
            streamlog_out(DEBUG0) << "Adding reco/sim clusters and relation to output collections" << std::endl;
            break;
          }
        }
      }
    }

    if (store_hit) {
      if (m_fillHistos) {
        m_clusterTheta_afterCut->Fill(incidentTheta);
        m_clusterLayer_afterCut->Fill(layerID);
        m_clusterSize_afterCut->Fill(cluster_size);
      }

      EVENT::LCRelation* rel = static_cast<EVENT::LCRelation*>(InRelationCollection->getElementAt(i));
      SimTrackerHit* simhit = dynamic_cast<SimTrackerHit*>(rel->getTo());

      OutSimTrackerHitCollection->addElement(simhit);
      OutTrackerHitCollection->addElement(trkhit);
      thitNav.addRelation(trkhit, simhit);
    } else {
      streamlog_out(DEBUG0) << "cluster rejected" << std::endl;
    }
  }

  LCCollection* OutRelationCollection = thitNav.createLCCollection();

  // Save output track collection
  evt->addCollection(OutTrackerHitCollection, _OutTrackerHitCollection);
  evt->addCollection(OutRelationCollection, _OutRelationCollection);
  evt->addCollection(OutSimTrackerHitCollection, _OutSimTrackerHitCollection);

  streamlog_out(MESSAGE) << " " << OutTrackerHitCollection->size()
                         << " reco clusters added to the collection: " << _OutTrackerHitCollection << ", "
                         << OutSimTrackerHitCollection->size()
                         << " sim hits added to the collection: " << _OutSimTrackerHitCollection << " and "
                         << OutRelationCollection->getNumberOfElements()
                         << " relations added to the collection: " << _OutRelationCollection << std::endl;
}

void FilterClusters::end() {}