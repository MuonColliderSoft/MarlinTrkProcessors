#include "FilterClusters.h"

#include <math.h>

#include <DD4hep/Detector.h>

#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCFlagImpl.h>

#include <marlin/AIDAProcessor.h>

#include <AIDA/ITree.h>

#include "marlin/VerbosityLevels.h"

FilterClusters aFilterClusters;

FilterClusters::FilterClusters()
    : Processor("FilterClusters")
{
  // modify processor description
  _description = "FilterClusters processor filters a collection of tracker hits based on cluster size and outputs a filtered collection";

  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("ThetaRanges",
                             "Divide theta into bins for different cluster size cuts",
                             _ThetaRanges,
                             {});
  registerProcessorParameter("ThetaBins",
                             "Number of bins in theta",
                             _ThetaBins,
                             {});
  registerProcessorParameter("ClusterSize",
                             "Maximum cluster size for each theta range",
                             _ClusterSize,
                             {});
  registerProcessorParameter("Layers",
                             "Layers to be filtered",
                             _Layers,
                             {});
  registerInputCollection(LCIO::SIMTRACKERHIT,
                          "InSimTrackerHitCollection",
                          "Name of the input sim tracker hit collection",
                          _InSimTrackerHitCollection,
                          _InSimTrackerHitCollection);

  registerInputCollection(LCIO::TRACKERHITPLANE,
                          "InTrackerHitCollection",
                          "Name of the input tracker hit collection",
                          _InTrackerHitCollection,
                          _InTrackerHitCollection);

  registerInputCollection(LCIO::LCRELATION,
                          "InRelationCollection",
                          "Name of the input relation collection",
                          _InRelationCollection,
                          _InRelationCollection);

  registerOutputCollection(LCIO::SIMTRACKERHIT,
                           "OutSimTrackerHitCollection",
                           "Name of output sim tracker hit collection",
                           _OutSimTrackerHitCollection,
                           _OutSimTrackerHitCollection);

  registerOutputCollection(LCIO::TRACKERHITPLANE,
                           "OutTrackerHitCollection",
                           "Name of output tracker hit collection",
                           _OutTrackerHitCollection,
                           _OutTrackerHitCollection);

  registerOutputCollection(LCIO::LCRELATION,
                           "OutRelationCollection",
                           "Name of output relation collection",
                           _OutRelationCollection,
                           _OutRelationCollection);

  registerProcessorParameter("FillHistograms",
                             "Flag to fill the diagnostic histograms",
                             m_fillHistos,
                             false);
}

void FilterClusters::init()
{
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

void FilterClusters::processRunHeader(LCRunHeader * /*run*/)
{
}

void FilterClusters::processEvent(LCEvent *evt)
{
  // Get input collection
  LCCollection *InTrackerHitCollection = evt->getCollection(_InTrackerHitCollection);
  LCCollection *InRelationCollection = evt->getCollection(_InRelationCollection);
  LCCollection *InSimTrackerHitCollection = evt->getCollection(_InSimTrackerHitCollection);

  streamlog_out(DEBUG0) << "Starting processing event." << std::endl;

  if (InTrackerHitCollection->getTypeName() != lcio::LCIO::TRACKERHITPLANE)
  {
    throw EVENT::Exception("Invalid collection type: " + InTrackerHitCollection->getTypeName());
    streamlog_out(DEBUG0) << "Wrong collection type for TrackerHitCollection. \n";
  }
  if (InRelationCollection->getTypeName() != lcio::LCIO::LCRELATION)
  {
    throw EVENT::Exception("Invalid collection type: " + InRelationCollection->getTypeName());
    streamlog_out(DEBUG0) << "Wrong collection type for InRelationCollection. \n";
  }
  if (InSimTrackerHitCollection->getTypeName() != lcio::LCIO::SIMTRACKERHIT)
  {
    throw EVENT::Exception("Invalid collection type: " + InSimTrackerHitCollection->getTypeName());
    streamlog_out(DEBUG0) << "Wrong collection type for SimTrackerHitCollection. \n";
  }

  streamlog_out(DEBUG1) << "Number of Elements in Tracker Hits Collection: " << InTrackerHitCollection->getNumberOfElements() << std::endl;

  // Make the output collections: reco hits, sim hits, reco-sim relationship
  std::string encoderString = InTrackerHitCollection->getParameters().getStringVal("CellIDEncoding");
  LCCollectionVec *OutTrackerHitCollection = new LCCollectionVec(LCIO::TRACKERHITPLANE);
  OutTrackerHitCollection->parameters().setValue("CellIDEncoding", encoderString);
  LCFlagImpl lcFlag(InTrackerHitCollection->getFlag());
  OutTrackerHitCollection->setFlag(lcFlag.getFlag());
  // OutTrackerHitCollection->setSubset(true);

  LCCollectionVec *OutSimTrackerHitCollection = new LCCollectionVec(LCIO::SIMTRACKERHIT);
  OutSimTrackerHitCollection->parameters().setValue("CellIDEncoding", encoderString);
  LCFlagImpl lcFlag_sim(InSimTrackerHitCollection->getFlag());
  OutSimTrackerHitCollection->setFlag(lcFlag_sim.getFlag());
  // OutSimTrackerHitCollection->setSubset(true);

  LCCollectionVec *OutRelationCollection = new LCCollectionVec(LCIO::LCRELATION);
  LCFlagImpl lcFlag_rel(InRelationCollection->getFlag());
  OutRelationCollection->setFlag(lcFlag_rel.getFlag());
  // OutRelationCollection->setSubset(true);

  // Filter
  for (int i = 0; i < InTrackerHitCollection->getNumberOfElements(); ++i) // loop through all hits
  {
    streamlog_out(DEBUG2) << "Loop over hits:" << i << std::endl;
    TrackerHitPlane *trkhit = static_cast<TrackerHitPlane *>(InTrackerHitCollection->getElementAt(i)); // define trkhit var, pointer to i'th element of tracker hits

    if (!trkhit)
    {
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
    const lcio::LCObjectVec &rawHits = trkhit->getRawHits();
    float ymax = -1000000;
    float xmax = -1000000;
    float ymin = 1000000;
    float xmin = 1000000;
    streamlog_out(DEBUG1) << "Looping over hits constituents." << std::endl;
    for (size_t j = 0; j < rawHits.size(); ++j)
    {
      lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit *>(rawHits[j]);
      const double *localPos = hitConstituent->getPosition();
      float x_local = localPos[0];
      float y_local = localPos[1];

      if (y_local < ymin)
      {
        ymin = y_local;
      }
      if (y_local > ymax)
      {
        ymax = y_local;
      }

      if (x_local < xmin)
      {
        xmin = x_local;
      }
      if (x_local > xmax)
      {
        xmax = x_local;
      }
    }
    float cluster_size_y = (ymax - ymin) + 1;
    float cluster_size_x = (xmax - xmin) + 1;
    float cluster_size = rawHits.size();

    streamlog_out(DEBUG2) << "Cluster size:" << cluster_size << std::endl;

    // Get hit subdetector/layer
    std::string _encoderString = lcio::LCTrackerCellID::encoding_string();
    UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
    uint32_t systemID = decoder(trkhit)["system"];
    uint32_t layerID = decoder(trkhit)["layer"];
    bool filter_layer = false;

    streamlog_out(DEBUG3) << "Decoded system/layer:" << systemID << "/" << layerID << std::endl;

    int rows = _Layers.size(), cols = std::stoi(_ThetaBins);
    if ((rows * cols != _ThetaRanges.size()) || (rows * (cols - 1) != _ClusterSize.size()))
    {
      std::cout << "Either theta cuts or cluster cuts not provided for each layer. Please change the config, exiting now..." << std::endl;
      return;
    }

    std::vector<std::vector<float>> _thetaCuts_byLayer;
    std::vector<std::vector<float>> _clusterSizeCuts_byLayer;

    for (int k = 0; k < rows; ++k)
    {
      std::vector<float> row;
      for (int j = 0; j < cols; ++j)
      {
        row.push_back(std::stof(_ThetaRanges[j]));
      }
      _thetaCuts_byLayer.push_back(row);
    }

    for (int k = 0; k < rows; ++k)
    {
      std::vector<float> row;
      for (int j = 0; j < cols; ++j)
      {
        row.push_back(std::stof(_ClusterSize[j]));
      }
      _clusterSizeCuts_byLayer.push_back(row);
    }

    for (size_t j = 0; j < _Layers.size(); ++j)
    {
      if (layerID == std::stof(_Layers[j]))
      {
        filter_layer = true;
        break;
      }
    }
    streamlog_out(DEBUG1) << "Filter layer: " << filter_layer << std::endl;

    if (m_fillHistos)
    {
      m_clusterTheta_beforeCut->Fill(incidentTheta);
      m_clusterLayer_beforeCut->Fill(layerID);
      m_clusterSize_beforeCut->Fill(cluster_size);
    }

    bool store_hit = true;
    if (filter_layer)
    {
      store_hit = false;
      for (size_t j = 0; j < _thetaCuts_byLayer[layerID].size() - 1; ++j)
      {
        streamlog_out(DEBUG0) << "theta: " << incidentTheta << std::endl;
        float min_theta = _thetaCuts_byLayer[layerID][j];
        float max_theta = _thetaCuts_byLayer[layerID][j + 1];
        streamlog_out(DEBUG0) << "theta range: " << min_theta << ", " << max_theta << std::endl;

        if (incidentTheta >= min_theta and incidentTheta <= max_theta and filter_layer)
        {
          store_hit = true;
          streamlog_out(DEBUG0) << "theta in range" << std::endl;
          streamlog_out(DEBUG0) << "cluster size cut off: " << _clusterSizeCuts_byLayer[layerID][j] << std::endl;
          streamlog_out(DEBUG0) << "cluster size: " << cluster_size << std::endl;
          if (cluster_size < _clusterSizeCuts_byLayer[layerID][j])
          {
            streamlog_out(DEBUG0) << "Adding reco/sim clusters and relation to output collections" << std::endl;
          }
        }
      }
    }

    if (store_hit)
    {
      if (m_fillHistos)
      {
        m_clusterTheta_afterCut->Fill(incidentTheta);
        m_clusterLayer_afterCut->Fill(layerID);
        m_clusterSize_afterCut->Fill(cluster_size);
      }

      EVENT::LCRelation *rel = static_cast<EVENT::LCRelation *>(InRelationCollection->getElementAt(i));
      SimTrackerHit *simhit = dynamic_cast<SimTrackerHit *>(rel->getTo());
      SimTrackerHitImpl *simhit_new = new SimTrackerHitImpl();
      simhit_new->setCellID0(simhit->getCellID0());
      simhit_new->setCellID1(simhit->getCellID1());
      simhit_new->setPosition(simhit->getPosition());
      simhit_new->setEDep(simhit->getEDep());
      simhit_new->setTime(simhit->getTime());
      simhit_new->setMCParticle(simhit->getMCParticle());
      simhit_new->setMomentum(simhit->getMomentum());
      simhit_new->setPathLength(simhit->getPathLength());
      simhit_new->setQuality(simhit->getQuality());
      simhit_new->setOverlay(simhit->isOverlay());
      simhit_new->setProducedBySecondary(simhit->isProducedBySecondary());
      OutSimTrackerHitCollection->addElement(simhit_new);

      TrackerHitPlaneImpl *hit_new = new TrackerHitPlaneImpl();
      hit_new->setCellID0(trkhit->getCellID0());
      hit_new->setCellID1(trkhit->getCellID1());
      hit_new->setType(trkhit->getType());
      hit_new->setPosition(trkhit->getPosition());
      hit_new->setU(trkhit->getU());
      hit_new->setV(trkhit->getV());
      hit_new->setdU(trkhit->getdU());
      hit_new->setdV(trkhit->getdV());
      hit_new->setEDep(trkhit->getEDep());
      hit_new->setEDepError(trkhit->getEDepError());
      hit_new->setTime(trkhit->getTime());
      hit_new->setQuality(trkhit->getQuality());
      const lcio::LCObjectVec &rawHits_new = trkhit->getRawHits();
      for (size_t k = 0; k < rawHits_new.size(); ++k)
      {
        lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit *>(rawHits_new[k]);
        hit_new->rawHits().push_back(hitConstituent);
      }
      OutTrackerHitCollection->addElement(hit_new);

      LCRelationImpl *rel_new = new LCRelationImpl();
      rel_new->setFrom(hit_new);
      rel_new->setTo(simhit_new);
      rel_new->setWeight(rel->getWeight());
      OutRelationCollection->addElement(rel_new);
    }
    else
    {
      streamlog_out(DEBUG0) << "cluster rejected" << std::endl;
    }
  }

  // Save output track collection
  evt->addCollection(OutTrackerHitCollection, _OutTrackerHitCollection);
  evt->addCollection(OutRelationCollection, _OutRelationCollection);
  evt->addCollection(OutSimTrackerHitCollection, _OutSimTrackerHitCollection);

  streamlog_out(MESSAGE) << " " << OutTrackerHitCollection->size() << " reco clusters added to the collection: " << _OutTrackerHitCollection << ", " << OutSimTrackerHitCollection->size() << " sim hits added to the collection: " << _OutSimTrackerHitCollection << " and " << OutRelationCollection->size() << " relations added to the collection: " << _OutRelationCollection << std::endl;
}

void FilterClusters::end()
{
}
