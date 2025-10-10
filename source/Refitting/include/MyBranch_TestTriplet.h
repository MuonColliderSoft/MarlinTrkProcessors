#ifndef MyBranch_TestTriplet_H
#define MyBranch_TestTriplet_H 1

#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <iostream>

class MCPStatus {
public:
  //////BasicInitializeで初期化されるもの//////
  int pdg;
  double energy;
  double mass;
  float charge;
  float time;
  double vx;
  double vy;
  double vz;
  double px;
  double py;
  double pz;
  double Pt;
  double Pabs;
  int isCreatedInSimulation;
  int isBackscatter;
  int vertexIsNotEndpointOfParent;
  int isDecayedInTracker;
  int isDecayedInCalorimeter;
  int hasLeftDetector;
  int isStopped;
  ///////BasicInitializeで初期化されないもの/////////////////
  float d0;
  float phi0;
  float omega;
  float z0;
  float tanL;
  int nhits;
  int nvxdhits;
  int nsithits;
  int nftdhits;
  int ntpchits;
  int nsethits;
  MCPStatus() {
    d0 = 0;
    phi0 = 0;
    omega = 0;
    z0 = 0;
    tanL = 0;
    nhits = 0;
    nvxdhits = 0;
    nsithits = 0;
    nftdhits = 0;
    ntpchits = 0;
    nsethits = 0;
  }
  void BasicInitialize(MCParticle* mcp);
};

void MCPStatus::BasicInitialize(MCParticle* mcp) {
  pdg = mcp->getPDG();
  energy = mcp->getEnergy();
  charge = mcp->getCharge();
  mass = mcp->getMass();
  time = mcp->getTime();
  px = mcp->getMomentum()[0];
  py = mcp->getMomentum()[1];
  pz = mcp->getMomentum()[2];
  Pt = sqrt(pow(px, 2) + pow(py, 2));
  Pabs = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2));
  vx = mcp->getVertex()[0];
  vy = mcp->getVertex()[1];
  vz = mcp->getVertex()[2];
  isCreatedInSimulation = mcp->isCreatedInSimulation();
  isBackscatter = mcp->isBackscatter();
  vertexIsNotEndpointOfParent = mcp->vertexIsNotEndpointOfParent();
  isDecayedInTracker = mcp->isDecayedInTracker();
  isDecayedInCalorimeter = mcp->isDecayedInCalorimeter();
  hasLeftDetector = mcp->hasLeftDetector();
  isStopped = mcp->isStopped();
}

class MyBranch_TestTriplet {
public:
  int id;
  int quality_code;
  int detID[3];
  int layer[3];
  int Comb;
  float chi2;
  int ndf;
  float prob;
  float purity;
  float d0;
  float z0;
  float phi0;
  float omega;
  float tanL;
  float cov[15]; // lcioと同じ詰め方
  MCPStatus mcp;
};

#endif
