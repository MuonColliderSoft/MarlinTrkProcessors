#ifndef GETPURITY_h
#define GETPURITY_h 1

#include "UTIL/LCTrackerConf.h"
#include <EVENT/LCCollection.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <UTIL/LCRelationNavigator.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include "UTIL/LCTrackerConf.h"
#include <UTIL/BitField64.h>

using namespace lcio;
using namespace marlin;

class purityMCP {
public:
  float purity;
  MCParticle* mcp;
};

class GetPurityUtil {
public:
  GetPurityUtil() { ; }
  typedef std::vector<LCRelationNavigator*> Navec;
  typedef std::vector<SimTrackerHit*> SimTrackerHitVec;
  typedef std::map<MCParticle*, int> mymap;
  typedef mymap::iterator myIter;
  typedef std::multimap<int, MCParticle*> revMCPMap;
  typedef revMCPMap::iterator revIter;
  // navecはsimtrackerhit,trackerhitの奴を使いそうなやつだけpush_backしてからGetPurityの引数に入れる。

  purityMCP GetPurity(TrackerHitVec tvec, Navec navec) {
    // pMCPにはdominantMCPのアドレスを指した状態で返ってくる。
    int ntrkhits = tvec.size();
    SimTrackerHitVec simvec;
    int bkgcount = 0;
    for (int i = 0; i < ntrkhits; i++) {
      LCObject* trk = dynamic_cast<LCObject*>(tvec[i]);
      const LCObjectVec* lcoVec = NULL;
      int navsize = navec.size();
      for (int j = 0; i < navsize; j++) {
        if (navec[j] != NULL) {
          lcoVec = &navec[j]->getRelatedToObjects(
              trk); // 見つからない場合はNULLを返すわけじゃないので、sizeが0かどうかで判定する。
          if (lcoVec->size() != 0) {
            break;
          }
        }
      }
      int size = lcoVec->size();
      simvec.reserve(size);
      if (size == 0) {
        bkgcount++;
      } else {
        for (int jlc = 0; jlc < size; jlc++)
          simvec.push_back(dynamic_cast<SimTrackerHit*>(lcoVec->at(jlc)));
      }
    }
    mymap map;
    myIter it;
    for (int i = 0; i < int(simvec.size()); i++) {
      it = map.find(simvec[i]->getMCParticle());
      //( it != map.end() )  ? ( it->second++ ) : ( map.insert(std::make_pair(simvec[i]->getMCParticle(),1)) );
      // 上の表記は何故かできない模様。
      if (it != map.end()) {
        it->second++;
      } else
        map.insert(std::make_pair(simvec[i]->getMCParticle(), 1));
    }

    float purity = 0;
    revMCPMap revmap;
    for (myIter itm = map.begin(); itm != map.end(); itm++) {
      revmap.insert(std::make_pair(itm->second, itm->first));
    }
    MCParticle* dominantmcp = NULL;
#if 0
      std::cout << "revmap.size() : " << revmap.size() << std::endl;
#endif
    if (revmap.size() != 0) {
      revIter maxit = std::max_element(revmap.begin(), revmap.end());
      dominantmcp = maxit->second;
#if 0
         std::cout << "dominantmcp->getPDG() : " << dominantmcp->getPDG() << std::endl;
#endif
      float ndmcp = float(maxit->first);
      purity = ndmcp / (float(simvec.size()) + float(bkgcount));
    }
    purityMCP pack;
    pack.mcp = dominantmcp;
    pack.purity = purity;
    return pack;
  }
};

#endif
