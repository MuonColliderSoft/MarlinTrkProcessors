#ifndef __MCPMAP__h
#define __MCPMAP__h 1

#include "lcio.h"
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <map>
#include <vector>

typedef std::vector<SimTrackerHit*> SimTrackerHitVec;
typedef std::multimap<MCParticle*, SimTrackerHit*> SegmentMap;
typedef SegmentMap::iterator segIter;
typedef std::map<MCParticle*, SimTrackerHitVec> MCPMap;
typedef std::pair<MCParticle*, SimTrackerHitVec> MCPMapPair;

// 以下のクラスはメソッドしかない。namespaceをクラス内で定義しちゃいけないらしく仕方なくこうした。
// 主な使用用途は、MCParticleとSimTrackerHitのLCRelationを作りたい場合などに有効である。
//
class moriUTIL {
public:
  moriUTIL() { ; }
  // もしこのクラスの外で以下のtypedefしてないならコメント外す

  /***以下は内部処理用の関数。ユーザーはMakeSegmentMap(Mが大文字なのに注意)を使えばそれで済む。*****/
  /***MakeSegmentMapの定義は一番下にある。**********************************************************/
  // 内部処理用
  SegmentMap makeSegmentMap(SimTrackerHitVec simthits) {
    int nsth = simthits.size();
    SegmentMap map;
    for (int i = 0; i < nsth; i++) {
      SimTrackerHit* sth = simthits[i];
      map.insert(std::make_pair(sth->getMCParticle(), sth));
    }
    return map;
  }

  // for inner process
  MCParticleVec getMCParticlesFromSegmentMap(SegmentMap map) {
    segIter it = map.begin();
    MCParticleVec mcps;
    mcps.push_back(it->first);
    it++;
    while (it != map.end()) {
      if (mcps.back() != it->first) {
        mcps.push_back(it->first);
      }
      it++;
    }
    return mcps;
  }

  // for inner process
  SimTrackerHitVec getSimTrackerHitVecFromSegmentMap(MCParticle* mcp, SegmentMap map) {
    SimTrackerHitVec simthits;
    simthits.reserve(10); // だいたい1MCPにつき10個以内のヒットがあると予想。
                          // push_backによるcapacityの自動調整を少なくする目的でつけた。
    typedef std::pair<segIter, segIter> EqualRange;
    EqualRange ret = map.equal_range(mcp);
    for (segIter it = ret.first; it != ret.second; it++) {
      simthits.push_back(it->second);
    }
    return simthits;
  }

  // for inner process
  MCPMap makeMCPMap(SegmentMap map) {
    MCParticleVec mcps = getMCParticlesFromSegmentMap(map);
    int size = mcps.size();
    MCPMap mcpmap;
    for (int i = 0; i < size; i++) {
      SimTrackerHitVec simthits = getSimTrackerHitVecFromSegmentMap(mcps[i], map);
      mcpmap.insert(std::make_pair(mcps[i], simthits));
    }
    return mcpmap;
  }

  //////////////////////////////////////
  // This is the tool for users.
  //////////////////////////////////////
  MCPMap MakeMCPMap(SimTrackerHitVec simthits) {
    if (simthits.size() == 0) {
      std::cout << "Fatal Error!! (MakeMCPMap cannot do with SimTrackerHitVec with size = 0 )" << std::endl;
      MCPMap empty;
      empty.clear();
      return empty;
    }
    SegmentMap segmap = makeSegmentMap(simthits);
    MCPMap mcpmap = makeMCPMap(segmap);
    return mcpmap;
  }
  // 例えば、SimTrackerHit Collectionにある全ヒットをここにぶち込めばMCPMapができあがる。
  // あるいは、Trackにアサインされたヒットのsimtrackerhitだけをぶち込めばtrackのヒットpurityや
  // picking efficiencyなんかも計算できるようになる。
  // あとは、LCRelationを作るのもいいし、ROOT解析のために使っても良い。
}; // moriUtil ends

#endif
