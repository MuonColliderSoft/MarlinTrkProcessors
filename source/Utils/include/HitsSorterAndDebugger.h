#ifndef HitsSorterAndDebugger_h
#define HitsSorterAndDebugger_h 1

#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>

// CxxUtils
#include "fpcompare.h"

// sorting by value of R(=x^2+y^2) in global coordinated so the hits are always
// sorted from close to the IP outward
inline bool sort_by_radius(EVENT::TrackerHit* hit1, EVENT::TrackerHit* hit2) {
  double radius1 =
      sqrt((hit1->getPosition()[0]) * (hit1->getPosition()[0]) + (hit1->getPosition()[1]) * (hit1->getPosition()[1]));
  double radius2 =
      sqrt((hit2->getPosition()[0]) * (hit2->getPosition()[0]) + (hit2->getPosition()[1]) * (hit2->getPosition()[1]));
  return CxxUtils::fpcompare::less(radius1, radius2);
}

// sorting by absolute value of Z so the hits are always sorted from close to
// the IP outward. This works as long as all hits are either in positive or
// negative side
inline bool sort_by_z(EVENT::TrackerHit* hit1, EVENT::TrackerHit* hit2) {
  const double z1 = fabs(hit1->getPosition()[2]);
  const double z2 = fabs(hit2->getPosition()[2]);
  return CxxUtils::fpcompare::less(z1, z2);
}

inline void printHits(const TrackerHitVec& hitVector) {
  int nHits = hitVector.size();
  for (int itHit = 0; itHit < nHits; itHit++) {
    // Get the tracker hit and print global coordinates of the hit
    TrackerHitPlane* hit = static_cast<TrackerHitPlane*>(hitVector.at(itHit));
    streamlog_out(DEBUG5) << " Hit #" << itHit << ", (x,y,z) = (" << hit->getPosition()[0] << ","
                          << hit->getPosition()[1] << "," << hit->getPosition()[2] << ")" << std::endl;
  }
  return;
}

// Print out the hits belonging to the track
inline void printHits(const Track* track) {
  const TrackerHitVec& hitVector = track->getTrackerHits();
  printHits(hitVector);
  return;
}

#endif
