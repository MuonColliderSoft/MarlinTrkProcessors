#ifndef DDSpacePointBuilder_h
#define DDSpacePointBuilder_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <IMPL/TrackerHitImpl.h>
#include <map>
#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

// DD4HEP classes
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/ISurface.h"
#include "DDRec/SurfaceHelper.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Vector2D.h"
#include "DDRec/Vector3D.h"

using namespace lcio;
using namespace marlin;

/** ================= FTD Space Point Builder =================
 *
 * Builds space points for pairs of silicon strip detectors.
 * The digitisers create TrackerHitPlanars for the front and the back strips.
 * In order to get a spacepoint (as is needed by track reconstruction) those strip measurements need to be combined into
 * one space point.
 *
 * This is done by this processor.
 *
 *  <h4>Input - Prerequisites</h4>
 *
 * The TrackerHitPlanars as created by the Digitisers of FTD, SET or SIT.
 * This could of course be used for other detectors as well, but as information about the detector
 * is acquired from DDRec, and as different detectors (at the moment) are stored differently
 * in DDRec, used detectors have to be taken care of in the code.
 *
 *  <h4>Output</h4>
 *
 * A collection of TrackerHits containing the constructed space points.
 * The created hits will store the original strip hits in their rawHits.
 *
 * @param TrackerHitCollection The name of the input collection of TrackerHits coming from strip detectors on the FTD,
 * SIT or SET <br> (default name FTDTrackerHits) <br>
 *
 * @param SpacePointsCollection The name of the output collection of the created spacepoints <br>
 * (default name FTDSpacePoints) <br>
 *
 * @param TrackerHitSimHitRelCollection The name of the input collection of the relations of the TrackerHits to
 * SimHits<br> (default name FTDTrackerHitRelations)<br>
 *
 * @param SimHitSpacePointRelCollection The name of the SpacePoint SimTrackerHit relation output collection <br>
 * (default name VTXTrackerHitRelations) <br>
 *
 * @author Robin Glattauer HEPHY, Vienna
 *
 */
class DDSpacePointBuilder : public Processor {
public:
  virtual Processor* newProcessor() { return new DDSpacePointBuilder; }

  DDSpacePointBuilder();

  DDSpacePointBuilder(const DDSpacePointBuilder&) = delete;
  DDSpacePointBuilder& operator=(const DDSpacePointBuilder&) = delete;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

protected:
  /** Input collection name.
   */
  std::string _TrackerHitCollection{};

  /** Input relation collection name.
   */
  std::string _TrackerHitSimHitRelCollection{};

  /** Output collection name.
   */
  std::string _SpacePointsCollection{};

  /** Output relations collection name.
   */
  std::string _relColName{};

  /** Calculates the 2 dimensional crossing point of two lines.
   * Each line is specified by a point (x,y) and a direction vector (ex,ey).
   *
   * @return 0, if the calculation has been successful, another number if there was a problem.
   */
  static int calculateXingPoint(double x1, double y1, float ex1, float ey1, double x2, double y2, float ex2, float ey2,
                                double& x, double& y);

  /** Calculates a point between two lines in a way that it is equidistant from both lines and as close to them as
   * possible.
   * So when you calculate the distance between two lines, and you actually draw this distance as a short connecting
   * line, the point in the middle of this line, is what this method finds.
   * \verbatim
                  /
                 /
                /|___here!!!
               / |
     ---------/-------------- (line a)
             /
            /(line b)

     \endverbatim
   * (If this sketch looks terrible, you probably read this in doxygen and it still needs to be formatted in a
   * way, that it still looks sensible in doxygen.)
   *
   * @param P1 a point on the first line
   * @param V1 the direction vector of the first line
   *
   * @param point the reference to a point where the result will be stored
   */
  static int calculatePointBetweenTwoLines(const CLHEP::Hep3Vector& P1, const CLHEP::Hep3Vector& V1,
                                           const CLHEP::Hep3Vector& P2, const CLHEP::Hep3Vector& V2,
                                           CLHEP::Hep3Vector& point);

  /** Calculates the intersection of a line L and line PAPB ( strip on the first sensor ).
   *  L is constrained to pass through the point "Vertex" and bisect both lines PAPB ( strip on the first sensor ) and
   * PCPD ( strip on the second sensor )
   * @param PA start of the first line
   * @param PB end of the first line
   * @param PC start of the second line
   * @param PD end of the second line
   * @param Vertex the position
   *
   * @param point the reference to a point where the result will be stored
   */

  static int calculatePointBetweenTwoLines_UsingVertex(const CLHEP::Hep3Vector& PA, const CLHEP::Hep3Vector& PB,
                                                       const CLHEP::Hep3Vector& PC, const CLHEP::Hep3Vector& PD,
                                                       const CLHEP::Hep3Vector& Vertex, CLHEP::Hep3Vector& point);

  /** @return a spacepoint (in the form of a TrackerHitImpl* ) created from two TrackerHitPlane* which stand for
   * si-strips */

  TrackerHitImpl* createSpacePoint(TrackerHitPlane* a, TrackerHitPlane* b, double stripLength);

  /** @return the CellID0s of the sensors that are back to back to a given front sensor. If the given sensor
   * is in the back itself or has no corresponding sensor(s) on the back the vector will be empty.
   *
   * @param cellID0 a CellID0 corresponding to a sensor
   */
  std::vector<int> getCellID0sAtBack(int cellID0);

  std::vector<int> getCellID0sAtBackOfFTD(int cellID0);

  std::vector<int> getCellID0sAtBackOfSET(int cellID0);

  std::vector<int> getCellID0sAtBackOfSIT(int cellID0);

  /** @return information about the contents of the passed CellID0 */
  std::string getCellID0Info(int cellID0);

  int _nRun{};
  int _nEvt{};

  unsigned _nOutOfBoundary{};
  unsigned _nStripsTooParallel{};
  unsigned _nPlanesNotParallel{};

  float _nominal_vertex_x{};
  float _nominal_vertex_y{};
  float _nominal_vertex_z{};

  CLHEP::Hep3Vector _nominal_vertex{};

  float _striplength_tolerance{};

  double _striplength{};
  std::string _subDetName{};

  // dd4hep::Detector& lcdd;
  const dd4hep::rec::SurfaceMap* surfMap{nullptr};
};

#endif
