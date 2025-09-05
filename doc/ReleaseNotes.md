# v02-12-07

* 2025-03-21 jmcarcell ([PR#73](https://github.com/iLCSoft/MarlinTrkProcessors/pull/73))
  - Set the subdetectorHitNumbers in ClonesAndSplitTracksFinder that was not being set (defaulted to 0) when `mergeSplitTracks` is false.

# v02-12-06

* 2024-08-28 Victor Schwan ([PR#72](https://github.com/iLCSoft/MarlinTrkProcessors/pull/72))
  - Remove the hardcoded (ILD) detector names from the FullLDC tracking to allow its usage with the `v11` model of ILD.
    - Keep defaults compatible with the non-hybrid ILD models

* 2024-08-23 tmadlener ([PR#70](https://github.com/iLCSoft/MarlinTrkProcessors/pull/70))
  - Add the `.clang-format` default file that is also used by other Key4hep repositories
  - Add a basic `pre-commit` configuration to run formatting checks
  - Run `clang-format` to fix all formatting issues

* 2024-08-22 tmadlener ([PR#71](https://github.com/iLCSoft/MarlinTrkProcessors/pull/71))
  - Fix warnings to make Key4hep based CI workflows pass again
  - Remove CentOS7 from Key4hep based workflows since it is no longer supported

# v02-12-05

* 2024-02-16 tmadlener ([PR#69](https://github.com/iLCSoft/MarlinTrkProcessors/pull/69))
  - Use the `LCRelationNavigator` to create output relation collections with the necessary type information.
    - This is necessary for automatic conversion from LCIO to EDM4hep when used within the MarlinWrapper and Gaudi.

* 2024-02-16 tmadlener ([PR#68](https://github.com/iLCSoft/MarlinTrkProcessors/pull/68))
  - Switch to the `key4hep-build` github action for building against key4hep stacks
  - Update the clicdp nightlies based workflow to the latest available version of the nightlies

# v02-12-04

* 2023-06-12 tmadlener ([PR#66](https://github.com/iLCSoft/MarlinTrkProcessors/pull/66))
  - Switch from `dd4hep::long64` to the more appropriate `dd4hep::CellID` after [AIDASoft/DD4hep#1125](https://github.com/AIDASoft/DD4hep/pull/1125)
  - Update CI
    - Switch to latest versions of github actions
    - Remove gcc8 based workflows
    - Add key4hep based workflow
  - Remove unused and deprecated usage of `std::binary_function`

# v02-12-03

* 2022-12-02 Thomas Madlener ([PR#64](https://github.com/iLCSoft/MarlinTrkProcessors/pull/64))
  - Remove the obsolote gcc8 based workflow and update github actions to the latest versions

# v02-12-02

* 2022-11-23 Thomas Madlener ([PR#62](https://github.com/iLCSoft/MarlinTrkProcessors/pull/62))
  - Make sure that the `_histos` pointer is at least initialized to a `nullptr` to avoid a spurious seg fault when trying to delete it uninitialized.

# v02-12-01

* 2022-06-28 Thomas Madlener ([PR#61](https://github.com/iLCSoft/MarlinTrkProcessors/pull/61))
  - Make doxygen cmake config work with newer cmakes

# v02-12

* 2021-11-05 Bohdan Dudar ([PR#55](https://github.com/iLCSoft/MarlinTrkProcessors/pull/55))
  - Center of the SET strips is now properly in the middle of the sensor modules.
  - Removed unused `vertex_x`, `vertex_y`, `vertex_z` steering parameters.

* 2021-09-22 Nazar Bartosik ([PR#53](https://github.com/iLCSoft/MarlinTrkProcessors/pull/53))
  - Added time smearing to the DDPlanarDigiProcessor 
  - Introduced a time window for digitized hits in the DDPlanarDigiProcessor 
  - Introduced TOF compensation for digitized hits in the DDPlanarDigiProcessor

* 2021-09-02 Placido Fernandez Declara ([PR#49](https://github.com/iLCSoft/MarlinTrkProcessors/pull/49))
  - DDPlanarDigiProcessor:  Add setting of the fromType and toType parameters for the LCRelation collection between TrackerHit and SimTrackerHit

* 2021-08-23 Andre Sailer ([PR#51](https://github.com/iLCSoft/MarlinTrkProcessors/pull/51))
  - CI: build against LCG_99python2 gcc8 and LCG_100 gcc10, clang11

* 2020-06-02 Erica Brondolin ([PR#45](https://github.com/iLCSoft/MarlinTrkProcessors/pull/45))
  - ClonesAndSplitTracksFinder 
    - If mergeSplitTracks set to `true`, pairs of tracks above a minimum pT requirements are compared in terms of pt, theta and phi
    - For each pair of reconstructed tracks: calculate significance for each parameter, impose a cut on the significance for each parameter, remove possible common hits and sort them in r, fit the new merged track
    - Results show positive effect on both single muons and complex events, with no significant increase of the run time ([Link to the slides](https://indico.cern.ch/event/920887/))

* 2020-01-27 Erica Brondolin ([PR#44](https://github.com/iLCSoft/MarlinTrkProcessors/pull/44))
  - In case finaliseLCIOTrack fails in the RefitFinal, the track should not be saved

# v02-11

* 2019-08-21 Frank Gaede ([PR#43](https://github.com/iLCSoft/MarlinTrkProcessors/pull/43))
  - fix RefitProcessor:
         - do not drop tracks in RefitProcessor, even if the refit failed
         - ensures that the output collection has the same number of entries
            in the same order as the input collection
         - needed to identify and compare refitted tracks

* 2019-04-01 Emilia Leogrande ([PR#42](https://github.com/iLCSoft/MarlinTrkProcessors/pull/42))
  - RefitFinal: added processor parameter `MinClustersOnTrackAfterFit`
    - defines the minimum number of hits on track after the fit for the track to be accepted 
    - previously hard-coded = 3, but tracks with 3 hits are ~100% fakes
    - now implemented as in ConformalTracking ilcsoft/ConformalTracking#52

* 2019-03-20 Emilia Leogrande ([PR#41](https://github.com/iLCSoft/MarlinTrkProcessors/pull/41))
  ClonesAndSplitTracksFinder:
    - The clone treatment has been uniformized with the updated clone treatment in the seeding (conformal tracking)
    - In few cases it is just a matter of including the equal sign in the comparison
    - In few other cases, the combined length-chi2 check is replaced with a length check only: longer tracks always preferred

* 2019-02-18 Shaojun Lu ([PR#39](https://github.com/iLCSoft/MarlinTrkProcessors/pull/39))
  - Improved SiliconTracking_MarlinTrk - added two options:
       - to use simple updated triplets searching core bin,
       - to use simple "AttachHitToTrack" for merging split track segments.

* 2019-01-18 Frank Gaede ([PR#40](https://github.com/iLCSoft/MarlinTrkProcessors/pull/40))
  - fixed DDTPCDigiProcessor
      - fix B-field correction for the r-phi resolution
      - now use (4./B)^2 rather than 4./B

* 2018-04-23 Emilia Leogrande ([PR#38](https://github.com/iLCSoft/MarlinTrkProcessors/pull/38))
  - ClonesAndSplitTracksFinder: the check for toBeMerged was in the wrong place. Fixed.

* 2018-04-20 Emilia Leogrande ([PR#37](https://github.com/iLCSoft/MarlinTrkProcessors/pull/37))
  - ClonesAndSplitTracksFinder: added mergeSplitTracks parameter, that allows one to decide whether to run only the clones skimming or also the merging of split tracks
    - mergeSplitTracks = false by default
    - mergeSplitTracks = true would enable the treatment of split tracks, which is still work in progress

# v02-10

* 2018-03-02 Emilia Leogrande ([PR#34](https://github.com/iLCSoft/MarlinTrkProcessors/pull/34))
  - New Processor: ClonesAndSplitTracksFinder for running after ConformalTracking
    - This processor has a track collection as input, removes all possible clones (tracks sharing at least 2 hits), looks for split tracks (tracks in a certain space region with similar track parameters) and merges them. If more than two tracks happen to fulfil the requirements, no merging is applied at the moment.

* 2018-02-25 Ete Remi ([PR#33](https://github.com/iLCSoft/MarlinTrkProcessors/pull/33))
  - DDPlanarDigiProcessor
      - Added time information to tracker hit 
  - DDSpacePointBuilder
      - Added time information to space point as the earliest time of the two tracker hits

* 2018-03-28 Marko Petric ([PR#36](https://github.com/iLCSoft/MarlinTrkProcessors/pull/36))
  - Fix for the removal of DDSurfaces which have been merged into DDRec 
    -  includes from `DDSurfaces` -> `DDRec`
    - namespace `DDSurfaces` -> `dd4hep::rec`

# v02-09-01

* 2017-10-18 Frank Gaede ([PR#31](https://github.com/iLCSoft/MarlinTrkProcessors/pull/31))
  - fix DDTPCDigiProcessor
      - allow smearing beyond cathode boundaries
  - reduce verbosity in ExtrToSIT

# v02-09

* 2017-10-12 Frank Gaede ([PR#30](https://github.com/iLCSoft/MarlinTrkProcessors/pull/30))
  - set correct MarlinTrkSystem configuration for every event with `MarlinTrk::TrkSysConfig`
     in all tracking processors

* 2017-07-07 Andre Sailer ([PR#24](https://github.com/iLCSoft/MarlinTrkProcessors/pull/24))
  - SiliconTracking_MarlinTrk: initialise fastfitter in init function, fixes small memory leak if processor is not active
  
  - FPCCDFullLDCTracking_MarlinTrk: initalise helper objects in init and add cleanup in end. Fixes small memory leaks

* 2017-10-09 Frank Gaede ([PR#29](https://github.com/iLCSoft/MarlinTrkProcessors/pull/29))
  - fix DDTPCDigiProcessor
      - ensure hits are not smeared beyond cathode

* 2017-08-21 Andre Sailer ([PR#25](https://github.com/iLCSoft/MarlinTrkProcessors/pull/25))
  - RefitFinal: new processor for refitting track collections. Does not re-sort hits (like RefitProcessor), different way to setup the initial trackstate for the fit, and different MarlinTrkUtil function being called
  - DDPlanarDigi: add absolute values of u and v smearing histograms

* 2017-08-23 Emilia Leogrande ([PR#26](https://github.com/iLCSoft/MarlinTrkProcessors/pull/26))
  * TruthTrackFinder: m_maxChi2perHit changed from 1.e3 to 1.e2 for consistency with Conformal Tracking

* 2017-09-12 Shaojun Lu ([PR#27](https://github.com/iLCSoft/MarlinTrkProcessors/pull/27))
  - DDTPCDigiProcessor: Added one control for the input value of asin to avoid NAN.

* 2017-10-06 Andre Sailer ([PR#28](https://github.com/iLCSoft/MarlinTrkProcessors/pull/28))
  - Drop unused and no longer existing header includes AidaSoft/DD4hep#241

# v02-08-01

* 2017-07-06 Frank Gaede ([PR#23](https://github.com/iLCSoft/MarlinTrkProcessors/pull/23))
  - fix phi0 of module rings in TPCModularEndplate

# v02-08

* 2017-06-28 Frank Gaede ([PR#21](https://github.com/iLCSoft/MarlinTrkProcessors/pull/21))
  - remove hits that are in module gaps on the TPC endplate in DDTPCDigiProcessor
        - use the following parameters to define the endplate:
        - TPCEndPlateModuleNumbers
        - TPCEndPlateModulePhi0s
        - TPCEndPlateModuleGapPhi
  - rm obsolete parameter DontEncodeSide from DDTPCDigiProcessor

* 2017-06-15 Marko Petric ([PR#14](https://github.com/iLCSoft/MarlinTrkProcessors/pull/14))
  - Replace auto_ptr with unique_ptr
  - explicitly call copy constructor of base class

* 2017-06-15 Emilia Leogrande ([PR#13](https://github.com/iLCSoft/MarlinTrkProcessors/pull/13))
  * TruthTrackFinder: hits in the same layer of the same subdetectors are now filtered. Those with higher radius are removed, then the fit is performed. In case the fit fails, those with higher z are removed and a new fit attempt is made.

* 2017-06-14 Frank Gaede ([PR#12](https://github.com/iLCSoft/MarlinTrkProcessors/pull/12))
  - remove some obsolete processors:
          - PlanarDigiProcessor
          - SpacePointBuilder
          - SimpleCylinderDigiProcessor
          - SimpleDiscDigiProcessor
          - SimplePlanarDigiProcessor.h
  - replace Gear with DD4hep/DDRec in all remaining processors

* 2017-06-30 Andre Sailer ([PR#22](https://github.com/iLCSoft/MarlinTrkProcessors/pull/22))
  - ExtrToTracker: cleanup transient collection

* 2017-06-20 Frank Gaede ([PR#17](https://github.com/iLCSoft/MarlinTrkProcessors/pull/17))
  - fix leftover namespace changes in dd4hep

* 2017-06-20 Andre Sailer ([PR#16](https://github.com/iLCSoft/MarlinTrkProcessors/pull/16))
  - Adapt to changes in namespaces and LCDD -->  Detector

* 2017-06-27 Andre Sailer ([PR#20](https://github.com/iLCSoft/MarlinTrkProcessors/pull/20))
  - FPCCDSiliconTracking_MarlinTrk: move object creation to init, add delete
  - DDPlanarDigiProcessor: free random number generator at the end

* 2017-06-24 Frank Gaede ([PR#19](https://github.com/iLCSoft/MarlinTrkProcessors/pull/19))
  - add DDTPCDigiProcessor and helper classes
          - adaptation of MarlinReco/TPCDigiProcessor for DD4hep

* 2017-06-23 Andre Sailer ([PR#18](https://github.com/iLCSoft/MarlinTrkProcessors/pull/18))
  - DDCellsAutomatonMV: fix duplicate parameter name 
     * MVHitsThetaDifference --> MVHitsThetaDifference_Adjacent

# v02-07

* 2017-04-28 Frank Gaede ([PR#11](https://github.com/iLCSoft/MarlinTrkProcessors/pull/11))
  - TruthTrackFinder: Hits produced by secondaries are discarded

* 2017-04-28 Andre Sailer ([PR#10](https://github.com/iLCSoft/MarlinTrkProcessors/pull/10))
  - TruthTrackFinder::sort_by_z: sort by absolute z. Fixes an issue with wrongly assigned charge and angle for tracks in the backward direction if the fit with hits sorted by radius failed.

# v02-07

* 2017-04-28 Frank Gaede ([PR#11](https://github.com/iLCSoft/MarlinTrkProcessors/pull/11))
  - TruthTrackFinder: Hits produced by secondaries are discarded

* 2017-04-28 Andre Sailer ([PR#10](https://github.com/iLCSoft/MarlinTrkProcessors/pull/10))
  - TruthTrackFinder::sort_by_z: sort by absolute z. Fixes an issue with wrongly assigned charge and angle for tracks in the backward direction if the fit with hits sorted by radius failed.

# v02-06

* 2017-04-06 Emilia Leogrande ([PR#9](https://github.com/iLCSoft/MarlinTrkProcessors/pull/9))
  - Replaced ILDCellID0 with LCTrackerCellID
  - Replaced DDKalTestConf with LCTrackerConf

* 2017-04-07 Andre Sailer ([PR#8](https://github.com/iLCSoft/MarlinTrkProcessors/pull/8))
  - Fix many warnings for gcc
  - ExtrToSIT::SelectBestCandidateLimited: remove variable shadowing pointer parameter
  - ExtrToTracker: Fix bug when using performFinalRefit where the refitted track was not used after refit

#  v02-05

- DDPlanarDigitiser: Use pre-filled map of neighboring surfaces to speed up extrapolation of tracks (R. Simoniello)
- DDPlanarDigitiser: Add parameter MinimumEnergyPerHit to set a threshold for accepting hits (A. Sailer)
- TruthTrackFinder: Fix duplicate trackstate problem from truthTrackFinder (A. Sailer)
- cellID encoding string from singleton in DDKalTest (momentarily added depedence on DDKalTest)
- add TrkSystemName parameter to ExtrToSIT processor

#  v02-04-01 

- FPCCDSiliconTracking_MarlinTrk.cc:
	- patch from K.Fujii for macos (use std::isnan)

  
# v02-04
  

  F.Gaede
   - add parameter ParticleMass to RefitProcessor
      - can fit with different mass hypotheses
   

  Y. Voutsinas 
   - added DDSpacePointBuilder :  using DDRec rather than GEAR 
      -> to become the new default

   - improved CellsAutomatonMV


  R. Simoniello
   - cleaned up ExtrToTracker and DDCellsAutomatonMV ( uss DDRec rather than gear, ...)



  
# v02-03
  

  R. Simoniello
 
 - Add procedure to handle with spiral tracks. Fit direction and fit initialisation configurable by steering macro (default fit direction: backward)
  - improved speed, map of hits on each element introduced, find neighbours of element ID and extrapolate to them (to be improved to some problem in the geo datat structure)	

  F. Gaede
  - made compatible with c++11
	- removed -ansi -pedantic -Wno-long-long
	- fixed narrowing in initializer lists

  Y. Voutsinas
  - Minivector creation from 1D SIT hits added 


# v02-02
  
  R.Simoniello
   - RefitProcessor
     - for each subdetctor X, store nhits in fit in lcio::ILDDetID::X -1 and nhits 
       in track in lcio::ILDDetIDX -2. Structure needed for Pandora.
   - file added for correct comparison for float and double: fpcompare.h 
   -DDCellsAutomatonMV/ ExtrToTracker
     - criteria based on nhits ans QI for CA vertex tracks, remove the possibility of geo from gear file (option GEO1 in compilation) 
       for extrapolated tracks

  F.Gaede
  - DDPlanarDigiProcessor:
   - introduced parameter ForceHitsOntoSurface:
     - if true hits are projected onto the planar surface
       if they are putside the bounds
     - no check on distance is applied !
    - store the resolution from the length of the wafer for 1D hits
      ( could be treated as 2d measurements )

  - RefitProcessor:
     - no longer write incomplete tracks (w/o) track state at calo
     - rely on MarlinTrkTrack's implementation of initialise() if
       no already fitted TrackState is used



  
# v02-01
  

   R.Simoniello
   - updated Refitting/src/ExtrToTracker.cc
     - extrapolation also in endcap subdetectors + CA selecting longest (nhits) tracks

   F.Gaede
   - handle empty collections gracefully ( for lcgeo/ddsim ) 
     in Digitisers/src/PlanarDigiProcessor.cc

   - use abstract ISurface


  
# v02-00  
  

   - adapted to  DD4hep::Surfaces and DDKalTest


   - DDCellsAutomatonMV: adapted CA to DD4hep (R.Simoniello)

   - modified RefitProcessor (F.Gaede)
    - InitialTrackState: allows to select one of the track states (or use hits) 
    - FitDirection: allow to fit forwards or backwards (default)
    -> can use for CLIC SiBarrel tracker to refit TrackState AtIP 
       in forward direction

   - adapted to new MarlinTrk::Factory: no longer delete the IMarlinTrkSystem pointer
     at the end of the processing for all fitters

   - ExtrToTracker: add new processor that extracts an existing track to another 
     tracking subdetector (R.Simoniello)

   - added DDPlanarDigiProcessor which uses the DDRec::Surface
     for smearing the hits (no gear needed) (F.Gaede)

   - add option to use DDKalTest to all fitters (TrackSubsetProcessor, FullLDCTracking_MarlinTrk,
      TrackSubsetProcessor, RefittProcessor)    (F.Gaede)
      - set parameter "TrackSystemName" to "DDKalTest"

   - PlanarDigiProcessor: add new utility processor that allows to split a collection 
    of hits into several collections, based on the layerID
    of the hits (works for all lcio hit types) (F.Gaede)

   - ExtrToSIT: added processor for refitting and propagating track from TPC or VXD to SIT (Y.Voutsinas)

   - reverted changes for dEdx calculationin FullLDCMarlinTrk: moved to new processor (M.Kurata)


   - new processors: 
      Refitting/src/TruthTrackFinder.cc (D.Hynds)
         alternative track cheating
      Utils/src/ClicEfficiencyCalculator.cc (D.Hynds)
         tracking efficiency for CLIC (all silicon)

  
# v01-11  
  
 -  added processor for VXD tracking using a cellular automaton algorithm based on mini - vectors


  
# v01-10  
  
 -  added FPCCDTracking code by Tatsuya Mori

# v01-09-01  

  - updated calling attrinbutes of ced_hit_ID to new CED # v01-09

  
# v01-09
  
  - FullLDCTracking_MarlinTrk
    - Use the Chi2 probability to remove badly fitted silicon tracks so that the don't pollute the TPC Si track merging.
    - Reject combinations of TPC and Silicon Tracks if more that a certain number of Silicon Hits, default 2, get rejected by the fit of the combined track.

  
# v01-08 
  

 - General 
   - SimplePlanarTestDigiProcessor renamed as PlanarDigiProcessor
   - Fixed signed unsigned warning.
   - Fix memory leaks.

 - PlanarDigiProcessor
   - Made processor parameters ResolutionU(V) vectors
   	 in order to allow for different point resolution values per layer
   	 ( needed for the VXD in the DBD)

 - FullLDCTracking_MarlinTrk
   - Added SET hits. 
   - Make sure hits flagged as not to be used in the fit are not included.
   - ForceTPCSegmentsMerging off by default
   - Ensure that trying to add a false hits, cannot derail the whole track creation. 
     Test number of Outliers using MaxAllowedPercentageOfOutliersForTrackCombination. 
   - When not using ForceTPCSegmentsMerging, make sure that TPC tracks which do not get 
   	 merged with Silicon Tracks have their hits set to setUsedInFit.


  
# v01-07 
  

  General:

   - Debug output has been made more consistent throughout.	

  FullLDCTracking_MarlinTrk
  
   - Use existing tracks parameters when refitting. 
   - Updated for new lcio TrackState copy constructor taking const reference. 
     
  TruthTracker	

   - Updated for new lcio TrackState copy constructor taking const reference. 
   - Corrected Helix orientation for pre-fit.


  RefitProcessor

   - Corrected Helix orientation for pre-fit.			
   - Protect against missing truth relations.

  
# v01-06-01 
  

  General:			

   - RIADA dependency added

  FullLDCTracking_MarlinTrk:

 - Increased debug output
 - Use fit from TPC to initialise fit when doing SiTrack TPCTrack comparisons.
 	Added debug output which can we removed at a later stage
 - Use Tracks from the TPC which are composed only of the innermost segment, which contains the fit needed 
 	for matching to the Silicon Tracks. The Silicon Track, TPC Track, and the segments are added to the final track, 
	which is then itself added to the Final Track collection. 
 - Added type bits for the final Track collection. 

 SiliconTracking_MarlinTrk:
    
- Make sure that float[6] is used for getPointOnCircle. 
- Added diagnostic histograms. 
- Check error return of fast fitter. 
- Use iopt 2 instead of 3, as newtonian part of fast fitter is not performing well.     
- Added type bits for the final Track collection. Increased debug for FTD tracking.
 
  SpacePointBuilder:
   
    - Added extra check for valid line intersection when creating space points. This fixes problems from very low pt tracks.


  
# v01-06 
   

  General:			 
  - Removed use of sort predicate. Use sorted list of std pairs instead.
  - Use _maxChi2PerHit in fits. 

  FullLDCTracking_MarlinTrk:	 
  - Corrected uninitialised covariance matrix in combine tracks method.
  - Fixed sorting hits in r2 during TPCTrack SiTracks combination.
  - Improved hit and geometry handling. 
  - Removed unused parameters. 
  - Fixed problem in allocating left over silicon hits, track fit was incorrect after adding a hit. 
  - Corrected point_res_rphi for FTD SpacePoints. 
  - Removed unused code.

  TruthTracker:	 
 - Added new mode to fit tracks iteratively, off by default. Allow the CED event display to be used during track creation.	 
 - For the iterative mode split Tracks now have the individual segments added to a separate collection. 
 - These segments are used to produce a combined track which uses the innermost and outermost segments to 
 - provided the track states at the IP, FirstHit, LastHit and Calo respectively. 
 - The combined track is added to the main track collection, the NDF and chi2 are set with the values from the innermost segment.
 - The NDF and chi2 values for the other segments can be retrieved from the tracks added to the composite track.  
	 
 -	Hits with no MCParticle are discarded.

  SiliconTracking_MarlinTrk:	 
  - Added protection against being overwhelmed by a very large number of hits within one sector. 
  - The default is set to 100. For now an ERROR message is printed and the sector is dropped.
  - Quality parameter added to output collection set to "Poor"
  - This will be removed when the problem is more properly addressed.
  - Make sure sort predicate for chi2 fulfils strict weak ordering.
  - Fixed problem where using return instead of continue caused loss of all fitted tracks when single track fit failed.
  - Improved hit and geometry handling. 
  - Instrumented with CEDEventDisplay to improve debugging. 
  - Removed unused parameters.

  SpacePointBuilder:		 
  - Corrected SpacePoint creation bug. Endpoint of rear strip swapped with endpoint of front strip. Added debug output for checking whether a strip is created from a real combination or ghost combination.

  TrackSubsetProcessor:		 
  - Add hit numbers and additional trackstates. Use MarlinTrk utilities. 

  
# v01-05 
   

 General:			 Added dependencies on KiTrack and KiTrackMarlin to MarlinTrkProcessors.
 
 Refitting:			 Added the TrackSubsetProcessor, to get the best subset from a redundant set of tracks. 
 				 TruthTracker restructured to use new convenience methods from MarlinTrk to produce 
				 cannonical trackstates. The diagonal elements of the Initial error matrix can now be set individually.
				 Fixed r_comparison sort problem.

 SiliconTracking_MarlinTrk	 Use convenience functions for track fitting from MarlinTrk. Provide steering of initial covariance matrix,  				 
& FullLDCTracking_MarlinTrk	 and max delta chi2 for adding hits. 
  				 Removed the creation of relations for MC <-> Tracks.

 SpacePointBuilder:		 Corrected bias in position by using vertex constraint	

  
# v01-04 
   

  SimplePlanarDigiProcessor:	 When calling getLadderNumber the "ladderNumber" was given as input instead of the "layerNumber" (typo).
  				 In getLadderNumber some hits were not assigned to the right layer due to rounding errors 
				 (they got INT_MAX which caused an exception when writing it to the bit field). Introduced a (hard coded) 1 nm rounding 
				 acceptance and all hits were found to be in their layer. In the case this should fail additionally put in that 
				 the layer with the closest distance to the layer centre is returned.
  				 Fixes from M. Killenberg 
   				
  SimplePlanarTestDigiProcessor: Updated to be able to work for all planar sub detectors, VXD, SIT, SET and FTD.
  				 Will very likely supperseed SimplePlanarDigiProcessor in the next release. 

  SiliconTracking_MarlinTrk &	 Added Strip Hits for FTD. 
  FullLDCTracking_MarlinTrk:     FTD restored to track finding and fitting.
 				 

  RefitProcessor:		 Updated to deal with recent addition of Strip Hits in FTD.

  TruthTracker:			 Allow Helix parameters for fit to be take from MCParticle using UseMCParticleParametersFotInitOfFit == true. 
  				 Allow steering of values used for the initial diagonal elements of the trackfit via InitialTrackErrors.

  
# v01-03 
   

  General:			Make use of MeasurementSurface classes from GEAR, as well as new UTIL::BitSet32 and 
   				UTIL::ILDTrkHitTypeBit from ILDConf (LCIO). 			

  SiliconTracking_MarlinTrk: 	Adapted to use composite Space-points composed from strip hits. 
     			    	Added new parameter NHitsChi2 which controls the maximal number of hits for which a track 
				with n hits is aways better than one with n-1 hits.
				General Cleanup and merger from changes from CLIC CDR. 				

  FullLDCTracking_MarlinTrk:	Adapted to use composite Space-points composed from strip hits.
	
  SimplePlanarDigiProcessor:	Adapted to have option to produce strip hits.

  TruthTracker:			Flexible collection input and able to construct track including composite Space-points.

  SpacePointBuilder:		New - Builds 3D composite Space-points from strip hits.



  
# v01-02 
   

SiliconTracking_MarlinTrk &	First version which is able to reconstruct LOI data.
 FullLDCTracking_MarlinTrk:     Fixed situation where, although more that 3 hits are excepted for fitting, 
 				less than 3 actually are included in the fit, by the kalman filter.

 SimpleCylinderDigiProcessor:	New - Uses simple gaussian smearing for cylinders.

 TruthTracker:			Fixed bug in SimTrackerHitSortPredicate.

 SimpleDiscDigiProcessor:	For LOI data the z coordinate is now set to that of the disk and the tolerance for dz has been increased to 10 microns.


  
# v01-01 
   

 SiliconTracking_MarlinTrk &   Modified to use getHitsInFit from MarlinTrk to get the correct TrackState for the first and last hits. 
 FullLDCTracking_MarlinTrk:    SiliconTracking has had the work around applied for the optimization problem of compare_r as applied 
                               to FullLDCTracking. 

  
# v01-00 
   

  SimplePlanarDigiProcessor:  
  			      Creates TrackerHits from SimTrackerHits, smearing them according to the input parameters. 
 			      The plannar geometry should be either VXD, SIT or SET described using ZPlannarLayout
 			      The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 			      of SimTrackerHits perpendicular and along the ladder according to the specified point resolutions.

  SimpleDiscDigiProcessor: 
  			      Produces a TrackerHit collection from SimTrackerHit collection. Intended for use with the FTD.
			      The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
			      of SimTrackerHits in (x,y) plane according to the specified point resolution. 
	
  SiliconTracking_MarlinTrk:  Processor performing stand-alone pattern recognition in the VTX, FTD and SIT.
  			      Adaption of the of the original SiliconTracking from A. Raspereza, to use MarlinTrk.

  FullLDCTracking_MarlinTrk:  Processor performing track finding procedure in the entire ILD detector by linking track segments found 
                              by the SiliconTracking module in the silicon detectors and by Clupatra in the TPC. 
			      Adaption of the of the original SiliconTracking from A. Raspereza, to use MarlinTrk.

  RefitProcessor:             Track Refitter processor for marlin. Refits an input track collection using MarlinTrk, 
                              producing a new collection of tracks. Primarily to demonstrate the use of MarlinTrk.

  TruthTracker:     	      Track creation based on MC truth. Uses Relations between TrackerHits and SimTrackerHits to produce a track 
  			      collection. It can also fit the tracks using MarlinTrk or fill the track collection with the helix parameters
                              of the MCParticle.			   
