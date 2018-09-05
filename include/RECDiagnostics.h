#ifndef RECDiagnostics_h
#define RECDiagnostics_h 1

#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/LCIterator.h"
#include "UTIL/Operators.h"
#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>

#include <marlin/VerbosityLevels.h>
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include <marlin/AIDAProcessor.h>
//#include <marlinutil/HelixClass.h>

#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/ICloud2D.h>
#include <AIDA/IHistogram2D.h>

#include <DD4hep/Detector.h>
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetType.h>
#include <DDRec/DetectorData.h>
#include <DD4hep/DetectorSelector.h>

#include <TROOT.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TEfficiency.h>


using namespace lcio ;
using namespace marlin ;
using namespace std ;

class TH1 ;

/** Extract various level track information.
 *  Creates ROOT histograms.
 *  For checking/monitoring the ILD optimisation.
 */

class RECDiagnostics : public Processor {
  
 public:
 
  virtual Processor*  newProcessor() { return new RECDiagnostics ; }
  
  RECDiagnostics(const RECDiagnostics&) = delete;
  RECDiagnostics& operator=(const RECDiagnostics&) = delete;

  RECDiagnostics() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

  /** Initialise the ROOT TTRee.
   */
  void initTree(void) ;

  /** Initialise the ROOT TH1F and so on ...
   */
  void initHist(void) ;

  /** Initialise the all the vectors.
   */
  void clearVec(void) ;

  /** Initialise the all the parameters for marlin.
   */
  void initParameters(void) ;

  
  //truth 

  void TrackGetSource(std::vector<Track*> &source, std::vector<std::vector<MCParticle*> >  &to, LCRelationNavigator* &relation); 
  
  void TrackGetSource2(std::vector<ReconstructedParticle*> &source, std::vector<std::vector<MCParticle*> >  &to, LCRelationNavigator* &relation); 

 protected:

  /** Input collection name.
   */
  std::string _colNameMCP{};
  std::string _colNameTrack{};
  std::string _colNamerec{};
  std::string _colNamerelaMcpRec{};
  std::string _colNamerelaRecMcp{}; 
  std::string _colNamerelaMcpTrk{};
  std::string _colNamerelaTrkMcp{};
  int _nRun{};
  int _nEvt{};

  float _bField{};

  std::vector<TH1*> _hReco{};

  TEfficiency* eff = nullptr;
  TEfficiency* eff_rec = nullptr;
  TEfficiency* eff_rec_type = nullptr;
  TEfficiency* eff_trk = nullptr;
  TEfficiency* eff2D = nullptr;
  TEfficiency* eff2D_pan = nullptr;
  TEfficiency* eff2D_pan_type = nullptr;
 
} ;


#endif
