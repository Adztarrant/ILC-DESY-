#include "RECDiagnostics.h"
#include <math.h>
#include <GeometryUtil.h>


RECDiagnostics aRECDiagnostics ;

// helper enum defining histogram index in vector 
namespace RecoHistos {
  enum index{
    hMCPMomentum,
    hRec,
    hTrkMomentum, 
    hPFOMomentum, 
    Momentumdif,
    ptdif,
    //Good,
    //-----  keep Size as last :
    Size   
  };

  
  class Histograms{
  public:

    //quick set up of histogram format 
    Histograms(std::vector<TH1*>& v) : _h(&v) {}
  
    void create(int idx, const char* n, int nBin=100, double min=0., double max=0. ){
      create( idx , n , n , nBin, min , max ) ; 
    }

    void create(int idx, const char* n, const char* t,  int nBin=100, double min=0., double max=0.){

      _h->at( idx ) = new TH1D( n, t , nBin , min, max ) ;
 
      streamlog_out( DEBUG ) << " create histo " <<  n << " at index " << idx << std::endl ;
    }
  
    void create(int idx, const char* n, const char* t,  int nBin , double* bins ){

      _h->at( idx ) = new TH1D( n, t , nBin , bins ) ;

      streamlog_out( DEBUG ) << " create histo " <<  n << " at index " << idx << std::endl ;
    }

    void fill( int idx , double val, double weight=1.0 ){  _h->at( idx )->Fill( val , weight ) ; }

  protected:

    std::vector<TH1*>* _h;

  };
}

using namespace RecoHistos ;

void RECDiagnostics::initParameters(void) {

  // register steering parameters: name, description, class-variable, default value
 
  registerInputCollection( LCIO::LCRELATION,
			   "track reco to truth info" , 
			   "Name of track reconstruction truth input collection"  ,
			   _colNamerelaTrkMcp ,
			   std::string("MarlinTrkTracksMCTruthLink")
			   ) ;
  
  registerInputCollection( LCIO::LCRELATION,
			   "truth to track reco info" , 
			   "Name of track reconstruction truth input collection"  ,
			   _colNamerelaMcpTrk ,
			   std::string("MCTruthMarlinTrkTracksLink")
			   ) ;
  
  registerInputCollection( LCIO::LCRELATION,
			   "reconstucted Particle to truth info" , 
			   "Name truth input collection for recon"  ,
			   _colNamerelaRecMcp ,
			   std::string("RecoMCTruthLink")
			   ) ;
  
  registerInputCollection( LCIO::LCRELATION,
			   "truth to reconstucted Particle info" , 
			   "Name truth input collection for recon"  ,
			   _colNamerelaMcpRec ,
			   std::string("MCTruthRecoLink")
			   ) ;
  

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection" , 
			   "Name of the MCParticle input collection"  ,
			   _colNameMCP ,
			   std::string("MCParticle")
			   ) ;

  registerInputCollection( LCIO::TRACK,
			   "StudiedTracks" , 
			   "Name of the FullLDC track collection"  ,
			   _colNameTrack ,
			   std::string("MarlinTrkTracks")
			   ) ;  
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "PFOCollection" , 
			   "Name of the ReconstructedParticle collection"  ,
			   _colNamerec,
			   std::string("ReconstructedParticle")
			   );
}


RECDiagnostics::RECDiagnostics() : Processor("RECDiagnostics") {

  _description = "RECDiagnostics generates information for tracking performance monitor." ;

  initParameters() ;

}


void RECDiagnostics::init() { 

  streamlog_out(DEBUG) << "   RECDiagnostics::init() called  "  << std::endl ;

  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

}

void RECDiagnostics::processRunHeader( LCRunHeader* run) { 
  streamlog_out(MESSAGE) << " start processing run "<<run->getRunNumber() << std::endl;
  _nRun++ ;
} 

void RECDiagnostics::processEvent( LCEvent * evt ) {

  streamlog_out(MESSAGE) << " start processing event "<<evt->getEventNumber() 
			 << "   in run:  " << evt->getRunNumber() << std::endl ;

  // define a histogram pointer

  static AIDA::ICloud2D* hpMCPvspTrk ;   
  static AIDA::ICloud2D* hpMCPvspPFO ;   
  static AIDA::ICloud2D* hpMCPvspRec ;  
  static AIDA::ICloud2D* Momdif ;
  static AIDA::ICloud2D* Ptdif ;
  //static AIDA::ICloud2D* good ;

  RecoHistos::Histograms h(_hReco) ;


  if( isFirstEvent() ) { 
    streamlog_out(DEBUG4) << " This is the first event " << std::endl;
  
    // this creates a directory for this processor ....

    AIDAProcessor::histogramFactory( this ) ;
      
    _hReco.resize( RecoHistos::Size ) ;
      
    const int nBins = 100 ;
    const float pmin = 0.0 ;
    const float pmax = 180.0 ;

    //defines all histograms made 
    h.create(hMCPMomentum  , "hMCPTheta ", "Theta of the MCParticles, muon gun simulation at 20GeV; #theta; Frequency", nBins, pmin, pmax ) ; 
    h.create(hTrkMomentum  , "hTrkTheta ", "Theta of the MarlinTrkTrack, muon gun simulation at 20GeV ; #theta; Frequency", nBins, pmin, pmax) ; 
    h.create(hRec , "hRecTheta ", "Theta of the ReconstructedParticles, muon gun simulation at 20GeV; #theta; Frequency", nBins, pmin, pmax) ; 
    h.create(Momentumdif  , "momdif", "Difference in theta of marlin trk and mcparticles muon gun simulation at 20GeV; #theta; Frequency", nBins, -0.5, 0.5 ) ; 
    h.create(ptdif  , "p_tdif", "Difference in transverse momentum of marlin trk  and mcparticles muon gun simulation at 20GeV; pt(GeV); Frequency", nBins, -0.5, 0.5 ) ; 
    //h.create(Good  , "Goodness of type casted", 20, 0, 1) ;


    hpMCPvspTrk  = AIDAProcessor::histogramFactory(this)->
      createCloud2D( "hpMCPvspTrk", "MCP vs MarlinTrkTrack", nBins ) ; 

    hpMCPvspPFO  = AIDAProcessor::histogramFactory(this)->
      createCloud2D( "hpMCPvspPFO", "MCP vs chagedPFO", nBins ) ; 

    hpMCPvspRec  = AIDAProcessor::histogramFactory(this)->
      createCloud2D( "hpMCPvspRec", "recon", nBins ) ; 
   
    hpMCPvspRec  = AIDAProcessor::histogramFactory(this)->
      createCloud2D( "hpMCPvspRec", "recon", nBins ) ; 

    Momdif = AIDAProcessor::histogramFactory(this)->
      createCloud2D( "Momdif", "mdif", nBins ) ; 
    Ptdif  = AIDAProcessor::histogramFactory(this)->
      createCloud2D( "ptdif", "pdif", nBins ) ; 
    // good  = AIDAProcessor::histogramFactory(this)->
    // createCloud2D( "Good", "Good", nBins ) ; 

    // defines all the Efficiency plots
    //defines range of plot change as required
    eff = new TEfficiency("eff","Efficiency of theta in track to mcparticle relation, muon gun simulation at 20GeV;#theta;efficiency",100,0.0,180.0);

    eff_rec = new TEfficiency("eff_rec","Efficiency of theta in reconstructed (PFO) to mcparticle relation, muon gun simulation at 20GeV;#theta;efficiency",100,0.0,180.0);
     eff_rec_type = new TEfficiency("eff_rec_t","Efficiency of theta in reconstructed (PFO) by type to mcparticle relation, muon gun simulation at 20GeV;#theta;efficiency",100,0.0,180.0);
   
    eff2D = new TEfficiency("eff_2D","theta and transverse momentum, efficiency track to mcparticle relation,muon gun simulation at 20GeV;#theta;pt(GeV)",100,0.0, 180.0, 100, 0., 25.); 
    eff2D_pan = new TEfficiency("eff_2D_pan","theta and transverse momentum, efficiency reconstructed(PFO) to mcparticle relation by charge, muon gun simulation at 20GeV;#theta;pt(GeV)",100,0.0, 180.0, 100, 0., 25.);
    eff2D_pan_type = new TEfficiency("eff_2D_pan_t","theta and transverse momentum, efficiency reconstructed(PFO) to mcparticle relation by type, muon gun simulation at 20GeV;#theta;pt(GeV)",100,0.0, 180.0, 100, 0., 25.);
 }

  // pointers to the hist 
  LCCollection* colMCP = nullptr ;
  LCCollection* colTrk = nullptr ;
  LCCollection* colrec = nullptr ;
  LCCollection* colrelaMcpTrk = nullptr ;   
  LCCollection* colrelaTrkMcp = nullptr ;
  LCCollection* colrelaMcpRec = nullptr ;   
  LCCollection* colrelaRecMcp = nullptr ;
  
  try{
    colMCP = evt->getCollection( _colNameMCP ) ;
  } catch(lcio::Exception){
    streamlog_out( DEBUG5 ) << " collection " << _colNameMCP <<  " not found ! "   << std::endl ;
  }

  try{
    colTrk = evt->getCollection( _colNameTrack ) ;
  } catch(lcio::Exception){
    streamlog_out( DEBUG5 ) << " collection " << _colNameTrack <<  " not found ! "   << std::endl ;
  }

  try{
    colrec = evt->getCollection( _colNamerec ) ;
  } catch(lcio::Exception){
    streamlog_out( DEBUG5 ) << " collection " << _colNamerec<<  " not found ! "   << std::endl ;
  }

  try{
    colrelaMcpTrk = evt->getCollection( _colNamerelaMcpTrk ) ;
  } catch(lcio::Exception){
    streamlog_out( DEBUG5 ) << " collection " << _colNamerelaMcpTrk <<  " not found ! "   << std::endl ;
  }

  try{
    colrelaTrkMcp = evt->getCollection( _colNamerelaTrkMcp ) ;
  } catch(lcio::Exception){
    streamlog_out( DEBUG5 ) << " collection " << _colNamerelaTrkMcp <<  " not found ! "   << std::endl ;
  }

  try{
    colrelaMcpRec = evt->getCollection( _colNamerelaMcpRec ) ;
  } catch(lcio::Exception){
    streamlog_out( DEBUG5 ) << " collection " << _colNamerelaMcpRec <<  " not found ! "   << std::endl ;
  }

  try{
    colrelaRecMcp = evt->getCollection( _colNamerelaRecMcp ) ;
  } catch(lcio::Exception){
    streamlog_out( DEBUG5 ) << " collection " << _colNamerelaRecMcp <<  " not found ! "   << std::endl ;
  }

  
  LCRelationNavigator *_navMcpTrk = new LCRelationNavigator( evt->getCollection( _colNamerelaMcpTrk ) );    //used to get truths to MarlinTrkTracks 
  LCRelationNavigator *_navTrkMcp = new LCRelationNavigator( evt->getCollection( _colNamerelaTrkMcp ) );    //used to get truths from MarlinTrkTracks 
  LCRelationNavigator *_navMcpRec = new LCRelationNavigator( evt->getCollection( _colNamerelaMcpRec ) );    //used to get truths to PFOs
  LCRelationNavigator *_navRecMcp = new LCRelationNavigator( evt->getCollection( _colNamerelaRecMcp ) );    //used to get truths from PFOs


  //bools useed to make cuts 
  bool istable = false;
  bool trkFound = false;
  bool recFound = false;
  bool recFoundCharge = false;
  bool recFoundType = false;

 
  //vectors used to store truth info
  std::vector<MCParticle* > origin_source;
  std::vector<Track*> recon_source;
  std::vector<ReconstructedParticle*> pando_source;
  std::vector<std::vector<MCParticle*> > original_source; 

  ////////////////////////plots/////////////////////////

  //  if( colMCP != NULL && colMCP->getTypeName() == LCIO::MCPARTICLE ){
  //    if( colTrk != NULL && colTrk->getTypeName() == LCIO::TRACK ){
  //     if( colrec != NULL && colrec->getTypeName() == LCIO::RECONSTRUCTEDPARTICLE){
      

	////////////////////// MCPARTICLE /////////////////////////////
	//looking for MCParticles nominator, when they pass selections/cuts
	for( int i=0, n=colMCP->getNumberOfElements(); i<n ; i++ ){
	    
	  MCParticle* mcp = dynamic_cast<MCParticle*>( colMCP->getElementAt(i) ) ;

	  const int stable = mcp->getGeneratorStatus();
	  const bool isdecay = mcp->isDecayedInTracker(); 
	
	  //cuts on stability and charge 
	  if( std::fabs( mcp->getCharge() ) < 0.5 )
	    continue;

	  if( stable != 1)
	    continue;

	  if(isdecay == true)
	    continue;	    

	  const int mcpid = mcp->id() ;
	  const int mcpdg = mcp->getPDG() ;

	  const double* vertex =mcp->getVertex();
	  double vtx = sqrt(vertex[0]*vertex[0] + vertex[1]*vertex[1] + vertex[2]*vertex[2]);


	  if(vtx < 10.0 && mcpdg==13){
	    
	    const double* mom = mcp->getMomentum() ;
	    double momentum = sqrt( mom[0] * mom[0] +  mom[1] * mom[1] +  mom[2] * mom[2] ) ;
	    
	    double theta = acos(mom[2]/momentum)*180/M_PI;
	    double tran_mom =sqrt( mom[0] * mom[0] +  mom[1] * mom[1]);
	    
	  
	    h.fill(hMCPMomentum, theta) ;
	    // This is the MCParticle I want to check
	    origin_source.push_back(mcp); //fill mc vector containing the true information 

	  }
	}

	// Now I have one MCP, looking for it in MarlinTrkTrack and PFOs.
	////////////////////  TRACK  ///////////////////////
	// set to false for this MCP, until it has been found in MarlinTrkTrack, it will be set to true
	trkFound=false;
	for ( int i=0, m=origin_source.size(); i<m ; i++){
	  trkFound=false;

	  // This MCPatilce information
	  const double* mom =  dynamic_cast<MCParticle*>(origin_source[i])->getMomentum() ;
	  double momentum = sqrt( mom[0] * mom[0] +  mom[1] * mom[1] +  mom[2] * mom[2] ) ;
	    
	  double theta_Mcp = acos(mom[2]/momentum)*180/M_PI;
	  double tran_mom_Mcp =sqrt( mom[0] * mom[0] +  mom[1] * mom[1]);


	  std::vector<Track*> mcp_Trk;

	  LCObjectVec frompars = _navMcpTrk->getRelatedToObjects(origin_source[i]);

	  for( unsigned int j = 0; j < frompars.size(); j++ ){
	    Track* recTrk = dynamic_cast< Track* >( frompars[j] );
	    mcp_Trk.push_back(recTrk);
	  } // check every track has relatrion to this MCParticle

	  //std::cout << "mcp_Trk.size(): "<< mcp_Trk.size() <<std::endl;


	  // If tmp_Trk is empty, this MCParticle has not been reconstructed successfully
	  // If tmp_Trk is not empty, found something related to this MCPartilce, check the details.
	  for ( int k=0, ntrk=mcp_Trk.size(); k<ntrk ; k++){

	    // Found Track information
	    /// momentum and theta calaculation 
	    float omega = mcp_Trk[k]->getOmega();
	    _bField = MarlinUtil::getBzAtOrigin(); //3.5 magnetic field of detector 
	    float Trans_mom = (0.0003)*(_bField/omega);
	    double tran_mom_Trk;
	    tran_mom_Trk = static_cast<double>(Trans_mom);
	    
	    double trk_dif_p = tran_mom_Trk + tran_mom_Mcp;
	    

	    float ftl = mcp_Trk[k]->getTanLambda();
	    double tanlan = static_cast<double>(ftl); 
	    double trk_theta = acos(tanlan/sqrt(1+(tanlan*tanlan)))*(180/M_PI);
	    double trk_dif_a = trk_theta-theta_Mcp;  

	    //differnece graphs 
	    h.fill(hTrkMomentum, trk_theta);
	    h.fill(Momentumdif, trk_dif_a);   
	    
	    h.fill(ptdif, trk_dif_p);

	    if(trk_dif_p<0.15 && trk_dif_a<0.1){
	      trkFound = true;
	    } 

	  } // End of details Track check

	  // Fill trkFound final information
	  //efficency graphs 
	  eff->Fill(trkFound, theta_Mcp);
	  eff2D->Fill(trkFound, theta_Mcp, tran_mom_Mcp);

	} // loop over all MCParticles




	// Now I have one MCP, looking for it in  PFOs.
	////////////////////  PFOs  ///////////////////////
	// set to false for this MCP, until it has been found in MarlinTrkTrack, it will be set to true
	recFoundCharge=false;
	recFoundType=false;
	for ( int i=0, m=origin_source.size(); i<m ; i++){
	  recFoundCharge=false;
	  recFoundType=false;

	  // This MCPatilce information
	  const int mcpdg = dynamic_cast<MCParticle*>(origin_source[i])->getPDG();
	  const double* mom =  dynamic_cast<MCParticle*>(origin_source[i])->getMomentum() ;
	  double momentum = sqrt( mom[0] * mom[0] +  mom[1] * mom[1] +  mom[2] * mom[2] ) ;
	    
	  double theta_Mcp = acos(mom[2]/momentum)*180/M_PI;
	  double tran_mom_Mcp =sqrt( mom[0] * mom[0] +  mom[1] * mom[1]);
	  


	  std::vector<ReconstructedParticle*> mcp_pfo;

	  LCObjectVec frompars = _navMcpRec->getRelatedToObjects(origin_source[i]);

	  for( unsigned int j = 0; j < frompars.size(); j++ ){
	    ReconstructedParticle* recPfo = dynamic_cast< ReconstructedParticle* >( frompars[j] );
	    mcp_pfo.push_back(recPfo);
	  } // check every PFOs has relatrion to this MCParticle
	  
	  // If tmp_pfo is empty, this MCParticle has not been reconstructed successfully
	  // If tmp_pfo is not empty, found something related to this MCPartilce, check the details.
	  for ( int k=0, np=mcp_pfo.size(); k<np ; k++){
	    
	    if( std::fabs( mcp_pfo[k]->getCharge() ) < 0.5 ) //charge check 
	      continue ;

	    /// momentum and theta calaculation 
	    const double* mom_p = mcp_pfo[k]->getMomentum() ;
	    double momentum_p = sqrt( mom_p[0] * mom_p[0] +  mom_p[1] * mom_p[1] +  mom_p[2] * mom_p[2] ) ;
  	    double theta_p = acos(mom_p[2]/momentum_p)*180/M_PI;

	    int typep = mcp_pfo[k]->getType();
	    // float goodness = mcp_pfo[k]->getGoodnessOfPID(); 
	    //std::cout<<goodness<<std::endl;
	    //h.fill(Good, goodness);
	    double tran_mom_p =sqrt( mom_p[0] * mom_p[0] +  mom_p[1] * mom_p[1]);

	    double difference = theta_p-theta_Mcp;
	    double mom_dif = tran_mom_p-tran_mom_Mcp;

	    if(difference<=0.1 && mom_dif<=0.2 ){
	      recFoundCharge=true;
	    }

	    if(difference<=0.1 && mom_dif<=0.2 && typep == mcpdg ){
	      recFoundType=true;
	    }
	    
	    if(typep == mcpdg){
	      h.fill(hRec, theta_p);}
	      //h.fill(Good, goodness);}
	  }

	  //efficency graphs 
	  eff_rec->Fill(recFoundCharge, theta_Mcp) ; 
	  eff2D_pan->Fill(recFoundCharge, theta_Mcp, tran_mom_Mcp);

	  eff_rec_type->Fill(recFoundType, theta_Mcp) ; 
	  eff2D_pan_type->Fill(recFoundType, theta_Mcp, tran_mom_Mcp);

	  
	} // loop over all MCParticles

	/*
	////////////////////  TRACK  ///////////////////////
	// set to false for this MCP, until it has been found in MarlinTrkTrack, it will be set to true
	trkFound=false;
	  for( int j=0, m=colTrk->getNumberOfElements(); j<m ; j++ ){
    
	    Track* mcp_Trk = dynamic_cast<Track*>( colTrk->getElementAt(j) ) ;

	    /// momentum and theta calaculation 
	    float omega = mcp_Trk->getOmega();
	    _bField = MarlinUtil::getBzAtOrigin(); //3.5 magnetic field of detector 
	    float Trans_mom = (0.0003)*(_bField/omega);
	    double trans_mom;
	    trans_mom = static_cast<double>(Trans_mom);
	    
	    double trk_dif_p = trans_mom + tran_mom;
	    

	    float ftl = mcp_Trk->getTanLambda();
	    double tanlan = static_cast<double>(ftl); 
	    double trk_theta = acos(tanlan/sqrt(1+(tanlan*tanlan)))*(180/M_PI);
	    double trk_dif_a = trk_theta-theta;  
	    if(*vertex < 10.0 && mcpdg==13){
	      h.fill(hTrkMomentum, trk_theta);}
	    recon_source.push_back(mcp_Trk);

  
	    //differnece graphs 
	    
	    h.fill(Momentumdif, trk_dif_a);   
	    
	    h.fill(ptdif, trk_dif_p);
	    
     
	    TrackGetSource(recon_source, origin_source, _navpo);

	    //check that truth is true 
	    for(int i = 0; i<origin_source.size(); i++){
	      for(int j =0; j<origin_source[i].size(); j++){
		const int truthPDG = origin_source[i][j]->getPDG();
		if(trk_dif_p<0.15 && trk_dif_a<0.1 && truthPDG==13){
		  trkFound = true;
		} }
	    } 

	  }
	  

	  //efficency graphs 
	  eff->Fill(trkFound, theta);
	  eff2D->Fill(trkFound, theta, tran_mom);

	  ////////////////////////////////////recol///////////////////////////////////////////   
	  for( int k=0, z=colrec->getNumberOfElements(); k<z ; k++ ){
    
	    ReconstructedParticle* mcp_P = dynamic_cast<ReconstructedParticle*>( colrec->getElementAt(k) ) ;

	    /// momentum and theta calaculation 
	    const double* mom_p = mcp_P->getMomentum() ;
	    double momentum_p = sqrt( mom_p[0] * mom_p[0] +  mom_p[1] * mom_p[1] +  mom_p[2] * mom_p[2] ) ;
  	    double theta_p = acos(mom_p[2]/momentum_p)*180/M_PI;

	    int type = mcp_P->getType();

	    if( std::fabs( mcp_P->getCharge() ) < 0.5 ) //charge check 
	      continue ;

	    double tran_mom_p =sqrt( mom_p[0] * mom_p[0] +  mom_p[1] * mom_p[1]);


	    pando_source.push_back(mcp_P);

	    if(*vertex < 10.0 && type == 13){
	      h.fill(hRec, theta_p);  }
	    TrackGetSource2(pando_source, original_source, _navpol);

	    double difference = theta_p-theta;


	    //get truth 
	    double mom_dif = tran_mom_p-tran_mom;
	    for(int i = 0; i<original_source.size(); i++){
	      for(int j =0; j<original_source[i].size(); j++){
		const int truthPDG = original_source[i][j]->getPDG();
		if(difference<=0.1 && mom_dif<=0.2 && type == 13 ){
		  recFounded=true;
		}
	      }
	    } }

	//efficency graphs 
	eff_rec->Fill(recFounded, theta) ; 
	eff2D_pan->Fill(recFounded, theta, tran_mom);

	
	  }}
	*/

	//      }}}
	
}


void RECDiagnostics::check( LCEvent * evt ) { 
  streamlog_out(DEBUG4) << " RECDiagnostics::check event "
			<< evt->getEventNumber() <<std::endl;
}
void RECDiagnostics::end(){ 

  streamlog_out(MESSAGE) << "RECDiagnostics::end()  " << name() 
			 << " processed " << _nEvt << " events in " << _nRun << " runs "
			 << std::endl ;
}

//get information from the LCRelation
void  RECDiagnostics::TrackGetSource(std::vector<Track*> &source, std::vector<std::vector<MCParticle*> >  &to, LCRelationNavigator* &relation)
{
  for( unsigned int i = 0; i < source.size(); i++ ){
    std::vector<MCParticle*> to_tmp;

    LCObjectVec frompars = relation->getRelatedToObjects(source[i]);

    for( unsigned int j = 0; j < frompars.size(); j++ ){
      MCParticle* recpfo = dynamic_cast< MCParticle* >( frompars[j] );
      to_tmp.push_back(recpfo);
    }
    to.push_back(to_tmp);
    to_tmp.clear();
  }

  return;
}


void  RECDiagnostics::TrackGetSource2(std::vector<ReconstructedParticle*> &source, std::vector<std::vector<MCParticle*> >  &to, LCRelationNavigator* &relation)
{
  for( unsigned int i = 0; i < source.size(); i++ ){
    std::vector<MCParticle*> to_tmp;

    LCObjectVec frompars = relation->getRelatedToObjects(source[i]);

    for( unsigned int j = 0; j < frompars.size(); j++ ){
      MCParticle* recpfo = dynamic_cast< MCParticle* >( frompars[j] );
      to_tmp.push_back(recpfo);
    }
    to.push_back(to_tmp);
    to_tmp.clear();
  }

  return;
}
