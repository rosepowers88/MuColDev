#include "AnaProcessors/RecoHistProcessor.hxx"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>

#include <marlin/AIDAProcessor.h>

#include <cmath>
#include <fstream>
#include <iostream>

RecoHistProcessor aRecoHistProcessor ;

RecoHistProcessor::RecoHistProcessor()
  : Processor("RecoHistProcessor") {
  // modify processor description
  _description = "RecoHistProcessor makes histograms at the reco level, analyzing clusters, tracks, and reconstructed particles." ;


  /***************************************************
   * Rose Powers (2023)
   * Use this processor to:
   * - Make root histograms at the reconstruction level
   * - Analyze clusters, PFOs, and tracks
   * - Isolate a certain reco particle (taus here) to analyze
   *
   *
   *
   *
   *
   *
   *
   ***************************************************/
  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("MinPt",
                             "Minimum particle pT.",
                             _minPt,
                             _minPt);
  registerProcessorParameter("MinTheta",
			     "Minimum polar angle",
			     _minTheta,
			     _minTheta);

  registerInputCollection(LCIO::SIMTRACKERHIT,
			  "InputCollectionSimHits",
			  "Name of STH Collection",
			  _inputCollectionSimHits,
			  _inputCollectionSimHits
			  );

  registerInputCollection(LCIO::MCPARTICLE,
			  "InputCollectionNameMCP",
			  "Name of MCP collection",
			  _inputCollectionNameMCP,
			  _inputCollectionNameMCP
			  );

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "InputCollectionNameReco",
			  "Name of RecoParticle collection.",
			  _inputCollectionNameR,
			  _inputCollectionNameR
			  );
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "TauCollection",
			  "Name of tau reco collection",
			  _TauCollection,
			  _TauCollection
			  );
  registerInputCollection(LCIO::CLUSTER,
			  "InputCollectionNameCluster",
			  "Name of cluster collection.",
			  _inputCollectionNameC,
			  _inputCollectionNameC
			  );
  registerInputCollection(LCIO::TRACK,
			  "InputCollectionNameTrack",
			  "Name of track collection.",
			  _inputCollectionNameT,
			  _inputCollectionNameT
			  ); 
  registerProcessorParameter("PDG",
			     "PDG for PFOs",
			     _pfopdg,
			     _pfopdg);

}


void RecoHistProcessor::book_pfo_histograms(){
  marlin::AIDAProcessor::histogramFactory(this);  
  _h_ptpf   = new TH1F("PFO_pT", ";pT [GeV]", 1000,0., 300);
  _h_ppf    = new TH1F("PFO_p", ";p [GeV]", 1000, 0, 300);
  _h_pdgpf  = new TH1F("PFO_PDG_ID", ";PDG ID", 10000, -500, 2500);
  _h_Npf    = new TH1F("PFO_nObjectsFound", ";N", 50, 0, 10);
  _h_Epf    = new TH1F("PFO_Energy", ";E [GeV]", 1000,0,300);
  _h_phipf  = new TH1F("PFO_phi", ";#phi [rad]", 100,-3, 3);
  _h_thetapf= new TH1F("PFO_theta", ";#theta [rad]", 100,0,3.14);
  _h_nClusters = new TH1F("PFO_ncl", ";NCl", 50,0,10);
 _h_inv_E = new TH1F("Invariant Energy - Set Energy", ";#Delta E [GeV]", 200,-500.,500.);
  _h_inv_msqr = new TH1F("Invariant Mass Squared", ";m^2 [GeV^2]", 200,-100000,50.);
  _h_mass = new TH1F("Set Mass", ";m [GeV]", 100,0,10);
  _h_inv_p = new TH1F("Invariant Momentum - Set Momentum", ";#Delta p [GeV]", 200,-500,100);
  _h_cl_res = new TH1F("Cluster Energy - PFO Energy", ";#Delta E [GeV]", 100,-500,500);
  _h_PDG_v_E = new TH2F("PDG against E", "PDG;E [GeV]", 2500,-250,2200, 250,0,300);


}
void RecoHistProcessor::book_cluster_histograms(){
  marlin::AIDAProcessor::histogramFactory(this);
  _h_Ecl = new TH1F("Cluster Energy", ";E [GeV]", 1000,0,300);
  _h_Xcl = new TH1F("Cluster X-position", ";x [mm]", 4000,-2000,2000);
  _h_Ycl = new TH1F("Cluster Y-position", ";y [mm]", 4000,-2000,2000);
  _h_Zcl = new TH1F("Cluster Z-position", ";z [mm]", 8000,-4000,4000);
  _h_iPhicl=new TH1F("Cluster Intrinsic Phi", ";#phi [rad]", 10,-3.2,3.2);
  _h_iThetacl=new TH1F("Cluster Intrinsic Theta", ";#theta [rad]",100,-6.4,6.4);
  _h_PIDcl = new TH1F("Cluster Most Likely PID", "; PID", 10000, -300, 2500);
  _h_Ncl    = new TH1F("nClustersFound", ";N", 50, 0, 10);
  _h_CLpdg = new TH1F("PDG of clust", ";PDG ID", 5000,-1000,2500);
  _h_Eres = new TH1F("Energy Resolution", "; (E_{truth}-E_{clust})/E_{truth}", 1000,-50,50);
  _h_nhits = new TH1F("NHits per cluster", ";NHits", 100,0,100);
  _h_nHcalHits = new TH1F("NHCal hits per event", ";NHCalHits", 500,0,1000);
  _h_nEcalHits = new TH1F("NECal hits per event", ";NECalHits", 500,0,1000);
  _h_E_v_theta = new TH2F("E against Theta", "E [GeV];#theta [rad]", 50,0,300, 50,-6.4,6.4);

}
void RecoHistProcessor::book_track_histograms(){
  marlin::AIDAProcessor::histogramFactory(this);
  _h_D0trk = new TH1F("Track Impact Parameter (r-phi)", ";D0 [mm]", 100,0,200);
  _h_phitrk = new TH1F("Track Phi", ";#phi [rad]", 100,-3.2,3.2);
  _h_omegatrk = new TH1F("Track Signed Curvature", ";#Omega [1/mm]", 1000,-200,200);
  _h_Z0trk = new TH1F("Track Impact Parameter (r-z)", ";Z0 [mm]", 100,0,200);
  _h_lamtrk =new TH1F("Track Dip Angle (r-z)", ";#lambda [rad]", 10,-3.2,3.2);
  _h_dEdxtrk=new TH1F("Track dEdx", ";dEdx [Gev/mm]",100,0,50);
  _h_Ntrk    = new TH1F("nTracksFound", ";N", 50, 0, 10);
  _h_trkpdg = new TH1F("pdg of Track", ";PDG ID", 10000, -500, 2500);
  _h_ntrkhits = new TH1F("number of track hits", ";nhits", 100,0,50);
  _h_deltaR= new TH1F("deltaR", ";#Delta r [mm]", 100,0,5000);

}

void RecoHistProcessor::book_trkzero_histograms(){ //to analyze events with no tracks
  marlin::AIDAProcessor::histogramFactory(this);
  _h_p0 = new TH1F("Momentum, 0trk",";p [GeV]", 100,0,200);
  _h_pt0 = new TH1F("pT, 0trk", ";pT [GeV]", 100, 0, 200);
  _h_E0 = new TH1F("Energy, 0trk", ";E [GeV]", 100, 0, 200);
  _h_pdg0 = new TH1F("PDG ID, 0trk", ";PDG ID", 10000, -250, 2500);
  _h_NMCP = new TH1F("Number of MCP, 0trk", ";NParticles", 100, 0, 500);
}

void RecoHistProcessor::book_tau_histograms(){
   marlin::AIDAProcessor::histogramFactory(this);  
  _h_ptTau   = new TH1F("Tau_pT", ";pT [GeV]", 1000,0., 200);
  _h_pTau    = new TH1F("Tau_p", ";p [GeV]", 1000, 0, 200);
  _h_pdgTau  = new TH1F("Tau_PDG_ID", ";PDG ID", 10000, -500, 2500);
  _h_NTau    = new TH1F("Tau_nObjectsFound", ";N", 50, 0, 10);
  _h_ETau    = new TH1F("Tau_Energy", ";E [GeV]", 1000,0,200);
}

void RecoHistProcessor::init() {
  // Print the initial parameters
  printParameters() ;

  //Basic histograms
  book_pfo_histograms();
  book_cluster_histograms();
  book_track_histograms();
  book_trkzero_histograms();
  book_tau_histograms();
}

bool _has_valid_pt=false;

void RecoHistProcessor::processRunHeader( LCRunHeader* /*run*/) {
}  

//Get a reference vector of sim hits and their PDGS for PDG tagging of tracks 
std::vector<std::vector<int>> RecoHistProcessor::getSimTrackerHitIDs(LCCollection* trackhitCol){
  std::vector<std::vector<int>> simhitIDs;
  const int nHits = trackhitCol->getNumberOfElements();
  for(int i=0; i<nHits; i++){
    std::vector<int> ID_PDG;
    const EVENT::SimTrackerHit *sth = static_cast<const EVENT::SimTrackerHit*>(trackhitCol->getElementAt(i));
    int ID=sth->getCellID0();
    const EVENT::MCParticle *mcp = static_cast<const EVENT::MCParticle*>(sth->getMCParticle());
    int pdg=mcp->getPDG();
    ID_PDG.push_back(ID);
    ID_PDG.push_back(pdg);
    simhitIDs.push_back(ID_PDG);
    
  }
  return(simhitIDs);
}
void RecoHistProcessor::fill_pfo_histograms(LCCollection* inputCol){

  //get number of elements
  const int nEl = inputCol->getNumberOfElements();
  _h_Npf->Fill(nEl);
  _has_valid_pt=false;
  int nPF = 0;
  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::ReconstructedParticle *pfo=static_cast<const EVENT::ReconstructedParticle*>(inputCol->getElementAt(i));
    // get pdg
    int pdg = pfo->getType();
    nPF++;
    _h_pdgpf->Fill(pdg);
    //get momentum and energy
    const double* mom = pfo->getMomentum();

    double E=pfo->getEnergy();
    double pt=std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    double p =std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2)+std::pow(mom[2],2));
    double phi = std::atan(mom[1]/mom[0]);
    double theta = std::atan(mom[2]/p);
    double mass = pfo->getMass();
    _h_mass->Fill(mass);
 double inv_E = std::sqrt(std::pow(p,2)+std::pow(mass,2));
    _h_inv_E->Fill(inv_E-E);
    _h_PDG_v_E->Fill(pdg, E);

    double inv_p = std::sqrt(std::pow(E,2)-std::pow(mass,2));
    _h_inv_p->Fill(inv_p-p);
    
    double invm = (std::pow(E,2)-std::pow(p,2));
    _h_inv_msqr->Fill(invm);
    double Eclust = 0;
    const EVENT::ClusterVec clusts = pfo->getClusters();
    if(clusts.size()>0){
      for(int cl_iter = 0; cl_iter < clusts.size(); cl_iter++){
	const EVENT::Cluster* cl = static_cast<const EVENT::Cluster*>(clusts[cl_iter]);
	Eclust+=cl->getEnergy();
      }
      _h_cl_res->Fill(Eclust-E);
    }
    
    

    //momentum check
    if(pt < _minPt){
      continue;
    }

    //theta check
    if(theta < _minTheta){
      continue;
    }
    EVENT::ClusterVec clvec = pfo->getClusters();
    const int cl_in_pf = clvec.size();
    _h_nClusters->Fill(cl_in_pf);
    _has_valid_pt=true;
    _h_ppf->Fill(p);
    _h_ptpf->Fill(pt);
    _h_Epf->Fill(E);
    _h_phipf->Fill(phi);
    _h_thetapf->Fill(theta);


  }
  if(_has_valid_pt==true){
    //_h_Npf->Fill(nPF);
  }

      
}


/*THIS IS THE FUNCTION THAT FILLS THE TAU HISTOGRAMS*/

void RecoHistProcessor::fill_tau_histograms(LCCollection* inputCol){
  //get number of elements
  const int nEl = inputCol->getNumberOfElements();
  // _h_Npf->Fill(nEl);
  _has_valid_pt=false;

  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::ReconstructedParticle *tau=static_cast<const EVENT::ReconstructedParticle*>(inputCol->getElementAt(i));
    //get momentum and energy
    const double* mom = tau->getMomentum();
    double E=tau->getEnergy();
    double pt=std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    double p =std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2))+std::sqrt(std::pow(mom[2],2));
    //momentum check
    if(pt < _minPt){
      continue;
    }

    //theta check
    double polar = std::acos(mom[2]/p);
    if(polar < _minTheta){
      continue;
    }

    _has_valid_pt=true;
    _h_pTau->Fill(p);
    _h_ptTau->Fill(pt);
    _h_ETau->Fill(E);

    //finally, get pdg ... should all be +/- 15, but this is on as a check
    int pdg = tau->getType();
    _h_pdgTau->Fill(pdg);
  }
  if(_has_valid_pt==true){
    _h_NTau->Fill(nEl);
  }

      
}
  
void RecoHistProcessor::fill_cluster_histograms(LCCollection* inputCol, LCEvent *evt){

  //get number of elements
  const int nEl = inputCol->getNumberOfElements();
  _h_Ncl->Fill(nEl);
 
 if(nEl!=0){
    LCCollection *ecb;
    LCCollection *ece;
    LCCollection *hcb;
    LCCollection *hce;
    int nhcalhits = 0;
    int necalhits = 0;
    bool hasECB = true;
    bool hasECE = true;
    bool hasHCB = true;
    bool hasHCE = true;
    try{
      ecb = evt->getCollection("EcalBarrelCollectionRec");
    }
    catch(Exception& e){
      hasECB = false;
    }
    try{
      ece = evt->getCollection("EcalEndcapCollectionRec");
    }
    catch(Exception& e){
      hasECE = false;
    }
    try{
      hcb = evt->getCollection("HcalBarrelsCollectionRec");
    }
    catch(Exception& e){
      hasHCB = false;
    }
    try{
      hce = evt->getCollection("HcalEndcapsCollectionRec");
    }
    catch(Exception& e){
      hasHCE = false;
    }
    if(hasECB) necalhits += ecb->getNumberOfElements();
    if(hasECE) necalhits += ece->getNumberOfElements();
    if(hasHCB){
      nhcalhits += hcb->getNumberOfElements();
    }
    if(hasHCE) nhcalhits += hce->getNumberOfElements();
    //if(!nhcalhits) std::cout<<necalhits<<std::endl;

    _h_nEcalHits->Fill(necalhits);
    _h_nHcalHits->Fill(nhcalhits);
  }

  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::Cluster *cl=static_cast<const EVENT::Cluster*>(inputCol->getElementAt(i));
    //get nhits
    const EVENT::CalorimeterHitVec hits = cl->getCalorimeterHits();
    _h_nhits->Fill(hits.size());
    //get energy and position
    double Ecl = cl->getEnergy();
    _h_Ecl->Fill(Ecl);
    const float* pos = cl->getPosition();
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double rad = std::sqrt(x*x + y*y + z*z);
    _h_Xcl->Fill(rad);
    _h_Ycl->Fill(y);
    _h_Zcl->Fill(z);

    //Get phi and theta
    double phi = cl->getIPhi();
    _h_iPhicl->Fill(phi);
    double theta = cl->getITheta();
    _h_iThetacl->Fill(theta);

    _h_E_v_theta->Fill(Ecl,theta);

    //Get most likely particle ID
    const ParticleIDVec& PID = cl->getParticleIDs();
    if(PID.size()!=0){
      double pid_pdg = PID[0]->getPDG();
      std::cout<<pid_pdg<<std::endl;
      _h_PIDcl->Fill(pid_pdg);
    }


  }
}
void RecoHistProcessor::zeroTrackAna(LCCollection* inputCol){
  //call if the event has zero tracks, take a look at the properties at the MC level to see where we are failing
  const int nEl = inputCol->getNumberOfElements();
  _h_NMCP->Fill(nEl);
  
  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::MCParticle *mcp=static_cast<const EVENT::MCParticle*>(inputCol->getElementAt(i));
    const double* mom = mcp->getMomentum();

    double E=mcp->getEnergy();
    double pt=std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    double p =std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2))+std::sqrt(std::pow(mom[2],2));

    _h_p0->Fill(p);
    _h_pt0->Fill(pt);
    _h_E0->Fill(E);

    //finally, get pdg
    int pdg = mcp->getPDG();
    _h_pdg0->Fill(pdg);
  }
}
void RecoHistProcessor::fill_track_histograms(LCCollection* inputCol, LCCollection* MCCol, LCCollection* trkCol){

  //get number of elements
  const int nEl = inputCol->getNumberOfElements();
  _h_Ntrk->Fill(nEl);
  /*if(nEl==0){
    zeroTrackAna(MCCol);
    }*/

  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::Track *trk=static_cast<const EVENT::Track*>(inputCol->getElementAt(i));
    const EVENT::TrackerHitVec trkhits = trk->getTrackerHits();
    int ntrkhits = trkhits.size();


    _h_ntrkhits->Fill(ntrkhits);
    
    //get D0
    double D0 = trk->getD0();
    _h_D0trk->Fill(D0);

    //get phi
    double phi = trk->getPhi();
    _h_phitrk->Fill(phi);

    //get curvature
    double omega = trk->getOmega();
    _h_omegatrk->Fill(omega);

    //get Z0
    double Z0 = trk->getZ0();
    _h_Z0trk->Fill(Z0);

    //get dip angle
    double tanlam = trk->getTanLambda();
    double lambda = std::atan(tanlam);
    _h_lamtrk->Fill(lambda);

    //get dE/dx
    double dEdx = trk->getdEdx();
    _h_dEdxtrk->Fill(dEdx);

    //get delta R
    double rmin=1000000;
    double rmax = -1;
    double dRmax = -1;
    double r_prev=0;
    for(uint32_t n=0; n<trkhits.size(); n++){
      const double * pos = trkhits[n]->getPosition();
      double R = std::sqrt(std::pow(pos[0],2)+std::pow(pos[1],2)+std::pow(pos[2],2));
      double dR=R-r_prev;
      //if(R < rmin) rmin = R;
      if(dR > dRmax) dRmax = dR;
      r_prev=R;
    }
    //DeltaR=rmax-rmin;
    _h_deltaR->Fill(dRmax);
  }
  
}



void RecoHistProcessor::processEvent( LCEvent * evt ) {
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.
  LCCollection* MCCol = evt->getCollection(_inputCollectionNameMCP);
  //LCCollection* TrkCol = evt->getCollection(_inputCollectionSimHits);
  LCCollection* tauCol = evt->getCollection(_TauCollection);

  /* LCCollection *SiTracks = evt->getCollection("MySiTracks");
  if(_has_valid_pt==true){
    fill_track_histograms(SiTracks);
    }*/

  if(_has_valid_pt==true){
    fill_tau_histograms(tauCol);
  }
  std::vector<std::string> inputnames {_inputCollectionNameR, _inputCollectionNameC, _inputCollectionNameT};
  for(int i=0; i<inputnames.size(); i++){
    LCCollection* inputCol;
    try{
      inputCol = evt->getCollection(inputnames[i]);
    }
    catch(Exception&e){
      continue;
    }

  
    if( inputCol->getTypeName() == lcio::LCIO::RECONSTRUCTEDPARTICLE ) {
      //fill the histograms for PFOs
      //fill_pfo_histograms(inputCol);
    }
    else if( inputCol->getTypeName() == lcio::LCIO::CLUSTER) {
      //fill the histograms for clusters
      //fill_cluster_histograms(inputCol,evt);
    }
    else if( inputCol->getTypeName() == lcio::LCIO::TRACK ) {
      //fill the histograms for tracks
      // fill_track_histograms(inputCol,MCCol,nullptr);
    }
    else{
      //throw EVENT::Exception( "Invalid collection type: " + inputCol->getTypeName() ) ;
      continue;
    }
  }
 }
  

//function that just finds the mode of a vector or tells you if it is multimodal
int RecoHistProcessor::findMode(std::vector<int> vec){
  std::vector<int> elements;
  std::vector<int> modes;
  for(int i=0; i<vec.size(); i++){
    if( std::find(elements.begin(), elements.end(), vec[i])== elements.end()){
      elements.push_back(vec[i]);
    }
    else continue;
  }
  for(int i=0; i<elements.size(); i++){
    int modeCounter=0;
    for(int j=0; j<vec.size(); j++){
      if(vec[j]==elements[i]){
	modeCounter++;
      }
    }
    modes.push_back(modeCounter);
  }
  int mode = elements[0];
  bool split = false;
  for(int t=1; t<modes.size(); t++){
    if(modes[t]>mode){
      mode=modes[t];
    }
    else if(modes[t]==mode){
      split=true;
    }
    else continue;
  }
  if(split) return(0);
  else return(mode);

  
}

void RecoHistProcessor::check( LCEvent * /*evt*/ )
{ }

void RecoHistProcessor::end()
{ }
