#include "AnaProcessors/CalAna.hxx"

#include <fstream>
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/TrackerHit.h>

#include <marlin/AIDAProcessor.h>
#include <marlin/Statusmonitor.h>

#include <cmath>

CalAna aCalAna;

CalAna::CalAna()
  : Processor("CalAna") {
  _description = "CalAna deals with the Calo cluster  analysis, unsurprisingly" ; 
  /********************************
   *
   * Rose Powers (2023)
   * Use this processor for doing specific studies on a calohit collection
   *
   ****************************************/

  // register steering params

  registerProcessorParameter("MinPt",
			     "Minimum transverse momentum",
			     _minPt,
			     _minPt);
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "PFOInput",
			  "Name of RecoParticle collection",
			  _inputcolPFO,
			  _inputcolPFO
			  );
  registerInputCollection(LCIO::MCPARTICLE,
			  "MCInput",
			  "Name of MCP collection",
			  _inputcolMC,
			  _inputcolMC
			  );
}

void CalAna::init(){
  printParameters();

  marlin::AIDAProcessor::histogramFactory(this);
  // _h_pion_misID = new TH1F("pion_ID", ";PDG ID", 10000, -500, 2500);
  _h_nocluster = new TH1F("ID of clusterless MCP", ";PDG ID", 10000, -500, 2500);


  _h_pass_CL1 = new TH1F("passCL1", "", 30,0,250);
  _h_all_CL1 = new TH1F("allCL1", "", 30,0,250);
  _h_pass_CL3 = new TH1F("passCL3", "", 30,0,250);
  _h_all_CL3 = new TH1F("allCL3", "", 30,0,250);

  _h_nCL = new TH1F("Cluster Occupancy", ";nCL", 100,0,50);
  _h_phi_CL=new TH1F("Cluster Azimuth", ";#phi", 50,-3.142,3.142);
  _h_theta_CL=new TH1F("Cluster Theta", ";#theta", 50,0,3.142);

  _h_pass_pdg = new TH1F("passpdg", "", 10000,-500,2500);
  _h_all_pdg = new TH1F("allpdg", "", 10000,-500,2500);

  _h_pass_phi = new TH1F("passphi", "", 50, -3.142, 3.142);
  _h_all_phi = new TH1F("allphi", "", 50, -3.142, 3.142);

  _h_pass_theta = new TH1F("passtheta", "", 50, 0, 3.142);
  _h_all_theta = new TH1F("alltheta", "", 50,0,3.142);

  _h_pass_E = new TH1F("passE", "",50, 5, 300);
  _h_all_E = new TH1F("allE", "", 50, 5, 300);

  _h_delta_E=new TH1F("#Delta E", ";E_{PFO}-E_{CL}", 100,-50,50);
  

  _h_nUnMatchedCL = new TH1F("Unmatched CL", ";nUMCL", 50,0,10);

  _h_pass_trks3=new TH1F("passtrks3p", "", 100,0,300);
  _h_all_trks3=new TH1F("alltrks3p", "", 100,0,300);

  _h_pass_trks1 = new TH1F("passtrks1p", "", 100,0,300);
  _h_all_trks1 = new TH1F("alltrks1p", "", 100,0,300);

  //_h_deltaBH = new TH1F("nsimBH - nBH", ";#Delta BH", 250, 0, 1000);
  //_h_deltaEH = new TH1F("nsimEH - nEH", ";#Delta EH", 250, 0, 1000);

  // _h_badBarrelContributers = new TH1F("contribution ID", ";PDG_cont", 10000, -500, 2500);

  //_h_noClpT = new TH1F("pT of clusterless PFO", ";pT", 100, 0, 500);

}

/*std::vector<std::vector<int>> CalAna::getSimTrackerHitIDs(LCCollection* trackhitCol){
  //makes a reference of simhits against which we can check track hits to get trackTruth info
  std::vector<std::vector<int>> simhitIDs;
  const int nHits = trackhitCol->getNumberOfElements();
  for(int i=0; i<nHits; i++){
    //get PDG, cell ID, and pT of MCP (useful for efficiency)
    std::vector<int> ID_PDG;
    const EVENT::SimTrackerHit *sth = static_cast<const EVENT::SimTrackerHit*>(trackhitCol->getElementAt(i));
    int ID=sth->getCellID0();
    const EVENT::MCParticle *mcp = static_cast<const EVENT::MCParticle*>(sth->getMCParticle());
    int pdg=mcp->getPDG();
    const double *mom = mcp->getMomentum();
    double pT = std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    double p = std::sqrt(std::pow(pT,2)+std::pow(mom[2],2));
    const EVENT::MCParticleVec Mvec = mcp->getParents();
    int MID=0;
    if(Mvec.size()!=0){
      MID = Mvec.back()->getPDG();
    }
    else{
      MID = -1;
    }
    ID_PDG.push_back(ID);
    ID_PDG.push_back(pdg);
    ID_PDG.push_back(pT);
    ID_PDG.push_back(MID);
    simhitIDs.push_back(ID_PDG);
    
  }
  return(simhitIDs);
}
double CalAna::anaPiEff( const LCObject * inputTrk, LCCollection* trkHitCol){
  //find the efficiency for charged pion FROM TAU reconstruction
  //trk part
  //first get the simhits and make a reference
  //find a way to associate to specific MCP tho
  
  std::vector<std::vector<int>> simHitRef = getSimTrackerHitIDs(trkHitCol);

  int npi = 0; //how many pions
  bool isPion = false; //the bool determining whether a track came from a tau-generated charged pi
  std::vector<std::vector<int>> usedHits; //for de-doubling fake tracks

// get the track 
  const EVENT::Track *trk=static_cast<const EVENT::Track*>(inputTrk);
  double trkpT = 0;

// get the hit collection
  const EVENT::TrackerHitVec trkhits = trk->getTrackerHits();

// get number of tracker hits in the track
  const int ntrkhits = trkhits.size();
// count how many of them came from a charged pion
  int pihits = 0;
  double ptavg = 0;
  double hitpT=0;

  for(uint32_t t=0; t<ntrkhits; t++){
    const EVENT::TrackerHit* th =static_cast<const EVENT::TrackerHit*>( trkhits[t]);
    int ID=th->getCellID0();
    int hitPDG = 0;
    int hitMID = 0;
    //now see if it matches any cell ID's in the reference, record PDG and pT
    for(uint32_t l=0; l<simHitRef.size(); l++){
      if(simHitRef[l][0]==ID){
	hitPDG = simHitRef[l][1];
	hitpT = simHitRef[l][2];
	hitMID = simHitRef[l][3];
      }
    }
    //check if it came from pion from tau
    if(abs(hitPDG)==211 && abs(hitMID)==15){
      pihits++;
    }
    ptavg+=hitpT;
  }
  trkpT = double(ptavg)/ntrkhits;

  //now we want to see what fraction of the hits that make up the track are pion-like
  //if over 50 percent, we have a pion
  double piFrac = double(pihits)/ntrkhits;
  if(piFrac >=0.5){
    isPion = true;
  }
  if(isPion==true){
    return(trkpT);
  }
  else{
    return(-1000);
  }   

*/
void CalAna::processRunHeader( LCRunHeader * /*run*/){
}

void CalAna::processEvent( LCEvent * evt) {
  LCCollection* MCPs = evt->getCollection(_inputcolMC);
  const int nMCP = MCPs->getNumberOfElements();
  LCCollection* taus;
  LCCollection* PFOs = evt->getCollection(_inputcolPFO);
  LCCollection* allCL = evt->getCollection("PandoraClusters");
  const int ncl = allCL->getNumberOfElements();
  LCCollection* trks = evt->getCollection("SiTracks_Refitted");
  const int ntrk = trks->getNumberOfElements();



  //isolate events that have at least one track or one cluster
  if(ntrk > 0 || ncl > 0){
    if(!PFOs->getNumberOfElements()){
      int en=evt->getEventNumber();
      std::cout<<ntrk<<", "<<ncl<<std::endl;
      std::cout<<en<<std::endl;
    }


  }

  //GOAL 9/27/23: compare energy of neutral CLUSTERS against gen-level PI-0 energy
  
  //cluster analysis
  const int nCL = allCL->getNumberOfElements();
  _h_nCL->Fill(nCL);
  for(uint32_t i=0; i<nCL; i++){
    const EVENT::Cluster *cl = static_cast<const EVENT::Cluster*>(allCL->getElementAt(i));
    float phi = cl->getIPhi();
    float theta = cl->getITheta();
    _h_phi_CL->Fill(phi);
    _h_theta_CL->Fill(theta);
  }

  //see how many unmatched CL there are
  //check eff for all pfo
  const int nPF = PFOs->getNumberOfElements();
  int nMatchedCl = 0;
  for(uint32_t i=0; i<nPF; i++){
    const EVENT::ReconstructedParticle *pf = static_cast<const EVENT::ReconstructedParticle*>(PFOs->getElementAt(i));
    EVENT::ClusterVec clvec = pf->getClusters();
    const int cl_in_pf = clvec.size();

    const int pdg = pf->getType();
    if(abs(pdg)!=211) continue;
    //_h_all_pdg->Fill(pdg);
    nMatchedCl+=cl_in_pf;
    
    const double* mom = pf->getMomentum();
    double pt=std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    double p =std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2))+std::sqrt(std::pow(mom[2],2));
    
    //find energy
    double E=pf->getEnergy();
    _h_all_E->Fill(E);

    //find and fill azimuth and polar angle

    double phi = std::atan(mom[1]/mom[0]);
    _h_all_phi->Fill(phi);

    double theta = std::acos(mom[2]/p);
    _h_all_theta->Fill(theta);
    if(cl_in_pf!=0){
      _h_pass_E->Fill(E);
      _h_pass_phi->Fill(phi);
      _h_pass_theta->Fill(theta);
      //_h_pass_pdg->Fill(pdg);
    }
      
  }
  _h_nUnMatchedCL->Fill(nCL-nMatchedCl);


 

  
  //before any PFO analysis, check the tau decay mode
  //if tau decay is not what we're looking for, skip the whole event
  bool isOneProng = false;
  bool isThreeProng = false;
  bool hasNeutrals = false;
  int ncp=0;
  int nnp=0;
  double npE=0;
  for(uint32_t i=0; i<nMCP; i++){
    const EVENT::MCParticle *mcp = static_cast<const EVENT::MCParticle*>(MCPs->getElementAt(i));
    const int mcpdg = mcp->getPDG();
    if(mcpdg != 15) continue; //keep going until we find tau
    //now get the daughter particles, and check the pdgs. 
    const EVENT::MCParticleVec dvec = mcp->getDaughters();
    ncp = 0;
    nnp = 0;
    for(uint32_t d=0; d<dvec.size(); d++){
      int did = dvec[d]->getPDG();
      double E = dvec[d]->getEnergy();
      if(abs(did)==211){
	ncp++;
      }
      if(did==111){
	 nnp++;
	 npE=E;
	 break;
      }
    }
  }

  if(ncp==1) isOneProng = true;
  else if(ncp==3) isThreeProng = true;
  if(nnp>0) hasNeutrals = true;

  int nneutral=0;
  if(hasNeutrals){
    //std::cout<<nnp<<std::endl;
  if(isOneProng || isThreeProng){
  int nPFO = PFOs->getNumberOfElements();
  //loop over the PFO collection and see how many of them have an associated cluster
  std::vector<std::vector<double>> PF_info;
  bool hasCl = false;
  double clE=0;
  for(uint32_t i=0; i<nPFO; i++){
    std::vector<double> infovec;
    const EVENT::ReconstructedParticle *pfo = static_cast<const EVENT::ReconstructedParticle*>(PFOs->getElementAt(i));
    int q=pfo->getCharge();
    EVENT::ClusterVec clvec = pfo->getClusters();
    if(q==0){
      for(uint32_t n=0; n<clvec.size(); n++){
	const EVENT::Cluster *cl = static_cast<const EVENT::Cluster*>(clvec[n]);
	clE+=cl->getEnergy();
    }
       _h_delta_E->Fill(npE-clE);
    }
   
    double pfE = pfo->getEnergy();
    EVENT::TrackVec trks = pfo->getTracks();
    const double* p = pfo->getMomentum();
    double pt = std::sqrt(std::pow(p[0],2)+std::pow(p[1],2));
    if(trks.size()!=0){

      if(isThreeProng){
      _h_pass_trks3->Fill(pt);
      }
      if(isOneProng){
	_h_pass_trks1->Fill(pt);
      }
    }
    if(isThreeProng)_h_all_trks3->Fill(pt);
    if(isOneProng) _h_all_trks1->Fill(pt);
    int PF = pfo->getType();
    if(PF==111) std::cout<< "Neutral PFO"<<std::endl;
    if(PF==2112 || PF == 22){
      nneutral++;
    }
    if(!clvec.size()) {
      if(abs(PF)!=211){
	std::cout<<PF<<std::endl;
      }
      //_h_noClpT->Fill(pt);
      _h_nocluster->Fill(PF);
      hasCl = false;
    }
    else hasCl = true;
    const double* mom = pfo->getMomentum();
    double pT = std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    if(hasCl){
      if(isOneProng) _h_pass_CL1->Fill(pT);
      if(isThreeProng) _h_pass_CL3->Fill(pT);
    }
    if(isOneProng) _h_all_CL1->Fill(pT);
    if(isThreeProng) _h_all_CL3->Fill(pT);
    const int PDG = pfo->getType();

  }

  //std::cout<<nnp<<" neutral pion MC, "<<nneutral<<" neutral PFO"<<std::endl;
  }
  }


  
}
void CalAna::check( LCEvent * /*evt*/ )
{}
void CalAna::end()
{}
