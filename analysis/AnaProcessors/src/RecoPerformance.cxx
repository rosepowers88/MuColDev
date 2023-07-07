#include "AnaProcessors/RecoPerformance.hxx"

#include <fstream>
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>

#include <marlin/AIDAProcessor.h>
#include <marlin/Statusmonitor.h>

#include <cmath>

RecoPerformance aRecoPerformance ;

RecoPerformance::RecoPerformance()
  : Processor("RecoPerformance") {
  // modify processor description
  _description = "RecoPerformance analyzes the efficiency of tracking, clustering, and PFO finding." ;

  // register steering parameters: name, description, class-variable, default value


  registerProcessorParameter("PDG",
			     "PDG code of particle type",
			     _PDG,
			     _PDG);
  registerProcessorParameter("MinPt",
			     "Minimum transverse momentum",
			     _minPt,
			     _minPt);
  registerProcessorParameter("MinTheta",
			     "Minimum polar angle",
			     _minTheta,
			     _minTheta);


  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
			  "InputCollectionNameReco",
			  "Name of RecoParticle collection.",
			  _inputCollectionNameR,
			  _inputCollectionNameR
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
  registerInputCollection(LCIO::MCPARTICLE,
			  "InputCollectionNameMCP",
			  "Name of MCP Collection.",
			  _inputCollectionNameMCP,
			  _inputCollectionNameMCP
			  );
  registerInputCollection(LCIO::SIMTRACKERHIT,
			  "InputCollectionSimHits",
			  "Name of STH Collection",
			  _inputCollectionSimHits,
			  _inputCollectionSimHits
			  );
  
}

void RecoPerformance::init() {
  // Print the initial parameters
  printParameters() ;

  //Basic histograms
  marlin::AIDAProcessor::histogramFactory(this);
  _h_effPFO  = new TH1F("PFO Efficiency", ";nPFO/nMCP", 100,0,1);
  _h_effCl = new TH1F("Cluster Efficiency", ";nCl/nMCP", 100,0,1); //keep on for a check to make sure mother id vetting works
  _h_effTrk = new TH1F("Track Efficiency", ";nTrk/nMCP", 100,0,1);

  //pi histograms
  _h_effPFO_cp = new TH1F("PFO Efficiency for Charged Pions",";nPFO/nCP", 100, 0, 1);
  _h_effPFO_np = new TH1F("PFO Efficiency for Neutral Pions", ";nPFO/nNP", 100, 0, 1);

  //other particle histograms
  _h_effPFO_n = new TH1F("PFO Efficiency for Neutrons", ";nPFO/nN", 100,0,1);
  _h_effPFO_gam = new TH1F("PFO Efficiency for Photons", ";nPFO/nGam", 100, 0 ,1);

}

std::vector<std::vector<int>> RecoPerformance::getSimTrackerHitIDs(LCCollection* trackhitCol){
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
    double phi = std::atan(mom[1]/mom[0]);
    double theta = std::acos(mom[2]/p);
    const EVENT::MCParticleVec Mvec = mcp->getParents();
    int MID = Mvec.back()->getPDG();
    ID_PDG.push_back(ID);
    ID_PDG.push_back(pdg);
    //ID_PDG.push_back(pT);
    ID_PDG.push_back(phi);
    ID_PDG.push_back(MID);
    simhitIDs.push_back(ID_PDG);
    
  }
  return(simhitIDs);
}
double RecoPerformance::anaPiEff( const LCObject * inputTrk, LCCollection* trkHitCol){
  //find the efficiency for charged pion FROM TAU reconstruction
  //trk part
  //first get the simhits and make a reference structure
  std::vector<std::vector<int>> simHitRef = getSimTrackerHitIDs(trkHitCol);
  
  int npi = 0; //how many pions
  bool isPion = false; //the bool determining whether a track came from a tau-generated charged pi
  std::vector<std::vector<int>> usedHits; //for de-doubling fake tracks

// get the track 
  const EVENT::Track *trk=static_cast<const EVENT::Track*>(inputTrk);
  double trkpT = 0;
  double trkphi=0;

// get the hit collection
  const EVENT::TrackerHitVec trkhits = trk->getTrackerHits();

// get number of tracker hits in the track
  const int ntrkhits = trkhits.size();
// count how many of them came from a charged pion
  int pihits = 0;
  //double ptavg = 0;
  //double hitpT=0;
  double phiavg=0;
  double hitphi=0;
  for(uint32_t t=0; t<ntrkhits; t++){
    const EVENT::TrackerHit* th =static_cast<const EVENT::TrackerHit*>( trkhits[t]);
    int ID=th->getCellID0();
    int hitPDG = 0;
    int hitMID = 0;
    //now see if it matches any cell ID's in the reference, record PDG and pT
    for(uint32_t l=0; l<simHitRef.size(); l++){
      if(simHitRef[l][0]==ID){
	hitPDG = simHitRef[l][1];
	hitphi = simHitRef[l][2];
	hitMID = simHitRef[l][3];
      }
    }
    //check if it came from pion from tau
    if(abs(hitPDG)==211 && abs(hitMID)==15){
      pihits++;
    }
    phiavg+=hitphi;
  }
  trkphi = phiavg/double(ntrkhits);
  //now we want to see what fraction of the hits that make up the track are pion-like
  double piFrac = double(pihits)/ntrkhits;
  if(piFrac >=0.5){
    isPion = true;
  }
  if(isPion==true){
    return(trkphi);
  }
  else{
    return(-1000);
  }   


}
//maybe we can de-double the tracks. Like maybe if they share the same hits? cuz that is a sure sign that they are fake doubles
//basically, after isPion is declared true, loop over the rest of the hits and make sure we are not sharing hits


void RecoPerformance::processRunHeader( LCRunHeader* /*run*/) {
}

void RecoPerformance::processEvent( LCEvent * evt ) {
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.

  LCCollection* inputColMC = evt->getCollection(_inputCollectionNameMCP);
  LCCollection* inputColPFO = evt->getCollection(_inputCollectionNameR);
  LCCollection* inputColCl = evt->getCollection(_inputCollectionNameC);
  LCCollection* inputColTrk = evt->getCollection("SiTracks");
  LCCollection* inputSimHits = evt->getCollection(_inputCollectionSimHits);
  
  const int nTrk = inputColTrk->getNumberOfElements();
  const int nMCP = inputColMC->getNumberOfElements();
  std::vector<double> ptvec;
  std::vector<double> phivec;
  int nPiTrks=0;
  
  for(uint32_t i=0; i<nTrk; i++){
    const EVENT::Track *trk=static_cast<const EVENT::Track*>(inputColTrk->getElementAt(i));
    double phi = anaPiEff(trk,inputSimHits );
    if(phi==-1000){
      continue;
    }
    nPiTrks++;
    phivec.push_back(phi);
  }

  int nPiMC=0;
  std::vector<double> MCptvec;
  std::vector<double> MCphivec;
  for(uint32_t i=0; i<nMCP; i++){
    const EVENT::MCParticle *mcp=static_cast<const EVENT::MCParticle*>(inputColMC->getElementAt(i));
    int pdg = mcp->getPDG();
    if(abs(pdg)!=211){
      continue;
    }
    const EVENT::MCParticleVec mvec = mcp->getParents();
    int MID = mvec.back()->getPDG();
    if(MID != 15){
      continue;
    }
    const double* mom = mcp->getMomentum();
    double mcpt = std::sqrt(std::pow(mom[0], 2)+std::pow(mom[1], 2));
    double mcP = std::sqrt(std::pow(mcpt,2)+std::pow(mom[2],2));
    double mcphi = std::atan(mom[1]/mom[0]);
    double mctheta=std::acos(mom[2]/mcP);
    MCptvec.push_back(mcpt);
    MCphivec.push_back(mcphi);
    nPiMC++;
  }
    
  

  LCCollection* inputColTaus = evt->getCollection("TauRec_PFO");
  LCCollection* simTrkHits = evt->getCollection(_inputCollectionSimHits);


  if( inputColMC->getTypeName() != lcio::LCIO::MCPARTICLE ) {
    throw EVENT::Exception( "Invalid collection type: " + inputColMC->getTypeName() ) ;
  }
  if( inputColPFO->getTypeName() != lcio::LCIO::RECONSTRUCTEDPARTICLE ) {
    throw EVENT::Exception( "Invalid collection type: " + inputColPFO->getTypeName() ) ;
  }
  if( inputColCl->getTypeName() != lcio::LCIO::CLUSTER ) {
    throw EVENT::Exception( "Invalid collection type: " + inputColCl->getTypeName() ) ;
  }
  if( inputColTrk->getTypeName() != lcio::LCIO::TRACK ) {
    throw EVENT::Exception( "Invalid collection type: " + inputColTrk->getTypeName() ) ;
  }

  //find how many particles in the event
  //const int nMCP = inputColMC->getNumberOfElements();
  const int nPFO = inputColPFO->getNumberOfElements();
  const int nCl = inputColCl->getNumberOfElements();
  //const int nTrk = inputColTrk->getNumberOfElements();
  const int nTau = inputColTaus->getNumberOfElements();


  double nMCTau=0;

  for(uint32_t i=0; i<nMCP; i++){
    const EVENT::MCParticle *tau=static_cast<const EVENT::MCParticle*>(inputColMC->getElementAt(i));
    int pdg = tau->getPDG();
    if(abs(pdg)!=15){
      continue;
    }
    
    nMCTau++;

    //get pT of the tau
    const double* mom = tau->getMomentum();
    double pt = std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    //ptvec.push_back(pt);
  }

  

  double tauEff = nTau/nMCTau;
    
  

  //basic efficiencies (not PDG specific)
  double PFO_eff = double(nPFO)/nMCP;
  double Cl_eff = double(nCl)/nMCP;
  double Trk_eff = double(nTrk)/nMCP;

  

  _h_effPFO->Fill(PFO_eff);
  _h_effCl->Fill(Cl_eff);
  _h_effTrk->Fill(Trk_eff);
  



  int evtnum = evt->getEventNumber();
  /* std::fstream EffFile;
  EffFile.open("PTTrk.txt", std::ios::app);
  for(int i=0; i<ptvec.size(); i++){
    EffFile<<ptvec[i]<<',';
  }
  EffFile.close();
  std::fstream EffFile2;
  EffFile2.open("PTMC.txt", std::ios::app);
  for(int i=0; i<MCptvec.size(); i++){
    EffFile2<<MCptvec[i]<<',';
  }

  EffFile2.close();
  if(ptvec.size()>MCptvec.size()){
    std::cout<<"Track pTs: ";
    for(int i=0; i<ptvec.size();i++){
      std::cout<<ptvec[i]<<' ';
    }
    std::cout<<'\n'<<"MC pTs: ";
    for(int i=0; i<MCptvec.size();i++){
      std::cout<<ptvec[i]<<' ';
    }
    std::cout<<'\n';
  }
  */

  std::fstream MCFile;
  MCFile.open("MCPhi.txt", std::ios::app);
  for(int i=0; i<MCphivec.size(); i++){
    MCFile<<MCphivec[i]<<',';
  }
  MCFile.close();
  std::fstream TrkFile;
  TrkFile.open("TrkPhi.txt",std::ios::app);
  for(int i=0; i<phivec.size(); i++){
    TrkFile<<phivec[i]<<',';
  }
  TrkFile.close();
  
  
}
void RecoPerformance::check( LCEvent * /*evt*/ )
{ }

void RecoPerformance::end()
{ }
