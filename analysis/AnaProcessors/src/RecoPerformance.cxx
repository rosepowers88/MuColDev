#include "AnaProcessors/RecoPerformance.hxx"

#include <fstream>
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Vertex.h>

#include <marlin/AIDAProcessor.h>
#include <marlin/Statusmonitor.h>

#include <cmath>

RecoPerformance aRecoPerformance ;

RecoPerformance::RecoPerformance()
  : Processor("RecoPerformance") {
  _description = "RecoPerformance analyzes the performance of tracking and PFO reco" ;
  /******************************************
   * Rose Powers (2023)
   * Use this processor for:
   * - Efficiency studies
   * - Writing results to separate files
   * - Getting track truth info w/o LCIO relations
   * 
   *
   *
   *
   *
   *
   *
   **************************************************/

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
double RecoPerformance::anaPiEff( const LCObject * inputTrk, LCCollection* trkHitCol){
  //find the efficiency for charged pion FROM TAU reconstruction
  //trk part
  //first get the simhits and make a reference
  
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


}



void RecoPerformance::processRunHeader( LCRunHeader* /*run*/) {
}

void RecoPerformance::processEvent( LCEvent * evt ) {
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.
  bool debug = false;
  if(debug){
    std::cout<<"Started processing event "<<evt->getEventNumber()<<std::endl;
  }

  LCCollection* inputColMC = evt->getCollection(_inputCollectionNameMCP);
  LCCollection* inputColPFO = evt->getCollection(_inputCollectionNameR);
  LCCollection* inputColPFONew = evt->getCollection("MyPandoraPFOs");
  LCCollection* inputColCl = evt->getCollection(_inputCollectionNameC);
  LCCollection* inputColTrk = evt->getCollection("SiTracks");
  LCCollection* inputColTrkNew = evt->getCollection("MySiTracks");
  LCCollection* inputSimHits = evt->getCollection(_inputCollectionSimHits);

  const int nTrk = inputColTrk->getNumberOfElements();
  const int nNTrk = inputColTrkNew->getNumberOfElements();
  const int nMCP = inputColMC->getNumberOfElements();
  std::vector<double> ptvec;
  std::vector<double> ptvecNew;
  std::vector<int> nHitsVec;
  int nPiTrks=0;
  int nPiTrksNew=0;
  int nPiFromTau=0;
  int evtnum = evt->getEventNumber();
  std::vector<double> mcptvec_temp;
  for(uint32_t i=0; i<nMCP; i++){
    const EVENT::MCParticle *mcp = static_cast<const EVENT::MCParticle*>(inputColMC->getElementAt(i));
    if(abs(mcp->getPDG())!=211) continue;
    const EVENT::MCParticleVec mvec = mcp->getParents();
    int MID = mvec.back()->getPDG();
    if(MID != 15) continue;
    const double* vert = mcp->getVertex();
    double radius = std::sqrt(std::pow(vert[0],2)+std::pow(vert[1], 2)+std::pow(vert[2],2));
    //std::cout<<radius<<std::endl;
    nPiFromTau++;
    const double* mom = mcp->getMomentum();
    double pt = std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    mcptvec_temp.push_back(pt);
  }
  std::vector<float> effs;
  std::vector<float> trkPTs;
  int nPiPFO=0;
  double pT=0;
  double checker=0;
  for(uint32_t i=0; i<inputColPFONew->getNumberOfElements(); i++){
    const EVENT::ReconstructedParticle *pfo = static_cast<const EVENT::ReconstructedParticle*>(inputColPFONew->getElementAt(i));
    if(abs(pfo->getType())!=211){ 
      continue; }

    //const float* pos = pfo->getReferencePoint();
    //double radius = std::sqrt(std::pow(pos[0],2)+std::pow(pos[1], 2)+std::pow(pos[2],2));
    //std::cout<<radius<<std::endl;
    pT=0;
    const double* mom = pfo->getMomentum();
    pT = std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    const EVENT::TrackVec tracks = pfo->getTracks();
    for(uint32_t j=0; j<tracks.size(); j++){
      const EVENT::Track *trk = static_cast<const EVENT::Track*>(tracks[j]);
      checker = anaPiEff(trk,inputSimHits);
      if(checker<0) continue;
    }
    if(checker<0) continue;
    nPiPFO++;
    trkPTs.push_back(pT);
    
  }

  // if(trkRadii.size()!=0){
   
  //}

  //want to get the pT that is most accurate for each, if there are multiples
  //so for each, make a candidate of which one it reconstructed and then minimize it
  std::vector<double> ptreco;
  std::vector<double> resolution;
  double res;
  if(trkPTs.size()==0&&nPiFromTau>0) effs.push_back(0);
  for(int i=0; i<trkPTs.size(); i++){
    double candidate=10000;
    double thePT=0;
    int erase=0;
    if(mcptvec_temp.size()==0) continue;
    for(int j=0; j<mcptvec_temp.size(); j++){
      double diff = trkPTs[i]-mcptvec_temp[j];
      if(abs(diff)<candidate){
	candidate=abs(diff);
	res=diff;
	thePT=mcptvec_temp[j];
	erase=j;
      }
    }
    ptreco.push_back(thePT);
    effs.push_back(1);
    resolution.push_back(res);
    mcptvec_temp.erase(mcptvec_temp.begin()+erase);
  }

  int iter = nPiFromTau-trkPTs.size();
  if(iter<0){
    std::cout<<"DUPLICATE: "<<evt->getEventNumber()<<std::endl;
  }
  else{
  for(int i=0; i<iter; i++){
    effs.push_back(0);
    ptreco.push_back(mcptvec_temp[i]);
  }
  }
  


  std::fstream pfoFile;
  pfoFile.open("piPFOs_new.txt", std::ios::app);
  for(int i=0; i<ptreco.size(); i++){
    pfoFile<<ptreco[i]<<' '<<effs[i]<<'\n';
  }
  pfoFile.close();

  std::fstream resFile;
  resFile.open("piRes_new.txt", std::ios::app);
  for(int i=0; i<resolution.size(); i++){
    resFile<<ptreco[i]<<' '<<resolution[i]<<'\n';
  }
  resFile.close();



    
    
      
  

  //get a vector of the pTs for charged pion tracks to do a pT efficiency study
  for(uint32_t i=0; i<nTrk; i++){
    const EVENT::Track *trk=static_cast<const EVENT::Track*>(inputColTrk->getElementAt(i));
    double pT = anaPiEff(trk,inputSimHits );
    if(pT==-1000){
      continue;
    }
    nPiTrks++;
    ptvec.push_back(pT);
    const EVENT::TrackerHitVec trkhits = trk->getTrackerHits();
    nHitsVec.push_back(trkhits.size());
  }

  for(uint32_t i=0; i<nNTrk; i++){
    const EVENT::Track *ntrk = static_cast<const EVENT::Track*>(inputColTrkNew->getElementAt(i));
    double pT = anaPiEff(ntrk, inputSimHits);
    if(pT==-1000){
      continue;
    }
    nPiTrksNew++;
    ptvecNew.push_back(pT);
  }

  //get a vector of truth pTs for charged pions to do a pT efficiency study
  int nPiMC=0;
  std::vector<double> MCptvec;
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
    nPiMC++;
  }

    
  //just a lil fix so we don't get errors if we're running over a set with no reco'd taus
  bool hastaus = true;
  bool _hastaus = true;
  LCCollection* oldTaus;
  LCCollection* newTaus;
  try{
    oldTaus = evt->getCollection("TauRec_PFO");
  }
  catch (Exception& e) {
    hastaus=false;
  }
  try{
    newTaus = evt->getCollection("MyTauRec_PFO");
  }
  catch (Exception& e){
    _hastaus=false;
  }
  

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
  const int nPFO = inputColPFO->getNumberOfElements();
  const int nCl = inputColCl->getNumberOfElements();
  if(hastaus&&_hastaus){
  const int nOTau = oldTaus->getNumberOfElements();
  const int nNTau = newTaus->getNumberOfElements();
  if(nOTau!=nNTau){
    std::cout<<"old: "<<nOTau<<std::endl;
    std::cout<<"new: "<<nNTau<<std::endl;
  }


  double nMCTau=0;
  double taupT;
  for(uint32_t i=0; i<nMCP; i++){
    const EVENT::MCParticle *tau=static_cast<const EVENT::MCParticle*>(inputColMC->getElementAt(i));
    int pdg = tau->getPDG();
    if(abs(pdg)!=15){
      continue;
    }
    const double* mom = tau->getMomentum();
    double pt = std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    taupT=pt;
    break; //if we got a tau, we're done bc there is only one per event
  }
  int nEff=0;
  int oEff=0;

  if(nOTau==0)oEff=0;
  else oEff=1;
  if(nNTau==0)nEff=0;
  else nEff=1;

  //stream out the results for later analysis (see the python files)
  std::fstream tauFileOld;
  tauFileOld.open("tauEff_old_test.txt", std::ios::app);
  tauFileOld << taupT <<' '<< oEff << '\n';
  tauFileOld.close(); 

  std::fstream tauFileNew;
  tauFileNew.open("tauEff_new_test.txt", std::ios::app);
  tauFileNew << taupT << ' '<<nEff << '\n';
  tauFileNew.close();
  

  /* //basic efficiencies (not PDG specific)
  double PFO_eff = double(nPFO)/nMCP;
  double Cl_eff = double(nCl)/nMCP;
  double Trk_eff = double(nTrk)/nMCP;
  */
  
  /*
  _h_effPFO->Fill(PFO_eff);
  _h_effCl->Fill(Cl_eff);
  _h_effTrk->Fill(Trk_eff);
  */


  }
  /*
  int evtnum=evt->getEventNumber();

  //stream out the pi information for further analysis
  std::fstream MCFile;
  MCFile.open("MCpT_pi.txt", std::ios::app);
  for(int i=0; i<MCptvec.size(); i++){
    MCFile<<MCptvec[i]<<',';
  }
  MCFile.close();

  std::fstream TrkFile;
  TrkFile.open("TrkpT_pi_old.txt",std::ios::app);
  for(int i=0; i<ptvec.size(); i++){
    TrkFile<<ptvec[i]<<',';
  }
  TrkFile.close();

  std::fstream TrkFileNew;
  TrkFileNew.open("TrkpT_pi_new.txt",std::ios::app);
  for(int i=0; i<ptvecNew.size(); i++){
    TrkFileNew<<ptvecNew[i]<<',';
  }
  TrkFileNew.close();
  */
  
  //stream out the duplicate track candidates for further analysis (just comment out if you don't want to do this)
  /*
  int run = evt->getRunNumber();
  std::fstream dupefile;
  dupefile.open("dupeTracks.txt", std::ios::app);
  std::sort(ptvec.begin(),ptvec.end());
  if(ptvec.size()>MCptvec.size()){
    int ndupes=0;
    dupefile<<"=========="<<run<<":"<<evtnum<<"============"<<'\n';
    dupefile<<"nMC: "<<MCptvec.size()<<" nTrk: "<<ptvec.size()<<'\n';
  for(int i=0; i<ptvec.size(); i++){
    for(int j=0; j<ptvec.size(); j++){
      if(i<j && abs(ptvec[i]-ptvec[j])<=0.1){
	dupefile<<ptvec[i]<<' '<<ptvec[j]<<'\n';
	ndupes++;
	if(nHitsVec[i]>=nHitsVec[j]){
	  dupefile<<nHitsVec[j]<<',';
	}
	else{
	  dupefile<<nHitsVec[i]<<','; //write the smallest number of hits
	}
      }
    }
  }
  dupefile<<'\n';
  dupefile.close();
  }
  */
 
}
void RecoPerformance::check( LCEvent * /*evt*/ )
{ }

void RecoPerformance::end()
{ }
