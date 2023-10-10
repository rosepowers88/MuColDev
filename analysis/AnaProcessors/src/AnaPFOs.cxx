#include "AnaProcessors/AnaPFOs.hxx"

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

AnaPFOs aAnaPFOs;

AnaPFOs::AnaPFOs()
  : Processor("AnaPFOs") {
  _description = "AnaPFOs deals with the PFO analysis, unsurprisingly" ; 
  /********************************
   *
   * Rose Powers (2023)
   * Use this processor for doing specific studies on a PFO collection
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

void AnaPFOs::init(){
  printParameters();

  marlin::AIDAProcessor::histogramFactory(this);
  _h_pion_misID = new TH1F("pion_ID", ";PDG ID", 10000, -500, 2500);
}

void AnaPFOs::processRunHeader( LCRunHeader * /*run*/){
}

void AnaPFOs::processEvent( LCEvent * evt) {
  LCCollection* MCPs = evt->getCollection(_inputcolMC);
  const int nMCP = MCPs->getNumberOfElements();
  LCCollection* taus;
  bool hastaus = true;
  try{
    taus = evt->getCollection("TauRec_PFO");
  }
  catch (Exception& e) {
    hastaus = false;
  }
  if(hastaus){
  //TAUS========================
  double nMCTau=0;
  double taupT=0;
  double tauTheta=0;
  for(uint32_t i=0; i<nMCP; i++){
    const EVENT::MCParticle *tau=static_cast<const EVENT::MCParticle*>(MCPs->getElementAt(i));
    int pdg = tau->getPDG();
    if(abs(pdg)!=15){
      continue;
    }
    const double* mom = tau->getMomentum();
    double pt = std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    taupT=pt;
    double theta = std::acos(mom[2]/std::sqrt(std::pow(pt,2)+std::pow(mom[2],2)));
    tauTheta = theta;
    break; //if we got a tau, we're done bc there is only one per event
  }
  const int nTau = taus->getNumberOfElements();
  std::fstream tauFile;
  tauFile.open("tauEff_Z0.txt", std::ios::app);
  tauFile<< taupT << ' ' << nTau << ' ' << tauTheta << '\n';
  tauFile.close();
  }
  bool badEvent = false;
  bool PFO = true;
  //get the collections
 
  LCCollection* PFOs;

  LCCollection* Trks = evt->getCollection("SiTracks_Refitted");

  try{
    PFOs=evt->getCollection(_inputcolPFO);
  }
  catch(Exception& e){
    PFO=false;
  }

  //get the sizes of the collections
 
  const int nTrk = Trks->getNumberOfElements();
  int goodMC=0;
  int goodPF=0;
  int goodTrk=0;

  //loop over MCP and get pTs,PDGs
  std::vector<std::vector<double>> MC_info;
  //also want to save all the hits and associate them with an MCP
  std::vector<std::vector<int>> mcpToSth;
  LCCollection* sthCol = evt->getCollection("AllSTH");
  const int nsth = sthCol->getNumberOfElements();

  
  std::vector<EVENT::MCParticle*> goodMCs;
  for(uint32_t i=0; i<nMCP; i++){
    std::vector<double> infovec; // each entry in MC_info is a vector doublet of (pT, PDG)
    EVENT::MCParticle *mcp = static_cast<EVENT::MCParticle*>(MCPs->getElementAt(i));
    const double* mom = mcp->getMomentum();
    double pT = std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    double theta = std::acos(mom[2]/std::sqrt(std::pow(pT,2)+std::pow(mom[2],2)));
    if(pT<_minPt){
      continue; //don't want a gazillion pieces of trash at MeV scale -- basically removing this particle
    }
    int q = mcp->getCharge();
    const int PDG = mcp->getPDG();
    //if(abs(PDG)!=211) continue;
    goodMC++;
    goodMCs.push_back(mcp);
    infovec.push_back(pT);

    infovec.push_back(PDG);
    infovec.push_back(q);
    infovec.push_back(theta);
    MC_info.push_back(infovec);
  }

  std::vector<int> simHitEffs;
  std::vector<int> pdgs;
  for(uint32_t i=0; i<goodMC; i++){
    EVENT::MCParticle *mcp = static_cast<EVENT::MCParticle*>(goodMCs[i]);
    pdgs.push_back(mcp->getPDG());
    std::vector<int> thIDs;
    for(uint32_t j=0; j<nsth; j++){
      EVENT::SimTrackerHit *sth = static_cast<EVENT::SimTrackerHit*>(sthCol->getElementAt(j));
      if(sth->getMCParticle() == mcp){
	int ID = sth->getCellID0();
	thIDs.push_back(ID);
      }
    }
    if(thIDs.size()==0){
      simHitEffs.push_back(0);
      //std::cout<<"MCP has no tracker hits!"<<std::endl;
      //std::cout<<"PDG: "<<mcp->getPDG()<<std::endl;
    }
    else{
      simHitEffs.push_back(1);
    }
    mcpToSth.push_back(thIDs);
  }

  /*
  std::fstream shfile;
  shfile.open("simhit.txt", std::ios::app);
  for(int i=0; i<goodMC; i++){
    shfile << MC_info[i][0] << ' '<< simHitEffs[i]<<' '<< pdgs[i]<<'\n';
  }
  shfile.close(); */


  std::vector<int> trkEffs;
  //okay, now get the track and check the hits against the simtrackerhits
  for(uint32_t i=0; i<goodMC; i++){
    int bestmatch=-1000;
    for(uint32_t l=0; l<nTrk; l++){
      int match=0;
      EVENT::Track *trk = static_cast< EVENT::Track*>(Trks->getElementAt(l));
      const EVENT::TrackerHitVec hits = trk->getTrackerHits();
      for(uint32_t j=0; j<hits.size(); j++){
	const EVENT::TrackerHit *th = static_cast<const EVENT::TrackerHit*>(hits[j]);
	int ID = th->getCellID0();
	for(int v=0; v<mcpToSth[i].size(); v++){
	  if(ID==mcpToSth[i][v]) match++;
	}
      }
      if(match>mcpToSth[i].size()*0.5){
	if(match>bestmatch) bestmatch=match;
      }
    }
    if(bestmatch>0) trkEffs.push_back(1);
    else trkEffs.push_back(0);
  }


      
  //okay so now we have trkeffs sorted by particle and pTs, so I feel like we can just write it out to a file now
  std::fstream trkFile;
  trkFile.open("trkEff_TEST.txt", std::ios::app);
  for(int i=0; i<goodMC; i++){
    double theta = MC_info[i][3];
    trkFile<<MC_info[i][0]<<' '<<trkEffs[i]<<' '<<pdgs[i]<<' '<<theta<<'\n';
  }
  trkFile.close();

  if(PFO){
    const int nPFO = PFOs->getNumberOfElements();
    //check for too many PFOs HERE (might fix the malloc issue)
    if(goodMC-nPFO<0){
      //std::cout<<"Too Many PFO: "<<evt->getEventNumber()<<std::endl;
      //std::cout<<nPFO<<", "<<goodMC<<std::endl;
      badEvent=true;
      return;
    }

    //now do the same for the PFOs
    std::vector<std::vector<double>> PF_info;
    std::vector<double> trkEff;
    std::vector<double> clEff;
    bool hasTrks=false;
    bool hasCl = false;
    for(uint32_t i=0; i<nPFO; i++){
      std::vector<double> infovec;
      const EVENT::ReconstructedParticle *pfo = static_cast<const EVENT::ReconstructedParticle*>(PFOs->getElementAt(i));
      int q = pfo->getCharge();
      //if(q == 0 ) continue;
      goodPF++;
      EVENT::TrackVec pfotrks = pfo->getTracks();
      EVENT::ClusterVec pfocl = pfo->getClusters();
      if(!pfotrks.size()){ hasTrks=false; }
      else{
	hasTrks=true;
      }
      if(!pfocl.size()){ hasCl = false; }
      else hasCl = true;
      trkEff.push_back(hasTrks);
      clEff.push_back(hasCl);
      const double* mom = pfo->getMomentum();
      double pT = std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
      infovec.push_back(pT);
      const int PDG = pfo->getType();
      infovec.push_back(PDG);
      PF_info.push_back(infovec);
    }

    //now we have all the info. We want to be equipped to make plots at different levels of intricacy
    //start with a very simple plot: calculate efficiency for each event as nPFO/nMCP. However we still need to match momenta
    //honestly let's do exactly what we did for the pion study, just with no PDG requirement. 

    //make an efficiency vector
    std::vector<double> effs;
    std::vector<double> truthpTs;
    std::vector<double> resolution; //while we're at it, might as well record resolution
    std::vector<double> truthPDGs; 
    std::vector<double> trkToPfo;
    std::vector<double> clToPfo;
    std::vector<double> truthQs;
    std::vector<double> truthThetas;
    std::vector<std::vector<double>> MC_temp = MC_info;
    if(goodPF==0 && goodMC>0) effs.push_back(0);
    for(int i=0; i<goodPF; i++){
      double candidate=10000;
      double thePT=0;
      double thePDG = 0;
      int erase=0;
      double res=0;
      double q = 0;
      double theta=0;
      for(int j=0; j<MC_temp.size(); j++){
	double diff = PF_info[i][0]-MC_temp[j][0];
	if(abs(diff)<candidate){
	  candidate=abs(diff);
	  res=diff;
	  thePT=MC_temp[j][0];
	  thePDG = MC_temp[j][1];
	  q=MC_temp[j][2];
	  theta=MC_temp[j][3];
	  erase=j;
	}
      }
      truthpTs.push_back(thePT);
      effs.push_back(1);
      resolution.push_back(res);
      truthPDGs.push_back(thePDG);
      trkToPfo.push_back(trkEff[i]);
      clToPfo.push_back(clEff[i]);
      truthQs.push_back(q);
      truthThetas.push_back(theta);
      MC_temp.erase(MC_temp.begin()+erase);
    }
 
    int iter = goodMC-goodPF;
    for(int i=0; i<iter; i++){
      effs.push_back(0);
      truthpTs.push_back(MC_temp[i][0]);
      truthPDGs.push_back(MC_temp[i][1]);
      truthQs.push_back(MC_temp[i][2]);
      truthThetas.push_back(MC_temp[i][3]);
      trkToPfo.push_back(0);
      clToPfo.push_back(0);
    }
  

    //we also want a separate diagnostic to see if the particle was misidentified. Check all the PFO's:
    //let's do a quick pion misidentification study.
    std::fstream misIDFile;
    misIDFile.open("misIDFile.txt", std::ios::app);
    std::vector<int> ID_true;
    for(int i=0; i<goodMC; i++){
      if(!effs[i]) ID_true.push_back(0);
      else{
	double reco_pdg = PF_info[i][1];
	if(reco_pdg!=truthPDGs[i]){
	  ID_true.push_back(0);
	  if(abs(truthPDGs[i])==211){
	    _h_pion_misID->Fill(reco_pdg);
	    misIDFile<<reco_pdg<<'\n';
	  }
	}
	else{
	  ID_true.push_back(1);
	}
      }
    }

    //finally, write out results to a file!
    if(!badEvent){ //only write out if we don't have the Problem (tm)
      std::fstream infoFile;
      infoFile.open("PFO_clStudy.txt", std::ios::app);
      for(int i=0; i<goodMC; i++){
	infoFile<<truthpTs[i]<<' '<<effs[i]<<' '<<ID_true[i]<<' '<<trkToPfo[i]<<' '<<truthQs[i]<<' '<<truthThetas[i]<<' '<<truthPDGs[i]<<' '<<clToPfo[i]<<'\n';
      }
      infoFile.close();
    }
    //write resolution to a file as well, but only if we we are misidentified
    std::fstream resFile;
    resFile.open("resolution_TEST.txt", std::ios::app);
    for(int i=0; i<goodMC; i++){
      if(!effs[i]) continue;
      if(ID_true[i]) continue; //only want the bad ones
      resFile<<truthpTs[i]<<' '<<resolution[i]<<'\n';
    }
    resFile.close();
  }
}
void AnaPFOs::check( LCEvent * /*evt*/ )
{}
void AnaPFOs::end()
{}
