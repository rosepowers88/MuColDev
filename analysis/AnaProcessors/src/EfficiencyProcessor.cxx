#include "AnaProcessors/EfficiencyProcessor.hxx"

#include <fstream>
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>

#include <marlin/AIDAProcessor.h>
#include <marlin/Statusmonitor.h>

#include <cmath>

EfficiencyProcessor aEfficiencyProcessor ;

EfficiencyProcessor::EfficiencyProcessor()
  : Processor("EfficiencyProcessor") {
  _description = "EfficiencyProcessor analyzes the performance of tracking and PFO reco and outputs root file for plotting" ;
  /******************************************
   * Rose Powers (2023)
   * Use this processor for:
   * - Efficiency studies
   * - Configured for taus but can change easily
   * - Writing out to a root file for plotmaking
   **************************************************/

  // register steering parameters

  registerProcessorParameter("MinPt",
			     "Minimum transverse momentum",
			     _minPt,
			     _minPt);
  registerProcessorParameter("MinTheta",
			     "Minimum polar angle",
			     _minTheta,
			     _minTheta);

//the relevant collections
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, //this will be PFO
			  "InputCollectionNameReco",
			  "Name of RecoParticle collection.",
			  _inputCollectionNameR,
			  _inputCollectionNameR
			  );
  registerInputCollection(LCIO::MCPARTICLE,
			  "InputCollectionNameMCP",
			  "Name of MCP Collection.",
			  _inputCollectionNameMCP,
			  _inputCollectionNameMCP
			  );
}

void EfficiencyProcessor::init() {
  // Print the initial parameters
  printParameters() ;

  //The histograms we want for efficiency. 
  marlin::AIDAProcessor::histogramFactory(this);
  //pT pass and all
  _h_pass_pT = new TH1F("pass, pT", "", 50,0,300);
  _h_all_pT = new TH1F("all, pT", "", 50,0,300);

  //theta pass and all
  _h_pass_theta = new TH1F("pass, theta", "", 50,-3.2,3.2);
  _h_all_theta = new TH1F("all, theta", "", 50, -3.2, 3.2);

}

void EfficiencyProcessor::processRunHeader( LCRunHeader* /*run*/) {
}

void EfficiencyProcessor::processEvent( LCEvent * evt ) {
  //


  LCCollection* inputColMC = evt->getCollection(_inputCollectionNameMCP);
  LCCollection* inputColPFO = evt->getCollection(_inputCollectionNameR);

  LCCollection* trkCol = evt->getCollection("SiTracks_Refitted");
  LCCollection* trk_to_MC;
  bool hasrel= true;
  try{
    trk_to_MC=evt->getCollection("MCParticle_SiTracks_Refitted");
  }
  catch(EVENT::Exception& e){
    hasrel=false;
  }
 
  //loop over charged MCPs and see if they have an associated track
  for(uint32_t i=0; i< inputColMC->getNumberOfElements(); i++){
    const EVENT::MCParticle *mc = static_cast<const EVENT::MCParticle*>(inputColMC->getElementAt(i));
    const int q = mc->getCharge();
    if(!q) continue;
    const int pdg = mc->getPDG();
    if(abs(pdg) != 211) continue; //select pions
    const double * p = mc->getMomentum();
    double pt = std::sqrt(std::pow(p[0],2)+std::pow(p[1],2));
    const double theta = std::acos(p[2]/std::sqrt(std::pow(pt,2)+std::pow(p[2],2)));
    //use LCRelations
    if(hasrel){
    for(uint32_t t=0; t < trk_to_MC->getNumberOfElements(); t++){
      const EVENT::LCRelation * rel = static_cast<const EVENT::LCRelation*>(trk_to_MC->getElementAt(t));
      const EVENT::MCParticle *relmc = static_cast<const EVENT::MCParticle*>(rel->getFrom());
      if(relmc=mc){
	//fill pass histograms
	_h_pass_pT->Fill(pt);
	_h_pass_theta->Fill(theta);
	break;
      }
    }
    }
    //fill "denominator" histograms
      _h_all_pT->Fill(pt);
      _h_all_theta->Fill(theta);
  }

  //end of track matching section, everything below for TAU REC


  /*
  
  const int nMCP = inputColMC->getNumberOfElements();
  std::vector<double> ptvec;
  int nPiTrks=0;
  //just a lil fix so we don't get errors if we're running over a set with no reco'd taus
  bool hastaus = true;
  LCCollection* inputColTaus;
  try{
    inputColTaus = evt->getCollection("TauRec_PFO");
  }
  catch (Exception& e) {
    hastaus=false;
  }



  //find how many particles in the event
  const int nPFO = inputColPFO->getNumberOfElements();
  if(hastaus){
  const int nTau = inputColTaus->getNumberOfElements();
  if(nTau > 1){
    std::cout<<"Extra taus"<<std::endl;
  }

  double nMCTau=0;
  double taupT;
  double tautheta;
  //get the MC tau and its pT and theta
  for(uint32_t i=0; i<nMCP; i++){
    const EVENT::MCParticle *tau=static_cast<const EVENT::MCParticle*>(inputColMC->getElementAt(i));
    int pdg = tau->getPDG();
    if(abs(pdg)!=15){
      continue;
    }
    const double* mom = tau->getMomentum();
    double pt = std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    taupT=pt;
    double mctheta=std::acos(mom[2]/std::sqrt(std::pow(pt,2)+std::pow(mom[2],2)));
    tautheta=mctheta;
    break; //if we got a tau, we're done bc there is only one per event
  }
    if(nTau>0){
      // _h_pass_pT->Fill(taupT);
      // _h_pass_theta->Fill(tautheta);
  }
    // _h_all_pT->Fill(taupT);
    //_h_all_theta->Fill(tautheta);



  


    } */
  
}
void EfficiencyProcessor::check( LCEvent * /*evt*/ )
{ }

void EfficiencyProcessor::end()
{ }
