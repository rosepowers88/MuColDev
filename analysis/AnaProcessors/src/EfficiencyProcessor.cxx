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
  _h_pass_theta = new TH1F("pass, theta", "", 100,-6.4,6.4);
  _h_all_theta = new TH1F("all, theta", "", 100, -6.4, 6.4);

  //phi pass and all
  _h_pass_phi = new TH1F("pass, phi", "", 100, -6.4, 6.4);
  _h_all_phi = new TH1F("all, phi", "", 100, -6.4, 6.4);

  //2d hists
  _h_pt_vs_theta = new TH2F("pT_vs_Theta", ";pT[GeV];#theta[rad]", 50,0,300,50,-6.4,6.4);


}

void EfficiencyProcessor::processRunHeader( LCRunHeader* /*run*/) {
}

void EfficiencyProcessor::processEvent( LCEvent * evt ) {
 
      const double pi = std::acos(-1);
      //begin tracking efficiency here
  LCCollection* inputColMC = evt->getCollection(_inputCollectionNameMCP);
  LCCollection* inputColPFO = evt->getCollection(_inputCollectionNameR);
  LCCollection* clusterCol = evt->getCollection("PandoraClusters");
  LCCollection* trkCol = evt->getCollection("SiTracks_Refitted");
 


  //loop over charged MCPs and see if they have an associated track
  for(uint32_t i=0; i< inputColMC->getNumberOfElements(); i++){
    const EVENT::MCParticle *mc = static_cast<const EVENT::MCParticle*>(inputColMC->getElementAt(i));
    const int q = mc->getCharge();
    if(!q) continue;
    const int pdg = mc->getPDG();
    if(abs(pdg) != 211) continue; //select pions
    //now select pions that specifically came from taus
    // const EVENT::MCParticleVec parentvec = mc->getParents();
    //if(parentvec.back()->getPDG() != 15) continue;
    //const EVENT::MCParticle *tau = static_cast<const EVENT::MCParticle*>(parentvec.back());
    //now cut on tau theta
    // const double * tau_p = tau->getMomentum();
    //double tau_theta = std::acos(tau_p[2]/std::sqrt(std::pow(tau_p[0],2)+std::pow(tau_p[1],2)+std::pow(tau_p[2],2)));
    //if(tau_theta > 0.75 ) continue;
    const double * p = mc->getMomentum();
    //cut on pt
    double pt = std::sqrt(std::pow(p[0],2)+std::pow(p[1],2));
    if(pt < _minPt) continue;
    //const double theta = std::acos(p[2]/std::sqrt(std::pow(pt,2)+std::pow(p[2],2)));
    //cut on theta
    const double theta = std::acos(p[2]/std::sqrt(std::pow(pt,2)+std::pow(p[2],2)));
    if(theta < 0.5 || theta > 3) continue;
    const double phi = std::acos(p[0]/pt);
    const double * vert = mc->getVertex();
    const double z = vert[2];
    //use geometrical matching
    double min_dR = 100;
    double d0=0;
    double z0=0;
    int matchID = 0;

    for(uint32_t t= 0; t<trkCol->getNumberOfElements(); t++){
      const EVENT::Track *trk = static_cast<const EVENT::Track*>(trkCol->getElementAt(t));
      //get trkphi
      const double trkphi = trk->getPhi();
      //get theta from dip angle -- Lambda = pi/2 - theta --> theta = pi/2-lambda

      const double trktheta = pi/2-std::atan(trk->getTanLambda());
      //find dR
      double dR = std::sqrt(std::pow(fabs(phi)-fabs(trkphi),2)+std::pow(fabs(theta)-fabs(trktheta),2));
      if(dR < min_dR){
	min_dR = dR;
	d0=trk->getD0();
	z0=trk->getZ0();
	matchID = t;
      }
      
    }
    if(min_dR < 0.005){
      //fill pass histograms (uncomment for track matching)
      //_h_pass_pT->Fill(pt);
      //_h_pass_theta->Fill(theta);
      // _h_pass_phi->Fill(phi);
      //_h_pass_z->Fill(z);

      }
    else {
      // fill here histograms of quantities of the particles that DON't get tracked
    }

    //fill "denominator" histograms (uncomment for track matching)
    //_h_all_pT->Fill(pt);
    //_h_all_theta->Fill(theta);
    ////_h_all_phi->Fill(phi);
    //_h_all_z->Fill(z);
  }

 
  //end of track matching section

  //Beginning of tau rec matching (only relevant if TauRecProcessor has been run)  
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
    //std::cout<<"no taus"<<std::endl;
    hastaus=false;
  }



  //find how many particles in the event
  const int nPFO = inputColPFO->getNumberOfElements();
  if(hastaus){
    // const int nTau = inputColTaus->getNumberOfElements();
    //if(nTau > 1){
    //std::cout<<"Extra taus"<<std::endl;
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
  //if(nTau>0){
      //_h_pass_pT->Fill(taupT);
      // _h_pass_theta->Fill(tautheta);
  // }
    //_h_all_pT->Fill(taupT);
    //_h_all_theta->Fill(tautheta)
 
  //end of tau rec matching 
  

  //beginning of CLUSTER MATCHING
  //perform cluster efficiency by doing dR matching for clusters, measure against theta and energy
  LCCollection* clusts = evt->getCollection("PandoraClusters");
  for(uint32_t mcIter = 0; mcIter < inputColMC->getNumberOfElements(); mcIter++){
    const EVENT::MCParticle *mc = static_cast<const EVENT::MCParticle*>(inputColMC->getElementAt(mcIter));
    if(abs(mc->getPDG()) != 211 && abs(mc->getPDG()) != 111) continue; // include neutral pions
    const double mcE = mc->getEnergy();
    const double * mcmom = mc->getMomentum();
    const double mcpt = std::sqrt(mcmom[0]*mcmom[0]+mcmom[1]*mcmom[1]);
    if(mcpt < 50) continue;
    const double mc_theta = std::acos(mcmom[2]/std::sqrt(mcpt*mcpt + mcmom[2]*mcmom[2]));
    const double mc_eta = eta(mc_theta);
    const double mc_phi = std::acos(mcmom[0]/std::sqrt(mcmom[0]*mcmom[0]+mcmom[1]*mcmom[1]));
    double mindR = 1000;
    for(uint32_t cl_Iter = 0; cl_Iter < clusts->getNumberOfElements(); cl_Iter++){
      const EVENT::Cluster *cl = static_cast<const EVENT::Cluster*>(clusts->getElementAt(cl_Iter));
      const double cl_theta = cl->getITheta();
      const double cl_eta = eta(cl_theta);
      const double cl_phi = cl->getIPhi();
      double dR = std::sqrt(std::pow(mc_phi-cl_phi,2)+std::pow(mc_eta-cl_eta,2));
      if(dR < mindR){
	mindR = dR;
      }

    }
    if(mindR < 0.1 ){
      _h_pass_pT->Fill(mcE);
      _h_pass_theta->Fill(mc_eta);
    }
    _h_all_pT->Fill(mcE);
    _h_all_theta->Fill(mc_eta);


  }
  //end of cluster matching

  //beginning of PFO MATCHING
  // LCCollection* newPFOs = evt->getCollection("PandoraPFOs_test_2"); //uncomment if comparing two PFO collections
  LCCollection* PFOs = evt->getCollection("PandoraPFOs_test"); //modified PFO collection
  for(uint32_t mcIter = 0; mcIter < inputColMC->getNumberOfElements(); mcIter++){
    const EVENT::MCParticle *mc = static_cast<const EVENT::MCParticle*>(inputColMC->getElementAt(mcIter));
    if(abs(mc->getPDG()) != 211) continue;
    const double mcE = mc->getEnergy();
    const double * mcmom = mc->getMomentum();
    const double mcpt = std::sqrt(std::pow(mcmom[0],2)+std::pow(mcmom[1],2));
    double dpt;
    if(mcpt < 7 || mcpt > 200) continue;
    const double mc_theta = std::acos(mcmom[2]/std::sqrt(std::pow(mcpt,2)+std::pow(mcmom[2],2)));
    const double mc_eta = eta(mc_theta);
    if(0.03 > mc_theta || mc_theta > 2.75) continue;
    const double mc_phi = std::acos(mcmom[0]/mcpt);
    double mindR = 1000;
    for(uint32_t pfIter= 0; pfIter <PFOs->getNumberOfElements(); pfIter++){
      const EVENT::ReconstructedParticle *pf_n = static_cast<const EVENT::ReconstructedParticle*>(PFOs->getElementAt(pfIter));
 if(abs(pf_n->getType()) != 211) continue;
      const EVENT::ClusterVec clv = pf_n->getClusters();
      double dE = 0;
      for(int i=0; i<clv.size(); i++){
	const EVENT::Cluster *clv_c = static_cast<const EVENT::Cluster*>(clv[i]);
	double Ediff = clv_c->getEnergy() - mcE;
	if(fabs(Ediff) > dE) dE = Ediff;
    }
      const double * pfmom = pf_n->getMomentum();
      const double pfpt = std::sqrt(std::pow(pfmom[0],2)+std::pow(pfmom[1],2));
      const double pf_theta = std::acos(pfmom[2]/std::sqrt(std::pow(pfpt,2)+std::pow(pfmom[2],2)));
      const double pf_eta = eta(pf_theta);
      const double pf_phi = std::acos(pfmom[0]/pfpt);
      double dR = std::sqrt(std::pow(mc_phi-pf_phi,2)+std::pow(mc_eta-pf_eta,2));
      if(dR < mindR){
	mindR = dR;
	dpt = fabs(pfpt-mcpt);
      }

    }
    if(mindR < 0.1 ){
      //fill pass hists, uncomment if doing PFO eff
      //  _h_pass_pT->Fill(mcpt);
      //_h_pass_theta->Fill(mc_eta);
    }
    //fill denominator hists, uncomment if doing PFO eff
    //_h_all_pT->Fill(mcpt);
    //_h_all_theta->Fill(mc_eta);

    // _h_pt_vs_theta->Fill(mcpt,mc_theta); //uncomment for 2d hist
  }
}


const double EfficiencyProcessor::eta(const double theta){
  //function to extract pseudorapidity from polar angle theta
  return(-std::log(std::tan(theta/2)));

}

void EfficiencyProcessor::check( LCEvent * /*evt*/ )
{ }

void EfficiencyProcessor::end()
{ }
