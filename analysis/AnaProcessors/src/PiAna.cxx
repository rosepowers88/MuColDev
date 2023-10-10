#include "AnaProcessors/PiAna.hxx"

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

PiAna aPiAna;

PiAna::PiAna()
  : Processor("PiAna") {
  _description = "PiAna deals with the pion analysis, unsurprisingly" ; 
  /********************************
   *
   * Rose Powers (2023)
   * Use this processor for doing specific studies on pion reco
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

void PiAna::init(){
  printParameters();

  marlin::AIDAProcessor::histogramFactory(this);
  // _h_pion_misID = new TH1F("pion_ID", ";PDG ID", 10000, -500, 2500);
  _h_E_neutron = new TH1F("Neutron energy MC", ";truth E [GeV]", 100,0,75);
  _h_E_PFneutron = new TH1F("Neutron energy PFO", ";reco E [GeV]", 100, 0, 75);
  _h_E_pi0 = new TH1F("Pi-0 energy MC", ";truth E [GeV]", 100,0,75);

}


void PiAna::processRunHeader( LCRunHeader * /*run*/){
}

void PiAna::processEvent( LCEvent * evt) {
  LCCollection* MCPs = evt->getCollection(_inputcolMC);
  const int nMCP = MCPs->getNumberOfElements();
  LCCollection* taus;
  LCCollection* PFOs = evt->getCollection(_inputcolPFO);
  const int nPFO = PFOs->getNumberOfElements();
  LCCollection* allCL = evt->getCollection("PandoraClusters");

 

  int nnp = 0;
  int ncp = 0;
  double npE=0;
  int nneutron = 0;
  for(uint32_t i=0; i<nMCP; i++){
    const EVENT::MCParticle *mcp = static_cast<const EVENT::MCParticle*>(MCPs->getElementAt(i));
    const int mcpdg = mcp->getPDG();
    if(mcpdg == 2112){
      nneutron++;
      _h_E_neutron->Fill(mcp->getEnergy());
    }
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
	 _h_E_pi0->Fill(E);
      }
    }
  }

 
  if(!nnp) return; //next event if no neutral pions in decay
  //std::cout<<"MCP: "<<nneutron<<std::endl;

  int nRecoNeutron = 0;
 for(uint32_t i=0; i<nPFO; i++){
    const EVENT::ReconstructedParticle *pfo = static_cast<const EVENT::ReconstructedParticle*>(PFOs->getElementAt(i));
    const int pfpdg = pfo->getType();
    if(pfpdg==2112){
      nRecoNeutron++;
      _h_E_PFneutron->Fill(pfo->getEnergy());
    }
 }
 //std::cout<<"Reco: "<<nRecoNeutron<<std::endl;

 if(nRecoNeutron>nneutron){
   //std::cout<<" delta neutron:"<< nRecoNeutron-nneutron<<std::endl;
   // std::cout<<"n pi 0: "<<nnp<<std::endl;
 }
}
void PiAna::check( LCEvent * /*evt*/ )
{}
void PiAna::end()
{}

