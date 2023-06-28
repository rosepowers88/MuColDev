#include "AnaProcessors/RecoPerformance.hxx"

#include <fstream>
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>

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

  //Graphs
  _gr_PFO_Eff = new TGraph();
  _gr_PFO_Eff->SetTitle("PFO Efficiency;Energy [GeV]; nPFO/nMCP");
  _gr_Cl_Eff = new TGraph();
  _gr_Trk_Eff = new TGraph();
  
}

void RecoPerformance::processRunHeader( LCRunHeader* /*run*/) {
}

void RecoPerformance::processEvent( LCEvent * evt ) {
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.

  LCCollection* inputColMC = evt->getCollection(_inputCollectionNameMCP);
  LCCollection* inputColPFO = evt->getCollection(_inputCollectionNameR);
  LCCollection* inputColCl = evt->getCollection(_inputCollectionNameC);
  LCCollection* inputColTrk = evt->getCollection(_inputCollectionNameT);


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
  const int nMCP = inputColMC->getNumberOfElements();
  const int nPFO = inputColPFO->getNumberOfElements();
  const int nCl = inputColCl->getNumberOfElements();
  const int nTrk = inputColTrk->getNumberOfElements();

  //basic efficiencies (not PDG specific)
  double PFO_eff = double(nPFO)/nMCP;
  double Cl_eff = double(nCl)/nMCP;
  double Trk_eff = double(nTrk)/nMCP;

  _h_effPFO->Fill(PFO_eff);
  _h_effCl->Fill(Cl_eff);
  _h_effTrk->Fill(Trk_eff);

  
  
  std::vector<int> pdglist = {0};
  //lets start off by getting particle-by-particle efficiency for just charged pions and neutral pions
  int nChargedPi = 0;
  int nNeutralPi = 0;
  int nNeutron = 0;
  int nGam = 0;
  std::vector<double> MCCPenergies;
  std::vector<double> MCNPenergies;

  double chargedEnergy = 0;
  double neutralEnergy = 0;

  for(uint32_t i=0;i<nMCP;i++) {
    const EVENT::MCParticle *mcp=static_cast<const EVENT::MCParticle*>(inputColMC->getElementAt(i));
    
    //get pdgid, continue if not desired pdg
    double pdg = mcp->getPDG();
    if(abs(pdg)==211){
      nChargedPi++;
      MCCPenergies.push_back(mcp->getEnergy());
    }

    if(abs(pdg)==111){
      nNeutralPi++;
      MCNPenergies.push_back(mcp->getEnergy());
    }
    else{ continue; }
    if(pdg==2112){
      nNeutron++;
    }
    if(pdg==22){
      nGam++;
    }
  }

  if(MCCPenergies.size()!=0){
    //chargedEnergy = std::reduce(MCCPenergies.begin(), MCCPenergies.end())/MCCPenergies.size();
  }
  if(MCNPenergies.size()!=0){
    //neutralEnergy = std::reduce(MCNPenergies.begin(), MCNPenergies.end())/MNCPenergies.size();
  }

  int n211PFOs = 0;
  int n111PFOs = 0;
  int n2112PFOs = 0;
  int n22PFOs = 0;

   for(uint32_t i=0;i<nPFO;i++) {
    const EVENT::ReconstructedParticle *pfo=static_cast<const EVENT::ReconstructedParticle*>(inputColPFO->getElementAt(i));

    double PFO_pdg = pfo->getType();
    if(abs(PFO_pdg)==211){
      n211PFOs++;
    }
    if(abs(PFO_pdg)==111){
      n111PFOs++;
    }
    if(abs(PFO_pdg)==2112){
      n2112PFOs++;
    }
    if(abs(PFO_pdg)==22){
      n22PFOs++;
    }
    
    else{continue;}
   }
   double effCP = 0;
   double effNP = 0;
   double effN = 0;
   double effGam = 0;

   if(nChargedPi!=0){
     effCP = double(n211PFOs)/nChargedPi;
     _h_effPFO_cp->Fill(effCP);
   }
   if(nNeutralPi!=0){
     effNP = double(n111PFOs)/nNeutralPi;
     _h_effPFO_np->Fill(effNP);
   }
   if(nNeutron!=0){
     effN = double(n2112PFOs)/nNeutron;
     _h_effPFO_n->Fill(effN);
   }
   if(nGam!=0){
     effGam = double(n22PFOs)/nGam;
     _h_effPFO_gam -> Fill(effGam);
   }

   //Fill efficiency and energy information to graphs (not perfect but a first go)
   //_gr_PFO_Eff->SetPoint(_gr_PFO_Eff->getN(),chargedEnergy, effC);
 
  
  
  
}
void RecoPerformance::check( LCEvent * /*evt*/ )
{ }

void RecoPerformance::end()
{ }
