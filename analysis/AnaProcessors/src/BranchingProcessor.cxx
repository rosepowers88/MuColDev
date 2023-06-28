#include "AnaProcessors/BranchingProcessor.hxx"

#include <fstream>
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>

#include <marlin/AIDAProcessor.h>
#include <marlin/Statusmonitor.h>

#include <cmath>

BranchingProcessor aBranchingProcessor ;

BranchingProcessor::BranchingProcessor()
  : Processor("BranchingProcessor") {
  // modify processor description
  _description = "BranchingProcessor makes histograms at the gen level." ;

  // register steering parameters: name, description, class-variable, default value


  registerProcessorParameter("PDG",
			     "PDG code of particle type",
			     _PDG,
			     _PDG);


  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputCollectionName" , 
			   "Name of an input collection.",
			   _inputCollectionName,
			   _inputCollectionName
			   );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "OutputCollectionName" ,
			    "Name of the output collection" ,
			    _outputCollectionName,
			    _outputCollectionName
			    );
  
}

void BranchingProcessor::init() {
  // Print the initial parameters
  printParameters() ;

  //Basic histograms
  marlin::AIDAProcessor::histogramFactory(this);
  _h_pdg  = new TH1F("PDG_ID", ";PDG ID", 1000, -500, 2500);
  _h_MID = new TH1F("Mother PDG", ";PDG ID", 1000, -500, 2500); //keep on for a check to make sure mother id vetting works
  _h_q = new TH1F("Charge", ";q [q_{e}]", 50, -2, 2);
}

void BranchingProcessor::processRunHeader( LCRunHeader* /*run*/) {
}

void BranchingProcessor::processEvent( LCEvent * evt ) {
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.
  // Loop over MCParticles
  LCCollection* inputCol = evt->getCollection(_inputCollectionName);


  if( inputCol->getTypeName() != lcio::LCIO::RECONSTRUCTEDPARTICLE ) {
    throw EVENT::Exception( "Invalid collection type: " + inputCol->getTypeName() ) ;
  }

  //find how many particles in the event
  const int nEl = inputCol->getNumberOfElements();
  //_h_N->Fill(nEl);
  
  std::vector<int> pdglist = {0};

  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::ReconstructedParticle *mcp=static_cast<const EVENT::ReconstructedParticle*>(inputCol->getElementAt(i));
    
    //get pdgid, continue if not desired pdg
    //double pdg = mcp->getPDG();
    double pdg = mcp->getType();
    /*if(abs(pdg) !=_PDG && pdg!=0){
      continue;
      }*/

    //Get the parent ID
    /* const EVENT::MCParticleVec parentvec = mcp -> getParents();
    if(parentvec.size()==0){
      continue;
    }
    int parentID = parentvec.back()->getPDG();
   
    // std::cout<<parentID<<std::endl;

    //If not direct tau decay, skip
    if(parentID != 15){
      continue;
      }*/
    

    //Fill charge
    int q = mcp->getCharge();
    _h_q->Fill(q);
    //_h_MID->Fill(parentID);
    _h_pdg->Fill(pdg);

    //if pdg not already in list, add
    bool pdg_in_list = false;
    for(int i; i<pdglist.size()-1; i++){
      if(pdglist[i]==pdg){
	pdg_in_list = true;
	break;
      }
    }
    if(pdg_in_list==false){
      pdglist.push_back(pdg);
      }

    
    
  }
  std::fstream pdgFile;
  pdgFile.open("pdgFile.txt", std::ios::app);
  for(int i=0; i<pdglist.size()-1; i++){
    pdgFile<<pdglist[i+1]<<'\n';
  }
  pdgFile.close();
  
  
  
  
}
void BranchingProcessor::check( LCEvent * /*evt*/ )
{ }

void BranchingProcessor::end()
{ }
