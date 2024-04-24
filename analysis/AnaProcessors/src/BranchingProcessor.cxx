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
#include <algorithm>

BranchingProcessor aBranchingProcessor ;

BranchingProcessor::BranchingProcessor()
  : Processor("BranchingProcessor") {
  // modify processor description
  _description = "BranchingProcessor finds decay fractions." ;

  // register steering parameters: name, description, class-variable, default value

  /*********************************************************
   * Rose Powers (2023)
   * Use this processor to:
   * - Find the branching rate for various decays of a simulated MC particle
   * - Stream out to a text file that can then be analyzed with a quick python "pdgreader"
   * - With a little modification, could use it for reco products as well
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *********************************************************/

  registerProcessorParameter("PDG",
			     "PDG code of particle type",
			     _PDG,
			     _PDG);


  registerInputCollection( LCIO::MCPARTICLE,
			   "InputCollectionName" , 
			   "Name of an input collection.",
			   _inputCollectionName,
			   _inputCollectionName
			   );

  registerOutputCollection( LCIO::MCPARTICLE,
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

  /********************************************************
   *Decay modes for tau:
   *0: -211,16, 111 (pi-, pi0, nu-tau) (25.49%)
   *1: -211, 16 (pi- nu-tau) (10.82%)
   *2: -211, 16, 111, 111 (pi-, pi0x2, nu-tau) (9.26%)
   *3: -211, -211, 16, 211 (3-prong plus nu-tau) (8.99%)
   *4: -211, -211, 16, 111, 211 (3-prong plus pi0 and nu-tau) (2.74%)
   *5: -211, 16, 111, 111, 111 (3 pi0, pi-, nu-tau) (1.04%)
   *6: -12, 11, 16 (nu-tau, e, nu-e-bar) (17.82%)
   *7: -14, 13, 16 (nu-tau, e, nu-mu-bar) (17.39%)
   *
   *
   *
   *
   *
   *
   *
   *
   *
   ********************************************************/
  
  LCCollection* inputCol = evt->getCollection(_inputCollectionName);

  std::vector<int>decaymodes={};
  
  if( inputCol->getTypeName() != lcio::LCIO::MCPARTICLE ) {
    throw EVENT::Exception( "Invalid collection type: " + inputCol->getTypeName() ) ;
  }

  //find how many particles in the event
  const int nEl = inputCol->getNumberOfElements();
  //_h_N->Fill(nEl);
  
  std::vector<int> pdglist={0};

  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::MCParticle *mcp=static_cast<const EVENT::MCParticle*>(inputCol->getElementAt(i));
    
    //get pdgid, continue if not desired pdg
    double pdg = mcp->getPDG();
    //double pdg = mcp->getType();
    if(abs(pdg) != 15){
      continue;
      }

    //get the daughter particles
    const EVENT::MCParticleVec daughtervec = mcp->getDaughters();
    std::vector<int> daughterpdgs={0};
    if(daughtervec.size()==0){continue;}
    for(int n=0; n<daughtervec.size(); n++){
      daughterpdgs.push_back(daughtervec[n]->getPDG());
    }
    std::sort(daughterpdgs.begin(), daughterpdgs.end());
    daughterpdgs.erase(std::remove(daughterpdgs.begin(), daughterpdgs.end(), 0), daughterpdgs.end());
    for(int o=0; o<daughterpdgs.size(); o++){
      //std::cout<<daughterpdgs[o]<<std::endl;
    }
    if(daughterpdgs[0]==-211){
      //hadronic decays
      if(daughterpdgs[1]==-211){
	//3-prongs
	if(daughterpdgs[3]==211){
	  // 3-prong, no neutral pion
	  decaymodes.push_back(3);
	}
	else if(daughterpdgs[3]==111){
	  //3-prong with neutral pion
	  decaymodes.push_back(4);
	}
      }
      else if(daughterpdgs[1]==16){
	//
	if(daughterpdgs.size()==2){
	  decaymodes.push_back(1);
	}
	else if(daughterpdgs[2]==111){
	  if(daughterpdgs.size()==3){
	    decaymodes.push_back(0);
	  }
	  else if(daughterpdgs.size()==4){
	    decaymodes.push_back(2);
	  }
	  else if(daughterpdgs.size()==5){
	    decaymodes.push_back(5);
	  }
	}
      }
    }
    else if(daughterpdgs[0]==-12){
      decaymodes.push_back(6);
    }
    else if(daughterpdgs[0]==-14){
      decaymodes.push_back(7);
    }
       
	
	
	
    
    //Get the parent ID
    
    const EVENT::MCParticleVec parentvec = mcp -> getParents();
    if(parentvec.size()==0){
      continue;
    }
    int parentID = parentvec.back()->getPDG();
   
    // std::cout<<parentID<<std::endl;

    //If not direct tau decay, skip
    if(parentID != 15){
      continue;
      }
    

    //Fill charge
    int q = mcp->getCharge();
    _h_q->Fill(q);
    _h_MID->Fill(parentID);
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

  //stream info out to files for further analysis
  std::fstream pdgFile;
  pdgFile.open("pdgFile.txt", std::ios::app);
  for(int i=0; i<pdglist.size()-1; i++){
    pdgFile<<pdglist[i+1]<<'\n';
  }
  pdgFile.close();

  std::fstream decayFile;
  decayFile.open("decayFile.txt", std::ios::app);
  for(int i=0; i<decaymodes.size();i++){
    decayFile<<decaymodes[i]<<'\n';
  }
  decayFile.close();
  
  
  
}
void BranchingProcessor::check( LCEvent * /*evt*/ )
{ }

void BranchingProcessor::end()
{ }
