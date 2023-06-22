#include "AnaProcessors/GenHistProcessor.hxx"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>

#include <marlin/AIDAProcessor.h>
#include <marlin/Statusmonitor.h>

#include <cmath>

GenHistProcessor aGenHistProcessor ;

GenHistProcessor::GenHistProcessor()
  : Processor("GenHistProcessor") {
  // modify processor description
  _description = "GenHistProcessor makes histograms at the gen level." ;

  // register steering parameters: name, description, class-variable, default value


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

void GenHistProcessor::init() {
  // Print the initial parameters
  printParameters() ;

  //Basic histograms
  marlin::AIDAProcessor::histogramFactory(this);
  _h_pt   = new TH1F("pT", ";pT [GeV]", 1000,0., 200);
  _h_p    = new TH1F("p", ";p [GeV]", 1000, 0, 500);
  _h_pdg  = new TH1F("PDG_ID", ";PDG ID", 1000, -500, 2500);
  _h_N    = new TH1F("nParticles", ";N", 50, 0, 1000);
  _h_E    = new TH1F("Energy", ";E [GeV]", 1000,0,500);

  _h_phi  = new TH1F("Azimuth", ";#phi [rad]", 50,-3.2,3.2); //defined between -pi and pi
  _h_theta = new TH1F("Polar", ";#theta [rad]", 50,0,3.2);
  _h_charge = new TH1F("Charge", ";q [q_{e}]", 20,-2,2);
  _h_vertex = new TH1F("Vertex distance", ";[mm]",1000,0,5000);
  _h_MID = new TH1F("Mother PDG", ";PDG ID", 1000, -500, 2500); //keep on for a check to make sure mother id vetting works
}

void GenHistProcessor::processRunHeader( LCRunHeader* /*run*/) {
}

void GenHistProcessor::processEvent( LCEvent * evt ) {
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.
  // Loop over MCParticles

  //map from pdg entry to actual pdg number for submitting in batch mode
  int pdgs[4]={0,15,211,111};
  
  LCCollection* inputCol = evt->getCollection(_inputCollectionName);


  if( inputCol->getTypeName() != lcio::LCIO::MCPARTICLE ) {
    throw EVENT::Exception( "Invalid collection type: " + inputCol->getTypeName() ) ;
  }

  //find how many particles in the event
  const int nEl = inputCol->getNumberOfElements();
  //_h_N->Fill(nEl);

  

  const double* mom;
  
  bool nELFilled = false;
  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::MCParticle *mcp=static_cast<const EVENT::MCParticle*>(inputCol->getElementAt(i));
    
    //get pdgid, continue if not desired pdg
    double pdg = mcp->getPDG();
    if(abs(pdg) !=pdgs[_PDG] && _PDG !=0 ){
      continue;
    }

    //Get the parent ID if we are not examining taus only or all
    if(pdgs[_PDG] != 15 && _PDG !=0){
      const EVENT::MCParticleVec parentvec = mcp -> getParents();
      int parentID = parentvec.back()->getPDG();

      //If not direct tau decay, skip
      if(parentID != 15){
	continue;
      }
      _h_MID->Fill(parentID);
    }

    //Fill charge and type
    int q = mcp->getCharge();
    _h_charge->Fill(q);
    _h_pdg->Fill(pdg);


    //fill momentum
    mom = mcp->getMomentum();
    double pt=std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    double p =std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2))+std::sqrt(std::pow(mom[2],2));
    _h_p->Fill(p);
    _h_pt->Fill(pt);
    
    //fill energy
    double E=mcp->getEnergy();
    _h_E->Fill(E);

    //find and fill azimuth and polar angle

    double phi = std::atan(mom[1]/mom[0]);
    _h_phi->Fill(phi);

    double theta = std::acos(mom[2]/p);
    _h_theta->Fill(theta);

    //fill distance from vertex

    
    const double *vertex = mcp->getVertex();
    double vert=std::sqrt(std::pow(vertex[0],2)+std::pow(vertex[1],2)+std::pow(vertex[2],2));
    _h_vertex->Fill(vert);

    //fill nEl once per cycle -- get accurate read on number of **relevant** elements in a collection (so after all the filter checks)
    if(nELFilled == false){
      _h_N->Fill(nEl);
      nELFilled = true;
    }
    
    
  }
}

void GenHistProcessor::check( LCEvent * /*evt*/ )
{ }

void GenHistProcessor::end()
{ }
