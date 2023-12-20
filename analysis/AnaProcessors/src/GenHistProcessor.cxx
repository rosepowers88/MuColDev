#include "AnaProcessors/GenHistProcessor.hxx"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>

#include <marlin/AIDAProcessor.h>
#include <marlin/Statusmonitor.h>

#include <cmath>

GenHistProcessor aGenHistProcessor ;

GenHistProcessor::GenHistProcessor()
  : Processor("GenHistProcessor") {
  // modify processor description
  _description = "GenHistProcessor makes histograms at the gen level." ;

  /******************************************************
   * Rose Powers (2023)
   * Use this processor to:
   * - get information about a sample at the MC truth (gen) level
   * - can restrict to a certain particle (i.e. taus, charged pions from taus, etc)
   *
   *
   *
   *
   *
   *
   *
   ******************************************************/

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
  _h_pt   = new TH1F("pT", ";pT [GeV]", 100,0., 250);
  _h_p    = new TH1F("p", ";p [GeV]", 100, 0, 250);
  _h_pdg  = new TH1F("PDG_ID", ";PDG ID", 1000, -500, 2500);
  _h_N    = new TH1F("nParticles", ";N", 50, 0, 1000);
  _h_E    = new TH1F("Energy", ";E [GeV]", 100,0,750000);
  _h_inv_E = new TH1F("Invariant Energy - Set Energy", ";#Delta E [GeV]", 200,-500.,500.);
  _h_inv_msqr = new TH1F("Invariant Mass Squared", ";m^2 [GeV^2]", 200,-100000,50.);
  _h_mass = new TH1F("Set Mass", ";m [GeV]", 100,0,10);
  _h_inv_p = new TH1F("Invariant Momentum - Set Momentum", ";#Delta p [GeV]", 200,-500,100);

  _h_phi  = new TH1F("Azimuth", ";#phi [rad]", 50,-6.4,6.4); //defined between -pi and pi
  _h_theta = new TH1F("Polar", ";#theta [rad]", 50,-6.4,6.4);
  _h_charge = new TH1F("Charge", ";q [q_{e}]", 20,-2,2);
  _h_vertex = new TH1F("Vertex distance", ";[mm]",1000,0,50);
  _h_vert_x = new TH1F("Vertex x", ";[mm]", 1000,-50,50);
  _h_vert_y = new TH1F("Vertex y", ";[mm]", 1000,-50,50);
  _h_vert_z = new TH1F("Vertex z", ";[mm]", 1000,-50,50);
  _h_endpoint = new TH1F("Endpoint distance", "; [mm]", 1000,0,100);
  _h_end_x = new TH1F("Endpoint x", ";[mm]", 1000,-100,100);
  _h_end_y = new TH1F("Endpoint y", ";[mm]", 1000,-100,100);
  _h_end_z = new TH1F("Endpoint z", ";[mm]", 1000,-100,100);
  _h_decaylength = new TH1F("Decay Length", ";[mm]", 1000,0,100);
  _h_MID = new TH1F("Mother PDG", ";PDG ID", 1000, -500, 2500); //keep on for a check to make sure mother id vetting works

  _h_rmax = new TH1F("RMax", ";rMax",1000,-5000,5000);
  _h_rmin = new TH1F("RMin", ";rMin", 1000, -500, 500);
}

void GenHistProcessor::processRunHeader( LCRunHeader* /*run*/) {
}

void GenHistProcessor::processEvent( LCEvent * evt ) {
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.
  // Loop over MCParticles

  //get the collection of tracker hits
  //sim tracker hits --> get them from the digitized sample tbh
  /*
  LCCollection *trkhits = evt->getCollection("AllSTH");
  const int nSTH = trkhits->getNumberOfElements();
  double rmax = -1;
  double rmin = 10000;
  for(uint32_t i=0; i<nSTH; i++){
    const EVENT::SimTrackerHit *sth = static_cast<const EVENT::SimTrackerHit*>(trkhits->getElementAt(i));
    //get the MCP it's associated with
    const EVENT::MCParticle *mcp = static_cast<const EVENT::MCParticle*>(sth->getMCParticle());
    int pdg = mcp->getPDG();
    if(abs(pdg)!=211) continue; //studying only pions
    //now get the maximum and minimum radius in the event
    const double *pos = sth->getPosition();
    double rad = std::sqrt(std::pow(pos[0],2)+std::pow(pos[1],2)+std::pow(pos[2],2));
    if(rad > rmax) rmax = rad;
    if(rad < rmin) rmin = rad;
  }
  if(rmax>0)
    _h_rmax->Fill(rmax);
  if(rmin<10000)
    _h_rmin->Fill(rmin);
  */
  
  //map from pdg entry to actual pdg number for submitting in batch mode
  int pdgs[7]={0,11,13,15,211,111,2112};
  
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
    //double p =std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2)+std::pow(mom[2],2));
    double p = std::sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
    _h_p->Fill(p);
    _h_pt->Fill(pt);
    // std::cout<<p<<std::endl;
    
    double mass = mcp->getMass();
    //std::cout<<mass<<std::endl;
    //fill energy
    double E=mcp->getEnergy();
    _h_E->Fill(E);
    //std::cout<<E<<std::endl;

    double inv_E = std::sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]+mass*mass);
    _h_inv_E->Fill(inv_E-E);

    double inv_p = std::sqrt(E*E-mass*mass);
    _h_inv_p->Fill(inv_p-p);
    //std::cout<<inv_p<<std::endl;
    
    double invm = (std::pow(E,2)-std::pow(p,2));
    _h_inv_msqr->Fill(invm);
    

    //find and fill azimuth and polar angle

    double phi = std::atan(mom[1]/mom[0]);
    _h_phi->Fill(phi);

    double theta = std::acos(mom[2]/p);
    _h_theta->Fill(theta);

    //fill distance of vertex

    
    const double *vertex = mcp->getVertex();
    double vert=std::sqrt(std::pow(vertex[0],2)+std::pow(vertex[1],2)+std::pow(vertex[2],2));
    _h_vertex->Fill(vert);
    _h_vert_x->Fill(vertex[0]);
    _h_vert_y->Fill(vertex[1]);
    _h_vert_z->Fill(vertex[2]);

    //fill endpoint info
    const double *endpoint = mcp->getEndpoint();
    double end = std::sqrt(std::pow(endpoint[0],2)+std::pow(endpoint[1],2)+std::pow(endpoint[2],2));
    _h_endpoint->Fill(end);
    _h_end_x->Fill(endpoint[0]);
    _h_end_y->Fill(endpoint[1]);
    _h_end_z->Fill(endpoint[2]);
    
    double decaylength = std::sqrt(std::pow(vertex[0]-endpoint[0],2)+std::pow(vertex[1]-endpoint[1],2)+std::pow(vertex[2]-endpoint[2],2));

    _h_decaylength->Fill(decaylength);

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
