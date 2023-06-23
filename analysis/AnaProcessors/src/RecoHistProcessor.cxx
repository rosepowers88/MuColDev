#include "AnaProcessors/RecoHistProcessor.hxx"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>

#include <marlin/AIDAProcessor.h>

#include <cmath>

RecoHistProcessor aRecoHistProcessor ;

RecoHistProcessor::RecoHistProcessor()
  : Processor("RecoHistProcessor") {
  // modify processor description
  _description = "RecoHistProcessor makes histograms at the reco level, analyzing clusters, tracks, and reconstructed particles." ;

  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("MinPt",
                             "Minimum particle pT.",
                             _minPt,
                             _minPt);

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

}


void RecoHistProcessor::book_pfo_histograms(){
  marlin::AIDAProcessor::histogramFactory(this);  
  _h_ptpf   = new TH1F("pT", ";pT [GeV]", 1000,0., 200);
  _h_ppf    = new TH1F("p", ";p [GeV]", 1000, 0, 200);
  _h_pdgpf  = new TH1F("PDG_ID", ";PDG ID", 10000, -500, 2500);
  _h_Npf    = new TH1F("nObjectsFound", ";N", 50, 0, 10);
  _h_Epf    = new TH1F("Energy", ";E [GeV]", 1000,0,200);

}
void RecoHistProcessor::book_cluster_histograms(){
  marlin::AIDAProcessor::histogramFactory(this);
  _h_Ecl = new TH1F("Cluster Energy", ";E [GeV]", 1000,0,200);
  _h_Xcl = new TH1F("X-position", ";x [mm]", 4000,-2000,2000);
  _h_Ycl = new TH1F("Y-position", ";y [mm]", 4000,-2000,2000);
  _h_Zcl = new TH1F("Z-position", ";z [mm]", 8000,-4000,4000);
  _h_iPhicl=new TH1F("Intrinsic Phi", ";#phi [rad]", 10,0,6.3);
  _h_iThetacl=new TH1F("Intrinsic Theta", ";#theta [rad]",10,0,3.2);
  _h_PIDcl = new TH1F("Most Likely PID", "; PID", 10000, -300, 2500);
  _h_Ncl    = new TH1F("nClustersFound", ";N", 50, 0, 10);

}
void RecoHistProcessor::book_track_histograms(){
  marlin::AIDAProcessor::histogramFactory(this);
  _h_D0trk = new TH1F("Impact Parameter (r-phi)", ";D0 [mm]", 100,0,200);
  _h_phitrk = new TH1F("Phi", ";#phi [rad]", 10,0,6.3);
  _h_omegatrk = new TH1F("Signed Curvature", ";#Omega [1/mm]", 1000,-200,200);
  _h_Z0trk = new TH1F("Impact Parameter (r-z)", ";Z0 [mm]", 100,0,200);
  _h_lamtrk =new TH1F("Dip Angle (r-z)", ";#lambda [rad]", 10,-3.2,3.2);
  _h_dEdxtrk=new TH1F("dEdx", ";dEdx [Gev/mm]",100,0,50);
  _h_Ntrk    = new TH1F("nTracksFound", ";N", 50, 0, 10);


}

void RecoHistProcessor::init() {
  // Print the initial parameters
  printParameters() ;

  //Basic histograms
  book_pfo_histograms();
  book_cluster_histograms();
  book_track_histograms();
}

void RecoHistProcessor::processRunHeader( LCRunHeader* /*run*/) {
}

void RecoHistProcessor::fill_pfo_histograms(LCCollection* inputCol){

  //get number of elements
  const int nEl = inputCol->getNumberOfElements();
  _h_Npf->Fill(nEl);

  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::ReconstructedParticle *pof=static_cast<const EVENT::ReconstructedParticle*>(inputCol->getElementAt(i));

    //get momentum and energy
    const double* mom = pof->getMomentum();

    double E=pof->getEnergy();
    double pt=std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    double p =std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2))+std::sqrt(std::pow(mom[2],2));
      
    _h_ppf->Fill(p);
    _h_ptpf->Fill(pt);
    _h_Epf->Fill(E);

    //finally, get pdg
    int pdg = pof->getType();
    _h_pdgpf->Fill(pdg);
  }

      
}
void RecoHistProcessor::fill_cluster_histograms(LCCollection* inputCol){

  //get number of elements
  const int nEl = inputCol->getNumberOfElements();
  _h_Ncl->Fill(nEl);

  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::Cluster *cl=static_cast<const EVENT::Cluster*>(inputCol->getElementAt(i));
    //get energy and position
    double E = cl->getEnergy();
    _h_Ecl->Fill(E);
    const float* pos = cl->getPosition();
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    _h_Xcl->Fill(x);
    _h_Ycl->Fill(y);
    _h_Zcl->Fill(z);

    //Get phi and theta
    double phi = cl->getIPhi();
    _h_iPhicl->Fill(phi);
    double theta = cl->getITheta();
    _h_iThetacl->Fill(theta);

    //Get most likely particle ID
    const ParticleIDVec& PID = cl->getParticleIDs();
    if(PID.size()!=0){
    double pid_pdg = PID[0]->getPDG();
    _h_PIDcl->Fill(pid_pdg);
    }



  }
}
void RecoHistProcessor::fill_track_histograms(LCCollection* inputCol){

  //get number of elements
  const int nEl = inputCol->getNumberOfElements();
  _h_Ntrk->Fill(nEl);

  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::Track *trk=static_cast<const EVENT::Track*>(inputCol->getElementAt(i));

    //get D0
    double D0 = trk->getD0();
    _h_D0trk->Fill(D0);

    //get phi
    double phi = trk->getPhi();
    _h_phitrk->Fill(phi);

    //get curvature
    double omega = trk->getOmega();
    _h_omegatrk->Fill(omega);

    //get Z0
    double Z0 = trk->getZ0();
    _h_Z0trk->Fill(Z0);

    //get dip angle
    double tanlam = trk->getTanLambda();
    double lambda = std::atan(tanlam);
    _h_lamtrk->Fill(lambda);

    //get dE/dx
    double dEdx = trk->getdEdx();
    _h_dEdxtrk->Fill(dEdx);
  }
  
}


void RecoHistProcessor::processEvent( LCEvent * evt ) {
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.
  std::vector<std::string> inputnames {_inputCollectionNameR, _inputCollectionNameC, _inputCollectionNameT};
  for(int i=0; i<inputnames.size(); i++){
  
  LCCollection* inputCol = evt->getCollection(inputnames[i]);
  
  if( inputCol->getTypeName() == lcio::LCIO::RECONSTRUCTEDPARTICLE ) {
    //fill the histograms for PFOs
    fill_pfo_histograms(inputCol);
  }
  else if( inputCol->getTypeName() == lcio::LCIO::CLUSTER) {
    //fill the histograms for clusters
    fill_cluster_histograms(inputCol);
  }
  else if( inputCol->getTypeName() == lcio::LCIO::TRACK) {
    //fill the histograms for tracks
    fill_track_histograms(inputCol);
  }
  else{
    throw EVENT::Exception( "Invalid collection type: " + inputCol->getTypeName() ) ;
  }
  }
}

void RecoHistProcessor::check( LCEvent * /*evt*/ )
{ }

void RecoHistProcessor::end()
{ }
