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
  _description = "RecoHistProcessor makes histograms at the reco level." ;

  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("MinPt",
                             "Minimum particle pT.",
                             _minPt,
                             _minPt);

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
  
  /*

    registerInputCollection( LCIO::SIMTRACKERHIT,
    "InputCollectionName",
    "Name of the input collection",
    _inputCollectionName,
    _inputCollectionName
    );
    registerOutputCollection(LCIO::SIMTRACKERHIT,
    "OutputCollectionName",
    "Name of the output collection",
    _outputCollectionName,
    _outputCollectionName
    );
  */
}

void RecoHistProcessor::init() {
  // Print the initial parameters
  printParameters() ;

  //Basic histograms
  marlin::AIDAProcessor::histogramFactory(this);
  _h_pt   = new TH1F("pT", ";pT [GeV]", 100,0., 20);
  _h_p    = new TH1F("p", ";p [GeV]", 100, 0, 20);
  _h_pdg  = new TH1F("PDG_ID", ";PDG ID", 10000, -500, 2500);
  _h_N    = new TH1F("nObjectsFound", ";N", 50, 0, 10);
  _h_E    = new TH1F("Energy", ";E [GeV]", 100,0,20);

  // _h_phi  = new TH1F("Phi", ";#phi [rad]", 50,-3.2,3.2);
  //_h_tanLam = new TH1F("tanLambda", ";tan(#lambda)", 5,-1,1);
  
  //restricted histograms
  _h_p_pi = new TH1F("p_{#pi}", ";p [GeV]",100,0,20);
  _h_pt_pi = new TH1F("pT_{#pi}", ";pT [GeV]", 100,0,20);
  _h_pdg_pi = new TH1F("PDG_ID_{#pi}", ";PDG_ID", 400, -220, 220);
  _h_N_pi  = new TH1F("nObjectsFound_{#pi}", ";N", 50,0,10);
  _h_E_pi    = new TH1F("E_{#pi}", ";E [GeV]", 100,0,20);
  _h_x_pi    = new TH1F("x_{#pi}", ";x", 100,-1000,1000);
  _h_y_pi    = new TH1F("y_{#pi}", ";y", 100,-1000,1000);
  _h_z_pi    = new TH1F("z_{#pi}", ";z", 100,-1000,1000);
  _h_r_pi    = new TH1F("r_{#pi}", ";z", 100,0,4000);
  
  
  //topology histograms
  //_h_theta = new TH1F("#theta of Cone", ";#theta [rad]", 50,0,3.142);

  


}

void RecoHistProcessor::processRunHeader( LCRunHeader* /*run*/) {
}

void RecoHistProcessor::processEvent( LCEvent * evt ) {
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.
  // Loop over MCParticles
  LCCollection* inputCol = evt->getCollection(_inputCollectionName);


  if( inputCol->getTypeName() != lcio::LCIO::RECONSTRUCTEDPARTICLE ) {
    throw EVENT::Exception( "Invalid collection type: " + inputCol->getTypeName() ) ;
  }



  const int nEl = inputCol->getNumberOfElements();
  if(nEl>0){
    std::cout<<nEl<<std::endl;
  }
  int nPi = 0;

  //std::cout << "Found " << nEl  << " particles." << std::endl;
  if(nEl==3){
    const int evtnum = evt->getEventNumber();
    //uncomment to get event numbers of 3 prong decays
    //std::cout<<"Three prong: "<<evtnum << std::endl;
  }
  
  _h_N->Fill(nEl);

  

  std::vector <double> xs;
  std::vector <double> ys;
  std::vector <double> zs;
  const double* mom;
  const double* mom1;
  const double* mom2;
  const double* mom3;
  double E;
  

  for(uint32_t i=0;i<nEl;i++) {
    const EVENT::ReconstructedParticle *mcp=static_cast<const EVENT::ReconstructedParticle*>(inputCol->getElementAt(i));
    //double phi = mcp->getPhi();
    //_h_phi->Fill(phi);
    //double tanlam = mcp->getTanLambda();
    //_h_tanLam->Fill(tanlam);
    mom = mcp->getMomentum();
    //Fill particle pT and p
    /*  if(i==0){
	mom1=mom;
      }
      else if(i==1){
	mom2=mom;
      }
      else{
	mom3=mom;
	} */
    E=mcp->getEnergy();
    double pt=std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
    double p =std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2))+std::sqrt(std::pow(mom[2],2));

    // if(pt<_minPt) {
    //continue; 
    //}
      
    _h_p->Fill(p);
    _h_pt->Fill(pt);
    _h_E->Fill(E);
    /*

    // get position, calculate cone topology
    if(nEl>=3){
    const float* pos = mcp-> getReferencePoint();
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];
    double r = std::sqrt(std::pow(x,2)+std::pow(y,2)+std::pow(z,2));
    _h_r_pi->Fill(r);
    _h_x_pi->Fill(x);
    _h_y_pi->Fill(y);
    _h_z_pi->Fill(z);
    //std::cout<<"x,y,z: "<<x<<", "<<y<<", "<<z<<std::endl; just a check to confirm working
    xs.push_back(x);
    ys.push_back(y);
    zs.push_back(z);
    }
    */

    //get Type
    int type = mcp -> getType();
    _h_pdg->Fill(type);
    //std::cout<<type<<std::endl;

    if(abs(type)==211){
    nPi++;
    }
    if(nPi>=3&&i==nEl-1){
    //std::cout<<evt->getEventNumber()<<": Cone candidate"<<std::endl;
    }

    if(nEl == 3 ){
    if(type==2112){
    int evtnum = evt->getEventNumber();
    //std::cout<<"neutron: "<<evtnum<<std::endl;
    }
    }

    //Examine only pions
    //if(abs(type)==211)&& nEl==3){
    if(abs(type)==211){
      _h_p_pi->Fill(p);
      _h_pt_pi->Fill(pt);
      _h_pdg_pi->Fill(type);
      _h_N_pi->Fill(nPi);
      _h_E->Fill(E);
    }
      

    /*}
      if(nEl==3){
      //define cross product
      auto norm=[&](std::vector<double> vec){
      std::vector<double> vecnorm;
      for(int i=0; i<3; i++){
      double coord = vec[i]/std::sqrt(std::pow(vec[0],2)+std::pow(vec[1],2)+std::pow(vec[2],2));
      vecnorm.push_back(coord);
      }
      return(vecnorm);
      };
      auto crossnorm=[&](std::vector<double> v1, std::vector<double> v2){
      std::vector<double> product;
      product.push_back( v1[1] * v2[2] - v1[2] * v2[1]);
      product.push_back( -(v1[0] * v2[2] - v1[2] * v2[0]));
      product.push_back(v1[0] * v2[1] - v1[1] * v2[0]);
      product=norm(product);
      return(product);
      };

      auto mid=[&](std::vector<double> pt1, std::vector<double> pt2){
      std::vector<double> midpoint;
      if(pt2[0]>pt1[0]){
      midpoint.push_back((pt2[0]-pt1[0])/2+pt1[0]);
      }
      else{
      midpoint.push_back((pt1[0]-pt2[0])/2+pt2[0]);
      }
      if(pt2[1]>pt1[1]){
      midpoint.push_back((pt2[1]-pt1[1])/2+pt1[1]);
      }
      else{
      midpoint.push_back((pt1[1]-pt2[1])/2+pt2[1]);
      }
      if(pt2[2]>pt1[2]){
      midpoint.push_back((pt2[2]-pt1[2])/2+pt1[2]);
      }
      else{
      midpoint.push_back((pt1[2]-pt2[2])/2+pt2[2]);
      }
      return(midpoint);
      };
      double l1 = std::sqrt(std::pow((xs[1]-xs[0]),2)+std::pow((ys[1]-ys[0]),2)+std::pow((zs[1]-zs[0]),2));
      double l2 = std::sqrt(std::pow((xs[2]-xs[0]),2)+std::pow((ys[2]-ys[0]),2)+std::pow((zs[2]-zs[0]),2));
      double l3 = std::sqrt(std::pow((xs[2]-xs[1]),2)+std::pow((ys[2]-ys[1]),2)+std::pow((zs[2]-zs[1]),2));
      //get the angle opposite l3 with law of cosines:
      double alpha = std::acos((std::pow(l1,2)+std::pow(l2,2)-std::pow(l3,2))/(2*l1*l2));
      //get the radius with law of sines:
      double R = l3/(2*sin(alpha));
      //now find the plane equation:
      std::vector<double> AB {xs[2]-xs[1], ys[2]-ys[1], zs[2]-zs[1]};
      std::vector<double> AC {xs[2]-xs[0], ys[2]-ys[0], zs[2]-zs[0]};
      std::vector<double> BC {xs[1]-xs[0], ys[1]-ys[0], zs[1]-zs[0]};
      std::vector<double> normal = crossnorm(AB, AC);
      //std::vector<double> normal;
      /* for(int i=0; i<3; i++){
      double coord = big_normal[i]/std::sqrt(std::pow(big_normal[0],2)+std::pow(big_normal[1],2)+std::pow(big_normal[2],2));
      normal.push_back(coord);
      }*/
    /*   std::vector<double> pt_1 {xs[0], ys[0], zs[0]}; std::vector<double> pt_2 {xs[1], ys[1], zs[1]}; std::vector<double> pt_3 {xs[2], ys[2], zs[2]};
	 std::vector<double> mid_l1 = mid(pt_1, pt_2);
	 std::vector<double> mid_l2 = mid(pt_3, pt_1);
	 std::vector<double> mid_l3 = mid(pt_3, pt_2);
	 std::vector<double> bisect_l1 = crossnorm(BC,normal);
	 std::vector<double> bisect_l2 = crossnorm(AC,normal);

	 std::vector<double> circumcenter;
	 for(int i=0; i<3; i++){
	 circumcenter.push_back((mid_l2[i]-mid_l1[i])/(bisect_l1[i]-bisect_l2[i]));
	 }

	 double p =std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2))+std::sqrt(std::pow(mom[2],2));
	 double beta = p/E;
	 std::vector<double> vel1;
	 std::vector<double> vel2;
	 std::vector<double> vel3;
	 for(int i=0; i<3; i++){
	 vel1.push_back(mom1[i]);
	 vel2.push_back(mom2[i]);
	 vel3.push_back(mom3[i]);
	 }
	 /*
	 std::vector<double> normvel1=norm(vel1);
	 std::vector<double> normvel2=norm(vel2);
	 std::vector<double> normvel3=norm(vel3);
	 std::cout<<normvel1[0]<<" : Made it past the normalization"<<std::endl;
    */
    /* vel1=norm(vel1);
       vel2=norm(vel2);
       vel3=norm(vel3);
       // for(int i=0; i<3; i++){
       //vel[i]=vel[i]*beta;
       //}

       //maybe we should leave it normed actually
       auto dot=[&](std::vector<double> vec1, std::vector<double> vec2){
       double product=0;
       for(int i=0; i<3; i++){
       double prod=vec1[i]*vec2[i];
       product+=prod;
       }
       return(product);
       };
       double cosThetaAvg = (dot(vel1,normal)+dot(vel2,normal)+dot(vel3,normal))/3;
       //std::cout<<cosThetaAvg<<std::endl;
       double theta = std::acos(cosThetaAvg);
       //std::cout<<theta<<std::endl;

       _h_theta->Fill(theta);
    
    
    
    
    
       //std::cout<<"circumcenter check: "<<circumcenter[0]<<" , "<<circumcenter[1]<<" , "<<circumcenter[2]<<std::endl;

       }*/
  }
}

void RecoHistProcessor::check( LCEvent * /*evt*/ )
{ }

void RecoHistProcessor::end()
{ }
