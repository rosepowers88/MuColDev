#pragma once

#include <TH1F.h>

#include <marlin/Processor.h>

class CalAna : public marlin::Processor
{
public:
  virtual Processor * newProcessor() {return new CalAna;}

  CalAna(const CalAna &) = delete;
  CalAna& operator = (const CalAna &) = delete;
  CalAna() ;

  virtual void init();
  virtual void processRunHeader(LCRunHeader * run);
  virtual void processEvent( LCEvent * evt);
  virtual void check( LCEvent * evt );
  virtual void end();

private:
  std::string _inputcolPFO {};
  std::string _inputcolMC {};

  float _minPt {};

  // TH1 *_h_pion_misID = nullptr;
   TH1 *_h_nocluster = nullptr;
  /*
  TH1 *_h_deltaBH = nullptr;
  TH1* _h_deltaEH = nullptr;
  TH1 *_h_badBarrelContributers = nullptr;
  TH1 *_h_noClpT = nullptr;
  */
  TH1 *_h_pass_CL1 = nullptr;
  TH1 *_h_all_CL1 = nullptr;
  TH1 *_h_pass_CL3 = nullptr;
  TH1 *_h_all_CL3 = nullptr;
  TH1 *_h_nCL = nullptr;
  TH1 *_h_phi_CL = nullptr;
  TH1 *_h_theta_CL = nullptr;
  TH1 *_h_nUnMatchedCL = nullptr;

  TH1 *_h_pass_pdg = nullptr;
  TH1 *_h_all_pdg = nullptr;
  TH1 *_h_pass_phi = nullptr;
  TH1 *_h_all_phi = nullptr;
  TH1 *_h_pass_theta = nullptr;
  TH1 *_h_all_theta = nullptr;
  TH1 *_h_pass_E = nullptr;
  TH1 *_h_all_E = nullptr;

  TH1 *_h_pass_trks3 = nullptr;
  TH1 *_h_all_trks3 = nullptr;
  TH1 *_h_pass_trks1 = nullptr;
  TH1 *_h_all_trks1 = nullptr;

  TH1 *_h_delta_E = nullptr;
  
};
