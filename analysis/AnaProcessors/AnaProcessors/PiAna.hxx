#pragma once

#include <TH1F.h>

#include <marlin/Processor.h>

class PiAna : public marlin::Processor
{
public:
  virtual Processor * newProcessor() {return new PiAna;}

  PiAna(const PiAna &) = delete;
  PiAna& operator = (const PiAna &) = delete;
  PiAna() ;

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
  TH1 *_h_E_neutron = nullptr;
  TH1 *_h_E_PFneutron = nullptr;
  TH1 *_h_E_pi0 = nullptr;

  
};
