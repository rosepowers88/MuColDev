#pragma once

#include <TH1F.h>

#include <marlin/Processor.h>

class AnaPFOs : public marlin::Processor
{
public:
  virtual Processor * newProcessor() {return new AnaPFOs;}

  AnaPFOs(const AnaPFOs &) = delete;
  AnaPFOs& operator = (const AnaPFOs &) = delete;
  AnaPFOs() ;

  virtual void init();
  virtual void processRunHeader(LCRunHeader * run);
  virtual void processEvent( LCEvent * evt);
  virtual void check( LCEvent * evt );
  virtual void end();

private:
  std::string _inputcolPFO {};
  std::string _inputcolMC {};

  float _minPt {};

  TH1 *_h_pion_misID = nullptr;
};
