#ifndef RecoLocalCalo_HcalRecAlgos_DoMahiAlgo_HH
#define RecoLocalCalo_HcalRecAlgos_DoMahiAlgo_HH

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/EigenMatrixTypes.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"


class DoMahiAlgo
{
 public:

  typedef BXVector::Index Index;
  
  DoMahiAlgo();
  ~DoMahiAlgo() { };

  void phase1Apply(const HBHEChannelInfo& channelData);//, float& reconstructedEnergy, float& chi2);
  
  bool DoFit(SampleVector amplitudes, SampleVector gains, std::vector<float> &correctedOutput);

  //void setPulseShapeTemplate();

 private:

  bool Minimize();
  bool UpdateCov();
  double CalculateChiSq();
  bool NNLS();
  SampleVector _amplitudes;
  SampleMatrix _invCovMat;
  
  FullSampleMatrix noiseCor;
  FullSampleMatrix pulseCov;

  SamplePulseMatrix _pulseMat;

  PulseVector _ampVec;
  PulseVector _ampVecMin;
  PulseVector _errVec;
  SampleVector pulseShape;
  SampleVector pulseShapeM;
  SampleVector pulseShapeP;
  SampleVector zeroShape;

  PulseVector ampvecpermtest;

  HcalDetId _detID;

  BXVector _bxs;
  BXVector _bxsMin;
  unsigned int _nPulseTot;
  unsigned int _nP;  

  double _chiSq;

  SamplePulseMatrix invcovp;
  PulseMatrix aTaMat; // A-transpose A (matrix)
  PulseVector aTbVec; // A-transpose b (vector)
  PulseVector wVec; // w (vector)
  PulseVector updateWork; // w (vector)

  SampleDecompLLT _covDecomp;
  SampleMatrix _covDecompLinv;
  PulseMatrix _topleft_work;
  PulseDecompLDLT _pulseDecomp;

  //PulseShapes pulseShapeObj;

}; 
#endif
