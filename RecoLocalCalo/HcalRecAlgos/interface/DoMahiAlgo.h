#ifndef RecoLocalCalo_HcalRecAlgos_DoMahiAlgo_HH
#define RecoLocalCalo_HcalRecAlgos_DoMahiAlgo_HH

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/EigenMatrixTypes.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"

#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/PulseShapeFitOOTPileupCorrection.h"

class DoMahiAlgo
{
 public:

  typedef BXVector::Index Index;
  
  DoMahiAlgo();
  ~DoMahiAlgo() { };

  void phase1Apply(const HBHEChannelInfo& channelData, float& reconstructedEnergy, float& chi2);
  
  bool DoFit(SampleVector amplitudes, std::vector<float> &correctedOutput);

  const HcalPulseShapes::Shape* currentPulseShape_=nullptr;

  void setPulseShapeTemplate  (const HcalPulseShapes::Shape& ps);
  void resetPulseShapeTemplate(const HcalPulseShapes::Shape& ps);

 private:

  //for pulse shapes
  int cntsetPulseShape;

  std::unique_ptr<FitterFuncs::PulseShapeFunctor> psfPtr_;
  std::unique_ptr<ROOT::Math::Functor> pfunctor_;

  bool Minimize();
  bool UpdateCov();
  double CalculateChiSq();
  bool NNLS();

  SampleVector _amplitudes;
  SampleMatrix _invCovMat;

  SampleVector _pedWidth;
  
  FullSampleMatrix noiseCor;
  FullSampleMatrix pulseCov;

  SamplePulseMatrix _pulseMat;

  PulseVector _ampVec;
  PulseVector _ampVecMin;
  PulseVector _errVec;
  FullSampleVector pulseShape;
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
