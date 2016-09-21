#ifndef HcalTesting_Mahi_DoMahiAlgo_HH
#define HcalTesting_Mahi_DoMahiAlgo_HH

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "HcalTesting/Mahi/interface/EigenMatrixTypes.h"
#include "HcalTesting/Mahi/interface/PulseShapes.h"

/*namespace HcalConst{

  constexpr int maxSamples = 10;
  constexpr int maxPSshapeBin = 256;
  constexpr int nsPerBX = 25;
  constexpr float iniTimeShift = 92.5f;
  constexpr double invertnsPerBx = 0.04;

  }*/

class DoMahiAlgo
{
 public:

  typedef BXVector::Index Index;

  DoMahiAlgo();
  ~DoMahiAlgo() { };

  void Apply(const CaloSamples & cs, const std::vector<int> & capidvec, const HcalDetId detID, const HcalCalibrations & calibs, std::vector<double> & correctedOutput);
  bool DoFit(SampleVector amplitudes, std::vector<double> &correctedOutput);

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

  PulseShapes pulseShapeObj;

  double deltaT;

};

#endif
