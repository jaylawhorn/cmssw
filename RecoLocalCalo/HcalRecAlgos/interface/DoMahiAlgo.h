#ifndef RecoLocalCalo_HcalRecAlgos_DoMahiAlgo_HH
#define RecoLocalCalo_HcalRecAlgos_DoMahiAlgo_HH

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/EigenMatrixTypes.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/PulseShapes.h"

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

  // FIXME: need to convert the correctedOutput to memory free
  void Apply(const CaloSamples & cs, const std::vector<int> & capidvec, const HcalDetId & detID, const HcalCalibrations & calibs, std::vector<double> & correctedOutput);
  
  bool DoFit(SampleVector amplitudes, SampleVector gains, std::vector<double> &correctedOutput);

  void setPulseShapeTemplate(bool useCSV, std::string filename);

 private:
  
  bool Minimize();
  bool UpdateCov();
  double CalculateChiSq();
  bool NNLS();

  void getPulseShape(float q, HcalDetId detID, float t, SampleVector &pulseShape, float sigma);

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

  bool _useCSV;

  float minCharge_[58];
  float maxCharge_[58];
  float pulseFrac_[58][10];
  float pulseFracDeriv_[58][10];

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
