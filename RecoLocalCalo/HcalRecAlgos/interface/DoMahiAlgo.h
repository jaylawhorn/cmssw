#ifndef RecoLocalCalo_HcalRecAlgos_DoMahiAlgo_HH
#define RecoLocalCalo_HcalRecAlgos_DoMahiAlgo_HH

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/EigenMatrixTypes.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"

#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/PulseShapeFitOOTPileupCorrection.h"

class DoMahiAlgo
{
 public:
  DoMahiAlgo();
  ~DoMahiAlgo() { };

  void phase1Apply(const HBHEChannelInfo& channelData, float& reconstructedEnergy, float& reconstructedTime, float& chi2);  
  bool DoFit(SampleVector amplitudes, std::vector<float> &correctedOutput, int nbx);

  void setParameters(double iTS4Thresh, bool iApplyTimeSlew, HcalTimeSlew::BiasSetting slewFlavor,
		     double iMeanTime, double iTimeSigmaHPD, double iTimeSigmaSiPM,
		     const std::vector <int> &iActiveBXs, int iNMaxItersMin, int iNMaxItersNNLS,
		     double iDeltaChiSqThresh, double iNnlsThresh);

  void setPulseShapeTemplate  (const HcalPulseShapes::Shape& ps);
  void resetPulseShapeTemplate(const HcalPulseShapes::Shape& ps);

  void setDebug(int val);

  typedef BXVector::Index Index;
  const HcalPulseShapes::Shape* currentPulseShape_=nullptr;

 private:

  // Python-configurables
  float TS4Thresh_; //0

  bool applyTimeSlew_; //true
  HcalTimeSlew::BiasSetting slewFlavor_; //medium

  float meanTime_; // 0
  float timeSigmaHPD_; // 5.0
  float timeSigmaSiPM_; //2.5

  std::vector <int> activeBXs_;

  int nMaxItersMin_; //500 
  int nMaxItersNNLS_; //500

  float deltaChiSqThresh_; //1e-3
  float nnlsThresh_; //1e-11

  int niterTot_;

  int doDebug;
  bool isHPD;
  //HcalTimeSlew::BiasSetting slewFlavor_;

  //for pulse shapes
  int cntsetPulseShape;

  //HcalDetId _detID;
  unsigned int _nPulseTot;

  //holds active bunch crossings
  BXVector _bxs;  

  BXVector _bxsMin;
  unsigned int _nP;
  double _chiSq;

  std::unique_ptr<FitterFuncs::PulseShapeFunctor> psfPtr_;
  std::unique_ptr<ROOT::Math::Functor> pfunctor_;

  bool Minimize();
  bool UpdateCov();
  bool UpdatePulseShape(double itQ, FullSampleVector &pulseShape, FullSampleMatrix &pulseCov);
  double CalculateChiSq();
  bool NNLS();

  //holds data samples
  SampleVector _amplitudes;
  //holds corrections per pulse for TS4->whole pulse charge
  SampleVector _fullPulseCorrection;
  //holds inverse covariance matrix
  SampleMatrix _invCovMat;

  //holds diagonal noise terms
  SampleVector _noiseTerms;
  //holds constant pedestal constraint
  double _pedConstraint;
  
  //holds full covariance matrix for a pulse shape 
  //varied in time
  FullSampleMatrix pulseCov;
  FullSampleMatrix pulseCovOOTM;
  FullSampleMatrix pulseCovOOTP;

  //holds full pulse shape template
  FullSampleVector pulseShape;
  FullSampleVector pulseShapeOOTM;
  FullSampleVector pulseShapeOOTP;

  //holds matrix of pulse shape templates for each BX
  SamplePulseMatrix _pulseMat;

  //for FNNLS algorithm
  PulseVector _ampVec;
  PulseVector _ampVecMin;
  PulseVector _errVec;
  PulseVector ampvecpermtest;

  SamplePulseMatrix invcovp;
  PulseMatrix aTaMat; // A-transpose A (matrix)
  PulseVector aTbVec; // A-transpose b (vector)
  PulseVector wVec; // w (vector)
  PulseVector updateWork; // w (vector)

  SampleDecompLLT _covDecomp;
  SampleMatrix _covDecompLinv;
  PulseMatrix _topleft_work;
  PulseDecompLDLT _pulseDecomp;

}; 
#endif
