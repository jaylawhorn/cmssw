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

  bool Minimize();
  bool UpdateCov();
  bool UpdatePulseShape(double itQ, FullSampleVector &pulseShape, FullSampleMatrix &pulseCov);
  double CalculateChiSq();
  bool NNLS();

  //hard coded in initializer
  const unsigned int FullTSSize_;
  const unsigned int FullTSofInterest_;

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

  unsigned int BXSize_;
  int BXOffset_;

  //from channelData
  float dt_;
  float darkCurrent_;
  float fcByPe_;

  unsigned int TSSize_;
  unsigned int TSOffset_;

  unsigned int FullTSOffset_;

  int niterTot_;

  int doDebug;
  //bool isHPD;

  //for pulse shapes
  int cntsetPulseShape_;

  //HcalDetId _detID;
  unsigned int nPulseTot_;

  //holds active bunch crossings
  BXVector bxs_;  
  BXVector bxsMin_;
  unsigned int nP_;
  double chiSq_;

  std::unique_ptr<FitterFuncs::PulseShapeFunctor> psfPtr_;
  std::unique_ptr<ROOT::Math::Functor> pfunctor_;

  //holds data samples
  SampleVector amplitudes_;
  //holds corrections per pulse for TS4->whole pulse charge
  SampleVector fullPulseCorrection_;
  //holds inverse covariance matrix
  SampleMatrix invCovMat_;

  //holds diagonal noise terms
  SampleVector noiseTerms_;
  //holds constant pedestal constraint
  double pedConstraint_;
  
  //holds full covariance matrix for a pulse shape 
  //varied in time
  FullSampleMatrix pulseCov_;
  FullSampleMatrix pulseCovOOTM_;
  FullSampleMatrix pulseCovOOTP_;

  //holds full pulse shape template
  FullSampleVector pulseShape_;
  FullSampleVector pulseShapeOOTM_;
  FullSampleVector pulseShapeOOTP_;

  //holds matrix of pulse shape templates for each BX
  SamplePulseMatrix pulseMat_;

  //for FNNLS algorithm
  PulseVector ampVec_;
  PulseVector ampVecMin_;
  PulseVector errVec_;
  PulseVector ampvecpermtest_;

  SamplePulseMatrix invcovp_;
  PulseMatrix aTaMat_; // A-transpose A (matrix)
  PulseVector aTbVec_; // A-transpose b (vector)
  PulseVector wVec_; // w (vector)
  PulseVector updateWork_; // w (vector)

  SampleDecompLLT covDecomp_;
  SampleMatrix covDecompLinv_;
  PulseMatrix topleft_work_;
  PulseDecompLDLT pulseDecomp_;

}; 
#endif
