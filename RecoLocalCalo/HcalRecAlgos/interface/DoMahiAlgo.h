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
#include "RecoLocalCalo/HcalRecAlgos/interface/PulseShapeFunctor.h"

#include <unordered_map>

class DoMahiAlgo
{
 public:
  DoMahiAlgo();
  ~DoMahiAlgo() { };

  void phase1Apply(const HBHEChannelInfo& channelData, float& reconstructedEnergy, float& chi2);  
  bool DoFit(SampleVector amplitudes, std::vector<float> &correctedOutput);

  void setParameters(bool iApplyTimeSlew, HcalTimeSlew::BiasSetting slewFlavor,
		     double iMeanTime, double iTimeSigmaHPD, double iTimeSigmaSiPM, 
		     const std::vector <int> &iActiveBXs, int iNMaxIters);

  void setPulseShapeTemplate  (const HcalPulseShapes::Shape& ps);
  void resetPulseShapeTemplate(const HcalPulseShapes::Shape& ps);

  typedef BXVector::Index Index;
  const HcalPulseShapes::Shape* currentPulseShape_=nullptr;

 private:

  bool applyTimeSlew_;

  double fcByPe_;

  double meanTime_;
  double timeSigmaHPD_;
  double timeSigmaSiPM_;

  std::vector <int> activeBXs_;
  int nMaxIters_;
  int nMaxItersNNLS_;

  double dt_;
  int pulseOffset_;

  HcalTimeSlew::BiasSetting slewFlavor_;

  //for pulse shapes
  int cntsetPulseShape;

  unsigned int nPulseTot_;

  //holds active bunch crossings
  BXVector bxs_;

  BXVector bxsMin_;
  unsigned int nP_;
  double chiSq_;

  std::unique_ptr<FitterFuncs::PulseShapeFunctor> psfPtr_;
  std::unique_ptr<ROOT::Math::Functor> pfunctor_;

  bool Minimize();
  bool UpdateCov();
  bool UpdatePulseShape(double itQ, FullSampleVector &pulseShape, FullSampleMatrix &pulseCov);
  double CalculateChiSq();
  bool NNLS();

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
  
  std::unordered_map<int, int> mapBXs;
  //holds full covariance matrix for a pulse shape 
  //varied in time
  std::vector<FullSampleMatrix> pulseCovArray_;

  //holds full pulse shape template
  std::vector<FullSampleVector> pulseShapeArray_;

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
