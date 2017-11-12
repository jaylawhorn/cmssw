#ifndef PulseShapeFunctor_h
#define PulseShapeFunctor_h 1

#include <typeinfo>

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"
#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"

#include <TH1F.h>
#include "Math/Functor.h"

namespace HcalConst{

   constexpr int maxSamples = 10;
   constexpr int maxPSshapeBin = 256;
   constexpr int nsPerBX = 25;
   constexpr float iniTimeShift = 92.5f;
   constexpr double invertnsPerBx = 0.04;
   constexpr int shiftTS = 4;

}

namespace FitterFuncs{
  
   class PulseShapeFunctor {
      public:
     PulseShapeFunctor(const HcalPulseShapes::Shape& pulse,bool iPedestalConstraint, bool iTimeConstraint,bool iAddPulseJitter,bool iAddTimeSlew,
		       double iPulseJitter,double iTimeMean,double iTimeSig,double iPedMean,double iPedSig,
		       double iNoise, unsigned int nSamplesToFit);
     ~PulseShapeFunctor();
     
     double EvalPulse(const double *pars, unsigned int nPar);
     
     void setDefaultcntNANinfit(){ cntNANinfit =0; }
     int getcntNANinfit(){ return cntNANinfit; }
     
     void setpsFitx(double *x ){ for(int i=0; i<HcalConst::maxSamples; ++i) psFit_x[i] = x[i]; }
     void setpsFity(double *y ){ for(int i=0; i<HcalConst::maxSamples; ++i) psFit_y[i] = y[i]; }
     void setpsFiterry (double *erry  ){ for(int i=0; i<HcalConst::maxSamples; ++i) psFit_erry  [i] = erry [i]; }
     void setpsFiterry2(double *erry2 ){ for(int i=0; i<HcalConst::maxSamples; ++i) psFit_erry2 [i] = erry2[i]; }
     void setpsFitslew (double *slew  ){ for(int i=0; i<HcalConst::maxSamples; ++i) {psFit_slew [i] = slew [i]; } }
     double sigmaHPDQIE8(double ifC);
     double sigmaSiPMQIE10(double ifC);
     double getSiPMDarkCurrent(double darkCurrent, double fcByPE, double lambda);
     void setinvertpedSig2(double x) { invertpedSig2_ = x; }

     double singlePulseShapeFunc( const double *x );
     double doublePulseShapeFunc( const double *x );
     double triplePulseShapeFunc( const double *x );

     double getPulseShape(int i) { 
       if (i>=0 && i<HcalConst::maxSamples)
	 return pulse_shape_[i]; 
       else
	 return 0;
     }
     
   private:
     std::array<float,HcalConst::maxPSshapeBin> pulse_hist;
     
     int cntNANinfit;
     std::vector<float> acc25nsVec, diff25nsItvlVec;
     std::vector<float> accVarLenIdxZEROVec, diffVarItvlIdxZEROVec;
     std::vector<float> accVarLenIdxMinusOneVec, diffVarItvlIdxMinusOneVec;
     void funcHPDShape(std::array<double,HcalConst::maxSamples> & ntmpbin, const double &pulseTime, const double &pulseHeight,const double &slew);
     double psFit_x[HcalConst::maxSamples], psFit_y[HcalConst::maxSamples], psFit_erry[HcalConst::maxSamples], psFit_erry2[HcalConst::maxSamples], psFit_slew[HcalConst::maxSamples];

     unsigned nSamplesToFit_;

     bool pedestalConstraint_;
     bool timeConstraint_;
     bool addPulseJitter_;
     bool unConstrainedFit_;
     double pulseJitter_;
     double timeMean_;
     double timeSig_;
     double pedMean_;
     double pedSig_;
     double noise_;
     double timeShift_;

     double inverttimeSig_, inverttimeSig2_;
     double invertpedSig_, invertpedSig2_;
     std::array<double,HcalConst::maxSamples> pulse_shape_;
     std::array<double,HcalConst::maxSamples> pulse_shape_sum_;

   };
   
}

#endif // PulseShapeFunctor_h
