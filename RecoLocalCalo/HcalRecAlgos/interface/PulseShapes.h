#ifndef HcalTesting_Mahi_PulseShapes_HH
#define HcalTesting_Mahi_PulseShapes_HH

#include "TMath.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/EigenMatrixTypes.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include <iostream>


class PulseShapes
{
 public:

  PulseShapes() { };
  ~PulseShapes() { };

 private:

  double landauWidth(double q, HcalDetId detID, double sigma=0) {
    double qt = (q<550 ? q : 550);
    if (detID.subdet()==HcalSubdetector::HcalBarrel) return exp( -2.69212 + 0.339142*sigma + (-0.0200032+0.00775682*sigma) * qt ) + 0.166766 + 0.00137787*sigma;
    else if (detID.subdet()==HcalSubdetector::HcalEndcap) return exp( -1.25575 + 0.0782487*sigma + (-0.033789+0.00217511*sigma) * qt ) + 0.138741 + 0.000636494*sigma;
    else {
      std::cout << "unknown subdetector!" << std::endl;
      return 0;
    }
  };

  double mostProb(double q, HcalDetId detID, double sigma=0) {
    double qt = (q<550 ? q : 550);
    if (detID.subdet()==HcalSubdetector::HcalBarrel) return 4.13324 + 0.00465576*sigma + (-0.000471641+ 3.44855e-05*sigma) * qt + (3.69809e-07 +5.40687e-08*sigma) * qt*qt;
    else if (detID.subdet()==HcalSubdetector::HcalEndcap) {
      double tmp=3.98878+0.00707374*sigma;
      tmp+=(0.00274516+0.000264674*sigma)*qt;
      tmp+=(-2.29749e-05+3.33247e-06*sigma)*qt*qt;
      tmp+=(9.044e-08+1.90669e-08*sigma)*qt*qt*qt;
      tmp+=(-1.92282e-10+5.43992e-11*sigma)*qt*qt*qt*qt;
      tmp+=(2.13577e-13+7.5391e-14*sigma)*qt*qt*qt*qt*qt;
      tmp+=(-9.72045e-17+4.04107e-17*sigma)*qt*qt*qt*qt*qt*qt;
      return tmp;
    }
    else {
      std::cout << "unknown subdetector!" << std::endl;
      return 0;
    }
  };

  double gausSigma(double q, HcalDetId detID, double sigma=0) {
    double qt = (q<550 ? q : 550);
    if (detID.subdet()==HcalSubdetector::HcalBarrel) return exp( -1.41872 + 0.573027*sigma + (-0.0558063 + 0.0186813*sigma) * qt ) + 0.315172 + 0.00127181*sigma;
    else if (detID.subdet()==HcalSubdetector::HcalEndcap) return exp( -1.9489 + 0.0675424*sigma + (-0.026614 +0.00164744*sigma) * qt ) + 0.316884 + 0.000473545*sigma;
    else {
      std::cout << "unknown subdetector!" << std::endl;
      return 0;
    }
  };
  double asym(double q, HcalDetId detID, double sigma=0) {
    double qt = (q<550 ? q : 550);
    if (detID.subdet()==HcalSubdetector::HcalBarrel) return exp( -1.85949 + 0.0476665*sigma + (-0.00946504+0.00072732*sigma) * qt ) + -0.0595725 + 0.00189319*sigma;
    else if (detID.subdet()==HcalSubdetector::HcalEndcap) return exp( -1.533 +0.0303588*sigma + (-0.00696893+0.000480861*sigma) * qt ) + -0.101749 + 0.00305867*sigma;
    else {
      std::cout << "unknown subdetector!" << std::endl;
      return 0;
    }
  };

 public:
  double integral(double q, HcalDetId detID, double t, double sigma=0) {
    SampleVector Shape;
    computeLAGShape(q, detID, t, Shape, sigma);

    double tmp=0;
    for (int i=0; i<10; i++) {
      tmp+=Shape.coeff(i);
    }
    return q*tmp;

  }

  void computeLAGShape(double q, HcalDetId detID, double t, SampleVector &Shape, double sigma=0) {
    double lwidth = landauWidth(q, detID, sigma);
    double MProb = mostProb(q, detID, sigma) + (t/25.0);
    double gsigma = gausSigma(q, detID, sigma);
    double asymmcoef = asym(q, detID, sigma);

    double timeslice[10] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0};
    double funcvalue[10];
    for (int m = 0; m<10; m++){ funcvalue[m] = 0.0;}

    double invsq2pi = 0.3989422804014;
    double mpshift  = -0.22278298;

    double np = 100.0;      // number of convolution steps
    double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
    
    for (int j=0; j<10 ;j++){
      double xx;
      double yy;
      double mpc = MProb - mpshift * lwidth;
      double fland;
      double summ = 0.0;
      double xlow = timeslice[j] - sc * gsigma;
      double xupp = timeslice[j] + sc * gsigma;
      double step = (xupp-xlow) / np;
      double asigma;

      for(double i=1.0; i<=np/2; i++) {
	xx = xlow + (i-.5) * step;
	yy = timeslice[j]-xx;
	asigma = gsigma+(yy>0.0)*asymmcoef*(yy-0.0);
	fland = TMath::Landau(xx,mpc,lwidth) / lwidth;
	summ += fland * exp(-0.5*pow((yy-0.0)/asigma,2)) / asigma;

	xx = xupp - (i-.5) * step;
	yy = timeslice[j]-xx;
	asigma = gsigma+(yy>0.0)*asymmcoef*(yy-0.0);
	fland = TMath::Landau(xx,mpc,lwidth) / lwidth;
	summ += fland * exp(-0.5*pow((yy-0.0)/asigma,2)) / asigma;
      }
      funcvalue[j] = step * summ * invsq2pi;
    }

    for (int i=0; i<10; i++) {
      if (funcvalue[i]>1e-6)
	Shape.coeffRef(i) = funcvalue[i];///funcvalue[4];
      else Shape.coeffRef(i) = 0;
    }
  };
 
};

#endif
