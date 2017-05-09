#include <iostream>
#include <cmath>
#include <climits>
#include "RecoLocalCalo/HcalRecAlgos/interface/NewPulseShapes.h"
#include <TMath.h>
#include <TF1.h>

using namespace std;

NewPulseShapes::NewPulseShapes() {
}

NewPulseShapes::~NewPulseShapes() { 
}

void NewPulseShapes::init() {
  configurePulseShapes();
}
void NewPulseShapes::configurePulseShapes() {

  loThresh=5;
  hiThresh=3000;

  flip[0]=-100; flip[1]=-100; flip[2]=-100; flip[3]=-100; flip[4]=-100;
  flip[5]=-100;
  flip[6]=-100; flip[7]=-100;
  flip[8]=-100;
  flip[9]=-100;

  for (int i=0; i<10; i++) {
    par0[i][0]=0; par0[i][1]=0;
    par1[i][0]=0; par1[i][1]=0;
    par2[i][0]=0; par2[i][1]=0;
    par3[i][0]=0; par3[i][1]=0;
    tpar0[i]=0;
    tpar1[i]=0;
    tpar2[i]=0;
  }

  par0[3][1] = 0.0960304;
  par1[3][1] = -0.0272066;
  par2[3][1] = 0.00230593;

  par0[4][1] = 0.323146;
  par1[4][1] = 0.067774;
  par2[4][1] = -0.00160776;

  par0[5][1] = 0.390316;
  par1[5][1] = -0.0107445;
  par2[5][1] = -0.00232944;

  par0[6][1] = 0.141146;
  par1[6][1] = -0.0120172;
  par2[6][1] = -5.27739e-05;

  tpar2[3] = -0.0437093;

  tpar0[4] = 0.565593;
  tpar1[4] = -0.00493221;
  tpar2[4] = -0.499619;

  tpar0[5] = -0.362077;
  tpar1[5] = -0.00286609;
  tpar2[5] = 0.402905;

  tpar0[6] = -0.135147;
  tpar1[6] = -0.00585069;
  tpar2[6] = 0.132601;

}

float NewPulseShapes::getPulseFrac(float fC, float time, int ts) const{
  float frac=0;

  if (ts<0 || ts>=10){
    cout << "wrong value for time slice!" << endl;
    return 0;
  }

  float tmpFC=fC;
  if (fC<10) { 
    tmpFC=10;
    //std::cout << "vv" << std::endl;
  }
  else if (fC>3000) {
    tmpFC=3000;
    //std::cout << "^^" << std::endl;
  }

  if (tmpFC < flip[ts]) {
    frac=par0[ts][0] + par1[ts][0]*log(tmpFC) + par2[ts][0]*log(tmpFC)*log(tmpFC) + par3[ts][0]*log(tmpFC)*log(tmpFC)*log(tmpFC);
  }
  else {
    frac=par0[ts][1] + par1[ts][1]*log(tmpFC) + par2[ts][1]*log(tmpFC)*log(tmpFC) + par3[ts][1]*log(tmpFC)*log(tmpFC)*log(tmpFC);
  }
  frac+=(tpar0[ts]*exp(tpar1[ts]*tmpFC)+tpar2[ts])*time; //hack b/c my derivatives are in fractional TS, not time [ns]

  if (frac>0.01) return frac;
  else return 0;
}

float NewPulseShapes::getPulseFracNorm(float fC, float time) const{
  float tempSum=0;
  //std::cout << "get pulse of " << fC << ", " << time << ": ";
  //std::cout << getPulseFrac(fC,0,4) << ", " << getPulseFrac(fC,0,5) << std::endl;
  for (int i=0; i<10; i++) {
    tempSum+=getPulseFrac(fC,time,i);
    //std::cout << getPulseFrac(fC,0,i) << ", ";
  }
  return tempSum;
}
