#include <iostream>
#include <cmath>
#include <climits>
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalDeterministicFit.h"
#include <TMath.h>
#include <TF1.h>

constexpr float HcalDeterministicFit::invGpar[3];
constexpr float HcalDeterministicFit::negThresh[2];
constexpr int HcalDeterministicFit::HcalRegion[2];
constexpr float HcalDeterministicFit::rCorr[2];
constexpr float HcalDeterministicFit::rCorrSiPM[2];

using namespace std;

HcalDeterministicFit::HcalDeterministicFit() {
}

HcalDeterministicFit::~HcalDeterministicFit() { 
}

void HcalDeterministicFit::init(HcalTimeSlew::ParaSource tsParam, HcalTimeSlew::BiasSetting bias, bool iApplyTimeSlew, PedestalSub pedSubFxn_, std::vector<double> pars, double respCorr) {
  for(int fi=0; fi<9; fi++){
	fpars[fi] = pars.at(fi);
  }

  applyTimeSlew_=iApplyTimeSlew;
  useExtPulse_=false;
  fTimeSlew=tsParam;
  fTimeSlewBias=bias;
  fPedestalSubFxn_=pedSubFxn_;
  frespCorr=respCorr;

  /*float flip[10];
  float par0[10][2];
  float par1[10][2];
  float par2[10][2];
  float par3[10][2];*/

  loThresh=10;

  flip[0]=-100; flip[1]=-100; flip[2]=-100; flip[3]=-100; flip[4]=-100; 
  flip[5]=39; 
  flip[6]=-100; flip[7]=-100; 
  flip[8]=110; 
  flip[9]=30;

  for (int i=0; i<10; i++) {
    par0[i][0]=0; par0[i][1]=0;
    par1[i][0]=0; par1[i][1]=0;
    par2[i][0]=0; par2[i][1]=0;
    par3[i][0]=0; par3[i][1]=0;
  }

  par0[3][1] = 0.0792124;
  par1[3][1] = -0.018064;
  par2[3][1] = 0.00135909;

  par0[4][1] = 0.121175;
  par1[4][1] = 0.11762;
  par2[4][1] = -0.00529703;

  par0[5][0] = -0.364612;
  par1[5][0] = 0.421838;
  par2[5][0] = -0.0629547;

  par0[5][1] = 0.488065;
  par1[5][1] = -0.0442216;
  par2[5][1] = 0.000826485;

  par0[6][1] =  0.147936;
  par1[6][1] = -0.000774686;
  par2[6][1] = -0.00353498;
  par3[6][1] = 0.000266972;

  par0[7][1] = 0.0582457;
  par1[7][1] = -0.0111091;
  par2[7][1] = 0.000724373;

  par0[8][0] = 0.0457377;
  par1[8][0] = -0.00809028;

  par0[9][0] = 0.123814;
  par1[9][0] = -0.0347304;

  for (int i=0; i<10; i++) {
    sumPars[0][0]+=par0[i][0];
    sumPars[1][0]+=par1[i][0];
    sumPars[2][0]+=par2[i][0];
    sumPars[3][0]+=par3[i][0];

    sumPars[0][1]+=par0[i][1];
    sumPars[1][1]+=par1[i][1];
    sumPars[2][1]+=par2[i][1];
    sumPars[3][1]+=par3[i][1];
  }

}

void HcalDeterministicFit::setExternalPulseShape(int shape) {

  shape_=shape;
  if(shape_>0) useExtPulse_=true;
}

float HcalDeterministicFit::getPulseFrac(float fC, int ts) const{
  double frac=0;
  if (ts<0 || ts>=10){
    cout << "wrong value for time slice!" << endl;
    return 0;
  }
  double tmpFC=fC;
  if (fC<loThresh) tmpFC=loThresh;

  if (tmpFC < flip[ts]) 
    frac=par0[ts][0] + par1[ts][0]*log(tmpFC) + par2[ts][0]*log(tmpFC)*log(tmpFC) + par3[ts][0]*log(tmpFC)*log(tmpFC)*log(tmpFC);
  else 
    frac=par0[ts][1] + par1[ts][1]*log(tmpFC) + par2[ts][1]*log(tmpFC)*log(tmpFC) + par3[ts][1]*log(tmpFC)*log(tmpFC)*log(tmpFC);

  if (frac>0.01) return frac;
  else return 0;
}

float HcalDeterministicFit::getPulseFracNorm(float fC) const{
  float tempSum=0;
  for (int i=0; i<10; i++) {
    tempSum+=getPulseFrac(fC,i);
  }
  return tempSum;
}

float HcalDeterministicFit::getNegativeEnergyCorr(float fC, float corrTS) const {

  float a4=0, a5=0, b4=0, b5=0, c4=0, c5=0, d4=0, d5=0;

  double tmpFC=fC;
  if (fC<loThresh) tmpFC=loThresh;

  if (tmpFC < flip[4]) {
    a4 = par0[4][0];
    b4 = par1[4][0];
    c4 = par2[4][0];
    d4 = par3[4][0];

    a5 = par0[5][0];
    b5 = par1[5][0];
    c5 = par2[5][0];
    d5 = par3[5][0];

  }
  else {
    a4 = par0[4][1];
    b4 = par1[4][1];
    c4 = par2[4][1];
    d4 = par3[4][1];

    a5 = par0[5][1];
    b5 = par1[5][1];
    c5 = par2[5][1];
    d5 = par3[5][1];

  }

  // set "Q5" pulse to lowest Q template
  /*  if (loThresh < flip[5]) {
  }
  else {
  }*/
  
  TF1 f4("f4","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)",0,3000);
  f4.SetParameter(0,a4);
  f4.SetParameter(1,b4);  
  f4.SetParameter(2,c4);
  f4.SetParameter(3,d4);

  TF1 f5("f5","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)",0,3000);
  f5.SetParameter(0,a5);
  f5.SetParameter(1,b5);  
  f5.SetParameter(2,c5);
  f5.SetParameter(3,d5);

  TF1 ts("ts","([0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x))/([4]+[5]*log(x)+[6]*log(x)*log(x)+[7]*log(x)*log(x)*log(x))",0,3000);
  ts.SetParameter(4,a4);
  ts.SetParameter(5,b4);  
  ts.SetParameter(6,c4);
  ts.SetParameter(7,d4);
  ts.SetParameter(0,a5);
  ts.SetParameter(1,b5);  
  ts.SetParameter(2,c5);
  ts.SetParameter(3,d5);

  if (ts.Eval(1500)>corrTS)
    return 2200;
  else if (ts.Eval(700)>corrTS) 
    return 1100;
  else if (ts.Eval(400)>corrTS) 
    return 800;
  else if (ts.Eval(200)>corrTS) 
    return 300;
  else if (ts.Eval(100)>corrTS) 
    return 50;

  //cout << ts.Eval(loThresh) << ", " << ts.Eval(3000) << endl;

  //cout << fC << ", "<< corrTS << ", " << ts.Eval(fC) << "; " << ts.GetX(ts.Eval(fC)) << endl;

  return 10;

}



constexpr float HcalDeterministicFit::landauFrac[];
// Landau function integrated in 1 ns intervals
//Landau pulse shape from https://indico.cern.ch/event/345283/contribution/3/material/slides/0.pdf
//Landau turn on by default at left edge of time slice 
// normalized to 1 on [0,10000]
void HcalDeterministicFit::getLandauFrac(float tStart, float tEnd, float &sum) const{

  if (std::abs(tStart-tEnd-tsWidth)<0.1) {
    sum=0;
    return;
  }
  sum= landauFrac[int(ceil(tStart+tsWidth))];
  return;
}

void HcalDeterministicFit::getLandauFrac(float fC, int offset, double fpar0, double fpar1, double fpar2, float &sum) const{

  //if (useExtPulse_ == true) {
    // hardcoded array :(
    int chargeBin = -1;
    for (int i=0; i<58; i++) {
      if (fC>minCharge_[i] && fC<maxCharge_[i]) chargeBin=i;
    }
    if (fC>maxCharge_[57]) chargeBin=57;

    if (chargeBin==-1) chargeBin=0;
    //if (chargeBin!=0 && offset==1) std::cout << "chargeBin = " << chargeBin << ", " << minCharge_[chargeBin] << ", " << maxCharge_[chargeBin] << std::endl;
    sum=pulseFrac_[chargeBin][offset+3];
    //std::cout << "sum = " << sum << std::endl;
    return;

    //}
}

void HcalDeterministicFit::phase1Apply(const HBHEChannelInfo& channelData,
				       float& reconstructedEnergy,
				       float& reconstructedTime) const
{

  std::vector<double> corrCharge;
  std::vector<double> inputCharge;
  std::vector<double> inputPedestal;
  double gainCorr = 0;
  double respCorr = 0;

  for(unsigned int ip=0; ip<channelData.nSamples(); ip++){

    double charge = channelData.tsRawCharge(ip);
    double ped = channelData.tsPedestal(ip); 
    double gain = channelData.tsGain(ip);

    gainCorr = gain;
    inputCharge.push_back(charge);
    inputPedestal.push_back(ped);

  }

  if (inputCharge[4] + inputCharge[5] <20) return;

  fPedestalSubFxn_.calculate(inputCharge, inputPedestal, corrCharge);

  const HcalDetId& cell = channelData.id();

  if (fTimeSlew==0)respCorr=1.0;
  else if (fTimeSlew==1) channelData.hasTimeInfo()?respCorr=rCorrSiPM[0]:respCorr=rCorr[0];
  else if (fTimeSlew==2) channelData.hasTimeInfo()?respCorr=rCorrSiPM[1]:respCorr=rCorr[1];
  else if (fTimeSlew==3)respCorr=frespCorr;

  /*  double fpar0, fpar1, fpar2;
  if(std::abs(cell.ieta())<HcalRegion[0]){
    fpar0 = fpars[0];
    fpar1 = fpars[1];
    fpar2 = fpars[2];
  }else if(std::abs(cell.ieta())==HcalRegion[0]||std::abs(cell.ieta())==HcalRegion[1]){
    fpar0 = fpars[3];
    fpar1 = fpars[4];
    fpar2 = fpars[5];
  }else{
    fpar0 = fpars[6];
    fpar1 = fpars[7];
    fpar2 = fpars[8];
    }*/

  float tsShift3=0;
  float tsShift4=0;
  float tsShift5=0;

  float ch3=0;
  float ch4=0;
  float ch5=0;

  //if(useExtPulse_ && shape_==1) {

  float sum3=getPulseFracNorm(inputCharge[3]);
  float sum4=getPulseFracNorm(inputCharge[4]);
  float sum5=getPulseFracNorm(inputCharge[5]);

  //cout << "sums: " << sum3 << ", " << sum4 << ", " << sum5 << endl;

  float i3=getPulseFrac(inputCharge[3],4)/sum3;
  float n3=getPulseFrac(inputCharge[3],5)/sum3;
  float nn3=getPulseFrac(inputCharge[3],6)/sum3;

  float i4=getPulseFrac(inputCharge[4],4)/sum4;
  float n4=getPulseFrac(inputCharge[4],5)/sum4;

  float i5=getPulseFrac(inputCharge[5],4)/sum5;
  float n5=getPulseFrac(inputCharge[5],5)/sum5;

  //cout << "i,n,nn 3: " << i3 << ", " << n3 << ", " << nn3 << endl;
  //cout << "i,n,nn 4: " << i4 << ", " << n4 << endl;
  //cout << "i,n,nn 5: " << i5 << ", " << n5 << endl;

  if (i3 != 0 && i4 != 0 && i5 != 0) {
    
    ch3=corrCharge[3]/i3;
    ch4=(i3*corrCharge[4]-n3*corrCharge[3])/(i3*i4);
    //ch5=(n3*n4*corrCharge[3]-i4*nn3*corrCharge[3]-i3*n4*corrCharge[4]+i3*i4*corrCharge[5])/(i3*i4*i5);
    
    if (ch3<negThresh[0]) {
      ch3=negThresh[0];
      ch4=corrCharge[4]/i4;
      //ch5=(i4*corrCharge[5]-n4*corrCharge[4])/(i4*i5);
    }
    if (ch5<negThresh[0] && ch4>negThresh[1]) {
      
      float newTS = (corrCharge[5]-negThresh[0]*i5)/(corrCharge[4]-ch3*i3);
      float newQ = getNegativeEnergyCorr(corrCharge[4], newTS);
      
      //cout << newTS << ", " << corrCharge[4] << "; " << getPulseFrac(newQ,5)/getPulseFrac(newQ,4) << ", " << newQ << endl;
      
      float i4_new = getPulseFrac(newQ,4)/getPulseFracNorm(newQ);
      
      if (i4_new!=0){
	ch5=negThresh[0];
	ch4=(corrCharge[4]-ch3*n3)/(i4_new);
      }
      ch4=0;
    }
    
  }
  
  if (ch4<1) ch4=0;
  
  /*
    tsShift3=HcalTimeSlew::delay(inputCharge[3], fTimeSlew, fTimeSlewBias, fpar0, fpar1 ,fpar2,!channelData.hasTimeInfo());
    tsShift4=HcalTimeSlew::delay(inputCharge[4], fTimeSlew, fTimeSlewBias, fpar0, fpar1 ,fpar2,!channelData.hasTimeInfo());
    tsShift5=HcalTimeSlew::delay(inputCharge[5], fTimeSlew, fTimeSlewBias, fpar0, fpar1 ,fpar2,!channelData.hasTimeInfo());

    //}

  float i32=0;
  getLandauFrac(-tsShift3,-tsShift3+tsWidth,i32);
  float n32=0;
  getLandauFrac(-tsShift3+tsWidth,-tsShift3+tsWidth*2,n32);
  float nn32=0;
  getLandauFrac(-tsShift3+tsWidth*2,-tsShift3+tsWidth*3,nn32);

  float i42=0;
  getLandauFrac(-tsShift4,-tsShift4+tsWidth,i42);
  float n42=0;
  getLandauFrac(-tsShift4+tsWidth,-tsShift4+tsWidth*2,n42);

  float i52=0;
  getLandauFrac(-tsShift5,-tsShift5+tsWidth,i52);
  float n52=0;
  getLandauFrac(-tsShift5+tsWidth,-tsShift5+tsWidth*2,n52);

  if (i3 != 0 && i4 != 0 && i5 != 0) {

    ch3=corrCharge[3]/i3;
    ch4=(i3*corrCharge[4]-n3*corrCharge[3])/(i3*i4);
    ch5=(n3*n4*corrCharge[3]-i4*nn3*corrCharge[3]-i3*n4*corrCharge[4]+i3*i4*corrCharge[5])/(i3*i4*i5);

    if (ch3<negThresh[0]) {
      ch3=negThresh[0];
      ch4=corrCharge[4]/i4;
      ch5=(i4*corrCharge[5]-n4*corrCharge[4])/(i4*i5);
    }
    if (ch5<negThresh[0] && ch4>negThresh[1]) {
      double ratio = (corrCharge[4]-ch3*i3)/(corrCharge[5]-negThresh[0]*i5);
      if (ratio < 5 && ratio > 0.5) {
        double invG = invGpar[0]+invGpar[1]*std::sqrt(2*std::log(invGpar[2]/ratio));
        float iG=0;
	getLandauFrac(-invG,-invG+tsWidth,iG);
	if (iG != 0 ) {
	  ch4=(corrCharge[4]-ch3*n3)/(iG);
	  tsShift4=invG;
	}
      }
    }
  }

  if (ch4<1) ch4=0;

  }
  */
  reconstructedEnergy=ch4*gainCorr*respCorr;
  reconstructedTime=tsShift4;

}
