#include <iostream>
#include <cmath>
#include <climits>
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalDeterministicFit.h"
#include <TMath.h>

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

  pulseFrac_[0][0]= 0; pulseFrac_[0][1]= 0; pulseFrac_[0][2]= 0; pulseFrac_[0][3]= 0.0297409; pulseFrac_[0][4]= 0.450687; pulseFrac_[0][5]= 0.341974; pulseFrac_[0][6]= 0.119703; pulseFrac_[0][7]= 0.0287843; pulseFrac_[0][8]= 0.0165336; pulseFrac_[0][9]= 0.0125778; 
  pulseFrac_[1][0]= 0; pulseFrac_[1][1]= 0; pulseFrac_[1][2]= 0; pulseFrac_[1][3]= 0.0272083; pulseFrac_[1][4]= 0.481118; pulseFrac_[1][5]= 0.337939; pulseFrac_[1][6]= 0.113655; pulseFrac_[1][7]= 0.0260069; pulseFrac_[1][8]= 0.0140733; pulseFrac_[1][9]= 0; 
  pulseFrac_[2][0]= 0; pulseFrac_[2][1]= 0; pulseFrac_[2][2]= 0; pulseFrac_[2][3]= 0.0265275; pulseFrac_[2][4]= 0.497913; pulseFrac_[2][5]= 0.329564; pulseFrac_[2][6]= 0.10861; pulseFrac_[2][7]= 0.0245201; pulseFrac_[2][8]= 0.0128656; pulseFrac_[2][9]= 0; 
  pulseFrac_[3][0]= 0; pulseFrac_[3][1]= 0; pulseFrac_[3][2]= 0; pulseFrac_[3][3]= 0.0264219; pulseFrac_[3][4]= 0.511151; pulseFrac_[3][5]= 0.322014; pulseFrac_[3][6]= 0.104761; pulseFrac_[3][7]= 0.0235931; pulseFrac_[3][8]= 0.0120597; pulseFrac_[3][9]= 0; 
  pulseFrac_[4][0]= 0; pulseFrac_[4][1]= 0; pulseFrac_[4][2]= 0; pulseFrac_[4][3]= 0.0264452; pulseFrac_[4][4]= 0.520779; pulseFrac_[4][5]= 0.316401; pulseFrac_[4][6]= 0.101922; pulseFrac_[4][7]= 0.0229341; pulseFrac_[4][8]= 0.0115185; pulseFrac_[4][9]= 0; 
  pulseFrac_[5][0]= 0; pulseFrac_[5][1]= 0; pulseFrac_[5][2]= 0; pulseFrac_[5][3]= 0.0265417; pulseFrac_[5][4]= 0.529211; pulseFrac_[5][5]= 0.311395; pulseFrac_[5][6]= 0.0993968; pulseFrac_[5][7]= 0.022417; pulseFrac_[5][8]= 0.0110384; pulseFrac_[5][9]= 0; 
  pulseFrac_[6][0]= 0; pulseFrac_[6][1]= 0; pulseFrac_[6][2]= 0; pulseFrac_[6][3]= 0.0265141; pulseFrac_[6][4]= 0.536232; pulseFrac_[6][5]= 0.30717; pulseFrac_[6][6]= 0.0973737; pulseFrac_[6][7]= 0.0220187; pulseFrac_[6][8]= 0.0106913; pulseFrac_[6][9]= 0; 
  pulseFrac_[7][0]= 0; pulseFrac_[7][1]= 0; pulseFrac_[7][2]= 0; pulseFrac_[7][3]= 0.0263896; pulseFrac_[7][4]= 0.542466; pulseFrac_[7][5]= 0.303487; pulseFrac_[7][6]= 0.0956399; pulseFrac_[7][7]= 0.0216533; pulseFrac_[7][8]= 0.0103645; pulseFrac_[7][9]= 0; 
  pulseFrac_[8][0]= 0; pulseFrac_[8][1]= 0; pulseFrac_[8][2]= 0; pulseFrac_[8][3]= 0.0262053; pulseFrac_[8][4]= 0.547516; pulseFrac_[8][5]= 0.300542; pulseFrac_[8][6]= 0.0942286; pulseFrac_[8][7]= 0.0213631; pulseFrac_[8][8]= 0.0101447; pulseFrac_[8][9]= 0; 
  pulseFrac_[9][0]= 0; pulseFrac_[9][1]= 0; pulseFrac_[9][2]= 0; pulseFrac_[9][3]= 0.0261446; pulseFrac_[9][4]= 0.55861; pulseFrac_[9][5]= 0.30027; pulseFrac_[9][6]= 0.0937033; pulseFrac_[9][7]= 0.0212717; pulseFrac_[9][8]= 0; pulseFrac_[9][9]= 0; 
  pulseFrac_[10][0]= 0; pulseFrac_[10][1]= 0; pulseFrac_[10][2]= 0; pulseFrac_[10][3]= 0.0257579; pulseFrac_[10][4]= 0.56418; pulseFrac_[10][5]= 0.296926; pulseFrac_[10][6]= 0.0921686; pulseFrac_[10][7]= 0.0209668; pulseFrac_[10][8]= 0; pulseFrac_[10][9]= 0; 
  pulseFrac_[11][0]= 0; pulseFrac_[11][1]= 0; pulseFrac_[11][2]= 0; pulseFrac_[11][3]= 0.0250805; pulseFrac_[11][4]= 0.569149; pulseFrac_[11][5]= 0.29421; pulseFrac_[11][6]= 0.0908511; pulseFrac_[11][7]= 0.0207095; pulseFrac_[11][8]= 0; pulseFrac_[11][9]= 0; 
  pulseFrac_[12][0]= 0; pulseFrac_[12][1]= 0; pulseFrac_[12][2]= 0; pulseFrac_[12][3]= 0.0240906; pulseFrac_[12][4]= 0.573843; pulseFrac_[12][5]= 0.291884; pulseFrac_[12][6]= 0.0897216; pulseFrac_[12][7]= 0.0204606; pulseFrac_[12][8]= 0; pulseFrac_[12][9]= 0; 
  pulseFrac_[13][0]= 0; pulseFrac_[13][1]= 0; pulseFrac_[13][2]= 0; pulseFrac_[13][3]= 0.0231029; pulseFrac_[13][4]= 0.578435; pulseFrac_[13][5]= 0.28956; pulseFrac_[13][6]= 0.0886047; pulseFrac_[13][7]= 0.0202976; pulseFrac_[13][8]= 0; pulseFrac_[13][9]= 0; 
  pulseFrac_[14][0]= 0; pulseFrac_[14][1]= 0; pulseFrac_[14][2]= 0; pulseFrac_[14][3]= 0.0222629; pulseFrac_[14][4]= 0.583002; pulseFrac_[14][5]= 0.287123; pulseFrac_[14][6]= 0.0874958; pulseFrac_[14][7]= 0.0201168; pulseFrac_[14][8]= 0; pulseFrac_[14][9]= 0; 
  pulseFrac_[15][0]= 0; pulseFrac_[15][1]= 0; pulseFrac_[15][2]= 0; pulseFrac_[15][3]= 0.0215465; pulseFrac_[15][4]= 0.586667; pulseFrac_[15][5]= 0.285279; pulseFrac_[15][6]= 0.0865954; pulseFrac_[15][7]= 0.0199122; pulseFrac_[15][8]= 0; pulseFrac_[15][9]= 0; 
  pulseFrac_[16][0]= 0; pulseFrac_[16][1]= 0; pulseFrac_[16][2]= 0; pulseFrac_[16][3]= 0.0209442; pulseFrac_[16][4]= 0.590551; pulseFrac_[16][5]= 0.283073; pulseFrac_[16][6]= 0.0856525; pulseFrac_[16][7]= 0.019779; pulseFrac_[16][8]= 0; pulseFrac_[16][9]= 0; 
  pulseFrac_[17][0]= 0; pulseFrac_[17][1]= 0; pulseFrac_[17][2]= 0; pulseFrac_[17][3]= 0.0204287; pulseFrac_[17][4]= 0.594137; pulseFrac_[17][5]= 0.281056; pulseFrac_[17][6]= 0.0847504; pulseFrac_[17][7]= 0.0196285; pulseFrac_[17][8]= 0; pulseFrac_[17][9]= 0; 
  pulseFrac_[18][0]= 0; pulseFrac_[18][1]= 0; pulseFrac_[18][2]= 0; pulseFrac_[18][3]= 0.0200044; pulseFrac_[18][4]= 0.597266; pulseFrac_[18][5]= 0.279254; pulseFrac_[18][6]= 0.0839624; pulseFrac_[18][7]= 0.0195131; pulseFrac_[18][8]= 0; pulseFrac_[18][9]= 0; 
  pulseFrac_[19][0]= 0; pulseFrac_[19][1]= 0; pulseFrac_[19][2]= 0; pulseFrac_[19][3]= 0.0196551; pulseFrac_[19][4]= 0.600902; pulseFrac_[19][5]= 0.277055; pulseFrac_[19][6]= 0.0830468; pulseFrac_[19][7]= 0.0193405; pulseFrac_[19][8]= 0; pulseFrac_[19][9]= 0; 
  pulseFrac_[20][0]= 0; pulseFrac_[20][1]= 0; pulseFrac_[20][2]= 0; pulseFrac_[20][3]= 0.0193426; pulseFrac_[20][4]= 0.603498; pulseFrac_[20][5]= 0.275623; pulseFrac_[20][6]= 0.0823089; pulseFrac_[20][7]= 0.0192276; pulseFrac_[20][8]= 0; pulseFrac_[20][9]= 0; 
  pulseFrac_[21][0]= 0; pulseFrac_[21][1]= 0; pulseFrac_[21][2]= 0; pulseFrac_[21][3]= 0.0191163; pulseFrac_[21][4]= 0.606443; pulseFrac_[21][5]= 0.273784; pulseFrac_[21][6]= 0.0815799; pulseFrac_[21][7]= 0.0190766; pulseFrac_[21][8]= 0; pulseFrac_[21][9]= 0; 
  pulseFrac_[22][0]= 0; pulseFrac_[22][1]= 0; pulseFrac_[22][2]= 0; pulseFrac_[22][3]= 0.0189606; pulseFrac_[22][4]= 0.608762; pulseFrac_[22][5]= 0.272322; pulseFrac_[22][6]= 0.0809561; pulseFrac_[22][7]= 0.0189994; pulseFrac_[22][8]= 0; pulseFrac_[22][9]= 0; 
  pulseFrac_[23][0]= 0; pulseFrac_[23][1]= 0; pulseFrac_[23][2]= 0; pulseFrac_[23][3]= 0.0187815; pulseFrac_[23][4]= 0.610771; pulseFrac_[23][5]= 0.271125; pulseFrac_[23][6]= 0.0804549; pulseFrac_[23][7]= 0.0188677; pulseFrac_[23][8]= 0; pulseFrac_[23][9]= 0; 
  pulseFrac_[24][0]= 0; pulseFrac_[24][1]= 0; pulseFrac_[24][2]= 0; pulseFrac_[24][3]= 0.0186165; pulseFrac_[24][4]= 0.612755; pulseFrac_[24][5]= 0.27; pulseFrac_[24][6]= 0.0798709; pulseFrac_[24][7]= 0.0187574; pulseFrac_[24][8]= 0; pulseFrac_[24][9]= 0; 
  pulseFrac_[25][0]= 0; pulseFrac_[25][1]= 0; pulseFrac_[25][2]= 0; pulseFrac_[25][3]= 0.0185591; pulseFrac_[25][4]= 0.614842; pulseFrac_[25][5]= 0.268555; pulseFrac_[25][6]= 0.0793308; pulseFrac_[25][7]= 0.0187124; pulseFrac_[25][8]= 0; pulseFrac_[25][9]= 0; 
  pulseFrac_[26][0]= 0; pulseFrac_[26][1]= 0; pulseFrac_[26][2]= 0; pulseFrac_[26][3]= 0.0185008; pulseFrac_[26][4]= 0.616771; pulseFrac_[26][5]= 0.267366; pulseFrac_[26][6]= 0.0787765; pulseFrac_[26][7]= 0.018585; pulseFrac_[26][8]= 0; pulseFrac_[26][9]= 0; 
  pulseFrac_[27][0]= 0; pulseFrac_[27][1]= 0; pulseFrac_[27][2]= 0; pulseFrac_[27][3]= 0.018486; pulseFrac_[27][4]= 0.618888; pulseFrac_[27][5]= 0.265864; pulseFrac_[27][6]= 0.0782841; pulseFrac_[27][7]= 0.0184778; pulseFrac_[27][8]= 0; pulseFrac_[27][9]= 0; 
  pulseFrac_[28][0]= 0; pulseFrac_[28][1]= 0; pulseFrac_[28][2]= 0; pulseFrac_[28][3]= 0.0184235; pulseFrac_[28][4]= 0.620561; pulseFrac_[28][5]= 0.264815; pulseFrac_[28][6]= 0.077793; pulseFrac_[28][7]= 0.0184071; pulseFrac_[28][8]= 0; pulseFrac_[28][9]= 0; 
  pulseFrac_[29][0]= 0; pulseFrac_[29][1]= 0; pulseFrac_[29][2]= 0; pulseFrac_[29][3]= 0.018404; pulseFrac_[29][4]= 0.622415; pulseFrac_[29][5]= 0.263508; pulseFrac_[29][6]= 0.0773295; pulseFrac_[29][7]= 0.0183434; pulseFrac_[29][8]= 0; pulseFrac_[29][9]= 0; 
  pulseFrac_[30][0]= 0; pulseFrac_[30][1]= 0; pulseFrac_[30][2]= 0; pulseFrac_[30][3]= 0.0184814; pulseFrac_[30][4]= 0.624184; pulseFrac_[30][5]= 0.262209; pulseFrac_[30][6]= 0.076844; pulseFrac_[30][7]= 0.0182813; pulseFrac_[30][8]= 0; pulseFrac_[30][9]= 0; 
  pulseFrac_[31][0]= 0; pulseFrac_[31][1]= 0; pulseFrac_[31][2]= 0; pulseFrac_[31][3]= 0.0184184; pulseFrac_[31][4]= 0.6264; pulseFrac_[31][5]= 0.260758; pulseFrac_[31][6]= 0.0762831; pulseFrac_[31][7]= 0.0181404; pulseFrac_[31][8]= 0; pulseFrac_[31][9]= 0; 
  pulseFrac_[32][0]= 0; pulseFrac_[32][1]= 0; pulseFrac_[32][2]= 0; pulseFrac_[32][3]= 0.0184492; pulseFrac_[32][4]= 0.628061; pulseFrac_[32][5]= 0.259537; pulseFrac_[32][6]= 0.0758762; pulseFrac_[32][7]= 0.0180771; pulseFrac_[32][8]= 0; pulseFrac_[32][9]= 0; 
  pulseFrac_[33][0]= 0; pulseFrac_[33][1]= 0; pulseFrac_[33][2]= 0; pulseFrac_[33][3]= 0.018456; pulseFrac_[33][4]= 0.629367; pulseFrac_[33][5]= 0.258605; pulseFrac_[33][6]= 0.0754965; pulseFrac_[33][7]= 0.0180747; pulseFrac_[33][8]= 0; pulseFrac_[33][9]= 0; 
  pulseFrac_[34][0]= 0; pulseFrac_[34][1]= 0; pulseFrac_[34][2]= 0; pulseFrac_[34][3]= 0.0184901; pulseFrac_[34][4]= 0.63077; pulseFrac_[34][5]= 0.257678; pulseFrac_[34][6]= 0.0751319; pulseFrac_[34][7]= 0.0179305; pulseFrac_[34][8]= 0; pulseFrac_[34][9]= 0; 
  pulseFrac_[35][0]= 0; pulseFrac_[35][1]= 0; pulseFrac_[35][2]= 0; pulseFrac_[35][3]= 0.0185222; pulseFrac_[35][4]= 0.632028; pulseFrac_[35][5]= 0.256833; pulseFrac_[35][6]= 0.0747378; pulseFrac_[35][7]= 0.0178787; pulseFrac_[35][8]= 0; pulseFrac_[35][9]= 0; 
  pulseFrac_[36][0]= 0; pulseFrac_[36][1]= 0; pulseFrac_[36][2]= 0; pulseFrac_[36][3]= 0.0185079; pulseFrac_[36][4]= 0.633404; pulseFrac_[36][5]= 0.255763; pulseFrac_[36][6]= 0.0744523; pulseFrac_[36][7]= 0.0178721; pulseFrac_[36][8]= 0; pulseFrac_[36][9]= 0; 
  pulseFrac_[37][0]= 0; pulseFrac_[37][1]= 0; pulseFrac_[37][2]= 0; pulseFrac_[37][3]= 0.0185818; pulseFrac_[37][4]= 0.634471; pulseFrac_[37][5]= 0.254988; pulseFrac_[37][6]= 0.0741514; pulseFrac_[37][7]= 0.0178078; pulseFrac_[37][8]= 0; pulseFrac_[37][9]= 0; 
  pulseFrac_[38][0]= 0; pulseFrac_[38][1]= 0; pulseFrac_[38][2]= 0; pulseFrac_[38][3]= 0.0186571; pulseFrac_[38][4]= 0.635789; pulseFrac_[38][5]= 0.254058; pulseFrac_[38][6]= 0.0737301; pulseFrac_[38][7]= 0.0177656; pulseFrac_[38][8]= 0; pulseFrac_[38][9]= 0; 
  pulseFrac_[39][0]= 0; pulseFrac_[39][1]= 0; pulseFrac_[39][2]= 0; pulseFrac_[39][3]= 0.0186282; pulseFrac_[39][4]= 0.63717; pulseFrac_[39][5]= 0.253115; pulseFrac_[39][6]= 0.0734096; pulseFrac_[39][7]= 0.0176766; pulseFrac_[39][8]= 0; pulseFrac_[39][9]= 0; 
  pulseFrac_[40][0]= 0; pulseFrac_[40][1]= 0; pulseFrac_[40][2]= 0; pulseFrac_[40][3]= 0.0187319; pulseFrac_[40][4]= 0.639093; pulseFrac_[40][5]= 0.251676; pulseFrac_[40][6]= 0.0728604; pulseFrac_[40][7]= 0.0176386; pulseFrac_[40][8]= 0; pulseFrac_[40][9]= 0; 
  pulseFrac_[41][0]= 0; pulseFrac_[41][1]= 0; pulseFrac_[41][2]= 0; pulseFrac_[41][3]= 0.0187251; pulseFrac_[41][4]= 0.640363; pulseFrac_[41][5]= 0.250851; pulseFrac_[41][6]= 0.072522; pulseFrac_[41][7]= 0.0175388; pulseFrac_[41][8]= 0; pulseFrac_[41][9]= 0; 
  pulseFrac_[42][0]= 0; pulseFrac_[42][1]= 0; pulseFrac_[42][2]= 0; pulseFrac_[42][3]= 0.0188316; pulseFrac_[42][4]= 0.642012; pulseFrac_[42][5]= 0.249589; pulseFrac_[42][6]= 0.0720674; pulseFrac_[42][7]= 0.0175001; pulseFrac_[42][8]= 0; pulseFrac_[42][9]= 0; 
  pulseFrac_[43][0]= 0; pulseFrac_[43][1]= 0; pulseFrac_[43][2]= 0; pulseFrac_[43][3]= 0.0187956; pulseFrac_[43][4]= 0.642613; pulseFrac_[43][5]= 0.249165; pulseFrac_[43][6]= 0.0719332; pulseFrac_[43][7]= 0.0174937; pulseFrac_[43][8]= 0; pulseFrac_[43][9]= 0; 
  pulseFrac_[44][0]= 0; pulseFrac_[44][1]= 0; pulseFrac_[44][2]= 0; pulseFrac_[44][3]= 0.0188741; pulseFrac_[44][4]= 0.643225; pulseFrac_[44][5]= 0.248612; pulseFrac_[44][6]= 0.0717685; pulseFrac_[44][7]= 0.0175203; pulseFrac_[44][8]= 0; pulseFrac_[44][9]= 0; 
  pulseFrac_[45][0]= 0; pulseFrac_[45][1]= 0; pulseFrac_[45][2]= 0; pulseFrac_[45][3]= 0.0188753; pulseFrac_[45][4]= 0.644537; pulseFrac_[45][5]= 0.247772; pulseFrac_[45][6]= 0.0713626; pulseFrac_[45][7]= 0.0174527; pulseFrac_[45][8]= 0; pulseFrac_[45][9]= 0; 
  pulseFrac_[46][0]= 0; pulseFrac_[46][1]= 0; pulseFrac_[46][2]= 0; pulseFrac_[46][3]= 0.0189965; pulseFrac_[46][4]= 0.645306; pulseFrac_[46][5]= 0.24719; pulseFrac_[46][6]= 0.0711349; pulseFrac_[46][7]= 0.0173724; pulseFrac_[46][8]= 0; pulseFrac_[46][9]= 0; 
  pulseFrac_[47][0]= 0; pulseFrac_[47][1]= 0; pulseFrac_[47][2]= 0; pulseFrac_[47][3]= 0.0189362; pulseFrac_[47][4]= 0.646786; pulseFrac_[47][5]= 0.246167; pulseFrac_[47][6]= 0.0707882; pulseFrac_[47][7]= 0.0173222; pulseFrac_[47][8]= 0; pulseFrac_[47][9]= 0; 
  pulseFrac_[48][0]= 0; pulseFrac_[48][1]= 0; pulseFrac_[48][2]= 0; pulseFrac_[48][3]= 0.0191089; pulseFrac_[48][4]= 0.647658; pulseFrac_[48][5]= 0.245393; pulseFrac_[48][6]= 0.0705472; pulseFrac_[48][7]= 0.0172931; pulseFrac_[48][8]= 0; pulseFrac_[48][9]= 0; 
  pulseFrac_[49][0]= 0; pulseFrac_[49][1]= 0; pulseFrac_[49][2]= 0; pulseFrac_[49][3]= 0.0190305; pulseFrac_[49][4]= 0.648874; pulseFrac_[49][5]= 0.244625; pulseFrac_[49][6]= 0.0702139; pulseFrac_[49][7]= 0.0172561; pulseFrac_[49][8]= 0; pulseFrac_[49][9]= 0; 
  pulseFrac_[50][0]= 0; pulseFrac_[50][1]= 0; pulseFrac_[50][2]= 0; pulseFrac_[50][3]= 0.0191553; pulseFrac_[50][4]= 0.649921; pulseFrac_[50][5]= 0.243874; pulseFrac_[50][6]= 0.0698875; pulseFrac_[50][7]= 0.0171628; pulseFrac_[50][8]= 0; pulseFrac_[50][9]= 0; 
  pulseFrac_[51][0]= 0; pulseFrac_[51][1]= 0; pulseFrac_[51][2]= 0; pulseFrac_[51][3]= 0.0191724; pulseFrac_[51][4]= 0.651095; pulseFrac_[51][5]= 0.242934; pulseFrac_[51][6]= 0.0696141; pulseFrac_[51][7]= 0.0171843; pulseFrac_[51][8]= 0; pulseFrac_[51][9]= 0; 
  pulseFrac_[52][0]= 0; pulseFrac_[52][1]= 0; pulseFrac_[52][2]= 0; pulseFrac_[52][3]= 0.0192112; pulseFrac_[52][4]= 0.651515; pulseFrac_[52][5]= 0.242665; pulseFrac_[52][6]= 0.0694841; pulseFrac_[52][7]= 0.0171251; pulseFrac_[52][8]= 0; pulseFrac_[52][9]= 0; 
  pulseFrac_[53][0]= 0; pulseFrac_[53][1]= 0; pulseFrac_[53][2]= 0; pulseFrac_[53][3]= 0.0191484; pulseFrac_[53][4]= 0.652073; pulseFrac_[53][5]= 0.24237; pulseFrac_[53][6]= 0.0693094; pulseFrac_[53][7]= 0.0170997; pulseFrac_[53][8]= 0; pulseFrac_[53][9]= 0; 
  pulseFrac_[54][0]= 0; pulseFrac_[54][1]= 0; pulseFrac_[54][2]= 0; pulseFrac_[54][3]= 0.0192645; pulseFrac_[54][4]= 0.651668; pulseFrac_[54][5]= 0.242639; pulseFrac_[54][6]= 0.0693267; pulseFrac_[54][7]= 0.0171015; pulseFrac_[54][8]= 0; pulseFrac_[54][9]= 0; 
  pulseFrac_[55][0]= 0; pulseFrac_[55][1]= 0; pulseFrac_[55][2]= 0; pulseFrac_[55][3]= 0.0192334; pulseFrac_[55][4]= 0.652532; pulseFrac_[55][5]= 0.241872; pulseFrac_[55][6]= 0.0692153; pulseFrac_[55][7]= 0.017147; pulseFrac_[55][8]= 0; pulseFrac_[55][9]= 0; 
  pulseFrac_[56][0]= 0; pulseFrac_[56][1]= 0; pulseFrac_[56][2]= 0; pulseFrac_[56][3]= 0.0192247; pulseFrac_[56][4]= 0.654486; pulseFrac_[56][5]= 0.240464; pulseFrac_[56][6]= 0.0687251; pulseFrac_[56][7]= 0.0170999; pulseFrac_[56][8]= 0; pulseFrac_[56][9]= 0; 
  pulseFrac_[57][0]= 0; pulseFrac_[57][1]= 0; pulseFrac_[57][2]= 0; pulseFrac_[57][3]= 0.0192489; pulseFrac_[57][4]= 0.654718; pulseFrac_[57][5]= 0.240403; pulseFrac_[57][6]= 0.0686233; pulseFrac_[57][7]= 0.0170065; pulseFrac_[57][8]= 0; pulseFrac_[57][9]= 0; 

  minCharge_[0] = 0;
  maxCharge_[0] = 30;
  for (int j = 0; j < 10; j++) pulseFracDeriv_[0][j] = 0;
  if (pulseFrac_[0][4] != 0)
    timeSlew_[0] = pulseFrac_[0][5]/pulseFrac_[0][4];
  else timeSlew_[0] = 0;

  int k = 1;
  for (int i=35; i<600; i+=10)
    {
      minCharge_[k] = i-5;
      maxCharge_[k] = i+5;
      for (int j = 0; j < 10; j++) pulseFracDeriv_[k][j] = 0;
      if (pulseFrac_[k][4] != 0)
        timeSlew_[k] = pulseFrac_[k][5]/pulseFrac_[k][4];
      else timeSlew_[k] = 0;
      k++;
    }

}

void HcalDeterministicFit::setExternalPulseShape(int shape) {

  shape_=shape;
  if(shape_>0) useExtPulse_=true;
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
    if (chargeBin!=0) std::cout << "chargeBin = " << chargeBin << std::endl;
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

  fPedestalSubFxn_.calculate(inputCharge, inputPedestal, corrCharge);

  const HcalDetId& cell = channelData.id();

  if (fTimeSlew==0)respCorr=1.0;
  else if (fTimeSlew==1) channelData.hasTimeInfo()?respCorr=rCorrSiPM[0]:respCorr=rCorr[0];
  else if (fTimeSlew==2) channelData.hasTimeInfo()?respCorr=rCorrSiPM[1]:respCorr=rCorr[1];
  else if (fTimeSlew==3)respCorr=frespCorr;

  double fpar0, fpar1, fpar2;
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
  }

  float tsShift3=0;
  float tsShift4=0;
  float tsShift5=0;

  float ch3=0;
  float ch4=0;
  float ch5=0;

  //if(useExtPulse_ && shape_==1) {


    float i3=0;
    getLandauFrac(inputCharge[3], 1, fpar0, fpar1, fpar2, i3);
    float n3=0;
    getLandauFrac(inputCharge[3], 2, fpar0, fpar1, fpar2, n3);
    float nn3=0;
    getLandauFrac(inputCharge[3], 3, fpar0, fpar1, fpar2, nn3);

    float i4=0;
    getLandauFrac(inputCharge[4], 1, fpar0, fpar1, fpar2, i4);
    float n4=0;
    getLandauFrac(inputCharge[4], 2, fpar0, fpar1, fpar2, n4);

    float i5=0;
    getLandauFrac(inputCharge[5], 1, fpar0, fpar1, fpar2, i5);
    float n5=0;
    getLandauFrac(inputCharge[5], 2, fpar0, fpar1, fpar2, n5);

    //std:: cout << inputCharge[4] << ", " << i4 << endl;

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
	float newTS = (corrCharge[5]-negThresh[0]*i5)/(corrCharge[4]-ch3*i3);
	int newBin=0;
	for (int k=0; k<58; k++) {
	  if (newTS < timeSlew_[k]) newBin=k;
	}
	float i4_new = pulseFrac_[newBin][4];

	if (i4_new!=0)
	  {
	    ch5=negThresh[0];
	    ch4=(corrCharge[4]-ch3*n3)/(i4_new);
	  }
      }
    }

    if (ch4<1) ch4=0;

    //std::cout << ch4 << std::endl;
    /*
  } else {
  //default Run2 like

  if(applyTimeSlew_) {

    tsShift3=HcalTimeSlew::delay(inputCharge[3], fTimeSlew, fTimeSlewBias, fpar0, fpar1 ,fpar2,!channelData.hasTimeInfo());
    tsShift4=HcalTimeSlew::delay(inputCharge[4], fTimeSlew, fTimeSlewBias, fpar0, fpar1 ,fpar2,!channelData.hasTimeInfo());
    tsShift5=HcalTimeSlew::delay(inputCharge[5], fTimeSlew, fTimeSlewBias, fpar0, fpar1 ,fpar2,!channelData.hasTimeInfo());

  }

  float i3=0;
  getLandauFrac(-tsShift3,-tsShift3+tsWidth,i3);
  float n3=0;
  getLandauFrac(-tsShift3+tsWidth,-tsShift3+tsWidth*2,n3);
  float nn3=0;
  getLandauFrac(-tsShift3+tsWidth*2,-tsShift3+tsWidth*3,nn3);

  float i4=0;
  getLandauFrac(-tsShift4,-tsShift4+tsWidth,i4);
  float n4=0;
  getLandauFrac(-tsShift4+tsWidth,-tsShift4+tsWidth*2,n4);

  float i5=0;
  getLandauFrac(-tsShift5,-tsShift5+tsWidth,i5);
  float n5=0;
  getLandauFrac(-tsShift5+tsWidth,-tsShift5+tsWidth*2,n5);

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
