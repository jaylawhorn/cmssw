#ifndef RecoLocalCalo_HcalRecAlgos_DoMahiAlgo_HH
#define RecoLocalCalo_HcalRecAlgos_DoMahiAlgo_HH

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/EigenMatrixTypes.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"


class DoMahiAlgo
{
 public:
  
  DoMahiAlgo();
  ~DoMahiAlgo() { };

  //  void phase1Apply(const HBHEChannelInfo& channelData, float& reconstructedEnergy, float& chi2);
  //  void setPulseShapeTemplate();

 private:
  
}; 
#endif
