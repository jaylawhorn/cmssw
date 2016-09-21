#ifndef HcalTesting_Mahi_Mahi_hh
#define HcalTesting_Mahi_Mahi_hh

#include <memory>
#include <string>
#include <map>
#include <iostream>

#include "TTree.h"
#include "TFile.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSimpleRecAlgo.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseContainmentManager.h"
#include "CondFormats/HcalObjects/interface/HcalRecoParams.h"
#include "CondFormats/HcalObjects/interface/HcalRecoParam.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalDbASCIIIO.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"

#include "CalibFormats/HcalObjects/interface/HcalCoder.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"

#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes.h"

#include "HcalTesting/Mahi/interface/DoMahiAlgo.h"

class Mahi : public edm::EDAnalyzer {
 public:
  explicit Mahi(const edm::ParameterSet&);
  ~Mahi();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------                                                                                              

  void ClearVariables();  

  std::map<int, double> hitEnergySumMap_;
  HcalSimParameterMap simParameterMap_;

  edm::EDGetTokenT<HBHERecHitCollection> tok_hbhe_;
  edm::EDGetTokenT<HBHEDigiCollection> tok_hbhe_digi_;

  bool FillHBHE;
  bool _IsData;
  double TotalChargeThreshold;
  double TS4ChargeThreshold;
  edm::Service<TFileService> FileService;

  // Basic event coordinates                                                   
  long long RunNumber;
  long long EventNumber;
  long long LumiSection;
  long long Bunch;
  long long Orbit;
  long long Time;

  // HBHE rechits and digis
  int PulseCount;
  double Charge[10];
  double Pedestal[10];
  double Gain[10];
  double SimHitEnergy;
  int IEta;
  int IPhi;
  int Depth;

  double Charge_M2;
  double Time_M2;
  double Charge_Mahi;

  SampleVector amplitudes;
  BXVector _bxs;

  DoMahiAlgo MahiAlgo;

  TTree *OutputTree;
  const CaloGeometry *Geometry;

  //double PulseCovArray[10][10];
  //double NoiseCovArray[10][10];
  //double TS[10];

};

#endif
