// -*- C++ -*-
//
// Package:    RecoLocalCalo/HcalRecAlgos
// Class:      HBHEReconstructionDebugger
// 
/**\class HBHEReconstructionDebugger HBHEReconstructionDebugger.cc RecoLocalCalo/HcalRecAlgos/plugins/HBHEReconstructionDebugger.cc

 Description: Tool to extract and store debugging information from various HBHE Reconstruction algorithms

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jay Mathew Lawhorn
//         Created:  Sat, 10 Feb 2018 10:02:38 GMT
//
//


// system include files
#include <utility>
#include <algorithm>
#include <memory>
#include <vector>
#include <iostream>
#include <fstream>

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

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

#include "RecoLocalCalo/HcalRecAlgos/interface/SimpleHBHEPhase1AlgoDebug.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseContainmentManager.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/PulseShapeFitOOTPileupCorrection.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalDeterministicFit.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/MahiFit.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/parseHBHEPhase1AlgoDescription.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/fetchHcalAlgoData.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class HBHEReconstructionDebugger : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit HBHEReconstructionDebugger(const edm::ParameterSet&);
      ~HBHEReconstructionDebugger();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  // Configuration parameters
//  std::string algoConfigClass_;
//  bool processQIE8_;
//  bool processQIE11_;
//  bool saveInfos_;
//  bool saveDroppedInfos_;
//  bool makeRecHits_;
//  bool dropZSmarkedPassed_;
//  bool tsFromDB_;
//  bool recoParamsFromDB_;
//  bool saveEffectivePedestal_;
//  int sipmQTSShift_;
//  int sipmQNTStoSum_;
//
//  // Parameters for turning status bit setters on/off
//  bool setNegativeFlagsQIE8_;
//  bool setNegativeFlagsQIE11_;
//  bool setNoiseFlagsQIE8_;
//  bool setNoiseFlagsQIE11_;
//  bool setPulseShapeFlagsQIE8_;
//  bool setPulseShapeFlagsQIE11_;
//
//  // Other members
//  edm::EDGetTokenT<HBHEDigiCollection> tok_qie8_;
//  edm::EDGetTokenT<QIE11DigiCollection> tok_qie11_;
//  std::unique_ptr<AbsHBHEPhase1Algo> reco_;
//  std::unique_ptr<AbsHcalAlgoData> recoConfig_;
//  std::unique_ptr<HcalRecoParams> paramTS_;

  std::unique_ptr<AbsHBHEPhase1Algo> reco_;
  edm::EDGetTokenT<HBHEChannelInfoCollection> token_ChannelInfo_;
  
  edm::Service<TFileService> FileService;
  TTree *outTree;
  
  //std::unique_ptr<const HcalTimeSlew> hcalTimeSlewDelay_;

  int ieta;
  int iphi;
  int depth;

  int   nSamples;
  int   soi;

  bool  use3;

  float mahiEnergy;//SOI charge
  float chiSq;
  float arrivalTime;

  float pEnergy; //SOI-1 charge
  float nEnergy; //SOI+1 charge
  float pedEnergy; //pedestal charge

  float count[10]; //TS value 0-9
  float inputTS[10];//input TS samples
  float itPulse[10];//SOI pulse shape
  float pPulse[10];//SOI-1 pulse shape
  float nPulse[10];//SOI+1 pulse shape


};

HBHEReconstructionDebugger::HBHEReconstructionDebugger(const edm::ParameterSet& iConfig)
  : reco_(parseHBHEPhase1AlgoDescription(iConfig.getParameter<edm::ParameterSet>("algorithm")))
{

  if (!reco_.get())
    throw cms::Exception("HBHEPhase1BadConfig")
      << "Invalid HBHEPhase1Algo algorithm configuration"
      << std::endl;
  
  usesResource("TFileService");
  token_ChannelInfo_ = consumes<HBHEChannelInfoCollection>(edm::InputTag("hbheprereco",""));
}


HBHEReconstructionDebugger::~HBHEReconstructionDebugger()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void HBHEReconstructionDebugger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //edm::ESHandle<HcalTimeSlew> delay;
   //iSetup.get<HcalTimeSlewRecord>().get("HBHE", delay);
   //hcalTimeSlewDelay_ = std::unique_ptr<const HcalTimeSlew>(&(*delay));

   //std::cout << hcalTimeSlewDelay_->delay(100, HcalTimeSlew::Medium) << std::endl;

   Handle<HBHEChannelInfoCollection> hChannelInfo;
   iEvent.getByToken(token_ChannelInfo_, hChannelInfo);

   for (HBHEChannelInfoCollection::const_iterator iter = hChannelInfo->begin(); 
	iter != hChannelInfo->end(); iter++) {

     const HBHEChannelInfo& hci(*iter);
     const HcalDetId detid=hci.id();

     ieta  = detid.ieta();
     iphi  = detid.iphi();
     depth = detid.depth();

     const bool isRealData = true;
     MahiDebugInfo mdi = static_cast<SimpleHBHEPhase1AlgoDebug*>(reco_.get())->recoDebug(hci, isRealData);
     nSamples = mdi.nSamples;
     soi = mdi.soi;
     use3 = mdi.use3;
     mahiEnergy = mdi.mahiEnergy;
     chiSq = mdi.chiSq;
     arrivalTime = mdi.arrivalTime;
     pEnergy=mdi.pEnergy;
     nEnergy=mdi.nEnergy;
     pedEnergy=mdi.pedEnergy;
     for (int i=0; i<nSamples; i++) {
       count[i]=mdi.count[i];
       inputTS[i]=mdi.inputTS[i];
       itPulse[i]=mdi.itPulse[i];
       pPulse[i]=mdi.pPulse[i];
       nPulse[i]=mdi.nPulse[i];
     }
     if (nSamples==8) {
       count[8]=8;
       count[9]=9;
     }

     outTree->Fill();
   }

   //hcalTimeSlewDelay_.reset(nullptr);

}


// ------------ method called once each job just before starting event loop  ------------
void 
HBHEReconstructionDebugger::beginJob()
{

  outTree = FileService->make<TTree>("HcalTree","HcalTree");
  
  outTree->Branch("ieta",  &ieta,  "ieta/I");
  outTree->Branch("iphi",  &iphi,  "iphi/I");
  outTree->Branch("depth", &depth, "depth/I");
  outTree->Branch("nSamples",   &nSamples,   "nSamples/I");
  outTree->Branch("soi",   &soi,   "soi/I");

  outTree->Branch("mahiEnergy",   &mahiEnergy,   "mahiEnergy/F");
  outTree->Branch("chiSq",   &chiSq,   "chiSq/F");
  outTree->Branch("arrivalTime",   &arrivalTime,   "arrivalTime/F");
  outTree->Branch("pEnergy",   &pEnergy,   "pEnergy/F");
  outTree->Branch("nEnergy",   &nEnergy,   "nEnergy/F");
  outTree->Branch("pedEnergy",   &pedEnergy,   "pedEnergy/F");
  outTree->Branch("count",   &count,   "count[10]/F");
  outTree->Branch("inputTS",   &inputTS,   "inputTS[10]/F");
  outTree->Branch("itPulse",   &itPulse,   "itPulse[10]/F");
  outTree->Branch("pPulse",   &pPulse,   "pPulse[10]/F");
  outTree->Branch("nPulse",   &nPulse,   "nPulse[10]/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HBHEReconstructionDebugger::endJob() 
{
}

#define add_param_set(name)	     \
  edm::ParameterSetDescription name; \
  name.setAllowAnything();           \
  desc.add<edm::ParameterSetDescription>(#name, name)


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HBHEReconstructionDebugger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  add_param_set(algorithm);
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(HBHEReconstructionDebugger);
