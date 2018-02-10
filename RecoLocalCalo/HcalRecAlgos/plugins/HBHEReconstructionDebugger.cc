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

  std::string algoConfigClass_;
  bool processQIE8_;
  bool processQIE11_;
  bool saveInfos_;
  bool saveDroppedInfos_;
  bool makeRecHits_;
  bool dropZSmarkedPassed_;
  bool tsFromDB_;
  bool recoParamsFromDB_;
  bool saveEffectivePedestal_;
  int sipmQTSShift_;
  int sipmQNTStoSum_;

  bool setNegativeFlagsQIE8_;
  bool setNegativeFlagsQIE11_;
  bool setNoiseFlagsQIE8_;
  bool setNoiseFlagsQIE11_;
  bool setPulseShapeFlagsQIE8_;
  bool setPulseShapeFlagsQIE11_;

  
  //std::map<int, double> hitEnergySumMap_;
  //HcalSimParameterMap simParameterMap_;

  edm::EDGetTokenT<HBHEChannelInfoCollection> token_ChannelInfo_;
  //edm::EDGetTokenT<HBHERecHitCollection> token_RecHit_;
  //edm::EDGetTokenT<edm::PCaloHitContainer> tok_hbhe_sim_;

  edm::Service<TFileService> FileService;

  std::unique_ptr<AbsHBHEPhase1Algo> reco_;

  //const HcalDDDRecConstants *hcons;
  //const CaloGeometry *Geometry;

  //TTree *outTree;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HBHEReconstructionDebugger::HBHEReconstructionDebugger(const edm::ParameterSet& iConfig)
  :  //algoConfigClass_(iConfig.getParameter<std::string>("algoConfigClass")),
     //processQIE8_(iConfig.getParameter<bool>("processQIE8")),
     //processQIE11_(iConfig.getParameter<bool>("processQIE11")),
     //saveInfos_(iConfig.getParameter<bool>("saveInfos")),
     //saveDroppedInfos_(iConfig.getParameter<bool>("saveDroppedInfos")),
     //makeRecHits_(iConfig.getParameter<bool>("makeRecHits")),
     //dropZSmarkedPassed_(iConfig.getParameter<bool>("dropZSmarkedPassed")),
     //tsFromDB_(iConfig.getParameter<bool>("tsFromDB")),
     //recoParamsFromDB_(iConfig.getParameter<bool>("recoParamsFromDB")),
     //saveEffectivePedestal_(iConfig.getParameter<bool>("saveEffectivePedestal")),
     //sipmQTSShift_(iConfig.getParameter<int>("sipmQTSShift")),
     //sipmQNTStoSum_(iConfig.getParameter<int>("sipmQNTStoSum")),
     //setNegativeFlagsQIE8_(iConfig.getParameter<bool>("setNegativeFlagsQIE8")),
     //setNegativeFlagsQIE11_(iConfig.getParameter<bool>("setNegativeFlagsQIE11")),
     //setNoiseFlagsQIE8_(iConfig.getParameter<bool>("setNoiseFlagsQIE8")),
     //setNoiseFlagsQIE11_(iConfig.getParameter<bool>("setNoiseFlagsQIE11")),
     //setPulseShapeFlagsQIE8_(iConfig.getParameter<bool>("setPulseShapeFlagsQIE8")),
     //setPulseShapeFlagsQIE11_(iConfig.getParameter<bool>("setPulseShapeFlagsQIE11")),
     reco_(parseHBHEPhase1AlgoDescription(iConfig.getParameter<edm::ParameterSet>("algorithm")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

   token_ChannelInfo_ = consumes<HBHEChannelInfoCollection>(edm::InputTag("hbheprereco",""));
   //token_RecHit_ = consumes<HBHERecHitCollection>(edm::InputTag("hbheprereco",""));
   //tok_hbhe_sim_ = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));

}


HBHEReconstructionDebugger::~HBHEReconstructionDebugger()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HBHEReconstructionDebugger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "hi?" << std::endl;

   using namespace edm;

   //ESHandle<HcalDDDRecConstants> pHRNDC;
   //iSetup.get<HcalRecNumberingRecord>().get( pHRNDC );
   //hcons = &(*pHRNDC);

   Handle<HBHEChannelInfoCollection> hChannelInfo;
   iEvent.getByToken(token_ChannelInfo_, hChannelInfo);

   for (HBHEChannelInfoCollection::const_iterator iter = hChannelInfo->begin(); 
	iter != hChannelInfo->end(); iter++) {

     const HBHEChannelInfo& hci(*iter);
     const bool isRealData = true;

     MahiDebugInfo mdi = static_cast<SimpleHBHEPhase1AlgoDebug*>(reco_.get())->recoDebug(hci, isRealData);
     std::cout << mdi.soi << std::endl;

   }

   //Handle<HBHERecHitCollection> hRecHit;
   //iEvent.getByToken(token_RecHit_, hRecHit);
   //
   //Handle<PCaloHitContainer> hSimHits;
   //iEvent.getByToken(tok_hbhe_sim_,hSimHits);
   //
   //ESHandle<HcalDbService> hConditions;
   //iSetup.get<HcalDbRecord>().get(hConditions);
   //
   //ESHandle<CaloGeometry> hGeometry;
   //iSetup.get<CaloGeometryRecord>().get(hGeometry);
   //Geometry = hGeometry.product();

}


// ------------ method called once each job just before starting event loop  ------------
void 
HBHEReconstructionDebugger::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HBHEReconstructionDebugger::endJob() 
{
}

#define add_param_set(name) /**/       \
  edm::ParameterSetDescription name; \
  name.setAllowAnything();           \
  desc.add<edm::ParameterSetDescription>(#name, name)


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HBHEReconstructionDebugger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;

  //desc.add<edm::InputTag>("digiLabelQIE8");
  //desc.add<edm::InputTag>("digiLabelQIE11");
  //desc.add<std::string>("algoConfigClass");
  //desc.add<bool>("processQIE8");
  //desc.add<bool>("processQIE11");
  //desc.add<bool>("saveInfos");
  //desc.add<bool>("saveDroppedInfos");
  //desc.add<bool>("makeRecHits");
  //desc.add<bool>("dropZSmarkedPassed");
  //desc.add<bool>("tsFromDB");
  //desc.add<bool>("recoParamsFromDB");
  //desc.add<bool>("saveEffectivePedestal", false);
  //desc.add<int>("sipmQTSShift", 0);
  //desc.add<int>("sipmQNTStoSum", 3);
  //desc.add<bool>("setNegativeFlagsQIE8");
  //desc.add<bool>("setNegativeFlagsQIE11");
  //desc.add<bool>("setNoiseFlagsQIE8");
  //desc.add<bool>("setNoiseFlagsQIE11");
  //desc.add<bool>("setPulseShapeFlagsQIE8");
  //desc.add<bool>("setPulseShapeFlagsQIE11");
  //desc.add<bool>("setLegacyFlagsQIE8");
  //desc.add<bool>("setLegacyFlagsQIE11");

  add_param_set(algorithm);
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(HBHEReconstructionDebugger);
