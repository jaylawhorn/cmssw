// -*- C++ -*-
//
// Package:    HcalTesting/Mahi
// Class:      Mahi
// 
/**\class Mahi Mahi.cc HcalTesting/Mahi/plugins/Mahi.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jay Mathew Lawhorn
//         Created:  Sun, 04 Sep 2016 13:49:49 GMT
//
//


#include "HcalTesting/Mahi/plugins/Mahi.h"

// system include files
#include <string>
#include <map>
#include <iostream>
#include <memory>
using namespace std;


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSimpleRecAlgo.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseContainmentManager.h"
#include "CondFormats/HcalObjects/interface/HcalRecoParams.h"
#include "CondFormats/HcalObjects/interface/HcalRecoParam.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalDbASCIIIO.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
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

//
// constructors and destructor
//
Mahi::Mahi(const edm::ParameterSet& iConfig)
{

  FillHBHE = iConfig.getUntrackedParameter<bool>("FillHBHE", true);
  TotalChargeThreshold = iConfig.getUntrackedParameter<double>("TotalChargeThreshold", 100);
  TS4ChargeThreshold = iConfig.getUntrackedParameter<double>("TS4ChargeThreshold", 20);

  tok_hbhe_ = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheInput"));
  tok_hbhe_digi_ = consumes<HBHEDigiCollection>(edm::InputTag("hcalDigis",""));

  _IsData = iConfig.getUntrackedParameter<bool>("IsData");

}


Mahi::~Mahi()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Mahi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   ClearVariables();

   Handle<HBHERecHitCollection> hRecHits;
   iEvent.getByToken(tok_hbhe_, hRecHits);

   Handle<HBHEDigiCollection> hHBHEDigis;
   iEvent.getByToken(tok_hbhe_digi_, hHBHEDigis);

   Handle<PCaloHitContainer> hSimHits;
   iEvent.getByLabel("g4SimHits","HcalHits",hSimHits);
   
   ESHandle<HcalDbService> hConditions;
   iSetup.get<HcalDbRecord>().get(hConditions);

   ESHandle<CaloGeometry> hGeometry;
   iSetup.get<CaloGeometryRecord>().get(hGeometry);
   Geometry = hGeometry.product();

   RunNumber = iEvent.id().run();
   EventNumber = iEvent.id().event();
   LumiSection = iEvent.luminosityBlock();
   Bunch = iEvent.bunchCrossing();
   Orbit = iEvent.orbitNumber();
   Time = iEvent.time().value();

   //ECAL does it classier 
   // https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit.cc#L151-L171
   _bxs.resize(3);
   _bxs << -1,0,1;

   //std::cout << RunNumber << ", " << EventNumber << ", " << LumiSection << std::endl;

   if(!_IsData) {
     // store the energy of each hit in a map
     hitEnergySumMap_.clear();
     PCaloHitContainer::const_iterator hitItr = hSimHits->begin();
     PCaloHitContainer::const_iterator last = hSimHits->end();
     
     for( ; hitItr != last; ++hitItr)
       {
	 HcalDetId hcalDetId(hitItr->id());
	  
	 if(hcalDetId.subdet()== HcalBarrel || hcalDetId.subdet() == HcalEndcap)
	   {
	     int id = hitItr->id();
	     double samplingFactor=1.;
	     if(hcalDetId.subdet()== HcalBarrel)
	       {
		 samplingFactor = simParameterMap_.hbParameters().samplingFactor(DetId(id));
	       }
	     else if(hcalDetId.subdet() == HcalEndcap)
	       {
		 samplingFactor = simParameterMap_.heParameters().samplingFactor(DetId(id));
	       }
	     double energy = hitItr->energy() * samplingFactor;
	     // add it to the map
	     std::map<int, double>::iterator mapItr = hitEnergySumMap_.find(id);
	     if(mapItr == hitEnergySumMap_.end()) {
	       hitEnergySumMap_[id] = energy;
	     }
	     else
	       {
		 hitEnergySumMap_[id] += energy;
	       }
	   }
       }
   }

   map<HcalDetId, int> RecHitIndex;

   for(int i = 0; i < (int)hRecHits->size(); i++)
     {
       HcalDetId id = (*hRecHits)[i].id();
       RecHitIndex.insert(pair<HcalDetId, int>(id, i));
     }

   // loop over digis
   for(HBHEDigiCollection::const_iterator iter = hHBHEDigis->begin(); iter != hHBHEDigis->end(); iter++)
     {
       HcalDetId id = iter->id();

       // First let's convert ADC to deposited charge
       const HcalCalibrations &Calibrations = hConditions->getHcalCalibrations(id);
       const HcalQIECoder *ChannelCoder = hConditions->getHcalCoder(id);
       const HcalQIEShape *Shape = hConditions->getHcalShape(ChannelCoder);
       HcalCoderDb Coder(*ChannelCoder, *Shape);
       CaloSamples Tool;
       Coder.adc2fC(*iter, Tool);

       std::vector<int> capidvec;
       for(int ip=0; ip<Tool.size(); ip++){
	 const int capid = iter->sample(ip).capid();
	 capidvec.push_back(capid);
       }

       std::vector<double> correctedOutput;
       MahiAlgo.Apply(Tool, capidvec, id, Calibrations, correctedOutput);

       for (uint i=0; i<correctedOutput.size(); i++) {
	 std::cout << correctedOutput.at(i) << ", ";
       }
       std::cout << std::endl;
     }

   
}


// ------------ method called once each job just before starting event loop  ------------
void 
Mahi::beginJob()
{

  MahiAlgo = DoMahiAlgo();

  OutputTree = FileService->make<TTree>("HcalTree", "Hcal tree");
  
  OutputTree->Branch("RunNumber", &RunNumber, "RunNumber/LL");
  OutputTree->Branch("EventNumber", &EventNumber, "EventNumber/LL");
  OutputTree->Branch("LumiSection", &LumiSection, "LumiSection/LL");
  OutputTree->Branch("Bunch", &Bunch, "Bunch/LL");
  OutputTree->Branch("Orbit", &Orbit, "Orbit/LL");
  OutputTree->Branch("Time", &Time, "Time/LL");
  
  if(FillHBHE == true) {
    OutputTree->Branch("Charge", &Charge, "Charge[10]/D");
    OutputTree->Branch("Pedestal", &Pedestal, "Pedestal[10]/D");
    OutputTree->Branch("Gain", &Gain, "Gain[10]/D");
    
    if(!_IsData) OutputTree->Branch("SimHitEnergy", &SimHitEnergy, "SimHitEnergy/D");
    
    OutputTree->Branch("IEta", &IEta, "IEta/I");
    OutputTree->Branch("IPhi", &IPhi, "IPhi/I");
    OutputTree->Branch("Depth", &Depth, "Depth/I");

    OutputTree->Branch("Charge_M2", &Charge_M2, "Charge_M2/D");
    OutputTree->Branch("Time_M2", &Time_M2, "Time_M2/D");
    OutputTree->Branch("Charge_Mahi", &Charge_Mahi, "Charge_Mahi/D");
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Mahi::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Mahi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
Mahi::ClearVariables() {

  RunNumber = 0;
  EventNumber = 0;
  LumiSection = 0;
  Bunch = 0;
  Orbit = 0;
  Time = 0;

  for (int i=0; i<10; i++) {
    Gain[i] = 0;
    Charge[i] = 0;
    Pedestal[i] = 0;
  }

  if(!_IsData) SimHitEnergy = 0;

  IEta = 0;
  IPhi = 0;
  Depth = 0;
  
  Charge_M2 = 0;
  Time_M2 = 0;
  Charge_Mahi = 0;

}

//define this as a plug-in
DEFINE_FWK_MODULE(Mahi);
