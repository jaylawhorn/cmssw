import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/work/j/jlawhorn/public/method2rnd2016/HIN_RAW.root'
        '/store/data/Run2016G/HLTPhysics0/RAW/v1/000/280/002/00000/580EACE6-B572-E611-99D2-02163E014733.root'
        #'/store/data/Run2016G/HLTPhysics0/RAW/v1/000/280/002/00000/10B719D7-B572-E611-991F-FA163E48A31D.root'
        )
                            )

from Configuration.StandardSequences.RawToDigi_cff import *
process.RawToDigi = cms.Sequence(hcalDigis)

from RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi import *
process.hcalOOTPileupESProducer = cms.ESProducer('OOTPileupDBCompatibilityESProducer')

from RecoLocalCalo.HcalRecProducers.HcalHitReconstructor_hbhe_cfi import *
hbheprereco.setNegativeFlags=False

process.hcalLocalRecoSequence = cms.Sequence(hbheprereco)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("HCALTree.root")
                                   )

process.demo = cms.EDAnalyzer('Mahi',
                              hbheInput = cms.InputTag('hbheprereco'), 
                              IsData = cms.untracked.bool(True),
                              TotalChargeThreshold = cms.untracked.double(-999),
                              TS4ChargeThreshold = cms.untracked.double(20)
                              )

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v8', '')

process.raw2digi_step = cms.Path(process.RawToDigi)
process.reco_step1 = cms.Path(process.hcalLocalRecoSequence)
process.reco_step2 = cms.Path(process.demo)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.schedule = cms.Schedule(process.raw2digi_step,process.reco_step1,process.reco_step2,process.endjob_step)
