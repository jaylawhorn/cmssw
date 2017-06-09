import FWCore.ParameterSet.Config as cms

from DQMOffline.Trigger.RazorMonitor_cfi import hltRazorMonitoring

# HLT_RsqMR300_Rsq0p09_MR200_v6
RsqMR300_Rsq0p09_MR200_RazorMonitoring = hltRazorMonitoring.clone()
RsqMR300_Rsq0p09_MR200_RazorMonitoring.FolderName = cms.string('HLT/SUSY/RsqMR300_Rsq0p09_MR200/')
RsqMR300_Rsq0p09_MR200_RazorMonitoring.numGenericTriggerEventPSet.hltPaths = cms.vstring("HLT_RsqMR300_Rsq0p09_MR200_v*")

susyHLTRazorMonitoring = cms.Sequence(
    RsqMR300_Rsq0p09_MR200_RazorMonitoring
)

