import FWCore.ParameterSet.Config as cms

# Configuration parameters for MAHI
mahiParameters = cms.PSet(

    doPrefit         = cms.bool(True),
    floatPedestal    = cms.bool(False),
    #activeBXs        = cms.vint32( -4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
    activeBXs        = cms.vint32(-1, 0, 1),
    #activeBXs        = cms.vint32(0),
    nMaxIters        = cms.int32(500),

    ###                                      units
    #meanTime              = cms.double(0.),   #ns 
    #timeSigmaHPD          = cms.double(5.),   #ns 
    #timeSigmaSiPM         = cms.double(2.5),  #ns

    ###
    #applyTimeSlew    = cms.bool(True),
)
