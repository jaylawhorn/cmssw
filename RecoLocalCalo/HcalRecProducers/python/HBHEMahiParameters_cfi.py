import FWCore.ParameterSet.Config as cms

# Configuration parameters for Method 2
mahiParameters = cms.PSet(

    doDebugMahi      = cms.int32(0),
    doPrefit         = cms.bool(True),
    floatPedestal    = cms.bool(True),
    activeBXs        = cms.vint32( -1, 0, 1),
    #activeBXs        = cms.vint32(0),
    nMaxIters        = cms.int32(500),

    ###                                      units
    #meanTime              = cms.double(0.),   #ns 
    #timeSigmaHPD          = cms.double(5.),   #ns 
    #timeSigmaSiPM         = cms.double(2.5),  #ns

    ###
    #applyTimeSlew    = cms.bool(True),
)
