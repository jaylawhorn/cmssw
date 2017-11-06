import FWCore.ParameterSet.Config as cms

# Configuration parameters for MAHI
mahiParameters = cms.PSet(

    activeBXs        = cms.vint32( -4, -3, -2, -1, 0, 1, 2, 3, 4),
    #activeBXs        = cms.vint32(-1, 0, 1),
    #activeBXs        = cms.vint32(0),
    nMaxIters        = cms.int32(500),

)
