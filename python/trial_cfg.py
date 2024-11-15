import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

# Add an empty source to allow cmsRun to run without actual event data
process.source = cms.Source("EmptySource")

process.trial = cms.EDProducer("trial",
    configString=cms.string("This is my configuration string")
)

process.p = cms.Path(process.trial)
