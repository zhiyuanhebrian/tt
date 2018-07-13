import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("USER")

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        
        'file:/nfs/dust/cms/user/lbenato/calo_jets/miniaod_calo.root'
    )
)

isData = ('/store/data/' in process.source.fileNames[0])

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag (process.GlobalTag, 'auto:run2_mc')#default option, but we have the best global tag manually

GT = ''
if isData:
    GT = '80X_dataRun2_2016SeptRepro_v7'
    print "data 2016, ReMiniaod GT"
elif not(isData):
    GT = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'#Moriond17 GT

process.GlobalTag = GlobalTag (process.GlobalTag, GT)

process.TFileService = cms.Service( "TFileService",
    fileName = cms.string('output_calo.root'),# if len(options.outputFile)==0 else options.outputFile),
    closeFileFast = cms.untracked.bool(True),
)

#################################################
## Remake jets
#################################################

fatjet_ptmin = 50.

##Dummy copy for having pf candidates into the process
#process.packed = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))


## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
## Define GenJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')

from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.ak8GenJetsNoNu = ak5GenJets.clone(src = 'packedGenParticlesForJetsNoNu', rParam = 0.8)

## Standard AK4 and AK8 jets without CHS
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4PFJets = ak4PFJets.clone(src = 'packedPFCandidates', doAreaFastjet = True)#HERE!!! add jetPtMin = (fatjet_ptmin)
process.ak8PFJets  = ak4PFJets.clone (src = 'packedPFCandidates', rParam = 0.8, doAreaFastjet = True, jetPtMin = fatjet_ptmin)

## Pruned AK8 without CHS
from RecoJets.JetProducers.ak4PFJetsPruned_cfi import ak4PFJetsPruned
process.ak8PFJetsPruned = ak4PFJetsPruned.clone(rParam = 0.8, doAreaFastjet = True, src = 'packedPFCandidates', jetPtMin = fatjet_ptmin)

## Softdrop AK8 without CHS
from RecoJets.JetProducers.ak4PFJetsSoftDrop_cfi import ak4PFJetsSoftDrop
process.ak8PFJetsSoftDrop = ak4PFJetsSoftDrop.clone(R0 = 0.8, rParam = 0.8, doAreaFastjet = True, src = 'packedPFCandidates', jetPtMin = fatjet_ptmin)

## PUPPI AK8
## Pre-requisites: re-build puppi jets

###Robin's way:

process.load('CommonTools/PileupAlgos/Puppi_cff')
## e.g. to run on miniAOD
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
process.puppi.clonePackedCands   = cms.bool(True)
process.puppi.useExistingWeights = cms.bool(True)

from RecoJets.JetProducers.ak4PFJetsPuppi_cfi import ak4PFJetsPuppi
process.ak8PFJetsPuppi = ak4PFJetsPuppi.clone(rParam = 0.8, doAreaFastjet = True, jetPtMin = fatjet_ptmin)

## Softdrop AK8 with Puppi
from RecoJets.Configuration.RecoPFJets_cff import ak8PFJetsPuppiSoftDrop
process.ak8PFJetsPuppiSoftDrop = ak8PFJetsPuppiSoftDrop.clone(jetPtMin = fatjet_ptmin)

###PERRCHE NON CONOSCE PIU' PACKED PF CANDIDATES??????????????
##print process.packedPFCandidates
#################
### Tentative n.2: association track-vtx



#################################################
## Remake PAT jets without CHS
#################################################

## b-tag discriminators
bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags'
]

from PhysicsTools.PatAlgos.tools.jetTools import *
## Add PAT jet collection based on the above-defined ak4PFJetsCHS
addJetCollection(
    process,
    labelName = 'AK4PF',
    jetSource = cms.InputTag('ak4PFJets'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak4GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    algo = 'AK',
    rParam = 0.4
)

addJetCollection(
    process,
    labelName = 'AK8PF',
    jetSource = cms.InputTag('ak8PFJets'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK8PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    algo = 'AK',
    rParam = 0.8
)

addJetCollection(
    process,
    labelName = 'AK8Puppi',
    jetSource = cms.InputTag('ak8PFJetsPuppi'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK8PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    algo = 'AK',
    rParam = 0.8
)


##
## PF AK8 groomed masses without CHS
process.patJetsAK8PF.userData.userFloats.src = []
from RecoJets.JetProducers.ak8PFJetsCHS_groomingValueMaps_cfi import ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass
process.ak8PFJetsSoftDropMass = ak8PFJetsCHSSoftDropMass.clone(src = 'ak8PFJets', matched = 'ak8PFJetsSoftDrop')
process.ak8PFJetsPrunedMass = ak8PFJetsCHSPrunedMass.clone(src = 'ak8PFJets', matched = 'ak8PFJetsPruned')
process.patJetsAK8PF.userData.userFloats.src += ['ak8PFJetsPrunedMass','ak8PFJetsSoftDropMass']

##
## PF AK8 Puppi groomed masses with Puppi
process.patJetsAK8Puppi.userData.userFloats.src = []
process.load("RecoJets.JetProducers.ak8PFJetsPuppi_groomingValueMaps_cfi")
process.patJetsAK8Puppi.userData.userFloats.src += ['ak8PFJetsPuppiSoftDropMass']

##
## N-subjettiness without CHS
from RecoJets.JetProducers.nJettinessAdder_cfi import *
process.NjettinessAK8 = Njettiness.clone(src='ak8PFJets')#src='ak8PFJets', cone=0.8)
process.NjettinessAK8.cone = cms.double(0.8)
process.patJetsAK8PF.userData.userFloats.src += ['NjettinessAK8:tau1','NjettinessAK8:tau2','NjettinessAK8:tau3']

##
## Puppi N-subjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import *
process.NjettinessAK8Puppi = Njettiness.clone(src='ak8PFJetsPuppi')
process.NjettinessAK8Puppi.cone = cms.double(0.8)
process.patJetsAK8Puppi.userData.userFloats.src += ['NjettinessAK8Puppi:tau1','NjettinessAK8Puppi:tau2','NjettinessAK8Puppi:tau3']


##THESE VALUE MAPS ARE PROBABLY NOT NEEDED
## PF AK8 matching to PF Puppi AK8
process.ak8PFJetsPuppiValueMap = cms.EDProducer("RecoJetToPatJetDeltaRValueMapProducer",
                                                src = cms.InputTag("ak8PFJets"),
                                                matched = cms.InputTag("patJetsAK8Puppi"),                                         
                                                distMax = cms.double(0.8),
                                                values = cms.vstring([
            'userFloat("NjettinessAK8Puppi:tau1")',
            'userFloat("NjettinessAK8Puppi:tau2")',
            'userFloat("NjettinessAK8Puppi:tau3")',
            'pt','eta','phi','mass'
            ]),
                                                valueLabels = cms.vstring( [
            'NjettinessAK8PuppiTau1',
            'NjettinessAK8PuppiTau2',
            'NjettinessAK8PuppiTau3',
            'pt','eta','phi','mass'
            ])
)


process.patJetsAK8PF.userData.userFloats.src += [cms.InputTag('ak8PFJetsPuppiValueMap','NjettinessAK8PuppiTau1'),
                                               cms.InputTag('ak8PFJetsPuppiValueMap','NjettinessAK8PuppiTau2'),
                                               cms.InputTag('ak8PFJetsPuppiValueMap','NjettinessAK8PuppiTau3'),
                                               cms.InputTag('ak8PFJetsPuppiValueMap','pt'),
                                               cms.InputTag('ak8PFJetsPuppiValueMap','eta'),
                                               cms.InputTag('ak8PFJetsPuppiValueMap','phi'),
                                               cms.InputTag('ak8PFJetsPuppiValueMap','mass'),
                                               ]

######



### here we try to have selectedPatJetsAK8PFCHSSoftDrop, slimmedJetsAK8PFCHSSoftDropSubjets

###LATER

###################################################

### here we need to pattify AK8PFSoftDrop jets

addJetCollection(
        process,
        labelName = 'AK8PFSoftDrop',
        jetSource = cms.InputTag('ak8PFJetsSoftDrop'),
        pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        btagDiscriminators = ['None'],
        jetCorrections = ('AK8PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
        getJetMCFlavour = False, # jet flavor disabled
        genParticles = cms.InputTag('prunedGenParticles'),
        genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
        algo = 'AK',
        rParam = 0.8
)

###try subjets, this will probably fail
addJetCollection(
        process,
        labelName = 'AK8PFSoftDropSubjets',
        jetSource = cms.InputTag('ak8PFJetsSoftDrop','SubJets'),
        algo = 'ak',  # needed for subjet flavor clustering
        rParam = 0.8, # needed for subjet flavor clustering
        btagDiscriminators = ['pfCombinedSecondaryVertexV2BJetTags', 'pfCombinedInclusiveSecondaryVertexV2BJetTags'],
        jetCorrections = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
        explicitJTA = True,  # needed for subjet b tagging
        svClustering = True, # needed for subjet b tagging
        pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        genParticles = cms.InputTag('prunedGenParticles'),
        genJetCollection = cms.InputTag('ak4GenJetsNoNu'), 
        fatJets=cms.InputTag('ak8PFJets'),             # needed for subjet flavor clustering
        groomedFatJets=cms.InputTag('ak8PFJetsSoftDrop') # needed for subjet flavor clustering
)

### here we need to pattify AK8PFPuppiSoftDrop jets
addJetCollection(
        process,
        labelName = 'AK8PFPuppiSoftDrop',
        jetSource = cms.InputTag('ak8PFJetsPuppiSoftDrop'),
        pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        btagDiscriminators = ['None'],
        #btagDiscriminators = bTagDiscriminators,
        jetCorrections = ('AK8PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
        getJetMCFlavour = False, # jet flavor disabled
        genParticles = cms.InputTag('prunedGenParticles'),
        genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
        algo = 'AK',
        rParam = 0.8
)

### here we need to pattify AK8PFPuppiSoftDropSubjets, that are AK4
addJetCollection(
        process,
        labelName = 'AK8PFPuppiSoftDropSubjets',
        jetSource = cms.InputTag('ak8PFJetsPuppiSoftDrop','SubJets'),
        algo = 'ak',  # needed for subjet flavor clustering
        rParam = 0.8, # needed for subjet flavor clustering
        btagDiscriminators = ['pfCombinedSecondaryVertexV2BJetTags', 'pfCombinedInclusiveSecondaryVertexV2BJetTags'],
        jetCorrections = ('AK4PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
        explicitJTA = True,  # needed for subjet b tagging
        svClustering = True, # needed for subjet b tagging
        pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        genParticles = cms.InputTag('prunedGenParticles'),
        genJetCollection = cms.InputTag('ak4GenJetsNoNu'), 
        fatJets=cms.InputTag('ak8PFJetsPuppi'),             # needed for subjet flavor clustering
        groomedFatJets=cms.InputTag('ak8PFJetsPuppiSoftDrop'), # needed for subjet flavor clustering
)


##
## here we try to have slimmedJetsAK8PFCHSSoftDropPacked
#process.slimmedJetsAK8PFCHSSoftDropPacked = cms.EDProducer("BoostedJetMerger",
#        jetSrc=cms.InputTag("selectedPatJetsAK8PFCHSSoftDrop"),
#        subjetSrc=cms.InputTag("slimmedJetsAK8PFCHSSoftDropSubjets")
#)

## here we try to have patJetsAK8PFPuppiSoftDropPacked and patJetsAK8PFPuppiSoftDropPacked

process.patJetsAK8PFSoftDropPacked = cms.EDProducer("BoostedJetMerger",
        jetSrc=cms.InputTag("patJetsAK8PFSoftDrop"),#selectedPat..
        subjetSrc=cms.InputTag("patJetsAK8PFSoftDropSubjets")#("slimmedJetsAK8PFPuppiSoftDropSubjets")
)

process.patJetsAK8PFPuppiSoftDropPacked = cms.EDProducer("BoostedJetMerger",
        jetSrc=cms.InputTag("patJetsAK8PFPuppiSoftDrop"),
        subjetSrc=cms.InputTag("patJetsAK8PFPuppiSoftDropSubjets")
)


#From miniaod dump ---> this should give all the proper subjets.

process.packedPatJetsAK8 = cms.EDProducer("JetSubstructurePacker",
    algoLabels = cms.vstring('SoftDrop',
                             'SoftDropPuppi'
                             ),
    algoTags = cms.VInputTag(cms.InputTag("patJetsAK8PFSoftDropPacked"), 
                             cms.InputTag("patJetsAK8PFPuppiSoftDropPacked")
                             ),#I need those, before
    #algoTags = cms.VInputTag(cms.InputTag("slimmedJetsAK8PFCHSSoftDropPacked"), cms.InputTag("slimmedJetsAK8PFPuppiSoftDropPacked")),#I need those, before
    distMax = cms.double(0.8),
    fixDaughters = cms.bool(False),
    jetSrc = cms.InputTag("patJetsAK8PF"),
    packedPFCandidates = cms.InputTag("packedPFCandidates"),
    #candidates = cms.InputTag("packedPFCandidates"),
)


##getattr(process,'selectedPatJetsAK4PFCHS').cut = cms.string('pt > 10')

process.seq = cms.Sequence(
    process.ak4PFJets *
    process.ak8PFJets
    #process.selectedPatJetsAK4PF
)


process.p = cms.Path(process.seq)

from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False), # while the timing of this is not reliable in unscheduled mode, it still helps understanding what was actually run
        allowUnscheduled = cms.untracked.bool(True)
)


#-----------------------#
#     DATA FLAGS        #
#-----------------------#
isData            = ('/store/data/' in process.source.fileNames[0])
isReHLT           = ('_reHLT_' in process.source.fileNames[0])
isReReco          = ('23Sep2016' in process.source.fileNames[0])
isReMiniAod       = ('03Feb2017' in process.source.fileNames[0])
isPromptReco      = ('PromptReco' in process.source.fileNames[0])
theRunBCD = ['Run2016B','Run2016C','Run2016D']
theRunEF  = ['Run2016E','Run2016F']
theRunG   = ['Run2016G']
theRunH   = ['Run2016H']

print 'isData',isData
print 'isReHLT',isReHLT
print 'isReReco',isReReco
print 'isReMiniAod',isReMiniAod
print 'isPromptReco',isPromptReco

#-----------------------#
#    VERTEX FILTER      #
#-----------------------#

import RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi
process.primaryVertexFilter = cms.EDFilter('GoodVertexFilter',
    vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),# if not isAOD else 'offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)

#-----------------------#
#     GLOBAL TAG        #
#-----------------------#

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
GT = ''
if isData:
    if isPromptReco: GT = "92X_dataRun2_Prompt_v4"
    print "data 2017, PromptReco"
#else:
#    GT = "90X_upgrade2017_realistic_v20"
elif not(isData):                                       
    GT = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'#Moriond17 GT

process.GlobalTag = GlobalTag(process.GlobalTag, GT)
print 'GlobalTag loaded: ', GT

#-----------------------#
#         JEC           #
#-----------------------#

JECstring = ''
if isData and (isReReco or isReMiniAod):
  if any(s in process.source.fileNames[0] for s in theRunBCD):
    JECstring = "Summer16_23Sep2016BCDV3_DATA" #if isReMiniAod else "Summer16_23Sep2016BCDV3_DATA"
  if any(s in process.source.fileNames[0] for s in theRunEF):
    JECstring = "Summer16_23Sep2016EFV3_DATA" #if isReMiniAod else "Summer16_23Sep2016EFV3_DATA"
  if any(s in process.source.fileNames[0] for s in theRunG):
    JECstring = "Summer16_23Sep2016GV3_DATA" #if isReMiniAod else "Summer16_23Sep2016GV3_DATA"
  if any(s in process.source.fileNames[0] for s in theRunH):
    JECstring = "Summer16_23Sep2016HV3_DATA" #if isReMiniAod else "Summer16_23Sep2016HV3_DATA"
elif isData and isPromptReco:
    JECstring = "Spring16_25nsV6_DATA"
elif not isData:
    JECstring = "Summer16_23Sep2016V3_MC"

print "JEC ->",JECstring

#-----------------------#
#        FILTERS        #
#-----------------------#

# JSON filter
if isData:
    import FWCore.PythonUtilities.LumiList as LumiList
    jsonName = "Cert_294927-305364_13TeV_PromptReco_Collisions17_JSON"#"Cert_294927-301567_13TeV_PromptReco_Collisions17_JSON" #golden json
    process.source.lumisToProcess = LumiList.LumiList(filename = 'data/JSON/'+jsonName+'.txt').getVLuminosityBlockRange()
    print "JSON file loaded: ", jsonName

# Trigger filter
triggerTag = 'HLT2' if isReHLT else 'HLT'

# MET filters
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag('slimmedMuons')# if not isAOD else 'muons')
process.BadPFMuonFilter.PFCandidates = cms.InputTag('packedPFCandidates')# if not isAOD else 'PFCandidates')

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag('slimmedMuons')# if not isAOD else 'muons')
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag('packedPFCandidates')# if not isAOD else 'PFCandidates')

if isData:
    filterString = "RECO"
else:
    filterString = "PAT"

#-----------------------#
#       ANALYZER        #
#-----------------------#

process.reconstruction = cms.EDAnalyzer('RecoStudiesCalo',
    genSet = cms.PSet(
        genProduct = cms.InputTag('generator'),
        lheProduct = cms.InputTag('externalLHEProducer'),
        genParticles = cms.InputTag('prunedGenParticles'),# if not isAOD else 'genParticles'),
        pdgId = cms.vint32(5,9000006),#(1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 21, 23, 24, 25, 36, 39, 1000022, 9100000, 9000001, 9000002, 9100012, 9100022, 9900032, 1023),
        status = cms.vint32(22,23),
        samplesDYJetsToLL = cms.vstring(),
        samplesZJetsToNuNu = cms.vstring(),
        samplesWJetsToLNu = cms.vstring(),
        samplesDir = cms.string('data/Stitch/'),
        sample = cms.string("" ), #( sample )
        ewkFile = cms.string('data/scalefactors_v4.root'),
        applyEWK = cms.bool(False),#(True if sample.startswith('DYJets') or sample.startswith('WJets') else False),
        applyTopPtReweigth = cms.bool(False),#(True if sample.startswith('TT_') else False),
        pythiaLOSample = cms.bool(False),#(True if isDibosonInclusive else False),
    ),
    pileupSet = cms.PSet(
        pileup = cms.InputTag('slimmedAddPileupInfo'),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
        dataFileName     = cms.string('data/PU_69200_ReReco.root'),#updated
        dataFileNameUp   = cms.string('data/PU_72380_ReReco.root'),#updated
        dataFileNameDown = cms.string('data/PU_66020_ReReco.root'),#updated
        mcFileName = cms.string('data/PU_MC_Moriond17.root'),#updated
        dataName = cms.string('pileup'),
        mcName = cms.string('2016_25ns_Moriond17MC_PoissonOOTPU'),#updated
    ),
    triggerSet = cms.PSet(
        trigger = cms.InputTag('TriggerResults', '', triggerTag),
        
        paths = cms.vstring(
*[
 'HLT_Ele27_WPTight_Gsf_v',  'HLT_Ele25_eta2p1_WPTight_Gsf_v',  'HLT_Ele27_eta2p1_WPTight_Gsf_v',  'HLT_Ele32_eta2p1_WPTight_Gsf_v',  'HLT_Ele27_WPLoose_Gsf_WHbbBoost_v',  'HLT_Ele115_CaloIdVT_GsfTrkIdT_v',  'HLT_IsoMu24_v',  'HLT_IsoMu22_eta2p1_v',  'HLT_IsoTkMu24_v',  'HLT_IsoMu27_v',  'HLT_IsoTkMu22_eta2p1_v',  'HLT_IsoTkMu27_v',  'HLT_Mu50_v',  'HLT_TkMu50_v',  'HLT_Mu55_v',  'HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v',  'HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v',  'HLT_Mu6_PFHT200_PFMET100_v',  'HLT_Mu15_IsoVVVL_PFHT400_v',  'HLT_Ele15_IsoVVVL_PFHT400_v',  'HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v',  'HLT_Mu50_IsoVVVL_PFHT400_v',  'HLT_Ele50_IsoVVVL_PFHT400_v',  'HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v',  'HLT_Mu30_eta2p1_PFJet150_PFJet50_v',  'HLT_DoubleMu3_PFMET50_v',  'HLT_Mu15_IsoVVVL_PFHT400_PFMET50_v',  'HLT_Ele15_IsoVVVL_PFHT400_PFMET50_v',  'HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v',  'HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose_v',  'HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v',  'HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v',  'HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v',  'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v',  'HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v',  'HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v',  'HLT_Ele27_WPTight_Gsf_L1JetTauSeeded_v',  'HLT_VLooseIsoPFTau140_Trk50_eta2p1_v',  'HLT_Mu30_TkMu11_v',  'HLT_Mu17_Mu8_SameSign_DZ_v',  'HLT_Mu40_TkMu11_v',  'HLT_Mu20_Mu10_SameSign_DZ_v',  'HLT_DoubleMu8_Mass8_PFHT300_v',  'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',  'HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v',  'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v',  'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v',  'HLT_AK8PFJet360_TrimMass30_v',  'HLT_AK8PFJet450_v',  'HLT_AK8PFJet500_v',  'HLT_BTagMu_AK8Jet300_Mu5_v',  'HLT_BTagMu_Jet300_Mu5_v',  'HLT_CaloJet500_NoJetID_v',  'HLT_DiCentralPFJet170_CFMax0p1_v',  'HLT_DiCentralPFJet330_CFMax0p5_v',  'HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v',  'HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v',  'HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6_v',  'HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172_v',  'HLT_DoubleMu3_PFMET50_v',  'HLT_DoubleMu8_Mass8_PFHT300_v',  'HLT_Ele115_CaloIdVT_GsfTrkIdT_v',  'HLT_Ele15_IsoVVVL_PFHT400_v',  'HLT_HT350_DisplacedDijet40_DisplacedTrack_v',  'HLT_HT350_DisplacedDijet80_DisplacedTrack_v',  'HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack_v',  'HLT_HT650_DisplacedDijet80_Inclusive_v',  'HLT_HT650_v',  'HLT_HT750_DisplacedDijet80_Inclusive_v',  'HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110_v',  'HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v',  'HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v',  'HLT_MET300_v',  'HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v',  'HLT_Mu15_IsoVVVL_PFHT400_PFMET50_v',  'HLT_Mu15_IsoVVVL_PFHT400_v',  'HLT_Mu15_IsoVVVL_PFHT600_v',  'HLT_Mu17_Mu8_SameSign_DZ_v',  'HLT_Mu25_TkMu0_dEta18_Onia_v',  'HLT_Mu30_TkMu11_v',  'HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose_v',  'HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose_v',  'HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight_v',  'HLT_Mu40_eta2p1_PFJet200_PFJet50_v',  'HLT_Mu6_PFHT200_PFMET100_v',  'HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v',  'HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v',  'HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53_v',  'HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52_v',  'HLT_PFHT400_SixJet30_DoubleBTagCSV_p056_v',  'HLT_PFHT450_SixJet40_BTagCSV_p056_v',  'HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v',  'HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v',  'HLT_PFHT900_v',  'HLT_PFJet450_v',  'HLT_PFJet500_v',  'HLT_PFMET400_v',  'HLT_PFMET500_v',  'HLT_PFMET600_v',  'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v',  'HLT_QuadJet45_TripleBTagCSV_p087_v',  'HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_v',  'HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500_v',  'HLT_Rsq0p25_v',  'HLT_Rsq0p30_v',  'HLT_RsqMR270_Rsq0p09_MR200_4jet_v',  'HLT_TkMu50_v',  'HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_v',  'HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5_v',  'HLT_VBF_DisplacedJet40_DisplacedTrack_v',  'HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack_v',  'HLT_VBF_DisplacedJet40_VTightID_Hadronic_v',  'HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack_v',  'HLT_VBF_DisplacedJet40_VVTightID_Hadronic_v',  'HLT_VLooseIsoPFTau140_Trk50_eta2p1_v',  'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v',  'HLT_PFMET110_PFMHT110_IDTight_v',  'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v',  'HLT_PFMET120_PFMHT120_IDTight_v',  'HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v',  'HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v',  'HLT_PFMET170_HBHECleaned_v',  'HLT_PFHT300_PFMET110_v',  'HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067_v',  'HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140_v',  'HLT_MET200_v',  'HLT_RsqMR270_Rsq0p09_MR200_v',  'HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55_v',  'HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63_v',  'HLT_MET250_v',  'HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58_v',   'HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54_v',  'HLT_PFMET300_v',  'HLT_MET75_IsoTrk50_v',  'HLT_MET90_IsoTrk50_v'
]
        ),
        metfilters = cms.InputTag('TriggerResults', '', filterString),
        metpaths = cms.vstring('Flag_HBHENoiseFilter', 'Flag_HBHENoiseIsoFilter', 'Flag_EcalDeadCellTriggerPrimitiveFilter', 'Flag_goodVertices', 'Flag_eeBadScFilter', 'Flag_globalTightHalo2016Filter','Flag_badMuons','Flag_duplicateMuons','Flag_noBadMuons') if isReMiniAod else cms.vstring('Flag_HBHENoiseFilter', 'Flag_HBHENoiseIsoFilter', 'Flag_EcalDeadCellTriggerPrimitiveFilter', 'Flag_goodVertices', 'Flag_eeBadScFilter', 'Flag_globalTightHalo2016Filter'),
        prescales = cms.InputTag('patTrigger','','PAT'),
        l1Minprescales = cms.InputTag('patTrigger','l1min','PAT'),
        l1Maxprescales = cms.InputTag('patTrigger','l1max','PAT'),
        badPFMuonFilter = cms.InputTag("BadPFMuonFilter"),
        badChCandFilter = cms.InputTag("BadChargedCandidateFilter"),
    ),
    chsJetSet = cms.PSet(
        #recoJets = cms.InputTag('ak4PFJets'),#('slimmedJetsAK8'), #selectedPatJetsAK8PFCHSPrunedPacked
        jets = cms.InputTag('slimmedJets'),#,#('slimmedJets'),#('slimmedJetsAK8'), #selectedPatJetsAK8PFCHSPrunedPacked
        jetid = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight                                                                                    
        jet1pt = cms.double(15.),
        jet2pt = cms.double(15.),
        jeteta = cms.double(2.4),
        addQGdiscriminator = cms.bool(False),
        recalibrateJets = cms.bool(True),
        recalibrateMass = cms.bool(False),
        recalibratePuppiMass = cms.bool(False),
        smearJets = cms.bool(False),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),# if not isAOD else 'offlinePrimaryVertices'),
        rho = cms.InputTag('fixedGridRhoFastjetAll'),
        jecUncertaintyDATA = cms.string('data/%s/%s_Uncertainty_AK4PFchs.txt' % (JECstring, JECstring)),#updating
        jecUncertaintyMC = cms.string('data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt'),#updating
        jecCorrectorDATA = cms.vstring(#updating
            'data/%s/%s_L1FastJet_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2Relative_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4PFchs.txt' % (JECstring, JECstring),
        ),
        jecCorrectorMC = cms.vstring(#updating!!!
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt',
        ),
        massCorrectorDATA = cms.vstring(#updating!!!
            'data/%s/%s_L2Relative_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4PFchs.txt' % (JECstring, JECstring),
        ),
        massCorrectorMC = cms.vstring(#updating!!!                                                                                                           
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt',
        ),
        massCorrectorPuppi = cms.string('data/puppiCorrSummer16.root'),#updating
        reshapeBTag = cms.bool(True),
        btag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
        btagDB = cms.string('data/CSVv2_Moriond17_B_H.csv'),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight                                                                                                       
        jet2btag = cms.int32(0),
        met = cms.InputTag('slimmedMETsMuEGClean', '', '') if isReMiniAod else cms.InputTag('slimmedMETs', '', 'ALPHA'),# if not isAOD else 'pfMet', '', 'ALPHA'),
        #recoMet = cms.InputTag('pfMet', '', 'ALPHA'),
        metRecoil = cms.bool(False),
        metRecoilMC = cms.string('data/recoilfit_gjetsMC_Zu1_pf_v5.root'),
        metRecoilData = cms.string('data/recoilfit_gjetsData_Zu1_pf_v5.root'),
        metTriggerFileName = cms.string('data/MET_trigger_eff_data_SingleMuRunBH.root'),
        jerNameRes = cms.string('data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt'),#v10 is the latest                                                                            
        jerNameSf = cms.string('data/JER/Spring16_25nsV10_MC_SF_AK4PFchs.txt'),#v10 is the latest
    ),
#PUPPI
    puppiJetSet = cms.PSet(
        #recoJets = cms.InputTag('ak4PFJets'),#('slimmedJetsAK8'), #selectedPatJetsAK8PFCHSPrunedPacked
        jets = cms.InputTag('slimmedJetsPuppi'),#,#('slimmedJets'),#('slimmedJetsAK8'), #selectedPatJetsAK8PFCHSPrunedPacked
        jetid = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight                                                                                    
        jet1pt = cms.double(15.),
        jet2pt = cms.double(15.),
        jeteta = cms.double(2.4),
        addQGdiscriminator = cms.bool(False),
        recalibrateJets = cms.bool(True),
        recalibrateMass = cms.bool(False),
        recalibratePuppiMass = cms.bool(False),
        smearJets = cms.bool(False),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),# if not isAOD else 'offlinePrimaryVertices'),
        rho = cms.InputTag('fixedGridRhoFastjetAll'),
        jecUncertaintyDATA = cms.string('data/%s/%s_Uncertainty_AK4PFPuppi.txt' % (JECstring, JECstring)),#updating
        jecUncertaintyMC = cms.string('data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFPuppi.txt'),#updating
        jecCorrectorDATA = cms.vstring(#updating
            'data/%s/%s_L1FastJet_AK4PFPuppi.txt' % (JECstring, JECstring),
            'data/%s/%s_L2Relative_AK4PFPuppi.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4PFPuppi.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4PFPuppi.txt' % (JECstring, JECstring),
        ),
        jecCorrectorMC = cms.vstring(#updating!!!
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFPuppi.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFPuppi.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFPuppi.txt',
        ),
        massCorrectorDATA = cms.vstring(#updating!!!
            'data/%s/%s_L2Relative_AK4PFPuppi.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4PFPuppi.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4PFPuppi.txt' % (JECstring, JECstring),
        ),
        massCorrectorMC = cms.vstring(#updating!!!                                                                                                           
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFPuppi.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFPuppi.txt',
        ),
        massCorrectorPuppi = cms.string('data/puppiCorrSummer16.root'),#updating
        reshapeBTag = cms.bool(True),
        btag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
        btagDB = cms.string('data/CSVv2_Moriond17_B_H.csv'),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight                                                                                                       
        jet2btag = cms.int32(0),
        met = cms.InputTag('slimmedMETsMuEGClean', '', '') if isReMiniAod else cms.InputTag('slimmedMETs', '', 'ALPHA'),# if not isAOD else 'pfMet', '', 'ALPHA'),
        #recoMet = cms.InputTag('pfMet', '', 'ALPHA'),
        metRecoil = cms.bool(False),
        metRecoilMC = cms.string('data/recoilfit_gjetsMC_Zu1_pf_v5.root'),
        metRecoilData = cms.string('data/recoilfit_gjetsData_Zu1_pf_v5.root'),
        metTriggerFileName = cms.string('data/MET_trigger_eff_data_SingleMuRunBH.root'),
        jerNameRes = cms.string('data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt'),#v10 is the latest                                                                            
        jerNameSf = cms.string('data/JER/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt'),#v10 is the latest
    ),
    jetSet = cms.PSet(
        #recoJets = cms.InputTag('ak4PFJets'),#('slimmedJetsAK8'), #selectedPatJetsAK8PFCHSPrunedPacked
        jets = cms.InputTag('patJetsAK4PF'),#('slimmedJets'),##('slimmedJetsAK8'), #selectedPatJetsAK8PFCHSPrunedPacked
        jetid = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight                                                                                    
        jet1pt = cms.double(15.),
        jet2pt = cms.double(15.),
        jeteta = cms.double(2.4),
        addQGdiscriminator = cms.bool(False),
        recalibrateJets = cms.bool(True),#(True),!!!!!!!!!!!!!!!!!!!
        recalibrateMass = cms.bool(False),
        recalibratePuppiMass = cms.bool(False),
        smearJets = cms.bool(False),#(True),!!!!!!!!!!!!!!!!!
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),# if not isAOD else 'offlinePrimaryVertices'),
        rho = cms.InputTag('fixedGridRhoFastjetAll'),
        jecUncertaintyDATA = cms.string('data/%s/%s_Uncertainty_AK4PF.txt' % (JECstring, JECstring)),#updating
        jecUncertaintyMC = cms.string('data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK4PF.txt'),#updating
        jecCorrectorDATA = cms.vstring(#updating
            'data/%s/%s_L1FastJet_AK4PF.txt' % (JECstring, JECstring),
            'data/%s/%s_L2Relative_AK4PF.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4PF.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4PF.txt' % (JECstring, JECstring),
        ),
        jecCorrectorMC = cms.vstring(#updating!!!
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PF.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PF.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PF.txt',
        ),
        massCorrectorDATA = cms.vstring(#updating!!!
            'data/%s/%s_L2Relative_AK4PF.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK4PF.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK4PF.txt' % (JECstring, JECstring),
        ),
        massCorrectorMC = cms.vstring(#updating!!!                                                                                                           
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PF.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PF.txt',
        ),
        massCorrectorPuppi = cms.string('data/puppiCorrSummer16.root'),#updating
        reshapeBTag = cms.bool(True),#(True),!!!!!
        btag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
        btagDB = cms.string('data/CSVv2_Moriond17_B_H.csv'),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight                                                                                                       
        jet2btag = cms.int32(0),
        met = cms.InputTag('slimmedMETsMuEGClean', '', '') if isReMiniAod else cms.InputTag('slimmedMETs', '', 'ALPHA'),# if not isAOD else 'pfMet', '', 'ALPHA'),
        #recoMet = cms.InputTag('pfMet', '', 'ALPHA'),
        metRecoil = cms.bool(False),
        metRecoilMC = cms.string('data/recoilfit_gjetsMC_Zu1_pf_v5.root'),
        metRecoilData = cms.string('data/recoilfit_gjetsData_Zu1_pf_v5.root'),
        metTriggerFileName = cms.string('data/MET_trigger_eff_data_SingleMuRunBH.root'),
        jerNameRes = cms.string('data/JER/Spring16_25nsV10_MC_PtResolution_AK4PF.txt'),#v10 is the latest                                                                            
        jerNameSf = cms.string('data/JER/Spring16_25nsV10_MC_SF_AK4PF.txt'),#v10 is the latest
    ),
    chsFatJetSet = cms.PSet(
        #recoJets = cms.InputTag('ak4PFJets'),#('slimmedJetsAK8'), #selectedPatJetsAK8PFCHSPrunedPacked
        jets = cms.InputTag('slimmedJetsAK8'), #selectedPatJetsAK8PFCHSPrunedPacked
        jetid = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight                                                                                    
        jet1pt = cms.double(170.),
        jet2pt = cms.double(170.),
        jeteta = cms.double(2.4),
        addQGdiscriminator = cms.bool(False),
        recalibrateJets = cms.bool(True),
        recalibrateMass = cms.bool(False),
        recalibratePuppiMass = cms.bool(False),
        smearJets = cms.bool(False),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),# if not isAOD else 'offlinePrimaryVertices'),
        rho = cms.InputTag('fixedGridRhoFastjetAll'),
        jecUncertaintyDATA = cms.string('data/%s/%s_Uncertainty_AK8PFchs.txt' % (JECstring, JECstring)),#updating
        jecUncertaintyMC = cms.string('data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK8PFchs.txt'),#updating
        jecCorrectorDATA = cms.vstring(#updating
            'data/%s/%s_L1FastJet_AK8PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2Relative_AK8PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK8PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK8PFchs.txt' % (JECstring, JECstring),
        ),
        jecCorrectorMC = cms.vstring(#updating!!!
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK8PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK8PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK8PFchs.txt',
        ),
        massCorrectorDATA = cms.vstring(#updating!!!
            'data/%s/%s_L2Relative_AK8PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK8PFchs.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK8PFchs.txt' % (JECstring, JECstring),
        ),
        massCorrectorMC = cms.vstring(#updating!!!                                                                                                           
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK8PFchs.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK8PFchs.txt',
        ),
        massCorrectorPuppi = cms.string('data/puppiCorrSummer16.root'),#updating
        reshapeBTag = cms.bool(True),
        btag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
        btagDB = cms.string('data/CSVv2_Moriond17_B_H.csv'),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight                                                                                                       
        jet2btag = cms.int32(0),
        met = cms.InputTag('slimmedMETsMuEGClean', '', '') if isReMiniAod else cms.InputTag('slimmedMETs', '', 'ALPHA'),# if not isAOD else 'pfMet', '', 'ALPHA'),
        #recoMet = cms.InputTag('pfMet', '', 'ALPHA'),
        metRecoil = cms.bool(False),
        metRecoilMC = cms.string('data/recoilfit_gjetsMC_Zu1_pf_v5.root'),
        metRecoilData = cms.string('data/recoilfit_gjetsData_Zu1_pf_v5.root'),
        metTriggerFileName = cms.string('data/MET_trigger_eff_data_SingleMuRunBH.root'),
        jerNameRes = cms.string('data/JER/Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt'),#v10 is the latest                                                                            
        jerNameSf = cms.string('data/JER/Spring16_25nsV10_MC_SF_AK8PFchs.txt'),#v10 is the latest
    ),
    fatJetSet = cms.PSet(
        jets = cms.InputTag('packedPatJetsAK8'),# '','USER'),#('patJetsAk8CHSJetsSoftDropPacked'),#('patJetsAK8PF', '','USER'),
        #jets = cms.InputTag('slimmedJetsAK8'),#('patJetsAK8PF'),#('slimmedJetsAK8'), #try AK8 first
        jetid = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight                                                                                    
        jet1pt = cms.double(50.),
        jet2pt = cms.double(50.),
        jeteta = cms.double(2.4),
        addQGdiscriminator = cms.bool(False),
        recalibrateJets = cms.bool(True),#(True),!!!!!!!!!!!!!!!!!!!
        recalibrateMass = cms.bool(False),
        recalibratePuppiMass = cms.bool(False),
        smearJets = cms.bool(False),#(True),!!!!!!!!!!!!!!!!!
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),# if not isAOD else 'offlinePrimaryVertices'),
        rho = cms.InputTag('fixedGridRhoFastjetAll'),
        jecUncertaintyDATA = cms.string('data/%s/%s_Uncertainty_AK8PF.txt' % (JECstring, JECstring)),#updating
        jecUncertaintyMC = cms.string('data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_Uncertainty_AK8PF.txt'),#updating
        jecCorrectorDATA = cms.vstring(#updating
            'data/%s/%s_L1FastJet_AK8PF.txt' % (JECstring, JECstring),
            'data/%s/%s_L2Relative_AK8PF.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK8PF.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK8PF.txt' % (JECstring, JECstring),
        ),
        jecCorrectorMC = cms.vstring(#updating!!!
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L1FastJet_AK8PF.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK8PF.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK8PF.txt',
        ),
        massCorrectorDATA = cms.vstring(#updating!!!
            'data/%s/%s_L2Relative_AK8PF.txt' % (JECstring, JECstring),
            'data/%s/%s_L3Absolute_AK8PF.txt' % (JECstring, JECstring),
            'data/%s/%s_L2L3Residual_AK8PF.txt' % (JECstring, JECstring),
        ),
        massCorrectorMC = cms.vstring(#updating!!!                                                                                                           
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L2Relative_AK8PF.txt',
            'data/Summer16_23Sep2016V3_MC/Summer16_23Sep2016V3_MC_L3Absolute_AK8PF.txt',
        ),
        massCorrectorPuppi = cms.string('data/puppiCorrSummer16.root'),#updating
        reshapeBTag = cms.bool(True),#(True),!!!!!
        btag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
        btagDB = cms.string('data/CSVv2_Moriond17_B_H.csv'),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight                                                                                                       
        jet2btag = cms.int32(0),
        met = cms.InputTag('slimmedMETsMuEGClean', '', '') if isReMiniAod else cms.InputTag('slimmedMETs', '', 'ALPHA'),# if not isAOD else 'pfMet', '', 'ALPHA'),
        #recoMet = cms.InputTag('pfMet', '', 'ALPHA'),
        metRecoil = cms.bool(False),
        metRecoilMC = cms.string('data/recoilfit_gjetsMC_Zu1_pf_v5.root'),
        metRecoilData = cms.string('data/recoilfit_gjetsData_Zu1_pf_v5.root'),
        metTriggerFileName = cms.string('data/MET_trigger_eff_data_SingleMuRunBH.root'),
        jerNameRes = cms.string('data/JER/Spring16_25nsV10_MC_PtResolution_AK8PF.txt'),#v10 is the latest                                                                            
        jerNameSf = cms.string('data/JER/Spring16_25nsV10_MC_SF_AK8PF.txt'),#v10 is the latest
    ),

    minGenBpt = cms.double(0.),
    maxGenBeta = cms.double(99999.),
    writeNJets = cms.int32(0),
    writeNFatJets = cms.int32(0),
    writeNGenBquarks = cms.int32(4),
    writeNGenLongLiveds = cms.int32(2),
    verbose = cms.bool(False),
    verboseTrigger  = cms.bool(False),
)

process.seq = cms.Sequence(
#    process.packedGenParticlesForJetsNoNu *
#    process.ak4GenJetsNoNu *
#    process.ak8GenJetsNoNu *
    process.ak4PFJets *
#    process.ak8PFJets *
#    process.ak8PFJetsPruned *
#    process.ak8PFJetsSoftDrop *
#    process.puppi *
#    process.ak8PFJetsPuppi *
#    process.ak8PFJetsPuppiSoftDrop *
#    process.patJetsAK4PF *
#    process.patJetsAK8PF *
#    process.patJetsAK8Puppi *
#    process.NjettinessAK8 *
#    process.NjettinessAK8Puppi *
#    process.ak8PFJetsPuppiValueMap *
#    process.patJetsAK8PFSoftDrop *
#    process.patJetsAK8PFSoftDropSubjets *
#    process.patJetsAK8PFPuppiSoftDrop *
#    process.patJetsAK8PFPuppiSoftDropSubjets *
#    process.patJetsAK8PFSoftDropPacked *
#    process.patJetsAK8PFPuppiSoftDropPacked *
#    process.packedPatJetsAK8 *
##unscheduled:
    process.BadPFMuonFilter *
    process.BadChargedCandidateFilter *
    process.reconstruction
)

process.p = cms.Path(process.seq)
print process.p


#open('pydump_recluster_and_reconstruct.py','w').write(process.dumpPython())
