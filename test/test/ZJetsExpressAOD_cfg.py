import FWCore.ParameterSet.Config as cms
process = cms.Process("ANA")
# ---- access the global tag (needed for the JEC) -----------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_42_V19::All'
# ---- load the reco jet configurations --------------------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
# ---- load the JEC services --------------------------------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
# ---- format the message service ---------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
# ---- maximum number of events to run over -----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(20962))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
# ---- define the source ------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#    '/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/873/80730A1C-85DD-E011-A7B2-BCAEC53296F9.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/916/B6C571DF-50DC-E011-88FC-003048D2BBF0.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/910/6E5E6A3C-EBDC-E011-BAD4-003048F11CF0.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/906/AC1E2DC1-E1DC-E011-81E8-003048F1BF66.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/906/6233A259-B2DE-E011-A8BA-003048D2C020.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/906/345F5965-F3DC-E011-808C-BCAEC5329718.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/896/8675A2DE-41DC-E011-B99B-001D09F24EE3.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/888/78F83C1C-A4DC-E011-A6CC-003048D2BE08.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/888/640F06D3-71DC-E011-8DD2-BCAEC5329708.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/888/481D452F-6CDC-E011-BF73-BCAEC532972C.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/888/3A73136A-70DC-E011-ACF1-BCAEC518FF6B.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/888/26127467-78DC-E011-9ADD-003048D2BB90.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/887/B80453D6-8CDC-E011-923E-001D09F251D1.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/887/2A8C72FC-56DD-E011-892B-001D09F23A34.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/886/F65FBB7B-B3DC-E011-AB7B-001D09F276CF.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/886/D4C0EAF2-4BDC-E011-9792-BCAEC532971C.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/886/C461BAAD-4EDC-E011-B69B-BCAEC53296FD.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/886/8626F8A1-4ADC-E011-BEDD-003048F1BF66.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/881/CC7549C5-36DC-E011-93A1-003048F1C424.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/881/82F18A2F-31DC-E011-8285-003048F11DE2.root'
#,'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/881/8243FB04-94DC-E011-8E79-BCAEC53296F4.root'

#'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/205/2208F0DF-E682-E011-BA61-000423D996C8.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/205/12A90001-CA82-E011-BC0A-001617C3B5E4.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/205/125B6282-C682-E011-B4C7-001D09F2305C.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/205/10D539F8-C982-E011-9C6E-001617C3B6C6.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/182/3C3A3D33-2382-E011-BAC5-003048F118DE.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/180/644FDD5C-1282-E011-A569-0030487CD6DA.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/168/F6345655-0182-E011-9A57-003048F024F6.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/165/3C61462B-F181-E011-B229-003048F118DE.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/163/52608B95-F581-E011-B0B8-003048F1110E.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/150/9418E87A-9B81-E011-8EB2-003048F024C2.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/129/7E0998CE-7B81-E011-8AD6-00304879FA4C.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/121/FAA66BF3-9681-E011-8FFD-003048F11CF0.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/121/E445B320-9481-E011-8293-003048F1183E.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/121/B4E24CD8-9481-E011-AA4F-003048F024FE.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/121/60EDD57D-1382-E011-AAE3-003048F024E0.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/121/5E8D272D-7A81-E011-8C0E-003048F11C28.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/121/509FEC54-9881-E011-9EF7-003048F1C424.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/121/42549B1C-9481-E011-A134-003048F1182E.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/121/267FBDCE-8D81-E011-86AC-003048F117EC.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/120/A26F9138-8481-E011-B7E4-0019B9F72BFF.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/103/86933D6B-EE80-E011-99F9-0030487CD716.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/102/4A4536AD-CF80-E011-964F-0030487CD7EA.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/099/2CA48355-D17F-E011-8703-0030487CD6E6.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/098/48FA4B21-DA7F-E011-8903-003048CFB40C.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/093/422AEE71-C27F-E011-B9EC-0030486730C6.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/088/F07DC738-E47F-E011-8C30-003048F118AC.root'
#,'/store/data/Run2011A/Photon/AOD/PromptReco-v4/000/165/071/5A2512AD-C37F-E011-AB7D-001D09F2932B.root'

# or MC
    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/06C8E521-62A4-E011-A5E9-0022198F5B1E.root',
    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/0669C32E-B69C-E011-8551-E0CB4E553643.root',
    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/06634156-5FA4-E011-BEA2-003048D293B4.root',
    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/0660AC04-B59C-E011-BEC7-E0CB4E19F9B4.root',
    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/062DC08D-AF9C-E011-B713-90E6BA0D09B9.root',
    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/061080FF-A99C-E011-90CD-0030487CDB24.root',
    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/04F65B00-A89C-E011-978E-E0CB4E5536EF.root',
    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/04B76CEA-AF9C-E011-8526-E0CB4E1A116A.root',
    '/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/047CEEC1-57A4-E011-86DC-0030487C6A1E.root'

    #'/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/06B46752-04CD-E011-B8B9-0015178C4994.root'
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/06883EDE-09CD-E011-899C-00A0D1EEE660.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/066F1232-05CD-E011-A14B-00A0D1EE8F34.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/065D3A98-0BCD-E011-8325-00A0D1EE95CC.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/0648740D-FECC-E011-A2C5-0024E876636C.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/063052DF-FCCC-E011-976C-00A0D1EE8A20.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/04920741-0DCD-E011-84FC-0024E8769B05.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/04437AFF-07CD-E011-8675-001D09676BAB.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/02A456DD-F8CC-E011-B0DF-00151796D680.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/02A3A00C-FECC-E011-82CC-0015178C0100.root',
   # '/store/mc/Summer11/GJets_TuneZ2_200_HT_inf_7TeV-madgraph/AODSIM/PU_S4_START42_V11-v1/0000/027A27F3-FFCC-E011-8210-00266CF97FF4.root'
    )
)

# ---- define the output file -------------------------------------------
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ZJetsExpress.root"),
    closeFileFast = cms.untracked.bool(True)
)
# ---- ZJetsExpress analyzer --------------------------------------------
process.accepted = cms.EDAnalyzer('ZJetsExpress',
    jets            = cms.InputTag('ak5PFJets'),
    srcRho          = cms.InputTag('kt6PFJets','rho'),
    srcRho25        = cms.InputTag('kt6PFJets25','rho'),
    minNjets        = cms.int32(1),
    jetLepIsoRadius = cms.double(0.4),
    jetLepPhoRadius = cms.double(0.4),
    minJetPt        = cms.double(30),
    maxJetEta       = cms.double(2.5),
    minPhoPt        = cms.double(70),
    maxPhoEta       = cms.double(3.0),
    minLepPt        = cms.double(20),
    maxLepEta       = cms.double(2.4),
    maxCombRelIso03 = cms.double(0.15),
    minLLMass       = cms.double(40),
    jecServiceDATA  = cms.string('ak5PFL1FastL2L3Residual'),
    jecServiceMC    = cms.string('ak5PFL1FastL2L3'),
    payload         = cms.string('AK5PF'),
    processName     = cms.string('HLT'),
    triggerName     = cms.vstring('HLT_DoubleMu6_v', ### IMPORTANT always put _v in the end of each bit 
                                  'HLT_DoubleMu7_v',
                                  'HLT_Mu13_Mu8_v',
                                  'HLT_Mu17_Mu8_v',
                                  'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v',
                                  'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v',
                                  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v',
				  #'HLT_Photon20_CaloIdVL_IsoL_v',
				  #'HLT_Photon30_CaloIdVL_v',
				  #'HLT_Photon30_CaloIdVL_IsoL_v', // we start offline at 70 GeV
				  'HLT_Photon50_CaloIdVL_IsoL_v',
				  'HLT_Photon75_CaloIdVL_v',
				  'HLT_Photon90_CaloIdVL_v',
				  'HLT_Photon90_CaloIdVL_IsoL_v',
				  'HLT_Photon125_v',
				  'HLT_Photon135_v',
         	                 ),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    triggerFamily1  = cms.vstring('HLT_DoubleMu6_v','HLT_DoubleMu7_v','HLT_Mu13_Mu8_v','HLT_Mu17_Mu8_v'),
    triggerFamily2  = cms.vstring('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v',
                                  'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v',
                                  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v'),
    triggerFamily3  = cms.vstring('HLT_Mu17_Ele8_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v'),
    triggerFamily4  = cms.vstring('HLT_Photon50_CaloIdVL_IsoL_v',
                                  'HLT_Photon75_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_IsoL_v',
                                  'HLT_Photon125_v',
                                  'HLT_Photon135_v'),

    prescaleDontAsk = cms.vstring('HLT_Mu17_Ele8_CaloIdL_v', # don't ask for L1 prescales for these bits
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v'),

)
# ---- duplicate the analyzer with different name -----------------------
process.rejected = process.accepted.clone()
# ---- filter the required HLT bits -------------------------------------
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring(
                                     'HLT_DoubleMu6_v*', # di-muon triggers
                                     'HLT_DoubleMu7_v*',
                                     'HLT_Mu13_Mu8_v*',
                                     'HLT_Mu17_Mu8*', 
                                     'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*', #di-electron trigger
                                     'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*',
                                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
				     'HLT_Mu17_Ele8_CaloIdL_v*', #di-emu trigger (for data-driven ttbar estimation)
                                     'HLT_Mu8_Ele17_CaloIdL_v*',
                                     'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*',
                         	     'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*'
),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)
############# turn-on the fastjet area calculation needed for the L1Fastjet ##############
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(5.0)
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(5.0)
############# turn-on the fastjet area in |eta|<2.5 needed for the photonISO #############
process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)
############# slimming the PFJet collection by raising the pt cut #################
process.ak5PFJets.jetPtMin = cms.double(15.0)

# ---- save all events for any trigger ---------------------
process.s0 = cms.Sequence(process.kt6PFJets + process.ak5PFJets + process.kt6PFJets25 + process.accepted)
process.p0 = cms.Path(process.s0)
process.schedule = cms.Schedule(process.p0)

# ---- first sequence: events that pass the trigger ---------------------
#process.s1 = cms.Sequence(process.hltFilter + process.kt6PFJets + process.ak5PFJets + process.accepted)
#process.p1 = cms.Path(process.s1)
# ---- second sequence: events that fail the trigger --------------------
#process.s2 = cms.Sequence(~process.hltFilter + process.kt6PFJets + process.ak5PFJets + process.rejected)
#process.p2 = cms.Path(process.s2)



# ---- schedule ---------------------------------------------------------
#process.schedule = cms.Schedule(process.p1,process.p2)