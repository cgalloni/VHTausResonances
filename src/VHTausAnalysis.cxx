// $Id: CycleCreators.py 344 2012-12-13 13:10:53Z krasznaa $

// Local include(s):
#include "../include/VHTausAnalysis.h"

// External include(s):
#include "../GoodRunsLists/include/TGoodRunsListReader.h"

#include <TMath.h>

ClassImp( VHTausAnalysis );

// define cut names

const std::string VHTausAnalysis::kCutName[ VHTausAnalysis::kNumCuts ] = {
  "BeforeCuts",            // C0
  "JSON",                  // C1
  "Trigger",               // C2
  "MetFilters",            // C3
  "Met",                   // C4
  "Jet",                   // C5
  "Tau",                   // C6
  "Lepton",                // C7
  "kTauIsolation",          // C8
  "MassWindow",               // C9
  // "HiggsWindow",           // C10
  "Tau21",                 // C11
  "SubjetSingleTag",       // C12
  "SubjetDoubleTag"        // C13
};

VHTausAnalysis::VHTausAnalysis()
   : SCycleBase()
   , m_jetAK4( this )
   , m_jetAK8( this )
   , m_eventInfo( this )
   , m_electron( this )
   , m_muon( this )
   , m_tau( this )
   , m_missingEt( this )
   , m_genParticle( this )
   , m_pileupReweightingTool( this )
   , m_bTaggingScaleTool( this )
{

   m_logger << INFO << "Hello!" << SLogger::endmsg;
   SetLogName( GetName() );
   
   // read configuration details from XML file
   // (defaults agree with July 2010 acceptance challenge)
   DeclareProperty( "RecoTreeName",             m_recoTreeName             = "physics" );
   DeclareProperty( "OutputTreeName",           m_outputTreeName           = "analysis" );
   DeclareProperty( "NtupleLevel",              m_ntupleLevel              = kTau );
   DeclareProperty( "JetAK4Name",               m_jetAK4Name               = "jetAK4" );
   DeclareProperty( "JetAK8Name",               m_jetAK8Name               = "jetAK8" );
   DeclareProperty( "ElectronName",             m_electronName             = "el" );
   DeclareProperty( "MuonName",                 m_muonName                 = "mu" );
   DeclareProperty( "TauName",                  m_tauName                  = "tau" );
   DeclareProperty( "MissingEtName",            m_missingEtName            = "MET" );
   DeclareProperty( "GenParticleName",          m_genParticleName          = "genParticle" );
   
   DeclareProperty( "IsData",                   m_isData                   = false );
   DeclareProperty( "IsSignal",                 m_isSignal                 = true );
   DeclareProperty( "ApplyMETFilters",          m_applyMETFilters          = true );
   
   DeclareProperty( "JetPtCut",                 m_jetPtCut           = 150. );
   DeclareProperty( "JetEtaCut",                m_jetEtaCut          = 2.4  );
   DeclareProperty( "MjjCut",                   m_mjjCut             = 600  );
   DeclareProperty( "JetDeltaEtaCut",           m_jetDeltaEtaCut     = 1.3  );

   DeclareProperty( "Tau21HPCut",                m_tau21HPCut             = 0.45 );
   DeclareProperty( "Tau21LPCut",                m_tau21LPCut             = 0.75 );
   DeclareProperty( "MVLowSidebandCut",          m_mVLowSidebandCut       = 40 );
   DeclareProperty( "MWLowerCut",                m_mWLowerCut             = 65 );
   DeclareProperty( "MWUpperCut",                m_mWUpperCut             = 85 );
   DeclareProperty( "MZLowerCut",                m_mZLowerCut             = 85 );
   DeclareProperty( "MZUpperCut",                m_mZUpperCut             = 105 );
   DeclareProperty( "MHLowCut",                 m_mHLowerCut             = 105 );
   DeclareProperty( "MHUpperCut",                m_mHUpperCut             = 135 );
   
   DeclareProperty( "CSVMin",                m_csvMin             = 0.605 );
   
   DeclareProperty( "ElectronPtCut",                 m_electronPtCut              = 10. );
   DeclareProperty( "ElectronEtaCut",                m_electronEtaCut             = 2.5 );
   
   DeclareProperty( "MuonPtCut",                 m_muonPtCut              = 10. );
   DeclareProperty( "MuonEtaCut",                m_muonEtaCut             = 2.4 );

   DeclareProperty( "TauPtCut",                 m_tauPtCut           = 20 );
   DeclareProperty( "TauEtaCut",                m_tauEtaCut          = 20.  );
   DeclareProperty( "LeptonPtCut",              m_leptonPtCut        = 10.   );
   DeclareProperty( "LeptonEtaCut",             m_leptonEtaCut       = 2.5  );
   DeclareProperty( "TauIsoCut",                m_tauIsoCut          = true ); 
   
   DeclareProperty( "METCut",                m_metCut             = 25 );
   
   DeclareProperty( "JSONName",                 m_jsonName             = std::string (std::getenv("SFRAME_DIR")) + "/../GoodRunsLists/JSON/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.txt" );
   
   
}

VHTausAnalysis::~VHTausAnalysis() {

  m_logger << INFO << "Tschoe!" << SLogger::endmsg;

}

void VHTausAnalysis::BeginCycle() throw( SError ) {

  m_logger << INFO << "Hello to you!" << SLogger::endmsg;
  
  //

  // Load the GRL:
  //
  if (m_isData) {
    m_logger << INFO << "Loading GoodRunsList from file: " << m_jsonName << SLogger::endmsg;
    Root::TGoodRunsListReader reader( TString( m_jsonName.c_str() ) );
    if( ! reader.Interpret() ) {
      m_logger << FATAL << "Couldn't read in the GRL!" << SLogger::endmsg;
      throw SError( ( "Couldn't read in file: " + m_jsonName ).c_str(), SError::SkipCycle );
    }
    m_grl = reader.GetMergedGoodRunsList();
    m_grl.Summary();
    m_grl.SetName( "MyGoodRunsList" );
  
    //
    // Add it as a configuration object, so that the worker nodes will receive it:
    //
    AddConfigObject( &m_grl );
		
    // m_logger << INFO << "Loading RunEventFilter from file: " << m_runEventFilterName << SLogger::endmsg;
    // m_runEventFilterReader.SetTextFile( TString( m_runEventFilterName.c_str() ) );
    // if( ! m_runEventFilterReader.Interpret() ) {
    //   m_logger << FATAL << "Couldn't read in the RunEventFilter file!" << SLogger::endmsg;
    //   throw SError( ( "Couldn't read in file: " + m_runEventFilterName ).c_str(), SError::SkipCycle );
    // }
		
  }
  
  m_triggerNames.clear();
	
  //Dijet triggers
  m_triggerNames.push_back("AK8PFJet360_TrimMass30") ;
  m_triggerNames.push_back("AK8PFHT700_TrimR0p1PT0p03Mass50") ;
  // trignames.push_back("AK8DiPFJet280_200_TrimMass30_BTagCSV0p45") ;
  m_triggerNames.push_back("PFHT650_WideJetMJJ950DEtaJJ1p5") ;
  m_triggerNames.push_back("PFHT650_WideJetMJJ900DEtaJJ1p5") ;
  m_triggerNames.push_back("PFHT800_v") ;
  
  // set names for various selections
  m_catNames.clear();

 
  // mutau channel
  m_catNames.push_back("mutau_NoWindow_NoTau21_SubjetPreTag");
 
  m_catNames.push_back("mutau_VWindow_NoTau21_SubjetPreTag");
  m_catNames.push_back("mutau_VWindow_NoTau21_SubjetNoTag");
  m_catNames.push_back("mutau_VWindow_NoTau21_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_VWindow_NoTau21_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_VWindow_NoTau21_SubjetDoubleTag");
  
  m_catNames.push_back("mutau_VWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("mutau_VWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("mutau_VWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_VWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_VWindow_Tau21LP_SubjetDoubleTag");

  m_catNames.push_back("mutau_WWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("mutau_WWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("mutau_WWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_WWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_WWindow_Tau21LP_SubjetDoubleTag");

  m_catNames.push_back("mutau_ZWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("mutau_ZWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("mutau_ZWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_ZWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_ZWindow_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("mutau_HWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("mutau_HWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("mutau_HWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_HWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_HWindow_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("mutau_VWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("mutau_VWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("mutau_VWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_VWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_VWindow_Tau21HP_SubjetDoubleTag");

  m_catNames.push_back("mutau_WWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("mutau_WWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("mutau_WWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_WWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_WWindow_Tau21HP_SubjetDoubleTag");

  m_catNames.push_back("mutau_ZWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("mutau_ZWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("mutau_ZWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_ZWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_ZWindow_Tau21HP_SubjetDoubleTag");

  m_catNames.push_back("mutau_HWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("mutau_HWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("mutau_HWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_HWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_HWindow_Tau21HP_SubjetDoubleTag");
  
  // low sideband for vector boson
  m_catNames.push_back("mutau_VLowSB_NoTau21_SubjetPreTag");
  m_catNames.push_back("mutau_VLowSB_NoTau21_SubjetNoTag");
  m_catNames.push_back("mutau_VLowSB_NoTau21_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_VLowSB_NoTau21_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_VLowSB_NoTau21_SubjetDoubleTag");
  
  m_catNames.push_back("mutau_VLowSB_Tau21LP_SubjetPreTag");
  m_catNames.push_back("mutau_VLowSB_Tau21LP_SubjetNoTag");
  m_catNames.push_back("mutau_VLowSB_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_VLowSB_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_VLowSB_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("mutau_VLowSB_Tau21HP_SubjetPreTag");
  m_catNames.push_back("mutau_VLowSB_Tau21HP_SubjetNoTag");
  m_catNames.push_back("mutau_VLowSB_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_VLowSB_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_VLowSB_Tau21HP_SubjetDoubleTag");
  
  // High sideband for Higgs boson
  m_catNames.push_back("mutau_HHighSB_NoTau21_SubjetPreTag");
  m_catNames.push_back("mutau_HHighSB_NoTau21_SubjetNoTag");
  m_catNames.push_back("mutau_HHighSB_NoTau21_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_HHighSB_NoTau21_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_HHighSB_NoTau21_SubjetDoubleTag");
  
  m_catNames.push_back("mutau_HHighSB_Tau21LP_SubjetPreTag");
  m_catNames.push_back("mutau_HHighSB_Tau21LP_SubjetNoTag");
  m_catNames.push_back("mutau_HHighSB_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_HHighSB_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_HHighSB_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("mutau_HHighSB_Tau21HP_SubjetPreTag");
  m_catNames.push_back("mutau_HHighSB_Tau21HP_SubjetNoTag");
  m_catNames.push_back("mutau_HHighSB_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("mutau_HHighSB_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("mutau_HHighSB_Tau21HP_SubjetDoubleTag");
 

  // eletau channel
  m_catNames.push_back("eletau_NoWindow_NoTau21_SubjetPreTag");
 

  m_catNames.push_back("eletau_VWindow_NoTau21_SubjetPreTag");
  m_catNames.push_back("eletau_VWindow_NoTau21_SubjetNoTag");
  m_catNames.push_back("eletau_VWindow_NoTau21_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_VWindow_NoTau21_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_VWindow_NoTau21_SubjetDoubleTag");
  
  m_catNames.push_back("eletau_VWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("eletau_VWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("eletau_VWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_VWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_VWindow_Tau21LP_SubjetDoubleTag");

  m_catNames.push_back("eletau_WWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("eletau_WWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("eletau_WWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_WWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_WWindow_Tau21LP_SubjetDoubleTag");

  m_catNames.push_back("eletau_ZWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("eletau_ZWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("eletau_ZWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_ZWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_ZWindow_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("eletau_HWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("eletau_HWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("eletau_HWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_HWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_HWindow_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("eletau_VWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("eletau_VWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("eletau_VWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_VWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_VWindow_Tau21HP_SubjetDoubleTag");

  m_catNames.push_back("eletau_WWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("eletau_WWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("eletau_WWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_WWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_WWindow_Tau21HP_SubjetDoubleTag");

  m_catNames.push_back("eletau_ZWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("eletau_ZWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("eletau_ZWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_ZWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_ZWindow_Tau21HP_SubjetDoubleTag");

  m_catNames.push_back("eletau_HWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("eletau_HWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("eletau_HWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_HWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_HWindow_Tau21HP_SubjetDoubleTag");
  
  // low sideband for vector boson
  m_catNames.push_back("eletau_VLowSB_NoTau21_SubjetPreTag");
  m_catNames.push_back("eletau_VLowSB_NoTau21_SubjetNoTag");
  m_catNames.push_back("eletau_VLowSB_NoTau21_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_VLowSB_NoTau21_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_VLowSB_NoTau21_SubjetDoubleTag");
  
  m_catNames.push_back("eletau_VLowSB_Tau21LP_SubjetPreTag");
  m_catNames.push_back("eletau_VLowSB_Tau21LP_SubjetNoTag");
  m_catNames.push_back("eletau_VLowSB_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_VLowSB_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_VLowSB_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("eletau_VLowSB_Tau21HP_SubjetPreTag");
  m_catNames.push_back("eletau_VLowSB_Tau21HP_SubjetNoTag");
  m_catNames.push_back("eletau_VLowSB_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_VLowSB_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_VLowSB_Tau21HP_SubjetDoubleTag");
  
  // High sideband for Higgs boson
  m_catNames.push_back("eletau_HHighSB_NoTau21_SubjetPreTag");
  m_catNames.push_back("eletau_HHighSB_NoTau21_SubjetNoTag");
  m_catNames.push_back("eletau_HHighSB_NoTau21_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_HHighSB_NoTau21_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_HHighSB_NoTau21_SubjetDoubleTag");
  
  m_catNames.push_back("eletau_HHighSB_Tau21LP_SubjetPreTag");
  m_catNames.push_back("eletau_HHighSB_Tau21LP_SubjetNoTag");
  m_catNames.push_back("eletau_HHighSB_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_HHighSB_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_HHighSB_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("eletau_HHighSB_Tau21HP_SubjetPreTag");
  m_catNames.push_back("eletau_HHighSB_Tau21HP_SubjetNoTag");
  m_catNames.push_back("eletau_HHighSB_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("eletau_HHighSB_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("eletau_HHighSB_Tau21HP_SubjetDoubleTag");
  
 // tautau channel

  m_catNames.push_back("tautau_NoWindow_NoTau21_SubjetPreTag");
 
  m_catNames.push_back("tautau_VWindow_NoTau21_SubjetPreTag");
  m_catNames.push_back("tautau_VWindow_NoTau21_SubjetNoTag");
  m_catNames.push_back("tautau_VWindow_NoTau21_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_VWindow_NoTau21_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_VWindow_NoTau21_SubjetDoubleTag");
  
  m_catNames.push_back("tautau_VWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("tautau_VWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("tautau_VWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_VWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_VWindow_Tau21LP_SubjetDoubleTag");

  m_catNames.push_back("tautau_WWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("tautau_WWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("tautau_WWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_WWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_WWindow_Tau21LP_SubjetDoubleTag");

  m_catNames.push_back("tautau_ZWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("tautau_ZWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("tautau_ZWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_ZWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_ZWindow_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("tautau_HWindow_Tau21LP_SubjetPreTag");
  m_catNames.push_back("tautau_HWindow_Tau21LP_SubjetNoTag");
  m_catNames.push_back("tautau_HWindow_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_HWindow_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_HWindow_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("tautau_VWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("tautau_VWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("tautau_VWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_VWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_VWindow_Tau21HP_SubjetDoubleTag");

  m_catNames.push_back("tautau_WWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("tautau_WWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("tautau_WWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_WWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_WWindow_Tau21HP_SubjetDoubleTag");

  m_catNames.push_back("tautau_ZWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("tautau_ZWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("tautau_ZWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_ZWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_ZWindow_Tau21HP_SubjetDoubleTag");

  m_catNames.push_back("tautau_HWindow_Tau21HP_SubjetPreTag");
  m_catNames.push_back("tautau_HWindow_Tau21HP_SubjetNoTag");
  m_catNames.push_back("tautau_HWindow_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_HWindow_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_HWindow_Tau21HP_SubjetDoubleTag");
  
  // low sideband for vector boson
  m_catNames.push_back("tautau_VLowSB_NoTau21_SubjetPreTag");
  m_catNames.push_back("tautau_VLowSB_NoTau21_SubjetNoTag");
  m_catNames.push_back("tautau_VLowSB_NoTau21_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_VLowSB_NoTau21_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_VLowSB_NoTau21_SubjetDoubleTag");
  
  m_catNames.push_back("tautau_VLowSB_Tau21LP_SubjetPreTag");
  m_catNames.push_back("tautau_VLowSB_Tau21LP_SubjetNoTag");
  m_catNames.push_back("tautau_VLowSB_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_VLowSB_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_VLowSB_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("tautau_VLowSB_Tau21HP_SubjetPreTag");
  m_catNames.push_back("tautau_VLowSB_Tau21HP_SubjetNoTag");
  m_catNames.push_back("tautau_VLowSB_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_VLowSB_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_VLowSB_Tau21HP_SubjetDoubleTag");
  
  // High sideband for Higgs boson
  m_catNames.push_back("tautau_HHighSB_NoTau21_SubjetPreTag");
  m_catNames.push_back("tautau_HHighSB_NoTau21_SubjetNoTag");
  m_catNames.push_back("tautau_HHighSB_NoTau21_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_HHighSB_NoTau21_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_HHighSB_NoTau21_SubjetDoubleTag");
  
  m_catNames.push_back("tautau_HHighSB_Tau21LP_SubjetPreTag");
  m_catNames.push_back("tautau_HHighSB_Tau21LP_SubjetNoTag");
  m_catNames.push_back("tautau_HHighSB_Tau21LP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_HHighSB_Tau21LP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_HHighSB_Tau21LP_SubjetDoubleTag");
  
  m_catNames.push_back("tautau_HHighSB_Tau21HP_SubjetPreTag");
  m_catNames.push_back("tautau_HHighSB_Tau21HP_SubjetNoTag");
  m_catNames.push_back("tautau_HHighSB_Tau21HP_SubjetSingleTagIncl");
  m_catNames.push_back("tautau_HHighSB_Tau21HP_SubjetSingleTagExcl");
  m_catNames.push_back("tautau_HHighSB_Tau21HP_SubjetDoubleTag");
  

   return;

}

void VHTausAnalysis::EndCycle() throw( SError ) {

   return;

}

void VHTausAnalysis::BeginInputData( const SInputData& id ) throw( SError ) {

  m_logger << INFO << "RecoTreeName:         " <<             m_recoTreeName << SLogger::endmsg;
  m_logger << INFO << "OutputTreeName:       " <<             m_outputTreeName << SLogger::endmsg;
  m_logger << INFO << "NtupleLevel:          " <<             m_ntupleLevel << SLogger::endmsg;
  m_logger << INFO << "JetAK4Name:           " <<             m_jetAK4Name << SLogger::endmsg;
  m_logger << INFO << "JetAK8Name:           " <<             m_jetAK8Name << SLogger::endmsg;
  m_logger << INFO << "ElectronName:         " <<             m_electronName << SLogger::endmsg;
  m_logger << INFO << "MuonName:             " <<             m_muonName << SLogger::endmsg;
  m_logger << INFO << "TauName:             " <<             m_tauName << SLogger::endmsg;
  m_logger << INFO << "GenParticleName:      " <<             m_genParticleName << SLogger::endmsg;
  
  m_logger << INFO << "IsData:           " <<                   (m_isData ? "TRUE" : "FALSE") << SLogger::endmsg;
  m_logger << INFO << "IsSignal:           " <<                 (m_isSignal ? "TRUE" : "FALSE") << SLogger::endmsg;
  m_logger << INFO << "ApplyMETFilters:           " <<          (m_applyMETFilters ? "TRUE" : "FALSE") << SLogger::endmsg;
  
  m_logger << INFO << "JetPtCut:           " <<                 m_jetPtCut << SLogger::endmsg;
  m_logger << INFO << "JetEtaCut:           " <<                m_jetEtaCut << SLogger::endmsg;
  m_logger << INFO << "MjjCut:           " <<                m_mjjCut << SLogger::endmsg;
  m_logger << INFO << "JetDeltaEtaCut:           " <<                m_jetDeltaEtaCut << SLogger::endmsg;

  m_logger << INFO << "Tau21HPCut:           " <<                m_tau21HPCut << SLogger::endmsg;
  m_logger << INFO << "Tau21LPCut:           " <<                m_tau21LPCut << SLogger::endmsg;
  
  m_logger << INFO << "MVLowSidebandCut:     " <<                m_mVLowSidebandCut << SLogger::endmsg;
  m_logger << INFO << "MWLowerCut:           " <<                m_mWLowerCut << SLogger::endmsg;
  m_logger << INFO << "MWUpperCut:           " <<                m_mWUpperCut << SLogger::endmsg;
  m_logger << INFO << "MZLowerCut:           " <<                m_mZLowerCut << SLogger::endmsg;
  m_logger << INFO << "MZUpperCut:           " <<                m_mZUpperCut << SLogger::endmsg;
  m_logger << INFO << "MHLowCut:             " <<                 m_mHLowerCut << SLogger::endmsg;
  m_logger << INFO << "MHUpperCut:           " <<                m_mHUpperCut << SLogger::endmsg;
  
  m_logger << INFO << "CSVMin:           " <<                m_csvMin << SLogger::endmsg;
  
  m_logger << INFO << "ElectronPtCut:           " <<                 m_electronPtCut << SLogger::endmsg;
  m_logger << INFO << "ElectronEtaCut:           " <<                m_electronEtaCut << SLogger::endmsg;
  
  m_logger << INFO << "MuonPtCut:           " <<                 m_muonPtCut << SLogger::endmsg;
  m_logger << INFO << "MuonEtaCut:           " <<                m_muonEtaCut << SLogger::endmsg;
  
  m_logger << INFO << "TauPtCut:           " <<                 m_tauPtCut << SLogger::endmsg;
  m_logger << INFO << "TauEtaCut:           " <<                m_tauEtaCut << SLogger::endmsg;
  m_logger << INFO << "TauIsoCut:           " <<                m_tauIsoCut << SLogger::endmsg;
  
  m_logger << INFO << "LeptonPtCut:           " <<                 m_leptonPtCut << SLogger::endmsg;
  m_logger << INFO << "LeptonEtaCut:           " <<                m_leptonEtaCut << SLogger::endmsg;
  
  m_logger << INFO << "METCut:           " <<                m_metCut << SLogger::endmsg;
  
  m_logger << INFO << "JSONName:           " <<                 m_jsonName << SLogger::endmsg;
  
  if (!m_isData) m_pileupReweightingTool.BeginInputData( id );
  
  if (m_isData) {
    TObject* grl;
    if( ! ( grl = GetConfigObject( "MyGoodRunsList" ) ) ) {
      m_logger << FATAL << "Can't access the GRL!" << SLogger::endmsg;
      throw SError( "Can't access the GRL!", SError::SkipCycle );
    }
    m_grl = *( dynamic_cast< Root::TGoodRunsList* >( grl ) );
  }
  
  //
  // output branches
  //
  DeclareVariable(b_weight              , "weight"                 );
  DeclareVariable(b_weightGen           , "weightGen"              );
  DeclareVariable(b_weightPU            , "weightPU"               );
  DeclareVariable(b_weightBtag          , "weightBtag"             );
  
  DeclareVariable( b_runNumber,           "b_runNumber"            );
  DeclareVariable( b_eventNumber,         "b_eventNumber"          );
  DeclareVariable( b_lumiBlock,           "b_lumiBlock"            );
  
  DeclareVariable( b_ak8jet0_pt,           "b_ak8jet0_pt"          );
  DeclareVariable( b_ak8jet0_phi,          "b_ak8jet0_phi"         );
  DeclareVariable( b_ak8jet0_eta,          "b_ak8jet0_eta"         );
  DeclareVariable( b_ak8jet0_e,            "b_ak8jet0_e"           );
  DeclareVariable( b_ak8jet0_tau21,        "b_ak8jet0_tau21"       );
  DeclareVariable( b_ak8jet0_m,            "b_ak8jet0_m"           );
  DeclareVariable( b_ak8jet0_mpruned,      "b_ak8jet0_mpruned"     );
  DeclareVariable( b_ak8jet0_csv,          "b_ak8jet0_csv"         );
  DeclareVariable( b_ak8jet1_pt,           "b_ak8jet1_pt"          );
  DeclareVariable( b_ak8jet1_phi,          "b_ak8jet1_phi"         );
  DeclareVariable( b_ak8jet1_eta,          "b_ak8jet1_eta"         );
  DeclareVariable( b_ak8jet1_e,            "b_ak8jet1_e"           );
  DeclareVariable( b_ak8jet1_tau21,        "b_ak8jet1_tau21"       );
  DeclareVariable( b_ak8jet1_m,            "b_ak8jet1_m"           );
  DeclareVariable( b_ak8jet1_mpruned,      "b_ak8jet1_mpruned"     );
  DeclareVariable( b_ak8jet1_csv,          "b_ak8jet1_csv"         );
  DeclareVariable( b_ak8jet0_subjet01_dr,   "b_ak8jet0_subjet01_dr"  );
  DeclareVariable( b_ak8jet0_subjet01_deta, "b_ak8jet0_subjet01_deta");
  DeclareVariable( b_ak8jet0_subjet01_dphi, "b_ak8jet0_subjet01_dphi");
  DeclareVariable( b_ak8jet0_subjet0_pt,   "b_ak8jet0_subjet0_pt"  );
  DeclareVariable( b_ak8jet0_subjet1_pt,   "b_ak8jet0_subjet1_pt"  );
  DeclareVariable( b_ak8jet0_subjet0_csv,  "b_ak8jet0_subjet0_csv" );
  DeclareVariable( b_ak8jet0_subjet1_csv,  "b_ak8jet0_subjet1_csv" );
  DeclareVariable( b_ak8jet1_subjet01_dr,   "b_ak8jet1_subjet01_dr"  );
  DeclareVariable( b_ak8jet1_subjet01_deta, "b_ak8jet1_subjet01_deta");
  DeclareVariable( b_ak8jet1_subjet01_dphi, "b_ak8jet1_subjet01_dphi");
  DeclareVariable( b_ak8jet1_subjet0_pt,   "b_ak8jet1_subjet0_pt"  );
  DeclareVariable( b_ak8jet1_subjet1_pt,   "b_ak8jet1_subjet1_pt"  );
  DeclareVariable( b_ak8jet1_subjet0_csv,  "b_ak8jet1_subjet0_csv" );
  DeclareVariable( b_ak8jet1_subjet1_csv,  "b_ak8jet1_subjet1_csv" );
  
  b_selection_bits.resize( m_catNames.size() );
  b_selection_lastcut.resize( m_catNames.size() );
  for (unsigned int s=0;s<m_catNames.size();++s) {
    DeclareVariable(b_selection_bits[s]    , Form("selection_bits_%s", m_catNames[s].c_str())    );
    DeclareVariable(b_selection_lastcut[s] , Form("selection_lastcut_%s", m_catNames[s].c_str()) );
  }
  
  //
  // Declare the output histograms:
  //
  Book( TH1F( "Events", "Events;weight", 10, -.5, 9.5 ) );
  Book( TH1F( "SumEvents", "SumEvents;weight", 10, -.5, 9.5 ) );
  Book( TH1F( "METFilters", "METFilters", 20, 0.5, 20.5 ));

  

  for (unsigned int s=0;s<m_catNames.size();++s) {
    TString directory = m_catNames[s].c_str();
    // cutflow
    Book( TH1F( "cutflow", "cutflow", 20, 0.5, 20.5 ), directory );  
    bookHistograms(directory);
    bookTriggerHistos(directory);
  }
  
  // b-tagging tool initialisation
  m_bTaggingScaleTool.BeginInputData( id );
  if (m_isSignal) {
    // b-tagging efficiencies
    m_bTaggingScaleTool.bookHistograms();
  }
  
 

  return;

}

void VHTausAnalysis::EndInputData( const SInputData& ) throw( SError ) {

  //
  // Final analysis of cut flow
  //
  
  TString defaultCutflow = m_catNames[m_catNames.size()-1]; //"VWindow_Tau21HP_SubjetDoubleTag";
  m_logger << INFO << "cut flow for " << defaultCutflow << SLogger::endmsg;
  m_logger << INFO << Form( "Cut\t%25.25s\tEvents\tRelEff\tAbsEff", "Name" ) << SLogger::endmsg;
  
  Double_t ntot = Hist( "cutflow", defaultCutflow )->GetBinContent( 1 );
  m_logger << INFO << Form( "\t%25.25s\t%6.0f", "start", ntot ) << SLogger::endmsg;
  for( Int_t ibin = 2; ibin <= kNumCuts; ++ibin ) {
    Int_t    icut    = ibin - 1;
    Double_t nevents = Hist( "cutflow", defaultCutflow )->GetBinContent( ibin );
    Double_t releff  = 100. * nevents / Hist( "cutflow", defaultCutflow )->GetBinContent( ibin-1 );
    Double_t abseff  = 100. * nevents / ntot;
    m_logger << INFO  << Form( "C%d\t%25.25s\t%6.0f\t%6.2f\t%6.2f", icut-1, kCutName[icut].c_str(), nevents, releff, abseff ) << SLogger::endmsg;
  }
   
   return;

}

void VHTausAnalysis::BeginInputFile( const SInputData& ) throw( SError ) {

  m_logger << INFO << "Connecting input variables" << SLogger::endmsg;
  if (m_isData) {
    // m_jetAK4.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis, (m_jetAK4Name + "_").c_str() );
    m_jetAK8.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetSubstructure|Ntuple::JetPrunedSubjets, (m_jetAK8Name + "_").c_str() );
    m_eventInfo.ConnectVariables(    m_recoTreeName.c_str(), Ntuple::EventInfoBasic|Ntuple::EventInfoTrigger|Ntuple::EventInfoMETFilters, "" );
  }
  else {
    // m_jetAK4.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetTruth, (m_jetAK4Name + "_").c_str() );
    m_jetAK8.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetSubstructure|Ntuple::JetTruth|Ntuple::JetPrunedSubjets|Ntuple::JetPrunedSubjetsTruth, (m_jetAK8Name + "_").c_str() );
    m_eventInfo.ConnectVariables(    m_recoTreeName.c_str(), Ntuple::EventInfoBasic|Ntuple::EventInfoTrigger|Ntuple::EventInfoMETFilters|Ntuple::EventInfoTruth, "" );
    m_genParticle.ConnectVariables(  m_recoTreeName.c_str(), Ntuple::GenParticleBasic, (m_genParticleName + "_").c_str() );
  }
  m_electron.ConnectVariables(     m_recoTreeName.c_str(), Ntuple::ElectronBasic|Ntuple::ElectronID|Ntuple::ElectronBoostedID, (m_electronName + "_").c_str() );
  m_muon.ConnectVariables(         m_recoTreeName.c_str(), Ntuple::MuonBasic|Ntuple::MuonID|Ntuple::MuonIsolation, (m_muonName + "_").c_str() );
  m_tau.ConnectVariables(         m_recoTreeName.c_str(), Ntuple::TauBasic|Ntuple::TauID|Ntuple::TauAdvancedID, (m_tauName + "_").c_str() );

  m_missingEt.ConnectVariables(    m_recoTreeName.c_str(), Ntuple::MissingEtBasic, (m_missingEtName + "_").c_str() );
  
  m_logger << INFO << "Connecting input variables completed" << SLogger::endmsg;

   return;

}

void VHTausAnalysis::ExecuteEvent( const SInputData&, Double_t ) throw( SError ) {

  m_logger << VERBOSE << "ExecuteEvent" << SLogger::endmsg;
  
  clearBranches();
  
  b_eventNumber = m_eventInfo.eventNumber;
  b_runNumber = m_eventInfo.runNumber;
  b_lumiBlock = m_eventInfo.lumiBlock;
  
  std::vector<TBits> selectionBits(m_catNames.size(), TBits(kNumCuts));
  for (unsigned int s=0;s<m_catNames.size();++s) {
    selectionBits[s].SetBitNumber( kBeforeCuts );
  }
  bool moveOn = true;
  
  // Cut 1: check for data if run/lumiblock in JSON
  if (m_isData) {
    if (isGoodEvent(m_eventInfo.runNumber, m_eventInfo.lumiBlock)) {
      for (unsigned int s=0;s<m_catNames.size();++s) {
        selectionBits[s].SetBitNumber( kJSON );
      }
    }
  }
  else {
    for (unsigned int s=0;s<m_catNames.size();++s) {
      selectionBits[s].SetBitNumber( kJSON );
    }
  }
  
  // Cut 2: pass trigger
 
  if (passTrigger()) {
    m_logger << VERBOSE << "Trigger pass" << SLogger::endmsg;
    for (unsigned int s=0;s<m_catNames.size();++s) {
      selectionBits[s].SetBitNumber( kTrigger );
    }
    // m_logger << INFO << "pass: " << selectionBits[0].TestBitNumber( kTrigger ) << SLogger::endmsg;
  }
  else {

    if (m_isSignal){
      for (unsigned int s=0;s<m_catNames.size();++s) {
	selectionBits[s].SetBitNumber( kTrigger );
      }
    }
    m_logger <<  VERBOSE  << "Trigger fail:" << selectionBits[0].TestBitNumber( kTrigger ) << SLogger::endmsg;
    for (unsigned int s=0;s<m_catNames.size();++s) {
      m_logger << VERBOSE  << selectionBits[s].TestBitNumber( kTrigger ) << SLogger::endmsg;
    }
  }
  
  // Cut 3: pass MET filters


  if (passMETFilters()) {
    m_logger << VERBOSE << "passMETFilters" << SLogger::endmsg;
    for (unsigned int s=0;s<m_catNames.size();++s) {
      selectionBits[s].SetBitNumber( kMetFilters );
    }
  }
  
  // Cut 4: select two fat jets
  std::vector<UZH::Jet> goodFatJets;
  for ( int i = 0; i < (m_jetAK8.N); ++i ) {
    UZH::Jet myjet( &m_jetAK8, i );
    if (myjet.pt() > m_jetPtCut) {
      if (fabs(myjet.eta()) < m_jetEtaCut) {

	// std::cout<<" myjet.pt() " << myjet.pt() <<" myjet.eta() " << myjet.eta() << " myjet.IDTight() " << myjet.IDTight() << 
	//   " myjet.pruned_massCorr() " << myjet.pruned_massCorr() <<std::endl;

        if (myjet.IDTight()) {
          goodFatJets.push_back(myjet);
          
          // fill output variables
          if (goodFatJets.size() == 1) {
            b_ak8jet0_pt = myjet.pt();
            b_ak8jet0_phi = myjet.phi();
            b_ak8jet0_eta = myjet.eta();
            b_ak8jet0_e = myjet.e();
            double tau21 = -1;
            if (myjet.tau1() != 0) {
              tau21 = myjet.tau2() / myjet.tau1();
            }
            b_ak8jet0_tau21 = tau21;
            b_ak8jet0_m = myjet.m();
            b_ak8jet0_mpruned = myjet.pruned_massCorr();
            b_ak8jet0_csv = myjet.csv();
            if (myjet.subjet_pruned_N() >= 2) {
              double deta = fabs(myjet.subjet_pruned_eta()[0] - myjet.subjet_pruned_eta()[1]);
              double dphi = fabs(myjet.subjet_pruned_phi()[0] - myjet.subjet_pruned_phi()[1]);
              double dr = sqrt(deta*deta + dphi*dphi);
              b_ak8jet0_subjet01_dr = dr;
              b_ak8jet0_subjet01_deta = deta;
              b_ak8jet0_subjet01_dphi = dphi;
              b_ak8jet0_subjet0_pt = myjet.subjet_pruned_pt()[0];
              b_ak8jet0_subjet1_pt = myjet.subjet_pruned_pt()[1];
              b_ak8jet0_subjet0_csv = myjet.subjet_pruned_csv()[0];
              b_ak8jet0_subjet1_csv = myjet.subjet_pruned_csv()[1];
            }
          }
          else if (goodFatJets.size() == 2) {
            b_ak8jet1_pt = myjet.pt();
            b_ak8jet1_phi = myjet.phi();
            b_ak8jet1_eta = myjet.eta();
            b_ak8jet1_e = myjet.e();
            double tau21 = -1;
            if (myjet.tau1() != 0) {
              tau21 = myjet.tau2() / myjet.tau1();
            }
            b_ak8jet1_tau21 = tau21;
            b_ak8jet1_m = myjet.m();
            b_ak8jet1_mpruned = myjet.pruned_massCorr();
            b_ak8jet1_csv = myjet.csv();
            if (myjet.subjet_pruned_N() >= 2) {
              double deta = fabs(myjet.subjet_pruned_eta()[0] - myjet.subjet_pruned_eta()[1]);
              double dphi = fabs(myjet.subjet_pruned_phi()[0] - myjet.subjet_pruned_phi()[1]);
              double dr = sqrt(deta*deta + dphi*dphi);
              b_ak8jet1_subjet01_dr = dr;
              b_ak8jet1_subjet01_deta = deta;
              b_ak8jet1_subjet01_dphi = dphi;
              b_ak8jet1_subjet0_pt = myjet.subjet_pruned_pt()[0];
              b_ak8jet1_subjet1_pt = myjet.subjet_pruned_pt()[1];
              b_ak8jet1_subjet0_csv = myjet.subjet_pruned_csv()[0];
              b_ak8jet1_subjet1_csv = myjet.subjet_pruned_csv()[1];
            }
          }
        }
      }
    }
  }


  if (goodFatJets.size() >= 1) {
    for (unsigned int s=0;s<m_catNames.size();++s) {
      selectionBits[s].SetBitNumber( kJet );
      if (m_catNames[s].find("NoWindow") != std::string::npos) {
	selectionBits[s].SetBitNumber( kMassWindow );
      }
    }
  }
  else {
    // can only continue with at least 1 selected fat jet
    moveOn = false;
  }
  
  m_logger << VERBOSE << "kJet" << SLogger::endmsg;
  
  int goodFatJet1Index = -1;

  
 
  UZH::Jet Jet;
  m_logger << VERBOSE << "kFatJetsDeltaEta" << SLogger::endmsg;
  
  if (moveOn) { // need to have selected the One candidate jet
    
    moveOn = false;
  
    m_logger << VERBOSE << "kDijetMass" << SLogger::endmsg;

    // Cut 7: require one of the jets to be in the V-boson mass window
    // make other selected jet the Higgs jet
    
    for (unsigned int i = 0; i < goodFatJets.size(); ++i) {
      moveOn=false;
      
      if ((goodFatJets[i].pruned_massCorr() > m_mWLowerCut) && (goodFatJets[i].pruned_massCorr() <= m_mZUpperCut)) {
	Jet = goodFatJets[i];
	goodFatJet1Index=i;
	moveOn = true;
      }
      
      if (moveOn) {
	for (unsigned int s=0;s<m_catNames.size();++s) {
	  if (m_catNames[s].find("VWindow") != std::string::npos) {
	    selectionBits[s].SetBitNumber( kMassWindow );
	  }
	  if ((m_catNames[s].find("WWindow") != std::string::npos) && (goodFatJets[i].pruned_massCorr() <= m_mWUpperCut)) {
	    selectionBits[s].SetBitNumber( kMassWindow );
	  }
	  if ((m_catNames[s].find("ZWindow") != std::string::npos) && (goodFatJets[i].pruned_massCorr() > m_mZLowerCut)) {
	    selectionBits[s].SetBitNumber( kMassWindow );
	  }
	}
      }
      else { // this means we're not in the V region of the pruned mass

	if ((goodFatJets[i].pruned_massCorr() > m_mVLowSidebandCut) && (goodFatJets[i].pruned_massCorr() <= m_mWLowerCut)) {
	  goodFatJet1Index=i;
	  Jet = goodFatJets[goodFatJet1Index];
        
	  moveOn = true;
	}
      
	if (moveOn) {
	  for (unsigned int s=0;s<m_catNames.size();++s) {
	    if (m_catNames[s].find("VLowSB") != std::string::npos) {
	      selectionBits[s].SetBitNumber( kMassWindow );
	    }
	  }      
	}
	
      }
      //  }// moveOn after having selected two candidate jets
     
      if (moveOn) { // move on only if V or Vlow sideband

	if (m_isSignal) {
	  std::vector<UZH::Jet> selectedJets;
	  selectedJets.push_back(Jet);
	  m_bTaggingScaleTool.fillEfficiencies(selectedJets);
	  m_bTaggingScaleTool.fillPrunedSubjetEfficiencies(selectedJets);
	}
      }
     
      m_logger << VERBOSE << "kMassWindow" << SLogger::endmsg;
    


      // Cut 8: check if Higgs candidate jet is in Higgs mass window
      if ((goodFatJets[i].pruned_massCorr() > m_mHLowerCut) && (goodFatJets[i].pruned_massCorr() <= m_mHUpperCut)) {
	//std::cout << "higgsJet mass SR "<<goodFatJets[i].pruned_massCorr()<<  std::endl;
	for (unsigned int s=0;s<m_catNames.size();++s) {
	  if (m_catNames[s].find("HWindow") != std::string::npos) { // HWindow category
	    selectionBits[s].SetBitNumber( kMassWindow );
	  }
	}
	goodFatJet1Index=i;
	Jet = goodFatJets[goodFatJet1Index];
	moveOn = true;
      }
      else if (goodFatJets[i].pruned_massCorr() > m_mHUpperCut ){

	for (unsigned int s=0;s<m_catNames.size();++s) {
	  if (m_catNames[s].find("HHighSB") != std::string::npos) { // HHighSB category
	    selectionBits[s].SetBitNumber( kMassWindow );
	  }
	}
	goodFatJet1Index=i;
	Jet = goodFatJets[goodFatJet1Index];
	moveOn = true;
      }
      if (moveOn) break; 
    }
      
    
    
    m_logger << VERBOSE << "kMassWindow" << SLogger::endmsg;

    // m_logger << VERBOSE << Jet << SLogger::endmsg;    
    
    if (moveOn==false) return;
    
    // Cut 9: pass tau21 for V-candidate jet
    //found at least 1 jet Vlow, v, H or HHigh
      
    if (Jet.tau1() != 0) {
      double tau21 = Jet.tau2() / Jet.tau1();
      m_logger << VERBOSE << tau21 << SLogger::endmsg;
      if (tau21 < m_tau21HPCut) {
        for (unsigned int s=0;s<m_catNames.size();++s) {
          if ((m_catNames[s].find("Tau21HP") != std::string::npos)) {
            m_logger << VERBOSE << m_catNames[s] << SLogger::endmsg;
            selectionBits[s].SetBitNumber( kTau21 );
          }
        }
      }
      else if ((tau21 >= m_tau21HPCut) && (tau21 < m_tau21LPCut)) {
        for (unsigned int s=0;s<m_catNames.size();++s) {
          if ((m_catNames[s].find("Tau21LP") != std::string::npos)) {
            m_logger << VERBOSE << m_catNames[s] << SLogger::endmsg;
            selectionBits[s].SetBitNumber( kTau21 );
          }
        }
      }
    }
    for (unsigned int s=0;s<m_catNames.size();++s) {
      if (m_catNames[s].find("NoTau21") != std::string::npos) {
        selectionBits[s].SetBitNumber( kTau21 );
      }
    }
    
    m_logger << VERBOSE << "kTau21" << SLogger::endmsg;
    
    // Cut 10: require at least one of the subjets from the Higgs jet to be b-tagged
    // Cut 11: require two subjets from the Higgs jet to be b-tagged
    // count number of b-tagged subjets
    int nTaggedSubjets = 0;
    for (int i = 0; i < Jet.subjet_pruned_N(); ++i) {
      if (m_bTaggingScaleTool.isTagged(Jet.subjet_pruned_csv()[i])) {
        ++nTaggedSubjets;
      }
    }
    for (unsigned int s=0;s<m_catNames.size();++s) {
      if (m_catNames[s].find("NoTag") != std::string::npos) {
        if (nTaggedSubjets == 0) {
          selectionBits[s].SetBitNumber( kSubjetSingleTag );
          selectionBits[s].SetBitNumber( kSubjetDoubleTag );
        }
      }
      else if (m_catNames[s].find("PreTag") != std::string::npos) {
        selectionBits[s].SetBitNumber( kSubjetSingleTag );
        selectionBits[s].SetBitNumber( kSubjetDoubleTag );
      }
      else if (m_catNames[s].find("SingleTagExcl") != std::string::npos) {
        if (nTaggedSubjets == 1) {
          selectionBits[s].SetBitNumber( kSubjetSingleTag );
          selectionBits[s].SetBitNumber( kSubjetDoubleTag );
        }
      }
      else if (m_catNames[s].find("SingleTagIncl") != std::string::npos) {
        if (nTaggedSubjets >= 1) {
          selectionBits[s].SetBitNumber( kSubjetSingleTag );
          selectionBits[s].SetBitNumber( kSubjetDoubleTag );
        }
      }
      else if (nTaggedSubjets >= 1) {
        // these will all be DoubleTag
        selectionBits[s].SetBitNumber( kSubjetSingleTag );
        if (nTaggedSubjets >= 2) {
          selectionBits[s].SetBitNumber( kSubjetDoubleTag );
        }
      }
    }
    TLorentzVector diJet = Jet.tlv() // + higgsJet.tlv()
      ;
   
    // selection done
    m_logger << VERBOSE << "selection done" << SLogger::endmsg;
    if (!m_isData) {
      b_weight = getEventWeight();
      std::vector<UZH::Jet> selectedJets;
      selectedJets.push_back(Jet);
    
      // b_weightBtag = m_bTaggingScaleTool.getScaleFactor(selectedJets); // event b-tag SF weight
      b_weightBtag = m_bTaggingScaleTool.getPrunedSubjetScaleFactor(selectedJets); // event b-tag SF weight
      
      b_weight *= b_weightBtag;
    }
  
    // std::vector<UZH::Jet> goodJetsAK4;
    // for ( int i = 0; i < (m_jetAK4.N); ++i ) {
    //   UZH::Jet myjet( &m_jetAK4, i );
    //   if (fabs(myjet.eta()) < m_jetEtaCut) {
    //     if (myjet.pt() > m_jetPtCut) {
    //       goodJetsAK4.push_back(myjet);
    //       // m_logger << INFO << myjet.pt() << SLogger::endmsg;
    //       Book( TH1F( "Jet_pt", "Jet p_{T} [GeV]", 100, 0.0, 1000 ), "AK4_low" )->Fill( myjet.pt() );
    //       Book( TH1F( "Jet_eta", "Jet #eta", 50, -2.5, 2.5 ), "AK4_low" )->Fill( myjet.eta() );
    //       Book( TH1F( "Jet_m", "Jet m [GeV]", 50, 0, 250 ), "AK4_low" )->Fill( myjet.m() );
    //       Book( TH1F( "Jet_phi", "Jet #phi", 50, -TMath::Pi(), TMath::Pi() ), "AK4_low" )->Fill( myjet.phi() );
    //       Book( TH1F( "Jet_e", "Jet e [GeV]", 100, 0.0, 1000 ), "AK4_low" )->Fill( myjet.e() );
    //
    //       Book( TH1F( "Jet_muf", "Jet muf", 50, 0.0, 1.0 ), "AK4_low" )->Fill( myjet.muf() );
    //       Book( TH1F( "Jet_phf", "Jet phf", 50, 0.0, 1.0 ), "AK4_low" )->Fill( myjet.phf() );
    //       Book( TH1F( "Jet_emf", "Jet emf", 50, 0.0, 1.0 ), "AK4_low" )->Fill( myjet.emf() );
    //       Book( TH1F( "Jet_nhf", "Jet nhf", 50, 0.0, 1.0 ), "AK4_low" )->Fill( myjet.nhf() );
    //       Book( TH1F( "Jet_chf", "Jet chf", 50, 0.0, 1.0 ), "AK4_low" )->Fill( myjet.chf() );
    //       Book( TH1F( "Jet_area", "Jet area", 40, 0.0, 4.0 ), "AK4_low" )->Fill( myjet.area() );
    //       Book( TH1F( "Jet_cm", "Jet cm", 50, 0.0, 50.0 ), "AK4_low" )->Fill( myjet.cm() );
    //       Book( TH1F( "Jet_nm", "Jet nm", 50, 0.0, 50.0 ), "AK4_low" )->Fill( myjet.nm() );
    //     }
    //   }
    // }
    //
  
  
  } // 1 fat jets
  
  // book-keeping
  Hist( "Events" )->Fill( 0., b_weightGen        ); // event with MC weight
  Hist( "Events" )->Fill( 1,  b_weight           ); // event total weight
  Hist( "Events" )->Fill( 2,  b_weightPU         ); // event pileup weight
  Hist( "Events" )->Fill( 3,  b_weightBtag ); // event b-tag SF weight
  // Hist( "Events" )->Fill( 3,  b_weight_elScale   ); // event electron SF weight
  // Hist( "Events" )->Fill( 4,  b_weight_muScale   ); // event muon SF weight
  Hist( "Events" )->Fill( 9,  1                  ); // event without MC weight
  Hist( "SumEvents" )->Fill( 0., fabs(b_weightGen)       ); // event with MC weight
  Hist( "SumEvents" )->Fill( 1,  fabs(b_weight)          ); // event total weight
  Hist( "SumEvents" )->Fill( 2,  fabs(b_weightPU)       ); // event pileup weight
  Hist( "SumEvents" )->Fill( 3,  fabs(b_weightBtag)); // event b-tag SF weight
  // Hist( "SumEvents" )->Fill( 3,  fabs(b_weight_elScale)  ); // event electron SF weight
  // Hist( "SumEvents" )->Fill( 4,  fabs(b_weight_muScale)  ); // event muon SF weight
  // Hist( "SumEvents" )->Fill( 6,  fabs(b_weight_jvfScale));  // event JVF SF weight
  Hist( "SumEvents" )->Fill( 9,  1                       ); // event without MC weight
  
 
  
  bool foundTau=false;
  int goodTauIndex=0;
  std::vector<UZH::Tau> goodTaus;
  for ( int i = 0; i <   (m_tau.N)
	  ; ++i ) {
    UZH::Tau mytau( &m_tau, i );
    if (mytau.pt() > m_tauPtCut){
      if (fabs(mytau.eta()) < m_tauEtaCut){

	 std::cout<<" mytau.pt() " << mytau.pt() <<" mytau.eta() " << mytau.eta() <<" mytau.TauType() " <<  mytau.TauType() <<std::endl;
	
	goodTaus.push_back(mytau);
          
	foundTau=true;
	std::cout<<" index i " << i <<std::endl;
	goodTauIndex=i;
	std::cout<<" address " << goodTaus[goodTauIndex] << std::endl;
	TLorentzVector tau_tlv = goodTaus[goodTauIndex].tlv();

	std::cout<<" goodTauIndex " <<goodTauIndex<<std::endl;
	std::cout<<"tau_tlv pt" << tau_tlv.Pt()<<std::endl;
      }
    }
  }
  
  if (!foundTau) return;
  for (unsigned int s=0;s<m_catNames.size();++s) {
    if (m_catNames[s].find("tau") != std::string::npos) { 
      selectionBits[s].SetBitNumber( kTau );
      selectionBits[s].SetBitNumber( kTauIsolation );
    }
  }
  bool foundMet=false;
 // Cut : select met
  UZH::MissingEt goodMet( &m_missingEt, 0 );
  if (goodMet.et() > m_metCut) {
    std::cout<<"MET " << goodMet.et() <<std::endl;
    foundMet=true;
  }
  for (unsigned int s=0;s<m_catNames.size();++s) {
    if (m_catNames[s].find("tau") != std::string::npos) { 
      selectionBits[s].SetBitNumber( kMet );
    }
  }
   if (!foundMet) return;

  bool foundSecondLepton=false;
  bool foundElectron=false;
  std::vector<UZH::Electron> goodElectrons;
  for ( int i = 0; i <   (m_electron.N)
	  ; ++i ) {
    UZH::Electron myelectron( &m_electron, i );
    if (myelectron.pt() > m_electronPtCut) {
      if (fabs(myelectron.eta()) < m_electronEtaCut) {

        if (myelectron.isLooseElectronBoosted()==1 || myelectron.isMediumElectronBoosted()==1 ){
	  
	   std::cout<<" myelectron.pt() " << myelectron.pt() <<" myelectron.eta() " << myelectron.eta() <<" myelectron.isLooseElectronBoosted() "<<  myelectron.isLooseElectronBoosted() <<" myelectron.isMediumElectronBoosted()"<<myelectron.isMediumElectronBoosted() <<std::endl;
	  
	  goodElectrons.push_back(myelectron);
          
	  foundElectron=true;
	}

      }
    }
  }
  if (foundElectron){ 
    for (unsigned int s=0;s<m_catNames.size();++s) {

      TLorentzVector ele_tlv = goodElectrons[0].tlv();
      TLorentzVector tau_tlv = goodTaus[0].tlv();
      if ( ele_tlv.DeltaR(tau_tlv) >0.1 && 
	   m_catNames[s].find("ele") != std::string::npos) {
	selectionBits[s].SetBitNumber( kLepton );
      }
    }
  }
  bool foundMuon=false;
  std::vector<UZH::Muon> goodMuons;
  for ( int i = 0; i <   (m_muon.N)
  	  ; ++i ) {
    UZH::Muon mymuon( &m_muon, i );
    if (mymuon.pt() > m_muonPtCut) {
      if (fabs(mymuon.eta()) < m_muonEtaCut) { 
	if (mymuon.isLooseMuon()==1){
	  std::cout<<" mymuon.pt() " << mymuon.pt() <<" mymuon.eta() " << mymuon.eta() <<" mymuon.isLooseMuon() "<<  mymuon.isLooseMuon() <<std::endl;
	  
	  goodMuons.push_back(mymuon);
          
	  foundMuon=true;

	}
      }
    }
  }
  for (unsigned int s=0;s<m_catNames.size();++s) {
    if (foundMuon){ 
      TLorentzVector muon_tlv = goodMuons[0].tlv();
      TLorentzVector tau_tlv = goodTaus[0].tlv();
      if( muon_tlv.DeltaR(tau_tlv) >0.1 && m_catNames[s].find("mu") != std::string::npos) { 
      selectionBits[s].SetBitNumber( kLepton );
      }
    }
  }
  for (unsigned int s=0;s<m_catNames.size();++s) {
    if (goodTaus.size()>1   // && goodTaus[goodTauIndex].tlv().DeltaR( goodTaus[goodTauIndex+1].tlv()) >0.1
			      && m_catNames[s].find("tautau") != std::string::npos) {
      selectionBits[s].SetBitNumber( kLepton );
    }
  }
  // if (foundMuon){
  //   if (goodTaus[goodTauIndex].tlv().DeltaR(goodMuons[0].tlv())>0.1) foundSecondLepton=true;
  // }
  // if (foundElectron){
  //   if (goodTaus[goodTauIndex].tlv().DeltaR(goodElectrons[0].tlv())>0.1) foundSecondLepton=true;
  // }
  // if (goodTaus.size()>1){
   
  //   if (goodTaus[goodTauIndex].tlv().DeltaR(goodTaus[goodTauIndex + 1].tlv())>0.1) foundSecondLepton=true;
  // }


  



  
  for (unsigned int s=0;s<m_catNames.size();++s) {
    // m_logger << INFO << selectionBits[s].TestBitNumber( kTrigger ) << SLogger::endmsg;
    fillCutflow("cutflow", m_catNames[s], selectionBits[s], b_weight);
  }
  
  bool doHistograms = false;
  // need vectorJet and higgsJet to be defined
  for (unsigned int s=0;s<m_catNames.size();++s) {
    if (selectionBits[s].TestBitNumber( kMassWindow )) {
      doHistograms = true;
    }
  }
  
  m_logger << VERBOSE << "before doHistograms" << SLogger::endmsg;
  
  if (doHistograms) {
    // calculate a few variables before filling histograms
    double vJet_tau21 = -1;
    double vJet_tau31 = -1;
    double vJet_tau32 = -1;
    if (Jet.tau1() != 0) {
      vJet_tau21 = Jet.tau2()/Jet.tau1();
      vJet_tau31 = Jet.tau3()/Jet.tau1();
    }
    if (Jet.tau2() != 0) {
      vJet_tau32 = Jet.tau3()/Jet.tau2();
    }
    int vJet_nTaggedSubjets = 0;
    double vJet_subjet0_csv = -99;
    double vJet_subjet1_csv = -99;

    
    for (int i = 0; i < Jet.subjet_pruned_N(); ++i) {
      switch(i) {
        case 0:
          vJet_subjet0_csv = Jet.subjet_pruned_csv()[i];
          break;
        case 1:
          vJet_subjet1_csv = Jet.subjet_pruned_csv()[i];
          break;
      }
     
      if ((m_bTaggingScaleTool.isTagged(Jet.subjet_pruned_csv()[i]))) {
        ++vJet_nTaggedSubjets;
      }
    }

    // fill cut info for ntuple
  std::vector<bool> passed_all(m_catNames.size(), true);
  for (unsigned int s=0;s<m_catNames.size();++s) {
    for( UInt_t icut = 0; icut < kNumCuts; ++icut ) {
      if( selectionBits[s].TestBitNumber( icut ) != kTRUE ){
        passed_all[s] = false;
      }
      else{
        if (icut) b_selection_bits[s]|=1<<icut; 
        if (icut-1==(unsigned)b_selection_lastcut[s])
          b_selection_lastcut[s]++;
      }
    }//cut loop
  }//category loop



    double deta = Jet.eta();
    double dphi = Jet.phi();
    double dr = sqrt(deta*deta + dphi*dphi);
    TLorentzVector diJet = Jet.tlv()   ;
   
    m_logger << VERBOSE << "category loopfillHistograms" << SLogger::endmsg;

    for (unsigned int s=0;s<m_catNames.size();++s) {
      if (passed_all[s]) {
        m_logger << VERBOSE << m_catNames[s] << SLogger::endmsg;
	if (m_catNames[s].find("NoWindow") != std::string::npos) {

	  checkTrigger(m_catNames[s]);
	}
	fillHistograms(m_catNames[s], Jet, Jet, diJet, vJet_tau21, vJet_tau31, vJet_tau32, vJet_nTaggedSubjets, vJet_subjet0_csv, vJet_subjet1_csv, vJet_tau21, vJet_tau31, vJet_tau32, vJet_nTaggedSubjets, vJet_subjet0_csv, vJet_subjet1_csv, deta, dphi, dr );
      }
    }
  }
  m_logger << VERBOSE << "return" << SLogger::endmsg;

  return;

}

bool VHTausAnalysis::isGoodEvent(int runNumber, int lumiSection) {
  
  bool isGood = true;
  if (m_isData) {
    isGood = m_grl.HasRunLumiBlock( runNumber, lumiSection );
    if( !isGood ) {
      m_logger << WARNING << "Bad event! Run: " << runNumber <<  " - Lumi Section: " << lumiSection << SLogger::endmsg;
      // throw SError( SError::SkipEvent );
    }
    else m_logger << VERBOSE << "Good event! Run: " << runNumber <<  " - Lumi Section: " << lumiSection << SLogger::endmsg;
  }
  return isGood;
  
}


bool VHTausAnalysis::passTrigger() {
  
  bool passTrigger = false;
  
  for (std::map<std::string,bool>::iterator it = (m_eventInfo.trigDecision)->begin(); it != (m_eventInfo.trigDecision)->end(); ++it){
    for (unsigned int t = 0; t < m_triggerNames.size(); ++t ){
      if ((it->first).find(m_triggerNames[t]) != std::string::npos) {
        if (it->second == true) {
          m_logger << VERBOSE << "Trigger pass: " << (it->first) << SLogger::endmsg;
          passTrigger = true;
          return passTrigger;
        }
      }
    }
  }
  
  return passTrigger;
  
}


bool VHTausAnalysis::passMETFilters() {
  
  bool passMetFilters = true;
  
  // using only what's recommended in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    
  if( !(m_eventInfo.PV_filter) ) {
    passMetFilters = false;
    m_logger << VERBOSE << "PV_filter" << SLogger::endmsg;
    Hist( "METFilters" )->Fill(1);
  }
  if( !(m_eventInfo.passFilter_CSCHalo) ) {
    passMetFilters = false;
    m_logger << VERBOSE << "passFilter_CSCHalo" << SLogger::endmsg;
    Hist( "METFilters" )->Fill(2);
  }
  if( !(m_eventInfo.passFilter_HBHELoose) ) {
    passMetFilters = false;
    m_logger << VERBOSE << "passFilter_HBHELoose" << SLogger::endmsg;
    Hist( "METFilters" )->Fill(3);
  }
  if( !(m_eventInfo.passFilter_HBHEIso) ) {
    passMetFilters = false;
    m_logger << VERBOSE << "passFilter_HBHEIso" << SLogger::endmsg;
    Hist( "METFilters" )->Fill(4);
  }
  if( !(m_eventInfo.passFilter_EEBadSc) ) {
    passMetFilters = false;
    m_logger << VERBOSE << "passFilter_EEBadSc" << SLogger::endmsg;
    Hist( "METFilters" )->Fill(5);
  }
  
  return passMetFilters;
  
}


double VHTausAnalysis::getEventWeight() {
  
  double weight = 1.;
  for( unsigned int v = 0; v < (m_eventInfo.actualIntPerXing)->size(); ++v ){
    
    if ( (*m_eventInfo.bunchCrossing)[v] == 0 ) {
      b_weightPU = m_pileupReweightingTool.getPileUpweight( (*m_eventInfo.actualIntPerXing)[v] );
      m_logger << VERBOSE << "Weight: " << b_weightPU << " for true: " << (*m_eventInfo.actualIntPerXing)[v] << SLogger::endmsg;
     
      break;
    }
  }
  b_weightGen = (m_eventInfo.genEventWeight < 0) ? -1 : 1; 
  weight *= b_weightPU*b_weightGen;
  
  return weight;
  
}

void VHTausAnalysis::clearBranches() {
  
  b_weight = 1.;
  b_weightGen = 1.;
  b_weightPU = 1.;
  b_weightBtag = 1.;
  
  b_runNumber = -99;;
  b_eventNumber = -99;
  b_lumiBlock = -99;
  
  b_ak8jet0_pt = -99;
  b_ak8jet0_phi = -99;
  b_ak8jet0_eta = -99;
  b_ak8jet0_e = -99;
  b_ak8jet0_tau21 = -99;
  b_ak8jet0_m = -99;
  b_ak8jet0_mpruned = -99;
  b_ak8jet0_csv = -99;
  b_ak8jet1_pt = -99;
  b_ak8jet1_phi = -99;
  b_ak8jet1_eta = -99;
  b_ak8jet1_e = -99;
  b_ak8jet1_tau21 = -99;
  b_ak8jet1_m = -99;
  b_ak8jet1_mpruned = -99;
  b_ak8jet1_csv = -99;
  b_ak8jet0_subjet01_dr = -99;
  b_ak8jet0_subjet01_deta = -99;
  b_ak8jet0_subjet01_dphi = -99;
  b_ak8jet0_subjet0_pt = -99;
  b_ak8jet0_subjet1_pt = -99;
  b_ak8jet0_subjet0_csv = -99;
  b_ak8jet0_subjet1_csv = -99;
  b_ak8jet1_subjet01_dr = -99;
  b_ak8jet1_subjet01_deta = -99;
  b_ak8jet1_subjet01_dphi = -99;
  b_ak8jet1_subjet0_pt = -99;
  b_ak8jet1_subjet1_pt = -99;
  b_ak8jet1_subjet0_csv = -99;
  b_ak8jet1_subjet1_csv = -99;
  
  b_selection_bits.clear();
  b_selection_lastcut.clear();
  
}

void VHTausAnalysis::fillCutflow( const std::string histName, const std::string dirName, const TBits& cutmap, const Double_t weight ) {

  // bool writeNtuple = false;
  // sequential cut flow -> stop at first failed cut
  // m_logger << INFO << histName << "\t" << dirName << SLogger::endmsg;
  for( UInt_t i = 0; i < cutmap.GetNbits(); ++i ) {
    // m_logger << INFO << i << ":\t" << cutmap.TestBitNumber( i ) << SLogger::endmsg;
    if( cutmap.TestBitNumber( i ) ) {
      Hist( histName.c_str(), dirName.c_str() )->Fill( i+1, weight );
      // if (i == (unsigned int) m_ntupleLevel) {
      //   writeNtuple = true;
      // }
    } else {
      break;
    }
  }
  // if (!writeNtuple) {
  // // this does something really bad...
  //   throw SError( SError::SkipEvent );
  // }
}
void VHTausAnalysis::bookTriggerHistos(const TString& directory  ){
  Book( TH1F( "Trigger", "Trigger",60,0,60  ), directory);
  // theAnalysis_->Book( TH1F( "TriggerMuTau", "TriggerMuTau", 45,0,45));
  // theAnalysis_->Book( TH1F( "TriggerEleTau", "TriggerEleTau", 45,0,45));
  Book( TH1F( "Trigger_Tau", "Trigger_Tau", 60,0,60), directory);
  Book( TH1F( "Trigger_Jet", "Trigger_Jet", 60,0,60), directory);
  // theAnalysis_->Book( TH1F( "Trigger_TauMuTau", "Trigger_TauMuTau", 45,0,45));
  // theAnalysis_->Book( TH1F( "Trigger_JetMuTau", "Trigger_JetMuTau", 45,0,45));
  // theAnalysis_->Book( TH1F( "Trigger_TauEleTau", "Trigger_TauEleTau", 45,0,45));
  // theAnalysis_->Book( TH1F( "Trigger_JetEleTau", "Trigger_JetEleTau", 45,0,45));
  Book( TH1F( "Trigger_TauOrJet", "Trigger_TauOrJet",7,-1.5,5.5  ) , directory);

  std::vector<std::string> trignames;

  trignames.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v");
  trignames.push_back("HLT_AK8PFJet360_TrimMass30_v");
  trignames.push_back("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v");
  trignames.push_back("HLT_AK8PFHT660_TrimR0p1PT0p03Mass50_BTagCSV_p20_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v");
  trignames.push_back("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_v");
  trignames.push_back("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v");
  trignames.push_back("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v");
  trignames.push_back("HLT_PFHT600_v");
  trignames.push_back("HLT_PFHT650_v") ;
  trignames.push_back("HLT_PFHT800_v") ;
  trignames.push_back("HLT_PFHT900_V") ;
  trignames.push_back("HLT_PFHT900_V") ;
  trignames.push_back("HLT_PFJet320_v") ;
  trignames.push_back("HLT_PFJet450_v") ;  
    
  trignames.push_back("HLT_IsoMu20_eta2p1_v") ;
  trignames.push_back("HLT_IsoMu24_eta2p1_v") ;
  trignames.push_back("HLT_Mu45_eta2p1_v") ;
  trignames.push_back("HLT_Ele27_eta2p1_WPLoose_v") ;
  trignames.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v") ;
  trignames.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v") ;
  trignames.push_back("HLT_Ele32_eta2p1_WP75_Gsf_v") ;
  trignames.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v") ;
  trignames.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v") ;
  trignames.push_back("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v") ;
  trignames.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau35_v") ;
  trignames.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v") ;
  trignames.push_back("HLT_IsoMu16_eta2p1_MET30_JetIDCleaned_LooseIsoPFTau50_v") ;
  //H->tautau triggers   
  trignames.push_back("HLT_LooseIsoPFTau50_v") ;
  trignames.push_back("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v") ;
  trignames.push_back("HLT_IsoMu17_eta2p1_v") ;
  trignames.push_back("HLT_IsoMu18_v") ;
  trignames.push_back("HLT_IsoMu27_v") ;
  trignames.push_back("HLT_IsoMu20_v") ;
  trignames.push_back("HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v") ;
  trignames.push_back("HLT_Ele22_eta2p1_WP75_Gsf_v") ;
  trignames.push_back("HLT_Ele32_eta2p1_WP75_Gsf_v") ;
  trignames.push_back("HLT_Ele23_WPLoose_Gsf_v") ;
  trignames.push_back("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v") ;
  trignames.push_back("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v") ;

  //Alternative triggers
  trignames.push_back("HLT_PFMET120_PFMHT120_IDTight_v") ;
  trignames.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v"); 
  trignames.push_back("HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v");

 
 
  std::vector<std::string> trignamesTau;
  //commenting out the prescaled ones
  // trignamesTau.push_back("HLT_IsoMu20_eta2p1_v");
  // trignamesTau.push_back("HLT_IsoMu24_eta2p1_v");
  // trignamesTau.push_back("HLT_Mu45_eta2p1_v");
  // trignamesTau.push_back("HLT_Ele27_eta2p1_WPLoose_v");
  // trignamesTau.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v");
  // trignamesTau.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  // trignamesTau.push_back("HLT_Ele32_eta2p1_WP75_Gsf_v");
  // trignamesTau.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");
  // trignamesTau.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
  // trignamesTau.push_back("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v");
  // trignamesTau.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau35_v");
  // trignamesTau.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v");
  // trignamesTau.push_back("HLT_IsoMu16_eta2p1_MET30_JetIDCleaned_LooseIsoPFTau50_v");  
  //H->tautau triggers   
  // trignamesTau.push_back("HLT_LooseIsoPFTau50_v");
  trignamesTau.push_back("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v");
  trignamesTau.push_back("HLT_IsoMu24_eta2p1_v");
  trignamesTau.push_back("HLT_IsoMu17_eta2p1_v");
  trignamesTau.push_back("HLT_IsoMu18_v");
  trignamesTau.push_back("HLT_IsoMu27_v");
  trignamesTau.push_back("HLT_IsoMu20_v");
  trignamesTau.push_back("HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v");
  trignamesTau.push_back("HLT_Ele22_eta2p1_WP75_Gsf_v");
  trignamesTau.push_back("HLT_Ele32_eta2p1_WP75_Gsf_v");
  trignamesTau.push_back("HLT_Ele23_WPLoose_Gsf_v");
  trignamesTau.push_back("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v");
  trignamesTau.push_back("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v");

  //Alternative triggers
  trignamesTau.push_back("HLT_PFMET120_PFMHT120_IDTight_v");
  trignamesTau.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v");
  trignamesTau.push_back("HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v");

 
  std::vector<std::string> trignamesJet;
  trignamesJet.push_back("HLT_AK8PFJet360_TrimMass30_v") ;
  trignamesJet.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v") ;
  trignamesJet.push_back("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v") ;
  trignamesJet.push_back("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v") ;
  trignamesJet.push_back("HLT_PFHT800_v") ;
  



  const char *TriggerNames[trignames.size()];
  const char *TriggerTauNames[trignamesTau.size()];
  const char *TriggerJetNames[trignamesJet.size()];




  Hist( "Trigger"  	, directory )->GetXaxis()->SetTitle( "Triggers" );
  Hist( "Trigger"  	, directory )->GetXaxis()->LabelsOption("v");
  Hist( "Trigger_Tau"  	, directory )->GetXaxis()->SetTitle( "Triggers" );
  Hist( "Trigger_Tau"  	, directory )->GetXaxis()->LabelsOption("v");
  Hist( "Trigger_Jet"  	, directory )->GetXaxis()->SetTitle( "Triggers" );
  Hist( "Trigger_Jet"   , directory )->GetXaxis()->LabelsOption("v");
 



  for ( unsigned int t=0; t< trignames.size(); t++){
    TriggerNames[t]= trignames[t].c_str();
    Hist( "Trigger"  , directory 	)->GetXaxis()->SetBinLabel(t+1, TriggerNames[t] );
  }

  for ( unsigned int t=0; t< trignamesTau.size(); t++){ 
    TriggerTauNames[t]= trignamesTau[t].c_str();
    Hist( "Trigger_Tau"  , directory 	)->GetXaxis()->SetBinLabel(t+1, TriggerTauNames[t] );
    }

  for ( unsigned int t=0; t< trignamesJet.size(); t++){ 
    TriggerJetNames[t]= trignamesJet[t].c_str();
    Hist( "Trigger_Jet"  , directory 	)->GetXaxis()->SetBinLabel(t+1, TriggerJetNames[t] );
  }

}


void VHTausAnalysis::bookHistograms( const TString& directory ) {
  
  // kinematics histograms
  Book( TH1F( "vjet_pt", "Vjet p_{T};Vjet p_{T} [GeV]", 200, 0, 2000 ), directory ); 
  Book( TH1F( "vjet_eta", "Vjet #eta;Vjet #eta", 50, -2.5, 2.5 ), directory ); 
  Book( TH1F( "vjet_phi", "Vjet #phi;Vjet #phi", 50, -3.15, 3.15 ), directory ); 
  Book( TH1F( "vjet_m", "Vjet m;Vjet m [GeV]", 40, 0, 200 ), directory ); 
  Book( TH1F( "vjet_mpruned", "Vjet pruned m;Vjet pruned m [GeV]", 40, 0, 200 ), directory ); 
  Book( TH1F( "vjet_tau1", "Vjet #tau_{1};Vjet #tau_{1}", 50, 0, 1 ), directory ); 
  Book( TH1F( "vjet_tau2", "Vjet #tau_{2};Vjet #tau_{2}", 50, 0, 1 ), directory ); 
  Book( TH1F( "vjet_tau3", "Vjet #tau_{3};Vjet #tau_{3}", 50, 0, 1 ), directory ); 
  Book( TH1F( "vjet_tau21", "Vjet #tau_{21};Vjet #tau_{21}", 50, 0, 1 ), directory ); 
  Book( TH1F( "vjet_tau31", "Vjet #tau_{31};Vjet #tau_{31}", 50, 0, 1 ), directory ); 
  Book( TH1F( "vjet_tau32", "Vjet #tau_{32};Vjet #tau_{32}", 50, 0, 1 ), directory ); 
  Book( TH1F( "vjet_nSubjets", "Vjet N subjets;Vjet N subjets", 10, -.5, 9.5 ), directory ); 
  Book( TH1F( "vjet_nTaggedSubjets", "Vjet N tagged subjets;Vjet N tagged subjets", 10, -.5, 9.5 ), directory ); 
  Book( TH1F( "vjet_subjet0_csv", "Vjet subjet 0 CSV;Vjet subjet0 CSV", 50, 0, 1 ), directory ); 
  Book( TH1F( "vjet_subjet1_csv", "Vjet subjet 1 CSV;Vjet subjet1 CSV", 50, 0, 1 ), directory ); 

  Book( TH1F( "hjet_pt", "Hjet p_{T};Hjet p_{T} [GeV]", 200, 0, 2000 ), directory ); 
  Book( TH1F( "hjet_eta", "Hjet #eta;Hjet #eta", 50, -2.5, 2.5 ), directory ); 
  Book( TH1F( "hjet_phi", "Hjet #phi;Hjet #phi", 50, -3.15, 3.15 ), directory ); 
  Book( TH1F( "hjet_m", "Hjet m;Hjet m [GeV]", 40, 0, 200 ), directory ); 
  Book( TH1F( "hjet_mpruned", "Hjet pruned m;Hjet pruned m [GeV]", 40, 0, 200 ), directory ); 
  Book( TH1F( "hjet_tau1", "Hjet #tau_{1};Hjet #tau_{1}", 50, 0, 1 ), directory ); 
  Book( TH1F( "hjet_tau2", "Hjet #tau_{2};Hjet #tau_{2}", 50, 0, 1 ), directory ); 
  Book( TH1F( "hjet_tau3", "Hjet #tau_{3};Hjet #tau_{3}", 50, 0, 1 ), directory ); 
  Book( TH1F( "hjet_tau21", "Hjet #tau_{21};Hjet #tau_{21}", 50, 0, 1 ), directory ); 
  Book( TH1F( "hjet_tau31", "Hjet #tau_{31};Hjet #tau_{31}", 50, 0, 1 ), directory ); 
  Book( TH1F( "hjet_tau32", "Hjet #tau_{32};Hjet #tau_{32}", 50, 0, 1 ), directory ); 
  Book( TH1F( "hjet_nSubjets", "Hjet N subjets;Hjet N subjets", 10, -.5, 9.5 ), directory ); 
  Book( TH1F( "hjet_nTaggedSubjets", "Hjet N tagged subjets;Hjet N tagged subjets", 10, -.5, 9.5 ), directory ); 
  Book( TH1F( "hjet_subjet0_csv", "Hjet subjet 0 CSV;Hjet subjet0 CSV", 100, -1, 1 ), directory ); 
  Book( TH1F( "hjet_subjet1_csv", "Hjet subjet 1 CSV;Hjet subjet1 CSV", 100, -1, 1 ), directory ); 

  Book( TH1F( "jets_deta", "jets #Delta #eta;jets #Delta #eta", 50, 0, 5 ), directory ); 
  Book( TH1F( "jets_dphi", "jets #Delta #phi;jets #Delta #phi", 50, 0, 6.3 ), directory ); 
  Book( TH1F( "jets_dr", "jets #Delta R;jets #Delta R", 50, 0, 5 ), directory ); 

  Book( TH1F( "dijet_pt", "dijet p_{T};dijet p_{T} [GeV]", 100, 0, 1000 ), directory ); 
  Book( TH1F( "dijet_eta", "dijet #eta;dijet #eta", 50, -2.5, 2.5 ), directory ); 
  Book( TH1F( "dijet_phi", "dijet #phi;Vdijet #phi", 50, -3.15, 3.15 ), directory ); 
  Book( TH1F( "dijet_m", "dijet m;dijet m [GeV]", 400, 0, 4000 ), directory );
  Book( TH1F( "dijet_template_m", "dijet m;dijet m [GeV]", 7000, 0, 7000 ), directory ); 
  
}



void VHTausAnalysis::fillHistograms( const TString& directory, const UZH::Jet& vectorJet, const UZH::Jet& higgsJet, const TLorentzVector& diJet, const double& vJet_tau21, const double& vJet_tau31, const double& vJet_tau32, const int& vJet_nTaggedSubjets, const double& vJet_subjet0_csv, const double& vJet_subjet1_csv, const double& hJet_tau21, const double& hJet_tau31, const double& hJet_tau32, const int& hJet_nTaggedSubjets, const double& hJet_subjet0_csv, const double& hJet_subjet1_csv, const double& deta, const double& dphi, const double& dr ) {
  
  // fill all histograms
  Hist( "vjet_pt", directory )->Fill( vectorJet.pt() , b_weight);
  Hist( "vjet_eta", directory )->Fill( vectorJet.eta() , b_weight);
  Hist( "vjet_phi", directory )->Fill( vectorJet.phi() , b_weight);
  Hist( "vjet_m", directory )->Fill( vectorJet.m() , b_weight);
  Hist( "vjet_mpruned", directory )->Fill( vectorJet.pruned_massCorr() , b_weight);
  Hist( "vjet_tau1", directory )->Fill( vectorJet.tau1() , b_weight);
  Hist( "vjet_tau2", directory )->Fill( vectorJet.tau2() , b_weight);
  Hist( "vjet_tau3", directory )->Fill( vectorJet.tau3() , b_weight);
  Hist( "vjet_tau21", directory )->Fill( vJet_tau21 , b_weight);
  Hist( "vjet_tau31", directory )->Fill( vJet_tau31 , b_weight);
  Hist( "vjet_tau32", directory )->Fill( vJet_tau32 , b_weight);
  Hist( "vjet_nSubjets", directory )->Fill( vectorJet.subjet_pruned_N() , b_weight);
  Hist( "vjet_nTaggedSubjets", directory )->Fill( vJet_nTaggedSubjets , b_weight);
  Hist( "vjet_subjet0_csv", directory )->Fill( vJet_subjet0_csv , b_weight);
  Hist( "vjet_subjet1_csv", directory )->Fill( vJet_subjet1_csv , b_weight);

  Hist( "hjet_pt", directory )->Fill( higgsJet.pt() , b_weight);
  Hist( "hjet_eta", directory )->Fill( higgsJet.eta() , b_weight);
  Hist( "hjet_phi", directory )->Fill( higgsJet.phi() , b_weight);
  Hist( "hjet_m", directory )->Fill( higgsJet.m() , b_weight);
  Hist( "hjet_mpruned", directory )->Fill( higgsJet.pruned_massCorr() , b_weight);
  Hist( "hjet_tau1", directory )->Fill( higgsJet.tau1() , b_weight);
  Hist( "hjet_tau2", directory )->Fill( higgsJet.tau2() , b_weight);
  Hist( "hjet_tau3", directory )->Fill( higgsJet.tau3() , b_weight);
  Hist( "hjet_tau21", directory )->Fill( hJet_tau21 , b_weight);
  Hist( "hjet_tau31", directory )->Fill( hJet_tau31 , b_weight);
  Hist( "hjet_tau32", directory )->Fill( hJet_tau32 , b_weight);
  Hist( "hjet_nSubjets", directory )->Fill( higgsJet.subjet_pruned_N() , b_weight);
  Hist( "hjet_nTaggedSubjets", directory )->Fill( hJet_nTaggedSubjets , b_weight);
  Hist( "hjet_subjet0_csv", directory )->Fill( hJet_subjet0_csv , b_weight);
  Hist( "hjet_subjet1_csv", directory )->Fill( hJet_subjet1_csv , b_weight);

  Hist( "jets_deta", directory )->Fill( deta , b_weight);
  Hist( "jets_dphi", directory )->Fill( dphi , b_weight);
  Hist( "jets_dr", directory )->Fill( dr , b_weight);

  Hist( "dijet_pt", directory )->Fill( diJet.Pt() , b_weight);
  Hist( "dijet_eta", directory )->Fill( diJet.Eta() , b_weight);
  Hist( "dijet_phi", directory )->Fill( diJet.Phi() , b_weight);
  Hist( "dijet_m", directory )->Fill( diJet.M() , b_weight);
  Hist( "dijet_template_m", directory )->Fill( diJet.M() , b_weight);
  
}


void VHTausAnalysis::checkTrigger(const TString& directory){
  
  std::vector<std::string> trignames;

  trignames.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v");
  trignames.push_back("HLT_AK8PFJet360_TrimMass30_v");
  trignames.push_back("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v");
  trignames.push_back("HLT_AK8PFHT660_TrimR0p1PT0p03Mass50_BTagCSV_p20_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v");
  trignames.push_back("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_v");

  trignames.push_back("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v");

  trignames.push_back("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v");
  trignames.push_back("HLT_PFHT600_v");
  trignames.push_back("HLT_PFHT650_v") ;
  trignames.push_back("HLT_PFHT800_v") ;
  trignames.push_back("HLT_PFHT900_V") ;
  trignames.push_back("HLT_PFHT900_V") ;
  trignames.push_back("HLT_PFJet320_v") ;
  trignames.push_back("HLT_PFJet450_v") ;  
    
  trignames.push_back("HLT_IsoMu20_eta2p1_v") ;
  trignames.push_back("HLT_IsoMu24_eta2p1_v") ;
  trignames.push_back("HLT_Mu45_eta2p1_v") ;
  trignames.push_back("HLT_Ele27_eta2p1_WPLoose_v") ;
  trignames.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v") ;
  trignames.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v") ;
  trignames.push_back("HLT_Ele32_eta2p1_WP75_Gsf_v") ;
  trignames.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v") ;
  trignames.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v") ;
  trignames.push_back("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v") ;
  trignames.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau35_v") ;
  trignames.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v") ;
  trignames.push_back("HLT_IsoMu16_eta2p1_MET30_JetIDCleaned_LooseIsoPFTau50_v") ;
  //H->tautau triggers   
  trignames.push_back("HLT_LooseIsoPFTau50_v") ;
  trignames.push_back("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v") ;
  trignames.push_back("HLT_IsoMu17_eta2p1_v") ;
  trignames.push_back("HLT_IsoMu18_v") ;
  trignames.push_back("HLT_IsoMu27_v") ;
  trignames.push_back("HLT_IsoMu20_v") ;
  trignames.push_back("HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v") ;
  trignames.push_back("HLT_Ele22_eta2p1_WP75_Gsf_v") ;
  trignames.push_back("HLT_Ele32_eta2p1_WP75_Gsf_v") ;
  trignames.push_back("HLT_Ele23_WPLoose_Gsf_v") ;
  trignames.push_back("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v") ;
  trignames.push_back("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v") ;

  //Alternative triggers
  trignames.push_back("HLT_PFMET120_PFMHT120_IDTight_v") ;
  trignames.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v"); 
  trignames.push_back("HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v");

 
 
  std::vector<std::string> trignamesTau;
  //commenting out the prescaled ones
  // trignamesTau.push_back("HLT_IsoMu20_eta2p1_v");
  // trignamesTau.push_back("HLT_IsoMu24_eta2p1_v");
  // trignamesTau.push_back("HLT_Mu45_eta2p1_v");
  // trignamesTau.push_back("HLT_Ele27_eta2p1_WPLoose_v");
  // trignamesTau.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v");
  // trignamesTau.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  // trignamesTau.push_back("HLT_Ele32_eta2p1_WP75_Gsf_v");
  // trignamesTau.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");
  // trignamesTau.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
  // trignamesTau.push_back("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v");
  // trignamesTau.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau35_v");
  // trignamesTau.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v");
  // trignamesTau.push_back("HLT_IsoMu16_eta2p1_MET30_JetIDCleaned_LooseIsoPFTau50_v");  
  //H->tautau triggers   
  // trignamesTau.push_back("HLT_LooseIsoPFTau50_v");
  trignamesTau.push_back("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v");
  trignamesTau.push_back("HLT_IsoMu24_eta2p1_v");
  trignamesTau.push_back("HLT_IsoMu17_eta2p1_v");
  trignamesTau.push_back("HLT_IsoMu18_v");
  trignamesTau.push_back("HLT_IsoMu27_v");
  trignamesTau.push_back("HLT_IsoMu20_v");
  trignamesTau.push_back("HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v");
  trignamesTau.push_back("HLT_Ele22_eta2p1_WP75_Gsf_v");
  trignamesTau.push_back("HLT_Ele32_eta2p1_WP75_Gsf_v");
  trignamesTau.push_back("HLT_Ele23_WPLoose_Gsf_v");
  trignamesTau.push_back("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v");
  trignamesTau.push_back("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v");

  //Alternative triggers
  trignamesTau.push_back("HLT_PFMET120_PFMHT120_IDTight_v");
  trignamesTau.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v");
  trignamesTau.push_back("HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v");

 
  std::vector<std::string> trignamesJet;
  trignamesJet.push_back("HLT_AK8PFJet360_TrimMass30_v") ;
  trignamesJet.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v") ;
  // trignames.push_back("AK8DiPFJet280_200_TrimMass30_BTagCSV0p45") ;
  trignamesJet.push_back("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v") ;
  trignamesJet.push_back("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v") ;
  trignamesJet.push_back("HLT_PFHT800_v") ;
  

  
  std::string GENChannel_="General";
  
  bool isfired = false;
 
  if (GENChannel_=="General"){
    bool isfiredTau = false;
    bool isfiredJet = false;
 
    for( std::map<std::string,bool>::iterator it = (m_eventInfo.trigDecision)->begin(); it != (m_eventInfo.trigDecision)->end(); ++it){
      //   std::cout <<"TriggerFired : "<< (it->first) <<std::endl; 

      for( unsigned int t = 0; t < trignames.size(); ++t ){
      
	if( (it->first).find(trignames[t]) != std::string::npos){
	  // std::cout <<"TriggerFired : "<< (it->first) <<" and number "<<t<< std::endl; 
	  if(it->second ){
	    isfired = true;
	    Hist( "Trigger",directory )->Fill( t ); 
	  }

	  // else std::cout <<"Trigger NOT Fired : "<< trignames[t] <<std::endl;  	
	}
     
      }

      for( unsigned int t = 0; t < trignamesJet.size(); ++t ){
	if( (it->first).find(trignamesJet[t]) != std::string::npos &&  it->second ){
	  isfiredJet = true;
	  Hist( "Trigger_Jet" ,directory)->Fill( t ); 
	}
      }
      for( unsigned int t = 0; t < trignamesTau.size(); ++t ){
	if( (it->first).find(trignamesTau[t]) != std::string::npos && it->second  // && !isfiredJet
	    ){
	  isfiredTau = true;
	  Hist( "Trigger_Tau",directory )->Fill( t ); 
	}
      }
    
    

    }
    Hist( "Trigger_TauOrJet" ,directory)->Fill( -1 ); 	 
    if (isfiredJet || isfiredTau)    Hist( "Trigger_TauOrJet" ,directory)->Fill( 0 ); 	 
    if (isfiredJet )    Hist( "Trigger_TauOrJet" ,directory)->Fill( 1 ); 
    if (isfiredTau)    Hist( "Trigger_TauOrJet" ,directory)->Fill( 2 ); 
    if (isfiredJet && !isfiredTau )    Hist( "Trigger_TauOrJet",directory )->Fill( 3 ); 
    if (isfiredTau && !isfiredJet)    Hist( "Trigger_TauOrJet",directory )->Fill( 4 ); 
    if (!isfiredTau && !isfiredJet)    Hist( "Trigger_TauOrJet",directory )->Fill( 5 ); 
  }
  
  
  
  trignames.clear();
  trignamesJet.clear();
  trignamesTau.clear();
}
