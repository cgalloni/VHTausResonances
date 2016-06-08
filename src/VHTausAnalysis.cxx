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

   /// Declare trees and varibles in trees
 
   channels_.push_back("mutau");
   channels_.push_back("eletau");
   channels_.push_back("tautau");
   
   
   
   for(unsigned int ch=0;ch< channels_.size();ch++){
     m_outputTreeName_ch_.push_back("");
   }   
     
 

   DeclareProperty( "OutputTreeName_mutau" , m_outputTreeName_ch_[0]);
   DeclareProperty( "OutputTreeName_eletau" , m_outputTreeName_ch_[1]);
   DeclareProperty( "OutputTreeName_tautau" , m_outputTreeName_ch_[2]);
 

   DeclareProperty( "NtupleLevel",              m_ntupleLevel              = kTau );
   DeclareProperty( "JetAK4Name",               m_jetAK4Name               = "jetAK4" );
   DeclareProperty( "JetAK8Name",               m_jetAK8Name               = "jetAK8" );
   DeclareProperty( "ElectronName",             m_electronName             = "el" );
   DeclareProperty( "MuonName",                 m_muonName                 = "mu" );
   DeclareProperty( "TauName",                  m_tauName                  = "tau" );
   DeclareProperty( "MissingEtName",            m_missingEtName            = "MET" );
   DeclareProperty( "GenParticleName",          m_genParticleName          = "genParticle" );
   
   DeclareProperty( "IsData",                   m_isData                   = false );
   DeclareProperty( "IsSignal",                 m_isSignal                 = false );
   DeclareProperty( "ApplyMETFilters",          m_applyMETFilters          = true );
   
   DeclareProperty( "JetPtCut",                 m_jetPtCut           = 120. );
   DeclareProperty( "JetEtaCut",                m_jetEtaCut          = 2.4  );
   DeclareProperty( "AK4JetPtCut",              m_AK4jetPtCut           = 40. );
   DeclareProperty( "AK4JetEtaCut",             m_AK4jetEtaCut          = 3.  );
   DeclareProperty( "MjjCut",                   m_mjjCut             = 600  );
   DeclareProperty( "JetDeltaEtaCut",           m_jetDeltaEtaCut     = 1.3  );

   DeclareProperty( "Tau21HPCut",                m_tau21HPCut             = 0.45 );
   DeclareProperty( "Tau21LPCut",                m_tau21LPCut             = 0.75 );
   DeclareProperty( "MVLowSidebandCut",          m_mVLowSidebandCut       = 0 );
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
   DeclareProperty( "TauEtaCut",                m_tauEtaCut          = 2.3  );
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
  // // trignames.push_back("AK8DiPFJet280_200_TrimMass30_BTagCSV0p45") ;
  // m_triggerNames.push_back("PFHT650_WideJetMJJ950DEtaJJ1p5") ;
  // m_triggerNames.push_back("PFHT650_WideJetMJJ900DEtaJJ1p5") ;
  m_triggerNames.push_back("PFHT800_v") ;
  m_triggerNames.push_back("PFHT900_v") ;
  m_triggerNames.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v") ;
  m_triggerNames.push_back("HLT_IsoMu20_v");
  m_triggerNames.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v") ;
  m_triggerNames.push_back("HLT_IsoTkMu20_v");
  m_triggerNames.push_back("HLT_IsoMu18_v");
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
  DeclareVariable(b_weight              , "weight"                 ,  m_outputTreeName.c_str());
  DeclareVariable(b_weightGen           , "weightGen"              ,  m_outputTreeName.c_str());
  DeclareVariable(b_weightPU            , "weightPU"               ,  m_outputTreeName.c_str());
  DeclareVariable(b_weightBtag          , "weightBtag"             ,  m_outputTreeName.c_str());
  
  DeclareVariable( b_runNumber,           "b_runNumber"            ,  m_outputTreeName.c_str());
  DeclareVariable( b_eventNumber,         "b_eventNumber"          ,  m_outputTreeName.c_str());
  DeclareVariable( b_lumiBlock,           "b_lumiBlock"            ,  m_outputTreeName.c_str());
  
  DeclareVariable( b_ak8jet0_pt,           "b_ak8jet0_pt"          ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_phi,          "b_ak8jet0_phi"         ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_eta,          "b_ak8jet0_eta"         ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_e,            "b_ak8jet0_e"           ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_tau21,        "b_ak8jet0_tau21"       ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_m,            "b_ak8jet0_m"           ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_mpruned,      "b_ak8jet0_mpruned"     ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_csv,          "b_ak8jet0_csv"         ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_pt,           "b_ak8jet1_pt"          ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_phi,          "b_ak8jet1_phi"         ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_eta,          "b_ak8jet1_eta"         ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_e,            "b_ak8jet1_e"           ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_tau21,        "b_ak8jet1_tau21"       ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_m,            "b_ak8jet1_m"           ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_mpruned,      "b_ak8jet1_mpruned"     ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_csv,          "b_ak8jet1_csv"         ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_subjet01_dr,   "b_ak8jet0_subjet01_dr"  ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_subjet01_deta, "b_ak8jet0_subjet01_deta",  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_subjet01_dphi, "b_ak8jet0_subjet01_dphi",  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_subjet0_pt,   "b_ak8jet0_subjet0_pt"  ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_subjet1_pt,   "b_ak8jet0_subjet1_pt"  ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_subjet0_csv,  "b_ak8jet0_subjet0_csv" ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet0_subjet1_csv,  "b_ak8jet0_subjet1_csv" ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_subjet01_dr,   "b_ak8jet1_subjet01_dr"  ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_subjet01_deta, "b_ak8jet1_subjet01_deta",  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_subjet01_dphi, "b_ak8jet1_subjet01_dphi",  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_subjet0_pt,   "b_ak8jet1_subjet0_pt"  ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_subjet1_pt,   "b_ak8jet1_subjet1_pt"  ,  m_outputTreeName.c_str());
  DeclareVariable( b_ak8jet1_subjet0_csv,  "b_ak8jet1_subjet0_csv",  m_outputTreeName.c_str() );
  DeclareVariable( b_ak8jet1_subjet1_csv,  "b_ak8jet1_subjet1_csv" ,  m_outputTreeName.c_str());
  
  b_selection_bits.resize( m_catNames.size() );
  b_selection_lastcut.resize( m_catNames.size() );
  for (unsigned int s=0;s<m_catNames.size();++s) {
    DeclareVariable(b_selection_bits[s]    , Form("selection_bits_%s", m_catNames[s].c_str())    ,  m_outputTreeName.c_str());
    DeclareVariable(b_selection_lastcut[s] , Form("selection_lastcut_%s", m_catNames[s].c_str()) ,  m_outputTreeName.c_str());
  }
  



 
  // std::cout << " channels_.size() "<< channels_.size() << std::endl;
    
  for(unsigned int ch=0;ch< channels_.size();ch++){
    // m_outputTreeName_ch_.push_back("tree_"+ channels_[ch]);

    std::cout << " declaring variables for channel " <<  channels_[ch] << " ch " <<ch << std::endl; 
    DeclareVariable(b_weight_[channels_[ch]]              , "weight"                 ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable(b_weightGen_[channels_[ch]]           , "weightGen"              ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable(b_weightPU_[channels_[ch]]            , "weightPU"               ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable(b_weightBtag_[channels_[ch]]          , "weightBtag"             ,  m_outputTreeName_ch_[ch].c_str()); 
    DeclareVariable(b_weightBtag_veto_[channels_[ch]]      , "weightBtag_veto"             ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_runNumber_[channels_[ch]],           "runNumber"            ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_eventNumber_[channels_[ch]],         "eventNumber"          ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_lumiBlock_[channels_[ch]],           "lumiBlock"            ,  m_outputTreeName_ch_[ch].c_str());
  


    DeclareVariable( b_ak8jet0_pt_[channels_[ch]],           "ak8jet0_pt"         ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_ak8jet0_phi_[channels_[ch]],          "ak8jet0_phi"         ,  m_outputTreeName_ch_[ch].c_str() );
    DeclareVariable( b_ak8jet0_eta_[channels_[ch]],          "ak8jet0_eta"          ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_ak8jet0_e_[channels_[ch]],            "ak8jet0_e"            ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_ak8jet0_tau21_[channels_[ch]],        "ak8jet0_tau21"        ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_ak8jet0_m_[channels_[ch]],            "ak8jet0_m"            ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_ak8jet0_mpruned_[channels_[ch]],      "ak8jet0_mpruned"     ,  m_outputTreeName_ch_[ch].c_str() );
    DeclareVariable( b_ak8jet0_csv_[channels_[ch]],          "ak8jet0_csv"         ,  m_outputTreeName_ch_[ch].c_str() );
    DeclareVariable( b_ak8jet0_subjet01_dr_[channels_[ch]],   "ak8jet0_subjet01_dr"  ,  m_outputTreeName_ch_[ch].c_str() );
    DeclareVariable( b_ak8jet0_subjet01_deta_[channels_[ch]], "ak8jet0_subjet01_deta" ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_ak8jet0_subjet01_dphi_[channels_[ch]], "ak8jet0_subjet01_dphi" ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_ak8jet0_subjet0_pt_[channels_[ch]],   "ak8jet0_subjet0_pt"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_ak8jet0_subjet1_pt_[channels_[ch]],   "ak8jet0_subjet1_pt"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_ak8jet0_subjet0_csv_[channels_[ch]],  "ak8jet0_subjet0_csv"  ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_ak8jet0_subjet1_csv_[channels_[ch]],  "ak8jet0_subjet1_csv"  ,  m_outputTreeName_ch_[ch].c_str());
  
    DeclareVariable( b_tau_pt_[channels_[ch]],   "tau_pt"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_tau_eta_[channels_[ch]],   "tau_eta"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_tau_phi_[channels_[ch]],   "tau_phi"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_tauIso_[channels_[ch]],   "tau_iso"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_tauIsoRaw_[channels_[ch]],   "tau_iso_raw"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_tau_againstEle_[channels_[ch]],   "tau_againstEle"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_tau_againstMu_[channels_[ch]],   "tau_againstMu"   ,  m_outputTreeName_ch_[ch].c_str());
  

    DeclareVariable( b_seclepton_pt_[channels_[ch]],   "lep2_pt"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_seclepton_eta_[channels_[ch]],   "lep2_eta"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_seclepton_phi_[channels_[ch]],   "lep2_phi"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_met_[channels_[ch]],   "met"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_met_phi_[channels_[ch]],   "met_phi"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_vis_mass_ll_[channels_[ch]],   "vis_mass_ll"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_dR_ll_[channels_[ch]],   "dR_ll"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_H_Mass_SVFit_[channels_[ch]],   "H_Mass_SVFit"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_H_Pt_SVFit_[channels_[ch]],   "H_Pt_SVFit"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_H_Eta_SVFit_[channels_[ch]],   "H_Eta_SVFit"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_H_Phi_SVFit_[channels_[ch]],   "H_Phi_SVFit"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_X_Mass_SVFit_[channels_[ch]],   "X_Mass_SVFit"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_X_Pt_SVFit_[channels_[ch]],   "X_Pt_SVFit"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_X_Eta_SVFit_[channels_[ch]],   "X_Eta_SVFit"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_X_Phi_SVFit_[channels_[ch]],   "X_Phi_SVFit"   ,  m_outputTreeName_ch_[ch].c_str());

    DeclareVariable( b_trig_HT800_[channels_[ch]], "trig_HT800"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_trig_ak8Jet360_Trim30_[channels_[ch]], "trig_ak8Jet360_Trim30"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_trig_Ele27_[channels_[ch]], "trig_Ele27"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_trig_Ele27_[channels_[ch]], "trig_Ele23"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_trig_Mu20_[channels_[ch]], "trig_Mu20"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_trig_Mu18_[channels_[ch]], "trig_Mu18"   ,  m_outputTreeName_ch_[ch].c_str());
    
    DeclareVariable( b_N_taggedAK4_L_[channels_[ch]], "N_taggedAK4_L"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_N_taggedAK4_M_[channels_[ch]], "N_taggedAK4_M"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_N_taggedAK4_T_[channels_[ch]], "N_taggedAK4_T"   ,  m_outputTreeName_ch_[ch].c_str());
    DeclareVariable( b_isBoostedTau_[channels_[ch]], "isBoostedTau"   ,  m_outputTreeName_ch_[ch].c_str());
 


  }


  //
  // Declare the output histograms:
  //
  Book( TH1F( "Events", "Events;weight", 10, -.5, 9.5 ) );
  Book( TH1F( "SumEvents", "SumEvents;weight", 10, -.5, 9.5 ) );
  Book( TH1F( "METFilters", "METFilters", 20, 0.5, 20.5 ));
  Book( TH1F( "JetPt", "JetPt", 2000, 0, 2000));
  Book( TH1F( "Jet_pruned_massCorr", "Jet_pruned_massCorr", 2000, 0, 2000));
  Book( TH1F( "SumEt","SumEt", 2000, 0, 10000));
  Book( TH1F( "HT","HT", 2000, 0, 10000));
  Book( TH1F( "MET","MET", 2000, 0, 10000));
  Book( TH1F( "N_vertex","N_vertex", 2000, 0, 100));
  Book( TH1F( "N_ak4","N_ak4", 50 , -0.5 ,49.5));
  Book( TH1F( "N_ak8","N_ak8", 50 , -0.5 ,49.5));

  for (unsigned int s=0;s<m_catNames.size();++s) {
    TString directory = m_catNames[s].c_str();
    // cutflow
    Book( TH1F( "cutflow", "cutflow", 20, 0.5, 20.5 ), directory );  
    bookHistograms(directory);
    bookTriggerHistos(directory);
  }
  // std::cout <<" beginning btagging  scale tool " <<std::endl;
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
    m_jetAK4.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis, (m_jetAK4Name + "_").c_str() );
    m_jetAK8.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetSubstructure|Ntuple::JetSoftdropSubjets, (m_jetAK8Name + "_").c_str() );
    m_eventInfo.ConnectVariables(    m_recoTreeName.c_str(), // Ntuple::EventInfoBasic|
				     Ntuple::EventInfoTrigger|Ntuple::EventInfoMETFilters, "" );
  }
  else {
    m_jetAK4.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetTruth, (m_jetAK4Name + "_").c_str() );
    m_jetAK8.ConnectVariables(       m_recoTreeName.c_str(), Ntuple::JetBasic|Ntuple::JetAnalysis|Ntuple::JetSubstructure|Ntuple::JetTruth|Ntuple::JetSoftdropSubjets|Ntuple::JetSoftdropSubjetsTruth, (m_jetAK8Name + "_").c_str() );
    m_eventInfo.ConnectVariables(    m_recoTreeName.c_str(), Ntuple::EventInfoBasic|Ntuple::EventInfoTrigger|Ntuple::EventInfoMETFilters|Ntuple::EventInfoTruth, "" );
    m_genParticle.ConnectVariables(  m_recoTreeName.c_str(), Ntuple::GenParticleBasic, (m_genParticleName + "_").c_str() );
  }
  m_electron.ConnectVariables(     m_recoTreeName.c_str(), Ntuple::ElectronBasic|Ntuple::ElectronID|Ntuple::ElectronBoostedID|Ntuple::ElectronBoostedIsolation, (m_electronName + "_").c_str() );
  m_muon.ConnectVariables(         m_recoTreeName.c_str(), Ntuple::MuonBasic|Ntuple::MuonID|Ntuple::MuonIsolation| Ntuple::MuonBoostedIsolation, (m_muonName + "_").c_str() );
  m_tau.ConnectVariables(         m_recoTreeName.c_str(), Ntuple::TauBasic|Ntuple::TauID|Ntuple::TauAdvancedID, (m_tauName + "_").c_str() );

  m_missingEt.ConnectVariables(    m_recoTreeName.c_str(), Ntuple::MissingEtBasic|Ntuple::MissingEtAnalysis , (m_missingEtName + "_").c_str() );
  
  m_logger << INFO << "Connecting input variables completed" << SLogger::endmsg;

   return;

}

void VHTausAnalysis::ExecuteEvent( const SInputData&, Double_t ) throw( SError ) {

  m_logger << VERBOSE << "ExecuteEvent" << SLogger::endmsg;
  // std::cout << "ExecuteEvent" << std::endl;

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
  /// gen Weigth for MC(

  if (!m_isData) {
    b_weight= getEventWeight();
  }

 
 // Rest of book-keeping
  Hist( "Events" )->Fill( 0., b_weightGen        ); // event with MC weight
  Hist( "Events" )->Fill( 2,  b_weightPU         ); // event pileup weight
  Hist( "Events" )->Fill( 9,  1                  ); // event without MC weight
  Hist( "SumEvents" )->Fill( 0., fabs(b_weightGen)       ); // event with MC weight
  Hist( "SumEvents" )->Fill( 2,  fabs(b_weightPU)       ); // event pileup weight
   Hist( "SumEvents" )->Fill( 9,  1                       ); // event without MC weight
  
 
   // Cut 2: pass trigger
   if (m_isData) checkTriggerTurnOn();

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
            if (myjet.subjet_softdrop_N() >= 2) {
              double deta = fabs(myjet.subjet_softdrop_eta()[0] - myjet.subjet_softdrop_eta()[1]);
              double dphi = fabs(myjet.subjet_softdrop_phi()[0] - myjet.subjet_softdrop_phi()[1]);
              double dr = sqrt(deta*deta + dphi*dphi);
              b_ak8jet0_subjet01_dr = dr;
              b_ak8jet0_subjet01_deta = deta;
              b_ak8jet0_subjet01_dphi = dphi;
              b_ak8jet0_subjet0_pt = myjet.subjet_softdrop_pt()[0];
              b_ak8jet0_subjet1_pt = myjet.subjet_softdrop_pt()[1];
              b_ak8jet0_subjet0_csv = myjet.subjet_softdrop_csv()[0];
              b_ak8jet0_subjet1_csv = myjet.subjet_softdrop_csv()[1];
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
            if (myjet.subjet_softdrop_N() >= 2) {
              double deta = fabs(myjet.subjet_softdrop_eta()[0] - myjet.subjet_softdrop_eta()[1]);
              double dphi = fabs(myjet.subjet_softdrop_phi()[0] - myjet.subjet_softdrop_phi()[1]);
              double dr = sqrt(deta*deta + dphi*dphi);
              b_ak8jet1_subjet01_dr = dr;
              b_ak8jet1_subjet01_deta = deta;
              b_ak8jet1_subjet01_dphi = dphi;
              b_ak8jet1_subjet0_pt = myjet.subjet_softdrop_pt()[0];
              b_ak8jet1_subjet1_pt = myjet.subjet_softdrop_pt()[1];
              b_ak8jet1_subjet0_csv = myjet.subjet_softdrop_csv()[0];
              b_ak8jet1_subjet1_csv = myjet.subjet_softdrop_csv()[1];
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

  
  bool foundJet=false;
  UZH::Jet Jet;
  m_logger << VERBOSE << "kFatJetsDeltaEta" << SLogger::endmsg;
  

  
  if(goodFatJets.size() >0)	 {
    Jet = goodFatJets[0];
    
    if (m_isSignal) {
      std::vector<UZH::Jet> selectedJets;
      selectedJets.push_back(Jet);
      m_bTaggingScaleTool.fillEfficiencies(selectedJets);
      m_bTaggingScaleTool.fillSoftdropSubjetEfficiencies(selectedJets);
    }
    // std::cout << " goodFatJets pt " << goodFatJets[0].pt()<< " ;goodFatJets[0].pruned_massCorr() " <<  goodFatJets[0].pruned_massCorr() << " ; eta " << goodFatJets[0].eta() <<" ; phi "<<  goodFatJets[0].phi() << std::endl;
    foundJet=true;
    for (unsigned int s=0;s<m_catNames.size();++s) {
      if (m_catNames[s].find("NoWindow") != std::string::npos) { 
	selectionBits[s].SetBitNumber( kMet );
      }
    }

  }
  else return;
  //} ///close loop goodFatJet
  
  
    
    m_logger << VERBOSE << "kMassWindow" << SLogger::endmsg;

    // m_logger << VERBOSE << Jet << SLogger::endmsg;    
    
    // if (moveOn==false) return;
    //commented 

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
    for (int i = 0; i < Jet.subjet_softdrop_N(); ++i) {
      if (m_bTaggingScaleTool.isTagged(Jet.subjet_softdrop_csv()[i])) {
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
      b_weight= getEventWeight();
      std::vector<UZH::Jet> selectedJets;
      selectedJets.push_back(Jet);
    
      // b_weightBtag = m_bTaggingScaleTool.getScaleFactor(selectedJets); // event b-tag SF weight
      b_weightBtag = m_bTaggingScaleTool.getSoftdropSubjetScaleFactor(selectedJets); // event b-tag SF weight
      
      b_weight *= b_weightBtag;
    }
  
     
  
  
    ///} // 1 fat jets
  
  // Rest of book-keeping
  // Hist( "Events" )->Fill( 0., b_weightGen        ); // event with MC weight
  Hist( "Events" )->Fill( 1,  b_weight           ); // event total weight
  // Hist( "Events" )->Fill( 2,  b_weightPU         ); // event pileup weight
  Hist( "Events" )->Fill( 3,  b_weightBtag ); // event b-tag SF weight
  // // Hist( "Events" )->Fill( 3,  b_weight_elScale   ); // event electron SF weight
  // // Hist( "Events" )->Fill( 4,  b_weight_muScale   ); // event muon SF weight
  // Hist( "Events" )->Fill( 9,  1                  ); // event without MC weight
  // Hist( "SumEvents" )->Fill( 0., fabs(b_weightGen)       ); // event with MC weight
  Hist( "SumEvents" )->Fill( 1,  fabs(b_weight)          ); // event total weight
  // Hist( "SumEvents" )->Fill( 2,  fabs(b_weightPU)       ); // event pileup weight
  Hist( "SumEvents" )->Fill( 3,  fabs(b_weightBtag)); // event b-tag SF weight
   Hist( "SumEvents" )->Fill( 3,  fabs(b_weightBtag_veto)); // event b-tag SF weight
  // // Hist( "SumEvents" )->Fill( 3,  fabs(b_weight_elScale)  ); // event electron SF weight
  // // Hist( "SumEvents" )->Fill( 4,  fabs(b_weight_muScale)  ); // event muon SF weight
  // // Hist( "SumEvents" )->Fill( 6,  fabs(b_weight_jvfScale));  // event JVF SF weight
  // Hist( "SumEvents" )->Fill( 9,  1                       ); // event without MC weight
  
 
  
  bool foundTau=false;
 
  std::vector<UZH::Tau> goodTaus;
  for ( int i = 0; i <   (m_tau.N)
	  ; ++i ) {
    UZH::Tau mytau( &m_tau, i );
    if (mytau.pt() > m_tauPtCut && mytau.TauType() ==2 ){
    
      if (fabs(mytau.eta()) < m_tauEtaCut){

	// std::cout<<" mytau.pt() " << mytau.pt() <<" mytau.eta() " << mytau.eta() <<" mytau.TauType() " <<  mytau.TauType() <<std::endl;
	if( mytau.byVLooseIsolationMVArun2v1PWnewDMwLT()>0.5   // && mytau.againstElectronVLooseMVA6>0.5
	    ){
	  // std::cout<<" mytau.pt() " << mytau.pt() <<" mytau.eta() " << mytau.eta() <<" mytau.TauType() " <<  mytau.TauType() <<std::endl;
	  goodTaus.push_back(mytau);
	  if (goodTaus.size()<2) isBoostedTau=true;
	  foundTau=true;
	}
      }
    }
  }
  
  for ( int i = 0; i <   (m_tau.N); ++i ) {
    UZH::Tau mytau( &m_tau, i );
    if (mytau.pt() > m_tauPtCut && mytau.TauType()==1 ){
      
      if (fabs(mytau.eta()) < m_tauEtaCut){
	
	// std::cout<<" mytau.pt() " << mytau.pt() <<" mytau.eta() " << mytau.eta() <<" mytau.TauType() " <<  mytau.TauType() <<std::endl;
	if( mytau.byVLooseIsolationMVArun2v1PWnewDMwLT()>0.5   // && mytau.againstElectronVLooseMVA5()>0.5
	    ){
	  // std::cout<<" mytau.pt() " << mytau.pt() <<" mytau.eta() " << mytau.eta() <<" mytau.TauType() " <<  mytau.TauType() <<std::endl;
	  goodTaus.push_back(mytau);
	   if (goodTaus.size()<2) isBoostedTau=false;
	  foundTau=true;
	}
	
      }
	
    }
  }


 
 
 
  for (unsigned int s=0;s<m_catNames.size();++s) {
    if (m_catNames[s].find("mutau") != std::string::npos || m_catNames[s].find("eletau") != std::string::npos|| m_catNames[s].find("tautau") != std::string::npos ) { 
      if  (goodTaus.size()>0){
	selectionBits[s].SetBitNumber( kTau );
	selectionBits[s].SetBitNumber( kTauIsolation );
      }
    }
  }
  bool foundMet=false;
 // Cut : select met
  UZH::MissingEt goodMet( &m_missingEt, 0 );
  if (goodMet.et() > m_metCut) {
    // std::cout<<"MET " << goodMet.et() <<std::endl;
    foundMet=true;
  }
  for (unsigned int s=0;s<m_catNames.size();++s) {
    if (m_catNames[s].find("mutau") != std::string::npos || m_catNames[s].find("eletau") != std::string::npos|| m_catNames[s].find("tautau") != std::string::npos) { 
      selectionBits[s].SetBitNumber( kMet );
    }
  }
  
  // std::cout << "before ak4 jets"<< std::endl;
 
  bool foundAK4Jets=false;
  HT=0.;
  int goodAK4Index=0;
  taggedAK4_L=0;
  taggedAK4_M=0;
  taggedAK4_T=0;
  // std::cout << "before ak4 jets: number of jets "<< (m_jetAK4.N)  <<std::endl;

  std::vector<UZH::Jet> goodJetsAK4;

  for ( int i = 0; i < (m_jetAK4.N); ++i ) {
    UZH::Jet myjetak4( &m_jetAK4, i );
    if (fabs(myjetak4.eta()) < m_AK4jetEtaCut) {
      if (myjetak4.pt() > m_AK4jetPtCut) {


	HT = HT + myjetak4.tlv().Et();
	// m_logger << INFO << myjetak4.pt() << SLogger::endmsg;
	Book( TH1F( "Jet_pt", "Jet p_{T} [GeV]", 100, 0.0, 1000 ), "AK4_low" )->Fill( myjetak4.pt() );
	Book( TH1F( "Jet_eta", "Jet #eta", 60, -3,3 ), "AK4_low" )->Fill( myjetak4.eta() );
	Book( TH1F( "Jet_m", "Jet m [GeV]", 50, 0, 250 ), "AK4_low" )->Fill( myjetak4.m() );
	Book( TH1F( "Jet_phi", "Jet #phi", 50, -TMath::Pi(), TMath::Pi() ), "AK4_low" )->Fill( myjetak4.phi() );
	Book( TH1F( "Jet_e", "Jet e [GeV]", 100, 0.0, 1000 ), "AK4_low" )->Fill( myjetak4.e() );
    
	Book( TH1F( "Jet_muf", "Jet muf", 50, 0.0, 1.0 ), "AK4_low" )->Fill( myjetak4.muf() );
	Book( TH1F( "Jet_phf", "Jet phf", 50, 0.0, 1.0 ), "AK4_low" )->Fill( myjetak4.phf() );
	Book( TH1F( "Jet_emf", "Jet emf", 50, 0.0, 1.0 ), "AK4_low" )->Fill( myjetak4.emf() );
	Book( TH1F( "Jet_nhf", "Jet nhf", 50, 0.0, 1.0 ), "AK4_low" )->Fill( myjetak4.nhf() );
	Book( TH1F( "Jet_chf", "Jet chf", 50, 0.0, 1.0 ), "AK4_low" )->Fill( myjetak4.chf() );
	Book( TH1F( "Jet_area", "Jet area", 40, 0.0, 4.0 ), "AK4_low" )->Fill( myjetak4.area() );
	Book( TH1F( "Jet_cm", "Jet cm", 50, 0.0, 50.0 ), "AK4_low" )->Fill( myjetak4.cm() );
	Book( TH1F( "Jet_nm", "Jet nm", 50, 0.0, 50.0 ), "AK4_low" )->Fill( myjetak4.nm() );
	if (myjetak4.IDTight()){
	  goodJetsAK4.push_back(myjetak4);
	  if ( myjetak4.csv() >  0.605 ) taggedAK4_L++;
	  if ( myjetak4.csv() >  0.89 ) taggedAK4_M++;
	  if ( myjetak4.csv() >  0.97 ) taggedAK4_T++;
	}
      }
    }
  }
  if (goodJetsAK4.size()>0){
     if (m_isSignal) {
         m_bTaggingScaleTool.fillEfficiencies_veto(goodJetsAK4);
     }    
     
    if (!m_isData) {
      b_weightBtag_veto = m_bTaggingScaleTool.getScaleFactor_veto(goodJetsAK4); // event b-tag SF weight
      b_weight *= b_weightBtag_veto;
    }
  }

  
  // std::cout << "after ak4 jets, HT " <<HT<< std::endl;
  
  //Jet for trigger turn on on the first cat, anyways just checking triggers and met filters doesn't change
  if   (selectionBits[0].TestBitNumber(  kTrigger )  && selectionBits[0].TestBitNumber( kMetFilters )) {
  
    //if (!foundJet) return; it already returns before if doesn't find a jet AK8
    if(Jet.pt()<20 ) std::cout << "Jet pt " <<Jet.pt()<< " ;Jet.pruned_massCorr() " << Jet.pruned_massCorr() << " ;Jet eta " <<Jet.eta() <<" ;Jet phi "<< Jet.phi() << std::endl;
    Hist( "JetPt" )->Fill(Jet.pt());
    Hist( "Jet_pruned_massCorr")->Fill(Jet.pruned_massCorr());
    Hist( "SumEt")->Fill(goodMet.sumEt());
    Hist( "MET")->Fill(goodMet.et());
    Hist( "HT")->Fill(HT);
    Hist( "N_ak4")->Fill(goodJetsAK4.size());
    Hist( "N_ak8")->Fill(goodFatJets.size());

    if (!m_isData) {
      for( unsigned int v = 0; v < (m_eventInfo.actualIntPerXing)->size(); ++v ){
      
	if ( (*m_eventInfo.bunchCrossing)[v] == 0 ) {
	Hist( "N_vertex")->Fill((*m_eventInfo.actualIntPerXing)[v] );
	}
    
      }
    }
  }
  // std::cout << "after hists for trigger efficiency "<< std::endl;

  if (!foundMet) return;
  if (!foundTau) return;
  bool foundSecondLepton=false;
  bool foundElectron=false;
  std::vector<UZH::Electron> goodElectrons;
  for ( int i = 0; i <   (m_electron.N)
	  ; ++i ) {
    UZH::Electron myelectron( &m_electron, i );
    if (myelectron.pt() > m_electronPtCut) {
      if (fabs(myelectron.eta()) < m_electronEtaCut) {

	if (myelectron.isLooseElectronBoosted()==1 || myelectron.isMediumElectronBoosted()==1   ){
	  //	if (myelectron.isLooseElectron()==1 || myelectron.isMediumElectron()==1   ){
	  // std::cout<<" myelectron.pt() " << myelectron.pt() <<" myelectron.eta() " << myelectron.eta() <<" myelectron.isLooseElectronBoosted() "<<  myelectron.isLooseElectronBoosted() <<" myelectron.isMediumElectronBoosted()"<<myelectron.isMediumElectronBoosted() <<" iso " << myelectron.SemileptonicPFIso()/myelectron.pt() << std::endl;
	  if(( fabs(myelectron.eta()) <=  1.479 &&  (myelectron.SemileptonicCorrPFIso()/myelectron.pt()) < 0.0983) ||
	     ( fabs(myelectron.eta()) >   1.479  &&  (myelectron.SemileptonicCorrPFIso()/myelectron.pt()) < 0.121 ) ){
	    goodElectrons.push_back(myelectron);
		 
	    foundElectron=true;
	  }
	}
      }
    }
  }
  if (foundElectron){ 
    for (unsigned int s=0;s<m_catNames.size();++s) {

      TLorentzVector ele_tlv = goodElectrons[0].tlv();
      TLorentzVector tau_tlv = goodTaus[0].tlv();
      if ( ele_tlv.DeltaR(tau_tlv) >0.1 // && ele_tlv.DeltaR(tau_tlv) < 1.5
	   &&
	   m_catNames[s].find("eletau") != std::string::npos) {
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
	if (mymuon.isLooseMuon()==1){//  mymuon.isLooseMuon()
	  // std::cout<<" mymuon.SemileptonicPFIso() "<<  mymuon.SemileptonicPFIso() <<std::endl;
	  // std::cout<<" mymuon.pt() " << mymuon.pt() <<" mymuon.eta() " << mymuon.eta() <<" mymuon.isLooseMuon() "<<  mymuon.isLooseMuon() << " iso "  <<mymuon.SemileptonicPFIso() /mymuon.pt() <<std::endl;
	  


	  if (mymuon.SemileptonicCorrPFIso() /mymuon.pt() <0.25){
	    goodMuons.push_back(mymuon);
	    
	    foundMuon=true;
	  }
	}
      }
    }
  }
  
  if (foundMuon){ 
    TLorentzVector muon_tlv = goodMuons[0].tlv();
    TLorentzVector tau_tlv = goodTaus[0].tlv();

    if( muon_tlv.DeltaR(tau_tlv) >0.1 ){// && muon_tlv.DeltaR(tau_tlv) <1.5
      // std::cout << "Checking DeltaR goodMuons[0].tlv().pt() " << goodMuons[0].tlv().Pt()<< " &&  goodTaus[0].tlv().Pt()" << goodTaus[0].tlv().Pt() <<std::endl;
      // std::cout <<"Checking DeltaR DR "<<  muon_tlv.DeltaR(tau_tlv) <<std::endl;
      for (unsigned int s=0;s<m_catNames.size();++s) {	
	if (m_catNames[s].find("mutau") != std::string::npos) { 
	  selectionBits[s].SetBitNumber( kLepton );
	}
      }
    }
  }
  
  for (unsigned int s=0;s<m_catNames.size();++s) {
    if (goodTaus.size()>1   && goodTaus[0].tlv().DeltaR( goodTaus[1].tlv()) >0.1
	&& m_catNames[s].find("tautau") != std::string::npos) {
      selectionBits[s].SetBitNumber( kLepton );
    }
  }
 


  
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
    double ak8jet0_subjet01_deta=-99;
    double ak8jet0_subjet01_dphi=-99;
    double ak8jet0_subjet01_dr = -99;
    
    if (Jet.subjet_softdrop_N() >= 2) {
      ak8jet0_subjet01_deta = fabs(Jet.subjet_softdrop_eta()[0] - Jet.subjet_softdrop_eta()[1]);
      ak8jet0_subjet01_dphi = fabs(Jet.subjet_softdrop_phi()[0] - Jet.subjet_softdrop_phi()[1]);
      ak8jet0_subjet01_dr = sqrt(ak8jet0_subjet01_deta*ak8jet0_subjet01_deta + ak8jet0_subjet01_dphi*ak8jet0_subjet01_dphi);
    }
    for (int i = 0; i < Jet.subjet_softdrop_N(); ++i) {
      switch(i) {
      case 0:
	vJet_subjet0_csv = Jet.subjet_softdrop_csv()[i];
	break;
      case 1:
	vJet_subjet1_csv = Jet.subjet_softdrop_csv()[i];
	break;
      }
      
      if ((m_bTaggingScaleTool.isTagged(Jet.subjet_softdrop_csv()[i]))) {
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
	  bool triggerStudy_=false;
	  if (triggerStudy_){
	    bool channel_alreadyassigned=false;

	    if (m_catNames[s].find("mutau")  != std::string::npos && selectionBits[s].TestBitNumber(kTau )&& selectionBits[s].TestBitNumber(kLepton ) && !channel_alreadyassigned ) 	{
	      // std::cout<<" checking out offline selection for mutau; mu iso = "<< goodMuons[0].SemileptonicPFIso() << std::endl;
	      checkTrigger(m_catNames[s],"mutau");
	    
	      checkTriggerOfflineSelection(m_catNames[s], Jet,goodTaus[0], goodMuons[0].tlv() , goodMuons[0].isTightMuon(), goodMuons[0].SemileptonicPFIso(),
					   goodMet, "mutau", HT);
	      channel_alreadyassigned=true;
	    }
	    if (m_catNames[s].find("eletau")  != std::string::npos && selectionBits[s].TestBitNumber(kTau )&& selectionBits[s].TestBitNumber(kLepton )  && !channel_alreadyassigned ) {
	      // std::cout<<" checking out offline selection for eletau; ele iso = " << goodElectrons[0].SemileptonicPFIso() <<std::endl;	
	      checkTrigger(m_catNames[s],"eletau");  
	      checkTriggerOfflineSelection(m_catNames[s], Jet,goodTaus[0], goodElectrons[0].tlv() , goodElectrons[0].isTightElectron(), goodElectrons[0].SemileptonicPFIso(),
					   goodMet, "eletau", HT);
	      channel_alreadyassigned=true; 
	    }
	    if (m_catNames[s].find("tautau")  != std::string::npos && selectionBits[s].TestBitNumber(kTau )&& selectionBits[s].TestBitNumber(kLepton ) && !channel_alreadyassigned ) { 
	      //std::cout<<" checking out offline selection for tautau"<< std::endl;
	      checkTrigger(m_catNames[s],"tautau");
	      checkTriggerOfflineSelection(m_catNames[s], Jet,goodTaus[0], goodTaus[1].tlv() , 0,0, goodMet, "tautau", HT);
	      channel_alreadyassigned=true;
	    }
	  }
	}
	fillHistograms(m_catNames[s], Jet, Jet, diJet, vJet_tau21, vJet_tau31, vJet_tau32, vJet_nTaggedSubjets, vJet_subjet0_csv, vJet_subjet1_csv, vJet_tau21, vJet_tau31, vJet_tau32, vJet_nTaggedSubjets, vJet_subjet0_csv, vJet_subjet1_csv, deta, dphi, dr );
      }
    }
    
    m_logger << VERBOSE << "after hist" << SLogger::endmsg;
    
 
    
    
    for (unsigned int s=0;s<m_catNames.size();++s) {
      
      if (m_catNames[s].find("NoWindow") != std::string::npos) {

	if (passed_all[s]) {
	   bool channel_alreadyassigned=false; 
	  if (m_catNames[s].find("mutau")  != std::string::npos && selectionBits[s].TestBitNumber(kTau )&& selectionBits[s].TestBitNumber(kLepton )   && !channel_alreadyassigned  ) {
	    // std::cout << "Bfore filling branch goodMuons[0].tlv().pt() " << goodMuons[0].tlv().Pt() <<  "&&  goodTaus[0].pt() " << goodTaus[0].pt()<<std::endl;
	    FillBranches( "mutau", Jet, vJet_tau21, vJet_tau31, vJet_tau32, vJet_nTaggedSubjets, vJet_subjet0_csv,vJet_subjet1_csv,  ak8jet0_subjet01_dr,   ak8jet0_subjet01_deta, ak8jet0_subjet01_dphi,  b_weight,  b_weightGen,  b_weightPU,  b_weightBtag, b_weightBtag_veto,   b_runNumber, b_eventNumber, b_lumiBlock,  goodTaus[0], goodMuons[0].tlv() , goodMet);
	      channel_alreadyassigned=true; 
	  }
	  
	  if (m_catNames[s].find("eletau")  != std::string::npos && selectionBits[s].TestBitNumber(kTau )&& selectionBits[s].TestBitNumber(kLepton ) && !channel_alreadyassigned ) {
	    FillBranches( "eletau", Jet, vJet_tau21, vJet_tau31, vJet_tau32, vJet_nTaggedSubjets, vJet_subjet0_csv,vJet_subjet1_csv,  ak8jet0_subjet01_dr,   ak8jet0_subjet01_deta, ak8jet0_subjet01_dphi,  b_weight,  b_weightGen,  b_weightPU,  b_weightBtag, b_weightBtag_veto,    b_runNumber, b_eventNumber, b_lumiBlock,  goodTaus[0], goodElectrons[0].tlv(), goodMet );
	    channel_alreadyassigned=true; 
	  }	
	  // if (m_catNames[s].find("tautau")  != std::string::npos && selectionBits[s].TestBitNumber(kTau )&& selectionBits[s].TestBitNumber(kLepton ) && !channel_alreadyassigned ) {
	  //   FillBranches( "tautau", Jet, vJet_tau21, vJet_tau31, vJet_tau32, vJet_nTaggedSubjets, vJet_subjet0_csv,vJet_subjet1_csv,  ak8jet0_subjet01_dr,   ak8jet0_subjet01_deta, ak8jet0_subjet01_dphi,  b_weight,  b_weightGen,  b_weightPU,  b_weightBtag,  b_weightBtag_veto,   b_runNumber, b_eventNumber, b_lumiBlock,  goodTaus[0], goodTaus[1].tlv() , goodMet  );
	  //   channel_alreadyassigned=true; 
	  // }
	}
      }
    }
    
    m_logger << VERBOSE << "after ntuple" << SLogger::endmsg;
  }    
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
  b_weightBtag_veto = 1.;
  
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
  for(unsigned int ch=0;ch< channels_.size();ch++){
    b_weight_[channels_[ch]] = 1.;
    b_weightGen_[channels_[ch]] = 1.;
    b_weightPU_[channels_[ch]] = 1.;
    b_weightBtag_[channels_[ch]] = 1.;
    b_weightBtag_veto_[channels_[ch]] = 1.;
  
    b_runNumber_[channels_[ch]] = -99;;
    b_eventNumber_[channels_[ch]] = -99;
    b_lumiBlock_[channels_[ch]] = -99;
  
    b_ak8jet0_pt_[channels_[ch]] = -99;
    b_ak8jet0_phi_[channels_[ch]] = -99;
    b_ak8jet0_eta_[channels_[ch]] = -99;
    b_ak8jet0_e_[channels_[ch]] = -99;
    b_ak8jet0_tau21_[channels_[ch]] = -99;
    b_ak8jet0_m_[channels_[ch]] = -99;
    b_ak8jet0_mpruned_[channels_[ch]] = -99;
    b_ak8jet0_csv_[channels_[ch]] = -99;
    b_ak8jet0_subjet01_dr_[channels_[ch]] = -99;
    b_ak8jet0_subjet01_deta_[channels_[ch]] = -99;
    b_ak8jet0_subjet01_dphi_[channels_[ch]] = -99;
    b_ak8jet0_subjet0_pt_[channels_[ch]] = -99;
    b_ak8jet0_subjet1_pt_[channels_[ch]] = -99;
    b_ak8jet0_subjet0_csv_[channels_[ch]] = -99;
    b_ak8jet0_subjet1_csv_[channels_[ch]] = -99;

    b_tau_pt_[channels_[ch]] =-99;
    b_tau_eta_[channels_[ch]] =-99;
    b_tau_phi_[channels_[ch]] =-99;
    b_seclepton_pt_[channels_[ch]] =-99;
    b_seclepton_eta_[channels_[ch]] =-99;
    b_seclepton_phi_[channels_[ch]] =-99;
    b_met_[channels_[ch]] =-99;
    b_met_phi_[channels_[ch]] =-99;
    b_vis_mass_ll_[channels_[ch]] =-99;
    b_dR_ll_[channels_[ch]] =-99;
    b_H_Mass_SVFit_[channels_[ch]] =-99;
    b_H_Pt_SVFit_[channels_[ch]] =-99;
    b_H_Eta_SVFit_[channels_[ch]] =-99;
    b_H_Phi_SVFit_[channels_[ch]] =-99;
    b_X_Mass_SVFit_[channels_[ch]] =-99;
    b_X_Pt_SVFit_[channels_[ch]] =-99;
    b_X_Eta_SVFit_[channels_[ch]] =-99;
    b_X_Phi_SVFit_[channels_[ch]] =-99;

    b_trig_HT800_[channels_[ch]] =false;
    b_trig_ak8Jet360_Trim30_[channels_[ch]] =false;
    b_trig_Ele27_[channels_[ch]] =false;
    b_trig_Ele23_[channels_[ch]] =false;
    b_trig_Mu20_[channels_[ch]] =false;
    b_trig_Mu18_[channels_[ch]] =false;

    b_selection_bits_.clear();
    b_selection_lastcut_.clear();
  }
  
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
  Book( TH1F( "Trigger_Tau", "Trigger_Tau", 60,0,60), directory);
  Book( TH1F( "Trigger_Jet", "Trigger_Jet", 60,0,60), directory);
  Book( TH1F( "Trigger_TauOrJet", "Trigger_TauOrJet",10,-1.5,8.5  ) , directory);

  Book( TH1F( "Trigger_Offline", "Trigger offline selection",60,0,60  ), directory);
  Book( TH1F( "Trigger_Tau_Offline", "Trigger_Tau offline selection", 60,0,60), directory);
  Book( TH1F( "Trigger_Jet_Offline", "Trigger_Jet offline selection", 60,0,60), directory);
  Book( TH1F( "Trigger_TauOrJet_Offline", "Trigger_TauOrJet  offline selection",10,-1.5,8.5  ) , directory);
  Book( TH1F( "Trigger_Attempts", "Trigger attempts fro 2016",30,0,30 ) , directory);

 

  std::vector<std::string> trignames;

  trignames.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v");
  trignames.push_back("HLT_AK8PFJet360_TrimMass30_v");
  trignames.push_back("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v");
  trignames.push_back("HLT_AK8PFHT660_TrimR0p1PT0p03Mass50_BTagCSV_p20_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p41_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v");
  trignames.push_back("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_v");
  trignames.push_back("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v");
  trignames.push_back("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v");
  trignames.push_back("HLT_PFHT600_v");
  trignames.push_back("HLT_PFHT650_v") ;
  trignames.push_back("HLT_PFHT800_v") ;
  trignames.push_back("HLT_PFHT900_v") ;
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

  trignames.push_back("HLT_Ele15_IsoVVVL_PFHT600");
  trignames.push_back("HLT_Mu15_IsoVVVL_PFHT600");
  trignames.push_back("HLT_Ele15_IsoVVVL_PFHT350_PFMET50");
  trignames.push_back("HLT_Mu15_IsoVVVL_PFHT350_PFMET50");
  trignames.push_back("HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400"); 
  trignames.push_back("HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400");
 
 
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
  trignamesJet.push_back("HLT_PFHT900_v") ;



  const char *TriggerNames[trignames.size()];
  const char *TriggerTauNames[trignamesTau.size()];
  const char *TriggerJetNames[trignamesJet.size()];




  Hist( "Trigger"  	, directory )->GetXaxis()->SetTitle( "Triggers" );
  Hist( "Trigger"  	, directory )->GetXaxis()->LabelsOption("v");
  Hist( "Trigger_Tau"  	, directory )->GetXaxis()->SetTitle( "Triggers" );
  Hist( "Trigger_Tau"  	, directory )->GetXaxis()->LabelsOption("v");
  Hist( "Trigger_Jet"  	, directory )->GetXaxis()->SetTitle( "Triggers" );
  Hist( "Trigger_Jet"   , directory )->GetXaxis()->LabelsOption("v");
 
  Hist( "Trigger_Offline"  	, directory )->GetXaxis()->SetTitle( "Triggers" );
  // Hist( "Trigger_Offline"  	, directory )->GetXaxis()->LabelsOption("v");
  Hist( "Trigger_Tau_Offline"  	, directory )->GetXaxis()->SetTitle( "Triggers" );
  //Hist( "Trigger_TaU_Offline"  	, directory )->GetXaxis()->LabelsOption("v");
  Hist( "Trigger_Jet_Offline"  	, directory )->GetXaxis()->SetTitle( "Triggers" );
  //Hist( "Trigger_Jet_Offline"   , directory )->GetXaxis()->LabelsOption("v");
 


  for ( unsigned int t=0; t< trignames.size(); t++){
    TriggerNames[t]= trignames[t].c_str();
    Hist( "Trigger"  , directory 	)->GetXaxis()->SetBinLabel(t+1, TriggerNames[t] );
    Hist( "Trigger_Offline"  , directory 	)->GetXaxis()->SetBinLabel(t+1, TriggerNames[t] );
  }

  for ( unsigned int t=0; t< trignamesTau.size(); t++){ 
    TriggerTauNames[t]= trignamesTau[t].c_str();
    Hist( "Trigger_Tau"  , directory 	)->GetXaxis()->SetBinLabel(t+1, TriggerTauNames[t] );
    Hist( "Trigger_Tau_Offline"  , directory 	)->GetXaxis()->SetBinLabel(t+1, TriggerTauNames[t] );
    }

  for ( unsigned int t=0; t< trignamesJet.size(); t++){ 
    TriggerJetNames[t]= trignamesJet[t].c_str();
    Hist( "Trigger_Jet"  , directory 	)->GetXaxis()->SetBinLabel(t+1, TriggerJetNames[t] );
    Hist( "Trigger_Jet_Offline"  , directory 	)->GetXaxis()->SetBinLabel(t+1, TriggerJetNames[t] );
  }

  std::vector<std::string> trignames_attempts;

  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet200_TrimMass20_v1");
  trignames_attempts.push_back("HLT_HT400_MET50_AK8PFJet200_TrimMass20_v1");
  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet200_TrimMass30_v1");
  trignames_attempts.push_back("HLT_HT350_MET50_AK8PFJet200_TrimMass30_v1"); 
  // trignames_attempts.push_back("HLT_HT200_MET50_AK8DiPFJet120_TrimMass20_v1");
  // trignames_attempts.push_back("HLT_HT250_MET50_AK8DiPFJet120_TrimMass20_v1");
  // trignames_attempts.push_back("HLT_HT300_MET50_AK8DiPFJet120_TrimMass20_v1");
  // trignames_attempts.push_back("HLT_HT300_MET50_AK8DiPFJet200_TrimMass20_v1");
  trignames_attempts.push_back("HLT_HT200_MET50_AK8PFJet120_TrimMass20_LooseIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT400_LooseIsoPFTau50_Trk30_v1"); 
  trignames_attempts.push_back("HLT_HT450_LooseIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT450_LooseIsoPFTau30_Trk30_v1");
  trignames_attempts.push_back("HLT_HT200_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT400_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet200_TrimMass20_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet200_TrimMass30_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT350_MET50_AK8PFJet200_TrimMass30_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT350_MET60_AK8PFJet200_TrimMass20_NotIsoPFTau50_Trk30_v1"); 
  trignames_attempts.push_back("HLT_HT350_MET50_AK8PFJet200_TrimMass20_NotIsoPFTau50_v1");
  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet120_TrimMass30_NotIsoPFTau50_v1");
  trignames_attempts.push_back("HLT_HT400_MET50_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT400_MET70_NotIsoPFTau50_Trk30_v1"); 
  trignames_attempts.push_back("HLT_HT450_MET50_NotIsoPFTau50_Trk30_v1"); 

trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet120_TrimMass30_NotIsoPFTau50_Trk30_v1"); 
trignames_attempts.push_back("HLT_HT300_MET60_AK8PFJet120_TrimMass30_NotIsoPFTau50_Trk30_v1");  
trignames_attempts.push_back("HLT_HT400_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30_v1");  

trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet150_TrimMass30_NotIsoPFTau50_Trk30_v1");   
trignames_attempts.push_back("HLT_HT300_MET60_AK8PFJet150_TrimMass30_NotIsoPFTau50_Trk30v1"); 
trignames_attempts.push_back("HLT_HT350_MET50_AK8PFJet150_TrimMass30_NotIsoPFTau50_Trk30v1"); 
trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet150_TrimMass20_NotIsoPFTau50_Trk30v1"); 


  const char *TriggerAttemptsNames[trignames_attempts.size()];
  for ( unsigned int t=0; t< trignames_attempts.size(); t++){ 
    TriggerAttemptsNames[t]= trignames_attempts[t].c_str();
    Hist( "Trigger_Attempts"  , directory 	)->GetXaxis()->SetBinLabel(t+1,   TriggerAttemptsNames[t] );
 
  }

}


void VHTausAnalysis::bookHistograms( const TString& directory ) {
  
  // kinematics histograms
  Book( TH1F( "vjet_pt", "Vjet p_{T};Vjet p_{T} [GeV]", 200, 0, 2000 ), directory ); 
  Book( TH1F( "vjet_eta", "Vjet #eta;Vjet #eta", 50, -2.5, 2.5 ), directory ); 
  Book( TH1F( "vjet_phi", "Vjet #phi;Vjet #phi", 50, -3.15, 3.15 ), directory ); 
  Book( TH1F( "vjet_m", "Vjet m;Vjet m [GeV]", 40, 0, 200 ), directory ); 
  Book( TH1F( "vjet_mpruned", "Vjet softdrop m;Vjet softdrop m [GeV]", 40, 0, 200 ), directory ); 
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
  Book( TH1F( "hjet_mpruned", "Hjet softdrop m;Hjet softdrop m [GeV]", 40, 0, 200 ), directory ); 
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
  Hist( "vjet_nSubjets", directory )->Fill( vectorJet.subjet_softdrop_N() , b_weight);
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
  Hist( "hjet_nSubjets", directory )->Fill( higgsJet.subjet_softdrop_N() , b_weight);
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


void VHTausAnalysis::checkTriggerOfflineSelection(const TString& directory,  const UZH::Jet& Jet,const UZH::Tau& tau,const  TLorentzVector& lepton, const  bool& ID,const  double& leptonIso, const UZH::MissingEt& met,const std::string& channel,const  double& HT  ){
  // std:: cout << "directory " <<  directory << " channel " << channel <<" jet pt " <<  Jet.pt() <<" tau tlv  pt " << tau.tlv().Pt()<< " lepton pt " << lepton.Pt() << std::endl;
  std::vector<std::string> trignames;

  trignames.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v");
  trignames.push_back("HLT_AK8PFJet360_TrimMass30_v");
  trignames.push_back("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v");
  trignames.push_back("HLT_AK8PFHT660_TrimR0p1PT0p03Mass50_BTagCSV_p20_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p41_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v");
  trignames.push_back("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_v");

  trignames.push_back("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v");

  trignames.push_back("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v");
  trignames.push_back("HLT_PFHT600_v");
  trignames.push_back("HLT_PFHT650_v") ;
  trignames.push_back("HLT_PFHT800_v") ;
  trignames.push_back("HLT_PFHT900_v") ;
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
  trignames.push_back("HLT_Ele15_IsoVVVL_PFHT600");
  trignames.push_back("HLT_Mu15_IsoVVVL_PFHT600");
  trignames.push_back("HLT_Ele15_IsoVVVL_PFHT350_PFMET50");
  trignames.push_back("HLT_Mu15_IsoVVVL_PFHT350_PFMET50");
  trignames.push_back("HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400"); 
  trignames.push_back("HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400");
 
 
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
  trignamesJet.push_back("HLT_PFHT900_v") ;

  
  std::string GENChannel_="General";
  bool isfired = false;
 
  bool isfired_Offline = false;
 
  if (GENChannel_=="General"){
    bool isfiredTau = false;
    bool isfiredJet = false;
    bool isfiredTau_Offline = false;
    bool isfiredJet_Offline = false;
 
    for( std::map<std::string,bool>::iterator it = (m_eventInfo.trigDecision)->begin(); it != (m_eventInfo.trigDecision)->end(); ++it){
      //   std::cout <<"TriggerFired : "<< (it->first) <<std::endl; 

      for( unsigned int t = 0; t < trignames.size(); ++t ){
      
	if( (it->first).find(trignames[t]) != std::string::npos){
	  // std::cout <<"TriggerFired : "<< (it->first) <<" and number "<<t<< std::endl; 
	  if(it->second  ){
	    isfired = true;
	    if (isfired &&( trignames[t].find("AK8") != std::string::npos|| trignames[t].find("HT") != std::string::npos )) {
	      if(Jet.pt()>500)	    Hist( "Trigger_Offline",directory )->Fill( t );
	    }
	    else   if (isfired && trignames[t].find("Ele") != std::string::npos ){
	      if( trignames[t].find("Ele105") == std::string::npos  && lepton.Pt()> 30  ){///not trigger ele 105, unisolated
		Hist( "Trigger_Offline",directory )->Fill( t );	
	      }
	      else if (isfired && lepton.Pt()>110){
		Hist( "Trigger_Offline",directory )->Fill( t );	
	      }
	    }
	    else if(isfired && trignames[t].find("Mu") != std::string::npos){
	      if( trignames[t].find("Mu45") == std::string::npos  && lepton.Pt()> 22){///not trigger ele 105, unisolated
		Hist( "Trigger_Offline",directory )->Fill( t );	
	      }
	      else if (isfired && lepton.Pt()>50){
		Hist( "Trigger_Offline",directory )->Fill( t );	
	      }
	    }else if(isfired && trignames[t].find("DoubleMediumIsoPFTau") != std::string::npos){
	      if (tau.tlv().Pt()>50 && lepton.Pt() >50)  Hist( "Trigger_Offline",directory )->Fill( t );	
	    }  
	    // else std::cout <<"Trigger NOT Fired : "<< trignames[t] <<std::endl;  	
	  }
	}
      }      
      
      
      for( unsigned int t = 0; t < trignamesJet.size(); ++t ){
	if( (it->first).find(trignamesJet[t]) != std::string::npos &&  it->second ){
	  isfiredJet = true;
	  if (isfiredJet &&( trignamesJet[t].find("AK8") != std::string::npos|| trignamesJet[t].find("HT") != std::string::npos )) {
	    if(Jet.pt()>500)	 {
	      Hist( "Trigger_Jet_Offline",directory )->Fill( t );
	      isfiredJet_Offline = true;
	    }
	  }
	}
      }
      for( unsigned int t = 0; t < trignamesTau.size(); ++t ){
	if( (it->first).find(trignamesTau[t]) != std::string::npos && it->second  // && !isfiredJet
	    ){
	  isfiredTau = true;
	  if (isfiredTau && trignamesTau[t].find("Ele") != std::string::npos  && channel=="eletau"){
	    if( trignamesTau[t].find("Ele105") == std::string::npos   && channel=="eletau" && ID && lepton.Pt()> 30){///not trigger ele 105, unisolated
	      Hist( "Trigger_Tau_Offline",directory )->Fill( t );	
	      isfiredTau_Offline = true;
	    }
	    else if (isfiredTau &&lepton.Pt()>110 && ID && channel=="eletau"){
	      Hist( "Trigger_Tau_Offline",directory )->Fill( t );	
	      isfiredTau_Offline = true;
	    }
	  }
	  else if(isfiredTau && trignamesTau[t].find("Mu") != std::string::npos && channel=="mutau"){
	    if(isfiredTau && trignamesTau[t].find("Mu45") == std::string::npos && channel=="mutau" && ID  && lepton.Pt()> 22 &&  leptonIso/lepton.Pt() <0.15 ){///not trigger mu45, unisolated
	      Hist( "Trigger_Tau_Offline",directory )->Fill( t );
	      isfiredTau_Offline = true;
	    }
	    else if (isfiredTau &&lepton.Pt()>50 && ID && channel=="mutau"){
	      Hist( "Trigger_Tau_Offline",directory )->Fill( t );	
	      isfiredTau_Offline = true;
	    }
	  }else if(isfiredTau && trignamesTau[t].find("DoubleMediumIsoPFTau") != std::string::npos && channel=="tautau"){
	    if (isfiredTau &&tau.tlv().Pt()>50 && lepton.Pt() >50 && channel=="tautau" ){
	      Hist( "Trigger_Tau_Offline",directory )->Fill( t );	
	      isfiredTau_Offline = true;
	    }  
	  }
	}
      }  
    }
    Hist( "Trigger_TauOrJet_Offline" ,directory)->Fill( -1 ); 	 
    // if ( (isfiredJet && isfiredJet_Offline) || ( isfiredTau  && isfiredTau_Offline )  )  Hist( "Trigger_TauOrJet_Offline" ,directory)->Fill( 0 ); 	 
    // if (isfiredJet && isfiredJet_Offline )    Hist( "Trigger_TauOrJet_Offline" ,directory)->Fill( 1 ); 
    // if (isfiredTau && isfiredTau_Offline )    Hist( "Trigger_TauOrJet_Offline" ,directory)->Fill( 2 ); 
    // if (isfiredJet && isfiredJet_Offline &&  !(isfiredTau_Offline))    Hist( "Trigger_TauOrJet_Offline",directory )->Fill( 3 ); 
    // if (isfiredTau && !( isfiredJet_Offline ))    Hist( "Trigger_TauOrJet_Offline",directory )->Fill( 4 ); 
    // if (!( isfiredTau_Offline ) && !( isfiredJet_Offline ))    Hist( "Trigger_TauOrJet_Offline",directory )->Fill( 5 ); 
   if ( isfiredJet_Offline ||isfiredTau_Offline   )  Hist( "Trigger_TauOrJet_Offline" ,directory)->Fill( 0 ); 	 
    if ( isfiredJet_Offline )    Hist( "Trigger_TauOrJet_Offline" ,directory)->Fill( 1 ); 
    if (isfiredTau_Offline )    Hist( "Trigger_TauOrJet_Offline" ,directory)->Fill( 2 ); 
    if ( isfiredJet_Offline &&  !(isfiredTau_Offline))    Hist( "Trigger_TauOrJet_Offline",directory )->Fill( 3 ); 
    if (isfiredTau_Offline && !( isfiredJet_Offline ))    Hist( "Trigger_TauOrJet_Offline",directory )->Fill( 4 ); 
    if (!( isfiredTau_Offline ) && !( isfiredJet_Offline ))    Hist( "Trigger_TauOrJet_Offline",directory )->Fill( 5 ); 
    
  }
  
  
  
  trignames.clear();
  trignamesJet.clear();
  trignamesTau.clear();

  std::vector<std::string> trignames_attempts;

  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet200_TrimMass20_v1");
  trignames_attempts.push_back("HLT_HT400_MET50_AK8PFJet200_TrimMass20_v1");
  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet200_TrimMass30_v1");
  trignames_attempts.push_back("HLT_HT350_MET50_AK8PFJet200_TrimMass30_v1"); 
  // trignames_attempts.push_back("HLT_HT200_MET50_AK8DiPFJet120_TrimMass20_v1");
  // trignames_attempts.push_back("HLT_HT250_MET50_AK8DiPFJet120_TrimMass20_v1");
  // trignames_attempts.push_back("HLT_HT300_MET50_AK8DiPFJet120_TrimMass20_v1");
  // trignames_attempts.push_back("HLT_HT300_MET50_AK8DiPFJet200_TrimMass20_v1");
  trignames_attempts.push_back("HLT_HT200_MET50_AK8PFJet120_TrimMass20_LooseIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT400_LooseIsoPFTau50_Trk30_v1"); 
  trignames_attempts.push_back("HLT_HT450_LooseIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT450_LooseIsoPFTau30_Trk30_v1");
  trignames_attempts.push_back("HLT_HT200_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT400_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet200_TrimMass20_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet200_TrimMass30_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT350_MET50_AK8PFJet200_TrimMass30_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT350_MET60_AK8PFJet200_TrimMass20_NotIsoPFTau50_Trk30_v1"); 
  trignames_attempts.push_back("HLT_HT350_MET50_AK8PFJet200_TrimMass20_NotIsoPFTau50_v1");
  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet120_TrimMass30_NotIsoPFTau50_v1");
  trignames_attempts.push_back("HLT_HT400_MET50_NotIsoPFTau50_Trk30_v1");
  trignames_attempts.push_back("HLT_HT400_MET70_NotIsoPFTau50_Trk30_v1"); 
  trignames_attempts.push_back("HLT_HT450_MET50_NotIsoPFTau50_Trk30_v1"); 

  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet120_TrimMass30_NotIsoPFTau50_Trk30_v1"); 
  trignames_attempts.push_back("HLT_HT300_MET60_AK8PFJet120_TrimMass30_NotIsoPFTau50_Trk30_v1");  
  trignames_attempts.push_back("HLT_HT400_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30_v1");  

  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet150_TrimMass30_NotIsoPFTau50_Trk30_v1");   
  trignames_attempts.push_back("HLT_HT300_MET60_AK8PFJet150_TrimMass30_NotIsoPFTau50_Trk30v1"); 
  trignames_attempts.push_back("HLT_HT350_MET50_AK8PFJet150_TrimMass30_NotIsoPFTau50_Trk30v1"); 
  trignames_attempts.push_back("HLT_HT300_MET50_AK8PFJet150_TrimMass20_NotIsoPFTau50_Trk30v1"); 





  if  (HT >350 && met.et() > 70 && Jet.pt()>300 && Jet.pruned_massCorr()>30)  Hist( "Trigger_Attempts"  , directory 	)->Fill(0);//"HLT_HT300_MET50_AK8PFJet200_TrimMass20_v1");
  if  (HT >450 && met.et() > 70 && Jet.pt()>300 && Jet.pruned_massCorr()>30)  Hist( "Trigger_Attempts"  , directory 	)->Fill(1);//"HLT_HT400_MET50_AK8PFJet200_TrimMass20_v1");
  if  (HT >350 && met.et() > 70 && Jet.pt()>300 && Jet.pruned_massCorr()>40)  Hist( "Trigger_Attempts"  , directory 	)->Fill(2);//"HLT_HT300_MET50_AK8PFJet200_TrimMass30_v1");
  if  (HT >400 && met.et() > 70 && Jet.pt()>300 && Jet.pruned_massCorr()>40)  Hist( "Trigger_Attempts"  , directory 	)->Fill(3);//"HLT_HT350_MET50_AK8PFJet200_TrimMass30_v1");


  if  (HT >250 && met.et() > 70 && Jet.pt()>180 && Jet.pruned_massCorr()>30 && tau.pt()>70 && tau.byLooseCombinedIsolationDeltaBetaCorr3Hits() >0.5)  Hist( "Trigger_Attempts"  , directory 	)->Fill(4);//"HLT_HT200_MET50_AK8PFJet120_TrimMass20_LooseIsoPFTau50_Trk30_v1");
  if  (HT >450 && met.et() > 70 && tau.pt()>70 && tau.byLooseCombinedIsolationDeltaBetaCorr3Hits() >0.5)  Hist( "Trigger_Attempts"  , directory 	)->Fill(5);//"HLT_HT400_LooseIsoPFtau50_Trk30_v1");
  if  (HT >500 && met.et() > 70 && tau.pt()>70 && tau.byLooseCombinedIsolationDeltaBetaCorr3Hits() >0.5)  Hist( "Trigger_Attempts"  , directory 	)->Fill(6);//"HLT_HT450_LooseIsoPFTau50_Trk30_v1");
  if  (HT >500 && met.et() > 70 && tau.pt()>50 && tau.byLooseCombinedIsolationDeltaBetaCorr3Hits() >0.5)  Hist( "Trigger_Attempts"  , directory 	)->Fill(7);//"HLT_HT450_LooseIsoPFTau30_Trk30_v1");

  if  (HT >250 && met.et() > 70 && Jet.pt()>180 && Jet.pruned_massCorr()>30 && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(8);//"HLT_HT200_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30_v1");
  if  (HT >350 && met.et() > 70 && Jet.pt()>180 && Jet.pruned_massCorr()>30 && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(9);//"HLT_HT300_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30_v1");
  if  (HT >450 && met.et() > 70 && Jet.pt()>180 && Jet.pruned_massCorr()>30 && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(10);//"HLT_HT400_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30_v1");
  if  (HT >350 && met.et() > 70 && Jet.pt()>300 && Jet.pruned_massCorr()>30 && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(11);//"HLT_HT300_MET50_AK8PFJet200_TrimMass20_NotIsoPFTau50_Trk30_v1");
  if  (HT >350 && met.et() > 70 && Jet.pt()>300 && Jet.pruned_massCorr()>40 && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(12);//"HLT_HT300_MET50_AK8PFJet200_TrimMass30_NotIsoPFTau50_Trk30_v1");
  if  (HT >400 && met.et() > 70 && Jet.pt()>300 && Jet.pruned_massCorr()>40 && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(13);//"HLT_HT350_MET50_AK8PFJet200_TrimMass30_NotIsoPFTau50_Trk30_v1");
  if  (HT >400 && met.et() > 80 && Jet.pt()>300 && Jet.pruned_massCorr()>30 && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(14);//"HLT_HT350_MET60_AK8PFJet200_TrimMass20_NotIsoPFTau50_Trk30_v1");
 if  (HT >400 && met.et() > 70 && Jet.pt()>300 && Jet.pruned_massCorr()>30 && tau.pt()>65)  Hist( "Trigger_Attempts"  , directory 	)->Fill(15);//"HLT_HT350_MET50_AK8PFJet200_TrimMass20_NotIsoPFTau50_v1");
 if  (HT >350 && met.et() > 70 && Jet.pt()>180 && Jet.pruned_massCorr()>40 && tau.pt()>65)  Hist( "Trigger_Attempts"  , directory 	)->Fill(16);//"HLT_HT300_MET50_AK8PFJet120_TrimMass30_NotIsoPFTau50_v1");

 if  (HT >450 && met.et() > 70  && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(17);//"HLT_HT400_MET50_NotIsoPFTau50_Trk30_v1");
 if  (HT >450 && met.et() > 90  && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(18);//"HLT_HT400_MET70_NotIsoPFTau50_Trk30_v1"); 
 if  (HT >500 && met.et() > 70  && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(19);//"HLT_HT450_MET50_NotIsoPFTau50_Trk30_v1");

 if  (HT >350 && met.et() > 70 && Jet.pt()>180 && Jet.pruned_massCorr()>45 && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(20); // HLT_HT300_MET50_AK8PFJet120_TrimMass30_NotIsoPFTau50_Trk30v1 
if  (HT >350 && met.et() > 80 && Jet.pt()>180 && Jet.pruned_massCorr()>45 && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(21);// HLT_HT300_MET60_AK8PFJet120_TrimMass30_NotIsoPFTau50_Trk30v1 
if  (HT >450 && met.et() > 70 && Jet.pt()>180 && Jet.pruned_massCorr()>40 && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(22);// HLT_HT400_MET50_AK8PFJet120_TrimMass20_NotIsoPFTau50_Trk30v1 

 if  (HT >350 && met.et() > 70 && Jet.pt()>210 && Jet.pruned_massCorr()>45 && tau.pt()>70)  Hist( "Trigger_Attempts"  , directory 	)->Fill(23);// HLT_HT300_MET50_AK8PFJet150_TrimMass30_NotIsoPFTau50_Trk30v1  
 if  (HT >350 && met.et() > 80 && Jet.pt()>210 && Jet.pruned_massCorr()>45 && tau.pt()>70) Hist( "Trigger_Attempts"  , directory 	)->Fill(24);// HLT_HT300_MET60_AK8PFJet150_TrimMass30_NotIsoPFTau50_Trk30v1
 if  (HT >400 && met.et() > 70 && Jet.pt()>210 && Jet.pruned_massCorr()>45 && tau.pt()>70) Hist( "Trigger_Attempts"  , directory 	)->Fill(25);// HLT_HT350_MET50_AK8PFJet150_TrimMass30_NotIsoPFTau50_Trk30v1
 if  (HT >350 && met.et() > 70 && Jet.pt()>210 && Jet.pruned_massCorr()>40 && tau.pt()>70) Hist( "Trigger_Attempts"  , directory 	)->Fill(26);// HLT_HT300_MET50_AK8PFJet150_TrimMass20_NotIsoPFTau50_Trk30v1



}

void VHTausAnalysis::FillBranches( const std::string& channel,  const UZH::Jet& Jet, const double& Jet_tau21, const double& Jet_tau31, const double& Jet_tau32, const int& Jet_nTaggedSubjets, const double& Jet_subjet0_csv, const double& Jet_subjet1_csv, const double& ak8jet0_subjet01_dr,  const double& ak8jet0_subjet01_deta, const double& ak8jet0_subjet01_dphi, const double& b_weight, const double& b_weightGen,  const double& b_weightPU, const double& b_weightBtag,  const double& b_weightBtag_veto,  const double& b_runNumber, const double& b_eventNumber,const double& b_lumiBlock,  const UZH::Tau& tau,const  TLorentzVector& lepton, const UZH::MissingEt& met ){
  // std::cout<< "Filling branches" <<std::endl;
  b_weight_[channel.c_str()] = b_weight;
  b_weightGen_[channel.c_str()] =  b_weightGen;
  b_weightPU_[channel.c_str()] =  b_weightPU;
  b_weightBtag_[channel.c_str()] =   b_weightBtag;
  b_weightBtag_veto_[channel.c_str()] =   b_weightBtag_veto;
  
  b_runNumber_[channel.c_str()] = b_runNumber;
  b_eventNumber_[channel.c_str()] = b_eventNumber;
  b_lumiBlock_[channel.c_str()] = b_lumiBlock;
  
  b_ak8jet0_pt_[channel.c_str()] = Jet.pt();
  b_ak8jet0_phi_[channel.c_str()] = Jet.phi() ;
  b_ak8jet0_eta_[channel.c_str()] = Jet.eta();
  b_ak8jet0_e_[channel.c_str()] = Jet.e() ;
  b_ak8jet0_tau21_[channel.c_str()] = Jet_tau21;
  b_ak8jet0_m_[channel.c_str()] =  Jet.m();
  b_ak8jet0_mpruned_[channel.c_str()] =  Jet.pruned_massCorr() ;
  b_ak8jet0_csv_[channel.c_str()] = Jet.csv();
  b_ak8jet0_subjet01_dr_[channel.c_str()] = ak8jet0_subjet01_dr;
  b_ak8jet0_subjet01_deta_[channel.c_str()] = ak8jet0_subjet01_deta;
  b_ak8jet0_subjet01_dphi_[channel.c_str()] = ak8jet0_subjet01_dphi;
  b_ak8jet0_subjet0_pt_[channel.c_str()] = b_ak8jet0_subjet0_pt;
  b_ak8jet0_subjet1_pt_[channel.c_str()] = b_ak8jet0_subjet1_pt;
  b_ak8jet0_subjet0_csv_[channel.c_str()] = Jet_subjet0_csv;
  b_ak8jet0_subjet1_csv_[channel.c_str()] = Jet_subjet1_csv ;


  b_tau_pt_[channel.c_str()]=tau.tlv().Pt();
  b_tau_eta_[channel.c_str()]=tau.tlv().Eta();
  b_tau_phi_[channel.c_str()]=tau.tlv().Phi();
  b_tauIsoRaw_[channel.c_str()]=tau.byIsolationMVArun2v1PWnewDMwLTraw();
  b_tauIso_[channel.c_str()]=tau.byVLooseIsolationMVArun2v1PWnewDMwLT();
  b_tau_againstEle_[channel.c_str()]=tau.againstElectronVLooseMVA6();
  b_tau_againstMu_[channel.c_str()]=tau.againstMuonLoose3();

  b_seclepton_pt_[channel.c_str()]=lepton.Pt();
  b_seclepton_eta_[channel.c_str()]=lepton.Phi();
  b_seclepton_phi_[channel.c_str()]=lepton.Eta();
  
  b_vis_mass_ll_[channel.c_str()]=(lepton +tau.tlv()).M();
  b_dR_ll_[channel.c_str()]= lepton.DeltaR(tau.tlv());
  // if (channel== "mutau") {
  //   std:: cout<<" Branch Dr mutau" <<  lepton.DeltaR(tau.tlv()) <<std::endl;
  //   std::cout << " Branch goodMuons[0].tlv().pt() " << lepton.Pt() <<std::endl;
  //   std::cout << " Branch Tau .pt() " << tau.tlv().Pt() <<std::endl;
  // }
  b_met_[channel.c_str()] =met.et();
  b_met_phi_[channel.c_str()] =met.phi();
  //std::cout << " exiting branch filling " << tau.tlv().Pt() <<std::endl;

  b_N_taggedAK4_L_[channel.c_str()] = taggedAK4_L ;
  b_N_taggedAK4_M_[channel.c_str()] = taggedAK4_M ;
  b_N_taggedAK4_T_[channel.c_str()] = taggedAK4_T ;

  b_isBoostedTau_[channel.c_str()] = isBoostedTau;

  passTrigger_branches(channel.c_str());


  
  bool doSVfit= true;
  if (doSVfit){
    TLorentzVector H_ll;
    H_ll.SetPtEtaPhiM(0,0,0,0);
  
    if (channel=="tautau"){
      H_ll=applySVFitHadronic(met.cov00(),met.cov10(),met.cov11(),met.et(),met.phi(),tau.tlv(),lepton);
    }
    else {
      H_ll=applySVFitSemileptonic(met.cov00(),met.cov10(),met.cov11(),met.et(),met.phi(),tau.tlv(),lepton);
    }

    b_H_Mass_SVFit_[channel.c_str()]=H_ll.M();
    b_H_Pt_SVFit_[channel.c_str()]=H_ll.Pt();
    b_H_Eta_SVFit_[channel.c_str()]=H_ll.Eta();
    b_H_Phi_SVFit_[channel.c_str()]=H_ll.Phi();
    b_X_Mass_SVFit_[channel.c_str()]=(H_ll+Jet.tlv()).M();
    b_X_Pt_SVFit_[channel.c_str()]=(H_ll+Jet.tlv()).Pt();
    b_X_Eta_SVFit_[channel.c_str()]=(H_ll+Jet.tlv()).Eta();
    b_X_Phi_SVFit_[channel.c_str()]=(H_ll+Jet.tlv()).Phi();
  }


  
  // b_selection_bits_[channel].clear();
  // b_selection_lastcut_[channel].clear();

}




TLorentzVector VHTausAnalysis::applySVFitHadronic( float cov00, float cov10, float cov11, float met, float met_phi, TLorentzVector lep1 , TLorentzVector lep2){
  // std::cout<<"inside applySVFitHadronic "<<std::endl;

  
  TLorentzVector   lBoson4;
  lBoson4.SetPtEtaPhiE(0,0,0,0);
  
  TMatrixD covMETraw(2,2);
  covMETraw[0][0]=  cov00;
  covMETraw[1][0]=  cov10;
  covMETraw[0][1]=  cov10;
  covMETraw[1][1]=  cov11;
  float lcov00 =  cov00;
  float lcov10 =  cov10;
  float lcov01 =  cov10;
  float lcov11 =  cov11;
  
  TLorentzVector   nullo;
  nullo.SetPtEtaPhiE(0,0,0,0);
  if(lcov00*lcov11-lcov01*lcov01 == 0) {
    std::cout<<"covMat det null "<<std::endl;
    return nullo;
  }

 
  lBoson4.SetPtEtaPhiE(0,0,0,0);
  NSVfitStandalone::Vector measuredMET(met *TMath::Cos(met_phi), met *TMath::Sin(met_phi), 0); 
  // setup measure tau lepton vectors 
  NSVfitStandalone::LorentzVector l1(lep1.Px(), lep1.Py(), lep1.Pz(), TMath::Sqrt(lep1.M()*lep1.M()+lep1.Px()*lep1.Px()+lep1.Py()*lep1.Py()+lep1.Pz()*lep1.Pz()));
  NSVfitStandalone::LorentzVector l2(lep2.Px(), lep2.Py(), lep2.Pz(), TMath::Sqrt(lep2.M()*lep2.M()+lep2.Px()*lep2.Px()+lep2.Py()*lep2.Py()+lep2.Pz()*lep2.Pz()));
  std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, l1));
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, l2));
  // construct the class object from the minimal necesarry information
  NSVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMET, covMETraw, 0);
  // apply customized configurations if wanted (examples are given below)
  //algo.maxObjFunctionCalls(5000);
  algo.addLogM(false);
  //algo.metPower(0.5)
  //algo.fit();
  // integration by VEGAS (same as function algo.integrate() that has been in use when markov chain integration had not yet been implemented)
  //algo.integrateVEGAS();
  // integration by markov chain MC
  algo.integrateMarkovChain();

  if(algo.isValidSolution()){
    lBoson4.SetPtEtaPhiM( algo.pt(), algo.eta(), algo.phi(), algo.getMass());
  }
  else{
   std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
   lBoson4.SetPtEtaPhiM( 0.,0,0,0);
  }
  
  measuredTauLeptons.clear(); 
  return   lBoson4 ;
}

TLorentzVector VHTausAnalysis::applySVFitSemileptonic(float cov00, float cov10, float cov11,  float met, float met_phi, TLorentzVector lep1 , TLorentzVector lep2){
  // std::cout<<"inside applySVFitSemileptonic "<<std::endl;
  
  TLorentzVector   lBoson4;
  lBoson4.SetPtEtaPhiE(0,0,0,0);

  TMatrixD covMETraw(2,2);
  covMETraw[0][0]=  cov00;
  covMETraw[1][0]=  cov10;
  covMETraw[0][1]=  cov10;
  covMETraw[1][1]=  cov11;

  float lcov00 =  cov00;
  float lcov10 =  cov10;
  float lcov01 =  cov10;
  float lcov11 =  cov11;
  
  
  TLorentzVector   nullo;
  nullo.SetPtEtaPhiE(0,0,0,0);
  if(lcov00*lcov11-lcov01*lcov01 == 0) {
    std::cout<<"covMat det null "<<std::endl;
    return nullo;
  }
  

  NSVfitStandalone::Vector measuredMET(met *TMath::Cos(met_phi), met *TMath::Sin(met_phi), 0); 
  // setup measure tau lepton vectors 
  NSVfitStandalone::LorentzVector l1(lep1.Px(), lep1.Py(), lep1.Pz(), TMath::Sqrt(lep1.M()*lep1.M()+lep1.Px()*lep1.Px()+lep1.Py()*lep1.Py()+lep1.Pz()*lep1.Pz()));
  NSVfitStandalone::LorentzVector l2(lep2.Px(), lep2.Py(), lep2.Pz(), TMath::Sqrt(lep2.M()*lep2.M()+lep2.Px()*lep2.Px()+lep2.Py()*lep2.Py()+lep2.Pz()*lep2.Pz()));
  std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, l1));
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, l2));
  // construct the class object from the minimal necesarry information
  NSVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMET, covMETraw, 0);
  // apply customized configurations if wanted (examples are given below)
  //algo.maxObjFunctionCalls(5000);
  algo.addLogM(false);
  //algo.metPower(0.5)
  //algo.fit();
  // integration by VEGAS (same as function algo.integrate() that has been in use when markov chain integration had not yet been implemented)
  //algo.integrateVEGAS();
  // integration by markov chain MC
  algo.integrateMarkovChain();
  if(algo.isValidSolution()){
   
    lBoson4.SetPtEtaPhiM( algo.pt(), algo.eta(), algo.phi(), algo.getMass());
  }
  else{
    std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
   lBoson4.SetPtEtaPhiM( 0.,0,0,0);
  }
  
  measuredTauLeptons.clear();
  return   lBoson4 ;
  
}
void VHTausAnalysis::checkTrigger(const TString& directory,const std::string& channel
				  ){
  
  // std::string channel;
  // if (directory.find("mutau") != std::string::npos) channel="mutau";
  // if (directory.find("eletau") != std::string::npos) channel="eletau";
  // if (directory.find("tautau") != std::string::npos) channel="tautau";

  std::vector<std::string> trignames;

  trignames.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v");
  trignames.push_back("HLT_AK8PFJet360_TrimMass30_v");
  trignames.push_back("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v");
  trignames.push_back("HLT_AK8PFHT660_TrimR0p1PT0p03Mass50_BTagCSV_p20_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p41_v");
  trignames.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v");
  trignames.push_back("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_v");

  trignames.push_back("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v");

  trignames.push_back("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v");
  trignames.push_back("HLT_PFHT600_v");
  trignames.push_back("HLT_PFHT650_v") ;
  trignames.push_back("HLT_PFHT800_v") ;
  trignames.push_back("HLT_PFHT900_v") ;
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
  trignames.push_back("HLT_Ele15_IsoVVVL_PFHT600");
  trignames.push_back("HLT_Mu15_IsoVVVL_PFHT600");
  trignames.push_back("HLT_Ele15_IsoVVVL_PFHT350_PFMET50");
  trignames.push_back("HLT_Mu15_IsoVVVL_PFHT350_PFMET50");
  trignames.push_back("HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400"); 
  trignames.push_back("HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400");
 
 
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
  trignamesJet.push_back("HLT_PFHT900_v") ;

  
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
      // 	if( (it->first).find(trignamesTau[t]) != std::string::npos && it->second   // && !isfiredJet
      // 	    ){
      // 	  isfiredTau = true;
      // 	  Hist( "Trigger_Tau",directory )->Fill( t ); 
      // 	}
	if( (it->first).find(trignamesTau[t]) != std::string::npos && it->second  // && !isfiredJet
	    ){
	
	  
	  if( trignamesTau[t].find("Ele") != std::string::npos   && channel=="eletau" ){///not trigger ele 105, unisolated
	    Hist( "Trigger_Tau",directory )->Fill( t );	
	    isfiredTau = true;
	    
	  }
	  
	  else if(trignamesTau[t].find("Mu") != std::string::npos && channel=="mutau"){
	    
	    Hist( "Trigger_Tau",directory )->Fill( t );
	    isfiredTau = true;
	  }
	  
	  else if(trignamesTau[t].find("DoubleMediumIsoPFTau") != std::string::npos && channel=="tautau"){
	    
	    Hist( "Trigger_Tau",directory )->Fill( t );	
	    isfiredTau = true;
	  }  
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


void VHTausAnalysis::checkTriggerTurnOn(){

  Book( TH1F( "HLT_AK8PFJet360_TrimMass30_pt","HLT_AK8PFJet360_TrimMass30_pt", 2000 , 0, 4000 ), "Trigger_eff") ;
  Book( TH1F( "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_pt","HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_pt", 2000 , 0, 4000 ), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_pt" ,"HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_pt" , 2000 , 0, 4000 ), "Trigger_eff");
  Book( TH1F( "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_pt","HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_pt", 2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT800_pt","HLT_PFHT800_pt",  2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT900_pt","HLT_PFHT900_pt", 2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT650_pt","HLT_PFHT650_pt", 2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_PFJet320_pt","HLT_PFJet320_pt", 2000 , 0, 4000), "Trigger_eff") ;

  Book( TH1F( "HLT_AK8PFJet360_TrimMass30_softdropMass","HLT_AK8PFJet360_TrimMass30_softdropMass", 200 , 0, 400 ), "Trigger_eff") ;
  Book( TH1F( "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_softdropMass","HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_softdropMass", 200 , 0, 400 ), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_softdropMass","HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_softdropMass" , 200 , 0, 400 ), "Trigger_eff");
  Book( TH1F( "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_softdropMass","HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_softdropMass", 200 , 0, 400), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT800_softdropMass","HLT_PFHT800_softdropMass", 200 , 0, 400), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT900_softdropMass","HLT_PFHT900_softdropMass", 200 , 0, 400), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT650_softdropMass","HLT_PFHT650_softdropMass", 200 , 0, 400), "Trigger_eff") ;
  Book( TH1F( "HLT_PFJet320_softdropMass","HLT_PFJet320_softdropMass", 200 , 0, 400), "Trigger_eff") ;

  Book( TH1F( "HLT_AK8PFJet360_TrimMass30_ht","HLT_AK8PFJet360_TrimMass30_ht", 2000 , 0, 4000 ), "Trigger_eff") ;
  Book( TH1F( "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_ht","HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_ht", 2000 , 0, 4000 ), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_ht","HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_ht" , 2000 , 0, 4000 ), "Trigger_eff");
  Book( TH1F( "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_ht","HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_ht", 2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT800_ht","HLT_PFHT800_ht", 2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT900_ht","HLT_PFHT900_ht", 2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_PFHT650_ht","HLT_PFHT650_ht", 2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_PFJet320_ht","HLT_PFJet320_ht", 2000 , 0, 4000), "Trigger_eff") ;
 
  Book( TH1F( "HLT_OR_pt","HLT_OR_pt", 2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_OR_ht","HLT_OR_ht", 2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_OR_softdropMass","HLT_OR_softdropMass", 200 , 0, 400), "Trigger_eff") ;
  
  Book( TH1F( "HLT_PreOR_pt","HLT_PreOR_pt", 2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_PreOR_ht","HLT_PreOR_ht", 2000 , 0, 4000), "Trigger_eff") ;
  Book( TH1F( "HLT_PreOR_softdropMass","HLT_PreOR_softdropMass", 200 , 0, 400), "Trigger_eff") ;

  std::vector<UZH::Jet> goodFatJets_trig;
  for ( int i = 0; i < (m_jetAK8.N); ++i ) {
    UZH::Jet myjet( &m_jetAK8, i );
    if (myjet.pt() > m_jetPtCut) {
      if (fabs(myjet.eta()) < m_jetEtaCut) {
	if (myjet.IDTight()) {
	  goodFatJets_trig.push_back(myjet);
	}
      }
    }
  }
  UZH::Jet Jet_trig;
  if(goodFatJets_trig.size()>=1) Jet_trig=goodFatJets_trig[0];
 

  float HT_trig=0.;


  for ( int i = 0; i < (m_jetAK4.N); ++i ) {
    UZH::Jet myjetak4( &m_jetAK4, i );
    if (fabs(myjetak4.eta()) < m_AK4jetEtaCut) {
      if (myjetak4.pt() > m_AK4jetPtCut) {


	HT = HT + myjetak4.tlv().Et();
      }
    }
  }


  bool looseTrigHT=false;
  for( std::map<std::string,bool>::iterator it = (m_eventInfo.trigDecision)->begin(); it != (m_eventInfo.trigDecision)->end(); ++it){
    if( (it->first).find("HLT_PFHT650_v") != std::string::npos  && it->second ){  
      if(goodFatJets_trig.size()>=1 && Jet_trig.pruned_massCorr()>50 )   Hist( "HLT_PFHT650_pt", "Trigger_eff")->Fill(Jet_trig.pt()); 
      Hist( "HLT_PFHT650_ht", "Trigger_eff")->Fill(HT); 
      if(goodFatJets_trig.size()>=1 && Jet_trig.pt()> 500  )  Hist( "HLT_PFHT650_softdropMass", "Trigger_eff")->Fill(Jet_trig.pruned_massCorr()) ;
      looseTrigHT=true;
      break; 
    }
  }

  bool looseTrigJet=false;
  
  for( std::map<std::string,bool>::iterator it = (m_eventInfo.trigDecision)->begin(); it != (m_eventInfo.trigDecision)->end(); ++it){
    if( (it->first).find("HLT_PFJet320_v") != std::string::npos  && it->second ){  
      if(goodFatJets_trig.size()>=1 && Jet_trig.pruned_massCorr()>50)   Hist( "HLT_PFJet320_pt", "Trigger_eff")->Fill(Jet_trig.pt()); 
      Hist( "HLT_PFJet320_ht", "Trigger_eff")->Fill(HT); 
      if(goodFatJets_trig.size()>=1  && Jet_trig.pt()> 500)  Hist( "HLT_PFJet320_softdropMass", "Trigger_eff")->Fill(Jet_trig.pruned_massCorr()) ;
      looseTrigJet=true;
      break;
    }
  }




  if (looseTrigJet || looseTrigHT){
    if(goodFatJets_trig.size()>=1 && Jet_trig.pruned_massCorr()>50)   Hist( "HLT_PreOR_pt", "Trigger_eff")->Fill(Jet_trig.pt()); 
    Hist( "HLT_PreOR_ht", "Trigger_eff")->Fill(HT); 
    if(goodFatJets_trig.size()>=1  && Jet_trig.pt()> 500)  Hist( "HLT_PreOR_softdropMass", "Trigger_eff")->Fill(Jet_trig.pruned_massCorr()) ;
  }

  if (  looseTrigJet){
   
    for( std::map<std::string,bool>::iterator it = (m_eventInfo.trigDecision)->begin(); it != (m_eventInfo.trigDecision)->end(); ++it){
      if( (it->first).find("HLT_AK8PFJet360_TrimMass30_v") != std::string::npos  &&  it->second){  
	if(goodFatJets_trig.size()>=1 && Jet_trig.pruned_massCorr()>50 )	Hist( "HLT_AK8PFJet360_TrimMass30_pt", "Trigger_eff")->Fill(Jet_trig.pt());
	Hist( "HLT_AK8PFJet360_TrimMass30_ht", "Trigger_eff")->Fill(HT);
	if(goodFatJets_trig.size()>=1 && Jet_trig.pt()> 500)   Hist( "HLT_AK8PFJet360_TrimMass30_softdropMass", "Trigger_eff")->Fill(Jet_trig.pruned_massCorr()) ;
      }
    } 
  }
    


  
  if ( looseTrigHT ){
    for( std::map<std::string,bool>::iterator it = (m_eventInfo.trigDecision)->begin(); it != (m_eventInfo.trigDecision)->end(); ++it){
     
      if ( (it->first).find("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v") != std::string::npos  &&  it->second){  
	if(goodFatJets_trig.size()>=1 && Jet_trig.pruned_massCorr()>50)	Hist( "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_pt", "Trigger_eff")->Fill(Jet_trig.pt());
	Hist( "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_ht", "Trigger_eff")->Fill(HT);
	if(goodFatJets_trig.size()>=1 && Jet_trig.pt()> 500)	Hist( "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_softdropMass", "Trigger_eff")->Fill(Jet_trig.pruned_massCorr()); 
      }
      if ( (it->first).find("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v") != std::string::npos  &&  it->second){  
	if(goodFatJets_trig.size()>=1 && Jet_trig.pruned_massCorr()>50)	Hist( "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_pt", "Trigger_eff")->Fill(Jet_trig.pt()); 
	Hist( "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_ht", "Trigger_eff")->Fill(HT);
	if(goodFatJets_trig.size()>=1 && Jet_trig.pt()> 500)	Hist( "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_softdropMass", "Trigger_eff")->Fill(Jet_trig.pruned_massCorr()) ;
      }
      if ( (it->first).find("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v") != std::string::npos  &&  it->second){  
	if(goodFatJets_trig.size()>=1 && Jet_trig.pruned_massCorr()>50)	Hist( "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_pt", "Trigger_eff")->Fill(Jet_trig.pt()); 
	Hist( "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_ht", "Trigger_eff")->Fill(HT);
	if(goodFatJets_trig.size()>=1 && Jet_trig.pt()> 500)	Hist( "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_softdropMass", "Trigger_eff")->Fill(Jet_trig.pruned_massCorr()) ;
      }
      if ( (it->first).find("HLT_PFHT800_v") != std::string::npos  &&  it->second){  
	if(goodFatJets_trig.size()>=1 && Jet_trig.pruned_massCorr()>50)	Hist( "HLT_PFHT800_pt", "Trigger_eff")->Fill(Jet_trig.pt()); 
	Hist( "HLT_PFHT800_ht", "Trigger_eff")->Fill(HT);
	if(goodFatJets_trig.size()>=1 && Jet_trig.pt()> 500)	Hist( "HLT_PFHT800_softdropMass", "Trigger_eff")->Fill(Jet_trig.pruned_massCorr()); 
      }
      if ( (it->first).find("HLT_PFHT900_v") != std::string::npos  &&  it->second){  
	if(goodFatJets_trig.size()>=1 && Jet_trig.pruned_massCorr()>50)	Hist( "HLT_PFHT900_pt", "Trigger_eff")->Fill(Jet_trig.pt());
	Hist( "HLT_PFHT900_ht", "Trigger_eff")->Fill(HT);
	if(goodFatJets_trig.size()>=1 && Jet_trig.pt()> 500)	Hist( "HLT_PFHT900_softdropMass", "Trigger_eff")->Fill(Jet_trig.pruned_massCorr());
      }
    }
  }
  
  if ( looseTrigHT || looseTrigJet ){

    for( std::map<std::string,bool>::iterator it = (m_eventInfo.trigDecision)->begin(); it != (m_eventInfo.trigDecision)->end(); ++it){
 
      if( ((it->first).find("HLT_PFHT800_v") != std::string::npos  &&  it->second   && looseTrigHT ) || ( (it->first).find("HLT_AK8PFJet360_TrimMass30_v") != std::string::npos  &&  it->second && looseTrigJet )){
      
	if(goodFatJets_trig.size()>=1 && Jet_trig.pruned_massCorr()>50)   Hist( "HLT_OR_pt", "Trigger_eff")->Fill(Jet_trig.pt()); 
	Hist( "HLT_OR_ht", "Trigger_eff")->Fill(HT); 
	if(goodFatJets_trig.size()>=1  && Jet_trig.pt()> 500)  Hist( "HLT_OR_softdropMass", "Trigger_eff")->Fill(Jet_trig.pruned_massCorr()) ;
	break;
      }

      
    }
  }

  
}


void VHTausAnalysis::passTrigger_branches(const std::string& channel) {
  
  
  for (std::map<std::string,bool>::iterator it = (m_eventInfo.trigDecision)->begin(); it != (m_eventInfo.trigDecision)->end(); ++it){
   
    if ((it->first).find("HLT_PFHT800_v") != std::string::npos) {
      if (it->second == true)  b_trig_HT800_[channel.c_str()]=true;
      else  b_trig_HT800_[channel.c_str()]=false;
    }
    if ((it->first).find("HLT_AK8PFJet360_TrimMass30_v") != std::string::npos) {
      if (it->second == true)  b_trig_ak8Jet360_Trim30_[channel.c_str()]=true;
      else b_trig_ak8Jet360_Trim30_[channel.c_str()]=false;
    }
    
    if ((it->first).find("HLT_Ele27_eta2p1_WP75_Gsf_v") != std::string::npos  || (it->first).find("HLT_Ele27_eta2p1_WPLoose_Gsf_v") != std::string::npos   ) {
      if (it->second == true)  b_trig_Ele27_[channel.c_str()]=true;
      else b_trig_Ele27_[channel.c_str()]=false;
    }  
    if ((it->first).find("HLT_Ele23_WPLoose_Gsf_v") ){
      if (it->second == true)  b_trig_Ele23_[channel.c_str()]=true;
      else b_trig_Ele23_[channel.c_str()]=false;
    }

    if ((it->first).find("HLT_IsoMu20_v") != std::string::npos || (it->first).find("HLT_IsoTkMu20_v") != std::string::npos) {
      if (it->second == true)  b_trig_Mu20_[channel.c_str()]=true;
      else b_trig_Mu20_[channel.c_str()]=false;
    }
    if ((it->first).find("HLT_IsoMu18_v") != std::string::npos ) {
      if (it->second == true)  b_trig_Mu18_[channel.c_str()]=true;
      else b_trig_Mu18_[channel.c_str()]=false;
    }
      
  }
}

