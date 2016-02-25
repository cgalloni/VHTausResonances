// Dear emacs, this is -*- c++ -*-
// $Id: CycleCreators.py 344 2012-12-13 13:10:53Z krasznaa $
#ifndef VHTausAnalysis_H
#define VHTausAnalysis_H

// SFrame include(s):
#include "core/include/SCycleBase.h"

// ROOT include(s):
#include <TBits.h>

// External include(s):
#include "../NtupleVariables/include/JetNtupleObject.h"
#include "../NtupleVariables/include/Jet.h"
#include "../NtupleVariables/include/EventInfoNtupleObject.h"
#include "../NtupleVariables/include/ElectronNtupleObject.h"
#include "../NtupleVariables/include/Electron.h"
#include "../NtupleVariables/include/MuonNtupleObject.h"
#include "../NtupleVariables/include/Muon.h"
#include "../NtupleVariables/include/TauNtupleObject.h"
#include "../NtupleVariables/include/Tau.h"
#include "../NtupleVariables/include/MissingEtNtupleObject.h"
#include "../NtupleVariables/include/MissingEt.h"
#include "../NtupleVariables/include/GenParticleNtupleObject.h"
#include "../NtupleVariables/include/GenParticle.h"
#include "../GoodRunsLists/include/TGoodRunsList.h"
#include "../PileupReweightingTool/include/PileupReweightingTool.h"
#include "../BTaggingTools/include/BTaggingScaleTool.h"

class TH1D;
class TH2D;
class TRandom3;
class TBits;
namespace UZH {
  class Jet;
  class Electron;
  class Muon;
  class Tau;
  class MissingEt;
  class GenParticle;
}

/**
 *   @short Put short description of class here
 *
 *          Put a longer description over here...
 *
 *  @author Put your name here
 * @version $Revision: 344 $
 */
class VHTausAnalysis : public SCycleBase {

public:
  
  /// Enumeration of all cut flags
  typedef enum {
    kBeforeCuts,            // C0
    kJSON,                  // C1
    kTrigger,               // C2
    kMetFilters,            // C3
    kMet,                   // C4
    kJet,                   // C5
    kTau,                   // C6
    kLepton,                // C7
    kTauIsolation,          // C8
    kMassWindow,               // C9
    // kMassWindow,           // C10
    kTau21,                 // C11
    kSubjetSingleTag,       // C12
    kSubjetDoubleTag,       // C13
    kNumCuts                // last!
  } SelectionCuts;
  
  /// static array of all cut names
  static const std::string kCutName[ kNumCuts ];
  
   /// Default constructor
   VHTausAnalysis();
   /// Default destructor
   ~VHTausAnalysis();

   /// Function called at the beginning of the cycle
   virtual void BeginCycle() throw( SError );
   /// Function called at the end of the cycle
   virtual void EndCycle() throw( SError );

   /// Function called at the beginning of a new input data
   virtual void BeginInputData( const SInputData& ) throw( SError );
   /// Function called after finishing to process an input data
   virtual void EndInputData  ( const SInputData& ) throw( SError );

   /// Function called after opening each new input file
   virtual void BeginInputFile( const SInputData& ) throw( SError );

   /// Function called for every event
   virtual void ExecuteEvent( const SInputData&, Double_t ) throw( SError );
   
   /// Function to check good lumi section
   virtual bool isGoodEvent(int runNumber, int lumiSection);
   
   /// Function to check for trigger pass
   virtual bool passTrigger();
   
   /// Function to check for MET filters pass
   virtual bool passMETFilters();
   
   /// Function to obtain event weights for MC
   virtual double getEventWeight();
   
   /// Function to clear/reset all output branches
   virtual void clearBranches();
   
   /// Function to fill cut flow
   virtual void fillCutflow( const std::string histName, const std::string dirName, const TBits& cutmap, const Double_t weight );
   
   /// Function to book histograms
   virtual void bookHistograms( const TString& directory );
   
  /// Function to book histograms for trigger choice
  virtual void bookTriggerHistos(const TString& directory );
  virtual void checkTrigger(const TString& directory );
 

   /// Function to fill histograms
   virtual void fillHistograms( const TString& directory, const UZH::Jet& vectorJet, const UZH::Jet& higgsJet, const TLorentzVector& diJet, const double& vJet_tau21, const double& vJet_tau31, const double& vJet_tau32, const int& vJet_nTaggedSubjets, const double& vJet_subjet0_csv, const double& vJet_subjet1_csv, const double& hJet_tau21, const double& hJet_tau31, const double& hJet_tau32, const int& hJet_nTaggedSubjets, const double& hJet_subjet0_csv, const double& hJet_subjet1_csv, const double& deta, const double& dphi, const double& dr );

private:
   //
   // Put all your private variables here
   //
  
  //
  // Input variable objects:
  //
  Ntuple::JetNtupleObject         m_jetAK4;            ///< jet container
  Ntuple::JetNtupleObject         m_jetAK8;            ///< jet container
  Ntuple::EventInfoNtupleObject   m_eventInfo; ///< event info container
  Ntuple::ElectronNtupleObject    m_electron;            ///< electron container
  Ntuple::MuonNtupleObject        m_muon;            ///< muon container
  Ntuple::TauNtupleObject         m_tau;            ///< tau container
  Ntuple::MissingEtNtupleObject   m_missingEt;            ///< missing E_T container
  Ntuple::GenParticleNtupleObject m_genParticle;            ///< gen particle container
  
  
  //
  // Further objects
  //
  Root::TGoodRunsList m_grl;
  PileupReweightingTool m_pileupReweightingTool;
  BTaggingScaleTool m_bTaggingScaleTool;
  
  //
  // XML settings for VHAnalysis
  //
  // naming
  std::string m_recoTreeName;       ///< name of tree with reconstructed objects in ntuple
  std::string m_outputTreeName;    ///< name of output tree
  int m_ntupleLevel;               ///< cut at which branches for ntuple are written out
  std::string m_jetAK4Name;            ///< name of AK4 jet collection in tree with reconstructed objects
  std::string m_jetAK8Name;       ///< name of AK8 jet collection in tree with reconstructed objects
  std::string m_electronName;       ///< name of electron collection in tree with reconstructed objects
  std::string m_muonName;       ///< name of muon collection in tree with reconstructed objects
  std::string m_tauName;       ///< name of tau collection in tree with reconstructed objects
  std::string m_missingEtName;       ///< name of missing E_T collection in tree with reconstructed objects
  std::string m_genParticleName;       ///< name of gen particle collection in tree with reconstructed objects
  
  // flags
  bool        m_isData;
  bool        m_isSignal;
  bool        m_applyMETFilters;
  
  // cuts
  // jets
  double      m_jetPtCut;         ///< cut on jet pT
  double      m_jetEtaCut;        ///< cut on jet eta
  double      m_mjjCut;           ///< cut on dijet mass
  double      m_jetDeltaEtaCut;   ///< cut on DeltaEta between selected jets
  // substructure cuts
  double      m_tau21HPCut;
  double      m_tau21LPCut;
  double      m_mVLowSidebandCut;
  double      m_mWLowerCut;
  double      m_mWUpperCut;
  double      m_mZLowerCut;
  double      m_mZUpperCut;
  double      m_mHLowerCut;
  double      m_mHUpperCut;
  // b-tagging
  double      m_csvMin;
  
  
  // electrons
  double         m_electronPtCut;
  double         m_electronEtaCut;

  // muons
  double         m_muonPtCut;
  double         m_muonEtaCut;
  
  // leptons=muons and electrons
  double         m_leptonPtCut;
  double         m_leptonEtaCut;

 // taus
  double         m_tauPtCut;
  double         m_tauEtaCut;
  bool           m_tauIsoCut;

  // MET
  double         m_metCut;
   
  // further settings
  std::string m_jsonName;
  
  // other variables needed
  std::vector<std::string> m_triggerNames;
  std::vector<std::string> m_catNames;
  
  ///
  /// branches
  ///
  double b_weight;
  double b_weightGen;
  double b_weightPU;
  double b_weightBtag;
  
  // event variables
  int b_runNumber;
  int b_eventNumber;
  int b_lumiBlock;
  
  // use flat variables for direct input in TMVA
  Double_t   b_ak8jet0_pt;
  Double_t   b_ak8jet0_phi;
  Double_t   b_ak8jet0_eta;
  Double_t   b_ak8jet0_e;
  Double_t   b_ak8jet0_tau21;
  Double_t   b_ak8jet0_m;
  Double_t   b_ak8jet0_mpruned;
  Double_t   b_ak8jet0_csv;
  Double_t   b_ak8jet1_pt;
  Double_t   b_ak8jet1_phi;
  Double_t   b_ak8jet1_eta;
  Double_t   b_ak8jet1_e;
  Double_t   b_ak8jet1_tau21;
  Double_t   b_ak8jet1_m;
  Double_t   b_ak8jet1_mpruned;
  Double_t   b_ak8jet1_csv;
  Double_t   b_ak8jet0_subjet01_dr;
  Double_t   b_ak8jet0_subjet01_deta;
  Double_t   b_ak8jet0_subjet01_dphi;
  Double_t   b_ak8jet0_subjet0_pt;
  Double_t   b_ak8jet0_subjet1_pt;
  Double_t   b_ak8jet0_subjet0_csv;
  Double_t   b_ak8jet0_subjet1_csv;
  Double_t   b_ak8jet1_subjet01_dr;
  Double_t   b_ak8jet1_subjet01_deta;
  Double_t   b_ak8jet1_subjet01_dphi;
  Double_t   b_ak8jet1_subjet0_pt;
  Double_t   b_ak8jet1_subjet1_pt;
  Double_t   b_ak8jet1_subjet0_csv;
  Double_t   b_ak8jet1_subjet1_csv;
  
  std::vector<Int_t>      b_selection_bits;
  std::vector<Int_t>      b_selection_lastcut;

   // Macro adding the functions for dictionary generation
   ClassDef( VHTausAnalysis, 0 );

}; // class VHAnalysis

#endif // VHAnalysis_H

