#ifndef HZAnalyzer_h
#define HZAnalyzer_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <set>
#include "TTree.h"
#include "TFile.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJetStructureBase.hh>
#include <fastjet/contrib/ValenciaPlugin.hh>
#include <fastjet/contrib/AxesDefinition.hh>
#include <fastjet/contrib/Nsubjettiness.hh>
#include <fastjet/contrib/EnergyCorrelator.hh>

typedef std::vector< fastjet::PseudoJet > PseudoJetList;

using namespace lcio ;
using namespace marlin ;

class FastJetUtil;


class HZAnalyzer : public Processor {
                
public:
        
    virtual Processor*  newProcessor() { return new HZAnalyzer ; }
    
    HZAnalyzer() ;

    HZAnalyzer(const HZAnalyzer&) = delete;
    HZAnalyzer& operator=(const HZAnalyzer&) = delete;
    
    // Initialisation - run at the beginning to start histograms, etc.
    virtual void init() ;
    
    // Called at the beginning of every run
    virtual void processRunHeader( LCRunHeader* run ) ;
    
    // Run over each event - the main algorithm
    virtual void processEvent( LCEvent * evt ) ;
    
    // Run at the end of each event
    virtual void check( LCEvent * evt ) ;
    
    // Called at the very end for cleanup, histogram saving, etc.
    virtual void end() ;
    
    friend class FastJetUtil; 
        
protected:
        
    //FastJetUtil* _fju=NULL;

    std::string m_inputMCParticleCollection="";
    std::string m_inputMCJetParticleCollection="";
    std::string m_inputRECOJetParticleCollection="";
    std::string m_inputRECOParticleCollection="";
    std::string m_rootFileName="";
    std::string m_inputTauCollection="";

    std::string m_recorefinedjet0ColName="";
    std::string m_recorefinedjet1ColName="";

    std::string m_vtx_rfj0ColName="";
    std::string m_vtx_rfj1ColName="";

    std::string m_recorefinedjet_4jets_ColName="";
    std::string m_vtx_rfj_4jets_ColName="";

    //isolated particles
    std::string m_recoIsoPartColName="";
    std::string m_genIsoPartColName="";
    std::string m_outputMCTrueLepPhParticleCollection="";

    float m_angleIso = 10.2;
    float m_genTrueLepPhEMin = 7.5;
    float m_R=0.7;
    float m_beta=1.0;
    float m_gamma=1.0;

    float m_cutPtTrackMin=-1.0;
 
    // Call to get collections
    void getCollection(LCCollection*&, std::string, LCEvent*);

    // Run and event counters
    int m_eventNumber=0;
    int m_runNumber=0;
    int m_eventCount=0;
    bool m_saveMEInfo=true;

    TFile * m_rootFile=NULL;
    TTree* m_outputTree=NULL;


    std::vector<float>* m_trueME_E=NULL;
    std::vector<float>* m_trueME_Px=NULL;
    std::vector<float>* m_trueME_Py=NULL;
    std::vector<float>* m_trueME_Pz=NULL;
    std::vector<int>* m_trueME_PDGID=NULL;
    //std::vector<int>* m_trueME_index=NULL;

    std::vector<float>* m_genTrueLepPh_E=NULL;
    std::vector<float>* m_genTrueLepPh_Px=NULL;
    std::vector<float>* m_genTrueLepPh_Py=NULL;
    std::vector<float>* m_genTrueLepPh_Pz=NULL;
    std::vector<int>* m_genTrueLepPh_PDGID=NULL;
    //actual ancestor, i.e. Higgs,W,tau, or incoming e+ or e-
    std::vector<int>* m_genTrueLepPh_ANC_PDGID=NULL;
    std::vector<float>* m_genTrueLepPh_relIso=NULL;
    std::vector<float>* m_genTrueLepPh_relIsoCH=NULL;
    std::vector<float>* m_genTrueLepPh_relIsoPh=NULL;
    std::vector<float>* m_genTrueLepPh_relIsoNH=NULL;
    std::vector<float>* m_genTrueLepPh_relIsoEl=NULL;
    std::vector<float>* m_genTrueLepPh_relIsoMu=NULL;

    //std::vector<float>* m_genjet_jetChargeE_kappa_0_25=NULL;
    //std::vector<float>* m_genjet_jetChargeE_kappa_0_50=NULL;
    //std::vector<float>* m_genjet_jetChargeE_kappa_0_75=NULL;
    //std::vector<float>* m_genjet_jetChargeE_kappa_1_00=NULL;
    //std::vector<float>* m_genjet_jetChargeE_kappa_0_10=NULL;
    //std::vector<float>* m_genjet_jetChargeE_kappa_0_15=NULL;
    //std::vector<float>* m_genjet_jetChargeE_kappa_0_20=NULL;
    //std::vector<float>* m_genjet_jetChargeE_kappa_0_30=NULL;

    //std::vector<float>* m_genjet_jetChargePt_kappa_0_25=NULL;
    //std::vector<float>* m_genjet_jetChargePt_kappa_0_50=NULL;
    //std::vector<float>* m_genjet_jetChargePt_kappa_0_75=NULL;
    //std::vector<float>* m_genjet_jetChargePt_kappa_1_00=NULL;
    //std::vector<float>* m_genjet_jetChargePt_kappa_0_10=NULL;
    //std::vector<float>* m_genjet_jetChargePt_kappa_0_15=NULL;
    //std::vector<float>* m_genjet_jetChargePt_kappa_0_20=NULL;
    //std::vector<float>* m_genjet_jetChargePt_kappa_0_30=NULL;

    //std::vector<float>* m_genjet_jetChargePProj_kappa_0_25=NULL;
    //std::vector<float>* m_genjet_jetChargePProj_kappa_0_50=NULL;
    //std::vector<float>* m_genjet_jetChargePProj_kappa_0_75=NULL;
    //std::vector<float>* m_genjet_jetChargePProj_kappa_1_00=NULL;
    //std::vector<float>* m_genjet_jetChargePProj_kappa_0_10=NULL;
    //std::vector<float>* m_genjet_jetChargePProj_kappa_0_15=NULL;
    //std::vector<float>* m_genjet_jetChargePProj_kappa_0_20=NULL;
    //std::vector<float>* m_genjet_jetChargePProj_kappa_0_30=NULL;

    std::vector<float>* m_genjet_E=NULL;
    std::vector<float>* m_genjet_Px=NULL;
    std::vector<float>* m_genjet_Py=NULL;
    std::vector<float>* m_genjet_Pz=NULL;
    std::vector<int>* m_genjet_Mult=NULL;
    std::vector<int>* m_genjet_NCH=NULL;
    std::vector<int>* m_genjet_NCH_trackPtMin=NULL;
    std::vector<int>* m_genjet_NPh=NULL;
    std::vector<int>* m_genjet_NNH=NULL;
    std::vector<float>* m_genjet_dij_21=NULL;
    std::vector<float>* m_genjet_dij_32=NULL;
    std::vector<float>* m_genjet_dij_43=NULL;
    std::vector<float>* m_genjet_dij_21_max=NULL;
    std::vector<float>* m_genjet_dij_32_max=NULL;
    std::vector<float>* m_genjet_dij_43_max=NULL;
    std::vector<float>* m_genjet_CHFraction=NULL;
    std::vector<float>* m_genjet_CHFraction_trackPtMin=NULL;
    std::vector<float>* m_genjet_PhFraction=NULL;
    std::vector<float>* m_genjet_ElFraction=NULL;
    std::vector<float>* m_genjet_MuFraction=NULL;
    std::vector<float>* m_genjet_NHFraction=NULL;
    std::vector<float>* m_genjet_nsubjettiness1=NULL;
    std::vector<float>* m_genjet_nsubjettiness2=NULL;
    std::vector<float>* m_genjet_nsubjettiness3=NULL;
    std::vector<float>* m_genjet_nsubjettiness1_lrz=NULL;
    std::vector<float>* m_genjet_nsubjettiness2_lrz=NULL;
    std::vector<float>* m_genjet_nsubjettiness3_lrz=NULL;

    std::vector<float>* m_genjet_beta1_ECorr2=NULL;
    std::vector<float>* m_genjet_beta1_ECorr3=NULL;
    std::vector<float>* m_genjet_beta1_N2=NULL;
    std::vector<float>* m_genjet_beta1_N3=NULL;
    std::vector<float>* m_genjet_beta1_C2=NULL;
    std::vector<float>* m_genjet_beta1_C3=NULL;
    std::vector<float>* m_genjet_beta1_D2=NULL;
    std::vector<float>* m_genjet_beta1_ECorr2_E_theta=NULL;
    std::vector<float>* m_genjet_beta1_ECorr3_E_theta=NULL;
    std::vector<float>* m_genjet_beta1_N2_E_theta=NULL;
    std::vector<float>* m_genjet_beta1_N3_E_theta=NULL;
    std::vector<float>* m_genjet_beta1_C2_E_theta=NULL;
    std::vector<float>* m_genjet_beta1_C3_E_theta=NULL;
    std::vector<float>* m_genjet_beta1_D2_E_theta=NULL;

    std::vector<float>* m_genjet_beta2_ECorr2=NULL;
    std::vector<float>* m_genjet_beta2_ECorr3=NULL;
    std::vector<float>* m_genjet_beta2_N2=NULL;
    std::vector<float>* m_genjet_beta2_N3=NULL;
    std::vector<float>* m_genjet_beta2_C2=NULL;
    std::vector<float>* m_genjet_beta2_C3=NULL;
    std::vector<float>* m_genjet_beta2_D2=NULL;
    std::vector<float>* m_genjet_beta2_ECorr2_E_theta=NULL;
    std::vector<float>* m_genjet_beta2_ECorr3_E_theta=NULL;
    std::vector<float>* m_genjet_beta2_N2_E_theta=NULL;
    std::vector<float>* m_genjet_beta2_N3_E_theta=NULL;
    std::vector<float>* m_genjet_beta2_C2_E_theta=NULL;
    std::vector<float>* m_genjet_beta2_C3_E_theta=NULL;
    std::vector<float>* m_genjet_beta2_D2_E_theta=NULL;

    std::vector<float>* m_genjet_beta0_5_ECorr2=NULL;
    std::vector<float>* m_genjet_beta0_5_ECorr3=NULL;
    std::vector<float>* m_genjet_beta0_5_N2=NULL;
    std::vector<float>* m_genjet_beta0_5_N3=NULL;
    std::vector<float>* m_genjet_beta0_5_C2=NULL;
    std::vector<float>* m_genjet_beta0_5_C3=NULL;
    std::vector<float>* m_genjet_beta0_5_D2=NULL;
    std::vector<float>* m_genjet_beta0_5_ECorr2_E_theta=NULL;
    std::vector<float>* m_genjet_beta0_5_ECorr3_E_theta=NULL;
    std::vector<float>* m_genjet_beta0_5_N2_E_theta=NULL;
    std::vector<float>* m_genjet_beta0_5_N3_E_theta=NULL;
    std::vector<float>* m_genjet_beta0_5_C2_E_theta=NULL;
    std::vector<float>* m_genjet_beta0_5_C3_E_theta=NULL;
    std::vector<float>* m_genjet_beta0_5_D2_E_theta=NULL;

    std::vector<float>* m_genjet_subjet_E=NULL;
    std::vector<float>* m_genjet_subjet_Px=NULL;
    std::vector<float>* m_genjet_subjet_Py=NULL;
    std::vector<float>* m_genjet_subjet_Pz=NULL;
    std::vector<int>* m_genjet_subjet_NCH=NULL;
    std::vector<int>* m_genjet_subjet_NCH_trackPtMin=NULL;
    std::vector<int>* m_genjet_subjet_NPh=NULL;
    std::vector<int>* m_genjet_subjet_NNH=NULL;
    std::vector<int>* m_genjet_subjet_jetindex=NULL;
    std::vector<float>* m_genjet_subjet_CHFraction=NULL;
    std::vector<float>* m_genjet_subjet_CHFraction_trackPtMin=NULL;
    std::vector<float>* m_genjet_subjet_PhFraction=NULL;
    std::vector<float>* m_genjet_subjet_ElFraction=NULL;
    std::vector<float>* m_genjet_subjet_MuFraction=NULL;
    std::vector<float>* m_genjet_subjet_NHFraction=NULL;

    std::vector<float>* m_genjet_subjet_jetChargeE_kappa_0_25=NULL;
    std::vector<float>* m_genjet_subjet_jetChargeE_kappa_0_50=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargeE_kappa_0_75=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargeE_kappa_1_00=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargeE_kappa_0_10=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargeE_kappa_0_15=NULL;
    std::vector<float>* m_genjet_subjet_jetChargeE_kappa_0_20=NULL;
    std::vector<float>* m_genjet_subjet_jetChargeE_kappa_0_30=NULL;

    std::vector<float>* m_genjet_subjet_jetChargePt_kappa_0_25=NULL;
    std::vector<float>* m_genjet_subjet_jetChargePt_kappa_0_50=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargePt_kappa_0_75=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargePt_kappa_1_00=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargePt_kappa_0_10=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargePt_kappa_0_15=NULL;
    std::vector<float>* m_genjet_subjet_jetChargePt_kappa_0_20=NULL;
    std::vector<float>* m_genjet_subjet_jetChargePt_kappa_0_30=NULL;

    //std::vector<float>* m_genjet_subjet_jetChargePProj_kappa_0_25=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargePProj_kappa_0_50=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargePProj_kappa_0_75=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargePProj_kappa_1_00=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargePProj_kappa_0_10=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargePProj_kappa_0_15=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargePProj_kappa_0_20=NULL;
    //std::vector<float>* m_genjet_subjet_jetChargePProj_kappa_0_30=NULL;

    std::vector<float>* m_isoPartGenDR10_E=NULL;
    std::vector<float>* m_isoPartGenDR10_Px=NULL;
    std::vector<float>* m_isoPartGenDR10_Py=NULL;
    std::vector<float>* m_isoPartGenDR10_Pz=NULL;
    std::vector<int>*   m_isoPartGenDR10_PDGID=NULL;
    std::vector<float>* m_isoPartGenDR10_relIso=NULL;
    std::vector<float>* m_isoPartGenDR10_relIsoCH=NULL;
    std::vector<float>* m_isoPartGenDR10_relIsoPh=NULL;
    std::vector<float>* m_isoPartGenDR10_relIsoNH=NULL;
    std::vector<float>* m_isoPartGenDR10_relIsoEl=NULL;
    std::vector<float>* m_isoPartGenDR10_relIsoMu=NULL;

    std::vector<float>* m_isoPartRecoDR10_E=NULL;
    std::vector<float>* m_isoPartRecoDR10_Px=NULL;
    std::vector<float>* m_isoPartRecoDR10_Py=NULL;
    std::vector<float>* m_isoPartRecoDR10_Pz=NULL;
    std::vector<int>*   m_isoPartRecoDR10_PDGID=NULL;
    std::vector<float>* m_isoPartRecoDR10_relIso=NULL;
    std::vector<float>* m_isoPartRecoDR10_relIsoCH=NULL;
    std::vector<float>* m_isoPartRecoDR10_relIsoPh=NULL;
    std::vector<float>* m_isoPartRecoDR10_relIsoNH=NULL;
    std::vector<float>* m_isoPartRecoDR10_relIsoEl=NULL;
    std::vector<float>* m_isoPartRecoDR10_relIsoMu=NULL;

    //std::vector<float>* m_recojet_jetChargeE_kappa_0_25=NULL;
    //std::vector<float>* m_recojet_jetChargeE_kappa_0_50=NULL;
    //std::vector<float>* m_recojet_jetChargeE_kappa_0_75=NULL;
    //std::vector<float>* m_recojet_jetChargeE_kappa_1_00=NULL;
    //std::vector<float>* m_recojet_jetChargeE_kappa_0_10=NULL;
    //std::vector<float>* m_recojet_jetChargeE_kappa_0_15=NULL;
    //std::vector<float>* m_recojet_jetChargeE_kappa_0_20=NULL;
    //std::vector<float>* m_recojet_jetChargeE_kappa_0_30=NULL;

    //std::vector<float>* m_recojet_jetChargePt_kappa_0_25=NULL;
    //std::vector<float>* m_recojet_jetChargePt_kappa_0_50=NULL;
    //std::vector<float>* m_recojet_jetChargePt_kappa_0_75=NULL;
    //std::vector<float>* m_recojet_jetChargePt_kappa_1_00=NULL;
    //std::vector<float>* m_recojet_jetChargePt_kappa_0_10=NULL;
    //std::vector<float>* m_recojet_jetChargePt_kappa_0_15=NULL;
    //std::vector<float>* m_recojet_jetChargePt_kappa_0_20=NULL;
    //std::vector<float>* m_recojet_jetChargePt_kappa_0_30=NULL;

    //std::vector<float>* m_recojet_jetChargePProj_kappa_0_25=NULL;
    //std::vector<float>* m_recojet_jetChargePProj_kappa_0_50=NULL;
    //std::vector<float>* m_recojet_jetChargePProj_kappa_0_75=NULL;
    //std::vector<float>* m_recojet_jetChargePProj_kappa_1_00=NULL;
    //std::vector<float>* m_recojet_jetChargePProj_kappa_0_10=NULL;
    //std::vector<float>* m_recojet_jetChargePProj_kappa_0_15=NULL;
    //std::vector<float>* m_recojet_jetChargePProj_kappa_0_20=NULL;
    //std::vector<float>* m_recojet_jetChargePProj_kappa_0_30=NULL;


    //std::vector<float>* m_recojet_rfj_4jets_E=NULL;
    //std::vector<float>* m_recojet_rfj_4jets_Px=NULL;
    //std::vector<float>* m_recojet_rfj_4jets_Py=NULL;
    //std::vector<float>* m_recojet_rfj_4jets_Pz=NULL;
    //std::vector<int>* m_recojet_rfj_4jets_Mult=NULL;
    //std::vector<int>* m_recojet_rfj_4jets_NCH=NULL;

    //std::vector<float>* m_recojet_rfj_4jets_CHFraction=NULL;
    //std::vector<float>* m_recojet_rfj_4jets_PhFraction=NULL;
    //std::vector<float>* m_recojet_rfj_4jets_ElFraction=NULL;
    //std::vector<float>* m_recojet_rfj_4jets_MuFraction=NULL;
    //std::vector<float>* m_recojet_rfj_4jets_NHFraction=NULL;

    //std::vector<float>* m_recojet_rfj_4jets_BTag=NULL;
    //std::vector<float>* m_recojet_rfj_4jets_CTag=NULL;
    //std::vector<float>* m_recojet_rfj_4jets_OTag=NULL;
    //std::vector<int>* m_recojet_rfj_4jets_cat=NULL;

    //std::vector<float>* m_recojet_rfj_4jets_svtx_r=NULL;
    //std::vector<float>* m_recojet_rfj_4jets_svtx_E=NULL;
    //std::vector<float>* m_recojet_rfj_4jets_svtx_Mass=NULL;
    //std::vector<int>* m_recojet_rfj_4jets_svtx_nTrack=NULL;
    //std::vector<int>* m_recojet_rfj_4jets_svtx_Charge=NULL;

    //std::vector<float>* m_recojet_fat2j_rfj4_E=NULL;
    //std::vector<float>* m_recojet_fat2j_rfj4_Px=NULL;
    //std::vector<float>* m_recojet_fat2j_rfj4_Py=NULL;
    //std::vector<float>* m_recojet_fat2j_rfj4_Pz=NULL;
    //std::vector<int>* m_recojet_fat2j_rfj4_Mult=NULL;

    std::vector<float>* m_recojet_E=NULL;
    std::vector<float>* m_recojet_Px=NULL;
    std::vector<float>* m_recojet_Py=NULL;
    std::vector<float>* m_recojet_Pz=NULL;
    std::vector<int>* m_recojet_Mult=NULL;
    std::vector<int>* m_recojet_NCH=NULL;
    std::vector<int>* m_recojet_NCH_trackPtMin=NULL;
    std::vector<int>* m_recojet_NPh=NULL;
    std::vector<int>* m_recojet_NNH=NULL;
    std::vector<float>* m_recojet_dij_21=NULL;
    std::vector<float>* m_recojet_dij_32=NULL;
    std::vector<float>* m_recojet_dij_43=NULL;
    std::vector<float>* m_recojet_dij_21_max=NULL;
    std::vector<float>* m_recojet_dij_32_max=NULL;
    std::vector<float>* m_recojet_dij_43_max=NULL;
    std::vector<float>* m_recojet_CHFraction=NULL;
    std::vector<float>* m_recojet_CHFraction_trackPtMin=NULL;
    std::vector<float>* m_recojet_PhFraction=NULL;
    std::vector<float>* m_recojet_ElFraction=NULL;
    std::vector<float>* m_recojet_MuFraction=NULL;
    std::vector<float>* m_recojet_NHFraction=NULL;
    std::vector<float>* m_recojet_nsubjettiness1=NULL;
    std::vector<float>* m_recojet_nsubjettiness2=NULL;
    std::vector<float>* m_recojet_nsubjettiness3=NULL;
    std::vector<float>* m_recojet_nsubjettiness1_lrz=NULL;
    std::vector<float>* m_recojet_nsubjettiness2_lrz=NULL;
    std::vector<float>* m_recojet_nsubjettiness3_lrz=NULL;
    std::vector<float>* m_recojet_beta1_ECorr2=NULL;
    std::vector<float>* m_recojet_beta1_ECorr3=NULL;
    std::vector<float>* m_recojet_beta1_N2=NULL;
    std::vector<float>* m_recojet_beta1_N3=NULL;
    std::vector<float>* m_recojet_beta1_C2=NULL;
    std::vector<float>* m_recojet_beta1_C3=NULL;
    std::vector<float>* m_recojet_beta1_D2=NULL;
    std::vector<float>* m_recojet_beta1_ECorr2_E_theta=NULL;
    std::vector<float>* m_recojet_beta1_ECorr3_E_theta=NULL;
    std::vector<float>* m_recojet_beta1_N2_E_theta=NULL;
    std::vector<float>* m_recojet_beta1_N3_E_theta=NULL;
    std::vector<float>* m_recojet_beta1_C2_E_theta=NULL;
    std::vector<float>* m_recojet_beta1_C3_E_theta=NULL;
    std::vector<float>* m_recojet_beta1_D2_E_theta=NULL;

    std::vector<float>* m_recojet_beta2_ECorr2=NULL;
    std::vector<float>* m_recojet_beta2_ECorr3=NULL;
    std::vector<float>* m_recojet_beta2_N2=NULL;
    std::vector<float>* m_recojet_beta2_N3=NULL;
    std::vector<float>* m_recojet_beta2_C2=NULL;
    std::vector<float>* m_recojet_beta2_C3=NULL;
    std::vector<float>* m_recojet_beta2_D2=NULL;
    std::vector<float>* m_recojet_beta2_ECorr2_E_theta=NULL;
    std::vector<float>* m_recojet_beta2_ECorr3_E_theta=NULL;
    std::vector<float>* m_recojet_beta2_N2_E_theta=NULL;
    std::vector<float>* m_recojet_beta2_N3_E_theta=NULL;
    std::vector<float>* m_recojet_beta2_C2_E_theta=NULL;
    std::vector<float>* m_recojet_beta2_C3_E_theta=NULL;
    std::vector<float>* m_recojet_beta2_D2_E_theta=NULL;

    std::vector<float>* m_recojet_beta0_5_ECorr2=NULL;
    std::vector<float>* m_recojet_beta0_5_ECorr3=NULL;
    std::vector<float>* m_recojet_beta0_5_N2=NULL;
    std::vector<float>* m_recojet_beta0_5_N3=NULL;
    std::vector<float>* m_recojet_beta0_5_C2=NULL;
    std::vector<float>* m_recojet_beta0_5_C3=NULL;
    std::vector<float>* m_recojet_beta0_5_D2=NULL;
    std::vector<float>* m_recojet_beta0_5_ECorr2_E_theta=NULL;
    std::vector<float>* m_recojet_beta0_5_ECorr3_E_theta=NULL;
    std::vector<float>* m_recojet_beta0_5_N2_E_theta=NULL;
    std::vector<float>* m_recojet_beta0_5_N3_E_theta=NULL;
    std::vector<float>* m_recojet_beta0_5_C2_E_theta=NULL;
    std::vector<float>* m_recojet_beta0_5_C3_E_theta=NULL;
    std::vector<float>* m_recojet_beta0_5_D2_E_theta=NULL;

    std::vector<float>* m_recojet_subjet_E=NULL;
    std::vector<float>* m_recojet_subjet_Px=NULL;
    std::vector<float>* m_recojet_subjet_Py=NULL;
    std::vector<float>* m_recojet_subjet_Pz=NULL;
    std::vector<int>* m_recojet_subjet_NCH=NULL;
    std::vector<int>* m_recojet_subjet_NCH_trackPtMin=NULL;
    std::vector<int>* m_recojet_subjet_NPh=NULL;
    std::vector<int>* m_recojet_subjet_NNH=NULL;
    std::vector<int>* m_recojet_subjet_jetindex=NULL;
    std::vector<float>* m_recojet_subjet_CHFraction=NULL;
    std::vector<float>* m_recojet_subjet_CHFraction_trackPtMin=NULL;
    std::vector<float>* m_recojet_subjet_PhFraction=NULL;
    std::vector<float>* m_recojet_subjet_ElFraction=NULL;
    std::vector<float>* m_recojet_subjet_MuFraction=NULL;
    std::vector<float>* m_recojet_subjet_NHFraction=NULL;

    std::vector<float>* m_recojet_subjet_jetChargeE_kappa_0_25=NULL;
    std::vector<float>* m_recojet_subjet_jetChargeE_kappa_0_50=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargeE_kappa_0_75=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargeE_kappa_1_00=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargeE_kappa_0_10=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargeE_kappa_0_15=NULL;
    std::vector<float>* m_recojet_subjet_jetChargeE_kappa_0_20=NULL;
    std::vector<float>* m_recojet_subjet_jetChargeE_kappa_0_30=NULL;

    std::vector<float>* m_recojet_subjet_jetChargePt_kappa_0_25=NULL;
    std::vector<float>* m_recojet_subjet_jetChargePt_kappa_0_50=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargePt_kappa_0_75=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargePt_kappa_1_00=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargePt_kappa_0_10=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargePt_kappa_0_15=NULL;
    std::vector<float>* m_recojet_subjet_jetChargePt_kappa_0_20=NULL;
    std::vector<float>* m_recojet_subjet_jetChargePt_kappa_0_30=NULL;

    //std::vector<float>* m_recojet_subjet_jetChargePProj_kappa_0_25=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargePProj_kappa_0_50=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargePProj_kappa_0_75=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargePProj_kappa_1_00=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargePProj_kappa_0_10=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargePProj_kappa_0_15=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargePProj_kappa_0_20=NULL;
    //std::vector<float>* m_recojet_subjet_jetChargePProj_kappa_0_30=NULL;

    std::vector<float>* m_recojet_subjet_rfj_j_E=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_Px=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_Py=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_Pz=NULL;
    std::vector<int>* m_recojet_subjet_rfj_j_NCH=NULL;
    std::vector<int>* m_recojet_subjet_rfj_j_NPh=NULL;
    std::vector<int>* m_recojet_subjet_rfj_j_NNH=NULL;
    std::vector<int>* m_recojet_subjet_rfj_j_jetindex=NULL;
    std::vector<int>* m_recojet_subjet_rfj_j_subjetindex=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_CHFraction=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_PhFraction=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_ElFraction=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_MuFraction=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_NHFraction=NULL;
    //std::vector<float>* m_recojet_subjet_rfj_j_jetChargeE_kappa_0_25=NULL;
    //std::vector<float>* m_recojet_subjet_rfj_j_jetChargeE_kappa_0_50=NULL;
    //std::vector<float>* m_recojet_subjet_rfj_j_jetChargeE_kappa_0_20=NULL;
    //std::vector<float>* m_recojet_subjet_rfj_j_jetChargeE_kappa_0_30=NULL;

    //std::vector<float>* m_recojet_subjet_rfj_j_jetChargePt_kappa_0_25=NULL;
    //std::vector<float>* m_recojet_subjet_rfj_j_jetChargePt_kappa_0_50=NULL;
    //std::vector<float>* m_recojet_subjet_rfj_j_jetChargePt_kappa_0_20=NULL;
    //std::vector<float>* m_recojet_subjet_rfj_j_jetChargePt_kappa_0_30=NULL;

    std::vector<float>* m_recojet_subjet_rfj_j_BTag=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_CTag=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_OTag=NULL;
    std::vector<int>* m_recojet_subjet_rfj_j_cat=NULL;

    std::vector<float>* m_recojet_subjet_rfj_j_svtx_r=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_svtx_E=NULL;
    std::vector<float>* m_recojet_subjet_rfj_j_svtx_Mass=NULL;
    std::vector<int>* m_recojet_subjet_rfj_j_svtx_nTrack=NULL;
    std::vector<int>* m_recojet_subjet_rfj_j_svtx_Charge=NULL;
    std::vector<int>* m_recojet_subjet_rfj_j_svtx_jetindex=NULL;

    std::vector<float>* m_reco_tau_E=NULL;
    std::vector<float>* m_reco_tau_Px=NULL;
    std::vector<float>* m_reco_tau_Py=NULL;
    std::vector<float>* m_reco_tau_Pz=NULL;
    std::vector<int>* m_reco_tau_Mult=NULL;
    std::vector<int>* m_reco_tau_NCH=NULL;
    std::vector<int>* m_reco_tau_NPh=NULL;
    std::vector<int>* m_reco_tau_NNH=NULL;
    std::vector<int>* m_reco_tau_Charge=NULL;
    std::vector<float>* m_reco_tau_CHFraction=NULL;
    std::vector<float>* m_reco_tau_PhFraction=NULL;
    std::vector<float>* m_reco_tau_ElFraction=NULL;
    std::vector<float>* m_reco_tau_MuFraction=NULL;
    std::vector<float>* m_reco_tau_NHFraction=NULL;

    //in gen tau we move the actual tau jet (or the tau lepton visible part summed up)
    //the MC tau details should all be listed in the MCLepPhTau collection
    std::vector<float>* m_gen_tau_E=NULL;
    std::vector<float>* m_gen_tau_Px=NULL;
    std::vector<float>* m_gen_tau_Py=NULL;
    std::vector<float>* m_gen_tau_Pz=NULL;
    std::vector<int>* m_gen_tau_Mult=NULL;
    std::vector<int>* m_gen_tau_NCH=NULL;
    std::vector<int>* m_gen_tau_NPh=NULL;
    std::vector<int>* m_gen_tau_NNH=NULL;
    std::vector<float>* m_gen_tau_CHFraction=NULL;
    std::vector<float>* m_gen_tau_PhFraction=NULL;
    std::vector<float>* m_gen_tau_ElFraction=NULL;
    std::vector<float>* m_gen_tau_MuFraction=NULL;
    std::vector<float>* m_gen_tau_NHFraction=NULL;
    std::vector<int>* m_gen_tau_MotherPDGID=NULL;
    std::vector<int>* m_gen_tau_Charge=NULL;//sum of charges of MC tau jet
    std::vector<int>* m_gen_tau_MCCharge=NULL;//actual charge of the tau

    float m_totPFO_E=0;
    float m_totPFO_Px=0;
    float m_totPFO_Py=0;
    float m_totPFO_Pz=0;
    int m_totPFO_Mult=0;
    int m_totPFO_NCH=0;
    int m_totPFO_NPh=0;
    int m_totPFO_NNH=0;
    float m_totPFO_CHFraction=0;
    float m_totPFO_PhFraction=0;
    float m_totPFO_ElFraction=0;
    float m_totPFO_MuFraction=0;
    float m_totPFO_NHFraction=0;


    float m_true_E=0;
    float m_true_Px=0;
    float m_true_Py=0;
    float m_true_Pz=0;
    int m_true_Mult=0;
    int m_true_NCH=0;
    int m_true_NPh=0;
    int m_true_NNH=0;
    float m_true_CHFraction=0;
    float m_true_PhFraction=0;
    float m_true_ElFraction=0;
    float m_true_MuFraction=0;
    float m_true_NHFraction=0;

    float m_true_inv_E=0;
    float m_true_inv_Px=0;
    float m_true_inv_Py=0;
    float m_true_inv_Pz=0;
    int m_true_inv_Mult=0;

    float m_gen_y21_max=0;
    float m_gen_y32_max=0;
    float m_gen_y43_max=0;
    
    float m_gen_y21=0;
    float m_gen_y32=0;
    float m_gen_y43=0;
    
    float m_reco_y21= 0;
    float m_reco_y32= 0;
    float m_reco_y43= 0;

    float m_reco_y21_max= 0;
    float m_reco_y32_max= 0;
    float m_reco_y43_max= 0;

    void fillStableDaughterSet(MCParticle*, std::set<MCParticle*> &);
     
    fastjet::contrib::ValenciaPlugin* vlcpl=NULL;
    fastjet::JetDefinition* _jetAlgoType = NULL;
    fastjet::contrib::ExclusiveJetAxes* _vlcAxes=NULL;
    fastjet::contrib::MeasureDefinition* _unnormalizedMeasure_ptR=NULL;
    fastjet::contrib::MeasureDefinition* _unnormalizedMeasure_Lorentz=NULL;

    fastjet::contrib:: EnergyCorrelator* _energycorr2_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr2_Etheta_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr3_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr3_Etheta_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr4_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr4_Etheta_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelatorC2* _energycorrC2_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelatorC2* _energycorrC2_Etheta_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelatorD2* _energycorrD2_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelatorD2* _energycorrD2_Etheta_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelatorN2* _energycorrN2_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelatorN2* _energycorrN2_Etheta_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelatorN3* _energycorrN3_beta0_5 = NULL;
    fastjet::contrib:: EnergyCorrelatorN3* _energycorrN3_Etheta_beta0_5 = NULL;

    fastjet::contrib:: EnergyCorrelator* _energycorr2_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr2_Etheta_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr3_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr3_Etheta_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr4_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr4_Etheta_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelatorC2* _energycorrC2_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelatorC2* _energycorrC2_Etheta_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelatorD2* _energycorrD2_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelatorD2* _energycorrD2_Etheta_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelatorN2* _energycorrN2_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelatorN2* _energycorrN2_Etheta_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelatorN3* _energycorrN3_beta1 = NULL;
    fastjet::contrib:: EnergyCorrelatorN3* _energycorrN3_Etheta_beta1 = NULL;

    fastjet::contrib:: EnergyCorrelator* _energycorr2_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr2_Etheta_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr3_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr3_Etheta_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr4_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelator* _energycorr4_Etheta_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelatorC2* _energycorrC2_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelatorC2* _energycorrC2_Etheta_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelatorD2* _energycorrD2_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelatorD2* _energycorrD2_Etheta_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelatorN2* _energycorrN2_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelatorN2* _energycorrN2_Etheta_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelatorN3* _energycorrN3_beta2 = NULL;
    fastjet::contrib:: EnergyCorrelatorN3* _energycorrN3_Etheta_beta2 = NULL;


    fastjet::contrib::Nsubjettiness* _nSubJettiness1_ptR = NULL;
    fastjet::contrib::Nsubjettiness* _nSubJettiness1_lorentz = NULL;

    fastjet::contrib::Nsubjettiness* _nSubJettiness2_ptR = NULL;
    fastjet::contrib::Nsubjettiness* _nSubJettiness2_lorentz = NULL;
    
    fastjet::contrib::Nsubjettiness* _nSubJettiness3_ptR = NULL;
    fastjet::contrib::Nsubjettiness* _nSubJettiness3_lorentz = NULL;




} ;

#endif



