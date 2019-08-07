#include "HZAnalyzer.h"
#include <EVENT/LCCollection.h>
#include <EVENT/Vertex.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include "TLorentzVector.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include <fastjet/contrib/MeasureDefinition.hh>
#include <fastjet/contrib/MeasureDefinition.hh>

#include "UTIL/PIDHandler.h"

using namespace lcio ;
using namespace marlin ;

//using namespace std;
//using namespace fastjet;
//using namespace fastjet::contrib;


HZAnalyzer aHZAnalyzer;

//HZAnalyzer::HZAnalyzer() : Processor("HZAnalyzer"),_fju(new FastJetUtil()) 
HZAnalyzer::HZAnalyzer() : Processor("HZAnalyzer") 
{


    // modify processor description
    _description = "HZAnalyzer calculates properties of calorimeter showers" ;
   
    registerProcessorParameter( "OutputRootFileName",
                                "ROOT File name to collect plots",
                                m_rootFileName,
                                std::string("HZAnalyzer.root"));

    registerInputCollection( LCIO::MCPARTICLE,
			     "MCParticleCollectionName",
			     "Name of the MCParticle input collection",
			     m_inputMCParticleCollection,
			     std::string("MCPhysicsParticles"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "MCForJetParticleCollectionName",
			     "Name of the MCParticle input collection",
			     m_inputMCJetParticleCollection,
			     std::string("MCParticlePandoraPFOsForJets"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoParticleCollectionName",
			     "Name of the ReconstructedParticle input collection",
			     m_inputRECOParticleCollection,
			     std::string("TightSelectedPandoraPFOs"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetParticleCollectionName",
			     "Name of the ReconstructedParticle jet input collection",
			     m_inputRECOJetParticleCollection,
			     std::string("PandoraPFOsInJets"));
 
   registerInputCollection( LCIO::MCPARTICLE,
			    "MCIsoLepPhParticleCollection",
			    "Name of the isolated gen particle collection",
			    m_genIsoPartColName,
			    std::string("MCIsoLepPhParticles"));

   registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "RecoParticleIsoLepPhCollectionName",
			    "Name of the isolated reco particle collection",
			    m_recoIsoPartColName,
			    std::string("PandoraPFOsIsoLepPhTau"));

   registerInputCollection( LCIO::MCPARTICLE,
			    "MCTrueLepPhParticleCollection",
			    "Name of the true MC lepton and ISR photon collection",
			    m_outputMCTrueLepPhParticleCollection,
			    std::string("MCTrueLepPhParticles"));
			    
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "TauJetCollection", 
			     "Name of the tauJet collection"  ,
			     m_inputTauCollection,
			     std::string("TauRec_PFO")
			     );

   registerProcessorParameter(
			      "saveMEInfo",
			      "safe Higgs ME Info",
			      m_saveMEInfo,
			      bool(true)
			      );

   registerProcessorParameter(
			       "IsoAngle", 
			       "isolation angle",
			       m_angleIso,
			       float(10.1)
			      );
   
   registerProcessorParameter(
			      "genTrueLepPhEMin", 
			      "minE to save true boson leptons and ISR photons",
			      m_genTrueLepPhEMin,
			      float(8.1)
			      );
   registerProcessorParameter(
			      "R", 
			      "jet radius",
			      m_R,
			      float(0.7)
			      );
   registerProcessorParameter(
			      "beta", 
			      "jet beta parameter",
			      m_beta,
			      float(1.0)
			      );

   registerProcessorParameter(
			      "gamma", 
			      "jet gamma parameter",
			      m_gamma,
			      float(1.0)
			      );
   registerProcessorParameter(
			      "TrackPtMin", 
			      "min track pt",
			      m_cutPtTrackMin,
			      float(-1.0)
			      );
   //registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
   //			    "RefinedJet4jCollectionName", 
   //			    "Name of the lcfiplus refined jet collection",
   //			    m_recorefinedjet_4jets_ColName,
   //			    std::string("RefinedVertexJets_4Jets")
   //			    );
   //registerInputCollection( LCIO::VERTEX,
   //			    "SVtxRfJ4jCollectionName", 
   //			    "Name of the lcfiplus secondary vertex collection belonging to refined jet collection",
   //			    m_vtx_rfj_4jets_ColName,
   //			    std::string("RefinedVertices_4Jets")
   //			    );

   registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "RefinedJet0CollectionName", 
			    "Name of the lcfiplus refined jet0 collection",
			    m_recorefinedjet0ColName,
			    std::string("RefinedVertexJets_Jet0")
			    );
   registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "RefinedJet1CollectionName", 
			    "Name of the lcfiplus refined jet1 collection",
			    m_recorefinedjet1ColName,
			    std::string("RefinedVertexJets_Jet1")
			    );


   registerInputCollection( LCIO::VERTEX,
			    "SVtxRfJ0CollectionName", 
			    "Name of the lcfiplus secondary vertex collection belonging to refined jet0 collection",
			    m_vtx_rfj0ColName,
			    std::string("RefinedVertices_Jet0")
			    );
   registerInputCollection( LCIO::VERTEX,
			    "SVtxRfJ1CollectionName", 
			    "Name of the lcfiplus secondary vertex collection belonging to refined jet1 collection",
			    m_vtx_rfj1ColName,
			    std::string("RefinedVertices_Jet1")
			    );
}  


void HZAnalyzer::init() {

  std::cout<<"R/beta/gamma parameter HZAnalyzer "<<m_R<<"/"<<m_beta<<"/"<<m_gamma<<std::endl;

  m_rootFile = new TFile(m_rootFileName.c_str(),"RECREATE");

  m_outputTree = new TTree("showerData","showerData");

  m_eventCount=0;
  
  // Print the initial parameters
  printParameters() ;
  
  // Reset counters
  m_runNumber = 0 ;
  m_eventNumber = 0 ;

  m_totPFO_E=0;
  m_totPFO_Px=0;
  m_totPFO_Py=0;
  m_totPFO_Pz=0;
  m_totPFO_Mult=0;
  m_totPFO_NCH=0;
  m_totPFO_NPh=0;
  m_totPFO_NNH=0;
  m_totPFO_CHFraction=0;
  m_totPFO_PhFraction=0;
  m_totPFO_ElFraction=0;
  m_totPFO_MuFraction=0;
  m_totPFO_NHFraction=0;
  
  m_true_E=0;
  m_true_Px=0;
  m_true_Py=0;
  m_true_Pz=0;
  m_true_Mult=0;
  m_true_NCH=0;
  m_true_NPh=0;
  m_true_NNH=0;
  m_true_CHFraction=0;
  m_true_PhFraction=0;
  m_true_ElFraction=0;
  m_true_MuFraction=0;
  m_true_NHFraction=0;
  
  m_true_inv_E=0;
  m_true_inv_Px=0;
  m_true_inv_Py=0;
  m_true_inv_Pz=0;
  m_true_inv_Mult=0;


  m_gen_y21_max=0;
  m_gen_y32_max=0;
  m_gen_y43_max=0;

  m_gen_y21=0;
  m_gen_y32=0;
  m_gen_y43=0;

  m_reco_y21= 0;
  m_reco_y32= 0;
  m_reco_y43= 0;
  m_reco_y21_max= 0;
  m_reco_y32_max= 0;
  m_reco_y43_max= 0;
  
  m_trueME_E= new std::vector<float>();
  m_trueME_Px= new std::vector<float>();
  m_trueME_Py= new std::vector<float>();
  m_trueME_Pz= new std::vector<float>();
  m_trueME_PDGID= new std::vector<int>();
  //m_trueME_index= new std::vector<int>();
  
  m_genTrueLepPh_E= new std::vector<float>();
  m_genTrueLepPh_Px= new std::vector<float>();
  m_genTrueLepPh_Py= new std::vector<float>();
  m_genTrueLepPh_Pz= new std::vector<float>();
  m_genTrueLepPh_PDGID= new std::vector<int>();
  //actual ancestor, i.e. Higgs,W,tau
  m_genTrueLepPh_ANC_PDGID= new std::vector<int>();
  m_genTrueLepPh_relIso= new std::vector<float>();
  m_genTrueLepPh_relIsoCH= new std::vector<float>();
  m_genTrueLepPh_relIsoPh= new std::vector<float>();
  m_genTrueLepPh_relIsoNH= new std::vector<float>();
  m_genTrueLepPh_relIsoEl= new std::vector<float>();
  m_genTrueLepPh_relIsoMu= new std::vector<float>();
  
  //m_genjet_jetChargeE_kappa_0_25= new std::vector<float>();
  //m_genjet_jetChargeE_kappa_0_50= new std::vector<float>();
  //m_genjet_jetChargeE_kappa_0_75= new std::vector<float>();
  //m_genjet_jetChargeE_kappa_1_00= new std::vector<float>();
  //m_genjet_jetChargeE_kappa_0_10= new std::vector<float>();
  //m_genjet_jetChargeE_kappa_0_15= new std::vector<float>();
  //m_genjet_jetChargeE_kappa_0_20= new std::vector<float>();
  //m_genjet_jetChargeE_kappa_0_30= new std::vector<float>();
  
  //m_genjet_jetChargePt_kappa_0_25= new std::vector<float>();
  //m_genjet_jetChargePt_kappa_0_50= new std::vector<float>();
  //m_genjet_jetChargePt_kappa_0_75= new std::vector<float>();
  //m_genjet_jetChargePt_kappa_1_00= new std::vector<float>();
  //m_genjet_jetChargePt_kappa_0_10= new std::vector<float>();
  //m_genjet_jetChargePt_kappa_0_15= new std::vector<float>();
  //m_genjet_jetChargePt_kappa_0_20= new std::vector<float>();
  //m_genjet_jetChargePt_kappa_0_30= new std::vector<float>();

  //m_genjet_jetChargePProj_kappa_0_25= new std::vector<float>();
  //m_genjet_jetChargePProj_kappa_0_50= new std::vector<float>();
  //m_genjet_jetChargePProj_kappa_0_75= new std::vector<float>();
  //m_genjet_jetChargePProj_kappa_1_00= new std::vector<float>();
  //m_genjet_jetChargePProj_kappa_0_10= new std::vector<float>();
  //m_genjet_jetChargePProj_kappa_0_15= new std::vector<float>();
  //m_genjet_jetChargePProj_kappa_0_20= new std::vector<float>();
  //m_genjet_jetChargePProj_kappa_0_30= new std::vector<float>();

  
  m_genjet_E= new std::vector<float>();
  m_genjet_Px= new std::vector<float>();
  m_genjet_Py= new std::vector<float>();
  m_genjet_Pz= new std::vector<float>();
  m_genjet_Mult= new std::vector<int>();
  m_genjet_NCH= new std::vector<int>();
  m_genjet_NCH_trackPtMin= new std::vector<int>();
  m_genjet_NPh= new std::vector<int>();
  m_genjet_NNH= new std::vector<int>();
  m_genjet_dij_21= new std::vector<float>();
  m_genjet_dij_32= new std::vector<float>();
  m_genjet_dij_43= new std::vector<float>();
  m_genjet_dij_21_max= new std::vector<float>();
  m_genjet_dij_32_max= new std::vector<float>();
  m_genjet_dij_43_max= new std::vector<float>();

  m_genjet_CHFraction= new std::vector<float>();
  m_genjet_CHFraction_trackPtMin= new std::vector<float>();
  m_genjet_PhFraction= new std::vector<float>();
  m_genjet_ElFraction= new std::vector<float>();
  m_genjet_MuFraction= new std::vector<float>();
  m_genjet_NHFraction= new std::vector<float>();
  m_genjet_nsubjettiness1 = new std::vector<float>();
  m_genjet_nsubjettiness2 = new std::vector<float>();
  m_genjet_nsubjettiness3 = new std::vector<float>();
  m_genjet_nsubjettiness1_lrz = new std::vector<float>();
  m_genjet_nsubjettiness2_lrz = new std::vector<float>();
  m_genjet_nsubjettiness3_lrz = new std::vector<float>();
  m_genjet_beta1_ECorr2 = new std::vector<float>();
  m_genjet_beta1_ECorr3 = new std::vector<float>();
  m_genjet_beta1_N2 = new std::vector<float>();
  m_genjet_beta1_N3 = new std::vector<float>();
  m_genjet_beta1_C2 = new std::vector<float>();
  m_genjet_beta1_C3 = new std::vector<float>();
  m_genjet_beta1_D2 = new std::vector<float>();
  m_genjet_beta1_ECorr2_E_theta = new std::vector<float>();
  m_genjet_beta1_ECorr3_E_theta = new std::vector<float>();
  m_genjet_beta1_N2_E_theta = new std::vector<float>();
  m_genjet_beta1_N3_E_theta = new std::vector<float>();
  m_genjet_beta1_C2_E_theta = new std::vector<float>();
  m_genjet_beta1_C3_E_theta = new std::vector<float>();
  m_genjet_beta1_D2_E_theta = new std::vector<float>();

  m_genjet_beta2_ECorr2 = new std::vector<float>();
  m_genjet_beta2_ECorr3 = new std::vector<float>();
  m_genjet_beta2_N2 = new std::vector<float>();
  m_genjet_beta2_N3 = new std::vector<float>();
  m_genjet_beta2_C2 = new std::vector<float>();
  m_genjet_beta2_C3 = new std::vector<float>();
  m_genjet_beta2_D2 = new std::vector<float>();
  m_genjet_beta2_ECorr2_E_theta = new std::vector<float>();
  m_genjet_beta2_ECorr3_E_theta = new std::vector<float>();
  m_genjet_beta2_N2_E_theta = new std::vector<float>();
  m_genjet_beta2_N3_E_theta = new std::vector<float>();
  m_genjet_beta2_C2_E_theta = new std::vector<float>();
  m_genjet_beta2_C3_E_theta = new std::vector<float>();
  m_genjet_beta2_D2_E_theta = new std::vector<float>();

  m_genjet_beta0_5_ECorr2 = new std::vector<float>();
  m_genjet_beta0_5_ECorr3 = new std::vector<float>();
  m_genjet_beta0_5_N2 = new std::vector<float>();
  m_genjet_beta0_5_N3 = new std::vector<float>();
  m_genjet_beta0_5_C2 = new std::vector<float>();
  m_genjet_beta0_5_C3 = new std::vector<float>();
  m_genjet_beta0_5_D2 = new std::vector<float>();
  m_genjet_beta0_5_ECorr2_E_theta = new std::vector<float>();
  m_genjet_beta0_5_ECorr3_E_theta = new std::vector<float>();
  m_genjet_beta0_5_N2_E_theta = new std::vector<float>();
  m_genjet_beta0_5_N3_E_theta = new std::vector<float>();
  m_genjet_beta0_5_C2_E_theta = new std::vector<float>();
  m_genjet_beta0_5_C3_E_theta = new std::vector<float>();
  m_genjet_beta0_5_D2_E_theta = new std::vector<float>();

  m_genjet_subjet_E= new std::vector<float>();
  m_genjet_subjet_Px= new std::vector<float>();
  m_genjet_subjet_Py= new std::vector<float>();
  m_genjet_subjet_Pz= new std::vector<float>();
  m_genjet_subjet_NCH= new std::vector<int>();
  m_genjet_subjet_NCH_trackPtMin= new std::vector<int>();
  m_genjet_subjet_NPh= new std::vector<int>();
  m_genjet_subjet_NNH= new std::vector<int>();
  m_genjet_subjet_jetindex= new std::vector<int>();
  m_genjet_subjet_CHFraction= new std::vector<float>();
  m_genjet_subjet_CHFraction_trackPtMin= new std::vector<float>();
  m_genjet_subjet_PhFraction= new std::vector<float>();
  m_genjet_subjet_ElFraction= new std::vector<float>();
  m_genjet_subjet_MuFraction= new std::vector<float>();
  m_genjet_subjet_NHFraction= new std::vector<float>();

  m_genjet_subjet_jetChargeE_kappa_0_25= new std::vector<float>();
  m_genjet_subjet_jetChargeE_kappa_0_50= new std::vector<float>();
  //m_genjet_subjet_jetChargeE_kappa_0_75= new std::vector<float>();
  //m_genjet_subjet_jetChargeE_kappa_1_00= new std::vector<float>();
  //m_genjet_subjet_jetChargeE_kappa_0_10= new std::vector<float>();
  //m_genjet_subjet_jetChargeE_kappa_0_15= new std::vector<float>();
  m_genjet_subjet_jetChargeE_kappa_0_20= new std::vector<float>();
  m_genjet_subjet_jetChargeE_kappa_0_30= new std::vector<float>();
  
  m_genjet_subjet_jetChargePt_kappa_0_25= new std::vector<float>();
  m_genjet_subjet_jetChargePt_kappa_0_50= new std::vector<float>();
  //m_genjet_subjet_jetChargePt_kappa_0_75= new std::vector<float>();
  //m_genjet_subjet_jetChargePt_kappa_1_00= new std::vector<float>();
  //m_genjet_subjet_jetChargePt_kappa_0_10= new std::vector<float>();
  //m_genjet_subjet_jetChargePt_kappa_0_15= new std::vector<float>();
  m_genjet_subjet_jetChargePt_kappa_0_20= new std::vector<float>();
  m_genjet_subjet_jetChargePt_kappa_0_30= new std::vector<float>();

  //m_genjet_subjet_jetChargePProj_kappa_0_25= new std::vector<float>();
  //m_genjet_subjet_jetChargePProj_kappa_0_50= new std::vector<float>();
  //m_genjet_subjet_jetChargePProj_kappa_0_75= new std::vector<float>();
  //m_genjet_subjet_jetChargePProj_kappa_1_00= new std::vector<float>();
  //m_genjet_subjet_jetChargePProj_kappa_0_10= new std::vector<float>();
  //m_genjet_subjet_jetChargePProj_kappa_0_15= new std::vector<float>();
  //m_genjet_subjet_jetChargePProj_kappa_0_20= new std::vector<float>();
  //m_genjet_subjet_jetChargePProj_kappa_0_30= new std::vector<float>();

  m_isoPartGenDR10_E= new std::vector<float>();
  m_isoPartGenDR10_Px= new std::vector<float>();
  m_isoPartGenDR10_Py= new std::vector<float>();
  m_isoPartGenDR10_Pz= new std::vector<float>();
  m_isoPartGenDR10_PDGID= new std::vector<int>();
  m_isoPartGenDR10_relIso= new std::vector<float>();
  m_isoPartGenDR10_relIsoCH= new std::vector<float>();
  m_isoPartGenDR10_relIsoPh= new std::vector<float>();
  m_isoPartGenDR10_relIsoNH= new std::vector<float>();
  m_isoPartGenDR10_relIsoEl= new std::vector<float>();
  m_isoPartGenDR10_relIsoMu= new std::vector<float>();
  
  m_isoPartRecoDR10_E= new std::vector<float>();
  m_isoPartRecoDR10_Px= new std::vector<float>();
  m_isoPartRecoDR10_Py= new std::vector<float>();
  m_isoPartRecoDR10_Pz= new std::vector<float>();
  m_isoPartRecoDR10_PDGID= new std::vector<int>();
  m_isoPartRecoDR10_relIso= new std::vector<float>();
  m_isoPartRecoDR10_relIsoCH= new std::vector<float>();
  m_isoPartRecoDR10_relIsoPh= new std::vector<float>();
  m_isoPartRecoDR10_relIsoNH= new std::vector<float>();
  m_isoPartRecoDR10_relIsoEl= new std::vector<float>();
  m_isoPartRecoDR10_relIsoMu= new std::vector<float>();
  
  //m_recojet_jetChargeE_kappa_0_25= new std::vector<float>();
  //m_recojet_jetChargeE_kappa_0_50= new std::vector<float>();
  //m_recojet_jetChargeE_kappa_0_75= new std::vector<float>();
  //m_recojet_jetChargeE_kappa_1_00= new std::vector<float>();
  //m_recojet_jetChargeE_kappa_0_10= new std::vector<float>();
  //m_recojet_jetChargeE_kappa_0_15= new std::vector<float>();
  //m_recojet_jetChargeE_kappa_0_20= new std::vector<float>();
  //m_recojet_jetChargeE_kappa_0_30= new std::vector<float>();
  
  //m_recojet_jetChargePt_kappa_0_25= new std::vector<float>();
  //m_recojet_jetChargePt_kappa_0_50= new std::vector<float>();
  //m_recojet_jetChargePt_kappa_0_75= new std::vector<float>();
  //m_recojet_jetChargePt_kappa_1_00= new std::vector<float>();
  //m_recojet_jetChargePt_kappa_0_10= new std::vector<float>();
  //m_recojet_jetChargePt_kappa_0_15= new std::vector<float>();
  //m_recojet_jetChargePt_kappa_0_20= new std::vector<float>();
  //m_recojet_jetChargePt_kappa_0_30= new std::vector<float>();

  //m_recojet_jetChargePProj_kappa_0_25= new std::vector<float>();
  //m_recojet_jetChargePProj_kappa_0_50= new std::vector<float>();
  //m_recojet_jetChargePProj_kappa_0_75= new std::vector<float>();
  //m_recojet_jetChargePProj_kappa_1_00= new std::vector<float>();
  //m_recojet_jetChargePProj_kappa_0_10= new std::vector<float>();
  //m_recojet_jetChargePProj_kappa_0_15= new std::vector<float>();
  //m_recojet_jetChargePProj_kappa_0_20= new std::vector<float>();
  //m_recojet_jetChargePProj_kappa_0_30= new std::vector<float>();

  /*
  m_recojet_rfj_4jets_E= new std::vector<float>();
  m_recojet_rfj_4jets_Px= new std::vector<float>();
  m_recojet_rfj_4jets_Py= new std::vector<float>();
  m_recojet_rfj_4jets_Pz= new std::vector<float>();
  m_recojet_rfj_4jets_Mult= new std::vector<int>();
  m_recojet_rfj_4jets_NCH= new std::vector<int>();
  m_recojet_rfj_4jets_CHFraction= new std::vector<float>();
  m_recojet_rfj_4jets_PhFraction= new std::vector<float>();
  m_recojet_rfj_4jets_ElFraction= new std::vector<float>();
  m_recojet_rfj_4jets_MuFraction= new std::vector<float>();
  m_recojet_rfj_4jets_NHFraction= new std::vector<float>();

  m_recojet_rfj_4jets_BTag = new std::vector<float>();
  m_recojet_rfj_4jets_CTag = new std::vector<float>();
  m_recojet_rfj_4jets_OTag = new std::vector<float>();
  m_recojet_rfj_4jets_cat = new std::vector<int>();
  
  m_recojet_rfj_4jets_svtx_r = new std::vector<float>();
  m_recojet_rfj_4jets_svtx_E = new std::vector<float>();
  m_recojet_rfj_4jets_svtx_Mass = new std::vector<float>();
  m_recojet_rfj_4jets_svtx_nTrack = new std::vector<int>();
  m_recojet_rfj_4jets_svtx_Charge = new std::vector<int>();
  
  m_recojet_fat2j_rfj4_E= new std::vector<float>();
  m_recojet_fat2j_rfj4_Px= new std::vector<float>();
  m_recojet_fat2j_rfj4_Py= new std::vector<float>();
  m_recojet_fat2j_rfj4_Pz= new std::vector<float>();
  m_recojet_fat2j_rfj4_Mult= new std::vector<int>();
  */

  m_recojet_E= new std::vector<float>();
  m_recojet_Px= new std::vector<float>();
  m_recojet_Py= new std::vector<float>();
  m_recojet_Pz= new std::vector<float>();
  m_recojet_Mult= new std::vector<int>();
  m_recojet_NCH= new std::vector<int>();
  m_recojet_NCH_trackPtMin= new std::vector<int>();
  m_recojet_NPh= new std::vector<int>();
  m_recojet_NNH= new std::vector<int>();
  m_recojet_dij_21= new std::vector<float>();
  m_recojet_dij_32= new std::vector<float>();
  m_recojet_dij_43= new std::vector<float>();
  m_recojet_dij_21_max= new std::vector<float>();
  m_recojet_dij_32_max= new std::vector<float>();
  m_recojet_dij_43_max= new std::vector<float>();

  m_recojet_CHFraction= new std::vector<float>();
  m_recojet_CHFraction_trackPtMin= new std::vector<float>();
  m_recojet_PhFraction= new std::vector<float>();
  m_recojet_ElFraction= new std::vector<float>();
  m_recojet_MuFraction= new std::vector<float>();
  m_recojet_NHFraction= new std::vector<float>();
  m_recojet_nsubjettiness1 = new std::vector<float>();
  m_recojet_nsubjettiness2 = new std::vector<float>();
  m_recojet_nsubjettiness3 = new std::vector<float>();
  m_recojet_nsubjettiness1_lrz = new std::vector<float>();
  m_recojet_nsubjettiness2_lrz = new std::vector<float>();
  m_recojet_nsubjettiness3_lrz = new std::vector<float>();

  m_recojet_beta1_ECorr2 = new std::vector<float>();
  m_recojet_beta1_ECorr3 = new std::vector<float>();
  m_recojet_beta1_N2 = new std::vector<float>();
  m_recojet_beta1_N3 = new std::vector<float>();
  m_recojet_beta1_C2 = new std::vector<float>();
  m_recojet_beta1_C3 = new std::vector<float>();
  m_recojet_beta1_D2 = new std::vector<float>();
  m_recojet_beta1_ECorr2_E_theta = new std::vector<float>();
  m_recojet_beta1_ECorr3_E_theta = new std::vector<float>();
  m_recojet_beta1_N2_E_theta = new std::vector<float>();
  m_recojet_beta1_N3_E_theta = new std::vector<float>();
  m_recojet_beta1_C2_E_theta = new std::vector<float>();
  m_recojet_beta1_C3_E_theta = new std::vector<float>();
  m_recojet_beta1_D2_E_theta = new std::vector<float>();

  m_recojet_beta2_ECorr2 = new std::vector<float>();
  m_recojet_beta2_ECorr3 = new std::vector<float>();
  m_recojet_beta2_N2 = new std::vector<float>();
  m_recojet_beta2_N3 = new std::vector<float>();
  m_recojet_beta2_C2 = new std::vector<float>();
  m_recojet_beta2_C3 = new std::vector<float>();
  m_recojet_beta2_D2 = new std::vector<float>();
  m_recojet_beta2_ECorr2_E_theta = new std::vector<float>();
  m_recojet_beta2_ECorr3_E_theta = new std::vector<float>();
  m_recojet_beta2_N2_E_theta = new std::vector<float>();
  m_recojet_beta2_N3_E_theta = new std::vector<float>();
  m_recojet_beta2_C2_E_theta = new std::vector<float>();
  m_recojet_beta2_C3_E_theta = new std::vector<float>();
  m_recojet_beta2_D2_E_theta = new std::vector<float>();

  m_recojet_beta0_5_ECorr2 = new std::vector<float>();
  m_recojet_beta0_5_ECorr3 = new std::vector<float>();
  m_recojet_beta0_5_N2 = new std::vector<float>();
  m_recojet_beta0_5_N3 = new std::vector<float>();
  m_recojet_beta0_5_C2 = new std::vector<float>();
  m_recojet_beta0_5_C3 = new std::vector<float>();
  m_recojet_beta0_5_D2 = new std::vector<float>();
  m_recojet_beta0_5_ECorr2_E_theta = new std::vector<float>();
  m_recojet_beta0_5_ECorr3_E_theta = new std::vector<float>();
  m_recojet_beta0_5_N2_E_theta = new std::vector<float>();
  m_recojet_beta0_5_N3_E_theta = new std::vector<float>();
  m_recojet_beta0_5_C2_E_theta = new std::vector<float>();
  m_recojet_beta0_5_C3_E_theta = new std::vector<float>();
  m_recojet_beta0_5_D2_E_theta = new std::vector<float>();

  m_recojet_subjet_E= new std::vector<float>();
  m_recojet_subjet_Px= new std::vector<float>();
  m_recojet_subjet_Py= new std::vector<float>();
  m_recojet_subjet_Pz= new std::vector<float>();
  m_recojet_subjet_NCH= new std::vector<int>();
  m_recojet_subjet_NCH_trackPtMin= new std::vector<int>();
  m_recojet_subjet_NPh= new std::vector<int>();
  m_recojet_subjet_NNH= new std::vector<int>();
  m_recojet_subjet_jetindex= new std::vector<int>();
  m_recojet_subjet_CHFraction= new std::vector<float>();
  m_recojet_subjet_CHFraction_trackPtMin= new std::vector<float>();
  m_recojet_subjet_PhFraction= new std::vector<float>();
  m_recojet_subjet_ElFraction= new std::vector<float>();
  m_recojet_subjet_MuFraction= new std::vector<float>();
  m_recojet_subjet_NHFraction= new std::vector<float>();
  
  m_recojet_subjet_jetChargeE_kappa_0_25= new std::vector<float>();
  m_recojet_subjet_jetChargeE_kappa_0_50= new std::vector<float>();
  //m_recojet_subjet_jetChargeE_kappa_0_75= new std::vector<float>();
  //m_recojet_subjet_jetChargeE_kappa_1_00= new std::vector<float>();
  //m_recojet_subjet_jetChargeE_kappa_0_10= new std::vector<float>();
  //m_recojet_subjet_jetChargeE_kappa_0_15= new std::vector<float>();
  m_recojet_subjet_jetChargeE_kappa_0_20= new std::vector<float>();
  m_recojet_subjet_jetChargeE_kappa_0_30= new std::vector<float>();
  
  m_recojet_subjet_jetChargePt_kappa_0_25= new std::vector<float>();
  m_recojet_subjet_jetChargePt_kappa_0_50= new std::vector<float>();
  //m_recojet_subjet_jetChargePt_kappa_0_75= new std::vector<float>();
  //m_recojet_subjet_jetChargePt_kappa_1_00= new std::vector<float>();
  //m_recojet_subjet_jetChargePt_kappa_0_10= new std::vector<float>();
  //m_recojet_subjet_jetChargePt_kappa_0_15= new std::vector<float>();
  m_recojet_subjet_jetChargePt_kappa_0_20= new std::vector<float>();
  m_recojet_subjet_jetChargePt_kappa_0_30= new std::vector<float>();

  //m_recojet_subjet_jetChargePProj_kappa_0_25= new std::vector<float>();
  //m_recojet_subjet_jetChargePProj_kappa_0_50= new std::vector<float>();
  //m_recojet_subjet_jetChargePProj_kappa_0_75= new std::vector<float>();
  //m_recojet_subjet_jetChargePProj_kappa_1_00= new std::vector<float>();
  //m_recojet_subjet_jetChargePProj_kappa_0_10= new std::vector<float>();
  //m_recojet_subjet_jetChargePProj_kappa_0_15= new std::vector<float>();
  //m_recojet_subjet_jetChargePProj_kappa_0_20= new std::vector<float>();
  //m_recojet_subjet_jetChargePProj_kappa_0_30= new std::vector<float>();
  
  m_recojet_subjet_rfj_j_E = new std::vector<float>();
  m_recojet_subjet_rfj_j_Px = new std::vector<float>();
  m_recojet_subjet_rfj_j_Py = new std::vector<float>();
  m_recojet_subjet_rfj_j_Pz = new std::vector<float>();
  m_recojet_subjet_rfj_j_NCH = new std::vector<int>();
  m_recojet_subjet_rfj_j_NPh = new std::vector<int>();
  m_recojet_subjet_rfj_j_NNH = new std::vector<int>();
  m_recojet_subjet_rfj_j_jetindex = new std::vector<int>();
  m_recojet_subjet_rfj_j_subjetindex = new std::vector<int>();
  m_recojet_subjet_rfj_j_CHFraction = new std::vector<float>();
  m_recojet_subjet_rfj_j_PhFraction = new std::vector<float>();
  m_recojet_subjet_rfj_j_ElFraction = new std::vector<float>();
  m_recojet_subjet_rfj_j_MuFraction = new std::vector<float>();
  m_recojet_subjet_rfj_j_NHFraction = new std::vector<float>();
  //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_25 = new std::vector<float>();
  //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_50 = new std::vector<float>();
  //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_20 = new std::vector<float>();
  //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_30 = new std::vector<float>();
  
  //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_25 = new std::vector<float>();
  //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_50 = new std::vector<float>();
  //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_20 = new std::vector<float>();
  //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_30 = new std::vector<float>();
  
  m_recojet_subjet_rfj_j_BTag = new std::vector<float>();
  m_recojet_subjet_rfj_j_CTag = new std::vector<float>();
  m_recojet_subjet_rfj_j_OTag = new std::vector<float>();
  m_recojet_subjet_rfj_j_cat = new std::vector<int>();
  
  m_recojet_subjet_rfj_j_svtx_r = new std::vector<float>();
  m_recojet_subjet_rfj_j_svtx_E = new std::vector<float>();
  m_recojet_subjet_rfj_j_svtx_Mass = new std::vector<float>();
  m_recojet_subjet_rfj_j_svtx_nTrack = new std::vector<int>();
  m_recojet_subjet_rfj_j_svtx_Charge = new std::vector<int>();
  m_recojet_subjet_rfj_j_svtx_jetindex = new std::vector<int>();

  m_reco_tau_E= new std::vector<float>();
  m_reco_tau_Px= new std::vector<float>();
  m_reco_tau_Py= new std::vector<float>();
  m_reco_tau_Pz= new std::vector<float>();
  m_reco_tau_Mult= new std::vector<int>();
  m_reco_tau_Charge= new std::vector<int>();
  m_reco_tau_NCH= new std::vector<int>();
  m_reco_tau_NPh= new std::vector<int>();
  m_reco_tau_NNH= new std::vector<int>();
  m_reco_tau_CHFraction= new std::vector<float>();
  m_reco_tau_PhFraction= new std::vector<float>();
  m_reco_tau_ElFraction= new std::vector<float>();
  m_reco_tau_MuFraction= new std::vector<float>();
  m_reco_tau_NHFraction= new std::vector<float>();

  m_gen_tau_E= new std::vector<float>();
  m_gen_tau_Px= new std::vector<float>();
  m_gen_tau_Py= new std::vector<float>();
  m_gen_tau_Pz= new std::vector<float>();
  m_gen_tau_Mult= new std::vector<int>();
  m_gen_tau_Charge= new std::vector<int>();
  m_gen_tau_MCCharge= new std::vector<int>();
  m_gen_tau_NCH= new std::vector<int>();
  m_gen_tau_NPh= new std::vector<int>();
  m_gen_tau_NNH= new std::vector<int>();
  m_gen_tau_MotherPDGID= new std::vector<int>();
  m_gen_tau_CHFraction= new std::vector<float>();
  m_gen_tau_PhFraction= new std::vector<float>();
  m_gen_tau_ElFraction= new std::vector<float>();
  m_gen_tau_MuFraction= new std::vector<float>();
  m_gen_tau_NHFraction= new std::vector<float>();
  
  std::cout<<"after creating the vectors "<<std::endl;

  m_trueME_E->clear();
  m_trueME_Px->clear();
  m_trueME_Py->clear();
  m_trueME_Pz->clear();
  m_trueME_PDGID->clear();
  //m_trueME_index->clear();
  
  m_genTrueLepPh_E->clear();
  m_genTrueLepPh_Px->clear();
  m_genTrueLepPh_Py->clear();
  m_genTrueLepPh_Pz->clear();
  m_genTrueLepPh_PDGID->clear();
  //actual ancestor, i.e. Higgs,W,tau
  m_genTrueLepPh_ANC_PDGID->clear();
  m_genTrueLepPh_relIso->clear();
  m_genTrueLepPh_relIsoCH->clear();
  m_genTrueLepPh_relIsoPh->clear();
  m_genTrueLepPh_relIsoNH->clear();
  m_genTrueLepPh_relIsoEl->clear();
  m_genTrueLepPh_relIsoMu->clear();
  
  //m_genjet_jetChargeE_kappa_0_25->clear();
  //m_genjet_jetChargeE_kappa_0_50->clear();
  //m_genjet_jetChargeE_kappa_0_75->clear();
  //m_genjet_jetChargeE_kappa_1_00->clear();
  //m_genjet_jetChargeE_kappa_0_10->clear();
  //m_genjet_jetChargeE_kappa_0_15->clear();
  //m_genjet_jetChargeE_kappa_0_20->clear();
  //m_genjet_jetChargeE_kappa_0_30->clear();
  
  //m_genjet_jetChargePt_kappa_0_25->clear();
  //m_genjet_jetChargePt_kappa_0_50->clear();
  //m_genjet_jetChargePt_kappa_0_75->clear();
  //m_genjet_jetChargePt_kappa_1_00->clear();
  //m_genjet_jetChargePt_kappa_0_10->clear();
  //m_genjet_jetChargePt_kappa_0_15->clear();
  //m_genjet_jetChargePt_kappa_0_20->clear();
  //m_genjet_jetChargePt_kappa_0_30->clear();


  //m_genjet_jetChargePProj_kappa_0_25->clear();
  //m_genjet_jetChargePProj_kappa_0_50->clear();
  //m_genjet_jetChargePProj_kappa_0_75->clear();
  //m_genjet_jetChargePProj_kappa_1_00->clear();
  //m_genjet_jetChargePProj_kappa_0_10->clear();
  //m_genjet_jetChargePProj_kappa_0_15->clear();
  //m_genjet_jetChargePProj_kappa_0_20->clear();
  //m_genjet_jetChargePProj_kappa_0_30->clear();
  
  m_genjet_E->clear();
  m_genjet_Px->clear();
  m_genjet_Py->clear();
  m_genjet_Pz->clear();
  m_genjet_Mult->clear();
  m_genjet_NCH->clear();
  m_genjet_NCH_trackPtMin->clear();
  m_genjet_NPh->clear();
  m_genjet_NNH->clear();
  m_genjet_dij_21->clear();
  m_genjet_dij_32->clear();
  m_genjet_dij_43->clear();
  m_genjet_dij_21_max->clear();
  m_genjet_dij_32_max->clear();
  m_genjet_dij_43_max->clear();
  m_genjet_CHFraction->clear();
  m_genjet_CHFraction_trackPtMin->clear();
  m_genjet_PhFraction->clear();
  m_genjet_ElFraction->clear();
  m_genjet_MuFraction->clear();
  m_genjet_NHFraction->clear();
  m_genjet_subjet_CHFraction->clear();
  m_genjet_subjet_CHFraction_trackPtMin->clear();
  m_genjet_subjet_PhFraction->clear();
  m_genjet_subjet_ElFraction->clear();
  m_genjet_subjet_MuFraction->clear();
  m_genjet_subjet_NHFraction->clear();
  m_genjet_nsubjettiness1->clear();
  m_genjet_nsubjettiness2->clear();
  m_genjet_nsubjettiness3->clear();
  m_genjet_nsubjettiness1_lrz->clear();
  m_genjet_nsubjettiness2_lrz->clear();
  m_genjet_nsubjettiness3_lrz->clear();
  m_genjet_beta1_ECorr2->clear();
  m_genjet_beta1_ECorr3->clear();
  m_genjet_beta1_N2->clear();
  m_genjet_beta1_N3->clear();
  m_genjet_beta1_C2->clear();
  m_genjet_beta1_C3->clear();
  m_genjet_beta1_D2->clear();
  m_genjet_beta1_ECorr2_E_theta->clear();
  m_genjet_beta1_ECorr3_E_theta->clear();
  m_genjet_beta1_N2_E_theta->clear();
  m_genjet_beta1_N3_E_theta->clear();
  m_genjet_beta1_C2_E_theta->clear();
  m_genjet_beta1_C3_E_theta->clear();
  m_genjet_beta1_D2_E_theta->clear();

  m_genjet_beta2_ECorr2->clear();
  m_genjet_beta2_ECorr3->clear();
  m_genjet_beta2_N2->clear();
  m_genjet_beta2_N3->clear();
  m_genjet_beta2_C2->clear();
  m_genjet_beta2_C3->clear();
  m_genjet_beta2_D2->clear();
  m_genjet_beta2_ECorr2_E_theta->clear();
  m_genjet_beta2_ECorr3_E_theta->clear();
  m_genjet_beta2_N2_E_theta->clear();
  m_genjet_beta2_N3_E_theta->clear();
  m_genjet_beta2_C2_E_theta->clear();
  m_genjet_beta2_C3_E_theta->clear();
  m_genjet_beta2_D2_E_theta->clear();

  m_genjet_beta0_5_ECorr2->clear();
  m_genjet_beta0_5_ECorr3->clear();
  m_genjet_beta0_5_N2->clear();
  m_genjet_beta0_5_N3->clear();
  m_genjet_beta0_5_C2->clear();
  m_genjet_beta0_5_C3->clear();
  m_genjet_beta0_5_D2->clear();
  m_genjet_beta0_5_ECorr2_E_theta->clear();
  m_genjet_beta0_5_ECorr3_E_theta->clear();
  m_genjet_beta0_5_N2_E_theta->clear();
  m_genjet_beta0_5_N3_E_theta->clear();
  m_genjet_beta0_5_C2_E_theta->clear();
  m_genjet_beta0_5_C3_E_theta->clear();
  m_genjet_beta0_5_D2_E_theta->clear();

  m_genjet_subjet_E->clear();
  m_genjet_subjet_Px->clear();
  m_genjet_subjet_Py->clear();
  m_genjet_subjet_Pz->clear();
  m_genjet_subjet_NCH->clear();
  m_genjet_subjet_NCH_trackPtMin->clear();
  m_genjet_subjet_NPh->clear();
  m_genjet_subjet_NNH->clear();
  m_genjet_subjet_jetindex->clear();

  m_genjet_subjet_jetChargeE_kappa_0_25->clear();
  m_genjet_subjet_jetChargeE_kappa_0_50->clear();
  //m_genjet_subjet_jetChargeE_kappa_0_75->clear();
  //m_genjet_subjet_jetChargeE_kappa_1_00->clear();
  //m_genjet_subjet_jetChargeE_kappa_0_10->clear();
  //m_genjet_subjet_jetChargeE_kappa_0_15->clear();
  m_genjet_subjet_jetChargeE_kappa_0_20->clear();
  m_genjet_subjet_jetChargeE_kappa_0_30->clear();
  
  m_genjet_subjet_jetChargePt_kappa_0_25->clear();
  m_genjet_subjet_jetChargePt_kappa_0_50->clear();
  //m_genjet_subjet_jetChargePt_kappa_0_75->clear();
  //m_genjet_subjet_jetChargePt_kappa_1_00->clear();
  //m_genjet_subjet_jetChargePt_kappa_0_10->clear();
  //m_genjet_subjet_jetChargePt_kappa_0_15->clear();
  m_genjet_subjet_jetChargePt_kappa_0_20->clear();
  m_genjet_subjet_jetChargePt_kappa_0_30->clear();

  //m_genjet_subjet_jetChargePProj_kappa_0_25->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_50->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_75->clear();
  //m_genjet_subjet_jetChargePProj_kappa_1_00->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_10->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_15->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_20->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_30->clear();

  m_isoPartGenDR10_E->clear();
  m_isoPartGenDR10_Px->clear();
  m_isoPartGenDR10_Py->clear();
  m_isoPartGenDR10_Pz->clear();
  m_isoPartGenDR10_PDGID->clear();
  m_isoPartGenDR10_relIso->clear();
  m_isoPartGenDR10_relIsoCH->clear();
  m_isoPartGenDR10_relIsoPh->clear();
  m_isoPartGenDR10_relIsoNH->clear();
  m_isoPartGenDR10_relIsoEl->clear();
  m_isoPartGenDR10_relIsoMu->clear();
  
  m_isoPartRecoDR10_E->clear();
  m_isoPartRecoDR10_Px->clear();
  m_isoPartRecoDR10_Py->clear();
  m_isoPartRecoDR10_Pz->clear();
  m_isoPartRecoDR10_PDGID->clear();
  m_isoPartRecoDR10_relIso->clear();
  m_isoPartRecoDR10_relIsoCH->clear();
  m_isoPartRecoDR10_relIsoPh->clear();
  m_isoPartRecoDR10_relIsoNH->clear();
  m_isoPartRecoDR10_relIsoEl->clear();
  m_isoPartRecoDR10_relIsoMu->clear();
  
  //m_recojet_jetChargeE_kappa_0_25->clear();
  //m_recojet_jetChargeE_kappa_0_50->clear();
  //m_recojet_jetChargeE_kappa_0_75->clear();
  //m_recojet_jetChargeE_kappa_1_00->clear();
  //m_recojet_jetChargeE_kappa_0_10->clear();
  //m_recojet_jetChargeE_kappa_0_15->clear();
  //m_recojet_jetChargeE_kappa_0_20->clear();
  //m_recojet_jetChargeE_kappa_0_30->clear();
  
  //m_recojet_jetChargePt_kappa_0_25->clear();
  //m_recojet_jetChargePt_kappa_0_50->clear();
  //m_recojet_jetChargePt_kappa_0_75->clear();
  //m_recojet_jetChargePt_kappa_1_00->clear();
  //m_recojet_jetChargePt_kappa_0_10->clear();
  //m_recojet_jetChargePt_kappa_0_15->clear();
  //m_recojet_jetChargePt_kappa_0_20->clear();
  //m_recojet_jetChargePt_kappa_0_30->clear();

  //m_recojet_jetChargePProj_kappa_0_25->clear();
  //m_recojet_jetChargePProj_kappa_0_50->clear();
  //m_recojet_jetChargePProj_kappa_0_75->clear();
  //m_recojet_jetChargePProj_kappa_1_00->clear();
  //m_recojet_jetChargePProj_kappa_0_10->clear();
  //m_recojet_jetChargePProj_kappa_0_15->clear();
  //m_recojet_jetChargePProj_kappa_0_20->clear();
  //m_recojet_jetChargePProj_kappa_0_30->clear();

  /*
  m_recojet_fat2j_rfj4_E->clear();
  m_recojet_fat2j_rfj4_Px->clear();
  m_recojet_fat2j_rfj4_Py->clear();
  m_recojet_fat2j_rfj4_Pz->clear();
  m_recojet_fat2j_rfj4_Mult->clear();
  */

  m_recojet_E->clear();
  m_recojet_Px->clear();
  m_recojet_Py->clear();
  m_recojet_Pz->clear();
  m_recojet_Mult->clear();
  m_recojet_NCH->clear();
  m_recojet_NCH_trackPtMin->clear();
  m_recojet_NPh->clear();
  m_recojet_NNH->clear();
  m_recojet_dij_21->clear();
  m_recojet_dij_32->clear();
  m_recojet_dij_43->clear();
  m_recojet_dij_21_max->clear();
  m_recojet_dij_32_max->clear();
  m_recojet_dij_43_max->clear();
  m_recojet_CHFraction->clear();
  m_recojet_CHFraction_trackPtMin->clear();
  m_recojet_PhFraction->clear();
  m_recojet_ElFraction->clear();
  m_recojet_MuFraction->clear();
  m_recojet_NHFraction->clear();
  m_recojet_nsubjettiness1->clear();
  m_recojet_nsubjettiness2->clear();
  m_recojet_nsubjettiness3->clear();
  m_recojet_nsubjettiness1_lrz->clear();
  m_recojet_nsubjettiness2_lrz->clear();
  m_recojet_nsubjettiness3_lrz->clear();
  m_recojet_beta1_ECorr2->clear();
  m_recojet_beta1_ECorr3->clear();
  m_recojet_beta1_N2->clear();
  m_recojet_beta1_N3->clear();
  m_recojet_beta1_C2->clear();
  m_recojet_beta1_C3->clear();
  m_recojet_beta1_D2->clear();
  m_recojet_beta1_ECorr2_E_theta->clear();
  m_recojet_beta1_ECorr3_E_theta->clear();
  m_recojet_beta1_N2_E_theta->clear();
  m_recojet_beta1_N3_E_theta->clear();
  m_recojet_beta1_C2_E_theta->clear();
  m_recojet_beta1_C3_E_theta->clear();
  m_recojet_beta1_D2_E_theta->clear();

  m_recojet_beta2_ECorr2->clear();
  m_recojet_beta2_ECorr3->clear();
  m_recojet_beta2_N2->clear();
  m_recojet_beta2_N3->clear();
  m_recojet_beta2_C2->clear();
  m_recojet_beta2_C3->clear();
  m_recojet_beta2_D2->clear();
  m_recojet_beta2_ECorr2_E_theta->clear();
  m_recojet_beta2_ECorr3_E_theta->clear();
  m_recojet_beta2_N2_E_theta->clear();
  m_recojet_beta2_N3_E_theta->clear();
  m_recojet_beta2_C2_E_theta->clear();
  m_recojet_beta2_C3_E_theta->clear();
  m_recojet_beta2_D2_E_theta->clear();

  m_recojet_beta0_5_ECorr2->clear();
  m_recojet_beta0_5_ECorr3->clear();
  m_recojet_beta0_5_N2->clear();
  m_recojet_beta0_5_N3->clear();
  m_recojet_beta0_5_C2->clear();
  m_recojet_beta0_5_C3->clear();
  m_recojet_beta0_5_D2->clear();
  m_recojet_beta0_5_ECorr2_E_theta->clear();
  m_recojet_beta0_5_ECorr3_E_theta->clear();
  m_recojet_beta0_5_N2_E_theta->clear();
  m_recojet_beta0_5_N3_E_theta->clear();
  m_recojet_beta0_5_C2_E_theta->clear();
  m_recojet_beta0_5_C3_E_theta->clear();
  m_recojet_beta0_5_D2_E_theta->clear();


  m_recojet_subjet_E->clear();
  m_recojet_subjet_Px->clear();
  m_recojet_subjet_Py->clear();
  m_recojet_subjet_Pz->clear();
  m_recojet_subjet_NCH->clear();
  m_recojet_subjet_NCH_trackPtMin->clear();
  m_recojet_subjet_NPh->clear();
  m_recojet_subjet_NNH->clear();
  m_recojet_subjet_jetindex->clear();
  m_recojet_subjet_CHFraction->clear();
  m_recojet_subjet_CHFraction_trackPtMin->clear();
  m_recojet_subjet_PhFraction->clear();
  m_recojet_subjet_ElFraction->clear();
  m_recojet_subjet_MuFraction->clear();
  m_recojet_subjet_NHFraction->clear();

  m_recojet_subjet_jetChargeE_kappa_0_25->clear();
  m_recojet_subjet_jetChargeE_kappa_0_50->clear();
  //m_recojet_subjet_jetChargeE_kappa_0_75->clear();
  //m_recojet_subjet_jetChargeE_kappa_1_00->clear();
  //m_recojet_subjet_jetChargeE_kappa_0_10->clear();
  //m_recojet_subjet_jetChargeE_kappa_0_15->clear();
  m_recojet_subjet_jetChargeE_kappa_0_20->clear();
  m_recojet_subjet_jetChargeE_kappa_0_30->clear();
  
  m_recojet_subjet_jetChargePt_kappa_0_25->clear();
  m_recojet_subjet_jetChargePt_kappa_0_50->clear();
  //m_recojet_subjet_jetChargePt_kappa_0_75->clear();
  //m_recojet_subjet_jetChargePt_kappa_1_00->clear();
  //m_recojet_subjet_jetChargePt_kappa_0_10->clear();
  //m_recojet_subjet_jetChargePt_kappa_0_15->clear();
  m_recojet_subjet_jetChargePt_kappa_0_20->clear();
  m_recojet_subjet_jetChargePt_kappa_0_30->clear();

  //m_recojet_subjet_jetChargePProj_kappa_0_25->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_50->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_75->clear();
  //m_recojet_subjet_jetChargePProj_kappa_1_00->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_10->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_15->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_20->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_30->clear();

  /*
  m_recojet_rfj_4jets_E->clear();
  m_recojet_rfj_4jets_Px->clear();
  m_recojet_rfj_4jets_Py->clear();
  m_recojet_rfj_4jets_Pz->clear();
  m_recojet_rfj_4jets_Mult->clear();
  m_recojet_rfj_4jets_NCH->clear();
  m_recojet_rfj_4jets_CHFraction->clear();
  m_recojet_rfj_4jets_PhFraction->clear();
  m_recojet_rfj_4jets_ElFraction->clear();
  m_recojet_rfj_4jets_MuFraction->clear();
  m_recojet_rfj_4jets_NHFraction->clear();

  m_recojet_rfj_4jets_BTag->clear();
  m_recojet_rfj_4jets_CTag->clear();
  m_recojet_rfj_4jets_OTag->clear();
  m_recojet_rfj_4jets_cat->clear();

  m_recojet_rfj_4jets_svtx_r->clear();
  m_recojet_rfj_4jets_svtx_E->clear();
  m_recojet_rfj_4jets_svtx_Mass->clear();
  m_recojet_rfj_4jets_svtx_nTrack->clear();
  m_recojet_rfj_4jets_svtx_Charge->clear();
  */
  
  m_recojet_subjet_rfj_j_E->clear();
  m_recojet_subjet_rfj_j_Px->clear();
  m_recojet_subjet_rfj_j_Py->clear();
  m_recojet_subjet_rfj_j_Pz->clear();
  m_recojet_subjet_rfj_j_NCH->clear();
  m_recojet_subjet_rfj_j_NPh->clear();
  m_recojet_subjet_rfj_j_NNH->clear();
  m_recojet_subjet_rfj_j_jetindex->clear();
  m_recojet_subjet_rfj_j_subjetindex->clear();
  m_recojet_subjet_rfj_j_CHFraction->clear();
  m_recojet_subjet_rfj_j_PhFraction->clear();
  m_recojet_subjet_rfj_j_ElFraction->clear();
  m_recojet_subjet_rfj_j_MuFraction->clear();
  m_recojet_subjet_rfj_j_NHFraction->clear();
  /*
  m_recojet_subjet_rfj_j_jetChargeE_kappa_0_25->clear();
  m_recojet_subjet_rfj_j_jetChargeE_kappa_0_50->clear();
  m_recojet_subjet_rfj_j_jetChargeE_kappa_0_20->clear();
  m_recojet_subjet_rfj_j_jetChargeE_kappa_0_30->clear();
  
  m_recojet_subjet_rfj_j_jetChargePt_kappa_0_25->clear();
  m_recojet_subjet_rfj_j_jetChargePt_kappa_0_50->clear();
  m_recojet_subjet_rfj_j_jetChargePt_kappa_0_20->clear();
  m_recojet_subjet_rfj_j_jetChargePt_kappa_0_30->clear();
  */
  m_recojet_subjet_rfj_j_BTag->clear();
  m_recojet_subjet_rfj_j_CTag->clear();
  m_recojet_subjet_rfj_j_OTag->clear();
  m_recojet_subjet_rfj_j_cat->clear();
  
  m_recojet_subjet_rfj_j_svtx_r->clear();
  m_recojet_subjet_rfj_j_svtx_E->clear();
  m_recojet_subjet_rfj_j_svtx_Mass->clear();
  m_recojet_subjet_rfj_j_svtx_nTrack->clear();
  m_recojet_subjet_rfj_j_svtx_Charge->clear();
  m_recojet_subjet_rfj_j_svtx_jetindex->clear();

  m_reco_tau_E->clear();
  m_reco_tau_Px->clear();
  m_reco_tau_Py->clear();
  m_reco_tau_Pz->clear();
  m_reco_tau_Mult->clear();
  m_reco_tau_Charge->clear();
  m_reco_tau_NCH->clear();
  m_reco_tau_NPh->clear();
  m_reco_tau_NNH->clear();
  m_reco_tau_CHFraction->clear();
  m_reco_tau_PhFraction->clear();
  m_reco_tau_ElFraction->clear();
  m_reco_tau_MuFraction->clear();
  m_reco_tau_NHFraction->clear();

  m_gen_tau_E->clear();
  m_gen_tau_Px->clear();
  m_gen_tau_Py->clear();
  m_gen_tau_Pz->clear();
  m_gen_tau_Mult->clear();
  m_gen_tau_Charge->clear();
  m_gen_tau_MCCharge->clear();
  m_gen_tau_NCH->clear();
  m_gen_tau_NPh->clear();
  m_gen_tau_NNH->clear();
  m_gen_tau_MotherPDGID->clear();
  m_gen_tau_CHFraction->clear();
  m_gen_tau_PhFraction->clear();
  m_gen_tau_ElFraction->clear();
  m_gen_tau_MuFraction->clear();
  m_gen_tau_NHFraction->clear();


  m_outputTree->Branch("runNumber", &m_runNumber,"runNumber/I");
  m_outputTree->Branch("eventNumber", &m_eventNumber,"eventNumber/I");
  m_outputTree->Branch("totPFO_E", &m_totPFO_E,"totPFO_E/F");
  m_outputTree->Branch("totPFO_Px", &m_totPFO_Px,"totPFO_Px/F");
  m_outputTree->Branch("totPFO_Py", &m_totPFO_Py,"totPFO_Py/F");
  m_outputTree->Branch("totPFO_Pz", &m_totPFO_Pz,"totPFO_Pz/F");
  m_outputTree->Branch("totPFO_Mult", &m_totPFO_Mult,"totPFO_Mult/I");
  m_outputTree->Branch("totPFO_NCH", &m_totPFO_NCH,"totPFO_NCH/I");
  m_outputTree->Branch("totPFO_NPh", &m_totPFO_NPh,"totPFO_NPh/I");
  m_outputTree->Branch("totPFO_NNH", &m_totPFO_NNH,"totPFO_NNH/I");
  m_outputTree->Branch("totPFO_CHFraction", &m_totPFO_CHFraction,"totPFO_CHFraction/F");
  m_outputTree->Branch("totPFO_PhFraction", &m_totPFO_PhFraction,"totPFO_PhFraction/F");
  m_outputTree->Branch("totPFO_ElFraction", &m_totPFO_ElFraction,"totPFO_ElFraction/F");
  m_outputTree->Branch("totPFO_MuFraction", &m_totPFO_MuFraction,"totPFO_MuFraction/F");
  m_outputTree->Branch("totPFO_NHFraction", &m_totPFO_NHFraction,"totPFO_NHFraction/F");

  m_outputTree->Branch("gen_y21",&m_gen_y21,"gen_y21/F");
  m_outputTree->Branch("gen_y32",&m_gen_y32,"gen_y21/F");
  m_outputTree->Branch("gen_y43",&m_gen_y43,"gen_y21/F");

  m_outputTree->Branch("gen_y21_max",&m_gen_y21_max,"gen_y21_max/F");
  m_outputTree->Branch("gen_y32_max",&m_gen_y32_max,"gen_y32_max/F");
  m_outputTree->Branch("gen_y43_max",&m_gen_y43_max,"gen_y43_max/F");

  m_outputTree->Branch("reco_y21",&m_reco_y21,"reco_y21/F");
  m_outputTree->Branch("reco_y32",&m_reco_y32,"reco_y21/F");
  m_outputTree->Branch("reco_y43",&m_reco_y43,"reco_y21/F");

  m_outputTree->Branch("reco_y21_max",&m_reco_y21_max,"reco_y21_max/F");
  m_outputTree->Branch("reco_y32_max",&m_reco_y32_max,"reco_y32_max/F");
  m_outputTree->Branch("reco_y43_max",&m_reco_y43_max,"reco_y43_max/F");

  m_outputTree->Branch("true_E", &m_true_E,"true_E/F");
  m_outputTree->Branch("true_Px", &m_true_Px,"true_Px/F");
  m_outputTree->Branch("true_Py", &m_true_Py,"true_Py/F");
  m_outputTree->Branch("true_Pz", &m_true_Pz,"true_Pz/F");
  m_outputTree->Branch("true_Mult", &m_true_Mult,"true_Mult/I");
  m_outputTree->Branch("true_NCH", &m_true_NCH,"true_NCH/I");
  m_outputTree->Branch("true_NPh", &m_true_NPh,"true_NPh/I");
  m_outputTree->Branch("true_NNH", &m_true_NNH,"true_NNH/I");
  m_outputTree->Branch("true_CHFraction", &m_true_CHFraction,"true_CHFraction/F");
  m_outputTree->Branch("true_PhFraction", &m_true_PhFraction,"true_PhFraction/F");
  m_outputTree->Branch("true_ElFraction", &m_true_ElFraction,"true_ElFraction/F");
  m_outputTree->Branch("true_MuFraction", &m_true_MuFraction,"true_MuFraction/F");
  m_outputTree->Branch("true_NHFraction", &m_true_NHFraction,"true_NHFraction/F");

  m_outputTree->Branch("true_inv_E", &m_true_inv_E,"true_inv_E/F");
  m_outputTree->Branch("true_inv_Px", &m_true_inv_Px,"true_inv_Px/F");
  m_outputTree->Branch("true_inv_Py", &m_true_inv_Py,"true_inv_Py/F");
  m_outputTree->Branch("true_inv_Pz", &m_true_inv_Pz,"true_inv_Pz/F");
  m_outputTree->Branch("true_inv_Mult", &m_true_inv_Mult,"true_inv_Mult/I");
  if(m_saveMEInfo){
    m_outputTree->Branch("trueME_Px", "std::vector< float >", &m_trueME_Px); 
    m_outputTree->Branch("trueME_Py", "std::vector< float >", &m_trueME_Py); 
    m_outputTree->Branch("trueME_Pz", "std::vector< float >", &m_trueME_Pz); 
    m_outputTree->Branch("trueME_E", "std::vector< float >", &m_trueME_E); 
    m_outputTree->Branch("trueME_PDGID", "std::vector< int >", &m_trueME_PDGID); 
    //m_outputTree->Branch("trueME_index", "std::vector< int >", &m_trueME_index);
  }
  m_outputTree->Branch("genTrueLepPh_E", "std::vector< float >", &m_genTrueLepPh_E);
  m_outputTree->Branch("genTrueLepPh_Px", "std::vector< float >", &m_genTrueLepPh_Px);
  m_outputTree->Branch("genTrueLepPh_Py", "std::vector< float >", &m_genTrueLepPh_Py);
  m_outputTree->Branch("genTrueLepPh_Pz", "std::vector< float >", &m_genTrueLepPh_Pz);
  m_outputTree->Branch("genTrueLepPh_PDGID", "std::vector< int >", &m_genTrueLepPh_PDGID);
  //actual ancestor, i.e. Higgs,W
  m_outputTree->Branch("genTrueLepPh_ANC_PDGID", "std::vector< int >", &m_genTrueLepPh_ANC_PDGID);
  m_outputTree->Branch("genTrueLepPh_relIso", "std::vector< float >", &m_genTrueLepPh_relIso);
  m_outputTree->Branch("genTrueLepPh_relIsoCH", "std::vector< float >", &m_genTrueLepPh_relIsoCH);
  m_outputTree->Branch("genTrueLepPh_relIsoPh", "std::vector< float >", &m_genTrueLepPh_relIsoPh);
  m_outputTree->Branch("genTrueLepPh_relIsoNH", "std::vector< float >", &m_genTrueLepPh_relIsoNH);
  m_outputTree->Branch("genTrueLepPh_relIsoEl", "std::vector< float >", &m_genTrueLepPh_relIsoEl);
  m_outputTree->Branch("genTrueLepPh_relIsoMu", "std::vector< float >", &m_genTrueLepPh_relIsoMu);
  
  //m_outputTree->Branch("genjet_jetChargeE_kappa_0_25", "std::vector< float >", &m_genjet_jetChargeE_kappa_0_25);
  //m_outputTree->Branch("genjet_jetChargeE_kappa_0_50", "std::vector< float >", &m_genjet_jetChargeE_kappa_0_50);
  //m_outputTree->Branch("genjet_jetChargeE_kappa_0_75", "std::vector< float >", &m_genjet_jetChargeE_kappa_0_75);
  //m_outputTree->Branch("genjet_jetChargeE_kappa_1_00", "std::vector< float >", &m_genjet_jetChargeE_kappa_1_00);
  //m_outputTree->Branch("genjet_jetChargeE_kappa_0_10", "std::vector< float >", &m_genjet_jetChargeE_kappa_0_10);
  //m_outputTree->Branch("genjet_jetChargeE_kappa_0_15", "std::vector< float >", &m_genjet_jetChargeE_kappa_0_15);
  //m_outputTree->Branch("genjet_jetChargeE_kappa_0_20", "std::vector< float >", &m_genjet_jetChargeE_kappa_0_20);
  //m_outputTree->Branch("genjet_jetChargeE_kappa_0_30", "std::vector< float >", &m_genjet_jetChargeE_kappa_0_30);

  //m_outputTree->Branch("genjet_jetChargePt_kappa_0_25", "std::vector< float >", &m_genjet_jetChargePt_kappa_0_25);
  //m_outputTree->Branch("genjet_jetChargePt_kappa_0_50", "std::vector< float >", &m_genjet_jetChargePt_kappa_0_50);
  //m_outputTree->Branch("genjet_jetChargePt_kappa_0_75", "std::vector< float >", &m_genjet_jetChargePt_kappa_0_75);
  //m_outputTree->Branch("genjet_jetChargePt_kappa_1_00", "std::vector< float >", &m_genjet_jetChargePt_kappa_1_00);
  //m_outputTree->Branch("genjet_jetChargePt_kappa_0_10", "std::vector< float >", &m_genjet_jetChargePt_kappa_0_10);
  //m_outputTree->Branch("genjet_jetChargePt_kappa_0_15", "std::vector< float >", &m_genjet_jetChargePt_kappa_0_15);
  //m_outputTree->Branch("genjet_jetChargePt_kappa_0_20", "std::vector< float >", &m_genjet_jetChargePt_kappa_0_20);
  //m_outputTree->Branch("genjet_jetChargePt_kappa_0_30", "std::vector< float >", &m_genjet_jetChargePt_kappa_0_30);

  //m_outputTree->Branch("genjet_jetChargePProj_kappa_0_25", "std::vector< float >", &m_genjet_jetChargePProj_kappa_0_25);
  //m_outputTree->Branch("genjet_jetChargePProj_kappa_0_50", "std::vector< float >", &m_genjet_jetChargePProj_kappa_0_50);
  //m_outputTree->Branch("genjet_jetChargePProj_kappa_0_75", "std::vector< float >", &m_genjet_jetChargePProj_kappa_0_75);
  //m_outputTree->Branch("genjet_jetChargePProj_kappa_1_00", "std::vector< float >", &m_genjet_jetChargePProj_kappa_1_00);
  //m_outputTree->Branch("genjet_jetChargePProj_kappa_0_10", "std::vector< float >", &m_genjet_jetChargePProj_kappa_0_10);
  //m_outputTree->Branch("genjet_jetChargePProj_kappa_0_15", "std::vector< float >", &m_genjet_jetChargePProj_kappa_0_15);
  //m_outputTree->Branch("genjet_jetChargePProj_kappa_0_20", "std::vector< float >", &m_genjet_jetChargePProj_kappa_0_20);
  //m_outputTree->Branch("genjet_jetChargePProj_kappa_0_30", "std::vector< float >", &m_genjet_jetChargePProj_kappa_0_30);

  m_outputTree->Branch("genjet_E", "std::vector< float >", &m_genjet_E);
  m_outputTree->Branch("genjet_Px", "std::vector< float >", &m_genjet_Px);
  m_outputTree->Branch("genjet_Py", "std::vector< float >", &m_genjet_Py);
  m_outputTree->Branch("genjet_Pz", "std::vector< float >", &m_genjet_Pz);
  m_outputTree->Branch("genjet_Mult", "std::vector< int >", &m_genjet_Mult);
  m_outputTree->Branch("genjet_NCH", "std::vector< int >", &m_genjet_NCH);
  m_outputTree->Branch("genjet_NCH_trackPtMin", "std::vector< int >", &m_genjet_NCH_trackPtMin);
  m_outputTree->Branch("genjet_NPh", "std::vector< int >", &m_genjet_NPh);
  m_outputTree->Branch("genjet_NNH", "std::vector< int >", &m_genjet_NNH);

  m_outputTree->Branch("genjet_dij_21","std::vector< float >", &m_genjet_dij_21);
  m_outputTree->Branch("genjet_dij_32","std::vector< float >", &m_genjet_dij_32);
  m_outputTree->Branch("genjet_dij_43","std::vector< float >", &m_genjet_dij_43);

  m_outputTree->Branch("genjet_dij_21_max","std::vector< float >", &m_genjet_dij_21_max);
  m_outputTree->Branch("genjet_dij_32_max","std::vector< float >", &m_genjet_dij_32_max);
  m_outputTree->Branch("genjet_dij_43_max","std::vector< float >", &m_genjet_dij_43_max);

  m_outputTree->Branch("genjet_CHFraction", "std::vector< float >", &m_genjet_CHFraction);
  m_outputTree->Branch("genjet_CHFraction_trackPtMin", "std::vector< float >", &m_genjet_CHFraction_trackPtMin);
  m_outputTree->Branch("genjet_PhFraction", "std::vector< float >", &m_genjet_PhFraction);
  m_outputTree->Branch("genjet_ElFraction", "std::vector< float >", &m_genjet_ElFraction);
  m_outputTree->Branch("genjet_MuFraction", "std::vector< float >", &m_genjet_MuFraction);
  m_outputTree->Branch("genjet_NHFraction", "std::vector< float >", &m_genjet_NHFraction);
  m_outputTree->Branch("genjet_nsubjettiness1", "std::vector< float >", m_genjet_nsubjettiness1);
  m_outputTree->Branch("genjet_nsubjettiness2", "std::vector< float >", m_genjet_nsubjettiness2);
  m_outputTree->Branch("genjet_nsubjettiness3", "std::vector< float >", m_genjet_nsubjettiness3); 
  m_outputTree->Branch("genjet_nsubjettiness1_lrz", "std::vector< float >", m_genjet_nsubjettiness1_lrz);
  m_outputTree->Branch("genjet_nsubjettiness2_lrz", "std::vector< float >", m_genjet_nsubjettiness2_lrz);
  m_outputTree->Branch("genjet_nsubjettiness3_lrz", "std::vector< float >", m_genjet_nsubjettiness3_lrz); 

  m_outputTree->Branch("genjet_beta1_ECorr2", "std::vector< float >", &m_genjet_beta1_ECorr2);
  m_outputTree->Branch("genjet_beta1_ECorr3", "std::vector< float >", &m_genjet_beta1_ECorr3);
  m_outputTree->Branch("genjet_beta1_N2", "std::vector< float >", &m_genjet_beta1_N2);
  m_outputTree->Branch("genjet_beta1_N3", "std::vector< float >", &m_genjet_beta1_N3);
  m_outputTree->Branch("genjet_beta1_C2", "std::vector< float >", &m_genjet_beta1_C2);
  m_outputTree->Branch("genjet_beta1_C3", "std::vector< float >", &m_genjet_beta1_C3);
  m_outputTree->Branch("genjet_beta1_D2", "std::vector< float >", &m_genjet_beta1_D2);
  m_outputTree->Branch("genjet_beta1_ECorr2_E_theta", "std::vector< float >", &m_genjet_beta1_ECorr2_E_theta);
  m_outputTree->Branch("genjet_beta1_ECorr3_E_theta", "std::vector< float >", &m_genjet_beta1_ECorr3_E_theta);
  m_outputTree->Branch("genjet_beta1_N2_E_theta", "std::vector< float >", &m_genjet_beta1_N2_E_theta);
  m_outputTree->Branch("genjet_beta1_N3_E_theta", "std::vector< float >", &m_genjet_beta1_N3_E_theta);
  m_outputTree->Branch("genjet_beta1_C2_E_theta", "std::vector< float >", &m_genjet_beta1_C2_E_theta);
  m_outputTree->Branch("genjet_beta1_C3_E_theta", "std::vector< float >", &m_genjet_beta1_C3_E_theta);
  m_outputTree->Branch("genjet_beta1_D2_E_theta", "std::vector< float >", &m_genjet_beta1_D2_E_theta);

  m_outputTree->Branch("genjet_beta2_ECorr2", "std::vector< float >", &m_genjet_beta2_ECorr2);
  m_outputTree->Branch("genjet_beta2_ECorr3", "std::vector< float >", &m_genjet_beta2_ECorr3);
  m_outputTree->Branch("genjet_beta2_N2", "std::vector< float >", &m_genjet_beta2_N2);
  m_outputTree->Branch("genjet_beta2_N3", "std::vector< float >", &m_genjet_beta2_N3);
  m_outputTree->Branch("genjet_beta2_C2", "std::vector< float >", &m_genjet_beta2_C2);
  m_outputTree->Branch("genjet_beta2_C3", "std::vector< float >", &m_genjet_beta2_C3);
  m_outputTree->Branch("genjet_beta2_D2", "std::vector< float >", &m_genjet_beta2_D2);
  m_outputTree->Branch("genjet_beta2_ECorr2_E_theta", "std::vector< float >", &m_genjet_beta2_ECorr2_E_theta);
  m_outputTree->Branch("genjet_beta2_ECorr3_E_theta", "std::vector< float >", &m_genjet_beta2_ECorr3_E_theta);
  m_outputTree->Branch("genjet_beta2_N2_E_theta", "std::vector< float >", &m_genjet_beta2_N2_E_theta);
  m_outputTree->Branch("genjet_beta2_N3_E_theta", "std::vector< float >", &m_genjet_beta2_N3_E_theta);
  m_outputTree->Branch("genjet_beta2_C2_E_theta", "std::vector< float >", &m_genjet_beta2_C2_E_theta);
  m_outputTree->Branch("genjet_beta2_C3_E_theta", "std::vector< float >", &m_genjet_beta2_C3_E_theta);
  m_outputTree->Branch("genjet_beta2_D2_E_theta", "std::vector< float >", &m_genjet_beta2_D2_E_theta);

  m_outputTree->Branch("genjet_beta0_5_ECorr2", "std::vector< float >", &m_genjet_beta0_5_ECorr2);
  m_outputTree->Branch("genjet_beta0_5_ECorr3", "std::vector< float >", &m_genjet_beta0_5_ECorr3);
  m_outputTree->Branch("genjet_beta0_5_N2", "std::vector< float >", &m_genjet_beta0_5_N2);
  m_outputTree->Branch("genjet_beta0_5_N3", "std::vector< float >", &m_genjet_beta0_5_N3);
  m_outputTree->Branch("genjet_beta0_5_C2", "std::vector< float >", &m_genjet_beta0_5_C2);
  m_outputTree->Branch("genjet_beta0_5_C3", "std::vector< float >", &m_genjet_beta0_5_C3);
  m_outputTree->Branch("genjet_beta0_5_D2", "std::vector< float >", &m_genjet_beta0_5_D2);
  m_outputTree->Branch("genjet_beta0_5_ECorr2_E_theta", "std::vector< float >", &m_genjet_beta0_5_ECorr2_E_theta);
  m_outputTree->Branch("genjet_beta0_5_ECorr3_E_theta", "std::vector< float >", &m_genjet_beta0_5_ECorr3_E_theta);
  m_outputTree->Branch("genjet_beta0_5_N2_E_theta", "std::vector< float >", &m_genjet_beta0_5_N2_E_theta);
  m_outputTree->Branch("genjet_beta0_5_N3_E_theta", "std::vector< float >", &m_genjet_beta0_5_N3_E_theta);
  m_outputTree->Branch("genjet_beta0_5_C2_E_theta", "std::vector< float >", &m_genjet_beta0_5_C2_E_theta);
  m_outputTree->Branch("genjet_beta0_5_C3_E_theta", "std::vector< float >", &m_genjet_beta0_5_C3_E_theta);
  m_outputTree->Branch("genjet_beta0_5_D2_E_theta", "std::vector< float >", &m_genjet_beta0_5_D2_E_theta);

  m_outputTree->Branch("genjet_subjet_E", "std::vector< float >", &m_genjet_subjet_E);
  m_outputTree->Branch("genjet_subjet_Px", "std::vector< float >", &m_genjet_subjet_Px);
  m_outputTree->Branch("genjet_subjet_Py", "std::vector< float >", &m_genjet_subjet_Py);
  m_outputTree->Branch("genjet_subjet_Pz", "std::vector< float >", &m_genjet_subjet_Pz);
  m_outputTree->Branch("genjet_subjet_NCH", "std::vector< int >", &m_genjet_subjet_NCH);
  m_outputTree->Branch("genjet_subjet_NCH_trackPtMin", "std::vector< int >", &m_genjet_subjet_NCH_trackPtMin);
  m_outputTree->Branch("genjet_subjet_NPh", "std::vector< int >", &m_genjet_subjet_NPh);
  m_outputTree->Branch("genjet_subjet_NNH", "std::vector< int >", &m_genjet_subjet_NNH);
  m_outputTree->Branch("genjet_subjet_jetindex", "std::vector< int >", &m_genjet_subjet_jetindex);
  m_outputTree->Branch("genjet_subjet_CHFraction", "std::vector< float >", &m_genjet_subjet_CHFraction);
  m_outputTree->Branch("genjet_subjet_CHFraction_trackPtMin", "std::vector< float >", &m_genjet_subjet_CHFraction_trackPtMin);
  m_outputTree->Branch("genjet_subjet_PhFraction", "std::vector< float >", &m_genjet_subjet_PhFraction);
  m_outputTree->Branch("genjet_subjet_ElFraction", "std::vector< float >", &m_genjet_subjet_ElFraction);
  m_outputTree->Branch("genjet_subjet_MuFraction", "std::vector< float >", &m_genjet_subjet_MuFraction);
  m_outputTree->Branch("genjet_subjet_NHFraction", "std::vector< float >", &m_genjet_subjet_NHFraction);

  m_outputTree->Branch("genjet_subjet_jetChargeE_kappa_0_25", "std::vector< float >", &m_genjet_subjet_jetChargeE_kappa_0_25);
  m_outputTree->Branch("genjet_subjet_jetChargeE_kappa_0_50", "std::vector< float >", &m_genjet_subjet_jetChargeE_kappa_0_50);
  //m_outputTree->Branch("genjet_subjet_jetChargeE_kappa_0_75", "std::vector< float >", &m_genjet_subjet_jetChargeE_kappa_0_75);
  //m_outputTree->Branch("genjet_subjet_jetChargeE_kappa_1_00", "std::vector< float >", &m_genjet_subjet_jetChargeE_kappa_1_00);
  //m_outputTree->Branch("genjet_subjet_jetChargeE_kappa_0_10", "std::vector< float >", &m_genjet_subjet_jetChargeE_kappa_0_10);
  //m_outputTree->Branch("genjet_subjet_jetChargeE_kappa_0_15", "std::vector< float >", &m_genjet_subjet_jetChargeE_kappa_0_15);
  m_outputTree->Branch("genjet_subjet_jetChargeE_kappa_0_20", "std::vector< float >", &m_genjet_subjet_jetChargeE_kappa_0_20);
  m_outputTree->Branch("genjet_subjet_jetChargeE_kappa_0_30", "std::vector< float >", &m_genjet_subjet_jetChargeE_kappa_0_30);

  m_outputTree->Branch("genjet_subjet_jetChargePt_kappa_0_25", "std::vector< float >", &m_genjet_subjet_jetChargePt_kappa_0_25);
  m_outputTree->Branch("genjet_subjet_jetChargePt_kappa_0_50", "std::vector< float >", &m_genjet_subjet_jetChargePt_kappa_0_50);
  //m_outputTree->Branch("genjet_subjet_jetChargePt_kappa_0_75", "std::vector< float >", &m_genjet_subjet_jetChargePt_kappa_0_75);
  //m_outputTree->Branch("genjet_subjet_jetChargePt_kappa_1_00", "std::vector< float >", &m_genjet_subjet_jetChargePt_kappa_1_00);
  //m_outputTree->Branch("genjet_subjet_jetChargePt_kappa_0_10", "std::vector< float >", &m_genjet_subjet_jetChargePt_kappa_0_10);
  //m_outputTree->Branch("genjet_subjet_jetChargePt_kappa_0_15", "std::vector< float >", &m_genjet_subjet_jetChargePt_kappa_0_15);
  m_outputTree->Branch("genjet_subjet_jetChargePt_kappa_0_20", "std::vector< float >", &m_genjet_subjet_jetChargePt_kappa_0_20);
  m_outputTree->Branch("genjet_subjet_jetChargePt_kappa_0_30", "std::vector< float >", &m_genjet_subjet_jetChargePt_kappa_0_30);

  //m_outputTree->Branch("genjet_subjet_jetChargePProj_kappa_0_25", "std::vector< float >", &m_genjet_subjet_jetChargePProj_kappa_0_25);
  //m_outputTree->Branch("genjet_subjet_jetChargePProj_kappa_0_50", "std::vector< float >", &m_genjet_subjet_jetChargePProj_kappa_0_50);
  //m_outputTree->Branch("genjet_subjet_jetChargePProj_kappa_0_75", "std::vector< float >", &m_genjet_subjet_jetChargePProj_kappa_0_75);
  //m_outputTree->Branch("genjet_subjet_jetChargePProj_kappa_1_00", "std::vector< float >", &m_genjet_subjet_jetChargePProj_kappa_1_00);
  //m_outputTree->Branch("genjet_subjet_jetChargePProj_kappa_0_10", "std::vector< float >", &m_genjet_subjet_jetChargePProj_kappa_0_10);
  //m_outputTree->Branch("genjet_subjet_jetChargePProj_kappa_0_15", "std::vector< float >", &m_genjet_subjet_jetChargePProj_kappa_0_15);
  //m_outputTree->Branch("genjet_subjet_jetChargePProj_kappa_0_20", "std::vector< float >", &m_genjet_subjet_jetChargePProj_kappa_0_20);
  //m_outputTree->Branch("genjet_subjet_jetChargePProj_kappa_0_30", "std::vector< float >", &m_genjet_subjet_jetChargePProj_kappa_0_30);
  
  m_outputTree->Branch("isoPartGenDR10_E", "std::vector< float >", &m_isoPartGenDR10_E);
  m_outputTree->Branch("isoPartGenDR10_Px", "std::vector< float >", &m_isoPartGenDR10_Px);
  m_outputTree->Branch("isoPartGenDR10_Py", "std::vector< float >", &m_isoPartGenDR10_Py);
  m_outputTree->Branch("isoPartGenDR10_Pz", "std::vector< float >", &m_isoPartGenDR10_Pz);
  m_outputTree->Branch("isoPartGenDR10_PDGID", "std::vector< int >", &m_isoPartGenDR10_PDGID);
  m_outputTree->Branch("isoPartGenDR10_relIso", "std::vector< float >", &m_isoPartGenDR10_relIso);
  m_outputTree->Branch("isoPartGenDR10_relIsoCH", "std::vector< float >", &m_isoPartGenDR10_relIsoCH);
  m_outputTree->Branch("isoPartGenDR10_relIsoPh", "std::vector< float >", &m_isoPartGenDR10_relIsoPh);
  m_outputTree->Branch("isoPartGenDR10_relIsoNH", "std::vector< float >", &m_isoPartGenDR10_relIsoNH);
  m_outputTree->Branch("isoPartGenDR10_relIsoEl", "std::vector< float >", &m_isoPartGenDR10_relIsoEl);
  m_outputTree->Branch("isoPartGenDR10_relIsoMu", "std::vector< float >", &m_isoPartGenDR10_relIsoMu);
  
  m_outputTree->Branch("isoPartRecoDR10_E", "std::vector< float >", &m_isoPartRecoDR10_E);
  m_outputTree->Branch("isoPartRecoDR10_Px", "std::vector< float >", &m_isoPartRecoDR10_Px);
  m_outputTree->Branch("isoPartRecoDR10_Py", "std::vector< float >", &m_isoPartRecoDR10_Py);
  m_outputTree->Branch("isoPartRecoDR10_Pz", "std::vector< float >", &m_isoPartRecoDR10_Pz);
  m_outputTree->Branch("isoPartRecoDR10_PDGID", "std::vector< int >", &m_isoPartRecoDR10_PDGID);
  m_outputTree->Branch("isoPartRecoDR10_relIso", "std::vector< float >", &m_isoPartRecoDR10_relIso);
  m_outputTree->Branch("isoPartRecoDR10_relIsoCH", "std::vector< float >", &m_isoPartRecoDR10_relIsoCH);
  m_outputTree->Branch("isoPartRecoDR10_relIsoPh", "std::vector< float >", &m_isoPartRecoDR10_relIsoPh);
  m_outputTree->Branch("isoPartRecoDR10_relIsoNH", "std::vector< float >", &m_isoPartRecoDR10_relIsoNH);
  m_outputTree->Branch("isoPartRecoDR10_relIsoEl", "std::vector< float >", &m_isoPartRecoDR10_relIsoEl);
  m_outputTree->Branch("isoPartRecoDR10_relIsoMu", "std::vector< float >", &m_isoPartRecoDR10_relIsoMu);
  
  //m_outputTree->Branch("recojet_jetChargeE_kappa_0_25", "std::vector< float >", &m_recojet_jetChargeE_kappa_0_25);
  //m_outputTree->Branch("recojet_jetChargeE_kappa_0_50", "std::vector< float >", &m_recojet_jetChargeE_kappa_0_50);
  //m_outputTree->Branch("recojet_jetChargeE_kappa_0_75", "std::vector< float >", &m_recojet_jetChargeE_kappa_0_75);
  //m_outputTree->Branch("recojet_jetChargeE_kappa_1_00", "std::vector< float >", &m_recojet_jetChargeE_kappa_1_00);
  //m_outputTree->Branch("recojet_jetChargeE_kappa_0_10", "std::vector< float >", &m_recojet_jetChargeE_kappa_0_10);
  //m_outputTree->Branch("recojet_jetChargeE_kappa_0_15", "std::vector< float >", &m_recojet_jetChargeE_kappa_0_15);
  //m_outputTree->Branch("recojet_jetChargeE_kappa_0_20", "std::vector< float >", &m_recojet_jetChargeE_kappa_0_20);
  //m_outputTree->Branch("recojet_jetChargeE_kappa_0_30", "std::vector< float >", &m_recojet_jetChargeE_kappa_0_30);

  //m_outputTree->Branch("recojet_jetChargePt_kappa_0_25", "std::vector< float >", &m_recojet_jetChargePt_kappa_0_25);
  //m_outputTree->Branch("recojet_jetChargePt_kappa_0_50", "std::vector< float >", &m_recojet_jetChargePt_kappa_0_50);
  //m_outputTree->Branch("recojet_jetChargePt_kappa_0_75", "std::vector< float >", &m_recojet_jetChargePt_kappa_0_75);
  //m_outputTree->Branch("recojet_jetChargePt_kappa_1_00", "std::vector< float >", &m_recojet_jetChargePt_kappa_1_00);
  //m_outputTree->Branch("recojet_jetChargePt_kappa_0_10", "std::vector< float >", &m_recojet_jetChargePt_kappa_0_10);
  //m_outputTree->Branch("recojet_jetChargePt_kappa_0_15", "std::vector< float >", &m_recojet_jetChargePt_kappa_0_15);
  //m_outputTree->Branch("recojet_jetChargePt_kappa_0_20", "std::vector< float >", &m_recojet_jetChargePt_kappa_0_20);
  //m_outputTree->Branch("recojet_jetChargePt_kappa_0_30", "std::vector< float >", &m_recojet_jetChargePt_kappa_0_30);

  //m_outputTree->Branch("recojet_jetChargePProj_kappa_0_25", "std::vector< float >", &m_recojet_jetChargePProj_kappa_0_25);
  //m_outputTree->Branch("recojet_jetChargePProj_kappa_0_50", "std::vector< float >", &m_recojet_jetChargePProj_kappa_0_50);
  //m_outputTree->Branch("recojet_jetChargePProj_kappa_0_75", "std::vector< float >", &m_recojet_jetChargePProj_kappa_0_75);
  //m_outputTree->Branch("recojet_jetChargePProj_kappa_1_00", "std::vector< float >", &m_recojet_jetChargePProj_kappa_1_00);
  //m_outputTree->Branch("recojet_jetChargePProj_kappa_0_10", "std::vector< float >", &m_recojet_jetChargePProj_kappa_0_10);
  //m_outputTree->Branch("recojet_jetChargePProj_kappa_0_15", "std::vector< float >", &m_recojet_jetChargePProj_kappa_0_15);
  //m_outputTree->Branch("recojet_jetChargePProj_kappa_0_20", "std::vector< float >", &m_recojet_jetChargePProj_kappa_0_20);
  //m_outputTree->Branch("recojet_jetChargePProj_kappa_0_30", "std::vector< float >", &m_recojet_jetChargePProj_kappa_0_30);



  m_outputTree->Branch("recojet_E", "std::vector< float >", &m_recojet_E);
  m_outputTree->Branch("recojet_Px", "std::vector< float >", &m_recojet_Px);
  m_outputTree->Branch("recojet_Py", "std::vector< float >", &m_recojet_Py);
  m_outputTree->Branch("recojet_Pz", "std::vector< float >", &m_recojet_Pz);
  m_outputTree->Branch("recojet_Mult", "std::vector< int >", &m_recojet_Mult);
  m_outputTree->Branch("recojet_NCH", "std::vector< int >", &m_recojet_NCH);
  m_outputTree->Branch("recojet_NCH_trackPtMin", "std::vector< int >", &m_recojet_NCH_trackPtMin);
  m_outputTree->Branch("recojet_NPh", "std::vector< int >", &m_recojet_NPh);
  m_outputTree->Branch("recojet_NNH", "std::vector< int >", &m_recojet_NNH);
  m_outputTree->Branch("recojet_dij_21","std::vector< float >", &m_recojet_dij_21);
  m_outputTree->Branch("recojet_dij_32","std::vector< float >", &m_recojet_dij_32);
  m_outputTree->Branch("recojet_dij_43","std::vector< float >", &m_recojet_dij_43);
  m_outputTree->Branch("recojet_dij_21_max","std::vector< float >", &m_recojet_dij_21_max);
  m_outputTree->Branch("recojet_dij_32_max","std::vector< float >", &m_recojet_dij_32_max);
  m_outputTree->Branch("recojet_dij_43_max","std::vector< float >", &m_recojet_dij_43_max);
  m_outputTree->Branch("recojet_CHFraction", "std::vector< float >", &m_recojet_CHFraction);
  m_outputTree->Branch("recojet_CHFraction_trackPtMin", "std::vector< float >", &m_recojet_CHFraction_trackPtMin);
  m_outputTree->Branch("recojet_PhFraction", "std::vector< float >", &m_recojet_PhFraction);
  m_outputTree->Branch("recojet_ElFraction", "std::vector< float >", &m_recojet_ElFraction);
  m_outputTree->Branch("recojet_MuFraction", "std::vector< float >", &m_recojet_MuFraction);
  m_outputTree->Branch("recojet_NHFraction", "std::vector< float >", &m_recojet_NHFraction);
  m_outputTree->Branch("recojet_nsubjettiness1", "std::vector< float >", m_recojet_nsubjettiness1);
  m_outputTree->Branch("recojet_nsubjettiness2", "std::vector< float >", m_recojet_nsubjettiness2);
  m_outputTree->Branch("recojet_nsubjettiness3", "std::vector< float >", m_recojet_nsubjettiness3); 
  m_outputTree->Branch("recojet_nsubjettiness1_lrz", "std::vector< float >", m_recojet_nsubjettiness1_lrz);
  m_outputTree->Branch("recojet_nsubjettiness2_lrz", "std::vector< float >", m_recojet_nsubjettiness2_lrz);
  m_outputTree->Branch("recojet_nsubjettiness3_lrz", "std::vector< float >", m_recojet_nsubjettiness3_lrz); 

  
  /*
  m_outputTree->Branch("recojet_rfj_4jets_E", "std::vector< float >", &m_recojet_rfj_4jets_E);
  m_outputTree->Branch("recojet_rfj_4jets_Px", "std::vector< float >", &m_recojet_rfj_4jets_Px);
  m_outputTree->Branch("recojet_rfj_4jets_Py", "std::vector< float >", &m_recojet_rfj_4jets_Py);
  m_outputTree->Branch("recojet_rfj_4jets_Pz", "std::vector< float >", &m_recojet_rfj_4jets_Pz);
  m_outputTree->Branch("recojet_rfj_4jets_NCH", "std::vector< int >", &m_recojet_rfj_4jets_NCH);
  m_outputTree->Branch("recojet_rfj_4jets_Mult", "std::vector< int >", &m_recojet_rfj_4jets_Mult);
  m_outputTree->Branch("recojet_rfj_4jets_CHFraction", "std::vector< float >", &m_recojet_rfj_4jets_CHFraction);
  m_outputTree->Branch("recojet_rfj_4jets_PhFraction", "std::vector< float >", &m_recojet_rfj_4jets_PhFraction);
  m_outputTree->Branch("recojet_rfj_4jets_ElFraction", "std::vector< float >", &m_recojet_rfj_4jets_ElFraction);
  m_outputTree->Branch("recojet_rfj_4jets_MuFraction", "std::vector< float >", &m_recojet_rfj_4jets_MuFraction);
  m_outputTree->Branch("recojet_rfj_4jets_NHFraction", "std::vector< float >", &m_recojet_rfj_4jets_NHFraction);
  m_outputTree->Branch("recojet_rfj_4jets_BTag","std::vector< float >", &m_recojet_rfj_4jets_BTag);
  m_outputTree->Branch("recojet_rfj_4jets_CTag","std::vector< float >", &m_recojet_rfj_4jets_CTag);
  m_outputTree->Branch("recojet_rfj_4jets_OTag","std::vector< float >", &m_recojet_rfj_4jets_OTag);
  m_outputTree->Branch("recojet_rfj_4jets_cat","std::vector< int >", &m_recojet_rfj_4jets_cat);

  m_outputTree->Branch("recojet_rfj_4jets_svtx_r","std::vector< float >", &m_recojet_rfj_4jets_svtx_r);
  m_outputTree->Branch("recojet_rfj_4jets_svtx_E","std::vector< float >", &m_recojet_rfj_4jets_svtx_E);
  m_outputTree->Branch("recojet_rfj_4jets_svtx_Mass","std::vector< float >", &m_recojet_rfj_4jets_svtx_Mass);
  m_outputTree->Branch("recojet_rfj_4jets_svtx_nTrack","std::vector< int >", &m_recojet_rfj_4jets_svtx_nTrack);
  m_outputTree->Branch("recojet_rfj_4jets_svtx_Charge","std::vector< int >", &m_recojet_rfj_4jets_svtx_Charge);

  m_outputTree->Branch("recojet_fat2j_rfj4_E", "std::vector< float >", &m_recojet_fat2j_rfj4_E);
  m_outputTree->Branch("recojet_fat2j_rfj4_Px", "std::vector< float >", &m_recojet_fat2j_rfj4_Px);
  m_outputTree->Branch("recojet_fat2j_rfj4_Py", "std::vector< float >", &m_recojet_fat2j_rfj4_Py);
  m_outputTree->Branch("recojet_fat2j_rfj4_Pz", "std::vector< float >", &m_recojet_fat2j_rfj4_Pz);
  m_outputTree->Branch("recojet_fat2j_rfj4_Mult", "std::vector< int >", &m_recojet_fat2j_rfj4_Mult);
  */

  m_outputTree->Branch("recojet_subjet_E", "std::vector< float >", &m_recojet_subjet_E);
  m_outputTree->Branch("recojet_subjet_Px", "std::vector< float >", &m_recojet_subjet_Px);
  m_outputTree->Branch("recojet_subjet_Py", "std::vector< float >", &m_recojet_subjet_Py);
  m_outputTree->Branch("recojet_subjet_Pz", "std::vector< float >", &m_recojet_subjet_Pz);
  m_outputTree->Branch("recojet_subjet_NCH", "std::vector< int >", &m_recojet_subjet_NCH);
  m_outputTree->Branch("recojet_subjet_NCH_trackPtMin", "std::vector< int >", &m_recojet_subjet_NCH_trackPtMin);
  m_outputTree->Branch("recojet_subjet_NPh", "std::vector< int >", &m_recojet_subjet_NPh);
  m_outputTree->Branch("recojet_subjet_NNH", "std::vector< int >", &m_recojet_subjet_NNH);
  m_outputTree->Branch("recojet_subjet_jetindex", "std::vector< int >", &m_recojet_subjet_jetindex);
  m_outputTree->Branch("recojet_subjet_CHFraction", "std::vector< float >", &m_recojet_subjet_CHFraction);
  m_outputTree->Branch("recojet_subjet_CHFraction_trackPtMin", "std::vector< float >", &m_recojet_subjet_CHFraction_trackPtMin);
  m_outputTree->Branch("recojet_subjet_PhFraction", "std::vector< float >", &m_recojet_subjet_PhFraction);
  m_outputTree->Branch("recojet_subjet_ElFraction", "std::vector< float >", &m_recojet_subjet_ElFraction);
  m_outputTree->Branch("recojet_subjet_MuFraction", "std::vector< float >", &m_recojet_subjet_MuFraction);
  m_outputTree->Branch("recojet_subjet_NHFraction", "std::vector< float >", &m_recojet_subjet_NHFraction);

  m_outputTree->Branch("recojet_subjet_jetChargeE_kappa_0_25", "std::vector< float >", &m_recojet_subjet_jetChargeE_kappa_0_25);
  m_outputTree->Branch("recojet_subjet_jetChargeE_kappa_0_50", "std::vector< float >", &m_recojet_subjet_jetChargeE_kappa_0_50);
  //m_outputTree->Branch("recojet_subjet_jetChargeE_kappa_0_75", "std::vector< float >", &m_recojet_subjet_jetChargeE_kappa_0_75);
  //m_outputTree->Branch("recojet_subjet_jetChargeE_kappa_1_00", "std::vector< float >", &m_recojet_subjet_jetChargeE_kappa_1_00);
  //m_outputTree->Branch("recojet_subjet_jetChargeE_kappa_0_10", "std::vector< float >", &m_recojet_subjet_jetChargeE_kappa_0_10);
  //m_outputTree->Branch("recojet_subjet_jetChargeE_kappa_0_15", "std::vector< float >", &m_recojet_subjet_jetChargeE_kappa_0_15);
  m_outputTree->Branch("recojet_subjet_jetChargeE_kappa_0_20", "std::vector< float >", &m_recojet_subjet_jetChargeE_kappa_0_20);
  m_outputTree->Branch("recojet_subjet_jetChargeE_kappa_0_30", "std::vector< float >", &m_recojet_subjet_jetChargeE_kappa_0_30);

  m_outputTree->Branch("recojet_subjet_jetChargePt_kappa_0_25", "std::vector< float >", &m_recojet_subjet_jetChargePt_kappa_0_25);
  m_outputTree->Branch("recojet_subjet_jetChargePt_kappa_0_50", "std::vector< float >", &m_recojet_subjet_jetChargePt_kappa_0_50);
  //m_outputTree->Branch("recojet_subjet_jetChargePt_kappa_0_75", "std::vector< float >", &m_recojet_subjet_jetChargePt_kappa_0_75);
  //m_outputTree->Branch("recojet_subjet_jetChargePt_kappa_1_00", "std::vector< float >", &m_recojet_subjet_jetChargePt_kappa_1_00);
  //m_outputTree->Branch("recojet_subjet_jetChargePt_kappa_0_10", "std::vector< float >", &m_recojet_subjet_jetChargePt_kappa_0_10);
  //m_outputTree->Branch("recojet_subjet_jetChargePt_kappa_0_15", "std::vector< float >", &m_recojet_subjet_jetChargePt_kappa_0_15);
  m_outputTree->Branch("recojet_subjet_jetChargePt_kappa_0_20", "std::vector< float >", &m_recojet_subjet_jetChargePt_kappa_0_20);
  m_outputTree->Branch("recojet_subjet_jetChargePt_kappa_0_30", "std::vector< float >", &m_recojet_subjet_jetChargePt_kappa_0_30);

  //m_outputTree->Branch("recojet_subjet_jetChargePProj_kappa_0_25", "std::vector< float >", &m_recojet_subjet_jetChargePProj_kappa_0_25);
  //m_outputTree->Branch("recojet_subjet_jetChargePProj_kappa_0_50", "std::vector< float >", &m_recojet_subjet_jetChargePProj_kappa_0_50);
  //m_outputTree->Branch("recojet_subjet_jetChargePProj_kappa_0_75", "std::vector< float >", &m_recojet_subjet_jetChargePProj_kappa_0_75);
  //m_outputTree->Branch("recojet_subjet_jetChargePProj_kappa_1_00", "std::vector< float >", &m_recojet_subjet_jetChargePProj_kappa_1_00);
  //m_outputTree->Branch("recojet_subjet_jetChargePProj_kappa_0_10", "std::vector< float >", &m_recojet_subjet_jetChargePProj_kappa_0_10);
  //m_outputTree->Branch("recojet_subjet_jetChargePProj_kappa_0_15", "std::vector< float >", &m_recojet_subjet_jetChargePProj_kappa_0_15);
  //m_outputTree->Branch("recojet_subjet_jetChargePProj_kappa_0_20", "std::vector< float >", &m_recojet_subjet_jetChargePProj_kappa_0_20);
  //m_outputTree->Branch("recojet_subjet_jetChargePProj_kappa_0_30", "std::vector< float >", &m_recojet_subjet_jetChargePProj_kappa_0_30);


  m_outputTree->Branch("recojet_subjet_rfj_j_E","std::vector< float >", &m_recojet_subjet_rfj_j_E);
  m_outputTree->Branch("recojet_subjet_rfj_j_Px","std::vector< float >", &m_recojet_subjet_rfj_j_Px);
  m_outputTree->Branch("recojet_subjet_rfj_j_Py","std::vector< float >", &m_recojet_subjet_rfj_j_Py);
  m_outputTree->Branch("recojet_subjet_rfj_j_Pz","std::vector< float >", &m_recojet_subjet_rfj_j_Pz);
  m_outputTree->Branch("recojet_subjet_rfj_j_NCH","std::vector< int >", &m_recojet_subjet_rfj_j_NCH);
  m_outputTree->Branch("recojet_subjet_rfj_j_NPh","std::vector< int >", &m_recojet_subjet_rfj_j_NPh);
  m_outputTree->Branch("recojet_subjet_rfj_j_NNH","std::vector< int >", &m_recojet_subjet_rfj_j_NNH);
  m_outputTree->Branch("recojet_subjet_rfj_j_jetindex","std::vector< int >", &m_recojet_subjet_rfj_j_jetindex);
  m_outputTree->Branch("recojet_subjet_rfj_j_subjetindex","std::vector< int >", &m_recojet_subjet_rfj_j_subjetindex);
  m_outputTree->Branch("recojet_subjet_rfj_j_CHFraction","std::vector< float >", &m_recojet_subjet_rfj_j_CHFraction);
  m_outputTree->Branch("recojet_subjet_rfj_j_PhFraction","std::vector< float >", &m_recojet_subjet_rfj_j_PhFraction);
  m_outputTree->Branch("recojet_subjet_rfj_j_ElFraction","std::vector< float >", &m_recojet_subjet_rfj_j_ElFraction);
  m_outputTree->Branch("recojet_subjet_rfj_j_MuFraction","std::vector< float >", &m_recojet_subjet_rfj_j_MuFraction);
  m_outputTree->Branch("recojet_subjet_rfj_j_NHFraction","std::vector< float >", &m_recojet_subjet_rfj_j_NHFraction);
  //m_outputTree->Branch("recojet_subjet_rfj_j_jetChargeE_kappa_0_25","std::vector< float >", &m_recojet_subjet_rfj_j_jetChargeE_kappa_0_25);
  //m_outputTree->Branch("recojet_subjet_rfj_j_jetChargeE_kappa_0_50","std::vector< float >", &m_recojet_subjet_rfj_j_jetChargeE_kappa_0_50);
  //m_outputTree->Branch("recojet_subjet_rfj_j_jetChargeE_kappa_0_20","std::vector< float >", &m_recojet_subjet_rfj_j_jetChargeE_kappa_0_20);
  //m_outputTree->Branch("recojet_subjet_rfj_j_jetChargeE_kappa_0_30","std::vector< float >", &m_recojet_subjet_rfj_j_jetChargeE_kappa_0_30);
  
  //m_outputTree->Branch("recojet_subjet_rfj_j_jetChargePt_kappa_0_25","std::vector< float >", &m_recojet_subjet_rfj_j_jetChargePt_kappa_0_25);
  //m_outputTree->Branch("recojet_subjet_rfj_j_jetChargePt_kappa_0_50","std::vector< float >", &m_recojet_subjet_rfj_j_jetChargePt_kappa_0_50);
  //m_outputTree->Branch("recojet_subjet_rfj_j_jetChargePt_kappa_0_20","std::vector< float >", &m_recojet_subjet_rfj_j_jetChargePt_kappa_0_20);
  //m_outputTree->Branch("recojet_subjet_rfj_j_jetChargePt_kappa_0_30","std::vector< float >", &m_recojet_subjet_rfj_j_jetChargePt_kappa_0_30);
  
  m_outputTree->Branch("recojet_subjet_rfj_j_BTag","std::vector< float >", &m_recojet_subjet_rfj_j_BTag);
  m_outputTree->Branch("recojet_subjet_rfj_j_CTag","std::vector< float >", &m_recojet_subjet_rfj_j_CTag);
  m_outputTree->Branch("recojet_subjet_rfj_j_OTag","std::vector< float >", &m_recojet_subjet_rfj_j_OTag);
  m_outputTree->Branch("recojet_subjet_rfj_j_cat","std::vector< int >", &m_recojet_subjet_rfj_j_cat);
  
  m_outputTree->Branch("recojet_subjet_rfj_j_svtx_r","std::vector< float >", &m_recojet_subjet_rfj_j_svtx_r);
  m_outputTree->Branch("recojet_subjet_rfj_j_svtx_E","std::vector< float >", &m_recojet_subjet_rfj_j_svtx_E);
  m_outputTree->Branch("recojet_subjet_rfj_j_svtx_Mass","std::vector< float >", &m_recojet_subjet_rfj_j_svtx_Mass);
  m_outputTree->Branch("recojet_subjet_rfj_j_svtx_nTrack","std::vector< int >", &m_recojet_subjet_rfj_j_svtx_nTrack);
  m_outputTree->Branch("recojet_subjet_rfj_j_svtx_Charge","std::vector< int >", &m_recojet_subjet_rfj_j_svtx_Charge);
  m_outputTree->Branch("recojet_subjet_rfj_j_svtx_jetindex","std::vector< int >", m_recojet_subjet_rfj_j_svtx_jetindex);

  m_outputTree->Branch("recojet_beta1_ECorr2", "std::vector< float >", &m_recojet_beta1_ECorr2);
  m_outputTree->Branch("recojet_beta1_ECorr3", "std::vector< float >", &m_recojet_beta1_ECorr3);
  m_outputTree->Branch("recojet_beta1_N2", "std::vector< float >", &m_recojet_beta1_N2);
  m_outputTree->Branch("recojet_beta1_N3", "std::vector< float >", &m_recojet_beta1_N3);
  m_outputTree->Branch("recojet_beta1_C2", "std::vector< float >", &m_recojet_beta1_C2);
  m_outputTree->Branch("recojet_beta1_C3", "std::vector< float >", &m_recojet_beta1_C3);
  m_outputTree->Branch("recojet_beta1_D2", "std::vector< float >", &m_recojet_beta1_D2);
  m_outputTree->Branch("recojet_beta1_ECorr2_E_theta", "std::vector< float >", &m_recojet_beta1_ECorr2_E_theta);
  m_outputTree->Branch("recojet_beta1_ECorr3_E_theta", "std::vector< float >", &m_recojet_beta1_ECorr3_E_theta);
  m_outputTree->Branch("recojet_beta1_N2_E_theta", "std::vector< float >", &m_recojet_beta1_N2_E_theta);
  m_outputTree->Branch("recojet_beta1_N3_E_theta", "std::vector< float >", &m_recojet_beta1_N3_E_theta);
  m_outputTree->Branch("recojet_beta1_C2_E_theta", "std::vector< float >", &m_recojet_beta1_C2_E_theta);
  m_outputTree->Branch("recojet_beta1_C3_E_theta", "std::vector< float >", &m_recojet_beta1_C3_E_theta);
  m_outputTree->Branch("recojet_beta1_D2_E_theta", "std::vector< float >", &m_recojet_beta1_D2_E_theta);

  m_outputTree->Branch("recojet_beta2_ECorr2", "std::vector< float >", &m_recojet_beta2_ECorr2);
  m_outputTree->Branch("recojet_beta2_ECorr3", "std::vector< float >", &m_recojet_beta2_ECorr3);
  m_outputTree->Branch("recojet_beta2_N2", "std::vector< float >", &m_recojet_beta2_N2);
  m_outputTree->Branch("recojet_beta2_N3", "std::vector< float >", &m_recojet_beta2_N3);
  m_outputTree->Branch("recojet_beta2_C2", "std::vector< float >", &m_recojet_beta2_C2);
  m_outputTree->Branch("recojet_beta2_C3", "std::vector< float >", &m_recojet_beta2_C3);
  m_outputTree->Branch("recojet_beta2_D2", "std::vector< float >", &m_recojet_beta2_D2);
  m_outputTree->Branch("recojet_beta2_ECorr2_E_theta", "std::vector< float >", &m_recojet_beta2_ECorr2_E_theta);
  m_outputTree->Branch("recojet_beta2_ECorr3_E_theta", "std::vector< float >", &m_recojet_beta2_ECorr3_E_theta);
  m_outputTree->Branch("recojet_beta2_N2_E_theta", "std::vector< float >", &m_recojet_beta2_N2_E_theta);
  m_outputTree->Branch("recojet_beta2_N3_E_theta", "std::vector< float >", &m_recojet_beta2_N3_E_theta);
  m_outputTree->Branch("recojet_beta2_C2_E_theta", "std::vector< float >", &m_recojet_beta2_C2_E_theta);
  m_outputTree->Branch("recojet_beta2_C3_E_theta", "std::vector< float >", &m_recojet_beta2_C3_E_theta);
  m_outputTree->Branch("recojet_beta2_D2_E_theta", "std::vector< float >", &m_recojet_beta2_D2_E_theta);

  m_outputTree->Branch("recojet_beta0_5_ECorr2", "std::vector< float >", &m_recojet_beta0_5_ECorr2);
  m_outputTree->Branch("recojet_beta0_5_ECorr3", "std::vector< float >", &m_recojet_beta0_5_ECorr3);
  m_outputTree->Branch("recojet_beta0_5_N2", "std::vector< float >", &m_recojet_beta0_5_N2);
  m_outputTree->Branch("recojet_beta0_5_N3", "std::vector< float >", &m_recojet_beta0_5_N3);
  m_outputTree->Branch("recojet_beta0_5_C2", "std::vector< float >", &m_recojet_beta0_5_C2);
  m_outputTree->Branch("recojet_beta0_5_C3", "std::vector< float >", &m_recojet_beta0_5_C3);
  m_outputTree->Branch("recojet_beta0_5_D2", "std::vector< float >", &m_recojet_beta0_5_D2);
  m_outputTree->Branch("recojet_beta0_5_ECorr2_E_theta", "std::vector< float >", &m_recojet_beta0_5_ECorr2_E_theta);
  m_outputTree->Branch("recojet_beta0_5_ECorr3_E_theta", "std::vector< float >", &m_recojet_beta0_5_ECorr3_E_theta);
  m_outputTree->Branch("recojet_beta0_5_N2_E_theta", "std::vector< float >", &m_recojet_beta0_5_N2_E_theta);
  m_outputTree->Branch("recojet_beta0_5_N3_E_theta", "std::vector< float >", &m_recojet_beta0_5_N3_E_theta);
  m_outputTree->Branch("recojet_beta0_5_C2_E_theta", "std::vector< float >", &m_recojet_beta0_5_C2_E_theta);
  m_outputTree->Branch("recojet_beta0_5_C3_E_theta", "std::vector< float >", &m_recojet_beta0_5_C3_E_theta);
  m_outputTree->Branch("recojet_beta0_5_D2_E_theta", "std::vector< float >", &m_recojet_beta0_5_D2_E_theta);
  
  m_outputTree->Branch("reco_tau_E", "std::vector< float >", &m_reco_tau_E);
  m_outputTree->Branch("reco_tau_Px", "std::vector< float >", &m_reco_tau_Px);
  m_outputTree->Branch("reco_tau_Py", "std::vector< float >", &m_reco_tau_Py);
  m_outputTree->Branch("reco_tau_Pz", "std::vector< float >", &m_reco_tau_Pz);
  m_outputTree->Branch("reco_tau_Mult", "std::vector< int >", &m_reco_tau_Mult);
  m_outputTree->Branch("reco_tau_Charge", "std::vector< int >", &m_reco_tau_Charge);
  m_outputTree->Branch("reco_tau_NCH", "std::vector< int >", &m_reco_tau_NCH);
  m_outputTree->Branch("reco_tau_NPh", "std::vector< int >", &m_reco_tau_NPh);
  m_outputTree->Branch("reco_tau_NNH", "std::vector< int >", &m_reco_tau_NNH);
  m_outputTree->Branch("reco_tau_CHFraction", "std::vector< float >", &m_reco_tau_CHFraction);
  m_outputTree->Branch("reco_tau_PhFraction", "std::vector< float >", &m_reco_tau_PhFraction);
  m_outputTree->Branch("reco_tau_ElFraction", "std::vector< float >", &m_reco_tau_ElFraction);
  m_outputTree->Branch("reco_tau_MuFraction", "std::vector< float >", &m_reco_tau_MuFraction);
  m_outputTree->Branch("reco_tau_NHFraction", "std::vector< float >", &m_reco_tau_NHFraction);

  m_outputTree->Branch("gen_tau_E", "std::vector< float >", &m_gen_tau_E);
  m_outputTree->Branch("gen_tau_Px", "std::vector< float >", &m_gen_tau_Px);
  m_outputTree->Branch("gen_tau_Py", "std::vector< float >", &m_gen_tau_Py);
  m_outputTree->Branch("gen_tau_Pz", "std::vector< float >", &m_gen_tau_Pz);
  m_outputTree->Branch("gen_tau_Mult", "std::vector< int >", &m_gen_tau_Mult);
  m_outputTree->Branch("gen_tau_Charge", "std::vector< int >", &m_gen_tau_Charge);
  m_outputTree->Branch("gen_tau_MCCharge", "std::vector< int >", &m_gen_tau_MCCharge);
  m_outputTree->Branch("gen_tau_NCH", "std::vector< int >", &m_gen_tau_NCH);
  m_outputTree->Branch("gen_tau_NPh", "std::vector< int >", &m_gen_tau_NPh);
  m_outputTree->Branch("gen_tau_NNH", "std::vector< int >", &m_gen_tau_NNH);
  m_outputTree->Branch("gen_tau_MotherPDGID", "std::vector< int >", &m_gen_tau_MotherPDGID);
  m_outputTree->Branch("gen_tau_CHFraction", "std::vector< float >", &m_gen_tau_CHFraction);
  m_outputTree->Branch("gen_tau_PhFraction", "std::vector< float >", &m_gen_tau_PhFraction);
  m_outputTree->Branch("gen_tau_ElFraction", "std::vector< float >", &m_gen_tau_ElFraction);
  m_outputTree->Branch("gen_tau_MuFraction", "std::vector< float >", &m_gen_tau_MuFraction);
  m_outputTree->Branch("gen_tau_NHFraction", "std::vector< float >", &m_gen_tau_NHFraction);

  //fastjet::contrib::ValenciaPlugin* vlcpl10;
  //R,beta,gamma
  //-->beta is weight of energy in distance measures, gamma gives weight to 
  vlcpl = new fastjet::contrib::ValenciaPlugin(m_R,m_beta,m_gamma);
  //					      atof(_jetAlgoNameAndParams[1].c_str()),  // R value
					      //atof(_jetAlgoNameAndParams[2].c_str()),  // beta value
					      //atof(_jetAlgoNameAndParams[3].c_str())   // gamma value
					      //);
  //fastjet::JetDefinition* 
  _jetAlgoType = new fastjet::JetDefinition(vlcpl);
 
  _vlcAxes = new fastjet::contrib::ExclusiveJetAxes((*_jetAlgoType));

  _unnormalizedMeasure_ptR = new fastjet::contrib::UnnormalizedMeasure(m_beta,fastjet::contrib::pt_R);
  _unnormalizedMeasure_Lorentz = new fastjet::contrib::UnnormalizedMeasure(m_beta,fastjet::contrib::perp_lorentz_dot);

  _nSubJettiness1_ptR = new fastjet::contrib::Nsubjettiness(1, (*_vlcAxes),(*_unnormalizedMeasure_ptR));
  _nSubJettiness1_lorentz= new fastjet::contrib::Nsubjettiness(1, (*_vlcAxes),(*_unnormalizedMeasure_Lorentz));
  _nSubJettiness2_ptR = new fastjet::contrib::Nsubjettiness(2, (*_vlcAxes),(*_unnormalizedMeasure_ptR));
  _nSubJettiness2_lorentz= new fastjet::contrib::Nsubjettiness(2, (*_vlcAxes),(*_unnormalizedMeasure_Lorentz));
  _nSubJettiness3_ptR = new fastjet::contrib::Nsubjettiness(3, (*_vlcAxes),(*_unnormalizedMeasure_ptR));
  _nSubJettiness3_lorentz= new fastjet::contrib::Nsubjettiness(3, (*_vlcAxes),(*_unnormalizedMeasure_Lorentz));


  //here beta is another weight, determines hw to weight the distance measure angles \theta_ij in correlations
  float m_beta_corr=1.0;
  _energycorr2_beta1= new fastjet::contrib::EnergyCorrelator(2,m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorr2_Etheta_beta1= new fastjet::contrib::EnergyCorrelator(2,m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorr3_beta1= new fastjet::contrib::EnergyCorrelator(3,m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorr3_Etheta_beta1= new fastjet::contrib::EnergyCorrelator(3,m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorr4_beta1= new fastjet::contrib::EnergyCorrelator(4,m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorr4_Etheta_beta1= new fastjet::contrib::EnergyCorrelator(4,m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  
  _energycorrC2_beta1= new fastjet::contrib::EnergyCorrelatorC2(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrC2_Etheta_beta1= new fastjet::contrib::EnergyCorrelatorC2(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorrD2_beta1= new fastjet::contrib::EnergyCorrelatorD2(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrD2_Etheta_beta1= new fastjet::contrib::EnergyCorrelatorD2(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorrN2_beta1= new fastjet::contrib::EnergyCorrelatorN2(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrN2_Etheta_beta1= new fastjet::contrib::EnergyCorrelatorN2(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorrN3_beta1= new fastjet::contrib::EnergyCorrelatorN3(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrN3_Etheta_beta1= new fastjet::contrib::EnergyCorrelatorN3(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);

  m_beta_corr=2.0;
  _energycorr2_beta2= new fastjet::contrib::EnergyCorrelator(2,m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorr2_Etheta_beta2= new fastjet::contrib::EnergyCorrelator(2,m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorr3_beta2= new fastjet::contrib::EnergyCorrelator(3,m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorr3_Etheta_beta2= new fastjet::contrib::EnergyCorrelator(3,m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorr4_beta2= new fastjet::contrib::EnergyCorrelator(4,m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorr4_Etheta_beta2= new fastjet::contrib::EnergyCorrelator(4,m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  
  _energycorrC2_beta2= new fastjet::contrib::EnergyCorrelatorC2(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrC2_Etheta_beta2= new fastjet::contrib::EnergyCorrelatorC2(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorrD2_beta2= new fastjet::contrib::EnergyCorrelatorD2(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrD2_Etheta_beta2= new fastjet::contrib::EnergyCorrelatorD2(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorrN2_beta2= new fastjet::contrib::EnergyCorrelatorN2(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrN2_Etheta_beta2= new fastjet::contrib::EnergyCorrelatorN2(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorrN3_beta2= new fastjet::contrib::EnergyCorrelatorN3(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrN3_Etheta_beta2= new fastjet::contrib::EnergyCorrelatorN3(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);

  m_beta_corr=0.5;
  _energycorr2_beta0_5= new fastjet::contrib::EnergyCorrelator(2,m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorr2_Etheta_beta0_5= new fastjet::contrib::EnergyCorrelator(2,m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorr3_beta0_5= new fastjet::contrib::EnergyCorrelator(3,m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorr3_Etheta_beta0_5= new fastjet::contrib::EnergyCorrelator(3,m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorr4_beta0_5= new fastjet::contrib::EnergyCorrelator(4,m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorr4_Etheta_beta0_5= new fastjet::contrib::EnergyCorrelator(4,m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  
  _energycorrC2_beta0_5= new fastjet::contrib::EnergyCorrelatorC2(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrC2_Etheta_beta0_5= new fastjet::contrib::EnergyCorrelatorC2(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorrD2_beta0_5= new fastjet::contrib::EnergyCorrelatorD2(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrD2_Etheta_beta0_5= new fastjet::contrib::EnergyCorrelatorD2(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorrN2_beta0_5= new fastjet::contrib::EnergyCorrelatorN2(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrN2_Etheta_beta0_5= new fastjet::contrib::EnergyCorrelatorN2(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  _energycorrN3_beta0_5= new fastjet::contrib::EnergyCorrelatorN3(m_beta_corr, fastjet::contrib::EnergyCorrelator::pt_R);
  _energycorrN3_Etheta_beta0_5= new fastjet::contrib::EnergyCorrelatorN3(m_beta_corr, fastjet::contrib::EnergyCorrelator::E_theta);
  
}


void HZAnalyzer::processRunHeader( LCRunHeader*) {

}

void HZAnalyzer::processEvent( LCEvent* evt ) { 
  m_eventCount+=1;

  if(m_eventCount==1){
    std::cout<<"params "<<m_angleIso<<"/"<<m_genTrueLepPhEMin<<"/"<<m_R<<"/"<<m_cutPtTrackMin<<std::endl;
  }
  //if(evt->getEventNumber()%50==0){
  //std::cout<<"run/evt "<<evt->getRunNumber()<<"/"<<evt->getEventNumber()<<"/"<<m_eventCount<<std::endl;
    //}

  m_runNumber = 0 ;
  m_eventNumber = 0 ;

  m_true_E=0;
  m_true_Px=0;
  m_true_Py=0;
  m_true_Pz=0;
  m_true_Mult=0;
  m_true_NCH=0;
  m_true_NPh=0;
  m_true_NNH=0;
  m_true_CHFraction=0;
  m_true_PhFraction=0;
  m_true_ElFraction=0;
  m_true_MuFraction=0;
  m_true_NHFraction=0;
  
  m_true_inv_E=0;
  m_true_inv_Px=0;
  m_true_inv_Py=0;
  m_true_inv_Pz=0;
  m_true_inv_Mult=0;

  m_totPFO_E=0;
  m_totPFO_Px=0;
  m_totPFO_Py=0;
  m_totPFO_Pz=0;
  m_totPFO_Mult=0;
  m_totPFO_NCH=0;
  m_totPFO_NPh=0;
  m_totPFO_NNH=0;
  m_totPFO_CHFraction=0;
  m_totPFO_PhFraction=0;
  m_totPFO_ElFraction=0;
  m_totPFO_MuFraction=0;
  m_totPFO_NHFraction=0;

  m_gen_y21_max=0;
  m_gen_y32_max=0;
  m_gen_y43_max=0;

  m_gen_y21=0;
  m_gen_y32=0;
  m_gen_y43=0;

  m_reco_y21= 0;
  m_reco_y32= 0;
  m_reco_y43= 0;
  m_reco_y21_max= 0;
  m_reco_y32_max= 0;
  m_reco_y43_max= 0;

  m_trueME_E->clear();
  m_trueME_Px->clear();
  m_trueME_Py->clear();
  m_trueME_Pz->clear();
  m_trueME_PDGID->clear();
  //m_trueME_index->clear();
  
  m_genTrueLepPh_E->clear();
  m_genTrueLepPh_Px->clear();
  m_genTrueLepPh_Py->clear();
  m_genTrueLepPh_Pz->clear();
  m_genTrueLepPh_PDGID->clear();
  //actual ancestor, i.e. Higgs,W,Z or e+e-
  m_genTrueLepPh_ANC_PDGID->clear();
  m_genTrueLepPh_relIso->clear();
  m_genTrueLepPh_relIsoCH->clear();
  m_genTrueLepPh_relIsoPh->clear();
  m_genTrueLepPh_relIsoNH->clear();
  m_genTrueLepPh_relIsoEl->clear();
  m_genTrueLepPh_relIsoMu->clear();
 
  //m_genjet_jetChargeE_kappa_0_25->clear();
  //m_genjet_jetChargeE_kappa_0_50->clear();
  //m_genjet_jetChargeE_kappa_0_75->clear();
  //m_genjet_jetChargeE_kappa_1_00->clear();
  //m_genjet_jetChargeE_kappa_0_10->clear();
  //m_genjet_jetChargeE_kappa_0_15->clear();
  //m_genjet_jetChargeE_kappa_0_20->clear();
  //m_genjet_jetChargeE_kappa_0_30->clear();
  
  //m_genjet_jetChargePt_kappa_0_25->clear();
  //m_genjet_jetChargePt_kappa_0_50->clear();
  //m_genjet_jetChargePt_kappa_0_75->clear();
  //m_genjet_jetChargePt_kappa_1_00->clear();
  //m_genjet_jetChargePt_kappa_0_10->clear();
  //m_genjet_jetChargePt_kappa_0_15->clear();
  //m_genjet_jetChargePt_kappa_0_20->clear();
  //m_genjet_jetChargePt_kappa_0_30->clear();

  //m_genjet_jetChargePProj_kappa_0_25->clear();
  //m_genjet_jetChargePProj_kappa_0_50->clear();
  //m_genjet_jetChargePProj_kappa_0_75->clear();
  //m_genjet_jetChargePProj_kappa_1_00->clear();
  //m_genjet_jetChargePProj_kappa_0_10->clear();
  //m_genjet_jetChargePProj_kappa_0_15->clear();
  //m_genjet_jetChargePProj_kappa_0_20->clear();
  //m_genjet_jetChargePProj_kappa_0_30->clear();
  
  m_genjet_E->clear();
  m_genjet_Px->clear();
  m_genjet_Py->clear();
  m_genjet_Pz->clear();
  m_genjet_Mult->clear();
  m_genjet_NCH->clear();
  m_genjet_NCH_trackPtMin->clear();
  m_genjet_NPh->clear();
  m_genjet_NNH->clear();
  m_genjet_dij_21->clear();
  m_genjet_dij_32->clear();
  m_genjet_dij_43->clear();
  m_genjet_dij_21_max->clear();
  m_genjet_dij_32_max->clear();
  m_genjet_dij_43_max->clear();
  m_genjet_CHFraction->clear();
  m_genjet_CHFraction_trackPtMin->clear();
  m_genjet_PhFraction->clear();
  m_genjet_ElFraction->clear();
  m_genjet_MuFraction->clear();
  m_genjet_NHFraction->clear();
  m_genjet_nsubjettiness1->clear();
  m_genjet_nsubjettiness2->clear();
  m_genjet_nsubjettiness3->clear();
  m_genjet_nsubjettiness1_lrz->clear();
  m_genjet_nsubjettiness2_lrz->clear();
  m_genjet_nsubjettiness3_lrz->clear();

  m_genjet_beta1_ECorr2->clear();
  m_genjet_beta1_ECorr3->clear();
  m_genjet_beta1_N2->clear();
  m_genjet_beta1_N3->clear();
  m_genjet_beta1_C2->clear();
  m_genjet_beta1_C3->clear();
  m_genjet_beta1_D2->clear();
  m_genjet_beta1_ECorr2_E_theta->clear();
  m_genjet_beta1_ECorr3_E_theta->clear();
  m_genjet_beta1_N2_E_theta->clear();
  m_genjet_beta1_N3_E_theta->clear();
  m_genjet_beta1_C2_E_theta->clear();
  m_genjet_beta1_C3_E_theta->clear();
  m_genjet_beta1_D2_E_theta->clear();

  m_genjet_beta2_ECorr2->clear();
  m_genjet_beta2_ECorr3->clear();
  m_genjet_beta2_N2->clear();
  m_genjet_beta2_N3->clear();
  m_genjet_beta2_C2->clear();
  m_genjet_beta2_C3->clear();
  m_genjet_beta2_D2->clear();
  m_genjet_beta2_ECorr2_E_theta->clear();
  m_genjet_beta2_ECorr3_E_theta->clear();
  m_genjet_beta2_N2_E_theta->clear();
  m_genjet_beta2_N3_E_theta->clear();
  m_genjet_beta2_C2_E_theta->clear();
  m_genjet_beta2_C3_E_theta->clear();
  m_genjet_beta2_D2_E_theta->clear();

  m_genjet_beta0_5_ECorr2->clear();
  m_genjet_beta0_5_ECorr3->clear();
  m_genjet_beta0_5_N2->clear();
  m_genjet_beta0_5_N3->clear();
  m_genjet_beta0_5_C2->clear();
  m_genjet_beta0_5_C3->clear();
  m_genjet_beta0_5_D2->clear();
  m_genjet_beta0_5_ECorr2_E_theta->clear();
  m_genjet_beta0_5_ECorr3_E_theta->clear();
  m_genjet_beta0_5_N2_E_theta->clear();
  m_genjet_beta0_5_N3_E_theta->clear();
  m_genjet_beta0_5_C2_E_theta->clear();
  m_genjet_beta0_5_C3_E_theta->clear();
  m_genjet_beta0_5_D2_E_theta->clear();

  
  m_genjet_subjet_E->clear();
  m_genjet_subjet_Px->clear();
  m_genjet_subjet_Py->clear();
  m_genjet_subjet_Pz->clear();
  m_genjet_subjet_NCH->clear();
  m_genjet_subjet_NCH_trackPtMin->clear();
  m_genjet_subjet_NPh->clear();
  m_genjet_subjet_NNH->clear();
  m_genjet_subjet_jetindex->clear();
  m_genjet_subjet_CHFraction->clear();
  m_genjet_subjet_CHFraction_trackPtMin->clear();
  m_genjet_subjet_PhFraction->clear();
  m_genjet_subjet_ElFraction->clear();
  m_genjet_subjet_MuFraction->clear();
  m_genjet_subjet_NHFraction->clear();

  m_genjet_subjet_jetChargeE_kappa_0_25->clear();
  m_genjet_subjet_jetChargeE_kappa_0_50->clear();
  //m_genjet_subjet_jetChargeE_kappa_0_75->clear();
  //m_genjet_subjet_jetChargeE_kappa_1_00->clear();
  //m_genjet_subjet_jetChargeE_kappa_0_10->clear();
  //m_genjet_subjet_jetChargeE_kappa_0_15->clear();
  m_genjet_subjet_jetChargeE_kappa_0_20->clear();
  m_genjet_subjet_jetChargeE_kappa_0_30->clear();
  
  m_genjet_subjet_jetChargePt_kappa_0_25->clear();
  m_genjet_subjet_jetChargePt_kappa_0_50->clear();
  //m_genjet_subjet_jetChargePt_kappa_0_75->clear();
  //m_genjet_subjet_jetChargePt_kappa_1_00->clear();
  //m_genjet_subjet_jetChargePt_kappa_0_10->clear();
  //m_genjet_subjet_jetChargePt_kappa_0_15->clear();
  m_genjet_subjet_jetChargePt_kappa_0_20->clear();
  m_genjet_subjet_jetChargePt_kappa_0_30->clear();

  //m_genjet_subjet_jetChargePProj_kappa_0_25->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_50->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_75->clear();
  //m_genjet_subjet_jetChargePProj_kappa_1_00->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_10->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_15->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_20->clear();
  //m_genjet_subjet_jetChargePProj_kappa_0_30->clear();
  
  m_isoPartGenDR10_E->clear();
  m_isoPartGenDR10_Px->clear();
  m_isoPartGenDR10_Py->clear();
  m_isoPartGenDR10_Pz->clear();
  m_isoPartGenDR10_PDGID->clear();
  m_isoPartGenDR10_relIso->clear();
  m_isoPartGenDR10_relIsoCH->clear();
  m_isoPartGenDR10_relIsoPh->clear();
  m_isoPartGenDR10_relIsoNH->clear();
  m_isoPartGenDR10_relIsoEl->clear();
  m_isoPartGenDR10_relIsoMu->clear();
  
  m_isoPartRecoDR10_E->clear();
  m_isoPartRecoDR10_Px->clear();
  m_isoPartRecoDR10_Py->clear();
  m_isoPartRecoDR10_Pz->clear();
  m_isoPartRecoDR10_PDGID->clear();
  m_isoPartRecoDR10_relIso->clear();
  m_isoPartRecoDR10_relIsoCH->clear();
  m_isoPartRecoDR10_relIsoPh->clear();
  m_isoPartRecoDR10_relIsoNH->clear();
  m_isoPartRecoDR10_relIsoEl->clear();
  m_isoPartRecoDR10_relIsoMu->clear();
  
  //m_recojet_jetChargeE_kappa_0_25->clear();
  //m_recojet_jetChargeE_kappa_0_50->clear();
  //m_recojet_jetChargeE_kappa_0_75->clear();
  //m_recojet_jetChargeE_kappa_1_00->clear();
  //m_recojet_jetChargeE_kappa_0_10->clear();
  //m_recojet_jetChargeE_kappa_0_15->clear();
  //m_recojet_jetChargeE_kappa_0_20->clear();
  //m_recojet_jetChargeE_kappa_0_30->clear();
  
  //m_recojet_jetChargePt_kappa_0_25->clear();
  //m_recojet_jetChargePt_kappa_0_50->clear();
  //m_recojet_jetChargePt_kappa_0_75->clear();
  //m_recojet_jetChargePt_kappa_1_00->clear();
  //m_recojet_jetChargePt_kappa_0_10->clear();
  //m_recojet_jetChargePt_kappa_0_15->clear();
  //m_recojet_jetChargePt_kappa_0_20->clear();
  //m_recojet_jetChargePt_kappa_0_30->clear();

  //m_recojet_jetChargePProj_kappa_0_25->clear();
  //m_recojet_jetChargePProj_kappa_0_50->clear();
  //m_recojet_jetChargePProj_kappa_0_75->clear();
  //m_recojet_jetChargePProj_kappa_1_00->clear();
  //m_recojet_jetChargePProj_kappa_0_10->clear();
  //m_recojet_jetChargePProj_kappa_0_15->clear();
  //m_recojet_jetChargePProj_kappa_0_20->clear();
  //m_recojet_jetChargePProj_kappa_0_30->clear();

  m_recojet_E->clear();
  m_recojet_Px->clear();
  m_recojet_Py->clear();
  m_recojet_Pz->clear();
  m_recojet_Mult->clear();
  m_recojet_NCH->clear();
  m_recojet_NCH_trackPtMin->clear();
  m_recojet_NPh->clear();
  m_recojet_NNH->clear();
  m_recojet_dij_21->clear();
  m_recojet_dij_32->clear();
  m_recojet_dij_43->clear();
  m_recojet_dij_21_max->clear();
  m_recojet_dij_32_max->clear();
  m_recojet_dij_43_max->clear();
  m_recojet_CHFraction->clear();
  m_recojet_CHFraction_trackPtMin->clear();
  m_recojet_PhFraction->clear();
  m_recojet_ElFraction->clear();
  m_recojet_MuFraction->clear();
  m_recojet_NHFraction->clear();
  m_recojet_nsubjettiness1->clear();
  m_recojet_nsubjettiness2->clear();
  m_recojet_nsubjettiness3->clear();
  m_recojet_nsubjettiness1_lrz->clear();
  m_recojet_nsubjettiness2_lrz->clear();
  m_recojet_nsubjettiness3_lrz->clear();

  m_recojet_beta1_ECorr2->clear();
  m_recojet_beta1_ECorr3->clear();
  m_recojet_beta1_N2->clear();
  m_recojet_beta1_N3->clear();
  m_recojet_beta1_C2->clear();
  m_recojet_beta1_C3->clear();
  m_recojet_beta1_D2->clear();
  m_recojet_beta1_ECorr2_E_theta->clear();
  m_recojet_beta1_ECorr3_E_theta->clear();
  m_recojet_beta1_N2_E_theta->clear();
  m_recojet_beta1_N3_E_theta->clear();
  m_recojet_beta1_C2_E_theta->clear();
  m_recojet_beta1_C3_E_theta->clear();
  m_recojet_beta1_D2_E_theta->clear();

  m_recojet_beta2_ECorr2->clear();
  m_recojet_beta2_ECorr3->clear();
  m_recojet_beta2_N2->clear();
  m_recojet_beta2_N3->clear();
  m_recojet_beta2_C2->clear();
  m_recojet_beta2_C3->clear();
  m_recojet_beta2_D2->clear();
  m_recojet_beta2_ECorr2_E_theta->clear();
  m_recojet_beta2_ECorr3_E_theta->clear();
  m_recojet_beta2_N2_E_theta->clear();
  m_recojet_beta2_N3_E_theta->clear();
  m_recojet_beta2_C2_E_theta->clear();
  m_recojet_beta2_C3_E_theta->clear();
  m_recojet_beta2_D2_E_theta->clear();

  m_recojet_beta0_5_ECorr2->clear();
  m_recojet_beta0_5_ECorr3->clear();
  m_recojet_beta0_5_N2->clear();
  m_recojet_beta0_5_N3->clear();
  m_recojet_beta0_5_C2->clear();
  m_recojet_beta0_5_C3->clear();
  m_recojet_beta0_5_D2->clear();
  m_recojet_beta0_5_ECorr2_E_theta->clear();
  m_recojet_beta0_5_ECorr3_E_theta->clear();
  m_recojet_beta0_5_N2_E_theta->clear();
  m_recojet_beta0_5_N3_E_theta->clear();
  m_recojet_beta0_5_C2_E_theta->clear();
  m_recojet_beta0_5_C3_E_theta->clear();
  m_recojet_beta0_5_D2_E_theta->clear();

  /*
  m_recojet_rfj_4jets_E->clear();
  m_recojet_rfj_4jets_Px->clear();
  m_recojet_rfj_4jets_Py->clear();
  m_recojet_rfj_4jets_Pz->clear();
  m_recojet_rfj_4jets_Mult->clear();
  m_recojet_rfj_4jets_NCH->clear();
  m_recojet_rfj_4jets_CHFraction->clear();
  m_recojet_rfj_4jets_PhFraction->clear();
  m_recojet_rfj_4jets_ElFraction->clear();
  m_recojet_rfj_4jets_MuFraction->clear();
  m_recojet_rfj_4jets_NHFraction->clear();

  m_recojet_rfj_4jets_BTag->clear();
  m_recojet_rfj_4jets_CTag->clear();
  m_recojet_rfj_4jets_OTag->clear();
  m_recojet_rfj_4jets_cat->clear();

  m_recojet_rfj_4jets_svtx_r->clear();
  m_recojet_rfj_4jets_svtx_E->clear();
  m_recojet_rfj_4jets_svtx_Mass->clear();
  m_recojet_rfj_4jets_svtx_nTrack->clear();
  m_recojet_rfj_4jets_svtx_Charge->clear();

  m_recojet_fat2j_rfj4_E->clear();
  m_recojet_fat2j_rfj4_Px->clear();
  m_recojet_fat2j_rfj4_Py->clear();
  m_recojet_fat2j_rfj4_Pz->clear();
  m_recojet_fat2j_rfj4_Mult->clear();
  */
  m_recojet_subjet_E->clear();
  m_recojet_subjet_Px->clear();
  m_recojet_subjet_Py->clear();
  m_recojet_subjet_Pz->clear();
  m_recojet_subjet_NCH->clear();
  m_recojet_subjet_NCH_trackPtMin->clear();
  m_recojet_subjet_NPh->clear();
  m_recojet_subjet_NNH->clear();
  m_recojet_subjet_jetindex->clear();
  m_recojet_subjet_CHFraction->clear();
  m_recojet_subjet_CHFraction_trackPtMin->clear();
  m_recojet_subjet_PhFraction->clear();
  m_recojet_subjet_ElFraction->clear();
  m_recojet_subjet_MuFraction->clear();
  m_recojet_subjet_NHFraction->clear();

  m_recojet_subjet_jetChargeE_kappa_0_25->clear();
  m_recojet_subjet_jetChargeE_kappa_0_50->clear();
  //m_recojet_subjet_jetChargeE_kappa_0_75->clear();
  //m_recojet_subjet_jetChargeE_kappa_1_00->clear();
  //m_recojet_subjet_jetChargeE_kappa_0_10->clear();
  //m_recojet_subjet_jetChargeE_kappa_0_15->clear();
  m_recojet_subjet_jetChargeE_kappa_0_20->clear();
  m_recojet_subjet_jetChargeE_kappa_0_30->clear();
  
  m_recojet_subjet_jetChargePt_kappa_0_25->clear();
  m_recojet_subjet_jetChargePt_kappa_0_50->clear();
  //m_recojet_subjet_jetChargePt_kappa_0_75->clear();
  //m_recojet_subjet_jetChargePt_kappa_1_00->clear();
  //m_recojet_subjet_jetChargePt_kappa_0_10->clear();
  //m_recojet_subjet_jetChargePt_kappa_0_15->clear();
  m_recojet_subjet_jetChargePt_kappa_0_20->clear();
  m_recojet_subjet_jetChargePt_kappa_0_30->clear();

  //m_recojet_subjet_jetChargePProj_kappa_0_25->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_50->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_75->clear();
  //m_recojet_subjet_jetChargePProj_kappa_1_00->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_10->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_15->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_20->clear();
  //m_recojet_subjet_jetChargePProj_kappa_0_30->clear();
   
  
  m_recojet_subjet_rfj_j_E->clear();
  m_recojet_subjet_rfj_j_Px->clear();
  m_recojet_subjet_rfj_j_Py->clear();
  m_recojet_subjet_rfj_j_Pz->clear();
  m_recojet_subjet_rfj_j_NCH->clear();
  m_recojet_subjet_rfj_j_NPh->clear();
  m_recojet_subjet_rfj_j_NNH->clear();
  m_recojet_subjet_rfj_j_jetindex->clear();
  m_recojet_subjet_rfj_j_subjetindex->clear();
  m_recojet_subjet_rfj_j_CHFraction->clear();
  m_recojet_subjet_rfj_j_PhFraction->clear();
  m_recojet_subjet_rfj_j_ElFraction->clear();
  m_recojet_subjet_rfj_j_MuFraction->clear();
  m_recojet_subjet_rfj_j_NHFraction->clear();
  //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_25->clear();
  //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_50->clear();
  //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_20->clear();
  //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_30->clear();
  
  //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_25->clear();
  //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_50->clear();
  //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_20->clear();
  //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_30->clear();
  
  m_recojet_subjet_rfj_j_BTag->clear();
  m_recojet_subjet_rfj_j_CTag->clear();
  m_recojet_subjet_rfj_j_OTag->clear();
  m_recojet_subjet_rfj_j_cat->clear();
  
  m_recojet_subjet_rfj_j_svtx_r->clear();
  m_recojet_subjet_rfj_j_svtx_E->clear();
  m_recojet_subjet_rfj_j_svtx_Mass->clear();
  m_recojet_subjet_rfj_j_svtx_nTrack->clear();
  m_recojet_subjet_rfj_j_svtx_Charge->clear();
  m_recojet_subjet_rfj_j_svtx_jetindex->clear();
 
  m_reco_tau_E->clear();
  m_reco_tau_Px->clear();
  m_reco_tau_Py->clear();
  m_reco_tau_Pz->clear();
  m_reco_tau_Mult->clear();
  m_reco_tau_Charge->clear();
  m_reco_tau_NCH->clear();
  m_reco_tau_NPh->clear();
  m_reco_tau_NNH->clear();
  m_reco_tau_CHFraction->clear();
  m_reco_tau_PhFraction->clear();
  m_reco_tau_ElFraction->clear();
  m_reco_tau_MuFraction->clear();
  m_reco_tau_NHFraction->clear();

  m_gen_tau_E->clear();
  m_gen_tau_Px->clear();
  m_gen_tau_Py->clear();
  m_gen_tau_Pz->clear();
  m_gen_tau_Mult->clear();
  m_gen_tau_Charge->clear();
  m_gen_tau_MCCharge->clear();
  m_gen_tau_NCH->clear();
  m_gen_tau_NPh->clear();
  m_gen_tau_NNH->clear();
  m_gen_tau_MotherPDGID->clear();
  m_gen_tau_CHFraction->clear();
  m_gen_tau_PhFraction->clear();
  m_gen_tau_ElFraction->clear();
  m_gen_tau_MuFraction->clear();
  m_gen_tau_NHFraction->clear();

  
  m_eventNumber=evt->getEventNumber();
  m_runNumber=evt->getRunNumber();


  m_totPFO_E=0;
  m_totPFO_Px=0;
  m_totPFO_Py=0;
  m_totPFO_Pz=0;
  m_totPFO_Mult=0;
  m_totPFO_NCH=0;
  m_totPFO_NPh=0;
  m_totPFO_NNH=0;
  m_totPFO_CHFraction=0;
  m_totPFO_PhFraction=0;
  m_totPFO_ElFraction=0;
  m_totPFO_MuFraction=0;
  m_totPFO_NHFraction=0;
  
  m_true_E=0;
  m_true_Px=0;
  m_true_Py=0;
  m_true_Pz=0;
  m_true_Mult=0;
  m_true_NCH=0;
  m_true_NPh=0;
  m_true_NNH=0;
  m_true_CHFraction=0;
  m_true_PhFraction=0;
  m_true_ElFraction=0;
  m_true_MuFraction=0;
  m_true_NHFraction=0;
  
  m_true_inv_E=0;
  m_true_inv_Px=0;
  m_true_inv_Py=0;
  m_true_inv_Pz=0;
  m_true_inv_Mult=0;

  //bool rfjet_has_duplicate=false;  

  LCCollection * mcColl =0;
  getCollection(mcColl,m_inputMCParticleCollection,evt);
  LCCollection * genIsoColl =0;
  getCollection(genIsoColl,m_genIsoPartColName,evt);
  
  if(genIsoColl != 0 ){
    for(int m=0;m<genIsoColl->getNumberOfElements();m++){
      MCParticle *genisop = static_cast<MCParticle*>(genIsoColl->getElementAt(m));
      m_isoPartGenDR10_E->push_back(genisop->getEnergy());
      m_isoPartGenDR10_Px->push_back(genisop->getMomentum()[0]);
      m_isoPartGenDR10_Py->push_back(genisop->getMomentum()[1]);
      m_isoPartGenDR10_Pz->push_back(genisop->getMomentum()[2]);
      m_isoPartGenDR10_PDGID->push_back(genisop->getPDG());
      TLorentzVector geniso(0,0,0,0);
      geniso.SetPxPyPzE(genisop->getMomentum()[0],genisop->getMomentum()[1],genisop->getMomentum()[2],genisop->getEnergy());
      float iso_absCH=0;
      float iso_absPh=0;
      float iso_absEl=0;
      float iso_absMu=0;
      float iso_absNH=0;
      float iso_absTot=0;
      if( mcColl!= 0 ){
	for(int g=0;g<mcColl->getNumberOfElements();g++){
	  MCParticle *genp = static_cast<MCParticle*>(mcColl->getElementAt(g));
	  //check only stable particles
	  if(genp->getGeneratorStatus()!=1){
	    continue;
	  }
	  //take out neutrinos
	  if(abs(genp->getPDG())==12 || abs(genp->getPDG())==14 || abs(genp->getPDG())==16){
	    continue;
	  }
	  TLorentzVector gen(0,0,0,0);
	  gen.SetPxPyPzE(genp->getMomentum()[0],genp->getMomentum()[1],genp->getMomentum()[2],genp->getEnergy());
	  if(genp==genisop){
	    std::cout<<"seems not the same here"<<std::endl;
	  }
	  if(gen!=geniso){
	    if((TMath::RadToDeg()*gen.Angle(geniso.Vect()))<m_angleIso){
	      iso_absTot+=genp->getEnergy();
	      //std::cout<<"gen iso check for "<<genisop->getEnergy()<< "/"<<genp->getEnergy()<<" type1/type2 "<<genisop->getPDG()<<"/"<<genp->getPDG()<<" status "<<genp->getGeneratorStatus()<<" angle "<<TMath::RadToDeg()*gen.Angle(geniso.Vect())<<" iso angle should be "<<m_angleIso <<" geniso after "<<iso_absTot/genisop->getEnergy()<<std::endl;
	      if(genp->getCharge()!=0){
		if(abs(genp->getPDG())==11){
		  iso_absEl+=genp->getEnergy();
		}else if(abs(genp->getPDG())==13){
		  iso_absMu+=genp->getEnergy();
		}else{
		  iso_absCH+=genp->getEnergy();
		}
	      }else{
		if(genp->getPDG()==22){
		  iso_absPh+=genp->getEnergy();
		}else{
		  iso_absNH+=genp->getEnergy();
		}
	      }
	    }
	  }//else{
	  //std::cout<<"particle iso/gen normal not equal "<<m<<"/"<<g<<std::endl;
	  //}
	}
	
	m_isoPartGenDR10_relIso->push_back(iso_absTot/genisop->getEnergy());
	m_isoPartGenDR10_relIsoCH->push_back(iso_absCH/genisop->getEnergy());
	m_isoPartGenDR10_relIsoPh->push_back(iso_absPh/genisop->getEnergy());
	m_isoPartGenDR10_relIsoNH->push_back(iso_absNH/genisop->getEnergy());
	m_isoPartGenDR10_relIsoEl->push_back(iso_absEl/genisop->getEnergy());
	m_isoPartGenDR10_relIsoMu->push_back(iso_absMu/genisop->getEnergy());
	if((iso_absTot/genisop->getEnergy())>0.10){
	  std::cout<<" HUGE genisop thing, what is wrong here "<<iso_absTot/genisop->getEnergy()<<std::endl;
	}
      }
    }

  }

  LCCollection * genTrueLepPhColl =0;
  getCollection(genTrueLepPhColl,m_outputMCTrueLepPhParticleCollection,evt);

  if(genTrueLepPhColl!=0){
  for(int m=0;m<genTrueLepPhColl->getNumberOfElements();m++){
      MCParticle *gentruelepphp = static_cast<MCParticle*>(genTrueLepPhColl->getElementAt(m));
      if(gentruelepphp->getEnergy()<m_genTrueLepPhEMin){
	continue;
      }
      m_genTrueLepPh_E->push_back(gentruelepphp->getEnergy());
      m_genTrueLepPh_Px->push_back(gentruelepphp->getMomentum()[0]);
      m_genTrueLepPh_Py->push_back(gentruelepphp->getMomentum()[1]);
      m_genTrueLepPh_Pz->push_back(gentruelepphp->getMomentum()[2]);
      m_genTrueLepPh_PDGID->push_back(gentruelepphp->getPDG());
      TLorentzVector gentruelepph(0,0,0,0);
      gentruelepph.SetPxPyPzE(gentruelepphp->getMomentum()[0],gentruelepphp->getMomentum()[1],gentruelepphp->getMomentum()[2],gentruelepphp->getEnergy());
      float iso_absCH=0;
      float iso_absPh=0;
      float iso_absEl=0;
      float iso_absMu=0;
      float iso_absNH=0;
      float iso_absTot=0;
      if( mcColl!= 0 ){
	for(int g=0;g<mcColl->getNumberOfElements();g++){
	  MCParticle *genp = static_cast<MCParticle*>(mcColl->getElementAt(g));
	  //check only stable particles
	  if(genp->getGeneratorStatus()!=1){
	    continue;
	  }
	  //take out neutrinos
	  if(abs(genp->getPDG())==12 || abs(genp->getPDG())==14 || abs(genp->getPDG())==16){
	    continue;
	  }
	  TLorentzVector gen(0,0,0,0);
	  gen.SetPxPyPzE(genp->getMomentum()[0],genp->getMomentum()[1],genp->getMomentum()[2],genp->getEnergy());
	  if(genp==gentruelepphp){
	    std::cout<<"seems not the same here"<<std::endl;
	  }
	  if(gen!=gentruelepph){
	    if((TMath::RadToDeg()*gen.Angle(gentruelepph.Vect()))<m_angleIso){
	      iso_absTot+=genp->getEnergy();
	      if(genp->getCharge()!=0){
		if(abs(genp->getPDG())==11){
		  iso_absEl+=genp->getEnergy();
		}else if(abs(genp->getPDG())==13){
		  iso_absMu+=genp->getEnergy();
		}else{
		  iso_absCH+=genp->getEnergy();
		}
	      }else{
		if(genp->getPDG()==22){
		  iso_absPh+=genp->getEnergy();
		}else{
		  iso_absNH+=genp->getEnergy();
		}
	      }
	    }
	  }else{ 
	    if(genp->getPDG()!=22){
	      if(genp->getParents().size()>0 && abs(genp->getParents()[0]->getPDG())>22 && abs(genp->getParents()[0]->getPDG())<26){
		m_genTrueLepPh_ANC_PDGID->push_back(genp->getParents()[0]->getPDG());
	      }else if(genp->getParents()[0]->getParents().size()>0 &&  abs(genp->getParents()[0]->getParents()[0]->getPDG())>22 && abs(genp->getParents()[0]->getParents()[0]->getPDG())<26){
		m_genTrueLepPh_ANC_PDGID->push_back(genp->getParents()[0]->getParents()[0]->getPDG());
	      }else if(genp->getParents()[0]->getParents()[0]->getParents().size()>0 && abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())>22 && abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())<26){
		m_genTrueLepPh_ANC_PDGID->push_back(genp->getParents()[0]->getParents()[0]->getParents()[0]->getPDG());
	      }else if(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents().size()>0 && abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())>22 && abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())<26){
		m_genTrueLepPh_ANC_PDGID->push_back(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG());
	      }else if(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents().size()>0 &&  abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())>22 && abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())<26){
		m_genTrueLepPh_ANC_PDGID->push_back(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG());
	      }else{
		m_genTrueLepPh_ANC_PDGID->push_back(-30);
	      }
	    }else{
	      //photon irradidated of H (or W, Z should not work, but well use the same code) --> or photon originated from the beam -->then parents of e+ or e- don't exis
	      if(genp->getParents().size()>0 && ( (abs(genp->getParents()[0]->getPDG())>22 && abs(genp->getParents()[0]->getPDG())<26) || (abs(genp->getParents()[0]->getPDG())==11 || genp->getParents()[0]->getParents().size()==0)) ){
		//for ISR photons both incoming beams appear as parents in the event history -->assign to 0 then
		if(abs(genp->getParents()[0]->getPDG())==11){
		  m_genTrueLepPh_ANC_PDGID->push_back(0);
		}else{
		  m_genTrueLepPh_ANC_PDGID->push_back(genp->getParents()[0]->getPDG());
		}
	      }else if (genp->getParents()[0]->getParents().size()>0 && ( (abs(genp->getParents()[0]->getParents()[0]->getPDG())>22 && abs(genp->getParents()[0]->getParents()[0]->getPDG())<26) || (abs(genp->getParents()[0]->getParents()[0]->getPDG())==11 || genp->getParents()[0]->getParents()[0]->getParents().size()==0)) ){

		if(abs(genp->getParents()[0]->getParents()[0]->getPDG())==11){
		  m_genTrueLepPh_ANC_PDGID->push_back(0);
		}else{
		m_genTrueLepPh_ANC_PDGID->push_back(genp->getParents()[0]->getParents()[0]->getPDG());
		}
	      }else if (genp->getParents()[0]->getParents()[0]->getParents().size()>0 && ( (abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())>22 && abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())<26) || (abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())==11 || genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents().size()==0)) ){

		if(abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())==11){
		  m_genTrueLepPh_ANC_PDGID->push_back(0);
		}else{
		  m_genTrueLepPh_ANC_PDGID->push_back(genp->getParents()[0]->getParents()[0]->getParents()[0]->getPDG());
		}
	      }else if (genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents().size()>0 && ( (abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())>22 && abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())<26) || (abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())==11 || genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents().size()==0)) ){

		if(abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())==11){
		  m_genTrueLepPh_ANC_PDGID->push_back(0);
		}else{
		m_genTrueLepPh_ANC_PDGID->push_back(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG());
		}
	      }else if (genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents().size()>0 && ( (abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())>22 && abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())<26) || (abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())==11 || genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents().size()==0)) ){

		if(abs(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG())==11){
		  m_genTrueLepPh_ANC_PDGID->push_back(0);
		}else{
		  m_genTrueLepPh_ANC_PDGID->push_back(genp->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getParents()[0]->getPDG());
		}
	      }else{
		m_genTrueLepPh_ANC_PDGID->push_back(-30);
	      }
	    }
	  }
	}	
	m_genTrueLepPh_relIso->push_back(iso_absTot/gentruelepphp->getEnergy());
	m_genTrueLepPh_relIsoCH->push_back(iso_absCH/gentruelepphp->getEnergy());
	m_genTrueLepPh_relIsoPh->push_back(iso_absPh/gentruelepphp->getEnergy());
	m_genTrueLepPh_relIsoNH->push_back(iso_absNH/gentruelepphp->getEnergy());
	m_genTrueLepPh_relIsoEl->push_back(iso_absEl/gentruelepphp->getEnergy());
	m_genTrueLepPh_relIsoMu->push_back(iso_absMu/gentruelepphp->getEnergy());
	//if((iso_absTot/gentruelepphp->getEnergy())>0.10){
	//std::cout<<m_eventCount<<" HUGE gentruelepphp thing, what is wrong here PDG/E/reliso "<<gentruelepphp->getPDG()<<"/"<<gentruelepphp->getEnergy()<<"/"<<iso_absTot/gentruelepphp->getEnergy()<<std::endl;
	//}
      }
    }
  }
  TLorentzVector temp_sqrtS(0,0,0,0);
  if( mcColl != 0 ){
    int nMCP = mcColl->getNumberOfElements();



    for(int m=0;m<nMCP;m++){
      MCParticle *mcp = static_cast<MCParticle*>(mcColl->getElementAt(m));
      //structure of events --> indices 0-1 incoming 4-5 outgoing electrons after ISR (i.e. the "real collision"
      //index 6 and 7 ISR photons of incoming electrons
      //index 8 is Z-->qq
      //index 9 is Higgs, but take daughter of that Higgs (gets boost to real system),
      //index 10 and 11 are qq from Z,
      //index 12 and 13 are decay products from H,
      //of that daughter choose the decay product (typically bbar)
      if(m>3 && m<14){
	m_trueME_E->push_back(mcp->getEnergy());
	m_trueME_Px->push_back(mcp->getMomentum()[0]);
	m_trueME_Py->push_back(mcp->getMomentum()[1]);
	m_trueME_Pz->push_back(mcp->getMomentum()[2]);
	m_trueME_PDGID->push_back(mcp->getPDG());
      }
      if(m==4 || m==5){
	TLorentzVector temp(0,0,0,0);
	temp.SetPxPyPzE(mcp->getMomentum()[0],mcp->getMomentum()[1],mcp->getMomentum()[2],mcp->getEnergy());
	temp_sqrtS+=temp;
      }
      /*
      if(m==9){
	//std::cout<<"should be higgs "<<mcp->getPDG()<<" no daugthers "<<mcp->getDaughters().size()<<"/"<<mcp->getDaughters()[0]->getEnergy()<<std::endl;
	if((mcp->getDaughters().size())>0){
	  m_trueME_E->push_back(mcp->getDaughters()[0]->getEnergy());
	  m_trueME_Px->push_back(mcp->getDaughters()[0]->getMomentum()[0]);
	  m_trueME_Py->push_back(mcp->getDaughters()[0]->getMomentum()[1]);
	  m_trueME_Pz->push_back(mcp->getDaughters()[0]->getMomentum()[2]);
	  m_trueME_PDGID->push_back(mcp->getDaughters()[0]->getPDG());
	  for(unsigned int i=0;i<mcp->getDaughters()[0]->getDaughters().size();i++){
	    if(mcp->getDaughters()[0]->getDaughters().size()!=2){
	      //std::cout<<"expected structure not in line with HZqq structure "<<mcp->getDaughters()[0]->getDaughters().size()<<std::endl;
	    }
	    m_trueME_E->push_back(mcp->getDaughters()[0]->getDaughters()[i]->getEnergy());
	    m_trueME_Px->push_back(mcp->getDaughters()[0]->getDaughters()[i]->getMomentum()[0]);
	    m_trueME_Py->push_back(mcp->getDaughters()[0]->getDaughters()[i]->getMomentum()[1]);
	    m_trueME_Pz->push_back(mcp->getDaughters()[0]->getDaughters()[i]->getMomentum()[2]);
	    m_trueME_PDGID->push_back(mcp->getDaughters()[0]->getDaughters()[i]->getPDG());
	  }
	}
      }
      */
      if((abs(mcp->getPDG())==24 || mcp->getPDG()==23  || mcp->getPDG()==25) && mcp->getDaughters().size()>1){
	for(unsigned int i=0;i<mcp->getDaughters().size();i++){
	  if(abs(mcp->getDaughters()[i]->getPDG())==15){
	    std::set<MCParticle*> tau_daughtersFunc;
	    fillStableDaughterSet(mcp->getDaughters()[i], tau_daughtersFunc);
	    float tau_genjet_E=0;
	    float tau_genjet_px=0;
	    float tau_genjet_py=0;
	    float tau_genjet_pz=0;
	    int tau_genjet_Charge=0;
	    int tau_genjet_Mult=0;
	    int tau_genjet_NCH=0;
	    int tau_genjet_NPh=0;
	    int tau_genjet_NNH=0;
	    float tau_genjet_CHF=0;
	    float tau_genjet_PhF=0;
	    float tau_genjet_ElF=0;
	    float tau_genjet_MuF=0;
	    float tau_genjet_NHF=0;
	    std::set<MCParticle*>::iterator tauDaughtIt;
	    for(tauDaughtIt=tau_daughtersFunc.begin();tauDaughtIt!=tau_daughtersFunc.end();tauDaughtIt++){
	      if(abs((*tauDaughtIt)->getPDG())==16 ||abs((*tauDaughtIt)->getPDG())==14 ||abs((*tauDaughtIt)->getPDG())==12 ){
		continue;
	      }
	      tau_genjet_E+=(*tauDaughtIt)->getEnergy();
	      tau_genjet_px+=(*tauDaughtIt)->getMomentum()[0];
	      tau_genjet_py+=(*tauDaughtIt)->getMomentum()[1];
	      tau_genjet_pz+=(*tauDaughtIt)->getMomentum()[2];
	      tau_genjet_Mult+=1;
	      if((*tauDaughtIt)->getCharge()==0){
		if((*tauDaughtIt)->getPDG()==22){
		  tau_genjet_NPh+=1;
		  tau_genjet_PhF+=(*tauDaughtIt)->getEnergy();
		}else{
		  tau_genjet_NNH+=1;
		  tau_genjet_NHF+=(*tauDaughtIt)->getEnergy();
		}
	      }else{
		tau_genjet_Charge+=(*tauDaughtIt)->getCharge();
		if(abs((*tauDaughtIt)->getPDG())==11){
		  tau_genjet_ElF+=(*tauDaughtIt)->getEnergy();
		}else if(abs((*tauDaughtIt)->getPDG())==13){
		  tau_genjet_MuF+=(*tauDaughtIt)->getEnergy();
		}else{
		  tau_genjet_NCH+=1;
		  tau_genjet_CHF+=(*tauDaughtIt)->getEnergy();
		}
	      }
	    }
	    m_gen_tau_E->push_back(tau_genjet_E);
	    m_gen_tau_Px->push_back(tau_genjet_px);
	    m_gen_tau_Py->push_back(tau_genjet_py);
	    m_gen_tau_Pz->push_back(tau_genjet_pz);
	    m_gen_tau_Mult->push_back(tau_genjet_Mult);
	    m_gen_tau_Charge->push_back(tau_genjet_Charge);
	    m_gen_tau_MCCharge->push_back(mcp->getDaughters()[i]->getCharge());
	    m_gen_tau_NCH->push_back(tau_genjet_NCH);
	    m_gen_tau_NPh->push_back(tau_genjet_NPh);
	    m_gen_tau_NNH->push_back(tau_genjet_NNH);
	    m_gen_tau_MotherPDGID->push_back(mcp->getPDG());
	    if(tau_genjet_E!=0){
	      m_gen_tau_CHFraction->push_back(tau_genjet_CHF/tau_genjet_E);
	      m_gen_tau_PhFraction->push_back(tau_genjet_PhF/tau_genjet_E);
	      m_gen_tau_ElFraction->push_back(tau_genjet_ElF/tau_genjet_E);
	      m_gen_tau_MuFraction->push_back(tau_genjet_MuF/tau_genjet_E);
	      m_gen_tau_NHFraction->push_back(tau_genjet_NHF/tau_genjet_E);
	    }
	  }
	}
      }


      if(mcp->getGeneratorStatus()==1){
	if(abs(mcp->getPDG())==12 || abs(mcp->getPDG())==14 || abs(mcp->getPDG())==16){
	  m_true_inv_E+=mcp->getEnergy();
	  m_true_inv_Px+=mcp->getMomentum()[0];
	  m_true_inv_Py+=mcp->getMomentum()[1];
	  m_true_inv_Pz+=mcp->getMomentum()[2];
	  m_true_inv_Mult+=1;
	}else{
	  m_true_E+=mcp->getEnergy();
	  m_true_Px+=mcp->getMomentum()[0];
	  m_true_Py+=mcp->getMomentum()[1];
	  m_true_Pz+=mcp->getMomentum()[2];
	  m_true_Mult+=1;
	  if(mcp->getCharge()==0){
	    if(mcp->getPDG()==22){
	      m_true_PhFraction+=mcp->getEnergy();
	      m_true_NPh+=1;
	    }else{
	      m_true_NHFraction+=mcp->getEnergy();
	      m_true_NNH+=1;
	    }
	  }else{
	    if(abs(mcp->getPDG())==11){
	      m_true_ElFraction+=mcp->getEnergy();
	    }else if(abs(mcp->getPDG())==13){
	      m_true_MuFraction+=mcp->getEnergy();
	    }else{
	      m_true_CHFraction+=mcp->getEnergy();
	      m_true_NCH+=1;
	    }
	  }
 	}
      }
    }
    if(m_true_E!=0){
      m_true_CHFraction=m_true_CHFraction/m_true_E;
      m_true_PhFraction=m_true_PhFraction/m_true_E;
      m_true_ElFraction=m_true_ElFraction/m_true_E;
      m_true_MuFraction=m_true_MuFraction/m_true_E;
      m_true_NHFraction=m_true_NHFraction/m_true_E;
    }
    //std::cout<<"sqrtS is "<<temp_sqrtS.M()<<std::endl;
  }
 
  
  LCCollection* particleMCJetIn(NULL);
  getCollection(particleMCJetIn,m_inputMCJetParticleCollection,evt);
  if(particleMCJetIn->getNumberOfElements()>1){
    PseudoJetList genpjList;
    for(int i = 0; i < particleMCJetIn->getNumberOfElements(); ++i){
      ReconstructedParticle* par = static_cast<ReconstructedParticle*> (particleMCJetIn->getElementAt(i));
      if(par->getEnergy()==0){
	std::cout<<"input particle "<<i<<" E/px/py/pz/type "<<par->getEnergy()<<"/"<<par->getMomentum()[0]<<"/"<<par->getMomentum()[1]<<"/"<<par->getMomentum()[2]<<"/"<<par->getType()<<std::endl;
      }
      genpjList.push_back( fastjet::PseudoJet( par->getMomentum()[0],
					       par->getMomentum()[1],
					       par->getMomentum()[2],
					       par->getEnergy() ) );
      genpjList.back().set_user_index(i);	// save the id of this recParticle
    }


    fastjet::ClusterSequence*_csgen = new fastjet::ClusterSequence(genpjList, *_jetAlgoType);


    PseudoJetList genjets = _csgen->exclusive_jets((int)(2));

    m_gen_y21_max=_csgen->exclusive_ymerge_max((int)(1));
    m_gen_y32_max=_csgen->exclusive_ymerge_max((int)(2));
    m_gen_y43_max=_csgen->exclusive_ymerge_max((int)(3));

    m_gen_y21=_csgen->exclusive_ymerge((int)(1));
    m_gen_y32=_csgen->exclusive_ymerge((int)(2));
    m_gen_y43=_csgen->exclusive_ymerge((int)(3));
    
    unsigned int genjet_index_=0;
    for (std::vector<fastjet::PseudoJet>::iterator genjetIt = genjets.begin(); genjetIt != genjets.end(); ++genjetIt ) {
      m_genjet_E->push_back(genjetIt->E());
      m_genjet_Px->push_back(genjetIt->px());
      m_genjet_Py->push_back(genjetIt->py());
      m_genjet_Pz->push_back(genjetIt->pz());
      m_genjet_Mult->push_back(_csgen->constituents((*genjetIt)).size());
      
      int num_pi=0;
      int num_pi_trackPtMin=0;
      int num_ph=0;
      int num_neutHad=0;
      
      float E_pi=0;
      float E_pi_trackPtMin=0;
      float E_ph=0;
      float E_el=0;
      float E_mu=0;
      float E_neutHad=0;

      //float genjet_jetChargeE_kappa_0_25=0;
      //float genjet_jetChargeE_kappa_0_50=0;
      //float genjet_jetChargeE_kappa_0_75=0;
      //float genjet_jetChargeE_kappa_1_00=0;
      //float genjet_jetChargeE_kappa_0_10=0;
      //float genjet_jetChargeE_kappa_0_15=0;
      //float genjet_jetChargeE_kappa_0_20=0;
      //float genjet_jetChargeE_kappa_0_30=0;

      //float genjet_jetChargePProj_kappa_0_25=0;
      //float genjet_jetChargePProj_kappa_0_50=0;
      //float genjet_jetChargePProj_kappa_0_75=0;
      //float genjet_jetChargePProj_kappa_1_00=0;
      //float genjet_jetChargePProj_kappa_0_10=0;
      //float genjet_jetChargePProj_kappa_0_15=0;
      //float genjet_jetChargePProj_kappa_0_20=0;
      //float genjet_jetChargePProj_kappa_0_30=0;

      //float genjet_jetPProj_kappa_0_25=0;
      //float genjet_jetPProj_kappa_0_50=0;
      //float genjet_jetPProj_kappa_0_75=0;
      //float genjet_jetPProj_kappa_1_00=0;
      //float genjet_jetPProj_kappa_0_10=0;
      //float genjet_jetPProj_kappa_0_15=0;
      //float genjet_jetPProj_kappa_0_20=0;
      //float genjet_jetPProj_kappa_0_30=0;
      
      //float genjet_jetChargePt_kappa_0_25=0;
      //float genjet_jetChargePt_kappa_0_50=0;
      //float genjet_jetChargePt_kappa_0_75=0;
      //float genjet_jetChargePt_kappa_1_00=0;
      //float genjet_jetChargePt_kappa_0_10=0;
      //float genjet_jetChargePt_kappa_0_15=0;
      //float genjet_jetChargePt_kappa_0_20=0;
      //float genjet_jetChargePt_kappa_0_30=0;
      
      for(unsigned int i=0;i<_csgen->constituents((*genjetIt)).size();i++){
	//genjet-particles are prepared as reconstructed particle
	//reco Type is the original PDGID, charge is filled as well
	ReconstructedParticle* mcgenreco = static_cast<ReconstructedParticle*> (particleMCJetIn->getElementAt(_csgen->constituents((*genjetIt))[i].user_index()));
	//if(sqrt(pow(mcgenreco->getMomentum()[0],2)+pow(mcgenreco->getMomentum()[1],2))>m_cutPtTrackMin){
	  //genjet_jetChargeE_kappa_0_25+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/genjetIt->E(),0.25);
	  //genjet_jetChargeE_kappa_0_50+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/genjetIt->E(),0.50);
	  //genjet_jetChargeE_kappa_0_75+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/genjetIt->E(),0.75);
	  //genjet_jetChargeE_kappa_1_00+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/genjetIt->E(),1.00);
	  //genjet_jetChargeE_kappa_0_10+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/genjetIt->E(),0.10);
	  //genjet_jetChargeE_kappa_0_15+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/genjetIt->E(),0.15);
	  //genjet_jetChargeE_kappa_0_20+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/genjetIt->E(),0.20);
	  //genjet_jetChargeE_kappa_0_30+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/genjetIt->E(),0.30);
	  //float genpartPt=sqrt(pow(mcgenreco->getMomentum()[0],2)+pow(mcgenreco->getMomentum()[1],2));
	  //genjet_jetChargePt_kappa_0_25+=mcgenreco->getCharge()*pow(genpartPt,0.25);
	  //genjet_jetChargePt_kappa_0_50+=mcgenreco->getCharge()*pow(genpartPt,0.50);
	  //genjet_jetChargePt_kappa_0_75+=mcgenreco->getCharge()*pow(genpartPt,0.75);
	  //genjet_jetChargePt_kappa_1_00+=mcgenreco->getCharge()*pow(genpartPt,1.00);
	  //genjet_jetChargePt_kappa_0_10+=mcgenreco->getCharge()*pow(genpartPt,0.10);
	  //genjet_jetChargePt_kappa_0_15+=mcgenreco->getCharge()*pow(genpartPt,0.15);
	  //genjet_jetChargePt_kappa_0_20+=mcgenreco->getCharge()*pow(genpartPt,0.20);
	  //genjet_jetChargePt_kappa_0_30+=mcgenreco->getCharge()*pow(genpartPt,0.30);
	  //genjet_jetChargePProj_kappa_0_25+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.25);
	  //genjet_jetChargePProj_kappa_0_50+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.50);
	  //genjet_jetChargePProj_kappa_0_75+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.75);
	  //genjet_jetChargePProj_kappa_1_00+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),1.00);
	  //genjet_jetChargePProj_kappa_0_10+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.10);
	  //genjet_jetChargePProj_kappa_0_15+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.15);
	  //genjet_jetChargePProj_kappa_0_20+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.20);
	  //genjet_jetChargePProj_kappa_0_30+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.30);
	//}
	//genjet_jetPProj_kappa_0_25+=pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.25);
	//genjet_jetPProj_kappa_0_50+=pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.50);
	//genjet_jetPProj_kappa_0_75+=pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.75);
	//genjet_jetPProj_kappa_1_00+=pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),1.00);
	//genjet_jetPProj_kappa_0_10+=pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.10);
	//genjet_jetPProj_kappa_0_15+=pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.15);
	//genjet_jetPProj_kappa_0_20+=pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.20);
	//genjet_jetPProj_kappa_0_30+=pow(fabs(mcgenreco->getMomentum()[0]*genjetIt->px()+mcgenreco->getMomentum()[1]*genjetIt->py()+mcgenreco->getMomentum()[2]*genjetIt->pz()),0.30);
	if(mcgenreco->getCharge()==0){
	  if(mcgenreco->getType()==22){
	    num_ph+=1;
	    E_ph+=mcgenreco->getEnergy();
	  }else{
	    num_neutHad+=1;
	    E_neutHad+=mcgenreco->getEnergy();
	  }
	}else{
	  if(abs(mcgenreco->getType())==11){
	    E_el+=mcgenreco->getEnergy();
	  }else if(abs(mcgenreco->getType())==13){
	    E_mu+=mcgenreco->getEnergy();
	  }else{//should be identified as pion on RECO level
	    num_pi+=1;
	    E_pi+=mcgenreco->getEnergy();
	    if(sqrt(pow(mcgenreco->getMomentum()[0],2)+pow(mcgenreco->getMomentum()[1],2))>m_cutPtTrackMin){
	      num_pi_trackPtMin+=1;
	      E_pi_trackPtMin+=mcgenreco->getEnergy();
	    }
	  }
	}
      }
      m_genjet_NCH->push_back(num_pi);
      m_genjet_NCH_trackPtMin->push_back(num_pi_trackPtMin);
      m_genjet_NPh->push_back(num_ph);
      m_genjet_NNH->push_back(num_neutHad);
      
      m_genjet_CHFraction->push_back(E_pi/genjetIt->E());
      m_genjet_CHFraction_trackPtMin->push_back(E_pi_trackPtMin/genjetIt->E());
      m_genjet_PhFraction->push_back(E_ph/genjetIt->E());
      m_genjet_ElFraction->push_back(E_el/genjetIt->E());
      m_genjet_MuFraction->push_back(E_mu/genjetIt->E());
      m_genjet_NHFraction->push_back(E_neutHad/genjetIt->E());

      //m_genjet_jetChargeE_kappa_0_25->push_back(genjet_jetChargeE_kappa_0_25);
      //m_genjet_jetChargeE_kappa_0_50->push_back(genjet_jetChargeE_kappa_0_50);
      //m_genjet_jetChargeE_kappa_0_75->push_back(genjet_jetChargeE_kappa_0_75);
      //m_genjet_jetChargeE_kappa_1_00->push_back(genjet_jetChargeE_kappa_1_00);
      //m_genjet_jetChargeE_kappa_0_10->push_back(genjet_jetChargeE_kappa_0_10);
      //m_genjet_jetChargeE_kappa_0_15->push_back(genjet_jetChargeE_kappa_0_15);
      //m_genjet_jetChargeE_kappa_0_20->push_back(genjet_jetChargeE_kappa_0_20);
      //m_genjet_jetChargeE_kappa_0_30->push_back(genjet_jetChargeE_kappa_0_30);
      
      //m_genjet_jetChargePt_kappa_0_25->push_back(genjet_jetChargePt_kappa_0_25/pow(genjetIt->pt(),0.25));
      //m_genjet_jetChargePt_kappa_0_50->push_back(genjet_jetChargePt_kappa_0_50/pow(genjetIt->pt(),0.50));
      //m_genjet_jetChargePt_kappa_0_75->push_back(genjet_jetChargePt_kappa_0_75/pow(genjetIt->pt(),0.75));
      //m_genjet_jetChargePt_kappa_1_00->push_back(genjet_jetChargePt_kappa_1_00/pow(genjetIt->pt(),1.00));
      //m_genjet_jetChargePt_kappa_0_10->push_back(genjet_jetChargePt_kappa_0_10/pow(genjetIt->pt(),0.10));
      //m_genjet_jetChargePt_kappa_0_15->push_back(genjet_jetChargePt_kappa_0_15/pow(genjetIt->pt(),0.15));
      //m_genjet_jetChargePt_kappa_0_20->push_back(genjet_jetChargePt_kappa_0_20/pow(genjetIt->pt(),0.20));
      //m_genjet_jetChargePt_kappa_0_30->push_back(genjet_jetChargePt_kappa_0_30/pow(genjetIt->pt(),0.30));

      //m_genjet_jetChargePProj_kappa_0_25->push_back(genjet_jetChargePProj_kappa_0_25/genjet_jetPProj_kappa_0_25);
      //m_genjet_jetChargePProj_kappa_0_50->push_back(genjet_jetChargePProj_kappa_0_50/genjet_jetPProj_kappa_0_50);
      //m_genjet_jetChargePProj_kappa_0_75->push_back(genjet_jetChargePProj_kappa_0_75/genjet_jetPProj_kappa_0_75);
      //m_genjet_jetChargePProj_kappa_1_00->push_back(genjet_jetChargePProj_kappa_1_00/genjet_jetPProj_kappa_1_00);
      //m_genjet_jetChargePProj_kappa_0_10->push_back(genjet_jetChargePProj_kappa_0_10/genjet_jetPProj_kappa_0_10);
      //m_genjet_jetChargePProj_kappa_0_15->push_back(genjet_jetChargePProj_kappa_0_15/genjet_jetPProj_kappa_0_15);
      //m_genjet_jetChargePProj_kappa_0_20->push_back(genjet_jetChargePProj_kappa_0_20/genjet_jetPProj_kappa_0_20);
      //m_genjet_jetChargePProj_kappa_0_30->push_back(genjet_jetChargePProj_kappa_0_30/genjet_jetPProj_kappa_0_30);

      m_genjet_nsubjettiness1->push_back((_nSubJettiness1_ptR->result(*genjetIt)));
      m_genjet_nsubjettiness1_lrz->push_back((_nSubJettiness1_lorentz->result(*genjetIt)));

      m_genjet_beta1_ECorr2->push_back(_energycorr2_beta1->result(*genjetIt));
      m_genjet_beta1_ECorr3->push_back(_energycorr3_beta1->result(*genjetIt));
      m_genjet_beta1_N2->push_back(_energycorrN2_beta1->result(*genjetIt));
      m_genjet_beta1_N3->push_back(_energycorrN3_beta1->result(*genjetIt));
      m_genjet_beta1_C2->push_back(_energycorrC2_beta1->result(*genjetIt));
      if(_energycorr3_beta1->result(*genjetIt)!=0){
	m_genjet_beta1_C3->push_back(_energycorr4_beta1->result(*genjetIt)*_energycorr2_beta1->result(*genjetIt)/pow(_energycorr3_beta1->result(*genjetIt),2));
      }else{
	m_genjet_beta1_C3->push_back(0.);
      }
      m_genjet_beta1_D2->push_back(_energycorrD2_beta1->result(*genjetIt));
      m_genjet_beta1_ECorr2_E_theta->push_back(_energycorr2_Etheta_beta1->result(*genjetIt));
      m_genjet_beta1_ECorr3_E_theta->push_back(_energycorr3_Etheta_beta1->result(*genjetIt));
      m_genjet_beta1_N2_E_theta->push_back(_energycorrN2_Etheta_beta1->result(*genjetIt));
      m_genjet_beta1_N3_E_theta->push_back(_energycorrN3_Etheta_beta1->result(*genjetIt));
      m_genjet_beta1_C2_E_theta->push_back(_energycorrC2_Etheta_beta1->result(*genjetIt));
      if(_energycorr3_Etheta_beta1->result(*genjetIt)!=0){
	m_genjet_beta1_C3_E_theta->push_back(_energycorr4_Etheta_beta1->result(*genjetIt)*_energycorr2_Etheta_beta1->result(*genjetIt)/pow(_energycorr3_Etheta_beta1->result(*genjetIt),2));
      }else{
	m_genjet_beta1_C3_E_theta->push_back(0.);
      }
      m_genjet_beta1_D2_E_theta->push_back(_energycorrD2_Etheta_beta1->result(*genjetIt));

      m_genjet_beta2_ECorr2->push_back(_energycorr2_beta2->result(*genjetIt));
      m_genjet_beta2_ECorr3->push_back(_energycorr3_beta2->result(*genjetIt));
      m_genjet_beta2_N2->push_back(_energycorrN2_beta2->result(*genjetIt));
      m_genjet_beta2_N3->push_back(_energycorrN3_beta2->result(*genjetIt));
      m_genjet_beta2_C2->push_back(_energycorrC2_beta2->result(*genjetIt));
      if(_energycorr3_beta2->result(*genjetIt)!=0){
	m_genjet_beta2_C3->push_back(_energycorr4_beta2->result(*genjetIt)*_energycorr2_beta2->result(*genjetIt)/pow(_energycorr3_beta2->result(*genjetIt),2));
      }else{
	m_genjet_beta2_C3->push_back(0.);
      }
      m_genjet_beta2_D2->push_back(_energycorrD2_beta2->result(*genjetIt));
      m_genjet_beta2_ECorr2_E_theta->push_back(_energycorr2_Etheta_beta2->result(*genjetIt));
      m_genjet_beta2_ECorr3_E_theta->push_back(_energycorr3_Etheta_beta2->result(*genjetIt));
      m_genjet_beta2_N2_E_theta->push_back(_energycorrN2_Etheta_beta2->result(*genjetIt));
      m_genjet_beta2_N3_E_theta->push_back(_energycorrN3_Etheta_beta2->result(*genjetIt));
      m_genjet_beta2_C2_E_theta->push_back(_energycorrC2_Etheta_beta2->result(*genjetIt));
      if(_energycorr3_Etheta_beta2->result(*genjetIt)!=0){
	m_genjet_beta2_C3_E_theta->push_back(_energycorr4_Etheta_beta2->result(*genjetIt)*_energycorr2_Etheta_beta2->result(*genjetIt)/pow(_energycorr3_Etheta_beta2->result(*genjetIt),2));
      }else{
	m_genjet_beta2_C3_E_theta->push_back(0.);
      }
      m_genjet_beta2_D2_E_theta->push_back(_energycorrD2_Etheta_beta2->result(*genjetIt));

      m_genjet_beta0_5_ECorr2->push_back(_energycorr2_beta0_5->result(*genjetIt));
      m_genjet_beta0_5_ECorr3->push_back(_energycorr3_beta0_5->result(*genjetIt));
      m_genjet_beta0_5_N2->push_back(_energycorrN2_beta0_5->result(*genjetIt));
      m_genjet_beta0_5_N3->push_back(_energycorrN3_beta0_5->result(*genjetIt));
      m_genjet_beta0_5_C2->push_back(_energycorrC2_beta0_5->result(*genjetIt));
      if(_energycorr3_beta0_5->result(*genjetIt)!=0){
	m_genjet_beta0_5_C3->push_back(_energycorr4_beta0_5->result(*genjetIt)*_energycorr2_beta0_5->result(*genjetIt)/pow(_energycorr3_beta0_5->result(*genjetIt),2));
      }else{
	m_genjet_beta0_5_C3->push_back(0.);
      }
      m_genjet_beta0_5_D2->push_back(_energycorrD2_beta0_5->result(*genjetIt));
      m_genjet_beta0_5_ECorr2_E_theta->push_back(_energycorr2_Etheta_beta0_5->result(*genjetIt));
      m_genjet_beta0_5_ECorr3_E_theta->push_back(_energycorr3_Etheta_beta0_5->result(*genjetIt));
      m_genjet_beta0_5_N2_E_theta->push_back(_energycorrN2_Etheta_beta0_5->result(*genjetIt));
      m_genjet_beta0_5_N3_E_theta->push_back(_energycorrN3_Etheta_beta0_5->result(*genjetIt));
      m_genjet_beta0_5_C2_E_theta->push_back(_energycorrC2_Etheta_beta0_5->result(*genjetIt));
      if(_energycorr3_Etheta_beta0_5->result(*genjetIt)!=0){
	m_genjet_beta0_5_C3_E_theta->push_back(_energycorr4_Etheta_beta0_5->result(*genjetIt)*_energycorr2_Etheta_beta0_5->result(*genjetIt)/pow(_energycorr3_Etheta_beta0_5->result(*genjetIt),2));
      }else{
	m_genjet_beta0_5_C3_E_theta->push_back(0.);
      }
      m_genjet_beta0_5_D2_E_theta->push_back(_energycorrD2_Etheta_beta0_5->result(*genjetIt));
      
      //decompose into requested number of subjets: 2 jets here, thus check if this can be done --> not possible for one 
      if(_csgen->constituents((*genjetIt)).size()>1){
	std::vector<fastjet::PseudoJet> subjets = _csgen->exclusive_subjets(*genjetIt, 2);
	//std::cout<<"has jet subjets "<<genjetIt->has_exclusive_subjets()<<std::endl;
	m_genjet_dij_21->push_back(genjetIt->exclusive_subdmerge(1));
	m_genjet_dij_32->push_back(genjetIt->exclusive_subdmerge(2));
	m_genjet_dij_43->push_back(genjetIt->exclusive_subdmerge(3));
 	m_genjet_dij_21_max->push_back(genjetIt->exclusive_subdmerge_max(1));
	m_genjet_dij_32_max->push_back(genjetIt->exclusive_subdmerge_max(2));
	m_genjet_dij_43_max->push_back(genjetIt->exclusive_subdmerge_max(3));

	//std::cout<<"input const "<<_csgen->constituents((*genjetIt)).size()<<" subt output list "<<_csgen->constituents(subjets[0]).size()<<"/"<<_csgen->constituents(subjets[1]).size()<<std::endl;
	int num_pi_sj=0;
	int num_pi_sj_trackPtMin=0;
	int num_ph_sj=0;
	int num_neutHad_sj=0;
	
	float E_pi_sj=0;
	float E_pi_sj_trackPtMin=0;
	float E_ph_sj=0;
	float E_el_sj=0;
	float E_mu_sj=0;
	float E_neutHad_sj=0;
	float genjet_subjet_jetChargeE_kappa_0_25=0;
	float genjet_subjet_jetChargeE_kappa_0_50=0;
	//float genjet_subjet_jetChargeE_kappa_0_75=0;
	//float genjet_subjet_jetChargeE_kappa_1_00=0;
	//float genjet_subjet_jetChargeE_kappa_0_10=0;
	//float genjet_subjet_jetChargeE_kappa_0_15=0;
	float genjet_subjet_jetChargeE_kappa_0_20=0;
	float genjet_subjet_jetChargeE_kappa_0_30=0;
	
	float genjet_subjet_jetChargePt_kappa_0_25=0;
	float genjet_subjet_jetChargePt_kappa_0_50=0;
	//float genjet_subjet_jetChargePt_kappa_0_75=0;
	//float genjet_subjet_jetChargePt_kappa_1_00=0;
	//float genjet_subjet_jetChargePt_kappa_0_10=0;
	//float genjet_subjet_jetChargePt_kappa_0_15=0;
	float genjet_subjet_jetChargePt_kappa_0_20=0;
	float genjet_subjet_jetChargePt_kappa_0_30=0;

	//float genjet_subjet_jetPProj_kappa_0_25=0;
	//float genjet_subjet_jetPProj_kappa_0_50=0;
	//float genjet_subjet_jetPProj_kappa_0_75=0;
	//float genjet_subjet_jetPProj_kappa_1_00=0;
	//float genjet_subjet_jetPProj_kappa_0_10=0;
	//float genjet_subjet_jetPProj_kappa_0_15=0;
	//float genjet_subjet_jetPProj_kappa_0_20=0;
	//float genjet_subjet_jetPProj_kappa_0_30=0;
	
	//float genjet_subjet_jetChargePProj_kappa_0_25=0;
	//float genjet_subjet_jetChargePProj_kappa_0_50=0;
	//float genjet_subjet_jetChargePProj_kappa_0_75=0;
	//float genjet_subjet_jetChargePProj_kappa_1_00=0;
	//float genjet_subjet_jetChargePProj_kappa_0_10=0;
	//loat genjet_subjet_jetChargePProj_kappa_0_15=0;
	//float genjet_subjet_jetChargePProj_kappa_0_20=0;
	//float genjet_subjet_jetChargePProj_kappa_0_30=0;
	
	m_genjet_subjet_E->push_back(subjets[0].E());
	m_genjet_subjet_Px->push_back(subjets[0].px());
	m_genjet_subjet_Py->push_back(subjets[0].py());
	m_genjet_subjet_Pz->push_back(subjets[0].pz());
	m_genjet_subjet_jetindex->push_back(genjet_index_);
	m_genjet_nsubjettiness2->push_back((_nSubJettiness2_ptR->result(*genjetIt)));
	//m_genjet_nsubjettiness2_lrz->push_back((_nSubJettiness2_lorentz->result(*genjetIt)));
	//double tau2_ptR = nSubJettiness2_ptR(*genjetIt);
	//double tau2_lorentz = nSubJettiness2_lorentz(*genjetIt);


	for(unsigned int i=0;i<_csgen->constituents(subjets[0]).size();i++){
	  ReconstructedParticle* mcgenreco = static_cast<ReconstructedParticle*> (particleMCJetIn->getElementAt(_csgen->constituents(subjets[0])[i].user_index()));
	  if(sqrt(pow(mcgenreco->getMomentum()[0],2)+pow(mcgenreco->getMomentum()[1],2))>m_cutPtTrackMin){
	    genjet_subjet_jetChargeE_kappa_0_25+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[0].E(),0.25);
	    genjet_subjet_jetChargeE_kappa_0_50+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[0].E(),0.50);
	    //genjet_subjet_jetChargeE_kappa_0_75+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[0].E(),0.75);
	    //genjet_subjet_jetChargeE_kappa_1_00+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[0].E(),1.00);
	    //genjet_subjet_jetChargeE_kappa_0_10+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[0].E(),0.10);
	    //genjet_subjet_jetChargeE_kappa_0_15+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[0].E(),0.15);
	    genjet_subjet_jetChargeE_kappa_0_20+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[0].E(),0.20);
	    genjet_subjet_jetChargeE_kappa_0_30+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[0].E(),0.30);
	    float genpartPt=sqrt(pow(mcgenreco->getMomentum()[0],2)+pow(mcgenreco->getMomentum()[1],2));
	    genjet_subjet_jetChargePt_kappa_0_25+=mcgenreco->getCharge()*pow(genpartPt,0.25);
	    genjet_subjet_jetChargePt_kappa_0_50+=mcgenreco->getCharge()*pow(genpartPt,0.50);
	    //genjet_subjet_jetChargePt_kappa_0_75+=mcgenreco->getCharge()*pow(genpartPt,0.75);
	    //genjet_subjet_jetChargePt_kappa_1_00+=mcgenreco->getCharge()*pow(genpartPt,1.00);
	    //genjet_subjet_jetChargePt_kappa_0_10+=mcgenreco->getCharge()*pow(genpartPt,0.10);
	    //genjet_subjet_jetChargePt_kappa_0_15+=mcgenreco->getCharge()*pow(genpartPt,0.15);
	    genjet_subjet_jetChargePt_kappa_0_20+=mcgenreco->getCharge()*pow(genpartPt,0.20);
	    genjet_subjet_jetChargePt_kappa_0_30+=mcgenreco->getCharge()*pow(genpartPt,0.30);
	    //genjet_subjet_jetChargePProj_kappa_0_25+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.25);
	    //genjet_subjet_jetChargePProj_kappa_0_50+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.50);
	    //genjet_subjet_jetChargePProj_kappa_0_75+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.75);
	    //genjet_subjet_jetChargePProj_kappa_1_00+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),1.00);
	    //genjet_subjet_jetChargePProj_kappa_0_10+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.10);
	    //genjet_subjet_jetChargePProj_kappa_0_15+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.15);
	    //genjet_subjet_jetChargePProj_kappa_0_20+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.20);
	    //genjet_subjet_jetChargePProj_kappa_0_30+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.30);
	  }
	  //genjet_subjet_jetPProj_kappa_0_25+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.25);
	  //genjet_subjet_jetPProj_kappa_0_50+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.50);
	  //genjet_subjet_jetPProj_kappa_0_75+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.75);
	  //genjet_subjet_jetPProj_kappa_1_00+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),1.00);
	  //genjet_subjet_jetPProj_kappa_0_10+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.10);
	  //genjet_subjet_jetPProj_kappa_0_15+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.15);
	  //genjet_subjet_jetPProj_kappa_0_20+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.20);
	  //genjet_subjet_jetPProj_kappa_0_30+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[0].px()+mcgenreco->getMomentum()[1]*subjets[0].py()+mcgenreco->getMomentum()[2]*subjets[0].pz()),0.30);
	  if(mcgenreco->getCharge()==0){
	    if(mcgenreco->getType()==22){
	      num_ph_sj+=1;
	      E_ph_sj+=mcgenreco->getEnergy();
	    }else{
	      num_neutHad_sj+=1;
	      E_neutHad_sj+=mcgenreco->getEnergy();
	    }
	  }else{
	    if(abs(mcgenreco->getType())==11){
	      E_el_sj+=mcgenreco->getEnergy();
	    }else if(abs(mcgenreco->getType())==13){
	      E_mu_sj+=mcgenreco->getEnergy();
	    }else{//should be identified as pion on RECO level
	      num_pi_sj+=1;
	      E_pi_sj+=mcgenreco->getEnergy();
	      if(sqrt(pow(mcgenreco->getMomentum()[0],2)+pow(mcgenreco->getMomentum()[1],2))>m_cutPtTrackMin){
		num_pi_sj_trackPtMin+=1;
		E_pi_sj_trackPtMin+=mcgenreco->getEnergy();
	      }
	    }
	  }
	}
	m_genjet_subjet_NCH->push_back(num_pi_sj);
	m_genjet_subjet_NCH_trackPtMin->push_back(num_pi_sj_trackPtMin);
	m_genjet_subjet_NPh->push_back(num_ph_sj);
	m_genjet_subjet_NNH->push_back(num_neutHad_sj);
	
	m_genjet_subjet_CHFraction->push_back(E_pi_sj/subjets[0].E());
	m_genjet_subjet_CHFraction_trackPtMin->push_back(E_pi_sj_trackPtMin/subjets[0].E());
	m_genjet_subjet_PhFraction->push_back(E_ph_sj/subjets[0].E());
	m_genjet_subjet_ElFraction->push_back(E_el_sj/subjets[0].E());
	m_genjet_subjet_MuFraction->push_back(E_mu_sj/subjets[0].E());
	m_genjet_subjet_NHFraction->push_back(E_neutHad_sj/subjets[0].E());

	m_genjet_subjet_jetChargeE_kappa_0_25->push_back(genjet_subjet_jetChargeE_kappa_0_25);
	m_genjet_subjet_jetChargeE_kappa_0_50->push_back(genjet_subjet_jetChargeE_kappa_0_50);
	//m_genjet_subjet_jetChargeE_kappa_0_75->push_back(genjet_subjet_jetChargeE_kappa_0_75);
	//m_genjet_subjet_jetChargeE_kappa_1_00->push_back(genjet_subjet_jetChargeE_kappa_1_00);
	//m_genjet_subjet_jetChargeE_kappa_0_10->push_back(genjet_subjet_jetChargeE_kappa_0_10);
	//m_genjet_subjet_jetChargeE_kappa_0_15->push_back(genjet_subjet_jetChargeE_kappa_0_15);
	m_genjet_subjet_jetChargeE_kappa_0_20->push_back(genjet_subjet_jetChargeE_kappa_0_20);
	m_genjet_subjet_jetChargeE_kappa_0_30->push_back(genjet_subjet_jetChargeE_kappa_0_30);
	m_genjet_subjet_jetChargePt_kappa_0_25->push_back(genjet_subjet_jetChargePt_kappa_0_25/pow(subjets[0].pt(),0.25));
	m_genjet_subjet_jetChargePt_kappa_0_50->push_back(genjet_subjet_jetChargePt_kappa_0_50/pow(subjets[0].pt(),0.50));
	//m_genjet_subjet_jetChargePt_kappa_0_75->push_back(genjet_subjet_jetChargePt_kappa_0_75/pow(subjets[0].pt(),0.75));
	//m_genjet_subjet_jetChargePt_kappa_1_00->push_back(genjet_subjet_jetChargePt_kappa_1_00/pow(subjets[0].pt(),1.00));
	//m_genjet_subjet_jetChargePt_kappa_0_10->push_back(genjet_subjet_jetChargePt_kappa_0_10/pow(subjets[0].pt(),0.10));
	//m_genjet_subjet_jetChargePt_kappa_0_15->push_back(genjet_subjet_jetChargePt_kappa_0_15/pow(subjets[0].pt(),0.15));
	m_genjet_subjet_jetChargePt_kappa_0_20->push_back(genjet_subjet_jetChargePt_kappa_0_20/pow(subjets[0].pt(),0.20));
	m_genjet_subjet_jetChargePt_kappa_0_30->push_back(genjet_subjet_jetChargePt_kappa_0_30/pow(subjets[0].pt(),0.30));

	//m_genjet_subjet_jetChargePProj_kappa_0_25->push_back(genjet_subjet_jetChargePProj_kappa_0_25/genjet_subjet_jetPProj_kappa_0_25);
	//m_genjet_subjet_jetChargePProj_kappa_0_50->push_back(genjet_subjet_jetChargePProj_kappa_0_50/genjet_subjet_jetPProj_kappa_0_50);
	//m_genjet_subjet_jetChargePProj_kappa_0_75->push_back(genjet_subjet_jetChargePProj_kappa_0_75/genjet_subjet_jetPProj_kappa_0_75);
	//m_genjet_subjet_jetChargePProj_kappa_1_00->push_back(genjet_subjet_jetChargePProj_kappa_1_00/genjet_subjet_jetPProj_kappa_1_00);
	//m_genjet_subjet_jetChargePProj_kappa_0_10->push_back(genjet_subjet_jetChargePProj_kappa_0_10/genjet_subjet_jetPProj_kappa_0_10);
	//m_genjet_subjet_jetChargePProj_kappa_0_15->push_back(genjet_subjet_jetChargePProj_kappa_0_15/genjet_subjet_jetPProj_kappa_0_15);
	//m_genjet_subjet_jetChargePProj_kappa_0_20->push_back(genjet_subjet_jetChargePProj_kappa_0_20/genjet_subjet_jetPProj_kappa_0_20);
	//m_genjet_subjet_jetChargePProj_kappa_0_30->push_back(genjet_subjet_jetChargePProj_kappa_0_30/genjet_subjet_jetPProj_kappa_0_30);

	
	num_pi_sj=0;
	num_pi_sj_trackPtMin=0;
	num_ph_sj=0;
	num_neutHad_sj=0;
	
	E_pi_sj=0;
	E_pi_sj_trackPtMin=0;
	E_ph_sj=0;
	E_el_sj=0;
	E_mu_sj=0;
	E_neutHad_sj=0;
	m_genjet_subjet_E->push_back(subjets[1].E());
	m_genjet_subjet_Px->push_back(subjets[1].px());
	m_genjet_subjet_Py->push_back(subjets[1].py());
	m_genjet_subjet_Pz->push_back(subjets[1].pz());
	m_genjet_subjet_jetindex->push_back(genjet_index_);


	genjet_subjet_jetChargeE_kappa_0_25=0;
	genjet_subjet_jetChargeE_kappa_0_50=0;
	//genjet_subjet_jetChargeE_kappa_0_75=0;
	//genjet_subjet_jetChargeE_kappa_1_00=0;
	//genjet_subjet_jetChargeE_kappa_0_10=0;
	//genjet_subjet_jetChargeE_kappa_0_15=0;
	genjet_subjet_jetChargeE_kappa_0_20=0;
	genjet_subjet_jetChargeE_kappa_0_30=0;
	
	genjet_subjet_jetChargePt_kappa_0_25=0;
	genjet_subjet_jetChargePt_kappa_0_50=0;
	//genjet_subjet_jetChargePt_kappa_0_75=0;
	//genjet_subjet_jetChargePt_kappa_1_00=0;
	//genjet_subjet_jetChargePt_kappa_0_10=0;
	//genjet_subjet_jetChargePt_kappa_0_15=0;
	genjet_subjet_jetChargePt_kappa_0_20=0;
	genjet_subjet_jetChargePt_kappa_0_30=0;

	//genjet_subjet_jetPProj_kappa_0_25=0;
	//genjet_subjet_jetPProj_kappa_0_50=0;
	//genjet_subjet_jetPProj_kappa_0_75=0;
	//genjet_subjet_jetPProj_kappa_1_00=0;
	//genjet_subjet_jetPProj_kappa_0_10=0;
	//genjet_subjet_jetPProj_kappa_0_15=0;
	//genjet_subjet_jetPProj_kappa_0_20=0;
	//genjet_subjet_jetPProj_kappa_0_30=0;
	
	//genjet_subjet_jetChargePProj_kappa_0_25=0;
	//genjet_subjet_jetChargePProj_kappa_0_50=0;
	//genjet_subjet_jetChargePProj_kappa_0_75=0;
	//genjet_subjet_jetChargePProj_kappa_1_00=0;
	//genjet_subjet_jetChargePProj_kappa_0_10=0;
	//genjet_subjet_jetChargePProj_kappa_0_15=0;
	//genjet_subjet_jetChargePProj_kappa_0_20=0;
	//genjet_subjet_jetChargePProj_kappa_0_30=0;

	for(unsigned int i=0;i<_csgen->constituents(subjets[1]).size();i++){
	  ReconstructedParticle* mcgenreco = static_cast<ReconstructedParticle*> (particleMCJetIn->getElementAt(_csgen->constituents(subjets[1])[i].user_index()));
	  if(sqrt(pow(mcgenreco->getMomentum()[0],2)+pow(mcgenreco->getMomentum()[1],2))>m_cutPtTrackMin){
	    genjet_subjet_jetChargeE_kappa_0_25+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[1].E(),0.25);
	    genjet_subjet_jetChargeE_kappa_0_50+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[1].E(),0.50);
	    //genjet_subjet_jetChargeE_kappa_0_75+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[1].E(),0.75);
	    //genjet_subjet_jetChargeE_kappa_1_00+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[1].E(),1.00);
	    //genjet_subjet_jetChargeE_kappa_0_10+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[1].E(),0.10);
	    //genjet_subjet_jetChargeE_kappa_0_15+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[1].E(),0.15);
	    genjet_subjet_jetChargeE_kappa_0_20+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[1].E(),0.20);
	    genjet_subjet_jetChargeE_kappa_0_30+=mcgenreco->getCharge()*pow(mcgenreco->getEnergy()/subjets[1].E(),0.30);
	    float genpartPt=sqrt(pow(mcgenreco->getMomentum()[0],2)+pow(mcgenreco->getMomentum()[1],2));
	    genjet_subjet_jetChargePt_kappa_0_25+=mcgenreco->getCharge()*pow(genpartPt,0.25);
	    genjet_subjet_jetChargePt_kappa_0_50+=mcgenreco->getCharge()*pow(genpartPt,0.50);
	    //genjet_subjet_jetChargePt_kappa_0_75+=mcgenreco->getCharge()*pow(genpartPt,0.75);
	    //genjet_subjet_jetChargePt_kappa_1_00+=mcgenreco->getCharge()*pow(genpartPt,1.00);
	    //genjet_subjet_jetChargePt_kappa_0_10+=mcgenreco->getCharge()*pow(genpartPt,0.10);
	    //genjet_subjet_jetChargePt_kappa_0_15+=mcgenreco->getCharge()*pow(genpartPt,0.15);
	    genjet_subjet_jetChargePt_kappa_0_20+=mcgenreco->getCharge()*pow(genpartPt,0.20);
	    genjet_subjet_jetChargePt_kappa_0_30+=mcgenreco->getCharge()*pow(genpartPt,0.30);
	    //genjet_subjet_jetChargePProj_kappa_0_25+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.25);
	    //genjet_subjet_jetChargePProj_kappa_0_50+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.50);
	    //genjet_subjet_jetChargePProj_kappa_0_75+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.75);
	    //genjet_subjet_jetChargePProj_kappa_1_00+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),1.00);
	    //genjet_subjet_jetChargePProj_kappa_0_10+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.10);
	    //genjet_subjet_jetChargePProj_kappa_0_15+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.15);
	    //genjet_subjet_jetChargePProj_kappa_0_20+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.20);
	    //genjet_subjet_jetChargePProj_kappa_0_30+=mcgenreco->getCharge()*pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.30);
	  }
	  //genjet_subjet_jetPProj_kappa_0_25+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.25);
	  //genjet_subjet_jetPProj_kappa_0_50+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.50);
	  //genjet_subjet_jetPProj_kappa_0_75+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.75);
	  //genjet_subjet_jetPProj_kappa_1_00+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),1.00);
	  //genjet_subjet_jetPProj_kappa_0_10+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.10);
	  //genjet_subjet_jetPProj_kappa_0_15+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.15);
	  //genjet_subjet_jetPProj_kappa_0_20+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.20);
	  //genjet_subjet_jetPProj_kappa_0_30+=pow(fabs(mcgenreco->getMomentum()[0]*subjets[1].px()+mcgenreco->getMomentum()[1]*subjets[1].py()+mcgenreco->getMomentum()[2]*subjets[1].pz()),0.30);
	  if(mcgenreco->getCharge()==0){
	    if(mcgenreco->getType()==22){
	      num_ph_sj+=1;
	      E_ph_sj+=mcgenreco->getEnergy();
	    }else{
	      num_neutHad_sj+=1;
	      E_neutHad_sj+=mcgenreco->getEnergy();
	    }
	  }else{
	    if(abs(mcgenreco->getType())==11){
	      E_el_sj+=mcgenreco->getEnergy();
	    }else if(abs(mcgenreco->getType())==13){
	      E_mu_sj+=mcgenreco->getEnergy();
	    }else{//should be identified as pion on RECO level
	      num_pi_sj+=1;
	      E_pi_sj+=mcgenreco->getEnergy();
	      if(sqrt(pow(mcgenreco->getMomentum()[0],2)+pow(mcgenreco->getMomentum()[1],2))>m_cutPtTrackMin){
		num_pi_sj_trackPtMin+=1;
		E_pi_sj_trackPtMin+=mcgenreco->getEnergy();
	      }
	    }
	  }
	}
	m_genjet_subjet_NCH_trackPtMin->push_back(num_pi_sj_trackPtMin);
	m_genjet_subjet_NCH->push_back(num_pi_sj);
	m_genjet_subjet_NPh->push_back(num_ph_sj);
	m_genjet_subjet_NNH->push_back(num_neutHad_sj);
	
	m_genjet_subjet_CHFraction_trackPtMin->push_back(E_pi_sj_trackPtMin/subjets[1].E());
	m_genjet_subjet_CHFraction->push_back(E_pi_sj/subjets[1].E());
	m_genjet_subjet_PhFraction->push_back(E_ph_sj/subjets[1].E());
	m_genjet_subjet_ElFraction->push_back(E_el_sj/subjets[1].E());
	m_genjet_subjet_MuFraction->push_back(E_mu_sj/subjets[1].E());
	m_genjet_subjet_NHFraction->push_back(E_neutHad_sj/subjets[1].E());

	m_genjet_subjet_jetChargeE_kappa_0_25->push_back(genjet_subjet_jetChargeE_kappa_0_25);
	m_genjet_subjet_jetChargeE_kappa_0_50->push_back(genjet_subjet_jetChargeE_kappa_0_50);
	//m_genjet_subjet_jetChargeE_kappa_0_75->push_back(genjet_subjet_jetChargeE_kappa_0_75);
	//m_genjet_subjet_jetChargeE_kappa_1_00->push_back(genjet_subjet_jetChargeE_kappa_1_00);
	//m_genjet_subjet_jetChargeE_kappa_0_10->push_back(genjet_subjet_jetChargeE_kappa_0_10);
	//m_genjet_subjet_jetChargeE_kappa_0_15->push_back(genjet_subjet_jetChargeE_kappa_0_15);
	m_genjet_subjet_jetChargeE_kappa_0_20->push_back(genjet_subjet_jetChargeE_kappa_0_20);
	m_genjet_subjet_jetChargeE_kappa_0_30->push_back(genjet_subjet_jetChargeE_kappa_0_30);
	m_genjet_subjet_jetChargePt_kappa_0_25->push_back(genjet_subjet_jetChargePt_kappa_0_25/pow(subjets[1].pt(),0.25));
	m_genjet_subjet_jetChargePt_kappa_0_50->push_back(genjet_subjet_jetChargePt_kappa_0_50/pow(subjets[1].pt(),0.50));
	//m_genjet_subjet_jetChargePt_kappa_0_75->push_back(genjet_subjet_jetChargePt_kappa_0_75/pow(subjets[1].pt(),0.75));
	//m_genjet_subjet_jetChargePt_kappa_1_00->push_back(genjet_subjet_jetChargePt_kappa_1_00/pow(subjets[1].pt(),1.00));
	//m_genjet_subjet_jetChargePt_kappa_0_10->push_back(genjet_subjet_jetChargePt_kappa_0_10/pow(subjets[1].pt(),0.10));
	//m_genjet_subjet_jetChargePt_kappa_0_15->push_back(genjet_subjet_jetChargePt_kappa_0_15/pow(subjets[1].pt(),0.15));
	m_genjet_subjet_jetChargePt_kappa_0_20->push_back(genjet_subjet_jetChargePt_kappa_0_20/pow(subjets[1].pt(),0.20));
	m_genjet_subjet_jetChargePt_kappa_0_30->push_back(genjet_subjet_jetChargePt_kappa_0_30/pow(subjets[1].pt(),0.30));

	//m_genjet_subjet_jetChargePProj_kappa_0_25->push_back(genjet_subjet_jetChargePProj_kappa_0_25/genjet_subjet_jetPProj_kappa_0_25);
	//m_genjet_subjet_jetChargePProj_kappa_0_50->push_back(genjet_subjet_jetChargePProj_kappa_0_50/genjet_subjet_jetPProj_kappa_0_50);
	//m_genjet_subjet_jetChargePProj_kappa_0_75->push_back(genjet_subjet_jetChargePProj_kappa_0_75/genjet_subjet_jetPProj_kappa_0_75);
	//m_genjet_subjet_jetChargePProj_kappa_1_00->push_back(genjet_subjet_jetChargePProj_kappa_1_00/genjet_subjet_jetPProj_kappa_1_00);
	//m_genjet_subjet_jetChargePProj_kappa_0_10->push_back(genjet_subjet_jetChargePProj_kappa_0_10/genjet_subjet_jetPProj_kappa_0_10);
	//m_genjet_subjet_jetChargePProj_kappa_0_15->push_back(genjet_subjet_jetChargePProj_kappa_0_15/genjet_subjet_jetPProj_kappa_0_15);
	//m_genjet_subjet_jetChargePProj_kappa_0_20->push_back(genjet_subjet_jetChargePProj_kappa_0_20/genjet_subjet_jetPProj_kappa_0_20);
	//m_genjet_subjet_jetChargePProj_kappa_0_30->push_back(genjet_subjet_jetChargePProj_kappa_0_30/genjet_subjet_jetPProj_kappa_0_30);

	//decompose into requested number of subjets: 2 jets here, thus check if this can be done --> not possible for one 
	if(_csgen->constituents((*genjetIt)).size()>2){
	  m_genjet_nsubjettiness3->push_back((_nSubJettiness3_ptR->result(*genjetIt)));
	  //m_genjet_nsubjettiness3_lrz->push_back((_nSubJettiness3_lorentz->result(*genjetIt)));
	}else{
	  m_genjet_nsubjettiness3->push_back(-1);
	  //m_genjet_nsubjettiness3_lrz->push_back(-1);
	}
      }else{//less than 2 particles per jet
	m_genjet_dij_21->push_back(0);
	m_genjet_dij_32->push_back(0);
	m_genjet_dij_43->push_back(0);
	m_genjet_dij_21_max->push_back(0);
	m_genjet_dij_32_max->push_back(0);
	m_genjet_dij_43_max->push_back(0);
	m_genjet_nsubjettiness2->push_back(-1);
	//m_genjet_nsubjettiness2_lrz->push_back(-1);
	m_genjet_nsubjettiness3->push_back(-1);
	//m_genjet_nsubjettiness3_lrz->push_back(-1);
      }
      genjet_index_+=1;
    }
    delete _csgen;
  }
  
  LCCollection * tauJetColl =0;
  getCollection(tauJetColl,m_inputTauCollection,evt);
  
  if(tauJetColl!=NULL){
    for(int i=0;i<tauJetColl->getNumberOfElements();i++){
      ReconstructedParticle* tauJet = dynamic_cast<ReconstructedParticle*>(tauJetColl->getElementAt(i));
      m_reco_tau_Charge->push_back(tauJet->getCharge()); 
      m_reco_tau_Px->push_back(tauJet->getMomentum()[0]); 
      m_reco_tau_Py->push_back(tauJet->getMomentum()[1]); 
      m_reco_tau_Pz->push_back(tauJet->getMomentum()[2]); 
      m_reco_tau_E->push_back(tauJet->getEnergy()); 
      m_reco_tau_Mult->push_back(tauJet->getParticles().size());
      int tau_recojet_NCH=0;
      int tau_recojet_NPh=0;
      int tau_recojet_NNH=0;
      float tau_recojet_CHF=0;
      float tau_recojet_PhF=0;
      float tau_recojet_ElF=0;
      float tau_recojet_MuF=0;
      float tau_recojet_NHF=0;
      for(unsigned int td=0;td<tauJet->getParticles().size();td++){
	if(tauJet->getParticles()[td]->getCharge()==0){
	  if(tauJet->getParticles()[td]->getType()==22){
	    tau_recojet_NPh+=1;
	    tau_recojet_PhF+=tauJet->getParticles()[td]->getEnergy();
	  }else{
	    tau_recojet_NNH+=1;
	    tau_recojet_NHF+=tauJet->getParticles()[td]->getEnergy();
	  }
	}else{
	  if(abs(tauJet->getParticles()[td]->getType())==11){
	    tau_recojet_ElF+=tauJet->getParticles()[td]->getEnergy();
	  }else if(abs(tauJet->getParticles()[td]->getType())==13){
	    tau_recojet_MuF+=tauJet->getParticles()[td]->getEnergy();
	  }else{
	    tau_recojet_NCH+=1;
	    tau_recojet_CHF+=tauJet->getParticles()[td]->getEnergy();
	  }
	}
      }
      m_reco_tau_NCH->push_back(tau_recojet_NCH);
      m_reco_tau_NPh->push_back(tau_recojet_NPh);
      m_reco_tau_NNH->push_back(tau_recojet_NNH);
      m_reco_tau_CHFraction->push_back(tau_recojet_CHF/tauJet->getEnergy());
      m_reco_tau_PhFraction->push_back(tau_recojet_PhF/tauJet->getEnergy());
      m_reco_tau_ElFraction->push_back(tau_recojet_ElF/tauJet->getEnergy());
      m_reco_tau_MuFraction->push_back(tau_recojet_MuF/tauJet->getEnergy());
      m_reco_tau_NHFraction->push_back(tau_recojet_NHF/tauJet->getEnergy());
    }
  }
  LCCollection * recoIsoColl =0;
  getCollection(recoIsoColl,m_recoIsoPartColName,evt);
  LCCollection * recoColl =0;
  getCollection(recoColl,m_inputRECOParticleCollection,evt);
  
  if( recoIsoColl != 0 ){
    for(int m=0;m<recoIsoColl->getNumberOfElements();m++){
      ReconstructedParticle *recoisop = static_cast<ReconstructedParticle*>(recoIsoColl->getElementAt(m));
      m_isoPartRecoDR10_E->push_back(recoisop->getEnergy());
      m_isoPartRecoDR10_Px->push_back(recoisop->getMomentum()[0]);
      m_isoPartRecoDR10_Py->push_back(recoisop->getMomentum()[1]);
      m_isoPartRecoDR10_Pz->push_back(recoisop->getMomentum()[2]);
      m_isoPartRecoDR10_PDGID->push_back(recoisop->getType());
      TLorentzVector recoiso(0,0,0,0);
      recoiso.SetPxPyPzE(recoisop->getMomentum()[0],recoisop->getMomentum()[1],recoisop->getMomentum()[2],recoisop->getEnergy());
      float iso_absCH=0;
      float iso_absPh=0;
      float iso_absEl=0;
      float iso_absMu=0;
      float iso_absNH=0;
      float iso_absTot=0;
      if( recoColl!= 0 ){
	for(int r=0;r<recoColl->getNumberOfElements();r++){
	  ReconstructedParticle *recop = static_cast<ReconstructedParticle*>(recoColl->getElementAt(r));
	  TLorentzVector reco(0,0,0,0);
	  reco.SetPxPyPzE(recop->getMomentum()[0],recop->getMomentum()[1],recop->getMomentum()[2],recop->getEnergy());
	  if(recop==recoisop){
	    std::cout<<"seems not the same here"<<std::endl;
	  }
	  if(reco!=recoiso){
	    if((TMath::RadToDeg()*reco.Angle(recoiso.Vect()))<m_angleIso){
	      iso_absTot+=recop->getEnergy();
	      if(recop->getCharge()!=0){
		if(abs(recop->getType())==11){
		  iso_absEl+=recop->getEnergy();
		}else if(abs(recop->getType())==13){
		  iso_absMu+=recop->getEnergy();
		}else{
		  iso_absCH+=recop->getEnergy();
		}
	      }else{
		if(recop->getType()==22){
		  iso_absPh+=recop->getEnergy();
		}else{
		  iso_absNH+=recop->getEnergy();
		}
	      }
	    }
	  }//else{
	  //std::cout<<"particle iso/reco normal not equal "<<m<<"/"<<r<<std::endl;
	  //}
	}
	
	m_isoPartRecoDR10_relIso->push_back(iso_absTot/recoisop->getEnergy());
	m_isoPartRecoDR10_relIsoCH->push_back(iso_absCH/recoisop->getEnergy());
	m_isoPartRecoDR10_relIsoPh->push_back(iso_absPh/recoisop->getEnergy());
	m_isoPartRecoDR10_relIsoNH->push_back(iso_absNH/recoisop->getEnergy());
	m_isoPartRecoDR10_relIsoEl->push_back(iso_absEl/recoisop->getEnergy());
	m_isoPartRecoDR10_relIsoMu->push_back(iso_absMu/recoisop->getEnergy());
	if((iso_absTot/recoisop->getEnergy())>0.10){
	  std::cout<<" HUGE recoisop thing, what is wrong here "<<iso_absTot/recoisop->getEnergy()<<std::endl;
	}
      }
    }

  }

  if( recoColl != 0 ){
    int nRECOP = recoColl->getNumberOfElements();
    for(int m=0;m<nRECOP;m++){
      ReconstructedParticle *recop = static_cast<ReconstructedParticle*>(recoColl->getElementAt(m));
      m_totPFO_E+=recop->getEnergy();
      m_totPFO_Px+=recop->getMomentum()[0];
      m_totPFO_Py+=recop->getMomentum()[1];
      m_totPFO_Pz+=recop->getMomentum()[2];
      m_totPFO_Mult+=1;
      if(recop->getCharge()==0){
	if(recop->getType()==22){
	  m_totPFO_NPh+=1;
	  m_totPFO_PhFraction+=recop->getEnergy();
	}else{
	  m_totPFO_NNH+=1;
	  m_totPFO_NHFraction+=recop->getEnergy();
	}
      }else{
	if(abs(recop->getType())==11){
	  m_totPFO_ElFraction+=recop->getEnergy();
	}else if(abs(recop->getType())==13){
	  m_totPFO_MuFraction+=recop->getEnergy();
	}else{
	  m_totPFO_NCH+=1;
	  m_totPFO_CHFraction+=recop->getEnergy();
	}
      }
    }
    if(m_totPFO_E>0){
      m_totPFO_CHFraction=m_totPFO_CHFraction/m_totPFO_E;
      m_totPFO_PhFraction=m_totPFO_PhFraction/m_totPFO_E;
      m_totPFO_ElFraction=m_totPFO_ElFraction/m_totPFO_E;
      m_totPFO_MuFraction=m_totPFO_MuFraction/m_totPFO_E;
      m_totPFO_NHFraction=m_totPFO_NHFraction/m_totPFO_E;
    }
  }
  
  LCCollection* particleRECOJetIn(NULL);
  getCollection(particleRECOJetIn,m_inputRECOJetParticleCollection,evt);
  
  LCCollection * recorefJet0Coll=0;
  getCollection(recorefJet0Coll,m_recorefinedjet0ColName,evt);

 if(recorefJet0Coll->getNumberOfElements()>0){
     PIDHandler pidh( recorefJet0Coll );
     // get algorithm ID associated with LCFIPlus
    int algo = pidh.getAlgorithmID( "lcfiplus" );
  // get index number for flavor tagging
    int ibtag = pidh.getParameterIndex(algo, "BTag");
    int ictag = pidh.getParameterIndex(algo, "CTag");
    int iotag = pidh.getParameterIndex(algo, "OTag");
    //int ibbtag = pidh.getParameterIndex(algo, "BBTag");
    //int icctag = pidh.getParameterIndex(algo, "CCTag");
    int icat = pidh.getParameterIndex(algo, "Category");
    // loop over jets to extract flavor tagging information
    for(int i=0; i < recorefJet0Coll->getNumberOfElements(); i++) {
      ReconstructedParticle *part = dynamic_cast<ReconstructedParticle*>( recorefJet0Coll->getElementAt( i ) );
      m_recojet_subjet_rfj_j_E->push_back(part->getEnergy());
      m_recojet_subjet_rfj_j_Px->push_back(part->getMomentum()[0]);
      m_recojet_subjet_rfj_j_Py->push_back(part->getMomentum()[1]);
      m_recojet_subjet_rfj_j_Pz->push_back(part->getMomentum()[2]);
      m_recojet_subjet_rfj_j_jetindex->push_back(0);
      m_recojet_subjet_rfj_j_subjetindex->push_back(i);

      //float subjet_pt=sqrt(pow(part->getMomentum()[0],2)+pow(part->getMomentum()[1],2));

      int num_pi_sj=0;
      int num_ph_sj=0;
      int num_neutHad_sj=0;
      
      float E_pi_sj=0;
      float E_ph_sj=0;
      float E_el_sj=0;
      float E_mu_sj=0;
      float E_neutHad_sj=0;

      //float recojet_subjet_jetChargeE_kappa_0_25=0;
      //float recojet_subjet_jetChargeE_kappa_0_50=0;
      //float recojet_subjet_jetChargeE_kappa_0_20=0;
      //float recojet_subjet_jetChargeE_kappa_0_30=0;
      
      //float recojet_subjet_jetChargePt_kappa_0_25=0;
      //float recojet_subjet_jetChargePt_kappa_0_50=0;
      //float recojet_subjet_jetChargePt_kappa_0_20=0;
      //float recojet_subjet_jetChargePt_kappa_0_30=0;
      for(unsigned int j=0; j < part->getParticles().size(); j++) {
	ReconstructedParticle* recopart = static_cast<ReconstructedParticle*> (part->getParticles()[j]);
	TLorentzVector temp(0,0,0,0);
	temp.SetPxPyPzE(recopart->getMomentum()[0],recopart->getMomentum()[1],recopart->getMomentum()[2],recopart->getEnergy());
	/*for(unsigned int j1=(j+1); j1 < part->getParticles().size(); j1++) {
	  ReconstructedParticle* recopart2 = static_cast<ReconstructedParticle*> (part->getParticles()[j1]);
	  TLorentzVector temp2(0,0,0,0);
	  temp2.SetPxPyPzE(recopart2->getMomentum()[0],recopart2->getMomentum()[1],recopart2->getMomentum()[2],recopart2->getEnergy());
	  if(temp.Angle(temp2.Vect())<1.e-3 && fabs(temp.E()-temp2.E())<1.e-3){
	    rfjet_has_duplicate=true;
	    std::cout<<"duplicate particles in jet "<<i<<" from rfj0 "<<j<<" E/px/py/pz "<<recopart->getEnergy()<<"/"<<recopart->getMomentum()[0]<<"/"<<recopart->getMomentum()[1]<<"/"<<recopart->getMomentum()[2]<<" PDG "<<recopart->getType()<<" duplicate to "<<j1<<" E/px/py/pz "<<recopart2->getEnergy()<<"/"<<recopart2->getMomentum()[0]<<"/"<<recopart2->getMomentum()[1]<<"/"<<recopart2->getMomentum()[2]<<" PDG "<<recopart2->getType()<<std::endl;
	  }
	  }*/


	//recojet_subjet_jetChargeE_kappa_0_25+=recopart->getCharge()*pow(recopart->getEnergy()/part->getEnergy(),0.25);
	//recojet_subjet_jetChargeE_kappa_0_50+=recopart->getCharge()*pow(recopart->getEnergy()/part->getEnergy(),0.50);
	//recojet_subjet_jetChargeE_kappa_0_20+=recopart->getCharge()*pow(recopart->getEnergy()/part->getEnergy(),0.20);
	//recojet_subjet_jetChargeE_kappa_0_30+=recopart->getCharge()*pow(recopart->getEnergy()/part->getEnergy(),0.30);
	//float recopartPt=sqrt(pow(recopart->getMomentum()[0],2)+pow(recopart->getMomentum()[1],2));
	//recojet_subjet_jetChargePt_kappa_0_25+=recopart->getCharge()*pow(recopartPt,0.25);
	//recojet_subjet_jetChargePt_kappa_0_50+=recopart->getCharge()*pow(recopartPt,0.50);
	//recojet_subjet_jetChargePt_kappa_0_20+=recopart->getCharge()*pow(recopartPt,0.20);
	//recojet_subjet_jetChargePt_kappa_0_30+=recopart->getCharge()*pow(recopartPt,0.30);
	if(recopart->getCharge()==0){
	  if(recopart->getType()==22){
	    num_ph_sj+=1;
	    E_ph_sj+=recopart->getEnergy();
	  }else{
	    num_neutHad_sj+=1;
	    E_neutHad_sj+=recopart->getEnergy();
	  }
	}else{
	  if(abs(recopart->getType())==11){
	    E_el_sj+=recopart->getEnergy();
	  }else if(abs(recopart->getType())==13){
	    E_mu_sj+=recopart->getEnergy();
	  }else{//should be identified as pion on RECO level
	    num_pi_sj+=1;
	    E_pi_sj+=recopart->getEnergy();
	  }
	}
      }
      /*if(rfjet_has_duplicate){
	for(unsigned int j=0; j < part->getParticles().size(); j++) {
	  ReconstructedParticle* recopart = static_cast<ReconstructedParticle*> (part->getParticles()[j]);
	  std::cout<<"duplicate particles in jet "<<i<<" from rfj0 "<<j<<" E/px/py/pz "<<recopart->getEnergy()<<"/"<<recopart->getMomentum()[0]<<"/"<<recopart->getMomentum()[1]<<"/"<<recopart->getMomentum()[2]<<" PDG "<<recopart->getType()<<std::endl;
	}
	}*/
      m_recojet_subjet_rfj_j_NCH->push_back(num_pi_sj);
      m_recojet_subjet_rfj_j_NPh->push_back(num_ph_sj);
      m_recojet_subjet_rfj_j_NNH->push_back(num_neutHad_sj);
      
      m_recojet_subjet_rfj_j_CHFraction->push_back(E_pi_sj/part->getEnergy());
      m_recojet_subjet_rfj_j_PhFraction->push_back(E_ph_sj/part->getEnergy());
      m_recojet_subjet_rfj_j_ElFraction->push_back(E_el_sj/part->getEnergy());
      m_recojet_subjet_rfj_j_MuFraction->push_back(E_mu_sj/part->getEnergy());
      m_recojet_subjet_rfj_j_NHFraction->push_back(E_neutHad_sj/part->getEnergy());
      
      //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_25->push_back(recojet_subjet_jetChargeE_kappa_0_25);
      //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_50->push_back(recojet_subjet_jetChargeE_kappa_0_50);
      //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_20->push_back(recojet_subjet_jetChargeE_kappa_0_20);
      //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_30->push_back(recojet_subjet_jetChargeE_kappa_0_30);
      //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_25->push_back(recojet_subjet_jetChargePt_kappa_0_25/pow(subjet_pt,0.25));
      //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_50->push_back(recojet_subjet_jetChargePt_kappa_0_50/pow(subjet_pt,0.50));
      //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_20->push_back(recojet_subjet_jetChargePt_kappa_0_20/pow(subjet_pt,0.20));
      //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_30->push_back(recojet_subjet_jetChargePt_kappa_0_30/pow(subjet_pt,0.30));
      const ParticleID &pid = pidh.getParticleID(part, algo);
      m_recojet_subjet_rfj_j_BTag->push_back(pid.getParameters()[ibtag]);
      m_recojet_subjet_rfj_j_CTag->push_back(pid.getParameters()[ictag]);
      m_recojet_subjet_rfj_j_OTag->push_back(pid.getParameters()[iotag]);
      m_recojet_subjet_rfj_j_cat->push_back(pid.getParameters()[icat]);
    }
 }
  LCCollection * SVtxrefJet0Coll=0;
  getCollection(SVtxrefJet0Coll,m_vtx_rfj0ColName,evt);

  for(int i=0; i < SVtxrefJet0Coll->getNumberOfElements(); i++) {
    Vertex *vtx = dynamic_cast<Vertex*>( SVtxrefJet0Coll->getElementAt( i ) );
    m_recojet_subjet_rfj_j_svtx_r->push_back(sqrt(pow(vtx->getPosition()[0],2)+pow(vtx->getPosition()[1],2)));
    ReconstructedParticle *part = vtx->getAssociatedParticle();
    m_recojet_subjet_rfj_j_svtx_E->push_back(part->getEnergy());
    m_recojet_subjet_rfj_j_svtx_Mass->push_back(part->getMass());
    m_recojet_subjet_rfj_j_svtx_nTrack->push_back(part->getParticles().size());
    m_recojet_subjet_rfj_j_svtx_Charge->push_back(part->getCharge());
    m_recojet_subjet_rfj_j_svtx_jetindex->push_back(0);
    /*if(rfjet_has_duplicate){
      std::cout<<"vtx of jet0 rfjets "<<i<<" has "<<part->getParticles().size()<<" particles "<<std::endl;
      for(unsigned int p=0;p<part->getParticles().size();p++){
	ReconstructedParticle *track_p = part->getParticles()[p];
	std::cout<<"vtx "<<i<<" of jet0 rfjets, track "<<p<<" E/px/py/pz "<<track_p->getEnergy()<<"/"<<track_p->getMomentum()[0]<<"/"<<track_p->getMomentum()[1]<<"/"<<track_p->getMomentum()[2]<<" PDG "<<track_p->getType()<<std::endl;
      }
      }*/
  }

  LCCollection * recorefJet1Coll=0;
  getCollection(recorefJet1Coll,m_recorefinedjet1ColName,evt);


 if(recorefJet1Coll->getNumberOfElements()>0){
    PIDHandler pidh( recorefJet1Coll );
    // get algorithm ID associated with LCFIPlus
    int algo = pidh.getAlgorithmID( "lcfiplus" );
  // get index number for flavor tagging
    int ibtag = pidh.getParameterIndex(algo, "BTag");
    int ictag = pidh.getParameterIndex(algo, "CTag");
    int iotag = pidh.getParameterIndex(algo, "OTag");
    //int ibbtag = pidh.getParameterIndex(algo, "BBTag");
    //int icctag = pidh.getParameterIndex(algo, "CCTag");
    int icat = pidh.getParameterIndex(algo, "Category");
    // loop over jets to extract flavor tagging information
    for(int i=0; i < recorefJet1Coll->getNumberOfElements(); i++) {
      ReconstructedParticle *part = dynamic_cast<ReconstructedParticle*>( recorefJet1Coll->getElementAt( i ) );
      m_recojet_subjet_rfj_j_E->push_back(part->getEnergy());
      m_recojet_subjet_rfj_j_Px->push_back(part->getMomentum()[0]);
      m_recojet_subjet_rfj_j_Py->push_back(part->getMomentum()[1]);
      m_recojet_subjet_rfj_j_Pz->push_back(part->getMomentum()[2]);
      m_recojet_subjet_rfj_j_jetindex->push_back(1);
      m_recojet_subjet_rfj_j_subjetindex->push_back(i);

      //float subjet_pt=sqrt(pow(part->getMomentum()[0],2)+pow(part->getMomentum()[1],2));

      int num_pi_sj=0;
      int num_ph_sj=0;
      int num_neutHad_sj=0;
      
      float E_pi_sj=0;
      float E_ph_sj=0;
      float E_el_sj=0;
      float E_mu_sj=0;
      float E_neutHad_sj=0;

      //float recojet_subjet_jetChargeE_kappa_0_25=0;
      //float recojet_subjet_jetChargeE_kappa_0_50=0;
      //float recojet_subjet_jetChargeE_kappa_0_20=0;
      //float recojet_subjet_jetChargeE_kappa_0_30=0;
      
      //float recojet_subjet_jetChargePt_kappa_0_25=0;
      //float recojet_subjet_jetChargePt_kappa_0_50=0;
      //float recojet_subjet_jetChargePt_kappa_0_20=0;
      //float recojet_subjet_jetChargePt_kappa_0_30=0;
      for(unsigned int j=0; j < part->getParticles().size(); j++) {
	ReconstructedParticle* recopart = static_cast<ReconstructedParticle*> (part->getParticles()[j]);
	TLorentzVector temp(0,0,0,0);
	temp.SetPxPyPzE(recopart->getMomentum()[0],recopart->getMomentum()[1],recopart->getMomentum()[2],recopart->getEnergy());
	/*for(unsigned int j1=(j+1); j1 < part->getParticles().size(); j1++) {
	  ReconstructedParticle* recopart2 = static_cast<ReconstructedParticle*> (part->getParticles()[j1]);
	  TLorentzVector temp2(0,0,0,0);
	  temp2.SetPxPyPzE(recopart2->getMomentum()[0],recopart2->getMomentum()[1],recopart2->getMomentum()[2],recopart2->getEnergy());
	  if(temp.Angle(temp2.Vect())<1.e-3 && fabs(temp.E()-temp2.E())<1.e-3){
	    rfjet_has_duplicate=true;
	    //std::cout<<"duplicate particles in jet "<<i<<" from rfj1 "<<j<<" E/px/py/pz "<<recopart->getEnergy()<<"/"<<recopart->getMomentum()[0]<<"/"<<recopart->getMomentum()[1]<<"/"<<recopart->getMomentum()[2]<<" PDG "<<recopart->getType()<<" duplicate to "<<j1<<" E/px/py/pz "<<recopart2->getEnergy()<<"/"<<recopart2->getMomentum()[0]<<"/"<<recopart2->getMomentum()[1]<<"/"<<recopart2->getMomentum()[2]<<" PDG "<<recopart2->getType()<<std::endl;
	  }
	  }*/
	//recojet_subjet_jetChargeE_kappa_0_25+=recopart->getCharge()*pow(recopart->getEnergy()/part->getEnergy(),0.25);
	//recojet_subjet_jetChargeE_kappa_0_50+=recopart->getCharge()*pow(recopart->getEnergy()/part->getEnergy(),0.50);
	//recojet_subjet_jetChargeE_kappa_0_20+=recopart->getCharge()*pow(recopart->getEnergy()/part->getEnergy(),0.20);
	//recojet_subjet_jetChargeE_kappa_0_30+=recopart->getCharge()*pow(recopart->getEnergy()/part->getEnergy(),0.30);
	//float recopartPt=sqrt(pow(recopart->getMomentum()[0],2)+pow(recopart->getMomentum()[1],2));
	//recojet_subjet_jetChargePt_kappa_0_25+=recopart->getCharge()*pow(recopartPt,0.25);
	//recojet_subjet_jetChargePt_kappa_0_50+=recopart->getCharge()*pow(recopartPt,0.50);
	//recojet_subjet_jetChargePt_kappa_0_20+=recopart->getCharge()*pow(recopartPt,0.20);
	//recojet_subjet_jetChargePt_kappa_0_30+=recopart->getCharge()*pow(recopartPt,0.30);
	if(recopart->getCharge()==0){
	  if(recopart->getType()==22){
	    num_ph_sj+=1;
	    E_ph_sj+=recopart->getEnergy();
	  }else{
	    num_neutHad_sj+=1;
	    E_neutHad_sj+=recopart->getEnergy();
	  }
	}else{
	  if(abs(recopart->getType())==11){
	    E_el_sj+=recopart->getEnergy();
	  }else if(abs(recopart->getType())==13){
	    E_mu_sj+=recopart->getEnergy();
	  }else{//should be identified as pion on RECO level
	    num_pi_sj+=1;
	    E_pi_sj+=recopart->getEnergy();
	  }
	}
      }
      //if(rfjet_has_duplicate){
	//for(unsigned int j=0; j < part->getParticles().size(); j++) {
	//ReconstructedParticle* recopart = static_cast<ReconstructedParticle*> (part->getParticles()[j]);
	  //std::cout<<"duplicate particles in jet "<<i<<" from rfj1 "<<j<<" E/px/py/pz "<<recopart->getEnergy()<<"/"<<recopart->getMomentum()[0]<<"/"<<recopart->getMomentum()[1]<<"/"<<recopart->getMomentum()[2]<<" PDG "<<recopart->getType()<<std::endl;
	//}
	//}
      m_recojet_subjet_rfj_j_NCH->push_back(num_pi_sj);
      m_recojet_subjet_rfj_j_NPh->push_back(num_ph_sj);
      m_recojet_subjet_rfj_j_NNH->push_back(num_neutHad_sj);
      
      m_recojet_subjet_rfj_j_CHFraction->push_back(E_pi_sj/part->getEnergy());
      m_recojet_subjet_rfj_j_PhFraction->push_back(E_ph_sj/part->getEnergy());
      m_recojet_subjet_rfj_j_ElFraction->push_back(E_el_sj/part->getEnergy());
      m_recojet_subjet_rfj_j_MuFraction->push_back(E_mu_sj/part->getEnergy());
      m_recojet_subjet_rfj_j_NHFraction->push_back(E_neutHad_sj/part->getEnergy());
      
      //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_25->push_back(recojet_subjet_jetChargeE_kappa_0_25);
      //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_50->push_back(recojet_subjet_jetChargeE_kappa_0_50);
      //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_20->push_back(recojet_subjet_jetChargeE_kappa_0_20);
      //m_recojet_subjet_rfj_j_jetChargeE_kappa_0_30->push_back(recojet_subjet_jetChargeE_kappa_0_30);
      //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_25->push_back(recojet_subjet_jetChargePt_kappa_0_25/pow(subjet_pt,0.25));
      //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_50->push_back(recojet_subjet_jetChargePt_kappa_0_50/pow(subjet_pt,0.50));
      //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_20->push_back(recojet_subjet_jetChargePt_kappa_0_20/pow(subjet_pt,0.20));
      //m_recojet_subjet_rfj_j_jetChargePt_kappa_0_30->push_back(recojet_subjet_jetChargePt_kappa_0_30/pow(subjet_pt,0.30));

      const ParticleID &pid = pidh.getParticleID(part, algo);
      m_recojet_subjet_rfj_j_BTag->push_back(pid.getParameters()[ibtag]);
      m_recojet_subjet_rfj_j_CTag->push_back(pid.getParameters()[ictag]);
      m_recojet_subjet_rfj_j_OTag->push_back(pid.getParameters()[iotag]);
      m_recojet_subjet_rfj_j_cat->push_back(pid.getParameters()[icat]);
    }
  }

 /*
  LCCollection * recorefJet_4jets_Coll=0;
  getCollection(recorefJet_4jets_Coll,m_recorefinedjet_4jets_ColName,evt);

 if(recorefJet_4jets_Coll->getNumberOfElements()>0){
    PIDHandler pidh( recorefJet_4jets_Coll );
    // get algorithm ID associated with LCFIPlus
    int algo = pidh.getAlgorithmID( "lcfiplus" );
  // get index number for flavor tagging
    int ibtag = pidh.getParameterIndex(algo, "BTag");
    int ictag = pidh.getParameterIndex(algo, "CTag");
    int iotag = pidh.getParameterIndex(algo, "OTag");
    //int ibbtag = pidh.getParameterIndex(algo, "BBTag");
    //int icctag = pidh.getParameterIndex(algo, "CCTag");
    int icat = pidh.getParameterIndex(algo, "Category");
    // loop over jets to extract flavor tagging information
    for(int i=0; i < recorefJet_4jets_Coll->getNumberOfElements(); i++) {
      ReconstructedParticle *part = dynamic_cast<ReconstructedParticle*>( recorefJet_4jets_Coll->getElementAt( i ) );
      m_recojet_rfj_4jets_E->push_back(part->getEnergy());
      m_recojet_rfj_4jets_Px->push_back(part->getMomentum()[0]);
      m_recojet_rfj_4jets_Py->push_back(part->getMomentum()[1]);
      m_recojet_rfj_4jets_Pz->push_back(part->getMomentum()[2]);

      int num_pi_sj=0;
      int num_ph_sj=0;
      int num_neutHad_sj=0;
      
      float E_pi_sj=0;
      float E_ph_sj=0;
      float E_el_sj=0;
      float E_mu_sj=0;
      float E_neutHad_sj=0;

      for(unsigned int j=0; j < part->getParticles().size(); j++) {
	ReconstructedParticle* recopart = static_cast<ReconstructedParticle*> (part->getParticles()[j]);
	if(recopart->getCharge()==0){
	  if(recopart->getType()==22){
	    num_ph_sj+=1;
	    E_ph_sj+=recopart->getEnergy();
	  }else{
	    num_neutHad_sj+=1;
	    E_neutHad_sj+=recopart->getEnergy();
	  }
	}else{
	  if(abs(recopart->getType())==11){
	    E_el_sj+=recopart->getEnergy();
	  }else if(abs(recopart->getType())==13){
	    E_mu_sj+=recopart->getEnergy();
	  }else{//should be identified as pion on RECO level
	    num_pi_sj+=1;
	    E_pi_sj+=recopart->getEnergy();
	  }
	}
      }
      m_recojet_rfj_4jets_Mult->push_back(part->getParticles().size());
      m_recojet_rfj_4jets_NCH->push_back(num_pi_sj);
      
      m_recojet_rfj_4jets_CHFraction->push_back(E_pi_sj/part->getEnergy());
      m_recojet_rfj_4jets_PhFraction->push_back(E_ph_sj/part->getEnergy());
      m_recojet_rfj_4jets_ElFraction->push_back(E_el_sj/part->getEnergy());
      m_recojet_rfj_4jets_MuFraction->push_back(E_mu_sj/part->getEnergy());
      m_recojet_rfj_4jets_NHFraction->push_back(E_neutHad_sj/part->getEnergy());

      const ParticleID &pid = pidh.getParticleID(part, algo);
      m_recojet_rfj_4jets_BTag->push_back(pid.getParameters()[ibtag]);
      m_recojet_rfj_4jets_CTag->push_back(pid.getParameters()[ictag]);
      m_recojet_rfj_4jets_OTag->push_back(pid.getParameters()[iotag]);
      m_recojet_rfj_4jets_cat->push_back(pid.getParameters()[icat]);
    }

    
    if(recorefJet_4jets_Coll->getNumberOfElements()>1){
      PseudoJetList recopj_rfj_4jetsList;
      for(int i = 0; i < recorefJet_4jets_Coll->getNumberOfElements(); ++i){
	ReconstructedParticle* par = static_cast<ReconstructedParticle*> (recorefJet_4jets_Coll->getElementAt(i));
	if(par->getEnergy()==0){
	  std::cout<<"input particle "<<i<<" E/px/py/pz/type "<<par->getEnergy()<<"/"<<par->getMomentum()[0]<<"/"<<par->getMomentum()[1]<<"/"<<par->getMomentum()[2]<<"/"<<par->getType()<<std::endl;
	}
	recopj_rfj_4jetsList.push_back( fastjet::PseudoJet( par->getMomentum()[0],
							    par->getMomentum()[1],
							    par->getMomentum()[2],
							    par->getEnergy() ) );
	recopj_rfj_4jetsList.back().set_user_index(i);	// save the id of this recParticle
      }

      
      fastjet::ClusterSequence*_csreco_rfj_4jets = new fastjet::ClusterSequence(recopj_rfj_4jetsList, *_jetAlgoType);
      PseudoJetList recojets = _csreco_rfj_4jets->exclusive_jets((int)(2));
      
      unsigned int recojet_index_=0;
      bool jet_has_2subjets=false;
      for (std::vector<fastjet::PseudoJet>::iterator recojetIt = recojets.begin(); recojetIt != recojets.end(); ++recojetIt ) {
	//std::cout<<"HZAnalzyer "<<particleRECOJetIn->getNumberOfElements()<<" size of jet "<<_csreco_rfj_4jets->constituents((*recojetIt)).size()<<std::endl;
	m_recojet_fat2j_rfj4_E->push_back(recojetIt->E());
	m_recojet_fat2j_rfj4_Px->push_back(recojetIt->px());
	m_recojet_fat2j_rfj4_Py->push_back(recojetIt->py());
	m_recojet_fat2j_rfj4_Pz->push_back(recojetIt->pz());
	m_recojet_fat2j_rfj4_Mult->push_back(_csreco_rfj_4jets->constituents((*recojetIt)).size());

	if(!jet_has_2subjets && _csreco_rfj_4jets->constituents((*recojetIt)).size()==2){
	  jet_has_2subjets=true;
	}
	if(jet_has_2subjets && _csreco_rfj_4jets->constituents((*recojetIt)).size()!=2){
	  std::cout<<"what is wrong, jet1 had 2 jets, and now jet 2 does not "<<recorefJet_4jets_Coll->getNumberOfElements()<<" "<<_csreco_rfj_4jets->constituents((*recojetIt)).size()<<std::endl;
	}
	recojet_index_+=1;
      }
    }
  }
 
  LCCollection * SVtxrefJet_4jets_Coll=0;
  getCollection(SVtxrefJet_4jets_Coll,m_vtx_rfj_4jets_ColName,evt);
  
  for(int i=0; i < SVtxrefJet_4jets_Coll->getNumberOfElements(); i++) {
    Vertex *vtx = dynamic_cast<Vertex*>( SVtxrefJet_4jets_Coll->getElementAt( i ) );
    m_recojet_rfj_4jets_svtx_r->push_back(sqrt(pow(vtx->getPosition()[0],2)+pow(vtx->getPosition()[1],2)));
    ReconstructedParticle *part = vtx->getAssociatedParticle();
    m_recojet_rfj_4jets_svtx_E->push_back(part->getEnergy());
    m_recojet_rfj_4jets_svtx_Mass->push_back(part->getMass());
    m_recojet_rfj_4jets_svtx_nTrack->push_back(part->getParticles().size());
    m_recojet_rfj_4jets_svtx_Charge->push_back(part->getCharge());
  }
 */
 
  LCCollection * SVtxrefJet1Coll=0;
  getCollection(SVtxrefJet1Coll,m_vtx_rfj1ColName,evt);


  for(int i=0; i < SVtxrefJet1Coll->getNumberOfElements(); i++) {
    Vertex *vtx = dynamic_cast<Vertex*>( SVtxrefJet1Coll->getElementAt( i ) );
    m_recojet_subjet_rfj_j_svtx_r->push_back(sqrt(pow(vtx->getPosition()[0],2)+pow(vtx->getPosition()[1],2)));
    ReconstructedParticle *part = vtx->getAssociatedParticle();
    m_recojet_subjet_rfj_j_svtx_E->push_back(part->getEnergy());
    m_recojet_subjet_rfj_j_svtx_Mass->push_back(part->getMass());
    m_recojet_subjet_rfj_j_svtx_nTrack->push_back(part->getParticles().size());
    m_recojet_subjet_rfj_j_svtx_Charge->push_back(part->getCharge());
    m_recojet_subjet_rfj_j_svtx_jetindex->push_back(1);
    //if(rfjet_has_duplicate){
    //std::cout<<"vtx of jet1 rfjets "<<i<<" has "<<part->getParticles().size()<<" particles "<<std::endl;
    //for(unsigned int p=0;p<part->getParticles().size();p++){
    //ReconstructedParticle *track_p = part->getParticles()[p];
    //std::cout<<"vtx "<<i<<" of jet1 rfjets, track "<<p<<" E/px/py/pz "<<track_p->getEnergy()<<"/"<<track_p->getMomentum()[0]<<"/"<<track_p->getMomentum()[1]<<"/"<<track_p->getMomentum()[2]<<" PDG "<<track_p->getType()<<std::endl;
    //}
    //}
  }


  if(particleRECOJetIn->getNumberOfElements()>1){
    PseudoJetList recopjList;
    for(int i = 0; i < particleRECOJetIn->getNumberOfElements(); ++i){
      ReconstructedParticle* par = static_cast<ReconstructedParticle*> (particleRECOJetIn->getElementAt(i));
      if(par->getEnergy()==0){
	std::cout<<"input particle "<<i<<" E/px/py/pz/type "<<par->getEnergy()<<"/"<<par->getMomentum()[0]<<"/"<<par->getMomentum()[1]<<"/"<<par->getMomentum()[2]<<"/"<<par->getType()<<std::endl;
      }
      //if(par->getEnergy()>20){
      //std::cout<<"input particle "<<i<<" E/px/py/pz/type "<<par->getEnergy()<<"/"<<par->getMomentum()[0]<<"/"<<par->getMomentum()[1]<<"/"<<par->getMomentum()[2]<<"/"<<par->getType()<<std::endl;
      //}
      recopjList.push_back( fastjet::PseudoJet( par->getMomentum()[0],
						par->getMomentum()[1],
						par->getMomentum()[2],
						par->getEnergy() ) );
      recopjList.back().set_user_index(i);	// save the id of this recParticle
    }


    fastjet::ClusterSequence*_csreco = new fastjet::ClusterSequence(recopjList, *_jetAlgoType);
    PseudoJetList recojets = _csreco->exclusive_jets((int)(2));

    m_reco_y21_max=_csreco->exclusive_ymerge_max((int)(1));
    m_reco_y32_max=_csreco->exclusive_ymerge_max((int)(2));
    m_reco_y43_max=_csreco->exclusive_ymerge_max((int)(3));

    m_reco_y21=_csreco->exclusive_ymerge((int)(1));
    m_reco_y32=_csreco->exclusive_ymerge((int)(2));
    m_reco_y43=_csreco->exclusive_ymerge((int)(3));
       
    unsigned int recojet_index_=0;
    for (std::vector<fastjet::PseudoJet>::iterator recojetIt = recojets.begin(); recojetIt != recojets.end(); ++recojetIt ) {
      //std::cout<<"HZAnalzyer "<<particleRECOJetIn->getNumberOfElements()<<" size of jet "<<_csreco->constituents((*recojetIt)).size()<<std::endl;
      m_recojet_E->push_back(recojetIt->E());
      m_recojet_Px->push_back(recojetIt->px());
      m_recojet_Py->push_back(recojetIt->py());
      m_recojet_Pz->push_back(recojetIt->pz());
      m_recojet_Mult->push_back(_csreco->constituents((*recojetIt)).size());
      //std::cout<<"jet "<<recojet_index_<<" mass "<<recojetIt->m()<<std::endl;
      
      int num_pi=0;
      int num_pi_trackPtMin=0;
      int num_ph=0;
      int num_neutHad=0;
      
      float E_pi=0;
      float E_pi_trackPtMin=0;
      float E_ph=0;
      float E_el=0;
      float E_mu=0;
      float E_neutHad=0;
      //float recojet_jetChargeE_kappa_0_25=0;
      //float recojet_jetChargeE_kappa_0_50=0;
      //float recojet_jetChargeE_kappa_0_75=0;
      //float recojet_jetChargeE_kappa_1_00=0;
      //float recojet_jetChargeE_kappa_0_10=0;
      //float recojet_jetChargeE_kappa_0_15=0;
      //float recojet_jetChargeE_kappa_0_20=0;
      //float recojet_jetChargeE_kappa_0_30=0;
      
      //float recojet_jetChargePt_kappa_0_25=0;
      //float recojet_jetChargePt_kappa_0_50=0;
      //float recojet_jetChargePt_kappa_0_75=0;
      //float recojet_jetChargePt_kappa_1_00=0;
      //float recojet_jetChargePt_kappa_0_10=0;
      //float recojet_jetChargePt_kappa_0_15=0;
      //float recojet_jetChargePt_kappa_0_20=0;
      //float recojet_jetChargePt_kappa_0_30=0;

      //float recojet_jetPProj_kappa_0_25=0;
      //float recojet_jetPProj_kappa_0_50=0;
      //float recojet_jetPProj_kappa_0_75=0;
      //float recojet_jetPProj_kappa_1_00=0;
      //float recojet_jetPProj_kappa_0_10=0;
      //float recojet_jetPProj_kappa_0_15=0;
      //float recojet_jetPProj_kappa_0_20=0;
      //float recojet_jetPProj_kappa_0_30=0;
      
      //float recojet_jetChargePProj_kappa_0_25=0;
      //float recojet_jetChargePProj_kappa_0_50=0;
      //float recojet_jetChargePProj_kappa_0_75=0;
      //float recojet_jetChargePProj_kappa_1_00=0;
      //float recojet_jetChargePProj_kappa_0_10=0;
      //float recojet_jetChargePProj_kappa_0_15=0;
      //float recojet_jetChargePProj_kappa_0_20=0;
      //float recojet_jetChargePProj_kappa_0_30=0;

      for(unsigned int i=0;i<_csreco->constituents((*recojetIt)).size();i++){
	//recojet-particles are prepared as reconstructed particle
	//reco Type is the original PDGID, charge is filled as well
	ReconstructedParticle* recopart = static_cast<ReconstructedParticle*> (particleRECOJetIn->getElementAt(_csreco->constituents((*recojetIt))[i].user_index()));
	//if(sqrt(pow(recopart->getMomentum()[0],2)+pow(recopart->getMomentum()[1],2))>m_cutPtTrackMin){
	//recojet_jetChargeE_kappa_0_25+=recopart->getCharge()*pow(recopart->getEnergy()/recojetIt->E(),0.25);
	//recojet_jetChargeE_kappa_0_50+=recopart->getCharge()*pow(recopart->getEnergy()/recojetIt->E(),0.50);
	  //recojet_jetChargeE_kappa_0_75+=recopart->getCharge()*pow(recopart->getEnergy()/recojetIt->E(),0.75);
	  //recojet_jetChargeE_kappa_1_00+=recopart->getCharge()*pow(recopart->getEnergy()/recojetIt->E(),1.00);
	  //recojet_jetChargeE_kappa_0_10+=recopart->getCharge()*pow(recopart->getEnergy()/recojetIt->E(),0.10);
	  //recojet_jetChargeE_kappa_0_15+=recopart->getCharge()*pow(recopart->getEnergy()/recojetIt->E(),0.15);
	  //recojet_jetChargeE_kappa_0_20+=recopart->getCharge()*pow(recopart->getEnergy()/recojetIt->E(),0.20);
	  //recojet_jetChargeE_kappa_0_30+=recopart->getCharge()*pow(recopart->getEnergy()/recojetIt->E(),0.30);
	  //float recopartPt=sqrt(pow(recopart->getMomentum()[0],2)+pow(recopart->getMomentum()[1],2));
	  //recojet_jetChargePt_kappa_0_25+=recopart->getCharge()*pow(recopartPt,0.25);
	  //recojet_jetChargePt_kappa_0_50+=recopart->getCharge()*pow(recopartPt,0.50);
	  //recojet_jetChargePt_kappa_0_75+=recopart->getCharge()*pow(recopartPt,0.75);
	  //recojet_jetChargePt_kappa_1_00+=recopart->getCharge()*pow(recopartPt,1.00);
	  //recojet_jetChargePt_kappa_0_10+=recopart->getCharge()*pow(recopartPt,0.10);
	  //recojet_jetChargePt_kappa_0_15+=recopart->getCharge()*pow(recopartPt,0.15);
	  //recojet_jetChargePt_kappa_0_20+=recopart->getCharge()*pow(recopartPt,0.20);
	  //recojet_jetChargePt_kappa_0_30+=recopart->getCharge()*pow(recopartPt,0.30);
	  //recojet_jetChargePProj_kappa_0_25+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.25);
	  //recojet_jetChargePProj_kappa_0_50+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.50);
	  //recojet_jetChargePProj_kappa_0_75+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.75);
	  //recojet_jetChargePProj_kappa_1_00+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),1.00);
	  //recojet_jetChargePProj_kappa_0_10+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.10);
	  //recojet_jetChargePProj_kappa_0_15+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.15);
	  //recojet_jetChargePProj_kappa_0_20+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.20);
	//recojet_jetChargePProj_kappa_0_30+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.30);
	//}
	//for denominator all particles are of relevance
	//recojet_jetPProj_kappa_0_25+=pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.25);
	//recojet_jetPProj_kappa_0_50+=pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.50);
	//recojet_jetPProj_kappa_0_75+=pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.75);
	//recojet_jetPProj_kappa_1_00+=pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),1.00);
	//recojet_jetPProj_kappa_0_10+=pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.10);
	//recojet_jetPProj_kappa_0_15+=pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.15);
	//recojet_jetPProj_kappa_0_20+=pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.20);
	//recojet_jetPProj_kappa_0_30+=pow(fabs(recopart->getMomentum()[0]*recojetIt->px()+recopart->getMomentum()[1]*recojetIt->py()+recopart->getMomentum()[2]*recojetIt->pz()),0.30);
	if(recopart->getCharge()==0){
	  if(recopart->getType()==22){
	    num_ph+=1;
	    E_ph+=recopart->getEnergy();
	  }else{
	    num_neutHad+=1;
	    E_neutHad+=recopart->getEnergy();
	  }
	}else{
	  if(abs(recopart->getType())==11){
	    E_el+=recopart->getEnergy();
	  }else if(abs(recopart->getType())==13){
	    E_mu+=recopart->getEnergy();
	  }else{//should be identified as pion on RECO level
	    num_pi+=1;
	    E_pi+=recopart->getEnergy();
	    if(sqrt(pow(recopart->getMomentum()[0],2)+pow(recopart->getMomentum()[1],2))>m_cutPtTrackMin){
	      num_pi_trackPtMin+=1;
	      E_pi_trackPtMin+=recopart->getEnergy();
	    }
	  }
	}
      }
      m_recojet_NCH->push_back(num_pi);
      m_recojet_NCH_trackPtMin->push_back(num_pi_trackPtMin);
      m_recojet_NPh->push_back(num_ph);
      m_recojet_NNH->push_back(num_neutHad);
      
      m_recojet_CHFraction->push_back(E_pi/recojetIt->E());
      m_recojet_CHFraction_trackPtMin->push_back(E_pi_trackPtMin/recojetIt->E());
      m_recojet_PhFraction->push_back(E_ph/recojetIt->E());
      m_recojet_ElFraction->push_back(E_el/recojetIt->E());
      m_recojet_MuFraction->push_back(E_mu/recojetIt->E());
      m_recojet_NHFraction->push_back(E_neutHad/recojetIt->E());
      //m_recojet_jetChargeE_kappa_0_25->push_back(recojet_jetChargeE_kappa_0_25);
      //m_recojet_jetChargeE_kappa_0_50->push_back(recojet_jetChargeE_kappa_0_50);
      //m_recojet_jetChargeE_kappa_0_75->push_back(recojet_jetChargeE_kappa_0_75);
      //m_recojet_jetChargeE_kappa_1_00->push_back(recojet_jetChargeE_kappa_1_00);
      //m_recojet_jetChargeE_kappa_0_10->push_back(recojet_jetChargeE_kappa_0_10);
      //m_recojet_jetChargeE_kappa_0_15->push_back(recojet_jetChargeE_kappa_0_15);
      //m_recojet_jetChargeE_kappa_0_20->push_back(recojet_jetChargeE_kappa_0_20);
      //m_recojet_jetChargeE_kappa_0_30->push_back(recojet_jetChargeE_kappa_0_30);
      //m_recojet_jetChargePt_kappa_0_25->push_back(recojet_jetChargePt_kappa_0_25/pow(recojetIt->pt(),0.25));
      //m_recojet_jetChargePt_kappa_0_50->push_back(recojet_jetChargePt_kappa_0_50/pow(recojetIt->pt(),0.50));
      //m_recojet_jetChargePt_kappa_0_75->push_back(recojet_jetChargePt_kappa_0_75/pow(recojetIt->pt(),0.75));
      //m_recojet_jetChargePt_kappa_1_00->push_back(recojet_jetChargePt_kappa_1_00/pow(recojetIt->pt(),1.00));
      //m_recojet_jetChargePt_kappa_0_10->push_back(recojet_jetChargePt_kappa_0_10/pow(recojetIt->pt(),0.10));
      //m_recojet_jetChargePt_kappa_0_15->push_back(recojet_jetChargePt_kappa_0_15/pow(recojetIt->pt(),0.15));
      //m_recojet_jetChargePt_kappa_0_20->push_back(recojet_jetChargePt_kappa_0_20/pow(recojetIt->pt(),0.20));
      //m_recojet_jetChargePt_kappa_0_30->push_back(recojet_jetChargePt_kappa_0_30/pow(recojetIt->pt(),0.30));

      //m_recojet_jetChargePProj_kappa_0_25->push_back(recojet_jetChargePProj_kappa_0_25/recojet_jetPProj_kappa_0_25);
      //m_recojet_jetChargePProj_kappa_0_50->push_back(recojet_jetChargePProj_kappa_0_50/recojet_jetPProj_kappa_0_50);
      //m_recojet_jetChargePProj_kappa_0_75->push_back(recojet_jetChargePProj_kappa_0_75/recojet_jetPProj_kappa_0_75);
      //m_recojet_jetChargePProj_kappa_1_00->push_back(recojet_jetChargePProj_kappa_1_00/recojet_jetPProj_kappa_1_00);
      //m_recojet_jetChargePProj_kappa_0_10->push_back(recojet_jetChargePProj_kappa_0_10/recojet_jetPProj_kappa_0_10);
      //m_recojet_jetChargePProj_kappa_0_15->push_back(recojet_jetChargePProj_kappa_0_15/recojet_jetPProj_kappa_0_15);
      //m_recojet_jetChargePProj_kappa_0_20->push_back(recojet_jetChargePProj_kappa_0_20/recojet_jetPProj_kappa_0_20);
      //m_recojet_jetChargePProj_kappa_0_30->push_back(recojet_jetChargePProj_kappa_0_30/recojet_jetPProj_kappa_0_30);

      m_recojet_nsubjettiness1->push_back((_nSubJettiness1_ptR->result(*recojetIt)));
      m_recojet_nsubjettiness1_lrz->push_back((_nSubJettiness1_lorentz->result(*recojetIt)));

      m_recojet_beta1_ECorr2->push_back(_energycorr2_beta1->result(*recojetIt));
      m_recojet_beta1_ECorr3->push_back(_energycorr3_beta1->result(*recojetIt));
      m_recojet_beta1_N2->push_back(_energycorrN2_beta1->result(*recojetIt));
      m_recojet_beta1_N3->push_back(_energycorrN3_beta1->result(*recojetIt));
      m_recojet_beta1_C2->push_back(_energycorrC2_beta1->result(*recojetIt));
      if(_energycorr3_beta1->result(*recojetIt)!=0){
	m_recojet_beta1_C3->push_back(_energycorr4_beta1->result(*recojetIt)*_energycorr2_beta1->result(*recojetIt)/pow(_energycorr3_beta1->result(*recojetIt),2));
      }else{
	m_recojet_beta1_C3->push_back(0.);
      }
      m_recojet_beta1_D2->push_back(_energycorrD2_beta1->result(*recojetIt));
      m_recojet_beta1_ECorr2_E_theta->push_back(_energycorr2_Etheta_beta1->result(*recojetIt));
      m_recojet_beta1_ECorr3_E_theta->push_back(_energycorr3_Etheta_beta1->result(*recojetIt));
      m_recojet_beta1_N2_E_theta->push_back(_energycorrN2_Etheta_beta1->result(*recojetIt));
      m_recojet_beta1_N3_E_theta->push_back(_energycorrN3_Etheta_beta1->result(*recojetIt));
      m_recojet_beta1_C2_E_theta->push_back(_energycorrC2_Etheta_beta1->result(*recojetIt));
      if(_energycorr3_Etheta_beta1->result(*recojetIt)!=0){
	m_recojet_beta1_C3_E_theta->push_back(_energycorr4_Etheta_beta1->result(*recojetIt)*_energycorr2_Etheta_beta1->result(*recojetIt)/pow(_energycorr3_Etheta_beta1->result(*recojetIt),2));
      }else{
	m_recojet_beta1_C3_E_theta->push_back(0.);
      }
      m_recojet_beta1_D2_E_theta->push_back(_energycorrD2_Etheta_beta1->result(*recojetIt));

      m_recojet_beta2_ECorr2->push_back(_energycorr2_beta2->result(*recojetIt));
      m_recojet_beta2_ECorr3->push_back(_energycorr3_beta2->result(*recojetIt));
      m_recojet_beta2_N2->push_back(_energycorrN2_beta2->result(*recojetIt));
      m_recojet_beta2_N3->push_back(_energycorrN3_beta2->result(*recojetIt));
      m_recojet_beta2_C2->push_back(_energycorrC2_beta2->result(*recojetIt));
      if(_energycorr3_beta2->result(*recojetIt)!=0){
	m_recojet_beta2_C3->push_back(_energycorr4_beta2->result(*recojetIt)*_energycorr2_beta2->result(*recojetIt)/pow(_energycorr3_beta2->result(*recojetIt),2));
      }else{
	m_recojet_beta2_C3->push_back(0.);
      }
      m_recojet_beta2_D2->push_back(_energycorrD2_beta2->result(*recojetIt));
      m_recojet_beta2_ECorr2_E_theta->push_back(_energycorr2_Etheta_beta2->result(*recojetIt));
      m_recojet_beta2_ECorr3_E_theta->push_back(_energycorr3_Etheta_beta2->result(*recojetIt));
      m_recojet_beta2_N2_E_theta->push_back(_energycorrN2_Etheta_beta2->result(*recojetIt));
      m_recojet_beta2_N3_E_theta->push_back(_energycorrN3_Etheta_beta2->result(*recojetIt));
      m_recojet_beta2_C2_E_theta->push_back(_energycorrC2_Etheta_beta2->result(*recojetIt));
      if(_energycorr3_Etheta_beta2->result(*recojetIt)!=0){
	m_recojet_beta2_C3_E_theta->push_back(_energycorr4_Etheta_beta2->result(*recojetIt)*_energycorr2_Etheta_beta2->result(*recojetIt)/pow(_energycorr3_Etheta_beta2->result(*recojetIt),2));
      }else{
	m_recojet_beta2_C3_E_theta->push_back(0.);
      }
      m_recojet_beta2_D2_E_theta->push_back(_energycorrD2_Etheta_beta2->result(*recojetIt));

      m_recojet_beta0_5_ECorr2->push_back(_energycorr2_beta0_5->result(*recojetIt));
      m_recojet_beta0_5_ECorr3->push_back(_energycorr3_beta0_5->result(*recojetIt));
      m_recojet_beta0_5_N2->push_back(_energycorrN2_beta0_5->result(*recojetIt));
      m_recojet_beta0_5_N3->push_back(_energycorrN3_beta0_5->result(*recojetIt));
      m_recojet_beta0_5_C2->push_back(_energycorrC2_beta0_5->result(*recojetIt));
      if(_energycorr3_beta0_5->result(*recojetIt)!=0){
	m_recojet_beta0_5_C3->push_back(_energycorr4_beta0_5->result(*recojetIt)*_energycorr2_beta0_5->result(*recojetIt)/pow(_energycorr3_beta0_5->result(*recojetIt),2));
      }else{
	m_recojet_beta0_5_C3->push_back(0.);
      }
      m_recojet_beta0_5_D2->push_back(_energycorrD2_beta0_5->result(*recojetIt));
      m_recojet_beta0_5_ECorr2_E_theta->push_back(_energycorr2_Etheta_beta0_5->result(*recojetIt));
      m_recojet_beta0_5_ECorr3_E_theta->push_back(_energycorr3_Etheta_beta0_5->result(*recojetIt));
      m_recojet_beta0_5_N2_E_theta->push_back(_energycorrN2_Etheta_beta0_5->result(*recojetIt));
      m_recojet_beta0_5_N3_E_theta->push_back(_energycorrN3_Etheta_beta0_5->result(*recojetIt));
      m_recojet_beta0_5_C2_E_theta->push_back(_energycorrC2_Etheta_beta0_5->result(*recojetIt));
      if(_energycorr3_Etheta_beta0_5->result(*recojetIt)!=0){
	m_recojet_beta0_5_C3_E_theta->push_back(_energycorr4_Etheta_beta0_5->result(*recojetIt)*_energycorr2_Etheta_beta0_5->result(*recojetIt)/pow(_energycorr3_Etheta_beta0_5->result(*recojetIt),2));
      }else{
	m_recojet_beta0_5_C3_E_theta->push_back(0.);
      }
      m_recojet_beta0_5_D2_E_theta->push_back(_energycorrD2_Etheta_beta0_5->result(*recojetIt));
      
      //decompose into requested number of subjets: 2 jets here, thus check if this can be done --> not possible for one 
      if(_csreco->constituents((*recojetIt)).size()>1){
	//std::cout<<"HZAnalyzer cluster subjet"<<std::endl;
	m_recojet_nsubjettiness2->push_back((_nSubJettiness2_ptR->result(*recojetIt)));
	m_recojet_nsubjettiness2_lrz->push_back((_nSubJettiness2_lorentz->result(*recojetIt)));
	std::vector<fastjet::PseudoJet> subjets = _csreco->exclusive_subjets(*recojetIt, 2);
	m_recojet_dij_21->push_back(recojetIt->exclusive_subdmerge(1));
	m_recojet_dij_32->push_back(recojetIt->exclusive_subdmerge(2));
	m_recojet_dij_43->push_back(recojetIt->exclusive_subdmerge(3));
	m_recojet_dij_21_max->push_back(recojetIt->exclusive_subdmerge_max(1));
	m_recojet_dij_32_max->push_back(recojetIt->exclusive_subdmerge_max(2));
	m_recojet_dij_43_max->push_back(recojetIt->exclusive_subdmerge_max(3));

	int num_pi_sj=0;
	int num_pi_sj_trackPtMin=0;
	int num_ph_sj=0;
	int num_neutHad_sj=0;
	
	float E_pi_sj=0;
	float E_pi_sj_trackPtMin=0;
	float E_ph_sj=0;
	float E_el_sj=0;
	float E_mu_sj=0;
	float E_neutHad_sj=0;
	m_recojet_subjet_E->push_back(subjets[0].E());
	m_recojet_subjet_Px->push_back(subjets[0].px());
	m_recojet_subjet_Py->push_back(subjets[0].py());
	m_recojet_subjet_Pz->push_back(subjets[0].pz());
	m_recojet_subjet_jetindex->push_back(recojet_index_);
	
	float recojet_subjet_jetChargeE_kappa_0_25=0;
	float recojet_subjet_jetChargeE_kappa_0_50=0;
	//float recojet_subjet_jetChargeE_kappa_0_75=0;
	//float recojet_subjet_jetChargeE_kappa_1_00=0;
	//float recojet_subjet_jetChargeE_kappa_0_10=0;
	//float recojet_subjet_jetChargeE_kappa_0_15=0;
	float recojet_subjet_jetChargeE_kappa_0_20=0;
	float recojet_subjet_jetChargeE_kappa_0_30=0;
	
	float recojet_subjet_jetChargePt_kappa_0_25=0;
	float recojet_subjet_jetChargePt_kappa_0_50=0;
	//float recojet_subjet_jetChargePt_kappa_0_75=0;
	//float recojet_subjet_jetChargePt_kappa_1_00=0;
	//float recojet_subjet_jetChargePt_kappa_0_10=0;
	//float recojet_subjet_jetChargePt_kappa_0_15=0;
	float recojet_subjet_jetChargePt_kappa_0_20=0;
	float recojet_subjet_jetChargePt_kappa_0_30=0;

	//float recojet_subjet_jetPProj_kappa_0_25=0;
	//float recojet_subjet_jetPProj_kappa_0_50=0;
	//float recojet_subjet_jetPProj_kappa_0_75=0;
	//float recojet_subjet_jetPProj_kappa_1_00=0;
	//float recojet_subjet_jetPProj_kappa_0_10=0;
	//float recojet_subjet_jetPProj_kappa_0_15=0;
	//float recojet_subjet_jetPProj_kappa_0_20=0;
	//float recojet_subjet_jetPProj_kappa_0_30=0;
 
	//float recojet_subjet_jetChargePProj_kappa_0_25=0;
	//float recojet_subjet_jetChargePProj_kappa_0_50=0;
	//float recojet_subjet_jetChargePProj_kappa_0_75=0;
	//float recojet_subjet_jetChargePProj_kappa_1_00=0;
	//float recojet_subjet_jetChargePProj_kappa_0_10=0;
	//float recojet_subjet_jetChargePProj_kappa_0_15=0;
	//float recojet_subjet_jetChargePProj_kappa_0_20=0;
	//float recojet_subjet_jetChargePProj_kappa_0_30=0;

	double E_total_subjet=0;

	for(unsigned int i=0;i<_csreco->constituents(subjets[0]).size();i++){
	  ReconstructedParticle* recopart = static_cast<ReconstructedParticle*> (particleRECOJetIn->getElementAt(_csreco->constituents(subjets[0])[i].user_index()));
	  //if(rfjet_has_duplicate){
	  //std::cout<<"duplicate particle "<<i<<" in sj0 of jet "<<recojet_index_<<" E/px/py/pz "<<recopart->getEnergy()<<"/"<<recopart->getMomentum()[0]<<"/"<<recopart->getMomentum()[1]<<"/"<<recopart->getMomentum()[2]<<" PDG "<<recopart->getType()<<std::endl;
	  //}
	  if(sqrt(pow(recopart->getMomentum()[0],2)+pow(recopart->getMomentum()[1],2))>m_cutPtTrackMin){
	    recojet_subjet_jetChargeE_kappa_0_25+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[0].E(),0.25);
	    recojet_subjet_jetChargeE_kappa_0_50+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[0].E(),0.50);
	    //recojet_subjet_jetChargeE_kappa_0_75+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[0].E(),0.75);
	    //recojet_subjet_jetChargeE_kappa_1_00+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[0].E(),1.00);
	    //recojet_subjet_jetChargeE_kappa_0_10+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[0].E(),0.10);
	    //recojet_subjet_jetChargeE_kappa_0_15+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[0].E(),0.15);
	    recojet_subjet_jetChargeE_kappa_0_20+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[0].E(),0.20);
	    recojet_subjet_jetChargeE_kappa_0_30+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[0].E(),0.30);
	    float recopartPt=sqrt(pow(recopart->getMomentum()[0],2)+pow(recopart->getMomentum()[1],2));
	    recojet_subjet_jetChargePt_kappa_0_25+=recopart->getCharge()*pow(recopartPt,0.25);
	    recojet_subjet_jetChargePt_kappa_0_50+=recopart->getCharge()*pow(recopartPt,0.50);
	    //recojet_subjet_jetChargePt_kappa_0_75+=recopart->getCharge()*pow(recopartPt,0.75);
	    //recojet_subjet_jetChargePt_kappa_1_00+=recopart->getCharge()*pow(recopartPt,1.00);
	    //recojet_subjet_jetChargePt_kappa_0_10+=recopart->getCharge()*pow(recopartPt,0.10);
	    //recojet_subjet_jetChargePt_kappa_0_15+=recopart->getCharge()*pow(recopartPt,0.15);
	    recojet_subjet_jetChargePt_kappa_0_20+=recopart->getCharge()*pow(recopartPt,0.20);
	    recojet_subjet_jetChargePt_kappa_0_30+=recopart->getCharge()*pow(recopartPt,0.30);
	    //recojet_subjet_jetChargePProj_kappa_0_25+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.25);
	    //recojet_subjet_jetChargePProj_kappa_0_50+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.50);
	    //recojet_subjet_jetChargePProj_kappa_0_75+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.75);
	    //recojet_subjet_jetChargePProj_kappa_1_00+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),1.00);
	    //recojet_subjet_jetChargePProj_kappa_0_10+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.10);
	    //recojet_subjet_jetChargePProj_kappa_0_15+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.15);
	    //recojet_subjet_jetChargePProj_kappa_0_20+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.20);
	    //recojet_subjet_jetChargePProj_kappa_0_30+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.30);
	  }
	  //recojet_subjet_jetPProj_kappa_0_25+=pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.25);
	  //recojet_subjet_jetPProj_kappa_0_50+=pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.50);
	  //recojet_subjet_jetPProj_kappa_0_75+=pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.75);
	  //recojet_subjet_jetPProj_kappa_1_00+=pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),1.00);
	  //recojet_subjet_jetPProj_kappa_0_10+=pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.10);
	  //recojet_subjet_jetPProj_kappa_0_15+=pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.15);
	  //recojet_subjet_jetPProj_kappa_0_20+=pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.20);
	  //recojet_subjet_jetPProj_kappa_0_30+=pow(fabs(recopart->getMomentum()[0]*subjets[0].px()+recopart->getMomentum()[1]*subjets[0].py()+recopart->getMomentum()[2]*subjets[0].pz()),0.30);
	  E_total_subjet+=recopart->getEnergy();
	  if(recopart->getCharge()==0){
	    if(recopart->getType()==22){
	      num_ph_sj+=1;
	      E_ph_sj+=recopart->getEnergy();
	    }else{
	      num_neutHad_sj+=1;
	      E_neutHad_sj+=recopart->getEnergy();
	    }
	  }else{
	    if(abs(recopart->getType())==11){
	      E_el_sj+=recopart->getEnergy();
	    }else if(abs(recopart->getType())==13){
	      E_mu_sj+=recopart->getEnergy();
	    }else{//should be identified as pion on RECO level
	      num_pi_sj+=1;
	      E_pi_sj+=recopart->getEnergy();
	      if(sqrt(pow(recopart->getMomentum()[0],2)+pow(recopart->getMomentum()[1],2))>m_cutPtTrackMin){
		num_pi_sj_trackPtMin+=1;
		E_pi_sj_trackPtMin+=recopart->getEnergy();
	      }
	    }
	  }
	}
	if(E_total_subjet!=subjets[0].E()){
	  std::cout<<"energy totally screwed up "<<E_total_subjet<<"/"<<subjets[0].E()<< "/"<< (E_total_subjet-subjets[0].E())/subjets[0].E()<<" calc from subjet ratio "<<std::endl;
	}
	m_recojet_subjet_NCH->push_back(num_pi_sj);
	m_recojet_subjet_NCH_trackPtMin->push_back(num_pi_sj_trackPtMin);
	m_recojet_subjet_NPh->push_back(num_ph_sj);
	m_recojet_subjet_NNH->push_back(num_neutHad_sj);
	
	m_recojet_subjet_CHFraction->push_back(E_pi_sj/subjets[0].E());
	m_recojet_subjet_CHFraction_trackPtMin->push_back(E_pi_sj_trackPtMin/subjets[0].E());
	m_recojet_subjet_PhFraction->push_back(E_ph_sj/subjets[0].E());
	m_recojet_subjet_ElFraction->push_back(E_el_sj/subjets[0].E());
	m_recojet_subjet_MuFraction->push_back(E_mu_sj/subjets[0].E());
	m_recojet_subjet_NHFraction->push_back(E_neutHad_sj/subjets[0].E());

	m_recojet_subjet_jetChargeE_kappa_0_25->push_back(recojet_subjet_jetChargeE_kappa_0_25);
	m_recojet_subjet_jetChargeE_kappa_0_50->push_back(recojet_subjet_jetChargeE_kappa_0_50);
	//m_recojet_subjet_jetChargeE_kappa_0_75->push_back(recojet_subjet_jetChargeE_kappa_0_75);
	//m_recojet_subjet_jetChargeE_kappa_1_00->push_back(recojet_subjet_jetChargeE_kappa_1_00);
	//m_recojet_subjet_jetChargeE_kappa_0_10->push_back(recojet_subjet_jetChargeE_kappa_0_10);
	//m_recojet_subjet_jetChargeE_kappa_0_15->push_back(recojet_subjet_jetChargeE_kappa_0_15);
	m_recojet_subjet_jetChargeE_kappa_0_20->push_back(recojet_subjet_jetChargeE_kappa_0_20);
	m_recojet_subjet_jetChargeE_kappa_0_30->push_back(recojet_subjet_jetChargeE_kappa_0_30);
	m_recojet_subjet_jetChargePt_kappa_0_25->push_back(recojet_subjet_jetChargePt_kappa_0_25/pow(subjets[0].pt(),0.25));
	m_recojet_subjet_jetChargePt_kappa_0_50->push_back(recojet_subjet_jetChargePt_kappa_0_50/pow(subjets[0].pt(),0.50));
	//m_recojet_subjet_jetChargePt_kappa_0_75->push_back(recojet_subjet_jetChargePt_kappa_0_75/pow(subjets[0].pt(),0.75));
	//m_recojet_subjet_jetChargePt_kappa_1_00->push_back(recojet_subjet_jetChargePt_kappa_1_00/pow(subjets[0].pt(),1.00));
	//m_recojet_subjet_jetChargePt_kappa_0_10->push_back(recojet_subjet_jetChargePt_kappa_0_10/pow(subjets[0].pt(),0.10));
	//m_recojet_subjet_jetChargePt_kappa_0_15->push_back(recojet_subjet_jetChargePt_kappa_0_15/pow(subjets[0].pt(),0.15));
	m_recojet_subjet_jetChargePt_kappa_0_20->push_back(recojet_subjet_jetChargePt_kappa_0_20/pow(subjets[0].pt(),0.20));
	m_recojet_subjet_jetChargePt_kappa_0_30->push_back(recojet_subjet_jetChargePt_kappa_0_30/pow(subjets[0].pt(),0.30));

	//m_recojet_subjet_jetChargePProj_kappa_0_25->push_back(recojet_subjet_jetChargePProj_kappa_0_25/recojet_subjet_jetPProj_kappa_0_25);
	//m_recojet_subjet_jetChargePProj_kappa_0_50->push_back(recojet_subjet_jetChargePProj_kappa_0_50/recojet_subjet_jetPProj_kappa_0_50);
	//m_recojet_subjet_jetChargePProj_kappa_0_75->push_back(recojet_subjet_jetChargePProj_kappa_0_75/recojet_subjet_jetPProj_kappa_0_75);
	//m_recojet_subjet_jetChargePProj_kappa_1_00->push_back(recojet_subjet_jetChargePProj_kappa_1_00/recojet_subjet_jetPProj_kappa_1_00);
	//m_recojet_subjet_jetChargePProj_kappa_0_10->push_back(recojet_subjet_jetChargePProj_kappa_0_10/recojet_subjet_jetPProj_kappa_0_10);
	//m_recojet_subjet_jetChargePProj_kappa_0_15->push_back(recojet_subjet_jetChargePProj_kappa_0_15/recojet_subjet_jetPProj_kappa_0_15);
	//m_recojet_subjet_jetChargePProj_kappa_0_20->push_back(recojet_subjet_jetChargePProj_kappa_0_20/recojet_subjet_jetPProj_kappa_0_20);
	//m_recojet_subjet_jetChargePProj_kappa_0_30->push_back(recojet_subjet_jetChargePProj_kappa_0_30/recojet_subjet_jetPProj_kappa_0_30);
	
	num_pi_sj=0;
	num_pi_sj_trackPtMin=0;
	num_ph_sj=0;
	num_neutHad_sj=0;

	E_pi_sj=0;
	E_pi_sj_trackPtMin=0;
	E_ph_sj=0;
	E_el_sj=0;
	E_mu_sj=0;
	E_neutHad_sj=0;
	m_recojet_subjet_E->push_back(subjets[1].E());
	m_recojet_subjet_Px->push_back(subjets[1].px());
	m_recojet_subjet_Py->push_back(subjets[1].py());
	m_recojet_subjet_Pz->push_back(subjets[1].pz());
	m_recojet_subjet_jetindex->push_back(recojet_index_);

	
	recojet_subjet_jetChargeE_kappa_0_25=0;
	recojet_subjet_jetChargeE_kappa_0_50=0;
	//recojet_subjet_jetChargeE_kappa_0_75=0;
	//recojet_subjet_jetChargeE_kappa_1_00=0;
	//recojet_subjet_jetChargeE_kappa_0_10=0;
	//recojet_subjet_jetChargeE_kappa_0_15=0;
	recojet_subjet_jetChargeE_kappa_0_20=0;
	recojet_subjet_jetChargeE_kappa_0_30=0;
	
	recojet_subjet_jetChargePt_kappa_0_25=0;
	recojet_subjet_jetChargePt_kappa_0_50=0;
	//recojet_subjet_jetChargePt_kappa_0_75=0;
	//recojet_subjet_jetChargePt_kappa_1_00=0;
	//recojet_subjet_jetChargePt_kappa_0_10=0;
	//recojet_subjet_jetChargePt_kappa_0_15=0;
	recojet_subjet_jetChargePt_kappa_0_20=0;
	recojet_subjet_jetChargePt_kappa_0_30=0;

	//recojet_subjet_jetPProj_kappa_0_25=0;
	//recojet_subjet_jetPProj_kappa_0_50=0;
	//recojet_subjet_jetPProj_kappa_0_75=0;
	//recojet_subjet_jetPProj_kappa_1_00=0;
	//recojet_subjet_jetPProj_kappa_0_10=0;
	//recojet_subjet_jetPProj_kappa_0_15=0;
	//recojet_subjet_jetPProj_kappa_0_20=0;
	//recojet_subjet_jetPProj_kappa_0_30=0;
	
	//recojet_subjet_jetChargePProj_kappa_0_25=0;
	//recojet_subjet_jetChargePProj_kappa_0_50=0;
	//recojet_subjet_jetChargePProj_kappa_0_75=0;
	//recojet_subjet_jetChargePProj_kappa_1_00=0;
	//recojet_subjet_jetChargePProj_kappa_0_10=0;
	//recojet_subjet_jetChargePProj_kappa_0_15=0;
	//recojet_subjet_jetChargePProj_kappa_0_20=0;
	//recojet_subjet_jetChargePProj_kappa_0_30=0;


	for(unsigned int i=0;i<_csreco->constituents(subjets[1]).size();i++){
	  ReconstructedParticle* recopart = static_cast<ReconstructedParticle*> (particleRECOJetIn->getElementAt(_csreco->constituents(subjets[1])[i].user_index()));
	  //if(rfjet_has_duplicate){
	  //std::cout<<"duplicate  particle "<<i<<" in sj1 of jet "<<recojet_index_<<" E/px/py/pz "<<recopart->getEnergy()<<"/"<<recopart->getMomentum()[0]<<"/"<<recopart->getMomentum()[1]<<"/"<<recopart->getMomentum()[2]<<" PDG "<<recopart->getType()<<std::endl;
	  //}
	  if(sqrt(pow(recopart->getMomentum()[0],2)+pow(recopart->getMomentum()[1],2))>m_cutPtTrackMin){
	    recojet_subjet_jetChargeE_kappa_0_25+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[1].E(),0.25);
	    recojet_subjet_jetChargeE_kappa_0_50+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[1].E(),0.50);
	    //recojet_subjet_jetChargeE_kappa_0_75+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[1].E(),0.75);
	    //recojet_subjet_jetChargeE_kappa_1_00+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[1].E(),1.00);
	    //recojet_subjet_jetChargeE_kappa_0_10+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[1].E(),0.10);
	    //recojet_subjet_jetChargeE_kappa_0_15+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[1].E(),0.15);
	    recojet_subjet_jetChargeE_kappa_0_20+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[1].E(),0.20);
	    recojet_subjet_jetChargeE_kappa_0_30+=recopart->getCharge()*pow(recopart->getEnergy()/subjets[1].E(),0.30);
	    float recopartPt=sqrt(pow(recopart->getMomentum()[0],2)+pow(recopart->getMomentum()[1],2));
	    recojet_subjet_jetChargePt_kappa_0_25+=recopart->getCharge()*pow(recopartPt,0.25);
	    recojet_subjet_jetChargePt_kappa_0_50+=recopart->getCharge()*pow(recopartPt,0.50);
	    //recojet_subjet_jetChargePt_kappa_0_75+=recopart->getCharge()*pow(recopartPt,0.75);
	    //recojet_subjet_jetChargePt_kappa_1_00+=recopart->getCharge()*pow(recopartPt,1.00);
	    //recojet_subjet_jetChargePt_kappa_0_10+=recopart->getCharge()*pow(recopartPt,0.10);
	    //recojet_subjet_jetChargePt_kappa_0_15+=recopart->getCharge()*pow(recopartPt,0.15);
	    recojet_subjet_jetChargePt_kappa_0_20+=recopart->getCharge()*pow(recopartPt,0.20);
	    recojet_subjet_jetChargePt_kappa_0_30+=recopart->getCharge()*pow(recopartPt,0.30);
	    //recojet_subjet_jetChargePProj_kappa_0_25+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.25);
	    //recojet_subjet_jetChargePProj_kappa_0_50+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.50);
	    //recojet_subjet_jetChargePProj_kappa_0_75+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.75);
	    //recojet_subjet_jetChargePProj_kappa_1_00+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),1.00);
	    //recojet_subjet_jetChargePProj_kappa_0_10+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.10);
	    //recojet_subjet_jetChargePProj_kappa_0_15+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.15);
	    //recojet_subjet_jetChargePProj_kappa_0_20+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.20);
	    //recojet_subjet_jetChargePProj_kappa_0_30+=recopart->getCharge()*pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.30);
	  }
	  //recojet_subjet_jetPProj_kappa_0_25+=pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.25);
	  //recojet_subjet_jetPProj_kappa_0_50+=pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.50);
	  //recojet_subjet_jetPProj_kappa_0_75+=pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.75);
	  //recojet_subjet_jetPProj_kappa_1_00+=pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),1.00);
	  //recojet_subjet_jetPProj_kappa_0_10+=pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.10);
	  //recojet_subjet_jetPProj_kappa_0_15+=pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.15);
	  //recojet_subjet_jetPProj_kappa_0_20+=pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.20);
	  //recojet_subjet_jetPProj_kappa_0_30+=pow(fabs(recopart->getMomentum()[0]*subjets[1].px()+recopart->getMomentum()[1]*subjets[1].py()+recopart->getMomentum()[2]*subjets[1].pz()),0.30);
	  if(recopart->getCharge()==0){
	    if(recopart->getType()==22){
	      num_ph_sj+=1;
	      E_ph_sj+=recopart->getEnergy();
	    }else{
	      num_neutHad_sj+=1;
	      E_neutHad_sj+=recopart->getEnergy();
	    }
	  }else{
	    if(abs(recopart->getType())==11){
	      E_el_sj+=recopart->getEnergy();
	    }else if(abs(recopart->getType())==13){
	      E_mu_sj+=recopart->getEnergy();
	    }else{//should be identified as pion on RECO level
	      num_pi_sj+=1;
	      E_pi_sj+=recopart->getEnergy();
	      if(sqrt(pow(recopart->getMomentum()[0],2)+pow(recopart->getMomentum()[1],2))>m_cutPtTrackMin){
		num_pi_sj_trackPtMin+=1;
		E_pi_sj_trackPtMin+=recopart->getEnergy();
	      }
	    }
	  }
	}
	m_recojet_subjet_NCH->push_back(num_pi_sj);
	m_recojet_subjet_NCH_trackPtMin->push_back(num_pi_sj_trackPtMin);
	m_recojet_subjet_NPh->push_back(num_ph_sj);
	m_recojet_subjet_NNH->push_back(num_neutHad_sj);
	
	m_recojet_subjet_CHFraction->push_back(E_pi_sj/subjets[1].E());
	m_recojet_subjet_CHFraction_trackPtMin->push_back(E_pi_sj_trackPtMin/subjets[1].E());
	m_recojet_subjet_PhFraction->push_back(E_ph_sj/subjets[1].E());
	m_recojet_subjet_ElFraction->push_back(E_el_sj/subjets[1].E());
	m_recojet_subjet_MuFraction->push_back(E_mu_sj/subjets[1].E());
	m_recojet_subjet_NHFraction->push_back(E_neutHad_sj/subjets[1].E());

	m_recojet_subjet_jetChargeE_kappa_0_25->push_back(recojet_subjet_jetChargeE_kappa_0_25);
	m_recojet_subjet_jetChargeE_kappa_0_50->push_back(recojet_subjet_jetChargeE_kappa_0_50);
	//m_recojet_subjet_jetChargeE_kappa_0_75->push_back(recojet_subjet_jetChargeE_kappa_0_75);
	//m_recojet_subjet_jetChargeE_kappa_1_00->push_back(recojet_subjet_jetChargeE_kappa_1_00);
	//m_recojet_subjet_jetChargeE_kappa_0_10->push_back(recojet_subjet_jetChargeE_kappa_0_10);
	//m_recojet_subjet_jetChargeE_kappa_0_15->push_back(recojet_subjet_jetChargeE_kappa_0_15);
	m_recojet_subjet_jetChargeE_kappa_0_20->push_back(recojet_subjet_jetChargeE_kappa_0_20);
	m_recojet_subjet_jetChargeE_kappa_0_30->push_back(recojet_subjet_jetChargeE_kappa_0_30);

	m_recojet_subjet_jetChargePt_kappa_0_25->push_back(recojet_subjet_jetChargePt_kappa_0_25/pow(subjets[1].pt(),0.25));
	m_recojet_subjet_jetChargePt_kappa_0_50->push_back(recojet_subjet_jetChargePt_kappa_0_50/pow(subjets[1].pt(),0.50));
	//m_recojet_subjet_jetChargePt_kappa_0_75->push_back(recojet_subjet_jetChargePt_kappa_0_75/pow(subjets[1].pt(),0.75));
	//m_recojet_subjet_jetChargePt_kappa_1_00->push_back(recojet_subjet_jetChargePt_kappa_1_00/pow(subjets[1].pt(),1.00));
	//m_recojet_subjet_jetChargePt_kappa_0_10->push_back(recojet_subjet_jetChargePt_kappa_0_10/pow(subjets[1].pt(),0.10));
	//m_recojet_subjet_jetChargePt_kappa_0_15->push_back(recojet_subjet_jetChargePt_kappa_0_15/pow(subjets[1].pt(),0.15));
	m_recojet_subjet_jetChargePt_kappa_0_20->push_back(recojet_subjet_jetChargePt_kappa_0_20/pow(subjets[1].pt(),0.20));
	m_recojet_subjet_jetChargePt_kappa_0_30->push_back(recojet_subjet_jetChargePt_kappa_0_30/pow(subjets[1].pt(),0.30));

	//m_recojet_subjet_jetChargePProj_kappa_0_25->push_back(recojet_subjet_jetChargePProj_kappa_0_25/recojet_subjet_jetPProj_kappa_0_25);
	//m_recojet_subjet_jetChargePProj_kappa_0_50->push_back(recojet_subjet_jetChargePProj_kappa_0_50/recojet_subjet_jetPProj_kappa_0_50);
	//m_recojet_subjet_jetChargePProj_kappa_0_75->push_back(recojet_subjet_jetChargePProj_kappa_0_75/recojet_subjet_jetPProj_kappa_0_75);
	//m_recojet_subjet_jetChargePProj_kappa_1_00->push_back(recojet_subjet_jetChargePProj_kappa_1_00/recojet_subjet_jetPProj_kappa_1_00);
	//m_recojet_subjet_jetChargePProj_kappa_0_10->push_back(recojet_subjet_jetChargePProj_kappa_0_10/recojet_subjet_jetPProj_kappa_0_10);
	//m_recojet_subjet_jetChargePProj_kappa_0_15->push_back(recojet_subjet_jetChargePProj_kappa_0_15/recojet_subjet_jetPProj_kappa_0_15);
	//m_recojet_subjet_jetChargePProj_kappa_0_20->push_back(recojet_subjet_jetChargePProj_kappa_0_20/recojet_subjet_jetPProj_kappa_0_20);
	//m_recojet_subjet_jetChargePProj_kappa_0_30->push_back(recojet_subjet_jetChargePProj_kappa_0_30/recojet_subjet_jetPProj_kappa_0_30);
	
	if(_csreco->constituents((*recojetIt)).size()>2){
	  m_recojet_nsubjettiness3->push_back((_nSubJettiness3_ptR->result(*recojetIt)));
	  m_recojet_nsubjettiness3_lrz->push_back((_nSubJettiness3_lorentz->result(*recojetIt)));
	}else{
	  m_recojet_nsubjettiness3->push_back(-1);
	  m_recojet_nsubjettiness3_lrz->push_back(-1);
	}
      }else{
	m_recojet_dij_21->push_back(0);
	m_recojet_dij_32->push_back(0);
	m_recojet_dij_43->push_back(0);
	m_recojet_dij_21_max->push_back(0);
	m_recojet_dij_32_max->push_back(0);
	m_recojet_dij_43_max->push_back(0);
	m_recojet_nsubjettiness2->push_back(-1);
	m_recojet_nsubjettiness2_lrz->push_back(-1);
	m_recojet_nsubjettiness3->push_back(-1);
	m_recojet_nsubjettiness3_lrz->push_back(-1);
      }
      recojet_index_+=1;
    }
    delete _csreco;
  }
  /*
  if(rfjet_has_duplicate){
    if( recoColl != 0 ){
      int nRECOP = recoColl->getNumberOfElements();
      for(int m=0;m<nRECOP;m++){
	ReconstructedParticle *recop = static_cast<ReconstructedParticle*>(recoColl->getElementAt(m));
	std::cout<<"evt "<<evt->getRunNumber()<<" particle "<<m<<" E/px/py/pz "<<recop->getEnergy()<<"/"<<recop->getMomentum()[0]<<"/"<<recop->getMomentum()[1]<<"/"<<recop->getMomentum()[2]<<" PDG "<<recop->getType()<<std::endl;
	if(recop->getTracks().size()>0){
	  for(unsigned int t=0;t<recop->getTracks().size();t++){
	    std::cout<<"track "<<t<<" z0/tanlambda "<<recop->getTracks()[t]->getZ0()<<"/"<<recop->getTracks()[t]->getTanLambda()<<std::endl;
	  }
	}
      }
    }
  }
  */
  //if((recorefSubJet0Jet0Coll->getNumberOfElements()+recorefSubJet1Jet0Coll->getNumberOfElements()+recorefSubJet0Jet1Coll->getNumberOfElements()+recorefSubJet1Jet1Coll->getNumberOfElements())!=(m_recojet_subjet_E->size()) || m_recojet_subjet_rfj_sj_E->size()!=(m_recojet_subjet_E->size()) ){
  //std::cout<<"rf subjets collection/subjet rfj size/not the same as subjet size "<<(recorefSubJet0Jet0Coll->getNumberOfElements()+recorefSubJet1Jet0Coll->getNumberOfElements()+recorefSubJet0Jet1Coll->getNumberOfElements()+recorefSubJet1Jet1Coll->getNumberOfElements())<<"/"<<m_recojet_subjet_rfj_sj_E->size()<<"/"<<m_recojet_subjet_E->size()<<std::endl;
  //}
  
  //missing are TAU's and isolated particles AND THEN WE ARE READY TO GO --> only two days I know sucks but oh well.
  m_outputTree->Fill();
}

void HZAnalyzer::check(LCEvent*){
}

void HZAnalyzer::end(){
  
  delete vlcpl;
  delete _jetAlgoType;
  delete _vlcAxes;
  delete _unnormalizedMeasure_ptR;
  delete _unnormalizedMeasure_Lorentz;
  delete _nSubJettiness1_ptR;
  delete _nSubJettiness1_lorentz;
  delete _nSubJettiness2_ptR;
  delete _nSubJettiness2_lorentz;
  delete _nSubJettiness3_ptR;
  delete _nSubJettiness3_lorentz;

  delete _energycorr2_beta1;
  delete _energycorr2_Etheta_beta1;
  delete _energycorr3_beta1;
  delete _energycorr3_Etheta_beta1;
  delete _energycorr4_beta1;
  delete _energycorr4_Etheta_beta1;
  delete _energycorrC2_beta1;
  delete _energycorrC2_Etheta_beta1;
  delete _energycorrD2_beta1;
  delete _energycorrD2_Etheta_beta1;
  delete _energycorrN2_beta1;
  delete _energycorrN2_Etheta_beta1;
  delete _energycorrN3_beta1;
  delete _energycorrN3_Etheta_beta1;

  delete _energycorr2_beta2;
  delete _energycorr2_Etheta_beta2;
  delete _energycorr3_beta2;
  delete _energycorr3_Etheta_beta2;
  delete _energycorr4_beta2;
  delete _energycorr4_Etheta_beta2;
  delete _energycorrC2_beta2;
  delete _energycorrC2_Etheta_beta2;
  delete _energycorrD2_beta2;
  delete _energycorrD2_Etheta_beta2;
  delete _energycorrN2_beta2;
  delete _energycorrN2_Etheta_beta2;
  delete _energycorrN3_beta2;
  delete _energycorrN3_Etheta_beta2;

  delete _energycorr2_beta0_5;
  delete _energycorr2_Etheta_beta0_5;
  delete _energycorr3_beta0_5;
  delete _energycorr3_Etheta_beta0_5;
  delete _energycorr4_beta0_5;
  delete _energycorr4_Etheta_beta0_5;
  delete _energycorrC2_beta0_5;
  delete _energycorrC2_Etheta_beta0_5;
  delete _energycorrD2_beta0_5;
  delete _energycorrD2_Etheta_beta0_5;
  delete _energycorrN2_beta0_5;
  delete _energycorrN2_Etheta_beta0_5;
  delete _energycorrN3_beta0_5;
  delete _energycorrN3_Etheta_beta0_5;
 

  m_rootFile->cd();  
 
  m_rootFile->Write();
  m_rootFile->Close();
}


void HZAnalyzer::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    
    streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

void HZAnalyzer::fillStableDaughterSet(MCParticle* mcp, std::set<MCParticle*> &stableDaughterSet){
  if(mcp->getGeneratorStatus()==1){
    stableDaughterSet.insert(mcp);
  }else if (mcp->getGeneratorStatus()==0){
    return;
  }
  for(unsigned int d=0;d<mcp->getDaughters().size();d++){
    fillStableDaughterSet(mcp->getDaughters()[d], stableDaughterSet);
  }
}
