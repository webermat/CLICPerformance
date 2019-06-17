#include "PhotonStudy.h"
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetectorSelector.h"

#include "UTIL/LCRelationNavigator.h"

#include "TFile.h"

using namespace lcio ;
using namespace marlin ;


PhotonStudy aPhotonStudy;

PhotonStudy::PhotonStudy() : Processor("PhotonStudy") {
    
    // modify processor description
    _description = "PhotonStudy calculates properties of calorimeter showers" ;
   

    registerInputCollection( LCIO::MCPARTICLE,
                            "MCParticleCollectionName",
                            "Name of the MCParticle input collection",
                            m_inputMCParticleCollection,
                            std::string("MCParticle"));
    
    registerProcessorParameter( "OutputRootFileName",
                                "ROOT File name to collect plots",
                                m_rootFileName,
                                std::string("PhotonStudy.root"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RECOParticleCollectionName",
                            "Name of the RECOParticle input collection",
                            m_inputRECOParticleCollection,
                            std::string("PandoraPFOs"));

    registerProcessorParameter(
			       "JetEMin" , 
			       "minimum Energy cut for RECO jets",
			       m_jetEMin ,
			       float(10.)
			       );

    registerProcessorParameter(
			       "ignoreGen" , 
			       "save generator information",
			       m_ignoreGen,
			       bool(false)			       
			       );

    registerProcessorParameter(
			       "fillMEInfo" , 
			       "save matrix element information",
			       m_fillMEInfo,
			       bool(false)			       
			       );


    registerProcessorParameter(
			       "reducedOutput" , 
			       "reduced output less energy bins etc",
			       m_reducedOutput,
			       bool(false)			       
			       );


    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoExclJetCollection" , 
			     "Name of the ReReco RecoJet collection"  ,
			     m_jetExclColName,
			     std::string("JetOut_excl")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoInclJetCollection" , 
			     "Name of the ReReco RecoJet collection"  ,
			     m_jetInclColName,
			     std::string("JetOut_incl")
			     );

}  


void PhotonStudy::init() {

  m_rootFile = new TFile(m_rootFileName.c_str(),"RECREATE");
  
    m_outputTree = new TTree("showerData","showerData");

    if(!m_ignoreGen){
      m_truePi0Energy = new std::vector<float>();
      m_truePi0_Px = new std::vector<float>();
      m_truePi0_Py = new std::vector<float>();
      m_truePi0_Pz = new std::vector<float>();
      m_truePi0_CosTheta = new std::vector<float>();
      m_truePi0_Phi = new std::vector<float>();
      m_truePi0_DRPh01 = new std::vector<float>();
      m_truePi0_CosAngPh01 = new std::vector<float>();
      m_truePi0_AngPh01 = new std::vector<float>();
      m_truePh0Energy = new std::vector<float>();
      m_truePh0_Px = new std::vector<float>();
      m_truePh0_Py = new std::vector<float>();
      m_truePh0_Pz = new std::vector<float>();
      m_truePh0_CosTheta = new std::vector<float>();
      m_truePh0_Phi = new std::vector<float>();
      m_truePh0_DRMin = new std::vector<float>();
      m_truePh0_DRMinPDGID = new std::vector<float>();
      m_truePh0_DRMinE = new std::vector<float>();
      m_truePh0_CosAngMax = new std::vector<float>();
      m_truePh0_CosAngMaxPDGID = new std::vector<float>();
      m_truePh0_CosAngMaxE = new std::vector<float>();
      m_truePh0_CosAng0995E = new std::vector<float>();
      m_truePh0_DR01E = new std::vector<float>();
      m_truePh1Energy = new std::vector<float>();
      m_truePh1_Px = new std::vector<float>();
      m_truePh1_Py = new std::vector<float>();
      m_truePh1_Pz = new std::vector<float>();
      m_truePh1_CosTheta = new std::vector<float>();
      m_truePh1_Phi = new std::vector<float>();
      m_truePh1_DRMin = new std::vector<float>();
      m_truePh1_DRMinPDGID = new std::vector<float>();
      m_truePh1_DRMinE = new std::vector<float>();
      m_truePh1_CosAngMax = new std::vector<float>();
      m_truePh1_CosAngMaxPDGID = new std::vector<float>();
      m_truePh1_CosAngMaxE = new std::vector<float>();
      m_truePh1_CosAng0995E = new std::vector<float>();
      m_truePh1_DR01E = new std::vector<float>();
    }
    if(m_fillMEInfo){
      m_trueME_E=new std::vector<float>();
      m_trueME_Px=new std::vector<float>();
      m_trueME_Py=new std::vector<float>();
      m_trueME_Pz=new std::vector<float>();
      m_trueME_PDGID=new std::vector<int>();
    }
    m_recoPhEnergy = new std::vector<float>();
    m_recoPh_Px = new std::vector<float>();
    m_recoPh_Py = new std::vector<float>();
    m_recoPh_Pz = new std::vector<float>();
    m_recoPh_CosTheta = new std::vector<float>();
    m_recoPh_Theta = new std::vector<float>();
    m_recoPh_Phi = new std::vector<float>();
    m_recoPh_PDGID = new std::vector<float>();
    m_recoPh_E_EB = new std::vector<float>();
    m_recoPh_E_EE = new std::vector<float>();
    m_recoPh_E_EO = new std::vector<float>();
    m_recoPh_E_HB = new std::vector<float>();
    m_recoPh_E_HE = new std::vector<float>();
    m_recoPh_E_HO = new std::vector<float>();
    m_recoPh_firstLayerECAL = new std::vector<int>();
    m_recoPh_lastLayerECAL = new std::vector<int>();
    m_recoPh_nhitsEB = new std::vector<int>();
    m_recoPh_nhitsEE = new std::vector<int>();
    m_recoPh_nhitsEO = new std::vector<int>();
    m_recoPh_firstLayerHCAL = new std::vector<int>();
    m_recoPh_lastLayerHCAL = new std::vector<int>();
    m_recoPh_nhitsHB = new std::vector<int>();
    m_recoPh_nhitsHE = new std::vector<int>();
    m_recoPh_nhitsHO = new std::vector<int>();

    m_recoPh_DRMin_Ph = new std::vector<float>();
    m_recoPh_DRMin_E_Ph = new std::vector<float>();
    m_recoPh_CosAngMax_Ph = new std::vector<float>();
    m_recoPh_CosAngMax_E_Ph = new std::vector<float>();

    m_recoPh_DRMin = new std::vector<float>();
    m_recoPh_DRMin_PDGID = new std::vector<float>();
    m_recoPh_DRMin_E = new std::vector<float>();
    m_recoPh_CosAngMax = new std::vector<float>();
    m_recoPh_CosAngMax_PDGID = new std::vector<float>();
    m_recoPh_CosAngMax_E = new std::vector<float>();
    m_recoPh_CosAng0995_E = new std::vector<float>();
    m_recoPh_DR01_E = new std::vector<float>();
    m_recoPi0Cand_DRMin_M = new std::vector<float>();
    m_recoPi0Cand_CosAngMax_M = new std::vector<float>();


    m_jet_excl_Px        = new std::vector<float>(); 
    m_jet_excl_Py        = new std::vector<float>(); 
    m_jet_excl_Pz        = new std::vector<float>(); 
    m_jet_excl_E         = new std::vector<float>(); 
    m_jet_excl_Phi       = new std::vector<float>(); 
    m_jet_excl_CosTheta  = new std::vector<float>(); 
    m_jet_excl_piE       = new std::vector<float>();
    m_jet_excl_phE       = new std::vector<float>();
    m_jet_excl_elE       = new std::vector<float>();
    m_jet_excl_muE       = new std::vector<float>();
    m_jet_excl_nE        = new std::vector<float>();
    m_jet_excl_neutElseE = new std::vector<float>();
    m_jet_excl_chElseE   = new std::vector<float>();
    m_jet_excl_neutMult  = new std::vector<int>();
    m_jet_excl_chMult    = new std::vector<int>();

    m_jet_incl_Px        = new std::vector<float>(); 
    m_jet_incl_Py        = new std::vector<float>(); 
    m_jet_incl_Pz        = new std::vector<float>(); 
    m_jet_incl_E         = new std::vector<float>(); 
    m_jet_incl_Phi       = new std::vector<float>(); 
    m_jet_incl_CosTheta  = new std::vector<float>(); 
    m_jet_incl_piE       = new std::vector<float>();
    m_jet_incl_phE       = new std::vector<float>();
    m_jet_incl_elE       = new std::vector<float>();
    m_jet_incl_muE       = new std::vector<float>();
    m_jet_incl_nE        = new std::vector<float>();
    m_jet_incl_neutElseE = new std::vector<float>();
    m_jet_incl_chElseE   = new std::vector<float>();
    m_jet_incl_neutMult  = new std::vector<int>();
    m_jet_incl_chMult    = new std::vector<int>();

    eventcount=0;
  
    // Print the initial parameters
    printParameters() ;

    // Reset counters
    m_runNumber = 0 ;
    m_eventNumber = 0 ;

    m_E_trueNeut=0;
    m_px_trueNeut=0;
    m_py_trueNeut=0;
    m_pz_trueNeut=0;

    m_E_trueZ2=0;
    m_px_trueZ2=0;
    m_py_trueZ2=0;
    m_pz_trueZ2=0;


    m_E_totPFO=0;
    m_E_totPi=0;
    m_E_totPh=0;
    m_E_totE=0;
    m_E_totMu=0;
    m_E_totN=0;

    m_px_totPFO=0;
    m_px_totPi=0;
    m_px_totPh=0;
    m_px_totE=0;
    m_px_totMu=0;
    m_px_totN=0;

    m_py_totPFO=0;
    m_py_totPi=0;
    m_py_totPh=0;
    m_py_totE=0;
    m_py_totMu=0;
    m_py_totN=0;

    m_pz_totPFO=0;
    m_pz_totPi=0;
    m_pz_totPh=0;
    m_pz_totE=0;
    m_pz_totMu=0;
    m_pz_totN=0;

    m_E_totPFO_0_95=0;
    m_E_totPi_0_95=0;
    m_E_totPh_0_95=0;
    m_E_totE_0_95=0;
    m_E_totMu_0_95=0;
    m_E_totN_0_95=0;
    m_E_totPFO_0_70=0;
    m_E_totPi_0_70=0;
    m_E_totPh_0_70=0;
    m_E_totE_0_70=0;
    m_E_totMu_0_70=0;
    m_E_totN_0_70=0;
 
    m_px_true_totAll=0;
    m_px_true_totInv=0;
    m_px_true_totPi=0;
    m_px_true_totPh=0;
    m_px_true_totK0L=0;
    m_px_true_totE=0;
    m_px_true_totMu=0;
    m_px_true_totN=0;
    m_px_true_totK=0;
    m_px_true_totP=0;

    m_py_true_totAll=0;
    m_py_true_totInv=0;
    m_py_true_totPi=0;
    m_py_true_totPh=0;
    m_py_true_totK0L=0;
    m_py_true_totE=0;
    m_py_true_totMu=0;
    m_py_true_totN=0;
    m_py_true_totK=0;
    m_py_true_totP=0;

    m_pz_true_totAll=0;
    m_pz_true_totInv=0;
    m_pz_true_totPi=0;
    m_pz_true_totPh=0;
    m_pz_true_totK0L=0;
    m_pz_true_totE=0;
    m_pz_true_totMu=0;
    m_pz_true_totN=0;
    m_pz_true_totK=0;
    m_pz_true_totP=0;

    m_E_true_totAll=0;
    m_E_true_totInv=0;
    m_E_true_totPi=0;
    m_E_true_totPh=0;
    m_E_true_totK0L=0;
    m_E_true_totE=0;
    m_E_true_totMu=0;
    m_E_true_totN=0;
    m_E_true_totK=0;
    m_E_true_totP=0;
    m_E_true_totOtherCH=0;
    m_E_true_totOtherNeut=0;
    m_E_true_totAll_0_95=0;
    m_E_true_totInv_0_95=0;
    m_E_true_totPi_0_95=0;
    m_E_true_totPh_0_95=0;
    m_E_true_totK0L_0_95=0;
    m_E_true_totE_0_95=0;
    m_E_true_totMu_0_95=0;
    m_E_true_totN_0_95=0;
    m_E_true_totK_0_95=0;
    m_E_true_totP_0_95=0;
    m_E_true_totOtherCH_0_95=0;
    m_E_true_totOtherNeut_0_95=0;
    m_E_true_totAll_0_70=0;
    m_E_true_totInv_0_70=0;
    m_E_true_totPi_0_70=0;
    m_E_true_totPh_0_70=0;
    m_E_true_totK0L_0_70=0;
    m_E_true_totE_0_70=0;
    m_E_true_totMu_0_70=0;
    m_E_true_totN_0_70=0;
    m_E_true_totK_0_70=0;
    m_E_true_totP_0_70=0;
    m_E_true_totOtherCH_0_70=0;
    m_E_true_totOtherNeut_0_70=0;

    m_E_0_2_totPFO=0;
    m_E_0_2_totPi=0;
    m_E_0_2_totPh=0;
    m_E_0_2_totE=0;
    m_E_0_2_totMu=0;
    m_E_0_2_totN=0;
    m_E_0_2_totPFO_0_95=0;
    m_E_0_2_totPi_0_95=0;
    m_E_0_2_totPh_0_95=0;
    m_E_0_2_totE_0_95=0;
    m_E_0_2_totMu_0_95=0;
    m_E_0_2_totN_0_95=0;
    m_E_0_2_totPFO_0_70=0;
    m_E_0_2_totPi_0_70=0;
    m_E_0_2_totPh_0_70=0;
    m_E_0_2_totE_0_70=0;
    m_E_0_2_totMu_0_70=0;
    m_E_0_2_totN_0_70=0;


    m_E_true_0_2_totAll=0;
    m_E_true_0_2_totInv=0;
    m_E_true_0_2_totPi=0;
    m_E_true_0_2_totPh=0;
    m_E_true_0_2_totK0L=0;
    m_E_true_0_2_totE=0;
    m_E_true_0_2_totMu=0;
    m_E_true_0_2_totN=0;
    m_E_true_0_2_totK=0;
    m_E_true_0_2_totP=0;
    m_E_true_0_2_totOtherCH=0;
    m_E_true_0_2_totOtherNeut=0;
    m_E_true_0_2_totAll_0_95=0;
    m_E_true_0_2_totInv_0_95=0;
    m_E_true_0_2_totPi_0_95=0;
    m_E_true_0_2_totPh_0_95=0;
    m_E_true_0_2_totK0L_0_95=0;
    m_E_true_0_2_totE_0_95=0;
    m_E_true_0_2_totMu_0_95=0;
    m_E_true_0_2_totN_0_95=0;
    m_E_true_0_2_totK_0_95=0;
    m_E_true_0_2_totP_0_95=0;
    m_E_true_0_2_totOtherCH_0_95=0;
    m_E_true_0_2_totOtherNeut_0_95=0;
    m_E_true_0_2_totAll_0_70=0;
    m_E_true_0_2_totInv_0_70=0;
    m_E_true_0_2_totPi_0_70=0;
    m_E_true_0_2_totPh_0_70=0;
    m_E_true_0_2_totK0L_0_70=0;
    m_E_true_0_2_totE_0_70=0;
    m_E_true_0_2_totMu_0_70=0;
    m_E_true_0_2_totN_0_70=0;
    m_E_true_0_2_totK_0_70=0;
    m_E_true_0_2_totP_0_70=0;
    m_E_true_0_2_totOtherCH_0_70=0;
    m_E_true_0_2_totOtherNeut_0_70=0;

    m_E_2_10_totPFO=0;
    m_E_2_10_totPi=0;
    m_E_2_10_totPh=0;
    m_E_2_10_totE=0;
    m_E_2_10_totMu=0;
    m_E_2_10_totN=0;
    m_E_2_10_totPFO_0_95=0;
    m_E_2_10_totPi_0_95=0;
    m_E_2_10_totPh_0_95=0;
    m_E_2_10_totE_0_95=0;
    m_E_2_10_totMu_0_95=0;
    m_E_2_10_totN_0_95=0;
    m_E_2_10_totPFO_0_70=0;
    m_E_2_10_totPi_0_70=0;
    m_E_2_10_totPh_0_70=0;
    m_E_2_10_totE_0_70=0;
    m_E_2_10_totMu_0_70=0;
    m_E_2_10_totN_0_70=0;

    m_E_true_2_10_totAll=0;
    m_E_true_2_10_totInv=0;
    m_E_true_2_10_totPi=0;
    m_E_true_2_10_totPh=0;
    m_E_true_2_10_totK0L=0;
    m_E_true_2_10_totE=0;
    m_E_true_2_10_totMu=0;
    m_E_true_2_10_totN=0;
    m_E_true_2_10_totK=0;
    m_E_true_2_10_totP=0;
    m_E_true_2_10_totOtherCH=0;
    m_E_true_2_10_totOtherNeut=0;
    m_E_true_2_10_totAll_0_95=0;
    m_E_true_2_10_totInv_0_95=0;
    m_E_true_2_10_totPi_0_95=0;
    m_E_true_2_10_totPh_0_95=0;
    m_E_true_2_10_totK0L_0_95=0;
    m_E_true_2_10_totE_0_95=0;
    m_E_true_2_10_totMu_0_95=0;
    m_E_true_2_10_totN_0_95=0;
    m_E_true_2_10_totK_0_95=0;
    m_E_true_2_10_totP_0_95=0;
    m_E_true_2_10_totOtherCH_0_95=0;
    m_E_true_2_10_totOtherNeut_0_95=0;
    m_E_true_2_10_totAll_0_70=0;
    m_E_true_2_10_totInv_0_70=0;
    m_E_true_2_10_totPi_0_70=0;
    m_E_true_2_10_totPh_0_70=0;
    m_E_true_2_10_totK0L_0_70=0;
    m_E_true_2_10_totE_0_70=0;
    m_E_true_2_10_totMu_0_70=0;
    m_E_true_2_10_totN_0_70=0;
    m_E_true_2_10_totK_0_70=0;
    m_E_true_2_10_totP_0_70=0;
    m_E_true_2_10_totOtherCH_0_70=0;
    m_E_true_2_10_totOtherNeut_0_70=0;

    m_E_10_50_totPFO=0;
    m_E_10_50_totPi=0;
    m_E_10_50_totPh=0;
    m_E_10_50_totE=0;
    m_E_10_50_totMu=0;
    m_E_10_50_totN=0;
    m_E_10_50_totPFO_0_95=0;
    m_E_10_50_totPi_0_95=0;
    m_E_10_50_totPh_0_95=0;
    m_E_10_50_totE_0_95=0;
    m_E_10_50_totMu_0_95=0;
    m_E_10_50_totN_0_95=0;
    m_E_10_50_totPFO_0_70=0;
    m_E_10_50_totPi_0_70=0;
    m_E_10_50_totPh_0_70=0;
    m_E_10_50_totE_0_70=0;
    m_E_10_50_totMu_0_70=0;
    m_E_10_50_totN_0_70=0;

    m_E_true_10_50_totAll=0;
    m_E_true_10_50_totInv=0;
    m_E_true_10_50_totPi=0;
    m_E_true_10_50_totPh=0;
    m_E_true_10_50_totK0L=0;
    m_E_true_10_50_totE=0;
    m_E_true_10_50_totMu=0;
    m_E_true_10_50_totN=0;
    m_E_true_10_50_totK=0;
    m_E_true_10_50_totP=0;
    m_E_true_10_50_totOtherCH=0;
    m_E_true_10_50_totOtherNeut=0;
    m_E_true_10_50_totAll_0_95=0;
    m_E_true_10_50_totInv_0_95=0;
    m_E_true_10_50_totPi_0_95=0;
    m_E_true_10_50_totPh_0_95=0;
    m_E_true_10_50_totK0L_0_95=0;
    m_E_true_10_50_totE_0_95=0;
    m_E_true_10_50_totMu_0_95=0;
    m_E_true_10_50_totN_0_95=0;
    m_E_true_10_50_totK_0_95=0;
    m_E_true_10_50_totP_0_95=0;
    m_E_true_10_50_totOtherCH_0_95=0;
    m_E_true_10_50_totOtherNeut_0_95=0;
    m_E_true_10_50_totAll_0_70=0;
    m_E_true_10_50_totInv_0_70=0;
    m_E_true_10_50_totPi_0_70=0;
    m_E_true_10_50_totPh_0_70=0;
    m_E_true_10_50_totK0L_0_70=0;
    m_E_true_10_50_totE_0_70=0;
    m_E_true_10_50_totMu_0_70=0;
    m_E_true_10_50_totN_0_70=0;
    m_E_true_10_50_totK_0_70=0;
    m_E_true_10_50_totP_0_70=0;
    m_E_true_10_50_totOtherCH_0_70=0;
    m_E_true_10_50_totOtherNeut_0_70=0;

    m_E_50_totPFO=0;
    m_E_50_totPi=0;
    m_E_50_totPh=0;
    m_E_50_totE=0;
    m_E_50_totMu=0;
    m_E_50_totN=0;
    m_E_50_totPFO_0_95=0;
    m_E_50_totPi_0_95=0;
    m_E_50_totPh_0_95=0;
    m_E_50_totE_0_95=0;
    m_E_50_totMu_0_95=0;
    m_E_50_totN_0_95=0;
    m_E_50_totPFO_0_70=0;
    m_E_50_totPi_0_70=0;
    m_E_50_totPh_0_70=0;
    m_E_50_totE_0_70=0;
    m_E_50_totMu_0_70=0;
    m_E_50_totN_0_70=0;

    m_E_true_50_totAll=0;
    m_E_true_50_totInv=0;
    m_E_true_50_totPi=0;
   m_E_true_50_totPh=0;
    m_E_true_50_totK0L=0;
    m_E_true_50_totE=0;
    m_E_true_50_totMu=0;
    m_E_true_50_totN=0;
    m_E_true_50_totK=0;
    m_E_true_50_totP=0;
    m_E_true_50_totOtherCH=0;
    m_E_true_50_totOtherNeut=0;
    m_E_true_50_totAll_0_95=0;
    m_E_true_50_totInv_0_95=0;
    m_E_true_50_totPi_0_95=0;
    m_E_true_50_totPh_0_95=0;
    m_E_true_50_totK0L_0_95=0;
    m_E_true_50_totE_0_95=0;
    m_E_true_50_totMu_0_95=0;
    m_E_true_50_totN_0_95=0;
    m_E_true_50_totK_0_95=0;
    m_E_true_50_totP_0_95=0;
    m_E_true_50_totOtherCH_0_95=0;
    m_E_true_50_totOtherNeut_0_95=0;
    m_E_true_50_totAll_0_70=0;
    m_E_true_50_totInv_0_70=0;
    m_E_true_50_totPi_0_70=0;
    m_E_true_50_totPh_0_70=0;
    m_E_true_50_totK0L_0_70=0;
    m_E_true_50_totE_0_70=0;
    m_E_true_50_totMu_0_70=0;
    m_E_true_50_totN_0_70=0;
    m_E_true_50_totK_0_70=0;
    m_E_true_50_totP_0_70=0;
    m_E_true_50_totOtherCH_0_70=0;
    m_E_true_50_totOtherNeut_0_70=0;

    //multiplicities
    m_n_totPFO=0;
    m_n_totPi=0;
    m_n_totPh=0;
    m_n_totE=0;
    m_n_totMu=0;
    m_n_totN=0;
    m_n_totPFO_0_95=0;
    m_n_totPi_0_95=0;
    m_n_totPh_0_95=0;
    m_n_totE_0_95=0;
    m_n_totMu_0_95=0;
    m_n_totN_0_95=0;
    m_n_totPFO_0_70=0;
    m_n_totPi_0_70=0;
    m_n_totPh_0_70=0;
    m_n_totE_0_70=0;
    m_n_totMu_0_70=0;
    m_n_totN_0_70=0;

    m_n_true_totAll=0;
    m_n_true_totInv=0;
    m_n_true_totPi=0;
    m_n_true_totPh=0;
    m_n_true_totK0L=0;
    m_n_true_totE=0;
    m_n_true_totMu=0;
    m_n_true_totN=0;
    m_n_true_totK=0;
    m_n_true_totP=0;
    m_n_true_totOtherCH=0;
    m_n_true_totOtherNeut=0;
    m_n_true_totAll_0_95=0;
    m_n_true_totInv_0_95=0;
    m_n_true_totPi_0_95=0;
    m_n_true_totPh_0_95=0;
    m_n_true_totK0L_0_95=0;
    m_n_true_totE_0_95=0;
    m_n_true_totMu_0_95=0;
    m_n_true_totN_0_95=0;
    m_n_true_totK_0_95=0;
    m_n_true_totP_0_95=0;
    m_n_true_totOtherCH_0_95=0;
    m_n_true_totOtherNeut_0_95=0;
    m_n_true_totAll_0_70=0;
    m_n_true_totInv_0_70=0;
    m_n_true_totPi_0_70=0;
    m_n_true_totPh_0_70=0;
    m_n_true_totK0L_0_70=0;
    m_n_true_totE_0_70=0;
    m_n_true_totMu_0_70=0;
    m_n_true_totN_0_70=0;
    m_n_true_totK_0_70=0;
    m_n_true_totP_0_70=0;
    m_n_true_totOtherCH_0_70=0;
    m_n_true_totOtherNeut_0_70=0;

    m_n_0_2_totPFO=0;
    m_n_0_2_totPi=0;
    m_n_0_2_totPh=0;
    m_n_0_2_totE=0;
    m_n_0_2_totMu=0;
    m_n_0_2_totN=0;
    m_n_0_2_totPFO_0_95=0;
    m_n_0_2_totPi_0_95=0;
    m_n_0_2_totPh_0_95=0;
    m_n_0_2_totE_0_95=0;
    m_n_0_2_totMu_0_95=0;
    m_n_0_2_totN_0_95=0;
    m_n_0_2_totPFO_0_70=0;
    m_n_0_2_totPi_0_70=0;
    m_n_0_2_totPh_0_70=0;
    m_n_0_2_totE_0_70=0;
    m_n_0_2_totMu_0_70=0;
    m_n_0_2_totN_0_70=0;

    m_n_true_0_2_totAll=0;
    m_n_true_0_2_totInv=0;
    m_n_true_0_2_totPi=0;
    m_n_true_0_2_totPh=0;
    m_n_true_0_2_totE=0;
    m_n_true_0_2_totMu=0;
    m_n_true_0_2_totN=0;
    m_n_true_0_2_totK=0;
    m_n_true_0_2_totP=0;
    m_n_true_0_2_totOtherCH=0;
    m_n_true_0_2_totOtherNeut=0;
    m_n_true_0_2_totAll_0_95=0;
    m_n_true_0_2_totInv_0_95=0;
    m_n_true_0_2_totPi_0_95=0;
    m_n_true_0_2_totPh_0_95=0;
    m_n_true_0_2_totK0L_0_95=0;
    m_n_true_0_2_totE_0_95=0;
    m_n_true_0_2_totMu_0_95=0;
    m_n_true_0_2_totN_0_95=0;
    m_n_true_0_2_totK_0_95=0;
    m_n_true_0_2_totP_0_95=0;
    m_n_true_0_2_totOtherCH_0_95=0;
    m_n_true_0_2_totOtherNeut_0_95=0;
    m_n_true_0_2_totAll_0_70=0;
    m_n_true_0_2_totInv_0_70=0;
    m_n_true_0_2_totPi_0_70=0;
    m_n_true_0_2_totPh_0_70=0;
    m_n_true_0_2_totK0L_0_70=0;
    m_n_true_0_2_totE_0_70=0;
    m_n_true_0_2_totMu_0_70=0;
    m_n_true_0_2_totN_0_70=0;
    m_n_true_0_2_totK_0_70=0;
    m_n_true_0_2_totP_0_70=0;
    m_n_true_0_2_totOtherCH_0_70=0;
    m_n_true_0_2_totOtherNeut_0_70=0;

    m_n_2_10_totPFO=0;
    m_n_2_10_totPi=0;
    m_n_2_10_totPh=0;
    m_n_2_10_totE=0;
    m_n_2_10_totMu=0;
    m_n_2_10_totN=0;
    m_n_2_10_totPFO_0_95=0;
    m_n_2_10_totPi_0_95=0;
    m_n_2_10_totPh_0_95=0;
    m_n_2_10_totE_0_95=0;
    m_n_2_10_totMu_0_95=0;
    m_n_2_10_totN_0_95=0;
    m_n_2_10_totPFO_0_70=0;
    m_n_2_10_totPi_0_70=0;
    m_n_2_10_totPh_0_70=0;
    m_n_2_10_totE_0_70=0;
    m_n_2_10_totMu_0_70=0;
    m_n_2_10_totN_0_70=0;

    m_n_true_2_10_totAll=0;
    m_n_true_2_10_totInv=0;
    m_n_true_2_10_totPi=0;
    m_n_true_2_10_totPh=0;
    m_n_true_2_10_totK0L=0;
    m_n_true_2_10_totE=0;
    m_n_true_2_10_totMu=0;
    m_n_true_2_10_totN=0;
    m_n_true_2_10_totK=0;
    m_n_true_2_10_totP=0;
    m_n_true_2_10_totOtherCH=0;
    m_n_true_2_10_totOtherNeut=0;
    m_n_true_2_10_totAll_0_95=0;
    m_n_true_2_10_totInv_0_95=0;
    m_n_true_2_10_totPi_0_95=0;
    m_n_true_2_10_totPh_0_95=0;
    m_n_true_2_10_totK0L_0_95=0;
    m_n_true_2_10_totE_0_95=0;
    m_n_true_2_10_totMu_0_95=0;
    m_n_true_2_10_totN_0_95=0;
    m_n_true_2_10_totK_0_95=0;
    m_n_true_2_10_totP_0_95=0;
    m_n_true_2_10_totOtherCH_0_95=0;
    m_n_true_2_10_totOtherNeut_0_95=0;
    m_n_true_2_10_totAll_0_70=0;
    m_n_true_2_10_totInv_0_70=0;
    m_n_true_2_10_totPi_0_70=0;
    m_n_true_2_10_totPh_0_70=0;
    m_n_true_2_10_totK0L_0_70=0;
    m_n_true_2_10_totE_0_70=0;
    m_n_true_2_10_totMu_0_70=0;
    m_n_true_2_10_totN_0_70=0;
    m_n_true_2_10_totK_0_70=0;
    m_n_true_2_10_totP_0_70=0;
    m_n_true_2_10_totOtherCH_0_70=0;
    m_n_true_2_10_totOtherNeut_0_70=0;

    m_n_10_50_totPFO=0;
    m_n_10_50_totPi=0;
    m_n_10_50_totPh=0;
    m_n_10_50_totE=0;
    m_n_10_50_totMu=0;
    m_n_10_50_totN=0;
    m_n_10_50_totPFO_0_95=0;
    m_n_10_50_totPi_0_95=0;
    m_n_10_50_totPh_0_95=0;
    m_n_10_50_totE_0_95=0;
    m_n_10_50_totMu_0_95=0;
    m_n_10_50_totN_0_95=0;
    m_n_10_50_totPFO_0_70=0;
    m_n_10_50_totPi_0_70=0;
    m_n_10_50_totPh_0_70=0;
    m_n_10_50_totE_0_70=0;
    m_n_10_50_totMu_0_70=0;
    m_n_10_50_totN_0_70=0;

    m_n_true_10_50_totAll=0;
    m_n_true_10_50_totInv=0;
    m_n_true_10_50_totPi=0;
    m_n_true_10_50_totPh=0;
    m_n_true_10_50_totK0L=0;
    m_n_true_10_50_totE=0;
    m_n_true_10_50_totMu=0;
    m_n_true_10_50_totN=0;
    m_n_true_10_50_totK=0;
    m_n_true_10_50_totP=0;
    m_n_true_10_50_totOtherCH=0;
    m_n_true_10_50_totOtherNeut=0;
    m_n_true_10_50_totAll_0_95=0;
    m_n_true_10_50_totInv_0_95=0;
    m_n_true_10_50_totPi_0_95=0;
    m_n_true_10_50_totPh_0_95=0;
    m_n_true_10_50_totK0L_0_95=0;
    m_n_true_10_50_totE_0_95=0;
    m_n_true_10_50_totMu_0_95=0;
    m_n_true_10_50_totN_0_95=0;
    m_n_true_10_50_totK_0_95=0;
    m_n_true_10_50_totP_0_95=0;
    m_n_true_10_50_totOtherCH_0_95=0;
    m_n_true_10_50_totOtherNeut_0_95=0;
    m_n_true_10_50_totAll_0_70=0;
    m_n_true_10_50_totInv_0_70=0;
    m_n_true_10_50_totPi_0_70=0;
    m_n_true_10_50_totPh_0_70=0;
    m_n_true_10_50_totK0L_0_70=0;
    m_n_true_10_50_totE_0_70=0;
    m_n_true_10_50_totMu_0_70=0;
    m_n_true_10_50_totN_0_70=0;
    m_n_true_10_50_totK_0_70=0;
    m_n_true_10_50_totP_0_70=0;
    m_n_true_10_50_totOtherCH_0_70=0;
    m_n_true_10_50_totOtherNeut_0_70=0;

    m_n_50_totPFO=0;
    m_n_50_totPi=0;
    m_n_50_totPh=0;
    m_n_50_totE=0;
    m_n_50_totMu=0;
    m_n_50_totN=0;
    m_n_50_totPFO_0_95=0;
    m_n_50_totPi_0_95=0;
    m_n_50_totPh_0_95=0;
    m_n_50_totE_0_95=0;
    m_n_50_totMu_0_95=0;
    m_n_50_totN_0_95=0;
    m_n_50_totPFO_0_70=0;
    m_n_50_totPi_0_70=0;
    m_n_50_totPh_0_70=0;
    m_n_50_totE_0_70=0;
    m_n_50_totMu_0_70=0;
    m_n_50_totN_0_70=0;

    m_n_true_50_totAll=0;
    m_n_true_50_totInv=0;
    m_n_true_50_totPi=0;
    m_n_true_50_totPh=0;
    m_n_true_50_totK0L=0;
    m_n_true_50_totE=0;
    m_n_true_50_totMu=0;
    m_n_true_50_totN=0;
    m_n_true_50_totK=0;
    m_n_true_50_totP=0;
    m_n_true_50_totOtherCH=0;
    m_n_true_50_totOtherNeut=0;
    m_n_true_50_totAll_0_95=0;
    m_n_true_50_totInv_0_95=0;
    m_n_true_50_totPi_0_95=0;
    m_n_true_50_totPh_0_95=0;
    m_n_true_50_totK0L_0_95=0;
    m_n_true_50_totE_0_95=0;
    m_n_true_50_totMu_0_95=0;
    m_n_true_50_totN_0_95=0;
    m_n_true_50_totK_0_95=0;
    m_n_true_50_totP_0_95=0;
    m_n_true_50_totOtherCH_0_95=0;
    m_n_true_50_totOtherNeut_0_95=0;
    m_n_true_50_totAll_0_70=0;
    m_n_true_50_totInv_0_70=0;
    m_n_true_50_totPi_0_70=0;
    m_n_true_50_totPh_0_70=0;
    m_n_true_50_totK0L_0_70=0;
    m_n_true_50_totE_0_70=0;
    m_n_true_50_totMu_0_70=0;
    m_n_true_50_totN_0_70=0;
    m_n_true_50_totK_0_70=0;
    m_n_true_50_totP_0_70=0;
    m_n_true_50_totOtherCH_0_70=0;
    m_n_true_50_totOtherNeut_0_70=0;

    m_Z_mcE=-10;
    m_Z_mcNDaughter=-10;

    m_d1_mcPDGID=-10;
    m_d1_mcE=-10;
    m_d1_mcPx=-10;
    m_d1_mcPy=-10;
    m_d1_mcPz=-10;
    m_d1_mcMass=-10;
    m_d1_mcPhi=-10;
    m_d1_mcTheta=-10; 
    m_d1_mcCosTheta=-10;

    m_d2_mcPDGID=-10;
    m_d2_mcE=-10;
    m_d2_mcPx=-10;
    m_d2_mcPy=-10;
    m_d2_mcPz=-10;
    m_d2_mcMass=-10;
    m_d2_mcPhi=-10;
    m_d2_mcTheta=-10; 
    m_d2_mcCosTheta=-10;

    if(!m_ignoreGen){
      m_truePi0Energy->clear();
      m_truePi0_Px->clear();
      m_truePi0_Py->clear();
      m_truePi0_Pz->clear();
      m_truePi0_CosTheta->clear();
      m_truePi0_Phi->clear();
      m_truePi0_DRPh01->clear();
      m_truePi0_CosAngPh01->clear();
      m_truePi0_AngPh01->clear();
      m_truePh0Energy->clear();
      m_truePh0_Px->clear();
      m_truePh0_Py->clear();
      m_truePh0_Pz->clear();
      m_truePh0_CosTheta->clear();
      m_truePh0_Phi->clear();
      m_truePh0_DRMin->clear();
      m_truePh0_DRMinPDGID->clear();
      m_truePh0_DRMinE->clear();
      m_truePh0_CosAngMax->clear();
      m_truePh0_CosAngMaxPDGID->clear();
      m_truePh0_CosAngMaxE->clear();
      m_truePh0_CosAng0995E->clear();
      m_truePh0_DR01E->clear();
      m_truePh1Energy->clear();
      m_truePh1_Px->clear();
      m_truePh1_Py->clear();
      m_truePh1_Pz->clear();
      m_truePh1_CosTheta->clear();
      m_truePh1_Phi->clear();
      m_truePh1_DRMin->clear();
      m_truePh1_DRMinPDGID->clear();
      m_truePh1_DRMinE->clear();
      m_truePh1_CosAngMax->clear();
      m_truePh1_CosAngMaxPDGID->clear();
      m_truePh1_CosAngMaxE->clear();
      m_truePh1_CosAng0995E->clear();
      m_truePh1_DR01E->clear();
    }
    m_recoPhEnergy->clear();
    m_recoPh_Px->clear();
    m_recoPh_Py->clear();
    m_recoPh_Pz->clear();
    m_recoPh_CosTheta->clear();
    m_recoPh_Theta->clear();
    m_recoPh_Phi->clear();
    m_recoPh_PDGID->clear();
    m_recoPh_E_EB->clear();
    m_recoPh_E_EE->clear();
    m_recoPh_E_EO->clear();
    m_recoPh_E_HB->clear();
    m_recoPh_E_HE->clear();
    m_recoPh_E_HO->clear();
    m_recoPh_firstLayerECAL->clear();
    m_recoPh_lastLayerECAL->clear();
    m_recoPh_nhitsEB->clear();
    m_recoPh_nhitsEE->clear();
    m_recoPh_nhitsEO->clear();
    m_recoPh_firstLayerHCAL->clear();
    m_recoPh_lastLayerHCAL->clear();
    m_recoPh_nhitsHB->clear();
    m_recoPh_nhitsHE->clear();
    m_recoPh_nhitsHO->clear();
    m_recoPh_DRMin_Ph->clear();
    m_recoPh_DRMin_E_Ph->clear();
    m_recoPh_CosAngMax_Ph->clear();
    m_recoPh_CosAngMax_E_Ph->clear();
    m_recoPh_DRMin->clear();
    m_recoPh_DRMin_PDGID->clear();
    m_recoPh_DRMin_E->clear();
    m_recoPh_CosAngMax->clear();
    m_recoPh_CosAngMax_PDGID->clear();
    m_recoPh_CosAngMax_E->clear();
    m_recoPh_CosAng0995_E->clear();
    m_recoPh_DR01_E->clear();
    m_recoPi0Cand_DRMin_M->clear();
    m_recoPi0Cand_CosAngMax_M->clear();

    m_jet_excl_Px->clear(); 
    m_jet_excl_Py->clear(); 
    m_jet_excl_Pz->clear(); 
    m_jet_excl_E->clear(); 
    m_jet_excl_Phi->clear(); 
    m_jet_excl_CosTheta->clear(); 
    m_jet_excl_piE->clear();
    m_jet_excl_phE->clear();
    m_jet_excl_elE->clear();
    m_jet_excl_muE->clear();
    m_jet_excl_nE->clear();
    m_jet_excl_neutElseE->clear();
    m_jet_excl_chElseE->clear();
    m_jet_excl_neutMult->clear();
    m_jet_excl_chMult->clear();

    m_jet_incl_Px->clear(); 
    m_jet_incl_Py->clear(); 
    m_jet_incl_Pz->clear(); 
    m_jet_incl_E->clear(); 
    m_jet_incl_Phi->clear(); 
    m_jet_incl_CosTheta->clear(); 
    m_jet_incl_piE->clear();
    m_jet_incl_phE->clear();
    m_jet_incl_elE->clear();
    m_jet_incl_muE->clear();
    m_jet_incl_nE->clear();
    m_jet_incl_neutElseE->clear();
    m_jet_incl_chElseE->clear();
    m_jet_incl_neutMult->clear();
    m_jet_incl_chMult->clear();

    if(m_fillMEInfo){
      m_trueME_E->clear();
      m_trueME_Px->clear();
      m_trueME_Py->clear();
      m_trueME_Pz->clear();
      m_trueME_PDGID->clear();
    }
    m_outputTree->Branch("Z_mcE",&m_Z_mcE,"Z_mcE/F");
    m_outputTree->Branch("Z_mcNDaughter",&m_Z_mcNDaughter,"Z_mcNDaugther/I");

    m_outputTree->Branch("d1_mcPDGID",&m_d1_mcPDGID,"d1_mcPDGID/I");
    m_outputTree->Branch("d1_mcE",&m_d1_mcE,"d1_mcE/F");
    m_outputTree->Branch("d1_mcPx",&m_d1_mcPx,"d1_mcPx/F");
    m_outputTree->Branch("d1_mcPy",&m_d1_mcPy,"d1_mcPy/F");
    m_outputTree->Branch("d1_mcPx",&m_d1_mcPz,"d1_mcPz/F");
    m_outputTree->Branch("d1_mcMass",&m_d1_mcMass,"d1_mcMass/F");
    m_outputTree->Branch("d1_mcPhi",&m_d1_mcPhi,"d1_mcPhi/F");
    m_outputTree->Branch("d1_mcTheta",&m_d1_mcTheta,"d1_mcTheta/F");
    m_outputTree->Branch("d1_mcCosTheta",&m_d1_mcCosTheta,"d1_mcCosTheta/F");

    m_outputTree->Branch("d2_mcPDGID",&m_d2_mcPDGID,"d2_mcPDGID/I");
    m_outputTree->Branch("d2_mcE",&m_d2_mcE,"d2_mcE/F");
    m_outputTree->Branch("d2_mcPx",&m_d2_mcPx,"d2_mcPx/F");
    m_outputTree->Branch("d2_mcPy",&m_d2_mcPy,"d2_mcPy/F");
    m_outputTree->Branch("d2_mcPx",&m_d2_mcPz,"d2_mcPz/F");
    m_outputTree->Branch("d2_mcMass",&m_d2_mcMass,"d2_mcMass/F");
    m_outputTree->Branch("d2_mcPhi",&m_d2_mcPhi,"d2_mcPhi/F");
    m_outputTree->Branch("d2_mcTheta",&m_d2_mcTheta,"d2_mcTheta/F");
    m_outputTree->Branch("d2_mcCosTheta",&m_d2_mcCosTheta,"d2_mcCosTheta/F");

    if(m_fillMEInfo){
      m_outputTree->Branch("trueME_Px", "std::vector< float >", &m_trueME_Px); 
      m_outputTree->Branch("trueME_Py", "std::vector< float >", &m_trueME_Py); 
      m_outputTree->Branch("trueME_Pz", "std::vector< float >", &m_trueME_Pz); 
      m_outputTree->Branch("trueME_E", "std::vector< float >", &m_trueME_E); 
      m_outputTree->Branch("trueME_PDGID", "std::vector< int >", &m_trueME_PDGID); 
    }


    m_outputTree->Branch("px_totPFO", &m_px_totPFO, "px_totPFO/F");
    m_outputTree->Branch("px_totPi", &m_px_totPi, "px_totPi/F");
    m_outputTree->Branch("px_totPh", &m_px_totPh, "px_totPh/F");
    m_outputTree->Branch("px_totE", &m_px_totE, "px_totE/F");
    m_outputTree->Branch("px_totMu", &m_px_totMu, "px_totMu/F");
    m_outputTree->Branch("px_totN", &m_px_totN, "px_totN/F");

    m_outputTree->Branch("py_totPFO", &m_py_totPFO, "py_totPFO/F");
    m_outputTree->Branch("py_totPi", &m_py_totPi, "py_totPi/F");
    m_outputTree->Branch("py_totPh", &m_py_totPh, "py_totPh/F");
    m_outputTree->Branch("py_totE", &m_py_totE, "py_totE/F");
    m_outputTree->Branch("py_totMu", &m_py_totMu, "py_totMu/F");
    m_outputTree->Branch("py_totN", &m_py_totN, "py_totN/F");

    m_outputTree->Branch("pz_totPFO", &m_pz_totPFO, "pz_totPFO/F");
    m_outputTree->Branch("pz_totPi", &m_pz_totPi, "pz_totPi/F");
    m_outputTree->Branch("pz_totPh", &m_pz_totPh, "pz_totPh/F");
    m_outputTree->Branch("pz_totE", &m_pz_totE, "pz_totE/F");
    m_outputTree->Branch("pz_totMu", &m_pz_totMu, "pz_totMu/F");
    m_outputTree->Branch("pz_totN", &m_pz_totN, "pz_totN/F");

    m_outputTree->Branch("E_totPFO", &m_E_totPFO, "E_totPFO/F");
    m_outputTree->Branch("E_totPi", &m_E_totPi, "E_totPi/F");
    m_outputTree->Branch("E_totPh", &m_E_totPh, "E_totPh/F");
    m_outputTree->Branch("E_totE", &m_E_totE, "E_totE/F");
    m_outputTree->Branch("E_totMu", &m_E_totMu, "E_totMu/F");
    m_outputTree->Branch("E_totN", &m_E_totN, "E_totN/F");
    m_outputTree->Branch("E_totPFO_0_95", &m_E_totPFO_0_95, "E_totPFO_0_95/F");
    m_outputTree->Branch("E_totPi_0_95", &m_E_totPi_0_95, "E_totPi_0_95/F");
    m_outputTree->Branch("E_totPh_0_95", &m_E_totPh_0_95, "E_totPh_0_95/F");
    m_outputTree->Branch("E_totE_0_95", &m_E_totE_0_95, "E_totE_0_95/F");
    m_outputTree->Branch("E_totMu_0_95", &m_E_totMu_0_95, "E_totMu_0_95/F");
    m_outputTree->Branch("E_totN_0_95", &m_E_totN_0_95, "E_totN_0_95/F");
    if( m_reducedOutput){
      m_outputTree->Branch("E_totPFO_0_70", &m_E_totPFO_0_70, "E_totPFO_0_70/F");
      m_outputTree->Branch("E_totPi_0_70", &m_E_totPi_0_70, "E_totPi_0_70/F");
      m_outputTree->Branch("E_totPh_0_70", &m_E_totPh_0_70, "E_totPh_0_70/F");
      m_outputTree->Branch("E_totE_0_70", &m_E_totE_0_70, "E_totE_0_70/F");
      m_outputTree->Branch("E_totMu_0_70", &m_E_totMu_0_70, "E_totMu_0_70/F");
      m_outputTree->Branch("E_totN_0_70", &m_E_totN_0_70, "E_totN_0_70/F");
    }
    
    if(!m_ignoreGen){
      m_outputTree->Branch("px_true_totAll", &m_px_true_totAll, "px_true_totAll/F");
      m_outputTree->Branch("px_true_totInv", &m_px_true_totInv, "px_true_totInv/F");
      m_outputTree->Branch("px_true_totPi", &m_px_true_totPi, "px_true_totPi/F");
      m_outputTree->Branch("px_true_totPh", &m_px_true_totPh, "px_true_totPh/F");
      m_outputTree->Branch("px_true_totK0L", &m_px_true_totK0L, "px_true_totK0L/F");
      m_outputTree->Branch("px_true_totE", &m_px_true_totE, "px_true_totE/F");
      m_outputTree->Branch("px_true_totMu", &m_px_true_totMu, "px_true_totMu/F");
      m_outputTree->Branch("px_true_totN", &m_px_true_totN, "px_true_totN/F");
      m_outputTree->Branch("px_true_totP", &m_px_true_totP, "px_true_totP/F");
      m_outputTree->Branch("px_true_totK", &m_px_true_totK, "px_true_totK/F");

      m_outputTree->Branch("py_true_totAll", &m_py_true_totAll, "py_true_totAll/F");
      m_outputTree->Branch("py_true_totInv", &m_py_true_totInv, "py_true_totInv/F");
      m_outputTree->Branch("py_true_totPi", &m_py_true_totPi, "py_true_totPi/F");
      m_outputTree->Branch("py_true_totPh", &m_py_true_totPh, "py_true_totPh/F");
      m_outputTree->Branch("py_true_totK0L", &m_py_true_totK0L, "py_true_totK0L/F");
      m_outputTree->Branch("py_true_totE", &m_py_true_totE, "py_true_totE/F");
      m_outputTree->Branch("py_true_totMu", &m_py_true_totMu, "py_true_totMu/F");
      m_outputTree->Branch("py_true_totN", &m_py_true_totN, "py_true_totN/F");
      m_outputTree->Branch("py_true_totP", &m_py_true_totP, "py_true_totP/F");
      m_outputTree->Branch("py_true_totK", &m_py_true_totK, "py_true_totK/F");

      m_outputTree->Branch("pz_true_totAll", &m_pz_true_totAll, "pz_true_totAll/F");
      m_outputTree->Branch("pz_true_totInv", &m_pz_true_totInv, "pz_true_totInv/F");
      m_outputTree->Branch("pz_true_totPi", &m_pz_true_totPi, "pz_true_totPi/F");
      m_outputTree->Branch("pz_true_totPh", &m_pz_true_totPh, "pz_true_totPh/F");
      m_outputTree->Branch("pz_true_totK0L", &m_pz_true_totK0L, "pz_true_totK0L/F");
      m_outputTree->Branch("pz_true_totE", &m_pz_true_totE, "pz_true_totE/F");
      m_outputTree->Branch("pz_true_totMu", &m_pz_true_totMu, "pz_true_totMu/F");
      m_outputTree->Branch("pz_true_totN", &m_pz_true_totN, "pz_true_totN/F");
      m_outputTree->Branch("pz_true_totP", &m_pz_true_totP, "pz_true_totP/F");
      m_outputTree->Branch("pz_true_totK", &m_pz_true_totK, "pz_true_totK/F");

      m_outputTree->Branch("E_true_totAll", &m_E_true_totAll, "E_true_totAll/F");
      m_outputTree->Branch("E_true_totInv", &m_E_true_totInv, "E_true_totInv/F");
      m_outputTree->Branch("E_true_totPi", &m_E_true_totPi, "E_true_totPi/F");
      m_outputTree->Branch("E_true_totPh", &m_E_true_totPh, "E_true_totPh/F");
      m_outputTree->Branch("E_true_totK0L", &m_E_true_totK0L, "E_true_totK0L/F");
      m_outputTree->Branch("E_true_totE", &m_E_true_totE, "E_true_totE/F");
      m_outputTree->Branch("E_true_totMu", &m_E_true_totMu, "E_true_totMu/F");
      m_outputTree->Branch("E_true_totN", &m_E_true_totN, "E_true_totN/F");
      m_outputTree->Branch("E_true_totP", &m_E_true_totP, "E_true_totP/F");
      m_outputTree->Branch("E_true_totK", &m_E_true_totK, "E_true_totK/F");
      m_outputTree->Branch("E_true_totOtherCH", &m_E_true_totOtherCH, "E_true_totOtherCH/F");
      m_outputTree->Branch("E_true_totOtherNeut", &m_E_true_totOtherNeut, "E_true_totOtherNeut/F");
      m_outputTree->Branch("E_true_totAll_0_95", &m_E_true_totAll_0_95, "E_true_totAll_0_95/F");
      m_outputTree->Branch("E_true_totInv_0_95", &m_E_true_totInv_0_95, "E_true_totInv_0_95/F");
      m_outputTree->Branch("E_true_totPi_0_95", &m_E_true_totPi_0_95, "E_true_totPi_0_95/F");
      m_outputTree->Branch("E_true_totPh_0_95", &m_E_true_totPh_0_95, "E_true_totPh_0_95/F");
      m_outputTree->Branch("E_true_totK0L_0_95", &m_E_true_totK0L_0_95, "E_true_totK0L_0_95/F");
      m_outputTree->Branch("E_true_totE_0_95", &m_E_true_totE_0_95, "E_true_totE_0_95/F");
      m_outputTree->Branch("E_true_totMu_0_95", &m_E_true_totMu_0_95, "E_true_totMu_0_95/F");
      m_outputTree->Branch("E_true_totN_0_95", &m_E_true_totN_0_95, "E_true_totN_0_95/F");
      m_outputTree->Branch("E_true_totP_0_95", &m_E_true_totP_0_95, "E_true_totP_0_95/F");
      m_outputTree->Branch("E_true_totK_0_95", &m_E_true_totK_0_95, "E_true_totK_0_95/F");
      m_outputTree->Branch("E_true_totOtherCH_0_95", &m_E_true_totOtherCH_0_95, "E_true_totOtherCH_0_95/F");
      m_outputTree->Branch("E_true_totOtherNeut_0_95", &m_E_true_totOtherNeut_0_95, "E_true_totOtherNeut_0_95/F");
      if(!m_reducedOutput){
	m_outputTree->Branch("E_true_totAll_0_70", &m_E_true_totAll_0_70, "E_true_totAll_0_70/F");
	m_outputTree->Branch("E_true_totInv_0_70", &m_E_true_totInv_0_70, "E_true_totInv_0_70/F");
	m_outputTree->Branch("E_true_totPi_0_70", &m_E_true_totPi_0_70, "E_true_totPi_0_70/F");
	m_outputTree->Branch("E_true_totPh_0_70", &m_E_true_totPh_0_70, "E_true_totPh_0_70/F");
	m_outputTree->Branch("E_true_totK0L_0_70", &m_E_true_totK0L_0_70, "E_true_totK0L_0_70/F");
	m_outputTree->Branch("E_true_totE_0_70", &m_E_true_totE_0_70, "E_true_totE_0_70/F");
	m_outputTree->Branch("E_true_totMu_0_70", &m_E_true_totMu_0_70, "E_true_totMu_0_70/F");
	m_outputTree->Branch("E_true_totN_0_70", &m_E_true_totN_0_70, "E_true_totN_0_70/F");
	m_outputTree->Branch("E_true_totP_0_70", &m_E_true_totP_0_70, "E_true_totP_0_70/F");
	m_outputTree->Branch("E_true_totK_0_70", &m_E_true_totK_0_70, "E_true_totK_0_70/F");
	m_outputTree->Branch("E_true_totOtherCH_0_70", &m_E_true_totOtherCH_0_70, "E_true_totOtherCH_0_70/F");
	m_outputTree->Branch("E_true_totOtherNeut_0_70", &m_E_true_totOtherNeut_0_70, "E_true_totOtherNeut_0_70/F");
      }
    }
    if(!m_reducedOutput){
      m_outputTree->Branch("E_0_2_totPFO", &m_E_0_2_totPFO, "E_0_2_totPFO/F");
      m_outputTree->Branch("E_0_2_totPi", &m_E_0_2_totPi, "E_0_2_totPi/F");
      m_outputTree->Branch("E_0_2_totPh", &m_E_0_2_totPh, "E_0_2_totPh/F");
      m_outputTree->Branch("E_0_2_totE", &m_E_0_2_totE, "E_0_2_totE/F");
      m_outputTree->Branch("E_0_2_totMu", &m_E_0_2_totMu, "E_0_2_totMu/F");
      m_outputTree->Branch("E_0_2_totN", &m_E_0_2_totN, "E_0_2_totN/F");
      m_outputTree->Branch("E_0_2_totPFO_0_95", &m_E_0_2_totPFO_0_95, "E_0_2_totPFO_0_95/F");
      m_outputTree->Branch("E_0_2_totPi_0_95", &m_E_0_2_totPi_0_95, "E_0_2_totPi_0_95/F");
      m_outputTree->Branch("E_0_2_totPh_0_95", &m_E_0_2_totPh_0_95, "E_0_2_totPh_0_95/F");
      m_outputTree->Branch("E_0_2_totE_0_95", &m_E_0_2_totE_0_95, "E_0_2_totE_0_95/F");
      m_outputTree->Branch("E_0_2_totMu_0_95", &m_E_0_2_totMu_0_95, "E_0_2_totMu_0_95/F");
      m_outputTree->Branch("E_0_2_totN_0_95", &m_E_0_2_totN_0_95, "E_0_2_totN_0_95/F");
      m_outputTree->Branch("E_0_2_totPFO_0_70", &m_E_0_2_totPFO_0_70, "E_0_2_totPFO_0_70/F");
      m_outputTree->Branch("E_0_2_totPi_0_70", &m_E_0_2_totPi_0_70, "E_0_2_totPi_0_70/F");
      m_outputTree->Branch("E_0_2_totPh_0_70", &m_E_0_2_totPh_0_70, "E_0_2_totPh_0_70/F");
      m_outputTree->Branch("E_0_2_totE_0_70", &m_E_0_2_totE_0_70, "E_0_2_totE_0_70/F");
      m_outputTree->Branch("E_0_2_totMu_0_70", &m_E_0_2_totMu_0_70, "E_0_2_totMu_0_70/F");
      m_outputTree->Branch("E_0_2_totN_0_70", &m_E_0_2_totN_0_70, "E_0_2_totN_0_70/F");
      if(!m_ignoreGen){
	m_outputTree->Branch("E_true_0_2_totAll", &m_E_true_0_2_totAll, "E_true_0_2_totAll/F");
	m_outputTree->Branch("E_true_0_2_totInv", &m_E_true_0_2_totInv, "E_true_0_2_totInv/F");
	m_outputTree->Branch("E_true_0_2_totPi", &m_E_true_0_2_totPi, "E_true_0_2_totPi/F");
	m_outputTree->Branch("E_true_0_2_totPh", &m_E_true_0_2_totPh, "E_true_0_2_totPh/F");
	m_outputTree->Branch("E_true_0_2_totK0L", &m_E_true_0_2_totK0L, "E_true_0_2_totK0L/F");
	m_outputTree->Branch("E_true_0_2_totE", &m_E_true_0_2_totE, "E_true_0_2_totE/F");
	m_outputTree->Branch("E_true_0_2_totMu", &m_E_true_0_2_totMu, "E_true_0_2_totMu/F");
	m_outputTree->Branch("E_true_0_2_totN", &m_E_true_0_2_totN, "E_true_0_2_totN/F");
	m_outputTree->Branch("E_true_0_2_totP", &m_E_true_0_2_totP, "E_true_0_2_totP/F");
	m_outputTree->Branch("E_true_0_2_totK", &m_E_true_0_2_totK, "E_true_0_2_totK/F");
	m_outputTree->Branch("E_true_0_2_totOtherCH", &m_E_true_0_2_totOtherCH, "E_true_0_2_totOtherCH/F");
	m_outputTree->Branch("E_true_0_2_totOtherNeut", &m_E_true_0_2_totOtherNeut, "E_true_0_2_totOtherNeut/F");
	m_outputTree->Branch("E_true_0_2_totAll_0_95", &m_E_true_0_2_totAll_0_95, "E_true_0_2_totAll_0_95/F");
	m_outputTree->Branch("E_true_0_2_totInv_0_95", &m_E_true_0_2_totInv_0_95, "E_true_0_2_totInv_0_95/F");
	m_outputTree->Branch("E_true_0_2_totPi_0_95", &m_E_true_0_2_totPi_0_95, "E_true_0_2_totPi_0_95/F");
	m_outputTree->Branch("E_true_0_2_totPh_0_95", &m_E_true_0_2_totPh_0_95, "E_true_0_2_totPh_0_95/F");
	m_outputTree->Branch("E_true_0_2_totK0L_0_95", &m_E_true_0_2_totK0L_0_95, "E_true_0_2_totK0L_0_95/F");
	m_outputTree->Branch("E_true_0_2_totE_0_95", &m_E_true_0_2_totE_0_95, "E_true_0_2_totE_0_95/F");
	m_outputTree->Branch("E_true_0_2_totMu_0_95", &m_E_true_0_2_totMu_0_95, "E_true_0_2_totMu_0_95/F");
	m_outputTree->Branch("E_true_0_2_totN_0_95", &m_E_true_0_2_totN_0_95, "E_true_0_2_totN_0_95/F");
	m_outputTree->Branch("E_true_0_2_totP_0_95", &m_E_true_0_2_totP_0_95, "E_true_0_2_totP_0_95/F");
	m_outputTree->Branch("E_true_0_2_totK_0_95", &m_E_true_0_2_totK_0_95, "E_true_0_2_totK_0_95/F");
	m_outputTree->Branch("E_true_0_2_totOtherCH_0_95", &m_E_true_0_2_totOtherCH_0_95, "E_true_0_2_totOtherCH_0_95/F");
	m_outputTree->Branch("E_true_0_2_totOtherNeut_0_95", &m_E_true_0_2_totOtherNeut_0_95, "E_true_0_2_totOtherNeut_0_95/F");
	m_outputTree->Branch("E_true_0_2_totAll_0_70", &m_E_true_0_2_totAll_0_70, "E_true_0_2_totAll_0_70/F");
	m_outputTree->Branch("E_true_0_2_totInv_0_70", &m_E_true_0_2_totInv_0_70, "E_true_0_2_totInv_0_70/F");
	m_outputTree->Branch("E_true_0_2_totPi_0_70", &m_E_true_0_2_totPi_0_70, "E_true_0_2_totPi_0_70/F");
	m_outputTree->Branch("E_true_0_2_totPh_0_70", &m_E_true_0_2_totPh_0_70, "E_true_0_2_totPh_0_70/F");
	m_outputTree->Branch("E_true_0_2_totK0L_0_70", &m_E_true_0_2_totK0L_0_70, "E_true_0_2_totK0L_0_70/F");
	m_outputTree->Branch("E_true_0_2_totE_0_70", &m_E_true_0_2_totE_0_70, "E_true_0_2_totE_0_70/F");
	m_outputTree->Branch("E_true_0_2_totMu_0_70", &m_E_true_0_2_totMu_0_70, "E_true_0_2_totMu_0_70/F");
	m_outputTree->Branch("E_true_0_2_totN_0_70", &m_E_true_0_2_totN_0_70, "E_true_0_2_totN_0_70/F");
	m_outputTree->Branch("E_true_0_2_totP_0_70", &m_E_true_0_2_totP_0_70, "E_true_0_2_totP_0_70/F");
	m_outputTree->Branch("E_true_0_2_totK_0_70", &m_E_true_0_2_totK_0_70, "E_true_0_2_totK_0_70/F");
	m_outputTree->Branch("E_true_0_2_totOtherCH_0_70", &m_E_true_0_2_totOtherCH_0_70, "E_true_0_2_totOtherCH_0_70/F");
	m_outputTree->Branch("E_true_0_2_totOtherNeut_0_70", &m_E_true_0_2_totOtherNeut_0_70, "E_true_0_2_totOtherNeut_0_70/F");
      }
      m_outputTree->Branch("E_2_10_totPFO", &m_E_2_10_totPFO, "E_2_10_totPFO/F");
      m_outputTree->Branch("E_2_10_totPi", &m_E_2_10_totPi, "E_2_10_totPi/F");
      m_outputTree->Branch("E_2_10_totPh", &m_E_2_10_totPh, "E_2_10_totPh/F");
      m_outputTree->Branch("E_2_10_totE", &m_E_2_10_totE, "E_2_10_totE/F");
      m_outputTree->Branch("E_2_10_totMu", &m_E_2_10_totMu, "E_2_10_totMu/F");
      m_outputTree->Branch("E_2_10_totN", &m_E_2_10_totN, "E_2_10_totN/F");
      m_outputTree->Branch("E_2_10_totPFO_0_95", &m_E_2_10_totPFO_0_95, "E_2_10_totPFO_0_95/F");
      m_outputTree->Branch("E_2_10_totPi_0_95", &m_E_2_10_totPi_0_95, "E_2_10_totPi_0_95/F");
      m_outputTree->Branch("E_2_10_totPh_0_95", &m_E_2_10_totPh_0_95, "E_2_10_totPh_0_95/F");
      m_outputTree->Branch("E_2_10_totE_0_95", &m_E_2_10_totE_0_95, "E_2_10_totE_0_95/F");
      m_outputTree->Branch("E_2_10_totMu_0_95", &m_E_2_10_totMu_0_95, "E_2_10_totMu_0_95/F");
      m_outputTree->Branch("E_2_10_totN_0_95", &m_E_2_10_totN_0_95, "E_2_10_totN_0_95/F");
      m_outputTree->Branch("E_2_10_totPFO_0_70", &m_E_2_10_totPFO_0_70, "E_2_10_totPFO_0_70/F");
      m_outputTree->Branch("E_2_10_totPi_0_70", &m_E_2_10_totPi_0_70, "E_2_10_totPi_0_70/F");
      m_outputTree->Branch("E_2_10_totPh_0_70", &m_E_2_10_totPh_0_70, "E_2_10_totPh_0_70/F");
      m_outputTree->Branch("E_2_10_totE_0_70", &m_E_2_10_totE_0_70, "E_2_10_totE_0_70/F");
      m_outputTree->Branch("E_2_10_totMu_0_70", &m_E_2_10_totMu_0_70, "E_2_10_totMu_0_70/F");
      m_outputTree->Branch("E_2_10_totN_0_70", &m_E_2_10_totN_0_70, "E_2_10_totN_0_70/F");
      
      if(!m_ignoreGen){
	m_outputTree->Branch("E_true_2_10_totAll", &m_E_true_2_10_totAll, "E_true_2_10_totAll/F");
	m_outputTree->Branch("E_true_2_10_totInv", &m_E_true_2_10_totInv, "E_true_2_10_totInv/F");
	m_outputTree->Branch("E_true_2_10_totPi", &m_E_true_2_10_totPi, "E_true_2_10_totPi/F");
	m_outputTree->Branch("E_true_2_10_totPh", &m_E_true_2_10_totPh, "E_true_2_10_totPh/F");
	m_outputTree->Branch("E_true_2_10_totK0L", &m_E_true_2_10_totK0L, "E_true_2_10_totK0L/F");
	m_outputTree->Branch("E_true_2_10_totE", &m_E_true_2_10_totE, "E_true_2_10_totE/F");
	m_outputTree->Branch("E_true_2_10_totMu", &m_E_true_2_10_totMu, "E_true_2_10_totMu/F");
	m_outputTree->Branch("E_true_2_10_totN", &m_E_true_2_10_totN, "E_true_2_10_totN/F");
	m_outputTree->Branch("E_true_2_10_totP", &m_E_true_2_10_totP, "E_true_2_10_totP/F");
	m_outputTree->Branch("E_true_2_10_totK", &m_E_true_2_10_totK, "E_true_2_10_totK/F");
	m_outputTree->Branch("E_true_2_10_totOtherCH", &m_E_true_2_10_totOtherCH, "E_true_2_10_totOtherCH/F");
	m_outputTree->Branch("E_true_2_10_totOtherNeut", &m_E_true_2_10_totOtherNeut, "E_true_2_10_totOtherNeut/F");
	m_outputTree->Branch("E_true_2_10_totAll_0_95", &m_E_true_2_10_totAll_0_95, "E_true_2_10_totAll_0_95/F");
	m_outputTree->Branch("E_true_2_10_totInv_0_95", &m_E_true_2_10_totInv_0_95, "E_true_2_10_totInv_0_95/F");
	m_outputTree->Branch("E_true_2_10_totPi_0_95", &m_E_true_2_10_totPi_0_95, "E_true_2_10_totPi_0_95/F");
	m_outputTree->Branch("E_true_2_10_totPh_0_95", &m_E_true_2_10_totPh_0_95, "E_true_2_10_totPh_0_95/F");
	m_outputTree->Branch("E_true_2_10_totK0L_0_95", &m_E_true_2_10_totK0L_0_95, "E_true_2_10_totK0L_0_95/F");
	m_outputTree->Branch("E_true_2_10_totE_0_95", &m_E_true_2_10_totE_0_95, "E_true_2_10_totE_0_95/F");
	m_outputTree->Branch("E_true_2_10_totMu_0_95", &m_E_true_2_10_totMu_0_95, "E_true_2_10_totMu_0_95/F");
	m_outputTree->Branch("E_true_2_10_totN_0_95", &m_E_true_2_10_totN_0_95, "E_true_2_10_totN_0_95/F");
	m_outputTree->Branch("E_true_2_10_totP_0_95", &m_E_true_2_10_totP_0_95, "E_true_2_10_totP_0_95/F");
	m_outputTree->Branch("E_true_2_10_totK_0_95", &m_E_true_2_10_totK_0_95, "E_true_2_10_totK_0_95/F");
	m_outputTree->Branch("E_true_2_10_totOtherCH_0_95", &m_E_true_2_10_totOtherCH_0_95, "E_true_2_10_totOtherCH_0_95/F");
	m_outputTree->Branch("E_true_2_10_totOtherNeut_0_95", &m_E_true_2_10_totOtherNeut_0_95, "E_true_2_10_totOtherNeut_0_95/F");
	m_outputTree->Branch("E_true_2_10_totAll_0_70", &m_E_true_2_10_totAll_0_70, "E_true_2_10_totAll_0_70/F");
	m_outputTree->Branch("E_true_2_10_totInv_0_70", &m_E_true_2_10_totInv_0_70, "E_true_2_10_totInv_0_70/F");
	m_outputTree->Branch("E_true_2_10_totPi_0_70", &m_E_true_2_10_totPi_0_70, "E_true_2_10_totPi_0_70/F");
	m_outputTree->Branch("E_true_2_10_totPh_0_70", &m_E_true_2_10_totPh_0_70, "E_true_2_10_totPh_0_70/F");
	m_outputTree->Branch("E_true_2_10_totK0L_0_70", &m_E_true_2_10_totK0L_0_70, "E_true_2_10_totK0L_0_70/F");
	m_outputTree->Branch("E_true_2_10_totE_0_70", &m_E_true_2_10_totE_0_70, "E_true_2_10_totE_0_70/F");
	m_outputTree->Branch("E_true_2_10_totMu_0_70", &m_E_true_2_10_totMu_0_70, "E_true_2_10_totMu_0_70/F");
	m_outputTree->Branch("E_true_2_10_totN_0_70", &m_E_true_2_10_totN_0_70, "E_true_2_10_totN_0_70/F");
	m_outputTree->Branch("E_true_2_10_totP_0_70", &m_E_true_2_10_totP_0_70, "E_true_2_10_totP_0_70/F");
	m_outputTree->Branch("E_true_2_10_totK_0_70", &m_E_true_2_10_totK_0_70, "E_true_2_10_totK_0_70/F");
	m_outputTree->Branch("E_true_2_10_totOtherCH_0_70", &m_E_true_2_10_totOtherCH_0_70, "E_true_2_10_totOtherCH_0_70/F");
	m_outputTree->Branch("E_true_2_10_totOtherNeut_0_70", &m_E_true_2_10_totOtherNeut_0_70, "E_true_2_10_totOtherNeut_0_70/F");
      }
      m_outputTree->Branch("E_10_50_totPFO", &m_E_10_50_totPFO, "E_10_50_totPFO/F");
      m_outputTree->Branch("E_10_50_totPi", &m_E_10_50_totPi, "E_10_50_totPi/F");
      m_outputTree->Branch("E_10_50_totPh", &m_E_10_50_totPh, "E_10_50_totPh/F");
      m_outputTree->Branch("E_10_50_totE", &m_E_10_50_totE, "E_10_50_totE/F");
      m_outputTree->Branch("E_10_50_totMu", &m_E_10_50_totMu, "E_10_50_totMu/F");
      m_outputTree->Branch("E_10_50_totN", &m_E_10_50_totN, "E_10_50_totN/F");
      m_outputTree->Branch("E_10_50_totPFO_0_95", &m_E_10_50_totPFO_0_95, "E_10_50_totPFO_0_95/F");
      m_outputTree->Branch("E_10_50_totPi_0_95", &m_E_10_50_totPi_0_95, "E_10_50_totPi_0_95/F");
      m_outputTree->Branch("E_10_50_totPh_0_95", &m_E_10_50_totPh_0_95, "E_10_50_totPh_0_95/F");
      m_outputTree->Branch("E_10_50_totE_0_95", &m_E_10_50_totE_0_95, "E_10_50_totE_0_95/F");
      m_outputTree->Branch("E_10_50_totMu_0_95", &m_E_10_50_totMu_0_95, "E_10_50_totMu_0_95/F");
      m_outputTree->Branch("E_10_50_totN_0_95", &m_E_10_50_totN_0_95, "E_10_50_totN_0_95/F");
      m_outputTree->Branch("E_10_50_totPFO_0_70", &m_E_10_50_totPFO_0_70, "E_10_50_totPFO_0_70/F");
      m_outputTree->Branch("E_10_50_totPi_0_70", &m_E_10_50_totPi_0_70, "E_10_50_totPi_0_70/F");
      m_outputTree->Branch("E_10_50_totPh_0_70", &m_E_10_50_totPh_0_70, "E_10_50_totPh_0_70/F");
      m_outputTree->Branch("E_10_50_totE_0_70", &m_E_10_50_totE_0_70, "E_10_50_totE_0_70/F");
      m_outputTree->Branch("E_10_50_totMu_0_70", &m_E_10_50_totMu_0_70, "E_10_50_totMu_0_70/F");
      m_outputTree->Branch("E_10_50_totN_0_70", &m_E_10_50_totN_0_70, "E_10_50_totN_0_70/F");
      
      if(!m_ignoreGen){
	m_outputTree->Branch("E_true_10_50_totAll", &m_E_true_10_50_totAll, "E_true_10_50_totAll/F");
	m_outputTree->Branch("E_true_10_50_totInv", &m_E_true_10_50_totInv, "E_true_10_50_totInv/F");
	m_outputTree->Branch("E_true_10_50_totPi", &m_E_true_10_50_totPi, "E_true_10_50_totPi/F");
	m_outputTree->Branch("E_true_10_50_totPh", &m_E_true_10_50_totPh, "E_true_10_50_totPh/F");
	m_outputTree->Branch("E_true_10_50_totK0L", &m_E_true_10_50_totK0L, "E_true_10_50_totK0L/F");
	m_outputTree->Branch("E_true_10_50_totE", &m_E_true_10_50_totE, "E_true_10_50_totE/F");
	m_outputTree->Branch("E_true_10_50_totMu", &m_E_true_10_50_totMu, "E_true_10_50_totMu/F");
	m_outputTree->Branch("E_true_10_50_totN", &m_E_true_10_50_totN, "E_true_10_50_totN/F");
	m_outputTree->Branch("E_true_10_50_totP", &m_E_true_10_50_totP, "E_true_10_50_totP/F");
	m_outputTree->Branch("E_true_10_50_totK", &m_E_true_10_50_totK, "E_true_10_50_totK/F");
	m_outputTree->Branch("E_true_10_50_totOtherCH", &m_E_true_10_50_totOtherCH, "E_true_10_50_totOtherCH/F");
	m_outputTree->Branch("E_true_10_50_totOtherNeut", &m_E_true_10_50_totOtherNeut, "E_true_10_50_totOtherNeut/F");
	m_outputTree->Branch("E_true_10_50_totAll_0_95", &m_E_true_10_50_totAll_0_95, "E_true_10_50_totAll_0_95/F");
	m_outputTree->Branch("E_true_10_50_totInv_0_95", &m_E_true_10_50_totInv_0_95, "E_true_10_50_totInv_0_95/F");
	m_outputTree->Branch("E_true_10_50_totPi_0_95", &m_E_true_10_50_totPi_0_95, "E_true_10_50_totPi_0_95/F");
	m_outputTree->Branch("E_true_10_50_totPh_0_95", &m_E_true_10_50_totPh_0_95, "E_true_10_50_totPh_0_95/F");
	m_outputTree->Branch("E_true_10_50_totK0L_0_95", &m_E_true_10_50_totK0L_0_95, "E_true_10_50_totK0L_0_95/F");
	m_outputTree->Branch("E_true_10_50_totE_0_95", &m_E_true_10_50_totE_0_95, "E_true_10_50_totE_0_95/F");
	m_outputTree->Branch("E_true_10_50_totMu_0_95", &m_E_true_10_50_totMu_0_95, "E_true_10_50_totMu_0_95/F");
	m_outputTree->Branch("E_true_10_50_totN_0_95", &m_E_true_10_50_totN_0_95, "E_true_10_50_totN_0_95/F");
	m_outputTree->Branch("E_true_10_50_totP_0_95", &m_E_true_10_50_totP_0_95, "E_true_10_50_totP_0_95/F");
	m_outputTree->Branch("E_true_10_50_totK_0_95", &m_E_true_10_50_totK_0_95, "E_true_10_50_totK_0_95/F");
	m_outputTree->Branch("E_true_10_50_totOtherCH_0_95", &m_E_true_10_50_totOtherCH_0_95, "E_true_10_50_totOtherCH_0_95/F");
	m_outputTree->Branch("E_true_10_50_totOtherNeut_0_95", &m_E_true_10_50_totOtherNeut_0_95, "E_true_10_50_totOtherNeut_0_95/F");
	m_outputTree->Branch("E_true_10_50_totAll_0_70", &m_E_true_10_50_totAll_0_70, "E_true_10_50_totAll_0_70/F");
	m_outputTree->Branch("E_true_10_50_totInv_0_70", &m_E_true_10_50_totInv_0_70, "E_true_10_50_totInv_0_70/F");
	m_outputTree->Branch("E_true_10_50_totPi_0_70", &m_E_true_10_50_totPi_0_70, "E_true_10_50_totPi_0_70/F");
	m_outputTree->Branch("E_true_10_50_totPh_0_70", &m_E_true_10_50_totPh_0_70, "E_true_10_50_totPh_0_70/F");
	m_outputTree->Branch("E_true_10_50_totK0L_0_70", &m_E_true_10_50_totK0L_0_70, "E_true_10_50_totK0L_0_70/F");
	m_outputTree->Branch("E_true_10_50_totE_0_70", &m_E_true_10_50_totE_0_70, "E_true_10_50_totE_0_70/F");
	m_outputTree->Branch("E_true_10_50_totMu_0_70", &m_E_true_10_50_totMu_0_70, "E_true_10_50_totMu_0_70/F");
	m_outputTree->Branch("E_true_10_50_totN_0_70", &m_E_true_10_50_totN_0_70, "E_true_10_50_totN_0_70/F");
	m_outputTree->Branch("E_true_10_50_totP_0_70", &m_E_true_10_50_totP_0_70, "E_true_10_50_totP_0_70/F");
	m_outputTree->Branch("E_true_10_50_totK_0_70", &m_E_true_10_50_totK_0_70, "E_true_10_50_totK_0_70/F");
	m_outputTree->Branch("E_true_10_50_totOtherCH_0_70", &m_E_true_10_50_totOtherCH_0_70, "E_true_10_50_totOtherCH_0_70/F");
	m_outputTree->Branch("E_true_10_50_totOtherNeut_0_70", &m_E_true_10_50_totOtherNeut_0_70, "E_true_10_50_totOtherNeut_0_70/F");
      }
      m_outputTree->Branch("E_50_totPFO", &m_E_50_totPFO, "E_50_totPFO/F");
      m_outputTree->Branch("E_50_totPi", &m_E_50_totPi, "E_50_totPi/F");
      m_outputTree->Branch("E_50_totPh", &m_E_50_totPh, "E_50_totPh/F");
      m_outputTree->Branch("E_50_totE", &m_E_50_totE, "E_50_totE/F");
      m_outputTree->Branch("E_50_totMu", &m_E_50_totMu, "E_50_totMu/F");
      m_outputTree->Branch("E_50_totN", &m_E_50_totN, "E_50_totN/F");
      m_outputTree->Branch("E_50_totPFO_0_95", &m_E_50_totPFO_0_95, "E_50_totPFO_0_95/F");
      m_outputTree->Branch("E_50_totPi_0_95", &m_E_50_totPi_0_95, "E_50_totPi_0_95/F");
      m_outputTree->Branch("E_50_totPh_0_95", &m_E_50_totPh_0_95, "E_50_totPh_0_95/F");
      m_outputTree->Branch("E_50_totE_0_95", &m_E_50_totE_0_95, "E_50_totE_0_95/F");
      m_outputTree->Branch("E_50_totMu_0_95", &m_E_50_totMu_0_95, "E_50_totMu_0_95/F");
      m_outputTree->Branch("E_50_totN_0_95", &m_E_50_totN_0_95, "E_50_totN_0_95/F");
      m_outputTree->Branch("E_50_totPFO_0_70", &m_E_50_totPFO_0_70, "E_50_totPFO_0_70/F");
      m_outputTree->Branch("E_50_totPi_0_70", &m_E_50_totPi_0_70, "E_50_totPi_0_70/F");
      m_outputTree->Branch("E_50_totPh_0_70", &m_E_50_totPh_0_70, "E_50_totPh_0_70/F");
      m_outputTree->Branch("E_50_totE_0_70", &m_E_50_totE_0_70, "E_50_totE_0_70/F");
      m_outputTree->Branch("E_50_totMu_0_70", &m_E_50_totMu_0_70, "E_50_totMu_0_70/F");
      m_outputTree->Branch("E_50_totN_0_70", &m_E_50_totN_0_70, "E_50_totN_0_70/F");
      if(!m_ignoreGen){
	m_outputTree->Branch("E_true_50_totAll", &m_E_true_50_totAll, "E_true_50_totAll/F");
	m_outputTree->Branch("E_true_50_totInv", &m_E_true_50_totInv, "E_true_50_totInv/F");
	m_outputTree->Branch("E_true_50_totPi", &m_E_true_50_totPi, "E_true_50_totPi/F");
	m_outputTree->Branch("E_true_50_totPh", &m_E_true_50_totPh, "E_true_50_totPh/F");
	m_outputTree->Branch("E_true_50_totK0L", &m_E_true_50_totK0L, "E_true_50_totK0L/F");
	m_outputTree->Branch("E_true_50_totE", &m_E_true_50_totE, "E_true_50_totE/F");
	m_outputTree->Branch("E_true_50_totMu", &m_E_true_50_totMu, "E_true_50_totMu/F");
	m_outputTree->Branch("E_true_50_totN", &m_E_true_50_totN, "E_true_50_totN/F");
	m_outputTree->Branch("E_true_50_totP", &m_E_true_50_totP, "E_true_50_totP/F");
	m_outputTree->Branch("E_true_50_totK", &m_E_true_50_totK, "E_true_50_totK/F");
	m_outputTree->Branch("E_true_50_totOtherCH", &m_E_true_50_totOtherCH, "E_true_50_totOtherCH/F");
	m_outputTree->Branch("E_true_50_totOtherNeut", &m_E_true_50_totOtherNeut, "E_true_50_totOtherNeut/F");
	m_outputTree->Branch("E_true_50_totAll_0_95", &m_E_true_50_totAll_0_95, "E_true_50_totAll_0_95/F");
	m_outputTree->Branch("E_true_50_totInv_0_95", &m_E_true_50_totInv_0_95, "E_true_50_totInv_0_95/F");
	m_outputTree->Branch("E_true_50_totPi_0_95", &m_E_true_50_totPi_0_95, "E_true_50_totPi_0_95/F");
	m_outputTree->Branch("E_true_50_totPh_0_95", &m_E_true_50_totPh_0_95, "E_true_50_totPh_0_95/F");
	m_outputTree->Branch("E_true_50_totK0L_0_95", &m_E_true_50_totK0L_0_95, "E_true_50_totK0L_0_95/F");
	m_outputTree->Branch("E_true_50_totE_0_95", &m_E_true_50_totE_0_95, "E_true_50_totE_0_95/F");
	m_outputTree->Branch("E_true_50_totMu_0_95", &m_E_true_50_totMu_0_95, "E_true_50_totMu_0_95/F");
	m_outputTree->Branch("E_true_50_totN_0_95", &m_E_true_50_totN_0_95, "E_true_50_totN_0_95/F");
	m_outputTree->Branch("E_true_50_totP_0_95", &m_E_true_50_totP_0_95, "E_true_50_totP_0_95/F");
	m_outputTree->Branch("E_true_50_totK_0_95", &m_E_true_50_totK_0_95, "E_true_50_totK_0_95/F");
	m_outputTree->Branch("E_true_50_totOtherCH_0_95", &m_E_true_50_totOtherCH_0_95, "E_true_50_totOtherCH_0_95/F");
	m_outputTree->Branch("E_true_50_totOtherNeut_0_95", &m_E_true_50_totOtherNeut_0_95, "E_true_50_totOtherNeut_0_95/F");
	m_outputTree->Branch("E_true_50_totAll_0_70", &m_E_true_50_totAll_0_70, "E_true_50_totAll_0_70/F");
	m_outputTree->Branch("E_true_50_totInv_0_70", &m_E_true_50_totInv_0_70, "E_true_50_totInv_0_70/F");
	m_outputTree->Branch("E_true_50_totPi_0_70", &m_E_true_50_totPi_0_70, "E_true_50_totPi_0_70/F");
	m_outputTree->Branch("E_true_50_totPh_0_70", &m_E_true_50_totPh_0_70, "E_true_50_totPh_0_70/F");
	m_outputTree->Branch("E_true_50_totK0L_0_70", &m_E_true_50_totK0L_0_70, "E_true_50_totK0L_0_70/F");
	m_outputTree->Branch("E_true_50_totE_0_70", &m_E_true_50_totE_0_70, "E_true_50_totE_0_70/F");
	m_outputTree->Branch("E_true_50_totMu_0_70", &m_E_true_50_totMu_0_70, "E_true_50_totMu_0_70/F");
	m_outputTree->Branch("E_true_50_totN_0_70", &m_E_true_50_totN_0_70, "E_true_50_totN_0_70/F");
	m_outputTree->Branch("E_true_50_totP_0_70", &m_E_true_50_totP_0_70, "E_true_50_totP_0_70/F");
	m_outputTree->Branch("E_true_50_totK_0_70", &m_E_true_50_totK_0_70, "E_true_50_totK_0_70/F");
	m_outputTree->Branch("E_true_50_totOtherCH_0_70", &m_E_true_50_totOtherCH_0_70, "E_true_50_totOtherCH_0_70/F");
	m_outputTree->Branch("E_true_50_totOtherNeut_0_70", &m_E_true_50_totOtherNeut_0_70, "E_true_50_totOtherNeut_0_70/F");
      }
      //now multiplicities
      m_outputTree->Branch("n_totPFO", &m_n_totPFO, "n_totPFO/i");
      m_outputTree->Branch("n_totPi", &m_n_totPi, "n_totPi/i");
      m_outputTree->Branch("n_totPh", &m_n_totPh, "n_totPh/i");
      m_outputTree->Branch("n_totE", &m_n_totE, "n_totE/i");
      m_outputTree->Branch("n_totMu", &m_n_totMu, "n_totMu/i");
      m_outputTree->Branch("n_totN", &m_n_totN, "n_totN/i");
      m_outputTree->Branch("n_totPFO_0_95", &m_n_totPFO_0_95, "n_totPFO_0_95/i");
      m_outputTree->Branch("n_totPi_0_95", &m_n_totPi_0_95, "n_totPi_0_95/i");
      m_outputTree->Branch("n_totPh_0_95", &m_n_totPh_0_95, "n_totPh_0_95/i");
      m_outputTree->Branch("n_totE_0_95", &m_n_totE_0_95, "n_totE_0_95/i");
      m_outputTree->Branch("n_totMu_0_95", &m_n_totMu_0_95, "n_totMu_0_95/i");
      m_outputTree->Branch("n_totN_0_95", &m_n_totN_0_95, "n_totN_0_95/i");
      m_outputTree->Branch("n_totPFO_0_70", &m_n_totPFO_0_70, "n_totPFO_0_70/i");
      m_outputTree->Branch("n_totPi_0_70", &m_n_totPi_0_70, "n_totPi_0_70/i");
      m_outputTree->Branch("n_totPh_0_70", &m_n_totPh_0_70, "n_totPh_0_70/i");
      m_outputTree->Branch("n_totE_0_70", &m_n_totE_0_70, "n_totE_0_70/i");
      m_outputTree->Branch("n_totMu_0_70", &m_n_totMu_0_70, "n_totMu_0_70/i");
      m_outputTree->Branch("n_totN_0_70", &m_n_totN_0_70, "n_totN_0_70/i");
      
      if(!m_ignoreGen){
	m_outputTree->Branch("n_true_totAll", &m_n_true_totAll, "n_true_totAll/i");
	m_outputTree->Branch("n_true_totInv", &m_n_true_totInv, "n_true_totInv/i");
	m_outputTree->Branch("n_true_totPi", &m_n_true_totPi, "n_true_totPi/i");
	m_outputTree->Branch("n_true_totPh", &m_n_true_totPh, "n_true_totPh/i");
	m_outputTree->Branch("n_true_totK0L", &m_n_true_totK0L, "n_true_totK0L/i");
	m_outputTree->Branch("n_true_totE", &m_n_true_totE, "n_true_totE/i");
	m_outputTree->Branch("n_true_totMu", &m_n_true_totMu, "n_true_totMu/i");
	m_outputTree->Branch("n_true_totN", &m_n_true_totN, "n_true_totN/i");
	m_outputTree->Branch("n_true_totP", &m_n_true_totP, "n_true_totP/i");
	m_outputTree->Branch("n_true_totK", &m_n_true_totK, "n_true_totK/i");
	m_outputTree->Branch("n_true_totOtherCH", &m_n_true_totOtherCH, "n_true_totOtherCH/i");
	m_outputTree->Branch("n_true_totOtherNeut", &m_n_true_totOtherNeut, "n_true_totOtherNeut/i");
	m_outputTree->Branch("n_true_totAll_0_95", &m_n_true_totAll_0_95, "n_true_totAll_0_95/i");
	m_outputTree->Branch("n_true_totInv_0_95", &m_n_true_totInv_0_95, "n_true_totInv_0_95/i");
	m_outputTree->Branch("n_true_totPi_0_95", &m_n_true_totPi_0_95, "n_true_totPi_0_95/i");
	m_outputTree->Branch("n_true_totPh_0_95", &m_n_true_totPh_0_95, "n_true_totPh_0_95/i");
	m_outputTree->Branch("n_true_totK0L_0_95", &m_n_true_totK0L_0_95, "n_true_totK0L_0_95/i");
	m_outputTree->Branch("n_true_totE_0_95", &m_n_true_totE_0_95, "n_true_totE_0_95/i");
	m_outputTree->Branch("n_true_totMu_0_95", &m_n_true_totMu_0_95, "n_true_totMu_0_95/i");
	m_outputTree->Branch("n_true_totN_0_95", &m_n_true_totN_0_95, "n_true_totN_0_95/i");
	m_outputTree->Branch("n_true_totP_0_95", &m_n_true_totP_0_95, "n_true_totP_0_95/i");
	m_outputTree->Branch("n_true_totK_0_95", &m_n_true_totK_0_95, "n_true_totK_0_95/i");
	m_outputTree->Branch("n_true_totOtherCH_0_95", &m_n_true_totOtherCH_0_95, "n_true_totOtherCH_0_95/i");
	m_outputTree->Branch("n_true_totOtherNeut_0_95", &m_n_true_totOtherNeut_0_95, "n_true_totOtherNeut_0_95/i");
	m_outputTree->Branch("n_true_totAll_0_70", &m_n_true_totAll_0_70, "n_true_totAll_0_70/i");
	m_outputTree->Branch("n_true_totInv_0_70", &m_n_true_totInv_0_70, "n_true_totInv_0_70/i");
	m_outputTree->Branch("n_true_totPi_0_70", &m_n_true_totPi_0_70, "n_true_totPi_0_70/i");
	m_outputTree->Branch("n_true_totPh_0_70", &m_n_true_totPh_0_70, "n_true_totPh_0_70/i");
	m_outputTree->Branch("n_true_totK0L_0_70", &m_n_true_totK0L_0_70, "n_true_totK0L_0_70/i");
	m_outputTree->Branch("n_true_totE_0_70", &m_n_true_totE_0_70, "n_true_totE_0_70/i");
	m_outputTree->Branch("n_true_totMu_0_70", &m_n_true_totMu_0_70, "n_true_totMu_0_70/i");
	m_outputTree->Branch("n_true_totN_0_70", &m_n_true_totN_0_70, "n_true_totN_0_70/i");
	m_outputTree->Branch("n_true_totP_0_70", &m_n_true_totP_0_70, "n_true_totP_0_70/i");
	m_outputTree->Branch("n_true_totK_0_70", &m_n_true_totK_0_70, "n_true_totK_0_70/i");
	m_outputTree->Branch("n_true_totOtherCH_0_70", &m_n_true_totOtherCH_0_70, "n_true_totOtherCH_0_70/i");
	m_outputTree->Branch("n_true_totOtherNeut_0_70", &m_n_true_totOtherNeut_0_70, "n_true_totOtherNeut_0_70/i");
      }
      m_outputTree->Branch("n_0_2_totPFO", &m_n_0_2_totPFO, "n_0_2_totPFO/i");
      m_outputTree->Branch("n_0_2_totPi", &m_n_0_2_totPi, "n_0_2_totPi/i");
      m_outputTree->Branch("n_0_2_totPh", &m_n_0_2_totPh, "n_0_2_totPh/i");
      m_outputTree->Branch("n_0_2_totE", &m_n_0_2_totE, "n_0_2_totE/i");
      m_outputTree->Branch("n_0_2_totMu", &m_n_0_2_totMu, "n_0_2_totMu/i");
      m_outputTree->Branch("n_0_2_totN", &m_n_0_2_totN, "n_0_2_totN/i");
      m_outputTree->Branch("n_0_2_totPFO_0_95", &m_n_0_2_totPFO_0_95, "n_0_2_totPFO_0_95/i");
      m_outputTree->Branch("n_0_2_totPi_0_95", &m_n_0_2_totPi_0_95, "n_0_2_totPi_0_95/i");
      m_outputTree->Branch("n_0_2_totPh_0_95", &m_n_0_2_totPh_0_95, "n_0_2_totPh_0_95/i");
      m_outputTree->Branch("n_0_2_totE_0_95", &m_n_0_2_totE_0_95, "n_0_2_totE_0_95/i");
      m_outputTree->Branch("n_0_2_totMu_0_95", &m_n_0_2_totMu_0_95, "n_0_2_totMu_0_95/i");
      m_outputTree->Branch("n_0_2_totN_0_95", &m_n_0_2_totN_0_95, "n_0_2_totN_0_95/i");
      m_outputTree->Branch("n_0_2_totPFO_0_70", &m_n_0_2_totPFO_0_70, "n_0_2_totPFO_0_70/i");
      m_outputTree->Branch("n_0_2_totPi_0_70", &m_n_0_2_totPi_0_70, "n_0_2_totPi_0_70/i");
      m_outputTree->Branch("n_0_2_totPh_0_70", &m_n_0_2_totPh_0_70, "n_0_2_totPh_0_70/i");
      m_outputTree->Branch("n_0_2_totE_0_70", &m_n_0_2_totE_0_70, "n_0_2_totE_0_70/i");
      m_outputTree->Branch("n_0_2_totMu_0_70", &m_n_0_2_totMu_0_70, "n_0_2_totMu_0_70/i");
      m_outputTree->Branch("n_0_2_totN_0_70", &m_n_0_2_totN_0_70, "n_0_2_totN_0_70/i");
      
      if(!m_ignoreGen){
	m_outputTree->Branch("n_true_0_2_totAll", &m_n_true_0_2_totAll, "n_true_0_2_totAll/i");
	m_outputTree->Branch("n_true_0_2_totInv", &m_n_true_0_2_totInv, "n_true_0_2_totInv/i");
	m_outputTree->Branch("n_true_0_2_totPi", &m_n_true_0_2_totPi, "n_true_0_2_totPi/i");
	m_outputTree->Branch("n_true_0_2_totPh", &m_n_true_0_2_totPh, "n_true_0_2_totPh/i");
	m_outputTree->Branch("n_true_0_2_totK0L", &m_n_true_0_2_totK0L, "n_true_0_2_totK0L/i");
	m_outputTree->Branch("n_true_0_2_totE", &m_n_true_0_2_totE, "n_true_0_2_totE/i");
	m_outputTree->Branch("n_true_0_2_totMu", &m_n_true_0_2_totMu, "n_true_0_2_totMu/i");
	m_outputTree->Branch("n_true_0_2_totN", &m_n_true_0_2_totN, "n_true_0_2_totN/i");
	m_outputTree->Branch("n_true_0_2_totP", &m_n_true_0_2_totP, "n_true_0_2_totP/i");
	m_outputTree->Branch("n_true_0_2_totK", &m_n_true_0_2_totK, "n_true_0_2_totK/i");
	m_outputTree->Branch("n_true_0_2_totOtherCH", &m_n_true_0_2_totOtherCH, "n_true_0_2_totOtherCH/i");
	m_outputTree->Branch("n_true_0_2_totOtherNeut", &m_n_true_0_2_totOtherNeut, "n_true_0_2_totOtherNeut/i");
	m_outputTree->Branch("n_true_0_2_totAll_0_95", &m_n_true_0_2_totAll_0_95, "n_true_0_2_totAll_0_95/i");
	m_outputTree->Branch("n_true_0_2_totInv_0_95", &m_n_true_0_2_totInv_0_95, "n_true_0_2_totInv_0_95/i");
	m_outputTree->Branch("n_true_0_2_totPi_0_95", &m_n_true_0_2_totPi_0_95, "n_true_0_2_totPi_0_95/i");
	m_outputTree->Branch("n_true_0_2_totPh_0_95", &m_n_true_0_2_totPh_0_95, "n_true_0_2_totPh_0_95/i");
	m_outputTree->Branch("n_true_0_2_totK0L_0_95", &m_n_true_0_2_totK0L_0_95, "n_true_0_2_totK0L_0_95/i");
	m_outputTree->Branch("n_true_0_2_totE_0_95", &m_n_true_0_2_totE_0_95, "n_true_0_2_totE_0_95/i");
	m_outputTree->Branch("n_true_0_2_totMu_0_95", &m_n_true_0_2_totMu_0_95, "n_true_0_2_totMu_0_95/i");
	m_outputTree->Branch("n_true_0_2_totN_0_95", &m_n_true_0_2_totN_0_95, "n_true_0_2_totN_0_95/i");
	m_outputTree->Branch("n_true_0_2_totP_0_95", &m_n_true_0_2_totP_0_95, "n_true_0_2_totP_0_95/i");
	m_outputTree->Branch("n_true_0_2_totK_0_95", &m_n_true_0_2_totK_0_95, "n_true_0_2_totK_0_95/i");
	m_outputTree->Branch("n_true_0_2_totOtherCH_0_95", &m_n_true_0_2_totOtherCH_0_95, "n_true_0_2_totOtherCH_0_95/i");
	m_outputTree->Branch("n_true_0_2_totOtherNeut_0_95", &m_n_true_0_2_totOtherNeut_0_95, "n_true_0_2_totOtherNeut_0_95/i");
	m_outputTree->Branch("n_true_0_2_totAll_0_70", &m_n_true_0_2_totAll_0_70, "n_true_0_2_totAll_0_70/i");
	m_outputTree->Branch("n_true_0_2_totInv_0_70", &m_n_true_0_2_totInv_0_70, "n_true_0_2_totInv_0_70/i");
	m_outputTree->Branch("n_true_0_2_totPi_0_70", &m_n_true_0_2_totPi_0_70, "n_true_0_2_totPi_0_70/i");
	m_outputTree->Branch("n_true_0_2_totPh_0_70", &m_n_true_0_2_totPh_0_70, "n_true_0_2_totPh_0_70/i");
	m_outputTree->Branch("n_true_0_2_totK0L_0_70", &m_n_true_0_2_totK0L_0_70, "n_true_0_2_totK0L_0_70/i");
	m_outputTree->Branch("n_true_0_2_totE_0_70", &m_n_true_0_2_totE_0_70, "n_true_0_2_totE_0_70/i");
	m_outputTree->Branch("n_true_0_2_totMu_0_70", &m_n_true_0_2_totMu_0_70, "n_true_0_2_totMu_0_70/i");
	m_outputTree->Branch("n_true_0_2_totN_0_70", &m_n_true_0_2_totN_0_70, "n_true_0_2_totN_0_70/i");
	m_outputTree->Branch("n_true_0_2_totP_0_70", &m_n_true_0_2_totP_0_70, "n_true_0_2_totP_0_70/i");
	m_outputTree->Branch("n_true_0_2_totK_0_70", &m_n_true_0_2_totK_0_70, "n_true_0_2_totK_0_70/i");
	m_outputTree->Branch("n_true_0_2_totOtherCH_0_70", &m_n_true_0_2_totOtherCH_0_70, "n_true_0_2_totOtherCH_0_70/i");
	m_outputTree->Branch("n_true_0_2_totOtherNeut_0_70", &m_n_true_0_2_totOtherNeut_0_70, "n_true_0_2_totOtherNeut_0_70/i");
      }
      m_outputTree->Branch("n_2_10_totPFO", &m_n_2_10_totPFO, "n_2_10_totPFO/i");
      m_outputTree->Branch("n_2_10_totPi", &m_n_2_10_totPi, "n_2_10_totPi/i");
      m_outputTree->Branch("n_2_10_totPh", &m_n_2_10_totPh, "n_2_10_totPh/i");
      m_outputTree->Branch("n_2_10_totE", &m_n_2_10_totE, "n_2_10_totE/i");
      m_outputTree->Branch("n_2_10_totMu", &m_n_2_10_totMu, "n_2_10_totMu/i");
      m_outputTree->Branch("n_2_10_totN", &m_n_2_10_totN, "n_2_10_totN/i");
      m_outputTree->Branch("n_2_10_totPFO_0_95", &m_n_2_10_totPFO_0_95, "n_2_10_totPFO_0_95/i");
      m_outputTree->Branch("n_2_10_totPi_0_95", &m_n_2_10_totPi_0_95, "n_2_10_totPi_0_95/i");
      m_outputTree->Branch("n_2_10_totPh_0_95", &m_n_2_10_totPh_0_95, "n_2_10_totPh_0_95/i");
      m_outputTree->Branch("n_2_10_totE_0_95", &m_n_2_10_totE_0_95, "n_2_10_totE_0_95/i");
      m_outputTree->Branch("n_2_10_totMu_0_95", &m_n_2_10_totMu_0_95, "n_2_10_totMu_0_95/i");
      m_outputTree->Branch("n_2_10_totN_0_95", &m_n_2_10_totN_0_95, "n_2_10_totN_0_95/i");
      m_outputTree->Branch("n_2_10_totPFO_0_70", &m_n_2_10_totPFO_0_70, "n_2_10_totPFO_0_70/i");
      m_outputTree->Branch("n_2_10_totPi_0_70", &m_n_2_10_totPi_0_70, "n_2_10_totPi_0_70/i");
      m_outputTree->Branch("n_2_10_totPh_0_70", &m_n_2_10_totPh_0_70, "n_2_10_totPh_0_70/i");
      m_outputTree->Branch("n_2_10_totE_0_70", &m_n_2_10_totE_0_70, "n_2_10_totE_0_70/i");
      m_outputTree->Branch("n_2_10_totMu_0_70", &m_n_2_10_totMu_0_70, "n_2_10_totMu_0_70/i");
      m_outputTree->Branch("n_2_10_totN_0_70", &m_n_2_10_totN_0_70, "n_2_10_totN_0_70/i");
      if(!m_ignoreGen){
	m_outputTree->Branch("n_true_2_10_totAll", &m_n_true_2_10_totAll, "n_true_2_10_totAll/i");
	m_outputTree->Branch("n_true_2_10_totInv", &m_n_true_2_10_totInv, "n_true_2_10_totInv/i");
	m_outputTree->Branch("n_true_2_10_totPi", &m_n_true_2_10_totPi, "n_true_2_10_totPi/i");
	m_outputTree->Branch("n_true_2_10_totPh", &m_n_true_2_10_totPh, "n_true_2_10_totPh/i");
	m_outputTree->Branch("n_true_2_10_totK0L", &m_n_true_2_10_totK0L, "n_true_2_10_totK0L/i");
	m_outputTree->Branch("n_true_2_10_totE", &m_n_true_2_10_totE, "n_true_2_10_totE/i");
	m_outputTree->Branch("n_true_2_10_totMu", &m_n_true_2_10_totMu, "n_true_2_10_totMu/i");
	m_outputTree->Branch("n_true_2_10_totN", &m_n_true_2_10_totN, "n_true_2_10_totN/i");
	m_outputTree->Branch("n_true_2_10_totP", &m_n_true_2_10_totP, "n_true_2_10_totP/i");
	m_outputTree->Branch("n_true_2_10_totK", &m_n_true_2_10_totK, "n_true_2_10_totK/i");
	m_outputTree->Branch("n_true_2_10_totOtherCH", &m_n_true_2_10_totOtherCH, "n_true_2_10_totOtherCH/i");
	m_outputTree->Branch("n_true_2_10_totOtherNeut", &m_n_true_2_10_totOtherNeut, "n_true_2_10_totOtherNeut/i");
	m_outputTree->Branch("n_true_2_10_totAll_0_95", &m_n_true_2_10_totAll_0_95, "n_true_2_10_totAll_0_95/i");
	m_outputTree->Branch("n_true_2_10_totInv_0_95", &m_n_true_2_10_totInv_0_95, "n_true_2_10_totInv_0_95/i");
	m_outputTree->Branch("n_true_2_10_totPi_0_95", &m_n_true_2_10_totPi_0_95, "n_true_2_10_totPi_0_95/i");
	m_outputTree->Branch("n_true_2_10_totPh_0_95", &m_n_true_2_10_totPh_0_95, "n_true_2_10_totPh_0_95/i");
	m_outputTree->Branch("n_true_2_10_totK0L_0_95", &m_n_true_2_10_totK0L_0_95, "n_true_2_10_totK0L_0_95/i");
	m_outputTree->Branch("n_true_2_10_totE_0_95", &m_n_true_2_10_totE_0_95, "n_true_2_10_totE_0_95/i");
	m_outputTree->Branch("n_true_2_10_totMu_0_95", &m_n_true_2_10_totMu_0_95, "n_true_2_10_totMu_0_95/i");
	m_outputTree->Branch("n_true_2_10_totN_0_95", &m_n_true_2_10_totN_0_95, "n_true_2_10_totN_0_95/i");
	m_outputTree->Branch("n_true_2_10_totP_0_95", &m_n_true_2_10_totP_0_95, "n_true_2_10_totP_0_95/i");
	m_outputTree->Branch("n_true_2_10_totK_0_95", &m_n_true_2_10_totK_0_95, "n_true_2_10_totK_0_95/i");
	m_outputTree->Branch("n_true_2_10_totOtherCH_0_95", &m_n_true_2_10_totOtherCH_0_95, "n_true_2_10_totOtherCH_0_95/i");
	m_outputTree->Branch("n_true_2_10_totOtherNeut_0_95", &m_n_true_2_10_totOtherNeut_0_95, "n_true_2_10_totOtherNeut_0_95/i");
	m_outputTree->Branch("n_true_2_10_totAll_0_70", &m_n_true_2_10_totAll_0_70, "n_true_2_10_totAll_0_70/i");
	m_outputTree->Branch("n_true_2_10_totInv_0_70", &m_n_true_2_10_totInv_0_70, "n_true_2_10_totInv_0_70/i");
	m_outputTree->Branch("n_true_2_10_totPi_0_70", &m_n_true_2_10_totPi_0_70, "n_true_2_10_totPi_0_70/i");
	m_outputTree->Branch("n_true_2_10_totPh_0_70", &m_n_true_2_10_totPh_0_70, "n_true_2_10_totPh_0_70/i");
	m_outputTree->Branch("n_true_2_10_totK0L_0_70", &m_n_true_2_10_totK0L_0_70, "n_true_2_10_totK0L_0_70/i");
	m_outputTree->Branch("n_true_2_10_totE_0_70", &m_n_true_2_10_totE_0_70, "n_true_2_10_totE_0_70/i");
	m_outputTree->Branch("n_true_2_10_totMu_0_70", &m_n_true_2_10_totMu_0_70, "n_true_2_10_totMu_0_70/i");
	m_outputTree->Branch("n_true_2_10_totN_0_70", &m_n_true_2_10_totN_0_70, "n_true_2_10_totN_0_70/i");
	m_outputTree->Branch("n_true_2_10_totP_0_70", &m_n_true_2_10_totP_0_70, "n_true_2_10_totP_0_70/i");
	m_outputTree->Branch("n_true_2_10_totK_0_70", &m_n_true_2_10_totK_0_70, "n_true_2_10_totK_0_70/i");
	m_outputTree->Branch("n_true_2_10_totOtherCH_0_70", &m_n_true_2_10_totOtherCH_0_70, "n_true_2_10_totOtherCH_0_70/i");
	m_outputTree->Branch("n_true_2_10_totOtherNeut_0_70", &m_n_true_2_10_totOtherNeut_0_70, "n_true_2_10_totOtherNeut_0_70/i");
      }
      m_outputTree->Branch("n_10_50_totPFO", &m_n_10_50_totPFO, "n_10_50_totPFO/i");
      m_outputTree->Branch("n_10_50_totPi", &m_n_10_50_totPi, "n_10_50_totPi/i");
      m_outputTree->Branch("n_10_50_totPh", &m_n_10_50_totPh, "n_10_50_totPh/i");
      m_outputTree->Branch("n_10_50_totE", &m_n_10_50_totE, "n_10_50_totE/i");
      m_outputTree->Branch("n_10_50_totMu", &m_n_10_50_totMu, "n_10_50_totMu/i");
      m_outputTree->Branch("n_10_50_totN", &m_n_10_50_totN, "n_10_50_totN/i");
      m_outputTree->Branch("n_10_50_totPFO_0_95", &m_n_10_50_totPFO_0_95, "n_10_50_totPFO_0_95/i");
      m_outputTree->Branch("n_10_50_totPi_0_95", &m_n_10_50_totPi_0_95, "n_10_50_totPi_0_95/i");
      m_outputTree->Branch("n_10_50_totPh_0_95", &m_n_10_50_totPh_0_95, "n_10_50_totPh_0_95/i");
      m_outputTree->Branch("n_10_50_totE_0_95", &m_n_10_50_totE_0_95, "n_10_50_totE_0_95/i");
      m_outputTree->Branch("n_10_50_totMu_0_95", &m_n_10_50_totMu_0_95, "n_10_50_totMu_0_95/i");
      m_outputTree->Branch("n_10_50_totN_0_95", &m_n_10_50_totN_0_95, "n_10_50_totN_0_95/i");
      m_outputTree->Branch("n_10_50_totPFO_0_70", &m_n_10_50_totPFO_0_70, "n_10_50_totPFO_0_70/i");
      m_outputTree->Branch("n_10_50_totPi_0_70", &m_n_10_50_totPi_0_70, "n_10_50_totPi_0_70/i");
      m_outputTree->Branch("n_10_50_totPh_0_70", &m_n_10_50_totPh_0_70, "n_10_50_totPh_0_70/i");
      m_outputTree->Branch("n_10_50_totE_0_70", &m_n_10_50_totE_0_70, "n_10_50_totE_0_70/i");
      m_outputTree->Branch("n_10_50_totMu_0_70", &m_n_10_50_totMu_0_70, "n_10_50_totMu_0_70/i");
      m_outputTree->Branch("n_10_50_totN_0_70", &m_n_10_50_totN_0_70, "n_10_50_totN_0_70/i");
      if(!m_ignoreGen){
	m_outputTree->Branch("n_true_10_50_totAll", &m_n_true_10_50_totAll, "n_true_10_50_totAll/i");
	m_outputTree->Branch("n_true_10_50_totInv", &m_n_true_10_50_totInv, "n_true_10_50_totInv/i");
	m_outputTree->Branch("n_true_10_50_totPi", &m_n_true_10_50_totPi, "n_true_10_50_totPi/i");
	m_outputTree->Branch("n_true_10_50_totPh", &m_n_true_10_50_totPh, "n_true_10_50_totPh/i");
	m_outputTree->Branch("n_true_10_50_totK0L", &m_n_true_10_50_totK0L, "n_true_10_50_totK0L/i");
	m_outputTree->Branch("n_true_10_50_totE", &m_n_true_10_50_totE, "n_true_10_50_totE/i");
	m_outputTree->Branch("n_true_10_50_totMu", &m_n_true_10_50_totMu, "n_true_10_50_totMu/i");
	m_outputTree->Branch("n_true_10_50_totN", &m_n_true_10_50_totN, "n_true_10_50_totN/i");
	m_outputTree->Branch("n_true_10_50_totP", &m_n_true_10_50_totP, "n_true_10_50_totP/i");
	m_outputTree->Branch("n_true_10_50_totK", &m_n_true_10_50_totK, "n_true_10_50_totK/i");
	m_outputTree->Branch("n_true_10_50_totOtherCH", &m_n_true_10_50_totOtherCH, "n_true_10_50_totOtherCH/i");
	m_outputTree->Branch("n_true_10_50_totOtherNeut", &m_n_true_10_50_totOtherNeut, "n_true_10_50_totOtherNeut/i");
	m_outputTree->Branch("n_true_10_50_totAll_0_95", &m_n_true_10_50_totAll_0_95, "n_true_10_50_totAll_0_95/i");
	m_outputTree->Branch("n_true_10_50_totInv_0_95", &m_n_true_10_50_totInv_0_95, "n_true_10_50_totInv_0_95/i");
	m_outputTree->Branch("n_true_10_50_totPi_0_95", &m_n_true_10_50_totPi_0_95, "n_true_10_50_totPi_0_95/i");
	m_outputTree->Branch("n_true_10_50_totPh_0_95", &m_n_true_10_50_totPh_0_95, "n_true_10_50_totPh_0_95/i");
	m_outputTree->Branch("n_true_10_50_totK0L_0_95", &m_n_true_10_50_totK0L_0_95, "n_true_10_50_totK0L_0_95/i");
	m_outputTree->Branch("n_true_10_50_totE_0_95", &m_n_true_10_50_totE_0_95, "n_true_10_50_totE_0_95/i");
	m_outputTree->Branch("n_true_10_50_totMu_0_95", &m_n_true_10_50_totMu_0_95, "n_true_10_50_totMu_0_95/i");
	m_outputTree->Branch("n_true_10_50_totN_0_95", &m_n_true_10_50_totN_0_95, "n_true_10_50_totN_0_95/i");
	m_outputTree->Branch("n_true_10_50_totP_0_95", &m_n_true_10_50_totP_0_95, "n_true_10_50_totP_0_95/i");
	m_outputTree->Branch("n_true_10_50_totK_0_95", &m_n_true_10_50_totK_0_95, "n_true_10_50_totK_0_95/i");
	m_outputTree->Branch("n_true_10_50_totOtherCH_0_95", &m_n_true_10_50_totOtherCH_0_95, "n_true_10_50_totOtherCH_0_95/i");
	m_outputTree->Branch("n_true_10_50_totOtherNeut_0_95", &m_n_true_10_50_totOtherNeut_0_95, "n_true_10_50_totOtherNeut_0_95/i");
	m_outputTree->Branch("n_true_10_50_totAll_0_70", &m_n_true_10_50_totAll_0_70, "n_true_10_50_totAll_0_70/i");
	m_outputTree->Branch("n_true_10_50_totInv_0_70", &m_n_true_10_50_totInv_0_70, "n_true_10_50_totInv_0_70/i");
	m_outputTree->Branch("n_true_10_50_totPi_0_70", &m_n_true_10_50_totPi_0_70, "n_true_10_50_totPi_0_70/i");
	m_outputTree->Branch("n_true_10_50_totPh_0_70", &m_n_true_10_50_totPh_0_70, "n_true_10_50_totPh_0_70/i");
	m_outputTree->Branch("n_true_10_50_totK0L_0_70", &m_n_true_10_50_totK0L_0_70, "n_true_10_50_totK0L_0_70/i");
	m_outputTree->Branch("n_true_10_50_totE_0_70", &m_n_true_10_50_totE_0_70, "n_true_10_50_totE_0_70/i");
	m_outputTree->Branch("n_true_10_50_totMu_0_70", &m_n_true_10_50_totMu_0_70, "n_true_10_50_totMu_0_70/i");
	m_outputTree->Branch("n_true_10_50_totN_0_70", &m_n_true_10_50_totN_0_70, "n_true_10_50_totN_0_70/i");
	m_outputTree->Branch("n_true_10_50_totP_0_70", &m_n_true_10_50_totP_0_70, "n_true_10_50_totP_0_70/i");
	m_outputTree->Branch("n_true_10_50_totK_0_70", &m_n_true_10_50_totK_0_70, "n_true_10_50_totK_0_70/i");
	m_outputTree->Branch("n_true_10_50_totOtherCH_0_70", &m_n_true_10_50_totOtherCH_0_70, "n_true_10_50_totOtherCH_0_70/i");
	m_outputTree->Branch("n_true_10_50_totOtherNeut_0_70", &m_n_true_10_50_totOtherNeut_0_70, "n_true_10_50_totOtherNeut_0_70/i");
      }
      m_outputTree->Branch("n_50_totPFO", &m_n_50_totPFO, "n_50_totPFO/i");
      m_outputTree->Branch("n_50_totPi", &m_n_50_totPi, "n_50_totPi/i");
      m_outputTree->Branch("n_50_totPh", &m_n_50_totPh, "n_50_totPh/i");
      m_outputTree->Branch("n_50_totE", &m_n_50_totE, "n_50_totE/i");
      m_outputTree->Branch("n_50_totMu", &m_n_50_totMu, "n_50_totMu/i");
      m_outputTree->Branch("n_50_totN", &m_n_50_totN, "n_50_totN/i");
      m_outputTree->Branch("n_50_totPFO_0_95", &m_n_50_totPFO_0_95, "n_50_totPFO_0_95/i");
      m_outputTree->Branch("n_50_totPi_0_95", &m_n_50_totPi_0_95, "n_50_totPi_0_95/i");
      m_outputTree->Branch("n_50_totPh_0_95", &m_n_50_totPh_0_95, "n_50_totPh_0_95/i");
      m_outputTree->Branch("n_50_totE_0_95", &m_n_50_totE_0_95, "n_50_totE_0_95/i");
      m_outputTree->Branch("n_50_totMu_0_95", &m_n_50_totMu_0_95, "n_50_totMu_0_95/i");
      m_outputTree->Branch("n_50_totN_0_95", &m_n_50_totN_0_95, "n_50_totN_0_95/i");
      m_outputTree->Branch("n_50_totPFO_0_70", &m_n_50_totPFO_0_70, "n_50_totPFO_0_70/i");
      m_outputTree->Branch("n_50_totPi_0_70", &m_n_50_totPi_0_70, "n_50_totPi_0_70/i");
      m_outputTree->Branch("n_50_totPh_0_70", &m_n_50_totPh_0_70, "n_50_totPh_0_70/i");
      m_outputTree->Branch("n_50_totE_0_70", &m_n_50_totE_0_70, "n_50_totE_0_70/i");
      m_outputTree->Branch("n_50_totMu_0_70", &m_n_50_totMu_0_70, "n_50_totMu_0_70/i");
      m_outputTree->Branch("n_50_totN_0_70", &m_n_50_totN_0_70, "n_50_totN_0_70/i");
      if(!m_ignoreGen){
	m_outputTree->Branch("n_true_50_totAll", &m_n_true_50_totAll, "n_true_50_totAll/i");
	m_outputTree->Branch("n_true_50_totInv", &m_n_true_50_totInv, "n_true_50_totInv/i");
	m_outputTree->Branch("n_true_50_totPi", &m_n_true_50_totPi, "n_true_50_totPi/i");
	m_outputTree->Branch("n_true_50_totPh", &m_n_true_50_totPh, "n_true_50_totPh/i");
	m_outputTree->Branch("n_true_50_totK0L", &m_n_true_50_totK0L, "n_true_50_totK0L/i");
	m_outputTree->Branch("n_true_50_totE", &m_n_true_50_totE, "n_true_50_totE/i");
	m_outputTree->Branch("n_true_50_totMu", &m_n_true_50_totMu, "n_true_50_totMu/i");
	m_outputTree->Branch("n_true_50_totN", &m_n_true_50_totN, "n_true_50_totN/i");
	m_outputTree->Branch("n_true_50_totP", &m_n_true_50_totP, "n_true_50_totP/i");
	m_outputTree->Branch("n_true_50_totK", &m_n_true_50_totK, "n_true_50_totK/i");
	m_outputTree->Branch("n_true_50_totOtherCH", &m_n_true_50_totOtherCH, "n_true_50_totOtherCH/i");
	m_outputTree->Branch("n_true_50_totOtherNeut", &m_n_true_50_totOtherNeut, "n_true_50_totOtherNeut/i");
	m_outputTree->Branch("n_true_50_totAll_0_95", &m_n_true_50_totAll_0_95, "n_true_50_totAll_0_95/i");
	m_outputTree->Branch("n_true_50_totInv_0_95", &m_n_true_50_totInv_0_95, "n_true_50_totInv_0_95/i");
	m_outputTree->Branch("n_true_50_totPi_0_95", &m_n_true_50_totPi_0_95, "n_true_50_totPi_0_95/i");
	m_outputTree->Branch("n_true_50_totPh_0_95", &m_n_true_50_totPh_0_95, "n_true_50_totPh_0_95/i");
	m_outputTree->Branch("n_true_50_totK0L_0_95", &m_n_true_50_totK0L_0_95, "n_true_50_totK0L_0_95/i");
	m_outputTree->Branch("n_true_50_totE_0_95", &m_n_true_50_totE_0_95, "n_true_50_totE_0_95/i");
	m_outputTree->Branch("n_true_50_totMu_0_95", &m_n_true_50_totMu_0_95, "n_true_50_totMu_0_95/i");
	m_outputTree->Branch("n_true_50_totN_0_95", &m_n_true_50_totN_0_95, "n_true_50_totN_0_95/i");
	m_outputTree->Branch("n_true_50_totP_0_95", &m_n_true_50_totP_0_95, "n_true_50_totP_0_95/i");
	m_outputTree->Branch("n_true_50_totK_0_95", &m_n_true_50_totK_0_95, "n_true_50_totK_0_95/i");
	m_outputTree->Branch("n_true_50_totOtherCH_0_95", &m_n_true_50_totOtherCH_0_95, "n_true_50_totOtherCH_0_95/i");
	m_outputTree->Branch("n_true_50_totOtherNeut_0_95", &m_n_true_50_totOtherNeut_0_95, "n_true_50_totOtherNeut_0_95/i");
	m_outputTree->Branch("n_true_50_totAll_0_70", &m_n_true_50_totAll_0_70, "n_true_50_totAll_0_70/i");
	m_outputTree->Branch("n_true_50_totInv_0_70", &m_n_true_50_totInv_0_70, "n_true_50_totInv_0_70/i");
	m_outputTree->Branch("n_true_50_totPi_0_70", &m_n_true_50_totPi_0_70, "n_true_50_totPi_0_70/i");
	m_outputTree->Branch("n_true_50_totPh_0_70", &m_n_true_50_totPh_0_70, "n_true_50_totPh_0_70/i");
	m_outputTree->Branch("n_true_50_totK0L_0_70", &m_n_true_50_totK0L_0_70, "n_true_50_totK0L_0_70/i");
	m_outputTree->Branch("n_true_50_totE_0_70", &m_n_true_50_totE_0_70, "n_true_50_totE_0_70/i");
	m_outputTree->Branch("n_true_50_totMu_0_70", &m_n_true_50_totMu_0_70, "n_true_50_totMu_0_70/i");
	m_outputTree->Branch("n_true_50_totN_0_70", &m_n_true_50_totN_0_70, "n_true_50_totN_0_70/i");
	m_outputTree->Branch("n_true_50_totP_0_70", &m_n_true_50_totP_0_70, "n_true_50_totP_0_70/i");
	m_outputTree->Branch("n_true_50_totK_0_70", &m_n_true_50_totK_0_70, "n_true_50_totK_0_70/i");
	m_outputTree->Branch("n_true_50_totOtherCH_0_70", &m_n_true_50_totOtherCH_0_70, "n_true_50_totOtherCH_0_70/i");
	m_outputTree->Branch("n_true_50_totOtherNeut_0_70", &m_n_true_50_totOtherNeut_0_70, "n_true_50_totOtherNeut_0_70/i");
      }
    }
    /*
    m_outputTree->Branch("truePi0Energy", "std::vector< float >", &m_truePi0Energy);
    m_outputTree->Branch("truePi0_Px", "std::vector< float >", &m_truePi0_Px);
    m_outputTree->Branch("truePi0_Py", "std::vector< float >", &m_truePi0_Py);
    m_outputTree->Branch("truePi0_Pz", "std::vector< float >", &m_truePi0_Pz);
    m_outputTree->Branch("truePi0_CosTheta", "std::vector< float >", &m_truePi0_CosTheta);
    m_outputTree->Branch("truePi0_Phi", "std::vector< float >", &m_truePi0_Phi);
    m_outputTree->Branch("truePi0_DRPh01", "std::vector< float >", &m_truePi0_DRPh01);
    m_outputTree->Branch("truePi0_CosAngPh01", "std::vector< float >", &m_truePi0_CosAngPh01);
    m_outputTree->Branch("truePi0_AngPh01", "std::vector< float >", &m_truePi0_AngPh01);
    m_outputTree->Branch("truePh0Energy", "std::vector< float >", &m_truePh0Energy);
    m_outputTree->Branch("truePh0_Px", "std::vector< float >", &m_truePh0_Px);
    m_outputTree->Branch("truePh0_Py", "std::vector< float >", &m_truePh0_Py);
    m_outputTree->Branch("truePh0_Pz", "std::vector< float >", &m_truePh0_Pz);
    m_outputTree->Branch("truePh0_CosTheta", "std::vector< float >", &m_truePh0_CosTheta);
    m_outputTree->Branch("truePh0_Phi", "std::vector< float >", &m_truePh0_Phi);
    m_outputTree->Branch("truePh0_DRMin", "std::vector< float >", &m_truePh0_DRMin);
    m_outputTree->Branch("truePh0_DRMinPDGID", "std::vector< float >", &m_truePh0_DRMinPDGID);
    m_outputTree->Branch("truePh0_DRMinE", "std::vector< float >", &m_truePh0_DRMinE);
    m_outputTree->Branch("truePh0_CosAngMax", "std::vector< float >", &m_truePh0_CosAngMax);
    m_outputTree->Branch("truePh0_CosAngMaxPDGID", "std::vector< float >", &m_truePh0_CosAngMaxPDGID);
    m_outputTree->Branch("truePh0_CosAngMaxE", "std::vector< float >", &m_truePh0_CosAngMaxE);
    m_outputTree->Branch("truePh0_CosAng0995E", "std::vector< float >", &m_truePh0_CosAng0995E);
    m_outputTree->Branch("truePh0_DR01E", "std::vector< float >", &m_truePh0_DR01E);
    m_outputTree->Branch("truePh1Energy", "std::vector< float >", &m_truePh1Energy);
    m_outputTree->Branch("truePh1_Px", "std::vector< float >", &m_truePh1_Px);
    m_outputTree->Branch("truePh1_Py", "std::vector< float >", &m_truePh1_Py);
    m_outputTree->Branch("truePh1_Pz", "std::vector< float >", &m_truePh1_Pz);
    m_outputTree->Branch("truePh1_CosTheta", "std::vector< float >", &m_truePh1_CosTheta);
    m_outputTree->Branch("truePh1_Phi", "std::vector< float >", &m_truePh1_Phi);
    m_outputTree->Branch("truePh1_DRMin", "std::vector< float >", &m_truePh1_DRMin);
    m_outputTree->Branch("truePh1_DRMinPDGID", "std::vector< float >", &m_truePh1_DRMinPDGID);
    m_outputTree->Branch("truePh1_DRMinE", "std::vector< float >", &m_truePh1_DRMinE);
    m_outputTree->Branch("truePh1_CosAngMax", "std::vector< float >", &m_truePh1_CosAngMax);
    m_outputTree->Branch("truePh1_CosAngMaxPDGID", "std::vector< float >", &m_truePh1_CosAngMaxPDGID);
    m_outputTree->Branch("truePh1_CosAngMaxE", "std::vector< float >", &m_truePh1_CosAngMaxE);
    m_outputTree->Branch("truePh1_CosAng0995E", "std::vector< float >", &m_truePh1_CosAng0995E);
    m_outputTree->Branch("truePh1_DR01E", "std::vector< float >", &m_truePh1_DR01E);

    m_outputTree->Branch("recoPhEnergy", "std::vector< float >", &m_recoPhEnergy);
    m_outputTree->Branch("recoPh_Px", "std::vector< float >", &m_recoPh_Px);
    m_outputTree->Branch("recoPh_Py", "std::vector< float >", &m_recoPh_Py);
    m_outputTree->Branch("recoPh_Pz", "std::vector< float >", &m_recoPh_Pz);
    m_outputTree->Branch("recoPh_CosTheta", "std::vector< float >", &m_recoPh_CosTheta);
    m_outputTree->Branch("recoPh_Theta", "std::vector< float >", &m_recoPh_Theta);
    m_outputTree->Branch("recoPh_Phi", "std::vector< float >", &m_recoPh_Phi);
    m_outputTree->Branch("recoPh_PDGID", "std::vector< float >", &m_recoPh_PDGID);
    m_outputTree->Branch("recoPh_E_EB", "std::vector< float >", &m_recoPh_E_EB);
    m_outputTree->Branch("recoPh_E_EE", "std::vector< float >", &m_recoPh_E_EE);
    m_outputTree->Branch("recoPh_E_EO", "std::vector< float >", &m_recoPh_E_EO);
    m_outputTree->Branch("recoPh_E_HB", "std::vector< float >", &m_recoPh_E_HB);
    m_outputTree->Branch("recoPh_E_HE", "std::vector< float >", &m_recoPh_E_HE);
    m_outputTree->Branch("recoPh_E_HO", "std::vector< float >", &m_recoPh_E_HO);
    m_outputTree->Branch("recoPh_firstLayerECAL", "std::vector< int >", &m_recoPh_firstLayerECAL);
    m_outputTree->Branch("recoPh_lastLayerECAL", "std::vector< int >", &m_recoPh_lastLayerECAL);
    m_outputTree->Branch("recoPh_nhitsEB", "std::vector< int >", &m_recoPh_nhitsEB);
    m_outputTree->Branch("recoPh_nhitsEE", "std::vector< int >", &m_recoPh_nhitsEE);
    m_outputTree->Branch("recoPh_nhitsEO", "std::vector< int >", &m_recoPh_nhitsEO);
    m_outputTree->Branch("recoPh_firstLayerHCAL", "std::vector< int >", &m_recoPh_firstLayerHCAL);
    m_outputTree->Branch("recoPh_lastLayerHCAL", "std::vector< int >", &m_recoPh_lastLayerHCAL);
    m_outputTree->Branch("recoPh_nhitsHB", "std::vector< int >", &m_recoPh_nhitsHB);
    m_outputTree->Branch("recoPh_nhitsHE", "std::vector< int >", &m_recoPh_nhitsHE);
    m_outputTree->Branch("recoPhnhitsHO", "std::vector< int >", &m_recoPh_nhitsHO);
    m_outputTree->Branch("recoPh_DRMin", "std::vector< float >", &m_recoPh_DRMin);
    m_outputTree->Branch("recoPh_DRMin_PDGID", "std::vector< float >", &m_recoPh_DRMin_PDGID);
    m_outputTree->Branch("recoPh_DRMin_E", "std::vector< float >", &m_recoPh_DRMin_E);
    m_outputTree->Branch("recoPh_CosAngMax", "std::vector< float >", &m_recoPh_CosAngMax);
    m_outputTree->Branch("recoPh_CosAngMax_PDGID", "std::vector< float >", &m_recoPh_CosAngMax_PDGID);
    m_outputTree->Branch("recoPh_CosAngMax_E", "std::vector< float >", &m_recoPh_CosAngMax_E);
    m_outputTree->Branch("recoPh_CosAng0995_E", "std::vector< float >", &m_recoPh_CosAng0995_E);
    m_outputTree->Branch("recoPh_DR01_E", "std::vector< float >", &m_recoPh_DR01_E);
    m_outputTree->Branch("recoPh_DRMin_Ph", "std::vector< float >", &m_recoPh_DRMin_Ph);
    m_outputTree->Branch("recoPh_DRMin_E_Ph", "std::vector< float >", &m_recoPh_DRMin_E_Ph);
    m_outputTree->Branch("recoPh_CosAngMax_Ph", "std::vector< float >", &m_recoPh_CosAngMax_Ph);
    m_outputTree->Branch("recoPh_CosAngMax_E_Ph", "std::vector< float >", &m_recoPh_CosAngMax_E_Ph);
    m_outputTree->Branch("recoPi0Cand_DRMin_M", "std::vector< float >", &m_recoPi0Cand_DRMin_M);
    m_outputTree->Branch("recoPi0Cand_CosMax_M", "std::vector< float >", &m_recoPi0Cand_CosAngMax_M);
    */
    m_outputTree->Branch("jetExclPx", "std::vector< float >", &m_jet_excl_Px); 
    m_outputTree->Branch("jetExclPy", "std::vector< float >", &m_jet_excl_Py); 
    m_outputTree->Branch("jetExclPz", "std::vector< float >", &m_jet_excl_Pz); 
    m_outputTree->Branch("jetExclE", "std::vector< float >", &m_jet_excl_E); 
    m_outputTree->Branch("jetExclPhi", "std::vector< float >", &m_jet_excl_Phi); 
    m_outputTree->Branch("jetExclCosTheta", "std::vector< float >", &m_jet_excl_CosTheta); 
    m_outputTree->Branch("jetExclPiE", "std::vector< float >", &m_jet_excl_piE);
    m_outputTree->Branch("jetExclPhE", "std::vector< float >", &m_jet_excl_phE);
    m_outputTree->Branch("jetExclElE", "std::vector< float >", &m_jet_excl_elE);
    m_outputTree->Branch("jetExclMuE", "std::vector< float >", &m_jet_excl_muE);
    m_outputTree->Branch("jetExclNE", "std::vector< float >", &m_jet_excl_nE);
    m_outputTree->Branch("jetExclNeutElseE", "std::vector< float >", &m_jet_excl_neutElseE);
    m_outputTree->Branch("jetExclChElseE", "std::vector< float >", &m_jet_excl_chElseE);
    m_outputTree->Branch("jetExclNeutMult", "std::vector< int >", &m_jet_excl_neutMult);
    m_outputTree->Branch("jetExclChMult", "std::vector< int >", &m_jet_excl_chMult);

    m_outputTree->Branch("jetInclPx", "std::vector< float >", &m_jet_incl_Px); 
    m_outputTree->Branch("jetInclPy", "std::vector< float >", &m_jet_incl_Py); 
    m_outputTree->Branch("jetInclPz", "std::vector< float >", &m_jet_incl_Pz); 
    m_outputTree->Branch("jetInclE", "std::vector< float >", &m_jet_incl_E); 
    m_outputTree->Branch("jetInclPhi", "std::vector< float >", &m_jet_incl_Phi); 
    m_outputTree->Branch("jetInclCosTheta", "std::vector< float >", &m_jet_incl_CosTheta); 
    m_outputTree->Branch("jetInclPiE", "std::vector< float >", &m_jet_incl_piE);
    m_outputTree->Branch("jetInclPhE", "std::vector< float >", &m_jet_incl_phE);
    m_outputTree->Branch("jetInclElE", "std::vector< float >", &m_jet_incl_elE);
    m_outputTree->Branch("jetInclMuE", "std::vector< float >", &m_jet_incl_muE);
    m_outputTree->Branch("jetInclNE", "std::vector< float >", &m_jet_incl_nE);
    m_outputTree->Branch("jetInclNeutElseE", "std::vector< float >", &m_jet_incl_neutElseE);
    m_outputTree->Branch("jetInclChElseE", "std::vector< float >", &m_jet_incl_chElseE);
    m_outputTree->Branch("jetInclNeutMult", "std::vector< int >", &m_jet_incl_neutMult);
    m_outputTree->Branch("jetInclChMult", "std::vector< int >", &m_jet_incl_chMult);
    m_outputTree->Branch("eventNumber", &m_eventNumber, "eventNumber/I");
    m_outputTree->Branch("runNumber", &m_runNumber, "runNumber/I");
    m_outputTree->Branch("eventCount", &m_eventCount, "eventCount/I");


    m_outputTree->Branch("E_trueNeut", &m_E_trueNeut, "E_trueNeut/F");
    m_outputTree->Branch("px_trueNeut", &m_px_trueNeut, "px_trueNeut/F");
    m_outputTree->Branch("py_trueNeut", &m_py_trueNeut, "py_trueNeut/F");
    m_outputTree->Branch("pz_trueNeut", &m_pz_trueNeut, "pz_trueNeut/F");

    m_outputTree->Branch("E_trueZ2", &m_E_trueZ2, "E_trueZ2/F");
    m_outputTree->Branch("px_trueZ2", &m_px_trueZ2, "px_truZ2/F");
    m_outputTree->Branch("py_trueZ2", &m_py_trueZ2, "py_trueZ2/F");
    m_outputTree->Branch("pz_trueZ2", &m_pz_trueZ2, "pz_trueZ2/F");
    
}


void PhotonStudy::processRunHeader( LCRunHeader*) {
	++m_runNumber ;
}

void PhotonStudy::processEvent( LCEvent* evt ) {
  //std::cout<<"at process event"<<std::endl;
  eventcount+=1;
  if(evt->getEventNumber()%50==0){
    std::cout<<"run/evt "<<evt->getRunNumber()<<"/"<<evt->getEventNumber()<<"/"<<eventcount<<std::endl;
  }


  m_eventNumber=evt->getEventNumber();
  m_runNumber=evt->getRunNumber();
  m_eventCount=eventcount;

    //std::cout<<"are we somewhere here or not 0 "<<std::endl;

    m_Z_mcE=-10;
    m_Z_mcNDaughter=-10;

    m_E_trueNeut=0;  
    m_px_trueNeut=0;
    m_py_trueNeut=0;
    m_pz_trueNeut=0;                                                                                                                                                                                                                             
    m_E_trueZ2=0;                   
    m_px_trueZ2=0;  
    m_py_trueZ2=0;   
    m_pz_trueZ2=0;   
    
    m_d1_mcPDGID=-10;
    m_d1_mcE=-10;
    m_d1_mcPx=-10;
    m_d1_mcPy=-10;
    m_d1_mcPz=-10;
    m_d1_mcMass=-10;
    m_d1_mcPhi=-10;
    m_d1_mcTheta=-10; 
    m_d1_mcCosTheta=-10;

    m_d2_mcPDGID=-10;
    m_d2_mcE=-10;
    m_d2_mcPx=-10;
    m_d2_mcPy=-10;
    m_d2_mcPz=-10;
    m_d2_mcMass=-10;
    m_d2_mcPhi=-10;
    m_d2_mcTheta=-10; 
    m_d2_mcCosTheta=-10;

    m_px_totPFO=0;
    m_px_totPi=0;
    m_px_totPh=0;
    m_px_totE=0;
    m_px_totMu=0;
    m_px_totN=0;

    m_py_totPFO=0;
    m_py_totPi=0;
    m_py_totPh=0;
    m_py_totE=0;
    m_py_totMu=0;
    m_py_totN=0;

    m_pz_totPFO=0;
    m_pz_totPi=0;
    m_pz_totPh=0;
    m_pz_totE=0;
    m_pz_totMu=0;
    m_pz_totN=0;

    m_E_totPFO=0;
    m_E_totPi=0;
    m_E_totPh=0;
    m_E_totE=0;
    m_E_totMu=0;
    m_E_totN=0;
    m_E_totPFO_0_95=0;
    m_E_totPi_0_95=0;
    m_E_totPh_0_95=0;
    m_E_totE_0_95=0;
    m_E_totMu_0_95=0;
    m_E_totN_0_95=0;
    m_E_totPFO_0_70=0;
    m_E_totPi_0_70=0;
    m_E_totPh_0_70=0;
    m_E_totE_0_70=0;
    m_E_totMu_0_70=0;
    m_E_totN_0_70=0;

    m_px_true_totAll=0;
    m_px_true_totInv=0;
    m_px_true_totPi=0;
    m_px_true_totPh=0;
    m_px_true_totK0L=0;
    m_px_true_totE=0;
    m_px_true_totMu=0;
    m_px_true_totN=0;
    m_px_true_totK=0;
    m_px_true_totP=0;

    m_py_true_totAll=0;
    m_py_true_totInv=0;
    m_py_true_totPi=0;
    m_py_true_totPh=0;
    m_py_true_totK0L=0;
    m_py_true_totE=0;
    m_py_true_totMu=0;
    m_py_true_totN=0;
    m_py_true_totK=0;
    m_py_true_totP=0;

    m_pz_true_totAll=0;
    m_pz_true_totInv=0;
    m_pz_true_totPi=0;
    m_pz_true_totPh=0;
    m_pz_true_totK0L=0;
    m_pz_true_totE=0;
    m_pz_true_totMu=0;
    m_pz_true_totN=0;
    m_pz_true_totK=0;
    m_pz_true_totP=0;


    m_E_true_totAll=0;
    m_E_true_totInv=0;
    m_E_true_totPi=0;
    m_E_true_totPh=0;
    m_E_true_totK0L=0;
    m_E_true_totE=0;
    m_E_true_totMu=0;
    m_E_true_totN=0;
    m_E_true_totK=0;
    m_E_true_totP=0;
    m_E_true_totOtherCH=0;
    m_E_true_totOtherNeut=0;
    m_E_true_totAll_0_95=0;
    m_E_true_totInv_0_95=0;
    m_E_true_totPi_0_95=0;
    m_E_true_totPh_0_95=0;
    m_E_true_totK0L_0_95=0;
    m_E_true_totE_0_95=0;
    m_E_true_totMu_0_95=0;
    m_E_true_totN_0_95=0;
    m_E_true_totK_0_95=0;
    m_E_true_totP_0_95=0;
    m_E_true_totOtherCH_0_95=0;
    m_E_true_totOtherNeut_0_95=0;
    m_E_true_totAll_0_70=0;
    m_E_true_totInv_0_70=0;
    m_E_true_totPi_0_70=0;
    m_E_true_totPh_0_70=0;
    m_E_true_totK0L_0_70=0;
    m_E_true_totE_0_70=0;
    m_E_true_totMu_0_70=0;
    m_E_true_totN_0_70=0;
    m_E_true_totK_0_70=0;
    m_E_true_totP_0_70=0;
    m_E_true_totOtherCH_0_70=0;
    m_E_true_totOtherNeut_0_70=0;

    m_E_0_2_totPFO=0;
    m_E_0_2_totPi=0;
    m_E_0_2_totPh=0;
    m_E_0_2_totE=0;
    m_E_0_2_totMu=0;
    m_E_0_2_totN=0;
    m_E_0_2_totPFO_0_95=0;
    m_E_0_2_totPi_0_95=0;
    m_E_0_2_totPh_0_95=0;
    m_E_0_2_totE_0_95=0;
    m_E_0_2_totMu_0_95=0;
    m_E_0_2_totN_0_95=0;
    m_E_0_2_totPFO_0_70=0;
    m_E_0_2_totPi_0_70=0;
    m_E_0_2_totPh_0_70=0;
    m_E_0_2_totE_0_70=0;
    m_E_0_2_totMu_0_70=0;
    m_E_0_2_totN_0_70=0;

    m_E_true_0_2_totAll=0;
    m_E_true_0_2_totInv=0;
    m_E_true_0_2_totPi=0;
    m_E_true_0_2_totPh=0;
    m_E_true_0_2_totE=0;
    m_E_true_0_2_totMu=0;
    m_E_true_0_2_totN=0;
    m_E_true_0_2_totK=0;
    m_E_true_0_2_totP=0;
    m_E_true_0_2_totOtherCH=0;
    m_E_true_0_2_totOtherNeut=0;
    m_E_true_0_2_totAll_0_95=0;
    m_E_true_0_2_totInv_0_95=0;
    m_E_true_0_2_totPi_0_95=0;
    m_E_true_0_2_totPh_0_95=0;
    m_E_true_0_2_totK0L_0_95=0;
    m_E_true_0_2_totE_0_95=0;
    m_E_true_0_2_totMu_0_95=0;
    m_E_true_0_2_totN_0_95=0;
    m_E_true_0_2_totK_0_95=0;
    m_E_true_0_2_totP_0_95=0;
    m_E_true_0_2_totOtherCH_0_95=0;
    m_E_true_0_2_totOtherNeut_0_95=0;
    m_E_true_0_2_totAll_0_70=0;
    m_E_true_0_2_totInv_0_70=0;
    m_E_true_0_2_totPi_0_70=0;
    m_E_true_0_2_totPh_0_70=0;
    m_E_true_0_2_totK0L_0_70=0;
    m_E_true_0_2_totE_0_70=0;
    m_E_true_0_2_totMu_0_70=0;
    m_E_true_0_2_totN_0_70=0;
    m_E_true_0_2_totK_0_70=0;
    m_E_true_0_2_totP_0_70=0;
    m_E_true_0_2_totOtherCH_0_70=0;
    m_E_true_0_2_totOtherNeut_0_70=0;

    m_E_2_10_totPFO=0;
    m_E_2_10_totPi=0;
    m_E_2_10_totPh=0;
    m_E_2_10_totE=0;
    m_E_2_10_totMu=0;
    m_E_2_10_totN=0;
    m_E_2_10_totPFO_0_95=0;
    m_E_2_10_totPi_0_95=0;
    m_E_2_10_totPh_0_95=0;
    m_E_2_10_totE_0_95=0;
    m_E_2_10_totMu_0_95=0;
    m_E_2_10_totN_0_95=0;
    m_E_2_10_totPFO_0_70=0;
    m_E_2_10_totPi_0_70=0;
    m_E_2_10_totPh_0_70=0;
    m_E_2_10_totE_0_70=0;
    m_E_2_10_totMu_0_70=0;
    m_E_2_10_totN_0_70=0;

    m_E_true_2_10_totAll=0;
    m_E_true_2_10_totInv=0;
    m_E_true_2_10_totPi=0;
    m_E_true_2_10_totPh=0;
    m_E_true_2_10_totK0L=0;
    m_E_true_2_10_totE=0;
    m_E_true_2_10_totMu=0;
    m_E_true_2_10_totN=0;
    m_E_true_2_10_totK=0;
    m_E_true_2_10_totP=0;
    m_E_true_2_10_totOtherCH=0;
    m_E_true_2_10_totOtherNeut=0;
    m_E_true_2_10_totAll_0_95=0;
    m_E_true_2_10_totInv_0_95=0;
    m_E_true_2_10_totPi_0_95=0;
    m_E_true_2_10_totPh_0_95=0;
    m_E_true_2_10_totK0L_0_95=0;
    m_E_true_2_10_totE_0_95=0;
    m_E_true_2_10_totMu_0_95=0;
    m_E_true_2_10_totN_0_95=0;
    m_E_true_2_10_totK_0_95=0;
    m_E_true_2_10_totP_0_95=0;
    m_E_true_2_10_totOtherCH_0_95=0;
    m_E_true_2_10_totOtherNeut_0_95=0;
    m_E_true_2_10_totAll_0_70=0;
    m_E_true_2_10_totInv_0_70=0;
    m_E_true_2_10_totPi_0_70=0;
    m_E_true_2_10_totPh_0_70=0;
    m_E_true_2_10_totK0L_0_70=0;
    m_E_true_2_10_totE_0_70=0;
    m_E_true_2_10_totMu_0_70=0;
    m_E_true_2_10_totN_0_70=0;
    m_E_true_2_10_totK_0_70=0;
    m_E_true_2_10_totP_0_70=0;
    m_E_true_2_10_totOtherCH_0_70=0;
    m_E_true_2_10_totOtherNeut_0_70=0;

    m_E_10_50_totPFO=0;
    m_E_10_50_totPi=0;
    m_E_10_50_totPh=0;
    m_E_10_50_totE=0;
    m_E_10_50_totMu=0;
    m_E_10_50_totN=0;
    m_E_10_50_totPFO_0_95=0;
    m_E_10_50_totPi_0_95=0;
    m_E_10_50_totPh_0_95=0;
    m_E_10_50_totE_0_95=0;
    m_E_10_50_totMu_0_95=0;
    m_E_10_50_totN_0_95=0;
    m_E_10_50_totPFO_0_70=0;
    m_E_10_50_totPi_0_70=0;
    m_E_10_50_totPh_0_70=0;
    m_E_10_50_totE_0_70=0;
    m_E_10_50_totMu_0_70=0;
    m_E_10_50_totN_0_70=0;

    m_E_true_10_50_totAll=0;
    m_E_true_10_50_totInv=0;
    m_E_true_10_50_totPi=0;
    m_E_true_10_50_totPh=0;
    m_E_true_10_50_totK0L=0;
    m_E_true_10_50_totE=0;
    m_E_true_10_50_totMu=0;
    m_E_true_10_50_totN=0;
    m_E_true_10_50_totK=0;
    m_E_true_10_50_totP=0;
    m_E_true_10_50_totOtherCH=0;
    m_E_true_10_50_totOtherNeut=0;
    m_E_true_10_50_totAll_0_95=0;
    m_E_true_10_50_totInv_0_95=0;
    m_E_true_10_50_totPi_0_95=0;
    m_E_true_10_50_totPh_0_95=0;
    m_E_true_10_50_totK0L_0_95=0;
    m_E_true_10_50_totE_0_95=0;
    m_E_true_10_50_totMu_0_95=0;
    m_E_true_10_50_totN_0_95=0;
    m_E_true_10_50_totK_0_95=0;
    m_E_true_10_50_totP_0_95=0;
    m_E_true_10_50_totOtherCH_0_95=0;
    m_E_true_10_50_totOtherNeut_0_95=0;
    m_E_true_10_50_totAll_0_70=0;
    m_E_true_10_50_totInv_0_70=0;
    m_E_true_10_50_totPi_0_70=0;
    m_E_true_10_50_totPh_0_70=0;
    m_E_true_10_50_totK0L_0_70=0;
    m_E_true_10_50_totE_0_70=0;
    m_E_true_10_50_totMu_0_70=0;
    m_E_true_10_50_totN_0_70=0;
    m_E_true_10_50_totK_0_70=0;
    m_E_true_10_50_totP_0_70=0;
    m_E_true_10_50_totOtherCH_0_70=0;
    m_E_true_10_50_totOtherNeut_0_70=0;

    m_E_50_totPFO=0;
    m_E_50_totPi=0;
    m_E_50_totPh=0;
    m_E_50_totE=0;
    m_E_50_totMu=0;
    m_E_50_totN=0;
    m_E_50_totPFO_0_95=0;
    m_E_50_totPi_0_95=0;
    m_E_50_totPh_0_95=0;
    m_E_50_totE_0_95=0;
    m_E_50_totMu_0_95=0;
    m_E_50_totN_0_95=0;
    m_E_50_totPFO_0_70=0;
    m_E_50_totPi_0_70=0;
    m_E_50_totPh_0_70=0;
    m_E_50_totE_0_70=0;
    m_E_50_totMu_0_70=0;
    m_E_50_totN_0_70=0;

    m_E_true_50_totAll=0;
    m_E_true_50_totInv=0;
    m_E_true_50_totPi=0;
    m_E_true_50_totPh=0;
    m_E_true_50_totK0L=0;
    m_E_true_50_totE=0;
    m_E_true_50_totMu=0;
    m_E_true_50_totN=0;
    m_E_true_50_totK=0;
    m_E_true_50_totP=0;
    m_E_true_50_totOtherCH=0;
    m_E_true_50_totOtherNeut=0;
    m_E_true_50_totAll_0_95=0;
    m_E_true_50_totInv_0_95=0;
    m_E_true_50_totPi_0_95=0;
    m_E_true_50_totPh_0_95=0;
    m_E_true_50_totK0L_0_95=0;
    m_E_true_50_totE_0_95=0;
    m_E_true_50_totMu_0_95=0;
    m_E_true_50_totN_0_95=0;
    m_E_true_50_totK_0_95=0;
    m_E_true_50_totP_0_95=0;
    m_E_true_50_totOtherCH_0_95=0;
    m_E_true_50_totOtherNeut_0_95=0;
    m_E_true_50_totAll_0_70=0;
    m_E_true_50_totInv_0_70=0;
    m_E_true_50_totPi_0_70=0;
    m_E_true_50_totPh_0_70=0;
    m_E_true_50_totK0L_0_70=0;
    m_E_true_50_totE_0_70=0;
    m_E_true_50_totMu_0_70=0;
    m_E_true_50_totN_0_70=0;
    m_E_true_50_totK_0_70=0;
    m_E_true_50_totP_0_70=0;
    m_E_true_50_totOtherCH_0_70=0;
    m_E_true_50_totOtherNeut_0_70=0;


    //multiplicities
    m_n_totPFO=0;
    m_n_totPi=0;
    m_n_totPh=0;
    m_n_totE=0;
    m_n_totMu=0;
    m_n_totN=0;
    m_n_totPFO_0_95=0;
    m_n_totPi_0_95=0;
    m_n_totPh_0_95=0;
    m_n_totE_0_95=0;
    m_n_totMu_0_95=0;
    m_n_totN_0_95=0;
    m_n_totPFO_0_70=0;
    m_n_totPi_0_70=0;
    m_n_totPh_0_70=0;
    m_n_totE_0_70=0;
    m_n_totMu_0_70=0;
    m_n_totN_0_70=0;

    m_n_true_totAll=0;
    m_n_true_totInv=0;
    m_n_true_totPi=0;
    m_n_true_totPh=0;
    m_n_true_totK0L=0;
    m_n_true_totE=0;
    m_n_true_totMu=0;
    m_n_true_totN=0;
    m_n_true_totK=0;
    m_n_true_totP=0;
    m_n_true_totOtherCH=0;
    m_n_true_totOtherNeut=0;
    m_n_true_totAll_0_95=0;
    m_n_true_totInv_0_95=0;
    m_n_true_totPi_0_95=0;
    m_n_true_totPh_0_95=0;
    m_n_true_totK0L_0_95=0;
    m_n_true_totE_0_95=0;
    m_n_true_totMu_0_95=0;
    m_n_true_totN_0_95=0;
    m_n_true_totK_0_95=0;
    m_n_true_totP_0_95=0;
    m_n_true_totOtherCH_0_95=0;
    m_n_true_totOtherNeut_0_95=0;
    m_n_true_totAll_0_70=0;
    m_n_true_totInv_0_70=0;
    m_n_true_totPi_0_70=0;
    m_n_true_totPh_0_70=0;
    m_n_true_totK0L_0_70=0;
    m_n_true_totE_0_70=0;
    m_n_true_totMu_0_70=0;
    m_n_true_totN_0_70=0;
    m_n_true_totK_0_70=0;
    m_n_true_totP_0_70=0;
    m_n_true_totOtherCH_0_70=0;
    m_n_true_totOtherNeut_0_70=0;

    m_n_0_2_totPFO=0;
    m_n_0_2_totPi=0;
    m_n_0_2_totPh=0;
    m_n_0_2_totE=0;
    m_n_0_2_totMu=0;
    m_n_0_2_totN=0;
    m_n_0_2_totPFO_0_95=0;
    m_n_0_2_totPi_0_95=0;
    m_n_0_2_totPh_0_95=0;
    m_n_0_2_totE_0_95=0;
    m_n_0_2_totMu_0_95=0;
    m_n_0_2_totN_0_95=0;
    m_n_0_2_totPFO_0_70=0;
    m_n_0_2_totPi_0_70=0;
    m_n_0_2_totPh_0_70=0;
    m_n_0_2_totE_0_70=0;
    m_n_0_2_totMu_0_70=0;
    m_n_0_2_totN_0_70=0;

    m_n_true_0_2_totAll=0;
    m_n_true_0_2_totInv=0;
    m_n_true_0_2_totPi=0;
    m_n_true_0_2_totPh=0;
    m_n_true_0_2_totE=0;
    m_n_true_0_2_totMu=0;
    m_n_true_0_2_totN=0;
    m_n_true_0_2_totK=0;
    m_n_true_0_2_totP=0;
    m_n_true_0_2_totOtherCH=0;
    m_n_true_0_2_totOtherNeut=0;
    m_n_true_0_2_totAll_0_95=0;
    m_n_true_0_2_totInv_0_95=0;
    m_n_true_0_2_totPi_0_95=0;
    m_n_true_0_2_totPh_0_95=0;
    m_n_true_0_2_totK0L_0_95=0;
    m_n_true_0_2_totE_0_95=0;
    m_n_true_0_2_totMu_0_95=0;
    m_n_true_0_2_totN_0_95=0;
    m_n_true_0_2_totK_0_95=0;
    m_n_true_0_2_totP_0_95=0;
    m_n_true_0_2_totOtherCH_0_95=0;
    m_n_true_0_2_totOtherNeut_0_95=0;
    m_n_true_0_2_totAll_0_70=0;
    m_n_true_0_2_totInv_0_70=0;
    m_n_true_0_2_totPi_0_70=0;
    m_n_true_0_2_totPh_0_70=0;
    m_n_true_0_2_totK0L_0_70=0;
    m_n_true_0_2_totE_0_70=0;
    m_n_true_0_2_totMu_0_70=0;
    m_n_true_0_2_totN_0_70=0;
    m_n_true_0_2_totK_0_70=0;
    m_n_true_0_2_totP_0_70=0;
    m_n_true_0_2_totOtherCH_0_70=0;
    m_n_true_0_2_totOtherNeut_0_70=0;

    m_n_2_10_totPFO=0;
    m_n_2_10_totPi=0;
    m_n_2_10_totPh=0;
    m_n_2_10_totE=0;
    m_n_2_10_totMu=0;
    m_n_2_10_totN=0;
    m_n_2_10_totPFO_0_95=0;
    m_n_2_10_totPi_0_95=0;
    m_n_2_10_totPh_0_95=0;
    m_n_2_10_totE_0_95=0;
    m_n_2_10_totMu_0_95=0;
    m_n_2_10_totN_0_95=0;
    m_n_2_10_totPFO_0_70=0;
    m_n_2_10_totPi_0_70=0;
    m_n_2_10_totPh_0_70=0;
    m_n_2_10_totE_0_70=0;
    m_n_2_10_totMu_0_70=0;
    m_n_2_10_totN_0_70=0;

    m_n_true_2_10_totAll=0;
    m_n_true_2_10_totInv=0;
    m_n_true_2_10_totPi=0;
    m_n_true_2_10_totPh=0;
    m_n_true_2_10_totK0L=0;
    m_n_true_2_10_totE=0;
    m_n_true_2_10_totMu=0;
    m_n_true_2_10_totN=0;
    m_n_true_2_10_totK=0;
    m_n_true_2_10_totP=0;
    m_n_true_2_10_totOtherCH=0;
    m_n_true_2_10_totOtherNeut=0;
    m_n_true_2_10_totAll_0_95=0;
    m_n_true_2_10_totInv_0_95=0;
    m_n_true_2_10_totPi_0_95=0;
    m_n_true_2_10_totPh_0_95=0;
    m_n_true_2_10_totK0L_0_95=0;
    m_n_true_2_10_totE_0_95=0;
    m_n_true_2_10_totMu_0_95=0;
    m_n_true_2_10_totN_0_95=0;
    m_n_true_2_10_totK_0_95=0;
    m_n_true_2_10_totP_0_95=0;
    m_n_true_2_10_totOtherCH_0_95=0;
    m_n_true_2_10_totOtherNeut_0_95=0;
    m_n_true_2_10_totAll_0_70=0;
    m_n_true_2_10_totInv_0_70=0;
    m_n_true_2_10_totPi_0_70=0;
    m_n_true_2_10_totPh_0_70=0;
    m_n_true_2_10_totK0L_0_70=0;
    m_n_true_2_10_totE_0_70=0;
    m_n_true_2_10_totMu_0_70=0;
    m_n_true_2_10_totN_0_70=0;
    m_n_true_2_10_totK_0_70=0;
    m_n_true_2_10_totP_0_70=0;
    m_n_true_2_10_totOtherCH_0_70=0;
    m_n_true_2_10_totOtherNeut_0_70=0;

    m_n_10_50_totPFO=0;
    m_n_10_50_totPi=0;
    m_n_10_50_totPh=0;
    m_n_10_50_totE=0;
    m_n_10_50_totMu=0;
    m_n_10_50_totN=0;
    m_n_10_50_totPFO_0_95=0;
    m_n_10_50_totPi_0_95=0;
    m_n_10_50_totPh_0_95=0;
    m_n_10_50_totE_0_95=0;
    m_n_10_50_totMu_0_95=0;
    m_n_10_50_totN_0_95=0;
    m_n_10_50_totPFO_0_70=0;
    m_n_10_50_totPi_0_70=0;
    m_n_10_50_totPh_0_70=0;
    m_n_10_50_totE_0_70=0;
    m_n_10_50_totMu_0_70=0;
    m_n_10_50_totN_0_70=0;

    m_n_true_10_50_totAll=0;
    m_n_true_10_50_totInv=0;
    m_n_true_10_50_totPi=0;
    m_n_true_10_50_totPh=0;
    m_n_true_10_50_totK0L=0;
    m_n_true_10_50_totE=0;
    m_n_true_10_50_totMu=0;
    m_n_true_10_50_totN=0;
    m_n_true_10_50_totK=0;
    m_n_true_10_50_totP=0;
    m_n_true_10_50_totOtherCH=0;
    m_n_true_10_50_totOtherNeut=0;
    m_n_true_10_50_totAll_0_95=0;
    m_n_true_10_50_totInv_0_95=0;
    m_n_true_10_50_totPi_0_95=0;
    m_n_true_10_50_totPh_0_95=0;
    m_n_true_10_50_totK0L_0_95=0;
    m_n_true_10_50_totE_0_95=0;
    m_n_true_10_50_totMu_0_95=0;
    m_n_true_10_50_totN_0_95=0;
    m_n_true_10_50_totK_0_95=0;
    m_n_true_10_50_totP_0_95=0;
    m_n_true_10_50_totOtherCH_0_95=0;
    m_n_true_10_50_totOtherNeut_0_95=0;
    m_n_true_10_50_totAll_0_70=0;
    m_n_true_10_50_totInv_0_70=0;
    m_n_true_10_50_totPi_0_70=0;
    m_n_true_10_50_totPh_0_70=0;
    m_n_true_10_50_totK0L_0_70=0;
    m_n_true_10_50_totE_0_70=0;
    m_n_true_10_50_totMu_0_70=0;
    m_n_true_10_50_totN_0_70=0;
    m_n_true_10_50_totK_0_70=0;
    m_n_true_10_50_totP_0_70=0;
    m_n_true_10_50_totOtherCH_0_70=0;
    m_n_true_10_50_totOtherNeut_0_70=0;

    m_n_50_totPFO=0;
    m_n_50_totPi=0;
    m_n_50_totPh=0;
    m_n_50_totE=0;
    m_n_50_totMu=0;
    m_n_50_totN=0;
    m_n_50_totPFO_0_95=0;
    m_n_50_totPi_0_95=0;
    m_n_50_totPh_0_95=0;
    m_n_50_totE_0_95=0;
    m_n_50_totMu_0_95=0;
    m_n_50_totN_0_95=0;
    m_n_50_totPFO_0_70=0;
    m_n_50_totPi_0_70=0;
    m_n_50_totPh_0_70=0;
    m_n_50_totE_0_70=0;
    m_n_50_totMu_0_70=0;
    m_n_50_totN_0_70=0;

    m_n_true_50_totAll=0;
    m_n_true_50_totInv=0;
    m_n_true_50_totPi=0;
    m_n_true_50_totPh=0;
    m_n_true_50_totK0L=0;
    m_n_true_50_totE=0;
    m_n_true_50_totMu=0;
    m_n_true_50_totN=0;
    m_n_true_50_totK=0;
    m_n_true_50_totP=0;
    m_n_true_50_totOtherCH=0;
    m_n_true_50_totOtherNeut=0;
    m_n_true_50_totAll_0_95=0;
    m_n_true_50_totInv_0_95=0;
    m_n_true_50_totPi_0_95=0;
    m_n_true_50_totPh_0_95=0;
    m_n_true_50_totK0L_0_95=0;
    m_n_true_50_totE_0_95=0;
    m_n_true_50_totMu_0_95=0;
    m_n_true_50_totN_0_95=0;
    m_n_true_50_totK_0_95=0;
    m_n_true_50_totP_0_95=0;
    m_n_true_50_totOtherCH_0_95=0;
    m_n_true_50_totOtherNeut_0_95=0;
    m_n_true_50_totAll_0_70=0;
    m_n_true_50_totInv_0_70=0;
    m_n_true_50_totPi_0_70=0;
    m_n_true_50_totPh_0_70=0;
    m_n_true_50_totK0L_0_70=0;
    m_n_true_50_totE_0_70=0;
    m_n_true_50_totMu_0_70=0;
    m_n_true_50_totN_0_70=0;
    m_n_true_50_totK_0_70=0;
    m_n_true_50_totP_0_70=0;
    m_n_true_50_totOtherCH_0_70=0;
    m_n_true_50_totOtherNeut_0_70=0;

    if(!m_ignoreGen){
      m_truePi0Energy->clear();
      m_truePi0_Px->clear();
      m_truePi0_Py->clear();
      m_truePi0_Pz->clear();
      m_truePi0_CosTheta->clear();
      m_truePi0_Phi->clear();
      m_truePi0_DRPh01->clear();
      m_truePi0_CosAngPh01->clear();
      m_truePi0_AngPh01->clear();
      m_truePh0Energy->clear();
      m_truePh0_Px->clear();
      m_truePh0_Py->clear();
      m_truePh0_Pz->clear();
      m_truePh0_CosTheta->clear();
      m_truePh0_Phi->clear();
      m_truePh0_DRMin->clear();
      m_truePh0_DRMinPDGID->clear();
      m_truePh0_DRMinE->clear();
      m_truePh0_CosAngMax->clear();
      m_truePh0_CosAngMaxPDGID->clear();
      m_truePh0_CosAngMaxE->clear();
      m_truePh0_CosAng0995E->clear();
      m_truePh0_DR01E->clear();
      m_truePh1Energy->clear();
      m_truePh1_Px->clear();
      m_truePh1_Py->clear();
      m_truePh1_Pz->clear();
      m_truePh1_CosTheta->clear();
      m_truePh1_Phi->clear();
      m_truePh1_DRMin->clear();
      m_truePh1_DRMinPDGID->clear();
      m_truePh1_DRMinE->clear();
      m_truePh1_CosAngMax->clear();
      m_truePh1_CosAngMaxPDGID->clear();
      m_truePh1_CosAngMaxE->clear();
      m_truePh1_CosAng0995E->clear();
      m_truePh1_DR01E->clear();
    }
    m_recoPhEnergy->clear();
    m_recoPh_Px->clear();
    m_recoPh_Py->clear();
    m_recoPh_Pz->clear();
    m_recoPh_CosTheta->clear();
    m_recoPh_Theta->clear();
    m_recoPh_Phi->clear();
    m_recoPh_PDGID->clear();
    m_recoPh_E_EB->clear();
    m_recoPh_E_EE->clear();
    m_recoPh_E_EO->clear();
    m_recoPh_E_HB->clear();
    m_recoPh_E_HE->clear();
    m_recoPh_E_HO->clear();
    m_recoPh_firstLayerECAL->clear();
    m_recoPh_lastLayerECAL->clear();
    m_recoPh_nhitsEB->clear();
    m_recoPh_nhitsEE->clear();
    m_recoPh_nhitsEO->clear();
    m_recoPh_firstLayerHCAL->clear();
    m_recoPh_lastLayerHCAL->clear();
    m_recoPh_nhitsHB->clear();
    m_recoPh_nhitsHE->clear();
    m_recoPh_nhitsHO->clear();
    m_recoPh_DRMin_Ph->clear();
    m_recoPh_DRMin_E_Ph->clear();
    m_recoPh_CosAngMax_Ph->clear();
    m_recoPh_CosAngMax_E_Ph->clear();
    m_recoPh_DRMin->clear();
    m_recoPh_DRMin_PDGID->clear();
    m_recoPh_DRMin_E->clear();
    m_recoPh_CosAngMax->clear();
    m_recoPh_CosAngMax_PDGID->clear();
    m_recoPh_CosAngMax_E->clear();
    m_recoPh_CosAng0995_E->clear();
    m_recoPh_DR01_E->clear();
    m_recoPi0Cand_CosAngMax_M->clear();
    m_recoPi0Cand_DRMin_M->clear();

    m_jet_excl_Px->clear(); 
    m_jet_excl_Py->clear(); 
    m_jet_excl_Pz->clear(); 
    m_jet_excl_E->clear(); 
    m_jet_excl_Phi->clear(); 
    m_jet_excl_CosTheta->clear(); 
    m_jet_excl_piE->clear();
    m_jet_excl_phE->clear();
    m_jet_excl_elE->clear();
    m_jet_excl_muE->clear();
    m_jet_excl_nE->clear();
    m_jet_excl_neutElseE->clear();
    m_jet_excl_chElseE->clear();
    m_jet_excl_neutMult->clear();
    m_jet_excl_chMult->clear();

    m_jet_incl_Px->clear(); 
    m_jet_incl_Py->clear(); 
    m_jet_incl_Pz->clear(); 
    m_jet_incl_E->clear(); 
    m_jet_incl_Phi->clear(); 
    m_jet_incl_CosTheta->clear(); 
    m_jet_incl_piE->clear();
    m_jet_incl_phE->clear();
    m_jet_incl_elE->clear();
    m_jet_incl_muE->clear();
    m_jet_incl_nE->clear();
    m_jet_incl_neutElseE->clear();
    m_jet_incl_chElseE->clear();
    m_jet_incl_neutMult->clear();
    m_jet_incl_chMult->clear();

    if(m_fillMEInfo){
      m_trueME_E->clear();
      m_trueME_Px->clear();
      m_trueME_Py->clear();
      m_trueME_Pz->clear();
      m_trueME_PDGID->clear();
    }
    //std::cout<<"before mcColl"<<std::endl;
    LCCollection * mcColl =0;
    getCollection(mcColl,m_inputMCParticleCollection,evt);
    //std::cout<<"after mcColl"<<std::endl;
    if(mcColl!=NULL){
      //std::cout<<"in mcColl"<<std::endl;
    //Look for a photon
      int index_first_Z_daughter_two_neut=-1;
      bool found_first_hadZ=false;
    for(int m =0; m< mcColl->getNumberOfElements(); m++){
        MCParticle* mcp= dynamic_cast<MCParticle*>( mcColl->getElementAt(m) ) ;
	//fill first quark anti quark pair which appears in the history
	if(abs(mcp->getPDG())<7 && m_d1_mcE<0) {
	  //std::cout<<"get to daughter "<<i<<" "<<mcp->getPDG()<<"/"<<mcp->getEnergy()<<std::endl;
	  m_d1_mcPDGID=mcp->getPDG();
	  m_d1_mcE=mcp->getEnergy();
	  m_d1_mcPx=mcp->getMomentum()[0];
	  m_d1_mcPy=mcp->getMomentum()[1];
	  m_d1_mcPz=mcp->getMomentum()[2];
	  m_d1_mcMass=mcp->getMass();
	  m_d1_mcPhi=atan2(m_d1_mcPy,m_d1_mcPx);
	  m_d1_mcCosTheta= m_d1_mcPz/sqrt(m_d1_mcPx*m_d1_mcPx+m_d1_mcPy*m_d1_mcPy+m_d1_mcPz*m_d1_mcPz);
	  m_d1_mcTheta=acos(m_d1_mcCosTheta);
	}
	if(m_d2_mcE<0 && (abs(mcp->getPDG())<7 && mcp->getPDG()==(-m_d1_mcPDGID))){
	  //std::cout<<"get to daughter "<<i<<" "<<mcp->getPDG()<<"/"<<mcp->getEnergy()<<std::endl;
	  m_d2_mcPDGID=mcp->getPDG();
	  m_d2_mcE=mcp->getEnergy();
	  m_d2_mcPx=mcp->getMomentum()[0];
	  m_d2_mcPy=mcp->getMomentum()[1];
	  m_d2_mcPz=mcp->getMomentum()[2];
	  m_d2_mcMass=mcp->getMass();
	  m_d2_mcPhi=atan2(m_d2_mcPy,m_d2_mcPx);
	  m_d2_mcCosTheta= m_d2_mcPz/sqrt(m_d2_mcPx*m_d2_mcPx+m_d2_mcPy*m_d2_mcPy+m_d2_mcPz*m_d2_mcPz);
	  m_d2_mcTheta=acos(m_d2_mcCosTheta);
	}
	if((mcp->getDaughters().size ()>1)){
	  unsigned int_neut12=0;
	  unsigned int_neut14=0;
	  unsigned int_neut16=0;
	  for(unsigned int i=0;i<mcp->getDaughters().size ();i++){
	    if(index_first_Z_daughter_two_neut==-1){//neutrinos not found yet
	      if(abs(mcp->getDaughters()[i]->getPDG())==12){
		int_neut12+=1;
	      }else if(abs(mcp->getDaughters()[i]->getPDG())==14){
		int_neut14+=1;
	      }else if(abs(mcp->getDaughters()[i]->getPDG())==16){
		int_neut16+=1;
	      }
	    }
	  }
          if(int_neut12==2 && int_neut14==0 && int_neut16==0){
            index_first_Z_daughter_two_neut=m;
          }else if(int_neut12==0 && int_neut14==2 && int_neut16==0){
            index_first_Z_daughter_two_neut=m;
          }else if (int_neut12==0 && int_neut14==0 && int_neut16==2){
            index_first_Z_daughter_two_neut=m;
          }
          if(index_first_Z_daughter_two_neut!=-1 && m_E_trueNeut==0){
            for(unsigned int i=0;i<mcp->getDaughters().size ();i++){
              if (abs(mcp->getDaughters()[i]->getPDG())==12 || abs(mcp->getDaughters()[i]->getPDG())==14 ||  abs(mcp->getDaughters()[i]->getPDG())==16){
		if( m_fillMEInfo){
		  m_trueME_Px->push_back(mcp->getDaughters()[i]->getMomentum()[0]);
		  m_trueME_Py->push_back(mcp->getDaughters()[i]->getMomentum()[1]);
		  m_trueME_Pz->push_back(mcp->getDaughters()[i]->getMomentum()[2]);
		  m_trueME_E->push_back(mcp->getDaughters()[i]->getEnergy());
		  m_trueME_PDGID->push_back(mcp->getDaughters()[i]->getPDG());
		}
                m_E_trueNeut+=mcp->getDaughters()[i]->getEnergy();
                m_px_trueNeut+=mcp->getDaughters()[i]->getMomentum()[0];
                m_py_trueNeut+=mcp->getDaughters()[i]->getMomentum()[1];
                m_pz_trueNeut+=mcp->getDaughters()[i]->getMomentum()[2];
              }
            }
	  }     
	  if(!found_first_hadZ){
	    for(unsigned int i=0;i<mcp->getDaughters().size ();i++){
              if (abs(mcp->getDaughters()[i]->getPDG())<7){
		found_first_hadZ=true;
		if( m_fillMEInfo){
		  m_trueME_Px->push_back(mcp->getDaughters()[i]->getMomentum()[0]);
		  m_trueME_Py->push_back(mcp->getDaughters()[i]->getMomentum()[1]);
		  m_trueME_Pz->push_back(mcp->getDaughters()[i]->getMomentum()[2]);
		  m_trueME_E->push_back(mcp->getDaughters()[i]->getEnergy());
		  m_trueME_PDGID->push_back(mcp->getDaughters()[i]->getPDG());
		}
                m_E_trueZ2+=mcp->getDaughters()[i]->getEnergy();
                m_px_trueZ2+=mcp->getDaughters()[i]->getMomentum()[0];
                m_py_trueZ2+=mcp->getDaughters()[i]->getMomentum()[1];
                m_pz_trueZ2+=mcp->getDaughters()[i]->getMomentum()[2];
              }
            }
	  }
	}
	if(mcp->getPDG() == 23 && (mcp->getDaughters().size ()==2)){
	  m_Z_mcE=mcp->getEnergy();
	  m_Z_mcNDaughter=mcp->getDaughters().size ();
	  //std::cout<<"get to Z here, two daughters "<<mcp->getDaughters().size ()<<std::endl;
	  //for(unsigned int i=0;i<mcp->getDaughters().size ();i++){
	    /*
	    */
	  //}
	  //std::cout<<"after loop over daughter "<<m_d2_mcE+m_d1_mcE<<"/"<<m_Z_mcE<<std::endl;
	}
	if(m_ignoreGen){
	  continue;
	}
	if(mcp->getGeneratorStatus()==1){//visible sum of stable particles --> take neutrinos out
	  double cos_theta=mcp->getMomentum()[2]/sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]+mcp->getMomentum()[2]*mcp->getMomentum()[2]);
	  if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
	    m_E_true_totAll+=mcp->getEnergy();
	    m_px_true_totAll+=mcp->getMomentum()[0];
	    m_py_true_totAll+=mcp->getMomentum()[1];
	    m_pz_true_totAll+=mcp->getMomentum()[2];
	    m_n_true_totAll+=1;
	    if(abs(mcp->getPDG())==211){
	      m_E_true_totPi+=mcp->getEnergy();
	      m_px_true_totPi+=mcp->getMomentum()[0];
	      m_py_true_totPi+=mcp->getMomentum()[1];
	      m_pz_true_totPi+=mcp->getMomentum()[2];
	      m_n_true_totPi+=1;
	    }else if (mcp->getPDG()==22){
	      m_E_true_totPh+=mcp->getEnergy();
	      m_px_true_totPh+=mcp->getMomentum()[0];
	      m_py_true_totPh+=mcp->getMomentum()[1];
	      m_pz_true_totPh+=mcp->getMomentum()[2];
	      m_n_true_totPh+=1;
	    }else if (abs(mcp->getPDG())==130){
	      m_E_true_totK0L+=mcp->getEnergy();
	      m_px_true_totK0L+=mcp->getMomentum()[0];
	      m_py_true_totK0L+=mcp->getMomentum()[1];
	      m_pz_true_totK0L+=mcp->getMomentum()[2];
	      m_n_true_totK0L+=1;
	    }else if(abs(mcp->getPDG())==11){
	      m_E_true_totE+=mcp->getEnergy();
	      m_px_true_totE+=mcp->getMomentum()[0];
	      m_py_true_totE+=mcp->getMomentum()[1];
	      m_pz_true_totE+=mcp->getMomentum()[2];
	      m_n_true_totE+=1;
	    }else if(abs(mcp->getPDG())==13 ){
	      m_E_true_totMu+=mcp->getEnergy();
	      m_px_true_totMu+=mcp->getMomentum()[0];
	      m_py_true_totMu+=mcp->getMomentum()[1];
	      m_pz_true_totMu+=mcp->getMomentum()[2];
	      m_n_true_totMu+=1;
	    }else if(abs(mcp->getPDG())==2112 ){
	      m_E_true_totN+=mcp->getEnergy();
	      m_px_true_totN+=mcp->getMomentum()[0];
	      m_py_true_totN+=mcp->getMomentum()[1];
	      m_pz_true_totN+=mcp->getMomentum()[2];
	      m_n_true_totN+=1;
	    }else if(abs(mcp->getPDG())==2212 ){
	      m_E_true_totP+=mcp->getEnergy();
	      m_px_true_totP+=mcp->getMomentum()[0];
	      m_py_true_totP+=mcp->getMomentum()[1];
	      m_pz_true_totP+=mcp->getMomentum()[2];
	      m_n_true_totP+=1;
	    }else if(abs(mcp->getPDG())==321 ){
	      m_E_true_totK+=mcp->getEnergy();
	      m_px_true_totK+=mcp->getMomentum()[0];
	      m_py_true_totK+=mcp->getMomentum()[1];
	      m_pz_true_totK+=mcp->getMomentum()[2];
	      m_n_true_totK+=1;
	    }else if (mcp->getCharge()==0){
	      m_E_true_totOtherNeut+=mcp->getEnergy();
	      m_n_true_totOtherNeut+=1;
	      std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
	    }else{
	      m_E_true_totOtherCH+=mcp->getEnergy();
	      m_n_true_totOtherCH+=1;
	      std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
	    }
	  }else{
	    m_E_true_totInv+=mcp->getEnergy();
	    m_n_true_totInv+=mcp->getEnergy();
	    m_px_true_totInv+=mcp->getMomentum()[0];
	    m_py_true_totInv+=mcp->getMomentum()[1];
	    m_pz_true_totInv+=mcp->getMomentum()[2];
	  }
	  if(mcp->getEnergy()<2.0){
	    if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
	      m_E_true_0_2_totAll+=mcp->getEnergy();
	      m_n_true_0_2_totAll+=1;
	      if(abs(mcp->getPDG())==211){
		m_E_true_0_2_totPi+=mcp->getEnergy();
		m_n_true_0_2_totPi+=1;
	      }else if (mcp->getPDG()==22){
		m_E_true_0_2_totPh+=mcp->getEnergy();
		m_n_true_0_2_totPh+=1;
	      }else if (abs(mcp->getPDG())==130){
		m_E_true_0_2_totK0L+=mcp->getEnergy();
		m_n_true_0_2_totK0L+=1;
	      }else if(abs(mcp->getPDG())==11){
		m_E_true_0_2_totE+=mcp->getEnergy();
		m_n_true_0_2_totE+=1;
	      }else if(abs(mcp->getPDG())==13 ){
		m_E_true_0_2_totMu+=mcp->getEnergy();
		m_n_true_0_2_totMu+=1;
	      }else if(abs(mcp->getPDG())==2112 ){
		m_E_true_0_2_totN+=mcp->getEnergy();
		m_n_true_0_2_totN+=1;
	      }else if(abs(mcp->getPDG())==2212 ){
		m_E_true_0_2_totP+=mcp->getEnergy();
		m_n_true_0_2_totP+=1;
	      }else if(abs(mcp->getPDG())==321 ){
		m_E_true_0_2_totK+=mcp->getEnergy();
		m_n_true_0_2_totK+=1;
	      }else if (mcp->getCharge()==0){
		m_E_true_0_2_totOtherNeut+=mcp->getEnergy();
		m_n_true_0_2_totOtherNeut+=1;
		std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
	      }else{
		m_E_true_0_2_totOtherCH+=mcp->getEnergy();
		m_n_true_0_2_totOtherCH+=1;
		std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
	      }
	    }else{
	      m_E_true_0_2_totInv+=mcp->getEnergy();
	      m_n_true_0_2_totInv+=mcp->getEnergy();
	    }
	  }else if(mcp->getEnergy()<10.0){
	    if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
	      m_E_true_2_10_totAll+=mcp->getEnergy();
	      m_n_true_2_10_totAll+=1;
	      if(abs(mcp->getPDG())==211){
		m_E_true_2_10_totPi+=mcp->getEnergy();
		m_n_true_2_10_totPi+=1;
	      }else if (mcp->getPDG()==22){
		m_E_true_2_10_totPh+=mcp->getEnergy();
		m_n_true_2_10_totPh+=1;
	      }else if (abs(mcp->getPDG())==130){
		m_E_true_2_10_totK0L+=mcp->getEnergy();
		m_n_true_2_10_totK0L+=1;
	      }else if(abs(mcp->getPDG())==11){
		m_E_true_2_10_totE+=mcp->getEnergy();
		m_n_true_2_10_totE+=1;
	      }else if(abs(mcp->getPDG())==13 ){
		m_E_true_2_10_totMu+=mcp->getEnergy();
		m_n_true_2_10_totMu+=1;
	      }else if(abs(mcp->getPDG())==2112 ){
		m_E_true_2_10_totN+=mcp->getEnergy();
		m_n_true_2_10_totN+=1;
	      }else if(abs(mcp->getPDG())==2212 ){
		m_E_true_2_10_totP+=mcp->getEnergy();
		m_n_true_2_10_totP+=1;
	      }else if(abs(mcp->getPDG())==321 ){
		m_E_true_2_10_totK+=mcp->getEnergy();
		m_n_true_2_10_totK+=1;
	      }else if (mcp->getCharge()==0){
		m_E_true_2_10_totOtherNeut+=mcp->getEnergy();
		m_n_true_2_10_totOtherNeut+=1;
		std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
	      }else{
		m_E_true_2_10_totOtherCH+=mcp->getEnergy();
		m_n_true_2_10_totOtherCH+=1;
		std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
	      }
	    }else{
	      m_E_true_2_10_totInv+=mcp->getEnergy();
	      m_n_true_2_10_totInv+=mcp->getEnergy();
	    }
	  }else if(mcp->getEnergy()<50.0){
	    if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
	      m_E_true_10_50_totAll+=mcp->getEnergy();
	      m_n_true_10_50_totAll+=1;
	      if(abs(mcp->getPDG())==211){
		m_E_true_10_50_totPi+=mcp->getEnergy();
		m_n_true_10_50_totPi+=1;
	      }else if (mcp->getPDG()==22){
		m_E_true_10_50_totPh+=mcp->getEnergy();
		m_n_true_10_50_totPh+=1;
	      }else if (abs(mcp->getPDG())==130){
		m_E_true_10_50_totK0L+=mcp->getEnergy();
		m_n_true_10_50_totK0L+=1;
	      }else if(abs(mcp->getPDG())==11){
		m_E_true_10_50_totE+=mcp->getEnergy();
		m_n_true_10_50_totE+=1;
	      }else if(abs(mcp->getPDG())==13 ){
		m_E_true_10_50_totMu+=mcp->getEnergy();
		m_n_true_10_50_totMu+=1;
	      }else if(abs(mcp->getPDG())==2112 ){
		m_E_true_10_50_totN+=mcp->getEnergy();
		m_n_true_10_50_totN+=1;
	      }else if(abs(mcp->getPDG())==2212 ){
		m_E_true_10_50_totP+=mcp->getEnergy();
		m_n_true_10_50_totP+=1;
	      }else if(abs(mcp->getPDG())==321 ){
		m_E_true_10_50_totK+=mcp->getEnergy();
		m_n_true_10_50_totK+=1;
	      }else if (mcp->getCharge()==0){
		m_E_true_10_50_totOtherNeut+=mcp->getEnergy();
		m_n_true_10_50_totOtherNeut+=1;
		std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
	      }else{
		m_E_true_10_50_totOtherCH+=mcp->getEnergy();
		m_n_true_10_50_totOtherCH+=1;
		std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
	      }
	    }else{
	      m_E_true_10_50_totInv+=mcp->getEnergy();
	      m_n_true_10_50_totInv+=mcp->getEnergy();
	    }
	  }else{
	    if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
	      m_E_true_50_totAll+=mcp->getEnergy();
	      m_n_true_50_totAll+=1;
	      if(abs(mcp->getPDG())==211){
		m_E_true_50_totPi+=mcp->getEnergy();
		m_n_true_50_totPi+=1;
	      }else if (mcp->getPDG()==22){
		m_E_true_50_totPh+=mcp->getEnergy();
		m_n_true_50_totPh+=1;
	      }else if (abs(mcp->getPDG())==130){
		m_E_true_50_totK0L+=mcp->getEnergy();
		m_n_true_50_totK0L+=1;
	      }else if(abs(mcp->getPDG())==11){
		m_E_true_50_totE+=mcp->getEnergy();
		m_n_true_50_totE+=1;
	      }else if(abs(mcp->getPDG())==13 ){
		m_E_true_50_totMu+=mcp->getEnergy();
		m_n_true_50_totMu+=1;
	      }else if(abs(mcp->getPDG())==2112 ){
		m_E_true_50_totN+=mcp->getEnergy();
		m_n_true_50_totN+=1;
	      }else if(abs(mcp->getPDG())==2212 ){
		m_E_true_50_totP+=mcp->getEnergy();
		m_n_true_50_totP+=1;
	      }else if(abs(mcp->getPDG())==321 ){
		m_E_true_50_totK+=mcp->getEnergy();
		m_n_true_50_totK+=1;
	      }else if (mcp->getCharge()==0){
		m_E_true_50_totOtherNeut+=mcp->getEnergy();
		m_n_true_50_totOtherNeut+=1;
		std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
	      }else{
		m_E_true_50_totOtherCH+=mcp->getEnergy();
		m_n_true_50_totOtherCH+=1;
		std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
	      }
	    }else{
	      m_E_true_50_totInv+=mcp->getEnergy();
	      m_n_true_50_totInv+=mcp->getEnergy();
	    }
	  }//general true values
	  if(fabs(cos_theta)<0.95){
	    if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
	      m_E_true_totAll_0_95+=mcp->getEnergy();
	      m_n_true_totAll_0_95+=1;
	      if(abs(mcp->getPDG())==211){
		m_E_true_totPi_0_95+=mcp->getEnergy();
		m_n_true_totPi_0_95+=1;
	      }else if (mcp->getPDG()==22){
		m_E_true_totPh_0_95+=mcp->getEnergy();
		m_n_true_totPh_0_95+=1;
	      }else if (abs(mcp->getPDG())==130){
		m_E_true_totK0L_0_95+=mcp->getEnergy();
		m_n_true_totK0L_0_95+=1;
	      }else if(abs(mcp->getPDG())==11){
		m_E_true_totE_0_95+=mcp->getEnergy();
		m_n_true_totE_0_95+=1;
	      }else if(abs(mcp->getPDG())==13 ){
		m_E_true_totMu_0_95+=mcp->getEnergy();
		m_n_true_totMu_0_95+=1;
	      }else if(abs(mcp->getPDG())==2112 ){
		m_E_true_totN_0_95+=mcp->getEnergy();
		m_n_true_totN_0_95+=1;
	      }else if(abs(mcp->getPDG())==2212 ){
		m_E_true_totP_0_95+=mcp->getEnergy();
		m_n_true_totP_0_95+=1;
	      }else if(abs(mcp->getPDG())==321 ){
		m_E_true_totK_0_95+=mcp->getEnergy();
		m_n_true_totK_0_95+=1;
	      }else if (mcp->getCharge()==0){
		m_E_true_totOtherNeut_0_95+=mcp->getEnergy();
		m_n_true_totOtherNeut_0_95+=1;
		std::cout<<"other type neutral PDG 95"<<mcp->getPDG()<<std::endl;
	      }else{
		m_E_true_totOtherCH_0_95+=mcp->getEnergy();
		m_n_true_totOtherCH_0_95+=1;
		std::cout<<"other type charged PDG 95"<<mcp->getPDG()<<std::endl;
	      }
	    }
	    if(mcp->getEnergy()<2.0){
	      if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
		m_E_true_0_2_totAll_0_95+=mcp->getEnergy();
		m_n_true_0_2_totAll_0_95+=1;
		if(abs(mcp->getPDG())==211){
		  m_E_true_0_2_totPi_0_95+=mcp->getEnergy();
		  m_n_true_0_2_totPi_0_95+=1;
		}else if (mcp->getPDG()==22){
		  m_E_true_0_2_totPh_0_95+=mcp->getEnergy();
		  m_n_true_0_2_totPh_0_95+=1;
		}else if (abs(mcp->getPDG())==130){
		  m_E_true_0_2_totK0L_0_95+=mcp->getEnergy();
		  m_n_true_0_2_totK0L_0_95+=1;
		}else if(abs(mcp->getPDG())==11){
		  m_E_true_0_2_totE_0_95+=mcp->getEnergy();
		  m_n_true_0_2_totE_0_95+=1;
		}else if(abs(mcp->getPDG())==13 ){
		  m_E_true_0_2_totMu_0_95+=mcp->getEnergy();
		  m_n_true_0_2_totMu_0_95+=1;
		}else if(abs(mcp->getPDG())==2112 ){
		  m_E_true_0_2_totN_0_95+=mcp->getEnergy();
		  m_n_true_0_2_totN_0_95+=1;
		}else if(abs(mcp->getPDG())==2212 ){
		  m_E_true_0_2_totP_0_95+=mcp->getEnergy();
		  m_n_true_0_2_totP_0_95+=1;
		}else if(abs(mcp->getPDG())==321 ){
		  m_E_true_0_2_totK_0_95+=mcp->getEnergy();
		  m_n_true_0_2_totK_0_95+=1;
		}else if (mcp->getCharge()==0){
		  m_E_true_0_2_totOtherNeut_0_95+=mcp->getEnergy();
		  m_n_true_0_2_totOtherNeut_0_95+=1;
		  std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
		}else{
		  m_E_true_0_2_totOtherCH_0_95+=mcp->getEnergy();
		  m_n_true_0_2_totOtherCH_0_95+=1;
		  std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
		}
	      }else{
		m_E_true_0_2_totInv_0_95+=mcp->getEnergy();
		m_n_true_0_2_totInv_0_95+=mcp->getEnergy();
	      }
	    }else if(mcp->getEnergy()<10.0){
	      if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
		m_E_true_2_10_totAll_0_95+=mcp->getEnergy();
		m_n_true_2_10_totAll_0_95+=1;
		if(abs(mcp->getPDG())==211){
		  m_E_true_2_10_totPi_0_95+=mcp->getEnergy();
		  m_n_true_2_10_totPi_0_95+=1;
		}else if (mcp->getPDG()==22){
		  m_E_true_2_10_totPh_0_95+=mcp->getEnergy();
		  m_n_true_2_10_totPh_0_95+=1;
		}else if (abs(mcp->getPDG())==130){
		  m_E_true_2_10_totK0L_0_95+=mcp->getEnergy();
		  m_n_true_2_10_totK0L_0_95+=1;
		}else if(abs(mcp->getPDG())==11){
		  m_E_true_2_10_totE_0_95+=mcp->getEnergy();
		  m_n_true_2_10_totE_0_95+=1;
		}else if(abs(mcp->getPDG())==13 ){
		  m_E_true_2_10_totMu_0_95+=mcp->getEnergy();
		  m_n_true_2_10_totMu_0_95+=1;
		}else if(abs(mcp->getPDG())==2112 ){
		  m_E_true_2_10_totN_0_95+=mcp->getEnergy();
		  m_n_true_2_10_totN_0_95+=1;
		}else if(abs(mcp->getPDG())==2212 ){
		  m_E_true_2_10_totP_0_95+=mcp->getEnergy();
		  m_n_true_2_10_totP_0_95+=1;
		}else if(abs(mcp->getPDG())==321 ){
		  m_E_true_2_10_totK_0_95+=mcp->getEnergy();
		  m_n_true_2_10_totK_0_95+=1;
		}else if (mcp->getCharge()==0){
		  m_E_true_2_10_totOtherNeut_0_95+=mcp->getEnergy();
		  m_n_true_2_10_totOtherNeut_0_95+=1;
		  std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
		}else{
		m_E_true_2_10_totOtherCH_0_95+=mcp->getEnergy();
		m_n_true_2_10_totOtherCH_0_95+=1;
		std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
		}
	      }else{
		m_E_true_2_10_totInv_0_95+=mcp->getEnergy();
		m_n_true_2_10_totInv_0_95+=mcp->getEnergy();
	      }
	    }else if(mcp->getEnergy()<50.0){
	      if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
		m_E_true_10_50_totAll_0_95+=mcp->getEnergy();
		m_n_true_10_50_totAll_0_95+=1;
		if(abs(mcp->getPDG())==211){
		  m_E_true_10_50_totPi_0_95+=mcp->getEnergy();
		  m_n_true_10_50_totPi_0_95+=1;
		}else if (mcp->getPDG()==22){
		  m_E_true_10_50_totPh_0_95+=mcp->getEnergy();
		  m_n_true_10_50_totPh_0_95+=1;
		}else if (abs(mcp->getPDG())==130){
		  m_E_true_10_50_totK0L_0_95+=mcp->getEnergy();
		  m_n_true_10_50_totK0L_0_95+=1;
		}else if(abs(mcp->getPDG())==11){
		  m_E_true_10_50_totE_0_95+=mcp->getEnergy();
		  m_n_true_10_50_totE_0_95+=1;
		}else if(abs(mcp->getPDG())==13 ){
		  m_E_true_10_50_totMu_0_95+=mcp->getEnergy();
		  m_n_true_10_50_totMu_0_95+=1;
		}else if(abs(mcp->getPDG())==2112 ){
		  m_E_true_10_50_totN_0_95+=mcp->getEnergy();
		  m_n_true_10_50_totN_0_95+=1;
		}else if(abs(mcp->getPDG())==2212 ){
		  m_E_true_10_50_totP_0_95+=mcp->getEnergy();
		  m_n_true_10_50_totP_0_95+=1;
		}else if(abs(mcp->getPDG())==321 ){
		  m_E_true_10_50_totK_0_95+=mcp->getEnergy();
		  m_n_true_10_50_totK_0_95+=1;
		}else if (mcp->getCharge()==0){
		  m_E_true_10_50_totOtherNeut_0_95+=mcp->getEnergy();
		  m_n_true_10_50_totOtherNeut_0_95+=1;
		  std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
		}else{
		  m_E_true_10_50_totOtherCH_0_95+=mcp->getEnergy();
		  m_n_true_10_50_totOtherCH_0_95+=1;
		  std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
		}
	      }else{
		m_E_true_10_50_totInv_0_95+=mcp->getEnergy();
		m_n_true_10_50_totInv_0_95+=mcp->getEnergy();
	      }
	    }else{
	      if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
		m_E_true_50_totAll_0_95+=mcp->getEnergy();
		m_n_true_50_totAll_0_95+=1;
		if(abs(mcp->getPDG())==211){
		  m_E_true_50_totPi_0_95+=mcp->getEnergy();
		  m_n_true_50_totPi_0_95+=1;
		}else if (mcp->getPDG()==22){
		  m_E_true_50_totPh_0_95+=mcp->getEnergy();
		  m_n_true_50_totPh_0_95+=1;
		}else if (abs(mcp->getPDG())==130){
		  m_E_true_50_totK0L_0_95+=mcp->getEnergy();
		  m_n_true_50_totK0L_0_95+=1;
		}else if(abs(mcp->getPDG())==11){
		  m_E_true_50_totE_0_95+=mcp->getEnergy();
		  m_n_true_50_totE_0_95+=1;
		}else if(abs(mcp->getPDG())==13 ){
		  m_E_true_50_totMu_0_95+=mcp->getEnergy();
		  m_n_true_50_totMu_0_95+=1;
		}else if(abs(mcp->getPDG())==2112 ){
		  m_E_true_50_totN_0_95+=mcp->getEnergy();
		  m_n_true_50_totN_0_95+=1;
		}else if(abs(mcp->getPDG())==2212 ){
		  m_E_true_50_totP_0_95+=mcp->getEnergy();
		  m_n_true_50_totP_0_95+=1;
		}else if(abs(mcp->getPDG())==321 ){
		  m_E_true_50_totK_0_95+=mcp->getEnergy();
		  m_n_true_50_totK_0_95+=1;
		}else if (mcp->getCharge()==0){
		  m_E_true_50_totOtherNeut_0_95+=mcp->getEnergy();
		  m_n_true_50_totOtherNeut_0_95+=1;
		  std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
		}else{
		  m_E_true_50_totOtherCH_0_95+=mcp->getEnergy();
		  m_n_true_50_totOtherCH_0_95+=1;
		  std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
		}
	      }else{
		m_E_true_50_totInv_0_95+=mcp->getEnergy();
		m_n_true_50_totInv_0_95+=mcp->getEnergy();
	      }
	    }//0.95 true values for various energy bins
	    if(fabs(cos_theta)<0.70){
	      if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
		m_E_true_totAll_0_70+=mcp->getEnergy();
		m_n_true_totAll_0_70+=1;
		if(abs(mcp->getPDG())==211){
		  m_E_true_totPi_0_70+=mcp->getEnergy();
		  m_n_true_totPi_0_70+=1;
		}else if (mcp->getPDG()==22){
		  m_E_true_totPh_0_70+=mcp->getEnergy();
		  m_n_true_totPh_0_70+=1;
		}else if (abs(mcp->getPDG())==130){
		  m_E_true_totK0L_0_70+=mcp->getEnergy();
		  m_n_true_totK0L_0_70+=1;
		}else if(abs(mcp->getPDG())==11){
		  m_E_true_totE_0_70+=mcp->getEnergy();
		  m_n_true_totE_0_70+=1;
		}else if(abs(mcp->getPDG())==13 ){
		  m_E_true_totMu_0_70+=mcp->getEnergy();
		  m_n_true_totMu_0_70+=1;
		}else if(abs(mcp->getPDG())==2112 ){
		  m_E_true_totN_0_70+=mcp->getEnergy();
		  m_n_true_totN_0_70+=1;
		}else if(abs(mcp->getPDG())==2212 ){
		  m_E_true_totP_0_70+=mcp->getEnergy();
		  m_n_true_totP_0_70+=1;
		}else if(abs(mcp->getPDG())==321 ){
		  m_E_true_totK_0_70+=mcp->getEnergy();
		  m_n_true_totK_0_70+=1;
		}else if (mcp->getCharge()==0){
		  m_E_true_totOtherNeut_0_70+=mcp->getEnergy();
		  m_n_true_totOtherNeut_0_70+=1;
		  std::cout<<"other type neutral PDG 070"<<mcp->getPDG()<<std::endl;
		}else{
		  m_E_true_totOtherCH_0_70+=mcp->getEnergy();
		  m_n_true_totOtherCH_0_70+=1;
		  std::cout<<"other type charged PDG 070"<<mcp->getPDG()<<std::endl;
		}
	      }
	      if(mcp->getEnergy()<2.0){
		if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
		  m_E_true_0_2_totAll_0_70+=mcp->getEnergy();
		  m_n_true_0_2_totAll_0_70+=1;
		  if(abs(mcp->getPDG())==211){
		    m_E_true_0_2_totPi_0_70+=mcp->getEnergy();
		    m_n_true_0_2_totPi_0_70+=1;
		  }else if (mcp->getPDG()==22){
		    m_E_true_0_2_totPh_0_70+=mcp->getEnergy();
		    m_n_true_0_2_totPh_0_70+=1;
		  }else if (abs(mcp->getPDG())==130){
		    m_E_true_0_2_totK0L_0_70+=mcp->getEnergy();
		    m_n_true_0_2_totK0L_0_70+=1;
		  }else if(abs(mcp->getPDG())==11){
		    m_E_true_0_2_totE_0_70+=mcp->getEnergy();
		    m_n_true_0_2_totE_0_70+=1;
		  }else if(abs(mcp->getPDG())==13 ){
		    m_E_true_0_2_totMu_0_70+=mcp->getEnergy();
		    m_n_true_0_2_totMu_0_70+=1;
		  }else if(abs(mcp->getPDG())==2112 ){
		    m_E_true_0_2_totN_0_70+=mcp->getEnergy();
		    m_n_true_0_2_totN_0_70+=1;
		  }else if(abs(mcp->getPDG())==2212 ){
		    m_E_true_0_2_totP_0_70+=mcp->getEnergy();
		    m_n_true_0_2_totP_0_70+=1;
		  }else if(abs(mcp->getPDG())==321 ){
		    m_E_true_0_2_totK_0_70+=mcp->getEnergy();
		    m_n_true_0_2_totK_0_70+=1;
		  }else if (mcp->getCharge()==0){
		    m_E_true_0_2_totOtherNeut_0_70+=mcp->getEnergy();
		    m_n_true_0_2_totOtherNeut_0_70+=1;
		    std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
		  }else{
		    m_E_true_0_2_totOtherCH_0_70+=mcp->getEnergy();
		    m_n_true_0_2_totOtherCH_0_70+=1;
		    std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
		  }
		}else{
		  m_E_true_0_2_totInv_0_70+=mcp->getEnergy();
		  m_n_true_0_2_totInv_0_70+=mcp->getEnergy();
		}
	      }else if(mcp->getEnergy()<10.0){
		if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
		  m_E_true_2_10_totAll_0_70+=mcp->getEnergy();
		  m_n_true_2_10_totAll_0_70+=1;
		  if(abs(mcp->getPDG())==211){
		    m_E_true_2_10_totPi_0_70+=mcp->getEnergy();
		    m_n_true_2_10_totPi_0_70+=1;
		  }else if (mcp->getPDG()==22){
		    m_E_true_2_10_totPh_0_70+=mcp->getEnergy();
		    m_n_true_2_10_totPh_0_70+=1;
		  }else if (abs(mcp->getPDG())==130){
		    m_E_true_2_10_totK0L_0_70+=mcp->getEnergy();
		    m_n_true_2_10_totK0L_0_70+=1;
		  }else if(abs(mcp->getPDG())==11){
		    m_E_true_2_10_totE_0_70+=mcp->getEnergy();
		    m_n_true_2_10_totE_0_70+=1;
		  }else if(abs(mcp->getPDG())==13 ){
		    m_E_true_2_10_totMu_0_70+=mcp->getEnergy();
		    m_n_true_2_10_totMu_0_70+=1;
		  }else if(abs(mcp->getPDG())==2112 ){
		    m_E_true_2_10_totN_0_70+=mcp->getEnergy();
		    m_n_true_2_10_totN_0_70+=1;
		  }else if(abs(mcp->getPDG())==2212 ){
		    m_E_true_2_10_totP_0_70+=mcp->getEnergy();
		    m_n_true_2_10_totP_0_70+=1;
		  }else if(abs(mcp->getPDG())==321 ){
		    m_E_true_2_10_totK_0_70+=mcp->getEnergy();
		    m_n_true_2_10_totK_0_70+=1;
		  }else if (mcp->getCharge()==0){
		    m_E_true_2_10_totOtherNeut_0_70+=mcp->getEnergy();
		    m_n_true_2_10_totOtherNeut_0_70+=1;
		    std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
		  }else{
		    m_E_true_2_10_totOtherCH_0_70+=mcp->getEnergy();
		    m_n_true_2_10_totOtherCH_0_70+=1;
		    std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
		  }
		}else{
		  m_E_true_2_10_totInv_0_70+=mcp->getEnergy();
		  m_n_true_2_10_totInv_0_70+=mcp->getEnergy();
		}
	      }else if(mcp->getEnergy()<50.0){
		if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
		  m_E_true_10_50_totAll_0_70+=mcp->getEnergy();
		  m_n_true_10_50_totAll_0_70+=1;
		  if(abs(mcp->getPDG())==211){
		    m_E_true_10_50_totPi_0_70+=mcp->getEnergy();
		    m_n_true_10_50_totPi_0_70+=1;
		  }else if (mcp->getPDG()==22){
		    m_E_true_10_50_totPh_0_70+=mcp->getEnergy();
		    m_n_true_10_50_totPh_0_70+=1;
		  }else if (abs(mcp->getPDG())==130){
		    m_E_true_10_50_totK0L_0_70+=mcp->getEnergy();
		    m_n_true_10_50_totK0L_0_70+=1;
		  }else if(abs(mcp->getPDG())==11){
		    m_E_true_10_50_totE_0_70+=mcp->getEnergy();
		    m_n_true_10_50_totE_0_70+=1;
		  }else if(abs(mcp->getPDG())==13 ){
		    m_E_true_10_50_totMu_0_70+=mcp->getEnergy();
		    m_n_true_10_50_totMu_0_70+=1;
		  }else if(abs(mcp->getPDG())==2112 ){
		    m_E_true_10_50_totN_0_70+=mcp->getEnergy();
		    m_n_true_10_50_totN_0_70+=1;
		  }else if(abs(mcp->getPDG())==2212 ){
		    m_E_true_10_50_totP_0_70+=mcp->getEnergy();
		    m_n_true_10_50_totP_0_70+=1;
		  }else if(abs(mcp->getPDG())==321 ){
		    m_E_true_10_50_totK_0_70+=mcp->getEnergy();
		    m_n_true_10_50_totK_0_70+=1;
		  }else if (mcp->getCharge()==0){
		    m_E_true_10_50_totOtherNeut_0_70+=mcp->getEnergy();
		    m_n_true_10_50_totOtherNeut_0_70+=1;
		    std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
		  }else{
		    m_E_true_10_50_totOtherCH_0_70+=mcp->getEnergy();
		    m_n_true_10_50_totOtherCH_0_70+=1;
		    std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
		  }
		}else{
		  m_E_true_10_50_totInv_0_70+=mcp->getEnergy();
		  m_n_true_10_50_totInv_0_70+=mcp->getEnergy();
		}
	      }else{
		if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
		  m_E_true_50_totAll_0_70+=mcp->getEnergy();
		  m_n_true_50_totAll_0_70+=1;
		  if(abs(mcp->getPDG())==211){
		    m_E_true_50_totPi_0_70+=mcp->getEnergy();
		    m_n_true_50_totPi_0_70+=1;
		  }else if (mcp->getPDG()==22){
		    m_E_true_50_totPh_0_70+=mcp->getEnergy();
		    m_n_true_50_totPh_0_70+=1;
		  }else if (abs(mcp->getPDG())==130){
		    m_E_true_50_totK0L_0_70+=mcp->getEnergy();
		    m_n_true_50_totK0L_0_70+=1;
		  }else if(abs(mcp->getPDG())==11){
		    m_E_true_50_totE_0_70+=mcp->getEnergy();
		    m_n_true_50_totE_0_70+=1;
		  }else if(abs(mcp->getPDG())==13 ){
		    m_E_true_50_totMu_0_70+=mcp->getEnergy();
		    m_n_true_50_totMu_0_70+=1;
		  }else if(abs(mcp->getPDG())==2112 ){
		    m_E_true_50_totN_0_70+=mcp->getEnergy();
		    m_n_true_50_totN_0_70+=1;
		  }else if(abs(mcp->getPDG())==2212 ){
		    m_E_true_50_totP_0_70+=mcp->getEnergy();
		    m_n_true_50_totP_0_70+=1;
		  }else if(abs(mcp->getPDG())==321 ){
		    m_E_true_50_totK_0_70+=mcp->getEnergy();
		    m_n_true_50_totK_0_70+=1;
		  }else if (mcp->getCharge()==0){
		    m_E_true_50_totOtherNeut_0_70+=mcp->getEnergy();
		    m_n_true_50_totOtherNeut_0_70+=1;
		    std::cout<<"other type neutral PDG "<<mcp->getPDG()<<std::endl;
		  }else{
		    m_E_true_50_totOtherCH_0_70+=mcp->getEnergy();
		    m_n_true_50_totOtherCH_0_70+=1;
		    std::cout<<"other type charged PDG "<<mcp->getPDG()<<std::endl;
		  }
		}else{
		  m_E_true_50_totInv_0_70+=mcp->getEnergy();
		  m_n_true_50_totInv_0_70+=mcp->getEnergy();
		}
	      }//0.70 true values in energy bins
	    }
	  }//0.95 loop
	}//check for pi0 here - > put that now to another particle ID which we should never encounter
	//if(mcp->getPDG()==111){
	if(mcp->getPDG()==3000111){  
	if(mcp->getDaughters().size ()==2 && (mcp->getDaughters()[0]->getPDG() == 22 && mcp->getDaughters()[0]->getGeneratorStatus()==1) &&(mcp->getDaughters()[1]->getPDG() == 22 && mcp->getDaughters()[1]->getGeneratorStatus()==1)){
	    m_truePi0Energy->push_back(mcp->getEnergy());
	    m_truePi0_Px->push_back(mcp->getMomentum()[0]);
	    m_truePi0_Py->push_back(mcp->getMomentum()[1]);
	    m_truePi0_Pz->push_back(mcp->getMomentum()[2]);
	    m_truePi0_CosTheta->push_back(mcp->getMomentum()[2]/sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]+mcp->getMomentum()[2]*mcp->getMomentum()[2]));
	    m_truePi0_Phi->push_back(atan2(mcp->getMomentum()[1],mcp->getMomentum()[0]));
	    m_truePh0Energy->push_back(mcp->getDaughters()[0]->getEnergy());
	    m_truePh0_Px->push_back(mcp->getDaughters()[0]->getMomentum()[0]);
	    m_truePh0_Py->push_back(mcp->getDaughters()[0]->getMomentum()[1]);
	    m_truePh0_Pz->push_back(mcp->getDaughters()[0]->getMomentum()[2]);
	    TLorentzVector mctempPh0;
	    mctempPh0.SetPxPyPzE(mcp->getDaughters()[0]->getMomentum()[0],mcp->getDaughters()[0]->getMomentum()[1],mcp->getDaughters()[0]->getMomentum()[2],mcp->getDaughters()[0]->getEnergy());
	    m_truePh0_CosTheta->push_back(mctempPh0.CosTheta());
	    m_truePh0_Phi->push_back(mctempPh0.Phi());
	    m_truePh1Energy->push_back(mcp->getDaughters()[1]->getEnergy());
	    m_truePh1_Px->push_back(mcp->getDaughters()[1]->getMomentum()[0]);
	    m_truePh1_Py->push_back(mcp->getDaughters()[1]->getMomentum()[1]);
	    m_truePh1_Pz->push_back(mcp->getDaughters()[1]->getMomentum()[2]);
	    TLorentzVector mctempPh1;
	    mctempPh1.SetPxPyPzE(mcp->getDaughters()[1]->getMomentum()[0],mcp->getDaughters()[1]->getMomentum()[1],mcp->getDaughters()[1]->getMomentum()[2],mcp->getDaughters()[1]->getEnergy());
	    m_truePh1_CosTheta->push_back(mctempPh1.CosTheta());
	    m_truePh1_Phi->push_back(mctempPh1.Phi());

	    m_truePi0_DRPh01->push_back(mctempPh1.DeltaR(mctempPh0));
	    double cosAngPh01=((mctempPh0.Px()*mctempPh1.Px()+mctempPh0.Py()*mctempPh1.Py()+mctempPh1.Pz()*mctempPh0.Pz())/(mctempPh0.P()*mctempPh1.P()));
	    m_truePi0_CosAngPh01->push_back(cosAngPh01);
	    m_truePi0_AngPh01->push_back(acos(cosAngPh01));
	    int index_0DRMin=-1;
	    int index_0CosAngMax=-1;
	    int index_1DRMin=-1;
	    int index_1CosAngMax=-1;
	    double deltaR0_min=100;
	    double deltaR1_min=100;
	    double cosAng0_max=-1.1;
	    double cosAng1_max=-1.1;
	    double EisoDR01_ph0=0;
	    double EisoDR01_ph1=0;
	    double EisoCosAng0995_ph0=0;
	    double EisoCosAng0995_ph1=0;
	    int found_pi0_photons=0;
	    for(int m1 =0; m1< mcColl->getNumberOfElements(); m1++){
	      MCParticle* mcp2= dynamic_cast<MCParticle*>( mcColl->getElementAt(m1) ) ;
	      if(mcp2->getGeneratorStatus()==1){
		TLorentzVector mcptest;
		mcptest.SetPxPyPzE(mcp2->getMomentum()[0],mcp2->getMomentum()[1],mcp2->getMomentum()[2],mcp->getEnergy());
		double deltaRcheck0=mcptest.DeltaR(mctempPh1);
		double deltaRcheck1=mcptest.DeltaR(mctempPh0);
		double cosAngcheck0=((mctempPh0.Px()*mcptest.Px()+mctempPh0.Py()*mcptest.Py()+mcptest.Pz()*mctempPh0.Pz())/(mctempPh0.P()*mcptest.P()));
		double cosAngcheck1=((mctempPh1.Px()*mcptest.Px()+mctempPh1.Py()*mcptest.Py()+mcptest.Pz()*mctempPh1.Pz())/(mctempPh1.P()*mcptest.P()));
		if(deltaRcheck0==0){
		  found_pi0_photons+=1;
		}
		if(deltaRcheck1==0){
		  found_pi0_photons+=1;
		}
		if((deltaRcheck0!=0 && deltaRcheck1!=0)){
		  if(deltaRcheck0<deltaR0_min){
		    deltaR0_min=deltaRcheck0;
		    index_0DRMin=m1;
		  }
		  if(deltaRcheck0<0.1){
		    EisoDR01_ph0+=mcp->getEnergy();
		  }
		  if(deltaRcheck1<deltaR1_min){
		    deltaR1_min=deltaRcheck1;
		    index_1DRMin=m1;
		  } 
		  if(deltaRcheck1<0.1){
		    EisoDR01_ph1+=mcp->getEnergy();
		  }
		  //delta R overlap requirement should also cover the different angular criteria
		  if(cosAngcheck0>cosAng0_max){
		    cosAng0_max=cosAngcheck0;
		    index_0CosAngMax=m1;
		  }
		  if(cosAngcheck0>0.995){
		    EisoCosAng0995_ph0+=mcp->getEnergy();
		  }
		  if(cosAngcheck1>cosAng1_max){
		    cosAng1_max=cosAngcheck1;
		    index_1CosAngMax=m1;
		  }
		  if(cosAngcheck1>0.995){
		    EisoCosAng0995_ph1+=mcp->getEnergy();
		  }
		}
	      }
	    }
	    if(found_pi0_photons!=2){
	      //std::cout<<"assumed we should find indeed 2 photons not more not less, what did we find instead "<<found_pi0_photons<<std::endl;
	    }
	    m_truePh0_DR01E->push_back(EisoDR01_ph0);
	    m_truePh0_CosAng0995E->push_back(EisoCosAng0995_ph0);
	    if((index_0DRMin!=-1 && index_0CosAngMax ==-1)||(index_0DRMin==-1 && index_0CosAngMax !=-1) || (index_0DRMin!=-1 && index_1DRMin ==-1) || (index_0DRMin==-1 && index_1DRMin !=-1)){
	      //std::cout<<"something wrong in closest particle definition of true photons"<<index_0DRMin<<"/"<<index_0CosAngMax <<"/"<<index_1DRMin<<"/"<<index_1CosAngMax<<std::endl;
	    }
	    if(index_0DRMin!=-1){//if true photon 0 finds any other particle (other than photon 1), then this particle is found by photon 1 automatically too
	      MCParticle* mcp0_isoDR=dynamic_cast<MCParticle*>( mcColl->getElementAt(index_0DRMin)); 
	      m_truePh0_DRMin->push_back(deltaR0_min);
	      m_truePh0_DRMinPDGID->push_back(mcp0_isoDR->getPDG());
	      m_truePh0_DRMinE->push_back(mcp0_isoDR->getEnergy());
	      MCParticle* mcp1_isoDR=dynamic_cast<MCParticle*>( mcColl->getElementAt(index_1DRMin)); 
	      m_truePh1_DRMin->push_back(deltaR1_min);
	      m_truePh1_DRMinPDGID->push_back(mcp1_isoDR->getPDG());
	      m_truePh1_DRMinE->push_back(mcp1_isoDR->getEnergy());
	      MCParticle* mcp0_isoCosAng=dynamic_cast<MCParticle*>( mcColl->getElementAt(index_0CosAngMax)); 
	      m_truePh0_CosAngMax->push_back(cosAng0_max);
	      m_truePh0_CosAngMaxPDGID->push_back(mcp0_isoCosAng->getPDG());
	      m_truePh0_CosAngMaxE->push_back(mcp0_isoCosAng->getEnergy());
	      MCParticle* mcp1_isoCosAng=dynamic_cast<MCParticle*>( mcColl->getElementAt(index_1CosAngMax)); 
	      m_truePh1_CosAngMax->push_back(cosAng1_max);
	      m_truePh1_CosAngMaxPDGID->push_back(mcp1_isoCosAng->getPDG());
	      m_truePh1_CosAngMaxE->push_back(mcp1_isoCosAng->getEnergy());
	    }else{
	      m_truePh0_DRMin->push_back(-1);
	      m_truePh0_DRMinPDGID->push_back(10);
	      m_truePh0_DRMinE->push_back(-1);
	      m_truePh1_DRMin->push_back(-1);
	      m_truePh1_DRMinPDGID->push_back(10);
	      m_truePh1_DRMinE->push_back(-1);

	      m_truePh0_CosAngMax->push_back(-1);
	      m_truePh0_CosAngMaxPDGID->push_back(10);
	      m_truePh0_CosAngMaxE->push_back(-1);
	      m_truePh1_CosAngMax->push_back(-1);
	      m_truePh1_CosAngMaxPDGID->push_back(10);
	      m_truePh1_CosAngMaxE->push_back(-1);

	    }
	    m_truePh1_DR01E->push_back(EisoDR01_ph1);
	    m_truePh1_CosAng0995E->push_back(EisoCosAng0995_ph1);
	  }
	}
    }
    //std::cout<<"energy values here "<<   m_E_true_totAll<<"/"<<    m_E_true_totInv<<"/"<<    m_E_true_totPi<<"/"<<    m_E_true_totPh<<"/"<<    m_E_true_totK0L<<"/"<<    m_E_true_totE<<"/"<<    m_E_true_totMu<<"/"<<    m_E_true_totN<<"/"<<    m_E_true_totK<<"/"<<    m_E_true_totP<<"/"<<    m_E_true_totOtherCH<<"/"<<    m_E_true_totOtherNeut<<"/"<<std::endl;
    if((m_E_true_totAll-m_Z_mcE)>20){
      //std::cout<<"run/event/ great difference "<<evt->getRunNumber()<<"/"<<evt->getEventNumber()<<"/"<<m_E_true_totAll-m_Z_mcE<<std::endl;
    }
    }
    //up to here the gen step should have been done, check if the original Z had been found

    if(m_Z_mcE==-10){//Z had not been found, now check if we are e.g. in a bbar sample
      //std::cout<<"after mcColl"<<std::endl;
      if(mcColl!=NULL){
	//std::cout<<"in mcColl"<<std::endl;
	//Look for a photon
	bool found_bbar=false;
	for(int m =0; m< mcColl->getNumberOfElements(); m++){
	  MCParticle* mcp= dynamic_cast<MCParticle*>( mcColl->getElementAt(m) ) ;
	  //quarks originate from 11 or -11 which are saved as producing the bbar event
	  if(found_bbar){
	    break;
	  }
	  if(abs(mcp->getPDG())==11 && mcp->getDaughters().size ()==2){
	    if(abs(mcp->getDaughters()[0]->getPDG())==5 && abs(mcp->getDaughters()[1]->getPDG())==5){
	      found_bbar=true;
	      m_d1_mcPDGID=mcp->getDaughters()[0]->getPDG();
	      m_d1_mcE=mcp->getDaughters()[0]->getEnergy();
	      m_d1_mcPx=mcp->getDaughters()[0]->getMomentum()[0];
	      m_d1_mcPy=mcp->getDaughters()[0]->getMomentum()[1];
	      m_d1_mcPz=mcp->getDaughters()[0]->getMomentum()[2];
	      m_d1_mcMass=mcp->getDaughters()[0]->getMass();
	      m_d1_mcPhi=atan2(m_d1_mcPy,m_d1_mcPx);
	      m_d1_mcCosTheta= m_d1_mcPz/sqrt(m_d1_mcPx*m_d1_mcPx+m_d1_mcPy*m_d1_mcPy+m_d1_mcPz*m_d1_mcPz);
	      m_d1_mcTheta=acos(m_d1_mcCosTheta);
	      m_d2_mcPDGID=mcp->getDaughters()[1]->getPDG();
	      m_d2_mcE=mcp->getDaughters()[1]->getEnergy();
	      m_d2_mcPx=mcp->getDaughters()[1]->getMomentum()[0];
	      m_d2_mcPy=mcp->getDaughters()[1]->getMomentum()[1];
	      m_d2_mcPz=mcp->getDaughters()[1]->getMomentum()[2];
	      m_d2_mcMass=mcp->getDaughters()[1]->getMass();
	      m_d2_mcPhi=atan2(m_d2_mcPy,m_d2_mcPx);
	      m_d2_mcCosTheta= m_d2_mcPz/sqrt(m_d2_mcPx*m_d2_mcPx+m_d2_mcPy*m_d2_mcPy+m_d2_mcPz*m_d2_mcPz);
	      m_d2_mcTheta=acos(m_d2_mcCosTheta);
	    }
	  }
	}
      }
    }

    
    //now reco PFO loop
    LCCollection* recoparticlecol = NULL;
    // Alternativelly if you do not want Marlin to exit in case of a non-existing collection
    // use the following (commented out) code:
    //run on H to gamma gamma -> in case there are no tracks around linker will fail to produce output collection
    
    //std::cout<<"before recColl"<<std::endl;
    recoparticlecol = evt->getCollection(m_inputRECOParticleCollection) ;
    // std::cout<<"after mcColl"<<std::endl;
    if(recoparticlecol!=NULL){
      //std::cout<<"in recColl"<<std::endl;
      //PandoraCandidate loop
      for(int i=0;i<recoparticlecol->getNumberOfElements();i++){
	ReconstructedParticle* pandorapart = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(i));
	double cosTheta = pandorapart->getMomentum()[2]/sqrt(pandorapart->getMomentum()[0]*pandorapart->getMomentum()[0]+pandorapart->getMomentum()[1]*pandorapart->getMomentum()[1]+pandorapart->getMomentum()[2]*pandorapart->getMomentum()[2]);
	m_E_totPFO+=pandorapart->getEnergy();
	m_px_totPFO+=pandorapart->getMomentum()[0];
	m_py_totPFO+=pandorapart->getMomentum()[1];
	m_pz_totPFO+=pandorapart->getMomentum()[2];
	if(abs(pandorapart->getType())==211){
	  m_E_totPi+=pandorapart->getEnergy();
	  m_px_totPi+=pandorapart->getMomentum()[0];
	  m_py_totPi+=pandorapart->getMomentum()[1];
	  m_pz_totPi+=pandorapart->getMomentum()[2];
	  m_n_totPi+=1;
	}else if (pandorapart->getType()==22){
	  m_E_totPh+=pandorapart->getEnergy();
	  m_px_totPh+=pandorapart->getMomentum()[0];
	  m_py_totPh+=pandorapart->getMomentum()[1];
	  m_pz_totPh+=pandorapart->getMomentum()[2];
	  m_n_totPh+=1;
	}else if(abs(pandorapart->getType())==11){
	  m_E_totE+=pandorapart->getEnergy();
	  m_px_totE+=pandorapart->getMomentum()[0];
	  m_py_totE+=pandorapart->getMomentum()[1];
	  m_pz_totE+=pandorapart->getMomentum()[2];
	  m_n_totE+=1;
	}else if(abs(pandorapart->getType())==13){
	  m_E_totMu+=pandorapart->getEnergy();
	  m_px_totMu+=pandorapart->getMomentum()[0];
	  m_py_totMu+=pandorapart->getMomentum()[1];
	  m_pz_totMu+=pandorapart->getMomentum()[2];
	  m_n_totMu+=1;
	}else if(pandorapart->getType()==2112){
	  m_E_totN+=pandorapart->getEnergy();
	  m_px_totN+=pandorapart->getMomentum()[0];
	  m_py_totN+=pandorapart->getMomentum()[1];
	  m_pz_totN+=pandorapart->getMomentum()[2];
	  m_n_totN+=1;
	}else{
	  std::cout<<"did we maybe miss lambdas "<<pandorapart->getType()<<std::endl;
	}
	//std::cout<<"after all energies"<<std::endl;
	if(pandorapart->getEnergy()<2){
	  m_E_0_2_totPFO+=pandorapart->getEnergy();
	  if(abs(pandorapart->getType())==211){
	    m_E_0_2_totPi+=pandorapart->getEnergy();
	    m_n_0_2_totPi+=1;
	  }else if (pandorapart->getType()==22){
	    m_E_0_2_totPh+=pandorapart->getEnergy();
	    m_n_0_2_totPh+=1;
	  }else if(abs(pandorapart->getType())==11){
	    m_E_0_2_totE+=pandorapart->getEnergy();
	    m_n_0_2_totE+=1;
	  }else if(abs(pandorapart->getType())==13){
	    m_E_0_2_totMu+=pandorapart->getEnergy();
	    m_n_0_2_totMu+=1;
	  }else if(pandorapart->getType()==2112){
	    m_E_0_2_totN+=pandorapart->getEnergy();
	    m_n_0_2_totN+=1;
	  }
	  //std::cout<<"after 2 energies"<<std::endl;
	}else if(pandorapart->getEnergy()<10){
	  m_E_2_10_totPFO+=pandorapart->getEnergy();
	  if(abs(pandorapart->getType())==211){
	    m_E_2_10_totPi+=pandorapart->getEnergy();
	    m_n_2_10_totPi+=1;
	  }else if (pandorapart->getType()==22){
	    m_E_2_10_totPh+=pandorapart->getEnergy();
	    m_n_2_10_totPh+=1;
	  }else if(abs(pandorapart->getType())==11){
	    m_E_2_10_totE+=pandorapart->getEnergy();
	    m_n_2_10_totE+=1;
	  }else if(abs(pandorapart->getType())==13){
	    m_E_2_10_totMu+=pandorapart->getEnergy();
	    m_n_2_10_totMu+=1;
	  }else if(pandorapart->getType()==2112){
	    m_E_2_10_totN+=pandorapart->getEnergy();
	    m_n_2_10_totN+=1;
	  }
	  //std::cout<<"after 10 energies"<<std::endl;
	}else if(pandorapart->getEnergy()<50){
	  m_E_10_50_totPFO+=pandorapart->getEnergy();
	  if(abs(pandorapart->getType())==211){
	    m_E_10_50_totPi+=pandorapart->getEnergy();
	    m_n_10_50_totPi+=1;
	  }else if (pandorapart->getType()==22){
	    m_E_10_50_totPh+=pandorapart->getEnergy();
	    m_n_10_50_totPh+=1;
	  }else if(abs(pandorapart->getType())==11){
	    m_E_10_50_totE+=pandorapart->getEnergy();
	    m_n_10_50_totE+=1;
	  }else if(abs(pandorapart->getType())==13){
	    m_E_10_50_totMu+=pandorapart->getEnergy();
	    m_n_10_50_totMu+=1;
	  }else if(pandorapart->getType()==2112){
	    m_E_10_50_totN+=pandorapart->getEnergy();
	    m_n_10_50_totN+=1;
	  }
	  //std::cout<<"after 50 energies"<<std::endl;
	}else{
	  m_E_50_totPFO+=pandorapart->getEnergy();
	  if(abs(pandorapart->getType())==211){
	    m_E_50_totPi+=pandorapart->getEnergy();
	    m_n_50_totPi+=1;
	  }else if (pandorapart->getType()==22){
	    m_E_50_totPh+=pandorapart->getEnergy();
	    m_n_50_totPh+=1;
	  }else if(abs(pandorapart->getType())==11){
	    m_E_50_totE+=pandorapart->getEnergy();
	    m_n_50_totE+=1;
	  }else if(abs(pandorapart->getType())==13){
	    m_E_50_totMu+=pandorapart->getEnergy();
	    m_n_50_totMu+=1;
	  }else if(pandorapart->getType()==2112){
	    m_E_50_totN+=pandorapart->getEnergy();
	    m_n_50_totN+=1;
	  }
	  //std::cout<<"after else energies"<<std::endl;
	}
	if(fabs(cosTheta)<0.95){
	  m_E_totPFO_0_95+=pandorapart->getEnergy();
	  m_n_totPFO_0_95+=1;
	  if(abs(pandorapart->getType())==211){
	    m_E_totPi_0_95+=pandorapart->getEnergy();
	    m_n_totPi_0_95+=1;
	  }else if (pandorapart->getType()==22){
	    m_E_totPh_0_95+=pandorapart->getEnergy();
	    m_n_totPh_0_95+=1;
	  }else if(abs(pandorapart->getType())==11){
	    m_E_totE_0_95+=pandorapart->getEnergy();
	    m_n_totE_0_95+=1;
	  }else if(abs(pandorapart->getType())==13){
	    m_E_totMu_0_95+=pandorapart->getEnergy();
	    m_n_totMu_0_95+=1;
	  }else if(pandorapart->getType()==2112){
	    m_E_totN_0_95+=pandorapart->getEnergy();
	    m_n_totN_0_95+=1;
	  }
	  //std::cout<<"after 095 all energies"<<std::endl;
	  if(pandorapart->getEnergy()<2){
	    m_E_0_2_totPFO_0_95+=pandorapart->getEnergy();
	    if(abs(pandorapart->getType())==211){
	      m_E_0_2_totPi_0_95+=pandorapart->getEnergy();
	      m_n_0_2_totPi_0_95+=1;
	    }else if (pandorapart->getType()==22){
	      m_E_0_2_totPh_0_95+=pandorapart->getEnergy();
	      m_n_0_2_totPh_0_95+=1;
	    }else if(abs(pandorapart->getType())==11){
	      m_E_0_2_totE_0_95+=pandorapart->getEnergy();
	      m_n_0_2_totE_0_95+=1;
	    }else if(abs(pandorapart->getType())==13){
	      m_E_0_2_totMu_0_95+=pandorapart->getEnergy();
	      m_n_0_2_totMu_0_95+=1;
	    }else if(pandorapart->getType()==2112){
	      m_E_0_2_totN_0_95+=pandorapart->getEnergy();
	      m_n_0_2_totN_0_95+=1;
	    }
	    //std::cout<<"after 095 2 energies"<<std::endl;
	  }else if(pandorapart->getEnergy()<10){
	    m_E_2_10_totPFO_0_95+=pandorapart->getEnergy();
	    if(abs(pandorapart->getType())==211){
	      m_E_2_10_totPi_0_95+=pandorapart->getEnergy();
	      m_n_2_10_totPi_0_95+=1;
	    }else if (pandorapart->getType()==22){
	      m_E_2_10_totPh_0_95+=pandorapart->getEnergy();
	      m_n_2_10_totPh_0_95+=1;
	    }else if(abs(pandorapart->getType())==11){
	      m_E_2_10_totE_0_95+=pandorapart->getEnergy();
	      m_n_2_10_totE_0_95+=1;
	    }else if(abs(pandorapart->getType())==13){
	      m_E_2_10_totMu_0_95+=pandorapart->getEnergy();
	      m_n_2_10_totMu_0_95+=1;
	    }else if(pandorapart->getType()==2112){
	      m_E_2_10_totN_0_95+=pandorapart->getEnergy();
	      m_n_2_10_totN_0_95+=1;
	    }
	    //std::cout<<"after 095 a10 energies"<<std::endl;
	  }else if(pandorapart->getEnergy()<50){
	    m_E_10_50_totPFO_0_95+=pandorapart->getEnergy();
	    if(abs(pandorapart->getType())==211){
	      m_E_10_50_totPi_0_95+=pandorapart->getEnergy();
	      m_n_10_50_totPi_0_95+=1;
	    }else if (pandorapart->getType()==22){
	      m_E_10_50_totPh_0_95+=pandorapart->getEnergy();
	      m_n_10_50_totPh_0_95+=1;
	    }else if(abs(pandorapart->getType())==11){
	      m_E_10_50_totE_0_95+=pandorapart->getEnergy();
	      m_n_10_50_totE_0_95+=1;
	    }else if(abs(pandorapart->getType())==13){
	      m_E_10_50_totMu_0_95+=pandorapart->getEnergy();
	      m_n_10_50_totMu_0_95+=1;
	    }else if(pandorapart->getType()==2112){
	      m_E_10_50_totN_0_95+=pandorapart->getEnergy();
	      m_n_10_50_totN_0_95+=1;
	    }
	    //std::cout<<"after 09510 energies"<<std::endl;
	  }else{
	    m_E_50_totPFO_0_95+=pandorapart->getEnergy();
	    if(abs(pandorapart->getType())==211){
	      m_E_50_totPi_0_95+=pandorapart->getEnergy();
	      m_n_50_totPi_0_95+=1;
	    }else if (pandorapart->getType()==22){
	      m_E_50_totPh_0_95+=pandorapart->getEnergy();
	      m_n_50_totPh_0_95+=1;
	    }else if(abs(pandorapart->getType())==11){
	      m_E_50_totE_0_95+=pandorapart->getEnergy();
	      m_n_50_totE_0_95+=1;
	    }else if(abs(pandorapart->getType())==13){
	      m_E_50_totMu_0_95+=pandorapart->getEnergy();
	      m_n_50_totMu_0_95+=1;
	    }else if(pandorapart->getType()==2112){
	      m_E_50_totN_0_95+=pandorapart->getEnergy();
	      m_n_50_totN_0_95+=1;
	    }
	    //std::cout<<"after 095else energies"<<std::endl;
	  }//0 95 PFO energy and multiplicity values in energy bins
	  if(fabs(cosTheta)<0.70){
	    m_E_totPFO_0_70+=pandorapart->getEnergy();
	    m_n_totPFO_0_70+=1;
	    if(abs(pandorapart->getType())==211){
	      m_E_totPi_0_70+=pandorapart->getEnergy();
	      m_n_totPi_0_70+=1;
	    }else if (pandorapart->getType()==22){
	      m_E_totPh_0_70+=pandorapart->getEnergy();
	      m_n_totPh_0_70+=1;
	    }else if(abs(pandorapart->getType())==11){
	      m_E_totE_0_70+=pandorapart->getEnergy();
	      m_n_totE_0_70+=1;
	    }else if(abs(pandorapart->getType())==13){
	      m_E_totMu_0_70+=pandorapart->getEnergy();
	      m_E_totMu_0_70+=1;
	    }else if(pandorapart->getType()==2112){
	      m_E_totN_0_70+=pandorapart->getEnergy();
	      m_E_totN_0_70+=1;
	    }
	    if(pandorapart->getEnergy()<2){
	      m_E_0_2_totPFO_0_70+=pandorapart->getEnergy();
	      if(abs(pandorapart->getType())==211){
		m_E_0_2_totPi_0_70+=pandorapart->getEnergy();
		m_n_0_2_totPi_0_70+=1;
	      }else if (pandorapart->getType()==22){
		m_E_0_2_totPh_0_70+=pandorapart->getEnergy();
		m_n_0_2_totPh_0_70+=1;
	      }else if(abs(pandorapart->getType())==11){
		m_E_0_2_totE_0_70+=pandorapart->getEnergy();
		m_n_0_2_totE_0_70+=1;
	      }else if(abs(pandorapart->getType())==13){
		m_E_0_2_totMu_0_70+=pandorapart->getEnergy();
		m_n_0_2_totMu_0_70+=1;
	      }else if(pandorapart->getType()==2112){
		m_E_0_2_totN_0_70+=pandorapart->getEnergy();
		m_n_0_2_totN_0_70+=1;
	      }
	    }else if(pandorapart->getEnergy()<10){
	      m_E_2_10_totPFO_0_70+=pandorapart->getEnergy();
	      if(abs(pandorapart->getType())==211){
		m_E_2_10_totPi_0_70+=pandorapart->getEnergy();
		m_n_2_10_totPi_0_70+=1;
	      }else if (pandorapart->getType()==22){
		m_E_2_10_totPh_0_70+=pandorapart->getEnergy();
		m_n_2_10_totPh_0_70+=1;
	      }else if(abs(pandorapart->getType())==11){
		m_E_2_10_totE_0_70+=pandorapart->getEnergy();
		m_n_2_10_totE_0_70+=1;
	      }else if(abs(pandorapart->getType())==13){
		m_E_2_10_totMu_0_70+=pandorapart->getEnergy();
		m_n_2_10_totMu_0_70+=1;
	      }else if(pandorapart->getType()==2112){
		m_E_2_10_totN_0_70+=pandorapart->getEnergy();
		m_n_2_10_totN_0_70+=1;
	      }
	    }else if(pandorapart->getEnergy()<50){
	      m_E_10_50_totPFO_0_70+=pandorapart->getEnergy();
	      if(abs(pandorapart->getType())==211){
		m_E_10_50_totPi_0_70+=pandorapart->getEnergy();
		m_n_10_50_totPi_0_70+=1;
	      }else if (pandorapart->getType()==22){
		m_E_10_50_totPh_0_70+=pandorapart->getEnergy();
		m_n_10_50_totPh_0_70+=1;
	      }else if(abs(pandorapart->getType())==11){
		m_E_10_50_totE_0_70+=pandorapart->getEnergy();
		m_n_10_50_totE_0_70+=1;
	      }else if(abs(pandorapart->getType())==13){
		m_E_10_50_totMu_0_70+=pandorapart->getEnergy();
		m_n_10_50_totMu_0_70+=1;
	      }else if(pandorapart->getType()==2112){
		m_E_10_50_totN_0_70+=pandorapart->getEnergy();
		m_n_10_50_totN_0_70+=1;
	      }
	    }else{
	      m_E_50_totPFO_0_70+=pandorapart->getEnergy();
	      if(abs(pandorapart->getType())==211){
		m_E_50_totPi_0_70+=pandorapart->getEnergy();
		m_n_50_totPi_0_70+=1;
	      }else if (pandorapart->getType()==22){
		m_E_50_totPh_0_70+=pandorapart->getEnergy();
		m_n_50_totPh_0_70+=1;
	      }else if(abs(pandorapart->getType())==11){
		m_E_50_totE_0_70+=pandorapart->getEnergy();
		m_n_50_totE_0_70+=1;
	      }else if(abs(pandorapart->getType())==13){
		m_E_50_totMu_0_70+=pandorapart->getEnergy();
		m_n_50_totMu_0_70+=1;
	      }else if(pandorapart->getType()==2112){
		m_E_50_totN_0_70+=pandorapart->getEnergy();
		m_n_50_totN_0_70+=1;
	      }
	    }//0 70 PFO energy and multiplicity values in energy bins
	  }//0.70 loop
	}//0.95 loop
	//check for photon here - > put that now to another particle ID which we should never encounter                                                                                  
        //if(pandorapart->getType()==22){	
	if(pandorapart->getType()==3000111){
	    //std::cout<<"in photon checking"<<std::endl;
	  //Check if in barrel
	  m_recoPhEnergy->push_back(pandorapart->getEnergy());
	  m_recoPh_Px->push_back(pandorapart->getMomentum()[0]);
	  m_recoPh_Py->push_back(pandorapart->getMomentum()[1]);
	  m_recoPh_Pz->push_back( pandorapart->getMomentum()[2]);
	  m_recoPh_CosTheta->push_back(cosTheta);    
	  //m_true_Phi=atan2(m_true_Py,m_true_Px);
	  m_recoPh_PDGID->push_back(pandorapart->getType());
	  float reco_px=pandorapart->getMomentum()[0];
	  float reco_py=pandorapart->getMomentum()[1];
	  //m_recoPh_Phi->push_back(atan2(pandorapart->getMomentum()[1],pandorapart->getMomentum()[0]));
	  m_recoPh_Phi->push_back(atan2(reco_py,reco_px));
	  m_recoPh_Theta->push_back(acos(cosTheta));
	  float En_ECAL_Barrel=0;
	  float En_ECAL_Endcap=0;
	  float En_ECAL_else=0;
	  float En_HCAL_Barrel=0;
	  float En_HCAL_Endcap=0;
	  float En_HCAL_else=0;
	  int phFirstLayerHCAL=-1;
	  int phFirstLayerECAL=-1;
	  int phLastLayerHCAL=-1;
	  int phLastLayerECAL=-1;
	  int hitsEB=0;
	  int hitsEE=0;
	  int hitsEO=0;
	  int hitsHB=0;
	  int hitsHE=0;
	  int hitsHO=0;
	  //std::cout<<"in photon with clusters"<<std::endl;
	  //two clusters per particle doesn't seem to happen, but maybe in future
	  for(unsigned int j=0;j<pandorapart->getClusters().size();j++){
	    //std::cout<<"in photon with calohits"<<std::endl;
	    for(unsigned int l=0; l<pandorapart->getClusters()[j]->getCalorimeterHits().size();l++){
	      int types=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getType();
	      //static const int fCaloType =     1 ;
	      static const int fCaloID   =    10 ;
	      static const int fLayout   =  1000 ;
	      static const int fLayer    = 10000 ;
	      int caloLayer=types/fLayer;
	      int caloLayout=(types-caloLayer*fLayer)/fLayout;
	      int caloID=(types-caloLayer*fLayer-caloLayout*fLayout)/fCaloID;
	      //int caloType=(types-caloLayer*fLayer-caloLayout*fLayout-caloID*fCaloID)/fCaloType;		  
	      if(caloID==1){//ecal
		if(phFirstLayerECAL<0){
		  phFirstLayerECAL=caloLayer;
		}else if (caloLayer<phFirstLayerECAL){
		  phFirstLayerECAL=caloLayer;
		}
		if(phLastLayerECAL<0){
		  phLastLayerECAL=caloLayer;
		}else if (caloLayer>phLastLayerECAL){
		  phLastLayerECAL=caloLayer;
		}
		if(caloLayout==1){//id 1=ecal, 2 hcal, layout 1 barrel, 2 endcap
		  En_ECAL_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  hitsEB+=1;
		}else if(caloLayout==2){
		  En_ECAL_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  hitsEE+=1;
		}else{
		  En_ECAL_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  hitsEO+=1;
		}
	      }else if(caloID==2){//hcal
		if(phFirstLayerHCAL<0){
		  phFirstLayerHCAL=caloLayer;
		}else if (caloLayer<phFirstLayerHCAL){
		  phFirstLayerHCAL=caloLayer;
		}
		if(phLastLayerHCAL<0){
		  phLastLayerHCAL=caloLayer;
		}else if (caloLayer>phLastLayerHCAL){
		  phLastLayerHCAL=caloLayer;
		}
		if(caloLayout==1){//id 1=ecal, 2 hcal, layout 1 barrel, 2 endcap
		  En_HCAL_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  hitsHB+=1;
		}else if(caloLayout==2){
		  En_HCAL_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  hitsHE+=1;
		}else{
		  En_HCAL_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  hitsHO+=1;
		}
	      }
	    }//hit loop
	  }//cluster loop
	  //need three indices: the two closest photon 
	  //next closest particle other than the photon
	  int ind_DR_ph0=-1;
	  int ind_DR_ph1=-1;
	  int ind_DR_other0=-1;
	  float DR_ph0=100;
	  float DR_ph1=100;
	  float DR_other0=100;
	  int ind_CosAngMax_ph0=-1;
	  int ind_CosAngMax_ph1=-1;
	  int ind_CosAngMax_other0=-1;
	  float CosAngMax_ph0=-0.1;
	  float CosAngMax_ph1=-0.1;
	  float CosAngMax_other0=-0.1;
	  float DR01_iso_E=0;
	  float CosAngleMax0995_iso_E=0;
	  //std::cout<<"for photon index loop"<<std::endl;
	  TLorentzVector recoPh;
	  recoPh.SetPxPyPzE(pandorapart->getMomentum()[0],pandorapart->getMomentum()[1],pandorapart->getMomentum()[2],pandorapart->getEnergy());
	  for(int i1=0;i1<recoparticlecol->getNumberOfElements();i1++){
	    ReconstructedParticle* pandorapart_check = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(i1));
	    TLorentzVector recotest;
	    recotest.SetPxPyPzE(pandorapart_check->getMomentum()[0],pandorapart_check->getMomentum()[1],pandorapart_check->getMomentum()[2],pandorapart_check->getEnergy());
	    float cosAngle=fabs( (recotest.Px()*recoPh.Px()+recotest.Py()*recoPh.Py()+recotest.Pz()*recoPh.Pz())/(sqrt(recoPh.Px()*recoPh.Px()+recoPh.Py()*recoPh.Py()+recoPh.Pz()*recoPh.Pz())*sqrt(recotest.Px()*recotest.Px()+recotest.Py()*recotest.Py()+recotest.Pz()*recotest.Pz())));
	    if(i!=i1){
	      if(recotest.DeltaR(recoPh)<0.10){
		DR01_iso_E+=pandorapart_check->getEnergy();
	      }
	      if(cosAngle>0.995){
		CosAngleMax0995_iso_E+=pandorapart_check->getEnergy();
	      }
	      if(pandorapart->getType()==22){
		if(recotest.DeltaR(recoPh)<DR_ph0){
		  //now the DR AND CosTheta values to check for closest values
		  ind_DR_ph1=ind_DR_ph0;
		  ind_DR_ph0=i1;
		  DR_ph1=DR_ph0;
		  DR_ph0=recotest.DeltaR(recoPh);
		}else if (recotest.DeltaR(recoPh)<DR_ph1){
		  DR_ph1=recotest.DeltaR(recoPh);
		  ind_DR_ph1=i1;
		}
		if(cosAngle>CosAngMax_ph0){
		  //now the DR AND CosTheta values to check for closest values
		  ind_CosAngMax_ph1=ind_CosAngMax_ph0;
		  ind_CosAngMax_ph0=i1;
		  CosAngMax_ph1=CosAngMax_ph0;
		  CosAngMax_ph0=cosAngle;
		}else if (cosAngle>CosAngMax_ph1){
		  CosAngMax_ph1=cosAngle;
		  ind_CosAngMax_ph1=i1;
		}
	      }
	    }else{
	      if(cosAngle>CosAngMax_other0){
		ind_DR_other0=i1;
		CosAngMax_other0=cosAngle;
	      }
	      if(recotest.DeltaR(recoPh)<DR_other0){
		//now the DR AND CosTheta values to check for closest values
		ind_DR_other0=i1;
		DR_other0=recotest.DeltaR(recoPh);
	      }
	    }
	  }//isolation values, closeness of other particles
	  m_recoPh_DR01_E->push_back(DR01_iso_E);
	  m_recoPh_CosAng0995_E->push_back(CosAngleMax0995_iso_E);
	  if(ind_DR_ph0>0){
	    //std::cout<<"ind 1 fill"<<std::endl;
	    ReconstructedParticle* pandoraph0 = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(ind_DR_ph0));
	    TLorentzVector recoph0;
	    recoph0.SetPxPyPzE(pandoraph0->getMomentum()[0],pandoraph0->getMomentum()[1],pandoraph0->getMomentum()[2],pandoraph0->getEnergy());
	    m_recoPh_DRMin_Ph->push_back(DR_ph0);
	    m_recoPh_DRMin_E_Ph->push_back(pandoraph0->getEnergy());
	    m_recoPi0Cand_DRMin_M->push_back((recoph0+recoPh).M());
	  }else{
	    m_recoPh_DRMin_Ph->push_back(-1);
	    m_recoPh_DRMin_E_Ph->push_back(-1);
	    m_recoPi0Cand_DRMin_M->push_back(-1);
	  }
	  if(DR_ph1<DR_other0){
	    //std::cout<<"ind 1 ph1 fill"<<std::endl;
	    ReconstructedParticle* pandoraphph1 = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(ind_DR_ph1));
	    m_recoPh_DRMin->push_back(DR_ph1);
	    m_recoPh_DRMin_PDGID->push_back(22);
	    m_recoPh_DRMin_E->push_back(pandoraphph1->getEnergy());
	  }else if (ind_DR_other0>=0){
	    // std::cout<<"ind 1 o fill"<<std::endl;
	    ReconstructedParticle* pandoraphother = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(ind_DR_other0));
	    m_recoPh_DRMin->push_back(DR_other0);
	    m_recoPh_DRMin_PDGID->push_back(pandoraphother->getType());
	    m_recoPh_DRMin_E->push_back(pandoraphother->getEnergy());
	  }else{
	    // std::cout<<"ind 1 else fill"<<std::endl;
	    m_recoPh_DRMin->push_back(-1);
	    m_recoPh_DRMin_PDGID->push_back(-1);
	    m_recoPh_DRMin_E->push_back(-1);
	  }
 
	  if(ind_CosAngMax_ph0>0){
	    // std::cout<<"ind cos max fill"<<std::endl;
	    ReconstructedParticle* pandoraCAph0 = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(ind_CosAngMax_ph0));
	    TLorentzVector recoCAph0;
	    recoCAph0.SetPxPyPzE(pandoraCAph0->getMomentum()[0],pandoraCAph0->getMomentum()[1],pandoraCAph0->getMomentum()[2],pandoraCAph0->getEnergy());
	    m_recoPh_CosAngMax_Ph->push_back(CosAngMax_ph0);
	    m_recoPh_CosAngMax_E_Ph->push_back(pandoraCAph0->getEnergy());
	    m_recoPi0Cand_CosAngMax_M->push_back((recoCAph0+recoPh).M());
	  }else{
	    m_recoPh_CosAngMax_Ph->push_back(-1);
	    m_recoPh_CosAngMax_E_Ph->push_back(-1);
	    m_recoPi0Cand_CosAngMax_M->push_back(-1);
	  }
	  if(CosAngMax_ph1<CosAngMax_other0 && ind_CosAngMax_ph1>0){
	    // std::cout<<"ind cos max 1 fill"<<std::endl;
	    ReconstructedParticle* pandoraCAphph1 = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(ind_CosAngMax_ph1));
	    m_recoPh_CosAngMax->push_back(CosAngMax_ph1);
	    m_recoPh_CosAngMax_PDGID->push_back(22);
	    m_recoPh_CosAngMax_E->push_back(pandoraCAphph1->getEnergy());
	  }else if (ind_CosAngMax_other0>=0){
            //std::cout<<"ind cos max else fill"<<std::endl;
	    ReconstructedParticle* pandoraCAphother = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(ind_CosAngMax_other0));
	    m_recoPh_CosAngMax->push_back(CosAngMax_other0);
	    m_recoPh_CosAngMax_PDGID->push_back(pandoraCAphother->getType());
	    m_recoPh_CosAngMax_E->push_back(pandoraCAphother->getEnergy());
	  }else{
	    m_recoPh_CosAngMax->push_back(-1);
	    m_recoPh_CosAngMax_PDGID->push_back(-1);
	    m_recoPh_CosAngMax_E->push_back(-1);
	  }
	  // std::cout<<"photon fill" <<std::endl;
	  m_recoPh_E_EB->push_back(En_ECAL_Barrel);
	  m_recoPh_E_EE->push_back(En_ECAL_Endcap);
	  m_recoPh_E_EO->push_back(En_ECAL_else);
	  m_recoPh_E_HB->push_back(En_HCAL_Barrel);
	  m_recoPh_E_HE->push_back(En_HCAL_Endcap);
	  m_recoPh_E_HO->push_back(En_HCAL_else);
	  m_recoPh_firstLayerECAL->push_back(phFirstLayerECAL);
	  m_recoPh_lastLayerECAL->push_back(phLastLayerECAL);
	  m_recoPh_nhitsEB->push_back(hitsEB);
	  m_recoPh_nhitsEE->push_back(hitsEE);
	  m_recoPh_nhitsEO->push_back(hitsEO);
	  m_recoPh_firstLayerHCAL->push_back(phFirstLayerHCAL);
	  m_recoPh_lastLayerHCAL->push_back(phLastLayerHCAL);
	  m_recoPh_nhitsHB->push_back(hitsHB);
	  m_recoPh_nhitsHE->push_back(hitsHE);
	  m_recoPh_nhitsHO->push_back(hitsHO);
	}
      }
    }

    LCCollection * jetExclColl =0;
    //std::cout<<"before jetColl"<<std::endl;
    getCollection(jetExclColl,m_jetExclColName,evt);
    //std::cout<<"aftrer jet exclColl"<<std::endl;
    if(jetExclColl!=NULL){
      //std::cout<<"in jet ecxclColl"<<std::endl;
      if(jetExclColl->getNumberOfElements()!=2){
	std::cout<<"expected number of events always to be 2 "<<jetExclColl->getNumberOfElements()<<std::endl;
      }
      for(int i=0;i<jetExclColl->getNumberOfElements();i++){
	ReconstructedParticle* jetExcl = dynamic_cast<ReconstructedParticle*>(jetExclColl->getElementAt(i));
	m_jet_excl_Px->push_back(jetExcl->getMomentum()[0]); 
	m_jet_excl_Py->push_back(jetExcl->getMomentum()[1]); 
	m_jet_excl_Pz->push_back(jetExcl->getMomentum()[2]); 
	m_jet_excl_E->push_back(jetExcl->getEnergy()); 
	m_jet_excl_Phi->push_back(atan2(jetExcl->getMomentum()[1],jetExcl->getMomentum()[0])); 
	double cosTheta = jetExcl->getMomentum()[2]/sqrt(jetExcl->getMomentum()[0]*jetExcl->getMomentum()[0]+jetExcl->getMomentum()[1]*jetExcl->getMomentum()[1]+jetExcl->getMomentum()[2]*jetExcl->getMomentum()[2]);
	m_jet_excl_CosTheta->push_back(cosTheta);
	float jetExclPiE=0;
	float jetExclPhE=0;
	float jetExclElE=0;
	float jetExclMuE=0;
	float jetExclNE=0;
	float jetExclNeutElseE=0;
	float jetExclChElseE=0;
	int jetExclNeutMult=0;
	int jetExclChMult=0;
	for(unsigned int j=0;j<jetExcl->getParticles().size();j++){
	  if(jetExcl->getParticles()[j]->getCharge()==0){
	    jetExclNeutMult+=1;
	  }else{
	    jetExclChMult+=1;
	  }
	  if(abs(jetExcl->getParticles()[j]->getType())==211){ 
	    jetExclPiE+=jetExcl->getParticles()[j]->getEnergy();
	  }else if(jetExcl->getParticles()[j]->getType()==22){ 
	    jetExclPhE+=jetExcl->getParticles()[j]->getEnergy();
	  }else if(abs(jetExcl->getParticles()[j]->getType())==11){ 
	    jetExclElE+=jetExcl->getParticles()[j]->getEnergy();
	  }else if(abs(jetExcl->getParticles()[j]->getType())==13){ 
	    jetExclMuE+=jetExcl->getParticles()[j]->getEnergy();
	  }else if(jetExcl->getParticles()[j]->getType()==2112){ 
	    jetExclNE+=jetExcl->getParticles()[j]->getEnergy();
	  }else if(jetExcl->getParticles()[j]->getCharge()==0){
	    jetExclNeutElseE+=jetExcl->getParticles()[j]->getEnergy();
	  }else{
	    jetExclChElseE+=jetExcl->getParticles()[j]->getEnergy();
	  }
	}
	m_jet_excl_piE->push_back(jetExclPiE);
	m_jet_excl_phE->push_back(jetExclPhE);
	m_jet_excl_elE->push_back(jetExclElE);
	m_jet_excl_muE->push_back(jetExclMuE);
	m_jet_excl_nE->push_back(jetExclNE);
	m_jet_excl_neutElseE->push_back(jetExclNeutElseE);
	m_jet_excl_chElseE->push_back(jetExclChElseE);
	m_jet_excl_neutMult->push_back(jetExclNeutMult);
	m_jet_excl_chMult->push_back(jetExclChMult);
      }
    }

    LCCollection * jetInclColl =0;
    
    //std::cout<<"before jet incColl"<<std::endl;
    getCollection(jetInclColl,m_jetInclColName,evt);
    //std::cout<<"after jet incColl"<<std::endl;
    if(jetInclColl!=NULL){
      //std::cout<<"in jet incoll"<<std::endl;
      for(int i=0;i<jetInclColl->getNumberOfElements();i++){
	ReconstructedParticle* jetIncl = dynamic_cast<ReconstructedParticle*>(jetInclColl->getElementAt(i));
	if(jetIncl->getEnergy()<m_jetEMin){
	  std::cout<<"I asssumed we should never get here "<<jetIncl->getEnergy()<<"/"<<m_jetEMin<<std::endl;
	  continue;
	}
	m_jet_incl_Px->push_back(jetIncl->getMomentum()[0]); 
	m_jet_incl_Py->push_back(jetIncl->getMomentum()[1]); 
	m_jet_incl_Pz->push_back(jetIncl->getMomentum()[2]); 
	m_jet_incl_E->push_back(jetIncl->getEnergy()); 
	m_jet_incl_Phi->push_back(atan2(jetIncl->getMomentum()[1],jetIncl->getMomentum()[0])); 
	double cosTheta = jetIncl->getMomentum()[2]/sqrt(jetIncl->getMomentum()[0]*jetIncl->getMomentum()[0]+jetIncl->getMomentum()[1]*jetIncl->getMomentum()[1]+jetIncl->getMomentum()[2]*jetIncl->getMomentum()[2]);
	m_jet_incl_CosTheta->push_back(cosTheta);
	float jetInclPiE=0;
	float jetInclPhE=0;
	float jetInclElE=0;
	float jetInclMuE=0;
	float jetInclNE=0;
	float jetInclNeutElseE=0;
	float jetInclChElseE=0;
	int jetInclNeutMult=0;
	int jetInclChMult=0;
	for(unsigned int j=0;j<jetIncl->getParticles().size();j++){
	  if(jetIncl->getParticles()[j]->getCharge()==0){
	    jetInclNeutMult+=1;
	  }else{
	    jetInclChMult+=1;
	  }
	  if(abs(jetIncl->getParticles()[j]->getType())==211){ 
	    jetInclPiE+=jetIncl->getParticles()[j]->getEnergy();
	  }else if(jetIncl->getParticles()[j]->getType()==22){ 
	    jetInclPhE+=jetIncl->getParticles()[j]->getEnergy();
	  }else if(abs(jetIncl->getParticles()[j]->getType())==11){ 
	    jetInclElE+=jetIncl->getParticles()[j]->getEnergy();
	  }else if(abs(jetIncl->getParticles()[j]->getType())==13){ 
	    jetInclMuE+=jetIncl->getParticles()[j]->getEnergy();
	  }else if(jetIncl->getParticles()[j]->getType()==2112){ 
	    jetInclNE+=jetIncl->getParticles()[j]->getEnergy();
	  }else if(jetIncl->getParticles()[j]->getCharge()==0){
	    jetInclNeutElseE+=jetIncl->getParticles()[j]->getEnergy();
	  }else{
	    jetInclChElseE+=jetIncl->getParticles()[j]->getEnergy();
	  }
	}
	m_jet_incl_piE->push_back(jetInclPiE);
	m_jet_incl_phE->push_back(jetInclPhE);
	m_jet_incl_elE->push_back(jetInclElE);
	m_jet_incl_muE->push_back(jetInclMuE);
	m_jet_incl_nE->push_back(jetInclNE);
	m_jet_incl_neutElseE->push_back(jetInclNeutElseE);
	m_jet_incl_chElseE->push_back(jetInclChElseE);
	m_jet_incl_neutMult->push_back(jetInclNeutMult);
	m_jet_incl_chMult->push_back(jetInclChMult);
      }
    }

    if((m_E_trueNeut+m_E_trueZ2)>3050){
      std::cout<<"high E sum for true stuff "<<(m_E_trueNeut+m_E_trueZ2)<<"/"<<m_E_trueNeut<<"/"<<m_E_trueZ2<<std::endl;
    }
    
    //std::cout<<"before tree filling"<<std::endl;
    m_outputTree->Fill();
    //std::cout<<"after tree filling"<<std::endl;
}

void PhotonStudy::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

void PhotonStudy::check(LCEvent*){
}

void PhotonStudy::end(){
    
  m_rootFile->cd();  
 
  m_rootFile->Write();
  m_rootFile->Close();

}
