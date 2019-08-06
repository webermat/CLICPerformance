#include "ZudsHistos.h"
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


ZudsHistos aZudsHistos;

ZudsHistos::ZudsHistos() : Processor("ZudsHistos") {
    
    // modify processor description
    _description = "ZudsHistos calculates properties of calorimeter showers" ;
   

    registerInputCollection( LCIO::MCPARTICLE,
                            "MCParticleCollectionName",
                            "Name of the MCParticle input collection",
                            m_inputMCParticleCollection,
                            std::string("MCParticle"));
    
    registerProcessorParameter( "OutputRootFileName",
                                "ROOT File name to collect plots",
                                m_rootFileName,
                                std::string("ZudsHistos.root"));

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


void ZudsHistos::init() {

  m_rootFile = new TFile(m_rootFileName.c_str(),"RECREATE");
  /*
  h_truePh_E_0_2_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_0_2_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_0_2_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_0_2_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_0_2_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_0_2_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_0_2_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_0_2_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_0_2_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_2_5_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_2_5_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_2_5_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_2_5_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_2_5_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_2_5_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_2_5_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_2_5_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_2_5_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_5_10_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_5_10_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_5_10_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_5_10_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_5_10_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_5_10_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_5_10_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_5_10_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_5_10_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_10_15_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_10_15_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_10_15_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_10_15_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_10_15_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_10_15_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_10_15_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_10_15_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_10_15_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_15_25_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_15_25_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_15_25_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_15_25_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_15_25_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_15_25_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_15_25_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_15_25_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_15_25_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_25_50_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_25_50_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_25_50_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_25_50_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_25_50_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_25_50_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_25_50_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_25_50_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,15_75,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_N_E_Closest_vs_CosTheta","",100,15_75,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_CH_E_Closest_vs_CosTheta","",100,15_75,100,-1.0,1.0);

  h_truePh_E_25_50_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,15_75,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_N_E_Closest_vs_CosTheta","",100,15_75,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_CH_E_Closest_vs_CosTheta","",100,15_75,100,-1.0,1.0);

  h_truePh_E_50_100_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_50_100_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_50_100_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_50_100_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_50_100_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_50_100_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_50_100_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_50_100_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,30,125,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_N_E_Closest_vs_CosTheta","",100,30,125,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_CH_E_Closest_vs_CosTheta","",100,30,125,100,-1.0,1.0);

  h_truePh_E_50_100_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,30,125,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_N_E_Closest_vs_CosTheta","",100,30,125,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_CH_E_Closest_vs_CosTheta","",100,30,125,100,-1.0,1.0);

  h_truePh_E_100_Inf_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_100_Inf_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_100_Inf_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_100_Inf_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_100_Inf_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_100_Inf_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_100_Inf_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_100_Inf_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_Ph_E_Closest_vs_CosTheta","",150,75,400,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_N_E_Closest_vs_CosTheta","",150,75,400,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_CH_E_Closest_vs_CosTheta","",150,75,400,100,-1.0,1.0);

  h_truePh_E_100_Inf_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_Ph_E_Closest_vs_CosTheta","",150,75,400,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_N_E_Closest_vs_CosTheta","",150,75,400,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_CH_E_Closest_vs_CosTheta","",150,75,400,100,-1.0,1.0);
  */


  h_truePi_E  = new TH1D("h_truePi_E","",1500,0,1500);
  h_truePr_E  = new TH1D("h_truePr_E","",1500,0,1500);
  h_trueCH_E  = new TH1D("h_trueCH_E","",1500,0,1500);
  h_truePh_E  = new TH1D("h_truePh_E","",1500,0,1500);
  h_trueNe_E  = new TH1D("h_trueNe_E","",1500,0,1500);
  h_trueK0L_E = new TH1D("h_trueK0L_E","",1500,0,1500);
  h_trueNH_E  = new TH1D("h_trueNH_E","",1500,0,1500);
  h_trueEl_E  = new TH1D("h_trueEl_E","",1500,0,1500);
  h_trueMu_E  = new TH1D("h_trueMu_E","",1500,0,1500);

  h_truePi_Pt  = new TH1D("h_truePi_Pt","",1500,0,1500);
  h_truePr_Pt  = new TH1D("h_truePr_Pt","",1500,0,1500);
  h_trueCH_Pt  = new TH1D("h_trueCH_Pt","",1500,0,1500);

  h_recoPi_E  = new TH1D("h_recoPi_E","",1500,0,1500);
  h_recoNe_E  = new TH1D("h_recoNe_E","",1500,0,1500);
  h_recoPh_E  = new TH1D("h_recoPh_E","",1500,0,1500);
  h_recoEl_E  = new TH1D("h_recoEl_E","",1500,0,1500);
  h_recoMu_E  = new TH1D("h_recoMu_E","",1500,0,1500);

  h_recoPi_Pt  = new TH1D("h_recoPi_Pt","",1500,0,1500);
 

    eventcount=0;
  
    // Print the initial parameters
    printParameters() ;

    // Reset counters
    m_runNumber = 0 ;
    m_eventNumber = 0 ;

}


void ZudsHistos::processRunHeader( LCRunHeader*) {
	++m_runNumber ;
}

void ZudsHistos::processEvent( LCEvent* evt ) {
    
  eventcount+=1;
  if(evt->getEventNumber()%50==0){
    std::cout<<"run/evt "<<evt->getRunNumber()<<"/"<<evt->getEventNumber()<<"/"<<eventcount<<std::endl;
  }


  m_eventNumber=evt->getEventNumber();
  m_runNumber=evt->getRunNumber();
  m_eventCount=eventcount;

    LCCollection * mcColl =0;
    getCollection(mcColl,m_inputMCParticleCollection,evt);

    //Look for a photon
    for(int m =0; m< mcColl->getNumberOfElements(); m++){
      MCParticle* mcp= dynamic_cast<MCParticle*>( mcColl->getElementAt(m) ) ;
      if(mcp->getGeneratorStatus()==1){//visible sum of stable particles --> take neutrinos out
	double cos_theta=mcp->getMomentum()[2]/sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]+mcp->getMomentum()[2]*mcp->getMomentum()[2]);
	if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
	  if(abs(mcp->getPDG())==211){
	    h_truePi_E->Fill(mcp->getEnergy());
	    h_truePi_Pt->Fill(sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]));
	    h_trueCH_E->Fill(mcp->getEnergy());
	    h_trueCH_Pt->Fill(sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]));
	  }else if (mcp->getPDG()==22){
	    h_truePh_E->Fill(mcp->getEnergy());


	    LCCollection* recoparticlecolmatch = NULL;
	    // Alternativelly if you do not want Marlin to exit in case of a non-existing collection
	    // use the following (commented out) code:
	    //run on H to gamma gamma -> in case there are no tracks around linker will fail to produce output collection
	    recoparticlecolmatch = evt->getCollection(m_inputRECOParticleCollection) ;
	    if(recoparticlecolmatch!=NULL){
	      //PandoraCandidate loop
	      float angle_min=7.5;
	      float index_min=-1;;
	      for(int i=0;i<recoparticlecolmatch->getNumberOfElements();i++){
		ReconstructedParticle* pandorapart = dynamic_cast<ReconstructedParticle*>(recoparticlecolmatch->getElementAt(i));
		float angle = acos((mcp->getMomentum()[0]*pandorapart->getMomentum()[0]+mcp->getMomentum()[1]*pandorapart->getMomentum()[1]+mcp->getMomentum()[2]*pandorapart->getMomentum()[2])/(sqrt(pandorapart->getMomentum()[0]*pandorapart->getMomentum()[0]+pandorapart->getMomentum()[1]*pandorapart->getMomentum()[1]+pandorapart->getMomentum()[2]*pandorapart->getMomentum()[2])*sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]+mcp->getMomentum()[2]*mcp->getMomentum()[2])));
		if(angle<angle_min){
		  angle_min=angle;
		  index_min=i;
		}
	      }
	    }
	    /*
 h_truePh_E_0_2_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_0_2_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_0_2_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_0_2_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_0_2_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_0_2_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_0_2_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_0_2_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_0_2_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_0_2_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_0_2_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_2_5_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_2_5_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_2_5_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_2_5_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_2_5_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_2_5_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_2_5_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_2_5_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_2_5_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_2_5_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_2_5_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_5_10_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_5_10_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_5_10_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_5_10_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_5_10_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_5_10_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_5_10_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_5_10_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_5_10_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_5_10_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_5_10_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_10_15_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_10_15_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_10_15_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_10_15_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_10_15_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_10_15_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_10_15_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_10_15_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_10_15_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_10_15_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_10_15_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_15_25_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_15_25_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_15_25_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_15_25_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_15_25_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_15_25_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_15_25_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_15_25_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_15_25_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_N_E_Closest_vs_CosTheta","",100,0,10,-1.0,1.0);
  h_truePh_E_15_25_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_15_25_DA_0_01_CH_E_Closest_vs_CosTheta","",100,0,10,100,-1.0,1.0);

  h_truePh_E_25_50_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_25_50_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_25_50_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_25_50_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_25_50_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_25_50_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_25_50_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_25_50_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,15_75,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_N_E_Closest_vs_CosTheta","",100,15_75,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_CH_E_Closest_vs_CosTheta","",100,15_75,100,-1.0,1.0);

  h_truePh_E_25_50_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,15_75,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_N_E_Closest_vs_CosTheta","",100,15_75,-1.0,1.0);
  h_truePh_E_25_50_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_25_50_DA_0_01_CH_E_Closest_vs_CosTheta","",100,15_75,100,-1.0,1.0);

  h_truePh_E_50_100_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_50_100_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_50_100_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_50_100_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_50_100_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_50_100_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_50_100_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_50_100_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,30,125,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_N_E_Closest_vs_CosTheta","",100,30,125,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_CH_E_Closest_vs_CosTheta","",100,30,125,100,-1.0,1.0);

  h_truePh_E_50_100_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_Ph_E_Closest_vs_CosTheta","",100,30,125,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_N_E_Closest_vs_CosTheta","",100,30,125,-1.0,1.0);
  h_truePh_E_50_100_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_50_100_DA_0_01_CH_E_Closest_vs_CosTheta","",100,30,125,100,-1.0,1.0);

  h_truePh_E_100_Inf_DA_0_01_N_part_vs_CosTheta  = new TH1D("h_truePh_E_100_Inf_DA_0_01_N_part_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_Ph_closest_vs_CosTheta  = new TH1D("h_truePh_E_100_Inf_DA_0_01_Ph_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_N_closest_vs_CosTheta  = new TH1D("h_truePh_E_100_Inf_DA_0_01_N_Closest_vs_CosTheta","",100,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_CH_closest_vs_CosTheta  = new TH1D("h_truePh_E_100_Inf_DA_0_01_CH_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_100_Inf_DA_0_01_NoPart_closest_vs_CosTheta  = new TH1D("h_truePh_E_100_Inf_DA_0_01_NoPart_Closest_vs_CosTheta","",100,-1.0,1.0);

  h_truePh_E_100_Inf_DA_0_01_Ph_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_Ph_E_Closest_vs_CosTheta","",150,75,400,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_N_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_N_E_Closest_vs_CosTheta","",150,75,400,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_CH_E_closest_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_CH_E_Closest_vs_CosTheta","",150,75,400,100,-1.0,1.0);

  h_truePh_E_100_Inf_DA_0_01_Ph_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_Ph_E_Closest_vs_CosTheta","",150,75,400,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_N_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_N_E_Closest_vs_CosTheta","",150,75,400,-1.0,1.0);
  h_truePh_E_100_Inf_DA_0_01_CH_E_closestE_vs_CosTheta  = new TH2D("h_truePh_E_100_Inf_DA_0_01_CH_E_Closest_vs_CosTheta","",150,75,400,100,-1.0,1.0);
	    */
	  }else if (abs(mcp->getPDG())==130){
	    h_trueK0L_E->Fill(mcp->getEnergy());
	    h_trueNH_E->Fill(mcp->getEnergy());
	  }else if(abs(mcp->getPDG())==11){
	    h_trueEl_E->Fill(mcp->getEnergy());
	  }else if(abs(mcp->getPDG())==13 ){
	    h_trueMu_E->Fill(mcp->getEnergy());
	  }else if(abs(mcp->getPDG())==2112 ){
	    h_trueNe_E->Fill(mcp->getEnergy());
	    h_trueNH_E->Fill(mcp->getEnergy());
	  }else if(abs(mcp->getPDG())==2212 ){
	    h_truePr_E->Fill(mcp->getEnergy());
	    h_truePr_Pt->Fill(sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]));
	    h_trueCH_E->Fill(mcp->getEnergy());
	    h_trueCH_Pt->Fill(sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]));
	  }else if(abs(mcp->getPDG())==321 || mcp->getCharge()==0 ){
	    h_trueNH_E->Fill(mcp->getEnergy());
	  }else{
	    h_trueCH_E->Fill(mcp->getEnergy());
	    h_trueCH_Pt->Fill(sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]));
	  }
	}
      }
    }
    
    //now reco PFO loop
    LCCollection* recoparticlecol = NULL;
    // Alternativelly if you do not want Marlin to exit in case of a non-existing collection
    // use the following (commented out) code:
    //run on H to gamma gamma -> in case there are no tracks around linker will fail to produce output collection
    recoparticlecol = evt->getCollection(m_inputRECOParticleCollection) ;
    if(recoparticlecol!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<recoparticlecol->getNumberOfElements();i++){
	ReconstructedParticle* pandorapart = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(i));
	double cosTheta = pandorapart->getMomentum()[2]/sqrt(pandorapart->getMomentum()[0]*pandorapart->getMomentum()[0]+pandorapart->getMomentum()[1]*pandorapart->getMomentum()[1]+pandorapart->getMomentum()[2]*pandorapart->getMomentum()[2]);
	if(abs(pandorapart->getType())==211){
	  h_recoPi_E->Fill(pandorapart->getEnergy());
	  h_recoPi_Pt->Fill(sqrt(pandorapart->getMomentum()[0]*pandorapart->getMomentum()[0]+pandorapart->getMomentum()[1]*pandorapart->getMomentum()[1]));
	}else if (pandorapart->getType()==22){
	  h_recoPh_E->Fill(pandorapart->getEnergy());
	}else if(abs(pandorapart->getType())==11){
	  h_recoEl_E->Fill(pandorapart->getEnergy());
	}else if(abs(pandorapart->getType())==13){
	  h_recoMu_E->Fill(pandorapart->getEnergy());
	}else if(pandorapart->getType()==2112){
	  h_recoNe_E->Fill(pandorapart->getEnergy());
	}
      }
    }	
}

void ZudsHistos::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    
    streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

void ZudsHistos::check(LCEvent*){
}

void ZudsHistos::end(){
    
  m_rootFile->cd();  
 
  m_rootFile->Write();
  m_rootFile->Close();

}
