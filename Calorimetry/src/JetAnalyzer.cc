#include "JetAnalyzer.h"
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


JetAnalyzer aJetAnalyzer;

JetAnalyzer::JetAnalyzer() : Processor("JetAnalyzer") {
    
    // modify processor description
    _description = "JetAnalyzer calculates properties of calorimeter showers" ;
   

    registerInputCollection( LCIO::MCPARTICLE,
                            "MCParticleCollectionName",
                            "Name of the MCParticle input collection",
                            m_inputMCParticleCollection,
                            std::string("MCPhysicsParticles"));
    
    registerProcessorParameter( "OutputRootFileName",
                                "ROOT File name to collect plots",
                                m_rootFileName,
                                std::string("JetAnalyzer.root"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RECOParticleCollectionName",
                            "Name of the RECOParticle input collection",
                            m_inputRECOParticleCollection,
                            std::string("PandoraPFOs"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "GenJetDR03Collection" ,  
			     "Name of the DR03 GenJet collection"  ,
			     m_genjetDR03ColName,
			     std::string("GenJet_VLC_R03")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "GenJetDR04Collection" ,  
			     "Name of the DR04 GenJet collection"  ,
			     m_genjetDR04ColName,
			     std::string("GenJet_VLC_R04")
			     );
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "GenJetDR05Collection" ,  
			     "Name of the DR05 GenJet collection"  ,
			     m_genjetDR05ColName,
			     std::string("GenJet_VLC_R05")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "GenJetDR06Collection" ,  
			     "Name of the DR06 GenJet collection"  ,
			     m_genjetDR06ColName,
			     std::string("GenJet_VLC_R06")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "GenJetDR07Collection" ,  
			     "Name of the DR07 GenJet collection"  ,
			     m_genjetDR07ColName,
			     std::string("GenJet_VLC_R07")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "GenJetDR08Collection" ,  
			     "Name of the DR08 GenJet collection"  ,
			     m_genjetDR08ColName,
			     std::string("GenJet_VLC_R08")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "GenJetDR09Collection" ,  
			     "Name of the DR09 GenJet collection"  ,
			     m_genjetDR09ColName,
			     std::string("GenJet_VLC_R09")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "GenJetDR10Collection" ,  
			     "Name of the DR10 GenJet collection"  ,
			     m_genjetDR10ColName,
			     std::string("GenJet_VLC_R10")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetDR03Collection" ,  
			     "Name of the DR03 RecoJet collection"  ,
			     m_recojetDR03ColName,
			     std::string("RecoJet_VLC_R03")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetDR04Collection" ,  
			     "Name of the DR04 RecoJet collection"  ,
			     m_recojetDR04ColName,
			     std::string("RecoJet_VLC_R04")
			     );
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetDR05Collection" ,  
			     "Name of the DR05 RecoJet collection"  ,
			     m_recojetDR05ColName,
			     std::string("RecoJet_VLC_R05")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetDR06Collection" ,  
			     "Name of the DR06 RecoJet collection"  ,
			     m_recojetDR06ColName,
			     std::string("RecoJet_VLC_R06")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetDR07Collection" ,  
			     "Name of the DR07 RecoJet collection"  ,
			     m_recojetDR07ColName,
			     std::string("RecoJet_VLC_R07")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetDR08Collection" ,  
			     "Name of the DR08 RecoJet collection"  ,
			     m_recojetDR08ColName,
			     std::string("RecoJet_VLC_R08")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetDR09Collection" ,  
			     "Name of the DR09 RecoJet collection"  ,
			     m_recojetDR09ColName,
			     std::string("RecoJet_VLC_R09")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetDR10Collection" ,  
			     "Name of the DR10 RecoJet collection"  ,
			     m_recojetDR10ColName,
			     std::string("RecoJet_VLC_R10")
			     );

    registerProcessorParameter(
			       "fillMEInfo" , 
			       "save matrix element information",
			       m_fillMEInfo,
			       bool(false)			       
			       );

    registerProcessorParameter(
			       "fillAllJets" , 
			       "save all jets, or only one selected chosen as R=0.7",
			       m_fillAllJets,
			       bool(true)			       
			       );

}  


void JetAnalyzer::init() {

  m_rootFile = new TFile(m_rootFileName.c_str(),"RECREATE");

  m_outputTree = new TTree("showerData","showerData");

  eventcount=0;
  
  // Print the initial parameters
  printParameters() ;
  
  // Reset counters
  m_runNumber = 0 ;
  m_eventNumber = 0 ;

  if(m_fillMEInfo){
    m_trueME_E=new std::vector<float>();
    m_trueME_Px=new std::vector<float>();
    m_trueME_Py=new std::vector<float>();
    m_trueME_Pz=new std::vector<float>();
    m_trueME_PDGID=new std::vector<int>();
  }

  m_genJetDR03_E   = new std::vector<float>(); 
  m_genJetDR03_Px  = new std::vector<float>(); 
  m_genJetDR03_Py  = new std::vector<float>(); 
  m_genJetDR03_Pz  = new std::vector<float>(); 
  
  m_genJetDR04_E   = new std::vector<float>(); 
  m_genJetDR04_Px  = new std::vector<float>(); 
  m_genJetDR04_Py  = new std::vector<float>(); 
  m_genJetDR04_Pz  = new std::vector<float>(); 
  
  m_genJetDR05_E   = new std::vector<float>(); 
  m_genJetDR05_Px  = new std::vector<float>(); 
  m_genJetDR05_Py  = new std::vector<float>(); 
  m_genJetDR05_Pz  = new std::vector<float>();
  
  m_genJetDR06_E   = new std::vector<float>(); 
  m_genJetDR06_Px  = new std::vector<float>(); 
  m_genJetDR06_Py  = new std::vector<float>(); 
  m_genJetDR06_Pz  = new std::vector<float>(); 
  
  m_genJetDR08_E   = new std::vector<float>(); 
  m_genJetDR08_Px  = new std::vector<float>(); 
  m_genJetDR08_Py  = new std::vector<float>(); 
  m_genJetDR08_Pz  = new std::vector<float>(); 

  m_genJetDR09_E   = new std::vector<float>(); 
  m_genJetDR09_Px  = new std::vector<float>(); 
  m_genJetDR09_Py  = new std::vector<float>(); 
  m_genJetDR09_Pz  = new std::vector<float>(); 

  m_genJetDR10_E   = new std::vector<float>(); 
  m_genJetDR10_Px  = new std::vector<float>(); 
  m_genJetDR10_Py  = new std::vector<float>(); 
  m_genJetDR10_Pz  = new std::vector<float>(); 

  m_recoJetDR03_E   = new std::vector<float>(); 
  m_recoJetDR03_Px  = new std::vector<float>(); 
  m_recoJetDR03_Py  = new std::vector<float>(); 
  m_recoJetDR03_Pz  = new std::vector<float>(); 

  m_recoJetDR04_E   = new std::vector<float>(); 
  m_recoJetDR04_Px  = new std::vector<float>(); 
  m_recoJetDR04_Py  = new std::vector<float>(); 
  m_recoJetDR04_Pz  = new std::vector<float>(); 

  m_recoJetDR05_E   = new std::vector<float>(); 
  m_recoJetDR05_Px  = new std::vector<float>(); 
  m_recoJetDR05_Py  = new std::vector<float>(); 
  m_recoJetDR05_Pz  = new std::vector<float>();
 
  m_recoJetDR06_Pz  = new std::vector<float>(); 
  m_recoJetDR06_E   = new std::vector<float>(); 
  m_recoJetDR06_Px  = new std::vector<float>(); 
  m_recoJetDR06_Py  = new std::vector<float>(); 

  m_recoJetDR08_E   = new std::vector<float>(); 
  m_recoJetDR08_Px  = new std::vector<float>(); 
  m_recoJetDR08_Py  = new std::vector<float>(); 
  m_recoJetDR08_Pz  = new std::vector<float>(); 

  m_recoJetDR09_E   = new std::vector<float>(); 
  m_recoJetDR09_Px  = new std::vector<float>(); 
  m_recoJetDR09_Py  = new std::vector<float>(); 
  m_recoJetDR09_Pz  = new std::vector<float>(); 

  m_recoJetDR10_E   = new std::vector<float>(); 
  m_recoJetDR10_Px  = new std::vector<float>(); 
  m_recoJetDR10_Py  = new std::vector<float>(); 
  m_recoJetDR10_Pz  = new std::vector<float>(); 

  m_genJetDR07_E   = new std::vector<float>(); 
  m_genJetDR07_Px  = new std::vector<float>(); 
  m_genJetDR07_Py  = new std::vector<float>(); 
  m_genJetDR07_Pz  = new std::vector<float>(); 

  m_recoJetDR07_E   = new std::vector<float>(); 
  m_recoJetDR07_Px  = new std::vector<float>(); 
  m_recoJetDR07_Py  = new std::vector<float>(); 
  m_recoJetDR07_Pz  = new std::vector<float>(); 
 


  m_E_trueInv  = 0;
  m_px_trueInv = 0;
  m_py_trueInv = 0;
  m_pz_trueInv = 0;

  m_E_trueAll  = 0;
  m_px_trueAll = 0;
  m_py_trueAll = 0;
  m_pz_trueAll = 0;

  m_E_totPFO  = 0;
  m_px_totPFO = 0;
  m_py_totPFO = 0;
  m_pz_totPFO = 0;

  m_E_trueInv_cos09=0;
  m_px_trueInv_cos09=0;
  m_py_trueInv_cos09=0;
  m_pz_trueInv_cos09=0;
  m_E_trueAll_cos09=0;
  m_px_trueAll_cos09=0;
  m_py_trueAll_cos09=0;
  m_pz_trueAll_cos09=0;
  m_E_totPFO_cos09=0;
  m_px_totPFO_cos09=0;
  m_py_totPFO_cos09=0;
  m_pz_totPFO_cos09=0;

  m_d1_mcPDGID = 0;
  m_d1_mcE     = 0;
  m_d1_mcPx    = 0;
  m_d1_mcPy    = 0;
  m_d1_mcPz    = 0;
  
  m_d2_mcPDGID = 0;
  m_d2_mcE     = 0;
  m_d2_mcPx    = 0;
  m_d2_mcPy    = 0;
  m_d2_mcPz    = 0;

  m_genJetDR03_E ->clear(); 
  m_genJetDR03_Px ->clear(); 
  m_genJetDR03_Py ->clear(); 
  m_genJetDR03_Pz ->clear(); 

  m_genJetDR04_E ->clear(); 
  m_genJetDR04_Px ->clear(); 
  m_genJetDR04_Py ->clear(); 
  m_genJetDR04_Pz ->clear(); 

  m_genJetDR05_E ->clear(); 
  m_genJetDR05_Px ->clear(); 
  m_genJetDR05_Py ->clear(); 
  m_genJetDR05_Pz ->clear(); 

  m_genJetDR06_E ->clear(); 
  m_genJetDR06_Px ->clear(); 
  m_genJetDR06_Py ->clear(); 
  m_genJetDR06_Pz ->clear(); 

  m_genJetDR07_E ->clear(); 
  m_genJetDR07_Px ->clear(); 
  m_genJetDR07_Py ->clear(); 
  m_genJetDR07_Pz ->clear(); 

  m_genJetDR08_E ->clear(); 
  m_genJetDR08_Px ->clear(); 
  m_genJetDR08_Py ->clear(); 
  m_genJetDR08_Pz ->clear(); 

  m_genJetDR09_E ->clear(); 
  m_genJetDR09_Px ->clear(); 
  m_genJetDR09_Py ->clear(); 
  m_genJetDR09_Pz ->clear(); 

  m_genJetDR10_E ->clear(); 
  m_genJetDR10_Px ->clear(); 
  m_genJetDR10_Py ->clear(); 
  m_genJetDR10_Pz ->clear(); 

  m_recoJetDR03_E ->clear(); 
  m_recoJetDR03_Px ->clear(); 
  m_recoJetDR03_Py ->clear(); 
  m_recoJetDR03_Pz ->clear(); 

  m_recoJetDR04_E ->clear(); 
  m_recoJetDR04_Px ->clear(); 
  m_recoJetDR04_Py ->clear(); 
  m_recoJetDR04_Pz ->clear(); 

  m_recoJetDR05_E ->clear(); 
  m_recoJetDR05_Px ->clear(); 
  m_recoJetDR05_Py ->clear(); 
  m_recoJetDR05_Pz ->clear(); 

  m_recoJetDR06_E ->clear(); 
  m_recoJetDR06_Px ->clear(); 
  m_recoJetDR06_Py ->clear(); 
  m_recoJetDR06_Pz ->clear(); 

  m_recoJetDR07_E ->clear(); 
  m_recoJetDR07_Px ->clear(); 
  m_recoJetDR07_Py ->clear(); 
  m_recoJetDR07_Pz ->clear(); 

  m_recoJetDR08_E ->clear(); 
  m_recoJetDR08_Px ->clear(); 
  m_recoJetDR08_Py ->clear(); 
  m_recoJetDR08_Pz ->clear(); 

  m_recoJetDR09_E ->clear(); 
  m_recoJetDR09_Px ->clear(); 
  m_recoJetDR09_Py ->clear(); 
  m_recoJetDR09_Pz ->clear(); 

  m_recoJetDR10_E ->clear(); 
  m_recoJetDR10_Px ->clear(); 
  m_recoJetDR10_Py ->clear(); 
  m_recoJetDR10_Pz ->clear(); 

  if(m_fillMEInfo){
    m_trueME_E->clear();
    m_trueME_Px->clear();
    m_trueME_Py->clear();
    m_trueME_Pz->clear();
    m_trueME_PDGID->clear();
  }

  if(m_fillMEInfo){
    m_outputTree->Branch("trueME_Px", "std::vector< float >", &m_trueME_Px); 
    m_outputTree->Branch("trueME_Py", "std::vector< float >", &m_trueME_Py); 
    m_outputTree->Branch("trueME_Pz", "std::vector< float >", &m_trueME_Pz); 
    m_outputTree->Branch("trueME_E", "std::vector< float >", &m_trueME_E); 
    m_outputTree->Branch("trueME_PDGID", "std::vector< int >", &m_trueME_PDGID); 
  }

  m_outputTree->Branch("d1_mcPDGID",&m_d1_mcPDGID,"d1_mcPDGID/I");
  m_outputTree->Branch("d1_mcE",&m_d1_mcE,"d1_mcE/F");
  m_outputTree->Branch("d1_mcPx",&m_d1_mcPx,"d1_mcPx/F");
  m_outputTree->Branch("d1_mcPy",&m_d1_mcPy,"d1_mcPy/F");
  m_outputTree->Branch("d1_mcPz",&m_d1_mcPz,"d1_mcPz/F");
  
  m_outputTree->Branch("d2_mcPDGID",&m_d2_mcPDGID,"d2_mcPDGID/I");
  m_outputTree->Branch("d2_mcE",&m_d2_mcE,"d2_mcE/F");
  m_outputTree->Branch("d2_mcPx",&m_d2_mcPx,"d2_mcPx/F");
  m_outputTree->Branch("d2_mcPy",&m_d2_mcPy,"d2_mcPy/F");
  m_outputTree->Branch("d2_mcPz",&m_d2_mcPz,"d2_mcPz/F");

  //true particle level, exclude neutrinos
  m_outputTree->Branch("E_trueAll" ,&m_E_trueAll, "E_trueAll/F");
  m_outputTree->Branch("Px_trueAll",&m_px_trueAll,"Px_trueAll/F");
  m_outputTree->Branch("Py_trueAll",&m_py_trueAll,"Py_trueAll/F");
  m_outputTree->Branch("Pz_trueAll",&m_pz_trueAll,"Pz_trueAll/F");

  //true particle level, only neutrinos
  m_outputTree->Branch("E_trueInv" ,&m_E_trueInv, "E_trueInv/F");
  m_outputTree->Branch("Px_trueInv",&m_px_trueInv,"Px_trueInv/F");
  m_outputTree->Branch("Py_trueInv",&m_py_trueInv,"Py_trueInv/F");
  m_outputTree->Branch("Pz_trueInv",&m_pz_trueInv,"Pz_trueInv/F");
  
  //reconstructed leve
  m_outputTree->Branch("E_totPFO" ,&m_E_totPFO, "E_totPFO/F");
  m_outputTree->Branch("Px_totPFO",&m_px_totPFO,"Px_totPFO/F");
  m_outputTree->Branch("Py_totPFO",&m_py_totPFO,"Py_totPFO/F");
  m_outputTree->Branch("Pz_totPFO",&m_pz_totPFO,"Pz_totPFO/F");
  //exclude very forward particles in order to avoid adding too much background from forward jets
  //true particle level, exclude neutrinos
  m_outputTree->Branch("E_trueAll_cos09" ,&m_E_trueAll_cos09, "E_trueAll_cos09/F");
  m_outputTree->Branch("Px_trueAll_cos09",&m_px_trueAll_cos09,"Px_trueAll_cos09/F");
  m_outputTree->Branch("Py_trueAll_cos09",&m_py_trueAll_cos09,"Py_trueAll_cos09/F");
  m_outputTree->Branch("Pz_trueAll_cos09",&m_pz_trueAll_cos09,"Pz_trueAll_cos09/F");

  //true particle level, only neutrinos
  m_outputTree->Branch("E_trueInv_cos09" ,&m_E_trueInv_cos09, "E_trueInv_cos09/F");
  m_outputTree->Branch("Px_trueInv_cos09",&m_px_trueInv_cos09,"Px_trueInv_cos09/F");
  m_outputTree->Branch("Py_trueInv_cos09",&m_py_trueInv_cos09,"Py_trueInv_cos09/F");
  m_outputTree->Branch("Pz_trueInv_cos09",&m_pz_trueInv_cos09,"Pz_trueInv_cos09/F");
  
  //reconstructed leve
  m_outputTree->Branch("E_totPFO_cos09" ,&m_E_totPFO_cos09, "E_totPFO_cos09/F");
  m_outputTree->Branch("Px_totPFO_cos09",&m_px_totPFO_cos09,"Px_totPFO_cos09/F");
  m_outputTree->Branch("Py_totPFO_cos09",&m_py_totPFO_cos09,"Py_totPFO_cos09/F");
  m_outputTree->Branch("Pz_totPFO_cos09",&m_pz_totPFO_cos09,"Pz_totPFO_cos09/F");

  if(m_fillAllJets){
    m_outputTree->Branch("genJetR03E",  "std::vector< float >", &m_genJetDR03_E); 
    m_outputTree->Branch("genJetR03Px", "std::vector< float >", &m_genJetDR03_Px); 
    m_outputTree->Branch("genJetR03Py", "std::vector< float >", &m_genJetDR03_Py); 
    m_outputTree->Branch("genJetR03Pz", "std::vector< float >", &m_genJetDR03_Pz); 
    
    m_outputTree->Branch("genJetR04E",  "std::vector< float >", &m_genJetDR04_E); 
    m_outputTree->Branch("genJetR04Px", "std::vector< float >", &m_genJetDR04_Px); 
    m_outputTree->Branch("genJetR04Py", "std::vector< float >", &m_genJetDR04_Py); 
    m_outputTree->Branch("genJetR04Pz", "std::vector< float >", &m_genJetDR04_Pz); 
    
    m_outputTree->Branch("genJetR05E",  "std::vector< float >", &m_genJetDR05_E); 
    m_outputTree->Branch("genJetR05Px", "std::vector< float >", &m_genJetDR05_Px); 
    m_outputTree->Branch("genJetR05Py", "std::vector< float >", &m_genJetDR05_Py); 
    m_outputTree->Branch("genJetR05Pz", "std::vector< float >", &m_genJetDR05_Pz); 
    
    m_outputTree->Branch("genJetR06E",  "std::vector< float >", &m_genJetDR06_E); 
    m_outputTree->Branch("genJetR06Px", "std::vector< float >", &m_genJetDR06_Px); 
    m_outputTree->Branch("genJetR06Py", "std::vector< float >", &m_genJetDR06_Py); 
    m_outputTree->Branch("genJetR06Pz", "std::vector< float >", &m_genJetDR06_Pz); 

    m_outputTree->Branch("genJetR08E",  "std::vector< float >", &m_genJetDR08_E); 
    m_outputTree->Branch("genJetR08Px", "std::vector< float >", &m_genJetDR08_Px); 
    m_outputTree->Branch("genJetR08Py", "std::vector< float >", &m_genJetDR08_Py); 
    m_outputTree->Branch("genJetR08Pz", "std::vector< float >", &m_genJetDR08_Pz); 
    
    m_outputTree->Branch("genJetR09E",  "std::vector< float >", &m_genJetDR09_E); 
    m_outputTree->Branch("genJetR09Px", "std::vector< float >", &m_genJetDR09_Px); 
    m_outputTree->Branch("genJetR09Py", "std::vector< float >", &m_genJetDR09_Py); 
    m_outputTree->Branch("genJetR09Pz", "std::vector< float >", &m_genJetDR09_Pz); 
    
    m_outputTree->Branch("genJetR10E",  "std::vector< float >", &m_genJetDR10_E); 
    m_outputTree->Branch("genJetR10Px", "std::vector< float >", &m_genJetDR10_Px); 
    m_outputTree->Branch("genJetR10Py", "std::vector< float >", &m_genJetDR10_Py); 
    m_outputTree->Branch("genJetR10Pz", "std::vector< float >", &m_genJetDR10_Pz); 
    
    m_outputTree->Branch("recoJetR03E",  "std::vector< float >", &m_recoJetDR03_E); 
    m_outputTree->Branch("recoJetR03Px", "std::vector< float >", &m_recoJetDR03_Px); 
    m_outputTree->Branch("recoJetR03Py", "std::vector< float >", &m_recoJetDR03_Py); 
    m_outputTree->Branch("recoJetR03Pz", "std::vector< float >", &m_recoJetDR03_Pz); 
    
    m_outputTree->Branch("recoJetR04E",  "std::vector< float >", &m_recoJetDR04_E); 
    m_outputTree->Branch("recoJetR04Px", "std::vector< float >", &m_recoJetDR04_Px); 
    m_outputTree->Branch("recoJetR04Py", "std::vector< float >", &m_recoJetDR04_Py); 
    m_outputTree->Branch("recoJetR04Pz", "std::vector< float >", &m_recoJetDR04_Pz); 
    
    m_outputTree->Branch("recoJetR05E",  "std::vector< float >", &m_recoJetDR05_E); 
    m_outputTree->Branch("recoJetR05Px", "std::vector< float >", &m_recoJetDR05_Px); 
    m_outputTree->Branch("recoJetR05Py", "std::vector< float >", &m_recoJetDR05_Py); 
    m_outputTree->Branch("recoJetR05Pz", "std::vector< float >", &m_recoJetDR05_Pz); 
    
    m_outputTree->Branch("recoJetR06E",  "std::vector< float >", &m_recoJetDR06_E); 
    m_outputTree->Branch("recoJetR06Px", "std::vector< float >", &m_recoJetDR06_Px); 
    m_outputTree->Branch("recoJetR06Py", "std::vector< float >", &m_recoJetDR06_Py); 
    m_outputTree->Branch("recoJetR06Pz", "std::vector< float >", &m_recoJetDR06_Pz); 

    m_outputTree->Branch("recoJetR08E",  "std::vector< float >", &m_recoJetDR08_E); 
    m_outputTree->Branch("recoJetR08Px", "std::vector< float >", &m_recoJetDR08_Px); 
    m_outputTree->Branch("recoJetR08Py", "std::vector< float >", &m_recoJetDR08_Py); 
    m_outputTree->Branch("recoJetR08Pz", "std::vector< float >", &m_recoJetDR08_Pz); 
    
    m_outputTree->Branch("recoJetR09E",  "std::vector< float >", &m_recoJetDR09_E); 
    m_outputTree->Branch("recoJetR09Px", "std::vector< float >", &m_recoJetDR09_Px); 
    m_outputTree->Branch("recoJetR09Py", "std::vector< float >", &m_recoJetDR09_Py); 
    m_outputTree->Branch("recoJetR09Pz", "std::vector< float >", &m_recoJetDR09_Pz); 
    
    m_outputTree->Branch("recoJetR10E",  "std::vector< float >", &m_recoJetDR10_E); 
    m_outputTree->Branch("recoJetR10Px", "std::vector< float >", &m_recoJetDR10_Px); 
    m_outputTree->Branch("recoJetR10Py", "std::vector< float >", &m_recoJetDR10_Py); 
    m_outputTree->Branch("recoJetR10Pz", "std::vector< float >", &m_recoJetDR10_Pz); 
  }

  m_outputTree->Branch("genJetR07E",  "std::vector< float >", &m_genJetDR07_E); 
  m_outputTree->Branch("genJetR07Px", "std::vector< float >", &m_genJetDR07_Px); 
  m_outputTree->Branch("genJetR07Py", "std::vector< float >", &m_genJetDR07_Py); 
  m_outputTree->Branch("genJetR07Pz", "std::vector< float >", &m_genJetDR07_Pz); 


  m_outputTree->Branch("recoJetR07E",  "std::vector< float >", &m_recoJetDR07_E); 
  m_outputTree->Branch("recoJetR07Px", "std::vector< float >", &m_recoJetDR07_Px); 
  m_outputTree->Branch("recoJetR07Py", "std::vector< float >", &m_recoJetDR07_Py); 
  m_outputTree->Branch("recoJetR07Pz", "std::vector< float >", &m_recoJetDR07_Pz); 


}


void JetAnalyzer::processRunHeader( LCRunHeader*) {
	++m_runNumber ;
}

void JetAnalyzer::processEvent( LCEvent* evt ) {
  eventcount+=1;
  if(evt->getEventNumber()%50==0){
    std::cout<<"run/evt "<<evt->getRunNumber()<<"/"<<evt->getEventNumber()<<"/"<<eventcount<<std::endl;
  }


  m_eventNumber=evt->getEventNumber();
  m_runNumber=evt->getRunNumber();
  m_eventCount=eventcount;

  m_E_trueInv_cos09=0;
  m_px_trueInv_cos09=0;
  m_py_trueInv_cos09=0;
  m_pz_trueInv_cos09=0;
  m_E_trueAll_cos09=0;
  m_px_trueAll_cos09=0;
  m_py_trueAll_cos09=0;
  m_pz_trueAll_cos09=0;
  m_E_totPFO_cos09=0;
  m_px_totPFO_cos09=0;
  m_py_totPFO_cos09=0;
  m_pz_totPFO_cos09=0;

  m_E_trueInv  = 0;
  m_px_trueInv = 0;
  m_py_trueInv = 0;
  m_pz_trueInv = 0;

  m_E_trueAll  = 0;
  m_px_trueAll = 0;
  m_py_trueAll = 0;
  m_pz_trueAll = 0;

  m_E_totPFO  = 0;
  m_px_totPFO = 0;
  m_py_totPFO = 0;
  m_pz_totPFO = 0;

  m_d1_mcPDGID=-10;
  m_d1_mcE=-10;
  m_d1_mcPx=-10;
  m_d1_mcPy=-10;
  m_d1_mcPz=-10;
  
  m_d2_mcPDGID=-10;
  m_d2_mcE=-10;
  m_d2_mcPx=-10;
  m_d2_mcPy=-10;
  m_d2_mcPz=-10;

  m_genJetDR03_E ->clear(); 
  m_genJetDR03_Px ->clear(); 
  m_genJetDR03_Py ->clear(); 
  m_genJetDR03_Pz ->clear(); 

  m_genJetDR04_E ->clear(); 
  m_genJetDR04_Px ->clear(); 
  m_genJetDR04_Py ->clear(); 
  m_genJetDR04_Pz ->clear(); 

  m_genJetDR05_E ->clear(); 
  m_genJetDR05_Px ->clear(); 
  m_genJetDR05_Py ->clear(); 
  m_genJetDR05_Pz ->clear(); 

  m_genJetDR06_E ->clear(); 
  m_genJetDR06_Px ->clear(); 
  m_genJetDR06_Py ->clear(); 
  m_genJetDR06_Pz ->clear(); 

  m_genJetDR07_E ->clear(); 
  m_genJetDR07_Px ->clear(); 
  m_genJetDR07_Py ->clear(); 
  m_genJetDR07_Pz ->clear(); 

  m_genJetDR08_E ->clear(); 
  m_genJetDR08_Px ->clear(); 
  m_genJetDR08_Py ->clear(); 
  m_genJetDR08_Pz ->clear(); 

  m_genJetDR09_E ->clear(); 
  m_genJetDR09_Px ->clear(); 
  m_genJetDR09_Py ->clear(); 
  m_genJetDR09_Pz ->clear(); 

  m_genJetDR10_E ->clear(); 
  m_genJetDR10_Px ->clear(); 
  m_genJetDR10_Py ->clear(); 
  m_genJetDR10_Pz ->clear(); 

  m_recoJetDR03_E ->clear(); 
  m_recoJetDR03_Px ->clear(); 
  m_recoJetDR03_Py ->clear(); 
  m_recoJetDR03_Pz ->clear(); 

  m_recoJetDR04_E ->clear(); 
  m_recoJetDR04_Px ->clear(); 
  m_recoJetDR04_Py ->clear(); 
  m_recoJetDR04_Pz ->clear(); 

  m_recoJetDR05_E ->clear(); 
  m_recoJetDR05_Px ->clear(); 
  m_recoJetDR05_Py ->clear(); 
  m_recoJetDR05_Pz ->clear(); 

  m_recoJetDR06_E ->clear(); 
  m_recoJetDR06_Px ->clear(); 
  m_recoJetDR06_Py ->clear(); 
  m_recoJetDR06_Pz ->clear(); 

  m_recoJetDR07_E ->clear(); 
  m_recoJetDR07_Px ->clear(); 
  m_recoJetDR07_Py ->clear(); 
  m_recoJetDR07_Pz ->clear(); 

  m_recoJetDR08_E ->clear(); 
  m_recoJetDR08_Px ->clear(); 
  m_recoJetDR08_Py ->clear(); 
  m_recoJetDR08_Pz ->clear(); 

  m_recoJetDR09_E ->clear(); 
  m_recoJetDR09_Px ->clear(); 
  m_recoJetDR09_Py ->clear(); 
  m_recoJetDR09_Pz ->clear(); 

  m_recoJetDR10_E ->clear(); 
  m_recoJetDR10_Px ->clear(); 
  m_recoJetDR10_Py ->clear(); 
  m_recoJetDR10_Pz ->clear(); 

  if(m_fillMEInfo){
    m_trueME_E->clear();
    m_trueME_Px->clear();
    m_trueME_Py->clear();
    m_trueME_Pz->clear();
    m_trueME_PDGID->clear();
  }


  LCCollection * mcColl =0;
  getCollection(mcColl,m_inputMCParticleCollection,evt);

  bool pass_W_boson_mass=true;
  bool pass_Z_boson_mass=true;
  if(mcColl!=NULL){

    //bool found_diboson_event=false;
    int boson_counter=0;
    int ind_second_boson=-1;
    for(int m =0; m< mcColl->getNumberOfElements(); m++){
      MCParticle* mcp= dynamic_cast<MCParticle*>( mcColl->getElementAt(m) ) ;

      if(abs(mcp->getPDG())==24 || mcp->getPDG()==23){
	boson_counter+=1;
	if(abs(mcp->getPDG())==24 && boson_counter<3){
	  if(fabs(mcp->getMass()-80.4)>20.0){
	    pass_W_boson_mass=false;
	    //std::cout<<evt->getEventNumber()<<" off shell W "<<mcp->getMass()<<"/"<< mcp->getPDG()<<"/"<<boson_counter<<std::endl;
	  }
	}
	if(mcp->getPDG()==23 && boson_counter<3){
	  if(fabs(mcp->getMass()-91.2)>20.0){
	    pass_Z_boson_mass=false;
	    //std::cout<<evt->getEventNumber()<<" off shell Z "<<mcp->getMass()<<"/"<< mcp->getPDG()<<"/"<<boson_counter<<std::endl;
	  }
	}
      }
      //first boson is the index immediately in front of the first boson
      if(boson_counter==2 && pass_W_boson_mass && ind_second_boson==-1 && m_fillMEInfo){//for ZZ passed anyway, for WW or WZ passsed if real W
	ind_second_boson=m;
	MCParticle* mcp_1st= dynamic_cast<MCParticle*>( mcColl->getElementAt(ind_second_boson-1) ) ;
	m_trueME_E->push_back(mcp_1st->getEnergy());
	m_trueME_Px->push_back(mcp_1st->getMomentum()[0]);
	m_trueME_Py->push_back(mcp_1st->getMomentum()[1]);
	m_trueME_Pz->push_back(mcp_1st->getMomentum()[2]);
	m_trueME_PDGID->push_back(mcp->getPDG());
      }
      if(ind_second_boson>-1 && m<(ind_second_boson+5)){
	m_trueME_E->push_back(mcp->getEnergy());
	m_trueME_Px->push_back(mcp->getMomentum()[0]);
	m_trueME_Py->push_back(mcp->getMomentum()[1]);
	m_trueME_Pz->push_back(mcp->getMomentum()[2]);
	m_trueME_PDGID->push_back(mcp->getPDG());
      }
      //fill first quark anti quark pair which appears in the history
      //this IS the hadronic W or Z
      if(abs(mcp->getPDG())<7 && m_d1_mcE<0) {
	m_d1_mcPDGID=mcp->getPDG();
	m_d1_mcE=mcp->getEnergy();
	m_d1_mcPx=mcp->getMomentum()[0];
	m_d1_mcPy=mcp->getMomentum()[1];
	m_d1_mcPz=mcp->getMomentum()[2];
      }
      if(m_d2_mcE<0 && (abs(mcp->getPDG())<7 && mcp->getPDG()==(-m_d1_mcPDGID))){
	m_d2_mcPDGID=mcp->getPDG();
	m_d2_mcE=mcp->getEnergy();
	m_d2_mcPx=mcp->getMomentum()[0];
	m_d2_mcPy=mcp->getMomentum()[1];
	m_d2_mcPz=mcp->getMomentum()[2];
      }
      if(mcp->getGeneratorStatus()==1){//visible sum of stable particles --> take neutrinos out
	if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
	  m_E_trueAll+=mcp->getEnergy();
	  m_px_trueAll+=mcp->getMomentum()[0];
	  m_py_trueAll+=mcp->getMomentum()[1];
	  m_pz_trueAll+=mcp->getMomentum()[2];
	  if(fabs(mcp->getMomentum()[2]/sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]+mcp->getMomentum()[2]*mcp->getMomentum()[2]))<0.9){
	    m_E_trueAll_cos09+=mcp->getEnergy();
	    m_px_trueAll_cos09+=mcp->getMomentum()[0];
	    m_py_trueAll_cos09+=mcp->getMomentum()[1];
	    m_pz_trueAll_cos09+=mcp->getMomentum()[2];
	  }
	}else{
	  m_E_trueInv+=mcp->getEnergy();
	  m_px_trueInv+=mcp->getMomentum()[0];
	  m_py_trueInv+=mcp->getMomentum()[1];
	  m_pz_trueInv+=mcp->getMomentum()[2];
	  if(fabs(mcp->getMomentum()[2]/sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]+mcp->getMomentum()[2]*mcp->getMomentum()[2]))<0.9){
	    m_E_trueInv_cos09+=mcp->getEnergy();
	    m_px_trueInv_cos09+=mcp->getMomentum()[0];
	    m_py_trueInv_cos09+=mcp->getMomentum()[1];
	    m_pz_trueInv_cos09+=mcp->getMomentum()[2];
	  }
	}
      }
    }
  }

  LCCollection* recoparticlecol = NULL;
  recoparticlecol = evt->getCollection(m_inputRECOParticleCollection) ;
  if(recoparticlecol!=NULL){
    //PandoraCandidate loop
    for(int i=0;i<recoparticlecol->getNumberOfElements();i++){
      ReconstructedParticle* pandorapart = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(i));
      m_E_totPFO+=pandorapart->getEnergy();
      m_px_totPFO+=pandorapart->getMomentum()[0];
      m_py_totPFO+=pandorapart->getMomentum()[1];
      m_pz_totPFO+=pandorapart->getMomentum()[2];
      if(fabs(pandorapart->getMomentum()[2]/sqrt(pandorapart->getMomentum()[0]*pandorapart->getMomentum()[0]+pandorapart->getMomentum()[1]*pandorapart->getMomentum()[1]+pandorapart->getMomentum()[2]*pandorapart->getMomentum()[2]))<0.9){
	m_E_totPFO_cos09+=pandorapart->getEnergy();
	m_px_totPFO_cos09+=pandorapart->getMomentum()[0];
	m_py_totPFO_cos09+=pandorapart->getMomentum()[1];
	m_pz_totPFO_cos09+=pandorapart->getMomentum()[2];
      }
    }
  }

 LCCollection* recojetsDR07 = NULL;
  recojetsDR07 = evt->getCollection(m_recojetDR07ColName) ;
  if(recojetsDR07!=NULL){
    //PandoraCandidate loop
    for(int i=0;i<recojetsDR07->getNumberOfElements();i++){
      ReconstructedParticle* recojet = dynamic_cast<ReconstructedParticle*>(recojetsDR07->getElementAt(i));
      m_recoJetDR07_E->push_back(recojet->getEnergy());
      m_recoJetDR07_Px->push_back(recojet->getMomentum()[0]);
      m_recoJetDR07_Py->push_back(recojet->getMomentum()[1]);
      m_recoJetDR07_Pz->push_back(recojet->getMomentum()[2]);
    }
  }
  LCCollection* genjetsDR07 = NULL;
  genjetsDR07 = evt->getCollection(m_genjetDR07ColName) ;
  if(genjetsDR07!=NULL){
    //PandoraCandidate loop
    for(int i=0;i<genjetsDR07->getNumberOfElements();i++){
      ReconstructedParticle* genjet = dynamic_cast<ReconstructedParticle*>(genjetsDR07->getElementAt(i));
      m_genJetDR07_E->push_back(genjet->getEnergy());
      m_genJetDR07_Px->push_back(genjet->getMomentum()[0]);
      m_genJetDR07_Py->push_back(genjet->getMomentum()[1]);
      m_genJetDR07_Pz->push_back(genjet->getMomentum()[2]);
    }
  }


  if(m_fillAllJets){
    LCCollection* recojetsDR03 = NULL;
    recojetsDR03 = evt->getCollection(m_recojetDR03ColName) ;
    if(recojetsDR03!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<recojetsDR03->getNumberOfElements();i++){
	ReconstructedParticle* recojet = dynamic_cast<ReconstructedParticle*>(recojetsDR03->getElementAt(i));
	m_recoJetDR03_E->push_back(recojet->getEnergy());
	m_recoJetDR03_Px->push_back(recojet->getMomentum()[0]);
	m_recoJetDR03_Py->push_back(recojet->getMomentum()[1]);
	m_recoJetDR03_Pz->push_back(recojet->getMomentum()[2]);
      }
    }
    LCCollection* recojetsDR04 = NULL;
    recojetsDR04 = evt->getCollection(m_recojetDR04ColName) ;
    if(recojetsDR04!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<recojetsDR04->getNumberOfElements();i++){
	ReconstructedParticle* recojet = dynamic_cast<ReconstructedParticle*>(recojetsDR04->getElementAt(i));
	m_recoJetDR04_E->push_back(recojet->getEnergy());
	m_recoJetDR04_Px->push_back(recojet->getMomentum()[0]);
	m_recoJetDR04_Py->push_back(recojet->getMomentum()[1]);
	m_recoJetDR04_Pz->push_back(recojet->getMomentum()[2]);
      }
    }
    LCCollection* recojetsDR05 = NULL;
    recojetsDR05 = evt->getCollection(m_recojetDR05ColName) ;
    if(recojetsDR05!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<recojetsDR05->getNumberOfElements();i++){
	ReconstructedParticle* recojet = dynamic_cast<ReconstructedParticle*>(recojetsDR05->getElementAt(i));
	m_recoJetDR05_E->push_back(recojet->getEnergy());
	m_recoJetDR05_Px->push_back(recojet->getMomentum()[0]);
	m_recoJetDR05_Py->push_back(recojet->getMomentum()[1]);
	m_recoJetDR05_Pz->push_back(recojet->getMomentum()[2]);
      }
    }
    LCCollection* recojetsDR06 = NULL;
    recojetsDR06 = evt->getCollection(m_recojetDR06ColName) ;
    if(recojetsDR06!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<recojetsDR06->getNumberOfElements();i++){
	ReconstructedParticle* recojet = dynamic_cast<ReconstructedParticle*>(recojetsDR06->getElementAt(i));
	m_recoJetDR06_E->push_back(recojet->getEnergy());
	m_recoJetDR06_Px->push_back(recojet->getMomentum()[0]);
	m_recoJetDR06_Py->push_back(recojet->getMomentum()[1]);
	m_recoJetDR06_Pz->push_back(recojet->getMomentum()[2]);
      }
    }
 
    LCCollection* recojetsDR08 = NULL;
    recojetsDR08 = evt->getCollection(m_recojetDR08ColName) ;
    if(recojetsDR08!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<recojetsDR08->getNumberOfElements();i++){
	ReconstructedParticle* recojet = dynamic_cast<ReconstructedParticle*>(recojetsDR08->getElementAt(i));
	m_recoJetDR08_E->push_back(recojet->getEnergy());
	m_recoJetDR08_Px->push_back(recojet->getMomentum()[0]);
	m_recoJetDR08_Py->push_back(recojet->getMomentum()[1]);
	m_recoJetDR08_Pz->push_back(recojet->getMomentum()[2]);
      }
    }
    LCCollection* recojetsDR09 = NULL;
    recojetsDR09 = evt->getCollection(m_recojetDR09ColName) ;
    if(recojetsDR09!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<recojetsDR09->getNumberOfElements();i++){
	ReconstructedParticle* recojet = dynamic_cast<ReconstructedParticle*>(recojetsDR09->getElementAt(i));
	m_recoJetDR09_E->push_back(recojet->getEnergy());
	m_recoJetDR09_Px->push_back(recojet->getMomentum()[0]);
	m_recoJetDR09_Py->push_back(recojet->getMomentum()[1]);
	m_recoJetDR09_Pz->push_back(recojet->getMomentum()[2]);
      }
    } 
    LCCollection* recojetsDR10 = NULL;
    recojetsDR10 = evt->getCollection(m_recojetDR10ColName) ;
    if(recojetsDR10!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<recojetsDR10->getNumberOfElements();i++){
	ReconstructedParticle* recojet = dynamic_cast<ReconstructedParticle*>(recojetsDR10->getElementAt(i));
	m_recoJetDR10_E->push_back(recojet->getEnergy());
	m_recoJetDR10_Px->push_back(recojet->getMomentum()[0]);
	m_recoJetDR10_Py->push_back(recojet->getMomentum()[1]);
	m_recoJetDR10_Pz->push_back(recojet->getMomentum()[2]);
      }
    }
    LCCollection* genjetsDR03 = NULL;
    genjetsDR03 = evt->getCollection(m_genjetDR03ColName) ;
    if(genjetsDR03!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<genjetsDR03->getNumberOfElements();i++){
	ReconstructedParticle* genjet = dynamic_cast<ReconstructedParticle*>(genjetsDR03->getElementAt(i));
	m_genJetDR03_E->push_back(genjet->getEnergy());
	m_genJetDR03_Px->push_back(genjet->getMomentum()[0]);
	m_genJetDR03_Py->push_back(genjet->getMomentum()[1]);
	m_genJetDR03_Pz->push_back(genjet->getMomentum()[2]);
      }
    }
    LCCollection* genjetsDR04 = NULL;
    genjetsDR04 = evt->getCollection(m_genjetDR04ColName) ;
    if(genjetsDR04!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<genjetsDR04->getNumberOfElements();i++){
	ReconstructedParticle* genjet = dynamic_cast<ReconstructedParticle*>(genjetsDR04->getElementAt(i));
	m_genJetDR04_E->push_back(genjet->getEnergy());
	m_genJetDR04_Px->push_back(genjet->getMomentum()[0]);
	m_genJetDR04_Py->push_back(genjet->getMomentum()[1]);
	m_genJetDR04_Pz->push_back(genjet->getMomentum()[2]);
      }
    }
    LCCollection* genjetsDR05 = NULL;
    genjetsDR05 = evt->getCollection(m_genjetDR05ColName) ;
    if(genjetsDR05!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<genjetsDR05->getNumberOfElements();i++){
	ReconstructedParticle* genjet = dynamic_cast<ReconstructedParticle*>(genjetsDR05->getElementAt(i));
	m_genJetDR05_E->push_back(genjet->getEnergy());
	m_genJetDR05_Px->push_back(genjet->getMomentum()[0]);
	m_genJetDR05_Py->push_back(genjet->getMomentum()[1]);
	m_genJetDR05_Pz->push_back(genjet->getMomentum()[2]);
      }
    }
    LCCollection* genjetsDR06 = NULL;
    genjetsDR06 = evt->getCollection(m_genjetDR06ColName) ;
    if(genjetsDR06!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<genjetsDR06->getNumberOfElements();i++){
      ReconstructedParticle* genjet = dynamic_cast<ReconstructedParticle*>(genjetsDR06->getElementAt(i));
      m_genJetDR06_E->push_back(genjet->getEnergy());
      m_genJetDR06_Px->push_back(genjet->getMomentum()[0]);
      m_genJetDR06_Py->push_back(genjet->getMomentum()[1]);
      m_genJetDR06_Pz->push_back(genjet->getMomentum()[2]);
      }
    }
    
    LCCollection* genjetsDR08 = NULL;
    genjetsDR08 = evt->getCollection(m_genjetDR08ColName) ;
    if(genjetsDR08!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<genjetsDR08->getNumberOfElements();i++){
	ReconstructedParticle* genjet = dynamic_cast<ReconstructedParticle*>(genjetsDR08->getElementAt(i));
	m_genJetDR08_E->push_back(genjet->getEnergy());
	m_genJetDR08_Px->push_back(genjet->getMomentum()[0]);
	m_genJetDR08_Py->push_back(genjet->getMomentum()[1]);
	m_genJetDR08_Pz->push_back(genjet->getMomentum()[2]);
      }
    }
    LCCollection* genjetsDR09 = NULL;
    genjetsDR09 = evt->getCollection(m_genjetDR09ColName) ;
    if(genjetsDR09!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<genjetsDR09->getNumberOfElements();i++){
	ReconstructedParticle* genjet = dynamic_cast<ReconstructedParticle*>(genjetsDR09->getElementAt(i));
	m_genJetDR09_E->push_back(genjet->getEnergy());
	m_genJetDR09_Px->push_back(genjet->getMomentum()[0]);
	m_genJetDR09_Py->push_back(genjet->getMomentum()[1]);
	m_genJetDR09_Pz->push_back(genjet->getMomentum()[2]);
      }
    } 
    LCCollection* genjetsDR10 = NULL;
    genjetsDR10 = evt->getCollection(m_genjetDR10ColName) ;
    if(genjetsDR10!=NULL){
      //PandoraCandidate loop
      for(int i=0;i<genjetsDR10->getNumberOfElements();i++){
	ReconstructedParticle* genjet = dynamic_cast<ReconstructedParticle*>(genjetsDR10->getElementAt(i));
	m_genJetDR10_E->push_back(genjet->getEnergy());
	m_genJetDR10_Px->push_back(genjet->getMomentum()[0]);
	m_genJetDR10_Py->push_back(genjet->getMomentum()[1]);
	m_genJetDR10_Pz->push_back(genjet->getMomentum()[2]);
      }
    }  
  }
  if(pass_W_boson_mass/* && pass_Z_boson_mass*/ ){
    m_outputTree->Fill();
  }//else{
    //std::cout<<evt->getEventNumber()<<" not filled due to offshell boson"<<std::endl;
    //}
}

void JetAnalyzer::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    
    streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

void JetAnalyzer::check(LCEvent*){
}

void JetAnalyzer::end(){
    
  m_rootFile->cd();  
 
  m_rootFile->Write();
  m_rootFile->Close();

}

