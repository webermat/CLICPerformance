#include "TTBarStudy.h"
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHitPlane.h>

#include "UTIL/LCRelationNavigator.h"
#include "CalorimeterHitType.h"

#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace lcio ;
using namespace marlin ;

TTBarStudy aTTBarStudy;


TTBarStudy::TTBarStudy() : Processor("TTBarStudy") {
    
    // modify processor description
    _description = "TTBarStudy calculates properties of calorimeter showers" ;

    registerInputCollection( LCIO::MCPARTICLE,
                            "MCParticleCollectionName",
                            "Name of the MCParticle input collection",
                            m_inputMCParticleCollection,
                            std::string("MCParticle"));
    
    registerProcessorParameter( "OutputRootFileName",
                                "ROOT File name to collect plots",
                                m_rootFileName,
                                std::string("showerStudy.root"));
    
    registerProcessorParameter("runTauMode",
			       "tau specific running mode",
			       m_runTauMode,
			       bool("false"));
    
    registerProcessorParameter("runGenOnly",
			       "run generator only part",
			       m_runGenOnly,
			       bool("false"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "TauJetCollection" , 
			     "Name of the tauJet collection"  ,
			     m_inputTauCollection,
			     std::string("TaJets")
			     );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RECOParticleCollectionName" , 
			     "Name of the reconstructed particle collection"  ,
			     m_inputRECOParticleCollection,
			     std::string("PandoraPFOs")
			     );
}


void TTBarStudy::init() {


  dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();
  
  const double position[3]={0,0,0}; // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3]={0,0,0}; // initialise object to hold magnetic field
  mainDetector.field().magneticField(position,magneticFieldVector); // get the magnetic field vector from dd4hep
    
  m_innerBField = magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)
  //m_innerBField = 4.0;//magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)
  //constant takes into account that omega is in untis of mm
  m_const_a=2.99792e-4;

  eventcount=0;
  
    // Print the initial parameters
    printParameters() ;

    // Reset counters
    m_runNumber = 0 ;
    m_eventNumber = 0 ;
    m_ttbar_decay_mode = -1 ; 

    m_rootFile = new TFile(m_rootFileName.c_str(),"RECREATE");
    
    //implement cluster track distance using hits -> check pandoraPFA in order how to do it properly

    m_outputTree = new TTree("showerData","showerData");

    m_trueEnergy = new std::vector<float>();
    m_true_Eptx = new std::vector<float>();
    m_true_Epty = new std::vector<float>();
    m_true_Eptz = new std::vector<float>();
    m_true_Eptr = new std::vector<float>();
    m_true_Vtxx = new std::vector<float>();
    m_true_Vtxy = new std::vector<float>();
    m_true_Vtxz = new std::vector<float>();
    m_true_Vtxr = new std::vector<float>();
    m_true_Px = new std::vector<float>();
    m_true_Py = new std::vector<float>();
    m_true_Pz = new std::vector<float>();
    m_true_CosTheta = new std::vector<float>();
    m_true_Phi = new std::vector<float>();
    m_true_PDGID = new std::vector<int>();
    m_true_GenStatus = new std::vector<int>();
    m_true_numDaughters = new std::vector<int>();
    m_true_decayTrackerCalo = new std::vector<int>();//1 for tracker, 2 for calo. 3 for both, 4 for leaving dector, 0 for else
    m_true_numMothers = new std::vector<int>();
    m_true_m1_PDGID = new std::vector<int>();
    //m_true_m2_PDGID = new std::vector<int>();
    //m_true_m2_status = new std::vector<int>();
    m_true_m1_index = new std::vector<int>();
    //m_true_m2_index = new std::vector<int>();
    m_true_m1_E = new std::vector<float>();
    //m_true_m2_E = new std::vector<float>();
    m_true_index = new std::vector<int>();

    m_true_all_6_deg_E = new std::vector<float>();
    m_true_all_6_deg_angle = new std::vector<float>();
    m_true_all_6_deg_index = new std::vector<int>();
    m_true_all_6_deg_PDG = new std::vector<int>();
    
    m_true_tauDaughter_Energy = new std::vector<float>();
    m_true_tauDaughter_Px = new std::vector<float>();
    m_true_tauDaughter_Py = new std::vector<float>();
    m_true_tauDaughter_Pz = new std::vector<float>();
    m_true_tauDaughter_Vtxx = new std::vector<float>();
    m_true_tauDaughter_Vtxy = new std::vector<float>();
    m_true_tauDaughter_Vtxz = new std::vector<float>();
    m_true_tauDaughter_Vtxr = new std::vector<float>();
    m_true_tauDaughter_PDGID = new std::vector<int>();
    m_true_tauDaughter_Charge = new std::vector<int>();
    m_true_tauDaughter_tauIndex = new std::vector<int>();
    m_true_tauDaughter_status = new std::vector<int>();
    m_true_tauDaughter_motherPDGID  = new std::vector<int>();
    m_true_tauDaughter_motherEnergy = new std::vector<float>();

    m_true_tauDaughter_all_6_deg_E = new std::vector<float>();
    m_true_tauDaughter_all_6_deg_angle = new std::vector<float>();
    m_true_tauDaughter_all_6_deg_index = new std::vector<int>();
    m_true_tauDaughter_all_6_deg_PDG = new std::vector<int>();


    if(m_runTauMode && !m_runGenOnly){  
      m_tauJet_Px        = new std::vector<float>(); 
      m_tauJet_Py        = new std::vector<float>(); 
      m_tauJet_Pz        = new std::vector<float>(); 
      m_tauJet_E         = new std::vector<float>(); 
      m_tauJet_Phi       = new std::vector<float>(); 
      m_tauJet_CosTheta  = new std::vector<float>(); 
      m_tauJet_neutMult  = new std::vector<int>();
      m_tauJet_chMult    = new std::vector<int>();
      m_tauJet_charge    = new std::vector<int>();
      
      m_tauJet_Part_Px       = new std::vector<float>(); 
      m_tauJet_Part_Py       = new std::vector<float>(); 
      m_tauJet_Part_Pz       = new std::vector<float>(); 
      m_tauJet_Part_E        = new std::vector<float>(); 
      m_tauJet_Part_charge   = new std::vector<int>();
      m_tauJet_Part_PDGID    = new std::vector<int>();
      m_tauJet_Part_JetIndex = new std::vector<int>();
    }

    if(!m_runGenOnly){
      m_recoEnergy = new std::vector<float>();
      m_reco_Px = new std::vector<float>();
      m_reco_Py = new std::vector<float>();
      m_reco_Pz = new std::vector<float>();
      m_reco_CosTheta = new std::vector<float>();
      m_reco_Phi = new std::vector<float>();
      m_reco_Charge = new std::vector<int>();
      
      m_reco_nTracks = new std::vector<int>();
      m_reco_track0_pt = new std::vector<float>();
      m_reco_track0_p = new std::vector<float>();
      m_reco_track0_nHits = new std::vector<int>();
      m_reco_track0_chi2OverNdof = new std::vector<float>();
      m_reco_nClusters = new std::vector<int>();
      m_reco_clusters_energy = new std::vector<float>();
      m_reco_PDGID = new std::vector<int>();
      m_reco_E_EB = new std::vector<float>();
      m_reco_E_EE = new std::vector<float>();
      m_reco_E_EO = new std::vector<float>();
      m_reco_E_HB = new std::vector<float>();
      m_reco_E_HE = new std::vector<float>();
      m_reco_E_HO = new std::vector<float>();
      
      m_reco_all_6_deg_index = new std::vector<int>();
      m_reco_all_6_deg_E = new std::vector<float>();
      m_reco_all_6_deg_angle = new std::vector<float>();
      //m_reco_all_6_deg_type = new std::vector<int>();
      m_reco_all_6_deg_PDG = new std::vector<int>();
    }
    m_trueEnergy->clear();
    m_true_Px->clear();
    m_true_Py->clear();
    m_true_Pz->clear();
    m_true_Eptx->clear();
    m_true_Epty->clear();
    m_true_Eptz->clear();
    m_true_Eptr->clear();
    m_true_Vtxx->clear();
    m_true_Vtxy->clear();
    m_true_Vtxz->clear();
    m_true_Vtxr->clear();
    m_true_index->clear();
    m_true_CosTheta->clear();
    m_true_Phi->clear();
    m_true_PDGID->clear();
    m_true_GenStatus->clear();
    m_true_numDaughters->clear();
    m_true_decayTrackerCalo->clear();
    m_true_numMothers->clear();
    m_true_m1_PDGID->clear();
    m_true_m1_E->clear();
    m_true_m1_index->clear();

    m_true_all_6_deg_E->clear();
    m_true_all_6_deg_angle->clear();
    m_true_all_6_deg_index->clear();
    m_true_all_6_deg_PDG->clear();
 
    m_true_tauDaughter_Energy->clear();
    m_true_tauDaughter_Px->clear();
    m_true_tauDaughter_Py->clear();
    m_true_tauDaughter_Pz->clear();
    m_true_tauDaughter_Vtxx->clear();
    m_true_tauDaughter_Vtxy->clear();
    m_true_tauDaughter_Vtxz->clear();
    m_true_tauDaughter_Vtxr->clear();
    m_true_tauDaughter_PDGID->clear();
    m_true_tauDaughter_Charge->clear();
    m_true_tauDaughter_tauIndex->clear();
    m_true_tauDaughter_status->clear();
    m_true_tauDaughter_motherPDGID  ->clear();
    m_true_tauDaughter_motherEnergy ->clear();

    m_true_tauDaughter_all_6_deg_E->clear();
    m_true_tauDaughter_all_6_deg_angle->clear();
    m_true_tauDaughter_all_6_deg_index->clear();
    m_true_tauDaughter_all_6_deg_PDG->clear();

    if(m_runTauMode && !m_runGenOnly){   
      m_tauJet_Px->clear(); 
      m_tauJet_Py->clear(); 
      m_tauJet_Pz->clear(); 
      m_tauJet_E->clear(); 
      m_tauJet_Phi->clear(); 
      m_tauJet_CosTheta->clear(); 
      m_tauJet_neutMult->clear();
      m_tauJet_chMult->clear();
      m_tauJet_charge->clear();
      m_tauJet_Part_Px    ->clear(); 
      m_tauJet_Part_Py    ->clear(); 
      m_tauJet_Part_Pz    ->clear(); 
      m_tauJet_Part_E     ->clear(); 
      m_tauJet_Part_charge->clear();
      m_tauJet_Part_PDGID ->clear();
      m_tauJet_Part_PDGID ->clear();
      m_tauJet_Part_JetIndex->clear();
    }
    if(!m_runGenOnly){
      m_recoEnergy->clear();
      m_reco_Px->clear();
      m_reco_Py->clear();
      m_reco_Pz->clear();
      m_reco_CosTheta->clear();
      m_reco_Phi->clear();
      m_reco_Charge->clear();
      m_reco_PDGID->clear();
      
      m_reco_nTracks->clear();
      m_reco_track0_pt->clear();
      m_reco_track0_p->clear();
      m_reco_track0_nHits->clear();
      m_reco_track0_chi2OverNdof->clear();
      m_reco_nClusters->clear();
      m_reco_clusters_energy->clear();
      m_reco_E_EB->clear();
      m_reco_E_EE->clear();
      m_reco_E_EO->clear();
      m_reco_E_HB->clear();
      m_reco_E_HE->clear();
      m_reco_E_HO->clear();
      
      m_reco_all_6_deg_index->clear();
      m_reco_all_6_deg_E->clear();
      m_reco_all_6_deg_angle->clear();
      //m_reco_all_6_deg_type->clear(); //0 all MCParticles, 1 all loose/2 all selected, 3 all tight
      m_reco_all_6_deg_PDG->clear();
    }

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

    m_E_top1  = 0;
    m_px_top1 = 0;
    m_py_top1 = 0;
    m_pz_top1 = 0;
    m_PDGID_top1 = 0;

    m_E_top2  = 0;
    m_px_top2 = 0;
    m_py_top2 = 0;
    m_pz_top2 = 0;
    m_PDGID_top2 = 0;


    m_outputTree->Branch("runNumber",&m_runNumber,"runNumber/I");
    m_outputTree->Branch("eventNumber",&m_eventNumber,"eventNumber/I");
    m_outputTree->Branch("ttbarDecayMode",&m_ttbar_decay_mode,"ttbarDecayMode/I");

    //true particle level, exclude neutrinos
    m_outputTree->Branch("E_trueAll" ,&m_E_trueAll, "E_trueAll/F");
    m_outputTree->Branch("Px_trueAll",&m_px_trueAll,"Px_trueAll/F");
    m_outputTree->Branch("Py_trueAll",&m_py_trueAll,"Py_trueAll/F");
    m_outputTree->Branch("Pz_trueAll",&m_pz_trueAll,"Pz_trueAll/F");

    m_outputTree->Branch("E_top1" ,&m_E_top1, "E_top1/F");
    m_outputTree->Branch("Px_top1",&m_px_top1,"Px_top1/F");
    m_outputTree->Branch("Py_top1",&m_py_top1,"Py_top1/F");
    m_outputTree->Branch("Pz_top1",&m_pz_top1,"Pz_top1/F");
    m_outputTree->Branch("PDGID_top1",&m_PDGID_top1,"PDGID_top1/I");

    m_outputTree->Branch("E_top2" ,&m_E_top2, "E_top2/F");
    m_outputTree->Branch("Px_top2",&m_px_top2,"Px_top2/F");
    m_outputTree->Branch("Py_top2",&m_py_top2,"Py_top2/F");
    m_outputTree->Branch("Pz_top2",&m_pz_top2,"Pz_top2/F");
    m_outputTree->Branch("PDGID_top2",&m_PDGID_top2,"PDGID_top2/I");
    
    
    //true particle level, only neutrinos
    m_outputTree->Branch("E_trueInv" ,&m_E_trueInv, "E_trueInv/F");
    m_outputTree->Branch("Px_trueInv",&m_px_trueInv,"Px_trueInv/F");
    m_outputTree->Branch("Py_trueInv",&m_py_trueInv,"Py_trueInv/F");
    m_outputTree->Branch("Pz_trueInv",&m_pz_trueInv,"Pz_trueInv/F");
    
    //reconstructed level
    m_outputTree->Branch("E_totPFO" ,&m_E_totPFO, "E_totPFO/F");
    m_outputTree->Branch("Px_totPFO",&m_px_totPFO,"Px_totPFO/F");
    m_outputTree->Branch("Py_totPFO",&m_py_totPFO,"Py_totPFO/F");
    m_outputTree->Branch("Pz_totPFO",&m_pz_totPFO,"Pz_totPFO/F");

    m_outputTree->Branch("true_Energy","std::vector< float >",m_trueEnergy);
    m_outputTree->Branch("true_Px","std::vector< float >",m_true_Px);
    m_outputTree->Branch("true_Py","std::vector< float >",m_true_Py);
    m_outputTree->Branch("true_Pz","std::vector< float >",m_true_Pz);
    m_outputTree->Branch("true_Vtxx","std::vector< float >",m_true_Vtxx);
    m_outputTree->Branch("true_Vtxy","std::vector< float >",m_true_Vtxy);
    m_outputTree->Branch("true_Vtxz","std::vector< float >",m_true_Vtxz);
    m_outputTree->Branch("true_Vtxr","std::vector< float >",m_true_Vtxr);
    m_outputTree->Branch("true_Eptx","std::vector< float >",m_true_Eptx);
    m_outputTree->Branch("true_Epty","std::vector< float >",m_true_Epty);
    m_outputTree->Branch("true_Eptz","std::vector< float >",m_true_Eptz);
    m_outputTree->Branch("true_Eptr","std::vector< float >",m_true_Eptr);
    m_outputTree->Branch("true_index","std::vector< int >",m_true_index);
    m_outputTree->Branch("true_CosTheta","std::vector< float >",m_true_CosTheta);
    //m_outputTree->Branch("true_Theta","std::vector< float >",m_true_Theta);
    m_outputTree->Branch("true_Phi","std::vector< float >",m_true_Phi);
    m_outputTree->Branch("true_PDGID","std::vector< int >",m_true_PDGID); 
    m_outputTree->Branch("true_GenStatus","std::vector< int >",m_true_GenStatus); 
    m_outputTree->Branch("true_numDaughters", "std::vector< int >" ,m_true_numDaughters);
    m_outputTree->Branch("true_decayTrackerCalo", "std::vector< int >" , m_true_decayTrackerCalo);
    m_outputTree->Branch("true_m1_PDGID","std::vector< int >",m_true_m1_PDGID);
    m_outputTree->Branch("true_m1_Energy","std::vector< float >",m_true_m1_E);
    //m_outputTree->Branch("true_m1_index","std::vector< int >",m_true_m1_index);

    m_outputTree->Branch("true_all_6_deg_index","std::vector< int >",m_true_all_6_deg_index);
    m_outputTree->Branch("true_all_6_deg_E","std::vector< float >",m_true_all_6_deg_E);
    m_outputTree->Branch("true_all_6_deg_angle","std::vector< float >",m_true_all_6_deg_angle);
    m_outputTree->Branch("true_all_6_deg_PDG","std::vector< int >",m_true_all_6_deg_PDG);
 
    if(!m_runGenOnly){
      m_outputTree->Branch("reco_Energy","std::vector< float >",m_recoEnergy);
      m_outputTree->Branch("reco_Px","std::vector< float >",m_reco_Px);
      m_outputTree->Branch("reco_Py","std::vector< float >",m_reco_Py);
      m_outputTree->Branch("reco_Pz","std::vector< float >",m_reco_Pz);
      m_outputTree->Branch("reco_CosTheta","std::vector< float >",m_reco_CosTheta);
      m_outputTree->Branch("reco_Phi","std::vector< float >",m_reco_Phi);
      m_outputTree->Branch("reco_PDGID","std::vector< int >",m_reco_PDGID);
      m_outputTree->Branch("reco_Charge","std::vector< int >",m_reco_Charge);
      
      m_outputTree->Branch("reco_nTracks","std::vector< int >",m_reco_nTracks);
      m_outputTree->Branch("reco_track0_pt","std::vector< float >",m_reco_track0_pt);
      m_outputTree->Branch("reco_track0_p","std::vector< float >",m_reco_track0_p);
      m_outputTree->Branch("reco_track0_nHits","std::vector< int >",m_reco_track0_nHits);
      m_outputTree->Branch("reco_track0_chi2OverNdof","std::vector< float >",m_reco_track0_chi2OverNdof);
      m_outputTree->Branch("reco_nClusters","std::vector< int >",m_reco_nClusters);
      m_outputTree->Branch("reco_clusters_energy","std::vector< float >",m_reco_clusters_energy);
      m_outputTree->Branch("reco_E_EB","std::vector< float >",m_reco_E_EB);
      m_outputTree->Branch("reco_E_EE","std::vector< float >",m_reco_E_EE);
      m_outputTree->Branch("reco_E_EO","std::vector< float >",m_reco_E_EO);
      m_outputTree->Branch("reco_E_HB","std::vector< float >",m_reco_E_HB);
      m_outputTree->Branch("reco_E_HE","std::vector< float >",m_reco_E_HE);
      m_outputTree->Branch("reco_E_HO","std::vector< float >",m_reco_E_HO);
      
      m_outputTree->Branch("reco_all_6_deg_index","std::vector< int >",m_reco_all_6_deg_index);
      m_outputTree->Branch("reco_all_6_deg_E","std::vector< float >",m_reco_all_6_deg_E);
      m_outputTree->Branch("reco_all_6_deg_angle","std::vector< float >",m_reco_all_6_deg_angle);
      //m_outputTree->Branch("reco_all_6_deg_type","std::vector< int >",m_reco_all_6_deg_type);
      m_outputTree->Branch("reco_all_6_deg_PDG","std::vector< int >",m_reco_all_6_deg_PDG);
    }
    if(m_runTauMode){
      m_outputTree->Branch("trueTauDaughterE", "std::vector< float >",m_true_tauDaughter_Energy);
      m_outputTree->Branch("trueTauDaughterPx", "std::vector< float >",m_true_tauDaughter_Px);
      m_outputTree->Branch("trueTauDaughterPy", "std::vector< float >",m_true_tauDaughter_Py);
      m_outputTree->Branch("trueTauDaughterPz", "std::vector< float >",m_true_tauDaughter_Pz);
      m_outputTree->Branch("trueTauDaughterVtxx", "std::vector< float >",m_true_tauDaughter_Vtxx);
      m_outputTree->Branch("trueTauDaughterVtxy", "std::vector< float >",m_true_tauDaughter_Vtxy);
      m_outputTree->Branch("trueTauDaughterVtxz", "std::vector< float >",m_true_tauDaughter_Vtxz);
      m_outputTree->Branch("trueTauDaughterVtxr", "std::vector< float >",m_true_tauDaughter_Vtxr);
      m_outputTree->Branch("trueTauDaughterPDGID", "std::vector< int >",m_true_tauDaughter_PDGID);
      m_outputTree->Branch("trueTauDaughterCharge", "std::vector< int >",m_true_tauDaughter_Charge);
      m_outputTree->Branch("trueTauDaughterTauIndex", "std::vector< int >",m_true_tauDaughter_tauIndex);
      m_outputTree->Branch("trueTauDaughterStatus", "std::vector< int >",m_true_tauDaughter_status);
      m_outputTree->Branch("trueTauDaughterMotherPDGID", "std::vector< int >",m_true_tauDaughter_motherPDGID); //direct mother, e.g. pi0
      m_outputTree->Branch("trueTauDaughterMotherEnergy", "std::vector< float >",m_true_tauDaughter_motherEnergy);
      
      m_outputTree->Branch("trueTauDaughter_all_6_deg_index","std::vector< int >",m_true_tauDaughter_all_6_deg_index);
      m_outputTree->Branch("trueTauDaughter_all_6_deg_E","std::vector< float >",m_true_tauDaughter_all_6_deg_E);
      m_outputTree->Branch("trueTauDaughter_all_6_deg_angle","std::vector< float >",m_true_tauDaughter_all_6_deg_angle);
      m_outputTree->Branch("trueTauDaughter_all_6_deg_PDG","std::vector< int >",m_true_tauDaughter_all_6_deg_PDG);
    }

    if(m_runTauMode && !m_runGenOnly){
      m_outputTree->Branch("tauJetPx", "std::vector< float >", &m_tauJet_Px); 
      m_outputTree->Branch("tauJetPy", "std::vector< float >", &m_tauJet_Py); 
      m_outputTree->Branch("tauJetPz", "std::vector< float >", &m_tauJet_Pz); 
      m_outputTree->Branch("tauJetE", "std::vector< float >", &m_tauJet_E); 
      //m_outputTree->Branch("tauJetPhi", "std::vector< float >", &m_tauJet_Phi); 
      m_outputTree->Branch("tauJetCosTheta", "std::vector< float >", &m_tauJet_CosTheta); 
      m_outputTree->Branch("tauJetNeutMult", "std::vector< int >", &m_tauJet_neutMult);
      m_outputTree->Branch("tauJetChMult", "std::vector< int >", &m_tauJet_chMult);
      m_outputTree->Branch("tauJetcharge", "std::vector< int >", &m_tauJet_charge);
      
      m_outputTree->Branch("tauJet_Part_Px", "std::vector< float >", &m_tauJet_Part_Px); 
      m_outputTree->Branch("tauJet_Part_Py", "std::vector< float >", &m_tauJet_Part_Py); 
      m_outputTree->Branch("tauJet_Part_Pz", "std::vector< float >", &m_tauJet_Part_Pz); 
      m_outputTree->Branch("tauJet_Part_E", "std::vector< float >", &m_tauJet_Part_E); 
      m_outputTree->Branch("tauJet_Part_charge", "std::vector< int >", &m_tauJet_Part_charge);
      m_outputTree->Branch("tauJet_Part_PDGID", "std::vector< int >", &m_tauJet_Part_PDGID);
      m_outputTree->Branch("tauJet_Part_JetIndex", "std::vector< int >", &m_tauJet_Part_JetIndex);
    }
}


void TTBarStudy::processRunHeader( LCRunHeader*) {
  //++m_runNumber ;
}

void TTBarStudy::processEvent( LCEvent* evt ) {

  m_runNumber=evt->getRunNumber();
  m_eventNumber=evt->getEventNumber();

  eventcount+=1;
  if(evt->getEventNumber()%50==0){
    std::cout<<"run/evt "<<evt->getRunNumber()<<"/"<<evt->getEventNumber()<<"/"<<eventcount<<std::endl;
  }
    m_trueEnergy->clear();
    m_true_Px->clear();
    m_true_Py->clear();
    m_true_Pz->clear();
    m_true_Vtxx->clear();
    m_true_Vtxy->clear();
    m_true_Vtxz->clear();
    m_true_Vtxr->clear();
    m_true_Eptx->clear();
    m_true_Epty->clear();
    m_true_Eptz->clear();
    m_true_Eptr->clear();
    m_true_CosTheta->clear();
    m_true_Phi->clear();
    m_true_PDGID->clear();
    m_true_GenStatus->clear();
    m_true_numDaughters->clear();
    m_true_decayTrackerCalo->clear();
    m_true_numMothers->clear();
    m_true_m1_PDGID->clear();
    m_true_m1_E->clear();
    m_true_m1_index->clear();
    m_true_all_6_deg_index->clear();
    m_true_all_6_deg_E->clear();
    m_true_all_6_deg_angle->clear();
    m_true_all_6_deg_PDG->clear();

    m_true_tauDaughter_Energy->clear();
    m_true_tauDaughter_Px->clear();
    m_true_tauDaughter_Py->clear();
    m_true_tauDaughter_Pz->clear();
    m_true_tauDaughter_Vtxx->clear();
    m_true_tauDaughter_Vtxy->clear();
    m_true_tauDaughter_Vtxz->clear();
    m_true_tauDaughter_Vtxr->clear();
    m_true_tauDaughter_PDGID->clear();
    m_true_tauDaughter_Charge->clear();
    m_true_tauDaughter_tauIndex->clear();
    m_true_tauDaughter_status->clear();
    m_true_tauDaughter_motherPDGID  ->clear();
    m_true_tauDaughter_motherEnergy ->clear();

    m_true_tauDaughter_all_6_deg_E->clear();
    m_true_tauDaughter_all_6_deg_angle->clear();
    m_true_tauDaughter_all_6_deg_index->clear();
    m_true_tauDaughter_all_6_deg_PDG->clear();

    
    if(m_runTauMode && !m_runGenOnly){  
      m_tauJet_Px->clear(); 
      m_tauJet_Py->clear(); 
      m_tauJet_Pz->clear(); 
      m_tauJet_E->clear(); 
      m_tauJet_Phi->clear(); 
      m_tauJet_CosTheta->clear(); 
      m_tauJet_neutMult->clear();
      m_tauJet_chMult->clear();
      m_tauJet_charge->clear();
      m_tauJet_Part_Px    ->clear(); 
      m_tauJet_Part_Py    ->clear(); 
      m_tauJet_Part_Pz    ->clear(); 
      m_tauJet_Part_E     ->clear(); 
      m_tauJet_Part_charge->clear();
      m_tauJet_Part_PDGID ->clear();
      m_tauJet_Part_JetIndex->clear();
    }
    if(!m_runGenOnly){
      m_recoEnergy->clear();
      m_reco_Px->clear();
      m_reco_Py->clear();
      m_reco_Pz->clear();
      m_reco_CosTheta->clear();
      m_reco_Phi->clear();
      m_reco_Charge->clear();
      m_reco_PDGID->clear();
      m_reco_nTracks->clear();
      m_reco_track0_pt->clear();
      m_reco_track0_p->clear();
      m_reco_track0_nHits->clear();
      m_reco_track0_chi2OverNdof->clear();
      m_reco_nClusters->clear();
      m_reco_clusters_energy->clear();
      m_reco_E_EB->clear();
      m_reco_E_EE->clear();
      m_reco_E_EO->clear();
      m_reco_E_HB->clear();
      m_reco_E_HE->clear();
      m_reco_E_HO->clear();
      
      m_reco_all_6_deg_E->clear();
      m_reco_all_6_deg_angle->clear();
      //m_reco_all_6_deg_type->clear(); //0 all MCParticles, 1 all loose/2 all selected, 3 all tight
      m_reco_all_6_deg_index->clear(); 
      m_reco_all_6_deg_PDG->clear(); 
    }


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

    m_E_top1  = 0;
    m_px_top1 = 0;
    m_py_top1 = 0;
    m_pz_top1 = 0;
    m_PDGID_top1 = 0;

    m_E_top2  = 0;
    m_px_top2 = 0;
    m_py_top2 = 0;
    m_pz_top2 = 0;
    m_PDGID_top2 = 0;

    //only run tracks and clusters for non tau events ->>validation and efficiency studies
    if(m_runTauMode && !m_runGenOnly){
      //tau jet loop
      LCCollection * tauJetColl =0;
      getCollection(tauJetColl,m_inputTauCollection,evt);
      if(tauJetColl!=NULL){
	for(int i=0;i<tauJetColl->getNumberOfElements();i++){
	  ReconstructedParticle* tauJet = dynamic_cast<ReconstructedParticle*>(tauJetColl->getElementAt(i));
	  m_tauJet_charge->push_back(tauJet->getCharge()); 
	  m_tauJet_Px->push_back(tauJet->getMomentum()[0]); 
	  m_tauJet_Py->push_back(tauJet->getMomentum()[1]); 
	  m_tauJet_Pz->push_back(tauJet->getMomentum()[2]); 
	  m_tauJet_E->push_back(tauJet->getEnergy()); 
	  m_tauJet_Phi->push_back(atan2(tauJet->getMomentum()[1],tauJet->getMomentum()[0])); 
	  double cosTheta = tauJet->getMomentum()[2]/sqrt(tauJet->getMomentum()[0]*tauJet->getMomentum()[0]+tauJet->getMomentum()[1]*tauJet->getMomentum()[1]+tauJet->getMomentum()[2]*tauJet->getMomentum()[2]);
	  m_tauJet_CosTheta->push_back(cosTheta);
	  int tauJetNeutMult=0;
	  int tauJetChMult=0;
	  for(unsigned int j=0;j<tauJet->getParticles().size();j++){
	    m_tauJet_Part_Px->push_back(tauJet->getParticles()[j]->getMomentum()[0]); 
	    m_tauJet_Part_Py->push_back(tauJet->getParticles()[j]->getMomentum()[1]); 
	    m_tauJet_Part_Pz->push_back(tauJet->getParticles()[j]->getMomentum()[2]); 
	    m_tauJet_Part_E->push_back(tauJet->getParticles()[j]->getEnergy()); 
	    m_tauJet_Part_charge->push_back(tauJet->getParticles()[j]->getCharge());
	    m_tauJet_Part_PDGID->push_back(tauJet->getParticles()[j]->getType());
	    m_tauJet_Part_JetIndex->push_back(i);
	    if(tauJet->getParticles()[j]->getCharge()==0){
	      tauJetNeutMult+=1;
	    }else{
	      tauJetChMult+=1;
	    }
	  }
	  m_tauJet_neutMult->push_back(tauJetNeutMult);
	  m_tauJet_chMult->push_back(tauJetChMult);
	}
      }
      //tau jet particle loop       
    }
    //std::cout<<"before mc collection"<<std::endl;
      
      LCCollection * mcColl =0;
      getCollection(mcColl,m_inputMCParticleCollection,evt);
      std::vector<TLorentzVector>tau_veto_vector;
      unsigned int n_W_found=0;
      bool W1_is_mu_el=false;
      bool W1_is_tau_lep=false;
      bool W1_is_tau_had=false;
      bool W2_is_mu_el=false;
      bool W2_is_tau_lep=false;
      bool W2_is_tau_had=false;
      std::set<MCParticle*> muon_electron_Func;
      for(int m =0; m< mcColl->getNumberOfElements(); m++){
	MCParticle* mcp= dynamic_cast<MCParticle*>( mcColl->getElementAt(m) ) ;
	if(m==4){
	  m_E_top1  = mcp->getEnergy();
	  m_px_top1 = mcp->getMomentum()[0];
	  m_py_top1 = mcp->getMomentum()[1];
	  m_pz_top1 = mcp->getMomentum()[2];
	  m_PDGID_top1 = mcp->getPDG();
	}
	if(m==5){
	  m_E_top2  = mcp->getEnergy();
	  m_px_top2 = mcp->getMomentum()[0];
	  m_py_top2 = mcp->getMomentum()[1];
	  m_pz_top2 = mcp->getMomentum()[2];
	  m_PDGID_top2 = mcp->getPDG();
	}
	if(mcp->getGeneratorStatus()==1){//visible sum of stable particles --> take neutrinos out
	  if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 &&  abs(mcp->getPDG())!=16){
	    m_E_trueAll+=mcp->getEnergy();
	    m_px_trueAll+=mcp->getMomentum()[0];
	    m_py_trueAll+=mcp->getMomentum()[1];
	    m_pz_trueAll+=mcp->getMomentum()[2];
	  }else{
	    m_E_trueInv+=mcp->getEnergy();
	    m_px_trueInv+=mcp->getMomentum()[0];
	    m_py_trueInv+=mcp->getMomentum()[1];
	    m_pz_trueInv+=mcp->getMomentum()[2];
	  }
	}


	int index_true=m_trueEnergy->size(); //starts at 0 counting up
	//FOR BBAR EVENTS NOW ----------------------------------->>>>>>>>>>>>>>>>>>>>>>>
	/*
	bool fill_all_stable_daughters=false;
	if( ((abs(mcp->getPDG())>510 && abs(mcp->getPDG())<600) || (abs(mcp->getPDG())>5100 && abs(mcp->getPDG())<5600))|| ((abs(mcp->getPDG())>410 && abs(mcp->getPDG())<500) || (abs(mcp->getPDG())>4100 && abs(mcp->getPDG())<4600)) ){
	  //since D mesons should be decay product of B's, they should have been taken care off
	  //if((abs(mcp->getPDG())>510 && abs(mcp->getPDG())<600) || (abs(mcp->getPDG())>5100 && abs(mcp->getPDG())<5600)){
	    for(unsigned int wd=0;wd<mcp->getDaughters().size();wd++){
	      if(mcp->getDaughters()[wd]->getGeneratorStatus()==1){
		fill_all_stable_daughters=true;
	      }
	    }
	    //}
	  if(fill_all_stable_daughters){
	    MCParticle* tau=mcp;
	    double cosTheta = tau->getMomentum()[2]/sqrt(tau->getMomentum()[0]*tau->getMomentum()[0]+tau->getMomentum()[1]*tau->getMomentum()[1]+tau->getMomentum()[2]*tau->getMomentum()[2]);
	    m_trueEnergy->push_back(tau->getEnergy());
	    m_true_Px->push_back(tau->getMomentum()[0]);
	    m_true_Py->push_back(tau->getMomentum()[1]);
	    m_true_Pz->push_back(tau->getMomentum()[2]);
	    m_true_CosTheta->push_back(cosTheta);    
	    m_true_Phi->push_back(atan2(tau->getMomentum()[1],tau->getMomentum()[0]));
	    m_true_Vtxx->push_back(tau->getVertex()[0]);
	    m_true_Vtxy->push_back(tau->getVertex()[1]);
	    m_true_Vtxz->push_back(tau->getVertex()[2]);
	    m_true_Vtxr->push_back(sqrt(pow(tau->getVertex()[0],2)+pow(tau->getVertex()[1],2)));
	    //std::cout<<"epx/epy/epx/epr "<<tau->getEndpoint()[0]<<"/"<<tau->getEndpoint()[1]<<"/"<<tau->getEndpoint()[2]<<"/"<<sqrt(pow(tau->getEndpoint()[0],2)+pow(tau->getEndpoint()[1],2))<<std::endl;
	    m_true_Eptx->push_back(tau->getEndpoint()[0]);
	    m_true_Epty->push_back(tau->getEndpoint()[1]);
	    m_true_Eptz->push_back(tau->getEndpoint()[2]);
	    m_true_Eptr->push_back(sqrt(pow(tau->getEndpoint()[0],2)+pow(tau->getEndpoint()[1],2)));
	    m_true_PDGID->push_back(tau->getPDG());
	    m_true_index->push_back(index_true);
	    m_true_GenStatus->push_back(tau->getGeneratorStatus());
	    m_true_numDaughters->push_back(tau->getDaughters().size ());
	    m_true_numMothers->push_back(tau->getParents().size());
	    m_true_m1_PDGID->push_back(tau->getParents()[0]->getPDG());
	    m_true_m1_E->push_back(tau->getParents()[0]->getEnergy());
	    m_true_decayTrackerCalo->push_back(-1);
	    for(unsigned int wd=0;wd<mcp->getDaughters().size();wd++){
	      m_true_tauDaughter_Energy->push_back(mcp->getDaughters()[wd]->getEnergy());
	      m_true_tauDaughter_Px->push_back(mcp->getDaughters()[wd]->getMomentum()[0]);
	      m_true_tauDaughter_Py->push_back(mcp->getDaughters()[wd]->getMomentum()[1]);
	      m_true_tauDaughter_Pz->push_back(mcp->getDaughters()[wd]->getMomentum()[2]);
	      m_true_tauDaughter_Vtxx->push_back(mcp->getDaughters()[wd]->getVertex()[0]);
	      m_true_tauDaughter_Vtxy->push_back(mcp->getDaughters()[wd]->getVertex()[1]);
	      m_true_tauDaughter_Vtxz->push_back(mcp->getDaughters()[wd]->getVertex()[2]);
	      m_true_tauDaughter_Vtxr->push_back(sqrt(pow(mcp->getDaughters()[wd]->getVertex()[0],2)+pow(mcp->getDaughters()[wd]->getVertex()[1],2)));
	      m_true_tauDaughter_PDGID->push_back(mcp->getDaughters()[wd]->getPDG());
	      m_true_tauDaughter_Charge->push_back(mcp->getDaughters()[wd]->getCharge());
	      m_true_tauDaughter_tauIndex->push_back(m);
	      m_true_tauDaughter_status->push_back(mcp->getDaughters()[wd]->getGeneratorStatus());
	      m_true_tauDaughter_motherPDGID->push_back(mcp->getDaughters()[wd]->getParents()[0]->getPDG());
	      m_true_tauDaughter_motherEnergy->push_back(mcp->getDaughters()[wd]->getParents()[0]->getEnergy());
	    }
	  }
	}
	*/
	//BBAR DONE here -------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//--------------------------------------------- now the actual ttbar W lepton code runs
	
	if(abs(mcp->getPDG()) == 24 && mcp->getDaughters().size ()>1){ 
	  n_W_found+=1;
	  if(n_W_found<3 || n_W_found>4){//1 and 2 are whizard internal, the candidates are then given to pythia and the proper decays move on from there
	    //so use indices 3 and 4
	    continue;
	  }
	  //check for lepton decays, also check if it is a tau decay
	  if((abs(mcp->getDaughters()[0]->getPDG())>=11 && abs(mcp->getDaughters()[0]->getPDG())<=16)||(abs(mcp->getDaughters()[1]->getPDG())>=11 && abs(mcp->getDaughters()[1]->getPDG())<=16)){
	  //check taus first
	    //std::cout<<"gen part 1"<<std::endl;
	    if((abs(mcp->getDaughters()[0]->getPDG())==15 || abs(mcp->getDaughters()[0]->getPDG())==16)||(abs(mcp->getDaughters()[1]->getPDG())==15 || abs(mcp->getDaughters()[1]->getPDG())==16)){
	      for(unsigned int wd=0;wd<mcp->getDaughters().size();wd++){
		if(abs(mcp->getDaughters()[wd]->getPDG())==15){
		  if(!m_runTauMode){
		    continue;
		  }
		  MCParticle* tau=mcp->getDaughters()[wd];
		  double cosTheta = tau->getMomentum()[2]/sqrt(tau->getMomentum()[0]*tau->getMomentum()[0]+tau->getMomentum()[1]*tau->getMomentum()[1]+tau->getMomentum()[2]*tau->getMomentum()[2]);
		  m_trueEnergy->push_back(tau->getEnergy());
		  m_true_Px->push_back(tau->getMomentum()[0]);
		  m_true_Py->push_back(tau->getMomentum()[1]);
		  m_true_Pz->push_back(tau->getMomentum()[2]);
		  m_true_CosTheta->push_back(cosTheta);    
		  m_true_Phi->push_back(atan2(tau->getMomentum()[1],tau->getMomentum()[0]));
		  m_true_PDGID->push_back(tau->getPDG());
		  m_true_index->push_back(index_true);
		  m_true_GenStatus->push_back(tau->getGeneratorStatus());
		  m_true_numDaughters->push_back(tau->getDaughters().size ());
		  m_true_numMothers->push_back(tau->getParents().size());
		  m_true_m1_PDGID->push_back(tau->getParents()[0]->getPDG());
		  m_true_m1_E->push_back(tau->getParents()[0]->getEnergy());
		  m_true_decayTrackerCalo->push_back(-1);
		  std::set<MCParticle*> tau_daughtersFunc;
		  fillStableDaughterSet(tau, tau_daughtersFunc);
		  std::set<MCParticle*>::iterator tauDaughtIt;
		  bool tau_has_el_mu_daughter=false;
		  bool tau_has_hadron_daughter=false;
		  //std::cout<<"should be in tau decays "<<std::endl;
		  //set flags to decide if particle isolation for muons and leptons is needed
		  //leave in for now to check tau 1 prongs as well
		  for(tauDaughtIt=tau_daughtersFunc.begin();tauDaughtIt!=tau_daughtersFunc.end();tauDaughtIt++){
		    m_true_tauDaughter_Energy->push_back((*tauDaughtIt)->getEnergy());
		    m_true_tauDaughter_Px->push_back((*tauDaughtIt)->getMomentum()[0]);
		    m_true_tauDaughter_Py->push_back((*tauDaughtIt)->getMomentum()[1]);
		    m_true_tauDaughter_Pz->push_back((*tauDaughtIt)->getMomentum()[2]);
		    m_true_tauDaughter_PDGID->push_back((*tauDaughtIt)->getPDG());
		    if(abs((*tauDaughtIt)->getPDG())==11 || abs((*tauDaughtIt)->getPDG())==13){
		      //std::cout<<"tau it type "<<(*tauDaughtIt)->getPDG()<<" "<<n_W_found <<" "<<(*tauDaughtIt)->getEnergy() <<" stat "<<(*tauDaughtIt)->getGeneratorStatus()<<std::endl;
		      tau_has_el_mu_daughter=true;
		      muon_electron_Func.insert((*tauDaughtIt));
		      if(abs(((*tauDaughtIt)->getPDG())==11 || abs((*tauDaughtIt)->getPDG())==13)){
			//we check in the end for true electron and muon isolations
			TLorentzVector tauDLep(0,0,0,0);
			int ind_tauD=m_true_tauDaughter_Energy->size()-1;
			tauDLep.SetPxPyPzE((*tauDaughtIt)->getMomentum()[0],(*tauDaughtIt)->getMomentum()[1],(*tauDaughtIt)->getMomentum()[2],(*tauDaughtIt)->getEnergy());
			for(int m_R =0; m_R< mcColl->getNumberOfElements(); m_R++){
			  MCParticle* mcp_R= dynamic_cast<MCParticle*>( mcColl->getElementAt(m_R) ) ;
			  if(mcp_R!=(*tauDaughtIt)  &&  mcp_R->getGeneratorStatus()==1 &&((abs(mcp_R->getPDG())!=12)||(abs(mcp_R->getPDG())!=14)||(abs(mcp_R->getPDG())!=16)) ){
			    TLorentzVector temp(0,0,0,0);
			    temp.SetPxPyPzE(mcp_R->getMomentum()[0],mcp_R->getMomentum()[1],mcp_R->getMomentum()[2],mcp_R->getEnergy());
			    if((tauDLep.Angle(temp.Vect())/M_PI*180.)<6.0){
			      m_true_tauDaughter_all_6_deg_E->push_back(mcp_R->getEnergy());
			      m_true_tauDaughter_all_6_deg_angle->push_back(tauDLep.Angle(temp.Vect()));
			      m_true_tauDaughter_all_6_deg_index->push_back(ind_tauD);
			      m_true_tauDaughter_all_6_deg_PDG->push_back(mcp_R->getPDG());
			    }
			  }
			}
		      }
		    }else if(abs((*tauDaughtIt)->getPDG())>22){//we don't have BSM stuff in the sample
		      tau_has_hadron_daughter=true;
		    }
		    m_true_tauDaughter_Charge->push_back((*tauDaughtIt)->getCharge());
		    m_true_tauDaughter_tauIndex->push_back(index_true);
		    m_true_tauDaughter_status->push_back((*tauDaughtIt)->getGeneratorStatus());
		    m_true_tauDaughter_motherPDGID->push_back((*tauDaughtIt)->getParents()[0]->getPDG());
		    m_true_tauDaughter_motherEnergy->push_back((*tauDaughtIt)->getParents()[0]->getEnergy());
		  }
		  //in this case, there is always a lep+/lep- (typically e+e-) pair present
		  //assign such cases to hadronic taus
		  if(tau_has_hadron_daughter){
		    if(n_W_found==3){
		      W1_is_tau_had=true;
		    }else if (n_W_found==4){
		      W2_is_tau_had=true;
		    }
		  }else if (tau_has_el_mu_daughter){
		    if(n_W_found==3){
		      W1_is_tau_lep=true;
		    }else if (n_W_found==4){
		      W2_is_tau_lep=true;
		    }
		  }else{
		    std::cout<<"sth is wrong in tau decay classification, should have hadron or el/mu, but doesn't seem to be the case"<<std::endl;
		    for(tauDaughtIt=tau_daughtersFunc.begin();tauDaughtIt!=tau_daughtersFunc.end();tauDaughtIt++){
		      std::cout<<"tau it type "<<(*tauDaughtIt)->getPDG()<<std::endl;
		    }
		  }
		  if(!tau_daughtersFunc.empty()){
		    tau_daughtersFunc.clear();
		  }
		}//check if it is the tau for W decays to tau-tau neutrino
	      }
	    }else{//lepton decay, but no tau's so
	      //direct W electrons and muons decays
	      std::set<MCParticle*> w_daughtersFunc;
	      fillStableDaughterSet(mcp, w_daughtersFunc);
	      //std::cout<<"w should be lepton decays"<<std::endl;
	      //std::set<MCParticle*>::iterator tauDaughtIt;
	      for(std::set<MCParticle*>::iterator  wDaughtIt=w_daughtersFunc.begin();wDaughtIt!=w_daughtersFunc.end();wDaughtIt++){
		if(abs((*wDaughtIt)->getPDG())==11 || abs((*wDaughtIt)->getPDG())==13){
		  //std::cout<<"w it type "<<(*wDaughtIt)->getPDG()<<" "<<n_W_found <<" "<<(*wDaughtIt)->getEnergy() <<" stat "<<(*wDaughtIt)->getGeneratorStatus()<<std::endl;
		  if(n_W_found==3){
		    W1_is_mu_el=true;
		  }else if (n_W_found==4){
		    W2_is_mu_el=true;
		  }
		  muon_electron_Func.insert((*wDaughtIt));
		  double cosTheta = (*wDaughtIt)->getMomentum()[2]/sqrt((*wDaughtIt)->getMomentum()[0]*(*wDaughtIt)->getMomentum()[0]+(*wDaughtIt)->getMomentum()[1]*(*wDaughtIt)->getMomentum()[1]+(*wDaughtIt)->getMomentum()[2]*(*wDaughtIt)->getMomentum()[2]);
		  m_trueEnergy->push_back((*wDaughtIt)->getEnergy());
		  m_true_Px->push_back((*wDaughtIt)->getMomentum()[0]);
		  m_true_Py->push_back((*wDaughtIt)->getMomentum()[1]);
		  m_true_Pz->push_back((*wDaughtIt)->getMomentum()[2]);
		  m_true_CosTheta->push_back(cosTheta);    
		  m_true_Phi->push_back(atan2((*wDaughtIt)->getMomentum()[1],(*wDaughtIt)->getMomentum()[0]));
		  m_true_PDGID->push_back((*wDaughtIt)->getPDG());
		  m_true_index->push_back(index_true);
		  m_true_GenStatus->push_back((*wDaughtIt)->getGeneratorStatus());
		  m_true_numDaughters->push_back((*wDaughtIt)->getDaughters().size ());
		  m_true_numMothers->push_back((*wDaughtIt)->getParents().size());
		  m_true_m1_PDGID->push_back((*wDaughtIt)->getParents()[0]->getPDG());
		  m_true_m1_E->push_back((*wDaughtIt)->getParents()[0]->getEnergy());
		  //we check in the end for true electron and muon isolations
		  TLorentzVector wDLep(0,0,0,0);
		  int ind_wD=m_trueEnergy->size()-1;
		  wDLep.SetPxPyPzE((*wDaughtIt)->getMomentum()[0],(*wDaughtIt)->getMomentum()[1],(*wDaughtIt)->getMomentum()[2],(*wDaughtIt)->getEnergy());
		  for(int m_R =0; m_R< mcColl->getNumberOfElements(); m_R++){
		    MCParticle* mcp_R= dynamic_cast<MCParticle*>( mcColl->getElementAt(m_R) ) ;
		    if(mcp_R!=(*wDaughtIt)  &&  mcp_R->getGeneratorStatus()==1 &&((abs(mcp_R->getPDG())!=12)||(abs(mcp_R->getPDG())!=14)||(abs(mcp_R->getPDG())!=16)) ){
		      TLorentzVector temp(0,0,0,0);
		      temp.SetPxPyPzE(mcp_R->getMomentum()[0],mcp_R->getMomentum()[1],mcp_R->getMomentum()[2],mcp_R->getEnergy());
		      if(((wDLep.Angle(temp.Vect()))/M_PI*180.)<6.0){
			m_true_all_6_deg_E->push_back(mcp_R->getEnergy());
			m_true_all_6_deg_angle->push_back(wDLep.Angle(temp.Vect()));
			m_true_all_6_deg_index->push_back(ind_wD);
			m_true_all_6_deg_PDG->push_back(mcp_R->getPDG());
		      }
		    }
		  }
		}
	      }//loop over w daughter set
	      if(!W1_is_mu_el && !W2_is_mu_el){
		std::cout<<"none of the W's claims be leptonic WTF "<<W1_is_mu_el<<"/"<<W2_is_mu_el<<std::endl;
	      }
	    }//check if leptonic decays, including ALL tau (no check for decays yet)
	  }//W decays into leptons, tau,mu or electrons
	  //check taus first)	  
	}//should be a loop over W's
	
	//-------------------------------------------now the W loop for ttbar is done


	//after W's and tau decays from W's have been checked - turn to muons and electrons from hadron decays
	//check so for NOT for hadron decays into leptons, these are not isolated enough in ttbar to be of relevance
      }//loop over MCParticles is done
      if(W1_is_tau_had && W1_is_mu_el){
	std::cout<<"W el and tau had WTF"<<std::endl;
      }
      if(W1_is_tau_lep && W1_is_mu_el){
	std::cout<<"W el and tau lep WTF"<<std::endl;
      }
      if(W1_is_tau_lep && W1_is_tau_had){
	std::cout<<"W tau had and tau lep maybe"<<std::endl;
      }
      if(W2_is_tau_had && W2_is_mu_el){
	std::cout<<"W2 el and tau had WTF"<<std::endl;
      }
      if(W2_is_tau_lep && W2_is_mu_el){
	std::cout<<"W2 el and tau lep WTF"<<std::endl;
      }
      if(W2_is_tau_lep && W2_is_tau_had){
	std::cout<<"W2 tau had and tau lep maybe"<<std::endl;
      }
      if(W1_is_tau_had && W2_is_tau_had){
	m_ttbar_decay_mode=9;
      }else if (W1_is_tau_lep && W2_is_tau_lep){
	m_ttbar_decay_mode=7;
      }else if ((W1_is_tau_had && W2_is_tau_lep) || (W1_is_tau_lep && W2_is_tau_had)){
	m_ttbar_decay_mode=8;
      }else if (W1_is_mu_el && W2_is_mu_el){
	m_ttbar_decay_mode=2;
      }else if ((W1_is_tau_had && W2_is_mu_el) || (W1_is_mu_el && W2_is_tau_had)){
	m_ttbar_decay_mode=6;
      }else if ((W1_is_tau_lep && W2_is_mu_el) || (W1_is_mu_el && W2_is_tau_lep)){
	m_ttbar_decay_mode=5;
      }else if( W1_is_mu_el || W2_is_mu_el){
	m_ttbar_decay_mode=1;
      }else if (W1_is_tau_lep || W2_is_tau_lep){
	m_ttbar_decay_mode=3;
      }else if(W1_is_tau_had || W2_is_tau_had){
	m_ttbar_decay_mode=4;
      }else{
	m_ttbar_decay_mode=0;
      }
      //std::cout<<"ttbar decay case "<<m_ttbar_decay_mode<<" W1 l/t_l/t_h "<<W1_is_mu_el<<"/"<<W1_is_tau_lep<<"/"<<W1_is_tau_had<<" W2 l/t_l/t_h  " <<W2_is_mu_el<<"/"<<W2_is_tau_lep<<"/"<<W2_is_tau_had<<std::endl;
      //std::cout<<"before reco collection"<<std::endl;
      if(!m_runGenOnly){
	//std::cout<<"in reco collection"<<std::endl;
	LCCollection* recoparticlecol = NULL;
	// Alternativelly if you do not want Marlin to exit in case of a non-existing collection
	// use the following (commented out) code:
	//run on H to gamma gamma -> in case there are no tracks around linker will fail to produce output collection
	recoparticlecol = evt->getCollection(m_inputRECOParticleCollection) ;
	if(recoparticlecol!=NULL){
	  //PandoraCandidate loop
	  for(int i=0;i<recoparticlecol->getNumberOfElements();i++){
	    ReconstructedParticle* pandorapart = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(i));	
	    m_E_totPFO+=pandorapart->getEnergy();
	    m_px_totPFO+=pandorapart->getMomentum()[0];
	    m_py_totPFO+=pandorapart->getMomentum()[1];
	    m_pz_totPFO+=pandorapart->getMomentum()[2];
	    if(abs(pandorapart->getType())!=11 && abs(pandorapart->getType())!=13){
	      //std::cout<<"veto this particle "<<pandorapart->getType()<<"/"<<pandorapart->getEnergy()<<std::endl;
	      continue;
	    }
	    //std::cout<<"WOOHOOOO electron or muon found reco "<<pandorapart->getType()<<"/"<<pandorapart->getEnergy()<<std::endl;
	    m_reco_nTracks->push_back(pandorapart->getTracks().size());
	    float track_p=-1;
	    if(pandorapart->getTracks().size()>0){
	      m_reco_track0_pt->push_back(m_innerBField * m_const_a/fabs( pandorapart->getTracks()[0]->getOmega()));
	      track_p=m_innerBField * m_const_a/(fabs( pandorapart->getTracks()[0]->getOmega())*cos(atan(pandorapart->getTracks()[0]->getTanLambda())));
	      m_reco_track0_p->push_back(track_p);
	      m_reco_track0_chi2OverNdof->push_back(pandorapart->getTracks()[0]->getChi2()/(float)pandorapart->getTracks()[0]->getNdf()); 
	      m_reco_track0_nHits->push_back(pandorapart->getTracks()[0]->getTrackerHits().size());
	    }else{
	      m_reco_track0_pt->push_back(-1.);
	      m_reco_track0_p->push_back(-1.);
	      m_reco_track0_nHits->push_back(-1);
	      m_reco_track0_chi2OverNdof->push_back(-1.);
	    }
	    m_reco_nClusters->push_back(pandorapart->getClusters().size());
	    float energy_sum_clusters=0;
	    if(pandorapart->getClusters().size()>0){
	      for(unsigned int c=0;c<pandorapart->getClusters().size();c++){
		energy_sum_clusters+=pandorapart->getClusters()[c]->getEnergy();
	      }
	      m_reco_clusters_energy->push_back(energy_sum_clusters);
	    }	  	  
	    //Check if in barrel
	    double cosTheta = pandorapart->getMomentum()[2]/sqrt(pandorapart->getMomentum()[0]*pandorapart->getMomentum()[0]+pandorapart->getMomentum()[1]*pandorapart->getMomentum()[1]+pandorapart->getMomentum()[2]*pandorapart->getMomentum()[2]);
	    m_recoEnergy->push_back(pandorapart->getEnergy());
	    m_reco_Px->push_back(pandorapart->getMomentum()[0]);
	    m_reco_Py->push_back(pandorapart->getMomentum()[1]);
	    m_reco_Pz->push_back( pandorapart->getMomentum()[2]);
	    m_reco_CosTheta->push_back(cosTheta);    
	    //m_true_Phi=atan2(m_true_Py,m_true_Px);
	    m_reco_PDGID->push_back(pandorapart->getType());
	    float reco_px=pandorapart->getMomentum()[0];
	    float reco_py=pandorapart->getMomentum()[1];
	    //m_reco_Phi->push_back(atan2(pandorapart->getMomentum()[1],pandorapart->getMomentum()[0]));
	    m_reco_Phi->push_back(atan2(reco_py,reco_px));
	    m_reco_Charge->push_back(pandorapart->getCharge());
	    float En_ECAL_Barrel=0;
	    float En_ECAL_Endcap=0;
	    float En_ECAL_else=0;
	    float En_HCAL_Barrel=0;
	    float En_HCAL_Endcap=0;
	    float En_HCAL_else=0;
	    for(unsigned int j=0;j<pandorapart->getClusters().size();j++){
	      for(unsigned int l=0; l<pandorapart->getClusters()[j]->getCalorimeterHits().size();l++){
		//loop over hits for logarithmic reweighting
		// constant is changed by hand --> why where and how was that decided
		//double logWeight = std::max ( 0.0 , 5.5 + std::log(pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy()/totEnergy));
		const CHT cht=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getType();		
		//int types=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getType();
		//static const int fCaloType =     1 ;
		//static const int fCaloID   =    10 ;
		//static const int fLayout   =  1000 ;
		//static const int fLayer    = 10000 ;
		//int caloLayer=types/fLayer;	      
		if(cht.is(CHT::ecal)){//ecal		
		  if(cht.is(CHT::barrel)){
		    En_ECAL_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  }else if (cht.is(CHT::endcap)){
		    En_ECAL_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  }else{
		    En_ECAL_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  }
		}else if(cht.is(CHT::hcal)){//h-cal
		  if(cht.is(CHT::barrel)){
		    En_HCAL_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  }else if(cht.is(CHT::endcap)){
		    En_HCAL_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  }else{
		    En_HCAL_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		  }
		}	
	      }//hit loop
	    }//cluster loop
	    m_reco_E_EB->push_back(En_ECAL_Barrel);
	    m_reco_E_EE->push_back(En_ECAL_Endcap);
	    m_reco_E_EO->push_back(En_ECAL_else);
	    m_reco_E_HB->push_back(En_HCAL_Barrel);
	    m_reco_E_HE->push_back(En_HCAL_Endcap);
	    m_reco_E_HO->push_back(En_HCAL_else);	  
	    //pf particle loop done
	    //we check in the end for true electron and muon isolations
	    TLorentzVector rLep(0,0,0,0);
	    int ind_rLep=m_recoEnergy->size()-1;
	    rLep.SetPxPyPzE(pandorapart->getMomentum()[0],pandorapart->getMomentum()[1],pandorapart->getMomentum()[2],pandorapart->getEnergy());
	    //double isoE=0;
	    for(int reco_R =0; reco_R< recoparticlecol->getNumberOfElements(); reco_R++){
	      if(reco_R==i){
		//std::cout<<"particle checked against itself"<<std::endl;
		continue;
	      }
	      ReconstructedParticle* pand_R= dynamic_cast<ReconstructedParticle*>( recoparticlecol->getElementAt(reco_R) ) ;
	      if(pand_R!=pandorapart){
		TLorentzVector temp(0,0,0,0);
		temp.SetPxPyPzE(pand_R->getMomentum()[0],pand_R->getMomentum()[1],pand_R->getMomentum()[2],pand_R->getEnergy());
		if(((rLep.Angle(temp.Vect()))/M_PI*180.)<6.0){
		  //std::cout<<"particle "<<reco_R<<" "<<rLep.Angle(temp.Vect())<<"/"<<rLep.Angle(temp.Vect())/M_PI*180<<std::endl;
		  m_reco_all_6_deg_E->push_back(pand_R->getEnergy());
		  m_reco_all_6_deg_angle->push_back(rLep.Angle(temp.Vect()));
		  m_reco_all_6_deg_index->push_back(ind_rLep);
		  m_reco_all_6_deg_PDG->push_back(pand_R->getType());
		  //isoE+=pand_R->getEnergy();
		}
	      }//else{
	      //std::cout<<"particle checked against itself"<<std::endl;
	      //}
	    }
	    //std::cout<<"rel iso at end "<<isoE/pandorapart->getEnergy()<<std::endl;
	  }
	}//reco collection available
      }//do gen only
 
    m_outputTree->Fill();
	
}
      
void TTBarStudy::fillStableDaughterSet(MCParticle* mcp, std::set<MCParticle*> &stableDaughterSet){
  if(mcp->getGeneratorStatus()==1){
    stableDaughterSet.insert(mcp);
  }else if (mcp->getGeneratorStatus()==0){
    return;
  }
  for(unsigned int d=0;d<mcp->getDaughters().size();d++){
    fillStableDaughterSet(mcp->getDaughters()[d], stableDaughterSet);
  }
}




void TTBarStudy::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    
    streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

void TTBarStudy::check(LCEvent*){
}

void TTBarStudy::end(){
 
  m_rootFile->Write();
  m_rootFile->Close();

}
