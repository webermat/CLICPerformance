#include "PionStudy.h"
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHitPlane.h>

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetectorSelector.h"

#include "UTIL/LCRelationNavigator.h"
#include "UTIL/ILDConf.h"
#include "CalorimeterHitType.h"

#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <marlinutil/GeometryUtil.h>

#include "UTIL/LCRelationNavigator.h"

using namespace lcio ;
using namespace marlin ;

using dd4hep::DetType;

PionStudy aPionStudy;

//dd4hep::rec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag=0);

PionStudy::PionStudy() : Processor("PionStudy") {
    
    // modify processor description
    _description = "PionStudy calculates properties of calorimeter showers" ;
   
    registerInputCollection( LCIO::TRACK,
			     "TrackCollectionName",
			     "Name of th eTrack Collection",
			     m_inputTrackCollection,
			     std::string("SiTracks")
			     );


    registerInputCollection( LCIO::SIMCALORIMETERHIT,
			     "SimCaloHitCollectionName",
			     "Name of the SimCaloHitCollection (Muon Barrel)",
			     m_inputSimCaloHitCollectionName,
			     std::string("YokeBarrelCollection")
			     );


    registerInputCollection( LCIO::MCPARTICLE,
                            "MCParticleCollectionName",
                            "Name of the MCParticle input collection",
                            m_inputMCParticleCollection,
                            std::string("MCParticle"));

    registerInputCollection( LCIO::CLUSTER,
                            "ClusterCollectionName",
                            "Name of the Cluster input collection",
                            m_inputClusterCollection,
                            std::string("PandoraClusters"));
    
    registerProcessorParameter( "OutputRootFileName",
                                "ROOT File name to collect plots",
                                m_rootFileName,
                                std::string("showerStudy.root"));
  

    registerProcessorParameter("runTauMode",
			       "tau specific running mode",
			       m_runTauMode,
			       bool("false"));

   registerProcessorParameter("runPhotonMode",
			       "photon specific running mode",
			       m_runPhotonMode,
			       bool("false"));

    registerProcessorParameter("fillClusterHits",
			       "fill cluster hits (only in non tauMode",
			       m_fillClusterHits,
			       bool("false"));

    registerProcessorParameter("fillTrackBranches",
			       "fill track branches",
			       m_fillTrackBranches,
			       bool("true"));

    registerProcessorParameter("fillSimHits",
			       "fill sim hits",
			       m_fillSimHits,
			       bool("true"));

    registerProcessorParameter("fillClusterBranches",
			       "fill cluster branches",
			       m_fillClusterBranches,
			       bool("true"));

    
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

    registerInputCollection( LCIO::CALORIMETERHIT,
			     "MuonHitCollectionName" , 
			     "Name of the Muon hit collection"  ,
			     m_inputMuonHitCollectionName,
			     std::string("MUON")
			     );

}


void PionStudy::init() {

  const dd4hep::rec::LayeredCalorimeterData * eCalEndcapExtension=       MarlinUtil::getLayeredCalorimeterData( ( DetType::CALORIMETER | DetType::ELECTROMAGNETIC | DetType::ENDCAP),
                                             ( DetType::AUXILIARY | DetType::FORWARD ) );


  const dd4hep::rec::LayeredCalorimeterData * eCalBarrelExtension= MarlinUtil::getLayeredCalorimeterData( ( DetType::CALORIMETER | DetType::ELECTROMAGNETIC | DetType::BARREL),
                                             ( DetType::AUXILIARY | DetType::FORWARD ) );

  //extent 0/1/2.3 : inner/outer R inner/outerZ
    m_ECAL_endcapZ_min=static_cast<float>(eCalEndcapExtension->extent[2]/dd4hep::mm);
    m_ECAL_endcapR_min=static_cast<float>(eCalEndcapExtension->extent[0]/dd4hep::mm);
    m_ECAL_barrelR_min=static_cast<float>(eCalBarrelExtension->extent[0]/dd4hep::mm);


  const dd4hep::rec::LayeredCalorimeterData * eCalRingExtension= MarlinUtil::getLayeredCalorimeterData( ( DetType::CALORIMETER | DetType::ELECTROMAGNETIC | DetType::ENDCAP | dd4hep::DetType::AUXILIARY ),
                                             ( DetType::FORWARD ) );

  m_ECAL_barrelZ_max=eCalBarrelExtension->extent[3]/dd4hep::mm;
  m_ECAL_ringZ_min=eCalRingExtension->extent[2]/dd4hep::mm;
  m_ECAL_ringZ_max=eCalRingExtension->extent[3]/dd4hep::mm;
  m_ECAL_endcapZ_min=eCalEndcapExtension->extent[2]/dd4hep::mm;
  m_ECAL_endcapZ_max=eCalEndcapExtension->extent[3]/dd4hep::mm;


  const dd4hep::rec::LayeredCalorimeterData * hcalEndcapExtension= MarlinUtil::getLayeredCalorimeterData( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) );

  const dd4hep::rec::LayeredCalorimeterData * hcalBarrelExtension= MarlinUtil::getLayeredCalorimeterData( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) );

  const dd4hep::rec::LayeredCalorimeterData * hcalRingExtension= MarlinUtil::getLayeredCalorimeterData( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP | dd4hep::DetType::AUXILIARY ),( dd4hep::DetType::FORWARD ) );

  m_HCAL_barrelZ_max=hcalBarrelExtension->extent[3]/dd4hep::mm;
  m_HCAL_ringZ_min=hcalRingExtension->extent[2]/dd4hep::mm;
  m_HCAL_ringZ_max=hcalRingExtension->extent[3]/dd4hep::mm;
  m_HCAL_endcapZ_min=hcalEndcapExtension->extent[2]/dd4hep::mm;
  m_HCAL_endcapZ_max=hcalEndcapExtension->extent[3]/dd4hep::mm;


  const dd4hep::rec::LayeredCalorimeterData * muonEndcapExtension= MarlinUtil::getLayeredCalorimeterData( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::ENDCAP), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) );

  const dd4hep::rec::LayeredCalorimeterData * muonBarrelExtension= MarlinUtil::getLayeredCalorimeterData( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::BARREL), ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) );

  m_MUON_barrelZ_max=muonBarrelExtension->extent[3]/dd4hep::mm;
  m_MUON_endcapZ_min=muonEndcapExtension->extent[2]/dd4hep::mm;
  m_MUON_endcapZ_max=muonEndcapExtension->extent[3]/dd4hep::mm;

  std::cout<<"MUON numbers "<<m_MUON_barrelZ_max<<"/"<<m_MUON_endcapZ_min<<"/"<<m_MUON_endcapZ_max<<std::endl;
   

  dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();
  
  const double position[3]={0,0,0}; // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3]={0,0,0}; // initialise object to hold magnetic field
  mainDetector.field().magneticField(position,magneticFieldVector); // get the magnetic field vector from dd4hep
    
  m_innerBField = magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)
  //m_innerBField = 4.0;//magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)
  //constant takes into account that omega is in untis of mm
  m_const_a=2.99792e-4;

  m_trackerInnerR=0.1; ///FIXME! CLIC-specific: Inner radius was set to 0 for SiD-type detectors --> also for CLIC it is NOT set
  m_trackerOuterR=mainDetector.constantAsDouble("tracker_region_rmax")/dd4hep::mm;
  m_trackerZmax=mainDetector.constantAsDouble("tracker_region_zmax")/dd4hep::mm;

  const std::vector< dd4hep::DetElement>& barrelDets = dd4hep::DetectorSelector(mainDetector).detectors( ( dd4hep::DetType::TRACKER | dd4hep::DetType::BARREL )) ;
  
  m_barrelTrackerLayers=0;
  
  for (std::vector< dd4hep::DetElement>::const_iterator iter = barrelDets.begin(), iterEnd = barrelDets.end();iter != iterEnd; ++iter){
    try
      {
	dd4hep::rec::ZPlanarData * theExtension = 0;
	
	const dd4hep::DetElement& theDetector = *iter;
	theExtension = theDetector.extension<dd4hep::rec::ZPlanarData>();
	
	unsigned int N = theExtension->layers.size();
	m_barrelTrackerLayers=m_barrelTrackerLayers+N;
	
	streamlog_out( DEBUG2 ) << " Adding layers for barrel tracker from DD4hep for "<< theDetector.name()<< "- n layers: " << N<< " sum up to now: "<<m_barrelTrackerLayers<<std::endl;
      } catch (std::runtime_error &exception){
      
      streamlog_out(WARNING) << "DDTrackCreatorCLIC exception during Barrel Tracker layer sum for "<<const_cast<dd4hep::DetElement&>(*iter).name()<<" : " << exception.what() << std::endl;
    }
  }
  

  eventcount=0;
  
    // Print the initial parameters
    printParameters() ;

    // Reset counters
    m_runNumber = 0 ;
    m_eventNumber = 0 ;
    

    m_rootFile = new TFile(m_rootFileName.c_str(),"RECREATE");
    
    //implement cluster track distance using hits -> check pandoraPFA in order how to do it properly

    m_outputTree = new TTree("showerData","showerData");

    m_trueEnergy = new std::vector<float>();
    m_true_Px = new std::vector<float>();
    m_true_Py = new std::vector<float>();
    m_true_Pz = new std::vector<float>();
    m_true_x = new std::vector<float>();
    m_true_y = new std::vector<float>();
    m_true_z = new std::vector<float>();
    m_true_CosTheta = new std::vector<float>();
    m_true_Theta = new std::vector<float>();
    m_true_Phi = new std::vector<float>();
    m_true_PDGID = new std::vector<int>();
    m_true_GenStatus = new std::vector<int>();
    m_true_numDaughters = new std::vector<int>();
    m_true_decayTrackerCalo = new std::vector<int>();//1 for tracker, 2 for calo. 3 for both, 4 for leaving dector, 0 for else
    m_true_motherDecayTrackerCalo = new std::vector<int>();//1 for tracker, 2 for calo
    m_true_numMothers = new std::vector<int>();
    m_true_m1_PDGID = new std::vector<int>();
    m_true_m2_PDGID = new std::vector<int>();
    m_true_m1_status = new std::vector<int>();
    m_true_m2_status = new std::vector<int>();
    m_true_m1_E = new std::vector<float>();
    m_true_m2_E = new std::vector<float>();
    m_true_index = new std::vector<int>();
    if(m_runPhotonMode){
      m_true_conv_e1_Px = new std::vector<float>();
      m_true_conv_e1_Py = new std::vector<float>();
      m_true_conv_e1_Pz = new std::vector<float>();
      m_true_conv_e1_E =  new std::vector<float>();
      m_true_conv_e1_PDGID = new std::vector<int>();
      m_true_conv_e2_Px = new std::vector<float>();
      m_true_conv_e2_Py = new std::vector<float>();
      m_true_conv_e2_Pz = new std::vector<float>();
      m_true_conv_e2_E =  new std::vector<float>();
      m_true_conv_e2_PDGID = new std::vector<int>();
      m_true_conv_Ph_VecInd = new std::vector<int>();
      m_true_conv_Vtx_x = new std::vector<float>();
      m_true_conv_Vtx_y = new std::vector<float>();
      m_true_conv_Vtx_z = new std::vector<float>();
    }

    if(m_runTauMode){
      m_true_tauDaughter_Energy = new std::vector<float>();
      m_true_tauDaughter_Px = new std::vector<float>();
      m_true_tauDaughter_Py = new std::vector<float>();
      m_true_tauDaughter_Pz = new std::vector<float>();
      m_true_tauDaughter_PDGID = new std::vector<int>();
      m_true_tauDaughter_Charge = new std::vector<int>();
      m_true_tauDaughter_tauIndex = new std::vector<int>();
      m_true_tauDaughter_status = new std::vector<int>();
      m_true_tauDaughter_motherPDGID  = new std::vector<int>();
      m_true_tauDaughter_motherEnergy = new std::vector<float>();

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
    }else{
      m_cluster_x = new std::vector<float>();
      m_cluster_y = new std::vector<float>();
      m_cluster_z = new std::vector<float>();
      m_cluster_dir_x = new std::vector<float>();
      m_cluster_dir_y = new std::vector<float>();
      m_cluster_dir_z = new std::vector<float>();
      m_cluster_energy = new std::vector<float>();
      m_cluster_energyError = new std::vector<float>();
      m_cluster_ECAL_energy = new std::vector<float>();
      m_cluster_HCAL_energy = new std::vector<float>();
      m_cluster_nHits_ECAL = new std::vector<int>();
      m_cluster_nHits_HCAL = new std::vector<int>();
      m_cluster_nHits_MUON = new std::vector<int>();
      m_cluster_iPhi = new std::vector<float>();
      m_cluster_iTheta = new std::vector<float>();
      m_cluster_track_PosAngle_min = new std::vector<float>();
      m_cluster_track_DirAngle = new std::vector<float>();
      m_cluster_track_minPFADist_cluster_EOverP = new std::vector<float>();
      m_cluster_track_minPFADist_cluster_Angle = new std::vector<float>();
      m_cluster_track_minPFADist_cluster_TCDistance = new std::vector<float>();
      m_cluster_E_EB = new std::vector<float>();
      m_cluster_E_EE = new std::vector<float>();
      m_cluster_E_EO = new std::vector<float>();
      m_cluster_E_HB = new std::vector<float>();
      m_cluster_E_HE = new std::vector<float>();
      m_cluster_E_HO = new std::vector<float>();
      m_cluster_E_MU = new std::vector<float>();

      if(m_fillClusterHits){
	m_cluster_hit_x = new std::vector<float>();
	m_cluster_hit_y = new std::vector<float>();
	m_cluster_hit_z = new std::vector<float>();
	m_cluster_hit_E = new std::vector<float>();
        m_cluster_hit_index = new std::vector<int>();
	m_cluster_hit_type = new std::vector<int>();

	m_no_cluster_hit_x = new std::vector<float>();
	m_no_cluster_hit_y = new std::vector<float>();
	m_no_cluster_hit_z = new std::vector<float>();
	m_no_cluster_hit_E = new std::vector<float>();
        //m_no_cluster_hit_index = new std::vector<int>();
	//m_no_cluster_hit_type = new std::vector<int>();
      }

      if(m_fillSimHits){
	m_muon_b_simhit_x = new std::vector<float>();
	m_muon_b_simhit_y = new std::vector<float>();
	m_muon_b_simhit_z = new std::vector<float>();
	m_muon_b_simhit_E = new std::vector<float>();
      }
    }
     
    if(m_fillTrackBranches){
      m_track_d0 = new std::vector<float>();
      m_track_z0 = new std::vector<float>();
      m_track_phi0 = new std::vector<float>();
      m_track_ndf = new std::vector<int>();
      m_track_nHits = new std::vector<int>();
      m_track_nExpectedTrackerHits = new std::vector<float>();
      m_track_nTrackerBarrelHits = new std::vector<int>(); //tracker+vertex
      m_track_nTrackerEndcapHits = new std::vector<int>();//tracker+vertex
      m_track_nVertexBarrelHits = new std::vector<int>();//vertex only
      m_track_nVertexEndcapHits = new std::vector<int>();//vertex only
      m_track_chi2 = new std::vector<float>();
      m_track_sigmaPOverP = new std::vector<float>();
      m_track_x_atIP = new std::vector<float>();
      m_track_y_atIP = new std::vector<float>();
      m_track_z_atIP = new std::vector<float>();
      m_track_Phi_atIP = new std::vector<float>();
      m_track_Theta_atIP = new std::vector<float>();
      m_track_pt_atIP = new std::vector<float>();
      m_track_p_atIP = new std::vector<float>();
      m_track_x_atCalo = new std::vector<float>();
      m_track_y_atCalo = new std::vector<float>();
      m_track_z_atCalo = new std::vector<float>();
      m_track_zMin = new std::vector<float>();
      m_track_x_innermostHit = new std::vector<float>();
      m_track_y_innermostHit = new std::vector<float>();
      m_track_z_innermostHit = new std::vector<float>();
      m_track_pt_innermostHit = new std::vector<float>();
      m_track_p_innermostHit = new std::vector<float>();
      m_track_x_outermostRHit = new std::vector<float>();
      m_track_y_outermostRHit = new std::vector<float>();
      m_track_z_outermostRHit = new std::vector<float>();
      m_track_pt_outermostRHit = new std::vector<float>();
      m_track_p_outermostRHit = new std::vector<float>();
      m_track_x_outermostZHit = new std::vector<float>();
      m_track_y_outermostZHit = new std::vector<float>();
      m_track_z_outermostZHit = new std::vector<float>();
      m_track_px_atCalo = new std::vector<float>();
      m_track_py_atCalo = new std::vector<float>();
      m_track_pz_atCalo = new std::vector<float>();
      m_track_Phi_atCalo = new std::vector<float>();
      m_track_Theta_atCalo = new std::vector<float>();
      m_track_pt_atCalo = new std::vector<float>();
      m_track_p_atCalo = new std::vector<float>();
      m_track_cluster_PosAngle_min_atCalo = new std::vector<float>();
      m_track_minDist_cluster_EOverP = new std::vector<float>();
      m_track_cluster_minPFADist_clusterEMEnergy = new std::vector<float>();
      m_track_cluster_minPFADist_clusterHadEnergy = new std::vector<float>();
      m_track_cluster_minPFADist_EOverP = new std::vector<float>();
      m_track_cluster_minPFAdistance_atCalo = new std::vector<float>();
      m_track_cluster_min_parPFAdistance_atCalo = new std::vector<float>();
      m_track_cluster_DirAngle_minPFAdistance_atCalo = new std::vector<float>();
      m_track_cluster_minPFADist_EOverP_tinyHad = new std::vector<float>();
      m_track_cluster_minPFADist_clusterEMEnergy_tinyHad = new std::vector<float>();
      m_track_cluster_minPFADist_clusterHadEnergy_tinyHad = new std::vector<float>();
      m_track_cluster_minPFAdistance_atCalo_tinyHad = new std::vector<float>();
      m_track_cluster_min_parPFAdistance_atCalo_tinyHad = new std::vector<float>();
      m_track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad = new std::vector<float>();
    }

    m_recoEnergy = new std::vector<float>();
    m_recohitEnergyCellCorrection = new std::vector<float>();
    m_reco_Px = new std::vector<float>();
    m_reco_Py = new std::vector<float>();
    m_reco_Pz = new std::vector<float>();
    m_reco_CosTheta = new std::vector<float>();
    m_reco_Theta = new std::vector<float>();
    m_reco_Phi = new std::vector<float>();
    m_reco_Charge = new std::vector<int>();

    m_reco_CosTheta_logERW = new std::vector<float>();
    m_reco_Theta_logERW = new std::vector<float>();
    m_reco_Phi_logERW = new std::vector<float>();
    m_reco_logE_tot = new std::vector<float>();
    m_reco_E_totRW = new std::vector<float>();

    m_reco_x_logERW = new std::vector<float>();
    m_reco_y_logERW = new std::vector<float>();
    m_reco_z_logERW = new std::vector<float>();

    m_reco_nTracks = new std::vector<int>();
    m_reco_track0_pt = new std::vector<float>();
    m_reco_track0_p = new std::vector<float>();
    m_reco_track0_nHits = new std::vector<int>();
    m_reco_track0_chi2OverNdof = new std::vector<float>();
    m_reco_nClusters = new std::vector<int>();
    m_reco_clusters_energy = new std::vector<float>();
    m_reco_cluster0_energy = new std::vector<float>();
    m_reco_cluster0_iPhi = new std::vector<float>();
    m_reco_cluster0_iTheta = new std::vector<float>();
    m_reco_cluster0_energyError = new std::vector<float>();
    m_reco_EClustersOverPTrack = new std::vector<float>();
    m_reco_EOverP_PFA = new std::vector<float>();
    m_reco_PDGID = new std::vector<int>();
    m_reco_E_EB = new std::vector<float>();
    m_reco_E_EE = new std::vector<float>();
    m_reco_E_EO = new std::vector<float>();
    m_reco_E_HB = new std::vector<float>();
    m_reco_E_HE = new std::vector<float>();
    m_reco_E_HO = new std::vector<float>();
    m_reco_E_MB = new std::vector<float>();
    m_reco_E_ME = new std::vector<float>();
    m_reco_E_MO = new std::vector<float>();
    m_reco_firstLayerECAL = new std::vector<int>();
    m_reco_lastLayerECAL = new std::vector<int>();
    m_reco_nhitsEB = new std::vector<int>();
    m_reco_nhitsEE = new std::vector<int>();
    m_reco_nhitsEO = new std::vector<int>();
    m_reco_firstLayerHCAL = new std::vector<int>();
    m_reco_lastLayerHCAL = new std::vector<int>();
    m_reco_nhitsHB = new std::vector<int>();
    m_reco_nhitsHE = new std::vector<int>();
    m_reco_nhitsHO = new std::vector<int>();
    m_reco_nhitsMB = new std::vector<int>();
    m_reco_nhitsME = new std::vector<int>();
    m_reco_nhitsMO = new std::vector<int>();

    m_trueEnergy->clear();
    m_true_Px->clear();
    m_true_Py->clear();
    m_true_Pz->clear();
    m_true_x->clear();
    m_true_y->clear();
    m_true_z->clear();
    m_true_index->clear();
    m_true_CosTheta->clear();
    m_true_Theta->clear();
    m_true_Phi->clear();
    m_true_PDGID->clear();
    m_true_GenStatus->clear();
    m_true_numDaughters->clear();
    m_true_decayTrackerCalo->clear();
    m_true_motherDecayTrackerCalo->clear();
    m_true_numMothers->clear();
    m_true_m1_PDGID->clear();
    m_true_m2_PDGID->clear();
    m_true_m1_status->clear();
    m_true_m2_status->clear();
    m_true_m1_E->clear();
    m_true_m2_E->clear();

    if(m_runPhotonMode){
      m_true_conv_e1_Px->clear();
      m_true_conv_e1_Py->clear();
      m_true_conv_e1_Pz->clear();
      m_true_conv_e1_E->clear();
      m_true_conv_e1_PDGID->clear();
      m_true_conv_e2_Px->clear();
      m_true_conv_e2_Py->clear();
      m_true_conv_e2_Pz->clear();
      m_true_conv_e2_E->clear();
      m_true_conv_e2_PDGID->clear();
      m_true_conv_Ph_VecInd->clear();
      m_true_conv_Vtx_x->clear();
      m_true_conv_Vtx_y->clear();
      m_true_conv_Vtx_z->clear();
    }

    if(m_runTauMode){
      m_true_tauDaughter_Energy->clear();
      m_true_tauDaughter_Px->clear();
      m_true_tauDaughter_Py->clear();
      m_true_tauDaughter_Pz->clear();
      m_true_tauDaughter_PDGID->clear();
      m_true_tauDaughter_Charge->clear();
      m_true_tauDaughter_tauIndex->clear();
      m_true_tauDaughter_status->clear();
      m_true_tauDaughter_motherPDGID  ->clear();
      m_true_tauDaughter_motherEnergy ->clear();

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
    }else{
      m_cluster_x->clear();
      m_cluster_y->clear();
      m_cluster_z->clear();
      m_cluster_dir_x->clear();
      m_cluster_dir_y->clear();
      m_cluster_dir_z->clear();
      m_cluster_energy->clear();
      m_cluster_energyError->clear();
      m_cluster_ECAL_energy->clear();
      m_cluster_HCAL_energy->clear();
      m_cluster_nHits_ECAL->clear();
      m_cluster_nHits_HCAL->clear();
      m_cluster_nHits_MUON->clear();
      m_cluster_E_EB->clear();
      m_cluster_E_EE->clear();
      m_cluster_E_EO->clear();
      m_cluster_E_HB->clear();
      m_cluster_E_HE->clear();
      m_cluster_E_HO->clear();
      m_cluster_E_MU->clear();
      m_cluster_iPhi->clear();
      m_cluster_iTheta->clear();
      m_cluster_track_PosAngle_min->clear();
      m_cluster_track_DirAngle->clear();
      m_cluster_track_minPFADist_cluster_EOverP->clear();
      m_cluster_track_minPFADist_cluster_Angle->clear();
      m_cluster_track_minPFADist_cluster_TCDistance->clear();

      if(m_fillClusterHits){
	m_cluster_hit_x->clear();
	m_cluster_hit_y->clear();
	m_cluster_hit_z->clear();
	m_cluster_hit_E->clear();
        m_cluster_hit_index->clear();
	m_cluster_hit_type->clear();

	m_no_cluster_hit_x->clear();
	m_no_cluster_hit_y->clear();
	m_no_cluster_hit_z->clear();
	m_no_cluster_hit_E->clear();
        //m_no_cluster_hit_index->clear();
	//m_no_cluster_hit_type->clear();
      }

      if(m_fillSimHits){
	m_muon_b_simhit_x->clear();
	m_muon_b_simhit_y->clear();
	m_muon_b_simhit_z->clear();
	m_muon_b_simhit_E->clear();
      }
    }
    if( m_fillTrackBranches){
      m_track_d0->clear();
      m_track_z0->clear();
      m_track_phi0->clear();
      m_track_ndf->clear();
      m_track_nHits->clear();
      m_track_nExpectedTrackerHits->clear();
      m_track_nTrackerBarrelHits->clear();
      m_track_nTrackerEndcapHits->clear();
      m_track_nVertexBarrelHits->clear();
      m_track_nVertexEndcapHits->clear();
      m_track_chi2->clear();
      m_track_sigmaPOverP->clear();
      m_track_x_atIP->clear();
      m_track_y_atIP->clear();
      m_track_z_atIP->clear();
      m_track_Phi_atIP->clear();
      m_track_Theta_atIP->clear();
      m_track_pt_atIP->clear();
      m_track_p_atIP->clear();
      m_track_x_atCalo->clear();
      m_track_y_atCalo->clear();
      m_track_z_atCalo->clear();
      m_track_zMin->clear();
      m_track_x_innermostHit->clear();
      m_track_y_innermostHit->clear();
      m_track_z_innermostHit->clear();
      m_track_p_innermostHit->clear();
      m_track_pt_innermostHit->clear();
      m_track_x_outermostRHit->clear();
      m_track_y_outermostRHit->clear();
      m_track_z_outermostRHit->clear();
      m_track_pt_outermostRHit->clear();
      m_track_p_outermostRHit->clear();
      m_track_x_outermostZHit->clear();
      m_track_y_outermostZHit->clear();
      m_track_z_outermostZHit->clear();
      m_track_px_atCalo->clear();
      m_track_px_atCalo->clear();
      m_track_px_atCalo->clear();
      m_track_Phi_atCalo->clear();
      m_track_Theta_atCalo->clear();
      m_track_pt_atCalo->clear();
      m_track_p_atCalo->clear();
      m_track_cluster_PosAngle_min_atCalo->clear();
      m_track_minDist_cluster_EOverP->clear();
      m_track_cluster_minPFADist_EOverP->clear();
      m_track_cluster_minPFADist_clusterEMEnergy->clear();
      m_track_cluster_minPFADist_clusterHadEnergy->clear();
      m_track_cluster_minPFAdistance_atCalo->clear();
      m_track_cluster_min_parPFAdistance_atCalo->clear();
      m_track_cluster_DirAngle_minPFAdistance_atCalo->clear();
      m_track_cluster_minPFADist_EOverP_tinyHad->clear();
      m_track_cluster_minPFADist_clusterEMEnergy_tinyHad->clear();
      m_track_cluster_minPFADist_clusterHadEnergy_tinyHad->clear();
      m_track_cluster_minPFAdistance_atCalo_tinyHad->clear();
      m_track_cluster_min_parPFAdistance_atCalo_tinyHad->clear();
      m_track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad->clear();
    }
    m_recoEnergy->clear();
    m_recohitEnergyCellCorrection->clear();
    m_reco_Px->clear();
    m_reco_Py->clear();
    m_reco_Pz->clear();
    m_reco_CosTheta->clear();
    m_reco_Theta->clear();
    m_reco_Phi->clear();
    m_reco_Charge->clear();
    m_reco_PDGID->clear();
    m_reco_CosTheta_logERW->clear();
    m_reco_Theta_logERW->clear();
    m_reco_Phi_logERW->clear();
    m_reco_x_logERW->clear();
    m_reco_y_logERW->clear();
    m_reco_z_logERW->clear();
    m_reco_logE_tot->clear();
    m_reco_E_totRW->clear();

    m_reco_nTracks->clear();
    m_reco_track0_pt->clear();
    m_reco_track0_p->clear();
    m_reco_track0_nHits->clear();
    m_reco_track0_chi2OverNdof->clear();
    m_reco_nClusters->clear();
    m_reco_clusters_energy->clear();
    m_reco_cluster0_energy->clear();
    m_reco_cluster0_iPhi->clear();
    m_reco_cluster0_iTheta->clear();
    m_reco_cluster0_energyError->clear();
    m_reco_EClustersOverPTrack->clear();
    m_reco_EOverP_PFA->clear();
    m_reco_E_EB->clear();
    m_reco_E_EE->clear();
    m_reco_E_EO->clear();
    m_reco_E_HB->clear();
    m_reco_E_HE->clear();
    m_reco_E_HO->clear();
    m_reco_E_MB->clear();
    m_reco_E_ME->clear();
    m_reco_E_MO->clear();
    m_reco_firstLayerECAL->clear();
    m_reco_lastLayerECAL->clear();
    m_reco_nhitsEB->clear();
    m_reco_nhitsEE->clear();
    m_reco_nhitsEO->clear();
    m_reco_firstLayerHCAL->clear();
    m_reco_lastLayerHCAL->clear();
    m_reco_nhitsHB->clear();
    m_reco_nhitsHE->clear();
    m_reco_nhitsHO->clear();
    m_reco_nhitsMB->clear();
    m_reco_nhitsME->clear();
    m_reco_nhitsMO->clear();


    m_outputTree->Branch("Z_mcE",&m_Z_mcE,"Z_mcE/F");
    m_outputTree->Branch("runNumber",&m_runNumber,"runNumber/I");
    m_outputTree->Branch("eventNumber",&m_eventNumber,"eventNumber/I");
    m_outputTree->Branch("Z_mcNDaughter",&m_Z_mcNDaughter,"Z_mcNDaughter/I");

    m_outputTree->Branch("d1_mcPDGID",&m_d1_mcPDGID,"d1_mcPDGID/I");
    m_outputTree->Branch("d1_mcE",&m_d1_mcE,"d1_mcE/F");
    m_outputTree->Branch("d1_mcPx",&m_d1_mcPx,"d1_mcPx/F");
    m_outputTree->Branch("d1_mcPy",&m_d1_mcPy,"d1_mcPy/F");
    m_outputTree->Branch("d1_mcPx",&m_d1_mcPz,"d1_mcPz/F");
    //m_outputTree->Branch("d1_mcE",&m_d1_mcE,"d1_mcE/F");
    //m_outputTree->Branch("d1_mcMass",&m_d1_mcMass,"d1_mcMass/F");
    //m_outputTree->Branch("d1_mcPhi",&m_d1_mcPhi,"d1_mcPhi/F");
    //m_outputTree->Branch("d1_mcTheta",&m_d1_mcTheta,"d1_mcTheta/F");
    m_outputTree->Branch("d1_mcCosTheta",&m_d1_mcCosTheta,"d1_mcCosTheta/F");

    m_outputTree->Branch("d2_mcPDGID",&m_d2_mcPDGID,"d2_mcPDGID/I");
    m_outputTree->Branch("d2_mcE",&m_d2_mcE,"d2_mcE/F");
    m_outputTree->Branch("d2_mcPx",&m_d2_mcPx,"d2_mcPx/F");
    m_outputTree->Branch("d2_mcPy",&m_d2_mcPy,"d2_mcPy/F");
    m_outputTree->Branch("d2_mcPx",&m_d2_mcPz,"d2_mcPz/F");
    //m_outputTree->Branch("d2_mcE",&m_d2_mcE,"d2_mcE/F");
    m_outputTree->Branch("d2_mcMass",&m_d2_mcMass,"d2_mcMass/F");
    //m_outputTree->Branch("d2_mcPhi",&m_d2_mcPhi,"d2_mcPhi/F");
    //m_outputTree->Branch("d2_mcTheta",&m_d2_mcTheta,"d2_mcTheta/F");
    m_outputTree->Branch("d2_mcCosTheta",&m_d2_mcCosTheta,"d2_mcCosTheta/F");

    m_outputTree->Branch("true_Energy","std::vector< float >",m_trueEnergy);
    m_outputTree->Branch("true_Px","std::vector< float >",m_true_Px);
    m_outputTree->Branch("true_Py","std::vector< float >",m_true_Py);
    m_outputTree->Branch("true_Pz","std::vector< float >",m_true_Pz);
    m_outputTree->Branch("true_index","std::vector< int >",m_true_index);
    m_outputTree->Branch("true_x","std::vector< float >",m_true_x);
    m_outputTree->Branch("true_y","std::vector< float >",m_true_y);
    m_outputTree->Branch("true_z","std::vector< float >",m_true_z);
    m_outputTree->Branch("true_CosTheta","std::vector< float >",m_true_CosTheta);
    //m_outputTree->Branch("true_Theta","std::vector< float >",m_true_Theta);
    m_outputTree->Branch("true_Phi","std::vector< float >",m_true_Phi);
    m_outputTree->Branch("true_PDGID","std::vector< int >",m_true_PDGID); 
    m_outputTree->Branch("true_GenStatus","std::vector< int >",m_true_GenStatus); 
    m_outputTree->Branch("true_numDaughters", "std::vector< int >" ,m_true_numDaughters);
    m_outputTree->Branch("true_decayTrackerCalo", "std::vector< int >" , m_true_decayTrackerCalo);
    m_outputTree->Branch("true_motherDecayTrackerCalo", "std::vector< int >" ,  m_true_motherDecayTrackerCalo);
    //m_outputTree->Branch("true_numMothers","std::vector< int >",m_true_numMothers);
    m_outputTree->Branch("true_m1_PDGID","std::vector< int >",m_true_m1_PDGID);
    //m_outputTree->Branch("true_m2_PDGID","std::vector< int >",m_true_m2_PDGID);
    m_outputTree->Branch("true_m1_GenStatus","std::vector< int >",m_true_m1_status);
    //m_outputTree->Branch("true_m2_GenStatus","std::vector< int >",m_true_m2_status);
    m_outputTree->Branch("true_m1_Energy","std::vector< float >",m_true_m1_E);
    //m_outputTree->Branch("true_m2_Energy","std::vector< float >",m_true_m2_E);
 
    m_outputTree->Branch("reco_Energy","std::vector< float >",m_recoEnergy);
    m_outputTree->Branch("reco_Px","std::vector< float >",m_reco_Px);
    m_outputTree->Branch("reco_Py","std::vector< float >",m_reco_Py);
    m_outputTree->Branch("reco_Pz","std::vector< float >",m_reco_Pz);
    m_outputTree->Branch("reco_CosTheta","std::vector< float >",m_reco_CosTheta);
    //m_outputTree->Branch("reco_Theta","std::vector< float >",m_reco_Theta);
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
    m_outputTree->Branch("reco_cluster0_energy","std::vector< float >",m_reco_cluster0_energy);
    m_outputTree->Branch("reco_cluster0_iPhi","std::vector< float >",m_reco_cluster0_iPhi);
    m_outputTree->Branch("reco_cluster0_iTheta","std::vector< float >",m_reco_cluster0_iTheta);
    m_outputTree->Branch("reco_cluster0_energyError","std::vector< float >",m_reco_cluster0_energyError);
    m_outputTree->Branch("reco_EClustersOverPTrack","std::vector< float >",m_reco_EClustersOverPTrack);
    m_outputTree->Branch("reco_EOverP_PFA","std::vector< float >",m_reco_EOverP_PFA);

    //m_outputTree->Branch("reco_CosTheta_logERW","std::vector< float >",m_reco_CosTheta_logERW);
    //m_outputTree->Branch("reco_Theta_logERW","std::vector< float >",m_reco_Theta_logERW);
    //m_outputTree->Branch("reco_Phi_logERW","std::vector< float >",m_reco_Phi_logERW);
    //m_outputTree->Branch("reco_logE_tot","std::vector< float >",m_reco_logE_tot);
    //m_outputTree->Branch("reco_E_totRW","std::vector< float >",m_reco_E_totRW);
    //m_outputTree->Branch("reco_x_logERW","std::vector< float >", m_reco_x_logERW);
    //m_outputTree->Branch("reco_y_logERW","std::vector< float >", m_reco_y_logERW);
    //m_outputTree->Branch("reco_z_logERW","std::vector< float >", m_reco_z_logERW);
    m_outputTree->Branch("reco_E_EB","std::vector< float >",m_reco_E_EB);
    m_outputTree->Branch("reco_E_EE","std::vector< float >",m_reco_E_EE);
    m_outputTree->Branch("reco_E_EO","std::vector< float >",m_reco_E_EO);
    m_outputTree->Branch("reco_E_HB","std::vector< float >",m_reco_E_HB);
    m_outputTree->Branch("reco_E_HE","std::vector< float >",m_reco_E_HE);
    m_outputTree->Branch("reco_E_HO","std::vector< float >",m_reco_E_HO);
    m_outputTree->Branch("reco_E_MB","std::vector< float >",m_reco_E_MB);
    m_outputTree->Branch("reco_E_ME","std::vector< float >",m_reco_E_ME);
    //m_outputTree->Branch("reco_E_MO","std::vector< float >",m_reco_E_MO);
    m_outputTree->Branch("reco_firstLayerECAL","std::vector< int >",m_reco_firstLayerECAL);
    m_outputTree->Branch("reco_lastLayerECAL","std::vector< int >",m_reco_lastLayerECAL);
    m_outputTree->Branch("reco_nhitsEB","std::vector< int >",m_reco_nhitsEB);
    m_outputTree->Branch("reco_nhitsEE","std::vector< int >",m_reco_nhitsEE);
    m_outputTree->Branch("reco_nhitsEO","std::vector< int >",m_reco_nhitsEO);
    m_outputTree->Branch("reco_firstLayerHCAL","std::vector< int >",m_reco_firstLayerHCAL);
    m_outputTree->Branch("reco_lastLayerHCAL","std::vector< int >",m_reco_lastLayerHCAL);
    m_outputTree->Branch("reco_nhitsHB","std::vector< int >",m_reco_nhitsHB);
    m_outputTree->Branch("reco_nhitsHE","std::vector< int >",m_reco_nhitsHE);
    m_outputTree->Branch("reco_nhitsHO","std::vector< int >",m_reco_nhitsHO);
    m_outputTree->Branch("reco_nhitsMB","std::vector< int >",m_reco_nhitsMB);
    m_outputTree->Branch("reco_nhitsME","std::vector< int >",m_reco_nhitsME);


    if(m_fillTrackBranches){
      m_outputTree->Branch("track_d0", "std::vector< float >",m_track_d0);
      m_outputTree->Branch("track_z0", "std::vector< float >",m_track_z0);
      m_outputTree->Branch("track_phi0", "std::vector< float >",m_track_phi0);
      m_outputTree->Branch("track_ndf", "std::vector< int >",m_track_ndf);
      m_outputTree->Branch("track_nHits", "std::vector< int >",m_track_nHits);
      m_outputTree->Branch("track_nExpectedTrackerHits", "std::vector< float >",m_track_nExpectedTrackerHits);
      m_outputTree->Branch("track_nTrackerBarrelHits", "std::vector< int >",m_track_nTrackerBarrelHits);
      m_outputTree->Branch("track_nTrackerEndcapHits", "std::vector< int >",m_track_nTrackerEndcapHits);
      m_outputTree->Branch("track_nVertexBarrelHits", "std::vector< int >",m_track_nVertexBarrelHits);
      m_outputTree->Branch("track_nVertexEndcapHits", "std::vector< int >",m_track_nVertexEndcapHits);
      m_outputTree->Branch("track_chi2", "std::vector< float >",m_track_chi2);
      m_outputTree->Branch("track_sigmaPOverP", "std::vector< float >",m_track_sigmaPOverP);
      m_outputTree->Branch("track_x_atIP", "std::vector< float >",m_track_x_atIP);
      m_outputTree->Branch("track_y_atIP", "std::vector< float >",m_track_y_atIP);
      m_outputTree->Branch("track_z_atIP", "std::vector< float >",m_track_z_atIP);
      m_outputTree->Branch("track_Phi_atIP", "std::vector< float >",m_track_Phi_atIP);
      m_outputTree->Branch("track_Theta_atIP", "std::vector< float >",m_track_Theta_atIP);
      m_outputTree->Branch("track_pt_atIP", "std::vector< float >",m_track_pt_atIP);
      m_outputTree->Branch("track_p_atIP", "std::vector< float >",m_track_p_atIP);
      m_outputTree->Branch("track_x_atCalo", "std::vector< float >",m_track_x_atCalo);
      m_outputTree->Branch("track_y_atCalo", "std::vector< float >",m_track_y_atCalo);
      m_outputTree->Branch("track_z_atCalo", "std::vector< float >",m_track_z_atCalo);
      m_outputTree->Branch("track_zMin", "std::vector< float >", m_track_zMin);
      m_outputTree->Branch("track_x_innermostHit", "std::vector< float >",m_track_x_innermostHit);
      m_outputTree->Branch("track_y_innermostHit", "std::vector< float >",m_track_y_innermostHit);
      m_outputTree->Branch("track_z_innermostHit", "std::vector< float >",m_track_z_innermostHit);
      m_outputTree->Branch("track_pt_innermostHit", "std::vector< float >",m_track_pt_innermostHit);
      m_outputTree->Branch("track_p_innermostHit", "std::vector< float >",m_track_p_innermostHit);
      m_outputTree->Branch("track_x_outermostRHit", "std::vector< float >",m_track_x_outermostRHit);
      m_outputTree->Branch("track_y_outermostRHit", "std::vector< float >",m_track_y_outermostRHit);
      m_outputTree->Branch("track_z_outermostRHit", "std::vector< float >",m_track_z_outermostRHit);
      m_outputTree->Branch("track_pt_outermostRHit", "std::vector< float >",m_track_pt_outermostRHit);
      m_outputTree->Branch("track_p_outermostRHit", "std::vector< float >",m_track_p_outermostRHit);
      m_outputTree->Branch("track_x_outermostZHit", "std::vector< float >",m_track_x_outermostZHit);
      m_outputTree->Branch("track_y_outermostZHit", "std::vector< float >",m_track_y_outermostZHit);
      m_outputTree->Branch("track_z_outermostZHit", "std::vector< float >",m_track_z_outermostZHit);
      m_outputTree->Branch("track_px_atCalo", "std::vector< float >",m_track_px_atCalo);
      m_outputTree->Branch("track_py_atCalo", "std::vector< float >",m_track_py_atCalo);
      m_outputTree->Branch("track_pz_atCalo", "std::vector< float >",m_track_pz_atCalo);
      m_outputTree->Branch("track_Phi_atCalo", "std::vector< float >",m_track_Phi_atCalo);
      m_outputTree->Branch("track_Theta_atCalo", "std::vector< float >",m_track_Theta_atCalo);
      m_outputTree->Branch("track_pt_atCalo", "std::vector< float >",m_track_pt_atCalo);
      m_outputTree->Branch("track_p_atCalo", "std::vector< float >",m_track_p_atCalo);
      m_outputTree->Branch("track_cluster_PosAngle_min_atCalo", "std::vector< float >",m_track_cluster_PosAngle_min_atCalo);
      m_outputTree->Branch("track_minDist_cluster_EOverP","std::vector< float >", m_track_minDist_cluster_EOverP);
      m_outputTree->Branch("track_cluster_minPFAdistance_atCalo","std::vector< float >",m_track_cluster_minPFAdistance_atCalo);
      m_outputTree->Branch("track_cluster_minPFADist_EOverP","std::vector< float >",m_track_cluster_minPFADist_EOverP);
      m_outputTree->Branch("track_cluster_minPFADist_clusterEMEnergy","std::vector< float >",m_track_cluster_minPFADist_clusterEMEnergy);
      m_outputTree->Branch("track_cluster_minPFADist_clusterHadEnergy","std::vector< float >",m_track_cluster_minPFADist_clusterHadEnergy);
      m_outputTree->Branch("track_cluster_min_parPFAdistance_atCalo","std::vector< float >",m_track_cluster_min_parPFAdistance_atCalo);
      m_outputTree->Branch("track_cluster_DirAngle_minPFAdistance_atCalo","std::vector< float >",m_track_cluster_DirAngle_minPFAdistance_atCalo);
      m_outputTree->Branch("track_cluster_minPFAdistance_atCalo_tinyHad","std::vector< float >",m_track_cluster_minPFAdistance_atCalo_tinyHad);
      m_outputTree->Branch("track_cluster_minPFADist_EOverP_tinyHad","std::vector< float >",m_track_cluster_minPFADist_EOverP_tinyHad);
      m_outputTree->Branch("track_cluster_minPFADist_clusterEMEnergy_tinyHad","std::vector< float >",m_track_cluster_minPFADist_clusterEMEnergy_tinyHad);
      m_outputTree->Branch("track_cluster_minPFADist_clusterHadEnergy_tinyHad","std::vector< float >",m_track_cluster_minPFADist_clusterHadEnergy_tinyHad);
      m_outputTree->Branch("track_cluster_min_parPFAdistance_atCalo_tinyHad","std::vector< float >",m_track_cluster_min_parPFAdistance_atCalo_tinyHad);
      m_outputTree->Branch("track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad","std::vector< float >",m_track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad);
    }

    if(m_runPhotonMode){
      m_outputTree->Branch("trueConvE1Px", "std::vector< float >",m_true_conv_e1_Px);
      m_outputTree->Branch("trueConvE1Py", "std::vector< float >",m_true_conv_e1_Py);
      m_outputTree->Branch("trueConvE1Pz", "std::vector< float >",m_true_conv_e1_Pz);
      m_outputTree->Branch("trueConvE1E", "std::vector< float >",m_true_conv_e1_E);
      m_outputTree->Branch("trueConvE1PDGID", "std::vector< int >",m_true_conv_e1_PDGID);
      m_outputTree->Branch("trueConvE2Px", "std::vector< float >",m_true_conv_e2_Px);
      m_outputTree->Branch("trueConvE2Py", "std::vector< float >",m_true_conv_e2_Py);
      m_outputTree->Branch("trueConvE2Pz", "std::vector< float >",m_true_conv_e2_Pz);
      m_outputTree->Branch("trueConvE2E", "std::vector< float >",m_true_conv_e2_E);
      m_outputTree->Branch("trueConvE2PDGID", "std::vector< int >",m_true_conv_e2_PDGID);
      m_outputTree->Branch("trueConvPhVecInd", "std::vector< int >",m_true_conv_Ph_VecInd);
      m_outputTree->Branch("trueConvVtxx", "std::vector< float >",m_true_conv_Vtx_x);
      m_outputTree->Branch("trueConvVtxy", "std::vector< float >",m_true_conv_Vtx_y);
      m_outputTree->Branch("trueConvVtxz", "std::vector< float >",m_true_conv_Vtx_z);
    }


    //m_outputTree->Branch("reco_nhitsMO","std::vector< int >",m_reco_nhitsMO); 
    if(m_runTauMode){
      m_outputTree->Branch("trueTauDaughterE", "std::vector< float >",m_true_tauDaughter_Energy);
      m_outputTree->Branch("trueTauDaughterPx", "std::vector< float >",m_true_tauDaughter_Px);
      m_outputTree->Branch("trueTauDaughterPy", "std::vector< float >",m_true_tauDaughter_Py);
      m_outputTree->Branch("trueTauDaughterPz", "std::vector< float >",m_true_tauDaughter_Pz);
      m_outputTree->Branch("trueTauDaughterPDGID", "std::vector< int >",m_true_tauDaughter_PDGID);
      m_outputTree->Branch("trueTauDaughterCharge", "std::vector< int >",m_true_tauDaughter_Charge);
      m_outputTree->Branch("trueTauDaughterTauIndex", "std::vector< int >",m_true_tauDaughter_tauIndex);
      m_outputTree->Branch("trueTauDaughterStatus", "std::vector< int >",m_true_tauDaughter_status);
      m_outputTree->Branch("trueTauDaughterMotherPDGID", "std::vector< int >",m_true_tauDaughter_motherPDGID); //direct mother, e.g. pi0
      m_outputTree->Branch("trueTauDaughterMotherEnergy", "std::vector< float >",m_true_tauDaughter_motherEnergy);

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
    }else{
      //output to determine things for particle gun samples
      m_outputTree->Branch("true_stable_firstDaughter_PDGID", &m_true_stable_firstDaughter_PDGID,"true_stable_firstDaughter_PDGID/I");
      m_outputTree->Branch("true_stable_firstDaughter_E",&m_true_stable_firstDaughter_E,"true_stable_firstDaughter_E/F");
      m_outputTree->Branch("true_stable_firstDaughterVtx_x",&m_true_stable_firstDaughterVtx_x,"true_stable_firstDaughterVtx_x/F");
      m_outputTree->Branch("true_stable_firstDaughterVtx_y",&m_true_stable_firstDaughterVtx_y,"true_stable_firstDaughterVtx_y/F");
      m_outputTree->Branch("true_stable_firstDaughterVtx_z",&m_true_stable_firstDaughterVtx_z,"true_stable_firstDaughterVtx_z/F");
      
      m_outputTree->Branch("true_stable_10pc_E",&m_true_stable_10pc_E,"true_stable_10pc_E/F");
      m_outputTree->Branch("true_stable_10pc_Esum",&m_true_stable_10pc_Esum,"true_stable_10pc_Esum/F");
      m_outputTree->Branch("true_stable_10pcVtx_x",&m_true_stable_10pcVtx_x,"true_stable_10pcVtx_x/F");
      m_outputTree->Branch("true_stable_10pcVtx_y",&m_true_stable_10pcVtx_y,"true_stable_10pcVtx_y/F");
      m_outputTree->Branch("true_stable_10pcVtx_z",&m_true_stable_10pcVtx_z,"true_stable_10pcVtx_z/F");
      
      m_outputTree->Branch("true_stable_Endpoint_x",&m_true_stable_Endpoint_x,"true_stable_Endpoint_x/F");
      m_outputTree->Branch("true_stable_Endpoint_y",&m_true_stable_Endpoint_y,"true_stable_Endpoint_y/F");
      m_outputTree->Branch("true_stable_Endpoint_z",&m_true_stable_Endpoint_z,"true_stable_Endpoint_z/F");
      if(m_fillClusterBranches){
	m_outputTree->Branch("cluster_E_EB","std::vector< float >",m_cluster_E_EB);
	m_outputTree->Branch("cluster_E_EE","std::vector< float >",m_cluster_E_EE);
	m_outputTree->Branch("cluster_E_EP","std::vector< float >",m_cluster_E_EO);
	m_outputTree->Branch("cluster_E_HB","std::vector< float >",m_cluster_E_HB);
	m_outputTree->Branch("cluster_E_HE","std::vector< float >",m_cluster_E_HE);
	m_outputTree->Branch("cluster_E_HP","std::vector< float >",m_cluster_E_HO);
	m_outputTree->Branch("cluster_E_MU","std::vector< float >",m_cluster_E_MU);
	m_outputTree->Branch("cluster_x", "std::vector< float >",m_cluster_x);
	m_outputTree->Branch("cluster_y", "std::vector< float >",m_cluster_y);
	m_outputTree->Branch("cluster_z", "std::vector< float >",m_cluster_z);
	m_outputTree->Branch("cluster_dir_x", "std::vector< float >",m_cluster_dir_x);
	m_outputTree->Branch("cluster_dir_y", "std::vector< float >",m_cluster_dir_y);
	m_outputTree->Branch("cluster_dir_z", "std::vector< float >",m_cluster_dir_z);
	m_outputTree->Branch("cluster_track_PosAngle_min", "std::vector< float >",m_cluster_track_PosAngle_min);
	m_outputTree->Branch("cluster_track_DirAngle", "std::vector< float >",m_cluster_track_DirAngle);
	m_outputTree->Branch("cluster_energy", "std::vector< float >",m_cluster_energy);
	m_outputTree->Branch("cluster_energyError", "std::vector< float >",m_cluster_energyError);
	//m_outputTree->Branch("cluster_ECAL_energy", "std::vector< float >",m_cluster_ECAL_energy);
	//m_outputTree->Branch("cluster_HCAL_energy", "std::vector< float >",m_cluster_HCAL_energy);
	m_outputTree->Branch("cluster_nHits_ECAL", "std::vector< int >",m_cluster_nHits_ECAL);
	m_outputTree->Branch("cluster_nHits_HCAL", "std::vector< int >",m_cluster_nHits_HCAL);
	m_outputTree->Branch("cluster_nHits_MUON", "std::vector< int >",m_cluster_nHits_MUON);
	m_outputTree->Branch("cluster_iPhi", "std::vector< float >",m_cluster_iPhi);
	m_outputTree->Branch("cluster_iTheta", "std::vector< float >",m_cluster_iTheta);
	m_outputTree->Branch("cluster_track_minPFADist_cluster_EOverP", "std::vector< float >",m_cluster_track_minPFADist_cluster_EOverP);
	m_outputTree->Branch("cluster_track_minPFADist_cluster_Angle", "std::vector< float >",m_cluster_track_minPFADist_cluster_Angle);
	m_outputTree->Branch("cluster_track_minPFADist_TrackClusterDistance", "std::vector< float >", m_cluster_track_minPFADist_cluster_TCDistance);
	if(m_fillClusterHits){
	  m_outputTree->Branch("cluster_hit_x", "std::vector< float >",m_cluster_hit_x);
	  m_outputTree->Branch("cluster_hit_y", "std::vector< float >",m_cluster_hit_y);
	  m_outputTree->Branch("cluster_hit_z", "std::vector< float >",m_cluster_hit_z);
	  m_outputTree->Branch("cluster_hit_E", "std::vector< float >",m_cluster_hit_E);
	  m_outputTree->Branch("cluster_hit_index", "std::vector< int >",m_cluster_hit_index);
	  m_outputTree->Branch("cluster_hit_type", "std::vector< int >",m_cluster_hit_type);
	  
	  m_outputTree->Branch("no_cluster_hit_x", "std::vector< float >",m_no_cluster_hit_x);
	  m_outputTree->Branch("no_cluster_hit_y", "std::vector< float >",m_no_cluster_hit_y);
	  m_outputTree->Branch("no_cluster_hit_z", "std::vector< float >",m_no_cluster_hit_z);
	  m_outputTree->Branch("no_cluster_hit_E", "std::vector< float >",m_no_cluster_hit_E);
	  //m_outputTree->Branch("no_cluster_hit_index", "std::vector< int >",m_no_cluster_hit_index);
	  //m_outputTree->Branch("no_cluster_hit_type", "std::vector< int >",m_no_cluster_hit_type);	  
	}
      }
      if(m_fillSimHits){
	  m_outputTree->Branch("muon_b_simhit_x", "std::vector< float >",m_muon_b_simhit_x);
	  m_outputTree->Branch("muon_b_simhit_y", "std::vector< float >",m_muon_b_simhit_y);
	  m_outputTree->Branch("muon_b_simhit_z", "std::vector< float >",m_muon_b_simhit_z);
	  m_outputTree->Branch("muon_b_simhit_E", "std::vector< float >",m_muon_b_simhit_E);
      }
    }
}


void PionStudy::processRunHeader( LCRunHeader*) {
  //++m_runNumber ;
}

void PionStudy::processEvent( LCEvent* evt ) {

  m_runNumber=evt->getRunNumber();
  m_eventNumber=evt->getEventNumber();

  eventcount+=1;
  if(evt->getEventNumber()%50==0){
    std::cout<<"run/evt "<<evt->getRunNumber()<<"/"<<evt->getEventNumber()<<"/"<<eventcount<<std::endl;
  }
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

    m_trueEnergy->clear();
    m_true_Px->clear();
    m_true_Py->clear();
    m_true_Pz->clear();
    m_true_x->clear();
    m_true_y->clear();
    m_true_z->clear();
    m_true_index->clear();
    m_true_CosTheta->clear();
    m_true_Theta->clear();
    m_true_Phi->clear();
    m_true_PDGID->clear();
    m_true_GenStatus->clear();
    m_true_numDaughters->clear();
    m_true_decayTrackerCalo->clear();
    m_true_motherDecayTrackerCalo->clear();
    m_true_numMothers->clear();
    m_true_m1_PDGID->clear();
    m_true_m2_PDGID->clear();
    m_true_m1_status->clear();
    m_true_m2_status->clear();
    m_true_m1_E->clear();
    m_true_m2_E->clear();

    if(m_runPhotonMode){
      m_true_conv_e1_Px->clear();
      m_true_conv_e1_Py->clear();
      m_true_conv_e1_Pz->clear();
      m_true_conv_e1_E->clear();
      m_true_conv_e1_PDGID->clear();
      m_true_conv_e2_Px->clear();
      m_true_conv_e2_Py->clear();
      m_true_conv_e2_Pz->clear();
      m_true_conv_e2_E->clear();
      m_true_conv_e2_PDGID->clear();
      m_true_conv_Ph_VecInd->clear();
      m_true_conv_Vtx_x->clear();
      m_true_conv_Vtx_y->clear();
      m_true_conv_Vtx_z->clear();
    }

    if(m_runTauMode){
      m_true_tauDaughter_Energy->clear();
      m_true_tauDaughter_Px->clear();
      m_true_tauDaughter_Py->clear();
      m_true_tauDaughter_Pz->clear();
      m_true_tauDaughter_PDGID->clear();
      m_true_tauDaughter_Charge->clear();
      m_true_tauDaughter_tauIndex->clear();
      m_true_tauDaughter_status->clear();
      m_true_tauDaughter_motherPDGID  ->clear();
      m_true_tauDaughter_motherEnergy ->clear();

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
    }else{
      m_true_stable_firstDaughter_PDGID=0;
      m_true_stable_firstDaughter_E=0;
      m_true_stable_firstDaughterVtx_x=0;
      m_true_stable_firstDaughterVtx_y=0;
      m_true_stable_firstDaughterVtx_z=0;

      m_true_stable_10pc_E=0;//last daughter to be created, when over 10 % of energy is lost
      m_true_stable_10pc_Esum=0;//sum of energy taken away by daughters beyond the 10 % energy loss point
      m_true_stable_10pcVtx_x=0;//to determine true vertex, at which point the original particle loses more than 10 % of its energy
      m_true_stable_10pcVtx_y=0;
      m_true_stable_10pcVtx_z=0;

      m_true_stable_Endpoint_x=0;
      m_true_stable_Endpoint_y=0;
      m_true_stable_Endpoint_z=0;

      m_cluster_x->clear();
      m_cluster_y->clear();
      m_cluster_z->clear();
      m_cluster_dir_x->clear();
      m_cluster_dir_y->clear();
      m_cluster_dir_z->clear();
      m_cluster_energy->clear();
      m_cluster_energyError->clear();
      m_cluster_ECAL_energy->clear();
      m_cluster_HCAL_energy->clear();
      m_cluster_nHits_ECAL->clear();
      m_cluster_nHits_HCAL->clear();
      m_cluster_nHits_MUON->clear();
      m_cluster_E_EB->clear();
      m_cluster_E_EE->clear();
      m_cluster_E_EO->clear();
      m_cluster_E_HB->clear();
      m_cluster_E_HE->clear();
      m_cluster_E_HO->clear();
      m_cluster_E_MU->clear();
      m_cluster_iPhi->clear();
      m_cluster_iTheta->clear();
      m_cluster_track_PosAngle_min->clear();
      m_cluster_track_DirAngle->clear();
      m_cluster_track_minPFADist_cluster_EOverP->clear();
      m_cluster_track_minPFADist_cluster_Angle->clear();
      m_cluster_track_minPFADist_cluster_TCDistance->clear();

      if(m_fillSimHits){
	m_muon_b_simhit_x->clear();
	m_muon_b_simhit_y->clear();
	m_muon_b_simhit_z->clear();
	m_muon_b_simhit_E->clear();
      }

      if(m_fillClusterHits){
	m_cluster_hit_x->clear();
	m_cluster_hit_y->clear();
	m_cluster_hit_z->clear();
	m_cluster_hit_E->clear();
        m_cluster_hit_index->clear();
        m_cluster_hit_type->clear();

	m_no_cluster_hit_x->clear();
	m_no_cluster_hit_y->clear();
	m_no_cluster_hit_z->clear();
	m_no_cluster_hit_E->clear();
        //m_no_cluster_hit_index->clear();
        //m_no_cluster_hit_type->clear();
      }
    }
      
    if(m_fillTrackBranches){
      m_track_d0->clear();
      m_track_z0->clear();
      m_track_phi0->clear();
      m_track_ndf->clear();
      m_track_nHits->clear();
      m_track_nExpectedTrackerHits->clear();
      m_track_nTrackerBarrelHits->clear();
      m_track_nTrackerEndcapHits->clear();
      m_track_nVertexBarrelHits->clear();
      m_track_nVertexEndcapHits->clear();
      m_track_chi2->clear();
      m_track_sigmaPOverP->clear();
      m_track_x_atIP->clear();
      m_track_y_atIP->clear();
      m_track_z_atIP->clear();
      m_track_Phi_atIP->clear();
      m_track_Theta_atIP->clear();
      m_track_pt_atIP->clear();
      m_track_p_atIP->clear();
      m_track_x_atCalo->clear();
      m_track_y_atCalo->clear();
      m_track_z_atCalo->clear();
      m_track_zMin->clear();
      m_track_x_innermostHit->clear();
      m_track_y_innermostHit->clear();
      m_track_z_innermostHit->clear();
      m_track_pt_innermostHit->clear();
      m_track_p_innermostHit->clear();
      m_track_x_outermostRHit->clear();
      m_track_y_outermostRHit->clear();
      m_track_z_outermostRHit->clear();
      m_track_pt_outermostRHit->clear();
      m_track_p_outermostRHit->clear();
      m_track_x_outermostZHit->clear();
      m_track_y_outermostZHit->clear();
      m_track_z_outermostZHit->clear();
      m_track_px_atCalo->clear();
      m_track_px_atCalo->clear();
      m_track_px_atCalo->clear();
      m_track_Phi_atCalo->clear();
      m_track_Theta_atCalo->clear();
      m_track_pt_atCalo->clear();
      m_track_p_atCalo->clear();
      m_track_cluster_PosAngle_min_atCalo->clear();
      m_track_minDist_cluster_EOverP->clear();
      m_track_cluster_minPFADist_EOverP->clear();
      m_track_cluster_minPFADist_clusterEMEnergy->clear();
      m_track_cluster_minPFADist_clusterHadEnergy->clear();
      m_track_cluster_minPFAdistance_atCalo->clear();
      m_track_cluster_min_parPFAdistance_atCalo->clear();
      m_track_cluster_DirAngle_minPFAdistance_atCalo->clear();
      m_track_cluster_minPFADist_EOverP_tinyHad->clear();
      m_track_cluster_minPFADist_clusterEMEnergy_tinyHad->clear();
      m_track_cluster_minPFADist_clusterHadEnergy_tinyHad->clear();
      m_track_cluster_minPFAdistance_atCalo_tinyHad->clear();
      m_track_cluster_min_parPFAdistance_atCalo_tinyHad->clear();
      m_track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad->clear();
    }
    m_recoEnergy->clear();
    m_recohitEnergyCellCorrection->clear();
    m_reco_Px->clear();
    m_reco_Py->clear();
    m_reco_Pz->clear();
    m_reco_CosTheta->clear();
    m_reco_Theta->clear();
    m_reco_Phi->clear();
    m_reco_Charge->clear();
    m_reco_PDGID->clear();
    m_reco_CosTheta_logERW->clear();
    m_reco_Theta_logERW->clear();
    m_reco_Phi_logERW->clear();
    m_reco_x_logERW->clear();
    m_reco_y_logERW->clear();
    m_reco_z_logERW->clear();
    m_reco_logE_tot->clear();
    m_reco_E_totRW->clear();

    m_reco_nTracks->clear();
    m_reco_track0_pt->clear();
    m_reco_track0_p->clear();
    m_reco_track0_nHits->clear();
    m_reco_track0_chi2OverNdof->clear();
    m_reco_nClusters->clear();
    m_reco_clusters_energy->clear();
    m_reco_cluster0_energy->clear();
    m_reco_cluster0_iPhi->clear();
    m_reco_cluster0_iTheta->clear();
    m_reco_cluster0_energyError->clear();
    m_reco_EClustersOverPTrack->clear();
    m_reco_EOverP_PFA->clear();
    m_reco_E_EB->clear();
    m_reco_E_EE->clear();
    m_reco_E_EO->clear();
    m_reco_E_HB->clear();
    m_reco_E_HE->clear();
    m_reco_E_HO->clear();
    m_reco_E_MB->clear();
    m_reco_E_ME->clear();
    m_reco_E_MO->clear();
    m_reco_firstLayerECAL->clear();
    m_reco_lastLayerECAL->clear();
    m_reco_nhitsEB->clear();
    m_reco_nhitsEE->clear();
    m_reco_nhitsEO->clear();
    m_reco_firstLayerHCAL->clear();
    m_reco_lastLayerHCAL->clear();
    m_reco_nhitsHB->clear();
    m_reco_nhitsHE->clear();
    m_reco_nhitsHO->clear();
    m_reco_nhitsMB->clear();
    m_reco_nhitsME->clear();
    m_reco_nhitsMO->clear();

    bool pass_tau_gen=false;
    bool pass_tau_reco=false;

    bool true_in_transition = true;

    float trueCosTheta=-1.1;

    std::set<CalorimeterHit*>clustermuonhits;

    //std::cout<<"cut 0"<<std::endl;

    //only run tracks and clusters for non tau events ->>validation and efficiency studies
    if(m_runTauMode){
      //tau jet loop
      LCCollection * tauJetColl =0;
      getCollection(tauJetColl,m_inputTauCollection,evt);
      if(tauJetColl!=NULL){
	if(tauJetColl->getNumberOfElements()>0){
	  pass_tau_reco=true;
	  //std::cout<<"should get to tau mode"<<tauJetColl->getNumberOfElements()<<std::endl;
	}
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

    //std::cout<<"cut 1"<<std::endl;


    LCCollection * mcColl =0;
    getCollection(mcColl,m_inputMCParticleCollection,evt);
    bool skip_tau_gen_standalone=false;
    std::vector<TLorentzVector>tau_veto_vector;
    int index_Z=-1;
    if(m_runTauMode){
      skip_tau_gen_standalone=true;
      for(int m =0; m< mcColl->getNumberOfElements(); m++){
	MCParticle* mcp= dynamic_cast<MCParticle*>( mcColl->getElementAt(m) ) ;
       if(abs(mcp->getPDG())==15 && (mcp->getGeneratorStatus()==2 || mcp->getGeneratorStatus()==1)){
	 skip_tau_gen_standalone=false;
       }
       if(mcp->getPDG() == 23 && (mcp->getDaughters().size ()>1)){
	 index_Z=m;
       }
      }
    }
    //check if we can avoid running on gen as well for taus, since none on reco level, maybe none on tau gen level
    if(m_runTauMode && skip_tau_gen_standalone && index_Z>=0){
      //otherwise skip everything for skip tau settings
      MCParticle* mcp= dynamic_cast<MCParticle*>( mcColl->getElementAt(index_Z) ) ;
      m_Z_mcE=mcp->getEnergy();
      m_Z_mcNDaughter=mcp->getDaughters().size ();
      for(unsigned int i=0;i<mcp->getDaughters().size ();i++){
	if(abs(mcp->getDaughters()[i]->getPDG())<7 && m_d1_mcE<0) {
	  m_d1_mcPDGID=mcp->getDaughters()[i]->getPDG();
	  m_d1_mcE=mcp->getDaughters()[i]->getEnergy();
	  m_d1_mcPx=mcp->getDaughters()[i]->getMomentum()[0];
	  m_d1_mcPy=mcp->getDaughters()[i]->getMomentum()[1];
	  m_d1_mcPz=mcp->getDaughters()[i]->getMomentum()[2];
	  m_d1_mcMass=mcp->getDaughters()[i]->getMass();
	  m_d1_mcPhi=atan2(m_d1_mcPy,m_d1_mcPx);
	  m_d1_mcCosTheta= m_d1_mcPz/sqrt(m_d1_mcPx*m_d1_mcPx+m_d1_mcPy*m_d1_mcPy+m_d1_mcPz*m_d1_mcPz);
	  m_d1_mcTheta=acos(m_d1_mcCosTheta);
	}
	if(m_d2_mcE<0 && (abs(mcp->getDaughters()[i]->getPDG())<7 && mcp->getDaughters()[i]->getPDG()==(-m_d1_mcPDGID))){
	  m_d2_mcPDGID=mcp->getDaughters()[i]->getPDG();
	  m_d2_mcE=mcp->getDaughters()[i]->getEnergy();
	  m_d2_mcPx=mcp->getDaughters()[i]->getMomentum()[0];
	  m_d2_mcPy=mcp->getDaughters()[i]->getMomentum()[1];
	  m_d2_mcPz=mcp->getDaughters()[i]->getMomentum()[2];
	  m_d2_mcMass=mcp->getDaughters()[i]->getMass();
	  m_d2_mcPhi=atan2(m_d2_mcPy,m_d2_mcPx);
	  m_d2_mcCosTheta= m_d2_mcPz/sqrt(m_d2_mcPx*m_d2_mcPx+m_d2_mcPy*m_d2_mcPy+m_d2_mcPz*m_d2_mcPz);
	  m_d2_mcTheta=acos(m_d2_mcCosTheta);
	}
      }
    }
    for(int m =0; m< mcColl->getNumberOfElements(); m++){
      MCParticle* mcp= dynamic_cast<MCParticle*>( mcColl->getElementAt(m) ) ;
      if(skip_tau_gen_standalone){
	break;
      }
      if(mcp->getPDG() == 23 && (mcp->getDaughters().size ()>1)){
	m_Z_mcE=mcp->getEnergy();
	m_Z_mcNDaughter=mcp->getDaughters().size ();
	for(unsigned int i=0;i<mcp->getDaughters().size ();i++){
	  if(abs(mcp->getDaughters()[i]->getPDG())<7 && m_d1_mcE<0) {
	    m_d1_mcPDGID=mcp->getDaughters()[i]->getPDG();
	    m_d1_mcE=mcp->getDaughters()[i]->getEnergy();
	    m_d1_mcPx=mcp->getDaughters()[i]->getMomentum()[0];
	    m_d1_mcPy=mcp->getDaughters()[i]->getMomentum()[1];
	    m_d1_mcPz=mcp->getDaughters()[i]->getMomentum()[2];
	    m_d1_mcMass=mcp->getDaughters()[i]->getMass();
	    m_d1_mcPhi=atan2(m_d1_mcPy,m_d1_mcPx);
	    m_d1_mcCosTheta= m_d1_mcPz/sqrt(m_d1_mcPx*m_d1_mcPx+m_d1_mcPy*m_d1_mcPy+m_d1_mcPz*m_d1_mcPz);
	    m_d1_mcTheta=acos(m_d1_mcCosTheta);
	  }
	  if(m_d2_mcE<0 && (abs(mcp->getDaughters()[i]->getPDG())<7 && mcp->getDaughters()[i]->getPDG()==(-m_d1_mcPDGID))){
	    m_d2_mcPDGID=mcp->getDaughters()[i]->getPDG();
	    m_d2_mcE=mcp->getDaughters()[i]->getEnergy();
	    m_d2_mcPx=mcp->getDaughters()[i]->getMomentum()[0];
	    m_d2_mcPy=mcp->getDaughters()[i]->getMomentum()[1];
	    m_d2_mcPz=mcp->getDaughters()[i]->getMomentum()[2];
	    m_d2_mcMass=mcp->getDaughters()[i]->getMass();
	    m_d2_mcPhi=atan2(m_d2_mcPy,m_d2_mcPx);
	    m_d2_mcCosTheta= m_d2_mcPz/sqrt(m_d2_mcPx*m_d2_mcPx+m_d2_mcPy*m_d2_mcPy+m_d2_mcPz*m_d2_mcPz);
	    m_d2_mcTheta=acos(m_d2_mcCosTheta);
	  }
	}	  	
      }
      if(m_runTauMode && abs(mcp->getPDG())!=15){
	continue;
      }
      //save only tau's which have a stable particle (status==1) as daughter, else go on
      bool tau_has_stable_daughter=false;
      bool tau_veto_inhistory=false;//if tau is a result of another tau decay
      //take it out of history maybe to avoid double counting?
      //usually the case of a brem photon radiated of the tau
      if(!m_runTauMode){
	tau_has_stable_daughter=true;
      }else{
	if(abs(mcp->getPDG()) == 15 && (mcp->getDaughters().size ()>1)){
	  TLorentzVector tau_test(0,0,0,0);
	  tau_test.SetPxPyPzE(mcp->getMomentum()[0],mcp->getMomentum()[1],mcp->getMomentum()[2],mcp->getEnergy());
	  for(unsigned int tv=0;tv<tau_veto_vector.size();tv++){
	    if(tau_test.DeltaR(tau_veto_vector[tv])==0){
	      tau_veto_inhistory=true;
	      //std::cout<<"tau already in history it seems"<<std::endl;
	    }
	  }
	  for(unsigned int i=0;i<mcp->getDaughters().size ();i++){
	    if(mcp->getDaughters()[i]->getGeneratorStatus()==1){
	      tau_has_stable_daughter=true;
	    }
	  }
	}
      }
      if(!tau_has_stable_daughter || tau_veto_inhistory){
	continue;
      }
      pass_tau_gen=true;
      if(mcp->getGeneratorStatus()==0 && mcp->isDecayedInTracker()){
	continue;//don't save particle which are created in simulation, but decay inside the tracker volume
      }
      if(mcp->getParents().size ()>0){
	if(mcp->getParents()[0]->isDecayedInTracker()){
	  m_true_motherDecayTrackerCalo->push_back(1);
	}else if(mcp->getParents()[0]->isDecayedInCalorimeter()){
	  m_true_motherDecayTrackerCalo->push_back(2);
	}else if(mcp->getParents()[0]->isBackscatter()){
	  m_true_motherDecayTrackerCalo->push_back(3);
	}else{
	  m_true_motherDecayTrackerCalo->push_back(-1);
	}
      }else{
	m_true_motherDecayTrackerCalo->push_back(0);
      }
      if(mcp->isBackscatter()){
	m_true_decayTrackerCalo->push_back(3);
      }else if (mcp->isDecayedInTracker()){
	m_true_decayTrackerCalo->push_back(1);
      }else if (mcp->isDecayedInCalorimeter()){
	m_true_decayTrackerCalo->push_back(2);
      }else if(mcp->hasLeftDetector()){
	m_true_decayTrackerCalo->push_back(4);
      }else{
	m_true_decayTrackerCalo->push_back(0);
      } 
      double cosTheta = mcp->getMomentum()[2]/sqrt(mcp->getMomentum()[0]*mcp->getMomentum()[0]+mcp->getMomentum()[1]*mcp->getMomentum()[1]+mcp->getMomentum()[2]*mcp->getMomentum()[2]);

      if(mcp->getGeneratorStatus()==1 && abs(mcp->getPDG()) == 13){
	trueCosTheta=cosTheta;
      }

      //if( (cosTheta>0.845) || (cosTheta< -0.845) || (cosTheta> -0.79 && cosTheta<0.79)){
      if( (cosTheta<-0.075) || (cosTheta>0.075)){
	true_in_transition = false;
	//continue;
      }
      //for particle gun samples, aka not tau checks, check history of brems and when particles send of daugthers inside the detector etc
      if(!m_runTauMode && mcp->getGeneratorStatus()==1){

	//std::cout<<"before MC brem collection"<<std::endl;

	m_true_stable_Endpoint_x=mcp->getEndpoint()[0]; 
	m_true_stable_Endpoint_y=mcp->getEndpoint()[1]; 
	m_true_stable_Endpoint_z=mcp->getEndpoint()[2]; 
	if(mcp->getDaughters().size ()>0){
	  m_true_stable_firstDaughter_PDGID=mcp->getDaughters()[0]->getPDG();
	  m_true_stable_firstDaughter_E=mcp->getDaughters()[0]->getEnergy();
	  m_true_stable_firstDaughterVtx_x=mcp->getDaughters()[0]->getVertex()[0];
	  m_true_stable_firstDaughterVtx_y=mcp->getDaughters()[0]->getVertex()[1];
	  m_true_stable_firstDaughterVtx_z=mcp->getDaughters()[0]->getVertex()[2];
	  float energy_sum=0;
	  for(unsigned int i=0;i<mcp->getDaughters().size ();i++){
	    energy_sum+=mcp->getDaughters()[i]->getEnergy();
	    if(energy_sum>(0.10*mcp->getEnergy())){
	      m_true_stable_10pc_E=mcp->getDaughters()[i]->getEnergy();
	      m_true_stable_10pc_Esum=energy_sum;
	      m_true_stable_10pcVtx_x=mcp->getDaughters()[i]->getVertex()[0];
	      m_true_stable_10pcVtx_y=mcp->getDaughters()[i]->getVertex()[1];
	      m_true_stable_10pcVtx_z=mcp->getDaughters()[i]->getVertex()[2];
	      break;
	    }
	  }
	}
      }
      m_trueEnergy->push_back(mcp->getEnergy());
      m_true_Px->push_back(mcp->getMomentum()[0]);
      m_true_Py->push_back(mcp->getMomentum()[1]);
      m_true_Pz->push_back(mcp->getMomentum()[2]);
      m_true_x->push_back(mcp->getVertex()[0]);
      m_true_y->push_back(mcp->getVertex()[1]);
      m_true_z->push_back(mcp->getVertex()[2]);
      m_true_CosTheta->push_back(cosTheta);    
      m_true_Phi->push_back(atan2(mcp->getMomentum()[1],mcp->getMomentum()[0]));
      m_true_Theta->push_back(acos(cosTheta));
      m_true_PDGID->push_back(mcp->getPDG());
      m_true_index->push_back(m);
      m_true_GenStatus->push_back(mcp->getGeneratorStatus());
      m_true_numDaughters->push_back(mcp->getDaughters().size ());
      m_true_numMothers->push_back(mcp->getParents().size());
      if(m_runPhotonMode){
	if(mcp->getPDG()==22 && mcp->getDaughters().size () ==2 && ( (mcp->getDaughters()[0]->getPDG()*mcp->getDaughters()[1]->getPDG())==-121)){
	  unsigned int ind_e1=0;
	  unsigned int ind_e2=1;
	  if(mcp->getDaughters()[1]->getEnergy()>mcp->getDaughters()[0]->getEnergy()){
	    ind_e1=1;
	    ind_e2=0;
	  }
	  m_true_conv_e1_Px->push_back(mcp->getDaughters()[ind_e1]->getMomentum()[0]);
	  m_true_conv_e1_Py->push_back(mcp->getDaughters()[ind_e1]->getMomentum()[1]);
	  m_true_conv_e1_Pz->push_back(mcp->getDaughters()[ind_e1]->getMomentum()[2]);
	  m_true_conv_e1_E->push_back(mcp->getDaughters()[ind_e1]->getEnergy());
	  m_true_conv_e1_PDGID->push_back(mcp->getDaughters()[ind_e1]->getPDG());
	  m_true_conv_e2_Px->push_back(mcp->getDaughters()[ind_e2]->getMomentum()[0]);
	  m_true_conv_e2_Py->push_back(mcp->getDaughters()[ind_e2]->getMomentum()[1]);
	  m_true_conv_e2_Pz->push_back(mcp->getDaughters()[ind_e2]->getMomentum()[2]);
	  m_true_conv_e2_E->push_back(mcp->getDaughters()[ind_e2]->getEnergy());
	  m_true_conv_e2_PDGID->push_back(mcp->getDaughters()[ind_e2]->getPDG());
	  m_true_conv_Ph_VecInd->push_back(m_trueEnergy->size()-1);
	  m_true_conv_Vtx_x->push_back(mcp->getDaughters()[ind_e1]->getVertex()[0]);
	  m_true_conv_Vtx_y->push_back(mcp->getDaughters()[ind_e1]->getVertex()[1]);
	  m_true_conv_Vtx_z->push_back(mcp->getDaughters()[ind_e1]->getVertex()[2]);
	}
	if(mcp->getPDG()==22 && mcp->getDaughters().size () ==2 && ( (mcp->getDaughters()[0]->getPDG()*mcp->getDaughters()[1]->getPDG())==121)){
	  //std::cout<<"seems we found a e e thing "<<mcp->isDecayedInTracker()<<"/"<<mcp->isDecayedInCalorimeter()<<std::endl;
	}
	
      }


      if(m_runTauMode){
	std::set<MCParticle*> tau_daughtersFunc;
	fillStableDaughterSet(mcp, tau_daughtersFunc);
	TLorentzVector temp(0,0,0,0);
	//avoid looking at tau daughters which are tau's agin --> happens in first step of photos, where a photon is irradiated 
	//and the daughter tau continues on its way
	for(unsigned int d=0;d<mcp->getDaughters().size();d++){
	  if(abs(mcp->getDaughters()[d]->getPDG())==15){
	    temp.SetPxPyPzE(mcp->getDaughters()[d]->getMomentum()[0],mcp->getDaughters()[d]->getMomentum()[1],mcp->getDaughters()[d]->getMomentum()[2],mcp->getDaughters()[d]->getEnergy());
	    tau_veto_vector.push_back(temp);
	  }
	  if(mcp->getDaughters()[d]->getGeneratorStatus()==2){
	    for(unsigned int gd=0;gd<mcp->getDaughters()[d]->getDaughters().size();gd++){
	      if(abs(mcp->getDaughters()[d]->getDaughters()[gd]->getPDG())==15){
		temp.SetPxPyPzE(mcp->getDaughters()[d]->getDaughters()[gd]->getMomentum()[0],mcp->getDaughters()[d]->getDaughters()[gd]->getMomentum()[1],mcp->getDaughters()[d]->getDaughters()[gd]->getMomentum()[2],mcp->getDaughters()[d]->getDaughters()[gd]->getEnergy());
		tau_veto_vector.push_back(temp);
	      }
	    }
	  }
	}
	std::set<MCParticle*>::iterator tauDaughtIt;
	for(tauDaughtIt=tau_daughtersFunc.begin();tauDaughtIt!=tau_daughtersFunc.end();tauDaughtIt++){
	  m_true_tauDaughter_Energy->push_back((*tauDaughtIt)->getEnergy());
	  m_true_tauDaughter_Px->push_back((*tauDaughtIt)->getMomentum()[0]);
	  m_true_tauDaughter_Py->push_back((*tauDaughtIt)->getMomentum()[1]);
	  m_true_tauDaughter_Pz->push_back((*tauDaughtIt)->getMomentum()[2]);
	  m_true_tauDaughter_PDGID->push_back((*tauDaughtIt)->getPDG());
	  m_true_tauDaughter_Charge->push_back((*tauDaughtIt)->getCharge());
	  m_true_tauDaughter_tauIndex->push_back(m);
	  m_true_tauDaughter_status->push_back((*tauDaughtIt)->getGeneratorStatus());
	  m_true_tauDaughter_motherPDGID->push_back((*tauDaughtIt)->getParents()[0]->getPDG());
	  m_true_tauDaughter_motherEnergy->push_back((*tauDaughtIt)->getParents()[0]->getEnergy());
	}
	if(!tau_daughtersFunc.empty()){
	  tau_daughtersFunc.clear();
	}
      }
      if(mcp->getParents().size()>0){
	m_true_m1_PDGID->push_back(mcp->getParents()[0]->getPDG());
	m_true_m1_status->push_back(mcp->getParents()[0]->getGeneratorStatus());
	m_true_m1_E->push_back(mcp->getParents()[0]->getEnergy());
	if(mcp->getParents().size()>1){
	  m_true_m2_PDGID->push_back(mcp->getParents()[1]->getPDG());
	  m_true_m2_status->push_back(mcp->getParents()[1]->getGeneratorStatus());
	  m_true_m2_E->push_back(mcp->getParents()[1]->getEnergy());
	}else{
	  m_true_m2_PDGID->push_back(-1);
	  m_true_m2_status->push_back(-1);
	  m_true_m2_E->push_back(-1);
	}
      }else{
	m_true_m1_PDGID->push_back(-1);
	m_true_m2_PDGID->push_back(-1);
	m_true_m1_status->push_back(-1);
	m_true_m2_status->push_back(-1);
	m_true_m1_E->push_back(-1);
	m_true_m2_E->push_back(-1);
      }
    }
    //std::cout<<"track collection"<<std::endl;

    //std::cout<<"cut 3"<<std::endl;

    LCCollection* trackcol = NULL;
    trackcol= evt->getCollection(m_inputTrackCollection);
    
    //std::cout<<"cluster collection"<<std::endl;

    LCCollection* clucol = NULL;
    clucol= evt->getCollection(m_inputClusterCollection);
    
    //std::cout<<"after cluster collection"<<std::endl;

    if(trackcol!=NULL){
      //std::cout<<"do we get here"<<std::endl;
      if(m_runPhotonMode && m_true_conv_Vtx_x->size()>0){
	//std::cout<<"we seem to have conversion /vtx r/z /cos ph/track size "<<sqrt((*m_true_conv_Vtx_x)[0]*(*m_true_conv_Vtx_x)[0]+(*m_true_conv_Vtx_y)[0]*(*m_true_conv_Vtx_y)[0])<<"/"<<(*m_true_conv_Vtx_z)[0]<<"/"<<(*m_true_CosTheta)[0]<<"/"<<trackcol->getNumberOfElements()<<std::endl;
	if(trackcol->getNumberOfElements()==2){
	  //float deltaR_min=10000;
	  //float deltaZ_min=10000;
	  float delta3D_min=10000;
	  Track* track1 = dynamic_cast<EVENT::Track*>(trackcol->getElementAt(0));
	  Track* track2 = dynamic_cast<EVENT::Track*>(trackcol->getElementAt(1));
	  const EVENT::TrackState* track1AtFirstHit=track1->getTrackState(EVENT::TrackState::AtFirstHit);
	  const EVENT::TrackState* track2AtFirstHit=track2->getTrackState(EVENT::TrackState::AtFirstHit);
	  //for(unsigned int j=0;j<track1->getTrackerHits().size();j++){
	  //TrackerHit* trackerhit = dynamic_cast<TrackerHit*>(track1->getTrackerHits()[j]);
	  //std::cout<<"all hit of track1 "<<j<<": "<<trackerhit->getPosition()[0]<<"/"<< trackerhit->getPosition()[1]<<"/"<<trackerhit->getPosition()[2] <<std::endl;
	  //}
	  /*
	  for(unsigned int j=0;j<track2->getTrackerHits().size();j++){
	    TrackerHit* trackerhit = dynamic_cast<TrackerHit*>(track2->getTrackerHits()[j]);
	    //std::cout<<"all hit of track2 "<<j<<": "<<trackerhit->getPosition()[0]<<"/"<< trackerhit->getPosition()[1]<<"/"<<trackerhit->getPosition()[2] <<" quality "<<trackerhit->getQuality()<<std::endl;
	    if(trackerhit->getQuality()==UTIL::ILDTrkHitQualityBit::USED_IN_FIT){
	      //std::cout<<"track used in fit"<<std::endl;
	    }
	    if(trackerhit->getQuality()==UTIL::ILDTrkHitQualityBit::USED_IN_TRACK){
	      //std::cout<<"track used in track"<<std::endl;
	    }
	    if(trackerhit->getQuality()==UTIL::ILDTrkHitQualityBit::DOUBLE_HIT_CANDIDATE){
	      //std::cout<<"double hit candidate"<<std::endl;
	    }
	    if(trackerhit->getQuality()==UTIL::ILDTrkHitQualityBit::GOOD){
	      //std::cout<<"good candidate"<<std::endl;
	    }
	    }*/
	  for(unsigned int j1=0;j1<track1->getTrackerHits().size();j1++){
	    for(unsigned int j2=0;j2<track2->getTrackerHits().size();j2++){
	      if(sqrt(pow(track1->getTrackerHits()[j1]->getPosition()[0]-track2->getTrackerHits()[j2]->getPosition()[0],2)+pow(track1->getTrackerHits()[j1]->getPosition()[1]-track2->getTrackerHits()[j2]->getPosition()[1],2)+pow(track1->getTrackerHits()[j1]->getPosition()[3]-track2->getTrackerHits()[j2]->getPosition()[3],2))<delta3D_min){
		delta3D_min=sqrt(pow(track1->getTrackerHits()[j1]->getPosition()[0]-track2->getTrackerHits()[j2]->getPosition()[0],2)+pow(track1->getTrackerHits()[j1]->getPosition()[1]-track2->getTrackerHits()[j2]->getPosition()[1],2)+pow(track1->getTrackerHits()[j1]->getPosition()[3]-track2->getTrackerHits()[j2]->getPosition()[3],2));
	      }
	    }
	  }
	  //std::set<TrackerHit*> track1hitset;
	  for(int i=0;i<trackcol->getNumberOfElements();i++){
	    Track* track = dynamic_cast<EVENT::Track*>(trackcol->getElementAt(i));
	    //const EVENT::TrackState* trackAtCalo=track->getTrackState(EVENT::TrackState::AtCalorimeter);
	    const EVENT::TrackState* trackAtFirstHit=track->getTrackState(EVENT::TrackState::AtFirstHit);
	    //std::cout<<"R_hitIn "<<i<<"/"<<track->getRadiusOfInnermostHit()<<" TS Hit 1 x/y/z/r "<<trackAtFirstHit->getReferencePoint()[0]<<"/"<<trackAtFirstHit->getReferencePoint()[1]<<"/"<<trackAtFirstHit->getReferencePoint()[2]<<"/"<<sqrt(trackAtFirstHit->getReferencePoint()[0]*trackAtFirstHit->getReferencePoint()[0]+trackAtFirstHit->getReferencePoint()[1]*trackAtFirstHit->getReferencePoint()[1]) <<" min dist "<<delta3D_min<<" first hit distance "<<sqrt(pow(track1AtFirstHit->getReferencePoint()[0] -track2AtFirstHit->getReferencePoint()[0] ,2)+ pow(track1AtFirstHit->getReferencePoint()[1] -track2AtFirstHit->getReferencePoint()[1] ,2)+pow(track1AtFirstHit->getReferencePoint()[2] -track2AtFirstHit->getReferencePoint()[2] ,2))<<" track hits "<<track->getTrackerHits().size()<<" d0/z0/chi2/ndf/chi2ndf/sigma/p "<<track->getD0()<<"/"<<track->getZ0()<<"/"<<track->getChi2()<<"/"<<track->getNdf()<<"/"<<track->getChi2()/(float)track->getNdf()<<"/"<<track->getCovMatrix()[5]/std::fabs(track->getOmega())<<"/"<<m_innerBField * m_const_a/fabs(trackAtFirstHit->getOmega())/cos(atan(trackAtFirstHit->getTanLambda()))<<std::endl;
	  }
	  LCCollection* recoparticlecol = NULL;
	  // Alternativelly if you do not want Marlin to exit in case of a non-existing collection
	  //use the following (commented out) code:
	  recoparticlecol = evt->getCollection(m_inputRECOParticleCollection);
	  //if(recoparticlecol!=NULL){
	  ////PandoraCandidate loop
	  ///for(int i=0;i<recoparticlecol->getNumberOfElements();i++){
	  //R/econstructedParticle* pandorapart = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(i));
	  //  std::cout<<"reco "<<i<<" type/E/px/py/pz "<<pandorapart->getType()<<"/"<<pandorapart->getEnergy()<<"/"<<pandorapart->getMomentum()[0]<<"/"<<pandorapart->getMomentum()[1]<<"/"<<pandorapart->getMomentum()[2]<<std::endl;
	  //}
	  //}
	    /*
	    if(i==0){
	      for(unsigned int j=0;j<track->getTrackerHits().size();j++){
		TrackerHit* trackerhit = dynamic_cast<TrackerHit*>(track->getTrackerHits()[j]);
		//std::cout<<"all hit of track1 "<<j<<": "<<trackerhit->getPosition()[0]<<"/"<< trackerhit->getPosition()[1]<<"/"<<trackerhit->getPosition()[2] <<std::endl;
		track1hitset.insert(trackerhit);
	      }
	    }
	    if(i==1){
	      //std::set<TrackerHit*>::iterator t1hitsIt;
	      //for(t1hitsIt=track1hitset.begin();t1hitsIt!=track1hitset.end();t1hitsIt++){
	      for(unsigned int j=0;j<track->getTrackerHits().size();j++){
		TrackerHit* trackerhit = dynamic_cast<TrackerHit*>(track->getTrackerHits()[j]);
		//std::cout<<"all hit of track2 "<<j<<": "<<trackerhit->getPosition()[0]<<"/"<< trackerhit->getPosition()[1]<<"/"<<trackerhit->getPosition()[2] <<std::endl;
		//if((trackerhit)==(*t1hitsIt)){
		//std::cout<<"hit of track2 in hitset of track 1 "<<trackerhit->getPosition()[0]<<"/"<< trackerhit->getPosition()[1]<<"/"<<trackerhit->getPosition()[2] <<std::endl;
		//}
	      }
	      //}
	      }*/
	}
      }
    }

    //std::cout<<"cut 4"<<std::endl;

    if(m_fillTrackBranches){
      //std::cout<<"what goes on here"<<std::endl;
      if(trackcol!=NULL){
	for(int i=0;i<trackcol->getNumberOfElements();i++){
	  Track* track = dynamic_cast<EVENT::Track*>(trackcol->getElementAt(i));
	  
	  int nBarrelTrackerHits = 0;
	  dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
	  const std::vector< dd4hep::DetElement>& barrelDets = dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER | dd4hep::DetType::BARREL )) ;
	  for (std::vector< dd4hep::DetElement>::const_iterator iter = barrelDets.begin(), iterEnd = barrelDets.end();iter != iterEnd; ++iter){
	    const dd4hep::DetElement & theDetector = *iter;
	    int detId = theDetector.id();
	    nBarrelTrackerHits+=track->getSubdetectorHitNumbers()[2*detId-2]; //Offset is 2 for hits in fit
	    //std::cout<<"get to barrel tracker hits "<<track->getSubdetectorHitNumbers()[2*detId-2]<<std::endl;
	  }      
	  int nEndcapTrackerHits = 0;
	  const std::vector< dd4hep::DetElement>& endcapDets = dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER | dd4hep::DetType::ENDCAP )) ;
	  for (std::vector< dd4hep::DetElement>::const_iterator iter = endcapDets.begin(), iterEnd = endcapDets.end();iter != iterEnd; ++iter){
            const dd4hep::DetElement& theDetector = *iter;
            int detId = theDetector.id();
            nEndcapTrackerHits +=track->getSubdetectorHitNumbers()[2*detId-2]; //Offset is 2 for hits in fit
	    //std::cout<<"get to endcap tracker hits "<<track->getSubdetectorHitNumbers()[2*detId-2]<<std::endl;
	  }
	  
	  m_track_nTrackerBarrelHits->push_back(nBarrelTrackerHits);
	  m_track_nTrackerEndcapHits->push_back(nEndcapTrackerHits);
	  
	  int nBarrelVtxTrackerHits = 0;
	  const std::vector< dd4hep::DetElement>& barrelVtxDets = dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER | dd4hep::DetType::VERTEX | dd4hep::DetType::BARREL )) ;
	  for (std::vector< dd4hep::DetElement>::const_iterator iter = barrelVtxDets.begin(), iterEnd = barrelVtxDets.end();iter != iterEnd; ++iter){
	    const dd4hep::DetElement & theDetector = *iter;
	    int detId = theDetector.id();
	    nBarrelVtxTrackerHits+=track->getSubdetectorHitNumbers()[2*detId-2]; //Offset is 2 for hits in fit
	    //std::cout<<"get to barrel vertex hits "<<track->getSubdetectorHitNumbers()[2*detId-2]<<std::endl;
	  }      
	  int nEndcapVtxTrackerHits = 0;
	  const std::vector< dd4hep::DetElement>& endcapVtxDets = dd4hep::DetectorSelector(mainDetector).detectors(  ( dd4hep::DetType::TRACKER | dd4hep::DetType::VERTEX | dd4hep::DetType::ENDCAP )) ;
	  for (std::vector< dd4hep::DetElement>::const_iterator iter = endcapVtxDets.begin(), iterEnd = endcapVtxDets.end();iter != iterEnd; ++iter){
	    const dd4hep::DetElement& theDetector = *iter;
	    int detId = theDetector.id();
	    nEndcapVtxTrackerHits +=track->getSubdetectorHitNumbers()[2*detId-2]; //Offset is 2 for hits in fit
	    //std::cout<<"get to endcap vertex hits "<<track->getSubdetectorHitNumbers()[2*detId-2]<<std::endl;
	  }
	  
	  m_track_nVertexBarrelHits->push_back(nBarrelVtxTrackerHits);
	  m_track_nVertexEndcapHits->push_back(nEndcapVtxTrackerHits);
	  unsigned int index_innerR=0;
	  unsigned int index_outerR=0;
	  unsigned int index_outerZ=0;
	  if(track->getTrackerHits().size()>0){
	    double R_min=sqrt(track->getTrackerHits()[0]->getPosition()[0]*track->getTrackerHits()[0]->getPosition()[0]+ track->getTrackerHits()[0]->getPosition()[1]*track->getTrackerHits()[0]->getPosition()[1]);
	    double R_max=R_min;
	    double Z_max=track->getTrackerHits()[0]->getPosition()[2];
	    double Z_min=track->getTrackerHits()[0]->getPosition()[2];
	    for(unsigned int i_thit=1;i_thit<track->getTrackerHits().size();i_thit++){
	      if(sqrt(track->getTrackerHits()[i_thit]->getPosition()[0]*track->getTrackerHits()[i_thit]->getPosition()[0]+ track->getTrackerHits()[i_thit]->getPosition()[1]*track->getTrackerHits()[i_thit]->getPosition()[1])<R_min){
		index_innerR=i_thit;
		R_min=sqrt(track->getTrackerHits()[i_thit]->getPosition()[0]*track->getTrackerHits()[i_thit]->getPosition()[0]+ track->getTrackerHits()[i_thit]->getPosition()[1]*track->getTrackerHits()[i_thit]->getPosition()[1]);
	      }  
	      if(sqrt(track->getTrackerHits()[i_thit]->getPosition()[0]*track->getTrackerHits()[i_thit]->getPosition()[0]+ track->getTrackerHits()[i_thit]->getPosition()[1]*track->getTrackerHits()[i_thit]->getPosition()[1])>R_max){
		index_outerR=i_thit;
		R_max=sqrt(track->getTrackerHits()[i_thit]->getPosition()[0]*track->getTrackerHits()[i_thit]->getPosition()[0]+ track->getTrackerHits()[i_thit]->getPosition()[1]*track->getTrackerHits()[i_thit]->getPosition()[1]);
	      } 
	      if(fabs(track->getTrackerHits()[i_thit]->getPosition()[2])<Z_min){
		Z_min=fabs(track->getTrackerHits()[i_thit]->getPosition()[2]);
	      } 
	      if(fabs(track->getTrackerHits()[i_thit]->getPosition()[2])>Z_max){
		index_outerZ=i_thit;
		Z_max=fabs(track->getTrackerHits()[i_thit]->getPosition()[2]);
	      } 
	    }
	    const EVENT::TrackState* trackAtFirstHit=track->getTrackState(EVENT::TrackState::AtFirstHit);
	    const EVENT::TrackState* trackAtLastHit=track->getTrackState(EVENT::TrackState::AtLastHit);
	    m_track_zMin->push_back(Z_min);
	    m_track_x_innermostHit->push_back(track->getTrackerHits()[index_innerR]->getPosition()[0]);
	    m_track_y_innermostHit->push_back(track->getTrackerHits()[index_innerR]->getPosition()[1]);
	    m_track_z_innermostHit->push_back(track->getTrackerHits()[index_innerR]->getPosition()[2]);
	    //m_track_r_innermostHit->push_back(R_min/*track->getRadiusOfInnermostHit()*/);
	    if(trackAtFirstHit!=0){
	      m_track_pt_innermostHit->push_back(m_innerBField * m_const_a/fabs(trackAtFirstHit->getOmega()));//omega is in terms of mm
	      m_track_p_innermostHit->push_back(m_innerBField * m_const_a/fabs(trackAtFirstHit->getOmega())/cos(atan(trackAtFirstHit->getTanLambda())));
	    }else{
	      m_track_pt_innermostHit->push_back(0);
	      m_track_p_innermostHit->push_back(0);
	    }
	    m_track_x_outermostRHit->push_back(track->getTrackerHits()[index_outerR]->getPosition()[0]);
	    m_track_y_outermostRHit->push_back(track->getTrackerHits()[index_outerR]->getPosition()[1]);
	    m_track_z_outermostRHit->push_back(track->getTrackerHits()[index_outerR]->getPosition()[2]);
	    if(trackAtLastHit!=0){ 
	      m_track_pt_outermostRHit->push_back(m_innerBField * m_const_a/fabs(trackAtLastHit->getOmega()));//omega is in terms of mm
	      m_track_p_outermostRHit->push_back(m_innerBField * m_const_a/fabs(trackAtLastHit->getOmega())/cos(atan(trackAtLastHit->getTanLambda())));
	    }else{
	      m_track_pt_outermostRHit->push_back(0);
	      m_track_p_outermostRHit->push_back(0);
	    }
	    // m_track_r_outermostRHit->push_back(sqrt(track->getTrackerHits()[index_outerR]->getPosition()[0]*track->getTrackerHits()[index_outerR]->getPosition()[0]+ track->getTrackerHits()[index_outerR]->getPosition()[1]*track->getTrackerHits()[index_outerR]->getPosition()[1]));
	    m_track_x_outermostZHit->push_back(track->getTrackerHits()[index_outerZ]->getPosition()[0]);
	    m_track_y_outermostZHit->push_back(track->getTrackerHits()[index_outerZ]->getPosition()[1]);
	    m_track_z_outermostZHit->push_back(track->getTrackerHits()[index_outerZ]->getPosition()[2]);

	    //m_track_r_outermostZHit->push_back(sqrt(track->getTrackerHits()[index_outerZ]->getPosition()[0]*track->getTrackerHits()[index_outerZ]->getPosition()[0]+ track->getTrackerHits()[index_outerZ]->getPosition()[1]*track->getTrackerHits()[index_outerZ]->getPosition()[1]));
	  }else{
	    m_track_zMin->push_back(0);
	    m_track_x_innermostHit->push_back(0);
	    m_track_y_innermostHit->push_back(0);
	    m_track_z_innermostHit->push_back(0);
	    m_track_pt_innermostHit->push_back(0);
	    m_track_p_innermostHit->push_back(0);
	    m_track_x_outermostRHit->push_back(0);
	    m_track_y_outermostRHit->push_back(0);
	    m_track_z_outermostRHit->push_back(0);
	    m_track_pt_outermostRHit->push_back(0);
	    m_track_p_outermostRHit->push_back(0);
	    m_track_x_outermostZHit->push_back(0);
	    m_track_y_outermostZHit->push_back(0);
	    m_track_z_outermostZHit->push_back(0);
	  }
	  m_track_d0->push_back(track->getD0());
	  m_track_z0->push_back(track->getZ0());
	  m_track_phi0->push_back(track->getPhi());
	  m_track_ndf->push_back(track->getNdf());
	  m_track_nHits->push_back(track->getTrackerHits().size());
	  m_track_chi2->push_back(track->getChi2());
	  m_track_sigmaPOverP->push_back(track->getCovMatrix()[5]/std::fabs(track->getOmega()));
	  const EVENT::TrackState* trackAtCalo=track->getTrackState(EVENT::TrackState::AtCalorimeter);
	  const EVENT::TrackState* trackAtIP=track->getTrackState(EVENT::TrackState::AtIP);
	  //std::cout<<"track "<<i <<" d0/z0/phi0/hits/chi2/pt/p atIP "<<track->getD0()<<"/"<<track->getZ0()<<"/"<<track->getPhi()<<"/"<<track->getTrackerHits().size()<<"/"<<track->getChi2()<<"/"<<m_innerBField * m_const_a/fabs(trackAtIP->getOmega())<<"/"<<m_innerBField * m_const_a/fabs(trackAtIP->getOmega())/cos(atan(trackAtIP->getTanLambda()))<<std::endl;
	  if(trackAtIP!=0){
	    if(fabs(track->getOmega())<0.0012){
	      //std::cout<<"omega of track "<<track->getOmega()<<"/"<<TMath::Pi()/2-atan(trackAtIP->getTanLambda())<<"/"<<m_innerBField * m_const_a/fabs(trackAtIP->getOmega())<<std::endl;
	    }
	    m_track_x_atIP->push_back(trackAtIP->getReferencePoint()[0]);
	    m_track_y_atIP->push_back(trackAtIP->getReferencePoint()[1]);
	    m_track_z_atIP->push_back(trackAtIP->getReferencePoint()[2]);
	    m_track_Phi_atIP->push_back(trackAtIP->getPhi());
	    m_track_Theta_atIP->push_back(TMath::Pi()/2-atan(trackAtIP->getTanLambda()));
	    m_track_pt_atIP->push_back(m_innerBField * m_const_a/fabs(trackAtIP->getOmega()));//omega is in terms of mm
	    m_track_p_atIP->push_back(m_innerBField * m_const_a/fabs(trackAtIP->getOmega())/cos(atan(trackAtIP->getTanLambda())));
	    
	    float nExpectedTrackerHits=0;
	    float p=m_innerBField * m_const_a/fabs(trackAtIP->getOmega())/cos(atan(trackAtIP->getTanLambda()));
	    float pT=m_innerBField * m_const_a/fabs(trackAtIP->getOmega());
	    float pZ=sqrt(p*p-pT*pT);
	    
	    if (pZ < (m_trackerZmax / m_trackerOuterR * pT)){
	      const float innerExpectedHitRadius=std::max(m_trackerInnerR, track->getRadiusOfInnermostHit());
	      const float frac=(m_trackerOuterR - innerExpectedHitRadius) / (m_trackerOuterR - m_trackerInnerR);
	      nExpectedTrackerHits = m_barrelTrackerLayers * frac; 
	    } 
	    if ((pZ <= (m_trackerZmax / m_trackerInnerR * pT)) && (pZ >= (m_trackerZmax / m_trackerOuterR * pT))){
	      const float innerExpectedHitRadius=std::max(m_trackerInnerR, track->getRadiusOfInnermostHit());
	      const float frac=(m_trackerZmax * pT / pZ - innerExpectedHitRadius) / (m_trackerOuterR - innerExpectedHitRadius);
	      nExpectedTrackerHits = frac * m_barrelTrackerLayers;
	    }
	    m_track_nExpectedTrackerHits->push_back(nExpectedTrackerHits);
	    if( (pZ < (m_trackerZmax / m_trackerOuterR * pT))   && ((pZ <= (m_trackerZmax / m_trackerInnerR * pT)) && (pZ >= (m_trackerZmax / m_trackerOuterR * pT))) ){
	      std::cout<<"HONESTLY BOTH CONDITIONS "<<std::endl;
	    }
	  }else{
	    m_track_nExpectedTrackerHits->push_back(-1);
	    m_track_x_atIP->push_back(-1);
	    m_track_y_atIP->push_back(-1);
	    m_track_z_atIP->push_back(-1);
	    m_track_Phi_atIP->push_back(-1);
	    m_track_Theta_atIP->push_back(-1);
	    m_track_pt_atIP->push_back(-1);
	    m_track_p_atIP->push_back(-1);
	  }
	  if(trackAtCalo!=0){
	    m_track_x_atCalo->push_back(trackAtCalo->getReferencePoint()[0]);
	    m_track_y_atCalo->push_back(trackAtCalo->getReferencePoint()[1]);
	    m_track_z_atCalo->push_back(trackAtCalo->getReferencePoint()[2]);
	    m_track_Phi_atCalo->push_back(trackAtCalo->getPhi());
	    m_track_Theta_atCalo->push_back(TMath::Pi()/2-atan(trackAtCalo->getTanLambda()));
	    m_track_pt_atCalo->push_back(m_innerBField * m_const_a/fabs(trackAtCalo->getOmega()));
	    m_track_p_atCalo->push_back(m_innerBField * m_const_a/fabs(trackAtCalo->getOmega())/cos(atan(trackAtCalo->getTanLambda())));
	    
	    float track_ptatcalo=m_innerBField * m_const_a/fabs(trackAtCalo->getOmega())/cos(atan(trackAtCalo->getTanLambda()));
	    TVector3 trackMomentumAtCaloDir(track_ptatcalo*std::cos(trackAtCalo->getPhi()),track_ptatcalo*std::sin(trackAtCalo->getPhi()), track_ptatcalo*trackAtCalo->getTanLambda());
	    m_track_px_atCalo->push_back(trackMomentumAtCaloDir[0]);
	    m_track_py_atCalo->push_back(trackMomentumAtCaloDir[1]);
	    m_track_pz_atCalo->push_back(trackMomentumAtCaloDir[2]);
	    trackMomentumAtCaloDir.SetMag(1.0);//track direction at calo face is momentum unit vector at caloface
	    if(clucol!=NULL){
	      int index_cluster_track=-1;
	      float angle_min=10;
	      //float dirAngleForPosAngle_min=10;
	      float EOverP=-10;
	      float PFdistance_cluster_track_min=-1;
	      int index_cluster_track_tinyHad=-1;
	      float PFdistance_cluster_track_min_tinyHad=-1;
	      for(int c=0;c<clucol->getNumberOfElements();c++){
		Cluster* pandoraclus1 = dynamic_cast<EVENT::Cluster*>(clucol->getElementAt(c));
		//calculate at that point the direction of the cluster in the manner PFAPandora does -->will be used later on for every single track check	
		double HadEnergy=pandoraclus1->getSubdetectorEnergies()[1];//1 is HCAL, 0 is ECAL       
		if(HadEnergy>0.2){
		  //check if high cluster or the low cluster point -> use subdetectors for that
		  for(unsigned int i_calohit1=0;i_calohit1< pandoraclus1->getCalorimeterHits().size();i_calohit1++){
		    TVector3 distancemomentum(pandoraclus1->getCalorimeterHits()[i_calohit1]->getPosition()[0]-trackAtCalo->getReferencePoint()[0],pandoraclus1->getCalorimeterHits()[i_calohit1]->getPosition()[1]-trackAtCalo->getReferencePoint()[1],pandoraclus1->getCalorimeterHits()[i_calohit1]->getPosition()[2]-trackAtCalo->getReferencePoint()[2]);
		    if(PFdistance_cluster_track_min<0){//first time we get to that check
		      PFdistance_cluster_track_min=(trackMomentumAtCaloDir.Cross(distancemomentum)).Mag();
		      index_cluster_track=c;
		    }else if((trackMomentumAtCaloDir.Cross(distancemomentum)).Mag()<PFdistance_cluster_track_min){
		      PFdistance_cluster_track_min=(trackMomentumAtCaloDir.Cross(distancemomentum)).Mag();
		      index_cluster_track=c;
		    }
		  }	  
		}else{
		  for(unsigned int i_calohit1=0;i_calohit1< pandoraclus1->getCalorimeterHits().size();i_calohit1++){
		    TVector3 distancemomentum(pandoraclus1->getCalorimeterHits()[i_calohit1]->getPosition()[0]-trackAtCalo->getReferencePoint()[0],pandoraclus1->getCalorimeterHits()[i_calohit1]->getPosition()[1]-trackAtCalo->getReferencePoint()[1],pandoraclus1->getCalorimeterHits()[i_calohit1]->getPosition()[2]-trackAtCalo->getReferencePoint()[2]);
		    if(PFdistance_cluster_track_min_tinyHad<0){//first time we get to that check
		      PFdistance_cluster_track_min_tinyHad=(trackMomentumAtCaloDir.Cross(distancemomentum)).Mag();
		      index_cluster_track_tinyHad=c;
		    }else if((trackMomentumAtCaloDir.Cross(distancemomentum)).Mag()<PFdistance_cluster_track_min){
		      PFdistance_cluster_track_min_tinyHad=(trackMomentumAtCaloDir.Cross(distancemomentum)).Mag();
		      index_cluster_track_tinyHad=c;
		    }
		  }
		}
		float angle = acos((trackAtCalo->getReferencePoint()[0]*pandoraclus1->getPosition()[0]+trackAtCalo->getReferencePoint()[1]*pandoraclus1->getPosition()[1]+trackAtCalo->getReferencePoint()[2]*pandoraclus1->getPosition()[2])/(sqrt(pandoraclus1->getPosition()[0]*pandoraclus1->getPosition()[0]+pandoraclus1->getPosition()[1]*pandoraclus1->getPosition()[1]+pandoraclus1->getPosition()[2]*pandoraclus1->getPosition()[2])*sqrt(trackAtCalo->getReferencePoint()[0]*trackAtCalo->getReferencePoint()[0]+trackAtCalo->getReferencePoint()[1]*trackAtCalo->getReferencePoint()[1]+trackAtCalo->getReferencePoint()[2]*trackAtCalo->getReferencePoint()[2])));
		if(angle<angle_min){
		  angle_min=angle;
		  EOverP= pandoraclus1->getEnergy()/(m_innerBField * m_const_a/(cos(atan(trackAtCalo->getTanLambda()))*fabs(trackAtCalo->getOmega()))) ;
		}
	      }//loop over clusters done
	      if(index_cluster_track>-1){//found one cluster - seems collection can be empty
		Cluster* pandoraclus_sel = dynamic_cast<EVENT::Cluster*> (clucol->getElementAt(index_cluster_track));
		TVector3 selclusterdirection(0.,0.,0.);
		for(unsigned int i_calohit_dir=0;i_calohit_dir< pandoraclus_sel->getCalorimeterHits().size();i_calohit_dir++){
		  TVector3 calohitdir(pandoraclus_sel->getCalorimeterHits()[i_calohit_dir]->getPosition()[0],pandoraclus_sel->getCalorimeterHits()[i_calohit_dir]->getPosition()[1],pandoraclus_sel->getCalorimeterHits()[i_calohit_dir]->getPosition()[2]);
		  calohitdir.SetMag(1.0);
		  selclusterdirection+=calohitdir;
		}
		selclusterdirection.SetMag(1.0);
		//angle between selected (min distance) cluster and track momentum
		float min_parallel_distance = -1;
		//do minimum parallel cluster distance
		for(unsigned int i_calohit1=0;i_calohit1< pandoraclus_sel->getCalorimeterHits().size();i_calohit1++){
		  TVector3 distancemomentum(pandoraclus_sel->getCalorimeterHits()[i_calohit1]->getPosition()[0]-trackAtCalo->getReferencePoint()[0],pandoraclus_sel->getCalorimeterHits()[i_calohit1]->getPosition()[1]-trackAtCalo->getReferencePoint()[1],pandoraclus_sel->getCalorimeterHits()[i_calohit1]->getPosition()[2]-trackAtCalo->getReferencePoint()[2]);
		  if(min_parallel_distance<0){
		    min_parallel_distance=fabs(trackMomentumAtCaloDir.Dot(distancemomentum));
		  }else if(fabs(trackMomentumAtCaloDir.Dot(distancemomentum))<min_parallel_distance){
		    min_parallel_distance=fabs(trackMomentumAtCaloDir.Dot(distancemomentum));
		  }
		}
		m_track_cluster_minPFADist_clusterEMEnergy->push_back(pandoraclus_sel->getSubdetectorEnergies()[0]);
		m_track_cluster_minPFADist_clusterHadEnergy->push_back(pandoraclus_sel->getSubdetectorEnergies()[1]);
		m_track_cluster_minPFAdistance_atCalo->push_back(PFdistance_cluster_track_min);
		m_track_cluster_DirAngle_minPFAdistance_atCalo->push_back(acos(trackMomentumAtCaloDir.Dot(selclusterdirection)));
		m_track_cluster_minPFADist_EOverP->push_back(pandoraclus_sel->getEnergy()/(track_ptatcalo/cos(atan(trackAtCalo->getTanLambda()))));
		m_track_cluster_min_parPFAdistance_atCalo->push_back(min_parallel_distance);
		m_track_cluster_PosAngle_min_atCalo->push_back(angle_min);
		m_track_minDist_cluster_EOverP->push_back(EOverP);
	      }else{
		m_track_cluster_PosAngle_min_atCalo->push_back(-1);
		m_track_minDist_cluster_EOverP->push_back(-1);
		m_track_cluster_minPFADist_EOverP->push_back(-1);
		m_track_cluster_minPFADist_clusterEMEnergy->push_back(-1);
		m_track_cluster_minPFADist_clusterHadEnergy->push_back(-1);
		m_track_cluster_minPFAdistance_atCalo->push_back(-1);
		m_track_cluster_min_parPFAdistance_atCalo->push_back(-1);
		m_track_cluster_DirAngle_minPFAdistance_atCalo->push_back(-1);
	      }
	      if(index_cluster_track_tinyHad>-1){//found one cluster - seems collection can be empty
		Cluster* pandoraclus_sel = dynamic_cast<EVENT::Cluster*> (clucol->getElementAt(index_cluster_track_tinyHad));
		TVector3 selclusterdirection(0.,0.,0.);
		for(unsigned int i_calohit_dir=0;i_calohit_dir< pandoraclus_sel->getCalorimeterHits().size();i_calohit_dir++){
		  TVector3 calohitdir(pandoraclus_sel->getCalorimeterHits()[i_calohit_dir]->getPosition()[0],pandoraclus_sel->getCalorimeterHits()[i_calohit_dir]->getPosition()[1],pandoraclus_sel->getCalorimeterHits()[i_calohit_dir]->getPosition()[2]);
		  calohitdir.SetMag(1.0);
		  selclusterdirection+=calohitdir;
		}
		selclusterdirection.SetMag(1.0);
		//angle between selected (min distance) cluster and track momentum
		float min_parallel_distance = -1;
		//do minimum parallel cluster distance
		for(unsigned int i_calohit1=0;i_calohit1< pandoraclus_sel->getCalorimeterHits().size();i_calohit1++){
		  TVector3 distancemomentum(pandoraclus_sel->getCalorimeterHits()[i_calohit1]->getPosition()[0]-trackAtCalo->getReferencePoint()[0],pandoraclus_sel->getCalorimeterHits()[i_calohit1]->getPosition()[1]-trackAtCalo->getReferencePoint()[1],pandoraclus_sel->getCalorimeterHits()[i_calohit1]->getPosition()[2]-trackAtCalo->getReferencePoint()[2]);
		  if(min_parallel_distance<0){
		    min_parallel_distance=fabs(trackMomentumAtCaloDir.Dot(distancemomentum));
		  }else if(fabs(trackMomentumAtCaloDir.Dot(distancemomentum))<min_parallel_distance){
		    min_parallel_distance=fabs(trackMomentumAtCaloDir.Dot(distancemomentum));
		  }
		}
		m_track_cluster_minPFADist_clusterEMEnergy_tinyHad->push_back(pandoraclus_sel->getSubdetectorEnergies()[0]);
		m_track_cluster_minPFADist_clusterHadEnergy_tinyHad->push_back(pandoraclus_sel->getSubdetectorEnergies()[1]);
		m_track_cluster_minPFAdistance_atCalo_tinyHad->push_back(PFdistance_cluster_track_min_tinyHad);
		m_track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad->push_back(acos(trackMomentumAtCaloDir.Dot(selclusterdirection)));
		m_track_cluster_minPFADist_EOverP_tinyHad->push_back(pandoraclus_sel->getEnergy()/(track_ptatcalo/cos(atan(trackAtCalo->getTanLambda()))));
		m_track_cluster_min_parPFAdistance_atCalo_tinyHad->push_back(min_parallel_distance);
	      }else{//no low HadEnergy cluster found
		m_track_cluster_minPFADist_EOverP_tinyHad->push_back(-1);
		m_track_cluster_minPFADist_clusterEMEnergy_tinyHad->push_back(-1);
		m_track_cluster_minPFADist_clusterHadEnergy_tinyHad->push_back(-1);
		m_track_cluster_minPFAdistance_atCalo_tinyHad->push_back(-1);
		m_track_cluster_min_parPFAdistance_atCalo_tinyHad->push_back(-1);
		m_track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad->push_back(-1);
	      }
	    }//cluster exits
	  }else{
	    m_track_cluster_minPFADist_EOverP->push_back(-1);
	    m_track_cluster_minPFADist_clusterEMEnergy->push_back(-1);
	    m_track_cluster_minPFADist_clusterHadEnergy->push_back(-1);;
	    m_track_cluster_minPFAdistance_atCalo->push_back(-1);
	    m_track_cluster_min_parPFAdistance_atCalo->push_back(-1);
	    m_track_cluster_DirAngle_minPFAdistance_atCalo->push_back(-1);
	    m_track_cluster_minPFADist_EOverP_tinyHad->push_back(-1);
	    m_track_cluster_minPFADist_clusterEMEnergy_tinyHad->push_back(-1);
	    m_track_cluster_minPFADist_clusterHadEnergy_tinyHad->push_back(-1);;
	    m_track_cluster_minPFAdistance_atCalo_tinyHad->push_back(-1);
	    m_track_cluster_min_parPFAdistance_atCalo_tinyHad->push_back(-1);
	    m_track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad->push_back(-1);
	    
	    m_track_cluster_PosAngle_min_atCalo->push_back(-1);
	    m_track_minDist_cluster_EOverP->push_back(-1);
	    m_track_x_atCalo->push_back(-1);
	    m_track_y_atCalo->push_back(-1);
	    m_track_z_atCalo->push_back(-1);
	    m_track_px_atCalo->push_back(-1);
	    m_track_py_atCalo->push_back(-1);
	    m_track_pz_atCalo->push_back(-1);
	    m_track_Phi_atCalo->push_back(-1);
	    m_track_Theta_atCalo->push_back(-1);
	    m_track_pt_atCalo->push_back(-1);
	    m_track_p_atCalo->push_back(-1);
	  }
	}
      }else{
	std::cout<<"WHERE IS THE TRACK"<<std::endl;
      }
    }

    //std::cout<<"cut 5"<<std::endl;

    bool nohits_present = true;

    if(!m_runTauMode && !nohits_present){
      //std::cout<<"what goes on here"<<std::endl;
      if(clucol!=NULL){
	//std::cout<<"cut 5a"<<std::endl;
	//std::cout<<"try to fill cluster collection"<<std::endl;
	//std::cout<<"before the cluster collection things "<<std::endl;
	for(int i=0;i<clucol->getNumberOfElements();i++){
	  //std::cout<<"cut 5b"<<std::endl;
	  //std::cout<<"go into the cluster collection things "<<i<<std::endl;
	  //std::cout<<"cluster col "<<i<<std::endl;
	  Cluster* pandoraclus = dynamic_cast<EVENT::Cluster*>(clucol->getElementAt(i));
	  m_cluster_x->push_back(pandoraclus->getPosition()[0]);
	  m_cluster_y->push_back(pandoraclus->getPosition()[1]);
	  m_cluster_z->push_back(pandoraclus->getPosition()[2]);
	  TVector3 clusterdirection(0.,0.,0.);
	  //calculate at that point the direction of the cluster in the manner PFAPandora does -->will be used later on for every single track check

	  for(unsigned int i_calohit=0;i_calohit< pandoraclus->getCalorimeterHits().size();i_calohit++){
	    if(m_fillClusterHits){
	      m_cluster_hit_x->push_back(pandoraclus->getCalorimeterHits()[i_calohit]->getPosition()[0]);
	      m_cluster_hit_y->push_back(pandoraclus->getCalorimeterHits()[i_calohit]->getPosition()[1]);
	      m_cluster_hit_z->push_back(pandoraclus->getCalorimeterHits()[i_calohit]->getPosition()[2]);
	      m_cluster_hit_E->push_back(pandoraclus->getCalorimeterHits()[i_calohit]->getEnergy());
	      m_cluster_hit_index->push_back(i);
	      
	      const CHT cht=pandoraclus->getCalorimeterHits()[i_calohit]->getType();	
	      if(cht.is(CHT::ecal)){////id 1=ecal, 2 hcal,3 muon
		m_cluster_hit_type->push_back(1); 
	      }else if(cht.is(CHT::hcal)){//hcal
		m_cluster_hit_type->push_back(2);
	      }else if (cht.is(CHT::yoke)){
		clustermuonhits.insert(pandoraclus->getCalorimeterHits()[i_calohit]);
		m_cluster_hit_type->push_back(3);
	      }   
	    }
	    //vHitE[i_calohit]=(float)pandoraclus->getCalorimeterHits()[i_calohit]->getEnergy();
	    //vHitZ[i_calohit]=(float)pandoraclus->getCalorimeterHits()[i_calohit]->getPosition()[0];
	    //vHitY[i_calohit]=(float)pandoraclus->getCalorimeterHits()[i_calohit]->getPosition()[1];
	    //vHitZ[i_calohit]=(float)pandoraclus->getCalorimeterHits()[i_calohit]->getPosition()[2];
	    TVector3 calohitdirection(pandoraclus->getCalorimeterHits()[i_calohit]->getPosition()[0],pandoraclus->getCalorimeterHits()[i_calohit]->getPosition()[1],pandoraclus->getCalorimeterHits()[i_calohit]->getPosition()[2]);
	    calohitdirection.SetMag(1.0);
	    clusterdirection+=calohitdirection;
	  }
	  //after loop over clusters, check what is contained in clustermuonhits
	  //std::cout<<"size of cluster muon hits after cluster "<<i<<"/"<<clustermuonhits.size()<<std::endl;
	  //std::set<CalorimeterHit*>::iterator mhitIt;
	  //unsigned int counter=0;
	  //for(mhitIt=clustermuonhits.begin();mhitIt!=clustermuonhits.end();mhitIt++){
	  //std::cout<<"x/y/z/E "<<counter<<"/"<<(*mhitIt)->getPosition()[0]<<"/"<<(*mhitIt)->getPosition()[1]<<"/"<<(*mhitIt)->getPosition()[2]<<"/"<<(*mhitIt)->getEnergy()<<std::endl;
	  //counter+=1;
	  //}
	  //ClusterShapes *const pClusterShapes= new ClusterShapes((int)pandoraclus->getCalorimeterHits().size(), vHitE, vHitX, vHitY, vHitZ);
	  //delete vHitE;
	  //delete vHitX;
	  //delete vHitY;
	  //delete vHitZ;
	  m_cluster_dir_x->push_back(clusterdirection[0]);
	  m_cluster_dir_y->push_back(clusterdirection[1]);
	  m_cluster_dir_z->push_back(clusterdirection[2]);
	  m_cluster_energy->push_back(pandoraclus->getEnergy());
	  m_cluster_energyError->push_back(pandoraclus->getEnergyError());
	  m_cluster_iPhi->push_back(pandoraclus->getIPhi());
	  m_cluster_iTheta->push_back(pandoraclus->getITheta());
	  //std::cout<<"cut 5c"<<std::endl;
	  //std::cout<<"cluster before track collection thing "<<i<<std::endl;
	  if(trackcol!=NULL){
	    //std::cout<<"cut 5d"<<std::endl;
	    clusterdirection.SetMag(1.0);
	    //std::cout<<"cut 5d00 "<<pandoraclus->getCalorimeterHits().size()<<std::endl;
	    TVector3 clusterposition(pandoraclus->getPosition()[0],pandoraclus->getPosition()[1],pandoraclus->getPosition()[2]);
	    clusterposition.SetMag(1.0);
	    //std::cout<<"cut 5d01"<<std::endl;
	    float angle_min=10;
	    float dirAngleForPosAngle_min=10;
	    int index_track_PFminDist=-1;
	    float par_PFdistance_max=1000.;
	    float PFdistance_min=-1;// parallel distance 1000,lowenergy 02 mintrack CosAngle 0, maxSearchLayer 9, maxTrackClusterDisatnce 10
	    for(int t=0;t<trackcol->getNumberOfElements();t++){
	      //std::cout<<"cut 5d1"<<std::endl;
	      //std::cout<<"in track collection thing "<<t<<" "<<trackcol->getNumberOfElements()<<std::endl;
	      Track* track1 = dynamic_cast<EVENT::Track*>(trackcol->getElementAt(t));
	      //std::cout<<"after track get "<<t<<std::endl;
	      const EVENT::TrackState* trackAtCalo1=track1->getTrackState(EVENT::TrackState::AtCalorimeter);
	      //std::cout<<"after track at Calo "<<t<<std::endl;
	      if(trackAtCalo1!=NULL){
		//std::cout<<"in track at calo "<<t<<std::endl;
		float trackptatcalo1=m_innerBField * m_const_a/fabs(trackAtCalo1->getOmega());
		TVector3 trackMomentumAtCalo1Dir(trackptatcalo1*std::cos(trackAtCalo1->getPhi()),trackptatcalo1*std::sin(trackAtCalo1->getPhi()), trackptatcalo1*trackAtCalo1->getTanLambda());
		trackMomentumAtCalo1Dir.SetMag(1.0);
		//unit vector of track direction unitvector and vector of position difference vector
		if((trackMomentumAtCalo1Dir.Dot(clusterdirection))<0.0){
		  continue;
		}
		for(unsigned int i_calohit1=0;i_calohit1< pandoraclus->getCalorimeterHits().size();i_calohit1++){
		  //std::cout<<"loop over calohits "<<i_calohit1<<std::endl;
		  TVector3 distancemomentum(pandoraclus->getCalorimeterHits()[i_calohit1]->getPosition()[0]-trackAtCalo1->getReferencePoint()[0],pandoraclus->getCalorimeterHits()[i_calohit1]->getPosition()[1]-trackAtCalo1->getReferencePoint()[1],pandoraclus->getCalorimeterHits()[i_calohit1]->getPosition()[2]-trackAtCalo1->getReferencePoint()[2]);
		  if(fabs(trackMomentumAtCalo1Dir.Dot(distancemomentum))>par_PFdistance_max){
		    //std::cout<<"big parallel value distance value "<<fabs(trackMomentumAtCalo1Dir.Dot(distancemomentum))<<std::endl;
		    continue;
		  }
		  if(PFdistance_min<0){//first time we get to that check
		    PFdistance_min=(trackMomentumAtCalo1Dir.Cross(distancemomentum)).Mag();
		    index_track_PFminDist=t;
		  }else if((trackMomentumAtCalo1Dir.Cross(distancemomentum)).Mag()<PFdistance_min){
		    PFdistance_min=(trackMomentumAtCalo1Dir.Cross(distancemomentum)).Mag();
		    index_track_PFminDist=t;
		  }
		}
		float angle = acos((trackAtCalo1->getReferencePoint()[0]*pandoraclus->getPosition()[0]+trackAtCalo1->getReferencePoint()[1]*pandoraclus->getPosition()[1]+trackAtCalo1->getReferencePoint()[2]*pandoraclus->getPosition()[2])/(sqrt(pandoraclus->getPosition()[0]*pandoraclus->getPosition()[0]+pandoraclus->getPosition()[1]*pandoraclus->getPosition()[1]+pandoraclus->getPosition()[2]*pandoraclus->getPosition()[2])*sqrt(trackAtCalo1->getReferencePoint()[0]*trackAtCalo1->getReferencePoint()[0]+trackAtCalo1->getReferencePoint()[1]*trackAtCalo1->getReferencePoint()[1]+trackAtCalo1->getReferencePoint()[2]*trackAtCalo1->getReferencePoint()[2])));
		if(angle<angle_min){
		  //std::cout<<"in track angle stuff"<<std::endl;
		  angle_min=angle;
		  //std::cout<<"abs track "<<trackAtCalo1->getPhi()<<"/"<<atan2(trackAtCalo1->getReferencePoint()[1],trackAtCalo1->getReferencePoint()[0])<<std::endl;
		  dirAngleForPosAngle_min=acos(cos(pandoraclus->getIPhi())*sin(pandoraclus->getITheta()) * cos(trackAtCalo1->getPhi())*sin(TMath::Pi()/2-atan(trackAtCalo1->getTanLambda()))
					       + sin(pandoraclus->getIPhi())*sin(pandoraclus->getITheta()) * sin(trackAtCalo1->getPhi())*sin(TMath::Pi()/2-atan(trackAtCalo1->getTanLambda()))
					       + cos(pandoraclus->getITheta())*cos(TMath::Pi()/2-atan(trackAtCalo1->getTanLambda())));
		}	  
	      }
	    }
	    //std::cout<<"in track before distance thing "<<std::endl;
	    if(index_track_PFminDist>-1){
	      //std::cout<<"in track collection distance thing "<<std::endl;
	      Track* seltrack = dynamic_cast<EVENT::Track*>(trackcol->getElementAt(index_track_PFminDist));
	      const EVENT::TrackState* seltrackAtCalo=seltrack->getTrackState(EVENT::TrackState::AtCalorimeter);
	      float seltrackptatcalo=m_innerBField * m_const_a/fabs(seltrackAtCalo->getOmega());
	      float seltrackpatcalo=m_innerBField * m_const_a/fabs(seltrackAtCalo->getOmega())/cos(atan(seltrackAtCalo->getTanLambda()));
	      TVector3 seltrackMomentumAtCaloDir(seltrackptatcalo*std::cos(seltrackAtCalo->getPhi()),seltrackptatcalo*std::sin(seltrackAtCalo->getPhi()), seltrackptatcalo*seltrackAtCalo->getTanLambda());
	      seltrackMomentumAtCaloDir.SetMag(1.0);
	      m_cluster_track_minPFADist_cluster_EOverP->push_back(pandoraclus->getEnergy()/seltrackpatcalo);
	      m_cluster_track_minPFADist_cluster_Angle->push_back(acos(seltrackMomentumAtCaloDir.Dot(clusterdirection)));//for min PFA like distance cluster direction and track direction angle
	      m_cluster_track_minPFADist_cluster_TCDistance->push_back(PFdistance_min);
	      m_cluster_track_PosAngle_min->push_back(angle_min);
	      m_cluster_track_DirAngle->push_back(dirAngleForPosAngle_min);
	    }else{	    
	      //std::cout<<"in distance else distance thing "<<std::endl;
	      m_cluster_track_minPFADist_cluster_EOverP->push_back(-1.);
	      m_cluster_track_minPFADist_cluster_Angle->push_back(-1.);
	      m_cluster_track_minPFADist_cluster_TCDistance->push_back(-1.);
	      m_cluster_track_PosAngle_min->push_back(-1);
	      m_cluster_track_DirAngle->push_back(-1);
	    }
	    //std::cout<<"cluster "<<i<<" tracks "<<trackcol->getNumberOfElements()<< "pos min/angle "<<angle_min<<"/"<<dirAngleForPosAngle_min<<std::endl;
	  }
	    //std::cout<<"cut 5e"<<std::endl;
	  double ecal_energy=0;
	  int necal_hits=0;	  
	  double hcal_energy=0;
	  int nhcal_hits=0;
	  int nmuon_hits=0;
	  float En_ECAL_Barrel=0;
	  float En_ECAL_Endcap=0;
	  float En_ECAL_else=0;
	  float En_HCAL_Barrel=0;
	  float En_HCAL_Endcap=0;
	  float En_HCAL_else=0;
	  float En_MU=0;
	  //std::cout<<"in cluster hit calorimeter stuff "<<std::endl;
	  for(unsigned int l=0; l<pandoraclus->getCalorimeterHits().size();l++){
	    //std::cout<<"cut 5d3"<<std::endl;
	    const CHT cht=pandoraclus->getCalorimeterHits()[l]->getType();
	    //int types=pandoraclus->getCalorimeterHits()[l]->getType();
	    //static const int fCaloType =     1 ;
	    //static const int fCaloID   =    10 ;
	    //static const int fLayout   =  1000 ;
	    //static const int fLayer    = 10000 ;
	    //int caloLayer=types/fLayer;
	    //int caloLayout=(types-caloLayer*fLayer)/fLayout;
	    //int caloID=(types-caloLayer*fLayer-caloLayout*fLayout)/fCaloID;
	    //int caloType=(types-caloLayer*fLayer-caloLayout*fLayout-caloID*fCaloID)/fCaloType;	
	    if(cht.is(CHT::ecal)){////id 1=ecal, 2 hcal
	      ecal_energy+=pandoraclus->getCalorimeterHits()[l]->getEnergy();
	      necal_hits+=1;
	      if(cht.is(CHT::barrel)){
	      double calohitZ=pandoraclus->getCalorimeterHits()[l]->getPosition()[2];
	      if(fabs(calohitZ)>m_ECAL_barrelZ_max && (fabs(calohitZ)<m_ECAL_ringZ_min || fabs(calohitZ)>m_ECAL_ringZ_max) && (fabs(calohitZ)<m_ECAL_endcapZ_min || fabs(calohitZ)>m_ECAL_endcapZ_max)){
		//std::cout<<"E calohit Z "<<calohitZ<<" bZ "<<m_ECAL_barrelZ_max<<" rZs "<<m_ECAL_ringZ_min<<"/"<<m_ECAL_ringZ_max<<" eZs "<<m_ECAL_endcapZ_min<<"/"<<m_ECAL_endcapZ_max<<" b/r/e type "<<cht.is(CHT::barrel)<<"/"<<cht.is(CHT::plug)<<"/"<<cht.is(CHT::endcap)<<std::endl;
	      }
		En_ECAL_Barrel+=pandoraclus->getCalorimeterHits()[l]->getEnergy();	       
	      }else if (cht.is(CHT::endcap)){
		En_ECAL_Endcap+=pandoraclus->getCalorimeterHits()[l]->getEnergy();
	      }else{
		En_ECAL_else+=pandoraclus->getCalorimeterHits()[l]->getEnergy();
	      }
	    }else if(cht.is(CHT::hcal)){//hcal
	      //double hcalohitZ=pandoraclus->getCalorimeterHits()[l]->getPosition()[2];
	      hcal_energy+=pandoraclus->getCalorimeterHits()[l]->getEnergy();
	      nhcal_hits+=1;
	      if(cht.is(CHT::barrel)){
		En_HCAL_Barrel+=pandoraclus->getCalorimeterHits()[l]->getEnergy();    
	      }else if(cht.is(CHT::endcap)){
		En_HCAL_Endcap+=pandoraclus->getCalorimeterHits()[l]->getEnergy();
	      }else{
		En_HCAL_else+=pandoraclus->getCalorimeterHits()[l]->getEnergy();
	      }
	    }else if (cht.is(CHT::yoke)){
	      double mcalohitZ=pandoraclus->getCalorimeterHits()[l]->getPosition()[2];
	      if(cht.is(CHT::barrel) && fabs(mcalohitZ)>m_MUON_barrelZ_max){
		std::cout<<"barrel hit but outside of barrel volume"<<std::endl;
	      }
	      if(cht.is(CHT::endcap) && fabs(mcalohitZ)<m_MUON_barrelZ_max){
		std::cout<<"endcap hit but outside of endcap volume"<<std::endl;
	      }
	      En_MU+=pandoraclus->getCalorimeterHits()[l]->getEnergy();
	      nmuon_hits+=1;
	    }
	    //std::cout<<"cut 5e"<<std::endl;
	  }//hit loop
	  m_cluster_E_EB->push_back(En_ECAL_Barrel);
	  m_cluster_E_EE->push_back(En_ECAL_Endcap);
	  m_cluster_E_EO->push_back(En_ECAL_else);
	  m_cluster_E_HB->push_back(En_HCAL_Barrel);
	  m_cluster_E_HE->push_back(En_HCAL_Endcap);
	  m_cluster_E_HO->push_back(En_HCAL_else);
	  m_cluster_E_MU->push_back(En_MU);
	  m_cluster_ECAL_energy->push_back(ecal_energy);
	  m_cluster_nHits_ECAL->push_back(necal_hits);
	  m_cluster_HCAL_energy->push_back(hcal_energy);
	  m_cluster_nHits_HCAL->push_back(nhcal_hits);
	  m_cluster_nHits_MUON->push_back(nmuon_hits);
	  //std::cout<<"cut 5f"<<std::endl;
	}
      }
   
      //std::cout<<"end of cluster collection"<<std::endl;

      if(m_fillClusterHits){
	//std::cout<<"cut 5g"<<std::endl;
	LCCollection* muonhitcol = NULL;
	muonhitcol = evt->getCollection(m_inputMuonHitCollectionName);
	std::set<CalorimeterHit*> samuonhits;
	if(muonhitcol!=NULL){
	  for(int i=0;i<muonhitcol->getNumberOfElements();i++){
	    CalorimeterHit* muonhit = dynamic_cast<CalorimeterHit*>(muonhitcol->getElementAt(i));
	    samuonhits.insert(muonhit);
	  }
	  if(clustermuonhits.size()!=samuonhits.size()){
	    std::set<CalorimeterHit*>::iterator samhitIt;
	    std::set<CalorimeterHit*>::iterator mhitIt;
	    for(samhitIt=samuonhits.begin();samhitIt!=samuonhits.end();samhitIt++){
	      bool hit_in_set=false;
	      for(mhitIt=clustermuonhits.begin();mhitIt!=clustermuonhits.end();mhitIt++){
		if((*mhitIt)==(*samhitIt)){
		  hit_in_set=true;
		  break;
		}
	      }
	      if(!hit_in_set){
		m_no_cluster_hit_x->push_back((*samhitIt)->getPosition()[0]);
		m_no_cluster_hit_y->push_back((*samhitIt)->getPosition()[1]);
		m_no_cluster_hit_z->push_back((*samhitIt)->getPosition()[2]);
		m_no_cluster_hit_E->push_back((*samhitIt)->getEnergy());
	      //std::cout<<"SAHIT x/y/z/E "<<(*samhitIt)->getPosition()[0]<<"/"<<(*samhitIt)->getPosition()[1]<<"/"<<(*samhitIt)->getPosition()[2]<<"/"<<(*samhitIt)->getEnergy()<<std::endl;
	      //for(unsigned int i=0;i< m_true_CosTheta->size();i++){
	      //    std::cout<<"true particles stat/ID/costheta "<< (*m_true_GenStatus)[i]<<"/"<<(*m_true_PDGID)[i]<<"/"<< (*m_true_CosTheta)[i]<<std::endl;
	      //}
	      }
	    }
	  }
	  //std::cout<<"size of cluster muon hits after mu hits "<<clustermuonhits.size()<<std::endl;
	  //std::set<CalorimeterHit*>::iterator mhitIt;
	  //unsigned int counter=0;
	  //for(mhitIt=clustermuonhits.begin();mhitIt!=clustermuonhits.end();mhitIt++){
	  //std::cout<<"x/y/z/E "<<counter<<"/"<<(*mhitIt)->getPosition()[0]<<"/"<<(*mhitIt)->getPosition()[1]<<"/"<<(*mhitIt)->getPosition()[2]<<"/"<<(*mhitIt)->getEnergy()<<std::endl;
	  //counter+=1;
	  //}
	}
      }//if cluster hits should be filled
      if(m_fillSimHits){
	//std::cout<<"cut 5h"<<std::endl;
	LCCollection* muonsimhitcol = NULL;
	muonsimhitcol = evt->getCollection(m_inputSimCaloHitCollectionName);
	if(muonsimhitcol!=NULL){
	  for(int mu_sh=0;mu_sh<muonsimhitcol->getNumberOfElements();mu_sh++){	
	    SimCalorimeterHit* muonsimhit = dynamic_cast<SimCalorimeterHit*>(muonsimhitcol->getElementAt(mu_sh));
	    m_muon_b_simhit_x->push_back(muonsimhit->getPosition()[0]);
	    m_muon_b_simhit_y->push_back(muonsimhit->getPosition()[1]);
	    m_muon_b_simhit_z->push_back(muonsimhit->getPosition()[2]);
	    m_muon_b_simhit_E->push_back(muonsimhit->getEnergy());
	    //std::cout<<"muon sim hit stuff "<<(*m_muon_b_simhit_E)[m_muon_b_simhit_E->size()-1]<<std::endl;
	  }
	}//else{
	//std::cout<<"seems like we have nothing in muon sim hits"<<std::endl;
	//}
      }
    }//only run tracks and clusters for non tau events ->>validation and efficiency studies
 
    //std::cout<<"cut 6"<<std::endl;

      //std::cout<<"do we get anywhere to reco particle business"<<std::endl;


    bool fill_reco_particle = true;
    
    //if we neither have taus on reco and gen level don't fill information
    //for MC truth and reco extended entries, i.e. clear the corresponding vectors
    if(m_runTauMode==true && pass_tau_reco==false && pass_tau_gen==false){
      fill_reco_particle=false;
    }
    if(fill_reco_particle){

      //std::cout<<"do we get anywhere to reco particle business 2"<<std::endl;

      //std::cout<<"before reco collection"<<std::endl;

      LCCollection* recoparticlecol = NULL;
      // Alternativelly if you do not want Marlin to exit in case of a non-existing collection
      // use the following (commented out) code:
      //run on H to gamma gamma -> in case there are no tracks around linker will fail to produce output collection
      recoparticlecol = evt->getCollection(m_inputRECOParticleCollection) ;
      if(recoparticlecol!=NULL){
	//PandoraCandidate loop
	for(int i=0;i<recoparticlecol->getNumberOfElements();i++){
	  ReconstructedParticle* pandorapart = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(i));
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
	    m_reco_cluster0_energy->push_back(pandorapart->getClusters()[0]->getEnergy());
	    m_reco_cluster0_iPhi->push_back(pandorapart->getClusters()[0]->getIPhi());
	    m_reco_cluster0_iTheta->push_back(pandorapart->getClusters()[0]->getITheta());
	    m_reco_cluster0_energyError->push_back(pandorapart->getClusters()[0]->getEnergyError());
	  }else{
	    m_reco_clusters_energy->push_back(-1.);
	    m_reco_cluster0_energy->push_back(-1.);
	    m_reco_cluster0_iPhi->push_back(-10.);
	    m_reco_cluster0_iTheta->push_back(-10.);
	    m_reco_cluster0_energyError->push_back(-1.);
	  }
	  if(pandorapart->getTracks().size()>0 && pandorapart->getClusters().size()>0){
	    m_reco_EClustersOverPTrack->push_back(energy_sum_clusters/track_p);
	  }else{
	    m_reco_EClustersOverPTrack->push_back(-1.);
	  }
	  float cellEnergySum(0.f);
	  const EVENT::ClusterVec &clusterVec(pandorapart->getClusters());
	  for (EVENT::ClusterVec::const_iterator iter = clusterVec.begin(), iterEnd = clusterVec.end(); iter != iterEnd; ++iter)
	    {
	      const EVENT::CalorimeterHitVec &calorimeterHitVec((*iter)->getCalorimeterHits());
	      
	      for (EVENT::CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter)
		{
		  const EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;
		  cellEnergySum += pCalorimeterHit->getEnergy();
		}
	    }
	  
	  const float correctionFactor=((cellEnergySum < std::numeric_limits<float>::epsilon()) ? 0.f : pandorapart->getEnergy() / cellEnergySum);
	  //std::cout<<"correctionFactor is "<<correctionFactor<<std::endl;
	  
	  
	  //Check if in barrel
	  double cosTheta = pandorapart->getMomentum()[2]/sqrt(pandorapart->getMomentum()[0]*pandorapart->getMomentum()[0]+pandorapart->getMomentum()[1]*pandorapart->getMomentum()[1]+pandorapart->getMomentum()[2]*pandorapart->getMomentum()[2]);
	  m_reco_EOverP_PFA->push_back(pandorapart->getEnergy()/sqrt(pandorapart->getMomentum()[0]*pandorapart->getMomentum()[0]+pandorapart->getMomentum()[1]*pandorapart->getMomentum()[1]+pandorapart->getMomentum()[2]*pandorapart->getMomentum()[2]));
	  m_recoEnergy->push_back(pandorapart->getEnergy());
	  m_recohitEnergyCellCorrection->push_back(correctionFactor);
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
	  m_reco_Theta->push_back(acos(cosTheta));
	  float En_ECAL_Barrel=0;
	  float En_ECAL_Endcap=0;
	  float En_ECAL_else=0;
	  float En_HCAL_Barrel=0;
	  float En_HCAL_Endcap=0;
	  float En_HCAL_else=0;
	  float En_MU_Barrel=0;
	  float En_MU_Endcap=0;
	  float En_MU_else=0;
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
	  int hitsMB=0;
	  int hitsME=0;
	  int hitsMO=0;
	  double totEnergy=0;
	  double totLogEnergy=0;
	  double posX=0;
	  double posY=0;
	  double posZ=0;


	  //std::cout<<"before reco cluster collection"<<std::endl;

	  for(unsigned int j=0;j<pandorapart->getClusters().size();j++){
	    for(unsigned int l=0; l<pandorapart->getClusters()[j]->getCalorimeterHits().size();l++){
	      //first need to know the total hit energy sum to correct for log weights --> photon E corrected by PFO additional calibration factor
	      totEnergy+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	    }
	  }
	  //two clusters per particle doesn't seem to happen, but maybe in future
	  for(unsigned int j=0;j<pandorapart->getClusters().size();j++){
	    for(unsigned int l=0; l<pandorapart->getClusters()[j]->getCalorimeterHits().size();l++){
	      //loop over hits for logarithmic reweighting
	      // constant is changed by hand --> why where and how was that decided
	      //double logWeight = std::max ( 0.0 , 5.5 + std::log(pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy()/totEnergy));
	      double logWeight = std::max ( 0.0 , 6.5 + std::log(pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy()/totEnergy));
	      posX += pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[0]*logWeight;
	      posY += pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[1]*logWeight;
	      posZ += pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[2]*logWeight;
	      totLogEnergy += logWeight;
	      
	      const CHT cht=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getType();
	      
	      int types=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getType();
	      //static const int fCaloType =     1 ;
	      //static const int fCaloID   =    10 ;
	      //static const int fLayout   =  1000 ;
	      static const int fLayer    = 10000 ;
	      int caloLayer=types/fLayer;
	      //int caloLayout=(types-caloLayer*fLayer)/fLayout;
	      //int caloID=(types-caloLayer*fLayer-caloLayout*fLayout)/fCaloID;
	      //int caloType=(types-caloLayer*fLayer-caloLayout*fLayout-caloID*fCaloID)/fCaloType;		  
	      //if(caloID==1){//ecal
	      //if(cht.is(CHT::muon)){
	      //std::cout<<"ec/had type "<<cht.is(CHT::em)<<"/"<<cht.is(CHT::had)<<" hcal/ecal/yoke "<<cht.is(CHT::hcal)<<"/"<<cht.is(CHT::ecal)<<"/"<<cht.is(CHT::yoke)<<" barrel/endcap "<<cht.is(CHT::barrel)<<"/"<<cht.is(CHT::endcap)<<std::endl;
	      //}
	      if(cht.is(CHT::ecal)){//ecal
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
		  //if(caloLayout==1){//id 1=ecal, 2 hcal, layout 1 barrel, 2 endcap
		  if(cht.is(CHT::barrel)){
		    En_ECAL_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		    hitsEB+=1;
		    //}else if(caloLayout==2){
		  }else if (cht.is(CHT::endcap)){
		    En_ECAL_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		    hitsEE+=1;
		  }else{
		    En_ECAL_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		    hitsEO+=1;
		  }
		  //}else if(caloID==2){//hcal
		}else if(cht.is(CHT::hcal)){//h-cal
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
		  //if(caloLayout==1){//id 1=ecal, 2 hcal, layout 1 barrel, 2 endcap
		  if(cht.is(CHT::barrel)){
		    En_HCAL_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		    hitsHB+=1;
		    //}else if(caloLayout==2){
		  }else if(cht.is(CHT::endcap)){
		    En_HCAL_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		    hitsHE+=1;
		  }else{
		    En_HCAL_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		    hitsHO+=1;
		  }
		}else if (cht.is(CHT::yoke)){
		  if(cht.is(CHT::barrel)){
		    En_MU_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		    hitsMB+=1;
		  }else if(cht.is(CHT::endcap)){
		    En_MU_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		    hitsME+=1;
		  }else{
		    En_MU_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		    hitsMO+=1;
		  }
		}else{
		  std::cout<<"do we really have even other calohits not ecal/hcal/muon "<<std::endl;
		}
	      }//hit loop
	    }//cluster loop
	    m_reco_E_totRW->push_back(totEnergy);
	    posX=posX/totLogEnergy;
	    posY=posY/totLogEnergy;
	    posZ=posZ/totLogEnergy;
	    m_reco_x_logERW->push_back(posX);
	    m_reco_y_logERW->push_back(posY);
	    m_reco_z_logERW->push_back(posZ);
	    m_reco_logE_tot->push_back(totLogEnergy);
	    m_reco_Phi_logERW->push_back(atan2(posY,posX));
	    double cosTheta_logERW=posZ/sqrt(posX*posX+posY*posY+posZ*posZ);
	    m_reco_CosTheta_logERW->push_back(cosTheta_logERW);
	    m_reco_Theta_logERW->push_back(acos(cosTheta_logERW));
	    m_reco_E_EB->push_back(En_ECAL_Barrel);
	    m_reco_E_EE->push_back(En_ECAL_Endcap);
	    m_reco_E_EO->push_back(En_ECAL_else);
	    m_reco_E_HB->push_back(En_HCAL_Barrel);
	    m_reco_E_HE->push_back(En_HCAL_Endcap);
	    m_reco_E_HO->push_back(En_HCAL_else);
	    m_reco_E_MB->push_back(En_MU_Barrel);
	    m_reco_E_ME->push_back(En_MU_Endcap);
	    m_reco_E_MO->push_back(En_MU_else);
	    m_reco_firstLayerECAL->push_back(phFirstLayerECAL);
	    m_reco_lastLayerECAL->push_back(phLastLayerECAL);
	    m_reco_nhitsEB->push_back(hitsEB);
	    m_reco_nhitsEE->push_back(hitsEE);
	    m_reco_nhitsEO->push_back(hitsEO);
	    m_reco_firstLayerHCAL->push_back(phFirstLayerHCAL);
	    m_reco_lastLayerHCAL->push_back(phLastLayerHCAL);
	    m_reco_nhitsHB->push_back(hitsHB);
	    m_reco_nhitsHE->push_back(hitsHE);
	    m_reco_nhitsHO->push_back(hitsHO);
	    m_reco_nhitsMB->push_back(hitsMB);
	    m_reco_nhitsME->push_back(hitsME);
	    m_reco_nhitsMO->push_back(hitsMO);
	    //pf particle loop done
	  }
	}//reco collection available
      }//skip reco calculation if no tau around for tau running mode
      
      //if we neither have taus on reco and gen level don't fill information
      //for MC truth and reco extended entries, i.e. clear the corresponding vectors
      if(m_runTauMode==true && pass_tau_reco==false && pass_tau_gen==false){
	m_trueEnergy->clear();
	m_true_Px->clear();
	m_true_Py->clear();
	m_true_Pz->clear();
	m_true_x->clear();
	m_true_y->clear();
	m_true_z->clear();
	m_true_CosTheta->clear();
	m_true_Theta->clear();
	m_true_Phi->clear();
	m_true_PDGID->clear();
	m_true_GenStatus->clear();
	m_true_numDaughters->clear();
	m_true_decayTrackerCalo->clear();
	m_true_motherDecayTrackerCalo->clear();
	m_true_numMothers->clear();
	m_true_m1_PDGID->clear();
	m_true_m2_PDGID->clear();
	m_true_m1_status->clear();
	m_true_m2_status->clear();
	m_true_m1_E->clear();
	m_true_m2_E->clear();

	m_recoEnergy->clear();
	m_recohitEnergyCellCorrection->clear();
	m_reco_Px->clear();
	m_reco_Py->clear();
	m_reco_Pz->clear();
	m_reco_CosTheta->clear();
	m_reco_Theta->clear();
	m_reco_Phi->clear();
	m_reco_Charge->clear();
	m_reco_PDGID->clear();
	m_reco_CosTheta_logERW->clear();
	m_reco_Theta_logERW->clear();
	m_reco_Phi_logERW->clear();
	m_reco_x_logERW->clear();
	m_reco_y_logERW->clear();
	m_reco_z_logERW->clear();
	m_reco_logE_tot->clear();
	m_reco_E_totRW->clear();
	
	m_reco_nTracks->clear();
	m_reco_track0_pt->clear();
	m_reco_track0_p->clear();
	m_reco_track0_nHits->clear();
	m_reco_track0_chi2OverNdof->clear();
	m_reco_nClusters->clear();
	m_reco_clusters_energy->clear();
	m_reco_cluster0_energy->clear();
	m_reco_cluster0_iPhi->clear();
	m_reco_cluster0_iTheta->clear();
	m_reco_cluster0_energyError->clear();
	m_reco_EClustersOverPTrack->clear();
	m_reco_EOverP_PFA->clear();
	m_reco_E_EB->clear();
	m_reco_E_EE->clear();
	m_reco_E_EO->clear();
	m_reco_E_HB->clear();
	m_reco_E_HE->clear();
	m_reco_E_HO->clear();
	m_reco_E_MB->clear();
	m_reco_E_ME->clear();
	m_reco_E_MO->clear();
	m_reco_firstLayerECAL->clear();
	m_reco_lastLayerECAL->clear();
	m_reco_nhitsEB->clear();
	m_reco_nhitsEE->clear();
	m_reco_nhitsEO->clear();
	m_reco_firstLayerHCAL->clear();
	m_reco_lastLayerHCAL->clear();
	m_reco_nhitsHB->clear();
	m_reco_nhitsHE->clear();
	m_reco_nhitsHO->clear();
	m_reco_nhitsMB->clear();
	m_reco_nhitsME->clear();
	m_reco_nhitsMO->clear();
      }

      if(m_fillClusterHits){//checked for cluster hit filling flag runs
	bool fill_muonhits=false;
	/*
	if(true_in_transition || !m_no_cluster_hit_x->empty()){
	fill_muonhits=true;
	}
	if(!fill_muonhits){//if muon hits are from non interesting events
	  m_cluster_hit_x->clear();
	  m_cluster_hit_y->clear();
	  m_cluster_hit_z->clear();
	  m_cluster_hit_E->clear();
	  m_cluster_hit_index->clear();
	  m_cluster_hit_type->clear();
	  
	  m_no_cluster_hit_x->clear();
	  m_no_cluster_hit_y->clear();
	  m_no_cluster_hit_z->clear();
	  m_no_cluster_hit_E->clear();
	}
	*/
      }


    //std::cout<<"cut 7"<<std::endl;
 
      //if(true_in_transition){
	m_outputTree->Fill();
	//}
}

void PionStudy::fillStableDaughterSet(MCParticle* mcp, std::set<MCParticle*> &stableDaughterSet){
  if(mcp->getGeneratorStatus()==1){
    stableDaughterSet.insert(mcp);
  }else if (mcp->getGeneratorStatus()==0){
    return;
  }
  for(unsigned int d=0;d<mcp->getDaughters().size();d++){
    fillStableDaughterSet(mcp->getDaughters()[d], stableDaughterSet);
  }
}



void PionStudy::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    
    streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

void PionStudy::check(LCEvent*){
}

void PionStudy::end(){
    

  
 
  m_rootFile->Write();
  m_rootFile->Close();

}
