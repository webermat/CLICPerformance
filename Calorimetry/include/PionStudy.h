#ifndef PionStudy_h
#define PionStudy_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>

#include "TH2F.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TFile.h"
#include <map>


using namespace lcio ;
using namespace marlin ;


class PionStudy : public Processor {
                
public:
        
    virtual Processor*  newProcessor() { return new PionStudy ; }
    
    PionStudy() ;

    PionStudy(const PionStudy&) = delete;
    PionStudy& operator=(const PionStudy&) = delete;
    
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
        
protected:

    // Collection names for (in/out)put
    std::string m_inputCalorimeterHitCollection="";
    std::string m_inputMCParticleCollection="";
    std::string m_inputRECOParticleCollection="";
    std::string m_inputTrackCollection="";
    std::string m_inputClusterCollection="";
    std::string m_inputTauCollection="";
    std::string m_inputSimCaloHitCollectionName="";
    std::string m_inputMuonHitCollectionName="";
    bool m_runTauMode=false;
    bool m_runPhotonMode=false;//safe conversion information
    bool m_fillClusterHits=false;
    bool m_fillSimHits=false;
    bool m_fillClusterBranches=false;
    bool m_fillTrackBranches=false;

    std::string m_rootFileName="";
    
    // Run and event counters
    int m_eventNumber=0;
    int m_runNumber=0;

    float  m_innerBField=0;
    float  m_const_a=0;

    float  m_trackerOuterR=0;
    float  m_trackerInnerR=0;
    float  m_trackerZmax=0;

    int m_barrelTrackerLayers=0;

    int eventcount =0;
    
    // Call to get collections
    void getCollection(LCCollection*&, std::string, LCEvent*);

    void fillStableDaughterSet(MCParticle*, std::set<MCParticle*> &);
    
    TTree * m_outputTree=NULL;
    TFile * m_rootFile=NULL;
    //Tree branch variables
    
    std::vector<float> * m_trueEnergy=NULL;
    std::vector<float> * m_true_Px=NULL;
    std::vector<float> * m_true_Py=NULL;
    std::vector<float> * m_true_Pz=NULL;
    std::vector<float> * m_true_x=NULL;//vertex x of particle
    std::vector<float> * m_true_y=NULL;//vertex y of particle
    std::vector<float> * m_true_z=NULL;//vertex z of particle
    std::vector<float> * m_true_CosTheta=NULL;
    std::vector<float> * m_true_Theta=NULL;
    std::vector<float> * m_true_Phi=NULL;
    std::vector<int> * m_true_PDGID=NULL;
    std::vector<int> * m_true_GenStatus=NULL;
    std::vector<int> * m_true_numDaughters=NULL;
    std::vector<int> *m_true_decayTrackerCalo=NULL;
    std::vector<int> *m_true_motherDecayTrackerCalo=NULL;
    std::vector<int> * m_true_numMothers=NULL;
    std::vector<int> * m_true_m1_PDGID=NULL;
    std::vector<int> * m_true_m2_PDGID=NULL;
    std::vector<int> * m_true_m1_status=NULL;
    std::vector<int> * m_true_m2_status=NULL;
    std::vector<int> * m_true_index=NULL;
    std::vector<float> * m_true_m1_E=NULL;
    std::vector<float> * m_true_m2_E=NULL;

    std::vector<float>* m_true_conv_e1_Px = NULL;
    std::vector<float>* m_true_conv_e1_Py = NULL;
    std::vector<float>* m_true_conv_e1_Pz = NULL;
    std::vector<float>* m_true_conv_e1_E = NULL;
    std::vector<int>* m_true_conv_e1_PDGID = NULL;
    std::vector<float>* m_true_conv_e2_Px = NULL;
    std::vector<float>* m_true_conv_e2_Py = NULL;
    std::vector<float>* m_true_conv_e2_Pz = NULL;
    std::vector<float>* m_true_conv_e2_E = NULL;
    std::vector<int>* m_true_conv_e2_PDGID = NULL;
    std::vector<int>* m_true_conv_Ph_VecInd = NULL;//index of converted photon
    std::vector<float>* m_true_conv_Vtx_x = NULL;
    std::vector<float>* m_true_conv_Vtx_y = NULL;
    std::vector<float>* m_true_conv_Vtx_z = NULL;

    std::vector<float>* m_true_tauDaughter_Energy = NULL;
    std::vector<float>* m_true_tauDaughter_Px = NULL;
    std::vector<float>* m_true_tauDaughter_Py = NULL;
    std::vector<float>* m_true_tauDaughter_Pz = NULL;
    std::vector<int>* m_true_tauDaughter_PDGID = NULL;
    std::vector<int>* m_true_tauDaughter_Charge = NULL;
    std::vector<int>* m_true_tauDaughter_tauIndex = NULL;
    std::vector<int>* m_true_tauDaughter_status = NULL;
    std::vector<int>* m_true_tauDaughter_motherPDGID = NULL;
    std::vector<float>* m_true_tauDaughter_motherEnergy = NULL;

    int m_true_stable_firstDaughter_PDGID=0;
    float m_true_stable_firstDaughter_E=0;
    float m_true_stable_firstDaughterVtx_x=0;
    float m_true_stable_firstDaughterVtx_y=0;
    float m_true_stable_firstDaughterVtx_z=0;

    float m_true_stable_10pc_E=0;//last daughter to be created, when over 10 % of energy is lost
    float m_true_stable_10pc_Esum=0;//sum of energy taken away by daughters beyond the 10 % energy loss point
    float m_true_stable_10pcVtx_x=0;//to determine true vertex, at which point the original particle loses more than 10 % of its energy
    float m_true_stable_10pcVtx_y=0;
    float m_true_stable_10pcVtx_z=0;

    float m_true_stable_Endpoint_x=0;
    float m_true_stable_Endpoint_y=0;
    float m_true_stable_Endpoint_z=0;

    float m_ECAL_endcapZ_min=0;
    float m_ECAL_endcapR_min=0;
    float m_ECAL_barrelR_min=0;


    float m_ECAL_barrelZ_max=0;
    float m_ECAL_ringZ_min=0;
    float m_ECAL_ringZ_max=0;
    //float m_ECAL_endcapZ_min=0;
    float m_ECAL_endcapZ_max=0;

    float m_HCAL_barrelZ_max=0;
    float m_HCAL_ringZ_min=0;
    float m_HCAL_ringZ_max=0;
    float m_HCAL_endcapZ_min=0;
    float m_HCAL_endcapZ_max=0;

    float m_MUON_barrelZ_max=0;
    float m_MUON_endcapZ_min=0;
    float m_MUON_endcapZ_max=0;

    float m_Z_mcE=0;
    int m_Z_mcNDaughter=0;

    int m_d1_mcPDGID=0;
    float m_d1_mcE=0;
    float m_d1_mcPx=0;
    float m_d1_mcPy=0;
    float m_d1_mcPz=0;
    float m_d1_mcPhi=0;
    float m_d1_mcTheta=0;
    float m_d1_mcCosTheta=0;
    float m_d1_mcMass=0;

    int m_d2_mcPDGID=0;
    float m_d2_mcE=0;
    float m_d2_mcPx=0;
    float m_d2_mcPy=0;
    float m_d2_mcPz=0;
    float m_d2_mcPhi=0;
    float m_d2_mcTheta=0;
    float m_d2_mcCosTheta=0;
    float m_d2_mcMass=0;

    std::vector<float>* m_tauJet_Px=NULL; 
    std::vector<float>* m_tauJet_Py=NULL; 
    std::vector<float>* m_tauJet_Pz=NULL; 
    std::vector<float>* m_tauJet_E=NULL; 
    std::vector<float>* m_tauJet_Phi=NULL; 
    std::vector<float>* m_tauJet_CosTheta=NULL; 
    std::vector<int>* m_tauJet_charge=NULL; 
    std::vector<int>* m_tauJet_neutMult=NULL;
    std::vector<int>* m_tauJet_chMult=NULL;

    std::vector<float>* m_tauJet_Part_Px=NULL; 
    std::vector<float>* m_tauJet_Part_Py=NULL; 
    std::vector<float>* m_tauJet_Part_Pz=NULL; 
    std::vector<float>* m_tauJet_Part_E=NULL; 
    std::vector<int>* m_tauJet_Part_charge=NULL;
    std::vector<int>* m_tauJet_Part_PDGID=NULL;
    std::vector<int>* m_tauJet_Part_JetIndex=NULL;

    std::vector<float> *m_cluster_x=NULL;
    std::vector<float> *m_cluster_y=NULL;
    std::vector<float> *m_cluster_z=NULL;
    std::vector<float> *m_cluster_dir_x=NULL;
    std::vector<float> *m_cluster_dir_y=NULL;
    std::vector<float> *m_cluster_dir_z=NULL;
    std::vector<float> *m_cluster_track_PosAngle_min=NULL;//min distance between cluster and any track at CaloFace
    std::vector<float> *m_cluster_track_DirAngle=NULL;//angle between cluster and min distance track track
    //this time use the PandoraPFA implementation, in LCContent: src/LCHelpers/ClusterHelper.cc
    std::vector<float> *m_cluster_track_minPFADist_cluster_EOverP=NULL;//EOverP for min PFA like distance cluster -> E cluster, P track
    std::vector<float> *m_cluster_track_minPFADist_cluster_Angle=NULL;//for min PFA like distance cluster direction and track direction angle
    std::vector<float> *m_cluster_track_minPFADist_cluster_TCDistance=NULL;//TrackClusterDistance for min PFA like distance cluster

    std::vector<float> *m_cluster_energy=NULL;
    std::vector<float> *m_cluster_energyError=NULL;
    std::vector<float> *m_cluster_ECAL_energy=NULL;
    std::vector<float> *m_cluster_HCAL_energy=NULL;
     std::vector<float> *m_cluster_E_EB=NULL;
    std::vector<float> *m_cluster_E_EE=NULL;
    std::vector<float> *m_cluster_E_EO=NULL;
    std::vector<float> *m_cluster_E_HB=NULL;
    std::vector<float> *m_cluster_E_HE=NULL;
    std::vector<float> *m_cluster_E_HO=NULL;
    std::vector<float> *m_cluster_E_MU=NULL;
    std::vector<float> *m_cluster_firstShowerLayer=NULL;//first layer which has neighbours with hits, assume try to avoid random first noise hits 
    std::vector<int> *m_cluster_nHits_ECAL=NULL;
    std::vector<int> *m_cluster_nHits_HCAL=NULL;
    std::vector<int> *m_cluster_nHits_MUON=NULL;
    std::vector<float> *m_cluster_iPhi=NULL;
    std::vector<float> *m_cluster_iTheta=NULL;

    std::vector<float> *m_cluster_hit_x=NULL;
    std::vector<float> *m_cluster_hit_y=NULL;
    std::vector<float> *m_cluster_hit_z=NULL;
    std::vector<float> *m_cluster_hit_E=NULL;
    std::vector<int> *m_cluster_hit_index=NULL;
    std::vector<int> *m_cluster_hit_type=NULL;//1 for ECAL, 2 for HCAL, 3 for muon hits

    std::vector<float> *m_no_cluster_hit_x=NULL;
    std::vector<float> *m_no_cluster_hit_y=NULL;
    std::vector<float> *m_no_cluster_hit_z=NULL;
    std::vector<float> *m_no_cluster_hit_E=NULL;
    //std::vector<int> *m_no_cluster_hit_index=NULL;
    //std::vector<int> *m_no_cluster_hit_type=NULL;//1 for ECAL, 2 for HCAL, 3 for muon hits

    std::vector<float> *m_muon_b_simhit_x=NULL;
    std::vector<float> *m_muon_b_simhit_y=NULL;
    std::vector<float> *m_muon_b_simhit_z=NULL;
    std::vector<float> *m_muon_b_simhit_E=NULL;

    //normalized position values
    std::vector<float> *m_track_d0=NULL;
    std::vector<float> *m_track_z0=NULL;
    std::vector<float> *m_track_phi0=NULL;
    std::vector<int> *m_track_ndf=NULL;
    std::vector<int> *m_track_nHits=NULL;
    std::vector<float> *m_track_nExpectedTrackerHits=NULL;
    std::vector<int> *m_track_nTrackerBarrelHits=NULL;//tracker+vertex
    std::vector<int> *m_track_nTrackerEndcapHits=NULL;//tracker+vertex
    std::vector<int> *m_track_nVertexBarrelHits=NULL;
    std::vector<int> *m_track_nVertexEndcapHits=NULL;
    std::vector<float> *m_track_chi2=NULL;
    std::vector<float> *m_track_sigmaPOverP=NULL;
    //define innermost as one with lowest R
    std::vector<float> *m_track_x_innermostHit=NULL;
    std::vector<float> *m_track_y_innermostHit=NULL;
    std::vector<float> *m_track_z_innermostHit=NULL;
    std::vector<float> *m_track_pt_innermostHit=NULL;//inner R track
    std::vector<float> *m_track_p_innermostHit=NULL;//inner R track
    std::vector<float> *m_track_zMin=NULL;
    std::vector<float> *m_track_x_outermostRHit=NULL;
    std::vector<float> *m_track_y_outermostRHit=NULL;
    std::vector<float> *m_track_z_outermostRHit=NULL;
    std::vector<float> *m_track_pt_outermostRHit=NULL;
    std::vector<float> *m_track_p_outermostRHit=NULL;
    std::vector<float> *m_track_x_outermostZHit=NULL;
    std::vector<float> *m_track_y_outermostZHit=NULL;
    std::vector<float> *m_track_z_outermostZHit=NULL;
    std::vector<float> *m_track_x_atIP=NULL;
    std::vector<float> *m_track_y_atIP=NULL;
    std::vector<float> *m_track_z_atIP=NULL;
    std::vector<float> *m_track_Phi_atIP=NULL;
    std::vector<float> *m_track_Theta_atIP=NULL;
    std::vector<float> *m_track_pt_atIP=NULL;
    std::vector<float> *m_track_p_atIP=NULL;
    std::vector<float> *m_track_x_atCalo=NULL;
    std::vector<float> *m_track_y_atCalo=NULL;
    std::vector<float> *m_track_z_atCalo=NULL;
    std::vector<float> *m_track_px_atCalo=NULL;
    std::vector<float> *m_track_py_atCalo=NULL;
    std::vector<float> *m_track_pz_atCalo=NULL;
    std::vector<float> *m_track_Phi_atCalo=NULL;
    std::vector<float> *m_track_Theta_atCalo=NULL;
    std::vector<float> *m_track_pt_atCalo=NULL;
    std::vector<float> *m_track_p_atCalo=NULL;
    std::vector<float> *m_track_cluster_PosAngle_min_atCalo=NULL;
    std::vector<float> *m_track_minDist_cluster_EOverP=NULL;//EOverP between track and min distance at CaloFace cluster
    std::vector<float> *m_track_cluster_minPFADist_EOverP=NULL;//EOverP for min PFA like distance cluster -> E cluster, P track
    std::vector<float> *m_track_cluster_minPFADist_clusterEMEnergy=NULL;//EMEnergy for min PFA like distance cluster
    std::vector<float> *m_track_cluster_minPFADist_clusterHadEnergy=NULL;//EMEnergy for min PFA like distance cluster
    std::vector<float> *m_track_cluster_minPFAdistance_atCalo=NULL;//min PFA like distance between an cluster and track at CaloFace
    std::vector<float> *m_track_cluster_min_parPFAdistance_atCalo=NULL;//min PFA like parallel distance between an cluster and track at CaloFace (used in preselection of track-cluster association
    std::vector<float> *m_track_cluster_DirAngle_minPFAdistance_atCalo=NULL;//diraction angle between cluster and track at CaloFace, for cluster with minPFADistance

    std::vector<float> *m_track_cluster_minPFADist_EOverP_tinyHad=NULL;//EOverP for min PFA like distance cluster -> E cluster, P track
    std::vector<float> *m_track_cluster_minPFADist_clusterEMEnergy_tinyHad=NULL;//EMEnergy for min PFA like distance cluster
    std::vector<float> *m_track_cluster_minPFADist_clusterHadEnergy_tinyHad=NULL;//EMEnergy for min PFA like distance cluster
    std::vector<float> *m_track_cluster_minPFAdistance_atCalo_tinyHad=NULL;//min PFA like distance between an cluster and track at CaloFace
    std::vector<float> *m_track_cluster_min_parPFAdistance_atCalo_tinyHad=NULL;//min PFA like parallel distance between an cluster and track at CaloFace (used in preselection of track-cluster association
    std::vector<float> *m_track_cluster_DirAngle_minPFAdistance_atCalo_tinyHad=NULL;//diraction angle between cluster and track at CaloFace, for cluster with minPFADistance


    std::vector<float> *m_reco_CosTheta_logERW=NULL;
    std::vector<float> *m_reco_Theta_logERW=NULL;
    std::vector<float> *m_reco_Phi_logERW=NULL;
    std::vector<float> *m_reco_logE_tot=NULL;
    std::vector<float> *m_reco_E_totRW=NULL;
    //normalized position values
    std::vector<float> *m_reco_x_logERW=NULL;
    std::vector<float> *m_reco_y_logERW=NULL;
    std::vector<float> *m_reco_z_logERW=NULL;

    std::vector<float> *m_recoEnergy=NULL;
    std::vector<float> *m_recohitEnergyCellCorrection=NULL;
    std::vector<float> *m_reco_Px=NULL;
    std::vector<float> *m_reco_Py=NULL;
    std::vector<float> *m_reco_Pz=NULL;
    std::vector<float> *m_reco_CosTheta=NULL;
    std::vector<float> *m_reco_Theta=NULL;
    std::vector<float> *m_reco_Phi=NULL;
    std::vector<int> *m_reco_Charge=NULL;

    std::vector<int> *m_reco_nTracks=NULL;
    std::vector<float> *m_reco_track0_pt=NULL;//calculated from defaults/IP
    std::vector<float> *m_reco_track0_p=NULL;//calculated from defaults/IP
    std::vector<int> *m_reco_track0_nHits=NULL;
    std::vector<float> *m_reco_track0_chi2OverNdof=NULL;
    std::vector<int> *m_reco_PDGID=NULL;
    std::vector<int> *m_reco_nClusters=NULL;
    std::vector<float> *m_reco_clusters_energy=NULL;//sum of all cluster energies
    std::vector<float> *m_reco_cluster0_energy=NULL;//sum of all cluster energies
    std::vector<float> *m_reco_cluster0_iPhi=NULL;
    std::vector<float> *m_reco_cluster0_iTheta=NULL;
    std::vector<float> *m_reco_cluster0_energyError=NULL;
    std::vector<float> *m_reco_EClustersOverPTrack=NULL;
    std::vector<float> *m_reco_EOverP_PFA=NULL;
    std::vector<float> *m_reco_E_EB=NULL;
    std::vector<float> *m_reco_E_EE=NULL;
    std::vector<float> *m_reco_E_EO=NULL;
    std::vector<float> *m_reco_E_HB=NULL;
    std::vector<float> *m_reco_E_HE=NULL;
    std::vector<float> *m_reco_E_HO=NULL;
    std::vector<float> *m_reco_E_MB=NULL;
    std::vector<float> *m_reco_E_ME=NULL;
    std::vector<float> *m_reco_E_MO=NULL;
    std::vector<int> *m_reco_firstLayerECAL=NULL;
    std::vector<int> *m_reco_firstConsecutiveLayerECAL=NULL;//layer in ECAL which has neighbouring layer filled
    std::vector<int> *m_reco_lastLayerECAL=NULL;
    std::vector<int> *m_reco_nhitsEB=NULL;
    std::vector<int> *m_reco_nhitsEE=NULL;
    std::vector<int> *m_reco_nhitsEO=NULL;
    std::vector<int> *m_reco_firstLayerHCAL=NULL;
    std::vector<int> *m_reco_lastLayerHCAL=NULL;
    std::vector<int> *m_reco_nhitsHB=NULL;
    std::vector<int> *m_reco_nhitsHE=NULL;
    std::vector<int> *m_reco_nhitsHO=NULL;
    std::vector<int> *m_reco_nhitsMB=NULL;
    std::vector<int> *m_reco_nhitsME=NULL;
    std::vector<int> *m_reco_nhitsMO=NULL;

} ;

#endif



