#ifndef TTBarStudy_h
#define TTBarStudy_h 1

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

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetectorSelector.h"


using namespace lcio ;
using namespace marlin ;


class TTBarStudy : public Processor {
                
public:
        
    virtual Processor*  newProcessor() { return new TTBarStudy ; }
    
    TTBarStudy() ;

    TTBarStudy(const TTBarStudy&) = delete;
    TTBarStudy& operator=(const TTBarStudy&) = delete;
    
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
    std::string m_inputMCParticleCollection="";
    std::string m_inputRECOParticleCollection="";
    std::string m_inputTauCollection="";
    std::string m_rootFileName="";
    bool m_runTauMode=false;
    bool m_runGenOnly=false;
    
    // Run and event counters
    int m_eventNumber=0;
    int m_runNumber=0;

    float  m_innerBField=0;
    float  m_const_a=0;

    int eventcount =0;
    
    // Call to get collections
    void getCollection(LCCollection*&, std::string, LCEvent*);

    void fillStableDaughterSet(MCParticle*, std::set<MCParticle*> &);
    
    TTree * m_outputTree=NULL;
    TFile * m_rootFile=NULL;
    //Tree branch variables
    
    //particle level information sums --> visible and invisible energy
    //only neutrinos
    float m_E_trueInv=0;
    float m_px_trueInv=0;
    float m_py_trueInv=0;
    float m_pz_trueInv=0;

    //all visible particles
    float m_E_trueAll=0;
    float m_px_trueAll=0;
    float m_py_trueAll=0;
    float m_pz_trueAll=0;

    //all reconstructed particles
    float m_E_totPFO=0;
    float m_px_totPFO=0;
    float m_py_totPFO=0;
    float m_pz_totPFO=0;

    float m_E_top1  = 0;
    float m_px_top1 = 0;
    float m_py_top1 = 0;
    float m_pz_top1 = 0;
    int m_PDGID_top1 = 0;

    float m_E_top2  = 0;
    float m_px_top2 = 0;
    float m_py_top2 = 0;
    float m_pz_top2 = 0;
    int m_PDGID_top2 = 0;

    //check if we are also including the bbar running
    std::vector<float> * m_trueEnergy=NULL;
    std::vector<float> * m_true_Px=NULL;
    std::vector<float> * m_true_Py=NULL;
    std::vector<float> * m_true_Pz=NULL;
    std::vector<float> * m_true_Vtxx=NULL;
    std::vector<float> * m_true_Vtxy=NULL;
    std::vector<float> * m_true_Vtxz=NULL;
    std::vector<float> * m_true_Vtxr=NULL;
    std::vector<float> * m_true_Eptx=NULL;
    std::vector<float> * m_true_Epty=NULL;
    std::vector<float> * m_true_Eptz=NULL;
    std::vector<float> * m_true_Eptr=NULL;
    std::vector<float> * m_true_CosTheta=NULL;
    std::vector<float> * m_true_Phi=NULL;
    std::vector<int> * m_true_PDGID=NULL;//only save muons, electrons, taus
    std::vector<int> * m_true_GenStatus=NULL;//for taus safe ONLY the ones from W, status then not stable
    std::vector<int> * m_true_numDaughters=NULL;
    std::vector<int> *m_true_decayTrackerCalo=NULL;    //-1 if decayed in generator (e.g. tau's), 1 for tracker, 2 for calorimeter
    std::vector<int> * m_true_numMothers=NULL;//check for first mothers which are NOT muons or electrons, for taus from W's we should know from the previous ones
    std::vector<int> * m_true_m1_PDGID=NULL;
    std::vector<int> * m_true_m1_index=NULL;
    std::vector<int> * m_true_index=NULL;//index in truth vector
    std::vector<float> * m_true_m1_E=NULL;
    //std::vector<float> * m_true_m2_E=NULL;

    //safe true particle energies of particles in angle of 6 around lepton candidate --> DON'T CHECK NEUTRINOS
    std::vector<float> * m_true_all_6_deg_E=NULL;
    //safe the actual angular difference
    std::vector<float> * m_true_all_6_deg_angle=NULL;
    std::vector<int> * m_true_all_6_deg_index=NULL;
    std::vector<int> * m_true_all_6_deg_PDG=NULL;//safe index of true tauDaughter index

    //only for taus from W's
    std::vector<float>* m_true_tauDaughter_Energy = NULL;
    std::vector<float>* m_true_tauDaughter_Px = NULL;
    std::vector<float>* m_true_tauDaughter_Py = NULL;
    std::vector<float>* m_true_tauDaughter_Pz = NULL;
    std::vector<float> * m_true_tauDaughter_Vtxx=NULL;
    std::vector<float> * m_true_tauDaughter_Vtxy=NULL;
    std::vector<float> * m_true_tauDaughter_Vtxz=NULL;
    std::vector<float> * m_true_tauDaughter_Vtxr=NULL;
    std::vector<int>* m_true_tauDaughter_PDGID = NULL;
    std::vector<int>* m_true_tauDaughter_Charge = NULL;
    std::vector<int>* m_true_tauDaughter_tauIndex = NULL;
    std::vector<int>* m_true_tauDaughter_status = NULL;
    std::vector<int>* m_true_tauDaughter_motherPDGID = NULL;
    std::vector<float>* m_true_tauDaughter_motherEnergy = NULL;
    //safe true particle energies of particles in angle of 6 around lepton candidate --> DON'T CHECK NEUTRINOS
    std::vector<float> * m_true_tauDaughter_all_6_deg_E=NULL;
    //safe the actual angular difference
    std::vector<float> * m_true_tauDaughter_all_6_deg_angle=NULL;
    std::vector<int> * m_true_tauDaughter_all_6_deg_index=NULL;//safe index of true tauDaughter index
    std::vector<int> * m_true_tauDaughter_all_6_deg_PDG=NULL;//safe index of true tauDaughter index



    int m_ttbar_decay_mode=-1;//0 for all hadronic, 1 for semileptonic decays with 1 muon and electron, 2 for dileptonic decays with two (electron, muon) 3 semileptonic, but tau into muon,electron, 4 semileptonic but tau into hadrons, 5 dileptonic but one tau lepton decay + electron/muon, 6 dilepton, but one hadronic tau decay, 7 two leptonic tau decays, 8 two taus, one into leptons, one into hadrons, 9 two taus, both hadronic

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

    //here we safe ONLY electrons and muons
    std::vector<float> *m_recoEnergy=NULL;
    std::vector<float> *m_reco_Px=NULL;
    std::vector<float> *m_reco_Py=NULL;
    std::vector<float> *m_reco_Pz=NULL;
    std::vector<float> *m_reco_CosTheta=NULL;
    std::vector<float> *m_reco_Phi=NULL;
    std::vector<int> *m_reco_Charge=NULL;
    std::vector<int> *m_reco_PDGID=NULL;

    std::vector<int> *m_reco_nTracks=NULL;
    std::vector<float> *m_reco_track0_pt=NULL;//calculated from defaults/IP
    std::vector<float> *m_reco_track0_p=NULL;//calculated from defaults/IP
    std::vector<int> *m_reco_track0_nHits=NULL;
    std::vector<float> *m_reco_track0_chi2OverNdof=NULL;
    std::vector<int> *m_reco_nClusters=NULL;
    std::vector<float> *m_reco_clusters_energy=NULL;//sum of all cluster energies
 
    std::vector<float> *m_reco_E_EB=NULL;
    std::vector<float> *m_reco_E_EE=NULL;
    std::vector<float> *m_reco_E_EO=NULL;
    std::vector<float> *m_reco_E_HB=NULL;
    std::vector<float> *m_reco_E_HE=NULL;
    std::vector<float> *m_reco_E_HO=NULL;
    //for NON leptons (electron, muons)
    //safe energies of particles around in cone of 6 degrees
    std::vector<int> *m_reco_all_6_deg_index = NULL;
    std::vector<float> *m_reco_all_6_deg_E = NULL;
    std::vector<float> *m_reco_all_6_deg_angle = NULL;
    std::vector<int> *m_reco_all_6_deg_type = NULL; //0 all MCParticles, 1 all loose/2 all selected, 3 all tight
    std::vector<int> *m_reco_all_6_deg_PDG = NULL; //0 all MCParticles, 1 all loose/2 all selected, 3 all tight   

} ;

#endif



