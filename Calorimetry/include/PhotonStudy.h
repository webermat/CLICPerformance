#ifndef PhotonStudy_h
#define PhotonStudy_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>


#include "TH2F.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <map>

using namespace lcio ;
using namespace marlin ;

class TTree;


class PhotonStudy : public Processor {
                
public:
        
    virtual Processor*  newProcessor() { return new PhotonStudy ; }
    
    PhotonStudy() ;

    PhotonStudy(const PhotonStudy&) = delete;
    PhotonStudy& operator=(const PhotonStudy&) = delete;
    
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
        
    std::string m_inputMCParticleCollection="";
    std::string m_inputRECOParticleCollection="";
    std::string m_jetInclColName="";
    std::string m_jetExclColName="";

    std::string m_rootFileName="";
    
    // Run and event counters
    int m_eventNumber=0;
    int m_runNumber=0;
    int m_eventCount=0;

    bool m_ignoreGen=false;
    bool m_fillMEInfo=false;
    bool m_reducedOutput=false;

    int eventcount=0;
    float m_jetEMin=0;
    
    // Call to get collections
    void getCollection(LCCollection*&, std::string, LCEvent*);
    
    TTree * m_outputTree=NULL;
    TFile * m_rootFile=NULL;

 

    float m_Z_mcE=0;
    int m_Z_mcNDaughter=0;
    //only for a ZZ with Z->qq Z->nunu configuration right now
    float m_E_trueNeut=0;
    float m_px_trueNeut=0;
    float m_py_trueNeut=0;
    float m_pz_trueNeut=0;

    float m_E_trueZ2=0;
    float m_px_trueZ2=0;
    float m_py_trueZ2=0;
    float m_pz_trueZ2=0;

    float m_E_totPFO=0;
    float m_E_totPi=0;
    float m_E_totPh=0;
    float m_E_totE=0;
    float m_E_totMu=0;
    float m_E_totN=0;

    float m_px_totPFO=0;
    float m_px_totPi=0;
    float m_px_totPh=0;
    float m_px_totE=0;
    float m_px_totMu=0;
    float m_px_totN=0;

    float m_py_totPFO=0;
    float m_py_totPi=0;
    float m_py_totPh=0;
    float m_py_totE=0;
    float m_py_totMu=0;
    float m_py_totN=0;

    float m_pz_totPFO=0;
    float m_pz_totPi=0;
    float m_pz_totPh=0;
    float m_pz_totE=0;
    float m_pz_totMu=0;
    float m_pz_totN=0;


    float m_E_totPFO_0_95=0;
    float m_E_totPi_0_95=0;
    float m_E_totPh_0_95=0;
    float m_E_totE_0_95=0;
    float m_E_totMu_0_95=0;
    float m_E_totN_0_95=0;
    float m_E_totPFO_0_70=0;
    float m_E_totPi_0_70=0;
    float m_E_totPh_0_70=0;
    float m_E_totE_0_70=0;
    float m_E_totMu_0_70=0;
    float m_E_totN_0_70=0;

    float m_E_true_totAll=0;
    float m_E_true_totInv=0;
    float m_E_true_totPi=0;
    float m_E_true_totPh=0;
    float m_E_true_totK0L=0;
    float m_E_true_totE=0;
    float m_E_true_totMu=0;
    float m_E_true_totN=0;
    float m_E_true_totK=0;
    float m_E_true_totP=0;

    float m_px_true_totAll=0;//no neutrino
    float m_px_true_totInv=0;//no neutrino
    float m_px_true_totPi=0;
    float m_px_true_totP=0;
    float m_px_true_totPh=0;
    float m_px_true_totE=0;
    float m_px_true_totMu=0;
    float m_px_true_totN=0;
    float m_px_true_totK0L=0;
    float m_px_true_totK=0;

    float m_py_true_totAll=0;//no neutrino
    float m_py_true_totInv=0;//no neutrino
    float m_py_true_totPi=0;
    float m_py_true_totP=0;
    float m_py_true_totPh=0;
    float m_py_true_totE=0;
    float m_py_true_totMu=0;
    float m_py_true_totN=0;
    float m_py_true_totK0L=0;
    float m_py_true_totK=0;

    float m_pz_true_totAll=0;//no neutrino
    float m_pz_true_totInv=0;//no neutrino
    float m_pz_true_totPi=0;
    float m_pz_true_totP=0;
    float m_pz_true_totPh=0;
    float m_pz_true_totE=0;
    float m_pz_true_totMu=0;
    float m_pz_true_totN=0;
    float m_pz_true_totK0L=0;
    float m_pz_true_totK=0;

    float m_E_true_totOtherCH=0;
    float m_E_true_totOtherNeut=0;
    float m_E_true_totAll_0_95=0;
    float m_E_true_totInv_0_95=0;
    float m_E_true_totPi_0_95=0;
    float m_E_true_totPh_0_95=0;
    float m_E_true_totK0L_0_95=0;
    float m_E_true_totE_0_95=0;
    float m_E_true_totMu_0_95=0;
    float m_E_true_totN_0_95=0;
    float m_E_true_totK_0_95=0;
    float m_E_true_totP_0_95=0;
    float m_E_true_totOtherCH_0_95=0;
    float m_E_true_totOtherNeut_0_95=0;
    float m_E_true_totAll_0_70=0;
    float m_E_true_totInv_0_70=0;
    float m_E_true_totPi_0_70=0;
    float m_E_true_totPh_0_70=0;
    float m_E_true_totK0L_0_70=0;
    float m_E_true_totE_0_70=0;
    float m_E_true_totMu_0_70=0;
    float m_E_true_totN_0_70=0;
    float m_E_true_totK_0_70=0;
    float m_E_true_totP_0_70=0;
    float m_E_true_totOtherCH_0_70=0;
    float m_E_true_totOtherNeut_0_70=0;

    float m_E_0_2_totPFO=0;
    float m_E_0_2_totPi=0;
    float m_E_0_2_totPh=0;
    float m_E_0_2_totE=0;
    float m_E_0_2_totMu=0;
    float m_E_0_2_totN=0;
    float m_E_0_2_totPFO_0_95=0;
    float m_E_0_2_totPi_0_95=0;
    float m_E_0_2_totPh_0_95=0;
    float m_E_0_2_totE_0_95=0;
    float m_E_0_2_totMu_0_95=0;
    float m_E_0_2_totN_0_95=0;
    float m_E_0_2_totPFO_0_70=0;
    float m_E_0_2_totPi_0_70=0;
    float m_E_0_2_totPh_0_70=0;
    float m_E_0_2_totE_0_70=0;
    float m_E_0_2_totMu_0_70=0;
    float m_E_0_2_totN_0_70=0;

    float m_E_true_0_2_totAll=0;
    float m_E_true_0_2_totInv=0;
    float m_E_true_0_2_totPi=0;
    float m_E_true_0_2_totPh=0;
    float m_E_true_0_2_totK0L=0;
    float m_E_true_0_2_totE=0;
    float m_E_true_0_2_totMu=0;
    float m_E_true_0_2_totN=0;
    float m_E_true_0_2_totK=0;
    float m_E_true_0_2_totP=0;
    float m_E_true_0_2_totOtherCH=0;
    float m_E_true_0_2_totOtherNeut=0;
    float m_E_true_0_2_totAll_0_95=0;
    float m_E_true_0_2_totInv_0_95=0;
    float m_E_true_0_2_totPi_0_95=0;
    float m_E_true_0_2_totPh_0_95=0;
    float m_E_true_0_2_totK0L_0_95=0;
    float m_E_true_0_2_totE_0_95=0;
    float m_E_true_0_2_totMu_0_95=0;
    float m_E_true_0_2_totN_0_95=0;
    float m_E_true_0_2_totK_0_95=0;
    float m_E_true_0_2_totP_0_95=0;
    float m_E_true_0_2_totOtherCH_0_95=0;
    float m_E_true_0_2_totOtherNeut_0_95=0;
    float m_E_true_0_2_totAll_0_70=0;
    float m_E_true_0_2_totInv_0_70=0;
    float m_E_true_0_2_totPi_0_70=0;
    float m_E_true_0_2_totPh_0_70=0;
    float m_E_true_0_2_totK0L_0_70=0;
    float m_E_true_0_2_totE_0_70=0;
    float m_E_true_0_2_totMu_0_70=0;
    float m_E_true_0_2_totN_0_70=0;
    float m_E_true_0_2_totK_0_70=0;
    float m_E_true_0_2_totP_0_70=0;
    float m_E_true_0_2_totOtherCH_0_70=0;
    float m_E_true_0_2_totOtherNeut_0_70=0;

    float m_E_2_10_totPFO=0;
    float m_E_2_10_totPi=0;
    float m_E_2_10_totPh=0;
    float m_E_2_10_totE=0;
    float m_E_2_10_totMu=0;
    float m_E_2_10_totN=0;
    float m_E_2_10_totPFO_0_95=0;
    float m_E_2_10_totPi_0_95=0;
    float m_E_2_10_totPh_0_95=0;
    float m_E_2_10_totE_0_95=0;
    float m_E_2_10_totMu_0_95=0;
    float m_E_2_10_totN_0_95=0;
    float m_E_2_10_totPFO_0_70=0;
    float m_E_2_10_totPi_0_70=0;
    float m_E_2_10_totPh_0_70=0;
    float m_E_2_10_totE_0_70=0;
    float m_E_2_10_totMu_0_70=0;
    float m_E_2_10_totN_0_70=0;

    float m_E_true_2_10_totAll=0;
    float m_E_true_2_10_totInv=0;
    float m_E_true_2_10_totPi=0;
    float m_E_true_2_10_totPh=0;
    float m_E_true_2_10_totK0L=0;
    float m_E_true_2_10_totE=0;
    float m_E_true_2_10_totMu=0;
    float m_E_true_2_10_totN=0;
    float m_E_true_2_10_totK=0;
    float m_E_true_2_10_totP=0;
    float m_E_true_2_10_totOtherCH=0;
    float m_E_true_2_10_totOtherNeut=0;
    float m_E_true_2_10_totAll_0_95=0;
    float m_E_true_2_10_totInv_0_95=0;
    float m_E_true_2_10_totPi_0_95=0;
    float m_E_true_2_10_totPh_0_95=0;
    float m_E_true_2_10_totK0L_0_95=0;
    float m_E_true_2_10_totE_0_95=0;
    float m_E_true_2_10_totMu_0_95=0;
    float m_E_true_2_10_totN_0_95=0;
    float m_E_true_2_10_totK_0_95=0;
    float m_E_true_2_10_totP_0_95=0;
    float m_E_true_2_10_totOtherCH_0_95=0;
    float m_E_true_2_10_totOtherNeut_0_95=0;
    float m_E_true_2_10_totAll_0_70=0;
    float m_E_true_2_10_totInv_0_70=0;
    float m_E_true_2_10_totPi_0_70=0;
    float m_E_true_2_10_totPh_0_70=0;
    float m_E_true_2_10_totK0L_0_70=0;
    float m_E_true_2_10_totE_0_70=0;
    float m_E_true_2_10_totMu_0_70=0;
    float m_E_true_2_10_totN_0_70=0;
    float m_E_true_2_10_totK_0_70=0;
    float m_E_true_2_10_totP_0_70=0;
    float m_E_true_2_10_totOtherCH_0_70=0;
    float m_E_true_2_10_totOtherNeut_0_70=0;

    float m_E_10_50_totPFO=0;
    float m_E_10_50_totPi=0;
    float m_E_10_50_totPh=0;
    float m_E_10_50_totE=0;
    float m_E_10_50_totMu=0;
    float m_E_10_50_totN=0;
    float m_E_10_50_totPFO_0_95=0;
    float m_E_10_50_totPi_0_95=0;
    float m_E_10_50_totPh_0_95=0;
    float m_E_10_50_totE_0_95=0;
    float m_E_10_50_totMu_0_95=0;
    float m_E_10_50_totN_0_95=0;
    float m_E_10_50_totPFO_0_70=0;
    float m_E_10_50_totPi_0_70=0;
    float m_E_10_50_totPh_0_70=0;
    float m_E_10_50_totE_0_70=0;
    float m_E_10_50_totMu_0_70=0;
    float m_E_10_50_totN_0_70=0;

    float m_E_true_10_50_totAll=0;
    float m_E_true_10_50_totInv=0;
    float m_E_true_10_50_totPi=0;
    float m_E_true_10_50_totPh=0;
    float m_E_true_10_50_totK0L=0;
    float m_E_true_10_50_totE=0;
    float m_E_true_10_50_totMu=0;
    float m_E_true_10_50_totN=0;
    float m_E_true_10_50_totK=0;
    float m_E_true_10_50_totP=0;
    float m_E_true_10_50_totOtherCH=0;
    float m_E_true_10_50_totOtherNeut=0;
    float m_E_true_10_50_totAll_0_95=0;
    float m_E_true_10_50_totInv_0_95=0;
    float m_E_true_10_50_totPi_0_95=0;
    float m_E_true_10_50_totPh_0_95=0;
    float m_E_true_10_50_totK0L_0_95=0;
    float m_E_true_10_50_totE_0_95=0;
    float m_E_true_10_50_totMu_0_95=0;
    float m_E_true_10_50_totN_0_95=0;
    float m_E_true_10_50_totK_0_95=0;
    float m_E_true_10_50_totP_0_95=0;
    float m_E_true_10_50_totOtherCH_0_95=0;
    float m_E_true_10_50_totOtherNeut_0_95=0;
    float m_E_true_10_50_totAll_0_70=0;
    float m_E_true_10_50_totInv_0_70=0;
    float m_E_true_10_50_totPi_0_70=0;
    float m_E_true_10_50_totPh_0_70=0;
    float m_E_true_10_50_totK0L_0_70=0;
    float m_E_true_10_50_totE_0_70=0;
    float m_E_true_10_50_totMu_0_70=0;
    float m_E_true_10_50_totN_0_70=0;
    float m_E_true_10_50_totK_0_70=0;
    float m_E_true_10_50_totP_0_70=0;
    float m_E_true_10_50_totOtherCH_0_70=0;
    float m_E_true_10_50_totOtherNeut_0_70=0;

    float m_E_50_totPFO=0;
    float m_E_50_totPi=0;
    float m_E_50_totPh=0;
    float m_E_50_totE=0;
    float m_E_50_totMu=0;
    float m_E_50_totN=0;
    float m_E_50_totPFO_0_95=0;
    float m_E_50_totPi_0_95=0;
    float m_E_50_totPh_0_95=0;
    float m_E_50_totE_0_95=0;
    float m_E_50_totMu_0_95=0;
    float m_E_50_totN_0_95=0;
    float m_E_50_totPFO_0_70=0;
    float m_E_50_totPi_0_70=0;
    float m_E_50_totPh_0_70=0;
    float m_E_50_totE_0_70=0;
    float m_E_50_totMu_0_70=0;
    float m_E_50_totN_0_70=0;

    float m_E_true_50_totAll=0;
    float m_E_true_50_totInv=0;
    float m_E_true_50_totPi=0;
    float m_E_true_50_totPh=0;
    float m_E_true_50_totK0L=0;
    float m_E_true_50_totE=0;
    float m_E_true_50_totMu=0;
    float m_E_true_50_totN=0;
    float m_E_true_50_totK=0;
    float m_E_true_50_totP=0;
    float m_E_true_50_totOtherCH=0;
    float m_E_true_50_totOtherNeut=0;
    float m_E_true_50_totAll_0_95=0;
    float m_E_true_50_totInv_0_95=0;
    float m_E_true_50_totPi_0_95=0;
    float m_E_true_50_totPh_0_95=0;
    float m_E_true_50_totK0L_0_95=0;
    float m_E_true_50_totE_0_95=0;
    float m_E_true_50_totMu_0_95=0;
    float m_E_true_50_totN_0_95=0;
    float m_E_true_50_totK_0_95=0;
    float m_E_true_50_totP_0_95=0;
    float m_E_true_50_totOtherCH_0_95=0;
    float m_E_true_50_totOtherNeut_0_95=0;
    float m_E_true_50_totAll_0_70=0;
    float m_E_true_50_totInv_0_70=0;
    float m_E_true_50_totPi_0_70=0;
    float m_E_true_50_totPh_0_70=0;
    float m_E_true_50_totK0L_0_70=0;
    float m_E_true_50_totE_0_70=0;
    float m_E_true_50_totMu_0_70=0;
    float m_E_true_50_totN_0_70=0;
    float m_E_true_50_totK_0_70=0;
    float m_E_true_50_totP_0_70=0;
    float m_E_true_50_totOtherCH_0_70=0;
    float m_E_true_50_totOtherNeut_0_70=0;

    //multiplicities
    unsigned int m_n_totPFO=0;
    unsigned int m_n_totPi=0;
    unsigned int m_n_totPh=0;
    unsigned int m_n_totE=0;
    unsigned int m_n_totMu=0;
    unsigned int m_n_totN=0;
    unsigned int m_n_totPFO_0_95=0;
    unsigned int m_n_totPi_0_95=0;
    unsigned int m_n_totPh_0_95=0;
    unsigned int m_n_totE_0_95=0;
    unsigned int m_n_totMu_0_95=0;
    unsigned int m_n_totN_0_95=0;
    unsigned int m_n_totPFO_0_70=0;
    unsigned int m_n_totPi_0_70=0;
    unsigned int m_n_totPh_0_70=0;
    unsigned int m_n_totE_0_70=0;
    unsigned int m_n_totMu_0_70=0;
    unsigned int m_n_totN_0_70=0;

    unsigned int m_n_true_totAll=0;
    unsigned int m_n_true_totInv=0;
    unsigned int m_n_true_totPi=0;
    unsigned int m_n_true_totPh=0;
    unsigned int m_n_true_totK0L=0;
    unsigned int m_n_true_totE=0;
    unsigned int m_n_true_totMu=0;
    unsigned int m_n_true_totN=0;
    unsigned int m_n_true_totK=0;
    unsigned int m_n_true_totP=0;
    unsigned int m_n_true_totOtherCH=0;
    unsigned int m_n_true_totOtherNeut=0;
    unsigned int m_n_true_totAll_0_95=0;
    unsigned int m_n_true_totInv_0_95=0;
    unsigned int m_n_true_totPi_0_95=0;
    unsigned int m_n_true_totPh_0_95=0;
    unsigned int m_n_true_totK0L_0_95=0;
    unsigned int m_n_true_totE_0_95=0;
    unsigned int m_n_true_totMu_0_95=0;
    unsigned int m_n_true_totN_0_95=0;
    unsigned int m_n_true_totK_0_95=0;
    unsigned int m_n_true_totP_0_95=0;
    unsigned int m_n_true_totOtherCH_0_95=0;
    unsigned int m_n_true_totOtherNeut_0_95=0;
    unsigned int m_n_true_totAll_0_70=0;
    unsigned int m_n_true_totInv_0_70=0;
    unsigned int m_n_true_totPi_0_70=0;
    unsigned int m_n_true_totPh_0_70=0;
    unsigned int m_n_true_totK0L_0_70=0;
    unsigned int m_n_true_totE_0_70=0;
    unsigned int m_n_true_totMu_0_70=0;
    unsigned int m_n_true_totN_0_70=0;
    unsigned int m_n_true_totK_0_70=0;
    unsigned int m_n_true_totP_0_70=0;
    unsigned int m_n_true_totOtherCH_0_70=0;
    unsigned int m_n_true_totOtherNeut_0_70=0;

    unsigned int m_n_0_2_totPFO=0;
    unsigned int m_n_0_2_totPi=0;
    unsigned int m_n_0_2_totPh=0;
    unsigned int m_n_0_2_totE=0;
    unsigned int m_n_0_2_totMu=0;
    unsigned int m_n_0_2_totN=0;
    unsigned int m_n_0_2_totPFO_0_95=0;
    unsigned int m_n_0_2_totPi_0_95=0;
    unsigned int m_n_0_2_totPh_0_95=0;
    unsigned int m_n_0_2_totE_0_95=0;
    unsigned int m_n_0_2_totMu_0_95=0;
    unsigned int m_n_0_2_totN_0_95=0;
    unsigned int m_n_0_2_totPFO_0_70=0;
    unsigned int m_n_0_2_totPi_0_70=0;
    unsigned int m_n_0_2_totPh_0_70=0;
    unsigned int m_n_0_2_totE_0_70=0;
    unsigned int m_n_0_2_totMu_0_70=0;
    unsigned int m_n_0_2_totN_0_70=0;

    unsigned int m_n_true_0_2_totAll=0;
    unsigned int m_n_true_0_2_totInv=0;
    unsigned int m_n_true_0_2_totPi=0;
    unsigned int m_n_true_0_2_totPh=0;
    unsigned int m_n_true_0_2_totK0L=0;
    unsigned int m_n_true_0_2_totE=0;
    unsigned int m_n_true_0_2_totMu=0;
    unsigned int m_n_true_0_2_totN=0;
    unsigned int m_n_true_0_2_totK=0;
    unsigned int m_n_true_0_2_totP=0;
    unsigned int m_n_true_0_2_totOtherCH=0;
    unsigned int m_n_true_0_2_totOtherNeut=0;
    unsigned int m_n_true_0_2_totAll_0_95=0;
    unsigned int m_n_true_0_2_totInv_0_95=0;
    unsigned int m_n_true_0_2_totPi_0_95=0;
    unsigned int m_n_true_0_2_totPh_0_95=0;
    unsigned int m_n_true_0_2_totK0L_0_95=0;
    unsigned int m_n_true_0_2_totE_0_95=0;
    unsigned int m_n_true_0_2_totMu_0_95=0;
    unsigned int m_n_true_0_2_totN_0_95=0;
    unsigned int m_n_true_0_2_totK_0_95=0;
    unsigned int m_n_true_0_2_totP_0_95=0;
    unsigned int m_n_true_0_2_totOtherCH_0_95=0;
    unsigned int m_n_true_0_2_totOtherNeut_0_95=0;
    unsigned int m_n_true_0_2_totAll_0_70=0;
    unsigned int m_n_true_0_2_totInv_0_70=0;
    unsigned int m_n_true_0_2_totPi_0_70=0;
    unsigned int m_n_true_0_2_totPh_0_70=0;
    unsigned int m_n_true_0_2_totK0L_0_70=0;
    unsigned int m_n_true_0_2_totE_0_70=0;
    unsigned int m_n_true_0_2_totMu_0_70=0;
    unsigned int m_n_true_0_2_totN_0_70=0;
    unsigned int m_n_true_0_2_totK_0_70=0;
    unsigned int m_n_true_0_2_totP_0_70=0;
    unsigned int m_n_true_0_2_totOtherCH_0_70=0;
    unsigned int m_n_true_0_2_totOtherNeut_0_70=0;

    unsigned int m_n_2_10_totPFO=0;
    unsigned int m_n_2_10_totPi=0;
    unsigned int m_n_2_10_totPh=0;
    unsigned int m_n_2_10_totE=0;
    unsigned int m_n_2_10_totMu=0;
    unsigned int m_n_2_10_totN=0;
    unsigned int m_n_2_10_totPFO_0_95=0;
    unsigned int m_n_2_10_totPi_0_95=0;
    unsigned int m_n_2_10_totPh_0_95=0;
    unsigned int m_n_2_10_totE_0_95=0;
    unsigned int m_n_2_10_totMu_0_95=0;
    unsigned int m_n_2_10_totN_0_95=0;
    unsigned int m_n_2_10_totPFO_0_70=0;
    unsigned int m_n_2_10_totPi_0_70=0;
    unsigned int m_n_2_10_totPh_0_70=0;
    unsigned int m_n_2_10_totE_0_70=0;
    unsigned int m_n_2_10_totMu_0_70=0;
    unsigned int m_n_2_10_totN_0_70=0;

    unsigned int m_n_true_2_10_totAll=0;
    unsigned int m_n_true_2_10_totInv=0;
    unsigned int m_n_true_2_10_totPi=0;
    unsigned int m_n_true_2_10_totPh=0;
    unsigned int m_n_true_2_10_totK0L=0;
    unsigned int m_n_true_2_10_totE=0;
    unsigned int m_n_true_2_10_totMu=0;
    unsigned int m_n_true_2_10_totN=0;
    unsigned int m_n_true_2_10_totK=0;
    unsigned int m_n_true_2_10_totP=0;
    unsigned int m_n_true_2_10_totOtherCH=0;
    unsigned int m_n_true_2_10_totOtherNeut=0;
    unsigned int m_n_true_2_10_totAll_0_95=0;
    unsigned int m_n_true_2_10_totInv_0_95=0;
    unsigned int m_n_true_2_10_totPi_0_95=0;
    unsigned int m_n_true_2_10_totPh_0_95=0;
    unsigned int m_n_true_2_10_totK0L_0_95=0;
    unsigned int m_n_true_2_10_totE_0_95=0;
    unsigned int m_n_true_2_10_totMu_0_95=0;
    unsigned int m_n_true_2_10_totN_0_95=0;
    unsigned int m_n_true_2_10_totK_0_95=0;
    unsigned int m_n_true_2_10_totP_0_95=0;
    unsigned int m_n_true_2_10_totOtherCH_0_95=0;
    unsigned int m_n_true_2_10_totOtherNeut_0_95=0;
    unsigned int m_n_true_2_10_totAll_0_70=0;
    unsigned int m_n_true_2_10_totInv_0_70=0;
    unsigned int m_n_true_2_10_totPi_0_70=0;
    unsigned int m_n_true_2_10_totPh_0_70=0;
    unsigned int m_n_true_2_10_totK0L_0_70=0;
    unsigned int m_n_true_2_10_totE_0_70=0;
    unsigned int m_n_true_2_10_totMu_0_70=0;
    unsigned int m_n_true_2_10_totN_0_70=0;
    unsigned int m_n_true_2_10_totK_0_70=0;
    unsigned int m_n_true_2_10_totP_0_70=0;
    unsigned int m_n_true_2_10_totOtherCH_0_70=0;
    unsigned int m_n_true_2_10_totOtherNeut_0_70=0;

    unsigned int m_n_10_50_totPFO=0;
    unsigned int m_n_10_50_totPi=0;
    unsigned int m_n_10_50_totPh=0;
    unsigned int m_n_10_50_totE=0;
    unsigned int m_n_10_50_totMu=0;
    unsigned int m_n_10_50_totN=0;
    unsigned int m_n_10_50_totPFO_0_95=0;
    unsigned int m_n_10_50_totPi_0_95=0;
    unsigned int m_n_10_50_totPh_0_95=0;
    unsigned int m_n_10_50_totE_0_95=0;
    unsigned int m_n_10_50_totMu_0_95=0;
    unsigned int m_n_10_50_totN_0_95=0;
    unsigned int m_n_10_50_totPFO_0_70=0;
    unsigned int m_n_10_50_totPi_0_70=0;
    unsigned int m_n_10_50_totPh_0_70=0;
    unsigned int m_n_10_50_totE_0_70=0;
    unsigned int m_n_10_50_totMu_0_70=0;
    unsigned int m_n_10_50_totN_0_70=0;

    unsigned int m_n_true_10_50_totAll=0;
    unsigned int m_n_true_10_50_totInv=0;
    unsigned int m_n_true_10_50_totPi=0;
    unsigned int m_n_true_10_50_totPh=0;
    unsigned int m_n_true_10_50_totK0L=0;
    unsigned int m_n_true_10_50_totE=0;
    unsigned int m_n_true_10_50_totMu=0;
    unsigned int m_n_true_10_50_totN=0;
    unsigned int m_n_true_10_50_totK=0;
    unsigned int m_n_true_10_50_totP=0;
    unsigned int m_n_true_10_50_totOtherCH=0;
    unsigned int m_n_true_10_50_totOtherNeut=0;
    unsigned int m_n_true_10_50_totAll_0_95=0;
    unsigned int m_n_true_10_50_totInv_0_95=0;
    unsigned int m_n_true_10_50_totPi_0_95=0;
    unsigned int m_n_true_10_50_totPh_0_95=0;
    unsigned int m_n_true_10_50_totK0L_0_95=0;
    unsigned int m_n_true_10_50_totE_0_95=0;
    unsigned int m_n_true_10_50_totMu_0_95=0;
    unsigned int m_n_true_10_50_totN_0_95=0;
    unsigned int m_n_true_10_50_totK_0_95=0;
    unsigned int m_n_true_10_50_totP_0_95=0;
    unsigned int m_n_true_10_50_totOtherCH_0_95=0;
    unsigned int m_n_true_10_50_totOtherNeut_0_95=0;
    unsigned int m_n_true_10_50_totAll_0_70=0;
    unsigned int m_n_true_10_50_totInv_0_70=0;
    unsigned int m_n_true_10_50_totPi_0_70=0;
    unsigned int m_n_true_10_50_totPh_0_70=0;
    unsigned int m_n_true_10_50_totK0L_0_70=0;
    unsigned int m_n_true_10_50_totE_0_70=0;
    unsigned int m_n_true_10_50_totMu_0_70=0;
    unsigned int m_n_true_10_50_totN_0_70=0;
    unsigned int m_n_true_10_50_totK_0_70=0;
    unsigned int m_n_true_10_50_totP_0_70=0;
    unsigned int m_n_true_10_50_totOtherCH_0_70=0;
    unsigned int m_n_true_10_50_totOtherNeut_0_70=0;

    unsigned int m_n_50_totPFO=0;
    unsigned int m_n_50_totPi=0;
    unsigned int m_n_50_totPh=0;
    unsigned int m_n_50_totE=0;
    unsigned int m_n_50_totMu=0;
    unsigned int m_n_50_totN=0;
    unsigned int m_n_50_totPFO_0_95=0;
    unsigned int m_n_50_totPi_0_95=0;
    unsigned int m_n_50_totPh_0_95=0;
    unsigned int m_n_50_totE_0_95=0;
    unsigned int m_n_50_totMu_0_95=0;
    unsigned int m_n_50_totN_0_95=0;
    unsigned int m_n_50_totPFO_0_70=0;
    unsigned int m_n_50_totPi_0_70=0;
    unsigned int m_n_50_totPh_0_70=0;
    unsigned int m_n_50_totE_0_70=0;
    unsigned int m_n_50_totMu_0_70=0;
    unsigned int m_n_50_totN_0_70=0;

    unsigned int m_n_true_50_totAll=0;
    unsigned int m_n_true_50_totInv=0;
    unsigned int m_n_true_50_totPi=0;
    unsigned int m_n_true_50_totPh=0;
    unsigned int m_n_true_50_totK0L=0;
    unsigned int m_n_true_50_totE=0;
    unsigned int m_n_true_50_totMu=0;
    unsigned int m_n_true_50_totN=0;
    unsigned int m_n_true_50_totK=0;
    unsigned int m_n_true_50_totP=0;
    unsigned int m_n_true_50_totOtherCH=0;
    unsigned int m_n_true_50_totOtherNeut=0;
    unsigned int m_n_true_50_totAll_0_95=0;
    unsigned int m_n_true_50_totInv_0_95=0;
    unsigned int m_n_true_50_totPi_0_95=0;
    unsigned int m_n_true_50_totPh_0_95=0;
    unsigned int m_n_true_50_totK0L_0_95=0;
    unsigned int m_n_true_50_totE_0_95=0;
    unsigned int m_n_true_50_totMu_0_95=0;
    unsigned int m_n_true_50_totN_0_95=0;
    unsigned int m_n_true_50_totK_0_95=0;
    unsigned int m_n_true_50_totP_0_95=0;
    unsigned int m_n_true_50_totOtherCH_0_95=0;
    unsigned int m_n_true_50_totOtherNeut_0_95=0;
    unsigned int m_n_true_50_totAll_0_70=0;
    unsigned int m_n_true_50_totInv_0_70=0;
    unsigned int m_n_true_50_totPi_0_70=0;
    unsigned int m_n_true_50_totPh_0_70=0;
    unsigned int m_n_true_50_totK0L_0_70=0;
    unsigned int m_n_true_50_totE_0_70=0;
    unsigned int m_n_true_50_totMu_0_70=0;
    unsigned int m_n_true_50_totN_0_70=0;
    unsigned int m_n_true_50_totK_0_70=0;
    unsigned int m_n_true_50_totP_0_70=0;
    unsigned int m_n_true_50_totOtherCH_0_70=0;
    unsigned int m_n_true_50_totOtherNeut_0_70=0;

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

    std::vector<float>* m_trueME_E=NULL;
    std::vector<float>* m_trueME_Px=NULL;
    std::vector<float>* m_trueME_Py=NULL;
    std::vector<float>* m_trueME_Pz=NULL;
    std::vector<int>* m_trueME_PDGID=NULL;

    std::vector<float>* m_truePi0Energy=NULL;
    std::vector<float>* m_truePi0_Px=NULL;
    std::vector<float>* m_truePi0_Py=NULL;
    std::vector<float>* m_truePi0_Pz=NULL;
    std::vector<float>* m_truePi0_CosTheta=NULL;
    std::vector<float>* m_truePi0_Phi=NULL;
    std::vector<float>* m_truePi0_DRPh01=NULL;
    std::vector<float>* m_truePi0_CosAngPh01=NULL;
    std::vector<float>* m_truePi0_AngPh01=NULL;
    std::vector<float>* m_truePh0Energy=NULL;
    std::vector<float>* m_truePh0_Px=NULL;
    std::vector<float>* m_truePh0_Py=NULL;
    std::vector<float>* m_truePh0_Pz=NULL;
    std::vector<float>* m_truePh0_CosTheta=NULL;
    std::vector<float>* m_truePh0_Phi=NULL;
    std::vector<float>* m_truePh0_DRMin=NULL;
    std::vector<float>* m_truePh0_DRMinPDGID=NULL;
    std::vector<float>* m_truePh0_DRMinE=NULL;
    std::vector<float>* m_truePh0_CosAngMax=NULL;
    std::vector<float>* m_truePh0_CosAngMaxPDGID=NULL;
    std::vector<float>* m_truePh0_CosAngMaxE=NULL;
    std::vector<float>* m_truePh0_CosAng0995E=NULL;
    std::vector<float>* m_truePh0_DR01E=NULL;
    std::vector<float>* m_truePh1Energy=NULL;
    std::vector<float>* m_truePh1_Px=NULL;
    std::vector<float>* m_truePh1_Py=NULL;
    std::vector<float>* m_truePh1_Pz=NULL;
    std::vector<float>* m_truePh1_CosTheta=NULL;
    std::vector<float>* m_truePh1_Phi=NULL;
    std::vector<float>* m_truePh1_DRMin=NULL;
    std::vector<float>* m_truePh1_DRMinPDGID=NULL;
    std::vector<float>* m_truePh1_DRMinE=NULL;
    std::vector<float>* m_truePh1_CosAngMax=NULL;
    std::vector<float>* m_truePh1_CosAngMaxPDGID=NULL;
    std::vector<float>* m_truePh1_CosAngMaxE=NULL;
    std::vector<float>* m_truePh1_CosAng0995E=NULL;
    std::vector<float>* m_truePh1_DR01E=NULL;


    std::vector<float>* m_recoPhEnergy=NULL;
    std::vector<float>* m_recoPh_Px=NULL;
    std::vector<float>* m_recoPh_Py=NULL;
    std::vector<float>* m_recoPh_Pz=NULL;
    std::vector<float>* m_recoPh_CosTheta=NULL;
    std::vector<float>* m_recoPh_Theta=NULL;
    std::vector<float>* m_recoPh_Phi=NULL;
    std::vector<float>* m_recoPh_PDGID=NULL;
    std::vector<float>* m_recoPh_E_EB=NULL;
    std::vector<float>* m_recoPh_E_EE=NULL;
    std::vector<float>* m_recoPh_E_EO=NULL;
    std::vector<float>* m_recoPh_E_HB=NULL;
    std::vector<float>* m_recoPh_E_HE=NULL;
    std::vector<float>* m_recoPh_E_HO=NULL;
    std::vector<int>* m_recoPh_firstLayerECAL=NULL;
    std::vector<int>* m_recoPh_lastLayerECAL=NULL;
    std::vector<int>* m_recoPh_nhitsEB=NULL;
    std::vector<int>* m_recoPh_nhitsEE=NULL;
    std::vector<int>* m_recoPh_nhitsEO=NULL;
    std::vector<int>* m_recoPh_firstLayerHCAL=NULL;
    std::vector<int>* m_recoPh_lastLayerHCAL=NULL;
    std::vector<int>* m_recoPh_nhitsHB=NULL;
    std::vector<int>* m_recoPh_nhitsHE=NULL;
    std::vector<int>* m_recoPh_nhitsHO=NULL;
    std::vector<float>* m_recoPh_DR01_E=NULL;//isolation including all particles, also the next closest photon
    std::vector<float>* m_recoPh_CosAng0995_E=NULL;//isolation including ALL particles, also the next closest photon
    //values for the closest particle other than the next closest photon
    std::vector<float>* m_recoPh_DRMin=NULL;
    std::vector<float>* m_recoPh_DRMin_PDGID=NULL;
    std::vector<float>* m_recoPh_DRMin_E=NULL;
    std::vector<float>* m_recoPh_CosAngMax=NULL;
    std::vector<float>* m_recoPh_CosAngMax_PDGID=NULL;
    std::vector<float>* m_recoPh_CosAngMax_E=NULL;
    //values for the next closest photon
    std::vector<float>* m_recoPh_DRMin_Ph=NULL;
    std::vector<float>* m_recoPh_DRMin_E_Ph=NULL;
    std::vector<float>* m_recoPi0Cand_DRMin_M=NULL;
    std::vector<float>* m_recoPh_CosAngMax_Ph=NULL;
    std::vector<float>* m_recoPh_CosAngMax_E_Ph=NULL;
    std::vector<float>* m_recoPi0Cand_CosAngMax_M=NULL;

    
    std::vector<float>* m_jet_excl_Px=NULL; 
    std::vector<float>* m_jet_excl_Py=NULL; 
    std::vector<float>* m_jet_excl_Pz=NULL; 
    std::vector<float>* m_jet_excl_E=NULL; 
    std::vector<float>* m_jet_excl_Phi=NULL; 
    std::vector<float>* m_jet_excl_CosTheta=NULL; 
    std::vector<float>* m_jet_excl_piE=NULL;
    std::vector<float>* m_jet_excl_phE=NULL;
    std::vector<float>* m_jet_excl_elE=NULL;
    std::vector<float>* m_jet_excl_muE=NULL;
    std::vector<float>* m_jet_excl_nE=NULL;
    std::vector<float>* m_jet_excl_neutElseE=NULL;
    std::vector<float>* m_jet_excl_chElseE=NULL;
    std::vector<int>* m_jet_excl_neutMult=NULL;
    std::vector<int>* m_jet_excl_chMult=NULL;

    std::vector<float>* m_jet_incl_Px=NULL; 
    std::vector<float>* m_jet_incl_Py=NULL; 
    std::vector<float>* m_jet_incl_Pz=NULL; 
    std::vector<float>* m_jet_incl_E=NULL; 
    std::vector<float>* m_jet_incl_Phi=NULL; 
    std::vector<float>* m_jet_incl_CosTheta=NULL; 
    std::vector<float>* m_jet_incl_piE=NULL;
    std::vector<float>* m_jet_incl_phE=NULL;
    std::vector<float>* m_jet_incl_elE=NULL;
    std::vector<float>* m_jet_incl_muE=NULL;
    std::vector<float>* m_jet_incl_nE=NULL;
    std::vector<float>* m_jet_incl_neutElseE=NULL;
    std::vector<float>* m_jet_incl_chElseE=NULL;
    std::vector<int>* m_jet_incl_neutMult=NULL;
    std::vector<int>* m_jet_incl_chMult=NULL;



} ;

#endif



