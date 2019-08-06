#ifndef ZudsHistos_h
#define ZudsHistos_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>


#include "TH2F.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <map>

using namespace lcio ;
using namespace marlin ;

class TTree;


class ZudsHistos : public Processor {
                
public:
        
    virtual Processor*  newProcessor() { return new ZudsHistos ; }
    
    ZudsHistos() ;

    ZudsHistos(const ZudsHistos&) = delete;
    ZudsHistos& operator=(const ZudsHistos&) = delete;
    
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


    TH1D* h_truePi_E= NULL;
    TH1D* h_truePr_E = NULL;
    TH1D* h_trueCH_E = NULL;
    TH1D* h_truePh_E = NULL;
    TH1D* h_trueNe_E = NULL;
    TH1D* h_trueK0L_E  =NULL;
    TH1D* h_trueNH_E = NULL;
    TH1D* h_trueEl_E = NULL;
    TH1D* h_trueMu_E = NULL;
    
    TH1D* h_truePi_Pt = NULL;
    TH1D* h_truePr_Pt = NULL;
    TH1D* h_trueCH_Pt = NULL;
    
    TH1D* h_recoPi_E = NULL;
    TH1D* h_recoNe_E = NULL;
    TH1D* h_recoPh_E = NULL;
    TH1D* h_recoEl_E = NULL;
    TH1D* h_recoMu_E = NULL;
    
    TH1D* h_recoPi_Pt = NULL;

    TH1D* h_truePh_nextPart_PID = NULL;

    TH1D* h_truePh_E_0_2_DA_0_01_N_part_vs_CosTheta=NULL;
    TH1D* h_truePh_E_0_2_DA_0_01_Ph_closest_vs_CosTheta=NULL;
    TH1D* h_truePh_E_0_2_DA_0_01_N_closest_vs_CosTheta=NULL;
    TH1D* h_truePh_E_0_2_DA_0_01__closest_vs_CosTheta=NULL;

    int eventcount=0;
    float m_jetEMin=0;
    
    // Call to get collections
    void getCollection(LCCollection*&, std::string, LCEvent*);
    
    
    TFile * m_rootFile=NULL;

} ;

#endif



