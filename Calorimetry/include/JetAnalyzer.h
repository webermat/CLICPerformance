#ifndef JetAnalyzer_h
#define JetAnalyzer_h 1

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


class JetAnalyzer : public Processor {
                
public:
        
    virtual Processor*  newProcessor() { return new JetAnalyzer ; }
    
    JetAnalyzer() ;

    JetAnalyzer(const JetAnalyzer&) = delete;
    JetAnalyzer& operator=(const JetAnalyzer&) = delete;
    
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
    std::string m_genjetDR03ColName="";
    std::string m_genjetDR04ColName="";
    std::string m_genjetDR05ColName="";
    std::string m_genjetDR06ColName="";
    std::string m_genjetDR07ColName="";
    std::string m_genjetDR08ColName="";
    std::string m_genjetDR09ColName="";
    std::string m_genjetDR10ColName="";

    std::string m_recojetDR03ColName="";
    std::string m_recojetDR04ColName="";
    std::string m_recojetDR05ColName="";
    std::string m_recojetDR06ColName="";
    std::string m_recojetDR07ColName="";
    std::string m_recojetDR08ColName="";
    std::string m_recojetDR09ColName="";
    std::string m_recojetDR10ColName="";

    std::string m_rootFileName="";
    
    // Run and event counters
    int m_eventNumber=0;
    int m_runNumber=0;
    int m_eventCount=0;
    //parton level information, save two Z daughters
    int m_d1_mcPDGID=0;
    float m_d1_mcE=0;
    float m_d1_mcPx=0;
    float m_d1_mcPy=0;
    float m_d1_mcPz=0;

    int m_d2_mcPDGID=0;
    float m_d2_mcE=0;
    float m_d2_mcPx=0;
    float m_d2_mcPy=0;
    float m_d2_mcPz=0;
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

    //for later cleaning cut with overlay
    float m_E_trueInv_cos09=0;
    float m_px_trueInv_cos09=0;
    float m_py_trueInv_cos09=0;
    float m_pz_trueInv_cos09=0;
    float m_E_trueAll_cos09=0;
    float m_px_trueAll_cos09=0;
    float m_py_trueAll_cos09=0;
    float m_pz_trueAll_cos09=0;
    float m_E_totPFO_cos09=0;
    float m_px_totPFO_cos09=0;
    float m_py_totPFO_cos09=0;
    float m_pz_totPFO_cos09=0;

    std::vector<float>* m_trueME_E=NULL;
    std::vector<float>* m_trueME_Px=NULL;
    std::vector<float>* m_trueME_Py=NULL;
    std::vector<float>* m_trueME_Pz=NULL;
    std::vector<int>* m_trueME_PDGID=NULL;

    std::vector<float>* m_genJetDR03_E=NULL;
    std::vector<float>* m_genJetDR03_Px=NULL;
    std::vector<float>* m_genJetDR03_Py=NULL;
    std::vector<float>* m_genJetDR03_Pz=NULL;

    std::vector<float>* m_genJetDR04_E=NULL;
    std::vector<float>* m_genJetDR04_Px=NULL;
    std::vector<float>* m_genJetDR04_Py=NULL;
    std::vector<float>* m_genJetDR04_Pz=NULL;

    std::vector<float>* m_genJetDR05_E=NULL;
    std::vector<float>* m_genJetDR05_Px=NULL;
    std::vector<float>* m_genJetDR05_Py=NULL;
    std::vector<float>* m_genJetDR05_Pz=NULL;

    std::vector<float>* m_genJetDR06_E=NULL;
    std::vector<float>* m_genJetDR06_Px=NULL;
    std::vector<float>* m_genJetDR06_Py=NULL;
    std::vector<float>* m_genJetDR06_Pz=NULL;

    std::vector<float>* m_genJetDR07_E=NULL;
    std::vector<float>* m_genJetDR07_Px=NULL;
    std::vector<float>* m_genJetDR07_Py=NULL;
    std::vector<float>* m_genJetDR07_Pz=NULL;

    std::vector<float>* m_genJetDR08_E=NULL;
    std::vector<float>* m_genJetDR08_Px=NULL;
    std::vector<float>* m_genJetDR08_Py=NULL;
    std::vector<float>* m_genJetDR08_Pz=NULL;

    std::vector<float>* m_genJetDR09_E=NULL;
    std::vector<float>* m_genJetDR09_Px=NULL;
    std::vector<float>* m_genJetDR09_Py=NULL;
    std::vector<float>* m_genJetDR09_Pz=NULL;

    std::vector<float>* m_genJetDR10_E=NULL;
    std::vector<float>* m_genJetDR10_Px=NULL;
    std::vector<float>* m_genJetDR10_Py=NULL;
    std::vector<float>* m_genJetDR10_Pz=NULL;

    std::vector<float>* m_recoJetDR03_E=NULL;
    std::vector<float>* m_recoJetDR03_Px=NULL;
    std::vector<float>* m_recoJetDR03_Py=NULL;
    std::vector<float>* m_recoJetDR03_Pz=NULL;

    std::vector<float>* m_recoJetDR04_E=NULL;
    std::vector<float>* m_recoJetDR04_Px=NULL;
    std::vector<float>* m_recoJetDR04_Py=NULL;
    std::vector<float>* m_recoJetDR04_Pz=NULL;

    std::vector<float>* m_recoJetDR05_E=NULL;
    std::vector<float>* m_recoJetDR05_Px=NULL;
    std::vector<float>* m_recoJetDR05_Py=NULL;
    std::vector<float>* m_recoJetDR05_Pz=NULL;

    std::vector<float>* m_recoJetDR06_E=NULL;
    std::vector<float>* m_recoJetDR06_Px=NULL;
    std::vector<float>* m_recoJetDR06_Py=NULL;
    std::vector<float>* m_recoJetDR06_Pz=NULL;

    std::vector<float>* m_recoJetDR07_E=NULL;
    std::vector<float>* m_recoJetDR07_Px=NULL;
    std::vector<float>* m_recoJetDR07_Py=NULL;
    std::vector<float>* m_recoJetDR07_Pz=NULL;

    std::vector<float>* m_recoJetDR08_E=NULL;
    std::vector<float>* m_recoJetDR08_Px=NULL;
    std::vector<float>* m_recoJetDR08_Py=NULL;
    std::vector<float>* m_recoJetDR08_Pz=NULL;

    std::vector<float>* m_recoJetDR09_E=NULL;
    std::vector<float>* m_recoJetDR09_Px=NULL;
    std::vector<float>* m_recoJetDR09_Py=NULL;
    std::vector<float>* m_recoJetDR09_Pz=NULL;

    std::vector<float>* m_recoJetDR10_E=NULL;
    std::vector<float>* m_recoJetDR10_Px=NULL;
    std::vector<float>* m_recoJetDR10_Py=NULL;
    std::vector<float>* m_recoJetDR10_Pz=NULL;

    int eventcount=0;
    float m_jetEMin=0;
    
    // Call to get collections
    void getCollection(LCCollection*&, std::string, LCEvent*);
    
    
    TFile * m_rootFile=NULL;

    TTree* m_outputTree=NULL;

    bool m_fillMEInfo=false;
    bool m_fillAllJets=true;

} ;

#endif



