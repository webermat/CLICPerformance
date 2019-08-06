#ifndef TrueHZParticlesForJets_h
#define TrueHZParticlesForJets_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "IMPL/ReconstructedParticleImpl.h"
#include <set>

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJetStructureBase.hh>
#include <fastjet/contrib/ValenciaPlugin.hh>

using namespace lcio ;
using namespace marlin ;

typedef std::vector< fastjet::PseudoJet > PseudoJetList;


class TrueHZParticlesForJets : public Processor {
                
public:
        
    virtual Processor*  newProcessor() { return new TrueHZParticlesForJets ; }
    
    TrueHZParticlesForJets() ;

    TrueHZParticlesForJets(const TrueHZParticlesForJets&) = delete;
    TrueHZParticlesForJets& operator=(const TrueHZParticlesForJets&) = delete;
    
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
    //MC Jet input collection, remove here isolated leptons and photons, isolation determined by relative isolation in cone of 10 degrees
    std::string m_outputMCForJetParticleCollection="";
    //take out isolated leptons/photons --> don't consider tau's for now
    std::string m_outputMCIsoLepPhParticleCollection="";
    //here the "true" leptons and photons from ISR are collected, as well as tau decay particles
    std::string m_outputMCTrueLepPhParticleCollection="";
    //here the "true" leptons and photons from ISR are removed, as well as tau decay particles (from Z,W and H)
    std::string m_outputMCTrueJetParticleCollection="";

    std::string m_inputRecoParticleCollection="";
    std::string m_outputRecoJetParticleCollection="";
    std::string m_outputRecoIsoLepPhParticleCollection="";

    float m_R=1.2;
    float m_beta=1.0;
    float m_gamma=1.0;

    int m_NJets =2;

    std::string m_outputPFOsInJetsCollection="";
    std::string m_outputPFOJet0Collection="";
    std::string m_outputPFOJet1Collection="";


    std::string m_inputTauCollection="";
 
    // Call to get collections
    void getCollection(LCCollection*&, std::string, LCEvent*);

    void fillStableDaughterSet(MCParticle*, std::set<MCParticle*> &);

    float m_angleIso = 10.4;
    float m_relIso = 0.101;
    float m_minE =1.5;

    fastjet::contrib::ValenciaPlugin* vlcpl=NULL;
    fastjet::JetDefinition* _jetAlgoType = NULL;

    /** Replace missing copy constructor by hand */
    ReconstructedParticleImpl* CopyRecoParticle ( ReconstructedParticle* pfo ) ;

} ;

#endif



