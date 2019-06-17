#ifndef RunEventStatisticsHistogram_h
#define RunEventStatisticsHistogram_h 1

#include "marlin/Processor.h"

#include "TH1F.h"
#include "TFile.h"


using namespace lcio ;
using namespace marlin ;


class RunEventStatisticsHistogram : public Processor {
                
public:
        
    virtual Processor*  newProcessor() { return new RunEventStatisticsHistogram ; }
    
    RunEventStatisticsHistogram() ;

    RunEventStatisticsHistogram(const RunEventStatisticsHistogram&) = delete;
    RunEventStatisticsHistogram& operator=(const RunEventStatisticsHistogram&) = delete;
    
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
        
    TH1F* runstatistics;

    std::string m_rootFileName="";

    TFile * m_rootFile=NULL;
    float eventcount=0;

} ;

#endif



