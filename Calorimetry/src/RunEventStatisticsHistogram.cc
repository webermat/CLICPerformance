#include "RunEventStatisticsHistogram.h"
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


RunEventStatisticsHistogram aRunEventStatisticsHistogram;

RunEventStatisticsHistogram::RunEventStatisticsHistogram() : Processor("RunEventStatisticsHistogram") {
    
    // modify processor description
    _description = "runstatistics histograms" ;
      
    registerProcessorParameter( "OutputRootFileName",
                                "ROOT File name to save the overview histogram",
                                m_rootFileName,
                                std::string("RunEventStatisticsHistogram.root"));

}  


void RunEventStatisticsHistogram::init() {

  m_rootFile = new TFile(m_rootFileName.c_str(),"RECREATE");

  runstatistics = new TH1F("h_runstatistics","",11,-0.5,10.5);
  eventcount=0;

}


void RunEventStatisticsHistogram::processRunHeader( LCRunHeader*) {
}

void RunEventStatisticsHistogram::processEvent( LCEvent* evt ) {
    
  eventcount+=1;

}

void RunEventStatisticsHistogram::check(LCEvent*){
}

void RunEventStatisticsHistogram::end(){
    

  runstatistics->SetBinContent(1,eventcount);
  m_rootFile->cd();  
 
  m_rootFile->Write();
  m_rootFile->Close();

}
