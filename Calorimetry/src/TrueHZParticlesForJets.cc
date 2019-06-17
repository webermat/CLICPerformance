#include "TrueHZParticlesForJets.h"
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/MCParticleImpl.h>
#include <fastjet/contrib/MeasureDefinition.hh>
#include <fastjet/contrib/MeasureDefinition.hh>

#include "DDRec/DetectorData.h"

#include "TLorentzVector.h"

using namespace lcio ;
using namespace marlin ;


TrueHZParticlesForJets aTrueHZParticlesForJets;

TrueHZParticlesForJets::TrueHZParticlesForJets() : Processor("TrueHZParticlesForJets") {
    
    // modify processor description
    _description = "TrueHZParticlesForJets calculates properties of calorimeter showers" ;
   

    registerInputCollection( LCIO::MCPARTICLE,
                            "MCParticleInputCollectionName",
                            "Name of the MCParticle input collection",
                            m_inputMCParticleCollection,
                            std::string("MCPhysicsParticles"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "MCForJetParticleCollectionName",
                            "Name of the MCForJets output RecoParticle collection",
                            m_outputMCForJetParticleCollection,
                            std::string("MCParticlePandoraPFOsForJets"));

    registerInputCollection( LCIO::MCPARTICLE,
                            "MCIsoLepPhParticleCollection",
                            "Name of the MCParticle IsoLeptonAndPhoton output collection",
                            m_outputMCIsoLepPhParticleCollection,
                            std::string("MCIsoLepPhParticles"));

    registerInputCollection( LCIO::MCPARTICLE,
                            "MCTrueLepPhParticleCollection",
                            "Name of the MCParticle TrueLeptonAndPhoton output collection",
                            m_outputMCTrueLepPhParticleCollection,
                            std::string("MCTrueLepPhParticles"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "MCTrueJetParticleCollection",
                            "Name of the MC true particle jet output collection",
			     m_outputMCTrueJetParticleCollection ,
                            std::string("MCTrueForJetsPandoraPFOs"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RecoParticleInputCollectionName",
                            "Name of the RecoParticle input collection",
                            m_inputRecoParticleCollection,
                            std::string("TightSelectedPandoraPFOs"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetParticleCollectionName",
			     "Name of the RecoJetParticle output collection",
			     m_outputRecoJetParticleCollection,
			     std::string("PandoraPFOsForJets"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RecoParticleIsoLepPhCollectionName",
                            "Name of the RecoIsoLepPhTauParticle output collection",
			     m_outputRecoIsoLepPhParticleCollection,
                            std::string("PandoraPFOsIsoLepPh"));


   registerProcessorParameter(
			       "isoAngle" , 
			       "isolation angle",
			       m_angleIso,
			       float(10.1)			       
			       );

   registerProcessorParameter(
			       "RelIso" , 
			       "relative isolation value",
			       m_relIso,
			       float(0.11)			       
			       );
   registerProcessorParameter(
			       "Emin" , 
			       "minimum energy for isolation checks",
			       m_minE,
			       float(1.5)
			      );
   registerProcessorParameter(
			      "R", 
			      "jet radius",
			      m_R,
			      float(0.7)
			      );
   registerProcessorParameter(
			      "beta", 
			      "jet beta parameter",
			      m_beta,
			      float(1.0)
			      );

   registerProcessorParameter(
			      "gamma", 
			      "jet gamma parameter",
			      m_gamma,
			      float(1.0)
			      );
   registerProcessorParameter(
			       "NJets" , 
			       "number of fat jets in exclusive clustering",
			       m_NJets,
			       int(2)
			      );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RecoParticleJet0",
                            "Name of the Jet0 PFO constituents",
			     m_outputPFOJet0Collection,
                            std::string("PandoraPFOsJet0"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RecoParticleJet1",
                            "Name of the Jet1 PFO constituents",
			     m_outputPFOJet1Collection,
                            std::string("PandoraPFOsJet1"));


    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoParticlesInJets",
			     "Name PFOsInJets collection",
			     m_outputPFOsInJetsCollection,
			     std::string("PandoraPFOsInJets"));

}  


void TrueHZParticlesForJets::init() {

  std::cout<<"R/beta/gamma parameter MCparticles Analyzer "<<m_R<<"/"<<m_beta<<"/"<<m_gamma<<std::endl;

  //fastjet::contrib::ValenciaPlugin* vlcpl10;
  //R,beta,gamma
  //-->beta is weight of energy in distance measures, gamma gives weight to 
  vlcpl = new fastjet::contrib::ValenciaPlugin(m_R,m_beta,m_gamma);
  //					      atof(_jetAlgoNameAndParams[1].c_str()),  // R value
					      //atof(_jetAlgoNameAndParams[2].c_str()),  // beta value
					      //atof(_jetAlgoNameAndParams[3].c_str())   // gamma value
					      //);
  //fastjet::JetDefinition* 
  _jetAlgoType = new fastjet::JetDefinition(vlcpl);

}


void TrueHZParticlesForJets::processRunHeader( LCRunHeader*) {

}

void TrueHZParticlesForJets::processEvent( LCEvent* evt ) {
  
  //std::cout<<"we are in event "<<evt->getEventNumber()<<std::endl;

    LCCollection * mcColl =0;
    getCollection(mcColl,m_inputMCParticleCollection,evt);


    LCCollectionVec *trueMCLepPhTauCol = new LCCollectionVec(LCIO::MCPARTICLE);
    LCCollectionVec *trueIsoMCLepPhCol = new LCCollectionVec(LCIO::MCPARTICLE);
    //jet particles removing "true" MC leptons,ISR photons as well as tau decays
    LCCollectionVec *trueMCJetCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    //here all particles bar isolated leptons and photons (E>10 are collected)
    LCCollectionVec *trueMCJetPartNoIsoLepPh = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

    //bool H_decays_ccbar=false;
    //bool H_decays_bbbar=false;

    //bool Z_decays_ccbar=false;
    //bool Z_decays_bbbar=false;

    if( mcColl != 0 ){
      int nMCP = mcColl->getNumberOfElements();
      
      std::set<MCParticle*> boson_daughtersFunc;
      //int ind_MCLep=-1;



      for(int m=0;m<nMCP;m++){
	MCParticle *mcp = static_cast<MCParticle*>(mcColl->getElementAt(m));
	TLorentzVector trueLepPhCand(0,0,0,0);
	trueLepPhCand.SetPxPyPzE(mcp->getMomentum()[0],mcp->getMomentum()[1],mcp->getMomentum()[2],mcp->getEnergy());	
	/*if(m<10){
	  if(mcp->getPDG()==25 && (mcp->getDaughters().size())>0){
	    //std::cout<<"higgs "<<mcp->getPDG()<<" daughter PID "<<mcp->getDaughters()[0]->getPDG()<<" phi/theta "<<trueLepPhCand.Phi()*TMath::RadToDeg()<<"/"<<trueLepPhCand.Theta()*TMath::RadToDeg()<<std::endl;
	    if(abs(mcp->getDaughters()[0]->getPDG())==5){
	      H_decays_bbbar=true;
	    }else if(abs(mcp->getDaughters()[0]->getPDG())==4){
	      H_decays_ccbar=true;
	    }	  
	  }else if(mcp->getPDG()==23 && (mcp->getDaughters().size())>0){
	    //std::cout<<"Z boson "<<mcp->getPDG()<<" daughter PID "<<mcp->getDaughters()[0]->getPDG()<<" phi/theta "<<trueLepPhCand.Phi()*TMath::RadToDeg()<<"/"<<trueLepPhCand.Theta()*TMath::RadToDeg()<<std::endl;
	    if(abs(mcp->getDaughters()[0]->getPDG())==5){
	      Z_decays_bbbar=true;
	    }else if(abs(mcp->getDaughters()[0]->getPDG())==4){
	      Z_decays_ccbar=true;
	    }	  
	  }
	  }*/
	//if(m==8){
	//std::cout<<"should be Z "<<mcp->getPDG()<<" no daugthers "<<mcp->getDaughters().size()<<" d1 "<<mcp->getDaughters()[0]->getEnergy()<<"/"<<mcp->getDaughters()[0]->getPDG()<<" d2 "<<mcp->getDaughters()[1]->getEnergy()<<"/"<<mcp->getDaughters()[1]->getPDG()<<"/"<<trueLepPhCand.Phi()*TMath::RadToDeg()<<"/"<<trueLepPhCand.Theta()*TMath::RadToDeg()<<std::endl;
	//}
	//check for isolation for stable photons and leptons aka electrons and muons
	bool isIsolated=true;
	float iso_absCH=0;
	float iso_absPh=0;
	float iso_absEl=0;
	float iso_absMu=0;
	float iso_absNH=0;
	//find here isolated photons,electrons and muons
	if(mcp->getGeneratorStatus()==1 && mcp->getEnergy()>m_minE && (mcp->getPDG()==22 || abs(mcp->getPDG())==11 || abs(mcp->getPDG())==13)){
	  //charged hadrons

	  for(int m1=0;m1<nMCP;m1++){
	    if(m1!=m && isIsolated){
	      MCParticle *mcp1 = static_cast<MCParticle*>(mcColl->getElementAt(m1));
	      //consider stable particles, exclude neutrinos
	      if(mcp1->getGeneratorStatus()==1 && (fabs(mcp1->getPDG())!=12 && fabs(mcp1->getPDG())!=14 && fabs(mcp1->getPDG())!=16)){
		TLorentzVector temp(0,0,0,0);
		temp.SetPxPyPzE(mcp1->getMomentum()[0],mcp1->getMomentum()[1],mcp1->getMomentum()[2],mcp1->getEnergy());	
		if((trueLepPhCand.Angle(temp.Vect())*TMath::RadToDeg())<m_angleIso){
		  //std::cout<<"part "<<m<<" type/E "<<mcp->getPDG()<<"/"<<mcp->getEnergy()<<" check "<<m1<<" type/E/cone/ cut "<<mcp1->getPDG()<<"/"<<mcp1->getEnergy()<<"/"<<trueLepPhCand.Angle(temp.Vect())*TMath::RadToDeg()<<"/"<<m_angleIso<<std::endl;
		  //fill correspondent isolation value
		  if(mcp1->getCharge()!=0){
		    if(abs(mcp1->getPDG())==11){
		      iso_absEl+=mcp1->getEnergy();
		    }else if(abs(mcp1->getPDG())==13){
		      iso_absMu+=mcp1->getEnergy();
		    }else{
		      iso_absCH+=mcp1->getEnergy();
		    }
		  }else{
		    if(mcp1->getPDG()==22){
		      iso_absPh+=mcp1->getEnergy();
		    }else{
		      iso_absNH+=mcp1->getEnergy();
		    }
		  }
		  if((iso_absCH+iso_absPh+iso_absEl+iso_absMu+iso_absNH)>(mcp->getEnergy()*m_relIso)){
		    isIsolated=false;
		  }

		}
	      }
	    }
	  }
	}else{
	  isIsolated=false;
	}
	if(mcp->getGeneratorStatus()==1 && (abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 && abs(mcp->getPDG())!=16)){
	  if(isIsolated){
	    MCParticleImpl* truePartIso = new MCParticleImpl;
	    truePartIso->setMomentum(mcp->getMomentum());
	    truePartIso->setPDG(mcp->getPDG());
	    truePartIso->setMass(mcp->getMass());
	    truePartIso->setCharge(mcp->getCharge());
	    trueIsoMCLepPhCol->addElement(truePartIso); 
	    //if(((abs(mcp->getPDG())==11 || abs(mcp->getPDG())==13)) && (H_decays_bbbar || H_decays_ccbar )){
	    //std::cout<<" we have an isolated gen lepton in an H->bb or H->cc event "<<mcp->getPDG()<<"/"<<(iso_absCH+iso_absPh+iso_absEl+iso_absMu+iso_absNH)/mcp->getEnergy()<<" "<<mcp->getEnergy() <<"/"<<trueLepPhCand.Phi()*TMath::RadToDeg()<<"/"<<trueLepPhCand.Theta()*TMath::RadToDeg()<<std::endl;
	    //}
	  }else{
	    if(mcp->getEnergy()==0){//seems the ISR photon can be of so low energy that it is E==0, then the jet algorithm is screwed up completely
	      continue;
	    }
	    ReconstructedParticleImpl* truePartNonIso = new ReconstructedParticleImpl;
	    truePartNonIso->setMomentum(mcp->getMomentum());
	    truePartNonIso->setType(mcp->getPDG());
	    truePartNonIso->setEnergy(mcp->getEnergy());
	    truePartNonIso->setMass(mcp->getMass());
	    truePartNonIso->setCharge(mcp->getCharge());
	    trueMCJetPartNoIsoLepPh->addElement(truePartNonIso); 
	  }
	}
	//check for ISR photon
	if(mcp->getPDG()==22 && m<6){
	  fillStableDaughterSet(mcp, boson_daughtersFunc);
	}
	//in order to catch events where the lepton branches off an FSR photon
	//for Higgs check only additional photons, electrons and muons, tau jets etc checked on in main analyzer, if H decays further into Z and W's the code should catch it too
	//
	if((abs(mcp->getPDG())==24 || mcp->getPDG()==23  || mcp->getPDG()==25) && mcp->getDaughters().size()>1 && ((abs(mcp->getDaughters()[0]->getPDG())>6 && abs(mcp->getDaughters()[0]->getPDG())<15) || (mcp->getDaughters()[0]->getPDG()==22 && mcp->getDaughters()[1]->getPDG()<22))){
	  fillStableDaughterSet(mcp, boson_daughtersFunc);
	}
	if(mcp->getGeneratorStatus()!=1){
	  continue;
	}
	if (boson_daughtersFunc.count(mcp) != 0)
	{
	  //save visible particles from resonances, as well as ISR photons
	  if(fabs(mcp->getPDG())!=12 && fabs(mcp->getPDG())!=14 && fabs(mcp->getPDG())!=16){
	    MCParticleImpl* truePart = new MCParticleImpl;
	    truePart->setMomentum(mcp->getMomentum());
	    truePart->setPDG(mcp->getPDG());
	    truePart->setMass(mcp->getMass());
	    truePart->setCharge(mcp->getCharge());
	    trueMCLepPhTauCol->addElement(truePart); 
	  }
	  continue;
	}
	if(fabs(mcp->getPDG())==12 || fabs(mcp->getPDG())==14 || fabs(mcp->getPDG())==16){
	  continue;
	}
	//particles for jets, take out leptons from and their FSR photons from W,Z and H's, also taken out ISR photons of the beam
	//here we fill the MC particle collection, cleaned from leptons and ISR photons from incoming lepton beams, all cleaned based on MC information
	ReconstructedParticleImpl* truePartIntoReco = new ReconstructedParticleImpl;
	truePartIntoReco->setMomentum(mcp->getMomentum());
	truePartIntoReco->setType(mcp->getPDG());
	truePartIntoReco->setEnergy(mcp->getEnergy());
	truePartIntoReco->setMass(mcp->getMass());
	truePartIntoReco->setCharge(mcp->getCharge());
	trueMCJetCol->addElement(truePartIntoReco); 
      }
    }

    LCCollection * recoColl =0;
    getCollection(recoColl,m_inputRecoParticleCollection,evt);
    int nReco = recoColl->getNumberOfElements();
    LCCollectionVec *reccolNoIsoLepPh = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    LCCollectionVec *reccolIsoLepPh = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);


    for(int r=0;r<nReco;r++){
      ReconstructedParticle *reco = static_cast<ReconstructedParticle*>(recoColl->getElementAt(r));
      if(reco->getEnergy()>m_minE && (abs(reco->getType())==11 || abs(reco->getType())==13 || reco->getType()==22)){
	TLorentzVector candLepPh(0,0,0,0);
	candLepPh.SetPxPyPzE(reco->getMomentum()[0],reco->getMomentum()[1],reco->getMomentum()[2],reco->getEnergy());	  
	bool isIsolated=true; 
	float iso_absCH=0;
	float iso_absPh=0;
	float iso_absEl=0;
	float iso_absMu=0;
	float iso_absNH=0;


	for(int r1=0;r1<nReco;r1++){
	  if(r1!=r && isIsolated){
	    ReconstructedParticle *reco1 = static_cast<ReconstructedParticle*>(recoColl->getElementAt(r1));
	    TLorentzVector recoVec(0,0,0,0);
	    recoVec.SetPxPyPzE(reco1->getMomentum()[0],reco1->getMomentum()[1],reco1->getMomentum()[2],reco1->getEnergy());
	    if((TMath::RadToDeg()*recoVec.Angle(candLepPh.Vect()))<m_angleIso){
	      if(reco1->getCharge()!=0){
		if(abs(reco1->getType())==11){
		  iso_absEl+=reco1->getEnergy();
		}else if(abs(reco1->getType())==13){
		  iso_absMu+=reco1->getEnergy();
		}else{
		  iso_absCH+=reco1->getEnergy();
		}
	      }else{
		if(reco1->getType()==22){
		  iso_absPh+=reco1->getEnergy();
		}else{
		  iso_absNH+=reco1->getEnergy();
		}
	      }
	      if((iso_absCH+iso_absPh+iso_absEl+iso_absMu+iso_absNH)>(reco->getEnergy()*m_relIso)){
		isIsolated=false;
	      }
	    }
	  }
	}
	ReconstructedParticleImpl* recoPart= CopyRecoParticle(reco);
	//std::cout<<" recopart energy, type, size of tracks "<<recoPart->getEnergy()<<"/"<<recoPart->getTracks().size()<<std::endl;

	if(isIsolated){
	  reccolIsoLepPh->addElement(recoPart); 
	  //if(((abs(reco->getType())==11 || abs(reco->getType())==13)) && (H_decays_bbbar || H_decays_ccbar)){
	  //std::cout<<" we have an isolated reco lepton in an H->bb or H->cc event "<<reco->getType()<<" "<<(iso_absCH+iso_absPh+iso_absEl+iso_absMu+iso_absNH)/reco->getEnergy()<<"/"<<reco->getEnergy()<<"/"<<candLepPh.Phi()*TMath::RadToDeg()<<"/"<<candLepPh.Theta()*TMath::RadToDeg()<<std::endl;
	  //}
	}else{
	  reccolNoIsoLepPh->addElement(recoPart); 
	}
      }else{
	//here all hadrons and soft leptons/photons will end up, per default these will end up in the non isolation part
	ReconstructedParticleImpl* recoPart = new ReconstructedParticleImpl;
	recoPart->setMomentum(reco->getMomentum());
	recoPart->setEnergy(reco->getEnergy());
	recoPart->setType(reco->getType());
	recoPart->setCovMatrix(reco->getCovMatrix());
	recoPart->setMass(reco->getMass());
	recoPart->setCharge(reco->getCharge());
	recoPart->setParticleIDUsed(reco->getParticleIDUsed());
	recoPart->setGoodnessOfPID(reco->getGoodnessOfPID());
	recoPart->setStartVertex(reco->getStartVertex());
	for(unsigned int i=0;i<reco->getTracks().size();i++){
	  recoPart->addTrack(reco->getTracks()[i]);
	}
	for(unsigned int i=0;i<reco->getClusters().size();i++){
	  recoPart->addCluster(reco->getClusters()[i]);
	}
	reccolNoIsoLepPh->addElement(recoPart); 
      }
    }

    LCCollectionVec *reccolPFOInJets = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

    LCCollectionVec *reccol_jet0 = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    LCCollectionVec *reccol_jet1 = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);


    //at this point we should have saved the particles in a collection
    if(reccolNoIsoLepPh->getNumberOfElements()>1){
      PseudoJetList recopjList;
      for(int i = 0; i < reccolNoIsoLepPh->getNumberOfElements(); ++i){
	ReconstructedParticle* par = static_cast<ReconstructedParticle*> (reccolNoIsoLepPh->getElementAt(i));
	if(par->getEnergy()==0){
	  std::cout<<"input particle "<<i<<" E/px/py/pz/type "<<par->getEnergy()<<"/"<<par->getMomentum()[0]<<"/"<<par->getMomentum()[1]<<"/"<<par->getMomentum()[2]<<"/"<<par->getType()<<std::endl;
	}
	recopjList.push_back( fastjet::PseudoJet( par->getMomentum()[0],
						  par->getMomentum()[1],
						  par->getMomentum()[2],
						  par->getEnergy() ) );
	recopjList.back().set_user_index(i);	// save the id of this recParticle
      }

      fastjet::ClusterSequence*_csreco = new fastjet::ClusterSequence(recopjList, *_jetAlgoType);
      PseudoJetList recojets = _csreco->exclusive_jets(m_NJets);

      unsigned int index_jet=0;

      for (std::vector<fastjet::PseudoJet>::iterator recojetIt = recojets.begin(); recojetIt != recojets.end(); ++recojetIt ) {

	ReconstructedParticleImpl* reco = new ReconstructedParticleImpl();
	
	// save the jet's parameters                                                                                                                                                                                                                                       
	reco->setEnergy( recojetIt->E() );
	reco->setMass( recojetIt->m() );
	
	double mom[3] = {recojetIt->px(), recojetIt->py(), recojetIt->pz()};
	reco->setMomentum( mom );
	

	//std::cout<<"HZTrueParticles cluster "<<reccolNoIsoLepPh->getNumberOfElements()<<" size of jet "<<_csreco->constituents((*recojetIt)).size()<<std::endl;
	
	for(unsigned int i=0;i< _csreco->constituents((*recojetIt)).size();i++){
	  //recojet-particles are prepared as reconstructed particle
	  //reco Type is the original PDGID, charge is filled as well
	  ReconstructedParticle* recopart = static_cast<ReconstructedParticle*> (reccolNoIsoLepPh->getElementAt(_csreco->constituents((*recojetIt))[i].user_index()));
	  reco->addParticle( recopart );
	  ReconstructedParticleImpl* recop_impl=CopyRecoParticle (recopart );
	  //second copy needed, for whatever reason programme crashes otherwise
	  ReconstructedParticleImpl* recop_impl_1=CopyRecoParticle (recopart );
	  reccolPFOInJets->addElement(recop_impl_1);
	  if(index_jet==0){
	    reccol_jet0->addElement(recop_impl);
	  }else if (index_jet==1){
	    reccol_jet1->addElement(recop_impl);
	  }
	}
	//iterator over pseudojet list
	index_jet+=1;
      }   
    }
 
    //if(H_decays_bbbar || H_decays_ccbar || Z_decays_bbbar || Z_decays_ccbar){
    //std::cout<<"event has HF decay in H or Z "<<H_decays_bbbar<<"/"<<H_decays_ccbar<<"/"<<Z_decays_bbbar<<"/"<<Z_decays_ccbar<<std::endl;
    //}
    
    
    evt->addCollection(reccol_jet0,m_outputPFOJet0Collection);
    evt->addCollection(reccol_jet1,m_outputPFOJet1Collection);


    evt->addCollection(trueIsoMCLepPhCol,m_outputMCIsoLepPhParticleCollection);
    evt->addCollection(trueMCJetPartNoIsoLepPh,m_outputMCForJetParticleCollection);

    evt->addCollection(trueMCLepPhTauCol,m_outputMCTrueLepPhParticleCollection);
    evt->addCollection(trueMCJetCol,m_outputMCTrueJetParticleCollection);

    evt->addCollection(reccolNoIsoLepPh,m_outputRecoJetParticleCollection);
    evt->addCollection(reccolIsoLepPh,m_outputRecoIsoLepPhParticleCollection);
 
    evt->addCollection(reccolPFOInJets,m_outputPFOsInJetsCollection);

 

}

void TrueHZParticlesForJets::fillStableDaughterSet(MCParticle* mcp, std::set<MCParticle*> &stableDaughterSet){
  if(mcp->getGeneratorStatus()==1){
    stableDaughterSet.insert(mcp);
  }else if (mcp->getGeneratorStatus()==0){
    return;
  }
  for(unsigned int d=0;d<mcp->getDaughters().size();d++){
    fillStableDaughterSet(mcp->getDaughters()[d], stableDaughterSet);
  }
}

void TrueHZParticlesForJets::check(LCEvent*){
}

void TrueHZParticlesForJets::end(){

  delete vlcpl;
  delete _jetAlgoType;

}


void TrueHZParticlesForJets::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    
    streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

ReconstructedParticleImpl* TrueHZParticlesForJets::CopyRecoParticle ( ReconstructedParticle* pfo_orig ) {
	// copy this in an ugly fashion to be modifiable - a versatile copy constructor would be much better!
	ReconstructedParticleImpl* pfo = new ReconstructedParticleImpl();
	pfo->setMomentum(pfo_orig->getMomentum());
	pfo->setEnergy(pfo_orig->getEnergy());
	pfo->setType(pfo_orig->getType());
	pfo->setCovMatrix(pfo_orig->getCovMatrix());
	pfo->setMass(pfo_orig->getMass());
	pfo->setCharge(pfo_orig->getCharge());
	pfo->setParticleIDUsed(pfo_orig->getParticleIDUsed());
	pfo->setGoodnessOfPID(pfo_orig->getGoodnessOfPID());
	pfo->setStartVertex(pfo_orig->getStartVertex());
	for(unsigned int i=0;i<pfo_orig->getTracks().size();i++){
	  pfo->addTrack(pfo_orig->getTracks()[i]);
	}
	for(unsigned int i=0;i<pfo_orig->getClusters().size();i++){
	  pfo->addCluster(pfo_orig->getClusters()[i]);
	}
	return pfo;
}
