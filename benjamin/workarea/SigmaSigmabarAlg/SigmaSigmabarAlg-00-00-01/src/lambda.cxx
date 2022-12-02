//author Liu Liang
//version 0.3
//add MCtruth four momentum of n nbar pi- pi+
//compare the sum of this four particles with J/psi
//check the kinematic fit
//set a new kinematicfit, use emc/2 -Epi as nbar input parameter 
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/IHistogramSvc.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "VertexFit/SecondVertexFit.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "DstEvent/TofHitStatus.h"

#include "McTruth/McParticle.h"
#include "Identifier/TofID.h"
#include "ParticleID/ParticleID.h"
#include "TrigEvent/TrigEvent.h"
#include "TrigEvent/TrigData.h"

#include "TMath.h"
#include "TF1.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"

#include "LambdaAlg/lambda.h"
#include "TString.h"


#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McDecayModeSvc/IMcDecayModeSvc.h"

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

#ifndef ENABLE_BACKWARDS_COMPATIBITY
//#include "MCDataCorrectionAlg/mcdatacorrection.h"
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include <vector>

const double mpi = 0.13957;
const double mkaon = 0.493677;
const double mproton = 0.938272;
const double mpion0 = 0.1349766;
const double mpion = 0.139570;
const double mEta = 0.547853;
const double n_mass = 0.939565;
const double nbar_mass = 0.939565;
const double Sigmam_mass = 1.197436;
const double Sigma_mass = 1.18937;
const double lambda_mass = 1.115683;
const double xmass[5] = {0.000511, 0.105658, 0.139570, 0.493677, 0.938272};
const double velc = 299.792458;			// tof path unit in mm
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<WTrackParameter> VWTP;

static int Ncut_none=0, 	Ncut_ncharge=0, 
		   Ncut_pid=0, 		Ncut_nshower=0, 
		   Ncut_npi0=0,		Ncut_photon=0;
static int Ncut_pppi0pi0=0,	Ncut_ppi0=0,
		   Ncut_pbarpi0=0,	Ncut_nnbarpipi,
		   Ncut_mc = 0;
static int Ncut_kin=0,		Ncut_nbar=0,
		   Ncut_asigma=0, 	Ncut_kin_flag1=0,
		   Ncut_kin_flag2=0, Ncut_vertex = 0;
static int Ncut_mass=0,	Ncut_length=0,
		   Ncut_chisq=0,    Liu = 0;

//------------------------------------------------------------------------------------------------
Lambda::Lambda(const std::string& name, ISvcLocator* pSvcLocator):
		Algorithm(name, pSvcLocator){
				declareProperty("energyTreshold",m_energyLambda=0.025);
				declareProperty("energyTreshold1",m_energyLambda1=0.025);
				declareProperty("energyTreshold2",m_energyLambda2=0.050);
				declareProperty("gammaAngleCut",m_gammaAngleCut=20.0);
				declareProperty("rdmSeed", m_anglename=2);
				declareProperty("cms",cms=3.097);
				declareProperty("WriteMyDst",       m_writeMyDst=false);
		}

//-------------------------------------------------------------------------------


//MCDataCorrection *nbarcorr;
//MCDataCorrection *nbarcorr1;


StatusCode Lambda::initialize(){
		MsgStream log(msgSvc(), name());

		log << MSG::INFO << "in initialize()" << endmsg;

		StatusCode status;

		// event filter set up

		TString Mother  = "EventWriter";

		if(m_writeMyDst)
		{
				status =  createSubAlgorithm( "EventWriter", "WriteMyDst", m_subAlgA);
				if( status.isFailure() )
				{
						log << MSG::ERROR << "  Error creating Sub-Algorithm WriteDst" <<endreq;
						return status;
				}
		}




		NTuplePtr nt0(ntupleSvc(), "FILE1/mctruth");
		if ( nt0 ) m_tuple0 = nt0;
		else {
				m_tuple0 = ntupleSvc()->book ("FILE1/mctruth", CLID_ColumnWiseTuple, "ks N-Tuple example");
				if (m_tuple0 ){
						status = m_tuple0->addItem ("mclambdacos",m_mclambdacos);
						status = m_tuple0->addItem ("mclambdabarcos",m_mclambdabarcos);
						status = m_tuple0->addItem ("mclambdapx",m_mclambdapx);
						status = m_tuple0->addItem ("mclambdapy",m_mclambdapy);
						status = m_tuple0->addItem ("mclambdapz",m_mclambdapz);
						status = m_tuple0->addItem ("mclambdae",m_mclambdae);
						status = m_tuple0->addItem ("mclambdabarpx",m_mclambdabarpx);
						status = m_tuple0->addItem ("mclambdabarpy",m_mclambdabarpy);
						status = m_tuple0->addItem ("mclambdabarpz",m_mclambdabarpz);
						status = m_tuple0->addItem ("mclambdabare",m_mclambdabare);
						status = m_tuple0->addItem ("mcpxpim", m_mcpxpim);
						status = m_tuple0->addItem ("mcpypim", m_mcpypim);
						status = m_tuple0->addItem ("mcpzpim", m_mcpzpim);
						status = m_tuple0->addItem ("mcepim", m_mcepim);
						status = m_tuple0->addItem ("mcpxpip", m_mcpxpip);
						status = m_tuple0->addItem ("mcpypip", m_mcpypip);
						status = m_tuple0->addItem ("mcpzpip", m_mcpzpip);
						status = m_tuple0->addItem ("mcepip", m_mcepip);
						status = m_tuple0->addItem ("mcpxpp", m_mcpxpp);
						status = m_tuple0->addItem ("mcpypp", m_mcpypp);
						status = m_tuple0->addItem ("mcpzpp", m_mcpzpp);
						status = m_tuple0->addItem ("mcepp", m_mcepp);
						status = m_tuple0->addItem ("mcpxpm", m_mcpxpm);
						status = m_tuple0->addItem ("mcpypm", m_mcpypm);
						status = m_tuple0->addItem ("mcpzpm", m_mcpzpm);
						status = m_tuple0->addItem ("mcepm", m_mcepm);

						status  = m_tuple0->addItem ("nmc_pp", m_nmc_pp);
						status  = m_tuple0->addItem ("nmc_pm", m_nmc_pm);
						status  = m_tuple0->addItem ("nmc_pip", m_nmc_pip);
						status  = m_tuple0->addItem ("nmc_pim", m_nmc_pim);
						status  = m_tuple0->addItem ("nmc_lambda", m_nmc_lambda);
						status  = m_tuple0->addItem ("nmc_lambdabar", m_nmc_lambdabar);
				}
				else{
						log << MSG::ERROR << "Cannot book N-Tuple:" << long(m_tuple0) << endmsg;
				}
		}


		NTuplePtr nt1(ntupleSvc(), "FILE1/FIT");
		if ( nt1 ) m_tuple1 = nt1;
		else	{
				m_tuple1 = ntupleSvc()->book ("FILE1/FIT", CLID_ColumnWiseTuple, "ks N-Tuple example");
				if( m_tuple1 )	{
						status  = m_tuple1->addItem ("indexmc", m_idxmc, 0, 100 );
						status  = m_tuple1->addIndexedItem ("trkidx", m_idxmc, m_trkidx);
						status  = m_tuple1->addIndexedItem ("pdgid", m_idxmc, m_pdgid);
						status  = m_tuple1->addIndexedItem ("motheridx", m_idxmc, m_motheridx);
						status  = m_tuple1->addItem ("nGood", m_nGood);
						status  = m_tuple1->addItem ("nCharge", m_nCharge);

						status  = m_tuple1->addItem ("runNo", m_runNo);
						status  = m_tuple1->addItem ("event", m_event);


						status = m_tuple1->addItem ("npip", m_npip, 0, 100);
						status = m_tuple1->addIndexedItem ("pip_p", m_npip, m_pip_p);
						status = m_tuple1->addIndexedItem ("pip_pt", m_npip, m_pip_pt);
						status = m_tuple1->addIndexedItem ("pip_cos", m_npip, m_pip_cos);
						status = m_tuple1->addIndexedItem ("pipRvxy0", m_npip, m_pipRvxy0);
						status = m_tuple1->addIndexedItem ("pipRvz0", m_npip, m_pipRvz0);
						status = m_tuple1->addIndexedItem ("pip_match", m_npip, m_pip_match);
						status = m_tuple1->addIndexedItem ("pippid_flag", m_npip, m_pippid_flag);

						status = m_tuple1->addItem ("npim", m_npim, 0, 100);
						status = m_tuple1->addIndexedItem ("pim_p", m_npim, m_pim_p);
						status = m_tuple1->addIndexedItem ("pim_pt", m_npim, m_pim_pt);
						status = m_tuple1->addIndexedItem ("pim_cos", m_npim, m_pim_cos);
						status = m_tuple1->addIndexedItem ("pimRvxy0", m_npim, m_pimRvxy0);
						status = m_tuple1->addIndexedItem ("pimRvz0", m_npim, m_pimRvz0);
						status = m_tuple1->addIndexedItem ("pim_match", m_npim, m_pim_match);
						status = m_tuple1->addIndexedItem ("pimpid_flag", m_npim, m_pimpid_flag);

						status = m_tuple1->addItem ("npm", m_npm, 0, 100);
						status = m_tuple1->addIndexedItem ("pm_p", m_npm,  m_pm_p);
						status = m_tuple1->addIndexedItem ("pm_pt", m_npm,  m_pm_pt);
						status = m_tuple1->addIndexedItem ("pm_cos", m_npm, m_pm_cos);
						status = m_tuple1->addIndexedItem ("pmRvxy0", m_npm, m_pmRvxy0);
						status = m_tuple1->addIndexedItem ("pmRvz0", m_npm, m_pmRvz0);
						status = m_tuple1->addIndexedItem ("pm_match", m_npm, m_pm_match);
						status = m_tuple1->addIndexedItem ("pmpid_flag", m_npm, m_pmpid_flag);

						status = m_tuple1->addItem ("npp", m_npp, 0, 100);
						status = m_tuple1->addIndexedItem ("pp_p", m_npp, m_pp_p);
						status = m_tuple1->addIndexedItem ("pp_pt", m_npp, m_pp_pt);
						status = m_tuple1->addIndexedItem ("pp_cos", m_npp, m_pp_cos);
						status = m_tuple1->addIndexedItem ("ppRvxy0", m_npp, m_ppRvxy0);
						status = m_tuple1->addIndexedItem ("ppRvz0", m_npp, m_ppRvz0);
						status = m_tuple1->addIndexedItem ("pp_match", m_npp, m_pp_match);
						status = m_tuple1->addIndexedItem ("pppid_flag", m_npp, m_pppid_flag);

						status = m_tuple1->addItem ("lambda_mass", m_lambda_mass);
						status = m_tuple1->addItem ("lambda_cos", m_lambda_cos);
						status = m_tuple1->addItem ("decayL_lambda", m_decayL_lambda);
						status = m_tuple1->addItem ("decayLerr_lambda", m_decayLerr_lambda);
						status = m_tuple1->addItem ("chisq_lambda", m_chisq_lambda);
						status = m_tuple1->addItem ("lambdabar_mass", m_lambdabar_mass);
						status = m_tuple1->addItem ("lambdabar_cos", m_lambdabar_cos);
						status = m_tuple1->addItem ("decayL_lambdabar", m_decayL_lambdabar);
						status = m_tuple1->addItem ("decayLerr_lambdabar", m_decayLerr_lambdabar);
						status = m_tuple1->addItem ("chisq_lambdabar", m_chisq_lambdabar);

						status = m_tuple1->addItem ("chisq", m_chisq);
						status = m_tuple1->addItem ("jpsi_mass", m_jpsi_mass);
						status = m_tuple1->addItem ("p4lambda", 4, m_p4lambda);
						status = m_tuple1->addItem ("p4lambdabar", 4, m_p4lambdabar);

						status = m_tuple1->addItem ("p4pip", 4, m_p4pip);
						status = m_tuple1->addItem ("p4ap", 4, m_p4ap);
						status = m_tuple1->addItem ("p4pim", 4, m_p4pim);
						status = m_tuple1->addItem ("p4p", 4, m_p4p);
						status = m_tuple1->addItem ("npppid", m_npppid);
						status = m_tuple1->addItem ("npmpid", m_npmpid);
						status = m_tuple1->addItem ("npippid", m_npippid);
						status = m_tuple1->addItem ("npimpid", m_npimpid);
				}
				else {
						log << MSG::ERROR << "Connot book N-tuple : " << long (m_tuple1) << endmsg;
						return StatusCode::FAILURE;
				}
		}


		NTuplePtr nt2(ntupleSvc(), "FILE1/kmf");
		if ( nt2 ) m_tuple2 = nt2;
		else	{
				m_tuple2 = ntupleSvc()->book ("FILE1/kmf", CLID_ColumnWiseTuple, "ks N-Tuple example");
				if( m_tuple2 )	{

						status = m_tuple2->addItem ("pxpim", m_pxpim);
						status = m_tuple2->addItem ("pypim", m_pypim);
						status = m_tuple2->addItem ("pzpim", m_pzpim);
						status = m_tuple2->addItem ("epim", m_epim);
						status = m_tuple2->addItem ("pxpip", m_pxpip);
						status = m_tuple2->addItem ("pypip", m_pypip);
						status = m_tuple2->addItem ("pzpip", m_pzpip);
						status = m_tuple2->addItem ("epip", m_epip);
						status = m_tuple2->addItem ("pxnbar", m_pxnbar);
						status = m_tuple2->addItem ("pynbar", m_pynbar);
						status = m_tuple2->addItem ("pznbar", m_pznbar);
						status = m_tuple2->addItem ("enbar", m_enbar);
						status = m_tuple2->addItem ("pxn", m_pxn);
						status = m_tuple2->addItem ("pyn", m_pyn);
						status = m_tuple2->addItem ("pzn", m_pzn);
						status = m_tuple2->addItem ("en", m_en);
						status = m_tuple2->addItem ("sigma_cos", m_sigma_cos);
						status = m_tuple2->addItem ("sigma_mass", m_sigma_mass);
						status = m_tuple2->addItem ("nshower", m_nshower);
						status = m_tuple2->addItem ("niter", m_niter);
						status = m_tuple2->addItem ("chisq_kmf", m_chisq_kmf);

						status = m_tuple2->addItem ("lamkmf_mass", m_lamkmf_mass);
						status = m_tuple2->addItem ("lamkmf_cos", m_lamkmf_cos);
						status = m_tuple2->addItem ("decayL_lamkmf", m_decayL_lamkmf);
						status = m_tuple2->addItem ("decayLerr_lamkmf", m_decayLerr_lamkmf);
						status = m_tuple2->addItem ("chisq_lamkmf", m_chisq_lamkmf);
						status = m_tuple2->addItem ("lamkmfbar_mass", m_lamkmfbar_mass);
						status = m_tuple2->addItem ("lamkmfbar_cos", m_lamkmfbar_cos);
						status = m_tuple2->addItem ("decayL_lamkmfbar", m_decayL_lamkmfbar);
						status = m_tuple2->addItem ("decayLerr_lamkmfbar", m_decayLerr_lamkmfbar);
						status = m_tuple2->addItem ("chisq_lamkmfbar", m_chisq_lamkmfbar);





				}
				else {
						log << MSG::ERROR << "Connot book N-tuple : " << long (m_tuple2) << endmsg;
						return StatusCode::FAILURE;
				}
		}





		hcos_total = new TH1D("hcos_total", "hcos_total", 20, -1, 1);
		hcos_charged = new TH1D("hcos_charged", "hcos_charged", 20, -1, 1);
		hcos_vertex = new TH1D("hcos_vertex", "hcos_vertex", 20, -1, 1);
		hcos_kin = new TH1D("hcos_kin", "hcos_kin", 20, -1, 1);



		log << MSG::INFO << "successfully return from initialize()" <<endmsg;
		return StatusCode::SUCCESS;
}



//========================================================================================



StatusCode Lambda::execute()	{
		MsgStream log(msgSvc(),name());
		log << MSG::INFO << "in execute()" << endreq;
		SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
		int runNo=eventHeader->runNumber();
		int event=eventHeader->eventNumber();
		m_runNo = eventHeader->runNumber();
		m_event = eventHeader->eventNumber();
		log << MSG::DEBUG << "run, evtnum = "
				<< runNo << " , "
				<< event << endreq;
		Ncut_none++;
		//	cout << "**********Ncut_none :" << Ncut_none <<endl; 
		/*
		   int RadTrig = 0;
		   if ( runNo>0 )	{
		   SmartDataPtr<TrigData> trigData(eventSvc(),EventModel::Trig::TrigData);
		   if( !trigData )	{
		   cout << "Could not find Trigger Data for physics analysis" << endl;
		   return StatusCode::FAILURE;
		   }
		   RadTrig = trigData->getTrigChannel(9);
		   }
		   */


		IMcDecayModeSvc* i_svc;
		StatusCode sc_DecayModeSvc = service("McDecayModeSvc", i_svc);
		if ( sc_DecayModeSvc.isFailure() ){
				log << MSG::FATAL << "Could not load McDecayModeSvc!" << endreq;
				return sc_DecayModeSvc;
		}
		m_svc = dynamic_cast<McDecayModeSvc*>(i_svc);

		if(runNo < 0) {
				SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
				if(!mcParticleCol){
						std::cout << "Could not retrieve McParticleCol" << std::endl;
						return StatusCode::FAILURE;
				}
				else{
						std::vector<int> pdgid;
						std::vector<int> motheridx;
						pdgid.clear();
						motheridx.clear();
						bool JpsiDecay = false; 
						int rootIndex = -1;
						int m_numParticle = 0;
						Event::McParticleCol::iterator iter_mc=mcParticleCol->begin();
						for(;iter_mc!=mcParticleCol->end();iter_mc++){
								if ((*iter_mc)->primaryParticle()) continue;
								if(!(*iter_mc)->decayFromGenerator()) continue;
								if ((*iter_mc)->particleProperty()==443){
										int mode = m_svc->extract(*iter_mc, pdgid, motheridx);
										m_numParticle = pdgid.size();
										for (int i=0; i < pdgid.size(); i++){
												m_pdgid[i] = pdgid[i];
												m_motheridx[i] = motheridx[i];
										}
										m_idxmc = m_numParticle;
								}
						}
				}
		}





		HepLorentzVector mc_p4pip, mc_p4pim, mc_p4pp, mc_p4pm, mc_p4lambda, mc_p4lambdabar, mc_p4pi0, mc_p4n;
		if( runNo < 0 ) {
				SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
				if( mcParticleCol ){
						bool vphoDecay = false;
						int rootIndex = -1;
						int nmc_pip = 0, nmc_pim = 0, nmc_nbar = 0, nmc_pp = 0, nmc_lambda = 0, nmc_pi0 = 0, nmc_pm = 0, nmc_lambdabar = 0;
						Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
						for (; iter_mc != mcParticleCol->end(); iter_mc++)	{
								if (!(*iter_mc)->decayFromGenerator()) continue;
								vphoDecay = true;
								HepLorentzVector mctrue_track = (*iter_mc)->initialFourMomentum();
								if((*iter_mc)->particleProperty()==-3122){ // anti-sigma+  -3112
										mc_p4lambdabar = mctrue_track;
										m_mclambdabarcos = cos(mc_p4lambdabar.theta());
										m_mclambdabare = mc_p4lambdabar.e();
										m_mclambdabarpx = mc_p4lambdabar.px();
										m_mclambdabarpy = mc_p4lambdabar.py();
										m_mclambdabarpz = mc_p4lambdabar.pz();

										nmc_lambdabar++;
								}
								if((*iter_mc)->particleProperty()==-2212 && ((*iter_mc)->mother()).particleProperty()==-3122){  //anti-n0 -2112
										mc_p4pm = mctrue_track;
										m_mcpxpm = mc_p4pm.px(); 	m_mcpypm = mc_p4pm.py();	m_mcpzpm = mc_p4pm.pz();	m_mcepm = mc_p4pm.e();
										nmc_pm++;
								}
								if ((*iter_mc)->particleProperty()==2212 ) {
										mc_p4pp = mctrue_track;
										m_mcpxpp = mc_p4pp.px(); 	m_mcpypp = mc_p4pp.py();	m_mcpzpp = mc_p4pp.pz();	m_mcepp = mc_p4pp.e();
										nmc_pp++;
								}
								if((*iter_mc)->particleProperty()==211 )	{
										mc_p4pip = mctrue_track;
										m_mcpxpip = mc_p4pip.px(); 	m_mcpypip = mc_p4pip.py();	m_mcpzpip = mc_p4pip.pz();	m_mcepip = mc_p4pip.e();
										nmc_pip++;
								}
								if((*iter_mc)->particleProperty()==-211) {
										mc_p4pim = mctrue_track;
										m_mcpxpim = mc_p4pim.px(); 	m_mcpypim = mc_p4pim.py();	m_mcpzpim = mc_p4pim.pz();	m_mcepim = mc_p4pim.e();
										nmc_pim++;
								} 
								if((*iter_mc)->particleProperty()==111 && ((*iter_mc)->mother()).particleProperty()==3122) {
										mc_p4pi0 = mctrue_track;
								}
								if((*iter_mc)->particleProperty()==2112 && ((*iter_mc)->mother()).particleProperty()==3122) {
										mc_p4n = mctrue_track;
								}

								if ((*iter_mc)->particleProperty()==3122){
										mc_p4lambda = mctrue_track;
										//m_mclambda_p = mctrue_track.rho();
										m_mclambdae = mctrue_track.e();
										m_mclambdapx = mctrue_track.px();
										m_mclambdapy = mctrue_track.py();
										m_mclambdapz = mctrue_track.pz();
										m_mclambdacos = cos(mctrue_track.theta());
										nmc_lambda++;
								}

						}
						m_nmc_pp = nmc_pp;
						m_nmc_pm = nmc_pm;
						m_nmc_pip = nmc_pip;
						m_nmc_pim = nmc_pim;
						m_nmc_lambda = nmc_lambda;
						m_nmc_lambdabar = nmc_lambdabar; 
						m_tuple0->write();
						Ncut_mc++;


				}

		}


		SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
		log 	<< MSG::DEBUG 
				<< " nCharge  = " << evtRecEvent->totalCharged() << " , "
				<< " nNeutral = " << evtRecEvent->totalNeutral() << " , "
				<< " tottrack = " << evtRecEvent->totalTracks()  << endreq;
		/*	cout 	<< " nCharge  = " << evtRecEvent->totalCharged() << " , "
			<< " nneutral = " << evtRecEvent->totalNeutral() << " , "
			<< " tottrack = " << evtRecEvent->totalTracks()  << endl; */
		SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
		//
		//	check x0, y0, z0, r0;
		//	suggest cut: |z0|<5 && r0<1
		//
		Vint iGood, ipip,ipim, iKp, iKm, ipp, ipm;
		iGood.clear();
		ipip.clear();
		ipim.clear();
		iKp.clear();
		iKm.clear();
		ipp.clear();
		ipm.clear();
		Vint ipppid, ipmpid, ipippid, ipimpid;
		ipppid.clear();
		ipmpid.clear();
		ipippid.clear();
		ipimpid.clear();

		Vp4 ppip, ppim;
		ppip.clear();
		ppim.clear();

		int nCharge = 0;

		hcos_total->Fill(mc_p4lambda.cosTheta());

		Hep3Vector xorigin(0,0,0);
		HepPoint3D vx(0., 0., 0.);
		HepSymMatrix Evx(3, 0);
		IVertexDbSvc* vtxsvc;
		Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
		if(vtxsvc->isVertexValid()){
				double* dbv = vtxsvc->PrimaryVertex();
				double*  vv = vtxsvc->SigmaPrimaryVertex();
				xorigin.setX(dbv[0]);
				xorigin.setY(dbv[1]);
				xorigin.setZ(dbv[2]);
				vx.setX(dbv[0]);
				vx.setY(dbv[1]);
				vx.setZ(dbv[2]);
				Evx[0][0]=vv[0]*vv[0];
				Evx[1][1]=vv[1]*vv[1];
				Evx[2][2]=vv[2]*vv[2];
		}
		VertexParameter vx_db;
		vx_db.setVx(vx);
		vx_db.setEvx(Evx);

		int idx_pp=0, idx_pm=0, idx_pip=0, idx_pim=0;
		for ( int i = 0; i < evtRecEvent->totalCharged(); i++){
				EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
				if (!(*itTrk)->isMdcTrackValid()) continue;
				RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
				double pch  = mdcTrk->p();
				double x0	= mdcTrk->x();
				double y0 	= mdcTrk->y();
				double z0   = mdcTrk->z();
				double phi0 = mdcTrk->helix(1);
				double xv   = xorigin.x();
				double yv   = xorigin.y();
				double Rzy  = (x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);



				HepVector a = mdcTrk->helix();
				HepSymMatrix Ea = mdcTrk->err();
				HepPoint3D point0(0.,0.,0.);	// the initial point for MDC recosntruction
				HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
				VFHelix helixip(point0, a, Ea);
				helixip.pivot(IP);
				HepVector vecipa = helixip.a();
				double Rvxy0 = fabs(vecipa[0]);	// the nearest distance to IP in xy plane
				double Rvz0 = vecipa[3];
				double Rvphi0 = vecipa[1];
				double cost = cos(mdcTrk->theta());

				//		if(pch > 2.0) continue;

				if(fabs(Rvz0) >= 30.0 ) continue;
				if(fabs(Rvxy0) >= 10.0 ) continue;
				if(fabs(cost)>=0.93) continue;
				RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();

				if (!(*itTrk)->isMdcKalTrackValid()) continue;
				iGood.push_back(i);
				if(mdcTrk->charge() > 0 && mdcTrk->p() > 0.5){
						ipp.push_back(i);
						m_pp_pt[idx_pp]  = mdcTrk->pxy();
						m_pp_p[idx_pp]  = mdcTrk->p();
						m_pp_cos[idx_pp] = cos(mdcTrk->theta());
						m_ppRvxy0[idx_pp] = Rvxy0;
						m_ppRvz0[idx_pp] = Rvz0;
						m_pp_match[idx_pp] = mc_p4pp.angle(mdcTrk->p3());
						idx_pp++;
				}
				if(mdcTrk->charge() > 0 && mdcTrk->p() < 0.5){
						ipip.push_back(i);
						m_pip_pt[idx_pip]  = mdcTrk->pxy();
						m_pip_p[idx_pip]  = mdcTrk->p();
						m_pip_cos[idx_pip] = cos(mdcTrk->theta());
						m_pipRvxy0[idx_pip] = Rvxy0;
						m_pipRvz0[idx_pip] = Rvz0;
						m_pip_match[idx_pip] = mc_p4pip.angle(mdcTrk->p3());
						idx_pip++;
				}
				if(mdcTrk->charge() < 0 && mdcTrk->p() > 0.5){
						ipm.push_back(i);
						m_pm_pt[idx_pm]  = mdcTrk->pxy();
						m_pm_p[idx_pm]  = mdcTrk->p();
						m_pm_cos[idx_pm] = cos(mdcTrk->theta());
						m_pmRvxy0[idx_pm] = Rvxy0;
						m_pmRvz0[idx_pm] = Rvz0;
						m_pm_match[idx_pm] = mc_p4pm.angle(mdcTrk->p3());
						idx_pm++;
				}
				if(mdcTrk->charge() < 0 && mdcTrk->p() < 0.5){
						ipim.push_back(i);
						m_pim_pt[idx_pim]  = mdcTrk->pxy();
						m_pim_p[idx_pim]  = mdcTrk->p();
						m_pim_cos[idx_pim] = cos(mdcTrk->theta());
						m_pimRvxy0[idx_pim] = Rvxy0;
						m_pimRvz0[idx_pim] = Rvz0;
						m_pim_match[idx_pim] = mc_p4pim.angle(mdcTrk->p3());
						idx_pim++;
				}
				nCharge += mdcTrk->charge();

		}
		int nGood = iGood.size();
		m_nGood = nGood;
		m_nCharge = nCharge;

		int npip = ipip.size();
		int npim = ipim.size();
		int npp = ipp.size();
		int npm = ipm.size();
		m_npip = npip;
		m_npim = npim;
		m_npp = ipp.size();
		m_npm = ipm.size();

		if(nGood != 4) return SUCCESS;
		if(npip < 1) return SUCCESS;
		if(npm < 1) return SUCCESS;
		if(npim < 1) return SUCCESS;
		if(npp < 1) return SUCCESS;
		Ncut_ncharge++;
		hcos_charged->Fill(mc_p4lambda.cosTheta());

		//////////////////////////////////////////////////////
		//
		//				PID 
		//
		////////////////////////////////////////////////////////


		ParticleID *pid = ParticleID::instance();

		for(int i = 0; i < npm; i++) {
				EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ipm[i];
				RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
				pid->init();
				pid->setMethod(pid->methodProbability());
				pid->setChiMinCut(4);
				pid->setRecTrack(*itTrk);
	//			pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()); // use PID sub-system
				pid->usePidSys(pid->useDedx() | pid->useTofCorr()); // use PID sub-system
				pid->identify(pid->onlyPionKaonProton());
				pid->calculate();
				if(!(pid->IsPidInfoValid())) continue;
				//	if(!(*itTrk)->isMdcDedxValid())continue;
				if((pid->probProton() > pid->probPion()) && (pid->probProton() > pid->probKaon()) ) {
						ipmpid.push_back(ipm[i]);
				}
		}

		for(int i = 0; i < npp; i++) {
				EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ipp[i];
				RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
				pid->init();
				pid->setMethod(pid->methodProbability());
				pid->setChiMinCut(4);
				pid->setRecTrack(*itTrk);
		//		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()); // use PID sub-system
				pid->usePidSys(pid->useDedx() | pid->useTofCorr()); // use PID sub-system
				pid->identify(pid->onlyPionKaonProton());
				pid->calculate();
				if(!(pid->IsPidInfoValid())) continue;
				//	if(!(*itTrk)->isMdcDedxValid())continue;
				if((pid->probProton() > pid->probPion()) && (pid->probProton() > pid->probKaon()) ) {
						ipppid.push_back(ipp[i]);
				}
		}


		for(int i = 0; i < npip; i++) {
				EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ipip[i];
				//	if(!(*itTrk)->isMdcDedxValid())continue;
				//	if (!(*itTrk)->isMdcTrackValid()) continue;

				RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
				pid->init();
				pid->setMethod(pid->methodProbability());
				pid->setChiMinCut(4);
				pid->setRecTrack(*itTrk);
			//	pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()); // use PID sub-system
				pid->usePidSys(pid->useDedx() | pid->useTofCorr()); // use PID sub-system
				pid->identify(pid->onlyPionKaonProton());
				pid->calculate();
				if(!(pid->IsPidInfoValid())) continue;
				if((pid->probPion() > pid->probProton()) && (pid->probPion() > pid->probKaon()) ) {
						ipippid.push_back(ipip[i]);
				}
		}



		for(int i = 0; i < npim; i++) {
				EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ipim[i];
				RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
				pid->init();
				pid->setMethod(pid->methodProbability());
				pid->setChiMinCut(4);
				pid->setRecTrack(*itTrk);
			//	pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()); // use PID sub-system
				pid->usePidSys(pid->useDedx() | pid->useTofCorr()); // use PID sub-system
				pid->identify(pid->onlyPionKaonProton());
				pid->calculate();
				if(!(pid->IsPidInfoValid())) continue;
				//	if(!(*itTrk)->isMdcDedxValid())continue;
				if((pid->probPion() > pid->probProton()) && (pid->probPion() > pid->probKaon()) ) {
						ipimpid.push_back(ipim[i]);
				}
		}


		m_npippid = ipippid.size();
		m_npimpid = ipimpid.size();
		m_npmpid = ipmpid.size();
		m_npppid = ipppid.size();
		int npimpid = ipimpid.size();
		int npippid = ipippid.size();
		int npmpid = ipmpid.size();
		int npppid = ipppid.size();


		if(npimpid != 1) return SUCCESS;
		if(npippid != 1) return SUCCESS;
		if(npppid != 1) return SUCCESS;
		if(npmpid != 1) return SUCCESS;
		Ncut_pid++;




		////////////////////////////////////////////////////////////////////
		//
		//=======================Kalman Kinematic Fit=========================
		//
		////////////////////////////////////////////////////////////////////
		HepLorentzVector ecms;
		ecms.setPy(0);
		ecms.setPz(0);
		ecms.setE(cms);
		ecms.setPx(cms*sin(0.022*0.5));


		Vp4 p4lambdabarvtx, p4pvtx, p4apvtx, p4pimvtx, p4pipvtx;
		p4lambdabarvtx.clear();
		p4apvtx.clear();
		p4pvtx.clear();
		p4pimvtx.clear();
		p4pipvtx.clear();
		Vdouble  decayL_lambdabar, decayLerr_lambdabar, chisq_lambdabar;
		decayLerr_lambdabar.clear();
		decayL_lambdabar.clear();
		chisq_lambdabar.clear();
		VWTP wlambdabar_vertex;
		wlambdabar_vertex.clear();
		HepPoint3D cPiont;
		RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin()+ipippid[0]))->mdcKalTrack();
		RecMdcKalTrack *pmTrk = (*(evtRecTrkCol->begin()+ipmpid[0]))->mdcKalTrack();

		WTrackParameter wvpmTrk = WTrackParameter(mproton, pmTrk->getZHelixP(), pmTrk->getZErrorP());
		HepLorentzVector p4pm = wvpmTrk.p();
		WTrackParameter wvpipTrk = WTrackParameter(mpion, pipTrk->getZHelix(), pipTrk->getZError());
		HepLorentzVector p4pip = wvpipTrk.p();

		HepPoint3D vx_vtx(0., 0., 0.);
		HepSymMatrix Evx_vtx(3, 0);
		double bx = 1E+6;
		double by = 1E+6;
		double bz = 1E+6;
		Evx_vtx[0][0] = bx*bx;
		Evx_vtx[1][1] = by*by;
		Evx_vtx[2][2] = bz*bz;

		VertexParameter vxpar;
		vxpar.setVx(vx_vtx);
		vxpar.setEvx(Evx_vtx);
		VertexFit* vtxfit = VertexFit::instance();
		vtxfit->init();
		vtxfit->AddTrack(0,  wvpmTrk);
		vtxfit->AddTrack(1,  wvpipTrk);
		vtxfit->AddVertex(0, vxpar,0, 1);
		if(vtxfit->Fit(0)){
				vtxfit->Swim(0);
				vtxfit->BuildVirtualParticle(0);
				WTrackParameter wlambdabar = vtxfit->wVirtualTrack(0);
				VertexParameter vtxlambdabar = vtxfit->vpar(0);
				WTrackParameter wtrkproton = vtxfit->wtrk(0);
				p4apvtx.push_back(wtrkproton.p());
				WTrackParameter wtrkpion = vtxfit->wtrk(1);
				p4pipvtx.push_back(wtrkpion.p());
				SecondVertexFit *vtxfit1 = SecondVertexFit::instance();
				vtxfit1->init();
				vtxfit1->setPrimaryVertex(vx_db);
				vtxfit1->AddTrack(0,wlambdabar);
				vtxfit1->setVpar(vtxlambdabar);
				if(vtxfit1->Fit()) {
						HepLorentzVector p4lambdabar = vtxfit1->p4par();

						p4lambdabarvtx.push_back(p4lambdabar);
						decayL_lambdabar.push_back(vtxfit1->decayLength());
						decayLerr_lambdabar.push_back(vtxfit1->decayLengthError());
						chisq_lambdabar.push_back(vtxfit1->chisq());
						wlambdabar_vertex.push_back(vtxfit1->wpar());
						cPiont = vtxfit1->crossPoint();
				}
		}


		Vp4 p4lambdavtx;
		p4lambdavtx.clear();
		Vdouble  decayL_lambda, decayLerr_lambda, chisq_lambda;
		decayLerr_lambda.clear();
		decayL_lambda.clear();
		chisq_lambda.clear();
		VWTP wlambda_vertex;
		wlambda_vertex.clear();
		RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin()+ipimpid[0]))->mdcKalTrack();
		RecMdcKalTrack *ppTrk = (*(evtRecTrkCol->begin()+ipppid[0]))->mdcKalTrack();

		WTrackParameter wvppTrk = WTrackParameter(mproton, ppTrk->getZHelixP(), ppTrk->getZErrorP());
		HepLorentzVector p4pp = wvppTrk.p();
		WTrackParameter wvpimTrk = WTrackParameter(mpion, pimTrk->getZHelix(), pimTrk->getZError());
		HepLorentzVector p4pim = wvpimTrk.p();
/*
		HepPoint3D vx_vtx(0., 0., 0.);
		HepSymMatrix Evx_vtx(3, 0);
		double bx = 1E+6;
		double by = 1E+6;
		double bz = 1E+6;
		Evx_vtx[0][0] = bx*bx;
		Evx_vtx[1][1] = by*by;
		Evx_vtx[2][2] = bz*bz;
*/
	//	VertexParameter vxpar;
	//	vxpar.setVx(vx_vtx);
	//	vxpar.setEvx(Evx_vtx);
	//	VertexFit* vtxfit = VertexFit::instance();
		vtxfit->init();
		vtxfit->AddTrack(0,  wvppTrk);
		vtxfit->AddTrack(1,  wvpimTrk);
		vtxfit->AddVertex(0, vxpar,0, 1);
		if(vtxfit->Fit(0)){
				vtxfit->Swim(0);
				vtxfit->BuildVirtualParticle(0);
				WTrackParameter wlambda = vtxfit->wVirtualTrack(0);
				VertexParameter vtxlambda = vtxfit->vpar(0);
				WTrackParameter wtrkproton = vtxfit->wtrk(0);
				p4pvtx.push_back(wtrkproton.p());
				WTrackParameter wtrkpion = vtxfit->wtrk(1);
				p4pimvtx.push_back(wtrkpion.p());
				SecondVertexFit *vtxfit1 = SecondVertexFit::instance();
				vtxfit1->init();
				vtxfit1->setPrimaryVertex(vx_db);
				vtxfit1->AddTrack(0,wlambda);
				vtxfit1->setVpar(vtxlambda);
				if(vtxfit1->Fit()) {
						HepLorentzVector p4lambda = vtxfit1->p4par();

						p4lambdavtx.push_back(p4lambda);
						decayL_lambda.push_back(vtxfit1->decayLength());
						decayLerr_lambda.push_back(vtxfit1->decayLengthError());
						chisq_lambda.push_back(vtxfit1->chisq());
						wlambda_vertex.push_back(vtxfit1->wpar());
				}
		}


		int nlambda = p4lambdavtx.size();
		int nlambdabar = p4lambdabarvtx.size();
		if(nlambda == 0 || nlambdabar == 0) return SUCCESS;
		Ncut_vertex++;
		double pull_temp = 999.;
		HepLorentzVector p4_lambda, p4_lambdabar, p4_p, p4_ap, p4_pip, p4_pim;
		double decayL_lam, decayL_lambar, decayLerr_lam, decayLerr_lambar, chisq_lam, chisq_lambar;
		WTrackParameter wlam, wlambar;

		for(int i = 0; i < nlambda; i++){
				for(int j = 0; j < nlambdabar; j++){
						HepLorentzVector p4lambda = p4lambdavtx[i];
						HepLorentzVector p4lambdabar = p4lambdabarvtx[j];

						double pull = (p4lambda.m() - lambda_mass)*(p4lambda.m() - lambda_mass) + (p4lambdabar.m() - lambda_mass)*(p4lambdabar.m() - lambda_mass);
						if(pull_temp > pull){
								pull_temp = pull;
								p4_lambda = p4lambdavtx[i];
								p4_lambdabar = p4lambdabarvtx[j];
								decayL_lam = decayL_lambda[i];
								decayLerr_lam = decayLerr_lambda[i];
								decayL_lambar = decayL_lambdabar[j];
								decayLerr_lambar = decayLerr_lambdabar[j];
								chisq_lam = chisq_lambda[i];
								chisq_lambar = chisq_lambdabar[i];
								wlam = wlambda_vertex[i];
								wlambar = wlambdabar_vertex[j];
								p4_p = p4pvtx[i];
								p4_pim = p4pimvtx[j];
								p4_ap = p4apvtx[j];
								p4_pip = p4pipvtx[j];
						}
				}
		}
		if(pull_temp > 998.) return SUCCESS;
		double lam_mass = p4_lambda.m();
		double lambar_mass = p4_lambdabar.m();
		HepLorentzVector p4_jpsi = p4_lambda + p4_lambdabar;

		double jpsi_mass = p4_jpsi.m();

		if( fabs(lam_mass - lambda_mass) > 0.1) return SUCCESS;
		if( fabs(lambar_mass - lambda_mass) > 0.1) return SUCCESS;
		if( fabs(jpsi_mass - cms) > 0.2) return SUCCESS;
		m_jpsi_mass = jpsi_mass;
		hcos_vertex->Fill(mc_p4lambda.cosTheta());
		Ncut_mass++;


		m_lamkmf_mass = p4_lambda.m();
		m_lamkmf_cos  = p4_lambda.cosTheta();
		m_lamkmfbar_mass = p4_lambdabar.m();
		m_lamkmfbar_cos = p4_lambdabar.cosTheta();
		m_decayL_lamkmf  = decayL_lam;
		m_decayL_lamkmfbar = decayL_lambar;
		m_decayLerr_lamkmf  = decayLerr_lam;
		m_decayLerr_lamkmfbar = decayLerr_lambar;
		m_chisq_lamkmf = chisq_lam;
		m_chisq_lamkmfbar = chisq_lambar;

		m_lambda_mass = p4_lambda.m();
		m_lambda_cos  = p4_lambda.cosTheta();
		m_lambdabar_mass = p4_lambdabar.m();
		m_lambdabar_cos = p4_lambdabar.cosTheta();
		m_decayL_lambda  = decayL_lam;
		m_decayL_lambdabar = decayL_lambar;
		m_decayLerr_lambda  = decayLerr_lam;
		m_decayLerr_lambdabar = decayLerr_lambar;
		m_chisq_lambda = chisq_lam;
		m_chisq_lambdabar = chisq_lambar;

		m_p4lambda[0] = p4_lambda.e(); m_p4lambda[1] = p4_lambda.px(); m_p4lambda[2] = p4_lambda.py(); m_p4lambda[3] = p4_lambda.pz();
		m_p4lambdabar[0] = p4_lambdabar.e(); m_p4lambdabar[1] = p4_lambdabar.px(); m_p4lambdabar[2] = p4_lambdabar.py(); m_p4lambdabar[3] = p4_lambdabar.pz();


		m_p4p[0] = p4_p.e(); m_p4p[1] = p4_p.px(); m_p4p[2] = p4_p.py(); m_p4p[3] = p4_p.pz();
		m_p4ap[0] = p4_ap.e(); m_p4ap[1] = p4_ap.px(); m_p4ap[2] = p4_ap.py(); m_p4ap[3] = p4_ap.pz();
		m_p4pim[0] = p4_pim.e(); m_p4pim[1] = p4_pim.px(); m_p4pim[2] = p4_pim.py(); m_p4pim[3] = p4_pim.pz();
		m_p4pip[0] = p4_pip.e(); m_p4pip[1] = p4_pip.px(); m_p4pip[2] = p4_pip.py(); m_p4pip[3] = p4_pip.pz();

		//============================ Kinematic Fitting =======================================

		double chisq = 999.;
		KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
		kmfit->init();
		kmfit->AddTrack(0, wlam);
		kmfit->AddTrack(1, wlambar);
		kmfit->AddFourMomentum(0, ecms);
		bool oksq = kmfit->Fit();
		if( oksq ){
				chisq = kmfit->chisq();
		}
		//	if(chisq > 990) return SUCCESS;
		hcos_kin->Fill(mc_p4lambda.cosTheta());
		m_chisq = chisq;

		m_tuple1->write();

		HepLorentzVector mc_p4nbar = p4_ap;
		HepLorentzVector mc_posinbar, mc_posisigmabar;
		mc_posinbar.setX(cPiont.x());
		mc_posinbar.setY(cPiont.y());
		mc_posinbar.setZ(cPiont.z());
		mc_posinbar.setT(0);
		mc_posisigmabar.setX(vx[0]);
		mc_posisigmabar.setY(vx[1]);
		mc_posisigmabar.setZ(vx[2]);
		mc_posisigmabar.setT(0);


		Vint ishower;
		ishower.clear();

		RecEmcShower *nbar_Trk;
		IAntiNeutronCorrectionSvc* nbar_svc;
		StatusCode sc_AntiNeutronCorrectionSvc = service("AntiNeutronCorrectionSvc", nbar_svc);
		if ( sc_AntiNeutronCorrectionSvc.isFailure() ){
				log << MSG::FATAL << "Could not load AntiNeutronCorrectionSvc!" << endreq;
				return sc_AntiNeutronCorrectionSvc;
		}
		m_nbar_svc = dynamic_cast<AntiNeutronCorrectionSvc*>(nbar_svc);
		runNo = -abs(runNo);
		if(runNo < 0){
				m_nbar_svc->setAntiNeutronTrk(mc_p4nbar, mc_posinbar, mc_posisigmabar, 0);
				if(m_nbar_svc->isAntiNeutronCorrectionValid()){
						nbar_Trk = m_nbar_svc->getNbarShower();
					//	EvtRecTrack* aNewEvtRecTrack = new EvtRecTrack;
					//	aNewEvtRecTrack->setTrackId(evtRecEvent->totalTracks());
					//	aNewEvtRecTrack->setEmcShower(nbar_Trk);
					//	evtRecTrkCol->push_back(aNewEvtRecTrack);
					//	StatusCode sc = eventSvc()->registerObject(EventModel::EvtRec::EvtRecTrackCol, evtRecTrkCol);
					//	if ( sc.isFailure() ) {
					//			log << MSG::FATAL << "Could not register EvtRecTrackCol in TDS!" << endreq;
					//			return( StatusCode::FAILURE);
					//	}
						Hep3Vector emcpos(nbar_Trk->x(), nbar_Trk->y(), nbar_Trk->z());
						double dang = 200.;
						for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
								EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
								if(!(*jtTrk)->isExtTrackValid()) continue;
								RecExtTrack *extTrk = (*jtTrk)->extTrack();
								if(extTrk->emcVolumeNumber()==-1) continue;
								Hep3Vector extpos = extTrk->emcPosition();
								double angd = extpos.angle(emcpos);

								if (angd < dang) {
										dang = angd;
								}
						}
						if(dang != 200.){
								dang = dang * 180 / (CLHEP::pi);
						}
						//	if(fabs(dang) < m_gammaAngleCut) return SUCCESS;
						ishower.push_back(evtRecEvent->totalTracks());   // means a brand new simulation of anti-neutron. 

				}
		}
		//////////////////////////////////////////////////////////////////////


		int nshower= ishower.size();
		if(nshower < 1) return SUCCESS;
		m_nshower = nshower;

		RecEmcShower *nbarTrk;
		int nbar_idx = -1;
		double energy = 0;
		if(runNo < 0){
				if(!m_nbar_svc->isAntiNeutronCorrectionValid()) return SUCCESS;
				nbar_idx = nshower - 1;
				nbarTrk = nbar_Trk;
		}
		else {
				for (int k = 0; k < nshower; k++){
						EvtRecTrackIterator itTrk = evtRecTrkCol->begin()+ishower[k];
						if(!(*itTrk)->isEmcShowerValid()) continue;
						RecEmcShower *emcTrk = (*itTrk)->emcShower();
						double eraw = emcTrk->energy();
						if(eraw>2.0) continue;
						if(eraw>energy) {
								energy = eraw;
								nbar_idx = k;
						}
				}
				if(nbar_idx == -1) return StatusCode::SUCCESS;
				nbarTrk = (*(evtRecTrkCol->begin()+ishower[nbar_idx]))->emcShower();
				if(nbarTrk->energy() < 0.48 || nbarTrk->numHits()<20 || nbarTrk->secondMoment()<18 ) return SUCCESS;
				m_nbar_svc->setErrorMatrix(nbarTrk);
		}

		Ncut_nbar++;


		double nbarEnergy =    cms/2 - p4_pip.e();
		nbarTrk->setEnergy(nbarEnergy);
		HepLorentzVector p4nbarraw;
		double eraw =  nbarEnergy;
		if(eraw*eraw < nbar_mass*nbar_mass) return StatusCode::SUCCESS;
		double praw = sqrt(eraw*eraw - nbar_mass*nbar_mass);
		double the = nbarTrk->theta();
		double phi = nbarTrk->phi();
		p4nbarraw.setPx(praw*sin(the)*cos(phi));
		p4nbarraw.setPy(praw*sin(the)*sin(phi));
		p4nbarraw.setPz(praw*cos(the));
		p4nbarraw.setE(eraw);

		HepLorentzVector p4nraw = ecms - p4_pim - p4_pip - p4nbarraw;

		double chi2 = 999.;
	//	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
		kmfit->init();
		kmfit->setIterNumber(20);
		kmfit->setBeamPosition(vx);
		kmfit->setVBeamPosition(Evx);
		kmfit->AddTrack(0, wvpimTrk);
		kmfit->AddTrack(1, wvpipTrk);
		kmfit->AddMissTrack(2, nbar_mass, nbarTrk);
		kmfit->myAddMissTrack(3, n_mass, p4nraw);
		kmfit->AddResonance(0, lambda_mass, 1, 2);
		kmfit->AddFourMomentum(1, ecms);
		bool oksq1 = kmfit->Fit();
	//	m_niter = kmfit->getIterNumber();
		if(oksq1) {
				chi2 = kmfit->chisq();
				HepLorentzVector p4pim = kmfit->pfit(0);
				HepLorentzVector p4pip = kmfit->pfit(1);
				HepLorentzVector p4nbar = kmfit->pfit(2);
				HepLorentzVector p4n = kmfit->pfit(3);
				HepLorentzVector p4sigmabar = p4pip+p4nbar;
				HepLorentzVector p4sigma = ecms -p4sigmabar;

				m_sigma_mass = p4sigma.m();
				m_sigma_cos = p4sigma.cosTheta();
				m_pxpim = p4pim.px(); 	m_pypim = p4pim.py();	m_pzpim = p4pim.pz();	m_epim = p4pim.e();
				m_pxpip = p4pip.px(); 	m_pypip = p4pip.py();	m_pzpip = p4pip.pz();	m_epip = p4pip.e();
				m_pxnbar = p4nbar.px(); 	m_pynbar = p4nbar.py();	m_pznbar = p4nbar.pz();	m_enbar = p4nbar.e();
				m_pxn = p4n.px(); 	m_pyn = p4n.py();	m_pzn = p4n.pz();	m_en = p4n.e();
		}	

		m_chisq_kmf = chi2;
		m_tuple2->write();

		if ( m_writeMyDst )  m_subAlgA->execute(); //  execute event filter
		Ncut_kin++;

		return StatusCode::SUCCESS;

}

//---------------------------------------------------------------------------------
StatusCode Lambda::finalize()	{
		cout << "Selection criteria:	" << "Survived : Taglambdabar pbarnpipi  V 0 - 4"    << endl;
		cout << "Total number:			" << Ncut_none     << endl;
		cout << "cut mc:				" << Ncut_mc       << endl;
		cout << "Charged trk:			" << Ncut_ncharge  << endl;
		cout << "Vertex sel:				" << Ncut_vertex    << endl;
		cout << "Mass sel:				" << Ncut_mass      << endl;
		cout << "kin :					" << Ncut_kin  << endl;
		MsgStream log(msgSvc(), name());
		log << MSG::INFO << "in finalize()" << endmsg;
		return StatusCode::SUCCESS;
}



