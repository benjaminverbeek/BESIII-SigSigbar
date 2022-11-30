#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include <iostream>
#include "VertexFit/Helix.h" 
#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "McTruth/McParticle.h"
#include "VertexFit/IVertexDbSvc.h"
#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h" 
#include "TLorentzVector.h"
#include "TVector3.h"
using CLHEP::Hep3Vector;
#include "CLHEP/Vector/LorentzVector.h"
using CLHEP::HepLorentzVector;
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep2Vector;
#include "CLHEP/Geometry/Point3D.h"
//#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
//#endif
//#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/HTrackParameter.h"
#include "VertexFit/SecondVertexFit.h" 
#include "ParticleID/ParticleID.h"
#include <vector>
#include "../XiXibarAlg/XiXibar.h"
#include "../XiXibarAlg/AngDis.hh"
#include "TRandom.h"
#include "DstEvent/DstTofTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "TofRecEvent/RecTofTrack.h"

//	SmartRefVector<RecTofTrack
//	EvtRecTrackCol
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
const double m_Psi2S    = 3.686097;
const double m_Jpsi     = 3.0969;
const double m_Xich     = 1.32171;
const double m_Lambda   = 1.115683;
const double xmass[5] = {0.000511, 0.105658, 0.139570, 0.493677, 0.938272};
//const double ecms = 3.686;
const double ecms = m_Jpsi;
const double xmpi0 = 0.13497;
const double xmk0 = 0.49767;
double par[7];

int Ncut0,Ncut1,Ncut2,Ncut3,Ncut4,Ncut5,Ncut6,Ncut7,Ncut8,Ncut9,Ncut10;
int Neff;

////////////////////////////////////////////////////////////////////
XiXibar::XiXibar(const std::string& name, ISvcLocator* pSvcLocator) : Algorithm(name, pSvcLocator) 
{  
//	outputFile = NULL;
//  	tree = NULL;
	// declareProperty("CmsEnergy", m_ecms = 3.686);
    declareProperty("CmsEnergy", m_ecms = m_Jpsi);
	declareProperty("Vr0cut", m_vr0cut=1.0);
	declareProperty("Vz0cut", m_vz0cut=10.0);
	declareProperty("EnergyThreshold_b",m_energyThreshold_b=0.025);
	declareProperty("EnergyThreshold_e",m_energyThreshold_e=0.050);
	declareProperty("BarrelEmc_th",  m_Barrel_th    = 0.8);
	declareProperty("EndcapEmc_th_1", m_Endcap_th_1  = 0.84);
	declareProperty("EndcapEmc_th_2", m_Endcap_th_2  = 0.92);
	declareProperty("EnergyThreshold", m_energyThreshold=0.02);
	declareProperty("GammaPhiCut", m_gammaPhiCut=20.0);
	declareProperty("GammaThetaCut", m_gammaThetaCut=20.0);
	declareProperty("GammaAngCut", m_gammaAngCut=20.0);
	declareProperty("CheckMC", m_checkMC = false);
	declareProperty("idNo",m_idNo = 443);
	//declareProperty("idNo",m_idNo = 100443);
	//Declare the properties  
}


HepLorentzVector XiXibar::GetEcmsP4(){
	HepLorentzVector Ecms_jpsi(0.011*m_ecms,0.0,0.0,m_ecms);
	return Ecms_jpsi;
}


void XiXibar::InitializeReconstructedSpectra() {
	hp_pr_lam_vx			= new TH1D("hp_pr_lam_vx","; momentum p (GeV/c)", 100, 0., 1.0); 
	hp_pr_bar_lambar_vx		= new TH1D("hp_pr_bar_lambar_vx","; momentum #bar{p} (GeV/c)", 100, 0., 1.0);
	hp_pim_lam_vx			= new TH1D("hp_pim_lam_vx","; momentum #pi^{-} #Lambda (GeV/c)", 100, 0., 1.0); 
	hp_pip_lambar_vx		= new TH1D("hp_pip_lambar_vx","; momentum #pi^{+} #bar{#Lambda} (GeV/c)", 100, 0., 1.0);	
	hp_pim_Xim_vx			= new TH1D("hp_pim_Xim_vx","; momentum bachelor #pi^{-} #Xi^{-} (GeV/c)", 100, 0., 1.0); 
	hp_pip_Xipbar_vx		= new TH1D("hp_pip_Xipbar_vx","; momentum bachelor #pi^{+} #bar{#Xi}^{+} GeV/c", 100, 0., 1.0);

	hpvth_pr_lam_vx			= new TH2D("hpvth_pr_lam_vx", "momentum vs #theta_{lab} proton; p (GeV/c); #theta (^{o})", 100, 0., 1.0, 90, 0, 180);
	hpvth_pim_lam_vx		= new TH2D("hpvth_pim_lam_vx", "momentum vs #theta_{lab} #pi^{-} #Lambda; p (GeV/c); #theta (^{o})", 100, 0., 1.0, 90, 0, 180);
	hpvth_pim_Xim_vx		= new TH2D("hpvth_pim_Xim_vx", "momentum vs #theta_{lab} #pi^{-} #Xi^{-}; p (GeV/c); #theta (^{o})", 100, 0., 1.0, 90, 0, 180);

	hpvth_pr_bar_lambar_vx	= new TH2D("hpvth_pr_bar_lambar_vx", "momentum vs #theta_{lab} anti-proton; p (GeV/c); #theta (^{o})", 100, 0., 1.0, 90, 0, 180);
	hpvth_pip_lambar_vx		= new TH2D("hpvth_pip_lambar_vx", "momentum vs #theta_{lab} #pi^{+} #bar{#Lambda}; p (GeV/c); #theta (^{o})", 100, 0., 1.0, 90, 0, 180);
	hpvth_pip_Xip_vx		= new TH2D("hpvth_pip_Xip_vx", "momentum vs #theta_{lab} #pi^{+} #bar{#Xi}^{+}; p (GeV/c); #theta (^{o})", 100, 0., 1.0, 90, 0, 180);

	hlam_decaylength		= new TH1D("hlam_decaylength","; #Lambda Decay Length (mm)", 100, 0., 100.); 
	hlam_bar_decaylength	= new TH1D("hlam_bar_decaylength","; #bar{#Lambda} Decay Length (mm)", 100, 0., 100.); 

	h_chi2_prvx_lam			= new TH1D("h_chi2_prvx_lam","; #chi^{2} #Lambda pr.vx.", 500, 0, 500);
	h_chi2_prvx_lam_bar		= new TH1D("h_chi2_prvx_lam_bar","; #chi^{2} #Lambda sec.vx.", 500, 0, 500);
	h_chi2_scvx_lam			= new TH1D("h_chi2_scvx_lam","; #chi^{2} #bar{#Lambda} pr.vx.", 500, 0, 500);
	h_chi2_scvx_lam_bar		= new TH1D("h_chi2_scvx_lam_bar","; #chi^{2} #bar{#Lambda} sec.vx.", 500, 0, 500);

	hm_lam_vx 		= new TH1D("hm_lam_vx","; m(#Lambda) GeV/c^{2}",100, 1.05, 1.15);
	hm_lam_bar_vx 	= new TH1D("hm_lam_bar_vx","; m(#bar{#Lambda}) GeV/c",100, 1.05, 1.15);
	hm_xim_vx 		= new TH1D("hm_xim_vx","; m(#Xi^{-}) GeV/c^{2}", 150, 1.2, 1.5);
	hm_xip_bar_vx 	= new TH1D("hm_xip_bar_vx","; m(#bar{#Xi}^{+}) GeV/c^{2}",150, 1.2, 1.5);
	hm_jpsi_vx 	= new TH1D("hm_jpsi_vx","; J/#Psi (GeV/c^{2})",100, 2.8, 3.2);

	hm_tof_proton 	= new TH1D("hm_tof_proton","; proton ToF )",100, -200., 200.);
	hm_tof_pion 	= new TH1D("hm_tof_pion","; pion ToF )",100, -200., 200.);

	hm_tof_proton_best 	= new TH1D("hm_tof_proton_best","; proton ToF )",100, -200., 200.);
	hm_tof_pion_best1 	= new TH1D("hm_tof_pion_best1","; best pion 1 ToF )",100, -200., 200.);
	hm_tof_pion_best2 	= new TH1D("hm_tof_pion_best2","; best pion 2 ToF )",100, -200., 200.);

	hweight 	= new TH1D("hweight","; weight ",250, -5, 5.);

	efficiency			= new TH1D("efficiency","; efficiency of cuts", 10, 0., 10); 

	hm_lam_lambar_vx 		= new TH2D("hm_lam_lambar_vx","; m(#Lambda) GeV/c^{2} ; m(#bar{#Lambda}) GeV/c^{2}",50, 1.05, 1.15, 50, 1.05, 1.15);
	hm_lam_lambar_vx_final  = new TH2D("hm_lam_lambar_vx_final","; m(#Lambda) GeV/c^{2} ; m(#bar{#Lambda}) GeV/c^{2}",50, 1.05, 1.15, 50, 1.05, 1.15);

	hm_xim_xipbar_vx 	= new TH2D("hm_xim_xipbar_vx","; m(#Xi^{-}) GeV/c^{2} ; m(#bar{#Xi}^{+}) GeV/c^{2}",75, 1.2, 1.5, 75, 1.2, 1.5);
	hm_xim_xipbar_fit 	= new TH2D("hm_xim_xipbar_fit","Fitted masses; m(#Xi^{-}) GeV/c^{2} ; m(#bar{#Xi}^{+}) GeV/c^{2}",75, 1.2, 1.5, 75, 1.2, 1.5);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode XiXibar::initialize()
{
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;

	initialiezeNtuple();
	InitializeReconstructedSpectra();

	log << MSG::INFO << "successfully return from initialize()" <<endmsg;

	return StatusCode::SUCCESS;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode XiXibar::execute() 
{

	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;

	StatusCode sc = StatusCode::SUCCESS;


	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	if (!eventHeader)
	{
		log << MSG::FATAL << "Could not find Event Header" << endreq;
		return StatusCode::FAILURE;
	}
	m_run = eventHeader->runNumber();
	m_rec = eventHeader->eventNumber();
	// m_flag1 = eventHeader->flag1(); m_flag2 = eventHeader->flag2(); // what is flag1 and flag2 used for?
	int runNo = m_run;
	int event = m_rec;

	// set true if you want to generate a MC sample with a specific angular distributions. Used for x-check of XiXibar results in data
	bool mc_true = false;
	
	//TrackCorrection * calibrate = new TrackCorrection();

	bool TrackCorr = false;
	if(eventHeader->runNumber() > 0) TrackCorr = false;

	if(eventHeader->runNumber() < 0) TrackCorr = false;

	if(m_rec%2000==0)
		cout<<"Run   "<<m_run<<"     Event   "<<m_rec<<endl; 
	//MC information
	if (eventHeader->runNumber() < 0){
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), 
		"/Event/MC/McParticleCol");
		int m_numParticle = 0;
		if (!mcParticleCol){
			std::cout << "Could not retrieve McParticelCol" << std::endl;
			return StatusCode::FAILURE;
		}
		else{

			bool jpsiDecay = false;
			bool psipDecay = false;
			int rootIndex = -1;
			Event::McParticleCol::iterator iter_mc;
			int m_pim = 0;
			int m_pip = 0;
			int mc_numParticle = 0;

			//std::cout << "event begin" << std::endl;
			for(iter_mc=mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++){
			
				// for each TLorentzVector produce 4 n-tuples and assign them the values of e, px, py, pz

				Double_t init = -9999.9;

				HepXiTrue.set(init, init, init, init);
				HepXibarTrue.set(init, init, init, init);
				HepLambdaTrue.set(init, init, init, init);
				HepLambdabarTrue.set(init, init, init, init);
				HepProtonTrue.set(init, init, init, init);
				HepPimLamTrue.set(init, init, init, init);
				HepPimXiTrue.set(init, init, init, init);
				HepProtonbarTrue.set(init, init, init, init);
				HepPipLambarTrue.set(init, init, init, init);
				HepPipXibarTrue.set(init, init, init, init);	

				if ((*iter_mc)->particleProperty()==m_idNo){
					jpsiDecay = true;
					rootIndex = (*iter_mc)->trackIndex();
				}

				long mother = ((*iter_mc)->mother()).trackIndex() - rootIndex;
				int pdgid = (*iter_mc)->particleProperty();
				mc_pdgid[mc_numParticle] = pdgid;
				mc_motheridx[mc_numParticle] = mother;
				mc_numParticle +=1;

				if(pdgid == 2112){return StatusCode::SUCCESS;} 

				if(pdgid == -2112){return StatusCode::SUCCESS;}

				if(pdgid == 3312) {
					HepXiTrue 	= (*iter_mc)->initialFourMomentum();
					mc_Xim_e 	= HepXiTrue.e();mc_Xim_px 	= HepXiTrue.px(); 
					mc_Xim_py 	= HepXiTrue.py();mc_Xim_pz 	= HepXiTrue.pz();
				}
	            if(pdgid == -3312){
	                HepXibarTrue = (*iter_mc)->initialFourMomentum();
					mc_Xip_e 	= HepXibarTrue.e();mc_Xip_px 	= HepXibarTrue.px();
        			mc_Xip_py 	= HepXibarTrue.py();mc_Xip_pz 	= HepXibarTrue.pz();	            
        		}
				if(pdgid == 3122){
				 	HepLambdaTrue = (*iter_mc)->initialFourMomentum();
					mc_Lam_e 	= HepLambdaTrue.e();mc_Lam_px 	= HepLambdaTrue.px();
        			mc_Lam_py 	= HepLambdaTrue.py();mc_Lam_pz 	= HepLambdaTrue.pz();				
        		}
				if(pdgid == -3122){
				 	HepLambdabarTrue = (*iter_mc)->initialFourMomentum();
					mc_Lambar_e 	= HepLambdabarTrue.e();mc_Lambar_px 	= HepLambdabarTrue.px();
        			mc_Lambar_py 	= HepLambdabarTrue.py();mc_Lambar_pz 	= HepLambdabarTrue.pz();					
        		}
				if(pdgid==2212){
					HepProtonTrue = (*iter_mc)->initialFourMomentum();
					mc_pr_e 	= HepProtonTrue.e();mc_pr_px 	= HepProtonTrue.px();
        			mc_pr_py 	= HepProtonTrue.py();mc_pr_pz 	= HepProtonTrue.pz();	
				}
				if(pdgid==-2212){
					HepProtonbarTrue = (*iter_mc)->initialFourMomentum();
					mc_prbar_e 		= HepProtonbarTrue.e();mc_prbar_px 	= HepProtonbarTrue.px();
        			mc_prbar_py 	= HepProtonbarTrue.py();mc_prbar_pz 	= HepProtonbarTrue.pz();	
   				}	
				if(pdgid==211){ // pi+
					if(m_pip == 0){
						HepPipXibarTrue=(*iter_mc)->initialFourMomentum();
						mc_pipXip_e		= HepPipXibarTrue.e();mc_pipXip_px 	= HepPipXibarTrue.px();
        				mc_pipXip_py 	= HepPipXibarTrue.py();mc_pipXip_pz 	= HepPipXibarTrue.pz();
        				m_pip++;
        			}
        			else{
        				HepPipLambarTrue = (*iter_mc)->initialFourMomentum();
						mc_pipLambar_e		= HepPipLambarTrue.e();mc_pipLambar_px 	= HepPipLambarTrue.px();
        				mc_pipLambar_py 	= HepPipLambarTrue.py();mc_pipLambar_pz 	= HepPipLambarTrue.pz();
        			}
				}	
				if(pdgid==-211){ // pi-
					if(m_pim == 0){
						HepPimXiTrue=(*iter_mc)->initialFourMomentum();
						mc_pimXim_e		= HepPimXiTrue.e();mc_pimXim_px 	= HepPimXiTrue.px();
        				mc_pimXim_py 	= HepPimXiTrue.py();mc_pimXim_pz 	= HepPimXiTrue.pz();
        				m_pim++;	
        			}			 	
        			else{
        				HepPimLamTrue=(*iter_mc)->initialFourMomentum();
						mc_pimLam_e		= HepPimLamTrue.e();mc_pimLam_px 	= HepPimLamTrue.px();
        				mc_pimLam_py 	= HepPimLamTrue.py();mc_pimLam_pz 	= HepPimLamTrue.pz();	
        			}
				}
								
			}
			mc_idxmc = mc_numParticle;
		}
	}

	HepXiTrue.set(mc_Xim_px, mc_Xim_py, mc_Xim_pz, mc_Xim_e);
	HepXibarTrue.set(mc_Xip_px, mc_Xip_py, mc_Xip_pz, mc_Xip_e);
	HepLambdaTrue.set(mc_Lam_px, mc_Lam_py, mc_Lam_pz, mc_Lam_e);
	HepLambdabarTrue.set(mc_Lambar_px, mc_Lambar_py, mc_Lambar_pz, mc_Lambar_e);
	HepProtonTrue.set(mc_pr_px, mc_pr_py, mc_pr_pz, mc_pr_e);	
	HepProtonbarTrue.set(mc_prbar_px, mc_prbar_py, mc_prbar_pz, mc_prbar_e);
	HepPimLamTrue.set(mc_pimLam_px,mc_pimLam_py,mc_pimLam_pz,mc_pimLam_e);
    HepPimXiTrue.set(mc_pimXim_px,mc_pimXim_py,mc_pimXim_pz,mc_pimXim_e);
    HepPipLambarTrue.set(mc_pipLambar_px,mc_pipLambar_py,mc_pipLambar_pz,mc_pipLambar_e);
    HepPipXibarTrue.set(mc_pipXip_px,mc_pipXip_py,mc_pipXip_pz,mc_pipXip_e);  

	if(eventHeader->runNumber() < 0){
		if(mc_true){
			bool test = WeighDataTrue();
			if(!test)
				return StatusCode::SUCCESS;
		}
	}	

	m_tuple->write();	

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),EventModel::EvtRec::EvtRecTrackCol);

	int totchrg = evtRecEvent->totalCharged();
	int totneu  = evtRecEvent->totalNeutral(); 
	Ncut0++;
	Neff++;
	efficiency->Fill(Neff);

	// first cut removes obviously false events from nr of collected tracks
	if( totchrg > 10 || totchrg < 2 ) return StatusCode::SUCCESS;

	HepLorentzVector Ecms_jpsi(0.011*m_ecms,0.0,0.0,m_ecms);
	HepLorentzVector p4psipCM(0.0,0.0,0.0,m_ecms);
	Hep3Vector betaLab = Ecms_jpsi.boostVector();
	Vint iGood,iChrgp,iChrgn;
	iGood.clear();
	iChrgp.clear();
	iChrgn.clear();
	int nCharge = 0;
       //..................................................
       HepPoint3D vx(0., 0., 0.);  //ini vertex
       HepSymMatrix Evx(3, 0);
       double bx = 1E+6;
       double by = 1E+6;
       double bz = 1E+6;
       Evx[0][0] = bx*bx;
       Evx[1][1] = by*by;
       Evx[2][2] = bz*bz;
       VertexParameter vxpar;
       vxpar.setVx(vx);
       vxpar.setEvx(Evx);
       
       Hep3Vector xorigin(0,0,0);
       Hep3Vector eorigin(0,0,0);
       IVertexDbSvc*  vtxsvc;
       Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
       if(vtxsvc->isVertexValid()){
          double* dbv = vtxsvc->PrimaryVertex(); 
          double*  vv = vtxsvc->SigmaPrimaryVertex();
          xorigin.setX(dbv[0]);
          xorigin.setY(dbv[1]);
          xorigin.setZ(dbv[2]);
          eorigin.setX(vv[0]);
          eorigin.setY(vv[1]);
          eorigin.setZ(vv[2]);
       }
       HepPoint3D pvx=xorigin;    //for every run
       HepSymMatrix pEvx(3,0);
       pEvx[0][0]=eorigin[0]*eorigin[0];
       pEvx[1][1]=eorigin[1]*eorigin[1];
       pEvx[2][2]=eorigin[2]*eorigin[2];
       VertexParameter pvxpar;
       pvxpar.setVx(pvx);
       pvxpar.setEvx(pEvx);
   //.................................................
                                                                                                        
	for(int i = 0; i < totchrg; i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		if (!(*itTrk)->isMdcKalTrackValid()) continue;//MdcKalTrk
		RecMdcTrack *mdcTrk =(*itTrk)->mdcTrack();
		double x0 = mdcTrk->x();
		double y0 = mdcTrk->y();
		double z0 = mdcTrk->z();
		double r0 = mdcTrk->r();
		double phi0 = mdcTrk->helix(1);
		double ptrk = mdcTrk->p();
		double chi2 = mdcTrk->chi2();
		double xv=xorigin.x();
		double yv=xorigin.y();
		double zv=xorigin.z();
		double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);
		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
		VFHelix helixip(point0,a,Ea);
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0=vecipa[1];
		double Vctheta =cos(mdcTrk->theta());
		if(fabs(Vctheta) > 0.93)continue;  // standard cut cos(theta) < 0.93
		iGood.push_back(i);
		if(mdcTrk->charge()>0)
		//if(mdcTrk->charge()<0)
			iChrgp.push_back(i);
		else
			iChrgn.push_back(i);
	}

	int ntot=iGood.size();
	m_nTrkp=iChrgp.size(); // nr of positively charged tracks
	m_nTrkm=iChrgn.size(); // nr of negatively charged tracks

    if((iChrgp.size() < 3) || (iChrgn.size() < 3 )) return sc;
    Ncut1++;
    Neff++;
	efficiency->Fill(Neff);    

	double mpip = xmass[2], mpim = xmass[2], mpm = xmass[4], mpp = xmass[4]; // m_pi, m_pi, m_proton, m_proton

	Vp4 ppip,ppm,ppim,ppp;
    Vint ipip,ipm,ipim,ipp;
    ppip.clear();
    ppm.clear();
    ppim.clear();
    ppp.clear();
    ipip.clear();
    ipm.clear();
    ipim.clear();
    ipp.clear();

    for(int i = 0; i<iChrgn.size(); i++){
    	RecMdcKalTrack *negTrk = (*(evtRecTrkCol->begin()+iChrgn[i]))->mdcKalTrack();
    	double mom_cand = negTrk->p();
		if (mom_cand < 0.30) ppim.push_back(iChrgn[i]);
		if (mom_cand > 0.32) ppm.push_back(iChrgn[i]);
    }

    for(int i = 0; i<iChrgp.size(); i++){
    	RecMdcKalTrack *posTrk = (*(evtRecTrkCol->begin()+iChrgp[i]))->mdcKalTrack();
    	double mom_cand = posTrk->p();
		if (mom_cand < 0.30) ppip.push_back(iChrgp[i]);
		if (mom_cand > 0.32) ppp.push_back(iChrgp[i]);
    }

    Vint ppip_cand, ppm_cand, ppim_cand, ppp_cand;
    ppip_cand.clear();
    ppm_cand.clear();
    ppim_cand.clear();
    ppp_cand.clear();


	// First lambda is identified by looping through all tracks and considering the pair with mass closest to Lambda to be the candidate

	//vertex fit
	double chisum=9999.;
	WTrackParameter wlamb1,wlamb2,wlambda1,wlambda2;
	WTrackParameter wLambda, wLambdabar;

	double chi = 9999.;
    std::vector<int> Ximtracks, Xiptracks;		// in order proton, pi(Lambda), pi(Xi)
    Ximtracks.resize(0); Xiptracks.resize(0);
	double length1,length2;
	HepLorentzVector plam,plambar;   
	VertexParameter LambdaVertex, LambdabarVertex;
	// reconstruct Lambda and anti-Lambda via looping positive and negative track combinations
	double deltmin=999.;
	double deltlam_min = 999;
	for(int i = 0; i<iChrgn.size(); i++){
		for(int j = 0; j<iChrgp.size(); j++){
			RecMdcKalTrack *prTrk  = (*(evtRecTrkCol->begin()+iChrgp[j]))->mdcKalTrack();
			RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin()+iChrgn[i]))->mdcKalTrack();

			double mom_prcand = prTrk->p();
			if (mom_prcand < 0.32) continue;
			double mom_picand = pimTrk->p();
			if (mom_picand > 0.30) continue;

		//	double time_proton = TimeOfFlight( iChrgp[j] , 1);
		//	std::cout << "time_proton " << time_proton << std::endl;
		//	hm_tof_proton->Fill(time_proton);

		//	double time_pion = TimeOfFlight( iChrgn[i] , 0);
		//	std::cout << "time_pion " << time_pion << std::endl;
		//	hm_tof_pion->Fill(time_pion);

			RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);
			WTrackParameter  wvprTrk;
			wvprTrk = WTrackParameter(mpp, prTrk->helix(), prTrk->err());
			RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
			WTrackParameter  wvpimTrk;
			wvpimTrk = WTrackParameter(mpim, pimTrk->helix(), pimTrk->err());

			//**************** lambda0 combination  *********************************
			VertexFit *vtxfitlam 			= VertexFit::instance();
			SecondVertexFit *svtxfitlam 	= SecondVertexFit::instance();

			VertexFit *vtxfitXi 			= VertexFit::instance();
			SecondVertexFit *svtxfitXi 		= SecondVertexFit::instance();

			vtxfitlam->init();
			vtxfitlam->AddTrack(0,  wvprTrk);
			vtxfitlam->AddTrack(1,  wvpimTrk);
			vtxfitlam->AddVertex(0, vxpar, 0, 1);
			if(!vtxfitlam->Fit(0)) continue;
			vtxfitlam->Swim(0);
			vtxfitlam->BuildVirtualParticle(0);
			wLambda = vtxfitlam->wVirtualTrack(0); 
			VertexParameter vtxlambda1 = vtxfitlam->vpar(0);

			// new addition 00-00-03
			// HepLorentzVector temp_lam = vtxfitlam->pfit(0) + vtxfitlam->pfit(1);
			// if(TMath::Abs(temp_lam.m()-m_Lambda) > deltlam_min) continue;

			// deltlam_min = TMath::Abs(temp_lam.m()-m_Lambda);
			// end new addition 00-00-03


			for(int k = 0; k < iChrgn.size(); k++){
				if(i == k) continue; // 

				RecMdcKalTrack *pim2Trk = (*(evtRecTrkCol->begin()+iChrgn[k]))->mdcKalTrack();
				double mom_picand = pim2Trk->p();
				if (mom_picand > 0.30) continue;

				RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
				WTrackParameter  wvpim2Trk;
				wvpim2Trk = WTrackParameter(xmass[2], pim2Trk->helix(), pim2Trk->err());

                vtxfitXi->init();
                vtxfitXi->AddTrack(0, wLambda);
                vtxfitXi->AddTrack(1, wvpim2Trk);
                vtxfitXi->AddVertex(0, vxpar, 0, 1); 
                vtxfitXi->setChisqCut(200.);
                vtxfitXi->Fit(0);
                if(!vtxfitXi->Fit()) continue;
                if(vtxfitXi->chisq(0)>200) continue;
                vtxfitXi->BuildVirtualParticle(0);

                VertexParameter vparXim = vtxfitXi->vpar(0);
                WTrackParameter wparXim = vtxfitXi->wVirtualTrack(0);
                HepLorentzVector temp_lam_v4 = vtxfitXi->pfit(0);	// Is this it?
                wparXim.setCharge(-1);

                svtxfitXi->init();
                svtxfitXi->AddTrack(0,wparXim);
                svtxfitXi->setVpar(vparXim);

				    svtxfitXi->setPrimaryVertex(pvxpar);
                //svtxfitXi->setChisqCut(1000.);
                //svtxfitXi->setIterNumber(20); //default is 10
                svtxfitXi->Fit();
                if(!svtxfitXi->Fit()) continue;
                // if(svtxfitXi->chisq()>200) continue;

                
                HepLorentzVector temp_xim_v4 = svtxfitXi->p4par();
                double temp_xim_chi=svtxfitXi->chisq();
                double temp_xim_length=svtxfitXi->decayLength();
                double temp_xim_length_err=svtxfitXi->decayLengthError();
                WTrackParameter temp_wxim=svtxfitXi->wpar();
                VertexParameter temp_vxim=svtxfitXi->vpar();

                double dltm1=sqrt((temp_lam_v4.m()-m_Lambda)*(temp_lam_v4.m()-m_Lambda)+(temp_xim_v4.m()-m_Xich)*(temp_xim_v4.m()-m_Xich))*1000.;
                if(dltm1<deltmin){
                	deltmin=dltm1;
                	LambdaVertex = vparXim;
                	Ximtracks.resize(0);
                	Ximtracks.push_back(iChrgp[j]); 
                	Ximtracks.push_back(iChrgn[i]);
                	Ximtracks.push_back(iChrgn[k]);

                	// obtaining unique number of protons and pions
                	bool pr_cand = true;
                	for(int kk = 0; kk < ppp_cand.size(); kk++)
                		if(iChrgp[j] == ppp_cand[kk]) pr_cand = false;

                	if(pr_cand) ppp_cand.push_back(iChrgp[j]);

                	bool pim_cand = true;
                	for(int kk = 0; kk < ppim_cand.size(); kk++)
                		if(iChrgn[i] == ppim_cand[kk]) pim_cand = false;

                	if(pim_cand) ppim_cand.push_back(iChrgn[i]);

                	pim_cand = true;
                	for(int kk = 0; kk < ppim_cand.size(); kk++)
                		if(iChrgn[k] == ppim_cand[kk]) pim_cand = false;

                	if(pim_cand) ppim_cand.push_back(iChrgn[k]);
   
                }
			}
		}
	}
    if(deltmin>998.)return sc;
    hlam_decaylength->Fill(length1);
    Ncut2++;
    Neff++;
	 efficiency->Fill(Neff);    

    // In the second step one considers all other combinations, not being the lambda pair, forming the anti-lambda pair
    // the track candidates forming the Xim are excluded in the for-loop

    double dltmin2=999.;
    double deltlambar_min = 999.;
    for(int k = 0; k < iChrgn.size(); k++){
     	if( (iChrgn[k] == Ximtracks[1]) || (iChrgn[k] == Ximtracks[2]) )	continue;
        for(int m = 0; m < iChrgp.size(); m++ ){
         	if( iChrgp[m] == Ximtracks[0] )		continue;

            RecMdcKalTrack *prbarTrk = (*(evtRecTrkCol->begin()+iChrgn[k]))->mdcKalTrack();
            RecMdcKalTrack *pipTrk   = (*(evtRecTrkCol->begin()+iChrgp[m]))->mdcKalTrack();
            double mom_prcand = prbarTrk->p();
			if (mom_prcand < 0.32) continue;
			double mom_picand = pipTrk->p();
			if (mom_picand > 0.30) continue;

            RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);
			WTrackParameter  wvprbarTrk;
			wvprbarTrk = WTrackParameter(mpm, prbarTrk->helix(), prbarTrk->err());
            RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
            WTrackParameter  wvpipTrk;
			wvpipTrk = WTrackParameter(mpip, pipTrk->helix(), pipTrk->err());

            //**************** anti-lambda0 combination  *********************************
            VertexFit *vtxfitlambar 		   = VertexFit::instance();
            SecondVertexFit *svtxfitlambar 	= SecondVertexFit::instance();

            VertexFit *vtxfitXi 			      = VertexFit::instance();
			   SecondVertexFit *svtxfitXi 		= SecondVertexFit::instance();

            vtxfitlambar->init();
            vtxfitlambar->AddTrack(0,  wvprbarTrk);
            vtxfitlambar->AddTrack(1,  wvpipTrk);
            vtxfitlambar->AddVertex(0, vxpar,0,1);
            if(!vtxfitlambar->Fit(0)) 	continue;
            vtxfitlambar->Swim(0);
            vtxfitlambar->BuildVirtualParticle(0);
            wLambdabar = vtxfitlambar->wVirtualTrack(0); 
            VertexParameter vtxlambda2 = vtxfitlambar->vpar(0);			

            // new addition 00-00-03
            // HepLorentzVector temp_lambar = vtxfitlambar->pfit(0) + vtxfitlambar->pfit(1);
			// if(TMath::Abs(temp_lambar.m()-m_Lambda) > deltlambar_min) continue;
				
			// deltlambar_min = TMath::Abs(temp_lambar.m()-m_Lambda);
			// end new addition 00-00-03

			for(int n = 0; n<iChrgp.size(); n++){
				if( (iChrgp[n] == Ximtracks[0]) || (iChrgp[n] == iChrgp[m]))	continue;

				RecMdcKalTrack *pip2Trk = (*(evtRecTrkCol->begin()+iChrgp[n]))->mdcKalTrack();
				double mom_picand = pip2Trk->p();
				if (mom_picand > 0.30) continue;
				RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
				WTrackParameter  wpip2Trk;
				wpip2Trk = WTrackParameter(xmass[2], pip2Trk->helix(), pip2Trk->err());
              
                vtxfitXi->init();
                vtxfitXi->AddTrack(0, wLambdabar);
                vtxfitXi->AddTrack(1, wpip2Trk);
                vtxfitXi->AddVertex(0, vxpar, 0, 1);
                vtxfitXi->setChisqCut(200.); 
                vtxfitXi->Fit(0);
                if(!vtxfitXi->Fit()) continue;
                if(vtxfitXi->chisq(0)>200) continue;
                vtxfitXi->BuildVirtualParticle(0);
                VertexParameter vparXip = vtxfitXi->vpar(0);
                WTrackParameter wparXip = vtxfitXi->wVirtualTrack(0);
                HepLorentzVector temp_lambar_v4 = vtxfitXi->pfit(0);
                wparXip.setCharge(+1);

                svtxfitXi->init();
                svtxfitXi->AddTrack(0,wparXip);
                svtxfitXi->setVpar(vparXip);
				svtxfitXi->setPrimaryVertex(pvxpar);
                //svtxfitXi->setChisqCut(1000.);
                //svtxfitXi->setIterNumber(20);
                svtxfitXi->Fit();
                if(!svtxfitXi->Fit()) continue;
                // if(svtxfitXi->chisq()>200) continue;
               
                HepLorentzVector temp_xip_v4=svtxfitXi->p4par();
                double temp_xip_chi=svtxfitXi->chisq();
                double temp_xip_length=svtxfitXi->decayLength();
                double temp_xip_length_err=svtxfitXi->decayLengthError();
                WTrackParameter temp_wxip=svtxfitXi->wpar();
                VertexParameter temp_vxip=svtxfitXi->vpar();

                double dltm2=sqrt((temp_lambar_v4.m()-m_Lambda)*(temp_lambar_v4.m()-m_Lambda)+(temp_xip_v4.m()-m_Xich)*(temp_xip_v4.m()-m_Xich))*1000.;
                if(dltm2<dltmin2){
                 	dltmin2=dltm2;
                 	LambdabarVertex = vparXip;
                 	Xiptracks.resize(0);
                 	Xiptracks.push_back(iChrgn[k]); 
                 	Xiptracks.push_back(iChrgp[m]);
                 	Xiptracks.push_back(iChrgp[n]);


                 	// obtaining unique number of antiprotons and pions
                	bool prbar_cand = true;
                	for(int kk = 0; kk < ppm_cand.size(); kk++)
                		if(iChrgn[k] == ppm_cand[kk]) prbar_cand = false;

                	if(prbar_cand) ppm_cand.push_back(iChrgn[k]);

                	bool pip_cand = true;
                	for(int kk = 0; kk < ppip_cand.size(); kk++)
                		if(iChrgp[m] == ppip_cand[kk]) pip_cand = false;

                	if(pip_cand) ppip_cand.push_back(iChrgp[m]);

                	pip_cand = true;
                	for(int kk = 0; kk < ppip_cand.size(); kk++)
                		if(iChrgp[n] == ppip_cand[kk]) pip_cand = false;

                	if(pip_cand) ppip_cand.push_back(iChrgp[n]);
                }
			}
        }
    }
    if( dltmin2 > 998. )	return sc;
    Ncut3++;
    Neff++;
	efficiency->Fill(Neff); 


   m_npr_cand    = ppp_cand.size();
   m_nprbar_cand = ppm_cand.size();
   m_npim_cand   = ppim_cand.size();
   m_npip_cand   = ppip_cand.size();


    // go through the analysis once more with the best combination and top it off with a kinfit_4C using the XiXibar tracks.

    double declength_lam, declength_lambar, declength_xim, declength_xip;
	HepLorentzVector V4Pr_lam(0., 0., 0., 0.);
	HepLorentzVector V4Pim_lam(0., 0., 0., 0.);	
	HepLorentzVector V4Pim_xim(0., 0., 0., 0.);
	HepLorentzVector V4Lam(0., 0., 0., 0.);
	HepLorentzVector V4Xim(0., 0., 0., 0.);
	HepLorentzVector V4Prbar_lambar(0., 0., 0., 0.);
	HepLorentzVector V4Pip_lambar(0., 0., 0., 0.);
	HepLorentzVector V4Pip_xibarp(0., 0., 0., 0.);
	HepLorentzVector V4Lambar(0., 0., 0., 0.);
	HepLorentzVector V4Xipbar(0., 0., 0., 0.);

	//**************** best Xi- combination  *******************************

    RecMdcKalTrack *prfinalTrk  = (*(evtRecTrkCol->begin()+Ximtracks[0]))->mdcKalTrack();
	RecMdcKalTrack *pimfinalTrk = (*(evtRecTrkCol->begin()+Ximtracks[1]))->mdcKalTrack();

	RecMdcTrack *prvxTrk  = (*(evtRecTrkCol->begin()+Ximtracks[0]))->mdcTrack();
	RecMdcTrack *pimvxTrk  = (*(evtRecTrkCol->begin()+Ximtracks[1]))->mdcTrack();

	double time_best_proton = TimeOfFlight( Ximtracks[0] , 1);
	//std::cout << "time_best_proton " << time_best_proton << std::endl;
	hm_tof_proton_best->Fill(time_best_proton);

	double time_best_pion_1 = TimeOfFlight( Ximtracks[1] , 0);
	//std::cout << "time_best_pion_1 " << time_best_pion_1 << std::endl;
	hm_tof_pion_best1->Fill(time_best_pion_1);


	pr_rec_px = prfinalTrk->px(); pr_rec_py = prfinalTrk->py(); 
	pr_rec_pz = prfinalTrk->pz(); pr_rec_e  = TMath::Sqrt(prfinalTrk->p()*prfinalTrk->p() + xmass[4]*xmass[4]);

	pimLam_rec_px = pimfinalTrk->px(); pimLam_rec_py = pimfinalTrk->py(); 
	pimLam_rec_pz = pimfinalTrk->pz(); pimLam_rec_e  = TMath::Sqrt(pimfinalTrk->p()*pimfinalTrk->p() + xmass[2]*xmass[2]);	

	Rxy_pr 		= GetRxy(Ximtracks[0]);
	Rxy_pimlam 	= GetRxy(Ximtracks[1]);

	RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);
	WTrackParameter  wvprTrk;
	wvprTrk = WTrackParameter(mpp, prfinalTrk->helix(), prfinalTrk->err());

	HTrackParameter hvprTrk = HTrackParameter(wvprTrk);
	HepVector prv = hvprTrk.hel();
	pr_drho = prv[0]; pr_phi0 = prv[1]; pr_kappa = prv[2]; pr_lambda = prv[3];

	RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
	WTrackParameter  wvpimTrk;
	wvpimTrk = WTrackParameter(mpim, pimfinalTrk->helix(), pimfinalTrk->err());

	HTrackParameter hvpimTrk = HTrackParameter(wvpimTrk);
	HepVector pimlamv = hvpimTrk.hel();
	pimlam_drho = pimlamv[0]; pimlam_phi0 = pimlamv[1]; pimlam_kappa = pimlamv[2]; pimlam_lambda = pimlamv[3];

	//**************** lambda0 combination  *********************************
	VertexFit *vtxfitlam 			= VertexFit::instance();
	SecondVertexFit *svtxfitlam 	= SecondVertexFit::instance();

	vtxfitlam->init();
	vtxfitlam->AddTrack(0,  wvprTrk);
	vtxfitlam->AddTrack(1,  wvpimTrk);
	vtxfitlam->AddVertex(0, vxpar, 0, 1);
	if(!vtxfitlam->Fit(0)) return sc;
	lam_prvxfitchi2 = vtxfitlam->chisq(0);
	pr_vx_x = prvxTrk->x();   pr_vx_y  = prvxTrk->y();  pr_vx_z = prvxTrk->z();
	pim_vx_x = pimvxTrk->x(); pim_vx_y = pimvxTrk->y(); pim_vx_z = pimvxTrk->z();


	vtxfitlam->Swim(0);
	vtxfitlam->BuildVirtualParticle(0);
	wLambda = vtxfitlam->wVirtualTrack(0); 
	VertexParameter vtxlambda1 = vtxfitlam->vpar(0);
	HepPoint3D vtxlam = vtxlambda1.vx();
	Lam_vx_x = vtxlam.x(); Lam_vx_y = vtxlam.y(); Lam_vx_z = vtxlam.z(); 

	svtxfitlam->init();
	svtxfitlam->AddTrack(0, wLambda);
	svtxfitlam->setVpar(vtxlambda1);
	svtxfitlam->setPrimaryVertex(LambdaVertex);
	//svtxfitlam->setIterNumber(20); //default is 10 
	//svtxfitlam->setChisqCut(1000.); 

	if(!svtxfitlam->Fit()) return sc;
	lam_decl = svtxfitlam->decayLength();
	lam_decl_err = svtxfitlam->decayLengthError();
	lam_vxfitchi2= svtxfitlam->chisq();
	svtxfitlam->ctau();

	WTrackParameter	wvprTrk_fit = vtxfitlam->wtrk(0);
	HTrackParameter hvprTrk_fit = HTrackParameter(wvprTrk_fit);
	HepVector prv_fit = hvprTrk_fit.hel();
	prvfit_drho = prv_fit[0]; prvfit_phi0 = prv_fit[1]; prvfit_kappa = prv_fit[2]; prvfit_lambda = prv_fit[3];

	WTrackParameter	wvpimTrk_fit = vtxfitlam->wtrk(1);
	HTrackParameter hvpimTrk_fit = HTrackParameter(wvpimTrk_fit);
	HepVector pimv_fit = hvpimTrk_fit.hel();
	pimlamvfit_drho = pimv_fit[0]; pimlamvfit_phi0 = pimv_fit[1]; pimlamvfit_kappa = pimv_fit[2]; pimlamvfit_lambda = pimv_fit[3];


	//wlambda1  = svtxfitlam->wpar();
	// double chisq2 = svtxfitlam->chisq();

	// four vectors after secondary vertex fit. Also 4-vec of proton and pion may be updated
	V4Pr_lam  = vtxfitlam->pfit(0);
	V4Pim_lam = vtxfitlam->pfit(1);

	pr_vx_e 	= V4Pr_lam.e();		pr_vx_px = V4Pr_lam.px();	
	pr_vx_py 	= V4Pr_lam.py();	pr_vx_pz = V4Pr_lam.pz();
	// pi- (lam) final track
	pimLam_vx_e  = V4Pim_lam.e();	pimLam_vx_px = V4Pim_lam.px();	
	pimLam_vx_py = V4Pim_lam.py();	pimLam_vx_pz = V4Pim_lam.pz();

    //**************** Xi- vertex fitting  *********************************

	RecMdcKalTrack *pim2finalTrk = (*(evtRecTrkCol->begin()+Ximtracks[2]))->mdcKalTrack();
	Rxy_pimxim 	= GetRxy(Ximtracks[2]);

	double time_best_pion_2 = TimeOfFlight( Ximtracks[2] , 0);
	//std::cout << "time_best_pion_2 " << time_best_pion_2 << std::endl;
	hm_tof_pion_best2->Fill(time_best_pion_2);

	pimXim_rec_px = pim2finalTrk->px(); pimXim_rec_py = pim2finalTrk->py(); 
	pimXim_rec_pz = pim2finalTrk->pz(); pimXim_rec_e  = TMath::Sqrt(pim2finalTrk->p()*pim2finalTrk->p() + xmass[2]*xmass[2]);

	RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
	WTrackParameter wvpim2Trk;
	wvpim2Trk = WTrackParameter(mpim, pim2finalTrk->helix(), pim2finalTrk->err());

	HTrackParameter hvpim2Trk = HTrackParameter(wvpim2Trk);
	HepVector pimximv = hvpim2Trk.hel();
	pimxim_drho = pimximv[0]; pimxim_phi0 = pimximv[1]; pimxim_kappa = pimximv[2]; pimxim_lambda = pimximv[3];

	VertexFit *vtxfitXi 			= VertexFit::instance();
	SecondVertexFit *svtxfitXi 		= SecondVertexFit::instance();

	vtxfitXi->init();
    vtxfitXi->AddTrack(0, wLambda);
    vtxfitXi->AddTrack(1, wvpim2Trk);
    vtxfitXi->AddVertex(0, vxpar, 0, 1); 
    vtxfitXi->Fit(0);
    if(!vtxfitXi->Fit()) return sc;
    if(vtxfitXi->chisq(0)>200) return sc;

    vtxfitXi->BuildVirtualParticle(0);
    VertexParameter vparXim = vtxfitXi->vpar(0);
	HepPoint3D vtxxim = vparXim.vx();
	Xim_vx_x = vtxxim.x(); Xim_vx_y = vtxxim.y(); Xim_vx_z = vtxxim.z();
    WTrackParameter wparXim = vtxfitXi->wVirtualTrack(0);
    wparXim.setCharge(-1);

    svtxfitXi->init();
    svtxfitXi->AddTrack(0,wparXim);
    svtxfitXi->setVpar(vparXim);
    //svtxfitXi->setChisqCut(1000.);
    svtxfitXi->setPrimaryVertex(pvxpar);
    //svtxfitXi->setIterNumber(20); //default is 10 
    svtxfitXi->Fit();
    if(!svtxfitXi->Fit()) return sc;
    // if(svtxfitXi->chisq()>200) continue;

    xim_decl     = svtxfitXi->decayLength();
    xim_decl_err = svtxfitXi->decayLengthError();
    xim_vxfitchi2 = svtxfitlam->chisq();
    svtxfitXi->ctau();

    WTrackParameter	wvpimximTrk_fit = vtxfitXi->wtrk(1);
	HTrackParameter hvpimximTrk_fit = HTrackParameter(wvpimximTrk_fit);
	HepVector pimximv_fit = hvpimximTrk_fit.hel();
	pimximvfit_drho = pimximv_fit[0]; pimximvfit_phi0 = pimximv_fit[1]; pimximvfit_kappa = pimximv_fit[2]; pimximvfit_lambda = pimximv_fit[3];
    
	// four vectors after secondary vertex fit. Also 4-vec of Lambda and pion may be updated
	V4Lam     = vtxfitXi->pfit(0);
	V4Pim_xim = vtxfitXi->pfit(1);
	V4Xim  	  = svtxfitXi->p4par();

	Lam_vx_e  = V4Lam.e();		Lam_vx_px = V4Lam.px();
	Lam_vx_py = V4Lam.py();		Lam_vx_pz = V4Lam.pz();

	pimXim_vx_e  = V4Pim_xim.e();	pimXim_vx_px = V4Pim_xim.px();
	pimXim_vx_py = V4Pim_xim.py();	pimXim_vx_pz = V4Pim_xim.pz();

	Xim_vx_e  = V4Xim.e();		Xim_vx_px = V4Xim.px();
	Xim_vx_py = V4Xim.py();		Xim_vx_pz = V4Xim.pz();


	//**************** best Xi+ combination  *******************************

    RecMdcKalTrack *prbarfinalTrk   = (*(evtRecTrkCol->begin()+Xiptracks[0]))->mdcKalTrack();
	RecMdcKalTrack *pipfinalTrk 	= (*(evtRecTrkCol->begin()+Xiptracks[1]))->mdcKalTrack();

	RecMdcTrack *prbarvxTrk  = (*(evtRecTrkCol->begin()+Ximtracks[0]))->mdcTrack();
	RecMdcTrack *pipvxTrk    = (*(evtRecTrkCol->begin()+Ximtracks[1]))->mdcTrack();

	prbar_rec_px = prbarfinalTrk->px(); prbar_rec_py = prbarfinalTrk->py(); 
	prbar_rec_pz = prbarfinalTrk->pz(); prbar_rec_e  = TMath::Sqrt(prbarfinalTrk->p()*prbarfinalTrk->p() + xmass[4]*xmass[4]);

	pipLambar_rec_px = pipfinalTrk->px(); pipLambar_rec_py = pipfinalTrk->py(); 
	pipLambar_rec_pz = pipfinalTrk->pz(); pipLambar_rec_e  = TMath::Sqrt(pipfinalTrk->p()*pipfinalTrk->p() + xmass[2]*xmass[2]);	

	Rxy_prbar 			= GetRxy(Xiptracks[0]);
	Rxy_piplambar 		= GetRxy(Xiptracks[1]);

    RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);
	WTrackParameter  wvprbarTrk;
	wvprbarTrk = WTrackParameter(mpm, prbarfinalTrk->helix(), prbarfinalTrk->err());

	HTrackParameter hvprbarTrk = HTrackParameter(wvprbarTrk);
	HepVector prbarv = hvprbarTrk.hel();
	prbar_drho = prbarv[0]; prbar_phi0 = prbarv[1]; prbar_kappa = prbarv[2]; prbar_lambda = prbarv[3];
	
	RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
	WTrackParameter  wvpipTrk;
	wvpipTrk = WTrackParameter(mpip, pipfinalTrk->helix(), pipfinalTrk->err());

	HTrackParameter hvpipTrk = HTrackParameter(wvpipTrk);
	HepVector piplambarv = hvpipTrk.hel();
	piplambar_drho = piplambarv[0]; piplambar_phi0 = piplambarv[1]; piplambar_kappa = piplambarv[2]; piplambar_lambda = piplambarv[3];
	
	//**************** lambdabar0 combination  *********************************
	VertexFit 		*vtxfitlambar 		= VertexFit::instance();
	SecondVertexFit *svtxfitlambar 		= SecondVertexFit::instance();

	vtxfitlambar->init();
	vtxfitlambar->AddTrack(0,  wvprbarTrk);
	vtxfitlambar->AddTrack(1,  wvpipTrk);
	vtxfitlambar->AddVertex(0, vxpar, 0, 1);
	if(!vtxfitlambar->Fit(0)) return sc;

	lambar_prvxfitchi2= vtxfitlambar->chisq(0);
	prbar_vx_x = prbarvxTrk->x(); prbar_vx_y = prbarvxTrk->y(); prbar_vx_z = prbarvxTrk->z();
	pip_vx_x = pipvxTrk->x(); pip_vx_y = pipvxTrk->y(); pip_vx_z = pipvxTrk->z();

	vtxfitlambar->Swim(0);
	vtxfitlambar->BuildVirtualParticle(0);
	wLambdabar = vtxfitlambar->wVirtualTrack(0); 
	VertexParameter vtxlambdabar1 = vtxfitlambar->vpar(0);
	HepPoint3D vtxlambar = vtxlambdabar1.vx();
	Lambar_vx_x = vtxlambar.x(); Lambar_vx_y = vtxlambar.y(); Lambar_vx_z = vtxlambar.z(); 

	svtxfitlambar->init();
	svtxfitlambar->AddTrack(0, wLambdabar);
	svtxfitlambar->setVpar(vtxlambdabar1);
	svtxfitlambar->setPrimaryVertex(LambdabarVertex);
    //svtxfitlambar->setIterNumber(20); //default is 10 
	if(!svtxfitlambar->Fit()) return sc;
	lambar_decl = svtxfitlambar->decayLength();
	lambar_decl_err = svtxfitlam->decayLengthError();
	lambar_vxfitchi2 = svtxfitlam->chisq();
	svtxfitlambar->ctau();

	WTrackParameter	wvprbarTrk_fit = vtxfitlambar->wtrk(0);
	HTrackParameter hvprbarTrk_fit = HTrackParameter(wvprbarTrk_fit);
	HepVector prbarv_fit = hvprbarTrk_fit.hel();
	prbarvfit_drho = prbarv_fit[0]; prbarvfit_phi0 = prbarv_fit[1]; prbarvfit_kappa = prbarv_fit[2]; prbarvfit_lambda = prbarv_fit[3];

	WTrackParameter	wvpipTrk_fit = vtxfitlambar->wtrk(1);
	HTrackParameter hvpipTrk_fit = HTrackParameter(wvpipTrk_fit);
	HepVector pipv_fit = hvpipTrk_fit.hel();
	piplambarvfit_drho = pipv_fit[0]; piplambarvfit_phi0 = pipv_fit[1]; piplambarvfit_kappa = pipv_fit[2]; piplambarvfit_lambda = pipv_fit[3];

	// four vectors after secondary vertex fit. Also 4-vec of proton and pion may be updated
	V4Prbar_lambar  = vtxfitlambar->pfit(0);
	V4Pip_lambar 	= vtxfitlambar->pfit(1);
	
	prbar_vx_e 		= V4Prbar_lambar.e();	prbar_vx_px = V4Prbar_lambar.px();	
	prbar_vx_py 	= V4Prbar_lambar.py();	prbar_vx_pz = V4Prbar_lambar.pz();
	// pi+ (lambar) final track
	pipLambar_vx_e  = V4Pip_lambar.e();		pipLambar_vx_px = V4Pip_lambar.px();	
	pipLambar_vx_py = V4Pip_lambar.py();	pipLambar_vx_pz = V4Pip_lambar.pz();

    //**************** Xi+ vertex fitting  *********************************
	RecMdcKalTrack *pip2finalTrk = (*(evtRecTrkCol->begin()+Xiptracks[2]))->mdcKalTrack();
	Rxy_pipxip 		= GetRxy(Xiptracks[2]);

	pipXip_rec_px = pip2finalTrk->px(); pipXip_rec_py = pip2finalTrk->py(); 
	pipXip_rec_pz = pip2finalTrk->pz(); pipXip_rec_e  = TMath::Sqrt(pip2finalTrk->p()*pip2finalTrk->p() + xmass[2]*xmass[2]);

	RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
	WTrackParameter  wvpip2Trk;
	wvpip2Trk = WTrackParameter(xmass[2], pip2finalTrk->helix(), pip2finalTrk->err());

	HTrackParameter hvpip2Trk = HTrackParameter(wvpip2Trk);
	HepVector pipxipv = hvpip2Trk.hel();
	pipxip_drho = pipxipv[0]; pipxip_phi0 = pipxipv[1]; pipxip_kappa = pipxipv[2]; pipxip_lambda = pipxipv[3];
    
	VertexFit *vtxfitXibar 				= VertexFit::instance();
	SecondVertexFit *svtxfitXibar 		= SecondVertexFit::instance();

	vtxfitXibar->init();
    vtxfitXibar->AddTrack(0, wLambdabar);
    vtxfitXibar->AddTrack(1, wvpip2Trk);
    vtxfitXibar->AddVertex(0, vxpar, 0, 1); 
    vtxfitXibar->Fit(0);
    if(!vtxfitXibar->Fit()) return sc;
    if(vtxfitXibar->chisq(0)>200) return sc;

    vtxfitXibar->BuildVirtualParticle(0);
    VertexParameter vparXip = vtxfitXibar->vpar(0);
    HepPoint3D vtxxip = vparXip.vx();
    Xip_vx_x = vtxxip.x(); Xip_vx_y = vtxxip.y(); Xip_vx_z = vtxxip.z();
    WTrackParameter wparXip = vtxfitXibar->wVirtualTrack(0);
    wparXip.setCharge(+1);

    svtxfitXibar->init();
    svtxfitXibar->AddTrack(0,wparXip);
    svtxfitXibar->setVpar(vparXip);
    //svtxfitXibar->setChisqCut(1000.);
    svtxfitXibar->setPrimaryVertex(pvxpar);
    //svtxfitXibar->setIterNumber(20); //default is 10 
    svtxfitXibar->Fit();
    if(!svtxfitXibar->Fit()) return sc;
    // if(svtxfitXi->chisq()>200) continue;

    xip_decl     = svtxfitXibar->decayLength();
    xip_decl_err = svtxfitXibar->decayLengthError(); 
	xip_vxfitchi2 = svtxfitXibar->chisq();
    svtxfitXibar->ctau();

    WTrackParameter	wvpipxipTrk_fit = vtxfitXibar->wtrk(1);
	HTrackParameter hvpipxipTrk_fit = HTrackParameter(wvpipxipTrk_fit);
	HepVector pipxipv_fit = hvpipxipTrk_fit.hel();
	pipxipvfit_drho = pipxipv_fit[0]; pipxipvfit_phi0 = pipxipv_fit[1]; pipxipvfit_kappa = pipxipv_fit[2]; pipxipvfit_lambda = pipxipv_fit[3];

	// four vectors after secondary vertex fit. Also 4-vec of Lambda and pion may be updated
	V4Lambar     = vtxfitXibar->pfit(0);
	V4Pip_xibarp = vtxfitXibar->pfit(1);
	V4Xipbar  	 = svtxfitXibar->p4par();

	Lambar_vx_e  = V4Lambar.e();		Lambar_vx_px = V4Lambar.px();
	Lambar_vx_py = V4Lambar.py();		Lambar_vx_pz = V4Lambar.pz();

	pipXip_vx_e  = V4Pip_xibarp.e();	pipXip_vx_px = V4Pip_xibarp.px();
	pipXip_vx_py = V4Pip_xibarp.py();	pipXip_vx_pz = V4Pip_xibarp.pz();

	Xip_vx_e  = V4Xipbar.e();		Xip_vx_px = V4Xipbar.px();
	Xip_vx_py = V4Xipbar.py();		Xip_vx_pz = V4Xipbar.pz();


    // particle 
    double mom_prcand 		= TMath::Sqrt(V4Pr_lam.px()*V4Pr_lam.px()   + V4Pr_lam.py()*V4Pr_lam.py()   + V4Pr_lam.pz()*V4Pr_lam.pz());
    double mom_pim_lam_cand = TMath::Sqrt(V4Pim_lam.px()*V4Pim_lam.px() + V4Pim_lam.py()*V4Pim_lam.py() + V4Pim_lam.pz()*V4Pim_lam.pz());
    double mom_pim_Xim_cand = TMath::Sqrt(V4Pim_xim.px()*V4Pim_xim.px() + V4Pim_xim.py()*V4Pim_xim.py() + V4Pim_xim.pz()*V4Pim_xim.pz());

	hp_pr_lam_vx->Fill(mom_prcand);
    hp_pim_lam_vx->Fill(mom_pim_lam_cand);
    hp_pim_Xim_vx->Fill(mom_pim_Xim_cand);

	hpvth_pr_lam_vx->Fill(mom_prcand, V4Pr_lam.theta()*180./3.14);
    hpvth_pim_lam_vx->Fill(mom_pim_lam_cand, V4Pim_lam.theta()*180./3.14);
    hpvth_pim_Xim_vx->Fill(mom_pim_Xim_cand, V4Pim_xim.theta()*180./3.14);

    // anti-particle 
    double mom_prbarcand 	= TMath::Sqrt(V4Prbar_lambar.px()*V4Prbar_lambar.px() + V4Prbar_lambar.py()*V4Prbar_lambar.py() + V4Prbar_lambar.pz()*V4Prbar_lambar.pz());
    double mom_pip_lam_cand = TMath::Sqrt(V4Pip_lambar.px()*V4Pip_lambar.px() + V4Pip_lambar.py()*V4Pip_lambar.py() + V4Pip_lambar.pz()*V4Pip_lambar.pz());
	double mom_pip_Xip_cand = TMath::Sqrt(V4Pip_xibarp.px()*V4Pip_xibarp.px() + V4Pip_xibarp.py()*V4Pip_xibarp.py() + V4Pip_xibarp.pz()*V4Pip_xibarp.pz());

	hp_pr_bar_lambar_vx->Fill(mom_prbarcand);
	hp_pip_lambar_vx->Fill(mom_pip_lam_cand);
    hp_pip_Xipbar_vx->Fill(mom_pip_Xip_cand);

    hpvth_pr_bar_lambar_vx->Fill(mom_prbarcand, V4Prbar_lambar.theta()*180./3.14);
    hpvth_pip_lambar_vx->Fill(mom_pip_lam_cand, V4Pip_lambar.theta()*180./3.14);
    hpvth_pip_Xip_vx->Fill(mom_pip_Xip_cand, V4Pip_xibarp.theta()*180./3.14);


    HepLorentzVector jpsi_vx;
    jpsi_vx = V4Xim + V4Xipbar;


    hm_lam_vx->Fill(V4Lam.m());
    hm_lam_bar_vx->Fill(V4Lambar.m());
    hm_xim_vx->Fill(V4Xim.m());
    hm_xip_bar_vx->Fill(V4Xipbar.m());
    hm_jpsi_vx->Fill(jpsi_vx.m());

    hm_lam_lambar_vx_final->Fill(V4Lam.m(), V4Lambar.m());
    hm_xim_xipbar_vx->Fill(V4Xim.m(), V4Xipbar.m());


    std::vector<WTrackParameter> tracks;
    tracks.push_back(wparXim);
    tracks.push_back(wparXip);

    TString constraint = "4C";
    KinematicFit(tracks, constraint);
    Neff++;
    efficiency->Fill(Neff);

    if(eventHeader->runNumber() < 0){

        true_Xim_e 	= HepXiTrue.e();true_Xim_px 	= HepXiTrue.px(); 
		true_Xim_py = HepXiTrue.py();true_Xim_pz 	= HepXiTrue.pz();

		true_Xip_e 	= HepXibarTrue.e();true_Xip_px 	= HepXibarTrue.px();
		true_Xip_py = HepXibarTrue.py();true_Xip_pz = HepXibarTrue.pz();	            

		true_Lam_e 	= HepLambdaTrue.e();true_Lam_px 	= HepLambdaTrue.px();
		true_Lam_py = HepLambdaTrue.py();true_Lam_pz 	= HepLambdaTrue.pz();				

		true_Lambar_e 	= HepLambdabarTrue.e();true_Lambar_px 	= HepLambdabarTrue.px();
		true_Lambar_py 	= HepLambdabarTrue.py();true_Lambar_pz 	= HepLambdabarTrue.pz();					

		true_pr_e 	= HepProtonTrue.e();true_pr_px 	= HepProtonTrue.px();
		true_pr_py 	= HepProtonTrue.py();true_pr_pz = HepProtonTrue.pz();	

		true_prbar_e 	= HepProtonbarTrue.e();true_prbar_px 	= HepProtonbarTrue.px();
		true_prbar_py 	= HepProtonbarTrue.py();true_prbar_pz 	= HepProtonbarTrue.pz();	

		true_pipXip_e	= HepPipXibarTrue.e();true_pipXip_px 	= HepPipXibarTrue.px();
		true_pipXip_py 	= HepPipXibarTrue.py();true_pipXip_pz 	= HepPipXibarTrue.pz();

		true_pipLambar_e	= HepPipLambarTrue.e();true_pipLambar_px 	= HepPipLambarTrue.px();
		true_pipLambar_py 	= HepPipLambarTrue.py();true_pipLambar_pz 	= HepPipLambarTrue.pz();

		true_pimXim_e	= HepPimXiTrue.e();true_pimXim_px 	= HepPimXiTrue.px();
		true_pimXim_py 	= HepPimXiTrue.py();true_pimXim_pz 	= HepPimXiTrue.pz();

		true_pimLam_e	= HepPimLamTrue.e();true_pimLam_px 	= HepPimLamTrue.px();
		true_pimLam_py 	= HepPimLamTrue.py();true_pimLam_pz 	= HepPimLamTrue.pz();			

    }
    Ncut7++;
    if(chi2_4C  > 0)

    m_tuple1->write();

	return StatusCode::SUCCESS;
} 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode XiXibar::finalize() 
{

	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	cout<<"total number:         "<<Ncut0<<endl;
	cout<<"iChrgp.size()>2 && iChrgn.size()>2:   "<<Ncut1<<endl;
	cout<<"Lambda reconstruction:"<<Ncut3<<endl;
	cout<<"pr/anti-pr mom. cut:  "<<Ncut4<<endl;
	cout<<"Xi pi-/pi+ mom. cut:  "<<Ncut5<<endl;
	cout<<"Kinfit 4C:  			 "<<Ncut6<<endl;
	cout<<"No Kinfit 4C (but written):  		 "<<Ncut7<<endl;
	

	//--------------------------------------------------------------
	return StatusCode::SUCCESS;
}

void XiXibar::KinematicFit(std::vector<WTrackParameter> tracks, TString constraint){
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();

	kmfit->init();

	kmfit->AddTrack(0, tracks[0]);
	kmfit->AddTrack(1, tracks[1]);

	kmfit->AddFourMomentum(0, GetEcmsP4());
	kmfit->setChisqCut(500.);
	if(constraint == "4C"){
		chi2_4C = -9999;
	}
	else{
		std::cout << "No correct constraint in kinematic fit given!" << std::endl;
		return;
	}

	bool oksq = kmfit->Fit();
	if(oksq){

		HepLorentzVector kfit_xim 		= kmfit->pfit(0);
		HepLorentzVector kfit_xibarp 	= kmfit->pfit(1);

		Xim_4C_e  = kfit_xim.e();		Xim_4C_px = kfit_xim.px();
		Xim_4C_py = kfit_xim.py();		Xim_4C_pz = kfit_xim.pz();

		Xip_4C_e  = kfit_xibarp.e();	Xip_4C_px = kfit_xibarp.px();
		Xip_4C_py = kfit_xibarp.py();	Xip_4C_pz = kfit_xibarp.pz();


		HepVector kfit_xim_pull 		= kmfit->pull(0);
		HepVector kfit_xibarp_pull 		= kmfit->pull(1);

		Xim_pull_x = kfit_xim_pull(0);
		Xim_pull_y = kfit_xim_pull(1);		Xim_pull_z = kfit_xim_pull(2);

		Xip_pull_x = kfit_xibarp_pull(0);
		Xip_pull_y = kfit_xibarp_pull(1);	Xip_pull_z = kfit_xibarp_pull(2);

		if(constraint == "4C"){
			Ncut6++;
			chi2_4C = kmfit->chisq();

         	hm_xim_xipbar_fit->Fill(kfit_xim.m(), kfit_xibarp.m());
		}
	}
	return;
}

Double_t XiXibar::GetRxy(int itr){

   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),EventModel::EvtRec::EvtRecTrackCol);

   Hep3Vector xorigin(0,0,0);
   Hep3Vector eorigin(0,0,0);
   IVertexDbSvc*  vtxsvc;
   Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
   if(vtxsvc->isVertexValid()){
      double* dbv = vtxsvc->PrimaryVertex(); 
      double*  vv = vtxsvc->SigmaPrimaryVertex();
      xorigin.setX(dbv[0]);
      xorigin.setY(dbv[1]);
      xorigin.setZ(dbv[2]);
      eorigin.setX(vv[0]);
      eorigin.setY(vv[1]);
      eorigin.setZ(vv[2]);
    }

	RecMdcTrack *Trk  = (*(evtRecTrkCol->begin()+itr))->mdcTrack();
	double x0 = Trk->x();
	double y0 = Trk->y();
	double z0 = Trk->z();
	double r0 = Trk->r();
	double phi0 = Trk->helix(1);

	double xv=xorigin.x();
	double yv=xorigin.y();
	double zv=xorigin.z();
	double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);

	return Rxy;
}

double XiXibar::TimeOfFlight(int itrack, int itype ){
	//if (itype == 0) pion
	//if (itype == 1) proton

	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),EventModel::EvtRec::EvtRecTrackCol);

	double tof_time = -100.;
	double m_tpi_etof, m_tp_etof, m_tpi_etof1, m_tpi_etof2, m_tp_etof1, m_tp_etof2;
	double velc = 299.792458;
    
	EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + itrack;
	if(!(*itTrk)->isMdcTrackValid()) return tof_time;
	if(!(*itTrk)->isTofTrackValid()) return tof_time;
	RecMdcTrack * mdcTrk = (*itTrk)->mdcTrack();
	SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
	
	double ptrk = mdcTrk->p();
	
	for (SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin(); iter_tof != tofTrkCol.end(); iter_tof++) {
	  TofHitStatus *status = new TofHitStatus;
	  status->setStatus((*iter_tof)->status());
	  if(!(status->is_barrel())){
		  if(!(status->is_counter())) continue;
		  if(status->layer()!=0) continue;
		  double path = (*iter_tof)->path();
		  double tof = (*iter_tof)->tof();
		  double ph = (*iter_tof)->ph();
		  double rhit = (*iter_tof)->zrhit();
		  double qual = 0.0 + (*iter_tof)->quality();
		  double cntr = 0.0 + (*iter_tof)->tofID();
		  double texp[5];
		  for (int j= 0; j <5; j++) {
			  double gb = ptrk/xmass[j];
			  double beta = gb/sqrt(1+gb*gb);
			  texp[j] = 10*path/beta/velc;
		  }
		  // m_cntr_etof = cntr;
		  // m_ptot_etof = ptrk;
		  // m_ph_etof = ph;
		  // m_rhit_etof = rhit;
		  // m_qual_etof = qual;
		  // m_te_etof   = tof - texp[0];
		  // m_tmu_etof  = tof - texp[1];
		  m_tpi_etof  = tof - texp[2];
		  // m_tk_etof   = tof - texp[3];
		  m_tp_etof   = tof - texp[4];

		  if(itype == 0)
		  	tof_time = m_tpi_etof;
		  else if(itype == 1)
		  	tof_time = m_tp_etof;
	  }
	  else {
		  if ( !(status->is_counter())) continue;
		  if (status->layer() ==1) {
			  double path = (*iter_tof)->path();
			  double tof = (*iter_tof)->tof();
			  double ph = (*iter_tof)->ph();
			  double rhit = (*iter_tof)->zrhit();
			  double qual = 0.0 + (*iter_tof)->quality();
			  double cntr = 0.0 + (*iter_tof)->tofID();
			  double texp[5];
			  for (int j= 0; j <5; j++) {
				  double gb = ptrk/xmass[j];
				  double beta = gb/sqrt(1+gb*gb);
				  texp[j] = 10*path/beta/velc;
			  }
			  // m_cntr_etof1 = cntr;
			  // m_ptot_etof1 = ptrk;
			  // m_ph_etof1 = ph;
			  // m_rhit_etof1 = rhit;
			  // m_qual_etof1 = qual;
			  // m_te_etof1   = tof - texp[0];
			  // m_tmu_etof1  = tof - texp[1];
			  m_tpi_etof1  = tof - texp[2];
			  // m_tk_etof1   = tof - texp[3];
			  m_tp_etof1   = tof - texp[4];

			  if(itype == 0)
		  		tof_time = m_tpi_etof1;
		  	  else if(itype == 1)
		  		tof_time = m_tp_etof1;
		  }
	  }
	  if(status->layer()==2) {
		  double path = (*iter_tof)->path();
		  double tof = (*iter_tof)->tof();
		  double ph = (*iter_tof)->ph();
		  double rhit = (*iter_tof)->zrhit();
		  double qual = 0.0 + (*iter_tof)->quality();
		  double cntr = 0.0 + (*iter_tof)->tofID();
		  double texp[5];
		  for (int j= 0; j <5; j++) {
			  double gb = ptrk/xmass[j];
			  double beta = gb/sqrt(1+gb*gb);
			  texp[j] = 10*path/beta/velc;
		  }
		  // m_cntr_etof2 = cntr;
		  // m_ptot_etof2 = ptrk;
		  // m_ph_etof2 = ph;
		  // m_rhit_etof2 = rhit;
		  // m_qual_etof2 = qual;
		  // m_te_etof2 = tof - texp[0];
		  // m_tmu_etof2 = tof - texp[1];
		  m_tpi_etof2  = tof - texp[2];
		  // m_tk_etof2 = tof - texp[3];
		  m_tp_etof2 = tof - texp[4];

		  if(itype == 0)
	  		tof_time = m_tpi_etof2;
	  	  else if(itype == 1)
	  		tof_time = m_tp_etof2;

	  }
	  delete status;
	}

	return tof_time;


}


StatusCode XiXibar::initialiezeNtuple(){
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;

	NTuplePtr nt0(ntupleSvc(), "FILE1/mctrue");
	
	if ( nt0 ) m_tuple = nt0;
	else{
		m_tuple = ntupleSvc()->book ("FILE1/mctrue", CLID_ColumnWiseTuple,
				"mctrue");
		if( m_tuple ){   
       		status = m_tuple->addItem ("mc_idxmc",  	mc_idxmc, 0, 100 );
        	status = m_tuple->addItem ("mc_pdgid",  	mc_idxmc, mc_pdgid  );
        	status = m_tuple->addItem ("mc_motheridx", 	mc_idxmc, mc_motheridx  );
			status = m_tuple->addItem ("mc_Xim_e",  	mc_Xim_e );
			status = m_tuple->addItem ("mc_Xim_px",  	mc_Xim_px );
			status = m_tuple->addItem ("mc_Xim_py",  	mc_Xim_py );
			status = m_tuple->addItem ("mc_Xim_pz",  	mc_Xim_pz );
			status = m_tuple->addItem ("mc_Xip_e",  	mc_Xip_e );
			status = m_tuple->addItem ("mc_Xip_px",  	mc_Xip_px );
			status = m_tuple->addItem ("mc_Xip_py",  	mc_Xip_py );
			status = m_tuple->addItem ("mc_Xip_pz",  	mc_Xip_pz );
			status = m_tuple->addItem ("mc_Lam_e",  	mc_Lam_e );
			status = m_tuple->addItem ("mc_Lam_px",  	mc_Lam_px );
			status = m_tuple->addItem ("mc_Lam_py",  	mc_Lam_py );
			status = m_tuple->addItem ("mc_Lam_pz",  	mc_Lam_pz );
			status = m_tuple->addItem ("mc_Lambar_e",  	mc_Lambar_e );
			status = m_tuple->addItem ("mc_Lambar_px",  mc_Lambar_px );
			status = m_tuple->addItem ("mc_Lambar_py",  mc_Lambar_py );
			status = m_tuple->addItem ("mc_Lambar_pz",  mc_Lambar_pz );
			status = m_tuple->addItem ("mc_pr_e",  		mc_pr_e );
			status = m_tuple->addItem ("mc_pr_px",  	mc_pr_px );
			status = m_tuple->addItem ("mc_pr_py",  	mc_pr_py );
			status = m_tuple->addItem ("mc_pr_pz",  	mc_pr_pz );	
			status = m_tuple->addItem ("mc_prbar_e",  	mc_prbar_e );
			status = m_tuple->addItem ("mc_prbar_px",  	mc_prbar_px );
			status = m_tuple->addItem ("mc_prbar_py",  	mc_prbar_py );
			status = m_tuple->addItem ("mc_prbar_pz",  	mc_prbar_pz );
			status = m_tuple->addItem ("mc_pimLam_e",  	mc_pimLam_e );
			status = m_tuple->addItem ("mc_pimLam_px",  mc_pimLam_px );
			status = m_tuple->addItem ("mc_pimLam_py",  mc_pimLam_py );
			status = m_tuple->addItem ("mc_pimLam_pz",  mc_pimLam_pz );
			status = m_tuple->addItem ("mc_pipLambar_e",  mc_pipLambar_e );
			status = m_tuple->addItem ("mc_pipLambar_px", mc_pipLambar_px );
			status = m_tuple->addItem ("mc_pipLambar_py", mc_pipLambar_py );
			status = m_tuple->addItem ("mc_pipLambar_pz", mc_pipLambar_pz );
			status = m_tuple->addItem ("mc_pimXim_e",  	mc_pimXim_e );
			status = m_tuple->addItem ("mc_pimXim_px",  mc_pimXim_px );
			status = m_tuple->addItem ("mc_pimXim_py",  mc_pimXim_py );
			status = m_tuple->addItem ("mc_pimXim_pz",  mc_pimXim_pz );
			status = m_tuple->addItem ("mc_pipXip_e",  	mc_pipXip_e );
			status = m_tuple->addItem ("mc_pipXip_px",  mc_pipXip_px );
			status = m_tuple->addItem ("mc_pipXip_py",  mc_pipXip_py );
			status = m_tuple->addItem ("mc_pipXip_pz",  mc_pipXip_pz );								
		}
		else{ 
			log << MSG::ERROR << "Cannot book N-tuple:" << long(m_tuple) <<endmsg;
//			return StatusCode::FAILURE;
		}
	}

	NTuplePtr nt1(ntupleSvc(), "FILE1/recobs");
	if ( nt1 ) m_tuple1 = nt1;
	else {
		m_tuple1 = ntupleSvc()->book ("FILE1/recobs", CLID_ColumnWiseTuple, "Reconstructed observables");
		if( m_tuple1 ){

			status = m_tuple1->addItem("indexmc",m_idxmc, 0, 100);
			status = m_tuple1->addIndexedItem("pdgid", m_idxmc, m_pdgid );
			status = m_tuple1->addIndexedItem("motheridx",m_idxmc, m_motheridx );
			status = m_tuple1->addItem ("run",  m_run );
			status = m_tuple1->addItem ("rec",  m_rec );
			status = m_tuple1->addItem ("nTrkm", m_nTrkm, 2, 10);
			status = m_tuple1->addItem ("nTrkp", m_nTrkp, 2, 10);

			status = m_tuple1->addItem ("npr_cand", m_npr_cand, 2, 10);
			status = m_tuple1->addItem ("nprbar_cand", m_nprbar_cand, 2, 10);
			status = m_tuple1->addItem ("npim_cand", m_npim_cand, 2, 10);
			status = m_tuple1->addItem ("npip_cand", m_npip_cand, 2, 10);

			status = m_tuple1->addItem("lam_decl", lam_decl);
        	status = m_tuple1->addItem("lambar_decl", lambar_decl);
        	status = m_tuple1->addItem("xim_decl", xim_decl);
        	status = m_tuple1->addItem("xip_decl", xip_decl);
        	status = m_tuple1->addItem("lam_decl_err", lam_decl_err);
        	status = m_tuple1->addItem("lambar_decl_err", lambar_decl_err);
        	status = m_tuple1->addItem("xim_decl_err", xim_decl_err);
        	status = m_tuple1->addItem("xip_decl_err", xip_decl_err);
        	status = m_tuple1->addItem("lam_prvxfitchi2", lam_prvxfitchi2);
        	status = m_tuple1->addItem("lambar_prvxfitchi2", lambar_prvxfitchi2);        	
        	status = m_tuple1->addItem("lam_vxfitchi2", lam_vxfitchi2);
        	status = m_tuple1->addItem("lambar_vxfitchi2", lambar_vxfitchi2);
        	status = m_tuple1->addItem("xim_vxfitchi2", xim_vxfitchi2);
        	status = m_tuple1->addItem("xip_vxfitchi2", xip_vxfitchi2);
			
			status = m_tuple1->addItem("pr_vx_e", pr_vx_e);
        	status = m_tuple1->addItem("pr_vx_px",pr_vx_px);
        	status = m_tuple1->addItem("pr_vx_py", pr_vx_py);
        	status = m_tuple1->addItem("pr_vx_pz", pr_vx_pz);
        	
        	status = m_tuple1->addItem("prbar_vx_e", prbar_vx_e);
        	status = m_tuple1->addItem("prbar_vx_px", prbar_vx_px);
        	status = m_tuple1->addItem("prbar_vx_py", prbar_vx_py);
        	status = m_tuple1->addItem("prbar_vx_pz", prbar_vx_pz);
        	
        	status = m_tuple1->addItem("pimLam_vx_e", pimLam_vx_e);
        	status = m_tuple1->addItem("pimLam_vx_px", pimLam_vx_px);
        	status = m_tuple1->addItem("pimLam_vx_py", pimLam_vx_py);
        	status = m_tuple1->addItem("pimLam_vx_pz", pimLam_vx_pz);

        	status = m_tuple1->addItem("pipLambar_vx_e", pipLambar_vx_e);
        	status = m_tuple1->addItem("pipLambar_vx_px", pipLambar_vx_px);
        	status = m_tuple1->addItem("pipLambar_vx_py", pipLambar_vx_py);
        	status = m_tuple1->addItem("pipLambar_vx_pz", pipLambar_vx_pz);

        	status = m_tuple1->addItem("pimXim_vx_e", pimXim_vx_e);
        	status = m_tuple1->addItem("pimXim_vx_px", pimXim_vx_px);
        	status = m_tuple1->addItem("pimXim_vx_py", pimXim_vx_py);
        	status = m_tuple1->addItem("pimXim_vx_pz", pimXim_vx_pz);

        	status = m_tuple1->addItem("pipXip_vx_e", pipXip_vx_e);
        	status = m_tuple1->addItem("pipXip_vx_px", pipXip_vx_px);
        	status = m_tuple1->addItem("pipXip_vx_py", pipXip_vx_py);
        	status = m_tuple1->addItem("pipXip_vx_pz", pipXip_vx_pz);

        	status = m_tuple1->addItem("Lam_vx_e", Lam_vx_e);
        	status = m_tuple1->addItem("Lam_vx_px", Lam_vx_px);
        	status = m_tuple1->addItem("Lam_vx_py", Lam_vx_py);
        	status = m_tuple1->addItem("Lam_vx_pz", Lam_vx_pz);

        	status = m_tuple1->addItem("Lambar_vx_e", Lambar_vx_e);
        	status = m_tuple1->addItem("Lambar_vx_px", Lambar_vx_px);
        	status = m_tuple1->addItem("Lambar_vx_py", Lambar_vx_py);
        	status = m_tuple1->addItem("Lambar_vx_pz", Lambar_vx_pz);

        	status = m_tuple1->addItem("Xim_vx_e", Xim_vx_e);
        	status = m_tuple1->addItem("Xim_vx_px", Xim_vx_px);
        	status = m_tuple1->addItem("Xim_vx_py", Xim_vx_py);
        	status = m_tuple1->addItem("Xim_vx_pz", Xim_vx_pz);

        	status = m_tuple1->addItem("Xip_vx_e", Xip_vx_e);
        	status = m_tuple1->addItem("Xip_vx_px", Xip_vx_px);
        	status = m_tuple1->addItem("Xip_vx_py", Xip_vx_py);
        	status = m_tuple1->addItem("Xip_vx_pz", Xip_vx_pz);

        	status = m_tuple1->addItem("Xim_4C_e", Xim_4C_e);
        	status = m_tuple1->addItem("Xim_4C_px", Xim_4C_px);
        	status = m_tuple1->addItem("Xim_4C_py", Xim_4C_py);
        	status = m_tuple1->addItem("Xim_4C_pz", Xim_4C_pz);

        	status = m_tuple1->addItem("Xip_4C_e", Xip_4C_e);
        	status = m_tuple1->addItem("Xip_4C_px", Xip_4C_px);
        	status = m_tuple1->addItem("Xip_4C_py", Xip_4C_py);
        	status = m_tuple1->addItem("Xip_4C_pz", Xip_4C_pz); 

        	status = m_tuple1->addItem("pr_vx_x", pr_vx_x);
        	status = m_tuple1->addItem("pr_vx_y", pr_vx_y);
        	status = m_tuple1->addItem("pr_vx_z", pr_vx_z);

        	status = m_tuple1->addItem("prbar_vx_x", prbar_vx_x);
        	status = m_tuple1->addItem("prbar_vx_y", prbar_vx_y);
        	status = m_tuple1->addItem("prbar_vx_z", prbar_vx_z);

        	status = m_tuple1->addItem("pim_vx_x", pim_vx_x);
        	status = m_tuple1->addItem("pim_vx_y", pim_vx_y);
        	status = m_tuple1->addItem("pim_vx_z", pim_vx_z);

        	status = m_tuple1->addItem("pip_vx_x", pip_vx_x);
        	status = m_tuple1->addItem("pip_vx_y", pip_vx_y);
        	status = m_tuple1->addItem("pip_vx_z", pip_vx_z);

        	status = m_tuple1->addItem("pr_rec_e", pr_rec_e);
        	status = m_tuple1->addItem("pr_rec_px", pr_rec_px);
	        status = m_tuple1->addItem("pr_rec_py", pr_rec_py);
	        status = m_tuple1->addItem("pr_rec_pz", pr_rec_pz);

	        status = m_tuple1->addItem("prbar_rec_e", prbar_rec_e);
	        status = m_tuple1->addItem("prbar_rec_px", prbar_rec_px);
	        status = m_tuple1->addItem("prbar_rec_py", prbar_rec_py);
	        status = m_tuple1->addItem("prbar_rec_pz", prbar_rec_pz);

	        status = m_tuple1->addItem("pimLam_rec_e", pimLam_rec_e);
	        status = m_tuple1->addItem("pimLam_rec_px", pimLam_rec_px);
	        status = m_tuple1->addItem("pimLam_rec_py", pimLam_rec_py);
	        status = m_tuple1->addItem("pimLam_rec_pz", pimLam_rec_pz);

	        status = m_tuple1->addItem("pipLambar_rec_e", pipLambar_rec_e);
	        status = m_tuple1->addItem("pipLambar_rec_px", pipLambar_rec_px);
	        status = m_tuple1->addItem("pipLambar_rec_py", pipLambar_rec_py);
	        status = m_tuple1->addItem("pipLambar_rec_pz", pipLambar_rec_pz);

	        status = m_tuple1->addItem("pimXim_rec_e", pimXim_rec_e);
	        status = m_tuple1->addItem("pimXim_rec_px", pimXim_rec_px);
	        status = m_tuple1->addItem("pimXim_rec_py", pimXim_rec_py);
	        status = m_tuple1->addItem("pimXim_rec_pz", pimXim_rec_pz);

	        status = m_tuple1->addItem("pipXip_rec_e", pipXip_rec_e);
	        status = m_tuple1->addItem("pipXip_rec_px", pipXip_rec_px);
	        status = m_tuple1->addItem("pipXip_rec_py", pipXip_rec_py);
	        status = m_tuple1->addItem("pipXip_rec_pz", pipXip_rec_pz);

        	status = m_tuple1->addItem("Lam_vx_x", Lam_vx_x);
        	status = m_tuple1->addItem("Lam_vx_y", Lam_vx_y);
        	status = m_tuple1->addItem("Lam_vx_z", Lam_vx_z);

        	status = m_tuple1->addItem("Lambar_vx_x", Lambar_vx_x);
        	status = m_tuple1->addItem("Lambar_vx_y", Lambar_vx_y);
        	status = m_tuple1->addItem("Lambar_vx_z", Lambar_vx_z);

        	status = m_tuple1->addItem("Xim_vx_x", Xim_vx_x);
        	status = m_tuple1->addItem("Xim_vx_y", Xim_vx_y);
        	status = m_tuple1->addItem("Xim_vx_z", Xim_vx_z);

        	status = m_tuple1->addItem("Xip_vx_x", Xip_vx_x);
        	status = m_tuple1->addItem("Xip_vx_y", Xip_vx_y);
        	status = m_tuple1->addItem("Xip_vx_z", Xip_vx_z);

        	status = m_tuple1->addItem("Xim_pull_x", Xim_pull_x);
        	status = m_tuple1->addItem("Xim_pull_y", Xim_pull_y);
        	status = m_tuple1->addItem("Xim_pull_z", Xim_pull_z);

        	status = m_tuple1->addItem("Xip_pull_x", Xip_pull_x);
        	status = m_tuple1->addItem("Xip_pull_y", Xip_pull_y);
        	status = m_tuple1->addItem("Xip_pull_z", Xip_pull_z);

        	status = m_tuple1->addItem("chi2_4C", chi2_4C); 

        	status = m_tuple1->addItem("Rxy_pr", Rxy_pr);
        	status = m_tuple1->addItem("Rxy_prbar", Rxy_prbar);
        	status = m_tuple1->addItem("Rxy_pimxim", Rxy_pimxim);
        	status = m_tuple1->addItem("Rxy_pipxip", Rxy_pipxip);
        	status = m_tuple1->addItem("Rxy_pimlam", Rxy_pimlam);
        	status = m_tuple1->addItem("Rxy_piplambar", Rxy_piplambar);

        	status = m_tuple1->addItem("pr_drho", pr_drho);
        	status = m_tuple1->addItem("pr_phi0", pr_phi0);
        	status = m_tuple1->addItem("pr_kappa", pr_kappa);
        	status = m_tuple1->addItem("pr_lambda", pr_lambda);

        	status = m_tuple1->addItem("pimlam_drho", pimlam_drho);
        	status = m_tuple1->addItem("pimlam_phi0", pimlam_phi0);
        	status = m_tuple1->addItem("pimlam_kappa", pimlam_kappa);
        	status = m_tuple1->addItem("pimlam_lambda", pimlam_lambda);

        	status = m_tuple1->addItem("pimxim_drho", pimxim_drho);
        	status = m_tuple1->addItem("pimxim_phi0", pimxim_phi0);
        	status = m_tuple1->addItem("pimxim_kappa", pimxim_kappa);
        	status = m_tuple1->addItem("pimxim_lambda", pimxim_lambda);

        	status = m_tuple1->addItem("prbar_drho", prbar_drho);
        	status = m_tuple1->addItem("prbar_phi0", prbar_phi0);
        	status = m_tuple1->addItem("prbar_kappa", prbar_kappa);
        	status = m_tuple1->addItem("prbar_lambda", prbar_lambda);

        	status = m_tuple1->addItem("piplambar_drho", piplambar_drho);
        	status = m_tuple1->addItem("piplambar_phi0", piplambar_phi0);
        	status = m_tuple1->addItem("piplambar_kappa", piplambar_kappa);
        	status = m_tuple1->addItem("piplambar_lambda", piplambar_lambda);

        	status = m_tuple1->addItem("pipxip_drho", pipxip_drho);
        	status = m_tuple1->addItem("pipxip_phi0", pipxip_phi0);
        	status = m_tuple1->addItem("pipxip_kappa", pipxip_kappa);
        	status = m_tuple1->addItem("pipxip_lambda", pipxip_lambda);  	

        	status = m_tuple1->addItem("prvfit_drho", prvfit_drho);
        	status = m_tuple1->addItem("prvfit_phi0", prvfit_phi0);
        	status = m_tuple1->addItem("prvfit_kappa", prvfit_kappa);
        	status = m_tuple1->addItem("prvfit_lambda", prvfit_lambda);

        	status = m_tuple1->addItem("pimlamvfit_drho", pimlamvfit_drho);
        	status = m_tuple1->addItem("pimlamvfit_phi0", pimlamvfit_phi0);
        	status = m_tuple1->addItem("pimlamvfit_kappa", pimlamvfit_kappa);
        	status = m_tuple1->addItem("pimlamvfit_lambda", pimlamvfit_lambda);     

        	status = m_tuple1->addItem("pimximvfit_drho", pimximvfit_drho);
        	status = m_tuple1->addItem("pimximvfit_phi0", pimximvfit_phi0);
        	status = m_tuple1->addItem("pimximvfit_kappa", pimximvfit_kappa);
        	status = m_tuple1->addItem("pimximvfit_lambda", pimximvfit_lambda);    

        	status = m_tuple1->addItem("prbarvfit_drho", prbarvfit_drho);
        	status = m_tuple1->addItem("prbarvfit_phi0", prbarvfit_phi0);
        	status = m_tuple1->addItem("prbarvfit_kappa", prbarvfit_kappa);
        	status = m_tuple1->addItem("prbarvfit_lambda", prbarvfit_lambda);

        	status = m_tuple1->addItem("piplambarvfit_drho", piplambarvfit_drho);
        	status = m_tuple1->addItem("piplambarvfit_phi0", piplambarvfit_phi0);
        	status = m_tuple1->addItem("piplambarvfit_kappa", piplambarvfit_kappa);
        	status = m_tuple1->addItem("piplambarvfit_lambda", piplambarvfit_lambda);     

        	status = m_tuple1->addItem("pipxipvfit_drho", pipxipvfit_drho);
        	status = m_tuple1->addItem("pipxipvfit_phi0", pipxipvfit_phi0);
        	status = m_tuple1->addItem("pipxipvfit_kappa", pipxipvfit_kappa);
        	status = m_tuple1->addItem("pipxipvfit_lambda", pipxipvfit_lambda);  

			status = m_tuple1->addItem ("true_Xim_e", true_Xim_e );
			status = m_tuple1->addItem ("true_Xim_px",true_Xim_px );
			status = m_tuple1->addItem ("true_Xim_py",true_Xim_py );
			status = m_tuple1->addItem ("true_Xim_pz",true_Xim_pz );
			status = m_tuple1->addItem ("true_Xip_e", true_Xip_e );
			status = m_tuple1->addItem ("true_Xip_px",true_Xip_px );
			status = m_tuple1->addItem ("true_Xip_py",true_Xip_py );
			status = m_tuple1->addItem ("true_Xip_pz",true_Xip_pz );
			status = m_tuple1->addItem ("true_Lam_e", true_Lam_e );
			status = m_tuple1->addItem ("true_Lam_px",true_Lam_px );
			status = m_tuple1->addItem ("true_Lam_py",true_Lam_py );
			status = m_tuple1->addItem ("true_Lam_pz",true_Lam_pz );
			status = m_tuple1->addItem ("true_Lambar_e",true_Lambar_e );
			status = m_tuple1->addItem ("true_Lambar_px",true_Lambar_px );
			status = m_tuple1->addItem ("true_Lambar_py",true_Lambar_py );
			status = m_tuple1->addItem ("true_Lambar_pz",true_Lambar_pz );
			status = m_tuple1->addItem ("true_pr_e",true_pr_e );
			status = m_tuple1->addItem ("true_pr_px",true_pr_px );
			status = m_tuple1->addItem ("true_pr_py",true_pr_py );
			status = m_tuple1->addItem ("true_pr_pz",true_pr_pz );	
			status = m_tuple1->addItem ("true_prbar_e",true_prbar_e );
			status = m_tuple1->addItem ("true_prbar_px",true_prbar_px );
			status = m_tuple1->addItem ("true_prbar_py",true_prbar_py );
			status = m_tuple1->addItem ("true_prbar_pz",true_prbar_pz );
			status = m_tuple1->addItem ("true_pimLam_e",true_pimLam_e );
			status = m_tuple1->addItem ("true_pimLam_px",true_pimLam_px );
			status = m_tuple1->addItem ("true_pimLam_py",true_pimLam_py );
			status = m_tuple1->addItem ("true_pimLam_pz",true_pimLam_pz );
			status = m_tuple1->addItem ("true_pipLambar_e",true_pipLambar_e );
			status = m_tuple1->addItem ("true_pipLambar_px",true_pipLambar_px );
			status = m_tuple1->addItem ("true_pipLambar_py",true_pipLambar_py );
			status = m_tuple1->addItem ("true_pipLambar_pz",true_pipLambar_pz );
			status = m_tuple1->addItem ("true_pimXim_e",true_pimXim_e );
			status = m_tuple1->addItem ("true_pimXim_px",true_pimXim_px );
			status = m_tuple1->addItem ("true_pimXim_py",true_pimXim_py );
			status = m_tuple1->addItem ("true_pimXim_pz",true_pimXim_pz );
			status = m_tuple1->addItem ("true_pipXip_e",true_pipXip_e );
			status = m_tuple1->addItem ("true_pipXip_px",true_pipXip_px );
			status = m_tuple1->addItem ("true_pipXip_py",true_pipXip_py );
			status = m_tuple1->addItem ("true_pipXip_pz",true_pipXip_pz );		
		}
		else{ 
			log << MSG::ERROR << "Cannot book N-tuple:" << long(m_tuple1) <<endmsg;
//			return StatusCode::FAILURE;
		}
	}
	log << MSG::INFO << "Successfully return from initialize()" << endmsg;
//	cout << "I made it this far" << endl;
//	status = StatusCode::SUCCESS;
	return StatusCode::SUCCESS; // status;
}

bool XiXibar::WeighDataTrue(){
	bool keep = 0;

	// weight max has to be calculated from a large data sample. This is used to compare a randomly drawn number between (0-1)*weightmax
	// to the calculated weight from the input angles.

	Double_t weight_max = 4.13707;


	for(int i = 0; i < 7; i++) par[i] = 0.;

	// par[0] =  0.57;	// alpha J/Psi
	// par[1] =  1.20;	// relative phase Dphi J/Psi
	// par[2] = -0.37; // alpha Xi--> Lambda pi-
	// par[3] =  0.020;	// Xi phi
	// par[4] =  0.75;	// Lambda -->ppi-
	// par[5] =  0.37; // alpha Xi--> Lambda pi-
	// par[6] = -0.75; // Lambda -->ppi-

	par[0] =  0.585; // alpha J/Psi
	par[1] =  1.25; // relative phase Dphi J/Psi
	par[2] = -0.375; // alpha Xi--> Lambda pi-
	par[3] =  0.020;  // Xi phi
	par[4] =  0.75; // Lambda -->ppi-
	par[5] =  0.375; // alpha Xi--> Lambda pi-
	par[6] = -0.75; // Lambda -->ppi-

	AngDis angdis(par[0], par[1], par[2], par[3], par[4], par[5], par[6]);

	TLorentzVector xi, xib, lam, lamb, pr, prb;

	xi.SetPxPyPzE(mc_Xim_px, mc_Xim_py, mc_Xim_pz, mc_Xim_e);
	xib.SetPxPyPzE(mc_Xip_px, mc_Xip_py, mc_Xip_pz, mc_Xip_e);

	// std::cout << xi.Px() << " " << xi.Py() << " " << xi.Pz() << " " << xi.E() << " " << (xi+xib).M() << endl;


	lam.SetPxPyPzE(mc_Lam_px, mc_Lam_py, mc_Lam_pz, mc_Lam_e);
	lamb.SetPxPyPzE(mc_Lambar_px, mc_Lambar_py, mc_Lambar_pz, mc_Lambar_e);

	pr.SetPxPyPzE(mc_pr_px, mc_pr_py, mc_pr_pz, mc_pr_e);
	prb.SetPxPyPzE(mc_prbar_px, mc_prbar_py, mc_prbar_pz, mc_prbar_e);

	std::vector<TLorentzVector> fourvec_true;

	fourvec_true.resize(0);
	fourvec_true.push_back(xi); 	fourvec_true.push_back(xib); 
	fourvec_true.push_back(lam); 	fourvec_true.push_back(lamb);
	fourvec_true.push_back(pr); 	fourvec_true.push_back(prb); 

	std::vector<double> angles = HelicityAnglesXi(fourvec_true);

	Double_t w = angdis(angles.at(0), angles.at(1), angles.at(2), angles.at(3), angles.at(4), angles.at(5), angles.at(6), angles.at(7), angles.at(8));
	hweight->Fill(w);

	//std::cout << "w is " << w <<endl;

	//double test = gRandom->Rndm()*weight_max;

	//std::cout << "test is " << test <<endl;

	if((gRandom->Rndm()*weight_max) < w){
		keep = 1;
	}

	return keep;
}

double XiXibar::Angles(CLHEP::HepLorentzVector& cm, int i){

	double rxy = sqrt( pow(cm.px(),2.)+pow(cm.py(),2.) );
	if( cm.rho()<=1e-10 ){
		if(i==0) return 0.;
		else if(i==1) return 0.;
		else if(i==2) return 0.;
		else if(i==3) return 0.;
		else {cout << "Angles(i): i<=2" << endl; abort();}
	}
	else if( rxy<=1e-10 ){
		if(i==0) return cm.rho();
		else if(i==1) return 0.;
		else if(i==2) return 0.;
		else {cout << "Angles(i): i<=2" << endl; abort();}
	}
	else{
		double theta = acos(cm.pz()/cm.rho());
		double csphi = cm.px()/cm.rho()/sin(theta);
		if( fabs(csphi)>1.0 ) csphi = csphi/fabs(csphi);
		double phi = acos(csphi);
		if(cm.py()<0.0) phi = 2*3.1415926-phi;

		if(i==0) return cm.rho();
		else if(i==1) return theta;
		else if(i==2) return phi;
		else {cout << "Angles(i): i<=2" << endl; abort();}
	}
}


CLHEP::HepLorentzVector XiXibar::Helrotate(CLHEP::HepLorentzVector& p1, double phi, double theta){

	HepLorentzVector Rp;
	double cosphi = cos(phi);
	double sinphi = sin(phi);
	double costheta = cos(theta);
	double sintheta = sin(theta);
	double t = p1.e(), x = p1.px(), y = p1.py(), z = p1.pz();
	
    double xp = x*cosphi*costheta+y*sinphi*costheta-z*sintheta; 
    double yp = -x*sinphi+y*cosphi;
    double zp = x*cosphi*sintheta+y*sinphi*sintheta+z*costheta;

	Rp.setE(t);
	Rp.setPx(xp);
	Rp.setPy(yp);
	Rp.setPz(zp);

	return Rp;
}

std::vector<double> XiXibar::HelicityAnglesXi(std::vector<TLorentzVector> fourvec){
    std::vector<double> HelicityAngles;
    HelicityAngles.resize(0);

    double xith, xiphi, xibth, xibphi, lamth, lamphi, lambth, lambphi;

    TLorentzVector xi, xibar, lam, lambar, pr, prbar;
    xi = fourvec.at(0); xibar = fourvec.at(1); lam = fourvec.at(2); lambar = fourvec.at(3); pr = fourvec.at(4); prbar = fourvec.at(5); 

    TVector3 bCMS = (xi+xibar).Vect(); bCMS *= 1/((xi+xibar).E());

    // Boost Xis in CM frame
    xi.Boost(-bCMS);        lam.Boost(-bCMS);       pr.Boost(-bCMS);
    xibar.Boost(-bCMS);     lambar.Boost(-bCMS);    prbar.Boost(-bCMS);

    HelicityAngles.push_back(xi.Theta());
    xith = xi.Theta();  xiphi = xi.Phi();

     // Velocity of Xim and Xip
    TVector3 bXi(xi.Vect());        bXi*=1/(xi.E());

    lam.Boost(-bXi);    lam.RotateZ(-xiphi);     lam.RotateY(-xith); 
    pr.Boost(-bXi);     pr.RotateZ(-xiphi);      pr.RotateY(-xith); 

    HelicityAngles.push_back(lam.Theta());  HelicityAngles.push_back(lam.Phi());
    lamth = lam.Theta();    lamphi = lam.Phi();

    TVector3 bLam(lam.Vect());  bLam*=1/(lam.E());

    pr.Boost(-bLam);
    pr.RotateZ(-lamphi); pr.RotateY(-lamth); 

    HelicityAngles.push_back(pr.Theta());  HelicityAngles.push_back(pr.Phi());
              
    // Doing the same for Xibar
    xibth = xibar.Theta();  xibphi = xibar.Phi();
    TVector3 bXibar(xibar.Vect());  bXibar*=1/(xibar.E());

    lambar.Boost(-bXibar); 
    lambar.RotateZ(-xibphi);    lambar.RotateY(-xibth);
    prbar.Boost(-bXibar);
    prbar.RotateZ(-xibphi);     prbar.RotateY(-xibth);

    HelicityAngles.push_back(lambar.Theta());  HelicityAngles.push_back(lambar.Phi());
    lambth = lambar.Theta();    lambphi = lambar.Phi(); 
            
    TVector3 bLambar(lambar.Vect());  bLambar*=1/(lambar.E());
    prbar.Boost(-bLambar);
    prbar.RotateZ(-lambphi); prbar.RotateY(-lambth);

    HelicityAngles.push_back(prbar.Theta());  HelicityAngles.push_back(prbar.Phi());
           
    return HelicityAngles;

}
