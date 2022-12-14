// This is v0.3, where nCharged is changed to 4, required IP-region is made larger
// photon cut unchanged. 
// And editing PID: it identifies pi+, pi-, p, pbar
// Seems to work. Time to check some plots. Cleaning up PID-part.
// Also added primary and secondary vetex fit!
//
// Benjamin Verbeek, Hefei, 2022-11-29

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"



#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "SigmaSigmabarAlg/SigmaSigmabar.h"

//#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "VertexFit/SecondVertexFit.h"	// added
#include "ParticleID/ParticleID.h"

#include <vector>
//const double twopi = 6.2831853;
//const double pi = 3.1415927;
const double mpi = 0.139570;  // pion mass [GeV]
const double mp  = 0.938272;  // Proton mass, from PDG 2022 [GeV]
// Lambda mass

// TODO: Maybe lambda too? Sigma?
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
//const double velc = 29.9792458;  tof_path unit in cm.
const double velc = 299.792458;   // tof path unit in mm   <--- WHAT IS THIS? /BV
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<WTrackParameter> VWTP;

int Ncut0,Ncut1,Ncut2,Ncut3,Ncut4,Ncut5,Ncut6,Ncut30, Ncut31;  // Initializing cut-counters. OK. /BV
int Npi, Nproton; // Initializing counters for pions and protons

/////////////////////////////////////////////////////////////////////////////

SigmaSigmabar::SigmaSigmabar(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator) {
  
  //Declare the properties  
  declareProperty("Vr0cut", m_vr0cut=1.0);  // Cylinder cut condition
  declareProperty("Vz0cut", m_vz0cut=20.0); // From memo /BV 221129
  // TODO: polar angle cut instead

  declareProperty("EnergyThreshold", m_energyThreshold=0.04); // Fake gamma cut
  declareProperty("GammaPhiCut", m_gammaPhiCut=20.0);
  declareProperty("GammaThetaCut", m_gammaThetaCut=20.0);
  declareProperty("GammaAngleCut", m_gammaAngleCut=20.0);

  // Flag: use or not (1/0)
  declareProperty("Test4C", m_test4C = 1);
  declareProperty("Test5C", m_test5C = 1);
  declareProperty("CheckDedx", m_checkDedx = 1);
  declareProperty("CheckTof",  m_checkTof = 1);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode SigmaSigmabar::initialize(){
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "in initialize()" << endmsg;
  
  //
  //  BOOK ALL OUTPUT NTUPLES
  //
  StatusCode status;
  NTuplePtr nt1(ntupleSvc(), "FILE1/vxyz");
  if ( nt1 ) m_tuple1 = nt1;
  else {
    m_tuple1 = ntupleSvc()->book ("FILE1/vxyz", CLID_ColumnWiseTuple, "ks N-Tuple example");
    if ( m_tuple1 )    {
      status = m_tuple1->addItem ("vx0",   m_vx0); // CHARGED TRACK VERTEX
      status = m_tuple1->addItem ("vy0",   m_vy0);
      status = m_tuple1->addItem ("vz0",   m_vz0);
      status = m_tuple1->addItem ("vr0",   m_vr0);
      status = m_tuple1->addItem ("rvxy0",  m_rvxy0); // what is this? Angle stuff? /BV
      status = m_tuple1->addItem ("rvz0",   m_rvz0);
      status = m_tuple1->addItem ("rvphi0", m_rvphi0);
    }
    else    { 
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
      return StatusCode::FAILURE;
    }
  }

  NTuplePtr nt2(ntupleSvc(), "FILE1/photon");
  if ( nt2 ) m_tuple2 = nt2;
  else {
    m_tuple2 = ntupleSvc()->book ("FILE1/photon", CLID_ColumnWiseTuple, "ks N-Tuple example");
    if ( m_tuple2 )    {
      status = m_tuple2->addItem ("dthe",   m_dthe);  // photon stuff
      status = m_tuple2->addItem ("dphi",   m_dphi);
      status = m_tuple2->addItem ("dang",   m_dang);
      status = m_tuple2->addItem ("eraw",   m_eraw);
    }
    else    { 
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple2) << endmsg;
      return StatusCode::FAILURE;
    }
  }


  NTuplePtr nt3(ntupleSvc(), "FILE1/etot");
  if ( nt3 ) m_tuple3 = nt3;
  else {
    m_tuple3 = ntupleSvc()->book ("FILE1/etot", CLID_ColumnWiseTuple, "ks N-Tuple example");
    if ( m_tuple3 )    {
      status = m_tuple3->addItem ("m2gg",   m_m2gg);  // 2 gamma mass
      status = m_tuple3->addItem ("etot",   m_etot);
    }
    else    { 
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple3) << endmsg;
      return StatusCode::FAILURE;
    }
  }
  if(m_test4C==1) {
    NTuplePtr nt4(ntupleSvc(), "FILE1/fit4c");
    if ( nt4 ) m_tuple4 = nt4;
    else {
      m_tuple4 = ntupleSvc()->book ("FILE1/fit4c", CLID_ColumnWiseTuple, "ks N-Tuple example");
      if ( m_tuple4 )    {
	status = m_tuple4->addItem ("chi2",   m_chi1);
	status = m_tuple4->addItem ("mpi0",   m_mpi0);  // change for my particles for kinematic fit? /BV
  status = m_tuple4->addItem ("mpip",   m_mpip);  // TODO
  status = m_tuple4->addItem ("mpim",   m_mpim);  // Added by me before
  status = m_tuple4->addItem ("mrho0",   m_mrho0);
      }
      else    { 
	log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple4) << endmsg;
	return StatusCode::FAILURE;
      }
    }
  } // test 4C


  if(m_test5C==1) {
    NTuplePtr nt5(ntupleSvc(), "FILE1/fit5c");
    if ( nt5 ) m_tuple5 = nt5;
    else {
      m_tuple5 = ntupleSvc()->book ("FILE1/fit5c", CLID_ColumnWiseTuple, "ks N-Tuple example");
      if ( m_tuple5 )    {
	status = m_tuple5->addItem ("chi2",   m_chi2);
	status = m_tuple5->addItem ("mrh0",   m_mrh0);  // Change to my particles? /BV
	status = m_tuple5->addItem ("mrhp",   m_mrhp);  // TODO
	status = m_tuple5->addItem ("mrhm",   m_mrhm);
      }
      else    { 
	log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple5) << endmsg;
	return StatusCode::FAILURE;
      }
    }
 
    NTuplePtr nt6(ntupleSvc(), "FILE1/geff");
    if ( nt6 ) m_tuple6 = nt6;
    else {
      m_tuple6 = ntupleSvc()->book ("FILE1/geff", CLID_ColumnWiseTuple, "ks N-Tuple example");
      if ( m_tuple6 )    {
	status = m_tuple6->addItem ("fcos",   m_fcos);
	status = m_tuple6->addItem ("elow",   m_elow);
      }
      else    { 
	log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple6) << endmsg;
	return StatusCode::FAILURE;
      }
    }
  } // test 5c

  if(m_checkDedx == 1) {
    NTuplePtr nt7(ntupleSvc(), "FILE1/dedx");
    if ( nt7 ) m_tuple7 = nt7;
    else {
      m_tuple7 = ntupleSvc()->book ("FILE1/dedx", CLID_ColumnWiseTuple, "ks N-Tuple example");
      if ( m_tuple7 )    {
	status = m_tuple7->addItem ("ptrk",   m_ptrk);
	status = m_tuple7->addItem ("chie",   m_chie);
	status = m_tuple7->addItem ("chimu",   m_chimu);  // Other particles? /BV
	status = m_tuple7->addItem ("chipi",   m_chipi);  // TODO
	status = m_tuple7->addItem ("chik",   m_chik);
	status = m_tuple7->addItem ("chip",   m_chip);
	status = m_tuple7->addItem ("probPH",   m_probPH);
	status = m_tuple7->addItem ("normPH",   m_normPH);
	status = m_tuple7->addItem ("ghit",   m_ghit);
	status = m_tuple7->addItem ("thit",   m_thit);
      }
      else    { 
	log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple7) << endmsg;
	return StatusCode::FAILURE;
      }
    }
  } // check dE/dx

  if(m_checkTof == 1) {
    NTuplePtr nt8(ntupleSvc(), "FILE1/tofe"); // ENDCAP TOF
    if ( nt8 ) m_tuple8 = nt8;
    else {
      m_tuple8 = ntupleSvc()->book ("FILE1/tofe",CLID_ColumnWiseTuple, "ks N-Tuple example");
      if ( m_tuple8 )    {
	status = m_tuple8->addItem ("ptrk",   m_ptot_etof);
	status = m_tuple8->addItem ("cntr",   m_cntr_etof);
	status = m_tuple8->addItem ("ph",  m_ph_etof);      // Other particles? /BV
	status = m_tuple8->addItem ("rhit", m_rhit_etof);   // TODO
	status = m_tuple8->addItem ("qual", m_qual_etof);
	status = m_tuple8->addItem ("te",   m_te_etof);
	status = m_tuple8->addItem ("tmu",   m_tmu_etof);
	status = m_tuple8->addItem ("tpi",   m_tpi_etof);
	status = m_tuple8->addItem ("tk",   m_tk_etof);
	status = m_tuple8->addItem ("tp",   m_tp_etof);
      }
      else    { 
	log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple8) << endmsg;
	return StatusCode::FAILURE;
      }
    }
  } // check Tof:endcap



  if(m_checkTof == 1) {
    NTuplePtr nt9(ntupleSvc(), "FILE1/tof1"); // BARREL INNER TOF
    if ( nt9 ) m_tuple9 = nt9;
    else {
      m_tuple9 = ntupleSvc()->book ("FILE1/tof1", CLID_ColumnWiseTuple, "ks N-Tuple example");
      if ( m_tuple9 )    {
	status = m_tuple9->addItem ("ptrk",   m_ptot_btof1);
	status = m_tuple9->addItem ("cntr",   m_cntr_btof1);
	status = m_tuple9->addItem ("ph",  m_ph_btof1);
	status = m_tuple9->addItem ("zhit", m_zhit_btof1);
	status = m_tuple9->addItem ("qual", m_qual_btof1);
	status = m_tuple9->addItem ("te",   m_te_btof1);
	status = m_tuple9->addItem ("tmu",   m_tmu_btof1);
	status = m_tuple9->addItem ("tpi",   m_tpi_btof1);
	status = m_tuple9->addItem ("tk",   m_tk_btof1);
	status = m_tuple9->addItem ("tp",   m_tp_btof1);
      }
      else    { 
	log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple9) << endmsg;
	return StatusCode::FAILURE;
      }
    }
  } // check Tof:barrel inner Tof 


  if(m_checkTof == 1) {
    NTuplePtr nt10(ntupleSvc(), "FILE1/tof2"); // BARREL OUTER TOF
    if ( nt10 ) m_tuple10 = nt10;
    else {
      m_tuple10 = ntupleSvc()->book ("FILE1/tof2", CLID_ColumnWiseTuple, "ks N-Tuple example");
      if ( m_tuple10 )    {
	status = m_tuple10->addItem ("ptrk",   m_ptot_btof2);
	status = m_tuple10->addItem ("cntr",   m_cntr_btof2);
	status = m_tuple10->addItem ("ph",  m_ph_btof2);
	status = m_tuple10->addItem ("zhit", m_zhit_btof2);
	status = m_tuple10->addItem ("qual", m_qual_btof2);
	status = m_tuple10->addItem ("te",   m_te_btof2);
	status = m_tuple10->addItem ("tmu",   m_tmu_btof2);
	status = m_tuple10->addItem ("tpi",   m_tpi_btof2);
	status = m_tuple10->addItem ("tk",   m_tk_btof2);
	status = m_tuple10->addItem ("tp",   m_tp_btof2);
      }
      else    { 
	log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple10) << endmsg;
	return StatusCode::FAILURE;
      }
    }
  } // check Tof:barrel outter Tof


  NTuplePtr nt11(ntupleSvc(), "FILE1/pid"); // PARTICLE ID INFO
  if ( nt11 ) m_tuple11 = nt11;
  else {
    m_tuple11 = ntupleSvc()->book ("FILE1/pid", CLID_ColumnWiseTuple, "ks N-Tuple example");
    if ( m_tuple11 )    {
      status = m_tuple11->addItem ("ptrk",   m_ptrk_pid);
      status = m_tuple11->addItem ("cost",   m_cost_pid);
      status = m_tuple11->addItem ("dedx",   m_dedx_pid);
      status = m_tuple11->addItem ("tof1",   m_tof1_pid);
      status = m_tuple11->addItem ("tof2",   m_tof2_pid);
      status = m_tuple11->addItem ("prob",   m_prob_pid);
      status = m_tuple11->addItem ("pproton",m_pproton_pid);
      status = m_tuple11->addItem ("ppion",  m_ppion_pid);
    }
    else    { 
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple11) << endmsg;
      return StatusCode::FAILURE;
    }
  }
  hist = new TH1D("h1","h1",100,0,10);


  //
  //--------end of book--------
  //

  log << MSG::INFO << "successfully return from initialize()" <<endmsg;
  return StatusCode::SUCCESS;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

// * * * * * * * * * * * *  END INITIALIZE * * * * * * * * * * * * * * * * * 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

// * * * * * * * * * * * * * BEGIN EXECUTE * * * * * * * * * * * * * * * * * 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

StatusCode SigmaSigmabar::execute() {
  
  // Messages
  std::cout << "execute()" << std::endl;

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  int runNo=eventHeader->runNumber();
  int event=eventHeader->eventNumber();
  log << MSG::DEBUG <<"run, evtnum = "
      << runNo << " , "
      << event <<endreq;
  cout<<"event "<<event<<endl;
  Ncut0++;  // Counts total number of events

  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  //  log << MSG::INFO << "get event tag OK" << endreq;
    log << MSG::DEBUG <<"ncharg, nneu, tottks = " 
      << evtRecEvent->totalCharged() << " , "
      << evtRecEvent->totalNeutral() << " , "
      << evtRecEvent->totalTracks() <<endreq;
  //
  //  ######### ########## ########## ########## ########## ########## #########
  // ######## START FIRST CUT: NUMBER OF CHARGED TRACKS AND TOTAL CHARGE + REGION OF ORGIN ########
  // ########## Here, we look for charge conservation (sum 0), 4 charged tracks and only tracks
  // ########## from near the interaction point
  //  ######### ########## ########## ########## ########## ########## #########
  //

  // Some initializations:
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
  //
  // check x0, y0, z0, r0
  // suggest cut: |z0|<5 && r0<1
  //
  Vint iGood, ipip, ipim, ip, ipbar; // CHNGD: add also for protons: ip, ipbar
  iGood.clear();
  ipip.clear();
  ipim.clear();
  ip.clear();
  ipbar.clear();
  Vp4 ppip, ppim, pp, ppbar; // CHNGD: add also for protons. pp, ppbar
  ppip.clear();
  ppim.clear();
  pp.clear();
  ppbar.clear();
  
  // Start cut 1: nCharge = 4, totCharge = 0
  int nCharge = 0;

  Hep3Vector xorigin(0,0,0);
  HepPoint3D vx(0., 0., 0.);
	HepSymMatrix Evx(3, 0);
  //if (m_reader.isRunNumberValid(runNo)) {
   IVertexDbSvc*  vtxsvc;
  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if(vtxsvc->isVertexValid()){
    double* dbv = vtxsvc->PrimaryVertex();  // Get primary vertex.
    double*  vv = vtxsvc->SigmaPrimaryVertex();  
  //    HepVector dbv = m_reader.PrimaryVertex(runNo);
  //    HepVector vv = m_reader.SigmaPrimaryVertex(runNo);
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

  cout << "Hello SigmaSigmabar world! v0.3.2. VX FIT" << endl; // edited!
  for(int i = 0; i < evtRecEvent->totalCharged(); i++){ // loop over charged tracks
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;  // get track
    if(!(*itTrk)->isMdcTrackValid()) continue;  // check if MDC track is valid
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack(); // get MDC track
    double pch=mdcTrk->p();
    double x0=mdcTrk->x();
    double y0=mdcTrk->y();
    double z0=mdcTrk->z();
    double phi0=mdcTrk->helix(1);
    double xv=xorigin.x();
    double yv=xorigin.y();
    double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);
    m_vx0 = x0; // Some stuff saved to output .root m_tuple1
    m_vy0 = y0;
    m_vz0 = z0;
    m_vr0 = Rxy;

    HepVector a = mdcTrk->helix();
    HepSymMatrix Ea = mdcTrk->err();
    HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
    HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
    VFHelix helixip(point0,a,Ea); 
    helixip.pivot(IP);
    HepVector vecipa = helixip.a();
    double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
    double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
    double  Rvphi0=vecipa[1];       // WHERE is cos theta here?
    m_rvxy0=Rvxy0;      // Some stuff saved to output .root m_tuple1
    m_rvz0=Rvz0;
    m_rvphi0=Rvphi0;

    m_tuple1->write();
//    if(fabs(z0) >= m_vz0cut) continue;  // Why this change? /BV 11-28
//    if(fabs(Rxy) >= m_vr0cut) continue; 

  // THROW AWAY TRACKS NOT ORIGINATING FROM INTERACTION REGION  
    //if(fabs(Rvz0) >= m_vz0cut) continue;  // CHNGD: from 10.0 to m_vz0cut
    //if(fabs(Rvxy0) >= 10.0) continue;  // TODO: Here? Angle? // NOTE: Xiaorong claims angle cut is superfluous anyway. /BV 11-30

    
    iGood.push_back(i); // OK, add to list of good charged tracks
    nCharge += mdcTrk->charge();  // add the charge (+1 or -1)
  }
  
  hist->Fill(iGood.size()); // Fill histogram with number of charged tracks

  //
  // Finish Good Charged Track Selection (we expect pi+, pi-, p, pbar)
  //
  int nGood = iGood.size();
  log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq;
  if((nGood != 4)||(nCharge!=0)){ // CHNGD: Here, we expect 4 charged tracks with total charge = 0
    return StatusCode::SUCCESS;   // Failed test? Return ok and move to next event
  }

  Ncut1++;  // Programme didn't return? Cool! Add to counter of items that came this far.
  //
  //  ######### ########## ########## ########## ########## ########## #########
  // ######### END FIRST CUT: NUMBER OF CHARGED TRACKS AND TOTAL CHARGE #########
  //  ######### ########## ########## ########## ########## ########## #########
  // ########## START SECOND CUT: GAMMA
  // ########## Here, we look for two or more gammas.
  // ######### ########## ########## ########## ########## ########## #########
  
  // Start gamma finder, nGam >= 2
  Vint iGam;
  iGam.clear();
  for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) { // Iterate over unchaarged tracks
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;  // Uncharged are appended after charged.
    if(!(*itTrk)->isEmcShowerValid()) continue; // Black-box. (Keep) If not valid, continue.
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
    // find the nearest charged track
    double dthe = 200.;
    double dphi = 200.;
    double dang = 200.; 
    for(int j = 0; j < evtRecEvent->totalCharged(); j++) { // loop over charged tracks
      EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
      if(!(*jtTrk)->isExtTrackValid()) continue;
      RecExtTrack *extTrk = (*jtTrk)->extTrack();
      if(extTrk->emcVolumeNumber() == -1) continue;
      Hep3Vector extpos = extTrk->emcPosition();
      //      double ctht = extpos.cosTheta(emcpos);
      double angd = extpos.angle(emcpos);
      double thed = extpos.theta() - emcpos.theta();
      double phid = extpos.deltaPhi(emcpos);
      thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      if(angd < dang){ // Update if closer
        dang = angd;
        dthe = thed;
        dphi = phid;
      }
    }
    if(dang>=200) continue; // CUT: Angle between charged track and gamma must be less than 200
    double eraw = emcTrk->energy(); // ^ TODO: ok limit?
    dthe = dthe * 180 / (CLHEP::pi);
    dphi = dphi * 180 / (CLHEP::pi);
    dang = dang * 180 / (CLHEP::pi);
    m_dthe = dthe;
    m_dphi = dphi;
    m_dang = dang;
    m_eraw = eraw;
    m_tuple2->write(); // Save to output .root m_tuple2
    if(eraw < m_energyThreshold) continue; // CUT: Energy of gamma must be below threshold
//    if((fabs(dthe) < m_gammaThetaCut) && (fabs(dphi)<m_gammaPhiCut) ) continue;
    if(fabs(dang) < m_gammaAngleCut) continue; // CUT: Angle between charged track and gamma must be less than m_gammaAngleCut
    //                                                ^ TODO ok limit?
    // good photon cut will be set here
    //
    iGam.push_back(i);  // Save gamma to list of good gammas
  }
  
  //
  // Finish Good Photon Selection
  //
  int nGam = iGam.size();

  log << MSG::DEBUG << "num Good Photon " << nGam  << " , " <<evtRecEvent->totalNeutral()<<endreq;
  if(nGam<2){
    return StatusCode::SUCCESS;
  }
  Ncut2++;
  //  ######### ########## ########## ########## ########## ########## #########
  // ######### END SCOND CUT: NUMBER OF OK GAMMAS >= 2 #########
  //  ######### ########## ########## ########## ########## ########## #########
  // ########## START PID PREP. #########
  // ######### ########## ########## ########## ########## ########## #########

  //
  //
  // check dedx infomation
  //
  //
  
  // Some validity checks? 
  // TODO: do I want this?
  if(m_checkDedx == 1) {
    for(int i = 0; i < nGood; i++) {
      EvtRecTrackIterator  itTrk = evtRecTrkCol->begin() + iGood[i];
      if(!(*itTrk)->isMdcTrackValid()) continue; // Black-box. (Keep) If not valid, continue.
      if(!(*itTrk)->isMdcDedxValid())continue;
      RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
      RecMdcDedx* dedxTrk = (*itTrk)->mdcDedx();
      m_ptrk = mdcTrk->p();   // This momentum can be of use! Clear split here.
      // Can make a cut here perhaps? What do I do with the ToF/dEdx then?
      // I don't know what to use when

      m_chie = dedxTrk->chiE();
      m_chimu = dedxTrk->chiMu();
      m_chipi = dedxTrk->chiPi();
      m_chik = dedxTrk->chiK();
      m_chip = dedxTrk->chiP();
      m_ghit = dedxTrk->numGoodHits();
      m_thit = dedxTrk->numTotalHits();
      m_probPH = dedxTrk->probPH();
      m_normPH = dedxTrk->normPH();
      m_tuple7->write();
    }
  }

  //
  // check TOF infomation
  //

// TODO: do I want this?
  if(m_checkTof == 1) {
    for(int i = 0; i < nGood; i++) {
      EvtRecTrackIterator  itTrk = evtRecTrkCol->begin() + iGood[i]; // Check valid tracks
      if(!(*itTrk)->isMdcTrackValid()) continue;
      if(!(*itTrk)->isTofTrackValid()) continue;

      RecMdcTrack * mdcTrk = (*itTrk)->mdcTrack();
      SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();

      double ptrk = mdcTrk->p();

      SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
      for(;iter_tof != tofTrkCol.end(); iter_tof++ ) { 
        TofHitStatus *status = new TofHitStatus; 
        status->setStatus((*iter_tof)->status());
        if(!(status->is_barrel())){//endcap
          if( !(status->is_counter()) ) continue; // ? 
          if( status->layer()!=0 ) continue;//layer1
          double path=(*iter_tof)->path(); // ? 
          double tof  = (*iter_tof)->tof();
          double ph   = (*iter_tof)->ph();
          double rhit = (*iter_tof)->zrhit();
          double qual = 0.0 + (*iter_tof)->quality();
          double cntr = 0.0 + (*iter_tof)->tofID();
          double texp[5];
          for(int j = 0; j < 5; j++) {
            double gb = ptrk/xmass[j];
            double beta = gb/sqrt(1+gb*gb);
            texp[j] = 10 * path /beta/velc;
          }
          m_cntr_etof  = cntr;
          m_ptot_etof = ptrk;
          m_ph_etof   = ph;
          m_rhit_etof  = rhit;
          m_qual_etof  = qual;
          m_te_etof    = tof - texp[0];
          m_tmu_etof   = tof - texp[1];
          m_tpi_etof   = tof - texp[2];
          m_tk_etof    = tof - texp[3];
          m_tp_etof    = tof - texp[4];
          m_tuple8->write();
        }
        else {//barrel
          if( !(status->is_counter()) ) continue; // ? 
          if(status->layer()==1){ //layer1
            double path=(*iter_tof)->path(); // ? 
            double tof  = (*iter_tof)->tof();
            double ph   = (*iter_tof)->ph();
            double rhit = (*iter_tof)->zrhit();
            double qual = 0.0 + (*iter_tof)->quality();
            double cntr = 0.0 + (*iter_tof)->tofID();
            double texp[5];
            for(int j = 0; j < 5; j++) {
              double gb = ptrk/xmass[j];
              double beta = gb/sqrt(1+gb*gb);
              texp[j] = 10 * path /beta/velc;
            }
 
            m_cntr_btof1  = cntr;
            m_ptot_btof1 = ptrk;
            m_ph_btof1   = ph;
            m_zhit_btof1  = rhit;
            m_qual_btof1  = qual;
            m_te_btof1    = tof - texp[0];
            m_tmu_btof1   = tof - texp[1];
            m_tpi_btof1   = tof - texp[2];
            m_tk_btof1    = tof - texp[3];
            m_tp_btof1    = tof - texp[4];
            m_tuple9->write();
          }

          if(status->layer()==2){//layer2
            double path=(*iter_tof)->path(); // ? 
            double tof  = (*iter_tof)->tof();
            double ph   = (*iter_tof)->ph();
            double rhit = (*iter_tof)->zrhit();
            double qual = 0.0 + (*iter_tof)->quality();
            double cntr = 0.0 + (*iter_tof)->tofID();
            double texp[5];
            for(int j = 0; j < 5; j++) {
              double gb = ptrk/xmass[j];
              double beta = gb/sqrt(1+gb*gb);
              texp[j] = 10 * path /beta/velc;
            }
 
            m_cntr_btof2  = cntr;
            m_ptot_btof2 = ptrk;
            m_ph_btof2   = ph;
            m_zhit_btof2  = rhit;
            m_qual_btof2  = qual;
            m_te_btof2    = tof - texp[0];
            m_tmu_btof2   = tof - texp[1];
            m_tpi_btof2   = tof - texp[2];
            m_tk_btof2    = tof - texp[3];
            m_tp_btof2    = tof - texp[4];
            m_tuple10->write();
          } 
        }

        delete status; 
      } 
    } // loop all charged track
  }  // check tof


// TODO: What is going on above? ^^^ Used below? vvv in PID

  //
  // Assign 4-momentum to each photon
  // 

  Vp4 pGam;
  pGam.clear();
  for(int i = 0; i < nGam; i++) { // loop all passed gamma
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i]; 
    RecEmcShower* emcTrk = (*itTrk)->emcShower();
    double eraw = emcTrk->energy();
    double phi = emcTrk->phi();
    double the = emcTrk->theta();
    HepLorentzVector ptrk;
    ptrk.setPx(eraw*sin(the)*cos(phi));
    ptrk.setPy(eraw*sin(the)*sin(phi));
    ptrk.setPz(eraw*cos(the));
    ptrk.setE(eraw);

//    ptrk = ptrk.boost(-0.011,0,0);// boost to cms //TODO: what is this? Maybe I want to?

    pGam.push_back(ptrk); // 4-momentum of each photon added to vector
  }
  cout<<"before pid"<<endl;

  //  ######### ########## ########## ########## ########## ########## #########
  // ######### Particle Identification (PID) and momentum assignments #########
  //  ######### ########## ########## ########## ########## ########## #########
//  ######### ########## ########## ########## ########## ########## #########
  // ######### START THIRD CUT: Exact check of decay products (PID) #########
  //  ######### ########## ########## ########## ########## ########## #########

  //
  // Assign 4-momentum to each charged track
  //
  ParticleID *pid = ParticleID::instance(); // Black-box to me but okay.
  for(int i = 0; i < nGood; i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];   // iterate good tracks
    //    if(pid) delete pid;
    pid->init();
    pid->setMethod(pid->methodProbability()); // TODO: Okay?
//    pid->setMethod(pid->methodLikelihood());  //for Likelihood Method  

    pid->setChiMinCut(4);
    pid->setRecTrack(*itTrk);
    pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()); // use PID sub-system

    // TODO: Branch up here maybe, into proton and pion?
    // pid->identify(pid->onlyPion() | pid->onlyKaon());    // seperater Pion/Kaon // TODO: what here?
    pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());    // Edited 11-30, remove Kaon? Or check it also?

    pid->calculate();
    if(!(pid->IsPidInfoValid())) continue;
    RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
    m_ptrk_pid = mdcTrk->p();   // // This shows a nice split into pions and protons. Cut here?
    // Possibly: discern proton/pion by momentum here. Skip the rest

    m_cost_pid = cos(mdcTrk->theta());
    m_dedx_pid = pid->chiDedx(2);
    m_tof1_pid = pid->chiTof1(2);
    m_tof2_pid = pid->chiTof2(2);
    // get max of pid probability for proton/pion... replace this by: (??)
    bool isPion = (m_ptrk_pid < 0.45) ? true : false; // Value taken from plot/memo // Seems to work fine. One switches.
    /*bool isPion = pid->probPion() > pid->probProton();  // added 11-30
    if (isPion) {  // added
      m_prob_pid = pid->probPion();
    } else {
      m_prob_pid = pid->probProton();
    }*/

    // If change, keep prob-checks here? 
    //if((pid->probPion() < 0.001) && (pid->probProton() < 0.001)) continue; // Edited. If too unlikely, throw it away.
    isPion ? Npi++ : Nproton++; // Count pions and protons
    isPion ? m_ppion_pid = (mdcTrk->p()) : m_pproton_pid = (mdcTrk->p()); // Assign momentum to pion/proton, can plot later
    
    m_tuple11->write();  


    RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
    RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);//PID can set to electron, muon, pion, kaon and proton;The default setting is pion

    // moved out from below if-else
    HepLorentzVector ptrk;
    ptrk.setPx(mdcKalTrk->px());
    ptrk.setPy(mdcKalTrk->py());
    ptrk.setPz(mdcKalTrk->pz());
    double p3 = ptrk.mag();
    // Add this momentum to output... curious is two separate peaks:


    // Characterize!
    if((mdcKalTrk->charge() > 0) && isPion ) { // if positive & pion, its pi+
      ipip.push_back(iGood[i]); // dep
      ptrk.setE(sqrt(p3*p3+mpi*mpi)); // This is why proton mass is needed. TODO: define in intialize
      ppip.push_back(ptrk); // 4-momentum of each pion added to vector

    } else if ((mdcKalTrk->charge() < 0) && isPion) {  // if not positive and is pion, its pi-
      ipim.push_back(iGood[i]);
      ptrk.setE(sqrt(p3*p3+mpi*mpi));
      ppim.push_back(ptrk); // 4-momentum of each pion added to vector

      // Added: now check for proton in the same way:
    } else if ((mdcKalTrk->charge() > 0) && !isPion) { // if positive & proton, its p+
      ip.push_back(iGood[i]);
      ptrk.setE(sqrt(p3*p3+mp*mp));
      pp.push_back(ptrk); // 4-momentum of each proton added to vector

    // and for pbar
    } else if ((mdcKalTrk->charge() < 0) && !isPion) {  // if not positive and is proton, its pbar-
      ipbar.push_back(iGood[i]);
      ptrk.setE(sqrt(p3*p3+mp*mp)); 
      ppbar.push_back(ptrk); // 4-momentum of each proton added to vector
    }
//      ptrk = ptrk.boost(-0.011,0,0);//boost to cms
  } // characterized tracks and momenta now stored to respective lists.


  int npip = ipip.size();
  int npim = ipim.size(); // Add for proton/-bar. TODO.
  int np = ip.size();
  int npbar = ipbar.size();
  cout << "Here come npip, npim, np, npbar: " << npip << npim << np << npbar << endl; // edited!

  if(npip*npim*np*npbar != 1) return SUCCESS;  // Here we check if there is only one pion and one anti-pion. TODO: Add proton/-bar.
          // Benefit of writing product is unclear to me. Harder to read and more computationally 
          // expensive, is my guess. TODO.  We could probably just use momentum split.

  Ncut3++;

//  ######### ########## ########## ########## ########## ########## #########
  // ######### END THIRD CUT: Exact check of decay products #########
  //  ######### ########## ########## ########## ########## ########## #########
  // ########## START KINEMATIC FIT. #########
  //######### ########## ########## ########## ########## ########## #########
  // TODO: What now? Secondary vx fit to get pion + proton to virtual lambda?


  // Construct virtual lambda. Check if pi+ and pbar origin from same vertex
  // Check if pi- and p origin from same vertex

  //
  // Loop each gamma pair, check ppi0 and pTot 
  //

  HepLorentzVector pTot;
  for(int i = 0; i < nGam - 1; i++){
    for(int j = i+1; j < nGam; j++) {
      HepLorentzVector p2g = pGam[i] + pGam[j];
      pTot = ppip[0] + ppim[0];
      pTot += p2g;
      m_m2gg = p2g.m();
      m_etot = pTot.e();
      m_tuple3 -> write();
    }
  }
  
  Vp4 p4Lambdavtx, p4pvtx, p4pbarvtx, p4pimvtx, p4pipvtx;
  p4Lambdavtx.clear();
  p4pbarvtx.clear();
  p4pvtx.clear();
  p4pimvtx.clear();
  p4pipvtx.clear();
  //Vdouble  decayL_lambdabar, decayLerr_lambdabar, chisq_lambdabar;
  //decayLerr_lambdabar.clear();
  //decayL_lambdabar.clear();
  //chisq_lambdabar.clear();
  //VWTP wlambdabar_vertex;
  //wlambdabar_vertex.clear();
  Vdouble  decayL_Lambda, decayLerr_Lambda, chisq_Lambda;
  decayLerr_Lambda.clear();
  decayL_Lambda.clear();
  chisq_Lambda.clear();
  VWTP wLambda_vertex;
  wLambda_vertex.clear();
  HepPoint3D cPoint;

  // point to start of pions. TODO: also for protons?
  RecMdcKalTrack *pipTrk  = (*(evtRecTrkCol->begin()+ipip[0]))->mdcKalTrack();
  RecMdcKalTrack *pimTrk  = (*(evtRecTrkCol->begin()+ipim[0]))->mdcKalTrack();
  RecMdcKalTrack *pTrk    = (*(evtRecTrkCol->begin()+ip[0]))->mdcKalTrack();
  RecMdcKalTrack *pbarTrk = (*(evtRecTrkCol->begin()+ipbar[0]))->mdcKalTrack();

  WTrackParameter wvpipTrk, wvpimTrk, wvpTrk, wvpbarTrk; // TODO: also for protons? see below
  wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());
  wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());
  wvpTrk   = WTrackParameter(mp,  pTrk->getZHelixP(),  pTrk->getZErrorP());
  wvpbarTrk= WTrackParameter(mp, pbarTrk->getZHelixP(), pbarTrk->getZErrorP());

/* Default is pion, for other particles:
  wvppTrk = WTrackParameter(mp, pipTrk->getZHelixP(), pipTrk->getZErrorP());//proton
  wvmupTrk = WTrackParameter(mmu, pipTrk->getZHelixMu(), pipTrk->getZErrorMu());//muon
  wvepTrk = WTrackParameter(me, pipTrk->getZHelixE(), pipTrk->getZErrorE());//electron
  wvkpTrk = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK());//kaon
*/
  //
  //    Test vertex fit
  //    TODO: Do I want to vertex fit, how? 2nd vx-fit. Lambda reconstruction?
  //    Reconstruct two oppositely charged particle tracks, converging.
  //

  // Initialize (origin + large error)  DEFINED ABOVE
  // HepPoint3D vx(0., 0., 0.);
  HepSymMatrix Evx2(3, 0);
  double bx = 1E+6;
  double by = 1E+6;
  double bz = 1E+6;
  Evx2[0][0] = bx*bx;
  Evx2[1][1] = by*by;
  Evx2[2][2] = bz*bz;

  // Set these initial values to the vertex object
  VertexParameter vxpar;
  vxpar.setVx(vx);

  vxpar.setEvx(Evx2);  // Previously wrong error
  // fit for Lambda from pi- and p
  VertexFit* vtxfit = VertexFit::instance();
  vtxfit->init();
  vtxfit->AddTrack(0,  wvpTrk); // TODO: proton/pi- and pi+/proton-bar instead.
  vtxfit->AddTrack(1,  wvpimTrk);
  vtxfit->AddVertex(0, vxpar,0, 1);
  //if(!vtxfit->Fit(0)) return SUCCESS; // If the fit fails, I guess? Is it assuming from origin?
  // Primary VERTEX FIT
  if(!vtxfit->Fit(0)) return SUCCESS;	// else, ... (to keep in scope)
    vtxfit->Swim(0);
    vtxfit->BuildVirtualParticle(0);
    WTrackParameter wLambda = vtxfit->wVirtualTrack(0);
    VertexParameter vtxLambda = vtxfit->vpar(0);
    WTrackParameter wtrkproton = vtxfit->wtrk(0);
    p4pvtx.push_back(wtrkproton.p());
    WTrackParameter wtrkpion = vtxfit->wtrk(1);
    p4pimvtx.push_back(wtrkpion.p());
    //HepLorentzVector pLambda = wLambda.p();
    //m_mLambda = pLambda.m();
    //m_chisq_vf = vtxfit->chisq(0);

    // Secondary VERTEX FIT
    SecondVertexFit *vtxfit2 = SecondVertexFit::instance();
    vtxfit2->init();
    vtxfit2->setPrimaryVertex(vx_db); // What is this? average interaction point from database
    vtxfit2->AddTrack(0, wLambda);
    vtxfit2->setVpar(vtxLambda);

    if(vtxfit2->Fit()) {
      HepLorentzVector p4Lambda = vtxfit2->p4par();

      p4Lambdavtx.push_back(p4Lambda);
      decayL_Lambda.push_back(vtxfit2->decayLength());
      decayLerr_Lambda.push_back(vtxfit2->decayLengthError());
      chisq_Lambda.push_back(vtxfit2->chisq());
      wLambda_vertex.push_back(vtxfit2->wpar());
      cPoint = vtxfit2->crossPoint();
    } else {
			return SUCCESS;
		}
  
  Ncut30++;
  
  Vp4 p4Lambdabarvtx;
	p4Lambdabarvtx.clear();
  Vdouble  decayL_Lambdabar, decayLerr_Lambdabar, chisq_Lambdabar;
  decayLerr_Lambdabar.clear();
  decayL_Lambdabar.clear();
  chisq_Lambdabar.clear();
  VWTP wLambdabar_vertex;
  wLambdabar_vertex.clear();

  // fit for Lambdabar from pi+ and pbar
  vtxfit->init();
  vtxfit->AddTrack(0,  wvpbarTrk);
  vtxfit->AddTrack(1,  wvpipTrk);
  vtxfit->AddVertex(0, vxpar,0, 1);

  if(!vtxfit->Fit(0)) return SUCCESS;
    vtxfit->Swim(0);
    vtxfit->BuildVirtualParticle(0);
    WTrackParameter wLambdabar = vtxfit->wVirtualTrack(0);
    VertexParameter vtxLambdabar = vtxfit->vpar(0);
    WTrackParameter wtrkprotonbar = vtxfit->wtrk(0);
    p4pbarvtx.push_back(wtrkprotonbar.p());
    WTrackParameter wtrkpionbar = vtxfit->wtrk(1);
    p4pipvtx.push_back(wtrkpionbar.p());

    // Secondary vtx fit
    //SecondVertexFit *vtxfit2 = SecondVertexFit::instance(); // Still in same scope.
    vtxfit2->init();
    vtxfit2->setPrimaryVertex(vx_db); // What is this?
    vtxfit2->AddTrack(0, wLambdabar);
    vtxfit2->setVpar(vtxLambdabar);
    if(vtxfit2->Fit()) {
      HepLorentzVector p4Lambdabar = vtxfit2->p4par();

      p4Lambdabarvtx.push_back(p4Lambdabar);
      decayL_Lambdabar.push_back(vtxfit2->decayLength());
      decayLerr_Lambdabar.push_back(vtxfit2->decayLengthError());
      chisq_Lambdabar.push_back(vtxfit2->chisq());
      wLambdabar_vertex.push_back(vtxfit2->wpar()); // push_back not needed
      //cPoint = vtxfit2->crossPoint(); // only done first time. Why?
    } else {
			return SUCCESS;
		}

  Ncut31++;

  // Now I have wLambdabar and wLambda (virtual tracks)
  // Kinematic fit with photons?

  //WTrackParameter wpip = vtxfit->wtrk(0);
  //WTrackParameter wpim = vtxfit->wtrk(1);

  //KinematicFit * kmfit = KinematicFit::instance();
  KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();

  //  
  //  Apply Kinematic 4C fit
  // 
  cout<<"before 4c"<<endl;
  if(m_test4C==1) {
//    double ecms = 3.097;
    HepLorentzVector ecms(0.034,0,0,3.097); // TODO: Check

    double chisq = 9999.;
    int ig1 = -1;
    int ig2 = -1;
    for(int i = 0; i < nGam-1; i++) { // iterate over pairs of gammas
      RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
      for(int j = i+1; j < nGam; j++) {
        RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();
        kmfit->init();
        kmfit->AddTrack(0, wLambda_vertex[0]);
        kmfit->AddTrack(1, wLambdabar_vertex[0]);
        kmfit->AddTrack(2, 0.0, g1Trk);
        kmfit->AddTrack(3, 0.0, g2Trk); // also protons here then...
        kmfit->AddFourMomentum(0, ecms);
        bool oksq = kmfit->Fit();
        if(oksq) {
          double chi2 = kmfit->chisq();
          if(chi2 < chisq) {  // better fit? Update gamma choice.
            chisq = chi2;
            ig1 = iGam[i];
            ig2 = iGam[j];
          }
        }
            }
    }
    
    // Get the best values
    // Else what? Crash?
    if(chisq < 200) {   // Memo suggests less than 100...

      RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+ig1))->emcShower();
      RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+ig2))->emcShower();
      kmfit->init();
      kmfit->AddTrack(0, wLambda);
      kmfit->AddTrack(1, wLambdabar);
      kmfit->AddTrack(2, 0.0, g1Trk); // <--- H??r l??ggs gamma1-sp??ret in
      kmfit->AddTrack(3, 0.0, g2Trk); // <--- H??r l??ggs gamma2-sp??ret in
      kmfit->AddFourMomentum(0, ecms);
      bool oksq = kmfit->Fit(); // Unnecessary re-check...?
      if(oksq) {
        // TODO: Also add protons here

				// For now, store ppip for pLambda, 
				// ppim for pLambdabar...
		// NOTE: TODO: These should be renamed.
	HepLorentzVector ppi0 = kmfit->pfit(2) + kmfit->pfit(3);  // <--- H??r konstrueras lorentzvektorn f??r pi0 fr??n g1, g2
  HepLorentzVector ppip = kmfit->pfit(0); // <-- Detta la jag till nu
  HepLorentzVector ppim = kmfit->pfit(1); // <-- Added 2022-11-24 21:30 /BV
  HepLorentzVector prho0 = kmfit->pfit(0) + kmfit->pfit(1); // For rho0? 
  // Lambda here.... ^

	// NOTE: TODO: These should be renamed.
	m_mpi0 = ppi0.m(); // <--- H??r tilldelas variabeln m_mpi0 som ??r definierad invarianta masssan f??r pi0
  m_mpip = ppip.m(); // <--- Detta la jag till nu
  m_mpim = ppim.m(); // <-- Added 2022-11-24 21:30 /BV	// should now be Lambda mass if ok.
  m_mrho0 = prho0.m(); // <-- Added 2022-11-24 21:30 /BV
  // TODO: Here for lambda also...
	m_chi1 = kmfit->chisq();
	m_tuple4->write();
        Ncut4++;    // What survived the reconstruction...
      }
    }
  }


 /*	// Commented out 5C fit. 2022-12-02 17:24 /BV
  //
  //  Apply Kinematic 5C Fit
  //

  // find the best combination over all possible pi+ pi- gamma gamma pair
  // TODO: Not relevant for me? Rather, 6C fit? How to know if sigma or just lambda + noise gamma?
  if(m_test5C==1) {
//    double ecms = 3.097;
    HepLorentzVector ecms(0.034,0,0,3.097);
    double chisq = 9999.;
    int ig1 = -1;
    int ig2 = -1;
    for(int i = 0; i < nGam-1; i++) {
      RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
      for(int j = i+1; j < nGam; j++) {
	RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();
	kmfit->init();
	kmfit->AddTrack(0, wpip);
	kmfit->AddTrack(1, wpim);
	kmfit->AddTrack(2, 0.0, g1Trk);
	kmfit->AddTrack(3, 0.0, g2Trk);
	kmfit->AddResonance(0, 0.135, 2, 3);
	kmfit->AddFourMomentum(1, ecms);
	if(!kmfit->Fit(0)) continue;
	if(!kmfit->Fit(1)) continue;
	bool oksq = kmfit->Fit();
	if(oksq) {
	  double chi2 = kmfit->chisq();
	  if(chi2 < chisq) {
	    chisq = chi2;
	    ig1 = iGam[i];
	    ig2 = iGam[j];
	  }
	}
      }
    }
  

    log << MSG::INFO << " chisq = " << chisq <<endreq;

    if(chisq < 200) {
      RecEmcShower* g1Trk = (*(evtRecTrkCol->begin()+ig1))->emcShower();
      RecEmcShower* g2Trk = (*(evtRecTrkCol->begin()+ig2))->emcShower();

      kmfit->init();
      kmfit->AddTrack(0, wpip);
      kmfit->AddTrack(1, wpim);
      kmfit->AddTrack(2, 0.0, g1Trk);
      kmfit->AddTrack(3, 0.0, g2Trk);
      kmfit->AddResonance(0, 0.135, 2, 3);
      kmfit->AddFourMomentum(1, ecms);
      bool oksq = kmfit->Fit();
      if(oksq){
	HepLorentzVector ppi0 = kmfit->pfit(2) + kmfit->pfit(3);
	HepLorentzVector prho0 = kmfit->pfit(0) + kmfit->pfit(1);
	HepLorentzVector prhop = ppi0 + kmfit->pfit(0);
	HepLorentzVector prhom = ppi0 + kmfit->pfit(1);
	
	m_chi2  = kmfit->chisq();
	m_mrh0 = prho0.m();
	m_mrhp = prhop.m();
	m_mrhm = prhom.m();
	double eg1 = (kmfit->pfit(2)).e();
	double eg2 = (kmfit->pfit(3)).e();
	double fcos = abs(eg1-eg2)/ppi0.rho();
	m_tuple5->write();
        Ncut5++;
	// 
	//  Measure the photon detection efficiences via
	//          J/psi -> rho0 pi0
	//
	if(fabs(prho0.m()-0.770)<0.150) {  
	  if(fabs(fcos)<0.99) {
	    m_fcos = (eg1-eg2)/ppi0.rho();
	    m_elow =  (eg1 < eg2) ? eg1 : eg2;
	    m_tuple6->write();
            Ncut6++;
	  }
	} // rho0 cut
      }  //oksq
    } 
  }
  return StatusCode::SUCCESS;
*/ 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode SigmaSigmabar::finalize() {
  cout<<"total number:         "<<Ncut0<<endl;
  cout<<"nGood==4, nCharge==0: "<<Ncut1<<endl;
  cout<<"nGam>=2:              "<<Ncut2<<endl;
  cout<<"Pass Pid:             "<<Ncut3<<endl;
  cout<<"Pass 4C:              "<<Ncut4<<endl;
  cout<<"Pass 5C:              "<<Ncut5<<endl;
  cout<<"J/psi->rho0 pi0:      "<<Ncut6<<endl;
  cout<<"Pass lambdavx:        "<<Ncut30<<endl;
  cout<<"Pass lambdabarvx:       "<<Ncut31<<endl;
  cout<<"Npi, Nproton, mom.: " << Npi << " " << Nproton << endl;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;
  return StatusCode::SUCCESS;
}

