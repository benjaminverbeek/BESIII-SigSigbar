#ifndef Physics_Analysis_Psip_SigLampi0_H
#define Physics_Analysis_Psip_SigLampi0_H

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "TH2.h"

class Psip_SigLampi0 : public Algorithm
{
    public:
        Psip_SigLampi0(const std::string& name, ISvcLocator* pSvcLocator);
        StatusCode initialize();
        StatusCode execute();
        StatusCode finalize();
    private:
        int m_evnt;
        double m_ecms;
        double m_vr0cut;
        double m_vz0cut;
        double m_pass[10];
        double m_chisqCut;
        double m_energyThreshold;
        double m_massRangeLower;
        double m_massRangeHigher;
        double m_gammaPhiCut;
        double m_gammaThetaCut;
        double m_gammaAngCut;
        double m_Endcap_th_1;
        double m_Endcap_th_2;
        double m_energyThreshold_b;
        double m_energyThreshold_e;
        double m_Barrel_th;
        bool m_checkMC;
        int m_idNo;
        //   McDecayModeSvc* m_svc;
        std::string m_log;
        //ReadBeamParFromDb m_reader;
        NTuple::Tuple*  m_tuple;

        NTuple::Item<long>   m_runnumber;
        NTuple::Item<double>   m_runsample;
        NTuple::Item<double>   m_runstatus;

        NTuple::Item<long>   m_run;
        NTuple::Item<long>   m_rec;
        NTuple::Item<long>   m_flag1;
        NTuple::Item<long>   m_flag2;

        NTuple::Item<long>   m_idxmc;
        NTuple::Array<long>  m_trkidx;//655
        NTuple::Array<long>  m_pdgid; //655
        NTuple::Array<long>  m_motheridx; //655

        //NTuple::Item<int>   m_idxmc;//664
        //NTuple::Array<int>  m_trkidx;//664
        //NTuple::Array<int>  m_pdgid;//664
        //NTuple::Array<int>  m_motheridx;//664
        NTuple::Item<double> m_cosThXi0_mc,m_cosThXi0b_mc,m_pxi0_mc,m_pxi0b_mc;
        NTuple::Item<double> m_cosThSig_mc,m_cosThSigbar_mc,m_pSig_mc,m_pSigbar_mc;
        NTuple::Item<double> m_mxi0_mc,m_mxi0rec_mc,m_mSig_mc,m_mSigrec_mc;
        NTuple::Item<double> m_mxi0b_mc,m_mxi0rec;
        NTuple::Item<double> m_plambda_mc;
        NTuple::Item<double> m_mlambda_mc;
        NTuple::Item<double> m_plamb_mc;
        NTuple::Item<long>  m_nTrkm;
        NTuple::Item<long>  m_nTrkp;
        NTuple::Item<long>  m_ntot;
        NTuple::Item<long>  m_ntot02;
        NTuple::Item<long>  m_ntot03;
        NTuple::Item<long>  m_nGam;
        NTuple::Array<double>  m_Eg;   

        NTuple::Item<double>  m_eg1,m_eg2,m_mgg,m_ntang01,m_ntang02;
        NTuple::Item<double>  m_pipirec01,m_pipirec02;
        NTuple::Item<double>  m_m2pi0rec,m_mpi01,m_mpi02;
        NTuple::Item<double>  m_pglam;
        NTuple::Item<double>  m_msigma0;


        NTuple::Item<long>     m_npi0;
        NTuple::Item<double>   m_ppi0;
        //NTuple::Array<double>  m_mxi0,m_mxi0b,m_pxi0,m_pxi0b,m_chi2pi0, m_mpi0,m_cosgg, m_Vig1, m_Vig2, m_delE_xi0b, m_delE_xi0, m_mbcXi0b, m_mbcXi0, m_mxi0Recoil, m_mxi0bRecoil, m_cosThXi0, m_cosThXi0b;
        NTuple::Item<double>  m_mxi0Recoil,m_mSigRecoil,m_mpi0;
        //NTuple::Array<double>  m_cosxi0, m_cosxi0Recoil,m_pxi0Recoil,m_cosxi0b;
        NTuple::Item<double>  m_cosxi0, m_cosxi0Recoil;
        NTuple::Item<double>  m_cosSig, m_cosSigRecoil;
        NTuple::Item<double>  m_mpip_rec;
        //NTuple::Array<double>  m_rhoxi0Recoil,m_rhoxi0bRecoil;
        NTuple::Item<double>  m_pxi0Recoil,m_pSigRecoil;
        NTuple::Item<double>  m_m4grec01,m_m4grec02;
        NTuple::Item<double>  m_ang,m_ang_xi0;
        //NTuple::Item<double>  m_mpi0;
        //NTuple::Item<double>  m_mpi02;
        NTuple::Item<double>  m_mlamrec,m_plamrec;
        NTuple::Item<double>  m_plambda;

        NTuple::Item<double>  m_mxi_star, m_mxi_star_rec,m_cos_star_rec,m_pxi_star, m_cos_star;
        NTuple::Item<double>  m_mlamb;
        NTuple::Item<double>  m_mlambda;
        NTuple::Item<double>  m_mLamLamb;

        //NTuple::Item<double>  m_mxi0b,m_pxi0brec,m_mxi0brec,m_pxi0b;
        NTuple::Item<double>   m_mlam01, m_mlam02,m_plam01,m_plam02;
        NTuple::Item<double>   m_costheta_lam01,m_costheta_lam02,m_lengthlam, m_lengthlamb;
        NTuple::Item<double>   m_mxi0,m_mSig,m_mxi0b,m_pxi0,m_pSig,m_pxi0b,m_costheta_mxi0,m_costheta_mxi0b;
        NTuple::Item<double>   m_dlt2pi0,m_dltt,m_chisq,m_chisq01;


        NTuple::Item<long>     m_npi0pi0;
        NTuple::Array<double>  m_mpi0pi0recoil;
        NTuple::Matrix<double>  m_ipi0pi0;


        // --- cosTh
        //NTuple::Tuple*  m_tuple;
        TH1D* myHis_pxi0_mc;
        TH1D* myHis_pxi0b_mc;
        TH1D* myHis_pSig_mc;
        TH1D* myHis_pSigbar_mc;
        TH1D* myHis_mSig_mc;
        TH1D* myHis_mSigbar_mc;
        TH1D* myHis_mSigrec_mc;
        TH1D* myHis_mxi0_mc;
        TH1D* myHis_mxi0rec_mc;
        TH1D* myHis_cosThXi0_mc;
        TH1D* myHis_cosThXi0b_mc;
        TH1D* myHis_cosThSig_mc;
        TH1D* myHis_cosThSigbar_mc;

        TH1D* hmxi0_mc;
        TH1D* hpxi0_mc;
        TH1D* hpxi0b_mc;
        TH1D* hplambda_mc;
        TH1D* hplamb_mc;
        TH1D* hcosThXi0_mc;
        TH1D* hcosThXi0b_mc;
        TH1D* hmxi0rec_mc;
        TH1D* hpxi0rec_mc;
        TH1D* hmxsim_mc;
        TH1D* hmxsibarp_mc;

};
#endif
