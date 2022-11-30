#ifndef Physics_Analysis_XiXibar_H
#define Physics_Analysis_XiXibar_H

#include "VertexFit/VertexFit.h"
#include "VertexFit/SecondVertexFit.h" 

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "GaudiKernel/AlgFactory.h"
#include "TClonesArray.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/IPartPropSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "TLorentzVector.h"
#include "TH2.h"
#include "TTree.h"
//#include "TrackCalibration/TrackCalibration.h"
//#include "TrackCorrection/TrackCorrection.h"

class XiXibar : public Algorithm
{
    public:
        XiXibar(const std::string& name, ISvcLocator* pSvcLocator);
        StatusCode initialize();
        StatusCode execute();
        StatusCode finalize();

        void InitializeReconstructedSpectra();
        void FinalizeReconstructedSpectra();
        StatusCode initialiezeNtuple();

        HepLorentzVector GetEcmsP4();
        void KinematicFit(std::vector<WTrackParameter> tracks, TString constraint);
        double GetRxy(int itr);
        double TimeOfFlight(int itrack, int itype );
        bool WeighDataTrue();
        double Angles(CLHEP::HepLorentzVector& cm, int i);
        //TLorentzVector Helrotate(TLorentzVector p1, double phi, double theta);
        CLHEP::HepLorentzVector Helrotate(CLHEP::HepLorentzVector& p1, double phi, double theta);
        std::vector<double> HelicityAnglesXi(std::vector<TLorentzVector> fourvec);

    private:

        HepLorentzVector HepXiTrue;
        HepLorentzVector HepXibarTrue;
        HepLorentzVector HepLambdaTrue;
        HepLorentzVector HepLambdabarTrue;
        HepLorentzVector HepProtonTrue;
        HepLorentzVector HepPimLamTrue;
        HepLorentzVector HepPimXiTrue;
        HepLorentzVector HepProtonbarTrue;
        HepLorentzVector HepPipLambarTrue;
        HepLorentzVector HepPipXibarTrue;  

        //  Vertex Database Service
        // IVertexDbSvc* vertexService;
        // IPartPropSvc* particlePropertiesService;

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
        NTuple::Tuple*  m_tuple1;

        // mc true

        NTuple::Item<int>   mc_idxmc;
        NTuple::Array<int>  mc_trkidx;
        NTuple::Array<int>  mc_pdgid;
        NTuple::Array<int>  mc_motheridx;

        NTuple::Item<double>    mc_Xim_e;
        NTuple::Item<double>    mc_Xim_px;
        NTuple::Item<double>    mc_Xim_py;
        NTuple::Item<double>    mc_Xim_pz;

        NTuple::Item<double>    mc_Xip_e;
        NTuple::Item<double>    mc_Xip_px;
        NTuple::Item<double>    mc_Xip_py;
        NTuple::Item<double>    mc_Xip_pz;

        NTuple::Item<double>    mc_Lam_e;
        NTuple::Item<double>    mc_Lam_px;
        NTuple::Item<double>    mc_Lam_py;
        NTuple::Item<double>    mc_Lam_pz;

        NTuple::Item<double>    mc_Lambar_e;
        NTuple::Item<double>    mc_Lambar_px;
        NTuple::Item<double>    mc_Lambar_py;
        NTuple::Item<double>    mc_Lambar_pz;

        NTuple::Item<double>    mc_pr_e;
        NTuple::Item<double>    mc_pr_px;
        NTuple::Item<double>    mc_pr_py;
        NTuple::Item<double>    mc_pr_pz;

        NTuple::Item<double>    mc_prbar_e;
        NTuple::Item<double>    mc_prbar_px;
        NTuple::Item<double>    mc_prbar_py;
        NTuple::Item<double>    mc_prbar_pz;

        NTuple::Item<double>    mc_pimLam_e;
        NTuple::Item<double>    mc_pimLam_px;
        NTuple::Item<double>    mc_pimLam_py;
        NTuple::Item<double>    mc_pimLam_pz;

        NTuple::Item<double>    mc_pipLambar_e;
        NTuple::Item<double>    mc_pipLambar_px;
        NTuple::Item<double>    mc_pipLambar_py;
        NTuple::Item<double>    mc_pipLambar_pz;

        NTuple::Item<double>    mc_pimXim_e;
        NTuple::Item<double>    mc_pimXim_px;
        NTuple::Item<double>    mc_pimXim_py;
        NTuple::Item<double>    mc_pimXim_pz;

        NTuple::Item<double>    mc_pipXip_e;
        NTuple::Item<double>    mc_pipXip_px;
        NTuple::Item<double>    mc_pipXip_py;
        NTuple::Item<double>    mc_pipXip_pz;

        // mc rec
        NTuple::Item<int>   m_idxmc;
        NTuple::Array<int>  m_trkidx;
        NTuple::Array<int>  m_pdgid;
        NTuple::Array<int>  m_motheridx;

        NTuple::Item<long>   m_nTrkm;
        NTuple::Item<long>   m_nTrkp;

        NTuple::Item<long>   m_npr_cand;
        NTuple::Item<long>   m_nprbar_cand;
        NTuple::Item<long>   m_npim_cand;
        NTuple::Item<long>   m_npip_cand;

        NTuple::Item<long>   m_run;
        NTuple::Item<long>   m_rec;

        NTuple::Item<double>    lam_decl;
        NTuple::Item<double>    lambar_decl;
        NTuple::Item<double>    xim_decl;
        NTuple::Item<double>    xip_decl;
        NTuple::Item<double>    lam_decl_err;
        NTuple::Item<double>    lambar_decl_err;
        NTuple::Item<double>    xim_decl_err;
        NTuple::Item<double>    xip_decl_err;
        NTuple::Item<double>    lam_prvxfitchi2;
        NTuple::Item<double>    lambar_prvxfitchi2;
        NTuple::Item<double>    lam_vxfitchi2;
        NTuple::Item<double>    lambar_vxfitchi2;
        NTuple::Item<double>    xim_vxfitchi2;
        NTuple::Item<double>    xip_vxfitchi2;

        NTuple::Item<double>    pr_vx_e;
        NTuple::Item<double>    pr_vx_px;
        NTuple::Item<double>    pr_vx_py;
        NTuple::Item<double>    pr_vx_pz;

        NTuple::Item<double>    prbar_vx_e;
        NTuple::Item<double>    prbar_vx_px;
        NTuple::Item<double>    prbar_vx_py;
        NTuple::Item<double>    prbar_vx_pz;

        NTuple::Item<double>    pimLam_vx_e;
        NTuple::Item<double>    pimLam_vx_px;
        NTuple::Item<double>    pimLam_vx_py;
        NTuple::Item<double>    pimLam_vx_pz;

        NTuple::Item<double>    pipLambar_vx_e;
        NTuple::Item<double>    pipLambar_vx_px;
        NTuple::Item<double>    pipLambar_vx_py;
        NTuple::Item<double>    pipLambar_vx_pz;

        NTuple::Item<double>    pimXim_vx_e;
        NTuple::Item<double>    pimXim_vx_px;
        NTuple::Item<double>    pimXim_vx_py;
        NTuple::Item<double>    pimXim_vx_pz;

        NTuple::Item<double>    pipXip_vx_e;
        NTuple::Item<double>    pipXip_vx_px;
        NTuple::Item<double>    pipXip_vx_py;
        NTuple::Item<double>    pipXip_vx_pz;

        NTuple::Item<double>    pr_rec_e;
        NTuple::Item<double>    pr_rec_px;
        NTuple::Item<double>    pr_rec_py;
        NTuple::Item<double>    pr_rec_pz;

        NTuple::Item<double>    prbar_rec_e;
        NTuple::Item<double>    prbar_rec_px;
        NTuple::Item<double>    prbar_rec_py;
        NTuple::Item<double>    prbar_rec_pz;

        NTuple::Item<double>    pimLam_rec_e;
        NTuple::Item<double>    pimLam_rec_px;
        NTuple::Item<double>    pimLam_rec_py;
        NTuple::Item<double>    pimLam_rec_pz;

        NTuple::Item<double>    pipLambar_rec_e;
        NTuple::Item<double>    pipLambar_rec_px;
        NTuple::Item<double>    pipLambar_rec_py;
        NTuple::Item<double>    pipLambar_rec_pz;

        NTuple::Item<double>    pimXim_rec_e;
        NTuple::Item<double>    pimXim_rec_px;
        NTuple::Item<double>    pimXim_rec_py;
        NTuple::Item<double>    pimXim_rec_pz;

        NTuple::Item<double>    pipXip_rec_e;
        NTuple::Item<double>    pipXip_rec_px;
        NTuple::Item<double>    pipXip_rec_py;
        NTuple::Item<double>    pipXip_rec_pz;

        NTuple::Item<double>    Lam_vx_e;
        NTuple::Item<double>    Lam_vx_px;
        NTuple::Item<double>    Lam_vx_py;
        NTuple::Item<double>    Lam_vx_pz;

        NTuple::Item<double>    Lambar_vx_e;
        NTuple::Item<double>    Lambar_vx_px;
        NTuple::Item<double>    Lambar_vx_py;
        NTuple::Item<double>    Lambar_vx_pz;

        NTuple::Item<double>    Xim_vx_e;
        NTuple::Item<double>    Xim_vx_px;
        NTuple::Item<double>    Xim_vx_py;
        NTuple::Item<double>    Xim_vx_pz;

        NTuple::Item<double>    Xip_vx_e;
        NTuple::Item<double>    Xip_vx_px;
        NTuple::Item<double>    Xip_vx_py;
        NTuple::Item<double>    Xip_vx_pz;

        NTuple::Item<double>    Xim_4C_e;
        NTuple::Item<double>    Xim_4C_px;
        NTuple::Item<double>    Xim_4C_py;
        NTuple::Item<double>    Xim_4C_pz;

        NTuple::Item<double>    Xip_4C_e;
        NTuple::Item<double>    Xip_4C_px;
        NTuple::Item<double>    Xip_4C_py;
        NTuple::Item<double>    Xip_4C_pz;

        NTuple::Item<double>    Xim_pull_x;
        NTuple::Item<double>    Xim_pull_y;
        NTuple::Item<double>    Xim_pull_z;
   
        NTuple::Item<double>    Xip_pull_x;
        NTuple::Item<double>    Xip_pull_y;
        NTuple::Item<double>    Xip_pull_z;

        NTuple::Item<double>    chi2_4C;

        NTuple::Item<double>    pr_vx_x;
        NTuple::Item<double>    pr_vx_y;
        NTuple::Item<double>    pr_vx_z;

        NTuple::Item<double>    prbar_vx_x;
        NTuple::Item<double>    prbar_vx_y;
        NTuple::Item<double>    prbar_vx_z;

        NTuple::Item<double>    pim_vx_x;
        NTuple::Item<double>    pim_vx_y;
        NTuple::Item<double>    pim_vx_z;

        NTuple::Item<double>    pip_vx_x;
        NTuple::Item<double>    pip_vx_y;
        NTuple::Item<double>    pip_vx_z;        

        NTuple::Item<double>    Lam_vx_x;
        NTuple::Item<double>    Lam_vx_y;
        NTuple::Item<double>    Lam_vx_z;

        NTuple::Item<double>    Lambar_vx_x;
        NTuple::Item<double>    Lambar_vx_y;
        NTuple::Item<double>    Lambar_vx_z;

        NTuple::Item<double>    Xim_vx_x;
        NTuple::Item<double>    Xim_vx_y;
        NTuple::Item<double>    Xim_vx_z;

        NTuple::Item<double>    Xip_vx_x;
        NTuple::Item<double>    Xip_vx_y;
        NTuple::Item<double>    Xip_vx_z;

        NTuple::Item<double>    Rxy_pr;
        NTuple::Item<double>    Rxy_prbar;
        NTuple::Item<double>    Rxy_pimxim;
        NTuple::Item<double>    Rxy_pipxip;
        NTuple::Item<double>    Rxy_pimlam;
        NTuple::Item<double>    Rxy_piplambar;

        NTuple::Item<double>    pr_drho;
        NTuple::Item<double>    pr_phi0;
        NTuple::Item<double>    pr_kappa;
        NTuple::Item<double>    pr_lambda;

        NTuple::Item<double>    pimlam_drho;
        NTuple::Item<double>    pimlam_phi0;
        NTuple::Item<double>    pimlam_kappa;
        NTuple::Item<double>    pimlam_lambda;

        NTuple::Item<double>    pimxim_drho;
        NTuple::Item<double>    pimxim_phi0;
        NTuple::Item<double>    pimxim_kappa;
        NTuple::Item<double>    pimxim_lambda;

        NTuple::Item<double>    prbar_drho;
        NTuple::Item<double>    prbar_phi0;
        NTuple::Item<double>    prbar_kappa;
        NTuple::Item<double>    prbar_lambda;

        NTuple::Item<double>    piplambar_drho;
        NTuple::Item<double>    piplambar_phi0;
        NTuple::Item<double>    piplambar_kappa;
        NTuple::Item<double>    piplambar_lambda;

        NTuple::Item<double>    pipxip_drho;
        NTuple::Item<double>    pipxip_phi0;
        NTuple::Item<double>    pipxip_kappa;
        NTuple::Item<double>    pipxip_lambda;      

        NTuple::Item<double>    prvfit_drho;
        NTuple::Item<double>    prvfit_phi0;
        NTuple::Item<double>    prvfit_kappa;
        NTuple::Item<double>    prvfit_lambda;  

        NTuple::Item<double>    pimlamvfit_drho;
        NTuple::Item<double>    pimlamvfit_phi0;
        NTuple::Item<double>    pimlamvfit_kappa;
        NTuple::Item<double>    pimlamvfit_lambda;

        NTuple::Item<double>    pimximvfit_drho;
        NTuple::Item<double>    pimximvfit_phi0;
        NTuple::Item<double>    pimximvfit_kappa;
        NTuple::Item<double>    pimximvfit_lambda;

        NTuple::Item<double>    prbarvfit_drho;
        NTuple::Item<double>    prbarvfit_phi0;
        NTuple::Item<double>    prbarvfit_kappa;
        NTuple::Item<double>    prbarvfit_lambda;  

        NTuple::Item<double>    piplambarvfit_drho;
        NTuple::Item<double>    piplambarvfit_phi0;
        NTuple::Item<double>    piplambarvfit_kappa;
        NTuple::Item<double>    piplambarvfit_lambda;

        NTuple::Item<double>    pipxipvfit_drho;
        NTuple::Item<double>    pipxipvfit_phi0;
        NTuple::Item<double>    pipxipvfit_kappa;
        NTuple::Item<double>    pipxipvfit_lambda;

        NTuple::Item<double>    true_Xim_e;
        NTuple::Item<double>    true_Xim_px;
        NTuple::Item<double>    true_Xim_py;
        NTuple::Item<double>    true_Xim_pz;

        NTuple::Item<double>    true_Xip_e;
        NTuple::Item<double>    true_Xip_px;
        NTuple::Item<double>    true_Xip_py;
        NTuple::Item<double>    true_Xip_pz;

        NTuple::Item<double>    true_Lam_e;
        NTuple::Item<double>    true_Lam_px;
        NTuple::Item<double>    true_Lam_py;
        NTuple::Item<double>    true_Lam_pz;

        NTuple::Item<double>    true_Lambar_e;
        NTuple::Item<double>    true_Lambar_px;
        NTuple::Item<double>    true_Lambar_py;
        NTuple::Item<double>    true_Lambar_pz;

        NTuple::Item<double>    true_pr_e;
        NTuple::Item<double>    true_pr_px;
        NTuple::Item<double>    true_pr_py;
        NTuple::Item<double>    true_pr_pz;

        NTuple::Item<double>    true_prbar_e;
        NTuple::Item<double>    true_prbar_px;
        NTuple::Item<double>    true_prbar_py;
        NTuple::Item<double>    true_prbar_pz;

        NTuple::Item<double>    true_pimLam_e;
        NTuple::Item<double>    true_pimLam_px;
        NTuple::Item<double>    true_pimLam_py;
        NTuple::Item<double>    true_pimLam_pz;

        NTuple::Item<double>    true_pipLambar_e;
        NTuple::Item<double>    true_pipLambar_px;
        NTuple::Item<double>    true_pipLambar_py;
        NTuple::Item<double>    true_pipLambar_pz;

        NTuple::Item<double>    true_pimXim_e;
        NTuple::Item<double>    true_pimXim_px;
        NTuple::Item<double>    true_pimXim_py;
        NTuple::Item<double>    true_pimXim_pz;

        NTuple::Item<double>    true_pipXip_e;
        NTuple::Item<double>    true_pipXip_px;
        NTuple::Item<double>    true_pipXip_py;
        NTuple::Item<double>    true_pipXip_pz;

        double ximtheta_tr, ximphi_tr;
        double xibarptheta_tr, xibarpphi_tr;

        HepLorentzVector rot_lam_tr, rot_lambar_tr;
        HepLorentzVector rot_pr_tr, rot_prbar_tr;

        double rot_lamtheta_tr, rot_lamphi_tr;
        double lamtheta_tr, lamphi_tr;
        double rot_prtheta_tr, rot_prphi_tr;

        double rot_lambartheta_tr, rot_lambarphi_tr;
        double lambartheta_tr, lambarphi_tr;
        double rot_prbartheta_tr, rot_prbarphi_tr;


        // --- cosTh
        //NTuple::Tuple*  m_tuple;
        TH1D* hmxi0_mc;
        TH1D* hpxi0_mc;
        TH1D* hpxi0b_mc;
        TH1D* hp_lambda_mc;
        TH1D* hp_lambda_bar_mc;
        TH1D* hcosThXi0_mc;
        TH1D* hcosThXi0b_mc;
        TH1D* hp_xim_mc;
        TH1D* hp_xipb_mc;
        TH1D* hmxsim_mc;
        TH1D* hmxsibarp_mc;

        TH1D* hp_pr_lam_vx;
        TH1D* hp_pr_bar_lambar_vx;
        TH1D* hp_pim_lam_vx;
        TH1D* hp_pip_lambar_vx;
        TH1D* hp_pim_Xim_vx;
        TH1D* hp_pip_Xipbar_vx; 

        TH1D* hlam_decaylength;
        TH1D* hlam_bar_decaylength; 

        TH1D* h_chi2_prvx_lam;
        TH1D* h_chi2_prvx_lam_bar;
        TH1D* h_chi2_scvx_lam;
        TH1D* h_chi2_scvx_lam_bar;       

        TH1D* hm_lam_vx;
        TH1D* hm_lam_bar_vx;
        TH1D* hm_xim_vx;
        TH1D* hm_xip_bar_vx;
        TH1D* hm_jpsi_vx;

        TH1D* hweight;
        TH1D* efficiency;

        TH2D* hpvth_pr_lam_vx;
        TH2D* hpvth_pim_lam_vx;
        TH2D* hpvth_pim_Xim_vx;
        TH2D* hpvth_pr_bar_lambar_vx;
        TH2D* hpvth_pip_lambar_vx;
        TH2D* hpvth_pip_Xip_vx;

        TH2D* hm_lam_lambar_vx;
        TH2D* hm_lam_lambar_vx_final;
        TH2D* hm_xim_xipbar_vx;

        TH2D* hm_xim_xipbar_fit;

        TH1D*     hm_tof_proton;  
        TH1D*     hm_tof_pion;

        TH1D*     hm_tof_proton_best;  
        TH1D*     hm_tof_pion_best1;
        TH1D*     hm_tof_pion_best2;


};
#endif
