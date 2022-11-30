#ifndef AA_Decay12_H_
#define AA_Decay12_H_
#include "TMath.h"
class AADecay12 {
  public:

  AADecay12(double aM, double pM, double th, double ph) {
    // decay 1/2 -> 1/2 + 0
    Double_t v=TMath::Sqrt(1-aM*aM);
    Double_t gM=v*TMath::Cos(pM);
    Double_t bM=v*TMath::Sin(pM);
    Double_t sp=TMath::Sin(ph);
    Double_t cp=TMath::Cos(ph);
    Double_t st=TMath::Sin(th);
    Double_t ct=TMath::Cos(th);
    for(int k=0;k<4;k++){
      for(int l=0;l<4;l++){
        Double_t a=0;

        if(k==0&&l==0)a= 1;
        if(k==0&&l==3)a= aM;
        if(k==1&&l==0)a= aM*cp*st;
        if(k==1&&l==1)a= gM*ct*cp-bM*sp;
        if(k==1&&l==2)a=-bM*ct*cp-gM*sp;
        if(k==1&&l==3)a= cp*st;
        if(k==2&&l==0)a= aM*sp*st;
        if(k==2&&l==1)a= bM*cp+gM*ct*sp;
        if(k==2&&l==2)a= gM*cp-bM*ct*sp;
        if(k==2&&l==3)a= sp*st;
        if(k==3&&l==0)a= aM*ct;
        if(k==3&&l==1)a=-gM*st;
        if(k==3&&l==2)a= bM*st;
        if(k==3&&l==3)a= ct;

        tHa[k][l]=a;
      }
    }
  }      
  ~AADecay12() {}
  Double_t A(Int_t i, Int_t j) const {return tHa[i][j];}
private:

  double tHa[4][4];
};
#endif
