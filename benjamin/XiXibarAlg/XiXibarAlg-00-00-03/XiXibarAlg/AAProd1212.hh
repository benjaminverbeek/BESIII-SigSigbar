#ifndef AA_Prod1212_H_
#define AA_Prod1212_H_
#include "TMath.h"
class AAProd1212 {
  public:

  Double_t C(Int_t i, Int_t j) const {return thC[i][j];}
  AAProd1212(double aP, double pP, double th) {
    // decay 1/2 -> 1/2 + 0
    Double_t v=TMath::Sqrt(1-aP*aP);
    Double_t gP=v*TMath::Cos(pP);
    Double_t bP=v*TMath::Sin(pP);
    Double_t st=TMath::Sin(th);
    Double_t ct=TMath::Cos(th);
    Double_t ct2=ct*ct;
    
    for(int k=0;k<4;k++){
      for(int l=0;l<4;l++){
         Double_t c=0;
         if(k==0&&l==0)c=1+aP*ct2;
         if(k==0&&l==2)c=bP*st*ct;
         if(k==2&&l==0)c=-bP*st*ct;
         if(k==1&&l==1)c=1-ct2;
         if(k==1&&l==3)c=gP*st*ct;
         if(k==3&&l==1)c=-gP*st*ct;
         if(k==3&&l==3)c=-aP-ct2;
         if(k==2&&l==2)c=aP*(1-ct2);
         thC[k][l]=c;
      }
    }
    
  }
  ~AAProd1212() {}
private:

  double thC[4][4];
};
#endif /* AA_Prod1212_H_ */
