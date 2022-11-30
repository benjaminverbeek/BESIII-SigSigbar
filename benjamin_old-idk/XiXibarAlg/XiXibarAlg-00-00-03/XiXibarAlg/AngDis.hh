#ifndef MN_AngDis_H_
#define MN_AngDis_H_
#include "AAProd1212.hh"
#include "AADecay12.hh"
#include <cmath>
#include "TMath.h"

using namespace std;

class AngDis {

  public:

    AngDis(double ap, double pp, double aXi, double pXi, double aL, double aXib, double aLb) :
      tHaPsi(ap),tHpPsi(pp),tHaXi(aXi), tHpXi(pXi), tHaL(aL), tHaXib(aXib), tHaLb(aLb)
  {}      
    ~AngDis() {}
    double  operator()(double th,double th1L,double ph1L,double th1p,
           double ph1p,double th2L,double ph2L,double th2p,double ph2p) 

    {
      AAProd1212 prod(tHaPsi,tHpPsi,th);
      AADecay12  dXi(tHaXi,tHpXi,th1L,ph1L);
      AADecay12  dXib(tHaXib,-tHpXi,th2L,ph2L);
      AADecay12  dL(tHaL,0,th1p,ph1p);
      AADecay12  dLb(tHaLb,0,th2p,ph2p);

      double tep=0;
      for(int i1=0;i1<4;i1++){// Xi loop
         for(int i2=0;i2<4;i2++){// Xi_bar loop
             for(int j1=0;j1<4;j1++){// L loop
                 for(int j2=0;j2<4;j2++){// L_bar loop
                   tep += prod.C(i1,i2)*dXi.A(i1,j1)*dXib.A(i2,j2)*
                          dL.A(j1,0)*dLb.A(j2,0);
      }
    }
  }
      }
      return tep;
    }



private:

double tHaPsi;
double tHpPsi;
double tHaXi;
double tHpXi;
double tHaL;
double tHaXib;
double tHaLb;
};

#endif // MN_AngDis_H_

