#include "TH1.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TCanvas.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;


#pragma link C++ function printNonZero(TH1*);
#pragma link C++ function stampAndSave(TCanvas *, const char*);
#pragma link C++ function shrinkPad(double, double, double, double);
#pragma link C++ function babar(double, double, double);
#pragma link C++ function zone(int, int, TCanvas *);
#pragma link C++ function wait();
#pragma link C++ function colors();

#pragma link C++ function setTitles(TH1*,const char*,const char*,float,float,float,float,int);
#pragma link C++ function setHist(TH1*,int,int,double,double);
#pragma link C++ function setGraph(TGraph*,int,int,double,double);
#pragma link C++ function setFilledHist(TH1 *,int,int,int,int);
#pragma link C++ function setMaximum(double, TH1 *, TH1 *);
#pragma link C++ function unmix(TH1D *, TH1D *, double);

// #pragma link C++ function scanParameter(int, int, TMinuit &, void (*func)(int &, double *, double &, double *, int)); 

#pragma link C++ function chi2Prob(double, double);
#pragma link C++ function chi2Test(TH1*, TH1*, double&, double&, int);
#pragma link C++ function chi2TestErr(TH1*, TH1*, double&, double&, int);

#pragma link C++ function dEff(int, int);
#pragma link C++ function dEff(int, double);
#pragma link C++ function dEff(double, int);
#pragma link C++ function dEff(double, double, double, double);

#pragma link C++ function dBinomial(double, double);
#pragma link C++ function dRatio(double, double);
#pragma link C++ function dRatio(double, double, double, double);
#pragma link C++ function dBF(double, double, double, double, double, double);

#pragma link C++ function getError(TH1*);
#pragma link C++ function getErrorRange(TH1*, int, int);


#endif

