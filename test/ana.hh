#ifndef ANA
#define ANA

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"


#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TVirtualPad.h"  // access to gPad

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>

class ana: public TObject {

public:

  ana(const char *files = "nada");

  // -- initialization and setup
  // ---------------------------
  void init(const char *files);
  void loadFiles(const char *filename);
  void loadDa(const char *filename, double vxs, const char *sign, const char *type);
  void loadMc(const char *filename, double vxs, const char *sign, const char *type);
  void loadSg(const char *filename, double vxs, const char *sign, const char *type);
  void dumpFiles();
  void dumpCuts();

  // -- main methods
  // --------------
  void makeAllPlots();

  void distributions();
  void showDistributions(int offset = 0, int wait = 0); 
  void showDistribution(const char *hname, int mode = 0, double x = -9999., double y = -9999.);

  void effTables();
  void effTable(TFile *f, const char *tag); 
  double massReduction(const char *s = "c030", const char *tag = "m0");

  // -- temporary 'working' methods
  // -------------------------------
  void loopOptimization(double pt = 4.); 
  void optimizerNumbers(const char *precuts, const char *line, double &eff, double &sg, double &bg);
  void runOptimization(const char *aoCuts, const char *extraVar, int nbin, double min, double max);
  void handOptimization(const char *aoCuts, const char *extraCuts); 

  // -- Utilities and helper methods
  // -------------------------------

  TH1* sumHistMC(const char *hist, int draw = 0, const char *sel = "c0");
  void sumHistMC_Add(int ch, const char *hist, TH1 *hist1, TH1 *hist2);

  void makeCanvas(int i = 3);
  void shrinkPad(double b = 0.1, double l = 0.1, double r = 0.1, double t = 0.1);
  void setTitles(TH1 *h, const char *sx, const char *sy, 
		 float size = 0.05, float xoff = 1.1, float yoff = 1.1, float lsize = 0.05, int font = 132);
  void setTitles2(TH2 *h, const char *sx, const char *sy, 
		 float size = 0.05, float xoff = 1.1, float yoff = 1.1, float lsize = 0.05, int font = 132);
  void setHist(TH1 *h, int color = kBlack, int symbol = 20, double size = 1., double width = 2.);
  void setFilledHist(TH1 *h, int lcol = kBlack, int fcol = kYellow, int fstyle = 1000, int width = 1, int style = 1);
  void emptyBinError(TH1 *h);

  TH1D* DivideHisto(TH1D *hist1, TH1D *hist2);
  TH2D* DivideHisto2(TH2D *hist1, TH2D *hist2);

  TString texForm(double e);
  TString texForm2(double e);
  TString texForm31(double e);
  TString formatTex(double n, const char *name, const char *tag);
  TString scaleFactor(int exp);

  double dBinomial(double n, double N);
  double dEff(int in, int iN);

  int  wait();

  // -- Files for Signal, Data, Monte Carlo, and Control Samples
  int nSg, nMc, nDa;
  TFile *fS[10], *fD[10], *fM[30], *fCS[10];

  int fShow; 
  TString fFile; 
  TString fNumbersFileName;


private:

  // -- for normalization
  double fvXsS[10],  fvXsD[10],  fvXsM[30];  // this is entered by hand from NSEL/NGEN
  double fLumiS[10], fLumiD[10], fLumiM[30]; // this is computed in loadFiles()
  double fNevtS[10], fNevtD[10], fNevtM[30]; // this is filled in loadFiles()
  double fNexpS[10], fNexpD[10], fNexpM[30]; // this is filled in loadFiles()
  double fMisIdM[30];                         // this is filled in muonMisId()
  TString fSignS[10], fSignD[10], fSignM[30]; // this is filled in loadFiles()
  TString fTypeS[10], fTypeD[10], fTypeM[30]; // this is filled in loadFiles()
  TString fSignTexS[10], fSignTexD[10], fSignTexM[30]; // this is filled in loadFiles()
  TString fSignLeggS[10], fSignLeggD[10], fSignLeggM[30]; // this is filled in loadFiles()


  int sgIndex, bgIndex, normSgIndex, normBgIndex; 

  // -- misc
  double fFom;

  double fMassBs, fMassBp;

  // -- functions
  TF1 *f0, *f1, *f2, *f3, *f4, *f5, *f6; 
  TF1 *f10, *f11;

  // -- Display utilities
  int fFont; 
  TCanvas *c0, *c1, *c2, *c3, *c4, *c5;
  TLatex *tl; 
  TBox *box;
  TArrow *pa; 
  TLine *pl; 
  TLegend *legg;
  TLegendEntry *legge;

  char inDir[200];
  char outDir[200];
  char line[200];

  char fSuffix[100]; 

  ClassDef(ana,1) //Testing ana
};


#endif

