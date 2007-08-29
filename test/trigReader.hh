#ifndef TRIGREADER_H
#define TRIGREADER_H

#include <iostream>

#include <TROOT.h>
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include "../rootio/TAna00Event.hh"
#include "../rootio/TGenCand.hh"
#include "../rootio/TAnaCand.hh"
#include "../rootio/TAnaTrack.hh"
#include "../rootio/TAnaVertex.hh"


using namespace std;

#define DR      57.29577951
#define PIPMASS 0.13957
#define ELMASS  0.000511
#define MUMASS  0.10567


class trigReader {
public:
  trigReader(TChain *tree, TString evtClassName);
  ~trigReader();
  void        init(TString evtClassName);

  void        openHistFile(TString filename);
  void        closeHistFile();
  void        bookHist();

  void        startAnalysis();
  int         loop(int nevents = 1, int start = -1);
  void        eventProcessing();

  void        dumpEvent();
  void        fillHist();
  int         matchingTrack(int genIndex, double &dR2, double &dpt);

private:

  TChain     *fpChain;        // pointer to the analyzed TTree or TChain
  TFile      *fpHistFile;     // for output histograms and reduced trees
  TString     fChainFileName; // the name of the chain file
  int         fNentries;      // number of events in chain
  int         fEvent;         // current event number

  TAna00Event*fpEvt; 

  // -- Variables
  int        fRun; 


  // -- Histogram pointers 
  TTree *fTree;

};

// ----------------------------------------------------------------------
inline void mk4Vector(TLorentzVector &p4, const Double_t p, const Double_t t, const Double_t f, const Double_t m) {
  p4.SetXYZM(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t), m);
}

#endif
