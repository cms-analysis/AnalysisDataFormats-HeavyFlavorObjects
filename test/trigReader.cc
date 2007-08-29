#include "trigReader.hh"

#include <cmath>

#include "TRandom.h"

#include "trigReader.icc"


// Run with: ./bmm -c chains/bg-test -D root ; ./bmm -c chains/sg-test -D root ;

// ----------------------------------------------------------------------
void trigReader::startAnalysis() {
  cout << "startAnalysis: ..." << endl;
}



// ----------------------------------------------------------------------
void trigReader::eventProcessing() {
  TAnaTrack *pTrack;

  // -- generated muons vs. reconstructed muons, trigger efficiency
  TGenCand *pGen; 
  double gpt, rpt; 
  double dR2, dpt;
  for (int ig = 0; ig < fpEvt->nGenCands(); ++ig) {
    pGen = fpEvt->getGenCand(ig); 
    if (TMath::Abs(pGen->fID) == 13) {
      gpt = pGen->fP.Perp();
      ((TH1D*)fpHistFile->Get("h100"))->Fill(gpt); 
      int isig = matchingTrack(ig, dR2, dpt);
      if (isig > -1) {
	pTrack = fpEvt->getSigTrack(isig);
	rpt = pTrack->fPlab.Perp(); 
	((TH1D*)fpHistFile->Get("h102"))->Fill(rpt - gpt); 
	((TH1D*)fpHistFile->Get("h101"))->Fill(rpt); 

	if (fpEvt->fTriggerWord1) {
	  ((TH1D*)fpHistFile->Get("h110"))->Fill(gpt); 
	}

	if (fpEvt->fTriggerWord2) {
	  ((TH1D*)fpHistFile->Get("h111"))->Fill(gpt); 
	}

	if (fpEvt->fTriggerWord3) {
	  ((TH1D*)fpHistFile->Get("h112"))->Fill(gpt); 
	}

// 	cout << "Match track to muon. gpt = " << gpt 
// 	     << " pt = " << rpt 
// 	     << " gphi = " << pGen->fP.Phi() 
// 	     << " phi = " << pTrack->fPlab.Phi() 
// 	     << " geta = " << pGen->fP.Eta() 
// 	     << " eta = " << pTrack->fPlab.Eta() 
// 	     << endl;
	  
      } else {
	((TH1D*)fpHistFile->Get("h90"))->Fill(dR2); 
	((TH1D*)fpHistFile->Get("h91"))->Fill(dpt); 
	((TH1D*)fpHistFile->Get("h92"))->Fill(gpt); 
	((TH1D*)fpHistFile->Get("h93"))->Fill(pGen->fP.Eta()); 
      }
    }	

  }
 

  fpHistFile->cd();
  fillHist(); 
  fTree->Fill();

}


// ----------------------------------------------------------------------
void trigReader::fillHist() {
  
  
}

// ----------------------------------------------------------------------
void trigReader::bookHist() {
  TH1 *h;
  cout << "-->bookHist> " << endl;

  h = new TH1D("h90", "dR2 for failed matches", 40, 0., 2.);
  h = new TH1D("h91", "dPt for failed matches", 40, 0., 2.);
  h = new TH1D("h92", "pt  for failed matches", 40, 0., 20.);
  h = new TH1D("h93", "eta for failed matches", 40, -5., 5.);

  int nbins(40); 
  double ptmax(20); 
  h = new TH1D("h100", "pt (muons) generator-level", nbins, 0., ptmax);
  h = new TH1D("h101", "pt (muons) reconstructed",   nbins, 0., ptmax);
  h = new TH1D("h102", "pt rec - pt gen (muons)",    nbins, -1., 1.0);
  h = new TH1D("h110", "pt (muons) relaxed dimuon",  nbins, 0., ptmax);
  h = new TH1D("h111", "pt (muons) displaced J/psi", nbins, 0., ptmax);
  h = new TH1D("h112", "pt (muons) displaced mmh",   nbins, 0., ptmax);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",    &fRun ,"run/I");

}


// ----------------------------------------------------------------------
int trigReader::matchingTrack(int genIndex, double &dR2, double &dpt) {

  static double dR2match = 0.2; 
  static double dptmatch = 0.3;

  dR2 = -99.;
  dpt = -99.;

  TGenCand *pGen = fpEvt->getGenCand(genIndex); 
  double gpt(pGen->fP.Perp()), gphi(pGen->fP.Phi()), geta(pGen->fP.Eta()); 
  double dphi, deta;
  TAnaTrack *pTrack; 
  for (int is = 0; is < fpEvt->nSigTracks(); ++is) {
    pTrack = fpEvt->getSigTrack(is);
    dpt  = TMath::Abs(gpt - pTrack->fPlab.Perp()); 
    dphi = gphi - pTrack->fPlab.Phi(); 
    while (dphi >= M_PI) dphi -= M_PI;
    while (dphi < -M_PI) dphi += 2.*M_PI; 

    deta = geta - pTrack->fPlab.Eta(); 
    dR2 = deta*deta + dphi*dphi; 
    if ((dR2 < dR2match) && (dpt < dptmatch)) {
      return is; 
    }
  }

  return -1; 

}
