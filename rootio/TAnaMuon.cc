#include "TAnaMuon.hh"
#include <iostream>
#include <bitset>

ClassImp(TAnaMuon)

using namespace std;

TAnaMuon::TAnaMuon(int index) { 
  fIndex = index;
}

// ----------------------------------------------------------------------
void TAnaMuon::clear() {
  fMuonChi2 = fMuonZ =  fMuonR = fCaloComp = fSegmComp = -999.;
  fTimeInOut = fTimeInOutE = fTimeOutIn = fTimeOutInE = -999.;
  fTimeNdof = -999;
  
  fNstTracks.clear();
  
  fMuonTrackPosAtM1.SetXYZ(-1e30,-1e30,-1e30);
  fMuonTrackPlabAtM1.SetXYZ(-1e30,-1e30,-1e30);

  fInnerPlab.SetXYZ(-1e30,-1e30,-1e30);
  fOuterPlab.SetXYZ(-1e30,-1e30,-1e30);
  fGlobalPlab.SetXYZ(-1e30,-1e30,-1e30);

  fPositionAtM1.SetXYZ(-1e30,-1e30,-1e30); 
  fPositionAtM2.SetXYZ(-1e30,-1e30,-1e30);
  
  fItrkValidFraction = fGtrkNormChi2 = fSegmentComp = fChi2LocalPosition = fChi2LocalMomentum = -999.;
  fNumberOfLostTrkHits = fNumberOfValidTrkHits = -99;

  fGtrkProb = -99.;

  fInt1 = fInt2 = fInt3 -99; 
  fDouble1 = fDouble2 = fDouble3 = -99.;

  for (unsigned int i = 0; i < NXPTRACKS; ++i) {
    fXpTracks[i].dist = 99.;
  }

}


// ----------------------------------------------------------------------
void TAnaMuon::dump() {
  //   cout << " q = " << fQ
  //        << " p = " << fPlab.Mag() 
  //        << " f = " << fPlab.Phi() 
  //        << " t = " << fPlab.Theta();

  cout << "Muon: idx = " << fIndex;
  cout << Form(" midx=%3d", fMuIndex) << " ID=" << static_cast<std::bitset<32> >(fMuID) << endl;

  if (fPlab.Perp() > 0) {
    cout << Form(" P:      q=%+2d pT=%6.2f f=%+4.3f eta=%+4.3f ", fQ, fPlab.Perp(), fPlab.Phi(), fPlab.Eta()) << endl;
  }
  if (fInnerPlab.Perp() > 0) {
    cout << Form(" inner:  q=%+2d pT=%6.2f f=%+4.3f eta=%+4.3f ", fQ, fInnerPlab.Perp(), fInnerPlab.Phi(), fInnerPlab.Eta()) << endl;
  } 
  if (fOuterPlab.Perp() > 0) {
    cout << Form(" outer:  q=%+2d pT=%6.2f f=%+4.3f eta=%+4.3f ", fQ, fOuterPlab.Perp(), fOuterPlab.Phi(), fOuterPlab.Eta()) << endl;
  } 
  if (fGlobalPlab.Perp() > 0) {
    cout << Form(" global: q=%+2d pT=%6.2f f=%+4.3f eta=%+4.3f ", fQ, fGlobalPlab.Perp(), fGlobalPlab.Phi(), fGlobalPlab.Eta()) << endl;
  } 
  cout << Form(" X=%+4.1f z=%4.1f r=%+4.1f cl=%4.1f c2=%4.1f", fMuonChi2, fMuonZ, fMuonR, fCaloComp, fSegmComp);
  cout << Form(" DT=%3d CSC=%3d RPC=%3d t=%4.1f+/-%4.1f (%3d)", fNhitsDT, fNhitsCSC, fNhitsRPC, fTimeInOut, fTimeInOutE, fTimeNdof);
  cout << endl;

  if (fPositionAtM1.Perp() > 0) {
    cout << " Position at M1: x1=" << fPositionAtM1.X() << " y1=" << fPositionAtM1.Y() << " z1=" << fPositionAtM1.Z() << endl;
  }
  if (fPositionAtM2.Perp() > 0) {
    cout << " Position at M2: x2=" << fPositionAtM2.X() << " y2=" << fPositionAtM2.Y() << " z2=" << fPositionAtM2.Z() << endl;
  }
  
  if (fMuonTrackPosAtM1.X() > -1e30 && fMuonTrackPosAtM1.Y() > -1e30 && fMuonTrackPosAtM1.Z() > -1e30) {
	  cout << " Position at M1 of muon track: x=" << fMuonTrackPosAtM1.X() << " y=" << fMuonTrackPosAtM1.Y() << " z=" << fMuonTrackPosAtM1.Z() << endl;
  }
  
  if (fMuonTrackPlabAtM1.X() > -1e30 && fMuonTrackPlabAtM1.Y() > -1e30 && fMuonTrackPlabAtM1.Z() > -1e30) {
	  cout << " Momentum at M1 of muon track: x=" << fMuonTrackPlabAtM1.X() << " y=" << fMuonTrackPlabAtM1.Y() << " z=" << fMuonTrackPlabAtM1.Z() << endl;
  }
  
  // std::map<int,float> fNstTracks; // usage: (trackIx, doca)
  cout << " closeby tracks (" << fNstTracks.size() << ") are" << endl;
  for (std::map<int,float>::const_iterator it = fNstTracks.begin(); it != fNstTracks.end(); ++it)
	  cout << "   (ix=" << it->first << ", doca=" << it->second << ")" << endl;

  //   if (fMCID != -99999) {
  //     cout << Form(" mcid =%+6d", fMCID);
  //   }
  
  //   if (fGenIndex > -1) {
  //     cout << " mcidx = " << fGenIndex;
  //   }
  
  //   if (fMuID > -1) {
  //     cout << Form(" mid=%3x", fMuID);
  //   }
  
  //   cout << endl;

}
  
