#include "TAna01Event.hh"
#include <iostream>

ClassImp(TAna01Event)
 
using namespace std;

// ----------------------------------------------------------------------
TAna01Event::TAna01Event(Int_t Option) {
  fGenCands        = new TClonesArray("TGenCand", 1000);
  fnGenCands       = 0;

  fRecTracks       = new TClonesArray("TAnaTrack", 1000);
  fnRecTracks      = 0;

  fMuons           = new TClonesArray("TAnaMuon", 1000);
  fnMuons          = 0;

  fTrgObj          = new TClonesArray("TTrgObj", 1000);
  fnTrgObj         = 0;

  fSigTracks       = new TClonesArray("TAnaTrack", 1000);
  fnSigTracks      = 0;

  fCaloJets        = new TClonesArray("TAnaJet", 1000);
  fnCaloJets       = 0;

  fGenJets         = new TClonesArray("TAnaJet", 1000);
  fnGenJets        = 0;

  fTrackJets       = new TClonesArray("TAnaJet", 1000);
  fnTrackJets      = 0;

  fCandidates      = new TClonesArray("TAnaCand", 1000);
  fnCandidates     = 0;

  fPV              = new TClonesArray("TAnaVertex", 100);
  fnPV             = 0;

}

// ----------------------------------------------------------------------
void TAna01Event::Clear(Option_t *option) {

  fLumi = fLumiSection = fBx = fOrbit = -9999; 
  fTimeLo = fTimeHi = 0; 
  
  TAnaTrack *pTrack;
  for (int i = 0; i < fnRecTracks; i++) {
    pTrack = getRecTrack(i);
    pTrack->clear();
  }
  fRecTracks->Clear(option);
  fnRecTracks = 0;

  TAnaMuon *pMuon;
  for (int i = 0; i < fnMuons; i++) {
    pMuon = getMuon(i);
    pMuon->clear();
  }
  fMuons->Clear(option);
  fnMuons = 0;

  TTrgObj *pTrgObj;
  for (int i = 0; i < fnTrgObj; i++) {
    pTrgObj = getTrgObj(i);
    pTrgObj->clear();
  }
  fTrgObj->Clear(option);
  fnTrgObj = 0;

  for (int i = 0; i < NL1T; ++i) {
    for (int j = 0; j < 32; ++j) {
      fL1TNames[i*32+j] = ""; 
    }
    
    fL1TWords[i]  = 0; 
    fL1TWasRun[i] = 0; 
  }

  for (int i = 0; i < NL1TT; ++i) {
    for (int j = 0; j < 32; ++j) {
      fL1TTNames[i*32+j] = ""; 
    }
    
    fL1TTWords[i]  = 0; 
    fL1TTWasRun[i] = 0; 
  }

  for (int i = 0; i < NHLT; ++i) {
    for (int j = 0; j < 32; ++j) {
      fHLTNames[i*32+j] = ""; 
    }
    fHLTWords[i] = 0; 
    fHLTWords[i] = 0; 
  }

  TAnaTrack *pSigTrack;
  for (int i = 0; i < fnSigTracks; i++) {
    pSigTrack = getSigTrack(i);
    pSigTrack->clear();
  }
  fSigTracks->Clear(option);
  fnSigTracks = 0;

  TAnaJet *pCaloJet;
  for (int i = 0; i < fnCaloJets; i++) {
    pCaloJet = getCaloJet(i);
    pCaloJet->clear();
  }
  fCaloJets->Clear(option);
  fnCaloJets = 0;

  TAnaJet *pGenJet;
  for (int i = 0; i < fnGenJets; i++) {
    pGenJet = getGenJet(i);
    pGenJet->clear();
  }
  fGenJets->Clear(option);
  fnGenJets = 0;

  TAnaJet *pTrackJet;
  for (int i = 0; i < fnTrackJets; i++) {
    pTrackJet = getTrackJet(i);
    pTrackJet->clear();
  }
  fTrackJets->Clear(option);
  fnTrackJets = 0;

  TAnaCand *pCand;
  for (int i = 0; i < fnCandidates; i++) {
    pCand = getCand(i);
    pCand->clear();
  }
  fCandidates->Clear(option);
  fnCandidates = 0;

  TGenCand *pGenCand;
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    pGenCand->clear();
  }
  fGenCands->Clear(option);
  fnGenCands = 0;

  TAnaVertex *pPV;
  for (int i = 0; i < fnPV; i++) {
    pPV = getPV(i);
    pPV->clear();
  }
  fPV->Clear(option);
  fnPV = 0;

}


// ----------------------------------------------------------------------
TGenCand* TAna01Event::getGenCand(Int_t n) { 
  return (TGenCand*)fGenCands->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TGenCand* TAna01Event::addGenCand() {
  TClonesArray& d = *fGenCands; 
  new(d[d.GetLast()+1]) TGenCand();
  ++fnGenCands;
  return (TGenCand*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
void TAna01Event::dumpGenBlock() {
  TGenCand *pGenCand;
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    pGenCand->dump();
  }
}

// ----------------------------------------------------------------------
int TAna01Event::getGenIndex(double px, double py, double pz, int charge, double precision) {
  TGenCand *pGenCand;
  int index(-1);

  //   cout << "Searching among " << fnGenCands << " cands for match to p = (" 
  //        << px << ", " << py << ", " << pz << ") with ID = " << id << endl;

  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    if (pGenCand->fQ * charge < 0) continue;
    if ((TMath::Abs((pGenCand->fP.X() - px))/px) > precision) continue;
    if ((TMath::Abs((pGenCand->fP.Y() - py))/py) > precision) continue;
    if ((TMath::Abs((pGenCand->fP.Z() - pz))/pz) > precision) continue;
    
    index = i; 
    //     cout << "I think this is a match at index= " << index 
    // 	 << "  " << px << " " << py << " " << pz << " " << endl;
    //     pGenCand->dump();
    break;
  }
  return index;
}


// ----------------------------------------------------------------------
int TAna01Event::getGenIndexWithDeltaR(double pt, double eta, double phi, double charge) {
  TGenCand *pGenCand;
  int index(-1);
  int tmp_index[30];
  TLorentzVector track, gencand;
  int count(0);
  track.SetPtEtaPhiM(pt, eta, phi, 0.);
				       
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    if (pGenCand->fStatus == 1) {
      if (charge * pGenCand->fQ > 0) {
	/*   
	     gencand.SetXYZM(pGenCand->fP.X(),
	     pGenCand->fP.Y(),
	     pGenCand->fP.Z(),
	     pGenCand->fMass);
	*/
	double dR = track.DeltaR(pGenCand->fP);
	double dpt = TMath::Abs((track.Perp() - pGenCand->fP.Perp())/track.Perp());
	
	if (dR < 0.12 && dpt < 0.3  && count < 30) {
	  tmp_index[count] = i;
	  count++; 
	}
      }
    }
    
  }
  
  double deltaR_min(0.12);
  int  i_min(-1); 	
  if (count == 1) index = tmp_index[0];	
  if (count > 1) { 
    for (int i = 0; i < count; i++) {
      pGenCand = getGenCand(tmp_index[i]);;
      /*	       
		      gencand.SetXYZM(pGenCand->fP.X(),
		      pGenCand->fP.Y(),
		      pGenCand->fP.Z(),
		      pGenCand->fMass);
      */
      if (track.DeltaR(pGenCand->fP) < deltaR_min) {
	deltaR_min = track.DeltaR(pGenCand->fP);			
	i_min = i;
      }
    }
    
    index = i_min;
    
  }
  
  return index;
}


// ----------------------------------------------------------------------
TAnaTrack* TAna01Event::getRecTrack(Int_t n) { 
  return (TAnaTrack*)fRecTracks->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaTrack* TAna01Event::addRecTrack() {
  TClonesArray& d = *fRecTracks; 
  new(d[d.GetLast()+1]) TAnaTrack(fnRecTracks);
  ++fnRecTracks;
  return (TAnaTrack*)d[d.GetLast()];
}


// ----------------------------------------------------------------------
TAnaMuon* TAna01Event::getMuon(Int_t n) { 
  return (TAnaMuon*)fMuons->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaMuon* TAna01Event::addMuon() {
  TClonesArray& d = *fMuons; 
  new(d[d.GetLast()+1]) TAnaMuon(fnMuons);
  ++fnMuons;
  return (TAnaMuon*)d[d.GetLast()];
}


// ----------------------------------------------------------------------
TTrgObj* TAna01Event::getTrgObj(Int_t n) { 
  return (TTrgObj*)fTrgObj->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TTrgObj* TAna01Event::addTrgObj() {
  TClonesArray& d = *fTrgObj; 
  new(d[d.GetLast()+1]) TTrgObj(fnTrgObj);
  ++fnTrgObj;
  return (TTrgObj*)d[d.GetLast()];
}


// ----------------------------------------------------------------------
TAnaTrack* TAna01Event::getSigTrack(Int_t n) { 
  return (TAnaTrack*)fSigTracks->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaTrack* TAna01Event::addSigTrack() {
  TClonesArray& d = *fSigTracks; 
  new(d[d.GetLast()+1]) TAnaTrack(fnSigTracks);
  ++fnSigTracks;
  return (TAnaTrack*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::getCaloJet(Int_t n) { 
  return (TAnaJet*)fCaloJets->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::addCaloJet() {
  TClonesArray& d = *fCaloJets; 
  new(d[d.GetLast()+1]) TAnaJet(fnCaloJets);
  ++fnCaloJets;
  return (TAnaJet*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::getTrackJet(Int_t n) {
 return (TAnaJet*)fTrackJets->UncheckedAt(n);
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::addTrackJet() {
 TClonesArray& d = *fTrackJets;
 new(d[d.GetLast()+1]) TAnaJet(fnTrackJets);
 ++fnTrackJets;
 return (TAnaJet*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::getGenJet(Int_t n) { 
  return (TAnaJet*)fGenJets->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::addGenJet() {
  TClonesArray& d = *fGenJets; 
  new(d[d.GetLast()+1]) TAnaJet(fnGenJets);
  ++fnGenJets;
  return (TAnaJet*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaCand* TAna01Event::getCand(Int_t n) { 
  return (TAnaCand*)fCandidates->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaCand* TAna01Event::addCand() {
  TClonesArray& d = *fCandidates; 
  new(d[d.GetLast()+1]) TAnaCand();
  ++fnCandidates;
  return (TAnaCand*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaVertex* TAna01Event::getPV(Int_t n) { 
  return (TAnaVertex*)fPV->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaVertex* TAna01Event::addPV() {
  TClonesArray& d = *fPV; 
  new(d[d.GetLast()+1]) TAnaVertex();
  ++fnPV;
  return (TAnaVertex*)d[d.GetLast()];
}



// ----------------------------------------------------------------------
void TAna01Event::dump() {
  cout << "Run: " << fRunNumber << ", Event: " << fEventNumber 
       << " nGenCand = " << nGenCands()
       << " nRecTrk = " << nRecTracks()
       << endl;

  cout << "Primary Vtx: "; 
  TAnaVertex *pPV;
  for (int i = 0; i < fnPV; i++) {
    pPV = getPV(i);
    if (i > 0) cout << "             ";
    pPV->dump();
  }
  
  for (int j = 0; j < 32; ++j) {
    for (int i = 0; i < NL1T; ++i) {
      cout << Form("%15s %2d ", fL1TNames[32*i+j].Data(), (fL1TWords[i] & (0x1<<j)? 1: 0)); 
    }
    cout << endl;

    for (int i = 0; i < NL1TT; ++i) {
      cout << Form("%15s %2d ", fL1TTNames[32*i+j].Data(), (fL1TTWords[i] & (0x1<<j)? 1: 0)); 
    }
    cout << endl;
  }    

  for (int j = 0; j < 32; ++j) {
    for (int i = 0; i < NHLT; ++i) {
      cout << Form("%15s %2d ", fHLTNames[32*i+j].Data(), (fHLTWords[i] & (0x1<<j)? 1: 0)); 
    }
    cout << endl;
  }    

  TGenCand *pGenCand;
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    cout << "genT: ";
    pGenCand->dump();
  }
  
  TAnaTrack *pTrk;
  TAnaCand  *pCand;
  for (int i = 0; i < nRecTracks(); ++i) {
    pTrk = getRecTrack(i);
    cout << "recT: ";
    pTrk->dump();
  }

  TAnaMuon  *pMuon;
  for (int i = 0; i < nMuons(); ++i) {
    pMuon = getMuon(i);
    cout << "muon: ";
    pMuon->dump();
  }
  
  TTrgObj  *pTrgObj;
  for (int i = 0; i < nTrgObj(); ++i) {
    pTrgObj = getTrgObj(i);
    cout << "trgObj: ";
    pTrgObj->dump();
  }
  
  for (int i = 0; i < nSigTracks(); ++i) {
    pTrk = getSigTrack(i);
    cout << "sigT: ";
    pTrk->dump();
  }

  for (int i = 0; i < nCands(); ++i) {
    pCand = getCand(i);
    cout << "cand: ";
    pCand->dump();
  }

  cout << endl;
}


