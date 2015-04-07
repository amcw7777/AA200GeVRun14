#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoD0AnaMaker.h"
#include "StPicoHFMaker/StHFCuts.h"
//////Refit include lib
#include "PhysicalConstants.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorD.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "StEvent/StDcaGeometry.h"
#include "KFParticle.h"
#include "TVector3.h"
#include "TRVector.h"
#include "TLorentzVector.h"
#include "MVertex.h"
#include "MTrack.h"
#include "StPrimaryVertex.h"

#include "StiMaker/StAnneling.h"
#include "StiMaker/StKFVertex.h"
//
#include "TSpectrum.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"
#include "TCanvas.h"
#include "TH1K.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "Sti/StiHit.h"


#include "StiMaker/StKFEvent.h"
#include "StiMaker/StKFTrack.h"
#include "StiMaker/StVertexP.h"
#include "StiMaker/StVertexT.h"
#include "TDirectory.h"
#include "Stypes.h"
#include "SystemOfUnits.h"
#include "StiMaker/StKFVertexMaker.h"
#include "StDetectorDbMaker/St_vertexSeedC.h"
#include "Sti/StiKalmanTrack.h"
#include "Sti/StiKalmanTrackNode.h"
#include <StEvent.h>
//#include "StiMaker/StiStEventFiller.h"
#include "TRMatrix.h"
#include "TRSymMatrix.h"

#include "Sti/StiToolkit.h"
#include "TArrayI.h"
#include "StiMaker/StiDefaultToolkit.h"
#include <vector>
//
#include <stdio.h>
#include <time.h>

ClassImp(StPicoD0AnaMaker)
StKFVerticesCollection *StPicoD0AnaMaker::fcVertices = 0;

StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name,char const * inputFilesList, 
    char const * outName,StPicoDstMaker* picoDstMaker): 
  StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL), mOutFileName(outName), mInputFileList(inputFilesList),
  mOutputFile(NULL), mChain(NULL), mEventCounter(0), mHFCuts(NULL)
{}

Int_t StPicoD0AnaMaker::Init()
{
   mPicoD0Event = new StPicoD0Event();

   mChain = new TChain("T");
   std::ifstream listOfFiles(mInputFileList.Data());
   if (listOfFiles.is_open())
   {
      std::string file;
      while (getline(listOfFiles, file))
      {
         LOG_INFO << "StPicoD0AnaMaker - Adding :" << file << endm;
         mChain->Add(file.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoD0AnaMaker - Could not open list of files. ABORT!" << endm;
      return kStErr;
   }

   mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
   mChain->SetBranchAddress("dEvent", &mPicoD0Event);

   mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
   mRefittuple = new TNtuple("mRefittuple","mRefittuple","refitx:refity:refitz:prmx:prmy:prmz:runId:eventId:ranking:time:mult");
   timemult = new TH2F("timemult","",100,0,500,100,0,5);
   mOutputFile->cd();

   if (!mHFCuts)
    mHFCuts = new StHFCuts;   

   // -------------- USER VARIABLES -------------------------
  //dcaG = new StDcaGeometry();
  dca = new StDcaGeometry();
  primV  = new StPrimaryVertex;
  fNzBins = 2500;
  fNPasses = 2;
  fSpectrum = 0;
  fzWindow = 2;
  fVtxM = 0;
  fVerticesPass = 0;
  fTempLog = 2;
  fminBrent = 0;
  func = 0;
  mBeamLine = kFALSE;
  fc1 = 0;

  //StKFVertexMaker constr.
  Int_t npeaks = 100;
  Double_t zmin = -250;
  Double_t zmax = 250;
  //  StKFVertex::_debug = 1;
  for (Int_t pass = 0; pass < fNPasses; pass++) {
      fVtxs[pass] = new TH1F(Form("Vtx%1i",pass),Form("z-dca distribution for pass = %1i",pass),fNzBins,zmin,zmax);
      fVtxs[pass]->SetDirectory(0);
      if (pass)  fVtxs[pass]->SetLineColor(5);
      fVtxs[pass]->SetDefaultSumw2();
      fVtxs[pass]->SetStats(0);
      fVtxKs[pass] = new TH1K(Form("VtxK%1i",pass),Form("z-dca distribution for pass = %1i",pass),fNzBins,zmin,zmax);
      fVtxKs[pass]->SetDirectory(0);
      fVtxKs[pass]->SetStats(0);
      fVtxKs[pass]->SetLineColor(2);
  }
  fVtxM = new TH1F("VtxM","MuDst reconstructed multiplicities versus Z",fNzBins,zmin,zmax);
  fVtxM->SetDirectory(0);
  fSpectrum = new TSpectrum(2*npeaks);
  func = new ROOT::Math::Functor1D(&StPicoD0AnaMaker::AnnelingFcn);
  fminBrent = new ROOT::Math::GSLMinimizer1D();
  fVerticesPass = new StKFVerticesCollection *[fNPasses+1];
  memset (fVerticesPass, 0, (fNPasses+1)*sizeof(StKFVerticesCollection *));
  fParticles = new TObjArray();
  fParticles->SetOwner(kTRUE);


   return kStOK;
}
//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
   /*  */
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
   LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
   mOutputFile->cd();
   // save user variables here
   mRefittuple->Write();
  // timemult-Write();
   mOutputFile->Close();
/////////Vertex Fitting
  SafeDelete(fVtxM);
  for (Int_t pass = 0; pass < fNPasses; pass++) {
    SafeDelete(fVtxs[pass]);
    SafeDelete(fVtxKs[pass]);
  }
  SafeDelete(fSpectrum);
  SafeDelete(func);
  SafeDelete(fminBrent);
  delete [] fVerticesPass; fVerticesPass = 0;
  SafeDelete(fParticles);
  SafeDelete(primV);
  SafeDelete(dca);
	
   return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Make()
{
   readNextEvent();
   time_t t_1,t_2;
   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   //StPicoDst const* picoDst = mPicoDstMaker->picoDst();
   picoDst = mPicoDstMaker->picoDst();

   if (!picoDst)
   {
      LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }
   t_1 = clock();
   if(mPicoD0Event->runId() != picoDst->event()->runId() ||
       mPicoD0Event->eventId() != picoDst->event()->eventId())
   {
     LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
     LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
     exit(1);
   }

   // -------------- USER ANALYSIS -------------------------
   TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();
   float refittuple_fill[20]; 
   
   StThreeVectorF pVtx(-999.,-999.,-999.);
   StThreeVectorF testVertex(-999.,-999.,-999.);
   StPicoEvent *event = (StPicoEvent *)picoDst->event();
   if(!(event->isMinBias()))
   {
     LOG_WARN << " Not Min Bias! Skip! " << endm;
     return kStWarn;
   }
   float mult = event->refMult();
   if(event) {
     pVtx = event->primaryVertex();
   }
   testVertex = pVtx;
   vector<int> daughter;
   daughter.clear();
   for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
   {
      // this is an example of how to get the kaonPion pairs and their corresponsing tracks
      StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
      if(!isGoodPair(kp)) continue;

      //StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
     // StPicoTrack const* pion = picoDst->track(kp->pionIdx());
      daughter.push_back(kp->kaonIdx());
      daughter.push_back(kp->pionIdx());

   }
   double mranking=0;
   mranking = primaryVertexRefit(&testVertex,daughter);
   cout<<"TIME for event with mult = "<<mult<<endl;
   refittuple_fill[0] = testVertex.x(); 
   refittuple_fill[1] = testVertex.y(); 
   refittuple_fill[2] = testVertex.z(); 
   refittuple_fill[3] = pVtx.x(); 
   refittuple_fill[4] = pVtx.y(); 
   refittuple_fill[5] = pVtx.z(); 
   refittuple_fill[6] = (double)mPicoD0Event->runId(); 
   refittuple_fill[7] = (double)mPicoD0Event->eventId(); 
   refittuple_fill[8] = mranking; 
   t_2 = clock();
   float dtime = 0.000001*difftime(t_2,t_1);
   refittuple_fill[9] = dtime; 
   refittuple_fill[10] = mult; 
   mRefittuple->Fill(refittuple_fill);
   cout<<"TIME to run is "<<dtime<<endl;
   timemult->Fill(mult,dtime);

   return kStOK;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodPair(StKaonPion const* const kp) const
{
  if(!kp) return false;

  StPicoTrack const* kaon = mPicoDstMaker->picoDst()->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDstMaker->picoDst()->track(kp->pionIdx());

  //  To be replaced by mHFCuts->isGoodSecondaryVertexPair(kp))
  bool pairCuts = kp->m() > mHFCuts->cutSecondaryPairMassMin() && 
    kp->m() < mHFCuts->cutSecondaryPairMassMax() &&
    std::cos(kp->pointingAngle()) > mHFCuts->cutSecondaryPairCosThetaMin() &&
    kp->decayLength()  > mHFCuts->cutSecondaryPairDecayLengthMin() && 
    kp->decayLength()  < mHFCuts->cutSecondaryPairDecayLengthMax() &&
    kp->dcaDaughters() < mHFCuts->cutSecondaryPairDcaDaughtersMax();

  return (mHFCuts->isGoodTrack(kaon) && mHFCuts->isGoodTrack(pion) &&
	  mHFCuts->isTPCKaon(kaon) && mHFCuts->isTPCPion(pion) && 
	  pairCuts);
}
/*
*/
double StPicoD0AnaMaker::primaryVertexRefit(StThreeVectorF *mRefitVertex, vector<int>& daughter) {
/*
  //StKFVertexMaker variables
  fNzBins = 2500;
  fNPasses = 2;
  fSpectrum = 0;
  fzWindow = 2;
  fVtxM = 0;
  fVerticesPass = 0;
  fTempLog = 2;
  fminBrent = 0;
  func = 0;
  mBeamLine = kFALSE;
  fc1 = 0;

  //StKFVertexMaker constr.
  Int_t npeaks = 100;
  Double_t zmin = -250;
  Double_t zmax = 250;
  //  StKFVertex::_debug = 1;
  for (Int_t pass = 0; pass < fNPasses; pass++) {
      fVtxs[pass] = new TH1F(Form("Vtx%1i",pass),Form("z-dca distribution for pass = %1i",pass),fNzBins,zmin,zmax);
      fVtxs[pass]->SetDirectory(0);
      if (pass)  fVtxs[pass]->SetLineColor(5);
      fVtxs[pass]->SetDefaultSumw2();
      fVtxs[pass]->SetStats(0);
      fVtxKs[pass] = new TH1K(Form("VtxK%1i",pass),Form("z-dca distribution for pass = %1i",pass),fNzBins,zmin,zmax);
      fVtxKs[pass]->SetDirectory(0);
      fVtxKs[pass]->SetStats(0);
      fVtxKs[pass]->SetLineColor(2);
  }
  fVtxM = new TH1F("VtxM","MuDst reconstructed multiplicities versus Z",fNzBins,zmin,zmax);
  fVtxM->SetDirectory(0);
  fSpectrum = new TSpectrum(2*npeaks);
  func = new ROOT::Math::Functor1D(&StPicoD0AnaMaker::AnnelingFcn);
  fminBrent = new ROOT::Math::GSLMinimizer1D();
  fVerticesPass = new StKFVerticesCollection *[fNPasses+1];
  memset (fVerticesPass, 0, (fNPasses+1)*sizeof(StKFVerticesCollection *));
  fParticles = new TObjArray();
  fParticles->SetOwner(kTRUE);

*/
  //StKFVertexMaker Make part
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  Double_t bField = event->bField();
  KFParticle::SetField(bField);
  Int_t NGoodGlobals = 0;
  int nTracks = picoDst->numberOfTracks();
  for (int i=0; i < nTracks; i++) {
     StPicoTrack *gTrack = (StPicoTrack*)picoDst->track(i);
     if (! gTrack) continue;
//     const StDcaGeometry* dca = gTrack->dcaGeometry();
//     StDcaGeometry *dca = new StDcaGeometry();
     dca->set(gTrack->params(),gTrack->errMatrix());
     if (! dca) continue;
//     if (gTrack->flag()     <   0) continue;     // Bad fit
//     if (gTrack->flag()     > 700) continue;     // FTPC
//     if (gTrack->flag()%100 == 11) continue;     // Short track pointing to EEMC
//     if ((gTrack->isWestTpcOnly() || gTrack->isEastTpcOnly()) && gTrack->isPostXTrack()) continue; // wrong TPC side track
     //Int_t kg = gTrack->key();
     Int_t kg = gTrack->id();
//     TrackNodeMap[kg] = node;
     bool flagDdaughterCand = false;
     for (vector<int>::size_type j=0; j < daughter.size(); j++) {
	if (daughter[j] == kg) 
        {
          flagDdaughterCand = 1;
        }
     }
     if (flagDdaughterCand == 1) continue;


     KFParticle *particle = AddTrackAt(dca,kg);
     NGoodGlobals++;
  }
  if (NGoodGlobals < 2 ) return false;
  Fit();
  if (! Vertices()) return false;
  //
  //  In case there are no tracks left we better quit
  //
//  StSPtrVecTrackDetectorInfo& detInfoVec = pEvent->trackDetectorInfo();
  Int_t Nvtx = Vertices()->NoVertices();
  //cout << "NGoodGlobals: " <<  NGoodGlobals << "   # Vertices: " << Nvtx << endl;
  double mrank=0;
  for (Int_t l = 0; l < Nvtx; l++) {
      const StKFVertex *V = Vertices()->Vertex(l);
      if (! V) continue;
      //if (Debug() > 2) 
//      V->PrintW();
      // Store vertex
//      StPrimaryVertex *primV  = new StPrimaryVertex;
      StThreeVectorF XVertex(&V->Vertex().X());
      primV->setPosition(XVertex);
      primV->setChiSquared(V->Vertex().Chi2()/V->Vertex().GetNDF());  
      primV->setProbChiSquared(TMath::Prob(V->Vertex().GetChi2(),V->Vertex().GetNDF()));
      Float_t cov[6];
      TCL::ucopy(&V->Vertex().Covariance(0),cov,6);
      primV->setCovariantMatrix(cov); 
      primV->setVertexFinderId(KFVertexFinder);
      primV->setFlag(1); // was not set earlier by this vertex finder ?? Jan
      primV->setRanking(333);
      primV->setNumTracksUsedInFinder(V->NoTracks());
      Bool_t beam = kFALSE;
      Double_t Pars[6];
      TCL::ucopy(&V->Vertex().X(), Pars, 6);
      Double_t Cov[21];
      TCL::ucopy(&V->Vertex().Covariance(0), Cov, 21);
      if (!StiToolkit::instance()) new StiDefaultToolkit;

///////?????????????????????
      StiHit *Vertex = StiToolkit::instance()->getHitFactory()->getInstance();
      Vertex->setGlobal(0, 0, V->Vertex().X(), V->Vertex().Y(), V->Vertex().Z(), 0);
      Vertex->setError(cov);
      Int_t NoTracks = V->NoTracks();
      //cout << "Vertex #: " <<  l << "   # tracks: " << NoTracks << endl;
      TArrayI indexT(NoTracks); Int_t *indexes = indexT.GetArray();
      TArrayI IdT(NoTracks);    Int_t *Ids     = IdT.GetArray();
      for (Int_t itk = 0; itk < NoTracks; itk++) {
	  Ids[itk] = 999999;
	  const StKFTrack*   track = V->Track(itk);
	  if (! track) continue;
	  const KFParticle   &P = track->Particle();
	  Int_t kg = P.GetID()%100000;
	  Ids[itk] = kg;
      }
      TMath::Sort(NoTracks,Ids,indexes,0);
      for (Int_t i = 0; i < NoTracks; i++) {
	  Int_t itk = indexes[i];
	  const StKFTrack*   track = V->Track(itk);
	  if (! track) continue;
	  const KFParticle   &P = track->Particle();
	  Int_t kg = P.GetID()%100000;
	  if (kg == 0) {
	     assert(!beam);
	     beam = kTRUE;
	     continue;
	  }
	  if (Debug() > 2) {
	     const KFParticle   *PO = track->OrigParticle();
	     const KFParticle *PS[2] = {PO, &P};
	     for (Int_t m = 0; m < 2; m++) {
		if (! m) cout << "Original";
		else     cout << "Fitted  ";
		static const Char_t *names[6] = {"x","y","z","px","py","pz"};
		for (Int_t j = 0; j < 6; j++) {
		   cout << Form(" %2s: %8.3f +/- %8.3f",names[j], PS[m]->GetParameter(j), PS[m]->GetCovariance(j,j) > 0 ? TMath::Sqrt(PS[m]->GetCovariance(j,j)) : -13);
		}
	//	cout << endl;
	     }
	  }
	 // node = TrackNodeMap[kg];
        }
//????????
      if (beam ) primV->setBeamConstrained();
      primV->setTrackNumbers();
      calculateRank(primV);
//      cout<<"checkpoint1 = "<<primV->ranking()<<endl;
//      cout<<"checkpoint2 = "<<primV->position()<<endl;
      if(primV->ranking()>mrank)
      {
        mrank = primV->ranking();
        *mRefitVertex = primV->position();
      }
		
   
  } 
  Clear();

  return mrank;

}


//StKFVertexMaker functions_______________________________________________________
void StPicoD0AnaMaker::calculateRank(StPrimaryVertex *primV) {    
  // Calculation of veretx ranks to select 'best' (i.e. triggered)  vertex
  // Simpilfied version (w/o weighting)
  Float_t rank = primV->probChiSquared();
  static Float_t Wveto = 1;
  static Float_t Wmatch = 4;
  if (primV->isBeamConstrained()) rank += Wmatch;
  rank -= Wveto*primV->numPostXTracks();
  rank += Wmatch*primV->numTracksWithPromptHit();
  rank += Wmatch*primV->numTracksCrossingCentralMembrane();
  rank += Wmatch*primV->numMatchesWithCTB()
    -     Wveto*primV->numNotMatchesWithCTB();
  rank += Wmatch*primV->numMatchesWithBTOF() 
    -     Wveto*primV->numNotMatchesWithBTOF();
  rank += Wmatch*(primV->numMatchesWithBEMC() + primV->numMatchesWithEEMC());
    -     Wveto*(primV->numNotMatchesWithBEMC() + primV->numNotMatchesWithEEMC());
  if (primV->numTracksTpcWestOnly() > 0 && primV->numTracksTpcEastOnly() > 0) 
    rank += Wmatch*TMath::Min(primV->numTracksTpcWestOnly(),primV->numTracksTpcEastOnly());
  rank += 100.0 + primV->numTracksUsedInFinder();
  primV->setRanking(rank); 
  //if (Debug()) primV->Print();
  //cout << "rank = " << rank << endl;
}

//________________________________________________________________________________
KFParticle *StPicoD0AnaMaker::AddTrackAt(const StDcaGeometry *dca, Int_t kg) {
  fParticles->AddAtAndExpand (0,kg);
  if (! dca) return 0;
  Double_t xyzp[6], CovXyzp[21];
  dca->GetXYZ(xyzp,CovXyzp);
  static MTrack track;
  track.SetParameters(xyzp);
  track.SetCovarianceMatrix(CovXyzp);
  track.SetNDF(1);
  //    track.SetChi2(GlobalTracks_mChiSqXY[k]);
  track.SetID(kg);
  Int_t q   = 1;
  Int_t pdg = 211;
  if (dca->charge() < 0) {
    q = -1;
    pdg = -211;
  } 
  track.SetCharge(q);
  KFParticle *particle = new KFParticle(track, pdg);
  particle->SetID(kg);
  fParticles->AddAt(particle,kg);
  return particle;
}
//________________________________________________________________________________
void StPicoD0AnaMaker::Fit() {
  if (Debug() != 2)  StKFVertex::SetDebug(Debug());
  fcVertices = 0;
  for (Int_t i = 0; i < fNPasses+1; i++) {
    SafeDelete(fVerticesPass[i]);
  }
  Int_t NGoodGlobals = Particles().GetLast();
  
  Double_t TempLog = fTempLog; // default Temperature Log
  for (Int_t pass = 0; pass < fNPasses; pass++) {
    Int_t nAccepted = 0;
    Double_t dZ = fVtxs[pass]->GetBinWidth(1);
    for (Int_t k = 0; k < NGoodGlobals; k++) {
      KFParticle *particle = (KFParticle *) Particles()[k];
      if (! particle) continue;
      Double_t pT;
      Double_t dpT;
      particle->GetPt(pT,dpT);
      Double_t offset = 0.5*particle->GetPz()/pT;
      Double_t SigmaZ = TMath::Sqrt(particle->Covariance(2,2) + offset*offset);
      SigmaZ += dZ;
      Double_t Z = particle->GetZ();
      fVtxKs[pass]->Fill(Z);
      Int_t bin1 = fVtxs[pass]->FindBin(Z - 5*SigmaZ);
      if (bin1 < 1) bin1 = 1;
      Int_t bin2 = fVtxs[pass]->FindBin(Z + 5*SigmaZ);
      if (bin2 > fNzBins) bin2 = fNzBins;
      Double_t z = fVtxs[pass]->GetBinCenter(bin1);
      for (Int_t bin = bin1; bin <= bin2; bin++, z += dZ) {
	fVtxs[pass]->Fill(z,(TMath::Erfc((z - Z - fzWindow)/SigmaZ) - TMath::Erfc((z - Z + fzWindow)/SigmaZ))/2.);
      }
      nAccepted++;
    }
    Double_t F = fVtxKs[pass]->GetEntries();
    if (F < 1) continue;
    fVtxKs[pass]->SetNormFactor(F/dZ);
    fVtx = fVtxs[0]; // << switch between types    Vtx = fVtxKs[0];
    TString opt("new");
    if (! Canvas()) opt = "goff";
    Int_t nfound = fSpectrum->Search(fVtx,3,opt,TMath::Min(0.1,5./NGoodGlobals));
    if (! nfound) continue;
    if (Canvas()) {
      Canvas()->cd();
      fVtxs[0]->Draw(); fVtxKs[0]->Draw("same");
      fVtxM->Draw("same");
      if (pass)    fVtx->Draw("same");
      Canvas()->Update();
    }
    if (StKFVertex::Debug() > 1) {
      LOG_INFO << "Found " << nfound 
	   << " candidate peaks to fit with " << NGoodGlobals
	   << " good globals from with " <<  nAccepted  << " accepted" << endm;
    }
    Double_t *zOfPeaks = new Double_t[nfound];
    Int_t npeaks = 0;
#if ROOT_VERSION_CODE > 336641 // ROOT_VERSION(5,35,1) 
    Double_t *xpeaks = fSpectrum->GetPositionX();
#else
    Float_t *xpeaks = fSpectrum->GetPositionX();
#endif
    for (Int_t p = 0; p < nfound; p++) {
#if ROOT_VERSION_CODE > 336641 // ROOT_VERSION(5,35,1) 
      Double_t xp = xpeaks[p];
#else
      Float_t xp = xpeaks[p];
#endif
      Int_t bin = fVtx->GetXaxis()->FindBin(xp);
      Double_t yp = fVtx->GetBinContent(bin);
      Double_t ep = fVtx->GetBinError(bin);
      if (yp-1.25*ep < 0) continue;
      zOfPeaks[npeaks] = xp;
      npeaks++;
    }
    if (StKFVertex::Debug() > 1) {
      LOG_INFO << "Found " << npeaks << " useful peaks to fit" << endm;
    }
    if (! npeaks) {delete [] zOfPeaks; break; }
    if (fVerticesPass[pass]) {delete fVerticesPass[pass]; fVerticesPass[pass] = 0;}
    fVerticesPass[pass] = new StKFVerticesCollection(npeaks, zOfPeaks);
    delete [] zOfPeaks;
    fcVertices = fVerticesPass[pass];
    fcVertices->DoTrack2VertexAssociation(Particles());
    if (! fcVertices->NoVertices())                         continue;
    if (AnnelingFcn(TMath::Exp(-TempLog)) <= 0) continue;
    if (! fcVertices->NoVertices())                         continue;
    fcVertices->UniqueTracks2VertexAssociation(); // Make track associated with only vertex
    //       fcVertices->PrintV(NoMuMcVertex,NoMuMcTrack,StMuMcVertex_time,
    // 		       StMuMcVertex_xyzV_mX1,StMuMcVertex_xyzV_mX2,StMuMcVertex_xyzV_mX3,
    // 		       StMuMcVertex_NoDaughters,StMuMcVertex_IdParTrk,StMuMcTrack_gePid);
  }
  if (! fVerticesPass[0]) return;
  if (fNPasses > 1 && Canvas()) {
    Canvas()->cd();
    fVtxs[1]->Draw("same"); 
    Canvas()->Update();
  }
  Int_t N1 = fVerticesPass[0]->NoVertices();
  if (! N1) return;
  if (fVerticesPass[1]) {
    *fVerticesPass[0] += *fVerticesPass[1];
  }
  fcVertices = fVerticesPass[0];
  fcVertices->MergeDuplicatedVertices();
  if (! fcVertices->NoVertices()) return;
  // Double_t Temperature = TMath::Exp(TempLog);
  TempLog = 5;
  Double_t Temperature = TMath::Exp(TempLog);
#if 1  
  // secondary vertices
  Int_t pass = fNPasses;
  if (fVerticesPass[pass]) {delete fVerticesPass[pass]; fVerticesPass[pass] = 0;}
  fVerticesPass[pass] = new StKFVerticesCollection();
  fcVertices = fVerticesPass[pass];
  StAnneling::SetTemperature(Temperature);
  for (Int_t k = 1; k < NGoodGlobals; k++) {
    KFParticle *particleK = (KFParticle *) Particles()[k];
    if (! particleK) continue;
    if (particleK->GetID() > 100000) continue;
    StKFVertex *vtx = 0;
    for (Int_t l = k+1; l < NGoodGlobals; l++) {
      KFParticle *particleL = (KFParticle *) Particles()[l];
      if (! particleL) continue;
      if (particleL->GetID() > 100000) continue;
      Double_t dist = particleK->GetDistanceFromParticle(*particleL);
      if (dist > 5.0) continue;
      if (! vtx) {
	vtx = new StKFVertex(fcVertices->NoVertices() + 1);
	vtx->AddTrack(new StKFTrack(k,particleK));
      }
      vtx->AddTrack(new StKFTrack(l,particleL));
    }
    if (! vtx) continue;
    vtx->Fit();
    Int_t N = vtx->NoTracks();
    if (! N) {delete vtx; vtx = 0; continue;}
    Double_t X = vtx->Vertex().X();
    Double_t Y = vtx->Vertex().Y();
    Double_t R = TMath::Sqrt(X*X + Y*Y);
    if (R > 200 ) {delete vtx; vtx = 0; continue;}
    Double_t prob = TMath::Prob(vtx->Vertex().GetChi2(),vtx->Vertex().GetNDF());
    if (N > 2 || prob > 1.e-3) {// Allow V2 to share tracks
      for (Int_t i = 0; i < N; i++) {
	KFParticle *particle = vtx->Track(i)->OrigParticle();;
	Int_t ID = particle->GetID()%100000 + 100000*vtx->ID();;
	particle->SetID(ID);
      }
    }
    fcVertices->AddVertex(vtx);
  }
  if (StKFVertex::Debug() > 1) {
    LOG_INFO << "Candidate for secondary vertices: " << fcVertices->NoVertices() << endm;
  }
  if ( fcVertices->NoVertices() ) {
    //       fcVertices->PrintV(NoMuMcVertex,NoMuMcTrack,StMuMcVertex_time,
    // 		       StMuMcVertex_xyzV_mX1,StMuMcVertex_xyzV_mX2,StMuMcVertex_xyzV_mX3,
    // 		       StMuMcVertex_NoDaughters,StMuMcVertex_IdParTrk,StMuMcTrack_gePid);
    *fVerticesPass[0] += *fVerticesPass[fNPasses];
  }
  // end of loop for secondary vertices
#endif
  fcVertices = fVerticesPass[0];
  fcVertices->Compress();
  if (! fcVertices->NoVertices()) return;
  fcVertices->MergeDuplicatedVertices();
  fminBrent->SetFunction(*func,TMath::Exp(-0.5*(TempLog)),TMath::Exp(-TempLog),1);
  if (! fminBrent->Minimize(10,0.1,0.1)) {
    LOG_WARN << "Temperature fit has failed" << endm;
    Temperature = 1;
  } else {
    Temperature = 1./fminBrent->XMinimum();
  }
  StAnneling::SetTemperature(Temperature);
  fcVertices->UniqueTracks2VertexAssociation(); // Make track associated with only vertex
  fcVertices->Fit(29,Canvas(),fVtx);
  if (Canvas()) Canvas()->Update();
}
//________________________________________________________________________________
Double_t StPicoD0AnaMaker::AnnelingFcn(Double_t TInv) {
  if (! fcVertices) return 0;
  Double_t Temperature = 1./TInv;
  StAnneling::SetTemperature(Temperature);
  Double_t Chi2 =  fcVertices->Fit();
  if (StKFVertex::Debug()) 
    LOG_INFO << "StKFVertexMaker::AnnelingFcn\tTemperature = " << Temperature << " Chi2 = " << Chi2 << endm;
  return Chi2;
}

void StPicoD0AnaMaker::Clear() {
  for (Int_t pass = 0; pass < fNPasses; pass++) {
    fVtxKs[pass]->Reset();
    fVtxKs[pass]->SetMaximum();
    fVtxs[pass]->Reset();
    fVtxs[pass]->SetMaximum();
  }
  fVtx = fVtxs[0]; // << switch between types    Vtx = fVtxKs[0];
  fVtxM->Reset();
  fcVertices = 0;
  fParticles->Clear("C");
}
