#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoD0AnaMaker.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StCuts.h"
#include "../StPicoPrescales/StPicoPrescales.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
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
#include "TMinuit.h"

#include "Sti/StiToolkit.h"
#include "TArrayI.h"
#include "StiMaker/StiDefaultToolkit.h"
#include <vector>
//
#include <stdio.h>
#include <time.h>

#include "StPicoKFVertexFitter/StPicoKFVertexFitter.h"
#include "TRandom.h"
vector<StDcaGeometry>     StPicoD0AnaMaker::mDCAs;
vector<StPhysicalHelixD>   StPicoD0AnaMaker::mHelices;
vector<UShort_t>           StPicoD0AnaMaker::mHelixFlags;
vector<Double_t >          StPicoD0AnaMaker::mSigma;
vector<Double_t >          StPicoD0AnaMaker::mZImpact;
Double_t                   StPicoD0AnaMaker::mWidthScale = 1;

ClassImp(StPicoD0AnaMaker)

  StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name,char const * inputFilesList, 
      char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil): 
    StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil), mOutFileName(outName), mInputFileList(inputFilesList),
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

  mPrescales = new StPicoPrescales(mycuts::prescalesFilesDirectoryName);
  int nRuns = mPrescales->numberOfRuns();

  mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
  mChain->SetBranchAddress("dEvent", &mPicoD0Event);

  mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
  mEventtuple = new TNtuple("mEventtuple","mEventtuple","kfx:kfy:kfz:kfSub1x:kfSub1y:kfSub1z:kfSub2x:kfSub2y:kfSub2z:minuitx:minuity:minuitz:minuitSub1x:minuitSub1y:minuitSub1z:minuitSub2x:minuitSub2y:minuitSub2z:grefMult:refMult:centrality:nHFTTracks:ZDCx:runId:eventId:nTrackFitting");
  mOutputFile->cd();

  //   if (!mHFCuts)
  //    mHFCuts = new StHFCuts;   

  // -------------- USER VARIABLES -------------------------
  dcaG = new StDcaGeometry();

  mMinuit = new TMinuit(3);         
  mMinuit->SetFCN(&StPicoD0AnaMaker::fcn);
  mMinuit->SetPrintLevel(-1);
  mMinuit->SetMaxIterations(1000);

  mGRefMultCorrUtil = new StRefMultCorr("grefmult");

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
  mEventtuple->Write();
  mOutputFile->Close();
  delete dcaG;
  delete mPrescales;
   delete mGRefMultCorrUtil;

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Make()
{
  readNextEvent();
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  //StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  picoDst = mPicoDstMaker->picoDst();
  clock_t t1,t2;

  if (!picoDst)
  {
    LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  if(mPicoD0Event->runId() != picoDst->event()->runId() ||
      mPicoD0Event->eventId() != picoDst->event()->eventId())
  {
    LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
    LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
    exit(1);
  }

  // -------------- USER ANALYSIS -------------------------
  TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();


  StThreeVectorF pVtx(-999.,-999.,-999.);
  StThreeVectorF minuit(-999.,-999.,-999.);
  StThreeVectorF minuitSub1(-999.,-999.,-999.);
  StThreeVectorF minuitSub2(-999.,-999.,-999.);
  StThreeVectorF kf(-999.,-999.,-999.);
  StThreeVectorF kfSub1(-999.,-999.,-999.);
  StThreeVectorF kfSub2(-999.,-999.,-999.);
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  if(!(isGoodEvent()) )
  {
    LOG_WARN << " Not Good Event! Skip! " << endm;
    return kStWarn;
  }
  float mult = event->grefMult();
  float const bField = event->bField();
  if(event) {
    pVtx = event->primaryVertex();
    minuit = event->primaryVertex();
  }
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }
  mGRefMultCorrUtil->init(picoDst->event()->runId());
  mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(),pVtx.z(),picoDst->event()->ZDCx()) ;

  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  const double reweight = mGRefMultCorrUtil->getWeight();
  const double refmultCor = mGRefMultCorrUtil->getRefMultCorr();
  kf = mPicoD0Event->kfVertex();


  vector<int> daughter;
  daughter.clear();//This daughter restore the index of tracking in fitting
  vector<int> daughterSub1;
  daughterSub1.clear();//This daughter restore the index of tracking in fitting
  vector<int> daughterSub2;
  daughterSub2.clear();//This daughter restore the index of tracking in fitting
  Int_t nTracks = picoDst->numberOfTracks();
  for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
  {
    StPicoTrack* trk = picoDst->track(iTrack);
    if (!trk ) continue;
    if(!isGoodForVertexFit(trk,pVtx)) continue;
    daughter.push_back(iTrack);
  } // .. end tracks loop
  for(int i=0;i<daughter.size();i++)
  {
    int seed = gRandom->Integer(2);
    if(seed==0)
    {
      if(daughterSub1.size()<=daughter.size())
        daughterSub1.push_back(daughter[i]);
      else
        daughterSub2.push_back(daughter[i]);
    }
    if(seed==1)
    {
      if(daughterSub2.size()<=daughter.size())
        daughterSub2.push_back(daughter[i]);
      else
        daughterSub1.push_back(daughter[i]);
    }
  }
 

  primaryVertexRefit(&kfSub1,daughterSub1);//Refit d0Vertex removing D0 daughters
  primaryVertexRefit(&kfSub2,daughterSub2);//Refit d0Vertex removing D0 daughters
  minuitSub1= vtxReFit(picoDst,daughterSub1);
  minuitSub2= vtxReFit(picoDst,daughterSub2);

  float refittuple_fill[100]; 
  refittuple_fill[0] = kf.x(); 
  refittuple_fill[1] = kf.y(); 
  refittuple_fill[2] = kf.z(); 
  refittuple_fill[3] = kfSub1.x(); 
  refittuple_fill[4] = kfSub1.y(); 
  refittuple_fill[5] = kfSub1.z(); 
  refittuple_fill[6] = kfSub2.x(); 
  refittuple_fill[7] = kfSub2.y(); 
  refittuple_fill[8] = kfSub2.z(); 
  refittuple_fill[9] = minuit.x(); 
  refittuple_fill[10] = minuit.y(); 
  refittuple_fill[11] = minuit.z(); 
  refittuple_fill[12] = minuitSub1.x(); 
  refittuple_fill[13] = minuitSub1.y(); 
  refittuple_fill[14] = minuitSub1.z(); 
  refittuple_fill[15] = minuitSub2.x(); 
  refittuple_fill[16] = minuitSub2.y(); 
  refittuple_fill[17] = minuitSub2.z(); 
  refittuple_fill[18] = event->grefMult();
  refittuple_fill[19] = event->refMult();
  refittuple_fill[20] = centrality;
  int nHFTTracks = event->numberOfPxlInnerHits()+event->numberOfPxlOuterHits()+event->numberOfIstHits()+event->numberOfSsdHits();
  refittuple_fill[21] = nHFTTracks; 
  refittuple_fill[22] = event->ZDCx(); 
  refittuple_fill[23] = event->runId();
  refittuple_fill[24] = event->eventId();
  refittuple_fill[25] = daughter.size();
  mEventtuple->Fill(refittuple_fill);
  return kStOK;
}
//-----------------------------------------------------------------------------


int StPicoD0AnaMaker::primaryVertexRefit(StThreeVectorF *mRefitVertex, vector<int>& daughter) 
{

  picoDst = mPicoDstMaker->picoDst();
  int nTracks = picoDst->numberOfTracks();
  int N = daughter.size();
  KFParticle *particles[N];
  Int_t NGoodGlobals = 0;
  for (int ii=0; ii < N; ii++) {
    int index = daughter[ii];
    StPicoTrack *gTrack = (StPicoTrack*)picoDst->track(index);
    if (! gTrack) continue;
    dcaG->set(gTrack->params(),gTrack->errMatrix());
    if (! dcaG) continue;
    Int_t kg = gTrack->id();
    Double_t xyzp[6], CovXyzp[21];
    dcaG->GetXYZ(xyzp,CovXyzp);
    static MTrack track;
    track.SetParameters(xyzp);
    track.SetCovarianceMatrix(CovXyzp);
    track.SetNDF(1);
    track.SetID(kg);
    Int_t q   = 1;
    Int_t pdg = 211;
    if (dcaG->charge() < 0) {
      q = -1;
      pdg = -211;
    } 
    track.SetCharge(q);
    particles[NGoodGlobals] = new KFParticle(track, pdg);
    particles[NGoodGlobals]->SetID(kg);
    NGoodGlobals++;
  }
  TArrayC Flag(N);
  KFVertex aVertex;
  aVertex.ConstructPrimaryVertex((const KFParticle **) particles, N,
      (Bool_t*) Flag.GetArray(),TMath::Sqrt(StAnneling::Chi2Cut()/2));
  //  delete [] particles;
  for(int i=0;i<N;++i)
  {
    delete particles[i];
  }

  if(aVertex.GetX()==0) return 0;
  mRefitVertex->set(aVertex.GetX(),aVertex.GetY(),aVertex.GetZ());
  return 1;
}


bool StPicoD0AnaMaker::isGoodEvent()
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  return (event->triggerWord() & mycuts::triggerWord) &&
    fabs(event->primaryVertex().z()) < mycuts::vz &&
    fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz;
  //return event->triggerWord() & mycuts::triggerWord;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
  // Require at least one hit on every layer of PXL and IST.
  // It is done here for tests on the preview II data.
  // The new StPicoTrack which is used in official production has a method to check this
  return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && trk->isHFTTrack();
}
//------------ Lines below are from Xin Dong's code to do Mnuit vertex fit -------------------------------
StThreeVectorF StPicoD0AnaMaker::vtxReFit(StPicoDst *picoDst, vector<int>& daughter)
{
  mNSeed = 0;
  mStatusMin = 0;
  StThreeVectorF XVertex(-999.,-999.,-999.);

  Double_t arglist[4];

  mDCAs.clear();
  mHelices.clear();
  mHelixFlags.clear();
  mSigma.clear();
  mZImpact.clear();

  Int_t nTracks = picoDst->numberOfTracks();
  Int_t nDaughter = daughter.size();
  for(Int_t ii=0;ii<nDaughter;ii++) {
    UInt_t k = daughter[ii];
    StPicoTrack* g = ( StPicoTrack*)picoDst->track(k);
    if (accept(g)) {
      mWidthScale = 0.1;// 1./TMath::Sqrt(5.);
      StDcaGeometry gDCA = g->dcaGeometry();
      // if (TMath::Abs(gDCA.impact()) >  mycuts::mRImpactMax) continue; // commented it out by Mustafa
      mDCAs.push_back(gDCA);
      StPhysicalHelixD helix = gDCA.helix();
      mHelices.push_back(helix);
      mHelixFlags.push_back(1);
      Double_t z_lin = gDCA.z();
      mZImpact.push_back(z_lin);
      mSigma.push_back(-1);
    }
  }

  if (mHelices.empty()) {
    LOG_WARN << "StPicoD0AnaMaker::fit: no tracks to fit." << endm;
    mStatusMin = -1;
    return XVertex;
  }
  LOG_INFO << "StPicoD0AnaMaker::fit size of helix vector: " << mHelices.size() << endm;

  //
  //  Reset and clear Minuit parameters
  // mStatusMin
  mMinuit->mnexcm("CLEar", 0, 0, mStatusMin);

  //
  //  Set parameters and start values. We do
  //  constrain the parameters since it harms
  //  the fit quality (see Minuit documentation).
  //
  // Initialize the seed with a z value which is not one of the discrete
  // values which it can tend to, implies zero not allowed.
  // Also need different initialization when vertex constraint.

  static Double_t step[3] = {0.03, 0.03, 0.03};

  //
  //  Scan z to find best seed for the actual fit.
  //  Skip this step if an external seed is given.
  //
  mNSeed = 1;  // use primary vertex as the seed

  mMinuit->mnexcm("CLEar", 0, 0, mStatusMin);

  Double_t seed_z = picoDst->event()->primaryVertex().z();

  mMinuit->mnparm(0, "x", 0, step[0], 0, 0, mStatusMin);
  mMinuit->mnparm(1, "y", 0, step[1], 0, 0, mStatusMin);
  mMinuit->mnparm(2, "z", seed_z, step[2], 0, 0, mStatusMin);

  Int_t done = 0;
  Int_t iter = 0;
  Double_t chisquare = 0;

  Int_t n_trk_vtx = 0;
  Int_t n_helix = mHelices.size();
  do {
    // For most vertices one pass is fine, but multiple passes
    // can be done
    n_trk_vtx = 0;
    for (Int_t i=0; i < n_helix; i++) {
      if (fabs(mZImpact[i]-seed_z) < mycuts::mDcaZMax) {
        mHelixFlags[i] |= kFlagDcaz;
        n_trk_vtx++;
      }
      else
        mHelixFlags[i] &= ~kFlagDcaz;
    }

    if (n_trk_vtx < mycuts::mMinTrack) {
      LOG_INFO << "Less than mMinTrack (=" << mycuts::mMinTrack << ") tracks, skipping vtx" << endm;
      continue;
    }
    mMinuit->mnexcm("MINImize", 0, 0, mStatusMin);
    done = 1;

    //
    //  Check fit result
    //

    if (mStatusMin) {
      LOG_WARN << "StPicoD0AnaMaker::fit: error in Minuit::mnexcm(), check status flag. ( iter=" << iter << endm;
      done = 0; // refit
    }

    Double_t fedm, errdef;
    Int_t npari, nparx;

    mMinuit->mnstat(chisquare, fedm, errdef, npari, nparx, mStatusMin);

    if (mStatusMin != 3) {
      LOG_INFO << "Warning: Minuit Status: " << mStatusMin << ", func val " << chisquare<< endm;
      done = 0;  // refit
    }
    mMinuit->mnhess();

    Double_t new_z, zerr;
    mMinuit->GetParameter(2, new_z, zerr);

    if (fabs(new_z - seed_z) > 1) // refit if vertex shifted
      done = 0;

    Int_t n_trk = 0;
    for (Int_t i=0; i < n_helix; i++) {
      if ( fabs(mZImpact[i] - new_z) < mycuts::mDcaZMax ) {
        n_trk++;
      }
    }
    if ( 10 * abs(n_trk - n_trk_vtx) >= n_trk_vtx ) // refit if number of track changed by more than 10%
      done = 0;

    iter++;
    seed_z = new_z; // seed for next iteration
  } while (!done && iter < 5 && n_trk_vtx >= mycuts::mMinTrack);


  if (n_trk_vtx < mycuts::mMinTrack)
    return XVertex;

  if (!done) {
    LOG_WARN << "Vertex unstable: no convergence after " << iter << " iterations. Skipping vertex " << endm;
    return XVertex;
  }

  // Store vertex
  Float_t cov[6];
  memset(cov,0,sizeof(cov));

  Double_t val, verr;
  XVertex = StThreeVectorF(mMinuit->fU[0],mMinuit->fU[1],mMinuit->fU[2]);
  Double_t emat[9];
  /* 0 1 2
     3 4 5
     6 7 8 */
  mMinuit->mnemat(emat,3);
  cov[0] = emat[0];
  cov[1] = emat[3];
  cov[2] = emat[4];
  cov[3] = emat[6];
  cov[4] = emat[7];
  cov[5] = emat[8];

  return XVertex;

}
//-----------------------------------------------------------------------------
Double_t StPicoD0AnaMaker::Chi2atVertex(StThreeVectorD &vtx) 
{
  static Int_t nCall=0; nCall++;
  Double_t f = 0;
  Double_t e;
  if (fabs(vtx.x())> 10) return 1e6;
  if (fabs(vtx.y())> 10) return 1e6;
  if (fabs(vtx.z())>300) return 1e6;
  for (UInt_t i=0; i<mDCAs.size(); i++) {
    if (mHelixFlags[i] & kFlagDcaz) {
      const StDcaGeometry gDCA = mDCAs[i];
      //      if (! gDCA) continue;
      const StPhysicalHelixD helix = gDCA.helix();
      e = helix.distance(vtx, kFALSE);  // false: don't do multiple loops
      //VP version
      //VP      Double_t chi2     = e*e/(errMatrix[0] + errMatrix[2]);
      static Int_t nCall=0;
      nCall++;
      Double_t err2;
      Double_t chi2 = gDCA.thelix().Dca(&(vtx.x()),&err2);
      chi2*=chi2/err2;
      //EndVP
      Double_t scale = 1./(mWidthScale*mWidthScale);
      f += scale*(1. - TMath::Exp(-chi2/scale)); // robust potential
      //        f -= scale*TMath::Exp(-chi2/scale); // robust potential
    }
  }
  return f;
}
//-----------------------------------------------------------------------------
void StPicoD0AnaMaker::fcn(int& npar, double* gin, double& f, double* par, Int_t iflag)
{
  StThreeVectorD vtx(par);
  f = Chi2atVertex(vtx);
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::accept(StPicoTrack* t) const
{
  return t;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodForVertexFit(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
  StPhysicalHelixD helix = trk->dcaGeometry().helix();
  float dca = (helix.at(helix.pathLength(vtx))-vtx).mag();

  if (dca > 3.0) return false;

  size_t numberOfFitPoints = popcount(trk->map0() >> 1); // drop first bit, primary vertex point
  numberOfFitPoints += popcount(trk->map1() & 0x1FFFFF); // only the first 21 bits are important, see StTrackTopologyMap.cxx

  return numberOfFitPoints >= 20;
}
