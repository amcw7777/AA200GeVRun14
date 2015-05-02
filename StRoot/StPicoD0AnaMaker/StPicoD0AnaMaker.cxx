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
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StCuts.h"
#include "../StPicoPrescales/StPicoPrescales.h"
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

vector<StDcaGeometry>     StMyAnalysisMaker::mDCAs;
vector<StPhysicalHelixD>   StMyAnalysisMaker::mHelices;
vector<UShort_t>           StMyAnalysisMaker::mHelixFlags;
vector<Double_t >          StMyAnalysisMaker::mSigma;
vector<Double_t >          StMyAnalysisMaker::mZImpact;
Double_t                   StMyAnalysisMaker::mWidthScale = 1;

ClassImp(StPicoD0AnaMaker)

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

   mPrescales = new StPicoPrescales(mycuts::prescalesFilesDirectoryName);
   int nRuns = mPrescales->numberOfRuns();

   mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
   mChain->SetBranchAddress("dEvent", &mPicoD0Event);

   mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
   mEventtuple = new TNtuple("mEventtuple","mEventtuple","testx:testy:testz:refitx:refity:refitz:prmx:prmy:prmz:mult");
   mOrigin= new TNtuple("mOrigin","mOrigin","mass:charge:mult");
   mTest= new TNtuple("mTest","mTest","mass:charge:mult");
   mRefit= new TNtuple("mRefit","mRefit","mass:charge:mult");
  mDmass_unlike = new TH1D("mDmass_unlike","",50,1.6,2.1);
  mDmass_like = new TH1D("mDmass_like","",500,1.6,2.1);
  mDmasscut_unlike = new TH1D("mDmasscut_unlike","",500,1.6,2.1);
  mDmasscut_like = new TH1D("mDmasscut_like","",500,1.6,2.1);
  mDmasstest_unlike = new TH1D("mDmasstest_unlike","",500,1.6,2.1);
  mDmasstest_like = new TH1D("mDmasstest_like","",500,1.6,2.1);
  mMult = new TH1D("mMult","",500,0,500);

   mOutputFile->cd();

//   if (!mHFCuts)
//    mHFCuts = new StHFCuts;   

   // -------------- USER VARIABLES -------------------------
  dcaG = new StDcaGeometry();
  
  mMinuit = new TMinuit(3);         
  mMinuit->SetFCN(&StMyAnalysisMaker::fcn);
  mMinuit->SetPrintLevel(-1);
  mMinuit->SetMaxIterations(1000);


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
   mOrigin->Write();
   mTest->Write();
   mRefit->Write();

  mMult->Write();
//  mDmass_unlike->Write();
//  mDmass_like->Write();
//  mDmasscut_unlike->Write();
//  mDmasscut_like->Write();
//  mDmasstest_unlike->Write();
//  mDmasstest_like->Write();
 

 
   mOutputFile->Close();
   delete dcaG;
   delete mPrescales;
	
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
   StThreeVectorF testVertex(-999.,-999.,-999.);
   StThreeVectorF d0Vertex(-999.,-999.,-999.);
   //StThreeVectorF Vertex(-999.,-999.,-999.);
   StPicoEvent *event = (StPicoEvent *)picoDst->event();
//   int aEventStat[mHFCuts->eventStatMax()];
//   if(!(mHFCuts->isGoodEvent(event,aEventStat)))
   if(!(isGoodEvent()))
   {
//     LOG_WARN << " Not Good Event! Skip! " << endm;
     return kStWarn;
   }
   float mult = event->refMult();
   float const bField = event->bField();
   mMult->Fill(event->refMult());
   if(event) {
     pVtx = event->primaryVertex();
   }
   testVertex = pVtx;
   d0Vertex = pVtx;

   vector<int> daughter;
   daughter.clear();
   primaryVertexRefit(&testVertex,daughter);//refit vertex using all tracks

   for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
   {
      StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
      if(!isGoodPair(kp)) continue;
      daughter.push_back(kp->kaonIdx());
      daughter.push_back(kp->pionIdx());
   }

   t1 = clock();
   primaryVertexRefit(&d0Vertex,daughter);//Refit d0Vertex removing D0 daughters
   t2 = clock();
   for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
   {
      StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
      StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
      StPicoTrack const* pion = picoDst->track(kp->pionIdx());
      StKaonPion testkp(kaon,pion,kp->kaonIdx(),kp->pionIdx(),testVertex,bField);
      StKaonPion d0kp(kaon,pion,kp->kaonIdx(),kp->pionIdx(),d0Vertex,bField);

      if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;
      if (!isTpcPion(pion)) continue;
      int charge=0;
      float mDmass_fill[3];
     
      if((charge=isD0Pair(kp))!=0 && isTpcKaon(kaon,&pVtx))
      {
	mDmass_fill[0]=kp->m();
	mDmass_fill[1]=charge;
	mDmass_fill[2]=mult;
        mOrigin->Fill(mDmass_fill);
        
      }
   
      if((charge=isD0Pair(&testkp))!=0 && isTpcKaon(kaon,&testVertex))
      {
	mDmass_fill[0]=testkp.m();
	mDmass_fill[1]=charge;
	mDmass_fill[2]=mult;
        mTest->Fill(mDmass_fill);
      }
      
      if((charge=isD0Pair(&d0kp))!=0 && isTpcKaon(kaon,&d0Vertex))
      {
	mDmass_fill[0]=d0kp.m();
	mDmass_fill[1]=charge;
	mDmass_fill[2]=mult;
        mRefit->Fill(mDmass_fill);
      }
   }
     
   float refittuple_fill[20]; 
   refittuple_fill[0] = testVertex.x(); 
   refittuple_fill[1] = testVertex.y(); 
   refittuple_fill[2] = testVertex.z(); 
   refittuple_fill[3] = d0Vertex.x(); 
   refittuple_fill[4] = d0Vertex.y(); 
   refittuple_fill[5] = d0Vertex.z(); 
   refittuple_fill[6] = pVtx.x(); 
   refittuple_fill[7] = pVtx.y(); 
   double dtime = (double)(t2- t1) / CLOCKS_PER_SEC;
   refittuple_fill[8] = pVtx.z(); 
   refittuple_fill[9] = mult; 
   mEventtuple->Fill(refittuple_fill);
   return kStOK;
}
//-----------------------------------------------------------------------------

bool StPicoD0AnaMaker::isGoodPair(StKaonPion const* const kp) const
{
  if(!kp) return false;

  StPicoTrack const* kaon = mPicoDstMaker->picoDst()->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDstMaker->picoDst()->track(kp->pionIdx());

  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  StThreeVectorF pVtx(-999.,-999.,-999.);
  pVtx = event->primaryVertex();
  int charge = kaon->charge() * pion->charge();
  bool pairCuts = kp->m()>1.6 && kp->m()<2.1 &&
                  charge==-1;

  return (isTpcKaon(kaon,&pVtx) && isTpcPion(pion) && 
	  pairCuts);
}


int StPicoD0AnaMaker::isD0Pair(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts =  cos(kp->pointingAngle()) > mycuts::cosTheta &&
          kp->pionDca() > mycuts::pDca && kp->kaonDca() > mycuts::kDca &&
          kp->dcaDaughters() < mycuts::dcaDaughters;
  int charge = kaon->charge() * pion->charge();
    

  if(pairCuts)
    return charge;
  else
    return 0;
}
/*
*/

int StPicoD0AnaMaker::primaryVertexRefit(StThreeVectorF *mRefitVertex, vector<int>& daughter) 
{

  picoDst = mPicoDstMaker->picoDst();
  int nTracks = picoDst->numberOfTracks();
  int N = 0;//nTracks - daughter.size();
  for (int i=0; i < nTracks; i++) {
     StPicoTrack *gTrack = (StPicoTrack*)picoDst->track(i);
     if (! gTrack) continue;
     Int_t kg = gTrack->id();

     bool flagDdaughterCand = false;
     for (vector<int>::size_type j=0; j < daughter.size(); j++) {
	if (daughter[j] == kg) 
        {
          flagDdaughterCand = 1;
        }
     }
     if (flagDdaughterCand == 1) continue;
     N++;
  }
  KFParticle *particles[N];
//  KFParticle **particles = new KFParticle*[N];
//  StKFVertexMaker fitter;
  Int_t NGoodGlobals = 0;
  for (int i=0; i < nTracks; i++) {
     StPicoTrack *gTrack = (StPicoTrack*)picoDst->track(i);
     if (! gTrack) continue;
     dcaG->set(gTrack->params(),gTrack->errMatrix());
     if (! dcaG) continue;
     Int_t kg = gTrack->id();

     bool flagDdaughterCand = false;
     for (vector<int>::size_type j=0; j < daughter.size(); j++) {
	if (daughter[j] == kg) 
        {
          flagDdaughterCand = 1;
        }
     }
     if (flagDdaughterCand == 1) continue;


    // particles[NGoodGlobals] = fitter.AddTrackAt(dcaG,kg);
    if (! dcaG) return 0;
   Double_t xyzp[6], CovXyzp[21];
   dcaG->GetXYZ(xyzp,CovXyzp);
   static MTrack track;
   track.SetParameters(xyzp);
   track.SetCovarianceMatrix(CovXyzp);
   track.SetNDF(1);
   //    track.SetChi2(GlobalTracks_mChiSqXY[k]);
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
//////////////////
     NGoodGlobals++;
  }
  TArrayC Flag(N);
  KFVertex aVertex;
  aVertex.ConstructPrimaryVertex((const KFParticle **) particles, N,
                                  (Bool_t*) Flag.GetArray(),TMath::Sqrt(StAnneling::Chi2Cut()/2));
//  delete [] particles;
  for(int i=0;i<N;++i)
  {
//    particles[i]->Clear();
    delete particles[i];
  }

  if(aVertex.GetX()==0) return 0;
  mRefitVertex->set(aVertex.GetX(),aVertex.GetY(),aVertex.GetZ());
  return 1;
}


bool StPicoD0AnaMaker::isGoodEvent()
{
   StPicoEvent *event = (StPicoEvent *)picoDst->event();
//   return (event->triggerWord() & mycuts::triggerWord) &&
//           fabs(event->primaryVertex().z()) < mycuts::vz &&
//         fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz;
  return event->triggerWord() & mycuts::triggerWord;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
   // Require at least one hit on every layer of PXL and IST.
   // It is done here for tests on the preview II data.
   // The new StPicoTrack which is used in official production has a method to check this
   return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcPion(StPicoTrack const * const trk) const
{
   return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcKaon(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
  float kBeta = getTofBeta(trk,pVtx);
  bool tofAvailable = kBeta>0;
  bool tofKaon = tofAvailable && isTofKaon(trk,kBeta);
   return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon
           || tofKaon;
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
   StPhysicalHelixD helix = trk->helix();

   int index2tof = trk->bTofPidTraitsIndex();

   float beta = std::numeric_limits<float>::quiet_NaN();

   if(index2tof >= 0)
   {
      StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

      if(tofPid)
      {
         beta = tofPid->btofBeta();

         if (beta < 1e-4)
         {
            StThreeVectorF const btofHitPos = tofPid->btofHitPos();

            float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();
         }
         else
         {
           beta = std::numeric_limits<float>::quiet_NaN();
         }
      }
   }

   return beta;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
   bool tofKaon = false;

   if(beta>0)
   {
     double ptot = trk->dcaGeometry().momentum().mag();
     float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);
     tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;
   }
   
   return tofKaon;
}
//------------ Lines below are from Xin Dong's code to do Mnuit vertex fit -------------------------------
StThreeVectorD StPicoD0AnaMaker::vtxReFit(StPicoDst *picoDst)
{
    mNSeed = 0;
    mStatusMin = 0;

    Double_t arglist[4];

    mDCAs.clear();
    mHelices.clear();
    mHelixFlags.clear();
    mSigma.clear();
    mZImpact.clear();

    Int_t nTracks = picoDst->numberOfTracks();
    for(UInt_t k=0;k<nTracks;k++) {
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
        return 0;
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

      if (n_trk_vtx < mMinTrack) {
        LOG_INFO << "Less than mMinTrack (=" << mMinTrack << ") tracks, skipping vtx" << endm;
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

    StThreeVectorD XVertex(0.,0.,0.);

    if (n_trk_vtx < mMinTrack)
      return XVertex;

    if (!done) {
      LOG_WARN << "Vertex unstable: no convergence after " << iter << " iterations. Skipping vertex " << endm;
      return XVertex;
    }

    // Store vertex
    Float_t cov[6];
    memset(cov,0,sizeof(cov));

    Double_t val, verr;
      XVertex = StThreeVectorD(mMinuit->fU[0],mMinuit->fU[1],mMinuit->fU[2]);
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
