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
//
#include <vector>
//
#include <stdio.h>
#include <time.h>

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

   mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
   mChain->SetBranchAddress("dEvent", &mPicoD0Event);

   mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
   mEventTuple = new TNtuple("mEventTuple","","v2Hadron:sumCosCond:sumPairCon:sumCosBkg:mult");
   mOutputFile->cd();

//   if (!mHFCuts)
//    mHFCuts = new StHFCuts;   

   // -------------- USER VARIABLES -------------------------

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
   mEventTuple->Write();
 
   mOutputFile->Close();
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
   //StThreeVectorF Vertex(-999.,-999.,-999.);
   StPicoEvent *event = (StPicoEvent *)picoDst->event();
//   int aEventStat[mHFCuts->eventStatMax()];
//   if(!(mHFCuts->isGoodEvent(event,aEventStat)))
   if(!(isGoodEvent()))
   {
//     LOG_WARN << " Not Good Event! Skip! " << endm;
     return kStWarn;
   }
   if(event) {
     pVtx = event->primaryVertex();
   }

   vector<const StKaonPion *> signal;
   vector<const StKaonPion *> bkg;
   signal.clear();
   bkg.clear();
   for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
   {
      StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
      StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
      StPicoTrack const* pion = picoDst->track(kp->pionIdx());

      if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;
      if (!isTpcPion(pion)) continue;
      int charge=0;
     
      if((charge=isD0Pair(kp))!=0 && isTpcKaon(kaon,&pVtx))
      {
        if(charge<0)
          signal.push_back(kp);
        if(charge>0)
          bkg.push_back(kp);
      }
   
   }
   float mEventFill[5];
   mEventFill[0] = getHadronCorV2();
   int sumCondPair = 0;
   mEventFill[1] = getD0CorV2(&sumCondPair,signal);
   mEventFill[2] = sumCondPair;
   mEventFill[3] = getD0CorV2(&sumCondPair,bkg);
   mEventFill[4] = event->refMult();
   mEventTuple->Fill(mEventFill);
     
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



bool StPicoD0AnaMaker::isGoodEvent()
{
   StPicoEvent *event = (StPicoEvent *)picoDst->event();
   return (event->triggerWord() & mycuts::triggerWord) &&
           fabs(event->primaryVertex().z()) < mycuts::vz &&
         fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz;
//  return event->triggerWord() & mycuts::triggerWord;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
   // Require at least one hit on every layer of PXL and IST.
   // It is done here for tests on the preview II data.
   // The new StPicoTrack which is used in official production has a method to check this
   //return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit;
   return  trk->nHitsFit() >= mycuts::nHitsFit;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodHadron(StPicoTrack const * const trk) const
{
   return trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= mycuts::nHitsFit && trk->charge()!=0;
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


bool StPicoD0AnaMaker::getCorHadron(float eta,vector<float> &hadronsPhi) 
{
  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* hadron = picoDst->track(i);
    if(!isGoodHadron(hadron)) continue;
    float dEta = fabs(hadron->pMom().pseudoRapidity() - eta);
    if(dEta< mycuts::corDetaMin || dEta>mycuts::corDetaMax)  continue;
    hadronsPhi.push_back(hadron->pMom().phi());
  }
  fixPhi(hadronsPhi);
  return true;
  
}
    
float StPicoD0AnaMaker::sumCos(float phi,vector<float> &hadronsPhi) 
{
  float sumOfCos = 0;
  for(unsigned int i=0;i<hadronsPhi.size();++i)
  {
    sumOfCos += cos(2*(phi-hadronsPhi[i]));
  }
  return sumOfCos;
}
    
bool StPicoD0AnaMaker::fixPhi(vector<float> &phi) 
{
  if(phi.size() == 0) return false;
  float sumPhi = 0;
  for(unsigned int i=0;i<phi.size();i++)
    sumPhi+=phi[i];
  float meanPhi = sumPhi/phi.size();
  for(unsigned int i=0;i<phi.size();i++)
    phi[i] = phi[i]-meanPhi;  
  return true;
}

  
float StPicoD0AnaMaker::getHadronCorV2()
{
  int sumHadronPair = 0;
  float sumCosPair = 0;
  vector<float> hadron1Phi;
  vector<float> hadron1Eta;
  vector<float> hadron2Phi;
  hadron1Phi.clear();
  hadron1Eta.clear();
  hadron2Phi.clear();
  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* hadron1 = picoDst->track(i);
    if(!isGoodHadron(hadron1)) continue;
    if(hadron1->pMom().pseudoRapidity()>-0.5) continue;//
    hadron1Phi.push_back(hadron1->pMom().phi());
    hadron1Eta.push_back(hadron1->pMom().pseudoRapidity());
  }
  if(hadron1Phi.size()==0)  return 0;
  fixPhi(hadron1Phi);
  for(unsigned int i=0;i<hadron1Phi.size();i++)
  { 
    float eta1 = hadron1Eta[i];
    float phi1 = hadron1Phi[i];
    getCorHadron(eta1,hadron2Phi);
    sumHadronPair += hadron2Phi.size();
    sumCosPair += sumCos(phi1,hadron2Phi);
  }
  return sqrt(sumCosPair/sumHadronPair);
}
  
   
float StPicoD0AnaMaker::getD0CorV2(int *sumPair, vector<const StKaonPion *> cand)
{
  float sumCosPair = 0;
  vector<float> hadron1Phi;
  vector<float> hadron1Eta;
  vector<float> hadron2Phi;
  hadron1Phi.clear();
  hadron1Eta.clear();
  hadron2Phi.clear();
  for(unsigned int i=0;i<cand.size();++i)
  {
    hadron1Phi.push_back(cand[i]->phi());
    hadron1Eta.push_back(cand[i]->eta());
  }
  if(hadron1Phi.size()==0)  return 0;
  fixPhi(hadron1Phi);
  for(unsigned int i=0;i<hadron1Phi.size();i++)
  { 
    float eta1 = hadron1Eta[i];
    float phi1 = hadron1Phi[i];
    getCorHadron(eta1,hadron2Phi);
    *sumPair += hadron2Phi.size();
    sumCosPair += sumCos(phi1,hadron2Phi);
  }
  return sumCosPair/(*sumPair);
}
  
    

























