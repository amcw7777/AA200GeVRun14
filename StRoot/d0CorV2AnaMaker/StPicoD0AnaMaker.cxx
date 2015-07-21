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
#include "TH3D.h"
#include "TFile.h"
#include "StEvent/StDcaGeometry.h"
//
#include <vector>
//
#include <stdio.h>
#include <time.h>

ClassImp(StPicoD0AnaMaker)

  StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name,char const * inputFilesList, 
      char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil): 
    StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil),
    mOutFileName(outName), mInputFileList(inputFilesList),mOutputFile(NULL), mChain(NULL), mEventCounter(0){}

Int_t StPicoD0AnaMaker::Init()
{
  mPicoD0Event = new StPicoD0Event();
  fout.open("check.txt");

  mChain = new TChain("T");
  std::ifstream listOfFiles(mInputFileList.Data());
  if (listOfFiles.is_open())
  {
    std::string file;
    while (getline(listOfFiles, file))
    {
      LOG_INFO << "StPicoD0AnaMaker - Adding :" << file <<endm;
      mChain->Add(file.c_str());
      LOG_INFO<<" Entries = "<<mChain->GetEntries()<< endm; 
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
  mHadronTuple = new TNtuple("mHadronTuple","","sum1:sin1:cos1:sum2:sin2:cos2:mult");
  etaPhi = new TH2F("etaPhi","",200,-6.29,6.29,100,0.5,1.5);
  etaPhi_D = new TH2F("etaPhiD","",100,-3.1416,3.1416,100,-2,2);
  etaPhi_Hadron = new TH2F("etaPhiHadron","",100,-3.1416,3.1416,100,-2,2);
  etaPhi_Hadron_all = new TH2F("etaPhiHadronAll","",100,-3.1416,3.1416,100,-2,2);
  dEtaDHadron = new TH1F("dEtaDHadron","",1000,0,10);
  hEtaD = new TH1F("hEtaD","",1000,0,10);
  hEtaHadron = new TH1F("hEtaHadron","",1000,0,10);
  TString flatten[5];
  flatten[0] = "v2";
  flatten[1] = "cosD";
  flatten[2] = "sinD";
  flatten[3] = "cosHadron";
  flatten[4] = "sinHadron";
  for(int i=0;i!=6;i++)
  {
    for(int j=0;j!=5;j++)
    {
      TString sign;
      //int gap = static_cast<int>i/2;
      int gap = (int)i/2;
      if(i%2) sign = "like";
      else sign = "unlike";
      TString name = flatten[j]+"_"+sign+Form("_%i",gap);
      profV2[i][j] = new TProfile(name.Data(),"",10,0,10);
    }
   }

  mOutputFile->cd();


  // -------------- USER VARIABLES -------------------------
  mGRefMultCorrUtil = new StRefMultCorr("grefmult");

  return kStOK;
}
//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
  /*  */
  delete mGRefMultCorrUtil;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
  LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
  fout.close();
  mOutputFile->cd();
  // save user variables here
  mHadronTuple->Write();
  dEtaDHadron->Write();
  hEtaD->Write();
  hEtaHadron->Write();
  for(int i=0;i!=6;i++)
  {
    for(int j=0;j!=5;j++)
    {
      profV2[i][j]->Write();
    }
  }

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
  if(!(isGoodEvent()))
  {
    //     LOG_WARN << " Not Good Event! Skip! " << endm;
    return kStWarn;
  }
  if(event) {
    pVtx = event->primaryVertex();
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
    //if(!isTpcKaon(kaon,&pVtx)) continue;  
    bool tpcKaon = isTpcKaon(kaon,&pVtx);
    float kBeta = getTofBeta(kaon,&pVtx);
    bool tofAvailable = kBeta>0;
    bool tofKaon = tofAvailable && isTofKaon(kaon,kBeta);
    bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
    if(!goodKaon) continue;
    int charge=0;
    if(kp->m()<1.85 || kp->m()>1.88)  continue;
    if((charge=isD0Pair(kp))!=0 )
    {
      float d0Fill[20] = {0}; 
      float d0Phi = kp->phi();
      float d0Eta = kp->eta();//hadron->pMom().pseudoRapidity();
      if(kp->pt()>10) continue;
      d0Fill[0] = d0Phi; 
      vector<float> hadronPhi[3];
      int index1 = kp->kaonIdx();
      int index2 = kp->pionIdx();
      float etaGap[3] = {0,0.1,0.75};
      vector<float>::iterator it; 
      for(int i=0;i<3;i++)
      {
        if(charge==1)
        {
          profV2[i*2][1]->Fill(kp->pt(),cos(d0Phi));
          profV2[i*2][2]->Fill(kp->pt(),sin(d0Phi));
        }
        if(charge==-1)
        {
          profV2[i*2+1][1]->Fill(kp->pt(),cos(d0Phi));
          profV2[i*2+1][2]->Fill(kp->pt(),sin(d0Phi));
        }
        getCorHadron(d0Eta,hadronPhi[i],index1,index2,d0Phi,etaGap[i]);
        for(it=hadronPhi[i].begin(); it!=hadronPhi[i].end(); it++)
        {
          double hadronPhi = *it;
          double cos2Phi = cos(2*(d0Phi-hadronPhi)); 
          if(charge==1)
          { 
            profV2[i*2][0]->Fill(kp->pt(),cos2Phi);
            profV2[i*2][3]->Fill(kp->pt(),cos(hadronPhi));
            profV2[i*2][4]->Fill(kp->pt(),sin(hadronPhi));
          }
          if(charge==-1)
          {
            profV2[i*2+1][0]->Fill(kp->pt(),cos2Phi);
            profV2[i*2+1][3]->Fill(kp->pt(),cos(hadronPhi));
            profV2[i*2+1][4]->Fill(kp->pt(),sin(hadronPhi));
          }
        }//Hadron loop
      }//different eta gap 
    }//D loop

  }
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
  //  bool pairCuts =  cos(kp->pointingAngle()) > mycuts::cosTheta &&
  //          kp->pionDca() > mycuts::pDca && kp->kaonDca() > mycuts::kDca &&
  //          kp->dcaDaughters() < mycuts::dcaDaughters;  
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0062 &&
      kp->pionDca() > 0.0109 && kp->kaonDca() > 0.0123 &&
      kp->dcaDaughters() < 0.0082 && kp->decayLength()>0.0149;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0047 &&
      kp->pionDca() > 0.0108 && kp->kaonDca() > 0.0097 &&
      kp->dcaDaughters() < 0.0070 && kp->decayLength()>0.0205;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0040 &&
      kp->pionDca() > 0.0100 && kp->kaonDca() > 0.0091 &&
      kp->dcaDaughters() < 0.0056 && kp->decayLength()>0.0216;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0041 &&
      kp->pionDca() > 0.0074 && kp->kaonDca() > 0.0075 &&
      kp->dcaDaughters() < 0.0065 && kp->decayLength()>0.0233;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0042 &&
      kp->pionDca() > 0.067 && kp->kaonDca() > 0.0053 &&
      kp->dcaDaughters() < 0.0065 && kp->decayLength()>0.0282;  
  }

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
  return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && trk->isHFTTrack();
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodHadron(StPicoTrack const * const trk) const
{
  return trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= mycuts::nHitsFit && trk->charge()!=0;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcPion(StPicoTrack const * const trk) const
{
  return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcKaon(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
  return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon;
  //      || tofKaon;
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{

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
        StPhysicalHelixD helix = trk->helix();

        float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
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


bool StPicoD0AnaMaker::getCorHadron(float eta,vector<float> &hadronsPhi, int index1,int index2, float phi, float etaCut) 
{
  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* hadron = picoDst->track(i);
    if(!hadron)  continue; 
    if(hadron->pMom().perp()<0.2) continue;
    etaPhi_Hadron_all->Fill(hadron->pMom().phi(),hadron->pMom().pseudoRapidity());
    if(i==index1 || i==index2) continue;
    if(!isGoodHadron(hadron)) continue;
    float dEta = fabs(hadron->pMom().pseudoRapidity() - eta);
    float dPhi = (hadron->pMom().phi() - phi);
    if(etaCut<0.001)
    {
      dEtaDHadron->Fill(dEta);
      hEtaD->Fill(eta);
      hEtaHadron->Fill(hadron->pMom().pseudoRapidity());
    }
    //if(dPhi>3.1416) dPhi = 2*3.1416-dPhi;
    if(dEta< etaCut|| dEta > mycuts::corDetaMax)  continue;
    etaPhi->Fill(dPhi,dEta);
    etaPhi_Hadron->Fill(hadron->pMom().phi(),hadron->pMom().pseudoRapidity());
    hadronsPhi.push_back(hadron->pMom().phi());
  }
  //  fixPhi(hadronsPhi);
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


bool StPicoD0AnaMaker::getHadronCorV2()
{
  float hadronFill[7] = {0};
  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* hadron = picoDst->track(i);
    if(hadron->pMom().perp()<0.2) continue;
    if(!isGoodHadron(hadron)) continue;
    float etaHadron = hadron->pMom().pseudoRapidity();
    float phiHadron = hadron->pMom().phi();
    if(etaHadron<-0.375)
    {
      hadronFill[0]++;
      hadronFill[1] += sin(2 * phiHadron);
      hadronFill[2] += cos(2 * phiHadron);
    }			
    if(etaHadron>0.375)
    {
      hadronFill[3]++;
      hadronFill[4] += sin(2 * phiHadron);
      hadronFill[5] += cos(2 * phiHadron);
    }			
    hadronFill[6] = picoDst->event()->grefMult();
  }
  mHadronTuple->Fill(hadronFill);
  return true;
}

/* 
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
   getCorHadron(eta1,hadron2Phi,-1,-1,0);
 *sumPair += hadron2Phi.size();
 sumCosPair += sumCos(phi1,hadron2Phi);
 }
 return sumCosPair/(*sumPair);
 }

*/ 

























