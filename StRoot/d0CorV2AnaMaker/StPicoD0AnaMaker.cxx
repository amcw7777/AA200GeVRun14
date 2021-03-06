#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TLorentzVector.h"

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
#include <algorithm>

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
  mHadronTuple = new TNtuple("mHadronTuple","","sum1:sin1:cos1:sum2:sin2:cos2:cent:reweight");
  etaPhi = new TH2F("etaPhi","",200,-6.29,6.29,100,0.5,1.5);
  etaPhi_D = new TH2F("etaPhiD","",100,-3.1416,3.1416,100,-2,2);
  etaPhi_Hadron = new TH2F("etaPhiHadron","",100,-3.1416,3.1416,100,-2,2);
  etaPhi_Hadron_all = new TH2F("etaPhiHadronAll","",100,-3.1416,3.1416,100,-2,2);
  dEtaDHadron = new TH1F("dEtaDHadron","",1000,0,10);
  hEtaD = new TH1F("hEtaD","",1000,0,10);
  hEtaHadron = new TH1F("hEtaHadron","",1000,0,10);
  hPhiHadron = new TH1F("hPhiHadron","",10000,-5,5);
  hPhiD= new TH1F("hPhiD","",10000,-5,5);
  vtxz = new TH1F("vtxz","",100,-10,10);
  TString flatten[5];
  flatten[0] = "v2";
  flatten[1] = "cosD";
  flatten[2] = "sinD";
  flatten[3] = "cosHadron";
  flatten[4] = "sinHadron";
  TString sb[8] = {"s1like","s3like","hSBlike","lSBlike","s1unlike","s3unlike","hSBunlike","lSBunlike"};
  float xbin[7] = {0,1,2,3,4,5,10};
  float binMass[2001];
  for(int i=0;i<2001;i++)
    binMass[i] = 0.01*i;
  massPt = new TH2D("massPt","",2000,binMass,6,xbin);
  massPtLike = new TH2D("massPtLike","",250,0,2.5,100,0,10);
  massLike = new TH2D("massLike1","",250,0,2.5,100,0,10);
  massLike2 = new TH2D("massLike2","",250,0,2.5,100,0,10);
  //massPtLike = new TH2D("massPtLike","",2000,binMass,6,xbin);
  //massLike1 = new TH2D("massLike1","",2000,binMass,6,xbin);
  //massLike2 = new TH2D("massLike2","",2000,binMass,6,xbin);
  massUnlike = new TH2D("massUnlike","",2000,binMass,6,xbin);
  for(int i=0;i!=8;i++)
  {
    for(int k=0;k!=3;k++)
    {
      for(int j=0;j!=5;j++)
      {
        TString name = sb[i]+flatten[j]+Form("_%i",k);
        profV2[i][j][k] = new TProfile(name.Data(),"",6,xbin);
      }
      TString weightName = sb[i]+Form("_%i_weigth",k);
      float xWeight[10];
      for(int ii=0;ii<10;ii++)
        xWeight[ii] = ii;
      v2Weight[i][k] = new TH2D(weightName.Data(),"",9,xWeight,6,xbin);

    }
  }
  float ptbin1[12] = {0.225,0.375,0.525,0.675,0.825,0.975,1.12,1.27,1.42,1.58,1.73,1.88};
  float ptbin2[11];
  for(int i=0;i<11;i++)
    ptbin2[i] = 0.5*(ptbin1[i]+ptbin1[i+1]);
  for(int i=0;i<5;i++)
  {
    hadronV2[i] = new TH1D(Form("hadron_%s",flatten[i].Data()),"",9,0,9);
    hadronV2[i]->Sumw2();
    hadronV2_sum[i] = new TH1D(Form("hadronsum_%s",flatten[i].Data()),"",9,0,9);
    hadronV2_sum[i]->Sumw2();
    for(int j=0;j<9;j++)
    {
      hadronV2_excl[i][j] = new TH1D(Form("hadron_%s_cent%i",flatten[i].Data(),j),"",10,ptbin2);
      hadronV2_excl[i][j]->Sumw2();
      hadronV2_excl_sum[i][j] = new TH1D(Form("hadronsum_%s_cent%i",flatten[i].Data(),j),"",10,ptbin2);
      hadronV2_excl_sum[i][j]->Sumw2();
    }
  }
 // fitmean = {1.8602,1.8626,1.8639,1.8645,1.8619,1.8619};
  //fitsigma = {0.0149,0.0146,0.0153,0.0167,0.0184,0.0184};
  fitmean = {1.85921,1.8633,1.86403,1.86475,1.86252,1.86534};
  fitsigma = {0.018139,0.0139476,0.0158346,0.0169282,0.0199567,0.0189131};
  ifstream ifs("efficiency.txt");
  for(int i=0; i<6; i++)
    for(int j=0; j<4; j++)
      ifs>>efficiency[j][i];
  //efficiency = {0.000199608, 0.000148673, 0.000225232, 0.000328405, 0.000273876, 0.000187721, 0.000371691, 0.000429313, 0.00066452, 0.00060398, 0.000917251, 0.00100015,0.00135655,0.00119798,0.00179572,0.00235767,0.00237739,0.00188794,0.00301619,0.00314772,0.0044066,0.00278407,0.00455717,0.00566198};
  for(int i=0;i<6;i++)//pt bin
  {
    for(int j=0;j<5;j++)//flatten
    {
      TString likename1 = Form("likeMass%i",i)+flatten[j];
      TString likename2 = Form("likeMass1%i",i)+flatten[j];
      TString unlikename = Form("unlikeMass%i",i)+flatten[j];
      likeV2Mass[i][j] = new TH2D(likename1.Data(),"",18,fitmean[i]-9*fitsigma[i],fitmean[i]+9*fitsigma[i],1000,-50,50);
      likeV2Mass2[i][j] = new TH2D(likename2.Data(),"",18,fitmean[i]-9*fitsigma[i],fitmean[i]+9*fitsigma[i],1000,-1,1);
      unlikeV2Mass[i][j] = new TH2D(unlikename.Data(),"",18,fitmean[i]-9*fitsigma[i],fitmean[i]+9*fitsigma[i],1000,-50,50);
    }
  }
  checkNew = new TNtuple("checkNew","","dphi:mass:pT:charge");
  checkOld = new TNtuple("checkOld","","dphi:mass:pT:charge");
  checkPeak = new TH2D("checkPeak","",18,fitmean[4]-9*fitsigma[4],fitmean[4]+9*fitsigma[4],100,-1,1);
  

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
  //mHadronTuple->Write();
//  dEtaDHadron->Write();
//  hEtaD->Write();
//  hEtaHadron->Write();
 // hPhiHadron->Write();
//  hPhiD->Write();
  massPt->Write();
  massPtLike->Write();
 // checkOld->Write();
//  checkNew->Write();
  vtxz->Write();
  for(int i=0;i!=8;i++)
  {
    for(int k=0;k!=3;k++)
    {
      for(int j=0;j!=5;j++)
      {
        profV2[i][j][k]->Write();
      }
      v2Weight[i][k]->Write();
    }
  }
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<5;j++)
    {
      likeV2Mass[i][j]->Write();
      likeV2Mass2[i][j]->Write();
      unlikeV2Mass[i][j]->Write();
    }
  }
    massLike->Write();
    massLike2->Write();
    massUnlike->Write();
  for(int i=0;i<5;i++)
  {
    hadronV2[i]->Write();
    hadronV2_sum[i]->Write();
    for(int j=0;j<9;j++)
    {
    //hadronV2_excl[i][j]->Write();
    //hadronV2_excl_sum[i][j]->Write();
    }
  }
    
  checkPeak->Write();
      

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
  vtxz->Fill(pVtx.z());
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }
  mGRefMultCorrUtil->init(picoDst->event()->runId());
  mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(),pVtx.z(),picoDst->event()->ZDCx()) ;

  getHadronCorV2();
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  int centBin = 0;
  if(centrality>=7) centBin=1;
  else if(centrality>=4)  centBin=2;
  else centBin=3;
  //if(centrality>1)
   // return kStWarn;
    
  double reweight = mGRefMultCorrUtil->getWeight();
  //const double refmultCor = mGRefMultCorrUtil->getRefMultCorr();

  vector<const StKaonPion *> signal;
  vector<const StKaonPion *> bkg;
  signal.clear();
  bkg.clear();
      
  vector<unsigned int> d0Daughter_l;
  vector<unsigned int> d0Daughter_u;
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
    if((charge=isD0Pair(kp))!=0 )
    {
      unsigned int index1 = kp->kaonIdx();
      unsigned int index2 = kp->pionIdx();
      if(charge==-1)
      {
        d0Daughter_u.push_back(index1);
        d0Daughter_u.push_back(index2);
      }
      if(charge>0)
      {
        d0Daughter_l.push_back(index1);
        d0Daughter_l.push_back(index2);
      }
    }
  }
  
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

    TLorentzVector rec_d[3];
    StThreeVectorF kaonMom = kaon->gMom(pVtx,event->bField());
    StThreeVectorF pionMom = pion->gMom(pVtx,event->bField());
    
    rec_d[0].SetXYZM(kaonMom.x(),kaonMom.y(),kaonMom.z(),0.493677);
    rec_d[1].SetXYZM(pionMom.x(),pionMom.y(),pionMom.z(),0.139570);
    int rec_seed = 0;//gRandom->Integer(2);
    //rec_d[rec_seed].SetPtEtaPhiE(rec_d[rec_seed].Pt(),rec_d[rec_seed].Eta(),rec_d[rec_seed].Phi()>3.1416 ? rec_d[rec_seed].Phi()-3.1416 : rec_d[rec_seed].Phi()+3.1416,rec_d[rec_seed].E());
   // rec_d[2] = rec_d[0]+rec_d[1];
    int charge=0;
    float d0Phi = kp->phi();
    float d0Eta = kp->eta();//hadron->pMom().pseudoRapidity();
    float d0Pt = kp->pt();
    double dMass = kp->m();

//    float d0Phi = rec_d[2].Phi();
//    float d0Eta = rec_d[2].Eta();
//    float d0Pt= rec_d[2].Pt();
//    float dMass= rec_d[2].M();
    if(d0Pt>10) continue;
    int fitindex = 5;
    if(d0Pt<5)
      fitindex = static_cast<int>(d0Pt);
    double mean = fitmean[fitindex];
    double sigma = fitsigma[fitindex];
    double hadronV2[9] = {0.0557348,0.0630083,0.0704378,0.0742203,0.0725183,0.0644258,0.0494283,0.0349057,0.024481};
    double weight_hadron = 1.0/hadronV2[centBin];
    double reweight_eff = (efficiency[0][fitindex]/efficiency[centBin][fitindex]);
     //cout<<"cent = "<<centrality<<"\t"<<"pt = "<<d0Pt<<"\teff = "<<(efficiency[0][fitindex]/efficiency[centBin][fitindex])<<"\tcentbin = "<<centBin<<endl;

    //cout<<"reweight_eff  = "<<reweight_eff<<endl;
    //cout<<"reweight  = "<<reweight<<endl;

    if((charge=isD0Pair(kp))!=0 )
    {
      float d0Fill[20] = {0}; 
      double etacut[3] = {0,0.2,0.5};
      d0Fill[0] = d0Phi; 
      bool fillSB[8];
      vector<float> hadronPhi[3];
      vector<float>::iterator it;
      for(int eta_i=0;eta_i<3;eta_i++)
      {
        if(charge<0)
          getCorHadron(d0Eta,hadronPhi[eta_i],d0Daughter_u,d0Phi,etacut[eta_i]);
        else 
          getCorHadron(d0Eta,hadronPhi[eta_i],d0Daughter_l,d0Phi,etacut[eta_i]);
      }
      if(charge==-1)
        massPt->Fill(dMass,d0Pt,reweight*reweight_eff);
      if(charge>0)
        massPtLike->Fill(dMass,d0Pt,reweight*reweight_eff);
      fillSB[0] =  (charge>0)&& (dMass>(mean-1*sigma)) &&  (dMass<(mean+1*sigma));
      fillSB[1] =  (charge>0)&& (dMass>(mean-3*sigma)) &&  (dMass<(mean+3*sigma));
      fillSB[2] =  (charge>0) && (((dMass>(mean+4*sigma)) &&  (dMass<(mean+9*sigma))) ||((dMass>(mean-9*sigma)) &&  (dMass<(mean-4*sigma))));
      fillSB[4] = (charge==-1)&& (((dMass>(mean+5.5*sigma)) &&  (dMass<(mean+7.5*sigma))) ||((dMass>(mean-7.5*sigma)) &&  (dMass<(mean-5.5*sigma))));
      fillSB[5] = (charge==-1)&& (dMass>(mean-3*sigma)) &&  (dMass<(mean+3*sigma));
      fillSB[6] = (charge==-1)&& (((dMass>(mean+4*sigma)) &&  (dMass<(mean+9*sigma))) ||((dMass>(mean-9*sigma)) &&  (dMass<(mean-4*sigma))));
      fillSB[3] = fillSB[1] || fillSB[2];
      fillSB[7] = fillSB[1] || fillSB[2] || fillSB[6];

      if(charge>0)
      {
        likeV2Mass[fitindex][1]->Fill(dMass,cos(2*d0Phi)*weight_hadron);
        likeV2Mass[fitindex][2]->Fill(dMass,sin(2*d0Phi)*weight_hadron);
        massLike->Fill(dMass,d0Pt,reweight*reweight_eff);
      }
      if(charge==2 )
      {
        likeV2Mass2[fitindex][1]->Fill(dMass,cos(2*d0Phi));
        likeV2Mass2[fitindex][2]->Fill(dMass,sin(2*d0Phi));
        massLike2->Fill(dMass,d0Pt,reweight*reweight_eff);
      }
      if(charge==-1)
      {
        unlikeV2Mass[fitindex][1]->Fill(dMass,cos(2*d0Phi)*weight_hadron);
        unlikeV2Mass[fitindex][2]->Fill(dMass,sin(2*d0Phi)*weight_hadron);
        massUnlike->Fill(dMass,d0Pt,reweight*reweight_eff);
      }
      //getCorHadron(d0Eta,hadronPhi[0],index1,index2,d0Phi,0);
      for(it=hadronPhi[0].begin(); it!=hadronPhi[0].end(); it++)
      {
        double mHadronPhi = *it;
        double deltaPhi = 2*(d0Phi - mHadronPhi);
        double cos2Phi = cos(deltaPhi); 
        if(d0Pt>10) continue;

        if(charge>0 )
        {
          likeV2Mass[fitindex][0]->Fill(dMass,cos2Phi*weight_hadron,reweight*reweight_eff);
          likeV2Mass[fitindex][3]->Fill(dMass,cos(2*mHadronPhi)*weight_hadron);
          likeV2Mass[fitindex][4]->Fill(dMass,sin(2*mHadronPhi)*weight_hadron);
      
          if(fitindex==4) 
            checkPeak->Fill(dMass,cos2Phi,reweight*reweight_eff);
        }
        if(charge==2 )
        {
          likeV2Mass2[fitindex][0]->Fill(dMass,cos2Phi,reweight*reweight_eff);
          likeV2Mass2[fitindex][3]->Fill(dMass,cos(2*mHadronPhi));
          likeV2Mass2[fitindex][4]->Fill(dMass,sin(2*mHadronPhi));
        }
        if(charge==-1)
        {
          unlikeV2Mass[fitindex][0]->Fill(dMass,cos2Phi*weight_hadron,reweight*reweight_eff);
          unlikeV2Mass[fitindex][3]->Fill(dMass,cos(2*mHadronPhi)*weight_hadron);
          unlikeV2Mass[fitindex][4]->Fill(dMass,sin(2*mHadronPhi)*weight_hadron);
        }
        if(d0Pt>4 && d0Pt<5)
        {
          float checkFill[4] = {(float)deltaPhi,(float)dMass,(float)d0Pt,(float)charge};
          checkNew->Fill(checkFill);
        }
      
      }

      for(int i=0;i<3;i++)//different eta gap
      {
        for(int j=0;j<8;j++)
        {
          if(fillSB[j])
          {
            profV2[j][1][i]->Fill(d0Pt,cos(2*d0Phi));
            profV2[j][2][i]->Fill(d0Pt,sin(2*d0Phi));
          }
        }
        for(it=hadronPhi[i].begin(); it!=hadronPhi[i].end(); it++)
        {
          double mHadronPhi = *it;
          double deltaPhi = d0Phi - mHadronPhi;
          double cos2Phi = cos(2*deltaPhi); 
          for(int j=0;j<8;j++)
          {
            if(fillSB[j])
            {
              profV2[j][0][i]->Fill(d0Pt,cos2Phi,reweight*reweight_eff);
              v2Weight[j][i]->Fill(centrality,d0Pt,reweight*reweight_eff);
              profV2[j][3][i]->Fill(d0Pt,cos(2*mHadronPhi));
              profV2[j][4][i]->Fill(d0Pt,sin(2*mHadronPhi));
            }
          }
        }//Hadron loop
      }//different eta gap 
/*
      int index1 = kp->kaonIdx();
      int index2 = kp->pionIdx();
      vector<unsigned int>::iterator it1;
      vector<unsigned int>::iterator it2;
      it1 = find(d0Daughter_u.begin(),d0Daughter_u.end(),index1);
      it2 = find(d0Daughter_u.begin(),d0Daughter_u.end(),index2);
      bool testLikesign = (it1==d0Daughter_u.end() && it2==d0Daughter_u.end());
      if(!testLikesign && charge>0) continue;
      for(it=hadronPhi[0].begin(); it!=hadronPhi[0].end(); it++)
      {
        double mHadronPhi = *it;
        double deltaPhi = 2*(d0Phi - mHadronPhi);
        double cos2Phi = cos(deltaPhi); 
        if(kp->pt()>5 ||kp->pt()<4) continue;

        if(charge==1 )
        {
          likeV2Mass1[5][0]->Fill(kp->m(),cos2Phi,reweight);
          likeV2Mass1[5][3]->Fill(kp->m(),cos(mHadronPhi));
          likeV2Mass1[5][4]->Fill(kp->m(),sin(mHadronPhi));
          if(5==4) 
            checkPeak->Fill(kp->m(),cos2Phi,reweight);
        }
        if(charge==2 )
        {
          likeV2Mass2[5][0]->Fill(kp->m(),cos2Phi,reweight);
          likeV2Mass2[5][3]->Fill(kp->m(),cos(mHadronPhi));
          likeV2Mass2[5][4]->Fill(kp->m(),sin(mHadronPhi));
        }
        if(charge==-1)
        {
          unlikeV2Mass[5][0]->Fill(kp->m(),cos2Phi,reweight);
          unlikeV2Mass[5][3]->Fill(kp->m(),cos(mHadronPhi));
          unlikeV2Mass[5][4]->Fill(kp->m(),sin(mHadronPhi));
        }
      }

      for(int i=1;i<2;i++)//different eta gap
      {
        for(int j=0;j<8;j++)
        {
          if(fillSB[j])
          {
            profV2[j][1][i]->Fill(kp->pt(),cos(d0Phi));
            profV2[j][2][i]->Fill(kp->pt(),sin(d0Phi));
          }
        }
        for(it=hadronPhi[i].begin(); it!=hadronPhi[i].end(); it++)
        {
          double mHadronPhi = *it;
          double deltaPhi = d0Phi - mHadronPhi;
          double cos2Phi = cos(2*deltaPhi); 
          for(int j=0;j<8;j++)
          {
            if(fillSB[j])
            {
              profV2[j][0][i]->Fill(kp->pt(),cos2Phi,reweight);
              v2Weight[j][i]->Fill(centrality,kp->pt(),reweight);
              profV2[j][3][i]->Fill(kp->pt(),cos(mHadronPhi));
              profV2[j][4][i]->Fill(kp->pt(),sin(mHadronPhi));
            }
          }
        }//Hadron loop
      }//different eta gap 
*/


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
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0061 &&
      kp->pionDca() > 0.0110 && kp->kaonDca() > 0.0103 &&
      kp->dcaDaughters() < 0.0084 && kp->decayLength()>0.0145;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0049 &&
      kp->pionDca() > 0.0111 && kp->kaonDca() > 0.0091 &&
      kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0181;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0095 &&
      kp->dcaDaughters() < 0.0057 && kp->decayLength()>0.0212;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0081 && kp->kaonDca() > 0.0079 &&
      kp->dcaDaughters() < 0.0050 && kp->decayLength()>0.0247;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0040 &&
      kp->pionDca() > 0.0062 && kp->kaonDca() > 0.0058 &&
      kp->dcaDaughters() < 0.0060 && kp->decayLength()>0.0259;  
  }

  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}
/*
*/

int StPicoD0AnaMaker::isD0PairOld(StKaonPion const* const kp) const
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
  return trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->pMom().pseudoRapidity())<1.&&fabs(trk->nSigmaElectron())>3 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
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


bool StPicoD0AnaMaker::getCorHadron(float eta,vector<float> &hadronsPhi, vector<unsigned int> index1, float phi, float etaCut) 
{
  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* hadron = picoDst->track(i);
    if(!hadron)  continue; 
    if(hadron->pMom().perp()<0.2) continue;
    etaPhi_Hadron_all->Fill(hadron->pMom().phi(),hadron->pMom().pseudoRapidity());
    vector<unsigned int>::iterator it_index;
    it_index = find(index1.begin(),index1.end(),i);
    if(it_index!=index1.end())  continue;
    if(!isGoodHadron(hadron)) continue;
    float dEta = fabs(hadron->pMom().pseudoRapidity() - eta);
    float dPhi = (hadron->pMom().phi() - phi);
    if(etaCut<0.001)
    {
      dEtaDHadron->Fill(dEta);
      hEtaD->Fill(eta);
      hPhiD->Fill(phi);
      hEtaHadron->Fill(hadron->pMom().pseudoRapidity());
      hPhiHadron->Fill(hadron->pMom().phi());
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
  const double reweight = mGRefMultCorrUtil->getWeight();
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* hadron = picoDst->track(i);
    if(hadron->pMom().perp()<0.2) continue;
    if(!isGoodHadron(hadron)) continue;
    float etaHadron = hadron->pMom().pseudoRapidity();
    float phiHadron = hadron->pMom().phi();
    if(etaHadron<-0.35)
    {
      hadronFill[0]++;
      hadronFill[1] += sin(2 * phiHadron);
      hadronFill[2] += cos(2 * phiHadron);
    }			
    if(etaHadron>0.35)
    {
      hadronFill[3]++;
      hadronFill[4] += sin(2 * phiHadron);
      hadronFill[5] += cos(2 * phiHadron);
    }			
  }
  hadronFill[6] = centrality;
  hadronFill[7] = reweight;
  //mHadronTuple->Fill(hadronFill);
  if(hadronFill[0]==0 || hadronFill[3]==0)
     return false; 
  double temp = (hadronFill[1]*hadronFill[4]+hadronFill[2]*hadronFill[5]);
  hadronV2[0]->Fill(centrality,temp*reweight);
  hadronV2[1]->Fill(centrality,hadronFill[2]*reweight);
  hadronV2[2]->Fill(centrality,hadronFill[1]*reweight);
  hadronV2[3]->Fill(centrality,hadronFill[5]*reweight);
  hadronV2[4]->Fill(centrality,hadronFill[4]*reweight);
  hadronV2_sum[0]->Fill(centrality,hadronFill[0]*hadronFill[3]*reweight);
  hadronV2_sum[1]->Fill(centrality,hadronFill[0]*reweight);
  hadronV2_sum[2]->Fill(centrality,hadronFill[0]*reweight);
  hadronV2_sum[3]->Fill(centrality,hadronFill[3]*reweight);
  hadronV2_sum[4]->Fill(centrality,hadronFill[3]*reweight);
//    StPicoTrack const* hadron = picoDst->track(i);
//  hadronV2_excl[0][centrality]->Fill(hadron->pMom().perp(),temp*reweight);
//  hadronV2_excl[1][centrality]->Fill(hadron->pMom().perp(),hadronFill[2]*reweight);
//  hadronV2_excl[2][centrality]->Fill(hadron->pMom().perp(),hadronFill[1]*reweight);
//  hadronV2_excl[3][centrality]->Fill(hadron->pMom().perp(),hadronFill[5]*reweight);
//  hadronV2_excl[4][centrality]->Fill(hadron->pMom().perp(),hadronFill[4]*reweight);
//  hadronV2_excl_sum[0][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*hadronFill[3]*reweight);
//  hadronV2_excl_sum[1][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*reweight);
//  hadronV2_excl_sum[2][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*reweight);
//  hadronV2_excl_sum[3][centrality]->Fill(hadron->pMom().perp(),hadronFill[3]*reweight);
//  hadronV2_excl_sum[4][centrality]->Fill(hadron->pMom().perp(),hadronFill[3]*reweight);
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

























