#include "StMyAnalysisMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoV0.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StPicoDstMaker/StPicoMtdPidTraits.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
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
ClassImp(StMyAnalysisMaker)
StKFVerticesCollection *StMyAnalysisMaker::fcVertices = 0;
//-----------------------------------------------------------------------------
StMyAnalysisMaker::StMyAnalysisMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mOutName = outName;
}

//----------------------------------------------------------------------------- 
StMyAnalysisMaker::~StMyAnalysisMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Init() {
  DeclareHistograms();

  dcaG = new StDcaGeometry();
  dca = new StDcaGeometry();
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
//  cout<<"checkpoint1"<<endl;
  func = new ROOT::Math::Functor1D(&StMyAnalysisMaker::AnnelingFcn);
  //func = new ROOT::Math::Functor1D(&StKFVertexMaker::AnnelingFcn);
  fminBrent = new ROOT::Math::GSLMinimizer1D();
  fVerticesPass = new StKFVerticesCollection *[fNPasses+1];
  memset (fVerticesPass, 0, (fNPasses+1)*sizeof(StKFVerticesCollection *));
  fParticles = new TObjArray();
  fParticles->SetOwner(kTRUE);
  mVertexOrderMethod = orderByRanking; // change ordering by ranking

  primV  = new StPrimaryVertex;

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Finish() {
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(),"RECREATE");
    fout->cd();
    WriteHistograms();
    fout->Close();
  }
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
  SafeDelete(dcaG);
	
  return kStOK;
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::DeclareHistograms() {
  mTracktuple = new TNtuple("mTracktuple","mTracktuple","chi2:chi2prob:gMomX:gMomY:gMomZ:pMomX:pMomY:pMomZ:originX:originY:originZ:dca:charge:nHitsFit:nHitsMax:nHitsDedx:dEdx:nSigmaPion:nSigmaKaon:nSigmaElectron:nSigmaProton:btofBeta:mtdBeta:mass2:pxl1hits:pxl2hits:isthits:triggerWord");
  mEventtuple = new TNtuple("mEventtuple","mEventtuple","bField:prmX:prmY:prmZ:refitprmX:refitprmY:refitprmZ:refMult:numOfGtrc:ranking:isMB:refitDX:refitDY:refitDZ");
  mDstartuple = new TNtuple("mDstartuple","mDstartuple","mass_d0:cos_prm:dR_prm:dR_kp:sign:kdca:pdca:softdca:triggerWord:thetastar:mass_kpp");
  mDtuple = new TNtuple("mDtuple","mDtuple","mass_d0:cos_prm:dR_prm:dR_kp:pt_d0:sign:kmass:kdca:kpxl2:kist:pmass:pdca:ppxl2:pist:isMinB:cos_prm_refit:dR_prm_refit:refMult:kdca_refit:pdca_refit");
  mCosTheta = new TH1D("mCosTheta","mCosTheta",1000,-1,1);
  mCosThetastar = new TH1D("mCosThetastar","mCosThetastar",1000,-1,1);
  mdRprm = new TH1D("mdRprm","mdRprm",1000,0,1);
  mdRKpi = new TH1D("mdRKpi","mdRKpi",1000,0,0.1);
  mPxlCosTheta = new TH1D("mPxlCosTheta","mPxlCosTheta",1000,-1,1);
  mPxlCosThetastar = new TH1D("mPxlCosThetastar","mPxlCosThetastar",1000,-1,1);
  mPxldRprm = new TH1D("mPxldRprm","mPxldRprm",1000,0,1);
  mPxldRKpi = new TH1D("mPxldRKpi","mPxldRKpi",1000,0,0.1);
  mDmass_unlike = new TH1D("mDmass_unlike","",2500,0.5,3);
  mDmass_like = new TH1D("mDmass_like","",2500,0.5,3);
  mDsmass_unlike = new TH1D("mDsmass_unlike","",200,1.4,1.6);
  mDsmass_like = new TH1D("mDsmass_like","",200,1.4,1.6);
  mDmasscut_unlike = new TH1D("mDmasscut_unlike","",2500,0.5,3);
  mDmasscut_like = new TH1D("mDmasscut_like","",2500,0.5,3);
  mDmasstest_unlike = new TH1D("mDmasstest_unlike","",2500,0.5,3);
  mDmasstest_like = new TH1D("mDmasstest_like","",2500,0.5,3);
  mDsmasscut_unlike = new TH1D("mDsmasscut_unlike","",200,1.4,1.6);
  mDsmasscut_like = new TH1D("mDsmasscut_like","",200,1.4,1.6);

  mDmasstrigger_unlike = new TH1D("mDmasstrigger_unlike","",2500,0.5,3);
  mDmasstrigger_like = new TH1D("mDmasstrigger_like","",2500,0.5,3);
  mDsmasstrigger_unlike = new TH1D("mDsmasstrigger_unlike","",200,1.4,1.6);
  mDsmasstrigger_like = new TH1D("mDsmasstrigger_like","",200,1.4,1.6);
  mDmasstriggercut_unlike = new TH1D("mDmasstriggercut_unlike","",2500,0.5,3);
  mDmasstriggercut_like = new TH1D("mDmasstriggercut_like","",2500,0.5,3);
  mDsmasstriggercut_unlike = new TH1D("mDsmasstriggercut_unlike","",200,1.4,1.6);
  mDsmasstriggercut_like = new TH1D("mDsmasstriggercut_like","",200,1.4,1.6);

}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::WriteHistograms() {
  mTracktuple->Write();
  mEventtuple->Write();
  mDtuple->Write();
  mDstartuple->Write();
  mCosTheta->Write();
  mCosThetastar->Write();
  mPxlCosThetastar->Write();
  mdRprm->Write();
  mdRKpi->Write();
  mPxlCosTheta->Write();
  mPxldRprm->Write();
  mPxldRKpi->Write();
  mDmass_unlike->Write();
  mDmass_like->Write();
  mDmasscut_unlike->Write();
  mDmasscut_like->Write();
  mDmasstest_unlike->Write();
  mDmasstest_like->Write();
  mDsmass_unlike->Write();
  mDsmass_like->Write();
  mDsmasscut_unlike->Write();
  mDsmasscut_like->Write();
  mDmasstrigger_unlike->Write();
  mDmasstrigger_like->Write();
  mDmasstriggercut_unlike->Write();
  mDmasstriggercut_like->Write();
  mDsmasstrigger_unlike->Write();
  mDsmasstrigger_like->Write();
  mDsmasstriggercut_unlike->Write();
  mDsmasstriggercut_like->Write();
}

//----------------------------------------------------------------------------- 
void StMyAnalysisMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Make() {
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }
//Added by Leon/////////////////////Vector has problem
  float track_fill[100];
  float event_fill[100];
  for(int i=0;i<100;i++)
    event_fill[i] = 0;
  StPhysicalHelixD v_k[10000];
  StPhysicalHelixD v_p[10000];
  StPhysicalHelixD v_soft_p[10000];
  double c_k[10000][10];
  double c_p[10000][10];
  double c_soft_p[10000][10];
  int count_k = 0;
  int count_p = 0;
  int count_soft_p = 0;
  
//////////////////////////////////
  vector<int> daughter;
  daughter.clear(); 
  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  if(!(event->isMinBias())||event->refMult()>100)
  { 
    LOG_WARN << " Not Min Bias! Skip! " << endm;
    return kStWarn;
  }
  StThreeVectorF pVtx(-999.,-999.,-999.);
  StThreeVectorF testVertex(-999.,-999.,-999.);
  if(event) {
     pVtx = event->primaryVertex();
 //    primaryVertexRefit(&testVertex,daughter);
  }
  testVertex = pVtx;

//Added by Leon////////////Event Information
//  double field = event->bField() * pow(10,-14);
  event_fill[0] = event->bField();
  event_fill[1] = pVtx.x();
  event_fill[2] = pVtx.y();
  event_fill[3] = pVtx.z();
  event_fill[7] = event->refMult();
  event_fill[8] = event->numberOfGlobalTracks();
  event_fill[9] = event->ranking();
  event_fill[10] = event->isMinBias();
//////////////////////////
  int nTracks = mPicoDst->numberOfTracks();
  //cout<<"Events loop"<<endl;
  for(int i=0;i<nTracks;i++) 
  {
    StPicoTrack *t = (StPicoTrack*)mPicoDst->track(i);
    if(!t) continue;
//    StDcaGeometry *dcaG = new StDcaGeometry();
    dcaG->set(t->params(),t->errMatrix());
    StPhysicalHelixD helix = dcaG->helix();
    double dca = helix.geometricSignedDistance(pVtx);

    // global helix
    StThreeVectorF origin;// = t->origin();
    StThreeVectorF mom = helix.momentum(event->bField()*kilogauss);
    StThreeVectorF pmom = t->pMom();

    int charge = t->charge();
//    StPhysicalHelixD gHelix(pmom, origin, event->bField()*kilogauss, charge);
    //Added by Leon////////////////Track Information
    track_fill[0] = t->chi2();
    track_fill[1] = 0;//t->chi2prob();
    track_fill[2] = mom.x();
    track_fill[3] = mom.y();
    track_fill[4] = mom.z();
    track_fill[5] = pmom.x();
    track_fill[6] = pmom.y();
    track_fill[7] = pmom.z();
    track_fill[8] = origin.x();
    track_fill[9] = origin.y();
    track_fill[10] = origin.z();
    track_fill[11] = dca;//t->dca();
    track_fill[12] = charge;
    track_fill[13] = t->nHitsFit();
    track_fill[14] = 0;//t->nHitsMax();
    track_fill[15] = t->nHitsDedx();
    track_fill[16] = t->dEdx();
    track_fill[17] = t->nSigmaPion();
    track_fill[18] = t->nSigmaKaon();
    track_fill[19] = t->nSigmaElectron();
    track_fill[20] = t->nSigmaProton();
// To get Beta/////////////
    int btofIndex = t->bTofPidTraitsIndex();
    if(btofIndex<0)
	track_fill[21] = -1;
    else
    {
	    StPicoBTofPidTraits *btofpid = (StPicoBTofPidTraits*)mPicoDst->btofPidTraits(btofIndex);
	    track_fill[21] = btofpid->btofBeta();
    }

    int mtdIndex = t->mtdPidTraitsIndex();
    if(mtdIndex<0)
	track_fill[22] = -1;
    else
    {
	    StPicoMtdPidTraits *mtdpid = (StPicoMtdPidTraits*)mPicoDst->mtdPidTraits(mtdIndex);
	    track_fill[22] = mtdpid->beta();
    }
    if(track_fill[21]>0.00000001)
      track_fill[23] = mom.mag2() * (1./track_fill[21]/track_fill[21] -1);
    else
      track_fill[23] = -999.;
/////////////////////////////////
//Hits store//////////////////////
    if((t->nHitsMapHFT()>>0 & 0x1)==1) track_fill[24]=1; 
    else track_fill[24] = 0;
    if((t->nHitsMapHFT()>>1 & 0x3)==3) track_fill[25]=2; 
    else if((t->nHitsMapHFT()>>1 & 0x3)<3&&(t->nHitsMapHFT()>>1 & 0x3)>0) track_fill[25]=1; 
      else track_fill[25] = 0;
    if((t->nHitsMapHFT()>>3 & 0x3)==3) track_fill[26]=2; 
    else if((t->nHitsMapHFT()>>3 & 0x3)<3 && (t->nHitsMapHFT()>>3 & 0x3)>0) track_fill[26]=1; 
      else track_fill[26] = 0;
/////////////////////////////////////

    track_fill[27] = event->isMinBias();
   
//
//
//Store K and Pi information//////////////////////////////////////
    if(t->nHitsFit()<20) continue;
    if(t->nHitsDedx()<15) continue;
    if(fabs(pVtx.z())>6) continue;
    StLorentzVectorD trk;
    trk.setPx(mom.x());
    trk.setPy(mom.y());
    trk.setPz(mom.z());
    trk.setE(mom.mag()/track_fill[21]);
    if((trk.pseudoRapidity())>1.) continue;
    if(!(t->nHitsMapHFT()>>3 & 0x3)) continue; //No IST hits;
    if(!(t->nHitsMapHFT()>>1 & 0x3)) continue; //No PXL2 hits;
    if(!(t->nHitsMapHFT()>>0 & 0x1)) continue; //No PXL1 hits;
////////Define parameters /////
//    vector<double> dPara;
//    dPara.push_back(charge);
//    dPara.push_back(t->id());
//    dPara.push_back(dca);
//    dPara.push_back(track_fill[25]);
//    dPara.push_back(track_fill[26]);
//    dPara.push_back(event->isMinBias());
////////////////////////////////////////////
 
////////////////////////Track cut///////////////////
//     if(fabs(t->nSigmaPion())<2 || (track_fill[23]>0.005 && track_fill[23]<0.04)){
//      v_soft_p.push_back(helix);
//      c_soft_p.push_back(dPara);
//    }
//////////////////Soft Pion for D0 star/////////////
    if(sqrt(mom.x()*mom.x()+mom.y()*mom.y())<1.0) continue;
/////////////////Add dca and pT track cuts///////////////
    if(fabs(t->nSigmaKaon())<2 || (track_fill[23]>0.2 && track_fill[23]<0.3))
    {
      v_k[count_k] = helix;
      c_k[count_k][0] = charge;
      c_k[count_k][1] = t->id();
      c_k[count_k][2] = dca;
      c_k[count_k][3] = track_fill[25];
      c_k[count_k][4] = track_fill[26];
      c_k[count_k][5] = event->isMinBias();
      count_k++;
    }
    if(fabs(t->nSigmaPion())<2 || (track_fill[23]>0.005 && track_fill[23]<0.04))
    {
      v_p[count_p] = helix;
      c_p[count_p][0] = charge;
      c_p[count_p][1] = t->id();
      c_p[count_p][2] = dca;
      c_p[count_p][3] = track_fill[25];
      c_p[count_p][4] = track_fill[26];
      c_p[count_p][5] = event->refMult();
      count_p++;
    }     
    
  /* 
    if(fabs(t->nSigmaKaon())<2 || (track_fill[23]>0.2 && track_fill[23]<0.3))
    {
      v_k.push_back(helix);
      c_k.push_back(dPara);
    }

    if(fabs(t->nSigmaPion())<2 || (track_fill[23]>0.005 && track_fill[23]<0.04))
    {
      v_p.push_back(helix);
      c_p.push_back(dPara);
    }
   
*/


//    mTracktuple->Fill(track_fill);
    

//    delete dcaG;
  }
  primaryVertexRefit(&testVertex,daughter);
  event_fill[4] = testVertex.x();
  event_fill[5] = testVertex.y();
  event_fill[6] = testVertex.z();
  int option = 1;//1:store d daughers; 2:store D information; 3:store Ds information
  Reco(v_k,v_p,v_soft_p,c_k,c_p,c_soft_p,pVtx,option,daughter,count_k,count_p);
  primaryVertexRefit(&testVertex,daughter);
  event_fill[11] = testVertex.x();
  event_fill[12] = testVertex.y();
  event_fill[13] = testVertex.z();
//  mEventtuple->Fill(event_fill);
  option = 2;
  Reco(v_k,v_p,v_soft_p,c_k,c_p,c_soft_p,testVertex,option,daughter,count_k,count_p);
//  cout<<"checkpoint4 = "<<testVertex<<endl;
//  cout<<"checkpoint5 = "<<pVtx<<endl;
//
//  vector<StPhysicalHelixD>(v_k).swap(v_k);
//  vector<StPhysicalHelixD>(v_p).swap(v_p);
//  vector<StPhysicalHelixD>(v_soft_p).swap(v_soft_p);
//  vector<vector<double> >(c_k).swap(c_k);
//  vector<vector<double> >(c_p).swap(c_p);
//  vector<vector<double> >(c_soft_p).swap(c_soft_p);
//    vector<double>(dPara).swap(dPara);
  return kStOK;
}
/*
*/
int StMyAnalysisMaker::primaryVertexRefit(StThreeVectorF *mRefitVertex, vector<int>& daughter) {
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
//  cout<<"checkpoint1"<<endl;
  func = new ROOT::Math::Functor1D(&StMyAnalysisMaker::AnnelingFcn);
  //func = new ROOT::Math::Functor1D(&StKFVertexMaker::AnnelingFcn);
  fminBrent = new ROOT::Math::GSLMinimizer1D();
  fVerticesPass = new StKFVerticesCollection *[fNPasses+1];
  memset (fVerticesPass, 0, (fNPasses+1)*sizeof(StKFVerticesCollection *));
  fParticles = new TObjArray();
  fParticles->SetOwner(kTRUE);
  mVertexOrderMethod = orderByRanking; // change ordering by ranking
*/

  //StKFVertexMaker Make part
  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  Double_t bField = event->bField();
  KFParticle::SetField(bField);
//  if (mBeamLine) {
//     St_vertexSeedC* vSeed = St_vertexSeedC::instance();
//     Double_t x0   = vSeed->x0()  ; Double_t err_x0   = vSeed->err_x0();  
//     Double_t y0   = vSeed->y0()  ; Double_t err_y0   = vSeed->err_y0();
//     Double_t dxdz = vSeed->dxdz(); Double_t err_dxdz = vSeed->err_dxdz();
//     Double_t dydz = vSeed->dydz(); Double_t err_dydz = vSeed->err_dydz();
//     Double_t weight = vSeed->weight();
//     if (err_x0 < 0.010) err_x0 = 0.010;
//     if (err_y0 < 0.010) err_y0 = 0.010;
//     static Bool_t firstTime = kTRUE;
//     if (firstTime) {
//	firstTime = kFALSE;
//	LOG_INFO << "BeamLine Constraint: weight =  " << weight << endm;
//	LOG_INFO << "x(z) = (" << x0 << " +- " << err_x0 << ") + (" << dxdz << " +- " << err_dxdz << ") * z" << endm;
//	LOG_INFO << "y(z) = (" << y0 << " +- " << err_y0 << ") + (" << dydz << " +- " << err_dydz << ") * z" << endm;
//     }
//     static Double_t pZ = 1000;
//     static MTrack track;
//     Double_t xyzP[6] = {     x0,      y0, 0.,
//			      pZ*dxdz, pZ*dydz, pZ};
//     Double_t CovXyzP[21] = {
//	      err_x0*err_x0,
//	      0            ,err_y0*err_y0,
//	      0            ,0              , 0,
//	      0            ,0              , 0, (err_dxdz*pZ)*(err_dxdz*pZ),
//	      0            ,0              , 0,                             0, (err_dydz*pZ)*(err_dydz*pZ)
//     };
//     track.SetParameters(xyzP);
//     track.SetCovarianceMatrix(CovXyzP);
//     track.SetNDF(1);
//     track.SetID(0);
//     track.SetCharge(1);
//     KFParticle *beam = new KFParticle(track, 321);
//     fParticles->AddAt(beam, 0);
//  }
//  StSPtrVecTrackNode& ptrackNode = pEvent->trackNodes();
//  UInt_t nTracks = ptrackNode.size();
//  StTrackNode *node=0;
  Int_t NGoodGlobals = 0;
//  map<Int_t,StTrackNode*> TrackNodeMap;
  int nTracks = mPicoDst->numberOfTracks();
  for (int i=0; i < nTracks; i++) {
//     node = ptrackNode[i]; 
//     if (!node) continue;
//     StGlobalTrack  *gTrack = static_cast<StGlobalTrack *>(node->track(global));
     StPicoTrack *gTrack = (StPicoTrack*)mPicoDst->track(i);
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
     if (Debug() > 1) {
	//cout << Form("particle: %4i/%4i ",NGoodGlobals,kg) << *particle << endl;
//	cout << "Add to map[" << kg << "] node = " << TrackNodeMap[kg] << endl;
     }
     NGoodGlobals++;
  //   cout << "Good Globals: " << i << endl;
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
      //V->PrintW();
      // Store vertex
//      StPrimaryVertex *primV  = new StPrimaryVertex;
      StThreeVectorF XVertex(&V->Vertex().X());
      //primV->Reset();
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
////////////////Free memory/////////////////////
  clear();
  return 1;

}


//StKFVertexMaker functions_______________________________________________________
void StMyAnalysisMaker::calculateRank(StPrimaryVertex *primV) {    
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
KFParticle *StMyAnalysisMaker::AddTrackAt(const StDcaGeometry *dca, Int_t kg) {
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
void StMyAnalysisMaker::Fit() {
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
Double_t StMyAnalysisMaker::AnnelingFcn(Double_t TInv) {
  if (! fcVertices) return 0;
  Double_t Temperature = 1./TInv;
  StAnneling::SetTemperature(Temperature);
  Double_t Chi2 =  fcVertices->Fit();
  if (StKFVertex::Debug()) 
    LOG_INFO << "StKFVertexMaker::AnnelingFcn\tTemperature = " << Temperature << " Chi2 = " << Chi2 << endm;
  return Chi2;
}
//________________________________________________________________________________
int StMyAnalysisMaker::Reco(StPhysicalHelixD v_kaon[10000],StPhysicalHelixD v_pion[10000],StPhysicalHelixD v_soft_pion[10000],double c_kaon[10000][10],double c_pion[10000][10],double c_soft_pion[10000][10],StThreeVectorF& pVtx,int opt,vector<int>& daughter,int count_k,int count_p)
//int StMyAnalysisMaker::Reco(vector<StPhysicalHelixD>& v_kaon, vector<StPhysicalHelixD>& v_pion,vector<StPhysicalHelixD>& v_soft_pion,vector<vector<double> >& c_kaon,vector<vector<double> >& c_pion,vector<vector<double> >& c_soft_pion,StThreeVectorF& pVtx,int opt,vector<int>& daughter)
{


  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  StThreeVectorF Vtx(-999,-999,-999);
  Vtx = event->primaryVertex();
/*
  vector<StPhysicalHelixD> v_pion;
  vector<StPhysicalHelixD> v_soft_pion;
  vector<vector<double> > c_kaon;
  vector<vector<double> > c_pion;
  vector<vector<double> > c_soft_pion;
  int opt;
  vector<int> daughter;
*/
//////////////Loop for K and Pi pairs//////////////////
//  for(vector<StPhysicalHelixD>::size_type i=0;i<v_kaon.size();i++){
//    for(vector<StPhysicalHelixD>::size_type j=0;j<v_pion.size();j++){
  for(int i=0;i<count_k;i++){
    for(int j=0;j<count_p;j++){
      float D_fill[100];
      float Dstar_fill[100];
      if(c_kaon[i][1] == c_pion[j][1]) continue;
      pair<double,double>pathLength = v_kaon[i].pathLengths(v_pion[j]);
      StThreeVectorD momDCAKaon = v_kaon[i].momentumAt(pathLength.first, event->bField()*kilogauss);
      StThreeVectorD momDCAPion = v_pion[j].momentumAt(pathLength.second, event->bField()*kilogauss);
      StThreeVectorD posDCAKaon = v_kaon[i].at(pathLength.first);
      StThreeVectorD posDCAPion = v_pion[j].at(pathLength.second);

      StLorentzVectorD LV1, LV2;
      LV1.setPx(momDCAKaon.x());
      LV1.setPy(momDCAKaon.y());
      LV1.setPz(momDCAKaon.z());
      LV1.setE(momDCAKaon.massHypothesis(0.493677));

      LV2.setPx(momDCAPion.x());
      LV2.setPy(momDCAPion.y());
      LV2.setPz(momDCAPion.z());
      LV2.setE(momDCAPion.massHypothesis(0.13957));
      StLorentzVectorD kstar = LV1;
      StLorentzVectorD rc_d0 = LV1+LV2;
      kstar.boost(-rc_d0);
      StThreeVectorF rc_d0_mom(rc_d0.x(),rc_d0.y(),rc_d0.z());
      float rc_d0_x = posDCAKaon.x() + (posDCAPion.x()-posDCAKaon.x())/2;
      float rc_d0_y = posDCAKaon.y() + (posDCAPion.y()-posDCAKaon.y())/2;
      float rc_d0_z = posDCAKaon.z() + (posDCAPion.z()-posDCAKaon.z())/2;
      StThreeVectorF rc_d0_pos(rc_d0_x,rc_d0_y,rc_d0_z);

      double DCA_pi_K = (posDCAPion-posDCAKaon).mag();//Distance of K and π
      double theta = rc_d0_mom.angle( (rc_d0_pos-Vtx) );
      double DCA_prm = (rc_d0_pos-Vtx).mag();//distance of D0 and primary vertex
      double theta_refit = rc_d0_mom.angle( (rc_d0_pos-pVtx) );
      double DCA_prm_refit = (rc_d0_pos-pVtx).mag();//distance of D0 and primary vertex
      //double thetastar = rc_d0.vect().angle(kstar.vect());
      double Dmass = (LV1+LV2).m();
      double dca_k = v_kaon[i].geometricSignedDistance(pVtx);
      double dca_p = v_pion[j].geometricSignedDistance(pVtx);
      D_fill[0] = Dmass;
      D_fill[1] = cos(theta);
      D_fill[2] = DCA_prm;
      D_fill[3] = DCA_pi_K;
      D_fill[4] = (LV1+LV2).perp();
      D_fill[5] = c_kaon[i][0] * c_pion[j][0];
      D_fill[6] = c_kaon[i][1];
      D_fill[7] = c_kaon[i][2];
      D_fill[8] = c_kaon[i][3];
      D_fill[9] = c_kaon[i][4];
      D_fill[10] = c_pion[j][1];
      D_fill[11] = c_pion[j][2];
      D_fill[12] = c_pion[j][3];
      D_fill[13] = c_pion[j][4];
      D_fill[14] = c_kaon[i][5];
      D_fill[15] = cos(theta_refit);
      D_fill[16] = DCA_prm_refit;
      D_fill[17] = c_pion[i][5];
      D_fill[18] = dca_k;
      D_fill[19] = dca_p;
     
//  To store the index of D0 daughters////
     // if(Dmass>1.8&&Dmass<1.92&&DCA_prm>0.01 && D_fill[5]==-1 && DCA_pi_K<0.005 && cos(theta)>0.995)
      if(Dmass>1.8&&Dmass<1.92)//&&DCA_prm>0.01 && D_fill[5]==-1 && DCA_pi_K<0.005 && cos(theta)>0.995)
      {
        daughter.push_back(c_kaon[i][1]);
        daughter.push_back(c_pion[j][1]);
      }
      if(dca_k>0.008&&dca_p>0.008&& cos(theta_refit)>0.995 && DCA_pi_K<0.005 && opt==1)
      {
        if(D_fill[5]==-1)
          mDmasstest_unlike->Fill(Dmass);
        if(D_fill[5]==1)
          mDmasstest_like->Fill(Dmass);
      }
      if(opt==1) continue;//To determine the D0 daughters
///////////////////////////////////////////
//    
//      mCosTheta->Fill(cos(theta));
//      mCosThetastar->Fill(D_fill[15]);
//      mdRprm->Fill(DCA_prm);
//      mdRKpi->Fill(D_fill[3]);
//      mPxlCosTheta->Fill(cos(theta_refit));
//      mPxldRprm->Fill(DCA_prm_refit);
      if(c_kaon[i][2]>0.008 && c_pion[j][2]>0.008 && cos(theta)>0.995 && DCA_pi_K<0.005)
      {
        if(D_fill[5]==-1)
          mDmass_unlike->Fill(Dmass);
        if(D_fill[5]==1)
          mDmass_like->Fill(Dmass);
      }
      if(dca_k>0.008&&dca_p>0.008&& cos(theta_refit)>0.995 && DCA_pi_K<0.005)
      {
        if(D_fill[5]==-1)
          mDmasscut_unlike->Fill(Dmass);
        if(D_fill[5]==1)
          mDmasscut_like->Fill(Dmass);
      }
      //if(DCA_prm<0.005) continue;
      //if(cos(theta)<0.95) continue;
//      mDtuple->Fill(D_fill);
     
///Reconstruct Dstar/////////
      if(opt==2) continue;
      if(Dmass<1.8||Dmass>1.92) continue;
/*
      for(vector<StPhysicalHelixD>::size_type k=0;k<v_soft_pion.size();k++){
        if(c_soft_pion[k][1] == c_pion[j][1]|| c_soft_pion[k][1] == c_kaon[i][1]) continue;
        StLorentzVectorD LV3;
        //StThreeVectorD softPion = v_pion[k].momentumAt(pVtx, event->bField()*kilogauss);
        StThreeVectorD softPion = v_soft_pion[k].momentum(event->bField()*kilogauss);
        LV3.setPx(softPion.x());
        LV3.setPy(softPion.y());
        LV3.setPz(softPion.z());
        LV3.setE(softPion.massHypothesis(0.13957));

        StLorentzVectorD rc_dstar = LV1+LV2+LV3;
	D_fill[16] = rc_dstar.m();
        if(D_fill[5]==1)
  	  mDsmass_like->Fill(rc_dstar.m());
        if(D_fill[5]==-1)
    	  mDsmass_unlike->Fill(rc_dstar.m());
        if(DCA_prm>0.01 && D_fill[5]==1 && DCA_pi_K<0.005 && cos(theta)>0.995 && D_fill[15]<0.8)
        //if(D_fill[5]==1 && DCA_pi_K<0.005 && cos(theta)>0.995 && D_fill[15]<0.8)
  	  mDsmasscut_like->Fill(rc_dstar.m());
        if(DCA_prm>0.01 && D_fill[5]==-1 && DCA_pi_K<0.005 && cos(theta)>0.995 && D_fill[15]<0.8)
        //if(D_fill[5]==-1 && DCA_pi_K<0.005 && cos(theta)>0.995 && D_fill[15]<0.8)
  	  mDsmasscut_unlike->Fill(rc_dstar.m());

        if(event->isMinBias())
	{
	  if(D_fill[5]==1)
	    mDsmasstrigger_like->Fill(rc_dstar.m());
	  if(D_fill[5]==-1)
	    mDsmasstrigger_unlike->Fill(rc_dstar.m());
          if(D_fill[5]==1 && DCA_pi_K<0.005 && cos(theta)>0.995 && D_fill[15]<0.8)
	    mDsmasstriggercut_like->Fill(rc_dstar.m());
          if(D_fill[5]==-1 && DCA_pi_K<0.005 && cos(theta)>0.995 && D_fill[15]<0.8)

	    mDsmasstriggercut_unlike->Fill(rc_dstar.m());
	}
	Dstar_fill[0] = D_fill[0];//D0mass
	Dstar_fill[1] = D_fill[1];//cos_prm
	Dstar_fill[2] = D_fill[2];//dR_prm
	Dstar_fill[3] = D_fill[3];//dR_kp
	Dstar_fill[4] = D_fill[5];//sign
	Dstar_fill[5] = c_kaon[i][2];//dca Kaon
	Dstar_fill[6] = c_pion[j][2];//dca pion
	Dstar_fill[7] = c_soft_pion[k][2];//dca soft pion
	Dstar_fill[8] = event->isMinBias();//trigger
	Dstar_fill[9] = D_fill[15];//theta star
	Dstar_fill[10] = rc_dstar.m();//Dstar mass


        mDstartuple->Fill(Dstar_fill);
      }
*/

    }
  }
  return 1;
}

void StMyAnalysisMaker::clear() {
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

