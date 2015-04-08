#ifndef StMyAnalysisMaker_h
#define StMyAnalysisMaker_h

#include "StMaker.h"
#include "TNtuple.h"


//StKFVertexMaker includes
//#include "TObjArray.h"
#include "StEnumerations.h"
#include "StThreeVectorF.hh"
//
#include "TSpectrum.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"
#include "TCanvas.h"
#include "TH1K.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StPhysicalHelixD.hh"
class StPrimaryVertex; 
class StEvent;
class StDcaGeometry; 
class KFParticle; 
class StKFVerticesCollection; 
class TSpectrum; 
class StAnneling;
//// 

class StPicoDst;
class StPicoDstMaker;
class TString;
class TH1F;
class TH2F;

class StMyAnalysisMaker : public StMaker {
  public:
    StMyAnalysisMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
    virtual ~StMyAnalysisMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();

    //void    AppKFParticle(StPicoEvent *evt, StPicoTrack *t1, StPicoTrack *t2);
    void    DeclareHistograms();
    void    WriteHistograms();
		//StKFVertexMaker public
    //int Reco(vector<StPhysicalHelixD>& v_kaon, vector<StPhysicalHelixD>& v_pion,vector<StPhysicalHelixD>& v_soft_pion,vector<vector<double> >& c_kaon,vector<vector<double> >& c_pion,vector<vector<double> >& c_soft_pion, StThreeVectorF& pVtx,int opt,vector<int>& daughter);
    int Reco(StPhysicalHelixD v_kaon[10000],StPhysicalHelixD v_pion[10000],StPhysicalHelixD v_soft_pion[10000],double c_kaon[10000][10],double c_pion[10000][10],double c_soft_pion[10000][10],StThreeVectorF& pVtx,int opt,vector<int>& daughter,int count_k,int count_p);
    int primaryVertexRefit(StThreeVectorF *, vector<int>& daughter);
    void clear();
    void Fit();
    TH1F *VtxM() {return fVtxM;}
    void SetZwindow(Double_t z = 2) {fzWindow = z;}
    void SetDefaultTempLog(Double_t tLog = 2) {fTempLog = tLog;}
    static Double_t AnnelingFcn(Double_t TInv=1);
  //  Double_t AnnelingFcn(Double_t TInv=1);
    TH1 *Vtx() {return fVtx;}
    StKFVerticesCollection* Vertices() {return fcVertices;}
    TObjArray &Particles() {return *fParticles;}
    KFParticle *AddTrackAt(const StDcaGeometry *dca,Int_t kg);
//    KFParticle *AddTrackAt(StDcaGeometry *dca,Int_t kg);
    void calculateRank(StPrimaryVertex *primV);
    void SetCanvas(TCanvas *c1) {fc1 = c1;}
    TCanvas *Canvas() {return fc1;}
    TH1F *GetVtxs(Int_t pass = 0) {return fVtxs[pass];}
    TH1K *GetVtxKs(Int_t pass = 0) {return fVtxKs[pass];}
    TH1F *GetVtxM() {return fVtxM;}
/*
*/
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    
    TString    mOutName;
    
    TNtuple*   mTracktuple;
    TNtuple*   mEventtuple;
    TNtuple*   mDtuple;
    TNtuple*   mDstartuple;
    TH1D* mCosTheta;
    TH1D* mCosThetastar;
    TH1D* mPxlCosTheta;
    TH1D* mPxlCosThetastar;
    TH1D* mdRprm;
    TH1D* mdRKpi;
    TH1D* mPxldRprm;
    TH1D* mPxldRKpi;
    TH1D*  mDmass_unlike;
    TH1D*  mDmass_like;
    TH1D*  mDsmass_unlike;
    TH1D*  mDsmass_like;
    TH1D*  mDmasscut_unlike;
    TH1D*  mDmasscut_like;
    TH1D*  mDmasstest_unlike;
    TH1D*  mDmasstest_like;
    TH1D*  mDsmasscut_unlike;
    TH1D*  mDsmasscut_like;

    TH1D*  mDmasstrigger_unlike;
    TH1D*  mDmasstrigger_like;
    TH1D*  mDmasstriggertest_unlike;
    TH1D*  mDmasstriggertest_like;
    TH1D*  mDmasstriggercut_unlike;
    TH1D*  mDmasstriggercut_like;

    TH1D*  mDsmasstrigger_unlike;
    TH1D*  mDsmasstrigger_like;
    TH1D*  mDsmasstriggercut_unlike;
    TH1D*  mDsmasstriggercut_like;
    TH1D*  mMult;
    StDcaGeometry *dca;
    StDcaGeometry *dcaG;
//
//	//StKFVertexMaker private
	TObjArray *fParticles; // KF particles
	Int_t fNzBins;
	Int_t fNPasses;
	TSpectrum *fSpectrum;
	Double_t fzWindow;
	TH1F  *fVtxM;
	StKFVerticesCollection **fVerticesPass;
	static StKFVerticesCollection *fcVertices;  // current vertex collection
	//StKFVerticesCollection *fcVertices;  // current vertex collection
	Double_t fTempLog;
	ROOT::Math::GSLMinimizer1D *fminBrent;
	ROOT::Math::Functor1D      *func;
	TH1F *fVtxs[2];
	TH1  *fVtx;
	TH1K *fVtxKs[2];
	Bool_t mBeamLine;
	StPrimaryVertexOrder     mVertexOrderMethod; // will default to 0 i.e. orderByNumberOfDaughters
	TCanvas                 *fc1;
	StPrimaryVertex *primV;
//
//                   
    ClassDef(StMyAnalysisMaker, 1)
};

#endif
