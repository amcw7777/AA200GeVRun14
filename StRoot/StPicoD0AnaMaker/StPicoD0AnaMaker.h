#ifndef StPicoD0AnaMaker_h
#define StPicoD0AnaMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoD0Event
 *  simultaneously and do analysis. 
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
//
#include "StThreeVectorF.hh"
#include "TSpectrum.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"
#include "TCanvas.h"
#include "TH1K.h"
#include "TH2D.h"
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


class TString;
class TFile;
class TNtuple;
class StPicoD0Event;
class StKaonPion;
class StPicoDstMaker;
class StHFCuts;
class StPicoDst;


class StPicoD0AnaMaker : public StMaker
{
  public:
    StPicoD0AnaMaker(char const * name, char const * inputFilesList, 
        char const * outName,StPicoDstMaker* picoDstMaker);
    virtual ~StPicoD0AnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;

    void setHFCuts(StHFCuts* cuts);    
/*
*/
/////Refit public functions//
    double primaryVertexRefit(StThreeVectorF *, vector<int>& daughter);
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
    void Clear();
    TCanvas *Canvas() {return fc1;}
    TH1F *GetVtxs(Int_t pass = 0) {return fVtxs[pass];}
    TH1K *GetVtxKs(Int_t pass = 0) {return fVtxKs[pass];}
    TH1F *GetVtxM() {return fVtxM;}

  private:
    StPicoD0AnaMaker() {}
    void readNextEvent();

    bool isGoodPair(StKaonPion const*) const;

    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;

    //StPicoDstMaker *
    StPicoDst *picoDst;

    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    TChain* mChain;
    int mEventCounter;

    StHFCuts* mHFCuts;

    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
//	//StKFVertexMaker private
        TNtuple *mRefittuple;
        TH2D *timemult;
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
	//StPrimaryVertexOrder     mVertexOrderMethod; // will default to 0 i.e. orderByNumberOfDaughters
	TCanvas                 *fc1;
//

    ClassDef(StPicoD0AnaMaker, 1)
};

inline int StPicoD0AnaMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoD0AnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}

inline void StPicoD0AnaMaker::setHFCuts(StHFCuts* cuts)   
{ 
  mHFCuts = cuts; 
}

#endif
