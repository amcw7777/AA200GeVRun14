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
#include <bitset>
#include <climits>
#include <vector>
#include "TChain.h"
#include "StMaker.h"
//
#include "StThreeVectorF.hh"
#include "TSpectrum.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"
#include "TCanvas.h"
#include "TH1K.h"
#include "TH1D.h"
#include "TH2F.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StPhysicalHelixD.hh"
class StPrimaryVertex; 
class StEvent;
class StDcaGeometry; 
class TMinuit;
class KFParticle; 
class StKFVerticesCollection; 
class StKFVertexMaker;
class TSpectrum; 
class StAnneling;
//// 


class TString;
class TFile;
class TNtuple;
class StPicoD0Event;
class StKaonPion;
class StPicoDstMaker;
class StPicoDst;
class StPicoTrack;
class StHFCuts;
class StPicoPrescales;
class StRefMultCorr;



class StPicoD0AnaMaker : public StMaker
{
  public:
    StPicoD0AnaMaker(char const * name, char const * inputFilesList, 
        char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil);
    virtual ~StPicoD0AnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;

    void setHFCuts(StHFCuts* cuts);    

  private:
    StPicoD0AnaMaker() {}
    void readNextEvent();

    bool isGoodPair(StKaonPion const*) const;
    int isD0Pair(StKaonPion const*) const;
    int isOptD0Pair(StKaonPion const*) const;
    int primaryVertexRefit(StThreeVectorF *, vector<int>& daughter);
    int D0Reco(StThreeVectorF *);
    bool isGoodEvent();
    bool  isGoodTrack(StPicoTrack const*) const;
    bool  isTpcPion(StPicoTrack const*) const;
    bool  isTpcKaon(StPicoTrack const*) const;
    bool  isKaon(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool isTofKaon(StPicoTrack const* const, float beta) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool  isGoodForVertexFit(StPicoTrack const*,StThreeVectorF const& vtx) const;
    size_t popcount(size_t) const;

    StRefMultCorr* mGRefMultCorrUtil;
    // Code from Xin Dong to do Minuit vertex fit
    StThreeVectorF vtxReFit(StPicoDst*);
    bool    accept(StPicoTrack*) const;   // track filter
    static void fcn(Int_t&, Double_t*, Double_t&, Double_t*, Int_t); // fit function
    static Double_t Chi2atVertex(StThreeVectorD &vtx);

    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
    StPicoPrescales* mPrescales;
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
     TNtuple *mEventtuple;
     TNtuple *mOrigin;
     TNtuple *mTest;
     TNtuple *mTest1;
     TNtuple *mRefit;
     TH2F *timemult;
     TH1D *mMult;
     StDcaGeometry *dcaG;
     TH1D *mDmass_like;
     TH1D *mDmasscut_like;
     TH1D *mDmasstest_like;
     TH1D *mDmass_unlike;
     TH1D *mDmasscut_unlike;
     TH1D *mDmasstest_unlike;
     TH1D *mCostheta;
     TH1D *mdRkp;
//
    
    // variables to do Minuit vertex fit
    enum  {kFlagDcaz = 1, kFlagCTBMatch = 2, kFlagBEMCMatch = 4, kFlagCrossMembrane = 8};
    static vector<StDcaGeometry>   mDCAs;
    static vector<StPhysicalHelixD> mHelices;
    static vector<UShort_t>         mHelixFlags;
    static vector<Double_t>         mSigma;
    static vector<Double_t>         mZImpact;
    Int_t                  mNSeed;
        
    static Double_t                 mWidthScale;
    Int_t                    mStatusMin;           // Minuit status flag 
    TMinuit*                 mMinuit;
 
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

inline size_t StPicoD0AnaMaker::popcount(size_t n) const
{
    std::bitset<sizeof(size_t) * CHAR_BIT> b(n);
    return b.count();
}


#endif
