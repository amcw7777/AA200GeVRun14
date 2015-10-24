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
#include "TH1D.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StPhysicalHelixD.hh"
class StPrimaryVertex; 
class StEvent;
class StDcaGeometry; 
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
    ofstream fout;

    bool isGoodPair(StKaonPion const*) const;
    int isD0Pair(StKaonPion const*) const;
    int isD0PairOld(StKaonPion const*) const;
    int D0Reco(StThreeVectorF *);
    bool isGoodEvent();
    bool  isGoodTrack(StPicoTrack const*) const;
    bool  isGoodHadron(StPicoTrack const*) const;
    bool  isTpcPion(StPicoTrack const*) const;
    bool  isTpcKaon(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool isTofKaon(StPicoTrack const* const, float beta) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool getCorHadron(float eta, vector<float> &hadronsPhi, const vector<unsigned int>, float, float);
    float sumCos(float phi, vector<float> &hadronsPhi);
    bool fixPhi(vector<float> &phi);
    bool getHadronCorV2();
    float getD0CorV2(int *sumPair, vector<const StKaonPion *> cand);

    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
    StPicoPrescales* mPrescales;
    StRefMultCorr* mGRefMultCorrUtil;
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
    TNtuple *mEventTuple;
		TNtuple *mDTuple;
		TNtuple *mHadronTuple;
    TH1F *dEtaDHadron;
    TH1F *hEtaD;
    TH1F *hEtaHadron;
    TH1F *hPhiHadron;
    TH1F *hPhiD;
    TH1F *vtxz;
		TH2F *etaPhi;
		TH2F *etaPhi_D;
		TH2F *etaPhi_Hadron;
		TH2F *etaPhi_Hadron_all;
    TProfile *profV2[8][5][3];//i.S or B; j.flatten; k. differetn etaGap
    TH1D *hadronV2[5];
    TH1D *hadronV2_sum[5];
    TH1D *hadronV2_excl[5][9];
    TH1D *hadronV2_excl_sum[5][9];
    TH2D *fitPhi[6];
    TH2D *massPt;
    TH2D *massPtLike;
    TH2D *massLike;
    TH2D *massLike2;
    TH2D *massUnlike;
    TH2D *v2Weight[8][3];
    TH2D *likeV2Mass[6][5];
    TH2D *likeV2Mass2[6][5];
    TH2D *unlikeV2Mass[6][5];
    TNtuple *checkNew;
    TNtuple *checkOld;

  TH2D *checkPeak;
    
    double fitsigma[6];
    double fitmean[6];
    double mHadronV2[9];
  double efficiency[4][6];
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
