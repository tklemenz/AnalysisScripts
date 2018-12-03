#include "TPCSimulation/Cluster.h"
#include "TPCBase/Mapper.h"
#include "TPCReconstruction/TrackTPC.h"
#include "TPCBase/CRU.h"
#include "TPCBase/ROC.h"
#include "TPCBase/CalDet.h"
#include "TPCBase/Painter.h"


#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TChain.h"

#include <math.h>


using namespace o2::TPC;


int compareGainMaps(TString GainMap1, TString GainMap2, const char *OutputPath="/home/tom/myana/Results")
{
  CalDet<float> *gm1S0 = nullptr, *gm1S1 = nullptr, *gm2S0 = nullptr, *gm2S1 = nullptr;

  TFile f(GainMap1);
  gROOT->cd();
  f.GetObject("GainMapS0", gm1S0);
  f.GetObject("GainMapS1", gm1S1);
  TFile g(GainMap2);
  gROOT->cd();
  g.GetObject("GainMapS0", gm2S0);
  g.GetObject("GainMapS1", gm2S1);

  Mapper& mapper = Mapper::instance();
  ROC roc(0);

  TH1F *ratioDistS0 = new TH1F ("ratio distribution S0","; ratio unbiased/biased; ;#counts",500,0,5);
  TH2F *ratioMapS0  = new TH2F ("ratio map S0","; row; pad",63,0,63,100,-50,50);
  TH1F *ratioDistS1 = new TH1F ("ratio distribution S1","; ratio unbiased/biased; ;#counts",500,0,5);
  TH2F *ratioMapS1  = new TH2F ("ratio map S1","; row; pad",63,0,63,100,-50,50);

  for (int irow=0; irow<63; ++irow) {
    int npads = mapper.getNumberOfPadsInRowROC(roc, irow);
    for (int ipad=0; ipad<npads; ++ipad) {
      float val1S0 = gm1S0->getValue(roc, irow, ipad);
      float val2S0 = gm2S0->getValue(roc, irow, ipad);
      float val1S1 = gm1S1->getValue(roc, irow, ipad);
      float val2S1 = gm2S1->getValue(roc, irow, ipad);
      float cpad = ipad - mapper.getNumberOfPadsInRowSector(irow)/2;

      if (val1S0 != 0 && val2S0 != 0) {
        ratioDistS0->Fill(val1S0/val2S0);
        ratioMapS0->Fill(irow,cpad,val1S0/val2S0);
      }

      if (val1S1 != 0 && val2S1 != 0) {
        ratioDistS1->Fill(val1S1/val2S1);
        ratioMapS1->Fill(irow,cpad,val1S1/val2S1);
      }
    }
  }
  TFile *OutFile = new TFile(Form("%s/compareGainMaps.root", OutputPath), "recreate");
  OutFile->cd();
  ratioDistS0->Write();
  ratioDistS1->Write();
  ratioMapS0->Write();
  ratioMapS1->Write();
  OutFile->Close();
  return 0;
}
