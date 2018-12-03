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

struct EventHeader
{
  int run;
  float cherenkovValue;
  int beamMomentum;
  int powerSupply;
  int HVSettings;
  int trigger;
  int dataType;
  int driftFieldStrength;
  int runType;
};

void makeHitMap(TString FileList, const char *OutputPath="/home/tom/myana/Results") {
  EventHeader Header;
  ROC roc(0);

  int runNr = 0;

  Mapper& mapper = Mapper::instance();

  CalDet<int> *HitMapS0 = new CalDet<int>(PadSubset::ROC);
  HitMapS0->setName("HitMapS0");
  CalArray<int>& HitMapArrayS0 = HitMapS0->getCalArray(roc);

  CalDet<int> *HitMapS1 = new CalDet<int>(PadSubset::ROC);
  HitMapS1->setName("HitMapS1");
  CalArray<int>& HitMapArrayS1 = HitMapS1->getCalArray(roc);

  TH2D *hHitMapS0 = new TH2D ("HitMapS0", "; Row; Pad", 63,0,63,100,0,100);
  TH2D *hHitMapS1 = new TH2D ("HitMapS1", "; Row; Pad", 63,0,63,100,0,100);
  //TH2D *TestHisto = new TH2D ("bla", "; row; pad", 63,0,63,100,0,100);

  for (int irow = 0; irow<63; ++irow) {
    for (int ipad = -1; ipad<101; ++ipad){
      HitMapArrayS0.setValue(irow, ipad, 0);
      HitMapArrayS1.setValue(irow, ipad, 0);
      hHitMapS0->Fill(irow, ipad, 0);
      hHitMapS1->Fill(irow, ipad, 0);
    }
  }

  std::vector<TrackTPC> *fChainEvent=0;

  TChain fChain("events");

  fChain.SetBranchAddress("Tracks",&fChainEvent);
  fChain.SetBranchAddress("header",&Header);

  TString allFiles;
  if (FileList.EndsWith(".txt")){
    allFiles=gSystem->GetFromPipe(Form("cat %s",FileList.Data()));
  }
  else {
    allFiles=gSystem->GetFromPipe(Form("ls %s",FileList.Data()));
  }
  TObjArray *arr = allFiles.Tokenize("\n");

  for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
    TString file=arr->At(ifile)->GetName();
    fChain.Add(file);
  }
  cout<<endl<<endl<<endl<<"Chain ready!!"<<endl<<endl;
  TObjArray *chainEntries = fChain.GetListOfFiles();
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
    cout<<chainEntries->At(ifile)->GetTitle()<<endl;
  }

  std::vector<TrackTPC> *vecEvent = 0;

  cout<<endl<<endl<<endl<<"Number of files to process: "<<chainEntries->GetEntriesFast()<<endl<<endl<<endl;
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
    TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
    cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;//<<endl;
    TTree *tree = (TTree*)TreeFile->Get("events");

    vecEvent=0;
    tree->SetBranchAddress("Tracks", &vecEvent);
    tree->SetBranchAddress("header", &Header);vecEvent=0;

    runNr = 0;

    for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
      tree->GetEntry(iev);
      runNr = Header.run;
      if (vecEvent->size() != 1) continue;
      for (auto trackObject : *vecEvent) {
        std::vector<Cluster> clCont;
        trackObject.getClusterVector(clCont);
        if (clCont.size() < 48) continue;
        for (auto &clusterObject : clCont) {
          DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
          float row = pos.getPadSecPos().getPadPos().getRow();
          float pad = pos.getPadSecPos().getPadPos().getPad();
          if (runNr <= 243) {
            HitMapArrayS0.setValue(row, pad, 1);
            hHitMapS0->SetBinContent(row+1, pad+1, 1);
          }
          else if (runNr > 255) {
            HitMapArrayS1.setValue(row, pad, 1);
            hHitMapS1->SetBinContent(row+1, pad+1, 1);
          }
        }
      }
    }
  }

//  for (int irow = 0; irow<63; ++irow) {
//    for (int ipad = 10; ipad<90; ++ipad) {
//      int HitS0val = HitMapS0->getValue(roc,irow,ipad);
//      if (HitS0val == 1) {
//        TestHisto->Fill(irow, ipad);
//      }
//    }
//  }

//  TestHisto->Draw("colz");

  TFile *f = TFile::Open(Form("%s/HitMap.root",OutputPath), "recreate");
  f->WriteObject(HitMapS0, "HitMapS0");
  f->WriteObject(HitMapS1, "HitMapS1");
  delete f;

  TFile *g = TFile::Open(Form("%s/HitMap_Histo.root",OutputPath), "recreate");
  g->WriteObject(hHitMapS0, "HitMapHistoS0");
  g->WriteObject(hHitMapS1, "HitMapHistoS1");
  delete g;
}
