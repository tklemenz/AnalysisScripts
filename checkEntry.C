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
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <math.h>
#include <boost/lambda/lambda.hpp>



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

void checkEntry(TString FileList, const char *OutputPath="/home/tom/myana/Results") {
  //gROOT->ProcessLine(".x /lustre/nyx/alice/users/tklemenz/rootlogon.C");
  //gStyle->SetOptStat(0);

  EventHeader Header;

  TFile *g = TFile::Open(Form("%s/checkEntry.root", OutputPath), "recreate");

  const Mapper& mapper = Mapper::instance();

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
    //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
    fChain.Add(file);
  }
  cout<<endl<<endl<<endl<<"Chain ready!!"<<endl<<endl;
  TObjArray *chainEntries = fChain.GetListOfFiles();
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
    cout<<chainEntries->At(ifile)->GetTitle()<<endl;
  }

  /* ============================================================
   * ============= Define histograms and graphs =================
   * ============================================================*/
  TH1F *hPiDistS0       = new TH1F("hPiDistS0","; Y; Counts",200,0,20);
  TH1F *hEleDistS0      = new TH1F("hEleDistS0","; Y; Counts",200,0,20);
  TH1F *hPiDistS1       = new TH1F("hPiDistS1","; Y; Counts",500,-250,250);
  TH1F *hEleDistS1      = new TH1F("hEleDistS1","; Y; Counts",500,-250,250);

  TH1F *hDistS0         = new TH1F("hDistS0","; Y; Counts",500,-250,250);
  TH1F *hDistS1         = new TH1F("hDistS1","; Y; Counts",500,-250,250);

  TH2D *hXYDist         = new TH2D("hXYDist","; X; Y",120,80,140,60,0,15);
  TH1F *hX              = new TH1F("hX","; x; Counts",150,0,150);
  TH1F *hY              = new TH1F("hY","; y; Counts",60,0,15);
  TH1F *hAlpha          = new TH1F("hAlpha","; Alpha; Counts",100,0,1);
  TH1F *hSnp            = new TH1F("hSnp","; Snp; Counts",100,0,0.2);
  TH1F *hTgl            = new TH1F("hTgl","; Tgl; Counts",200,-0.1,0.1);
  TH1F *hCurvature      = new TH1F("hCurvature","; Curvature; Counts",100,0,10);
  TH1F *hPhi            = new TH1F("hPhi","; Phi; Counts",100,0,1);

  int trackCounter = 0;
  float CherenkovValue = 0.;
  float CherCutLow = 0.009;
  float CherCutHigh = 0.011;
  int setting = 0;
  int runNr = 0;
  float yEntryPi = 0;
  float yEntryEle = 0;
  std::vector<TrackTPC> *vecEvent = 0;
  cout<<endl<<endl<<endl<<"Number of files to process: "<<chainEntries->GetEntriesFast()<<endl<<endl<<endl;

  /*=============================================== Loop over runs =========================================================*/
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
    TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
    cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;//<<endl;
    TTree *tree = (TTree*)TreeFile->Get("events");

    vecEvent=0;
    tree->SetBranchAddress("Tracks", &vecEvent);
    tree->SetBranchAddress("header", &Header);

    CherenkovValue = 0;
    runNr = 0;

    for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
      tree->GetEntry(iev);
      int NTracks = vecEvent->size();
      if (NTracks != 1) continue;
      CherenkovValue = Header.cherenkovValue;
      runNr = Header.run;
       setting = 0;
      if (runNr >= 255) {setting = 1;}
      for (auto trackObject : *vecEvent) {
        std::vector<Cluster> clCont;
        trackObject.getClusterVector(clCont);
        if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;
        ++trackCounter;
//        for (auto &clusterObject : clCont) {
//          DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
//          float row = pos.getPadSecPos().getPadPos().getRow();
//          const CRU cru(clusterObject.getCRU());

//          const PadRegionInfo& region = mapper.getPadRegionInfo(cru.region());
//          const int rowInSector       = clusterObject.getRow() + region.getGlobalRowOffset();
//          const GlobalPadNumber pad   = mapper.globalPadNumber(PadPos(rowInSector, clusterObject.getPad()));
//          const PadCentre& padCentre  = mapper.padCentre(pad);
//          const float localYfactor    = (cru.side()==Side::A)?-1.f:1.f;
//          float zPosition       = Digitizer::getZfromTimeBin(clusterObject.getTimeStamp(), cru.side());

//          LocalPosition3D posLoc(padCentre.X(), localYfactor*padCentre.Y(), zPosition);
//          GlobalPosition3D posGlob = Mapper::LocalToGlobal(posLoc, cru.sector());

//          float clusX = posGlob.X();
//          float clusY = posGlob.Y();
//          float clusZ = zPosition;

//          if (CherenkovValue < CherCutLow) {
//            if (setting == 0) {
//              hPiDistS0->Fill(clusY);
//              hDistS0->Fill(clusY);
//            }
//            else if (setting ==1) {
//              hPiDistS1->Fill(clusY);
//              hDistS1->Fill(clusY);
//            }
//          }

//          else if (CherenkovValue > CherCutHigh) {
//            if (setting == 0) {
//              hEleDistS0->Fill(clusY);
//              hDistS0->Fill(clusY);
//            }
//            else if (setting ==1) {
//             hEleDistS1->Fill(clusY);
//             hDistS1->Fill(clusY);
//            }
//          }
//        }

        float x          = trackObject.getX();
        float y          = trackObject.getY();
        float alpha      = trackObject.getAlpha();
        float Snp        = trackObject.getSnp();
        float Tgl        = trackObject.getTgl();
        //float curvature  = trackObject.getCurvature(b);
        float Phi        = trackObject.getPhi();

        if (CherenkovValue < CherCutLow) {
          yEntryPi   = y + TMath::Tan(TMath::ASin(Snp))*(131.725-x);
        }
        else if (CherenkovValue > CherCutHigh) {
          yEntryEle   = y + TMath::Tan(TMath::ASin(Snp))*(131.725-x);
        }

        hXYDist->Fill(x,y);
        hX->Fill(x);
        hY->Fill(y);
        hAlpha->Fill(alpha);
        hSnp->Fill(Snp);
        hTgl->Fill(Tgl);
        hPhi->Fill(Phi);
        hPiDistS0->Fill(yEntryPi);
        hEleDistS0->Fill(yEntryEle);

      }
    }
  }

  std::cout << "\n \n \n" << trackCounter << " Tracks used" << "\n \n \n";
  g->WriteObject(hDistS0, "hDistS0");
  g->WriteObject(hDistS1, "hDistS1");
  g->WriteObject(hPiDistS0, "hPiDistS0");
  g->WriteObject(hEleDistS0, "hEleDistS0");
  g->WriteObject(hPiDistS1, "hPiDistS1");
  g->WriteObject(hEleDistS1, "hEleDistS1");
  g->WriteObject(hXYDist, "hXYDist");
  g->WriteObject(hX, "hX");
  g->WriteObject(hY, "hY");
  g->WriteObject(hAlpha, "hAlpha");
  g->WriteObject(hSnp, "hSnp");
  g->WriteObject(hTgl, "hTgl");
  g->WriteObject(hPhi, "hPhi");

  delete g;
}

