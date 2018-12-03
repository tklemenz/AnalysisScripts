#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <vector>
#include <fstream>
#include <iostream>

#include "TSystem.h"
#include "TROOT.h"
#include "TLine.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TChain.h"
#include "TPCSimulation/Point.h"
#include "TPCBase/Digit.h"
#include "TPCSimulation/Digitizer.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsTPC/Cluster.h"
#include "TPCBase/Mapper.h"
#include "TPCSimulation/SAMPAProcessing.h"
#endif

using namespace o2::TPC;

void studyClusters(TString FileList, const char* Filename="studyClusters.root",
                   const char* OutputPath="/home/gu74yub/SimAna/Results/compareSimvsReco/clusters") {

  Mapper &mapper = Mapper::instance();

  std::vector<o2::TPC::Cluster> *ClusterVector;

  TH1F *hQmax       = new TH1F("hQmax","; Qmax; Counts",1000,0,1000);
  TH1F *hQtot       = new TH1F("hQtot","; Qtot; Counts",1000,0,1000);
  TH1F *hPadMean    = new TH1F("hPadMean","; PadMean; Counts",100,0,100);
  TH1F *hTimeMean   = new TH1F("hTimeMean","; TimeMean; Counts",300,0,30);
  TH1F *hPadSigma   = new TH1F("hPadSigma","; PadSigma; Counts",100,0,1);
  TH1F *hTimeSigma  = new TH1F("hTimeSigma","; TimeSigma; Counts",100,0,1);


  /// ======| get files ready |=============================================================

    std::vector<o2::TPC::Cluster> *fChainEvent=0;

    TChain fChain("o2sim");

    fChain.SetBranchAddress("TPCClusterHW",&fChainEvent);

    TString allFiles;
    if (FileList.EndsWith(".txt")){
      allFiles=gSystem->GetFromPipe(Form("cat %s",FileList.Data()));

      TObjArray *arr = allFiles.Tokenize("\n");

      for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
        TString file=arr->At(ifile)->GetName();
        fChain.Add(file);
      }
    }
    else if (FileList.EndsWith(".root")){
      allFiles=gSystem->GetFromPipe(Form("ls %s",FileList.Data()));

      TObjArray *arr = allFiles.Tokenize("\n");

      for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
        TString file=arr->At(ifile)->GetName();
        fChain.Add(file);
      }
    }
    else {
      allFiles=gSystem->GetFromPipe(Form("ls %s",FileList.Data()));

      TObjArray *arr = allFiles.Tokenize("\n");

      for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
        TString file=arr->At(ifile)->GetName();
        file = FileList+"/"+file;
        fChain.Add(file);
      }
    }
    cout<<endl<<endl<<endl<<"Chain ready!!"<<endl<<endl;
    TObjArray *chainEntries = fChain.GetListOfFiles();
    for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
      cout<<chainEntries->At(ifile)->GetTitle()<<endl;
    }

  ///===================================================================================
  ///======| process clusters |=========================================================

    for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile) {
      TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
      cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;
      TTree *clusTree = (TTree*)TreeFile->Get("o2sim");

      ClusterVector = 0;
      clusTree->SetBranchAddress("TPCClusterHW", &ClusterVector);

      for (int event=0; event<clusTree->GetEntriesFast(); ++event) {

        clusTree->GetEntry(event);

        for (auto& clusterObject : *ClusterVector) {

          const float Qmax      = clusterObject.getQmax();
          const float Qtot      = clusterObject.getQ();
          const float PadMean   = clusterObject.getPadMean();
          const float TimeMean  = clusterObject.getTimeMean();
          const float PadSigma  = clusterObject.getPadSigma();
          const float TimeSigma = clusterObject.getTimeSigma();

          hQtot->       Fill(Qtot);
          hQmax->       Fill(Qmax);
          hPadMean->    Fill(PadMean);
          hTimeMean->   Fill(TimeMean);
          hPadSigma->   Fill(PadSigma);
          hTimeSigma->  Fill(TimeSigma);

          const CRU cru(clusterObject.getCRU());

          const PadRegionInfo& region     = mapper.getPadRegionInfo(cru.region());
          const int rowInSector           = clusterObject.getRow() + region.getGlobalRowOffset();
          const GlobalPadNumber pad       = mapper.globalPadNumber(PadPos(rowInSector, clusterObject.getPadMean()));
          const PadCentre& padCentre      = mapper.padCentre(pad);
          const float localYfactor        = (cru.side() == Side::A) ? -1.f : 1.f;
          float zPosition                 = SAMPAProcessing::getZfromTimeBin(clusterObject.getTimeMean(), cru.side());

          LocalPosition3D posLoc(padCentre.X(), localYfactor * padCentre.Y(), zPosition);
          GlobalPosition3D posGlob = Mapper::LocalToGlobal(posLoc, cru.sector());

          DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));

          const float PosX            = posGlob.X();
          const float PosYCentreGlob  = posGlob.Y();
          const float Row             = pos.getPadSecPos().getPadPos().getRow();
          const float PadLoc          = pos.getPadSecPos().getPadPos().getPad();

          const float relPadPos       = PadMean - PadLoc;

          const float PosY            = posGlob.Y() - 0.2+relPadPos*0.4;

        }
      }

    }

  ///===================================================================================

    TFile *OutFile = new TFile(Form("%s/%s", OutputPath,Filename), "recreate");

    OutFile->WriteObject(hQmax,       "hQmax");
    OutFile->WriteObject(hQtot,       "hQtot");
    OutFile->WriteObject(hPadMean,    "hPadMean");
    OutFile->WriteObject(hTimeMean,   "hTimeMean");
    OutFile->WriteObject(hPadSigma,   "hPadSigma");
    OutFile->WriteObject(hTimeSigma,  "hTimeSigma");

    OutFile->Close();
    delete OutFile;

}
