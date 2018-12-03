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
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TPCSimulation/Point.h"
#include "TPCBase/Digit.h"
#include "TPCSimulation/Digitizer.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsTPC/Cluster.h"
#include "TPCBase/Mapper.h"
#endif

using namespace o2::TPC;

void studySamplingSim(TString FileList, const char* Filename="studySamplingSim.root", const char* OutputPath="/home/gu74yub/SimAna/Results/studySamplingSim")
{
  Mapper &mapper = Mapper::instance();

  TH2F *hChargeVsTime2D          = new TH2F("ChargeVsTime2D",";#Deltat_Digit;charge",40,-20,20,1100,0,1100);
  TH1F *hChargeVsTime1D          = new TH1F("ChargeVSTime1D",";#Deltat_Digit;charge",40,-20,20);
  TH2F *hPadOccupancy            = new TH2F("PadOccupancy", ";row;pad;counts", 63,0,63,100,-50,50);

  std::vector<Digit> *TPCDigits = 0;

  /// ======| get files ready |=============================================================

  std::vector<Digit> *fChainEvent=0;

  TChain fChain("o2sim");

  fChain.SetBranchAddress("TPCDigit",&fChainEvent);

  TString allFiles;
  if (FileList.EndsWith(".txt")){
    allFiles=gSystem->GetFromPipe(Form("cat %s",FileList.Data()));

    TObjArray *arr = allFiles.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      fChain.Add(file);
    }
  }
  else if (FileList.EndsWith(".root")){
    allFiles=gSystem->GetFromPipe(Form("ls %s",FileList.Data()));

    TObjArray *arr = allFiles.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      fChain.Add(file);
    }
  }
  else {
    allFiles=gSystem->GetFromPipe(Form("ls %s",FileList.Data()));

    TObjArray *arr = allFiles.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
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
  ///======| process digits |=============================================================


  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile) {
    TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
    cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;
    TTree *digitTree = (TTree*)TreeFile->Get("o2sim");

    TPCDigits = 0;
    digitTree->SetBranchAddress("TPCDigit", &TPCDigits);

    float maxCharge = -1;
    int maxTime = -1;
    float maxPad = -1;
    std::vector<float> maxTimeVector;
    int eventCounter = 0;

    /// find timeStamp of highest value per event => t0
    /// _____________________________________________________________________________________________________
    for (int event=0; event<digitTree->GetEntriesFast(); ++event) {
      digitTree->GetEntry(event);
      maxCharge = -1;
      maxTime = -1;
      maxPad = -1;
      for (auto& digit : *TPCDigits) {
        CRU cruObj = digit.getCRU();
        const int crurow = digit.getRow();
        const int crupad = digit.getPad();
        DigitPos pos(cruObj, PadPos(crurow, crupad));
        const float row = pos.getPadSecPos().getPadPos().getRow();
        const float pad = pos.getPadSecPos().getPadPos().getPad();
        const float cpad  = pad - mapper.getNumberOfPadsInRowSector(row)/2;
        const int timeStamp = digit.getTimeStamp();
        const float charge = digit.getChargeFloat();
        const int roc = cruObj.roc();
        const GlobalPadNumber padInROC = mapper.getPadNumberInROC(PadROCPos(roc, row, pad));
        if (charge > maxCharge && charge < 150) {
          maxCharge = charge;
          maxTime = timeStamp;
          maxPad = padInROC;
        }
      }
      maxTimeVector.push_back(maxTime);
      //std::cout << "maxCharge: " << maxCharge << "\t maxTime: " << maxTime << std::endl;
    }
    /// _____________________________________________________________________________________________________


    for (int event=0; event<digitTree->GetEntriesFast(); ++event) {
      digitTree->GetEntry(event);
      for (auto& digit : *TPCDigits) {
        CRU cruObj = digit.getCRU();
        const int crurow = digit.getRow();
        const int crupad = digit.getPad();
        DigitPos pos(cruObj, PadPos(crurow, crupad));
        const float row = pos.getPadSecPos().getPadPos().getRow();
        const float pad = pos.getPadSecPos().getPadPos().getPad();
        const float cpad  = pad - mapper.getNumberOfPadsInRowSector(row)/2;

        const int timeStamp = digit.getTimeStamp();
        const float charge = digit.getChargeFloat();
        const int roc = cruObj.roc();
        const GlobalPadNumber padInROC = mapper.getPadNumberInROC(PadROCPos(roc, row, pad));

        const float deltaT = timeStamp - maxTimeVector[eventCounter];
        hChargeVsTime2D->Fill(deltaT,charge);
        hChargeVsTime1D->Fill(deltaT,charge);
        hPadOccupancy->Fill(row,cpad,charge);

        //std::cout << "pad: " << padInROC << "\t " << "timeStamp: " << timeStamp << "\t" << "charge: " << charge
        //          << "      \t"<<std::endl;// << "padCounter: " << padCounter << "\t" << "nPads: " << nPads << "\t" << "timeStamp: "
                  //<< mChargeVectorVector[padCounter].size() << std::endl;
      }
      ++eventCounter;
    }
  }

  TFile Outputfile(Form("%s/%s", OutputPath, Filename), "recreate");
  Outputfile.WriteObject(hChargeVsTime2D, "hChargeVsTime_2D");
  Outputfile.WriteObject(hChargeVsTime1D, "hChargeVsTime_1D");
  Outputfile.WriteObject(hPadOccupancy, "hPadOccupancy");
  Outputfile.Close();
}
