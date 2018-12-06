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
#include "TF1.h"
#include "TChain.h"
#include "TPCSimulation/Point.h"
#include "TPCBase/Digit.h"
//#include "TPCSimulation/Digitizer.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsTPC/Cluster.h"
#include "TPCBase/Mapper.h"
#endif

using namespace o2::TPC;

void studySamplingSim(TString FileList, const char* Filename="studySamplingSim.root", const char* OutputPath="/scratch2/tklemenz/SimAna/Results/newO2/studySamplingSim")
{
  Mapper &mapper = Mapper::instance();

  TH2F *hChargeVsTime2D          = new TH2F("ChargeVsTime2D",";#Deltat_Digit;charge",60,-30,30,1100,0,1100);
  TH1F *hChargeVsTime1D          = new TH1F("ChargeVSTime1D",";#Deltat_Digit;charge",60,-30,30);
//  TH2F *hChargeVsTime2D          = new TH2F("ChargeVsTime2D",";#Deltat_Digit;charge",2000,0,2000,1100,0,1100);
//  TH1F *hChargeVsTime1D          = new TH1F("ChargeVSTime1D",";#Deltat_Digit;charge",2000,0,2000);
  TH2F *hTimeVsRow2DForFit       = new TH2F("hTimeVsRow2DForFit",";row;time;charge",63,0,63,10000,0,10000);
  TH2F *hSingleEventTrack		 = new TH2F("hSingleEventTrack", ";row;pad;counts", 63,0,63,100,-50,50);

  TH2F *hPadOccupancy            = new TH2F("PadOccupancy", ";row;pad;counts", 63,0,63,100,-50,50);

  TH1F *hT0FitDist               = new TH1F("hT0FitDist", ";T0fromFit; counts", 10000,0,10000);
  TH2F *hRawChargeVsTime2D       = new TH2F("RawChargeVsTime2D",";#Deltat_Digit;charge",60,-30,30,1100,0,1100);

  std::vector<Digit> *TPCDigits = 0;

  /// ======| get files ready |=============================================================

  std::vector<Digit> *fChainEvent=0;

  TChain fChain("o2sim");

  fChain.SetBranchAddress("TPCDigit_0",&fChainEvent);

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
    digitTree->SetBranchAddress("TPCDigit_0", &TPCDigits);

    float maxCharge = -1;
    int maxTime = -1;
    float maxPad = -1;
    std::vector<float> maxTimeVector;
    int eventCounter = 0;

    /// find timeStamp of highest value per event => t0
    /// _____________________________________________________________________________________________________
/*    for (int event=0; event<digitTree->GetEntriesFast(); ++event) {
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
*/      //std::cout << "maxCharge: " << maxCharge << "\t maxTime: " << maxTime << std::endl;
    //}
    /// _____________________________________________________________________________________________________


    /// fit function for t0 fit and "tracking"
    /// _____________________________________________________________________________________________________

    TF1* FitT0 = new TF1("FitT0", "[0]*x+[1]",0,63);
    FitT0->SetParName(0,"slopeT0");
    FitT0->SetParName(1,"constantT0");
    FitT0->SetParameter(0,1);
    FitT0->SetParameter(1,0);

    TF1* FitPadOcc = new TF1("FitPadOcc", "[0]*x+[1]",0,63);
    FitPadOcc->SetParName(0,"slopeTrack");
    FitPadOcc->SetParName(1,"constantTrack");
    FitPadOcc->SetParameter(0,1);
    FitPadOcc->SetParameter(1,0);

    /// _____________________________________________________________________________________________________


    for (int event=0; event<digitTree->GetEntriesFast(); ++event) {
      digitTree->GetEntry(event);
      hTimeVsRow2DForFit->Reset("ICESM");
      hSingleEventTrack->Reset("ICESM");
      for (auto& digit : *TPCDigits) {
        CRU cruObj = digit.getCRU();
        const int crurow = digit.getRow();
        const int crupad = digit.getPad();
        DigitPos pos(cruObj, PadPos(crurow, crupad));
        //const float row = pos.getPadSecPos().getPadPos().getRow();
        const int row = crurow;
        const int timeStamp = digit.getTimeStamp();
        const float charge = digit.getChargeFloat();
        const float pad = pos.getPadSecPos().getPadPos().getPad();
        const float cpad  = pad - mapper.getNumberOfPadsInRowSector(row)/2;

        if (row <= 4 || row >= 56) continue;
        //if (timeStamp <= 4 || timeStamp >= 14) continue;
        if (charge < 5) continue;

        hTimeVsRow2DForFit->Fill(row, timeStamp, charge);
        hSingleEventTrack->Fill(row,cpad,charge);
        hPadOccupancy->Fill(row,cpad,charge);

      }

      hTimeVsRow2DForFit->Fit(FitT0);
      hSingleEventTrack->Fit(FitPadOcc);

      for (int i = 0; i<63; i++) {
        if ((TMath::Abs(FitT0->Eval(i) - hTimeVsRow2DForFit->GetBinContent(i+1)) > 1.3) && hTimeVsRow2DForFit->GetBinContent(i+1) > 0)
        {
          continue;
        }
        else {
          hT0FitDist->Fill(FitT0->Eval(i));
        }
      }
      for (auto& digit : *TPCDigits) {
        CRU cruObj = digit.getCRU();
        const int crurow = digit.getRow();
        const int crupad = digit.getPad();
        DigitPos pos(cruObj, PadPos(crurow, crupad));
        const float row = crurow;//pos.getPadSecPos().getPadPos().getRow();
        const float pad = pos.getPadSecPos().getPadPos().getPad();

        const int timeStamp = digit.getTimeStamp();
        const float charge = digit.getChargeFloat();
        const int roc = cruObj.roc();
        const GlobalPadNumber padInROC = mapper.getPadNumberInROC(PadROCPos(roc, row, pad));

        if (row <= 4 || row >= 56) continue;
        //if (timeStamp <= 4 || timeStamp >= 14) continue;
        if (charge < 5) continue;

        const float deltaT = timeStamp - FitT0->Eval(row);
        hChargeVsTime2D->Fill(deltaT,charge);
        hChargeVsTime1D->Fill(deltaT,charge);
        hRawChargeVsTime2D->Fill(timeStamp,charge);

//        std::cout << "pad: " << padInROC << "\t " << "timeStamp: " << timeStamp << "\t" << "charge: " << charge
//                  << "      \t"<<std::endl;// << "padCounter: " << padCounter << "\t" << "nPads: " << nPads << "\t" << "timeStamp: "
                  //<< mChargeVectorVector[padCounter].size() << std::endl;
      }
      ++eventCounter;
    }
  }

  TFile Outputfile(Form("%s/%s", OutputPath, Filename), "recreate");
  Outputfile.WriteObject(hRawChargeVsTime2D, "hRawChargeVsTime_2D");
  Outputfile.WriteObject(hChargeVsTime2D, "hChargeVsTime_2D");
  Outputfile.WriteObject(hChargeVsTime1D, "hChargeVsTime_1D");
  Outputfile.WriteObject(hTimeVsRow2DForFit, "hTimeVsRow2DForFit_last_event");
  Outputfile.WriteObject(hSingleEventTrack, "hSingleEventTrack_last_event");
  Outputfile.WriteObject(hPadOccupancy, "hPadOccupancy");
  Outputfile.WriteObject(hT0FitDist, "hT0FitDist");
  Outputfile.Close();
}
