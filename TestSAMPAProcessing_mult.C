#include "TPCSimulation/SAMPAProcessing.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TString.h"

using namespace o2::TPC;

void TestSAMPAProcessing_mult(int nSampling, float TimeBin, TString particle="pi", const char *OutputPath="/home/tom/myana/Results/TestSAMPAProcessing/testSAMPAProcessing_mult.root") {

  /// NBins is the number of bins, one time bin is divided in
  /// TimeBin [ns] corresponds to the Sampling frequency, 200 -> 5 MHz, 50 -> 20 MHz

  const SAMPAProcessing& sampa = SAMPAProcessing::instance();

  TRandom3 *rRandom = new TRandom3();

  float ADCsignal = .0;
  float amplitude = 120;
  double timeShift = .0;
  float randomNumber = 0.;
  double Qmax = .0;
  int QmaxInd = 0;
  double Qtot = 0.;
  float TimeBinSize = TimeBin/1000;
///  int nbins = NBins;    NBins needs to be variable!!
  std::vector<float> signalArray = { 0.,0.,0.,0.,0.,0.,0.,0. };

  std::vector<float> inputSigma = { 9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14 };

  int nHist = inputSigma.size();

  std::vector<float> outputSigma = { 0,0,0,0,0,0,0,0,0,0,0 };

  TH2D *hResult[nHist];
  TH2D *hResultQmax[nHist];
  TH1F *hCharge[nHist];
  TH1F *hQmaxDist[nHist];
  TH1F *hQtotDist[nHist];
//  TF1  *gauss[nHist];

  TH1F *hQmaxQtotRatio = new TH1F("hQmaxQtotRatio","; Q_{max}/Q_{tot}; Counts",120,0,1.2);
  TH2D *hQmaxQtotRatioVsDt = new TH2D("hQmaxQtotRatioVsDt",";#Delta t; Q_{max}/Q_{tot}",500,0,200e-3f,120,0,1.2);

  TGraphErrors *gInSigVsOutSig = new TGraphErrors();
  gInSigVsOutSig->GetXaxis()->SetTitle("#sigma_{in} [%]");
  gInSigVsOutSig->GetYaxis()->SetTitle("#sigma_{out} [%]");
  gInSigVsOutSig->SetMarkerStyle(20);

  for (int i=0; i<nHist; ++i) {

    hResult[i] = new TH2D(Form("hResult_%f",inputSigma[i]),"; nr; Signal",8,0,8,15000,0,150);
    hResultQmax[i] = new TH2D(Form("hResultQmax_%f",inputSigma[i]),"; nr; Q_{max}",8,0,8,15000,0,150);
    hCharge[i] = new TH1F(Form("hCharge_%f",inputSigma[i]),"; Q; Counts",480,0,120);
    hQmaxDist[i] = new TH1F(Form("hQmaxDist_%f",inputSigma[i]),"; Q_{max}; Counts",1500,0,150);
    if (TimeBin == 50) {
      hQtotDist[i] = new TH1F(Form("hQtotDist_%f",inputSigma[i]),"; Q_{tot}; Counts",700,150,850);
    }
    else {
      hQtotDist[i] = new TH1F(Form("hQtotDist_%f",inputSigma[i]),"; Q_{tot}; Counts",300,0,300);
    }
//    gauss[i] = new TF1(Form("hQtotDistFit[%i]", i),"gaus");

//    TH1D *htimeShift = new TH1D("htimeShift","; #Delta t",500,0,200e-3f);
//    TH2D *hQmaxVstimeShift = new TH2D("hQmaxVstimeShift","; #Delta t; Q_{max}",500,0,200e-3f,360,0,120);
//    TH2D *hQtotVstimeShift = new TH2D("hQtotVstimeShift","; #Delta t; Q_{tot}",500,0,200e-3f,1200,0,400);
//    TH2D *hQ1VstimeShift = new TH2D("hQ1VstimeShift","; #Delta t; Q_{1}",500,0,200e-3f,360,0,120);
//    TH2D *hQ2VstimeShift = new TH2D("hQ2VstimeShift","; #Delta t; Q_{2}",500,0,200e-3f,360,0,120);
//    TH2D *hQ3VstimeShift = new TH2D("hQ3VstimeShift","; #Delta t; Q_{3}",500,0,200e-3f,360,0,120);
//    TH2D *hQ4VstimeShift = new TH2D("hQ4VstimeShift","; #Delta t; Q_{4}",500,0,200e-3f,360,0,120);
//    TH1F *hQmaxInd = new TH1F("QmaxInd","; Q_{max} Index; Counts",8,0,8);


///  TH1F *hQmaxDistTime[nbins];
///  TH1F *hQtotDistTime[nbins];
  /// Q1 is the sampled value of QmaxInd-2
  /// Q2 is the sampled value of QmaxInd-1
  /// Q3 is the sampled value of QmaxInd+1
  /// Q4 is the sampled value of QmaxInd+2
///  TH1F *hQ1DistTime[nbins];
///  TH1F *hQ2DistTime[nbins];
///  TH1F *hQ3DistTime[nbins];
///  TH1F *hQ4DistTime[nbins];
///  TH2D *hResultTime[nbins];
//  TH2D *hQ1VstimeShift[nbins];
//  TH2D *hQ2VstimeShift[nbins];
//  TH2D *hQ3VstimeShift[nbins];
//  TH2D *hQ4VstimeShift[nbins];
///  TH2D *hResultTimeQmax[nbins];

/*
  for (int i = 0; i < nbins; ++i) {
    hQmaxDistTime[i] = new TH1F(Form("hQmaxDistTime_%i",i),"; Q_{max}; Counts",12000,0,120);
    hQtotDistTime[i] = new TH1F(Form("hQtotDistTime_%i",i),"; Q_{tot}; Counts",40000,0,400);
    hQ1DistTime[i] = new TH1F(Form("hQ1totDistTime_%i",i),"; Q_{1}; Counts",12000,0,120);
    hQ2DistTime[i] = new TH1F(Form("hQ2totDistTime_%i",i),"; Q_{2}; Counts",12000,0,120);
    hQ3DistTime[i] = new TH1F(Form("hQ3totDistTime_%i",i),"; Q_{3}; Counts",12000,0,120);
    hQ4DistTime[i] = new TH1F(Form("hQ4totDistTime_%i",i),"; Q_{4}; Counts",12000,0,120);
    hResultTime[i] = new TH2D(Form("hResultTime_%i",i),"; nr; Signal",8,0,8,12000,0,120);
    hResultTimeQmax[i] = new TH2D(Form("hResultTimeQmax_%i",i),"; nr; Q_{max}",8,0,8,12000,0,120);
    if (TimeBin == 200) {
      hQtotDistTime[i]->GetXaxis()->SetRangeUser(70,130);
    }
  }
*/

/// START OF SAMPLING TEST

    std::cout << std::endl << "sampling signal " << nSampling << " times with input sigma " << inputSigma[i] << "% ...";

    for (int j = 0; j < nSampling; ++j) {


      signalArray = {0.,0.,0.,0.,0.,0.,0.,0.};
      Qmax = 0.;
      QmaxInd = 0;
      Qtot = 0.;
//    timeShift = rRandom->Gaus(0, 0.16);

      randomNumber = rRandom->Rndm();
      timeShift = TimeBinSize * randomNumber;

      ADCsignal = rRandom->Gaus(amplitude,amplitude*(inputSigma[i]/100));

      //htimeShift->Fill(timeShift);
      sampa.getShapedSignal(ADCsignal, timeShift, signalArray);

    ///find Qmax
      Qmax = *std::max_element(signalArray.begin(), signalArray.end());

    ///check index of Qmax
      QmaxInd = std::distance(signalArray.begin(), std::find(signalArray.begin(), signalArray.end(), Qmax));
      //hQmaxInd->Fill(QmaxInd);

      if (Qmax > 4) {
        //hQmaxVstimeShift->Fill(timeShift, Qmax);
        hQmaxDist[i]->Fill(Qmax);
        /*if (signalArray[QmaxInd-2] > 4) {
          hQ1VstimeShift->Fill(timeShift,signalArray[QmaxInd-2]);
        }
        if (signalArray[QmaxInd-1] > 4) {
          hQ2VstimeShift->Fill(timeShift,signalArray[QmaxInd-1]);
        }
        if (signalArray[QmaxInd+1] > 4) {
          hQ3VstimeShift->Fill(timeShift,signalArray[QmaxInd+1]);
        }
        if (signalArray[QmaxInd+2] > 4) {
          hQ4VstimeShift->Fill(timeShift,signalArray[QmaxInd+2]);
        }*/
      }






    ///calculate Qtot
      if (QmaxInd == 0) {
        if (signalArray[QmaxInd] > 4) {
          Qtot += signalArray[QmaxInd];
        }
        if (signalArray[QmaxInd +1] > 4) {
          Qtot += signalArray[QmaxInd+1];
        }
        if (signalArray[QmaxInd+2] > 4) {
          Qtot += signalArray[QmaxInd+2];
        }
      }
      else if (QmaxInd == 1) {
        if (signalArray[QmaxInd-1] > 4) {
          Qtot += signalArray[QmaxInd-1];
        }
        if (signalArray[QmaxInd] > 4) {
          Qtot += signalArray[QmaxInd];
        }
        if (signalArray[QmaxInd+1] > 4) {
          Qtot += signalArray[QmaxInd+1];
        }
        if (signalArray[QmaxInd+2] > 4) {
          Qtot += signalArray[QmaxInd+2];
        }
      }
      else if (QmaxInd == 6) {
        if (signalArray[QmaxInd-2] > 4){
          Qtot += signalArray[QmaxInd-2];
        }
        if (signalArray[QmaxInd-1] > 4) {
          Qtot += signalArray[QmaxInd-1];
        }
        if (signalArray[QmaxInd] > 4) {
          Qtot += signalArray[QmaxInd];
        }
        if (signalArray[QmaxInd+1] > 4) {
          Qtot += signalArray[QmaxInd+1];
        }
        if (signalArray[QmaxInd+2] > 4) {
          Qtot += signalArray[QmaxInd+2];
        }
     }
      else if (QmaxInd == 7) {
        if (signalArray[QmaxInd-2] > 4) {
          Qtot += signalArray[QmaxInd-2];
        }
        if (signalArray[QmaxInd-1] > 4) {
          Qtot += signalArray[QmaxInd-1];
        }
        if (signalArray[QmaxInd] > 4) {
          Qtot += signalArray[QmaxInd];
        }
      }
     else {
       if (signalArray[QmaxInd-2] > 4) {
          Qtot += signalArray[QmaxInd-2];
        }
        if (signalArray[QmaxInd-1] > 4) {
          Qtot += signalArray[QmaxInd-1];
        }
        if (signalArray[QmaxInd] > 4) {
          Qtot += signalArray[QmaxInd];
        }
        if (signalArray[QmaxInd+1] > 4) {
          Qtot += signalArray[QmaxInd+1];
        }
        if (signalArray[QmaxInd+2] > 4) {
          Qtot += signalArray[QmaxInd+2];
        }
    }
    if (Qtot != 0) {
      //hQtotVstimeShift->Fill(timeShift,Qtot);
      hQtotDist[i]->Fill(Qtot);
      hQmaxQtotRatio->Fill(Qmax/Qtot);
      hQmaxQtotRatioVsDt->Fill(timeShift,Qmax/Qtot);
    }


/*    float max = 0;
    int index;
    for (int i = 0; i < signalArray.size(); ++i) {
      float tmp = signalArray[i];
      if (tmp > max) {
        max = tmp;
      }
      else if (tmp < max) {
        index = i-1;
        break;
      }
    }
*/
      hResultQmax[i]->Fill(QmaxInd,Qmax);

    ///check sampled values
      for (int k = 0; k < signalArray.size(); ++k) {
        hResult[i]->Fill(k,signalArray[k]);
        hCharge[i]->Fill(signalArray[k]);
//      if (signalArray[i] > 4) {
//        hQmaxVstimeShift->Fill(timeShift, signalArray[i]);
//      }
      }
/*
    for (int i = 1; i <= nbins; ++i) {
      if (randomNumber >= (i-1)/float(nbins) && randomNumber < i/float(nbins)) {
        hQmaxDistTime[i-1]->Fill(Qmax);
        hQtotDistTime[i-1]->Fill(Qtot);
        for (int j = 0; j < signalArray.size(); ++j) {
          hResultTime[i-1]->Fill(j,signalArray[j]);
        }
        hResultTimeQmax[i-1]->Fill(QmaxInd,Qmax);
        hQ1DistTime[i-1]->Fill(signalArray[QmaxInd-2]);
        hQ2DistTime[i-1]->Fill(signalArray[QmaxInd-1]);
        hQ3DistTime[i-1]->Fill(signalArray[QmaxInd+1]);
        hQ4DistTime[i-1]->Fill(signalArray[QmaxInd+2]);
      }
    }*/

    }
/// END OF SAMPLING TEST
  }

  for (int i=0; i<nHist; ++i) {
    TF1 *QtotFit = new TF1("QtotFit","gaus");
    hQtotDist[i]->Fit("QtotFit");

    float mean = QtotFit->GetParameter(1);
    float meanErr = QtotFit->GetParError(1);
    float width = QtotFit->GetParameter(2);
    float widthErr = QtotFit->GetParError(2);
    float sigma = width/mean;
    float sigmaErr = sqrt(pow((widthErr/mean),2)+pow(((width*meanErr)/pow(mean,2)),2));

    gInSigVsOutSig->SetPoint(i,inputSigma[i],sigma*100);
    gInSigVsOutSig->SetPointError(i,0,sigmaErr*100);

    std::cout << std::endl << std::endl << "input: " << inputSigma[i] << "%\t" << "output: " << sigma*100 << "%" << std::endl;
  }


  TFile *g = TFile::Open(Form("%s",OutputPath), "recreate");
  g->WriteObject(gInSigVsOutSig, "InputVsOutput");
//  g->WriteObject(outputSigma, "outputSigma");
  g->WriteObject(hQmaxQtotRatio, "hQmaxQtotRatio");
  g->WriteObject(hQmaxQtotRatioVsDt, "hQmaxQtotRatioVsDt");
  for (int i=0; i<nHist; ++i) {
    g->WriteObject(hResult[i], Form("hResult_%f",inputSigma[i]));
    g->WriteObject(hResultQmax[i], Form("hResultQmax_%f",inputSigma[i]));
    g->WriteObject(hCharge[i], Form("hCharge_%f",inputSigma[i]));
    g->WriteObject(hQmaxDist[i], Form("hQmaxDist_%f",inputSigma[i]));
    g->WriteObject(hQtotDist[i], Form("hQtotDist_%f",inputSigma[i]));
  }

/*  g->WriteObject(htimeShift, "htimeShift");
  g->WriteObject(hQmaxVstimeShift, "hQmaxVstimeShift");
  g->WriteObject(hQtotVstimeShift, "hQtotVstimeShift");
  g->WriteObject(hQ1VstimeShift, "hQ1VstimeShift");
  g->WriteObject(hQ2VstimeShift, "hQ2VstimeShift");
  g->WriteObject(hQ3VstimeShift, "hQ3VstimeShift");
  g->WriteObject(hQ4VstimeShift, "hQ4VstimeShift");
  g->WriteObject(hQmaxInd, "hQmaxInd");*/
/*  for (int i = 0; i < nbins; ++i) {
    int timebin = (((i+1)/float(nbins))*TimeBinSize*1000);
    g->WriteObject(hResultTime[i], Form("hResultTime_%i_ns",timebin));
    g->WriteObject(hResultTimeQmax[i], Form("hResultTimeQmax_%i_ns",timebin));
    g->WriteObject(hQmaxDistTime[i], Form("hQmaxTime_%i_ns",timebin));
    g->WriteObject(hQtotDistTime[i], Form("hQtotTime_%i_ns",timebin));
    g->WriteObject(hQ1DistTime[i], Form("hQ1DistTime_%i_ns",timebin));
    g->WriteObject(hQ2DistTime[i], Form("hQ2DistTime_%i_ns",timebin));
    g->WriteObject(hQ3DistTime[i], Form("hQ3DistTime_%i_ns",timebin));
    g->WriteObject(hQ4DistTime[i], Form("hQ4DistTime_%i_ns",timebin));
  }*/
  delete g;
  delete hQmaxQtotRatio;
  delete hQmaxQtotRatioVsDt;
  for (int i=0; i<nHist; ++i) {
    delete hResult[i];
    delete hResultQmax[i];
    delete hCharge[i];
    delete hQmaxDist[i];
    delete hQtotDist[i];
  }

  std::cout << "====================================================================" << "\n \n" << "done" << "\n \n" << "==================================================================== \n";
}
