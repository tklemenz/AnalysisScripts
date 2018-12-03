#include "TPCSimulation/Detector.h"
#include "TPCBase/ParameterGas.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>

using std::cout;
using std::endl;

using namespace o2::TPC;

float RealGain(float posPiMeasured, float posPiExpected, float GainExpected)
{
  return GainExpected * posPiMeasured/posPiExpected;
}

float RealGainFromNoise(float posPiMeasured, float NElectrons)
{
  int ePerADC = 670;

  return posPiMeasured * ePerADC/NElectrons;
}

float BetheBlochAleph(float bg, float kp1, float kp2, float kp3, float kp4, float kp5)
{
  float beta = bg/std::sqrt(1.+ bg*bg);
  float aa = std::pow(beta,kp4);
  float bb = std::pow(1./bg,kp5);
  bb=std::log(kp3+bb);

  return (kp2-aa-bb)*kp1/aa;
}

float BetheBloch(std::string Particle = "", int BeamMomentum = 1, float gain = 2000)
{
  const static ParameterGas &gasParam = ParameterGas::defaultInstance();

  float Mpi = 139.57018*1e6;               // eV/c^2
  float Mele = 0.5109989461*1e6;           // eV/c^2
  float beamMomentum = BeamMomentum*1e9;   // eV/c

  float NElectronsPad = 0;
  float m = 0;
  float betaGamma = 0;
  float ConversionGain = 2*1e13;           // Volt/Coulomb
  float ChargeTot = 0;
  float ChargeEle = 1.6021766208*1e-19;    // Coulomb
  float ADCChannels = 1024;
  float ADCRange = 2.2;                    // Volt
  float ADCValue = 0;
  float NeCON = 36.1;

  if (Particle == "pi") {
    m = Mpi;
  }
  else if (Particle == "e") {
    m = Mele;
  }
  else {
    cout << endl << endl << "Choose either pi or e as particle" << endl << endl;
    return 0;
  }

  betaGamma = beamMomentum/m;
  NElectronsPad = gain * 0.7513 * NeCON * BetheBlochAleph(betaGamma, gasParam.getBetheBlochParam(0), gasParam.getBetheBlochParam(1), gasParam.getBetheBlochParam(2), gasParam.getBetheBlochParam(3), gasParam.getBetheBlochParam(4));
//  NElectronsPad = gain * 0.7513 * gasParam.getNprim() * BetheBlochAleph(betaGamma, gasParam.getBetheBlochParam(0), gasParam.getBetheBlochParam(1), gasParam.getBetheBlochParam(2), gasParam.getBetheBlochParam(3), gasParam.getBetheBlochParam(4));

  ChargeTot = NElectronsPad * ChargeEle;
  ADCValue = ConversionGain * ChargeTot * 0.64/(ADCRange/ADCChannels);

  return ADCValue;
}

float nElectronsNoGain(std::string Particle = "", int BeamMomentum = 1)
{
  const static ParameterGas &gasParam = ParameterGas::defaultInstance();

  float Mpi = 139.57018*1e6;               // eV/c^2
  float Mele = 0.5109989461*1e6;           // eV/c^2
  float beamMomentum = BeamMomentum*1e9;   // eV/c

  float NElectronsPad = 0;
  float m = 0;
  float betaGamma = 0;
  float NeCON = 36.1;
  float ArCO = 74.9;

  if (Particle == "pi") {
    m = Mpi;
  }
  else if (Particle == "e") {
    m = Mele;
  }
  else {
    cout << endl << endl << "Choose either pi or e as particle" << endl << endl;
    return 0;
  }

  betaGamma = beamMomentum/m;
  NElectronsPad = 0.7513 * NeCON * BetheBlochAleph(betaGamma, gasParam.getBetheBlochParam(0), gasParam.getBetheBlochParam(1), gasParam.getBetheBlochParam(2), gasParam.getBetheBlochParam(3), gasParam.getBetheBlochParam(4));

  return NElectronsPad;
}

void HVRealGain(TString filename)
{
  TFile *TreeFile = TFile::Open(filename.Data());
  TTree *tree = (TTree*)gDirectory->Get("TrackAna");

  float posPiQmax = 0;
  int beamMomentum = 0;
  int HVSetting = 0;

  TGraphErrors *gRealGain = new TGraphErrors();
  TGraph *gGainFrankfurt = new TGraph();
  float realGains[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
  float realGainsErr[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
  float gainFrankfurt[12] = {2000,2000,2000,1500,1000,3000,2000,2000,2000,2000,2000,2000};

  gRealGain->SetMarkerStyle(20);
  gRealGain->SetMarkerColor(kRed);
  gRealGain->SetLineColor(kRed);
  gRealGain->GetXaxis()->SetTitle("HV setting");
  gRealGain->GetYaxis()->SetTitle("gain");
  gGainFrankfurt->SetMarkerStyle(20);
  gGainFrankfurt->SetMarkerColor(kBlue);
  gGainFrankfurt->SetLineColor(kBlue);

  TH1F *hposPiQmax_01_2= new TH1F("hposPiQmax_01_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_0_2= new TH1F("hposPiQmax_0_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_1_2= new TH1F("hposPiQmax_1_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_2_2= new TH1F("hposPiQmax_2_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_3_2= new TH1F("hposPiQmax_3_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_4_2= new TH1F("hposPiQmax_4_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_5_2= new TH1F("hposPiQmax_5_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_6_2= new TH1F("hposPiQmax_6_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_7_2= new TH1F("hposPiQmax_7_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_8_2= new TH1F("hposPiQmax_8_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_9_2= new TH1F("hposPiQmax_9_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_10_1= new TH1F("hposPiQmax_10_1", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_10_2= new TH1F("hposPiQmax_10_2", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_10_3= new TH1F("hposPiQmax_10_3", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_10_4= new TH1F("hposPiQmax_10_4", "; posPiQmax; Counts", 600,0,120);
  TH1F *hposPiQmax_10_5= new TH1F("hposPiQmax_10_5", "; posPiQmax; Counts", 600,0,120);

  tree->SetBranchAddress("posPiQmax", &posPiQmax);
  tree->SetBranchAddress("beamMomentum", &beamMomentum);
  tree->SetBranchAddress("HVSetting", &HVSetting);

  for (int iRun = 0; iRun<tree->GetEntriesFast(); ++iRun) {
    tree->GetEntry(iRun);

    if (beamMomentum == -2 && HVSetting == -1) {
      hposPiQmax_01_2->Fill(posPiQmax);
    }
    else if(beamMomentum == -2 && HVSetting == 0) {
      hposPiQmax_0_2->Fill(posPiQmax);
    }
    else if(beamMomentum == -2 && HVSetting == 1) {
      hposPiQmax_1_2->Fill(posPiQmax);
    }
    else if(beamMomentum == -2 && HVSetting == 2) {
      hposPiQmax_2_2->Fill(posPiQmax);
    }
    else if(beamMomentum == -2 && HVSetting == 3) {
      hposPiQmax_3_2->Fill(posPiQmax);
    }
    else if(beamMomentum == -2 && HVSetting == 4) {
      hposPiQmax_4_2->Fill(posPiQmax);
    }
    else if(beamMomentum == -2 && HVSetting == 5) {
      hposPiQmax_5_2->Fill(posPiQmax);
    }
    else if(beamMomentum == -2 && HVSetting == 6) {
      hposPiQmax_6_2->Fill(posPiQmax);
    }
    else if(beamMomentum == -2 && HVSetting == 7) {
      hposPiQmax_7_2->Fill(posPiQmax);
    }
    else if(beamMomentum == -2 && HVSetting == 8) {
      hposPiQmax_8_2->Fill(posPiQmax);
    }
    else if(beamMomentum == -2 && HVSetting == 9) {
      hposPiQmax_9_2->Fill(posPiQmax);
    }
    else if(beamMomentum == -1 && HVSetting == 10) {
      hposPiQmax_10_1->Fill(posPiQmax);
    }
    else if(beamMomentum == -2 && HVSetting == 10) {
      hposPiQmax_10_2->Fill(posPiQmax);
    }
    //else if(beamMomentum == -2 && HVSetting == 10) {     /// do not treat HV 10 seperate from HVB
      //hposPiQmax_0_2->Fill(posPiQmax);
    //}
    else if(beamMomentum == -3 && HVSetting == 10) {
      hposPiQmax_10_3->Fill(posPiQmax);
    }
    else if(beamMomentum == -4 && HVSetting == 10) {
      hposPiQmax_10_4->Fill(posPiQmax);
    }
    else if(beamMomentum == -5 && HVSetting == 10) {
      hposPiQmax_10_5->Fill(posPiQmax);
    }
    else {
      cout << endl << endl << "Something went wrong!!" << endl << endl << endl;
    }

  } // end of run loop

  float gain_01_2, gain_0_2, gain_1_2, gain_2_2, gain_3_2, gain_4_2, gain_5_2, gain_6_2, gain_7_2,
      gain_8_2, gain_9_2, gain_10_1, gain_10_2, gain_10_3, gain_10_4, gain_10_5;

  gain_01_2 = gain_0_2 = gain_1_2 = gain_2_2 = gain_3_2 = gain_4_2 = gain_5_2 = gain_6_2 = gain_7_2 =
        gain_8_2 = gain_9_2 = gain_10_1 = gain_10_2 = gain_10_3 = gain_10_4 = gain_10_5 = 0.;

  TTree *RealGains = new TTree("RealGains", "RealGain");
  RealGains->Branch("HVADivider_2GeV", &gain_01_2, "gain_01_2/F");
  RealGains->Branch("HVBDivider_2GeV", &gain_0_2, "gain_0_2/F");
  RealGains->Branch("HV1_2GeV", &gain_1_2, "gain_1_2/F");
  RealGains->Branch("HV2_2GeV", &gain_2_2, "gain_2_2/F");
  RealGains->Branch("HV3_2GeV", &gain_3_2, "gain_3_2/F");
  RealGains->Branch("HV4_2GeV", &gain_4_2, "gain_4_2/F");
  RealGains->Branch("HV5_2GeV", &gain_5_2, "gain_5_2/F");
  RealGains->Branch("HV6_2GeV", &gain_6_2, "gain_6_2/F");
  RealGains->Branch("HV7_2GeV", &gain_7_2, "gain_7_2/F");
  RealGains->Branch("HV8_2GeV", &gain_8_2, "gain_8_2/F");
  RealGains->Branch("HV9_2GeV", &gain_9_2, "gain_9_2/F");
  RealGains->Branch("HVBCascaded_2GeV", &gain_10_2, "gain_10_2/F");
  RealGains->Branch("HVBCascaded_1GeV", &gain_10_1, "gain_10_1/F");
  RealGains->Branch("HVBCascaded_3GeV", &gain_10_3, "gain_10_3/F");
  RealGains->Branch("HVBCascaded_4GeV", &gain_10_4, "gain_10_4/F");
  RealGains->Branch("HVBCascaded_5GeV", &gain_10_5, "gain_10_5/F");

  gain_01_2 = RealGain(hposPiQmax_01_2->GetMean(), BetheBloch("pi",2,2000), 2000);
  gain_0_2 = RealGain(hposPiQmax_0_2->GetMean(), BetheBloch("pi",2,2000), 2000);
  gain_1_2 = RealGain(hposPiQmax_1_2->GetMean(), BetheBloch("pi",2,2000), 2000);
  gain_2_2 = RealGain(hposPiQmax_2_2->GetMean(), BetheBloch("pi",2,1500), 1500);
  gain_3_2 = RealGain(hposPiQmax_3_2->GetMean(), BetheBloch("pi",2,1000), 1000);
  gain_4_2 = RealGain(hposPiQmax_4_2->GetMean(), BetheBloch("pi",2,3000), 3000);
  gain_5_2 = RealGain(hposPiQmax_5_2->GetMean(), BetheBloch("pi",2,2000), 2000);
  gain_6_2 = RealGain(hposPiQmax_6_2->GetMean(), BetheBloch("pi",2,2000), 2000);
  gain_7_2 = RealGain(hposPiQmax_7_2->GetMean(), BetheBloch("pi",2,2000), 2000);
  gain_8_2 = RealGain(hposPiQmax_8_2->GetMean(), BetheBloch("pi",2,2000), 2000);
  gain_9_2 = RealGain(hposPiQmax_9_2->GetMean(), BetheBloch("pi",2,2000), 2000);
  gain_10_2 = RealGain(hposPiQmax_10_2->GetMean(), BetheBloch("pi",2,2000), 2000);
  gain_10_1 = RealGain(hposPiQmax_10_1->GetMean(), BetheBloch("pi",1,2000), 2000);
  gain_10_3 = RealGain(hposPiQmax_10_3->GetMean(), BetheBloch("pi",3,2000), 2000);
  gain_10_4 = RealGain(hposPiQmax_10_4->GetMean(), BetheBloch("pi",4,2000), 2000);
  gain_10_5 = RealGain(hposPiQmax_10_5->GetMean(), BetheBloch("pi",5,2000), 2000);

  realGains[0] = gain_01_2;
  realGains[1] = gain_0_2;
  realGains[2] = gain_1_2;
  realGains[3] = gain_2_2;
  realGains[4] = gain_3_2;
  realGains[5] = gain_4_2;
  realGains[6] = gain_5_2;
  realGains[7] = gain_6_2;
  realGains[8] = gain_7_2;
  realGains[9] = gain_8_2;
  realGains[10] = gain_9_2;
  realGains[11] = gain_10_2;

  realGainsErr[0] = hposPiQmax_01_2->GetRMS()*(2000/BetheBloch("pi",2,2000));
  realGainsErr[1] = hposPiQmax_0_2->GetRMS()*(2000/BetheBloch("pi",2,2000));
  realGainsErr[2] = hposPiQmax_1_2->GetRMS()*(2000/BetheBloch("pi",2,2000));
  realGainsErr[3] = hposPiQmax_2_2->GetRMS()*(1500/BetheBloch("pi",2,1500));
  realGainsErr[4] = hposPiQmax_3_2->GetRMS()*(1000/BetheBloch("pi",2,1000));
  realGainsErr[5] = hposPiQmax_4_2->GetRMS()*(3000/BetheBloch("pi",2,3000));
  realGainsErr[6] = hposPiQmax_5_2->GetRMS()*(2000/BetheBloch("pi",2,2000));
  realGainsErr[7] = hposPiQmax_6_2->GetRMS()*(2000/BetheBloch("pi",2,2000));
  realGainsErr[8] = hposPiQmax_7_2->GetRMS()*(2000/BetheBloch("pi",2,2000));
  realGainsErr[9] = hposPiQmax_8_2->GetRMS()*(2000/BetheBloch("pi",2,2000));
  realGainsErr[10] = hposPiQmax_9_2->GetRMS()*(2000/BetheBloch("pi",2,2000));
  realGainsErr[11] = hposPiQmax_10_2->GetRMS()*(2000/BetheBloch("pi",2,2000));

  for (int i = 0; i<12; ++i) {                      /// 11 if HV 10 not seperate from HV B
    gRealGain->SetPoint(i,-1+i,realGains[i]);
    gRealGain->SetPointError(i,0,realGainsErr[i]);
    gGainFrankfurt->SetPoint(i,-1+i,gainFrankfurt[i]);
  }



/*  gain_01_2 = RealGainFromNoise(hposPiQmax_01_2->GetMean(), nElectronsNoGain("pi",2));
  gain_0_2 = RealGainFromNoise(hposPiQmax_0_2->GetMean(), nElectronsNoGain("pi",2));
  gain_1_2 = RealGainFromNoise(hposPiQmax_1_2->GetMean(), nElectronsNoGain("pi",2));
  gain_2_2 = RealGainFromNoise(hposPiQmax_2_2->GetMean(), nElectronsNoGain("pi",2));
  gain_3_2 = RealGainFromNoise(hposPiQmax_3_2->GetMean(), nElectronsNoGain("pi",2));
  gain_4_2 = RealGainFromNoise(hposPiQmax_4_2->GetMean(), nElectronsNoGain("pi",2));
  gain_5_2 = RealGainFromNoise(hposPiQmax_5_2->GetMean(), nElectronsNoGain("pi",2));
  gain_6_2 = RealGainFromNoise(hposPiQmax_6_2->GetMean(), nElectronsNoGain("pi",2));
  gain_7_2 = RealGainFromNoise(hposPiQmax_7_2->GetMean(), nElectronsNoGain("pi",2));
  gain_8_2 = RealGainFromNoise(hposPiQmax_8_2->GetMean(), nElectronsNoGain("pi",2));
  gain_9_2 = RealGainFromNoise(hposPiQmax_9_2->GetMean(), nElectronsNoGain("pi",2));
  gain_10_2 = RealGainFromNoise(hposPiQmax_10_2->GetMean(), nElectronsNoGain("pi",2));
  gain_10_1 = RealGainFromNoise(hposPiQmax_10_1->GetMean(), nElectronsNoGain("pi",1));
  gain_10_3 = RealGainFromNoise(hposPiQmax_10_3->GetMean(), nElectronsNoGain("pi",3));
  gain_10_4 = RealGainFromNoise(hposPiQmax_10_4->GetMean(), nElectronsNoGain("pi",4));
  gain_10_5 = RealGainFromNoise(hposPiQmax_10_5->GetMean(), nElectronsNoGain("pi",5));
*/
  RealGains->Fill();

  TCanvas *gainComp = new TCanvas();
  gRealGain->GetYaxis()->SetTitleSize(24);
  gRealGain->GetYaxis()->SetLabelSize(21);
  gRealGain->GetYaxis()->SetTitleFont(43);
  gRealGain->GetYaxis()->SetLabelFont(43);
  gRealGain->GetYaxis()->SetTitleOffset(1.1);

  gRealGain->GetXaxis()->SetTitleSize(24);
  gRealGain->GetXaxis()->SetLabelSize(21);
  gRealGain->GetXaxis()->SetTitleFont(43);
  gRealGain->GetXaxis()->SetLabelFont(43);
  gRealGain->GetXaxis()->SetTitleOffset(1.1);

  gGainFrankfurt->GetYaxis()->SetTitleSize(24);
  gGainFrankfurt->GetYaxis()->SetLabelSize(21);
  gGainFrankfurt->GetYaxis()->SetTitleFont(43);
  gGainFrankfurt->GetYaxis()->SetLabelFont(43);
  gGainFrankfurt->GetYaxis()->SetTitleOffset(1.1);

  gGainFrankfurt->GetXaxis()->SetTitleSize(24);
  gGainFrankfurt->GetXaxis()->SetLabelSize(21);
  gGainFrankfurt->GetXaxis()->SetTitleFont(43);
  gGainFrankfurt->GetXaxis()->SetLabelFont(43);
  gGainFrankfurt->GetXaxis()->SetTitleOffset(1.1);

  gRealGain->Draw("ap");
  gGainFrankfurt->Draw("same p");
  TLegend *leg = new TLegend(0.75,0.75,0.9,0.9);
  leg->AddEntry(gRealGain,"test beam","p");
  leg->AddEntry(gGainFrankfurt,"nominal","p");
  leg->Draw("same");


  TFile *f = new TFile("/home/tom/myana/Results/HVRealGain/HVRealGain_exp_gain_correct.root", "recreate");

  f->WriteObject(RealGains, "RealGains");
  f->WriteObject(gainComp, "gainComparison");
  f->WriteObject(gRealGain, "RealGainVsHVSetting");
  f->WriteObject(gGainFrankfurt, "FrankfurtGainVsHVSetting");
  f->WriteObject(hposPiQmax_01_2, "hposPiQmax_HVA_2GeV");
  f->WriteObject(hposPiQmax_0_2, "hposPiQmax_HVB_2GeV");
  f->WriteObject(hposPiQmax_1_2, "hposPiQmax_HV1_2GeV");
  f->WriteObject(hposPiQmax_2_2, "hposPiQmax_HV2_2GeV");
  f->WriteObject(hposPiQmax_3_2, "hposPiQmax_HV3_2GeV");
  f->WriteObject(hposPiQmax_4_2, "hposPiQmax_HV4_2GeV");
  f->WriteObject(hposPiQmax_5_2, "hposPiQmax_HV5_2GeV");
  f->WriteObject(hposPiQmax_6_2, "hposPiQmax_HV6_2GeV");
  f->WriteObject(hposPiQmax_7_2, "hposPiQmax_HV7_2GeV");
  f->WriteObject(hposPiQmax_8_2, "hposPiQmax_HV8_2GeV");
  f->WriteObject(hposPiQmax_9_2, "hposPiQmax_HV9_2GeV");
  f->WriteObject(hposPiQmax_10_1, "hposPiQmax_HVBCas_1GeV");
  f->WriteObject(hposPiQmax_10_2, "hposPiQmax_HVBCas_2GeV");
  f->WriteObject(hposPiQmax_10_3, "hposPiQmax_HVBCas_3GeV");
  f->WriteObject(hposPiQmax_10_4, "hposPiQmax_HVBCas_4GeV");
  f->WriteObject(hposPiQmax_10_5, "hposPiQmax_HVBCas_5GeV");

  delete f;
}
