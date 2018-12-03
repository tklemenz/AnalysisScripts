//#include "TFile.h"
//#include "TGraphErrors.h"
//#include "TTree.h"

void analyseTree2017()

{
   TFile file ("/home/tom/myana/Results/GSI/anaOut_tree_true_true_true.root");

   TTree *dEdx = (TTree*)file.Get("TrackAna");

   dEdx->SetMarkerStyle(20);

   dEdx->SetAlias("relGain","(posPiQmax/59)*2000"); /// 54 for 1 GeV/c and 59 for 2 GeV/c


   const int nHV = 9;

//   int iHV[nHV] = {-1,0,1,2,3,4,5,6,7,8,9,10};
//   float isigmaFe[nHV] = {12.7,10.8,12.2,12.1,12.3,12.0,13.7,11.4,10.2,9.0,8.1,10.8};
//   float iIBF[nHV] = {0.65,1.12,0.98,1.12,1.6,0.92,0.52,0.68,0.82,1.16,2.54,1.12};

   int iHV[nHV] = {-1,0,1,5,6,7,8,9,10};
   float isigmaFe[nHV] = {12.7,10.8,12.2,13.7,11.4,10.2,9.0,8.1,10.8};
   float iIBF[nHV] = {0.65,1.12,0.98,0.52,0.68,0.82,1.16,2.54,1.12};


   int iBeamM[5] = {-1,-2,-3,-4,-5};
   float BetheBloch[6] = {0,1.04586,1.13406,1.19109,1.23241,1.26464};

   TH1F *dummyHistoBeam = new TH1F("hDummyBeam",";;;",100000,0,5);
   TH1F *dummyHistoSigma = new TH1F("hDummySigma",";;;",10000000,0,100);
   TH1F *dummyHistoIBF = new TH1F("hDummyIBF",";;;",10000000,0,100);

   dummyHistoBeam->GetXaxis()->SetTitle("beam momentum [GeV/#it{c}]");
   dummyHistoBeam->GetYaxis()->SetTitle("separation");
   dummyHistoBeam->GetYaxis()->SetTitleSize(30);
   dummyHistoBeam->GetYaxis()->SetLabelSize(27);
   dummyHistoBeam->GetYaxis()->SetTitleFont(43);
   dummyHistoBeam->GetYaxis()->SetLabelFont(43);
   dummyHistoBeam->GetXaxis()->SetTitleSize(30);
   dummyHistoBeam->GetXaxis()->SetLabelSize(27);
   dummyHistoBeam->GetXaxis()->SetTitleFont(43);
   dummyHistoBeam->GetXaxis()->SetLabelFont(43);

   dummyHistoSigma->GetXaxis()->SetTitle("#sigma(^{55}Fe) [%]");
   dummyHistoSigma->GetYaxis()->SetTitle("separation");
   dummyHistoSigma->GetYaxis()->SetTitleSize(30);
   dummyHistoSigma->GetYaxis()->SetLabelSize(27);
   dummyHistoSigma->GetYaxis()->SetTitleFont(43);
   dummyHistoSigma->GetYaxis()->SetLabelFont(43);
   dummyHistoSigma->GetXaxis()->SetTitleSize(30);
   dummyHistoSigma->GetXaxis()->SetLabelSize(27);
   dummyHistoSigma->GetXaxis()->SetTitleFont(43);
   dummyHistoSigma->GetXaxis()->SetLabelFont(43);

   dummyHistoIBF->GetXaxis()->SetTitle("IB [%]");
   dummyHistoIBF->GetYaxis()->SetTitle("separation");
   dummyHistoIBF->GetYaxis()->SetTitleSize(30);
   dummyHistoIBF->GetYaxis()->SetLabelSize(27);
   dummyHistoIBF->GetYaxis()->SetTitleFont(43);
   dummyHistoIBF->GetYaxis()->SetLabelFont(43);
   dummyHistoIBF->GetXaxis()->SetTitleSize(30);
   dummyHistoIBF->GetXaxis()->SetLabelSize(27);
   dummyHistoIBF->GetXaxis()->SetTitleFont(43);
   dummyHistoIBF->GetXaxis()->SetLabelFont(43);



   TGraphErrors *graphsSep[nHV][5];
   TGraphErrors *gSepVsBeamM = new TGraphErrors();
   gSepVsBeamM->SetMarkerStyle(20);
   gSepVsBeamM->SetMarkerSize(1.3);
   gSepVsBeamM->SetMarkerColor(kBlue+2);
   gSepVsBeamM->SetLineColor(kBlue+2);

   float SepHVBCasc[6] = {0,0,0,0,0,0};
   float SepHVBCascErr[6] = {0,0,0,0,0,0};
   double *SepHVBCascErrTemp1 = new double{0};
   double *SepHVBCascErrTemp2 = new double{0};
   double *SepHVBCascErrTemp3 = new double{0};
   double *SepHVBCascErrTemp4 = new double{0};
   double *SepHVBCascErrTemp5 = new double{0};

   TGraphErrors *gSepVsSigmaFe = new TGraphErrors();
   gSepVsSigmaFe->SetMarkerStyle(20);
   gSepVsSigmaFe->SetMarkerSize(1.3);
   gSepVsSigmaFe->SetMarkerColor(kBlue+2);
   gSepVsSigmaFe->SetLineColor(kBlue+2);

   TGraphErrors *gSepVsIBF = new TGraphErrors();
   gSepVsIBF->SetMarkerStyle(20);
   gSepVsIBF->SetMarkerSize(1.3);
   gSepVsIBF->SetMarkerColor(kBlue+2);
   gSepVsIBF->SetLineColor(kBlue+2);

//   float SepVsHV[nHV] = {0,0,0,0,0,0,0,0,0,0,0,0};
//   float SepVsHVErr[nHV] = {0,0,0,0,0,0,0,0,0,0,0,0};

   float SepVsHV[nHV] = {0,0,0,0,0,0,0,0,0};
   float SepVsHVErr[nHV] = {0,0,0,0,0,0,0,0,0};

   const double *ReadOut = new double{2000};
   const double *ReadOut0 = new double{2280};
   const double *ReadOut1 = new double{2000};
   const double *ReadOut2 = new double{1600};
//   const double *ReadOut3 = new double{1245};
//   const double *ReadOut4 = new double{765};
//   const double *ReadOut5 = new double{1995};
   const double *ReadOut3 = new double{2675};
   const double *ReadOut4 = new double{2960};
   const double *ReadOut5 = new double{2640};
   const double *ReadOut6 = new double{3040};
   const double *ReadOut7 = new double{3200};
   const double *ReadOut8 = new double{1700};

   double *SepVsHVErrTemp1 = new double{0};
   double *SepVsHVErrTemp2 = new double{0};
//   double *SepVsHVErrTemp3 = new double{0};
//   double *SepVsHVErrTemp4 = new double{0};
//   double *SepVsHVErrTemp5 = new double{0};
   double *SepVsHVErrTemp3 = new double{0};
   double *SepVsHVErrTemp4 = new double{0};
   double *SepVsHVErrTemp5 = new double{0};
   double *SepVsHVErrTemp6 = new double{0};
   double *SepVsHVErrTemp7 = new double{0};
   double *SepVsHVErrTemp8 = new double{0};
   double *SepVsHVErrTemp0 = new double{0};

   TF1 p1("p1","pol1");

   for (int iiHV = 0; iiHV<nHV; ++iiHV) {
     for (int iiBeamM = 0; iiBeamM<5; ++iiBeamM) {
       dEdx->Draw("SepPowerQtot:SepPowerQtotErr:relGain",Form("beamMomentum == %i && HVSetting == %i && abs(SepPowerQtotErr/SepPowerQtot)<0.25",iBeamM[iiBeamM], iHV[iiHV]));
       graphsSep[iiHV][iiBeamM] = new TGraphErrors(dEdx->GetSelectedRows(),dEdx->GetV3(), dEdx->GetV1(), 0, dEdx->GetV2());
       graphsSep[iiHV][iiBeamM]->SetTitle(Form("HV %i, %i GeV/#it{c};#LTd#it{E}/d#it{x}_{Qmax}#GT^{#pi};separation",iHV[iiHV],iBeamM[iiBeamM]));
       graphsSep[iiHV][iiBeamM]->SetMarkerStyle(20);
       if (graphsSep[iiHV][iiBeamM]->GetN() != 0) {
         TFitResultPtr graphsSepFitRes = graphsSep[iiHV][iiBeamM]->Fit(&p1, "S");


       /// Separation Vs BeamMommentum
       if (iHV[iiHV] == 10) {
         if (iBeamM[iiBeamM] == -1) {
           SepHVBCasc[1] = p1.Eval(1550);
           const double *GainTemp1 = new double{1550};
           graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, GainTemp1, SepHVBCascErrTemp1, 0.683, false);
           SepHVBCascErr[1] = SepHVBCascErrTemp1[0];
         }
         else if (iBeamM[iiBeamM] == -2) {
           SepHVBCasc[2] = p1.Eval(1550*(BetheBloch[2]/BetheBloch[1]));
           const double *GainTemp2 = new double{1550*(BetheBloch[2]/BetheBloch[1])};
           graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, GainTemp2, SepHVBCascErrTemp2, 0.683, false);
           SepHVBCascErr[2] = SepHVBCascErrTemp2[0];
         }
         else if (iBeamM[iiBeamM] == -3) {
           SepHVBCasc[3] = p1.Eval(1550*(BetheBloch[3]/BetheBloch[1]));
           //SepHVBCascErr[3] = sqrt(pow(1550*(BetheBloch[3]/BetheBloch[1])*p1.GetParError(0),2)+pow(p1.GetParError(1),2));
           const double *GainTemp3 = new double{1550*(BetheBloch[3]/BetheBloch[1])};
           graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, GainTemp3, SepHVBCascErrTemp3, 0.683, false);
           SepHVBCascErr[3] = SepHVBCascErrTemp3[0];
         }
         else if (iBeamM[iiBeamM] == -4) {
           SepHVBCasc[4] = p1.Eval(1550*(BetheBloch[4]/BetheBloch[1]));
           //SepHVBCascErr[4] = sqrt(pow(1550*(BetheBloch[4]/BetheBloch[1])*p1.GetParError(0),2)+pow(p1.GetParError(1),2));
           const double *GainTemp4 = new double{1550*(BetheBloch[4]/BetheBloch[1])};
           graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, GainTemp4, SepHVBCascErrTemp4, 0.683, false);
           SepHVBCascErr[4] = SepHVBCascErrTemp4[0];
         }
         else if (iBeamM[iiBeamM] == -5) {
           SepHVBCasc[5] = p1.Eval(1550*(BetheBloch[5]/BetheBloch[1]));
           //SepHVBCascErr[5] = sqrt(pow(1550*(BetheBloch[5]/BetheBloch[1])*p1.GetParError(0),2)+pow(p1.GetParError(1),2));
           const double *GainTemp5 = new double{1550*(BetheBloch[5]/BetheBloch[1])};
           graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, GainTemp5, SepHVBCascErrTemp5, 0.683, false);
           SepHVBCascErr[5] = SepHVBCascErrTemp5[0];
         }
       }
       else if (iHV[iiHV] == -1 && iBeamM[iiBeamM] == -2) {
         SepVsHV[0] = p1.Eval(2280);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut0, SepVsHVErrTemp0, 0.683, false);
         SepVsHVErr[0] = SepVsHVErrTemp0[0];
       }
       else if (iHV[iiHV] == 0 && iBeamM[iiBeamM] == -2) {
         SepVsHV[1] = p1.Eval(2000);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut1, SepVsHVErrTemp1, 0.683, false);
         SepVsHVErr[1] = SepVsHVErrTemp1[0];
       }
       else if (iHV[iiHV] == 1 && iBeamM[iiBeamM] == -2) {
         SepVsHV[2] = p1.Eval(1600);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut2, SepVsHVErrTemp2, 0.683, false);
         SepVsHVErr[2] = SepVsHVErrTemp2[0];
       }
/*       else if (iHV[iiHV] == 2 && iBeamM[iiBeamM] == -2) {
         SepVsHV[3] = p1.Eval(1245);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut3, SepVsHVErrTemp3, 0.683, false);
         SepVsHVErr[3] = SepVsHVErrTemp3[0];
       }
       else if (iHV[iiHV] == 3 && iBeamM[iiBeamM] == -2) {
         SepVsHV[4] = p1.Eval(765);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut4, SepVsHVErrTemp4, 0.683, false);
         SepVsHVErr[4] = SepVsHVErrTemp4[0];
       }
       else if (iHV[iiHV] == 4 && iBeamM[iiBeamM] == -2) {
         SepVsHV[5] = p1.Eval(1995);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut5, SepVsHVErrTemp5, 0.683, false);
         SepVsHVErr[5] = SepVsHVErrTemp5[0];
       } */
       else if (iHV[iiHV] == 5 && iBeamM[iiBeamM] == -2) {
         SepVsHV[3] = p1.Eval(2675);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut3, SepVsHVErrTemp3, 0.683, false);
         SepVsHVErr[3] = SepVsHVErrTemp3[0];
       }
       else if (iHV[iiHV] == 6 && iBeamM[iiBeamM] == -2) {
         SepVsHV[4] = p1.Eval(2960);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut4, SepVsHVErrTemp4, 0.683, false);
         SepVsHVErr[4] = SepVsHVErrTemp4[0];
       }
       else if (iHV[iiHV] == 7 && iBeamM[iiBeamM] == -2) {
         SepVsHV[5] = p1.Eval(2640);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut5, SepVsHVErrTemp5, 0.683, false);
         SepVsHVErr[5] = SepVsHVErrTemp5[0];
       }
       else if (iHV[iiHV] == 8 && iBeamM[iiBeamM] == -2) {
         SepVsHV[6] = p1.Eval(3040);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut6, SepVsHVErrTemp6, 0.683, false);
         SepVsHVErr[6] = SepVsHVErrTemp6[0];
       }
       else if (iHV[iiHV] == 9 && iBeamM[iiBeamM] == -2) {
         SepVsHV[7] = p1.Eval(3200);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut7, SepVsHVErrTemp7, 0.683, false);
         SepVsHVErr[7] = SepVsHVErrTemp7[0];
       }
       else if (iHV[iiHV] == 10 && iBeamM[iiBeamM] == -2) {
         SepVsHV[8] = p1.Eval(1700);
         graphsSepFitRes->GetConfidenceIntervals(1, 1, 1, ReadOut8, SepVsHVErrTemp8, 0.683, false);
         SepVsHVErr[8] = SepVsHVErrTemp8[0];
       }
     }
     }
   }

  for (int i = 0; i<2; ++i) {
    if (SepHVBCasc[i+1] != 0) {
      gSepVsBeamM->SetPoint(i,iBeamM[i]*(-1),SepHVBCasc[i+1]);
      gSepVsBeamM->SetPointError(i,0,SepHVBCascErr[i+1]);
    }
  }
  gSepVsBeamM->GetXaxis()->SetTitle("beam momentum [GeV/#it{c}]");
  gSepVsBeamM->GetYaxis()->SetTitle("separation");
  gSepVsBeamM->GetYaxis()->SetTitleSize(30);
  gSepVsBeamM->GetYaxis()->SetLabelSize(27);
  gSepVsBeamM->GetYaxis()->SetTitleFont(43);
  gSepVsBeamM->GetYaxis()->SetLabelFont(43);

  gSepVsBeamM->GetXaxis()->SetTitleSize(30);
  gSepVsBeamM->GetXaxis()->SetLabelSize(27);
  gSepVsBeamM->GetXaxis()->SetTitleFont(43);
  gSepVsBeamM->GetXaxis()->SetLabelFont(43);

  for (int i = 0; i<9; ++i) {
    if (SepVsHV[i] != 0) {
      gSepVsSigmaFe->SetPoint(i,isigmaFe[i],SepVsHV[i]);
      gSepVsSigmaFe->SetPointError(i,0,SepVsHVErr[i]);

      gSepVsIBF->SetPoint(i,iIBF[i],SepVsHV[i]);
      gSepVsIBF->SetPointError(i,0,SepVsHVErr[i]);
    }
  }

  gSepVsSigmaFe->GetXaxis()->SetTitle("#sigma(^{55}Fe) [%]");
  gSepVsSigmaFe->GetYaxis()->SetTitle("separation");
  gSepVsSigmaFe->GetYaxis()->SetTitleSize(30);
  gSepVsSigmaFe->GetYaxis()->SetLabelSize(27);
  gSepVsSigmaFe->GetYaxis()->SetTitleFont(43);
  gSepVsSigmaFe->GetYaxis()->SetLabelFont(43);

  gSepVsSigmaFe->GetXaxis()->SetTitleSize(30);
  gSepVsSigmaFe->GetXaxis()->SetLabelSize(27);
  gSepVsSigmaFe->GetXaxis()->SetTitleFont(43);
  gSepVsSigmaFe->GetXaxis()->SetLabelFont(43);

  gSepVsIBF->GetXaxis()->SetTitle("IB [%]");
  gSepVsIBF->GetYaxis()->SetTitle("separation");
  gSepVsIBF->GetYaxis()->SetTitleSize(30);
  gSepVsIBF->GetYaxis()->SetLabelSize(27);
  gSepVsIBF->GetYaxis()->SetTitleFont(43);
  gSepVsIBF->GetYaxis()->SetLabelFont(43);

  gSepVsIBF->GetXaxis()->SetTitleSize(30);
  gSepVsIBF->GetXaxis()->SetLabelSize(27);
  gSepVsIBF->GetXaxis()->SetTitleFont(43);
  gSepVsIBF->GetXaxis()->SetLabelFont(43);

   TFile *results = new TFile("/home/tom/myana/Results/analyseTree2017/resultsFile_only_gain2000.root","recreate");
   results->WriteObject(gSepVsBeamM, "SepPowerVsBeam[HV10]");
   results->WriteObject(gSepVsSigmaFe, "SepPowerVsSigmaFe");
   results->WriteObject(gSepVsIBF, "SepPowerVsIBF");
   results->WriteObject(dummyHistoBeam, "dummyHistoBeam");
   results->WriteObject(dummyHistoSigma, "dummyHistoSigma");
   results->WriteObject(dummyHistoIBF, "dummyHistoIBF");
   for (int i = 0; i<nHV; ++i) {
     for (int j = 0; j<5; ++j) {
       if (graphsSep[i][j]->GetN() != 0) {
         results->WriteObject(graphsSep[i][j],Form("graphsSep[HV%i][%i GeV/c]",iHV[i],iBeamM[j]));
       }
     }
   }
   delete results;
}
