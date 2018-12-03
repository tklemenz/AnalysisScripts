int makeniceplots(TString Filename, TString HistoName)
{
//  using namespace o2::TPC;

//  gStyle->SetOptStat(0);

  TFile f(Filename);

  TH2D *Histo = nullptr;
  f.GetObject(HistoName, Histo);
//  if (f.IsOpen()) {
//    std::cout<<"I bims 1 file"<<endl;
//  }

  gPiondEdxMean->GetYaxis()->SetTitleSize(24);
  gPiondEdxMean->GetYaxis()->SetLabelSize(21);
  gPiondEdxMean->GetYaxis()->SetTitleFont(43);
  gPiondEdxMean->GetYaxis()->SetLabelFont(43);
  gPiondEdxMean->GetYaxis()->SetTitleOffset(1.1);

  gPiondEdxMean->GetXaxis()->SetTitleSize(24);
  gPiondEdxMean->GetXaxis()->SetLabelSize(21);
  gPiondEdxMean->GetXaxis()->SetTitleFont(43);
  gPiondEdxMean->GetXaxis()->SetLabelFont(43);
  gPiondEdxMean->GetXaxis()->SetTitleOffset(1.1);

//  TCanvas *Canvas = new TCanvas();
//  Canvas->SetRightMargin(0.16);
  Histo->Draw("colz");
  return 0;
}
