void CalDetTest(TString File)
{
  using namespace o2::TPC;
  TFile f(File);
  gROOT->cd();

  ROC roc(0);

  CalDet<float> *gainmappi = nullptr, *gainmapele = nullptr;
  f.GetObject("GainMapPi", gainmappi);
  f.GetObject("GainMapEle", gainmapele);

  if (gainmappi == nullptr) {cout<<endl<<"BIG FAIL"<<endl;}

  CalDet<float> *GainMapPi = new CalDet<float>(PadSubset::ROC);
  CalDet<float> *GainMapEle = new CalDet<float>(PadSubset::ROC);

  CalArray<float>& GainMapEleArr = GainMapEle->getCalArray(roc);
  CalArray<float>& GainMapPiArr = GainMapPi->getCalArray(roc);

  GainMapEleArr = gainmapele->getCalArray(roc);
  GainMapPiArr = gainmappi->getCalArray(roc);

  /*(if (GainMapPi == nullptr) {cout<<endl<<"FAIL"<<endl;}
  else if (GainMapPi == gainmappi) {cout<<endl<<"SUCCESS"<<endl;}
  else {cout<<endl<<"WHAT HAPPENED?"<<endl;}

  auto real = Painter::getHistogram2D(gainmappi->getCalArray(roc));
  auto copy = Painter::getHistogram2D(GainMapPiArr);

  auto creal = new TCanvas("creal", "real");
  real->Draw("colz");
  auto ccopy = new TCanvas("ccopy", "copy");
  copy->Draw("colz");*/


  cout<<"done"<<endl;
}

