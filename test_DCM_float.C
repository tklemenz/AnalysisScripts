void testDCM(TString deadchannelmap)
{
  using namespace o2::TPC;
  TFile f(deadchannelmap);
  gROOT->cd();

  static const Mapper& mapper = Mapper::instance();
  ROC roc(0);

  CalDet<float> *dcm = nullptr;
  f.GetObject("DeadChannelMap", dcm);

  CalDet<float> *dcmVis = new CalDet<float>(PadSubset::ROC);
  CalArray<float> dcmVisualisation = dcmVis->getCalArray(roc);

  for (int irow=0; irow<63; ++irow) {
    int npads = mapper.getNumberOfPadsInRowROC(roc, irow);
    for (int ipad=0; ipad<npads; ++ipad) {
      bool val = dcm->getValue(roc, irow, ipad);
      if (val == 1) {
        dcmVisualisation.setValue(irow, ipad, 10);
      }
      else {
        dcmVisualisation.setValue(irow, ipad, 5);
      }
    }
  }

  auto visplot = Painter::getHistogram2D(dcmVisualisation);
  auto cdcmVis = new TCanvas ("cdcmVis","DCM Visualization");
  visplot->Draw("colz");
}

