void getDeadChannelMap(TString pedestalFile0, TString pedestalFile1)
{
  using namespace o2::TPC;
  TFile file0(pedestalFile0);
  TFile file1(pedestalFile1);
  gROOT->cd();
  static const Mapper& mapper = Mapper::instance();

  CalDet<float> *pedestal0=nullptr, *noise0=nullptr, *pedestal1=nullptr, *noise1=nullptr;
  file0.GetObject("Pedestals", pedestal0);
  file0.GetObject("Noise", noise0);
  file1.GetObject("Pedestals", pedestal1);
  file1.GetObject("Noise", noise1);

  //const auto& rocPedestal0 = pedestal0->getCalArray(0);
  //const auto& rocNoise0 = noise0->getCalArray(0);

  auto hPedestal2D = new TH2F("hPedestals2D",";pad row; pad",63,0,62,200,-100,100);

  ROC roc(0);
  //CalPad *mPedestal = new CalPad(roc);
  //CalDet<float> *mPedestal = new CalDet<float>(PadSubset::ROC);
  //CalDet<float> *mNoise = new CalDet<float>(PadSubset::ROC);


  //CalArray<float> calROCPedestal = mPedestal->getCalArray(roc);
  //CalArray<float> calROCNoise = mNoise->getCalArray(roc);

  CalDet<float> *DCM = new CalDet<float>(PadSubset::ROC);
  CalArray<float> DeadChannelMap = DCM->getCalArray(roc);

  for (int irow = 0; irow<63; ++irow) {
    int npads = mapper.getNumberOfPadsInRowROC(roc, irow);
    //cout<<endl<<"npads: "<<npads<<endl;
    for (int ipad = 0; ipad<npads; ++ipad) {
      const auto ped0 = pedestal0->getValue(roc, irow, ipad);
      const auto no0 = noise0->getValue(roc, irow, ipad);
      const auto ped1 = pedestal1->getValue(roc, irow, ipad);
      const auto no1 = noise1->getValue(roc, irow, ipad);
      int pad = ipad + npads/2;
      //calROCPedestal.setValue(irow, pad, ped0);
      //calROCNoise.setValue(irow, pad, no0);

      if (((irow >= 0 && irow <=16 && ipad >= 7 && ipad <= 12 && ped0 >= 200) || (irow >= 0 && irow <=16 && ipad >= 7 && ipad <= 12 && ped1 >= 200)) || ((irow >= 0 && irow <=16 && ipad >= 7 && ipad <= 12 && no0 <= .5) && (irow >= 0 && irow <=16 && ipad >= 7 && ipad <= 12 && no1 <= .5))) {
        cout<<endl<<"Warning: Bad channel in overlap region!\t row: "<<irow<<"\t"<<"pad: "<<ipad<<endl;
      }

      if (!(ped0 > 200 && ped1 > 200)) {
        DeadChannelMap.setValue(irow, ipad, 1);
      }
      if (no0 <= .5 && no1 <= .5) {
        DeadChannelMap.setValue(irow, ipad, 0);
      }
      else if (ped0 > 10 && ped0 <= 200 && ped1 > 10 && ped1 <= 200) {
        DeadChannelMap.setValue(irow, ipad, 2);
      }
      else {
        DeadChannelMap.setValue(irow,ipad, 1);
      }
      hPedestal2D->Fill(irow,ipad,ped0);
    }
  }

  //auto pedtest = Painter::getHistogram2D(calROCPedestal);
  //auto noisetest = Painter::getHistogram2D(calROCNoise);
  auto dcm = Painter::getHistogram2D(DeadChannelMap);


  auto cPedestal2D=new TCanvas("cPedestal2D","Pedestal2D");
  hPedestal2D->Draw("colz");

  //auto cPedtest=new TCanvas("cPedtest","TestPeds");
  //pedtest->Draw("colz");

  //auto cNoisetest=new TCanvas("cNoisetest","NoiseTest");
  //noisetest->Draw("colz");

  auto cDCM=new TCanvas("cDeadChannelMap","DeadChannelMap");
  dcm->Draw("colz");
}
