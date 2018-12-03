/// pedestalFile0 from first FEC setting
/// pedestalFile1 from second FEC setting

int getDeadChannelMap(TString pedestalFile0, TString pedestalFile1, TString filename="DeadChannelMap.root")
{
  using namespace o2::TPC;
  TFile f0(pedestalFile0);
  TFile f1(pedestalFile1);
  gROOT->cd();

  static const Mapper& mapper = Mapper::instance();
  ROC roc(0);

  CalDet<float> *pedestal0 = nullptr, *pedestal1 = nullptr;
  f0.GetObject("Pedestals", pedestal0);
  f1.GetObject("Pedestals", pedestal1);

  CalDet<bool> *DCM = new CalDet<bool>(PadSubset::ROC);
  DCM->setName("DeadChannelMap");
  CalArray<bool>& DeadChannelMap = DCM->getCalArray(roc);

  //CalDet<float> *dummy = new CalDet<float>(PadSubset::ROC);
  //CalArray<float>& dummyArr = dummy->getCalArray(roc);

  for (int irow=0; irow<63; ++irow) {
    int npads = mapper.getNumberOfPadsInRowROC(roc, irow);
    //cout<<irow<<"\t"<<npads<<endl;
    for (int ipad=0; ipad<npads; ++ipad) {
      //cout<<endl<<pad<<endl;
      const auto ped0 = pedestal0->getValue(roc, irow, ipad);
      const auto ped1 = pedestal1->getValue(roc, irow, ipad);
      if (ped0 != 0 || ped1 != 0) {
        //dummyArr.setValue(irow, ipad, 1);
        DeadChannelMap.setValue(irow, ipad, true);
        //std::cout << DeadChannelMap.getValue(irow, ipad) << "true \n";
        //cout<<irow<<"\t"<<ipad<<endl;
      }
      else {
        //dummyArr.setValue(irow, ipad, 0);
        DeadChannelMap.setValue(irow, ipad, false);
        //std::cout << DeadChannelMap.getValue(irow, ipad) << "false \n";
      }
    }
  }

  /*dummyArr.setValue(32, 51, 0);
  dummyArr.setValue(28, 60, 0);
  dummyArr.setValue(47, 57, 0);
  dummyArr.setValue(48, 71, 0);
  dummyArr.setValue(58, 63, 0);
  dummyArr.setValue(59, 61, 0);
  dummyArr.setValue(61, 69, 0);
  dummyArr.setValue(44, 57, 0);
  dummyArr.setValue(4, 42, 0);
  dummyArr.setValue(4, 33, 0);
  dummyArr.setValue(32, 33, 0);
  dummyArr.setValue(58, 35, 0);
  dummyArr.setValue(59, 34, 0);*/

  DeadChannelMap.setValue(32, 51, false);
  DeadChannelMap.setValue(28, 60, false);
  DeadChannelMap.setValue(47, 57, false);
  DeadChannelMap.setValue(48, 71, false);
  DeadChannelMap.setValue(58, 63, false);
  DeadChannelMap.setValue(59, 61, false);
  DeadChannelMap.setValue(61, 69, false);
  DeadChannelMap.setValue(44, 57, false);
  DeadChannelMap.setValue(4, 42, false);
  DeadChannelMap.setValue(4, 33, false);
  DeadChannelMap.setValue(32, 33, false);
  DeadChannelMap.setValue(58, 35, false);
  DeadChannelMap.setValue(59, 34, false);

  /*auto dummyplot = Painter::getHistogram2D(dummyArr);

  auto cDummyDCM = new TCanvas("cDummyDCM", "DummyDCM");
  dummyplot->Draw("colz");*/

  TFile *f = TFile::Open(filename, "recreate");
  f->WriteObject(DCM, "DeadChannelMap");
  f->Close();
  delete f;

  return 0;
}
