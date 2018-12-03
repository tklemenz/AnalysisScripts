#include "TPCSimulation/Cluster.h"
#include "TPCBase/Mapper.h"
#include "TPCReconstruction/TrackTPC.h"
#include "TPCBase/CRU.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"

#pragma link C++ class std::vector<o2::TPC::TrackTPC>+;

using namespace o2::TPC;


void GetBinMinMax(const TH1 *hist, const float frac, int &bin1, int &bin2)
{

  const int binMax=hist->GetMaximumBin();
  const float contMax=hist->GetBinContent(binMax);
  bin1=binMax;
  bin2=binMax;
  while ( (bin1--)>binMax/3. ) if (hist->GetBinContent(bin1)<frac*contMax) break;
  while ( (bin2++)<binMax*3. ) if (hist->GetBinContent(bin2)<frac*contMax) break;
}

float GetCPad(Cluster cl)
{
  Mapper &mapper = Mapper::instance();
  DigitPos pos(cl.getCRU(), PadPos(cl.getRow(), cl.getPadMean()));
  float row = pos.getPadSecPos().getPadPos().getRow();
  float pad = pos.getPadSecPos().getPadPos().getPad();
  return pad - mapper.getNumberOfPadsInRowSector(row)/2;
}



int getdEdxResolutionO2_Sectors(TString argv){

  //const char *sarg = argv;
  //const char *exRows = argv[2];

  TString filename(argv);
  ///int dEdxExcludeRows(exRows);
  Mapper &mapper = Mapper::instance();

  TFile *TreeFile = TFile::Open(filename.Data());
  cout<<endl<<endl<<filename<<endl<<endl;
  TTree *tree = (TTree*)gDirectory->Get("events");

  ///TH2D *fGainMap = 0x0;

  ///TFile *fGain = new TFile("~/GEMtest/code/myana/GainMap_FilesForGainMap.txt_2017-04-24_1759_0_07_trunc/GainMap_Trunc.root");

  ///fGainMap = (TH2D*)fGain->Get("GainMap");

  cout<<"DEBUG POINT 1"<<endl;




  /* ============================================================
   * ====================== Cut parameter =======================
   * ============================================================ */


  int TrPerEv = 1;
  int nclCut = 32;
  int CherCutLow = 40;
  int CherCutHigh = 50;
  int timeMeanLow = 25;
  int timeMeanHigh = 90;
  int cpadLow = 3;
  int cpadHigh = 17;
  float nclFracoutofCPad = .3;
  //int dEdxExcludeRows = 0;
  int excludeEdge = 1;//(dEdxExcludeRows/10)%10;
  //dEdxExcludeRows%=10;



  /* ============================================================
   * ============= Define histograms and graphs =================
   * ============================================================*/


  //TH1D *hNclusters   = new TH1D("hNclusters", ";Number of clusters; Counts", 65,0,65);
  TH1D *hNclustersUsed     = new TH1D("hNclustersUsed", ";Number of clusters per track; Counts", 300,0,300);


  TH1D *hNIROCtracks       = new TH1D("hNtracksIROC","; Number of tracks per event; # counts",20,0,20);
  TH1D *hNIROCtracksUsed   = new TH1D("hNtracksIROCUsed","; Number of tracks per event; # counts",50,0,50);

  TH1F *hdEdxEleTot        = new TH1F("hdEdxEleTot", "; d#it{E}/d#it{x} Q_{tot} [a.u.]; # counts", 250, 0, 500);
  TH1F *hdEdxPionTot       = new TH1F("hdEdxTot", "; d#it{E}/d#it{x} Q_{tot} [a.u.]; # counts", 250, 0, 500);

//  TH1F *hdEdxEleMax        = new TH1F("hdEdxEleMax", "; d#it{E}/d#it{x} Q_{tot} [a.u.]; # counts", 100, 0, 120);
//  TH1F *hdEdxPionMax       = new TH1F("hdEdxMax", "; d#it{E}/d#it{x} Q_{tot} [a.u.]; # counts", 100, 0, 120);

  TH1D *hQ              = new TH1D("hQ", "; Total cluster charge Q [ADC counts]; # counts", 600,0,600);
  TH1D *hQUsed          = new TH1D("hQUsed", "; Total cluster charge Q_{tot} [ADC counts]; # counts", 600,0,600);

  ///TH1D *hCherenkov         = new TH1D("hCherenkov", "; ADC signal; # counts", 90,20,196);
  ///TH1D *hCherenkovUsedPions     = new TH1D("hCherenkovUsed", "; ADC signal; # counts", 90,20,196);
  ///TH1D *hCherenkovUsedEle     = new TH1D("hCherenkovUsedEle", "; ADC signal; # counts", 90,20,196);

  TH2D *hPadOccupancy      = new TH2D("hOcc", "; Row; Pad", 250,-10,240,300,-100,200);
  TH2D *hPadOccupancyUsed  = new TH2D("hOccUsed", "; Row; Pad", 63,0,63,35,-5,30);

  TH2D *hQtimeMean       = new TH2D("hQTime","; TimeBin; Q_{tot}",1500,0,1500,1000,0,3000);
  TH2D *hQtimeMeanUsed   = new TH2D("hQTimeUsed","; TimeBin; Q_{tot}",105,0,105,200,0,600);

  const int nsectors = 36;
  TH2D *hOcc[nsectors];
  if (showSec){
    for (int isec = 0; isec<nsectors; ++isec){
        hOcc[isec] = new TH2D(Form("hOcc_Sec_%i", isec), "; row; pad", 170, -10, 160,300,-100,200);
    }
  }



/* =============================================================================
 * ================== Loop over events and apply all cuts ======================
 * ======================= meanwhile fill histograms ===========================
 * ============================================================================= */


  std::vector<TrackTPC> *fEvent = 0;
  //std::vector<TrackTPC> *pEvent = &Event;
  //std::vector<TrackTPC> *pEvent = 0;
  cout<<"DEBUG POINT 2"<<endl;

  tree->SetBranchAddress("Tracks", &fEvent);

  //GEMTrack::SetRemoveRowsdEdx("9,17,25,32,39,45,51,57,58,59,60,61,62");

  int usedTracks = 0;
  int Tracks = 0;
  int OneTrackEvents = 0;
  float pions = 0;
  float electrons = 0;
  int usedcl = 0;
  int onetrcl = 0;
  float pioncl = 0;
  float elecl = 0;

  cout<<"DEBUG POINT 3"<<endl;


  for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
    tree->GetEntry(iev);

    int NTracks = fEvent->size();
    ///int CherenkovValue = fEvent->GetCherenkovValue();

    Tracks += NTracks;

    hNIROCtracks->Fill(NTracks);



/*========================================== CUT ================================================*/
    ///if (NTracks != TrPerEv) continue;

//    ++OneTrackEvents;
//    cout<<"DEBUG POINT 4"<<endl;

  for (auto& trackObject : *fEvent) {
    std::vector<Cluster> clCont;
    trackObject.getClusterVector(clCont);


/*========================================== CUT ================================================*/
      ///if (track->GetROC() != 0) continue;															// only ROC 0


      int ncl = clCont.size();
      bool isok = true;
      int ncledge = 0;

      onetrcl += ncl;
///////////////////////////////////////////////////////////////////////
      for (auto& clusterObject : clCont) {
        DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
        float row = pos.getPadSecPos().getPadPos().getRow();
        float pad = pos.getPadSecPos().getPadPos().getPad();
        const int cpad = GetCPad(clusterObject);
        UShort_t timeMean = clusterObject.getTimeMean();                          /// find solution for timeMean problem
        double Q = clusterObject.getQ();
        hQtimeMean->Fill(timeMean,Q);
        hQ->Fill(Q);
        hPadOccupancy->Fill(row,cpad);
      }


      cout<<"DEBUG POINT 5"<<endl;


/*========================================== CUT ================================================*/
     /// if (ncl < nclCut) continue;																	// cut on number of clusters per track



      ///hCherenkov->Fill(CherenkovValue);
/*========================================== CUT ================================================*/
///	  if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;					// PID via Cherenkov
///	    if (CherenkovValue < CherCutLow) {
///		  hCherenkovUsedPions->Fill(CherenkovValue);
///		}
///		else{
///		  hCherenkovUsedEle->Fill(CherenkovValue);
///		}



    for (auto& clusterObject : clCont) {															// make cuts on cluster properties
      DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
      float row = pos.getPadSecPos().getPadPos().getRow();
      float pad = pos.getPadSecPos().getPadPos().getPad();
      const int cpad = GetCPad(clusterObject);
      ushort timeMean = clusterObject.getTimeMean();
      double gainMapCorr = 1.;

      //cout<<endl<<"Row: "<<double(row)<<"\t"<<"Pad: "<<double(pad)<<"\t"<<"CPad: "<<double(cpad)<<endl;


/*========================================== CUT ================================================*/
    ///if (timeMean < timeMeanLow || timeMean > timeMeanHigh){isok = false; break;}					// cut on time max, only clusters from a certain range on z-axis can come from a particle


    ///if (cpad < cpadLow || cpad > cpadHigh){++ncledge;}
      //cout<<"DEBUG POINT 3"<<endl;
    ///if (fGainMap){
       /// gainMapCorr = fGainMap->GetBinContent(row+6,cpad+60);
        //cout<<endl<<endl<<fGainMap2->GetBinContent(row+6,cpad+60)<<"\t\t"<<row<<"\t"<<pad<<endl;
        ///if (gainMapCorr<1e-4) gainMapCorr=1.;
      ///}
    ///cl->CorrectCharge(gainMapCorr);                                                              //has to be put back in at some point
    //cout<<"DEBUG POINT 5"<<endl;															// gain correction from gainmap
    }

/*========================================== CUT ================================================*/
    ///  if (((float(ncledge)/float(ncl)) > nclFracoutofCPad) && excludeEdge) continue;					// cut on detector geometry


    if (isok == true){																			// track accepted
    usedcl += ncl;
    float dEdxTot = trackObject.getTruncatedMean(.05,.7,1);
    float dEdxMax = trackObject.getTruncatedMean(.05,.7,0);

//		if (CherenkovValue < CherCutLow){
//          hdEdxPionTot->Fill(dEdxTot);
//          hdEdxPionMax->Fill(dEdxMax);
//          ++pions;
//          pioncl += ncl;
//		}
//		if (CherenkovValue > CherCutHigh){
//          hdEdxEleTot->Fill(dEdxTot);
//          hdEdxEleMax->Fill(dEdxMax);
//          ++electrons;
//          elecl += ncl;
//		}
    hdEdxPionTot->Fill(dEdxTot);
    ++usedTracks;
    hNIROCtracksUsed->Fill(NTracks);
    hNclustersUsed->Fill(ncl);
    for (auto& clusterObject : clCont) {															// loop over clusters
      DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
      float row = pos.getPadSecPos().getPadPos().getRow();
      float pad = pos.getPadSecPos().getPadPos().getPad();
      const int cpad = GetCPad(clusterObject);
      ushort timeMean = clusterObject.getTimeMean();
      double Q = clusterObject.getQ();

      hQtimeMeanUsed->Fill(timeMean,Q);
      QUsed->Fill(Q);
      hPadOccupancyUsed->Fill(row,cpad);

      }
    }
    }
  }



  /* ======================================================================
   * ===================== Plot everything ================================
   * ====================================================================== */

  hQ->SetLineColor(kBlue+2);

  TCanvas *c1 = new TCanvas();
  hNIROCtracksUsed->Draw();

  TCanvas *c2 = new TCanvas();
  hNIROCtracks->Draw();

  TCanvas *c3 = new TCanvas();
  hNclustersUsed->Draw();

  TCanvas *c4 = new TCanvas();
  //hQtimeMean->GetZaxis()->SetRangeUser(0,820);
  hQtimeMean->Draw("colz");

  TCanvas *c5 = new TCanvas();
  //hQtimeMeanUsed->GetZaxis()->SetRangeUser(0,820);
  hQtimeMeanUsed->Draw("colz");

  TCanvas *c6 = new TCanvas();
  hQ->Draw();

  TCanvas *c7 = new TCanvas();
  hQUsed->Draw();

  TCanvas *c8 = new TCanvas();
  //hPadOccupancy->GetZaxis()->SetRangeUser(0,3700);
  hPadOccupancy->Draw("colz");

  TCanvas *c9 = new TCanvas();
  //hPadOccupancyUsed->GetZaxis()->SetRangeUser(0,3700);
  hPadOccupancyUsed->Draw("colz");

  ///TCanvas *c10 = new TCanvas();
  ///c10->SetLogy();
  ///hCherenkov->SetLineColor(kBlue);
  ///hCherenkov->SetLineWidth(2);
  ///hCherenkov->Draw();

  ///TCanvas *c11 = new TCanvas();
  ///c11->SetLogy();
  ///hCherenkovUsedPions->GetYaxis()->SetRangeUser(0.5,21000);
  ///hCherenkovUsedPions->SetLineColor(kBlue+2);
  ///hCherenkovUsedEle->SetLineColor(kGreen+2);
  ///hCherenkovUsedPions->SetLineWidth(2);
  ///hCherenkovUsedEle->SetLineWidth(2);

  ///hCherenkovUsedPions->Draw();
  ///hCherenkovUsedEle->Draw("same");


  TCanvas *c15 = new TCanvas();
  TF1 *pionfit = new TF1("pionfit","gaus",hdEdxPionTot->GetXaxis()->GetXmin(),hdEdxPionTot->GetXaxis()->GetXmax());
  TF1 *electronfit = new TF1("electronfit","gaus",hdEdxEleTot->GetXaxis()->GetXmin(),hdEdxEleTot->GetXaxis()->GetXmax());
  hdEdxEleTot->SetLineColor(kGreen+2);
  hdEdxPionTot->SetLineColor(kBlue+2);

  const Float_t frac=0.2;
  Int_t bin1=0,bin2=0;

  GetBinMinMax(hdEdxPionTot,frac,bin1,bin2);
  hdEdxPionTot->Fit("pionfit","","",hdEdxPionTot->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionTot->GetXaxis()->GetBinUpEdge(bin2));
//  GetBinMinMax(hdEdxEleTot,frac,bin1,bin2);
//  hdEdxEleTot->Fit("electronfit","","",hdEdxEleTot->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleTot->GetXaxis()->GetBinUpEdge(bin2));

  //alternative fit
  //hdEdxPionTot->Fit("pionfit");
  //hdEdxEleTot->Fit("electronfit");


  //hdEdxEleTot->GetFunction("electronfit")->SetLineColor(kGreen+2);
  hdEdxPionTot->GetFunction("pionfit")->SetLineColor(kBlue+2);
  hdEdxPionTot->Draw();
  //hdEdxEleTot->Draw("same");
  float pionmeanTot = pionfit->GetParameter(1);
  float pionsigmaTot = pionfit->GetParameter(2);
//  float electronmeanTot = electronfit->GetParameter(1);
//  float electronsigmaTot = electronfit->GetParameter(2);

  float pionres = pionsigmaTot/pionmeanTot;
  //float electronres = electronsigmaTot/electronmeanTot;
  float pionsigmaerr = pionfit->GetParError(2);
  float pionmeanerr = pionfit->GetParError(1);
//  float electronsigmaerr = electronfit->GetParError(2);
//  float electronmeanerr = electronfit->GetParError(1);


//  float separationpower = 2*(electronmeanTot-pionmeanTot)/(pionsigmaTot+electronsigmaTot);


  float pionreserror = sqrt(pow((pionsigmaerr/pionmeanTot),2)+pow(((pionsigmaTot*pionmeanerr)/pow(pionmeanTot,2)),2));
//  float electronreserror = sqrt(pow((electronsigmaerr/electronmeanTot),2)+pow(((electronsigmaTot*electronmeanerr)/pow(electronmeanTot,2)),2));

//  float separationpowererr = sqrt(pow((2*electronmeanerr/(electronsigmaTot+pionsigmaTot)),2)+pow((2*pionmeanerr/(electronsigmaTot+pionsigmaTot)),2)+
//        pow(((2*electronsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2)+pow(((2*pionsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2));

//  float elechisquareTot = electronfit->GetChisquare();
//  float eleNDFTot = electronfit->GetNDF();
  float pionchisquareTot = pionfit->GetChisquare();
  float pionNDFTot = pionfit->GetNDF();



//  TCanvas *c16 = new TCanvas();
//  TF1 *pionfitMax = new TF1("pionfitMax","gaus",hdEdxPionMax->GetXaxis()->GetXmin(),hdEdxPionMax->GetXaxis()->GetXmax());
//  TF1 *electronfitMax = new TF1("electronfitMax","gaus",hdEdxEleMax->GetXaxis()->GetXmin(),hdEdxEleMax->GetXaxis()->GetXmax());
//  hdEdxEleMax->SetLineColor(kGreen+2);
//  hdEdxPionMax->SetLineColor(kBlue+2);

//  GetBinMinMax(hdEdxPionMax,frac,bin1,bin2);
//  hdEdxPionMax->Fit("pionfitMax","","",hdEdxPionMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionMax->GetXaxis()->GetBinUpEdge(bin2));
//  GetBinMinMax(hdEdxEleMax,frac,bin1,bin2);
//  hdEdxEleMax->Fit("electronfitMax","","",hdEdxEleMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleMax->GetXaxis()->GetBinUpEdge(bin2));



//  hdEdxEleMax->GetFunction("electronfitMax")->SetLineColor(kGreen+2);
//  hdEdxPionMax->GetFunction("pionfitMax")->SetLineColor(kBlue+2);
//  hdEdxPionMax->Draw();
//  hdEdxEleMax->Draw("same");
//  float pionmeanMax = pionfitMax->GetParameter(1);
//  float pionsigmaMax = pionfitMax->GetParameter(2);
//  float electronmeanMax = electronfitMax->GetParameter(1);
//  float electronsigmaMax = electronfitMax->GetParameter(2);

//  float pionresMax = pionsigmaMax/pionmeanMax;
//  float electronresMax = electronsigmaMax/electronmeanMax;
//  float pionsigmaerrMax = pionfitMax->GetParError(2);
//  float pionmeanerrMax = pionfitMax->GetParError(1);
//  float electronsigmaerrMax = electronfitMax->GetParError(2);
//  float electronmeanerrMax = electronfitMax->GetParError(1);


//  float separationpowerMax = 2*(electronmeanMax-pionmeanMax)/(pionsigmaMax+electronsigmaMax);


//  float pionreserrorMax = sqrt(pow((pionsigmaerrMax/pionmeanMax),2)+pow(((pionsigmaMax*pionmeanerrMax)/pow(pionmeanMax,2)),2));
//  float electronreserrorMax = sqrt(pow((electronsigmaerrMax/electronmeanMax),2)+pow(((electronsigmaMax*electronmeanerrMax)/(electronmeanMax,2)),2));

//  float separationpowererrMax = sqrt(pow((2*electronmeanerrMax/(electronsigmaMax+pionsigmaMax)),2)+pow((2*pionmeanerrMax/(electronsigmaMax+pionsigmaMax)),2)+
//        pow(((2*electronsigmaerrMax*(pionmeanMax-electronmeanMax))/(pow((electronsigmaMax+pionsigmaMax),2))),2)+pow(((2*pionsigmaerrMax*(pionmeanMax-electronmeanMax))/(pow((electronsigmaMax+pionsigmaMax),2))),2));

//  float elechisquareMax = electronfitMax->GetChisquare();
//  float eleNDFMax = electronfitMax->GetNDF();
//  float pionchisquareMax = pionfitMax->GetChisquare();
//  float pionNDFMax = pionfitMax->GetNDF();


  cout<<endl;cout<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<endl;
  cout<<"============================= Total Charge ======================================="<<endl;
  cout<<"Pion dE/dx resolution:          ("<<pionres*100<<" +- "<<pionreserror*100<<") %"<<endl;
  //cout<<"Electron dE/dx resolution:      ("<<electronres*100<<" +- "<<electronreserror*100<<") %"<<endl;
  //cout<<"Separation power:               ("<<separationpower<<" +- "<<separationpowererr<<") sigma"<<endl;
  cout<<"Pion Chisquare:                  "<<pionchisquareTot/pionNDFTot<<endl;
  //cout<<"Electron Chisquare:              "<<elechisquareTot/eleNDFTot<<endl;
  cout<<"Pion mean:                       "<<pionmeanTot<<" +- "<<pionsigmaTot<<endl;
  //cout<<"Electron mean:                   "<<electronmeanTot<<" +- "<<electronsigmaTot<<endl;
  cout<<endl;
  cout<<"============================== Max Charge ========================================"<<endl;
//  cout<<"Pion dE/dx resolution:          ("<<pionresMax*100<<" +- "<<pionreserrorMax*100<<") %"<<endl;
//  cout<<"Electron dE/dx resolution:      ("<<electronresMax*100<<" +- "<<electronreserrorMax*100<<") %"<<endl;
//  cout<<"Separation power:               ("<<separationpowerMax<<" +- "<<separationpowererrMax<<") sigma"<<endl;
//  cout<<"Pion Chisquare:                  "<<pionchisquareMax/pionNDFMax<<endl;
//  cout<<"Electron Chisquare:              "<<elechisquareMax/eleNDFMax<<endl;
//  cout<<"Pion mean:                       "<<pionmeanMax<<" +- "<<pionsigmaMax<<endl;
//  cout<<"Electron mean:                   "<<electronmeanMax<<" +- "<<electronsigmaMax<<endl;
//  cout<<endl;
  cout<<"============================== Statistics ========================================"<<endl;
  cout<<"Number of tracks:                "<<Tracks<<endl;
  cout<<"Number of one-track-events:      "<<OneTrackEvents<<endl;
  cout<<"Number of used tracks:           "<<usedTracks<<endl;
  cout<<"Number of pions:                 "<<pions<<endl;
  //cout<<"Number of electrons:             "<<electrons<<endl;
  cout<<endl;
  /*cout<<"Number of one-track-clusters:    "<<onetrcl<<endl;
  cout<<"Number of used clusters:         "<<usedcl<<endl;
  cout<<"Number of electron clusters:     "<<elecl<<endl;
  cout<<"Number of pion clusters:         "<<pioncl<<endl;*/
  cout<<"Average clusters per pion:       "<<pioncl/pions<<endl;
  //cout<<"Average clusters per electron:   "<<elecl/electrons<<endl;
  cout<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<endl;cout<<endl;

  return 0;
}


