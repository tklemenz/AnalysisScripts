/// dEdx resolution, separation power between electrons and pios, some QA


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
#include "TProfile.h"
#include "TPaveText.h"

#pragma link C++ class std::vector<o2::TPC::TrackTPC>+;

using namespace o2::TPC;

struct EventHeader
{
  int run;
  float cherenkovValue;

};


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

int dEdxQA(TString trackfile, bool showSectors=0){

  bool showSec = showSectors;

  float electronmeanTot, electronsigmaTot, pionmeanTot, pionsigmaTot, electronmeanerr, electronsigmaerr, pionmeanerr, pionsigmaerr, electronmeanMax,
      electronsigmaMax, pionmeanMax, pionsigmaMax, electronmeanerrMax, electronsigmaerrMax, pionmeanerrMax, pionsigmaerrMax, pionres, pionresMax, electronres, electronresMax, separationpower, separationpowerMax;
  int Pions, Electrons, PionClusters, ElectronClusters, Tracks, TracksUsed, EvWithCorrNTracks, DesiredNTrackClusters, ClustersUsed;

  electronmeanTot = electronsigmaTot = pionmeanTot = pionsigmaTot = electronmeanerr = electronsigmaerr = pionmeanerr = pionsigmaerr = electronmeanMax =
      electronsigmaMax = pionmeanMax = pionsigmaMax = electronmeanerrMax = electronsigmaerrMax = pionmeanerrMax = pionsigmaerrMax = pionres =
      pionresMax = electronres = electronresMax = separationpower = separationpowerMax = 0.;
  Pions = Electrons = PionClusters = ElectronClusters = Tracks = TracksUsed = EvWithCorrNTracks = DesiredNTrackClusters = ClustersUsed = 0;


  Mapper &mapper = Mapper::instance();

  ///read inputfile (e.g. tracks.root)
  TString filename = trackfile;

  TFile *TreeFile = TFile::Open(trackfile.Data());
  TTree *tree = (TTree*)gDirectory->Get("events");


  TFile *OutFile = new TFile("anaOut.root", "UPDATE");
  TTree *TrackAna = new TTree("TrackAna","dEdxAna");

  TrackAna->Branch("PidEdxResQ", &pionres, "pionres/F");
  TrackAna->Branch("PidEdxResQmax", &pionresMax, "pionresMax/F");
  TrackAna->Branch("EledEdxResQ", &electronres, "electronres/F");
  TrackAna->Branch("EledEdxResQmax", &electronresMax, "electronresMax/F");
  TrackAna->Branch("SepPowerQ", &separationpower, "separationpower/F");
  TrackAna->Branch("SepPowerQmax", &separationpowerMax);
  TrackAna->Branch("posEle", &electronmeanTot, "electronmeanTot/F");
  TrackAna->Branch("resEle", &electronsigmaTot, "electronsigmaTot/F");
  TrackAna->Branch("posPi", &pionmeanTot, "pionmeanTot/F");
  TrackAna->Branch("resPi", &pionsigmaTot, "pionsigmaTot/F");
  TrackAna->Branch("posEleErr", &electronmeanerr, "electronmeanerr/F");
  TrackAna->Branch("resEleErr", &electronsigmaerr, "electronsigmaerr/F");
  TrackAna->Branch("posPiErr", &pionmeanerr, "pionmeanerr/F");
  TrackAna->Branch("resPiErr", &pionsigmaerr, "pionsigmaerr/F");
  TrackAna->Branch("posEleQmax", &electronmeanMax, "electronmeanMax/F");
  TrackAna->Branch("resEleQmax", &electronsigmaMax, "electronsigmaMax/F");
  TrackAna->Branch("posPiQmax", &pionmeanMax, "pionmeanMax/F");
  TrackAna->Branch("resPiQmax", &pionsigmaMax, "pionsigmaMax/F");
  TrackAna->Branch("posEleQmaxErr", &electronmeanerrMax, "electronmeanerrMax/F");
  TrackAna->Branch("resEleQmaxErr", &electronsigmaerrMax, "electronsigmaerrMax/F");
  TrackAna->Branch("posPiQmaxErr", &pionmeanerrMax, "pionmeanerrMax/F");
  TrackAna->Branch("resPiQmaxErr", &pionsigmaerrMax, "pionsigmaerrMax/F");
  TrackAna->Branch("nPions", &Pions, "Pions/I");
  TrackAna->Branch("nElectrons", &Electrons, "Electrons/I");
  TrackAna->Branch("nPionCl", &PionClusters, "PionClusters/I");
  TrackAna->Branch("nEleCl", &ElectronClusters, "ElectronClusters/I");



  /// set cut parameters
  /// \param TrPerEv Tracks per Event
  /// \param nclCut minimum number of clusters in track
  /// \param CherCutLow upper limit of Cherenkov value from track for pions
  /// \param CherCutHigh lower limit of Cherenkov value from track for electrons
  /// \param timeMeanLow lower boundary for TimeMean value of cluster
  /// \param timeMeanHigh upper boundary for TimeMean value of cluster
  /// \param cpadLow lower cpad boundary of clusters
  /// \param cpadHigh upper cpad boundary of clusters
  /// \param nclFractionoutofCPad if fraction x of the clusters of a track are outside the above pad range, the track is cut
  int TrPerEv = 21;
  int nclCut = 132;
  int CherCutLow = 40;
  int CherCutHigh = 50;
  int timeMeanLow = 150;
  int timeMeanHigh = 600;
  int cpadLow = -20;
  int cpadHigh = 40;
  float nclFracoutofCPad = .5;
  //int dEdxExcludeRows = 0;
  int excludeEdge = 1;//(dEdxExcludeRows/10)%10;
  //dEdxExcludeRows%=10;



  /// define histograms and graphs
  TH1D *hNEvents           = new TH1D("hNEvents", ";Number of events; # counts", 500,0,500);

  TH1D *hNClusters         = new TH1D("hNClusters", ";Number of clusters per track; # counts", 170,0,170);
  TH1D *hNClustersUsed     = new TH1D("hNClustersUsed", ";Number of clusters per track; # counts", 170,0,170);


  TH1D *hNTracks           = new TH1D("hNTracks","; Number of tracks per event; # counts",50,0,50);
  //TH1D *hNTracksUsed       = new TH1D("hNTracksUsed","; Number of tracks per event; # counts",50,0,50);

  TH1F *hdEdxEleTot        = new TH1F("hdEdxEleTot", "; d#it{E}/d#it{x} Q [a.u.]; # counts", 500, 0, 600);
  TH1F *hdEdxPionTot       = new TH1F("hdEdxTot", "; d#it{E}/d#it{x} Q [a.u.]; # counts", 500, 0, 600);

  TH1F *hdEdxEleMax        = new TH1F("hdEdxEleMax", "; d#it{E}/d#it{x} Q_{max} [a.u.]; # counts", 100, 0, 120);
  TH1F *hdEdxPionMax       = new TH1F("hdEdxMax", "; d#it{E}/d#it{x} Q_{max} [a.u.]; # counts", 100, 0, 120);

  TH1D *hQ                 = new TH1D("hQ", "; Total cluster charge Q [ADC counts]; # counts",600,0,1200);
  TH1D *hQUsed             = new TH1D("hQUsed", "; Total cluster charge Q [ADC counts]; # counts", 600,0,1200);

  ///TH1D *hCherenkov         = new TH1D("hCherenkov", "; ADC signal; # counts", 90,20,196);
  ///TH1D *hCherenkovUsedPions     = new TH1D("hCherenkovUsed", "; ADC signal; # counts", 90,20,196);
  ///TH1D *hCherenkovUsedEle     = new TH1D("hCherenkovUsedEle", "; ADC signal; # counts", 90,20,196);

  TH2D *hPadOccupancy      = new TH2D("hOcc", "; Row; Pad", 250,-10,240,300,-100,200);
  TH2D *hPadOccupancyUsed  = new TH2D("hOccUsed", "; Row; Pad", 63,0,63,35,-5,30);

  TH2D *hQTimeMean         = new TH2D("hQTime","; TimeBin; Q",600,0,600,400,0,2000);
  TH2D *hQTimeMeanUsed     = new TH2D("hQTimeUsed","; TimeBin; Q",600,0,600,400,0,2000);
  TH2D *hQmaxTimeMean      = new TH2D("hQMaxTime","; TimeBin; Q_{max}",600,0,600,400,0,2000);

  TH2D *hQmaxRelPadPos     = new TH2D("hQmaxRelPadPos", "; relative position of the cluster; Q_{max}", 100,-.2,1.2,600,0,1200);

  const int nsectors = 36;
  TH2D *hOcc[nsectors];
  if (showSec){
    for (int isec = 0; isec<nsectors; ++isec){
        hOcc[isec] = new TH2D(Form("hOcc_Sec_%i", isec), "; row; pad", 1700, -10, 160,300,-100,200);
    }
  }

  //TH2D *PadSec17 = new TH2D("Sector 17 pads",";Row;Pad",1700, -10, 160, 300, -100, 200);

  //Histograms for cut testing
  TH1D *hdEdxTotNOCut = new TH1D("dEdxNOCut", "; d#it{E}/d#it{x} Q [a.u.]; # counts", 500,0,600);
  TH1D *hdEdxTotAfterNTrackCut = new TH1D("dEdxAfterNTrackCut", "; d#it{E}/d#it{x} Q [a.u.]; # counts", 500,0,600);
  TH1D *hdEdxTotAfterROCCut = new TH1D("dEdxAfterROCCut", "; d#it{E}/d#it{x} Q [a.u.]; # counts", 500,0,600);
  TH1D *hdEdxTotAfterNClCut = new TH1D("dEdxAfterNClCut", "; d#it{E}/d#it{x} Q [a.u.]; # counts", 500,0,600);
  //TH1D *hdEdxTotAfterCherenkovCut = new TH1D("dEdxAfterCherenkovCut", "; d#it{E}/d#it{x} Q [a.u.]; # counts", 500,0,600);
  TH1D *hdEdxTotAfterTimeCut = new TH1D("dEdxAfterTimeCut", "; d#it{E}/d#it{x} Q [a.u.]; # counts", 500,0,600);
  TH1D *hdEdxTotAfterNClEdgeCut = new TH1D("dEdxAfterNClEdgeCut", "; d#it{E}/d#it{x} Q [a.u.]; # counts", 500,0,600);

  TH1D *hdEdxMaxNOCut = new TH1D("dEdxMaxNoCut", "; d#it{E}/d#it{x} Q_{max} [a.u.]; # counts", 500,0,600);




  /// loop over events, apply cuts, fill histograms
  std::vector<TrackTPC> *vecEvent = nullptr;
  EventHeader Header;
  tree->SetBranchAddress("Tracks", &vecEvent);
  tree->SetBranchAddress("header", &Header);


  hNEvents->Fill(tree->GetEntriesFast());

  for (int iEv=0; iEv<tree->GetEntriesFast(); ++iEv){
    tree->GetEntry(iEv);

    int runNr = Header.run;
    float CherenkovValue = Header.cherenkovValue;


    int nTracks = vecEvent->size();
    hNTracks->Fill(nTracks);
    Tracks += nTracks;
    for (auto& trackObject : *vecEvent){
      std::vector<Cluster> clCont;
      trackObject.getClusterVector(clCont);
      int ncl = clCont.size();
      hNClusters->Fill(ncl);
      float dEdxTot = trackObject.getTruncatedMean(0.05,0.7,1);
      float dEdxMax = trackObject.getTruncatedMean(0.05,0.7,0);

      hdEdxTotNOCut->Fill(dEdxTot);
      hdEdxMaxNOCut->Fill(dEdxMax);

      for (auto& clusterObject : clCont){
        //const int row = clusterObject.getRow();
        const int cpad = GetCPad(clusterObject);
        ushort timeMean = clusterObject.getTimeMean();
        float Q = clusterObject.getQ();
        float QMax = clusterObject.getQmax();
        const CRU cru(clusterObject.getCRU());
        const PadRegionInfo& region = mapper.getPadRegionInfo(cru.region());
        const int row = clusterObject.getRow() + region.getGlobalRowOffset();



        hQmaxTimeMean->Fill(timeMean,QMax);
        hQTimeMean->Fill(timeMean,Q);
        hQ->Fill(Q);
        hPadOccupancy->Fill(row,cpad);
      }
    }

    //if (nTracks != TrPerEv) continue;
    ++EvWithCorrNTracks;
//________________________________________________________________________________
    for (auto& trackObject : *vecEvent){
      float dEdxTot = trackObject.getTruncatedMean(0.05,0.7,1);
      float dEdxMax = trackObject.getTruncatedMean(0.05,0.7,0);

      hdEdxTotAfterNTrackCut->Fill(dEdxTot);
      std::vector<Cluster> clCont;
      trackObject.getClusterVector(clCont);

      //if (trackObject->GetROC() != 0) continue;

      hdEdxTotAfterROCCut->Fill(dEdxTot);
      int ncl = clCont.size();
      bool isok = true;
      float ncledge = 0;
      DesiredNTrackClusters += ncl;



      ///if (ncl < nclCut) continue;
      //dEdxTot = trackObject.getTruncatedMean(0.05,0.7,1);

      hdEdxTotAfterNClCut->Fill(dEdxTot);
      //hCherenkov->Fill(CherenkovValue);
      //if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;
      //hdEdxTotAfterCherenkovCut->Fill(dEdx);
      for (auto& clusterObject : clCont) {															// make cuts on cluster properties
        DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
        float row = pos.getPadSecPos().getPadPos().getRow();
        float pad = pos.getPadSecPos().getPadPos().getPad();
        const int cpad = GetCPad(clusterObject);
        ushort timeMean = clusterObject.getTimeMean();
        double gainMapCorr = 1.;
        float Qmax = clusterObject.getQmax();

        const float padmean = clusterObject.getPadMean();
        float relPadPos = padmean-pad;

        hQmaxRelPadPos->Fill(relPadPos,Qmax);

        ///if (timeMean < timeMeanLow || timeMean > timeMeanHigh){isok =false; break;}
        ///if (cpad < cpadLow || cpad > cpadHigh){++ncledge;}
        ///\todo gain correction with gain map
      }
      if (isok==true){
        hdEdxTotAfterTimeCut->Fill(dEdxTot);
      }
        ///if (((ncledge/float(ncl)) > nclFracoutofCPad) && excludeEdge) continue;
      if (isok==true){
        hdEdxTotAfterNClEdgeCut->Fill(dEdxTot);
      }
        if (isok == true){
          ClustersUsed += ncl;
          if (CherenkovValue < CherCutLow){
            hdEdxPionTot->Fill(dEdxTot);
            hdEdxPionMax->Fill(dEdxMax);
            ++Pions;
            PionClusters += ncl;
          }
          if (CherenkovValue > CherCutHigh){
            hdEdxEleTot->Fill(dEdxTot);
            hdEdxEleMax->Fill(dEdxMax);
            ++Electrons;
            ElectronClusters += ncl;
          }
//          hdEdxPionTot->Fill(dEdxTot);    ///remove when cherenkov value is implemented
//          hdEdxPionMax->Fill(dEdxMax);
          ++TracksUsed;
          //hNTracksUsed->Fill(nTracks);  /// this method works only for events with one track
          hNClustersUsed->Fill(ncl);
          for (auto& clusterObject : clCont){
            //const int row = clusterObject.getRow();
            const int cpad = GetCPad(clusterObject);
            //const int pad = TMath::Nint(clusterObject.getPadMean());
            const ushort timeMean = clusterObject.getTimeMean();
            double Q = clusterObject.getQ();
            const CRU cru = clusterObject.getCRU();
            const Sector sector = cru.sector();
            int isector = sector.getSector();
//            const PadRegionInfo& region = mapper.getPadRegionInfo(cru.region());
//            const int row       = clusterObject.getRow() + region.getGlobalRowOffset();
            DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
            float row = pos.getPadSecPos().getPadPos().getRow();
            float pad = pos.getPadSecPos().getPadPos().getPad();


            hQTimeMeanUsed->Fill(timeMean,Q);
            hQUsed->Fill(Q);
            hPadOccupancyUsed->Fill(row,cpad);
            if (showSec){
              hOcc[isector]->Fill(row,cpad);
            }
//            if (isector == 17){
//              PadSec17->Fill(row,pad);
//            }
          }
        }
      }
   // __________________________________________________________________
    }
  ///create plots
  hQ->SetLineColor(kBlue+2);

  TCanvas *NEvents = new TCanvas();
  hNEvents->Draw();

//  TCanvas *TrUsed = new TCanvas();
//  hNTracksUsed->Draw();

  TCanvas *Tr = new TCanvas();
  hNTracks->Draw();

  TCanvas *ClusUsed = new TCanvas();
  hNClustersUsed->Draw();

  TCanvas *Clus = new TCanvas();
  hNClusters->Draw();

  TCanvas *TimeMean = new TCanvas();
  ///hQTimeMean->GetZaxis()->SetRangeUser(0,820);
  hQTimeMean->Draw("colz");

  TCanvas *TimeMeanUsed = new TCanvas();
  ///hQTimeMeanUsed->GetZaxis()->SetRangeUser(0,820);
  hQTimeMeanUsed->Draw("colz");

  TCanvas *ClusterCharge = new TCanvas();
  hQ->Draw();

  TCanvas *ClusterChargeUsed = new TCanvas();
  hQUsed->Draw();

  TCanvas *PadOcc = new TCanvas();
  ///hPadOccupancy->GetZaxis()->SetRangeUser(0,3700);
  hPadOccupancy->Draw("colz");

  TCanvas *PadOccUsed = new TCanvas();
  ///hPadOccupancyUsed->GetZaxis()->SetRangeUser(0,3700);
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

/* ===========================================================================================================
 * Fit the dEdx plot for the total cluster charge
*/

  TCanvas *c15 = new TCanvas();
  TF1 *pionfit = new TF1("pionfit","gaus",hdEdxPionTot->GetXaxis()->GetXmin(),hdEdxPionTot->GetXaxis()->GetXmax());
  ///TF1 *electronfit = new TF1("electronfit","gaus",hdEdxEleTot->GetXaxis()->GetXmin(),hdEdxEleTot->GetXaxis()->GetXmax());
  //hdEdxEleTot->SetLineColor(kGreen+2);
  //hdEdxPionTot->SetLineColor(kBlue+2);

  const float frac=0.2;
  int bin1=0,bin2=0;

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
  pionmeanTot = pionfit->GetParameter(1);
  pionsigmaTot = pionfit->GetParameter(2);
//  electronmeanTot = electronfit->GetParameter(1);
//  electronsigmaTot = electronfit->GetParameter(2);

  pionres = pionsigmaTot/pionmeanTot;
  //electronres = electronsigmaTot/electronmeanTot;
  pionsigmaerr = pionfit->GetParError(2);
  pionmeanerr = pionfit->GetParError(1);
//  electronsigmaerr = electronfit->GetParError(2);
//  electronmeanerr = electronfit->GetParError(1);


//  separationpower = 2*(electronmeanTot-pionmeanTot)/(pionsigmaTot+electronsigmaTot);


  float pionreserror = sqrt(pow((pionsigmaerr/pionmeanTot),2)+pow(((pionsigmaTot*pionmeanerr)/pow(pionmeanTot,2)),2));
//  float electronreserror = sqrt(pow((electronsigmaerr/electronmeanTot),2)+pow(((electronsigmaTot*electronmeanerr)/pow(electronmeanTot,2)),2));

//  float separationpowererr = sqrt(pow((2*electronmeanerr/(electronsigmaTot+pionsigmaTot)),2)+pow((2*pionmeanerr/(electronsigmaTot+pionsigmaTot)),2)+
//        pow(((2*electronsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2)+pow(((2*pionsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2));

//  float elechisquareTot = electronfit->GetChisquare();
//  float eleNDFTot = electronfit->GetNDF();
  float pionchisquareTot = pionfit->GetChisquare();
  float pionNDFTot = pionfit->GetNDF();

  TPaveText *pave1=new TPaveText(0.6,.7,.9,.9,"NDC");
  pave1->SetBorderSize(1);
  pave1->SetFillColor(10);
  pave1->AddText("electron resolution");
  //pave1->AddText(Form("e: %.2f #pm %.2f (%.2f%%)",electronmeanTot,electronsigmaTot, electronsigmaTot/electronmeanTot*100));
  pave1->AddText(Form("#pi: %.2f #pm %.2f (%.2f%%)",pionmeanTot,pionsigmaTot,pionsigmaTot/pionmeanTot*100));
  pave1->AddText("separation power");
  //pave1->AddText(Form("Separation: %.2f#sigma", TMath::Abs(electronmeanTot-pionmeanTot)/((electronsigmaTot+pionsigmaTot)/2.)));
  pave1->Draw("same");
/**/
/*=============================================================================================================
 *Fit the dEdx plot for the maximal cluster charge
*/
  TCanvas *c16 = new TCanvas();
  TF1 *pionfitMax = new TF1("pionfitMax","gaus",hdEdxPionMax->GetXaxis()->GetXmin(),hdEdxPionMax->GetXaxis()->GetXmax());
//  TF1 *electronfitMax = new TF1("electronfitMax","gaus",hdEdxEleMax->GetXaxis()->GetXmin(),hdEdxEleMax->GetXaxis()->GetXmax());
//  hdEdxEleMax->SetLineColor(kGreen+2);
  hdEdxPionMax->SetLineColor(kBlue+2);

  GetBinMinMax(hdEdxPionMax,frac,bin1,bin2);
  hdEdxPionMax->Fit("pionfitMax","","",hdEdxPionMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionMax->GetXaxis()->GetBinUpEdge(bin2));
//  GetBinMinMax(hdEdxEleMax,frac,bin1,bin2);
//  hdEdxEleMax->Fit("electronfitMax","","",hdEdxEleMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleMax->GetXaxis()->GetBinUpEdge(bin2));



//  hdEdxEleMax->GetFunction("electronfitMax")->SetLineColor(kGreen+2);
  hdEdxPionMax->GetFunction("pionfitMax")->SetLineColor(kBlue+2);
  hdEdxPionMax->Draw();
//  hdEdxEleMax->Draw("same");
  pionmeanMax = pionfitMax->GetParameter(1);
  pionsigmaMax = pionfitMax->GetParameter(2);
//  float electronmeanMax = electronfitMax->GetParameter(1);
//  float electronsigmaMax = electronfitMax->GetParameter(2);

  pionresMax = pionsigmaMax/pionmeanMax;
//  electronresMax = electronsigmaMax/electronmeanMax;
  pionsigmaerrMax = pionfitMax->GetParError(2);
  pionmeanerrMax = pionfitMax->GetParError(1);
//  electronsigmaerrMax = electronfitMax->GetParError(2);
//  electronmeanerrMax = electronfitMax->GetParError(1);


//  separationpowerMax = 2*(electronmeanMax-pionmeanMax)/(pionsigmaMax+electronsigmaMax);


  float pionreserrorMax = sqrt(pow((pionsigmaerrMax/pionmeanMax),2)+pow(((pionsigmaMax*pionmeanerrMax)/pow(pionmeanMax,2)),2));
//  float electronreserrorMax = sqrt(pow((electronsigmaerrMax/electronmeanMax),2)+pow(((electronsigmaMax*electronmeanerrMax)/(electronmeanMax,2)),2));

//  float separationpowererrMax = sqrt(pow((2*electronmeanerrMax/(electronsigmaMax+pionsigmaMax)),2)+pow((2*pionmeanerrMax/(electronsigmaMax+pionsigmaMax)),2)+
//        pow(((2*electronsigmaerrMax*(pionmeanMax-electronmeanMax))/(pow((electronsigmaMax+pionsigmaMax),2))),2)+pow(((2*pionsigmaerrMax*(pionmeanMax-electronmeanMax))/(pow((electronsigmaMax+pionsigmaMax),2))),2));

//  float elechisquareMax = electronfitMax->GetChisquare();
//  float eleNDFMax = electronfitMax->GetNDF();
  float pionchisquareMax = pionfitMax->GetChisquare();
  float pionNDFMax = pionfitMax->GetNDF();

  TPaveText *pave3=new TPaveText(0.6,.7,.9,.9,"NDC");
  pave3->SetBorderSize(1);
  pave3->SetFillColor(10);
  pave3->AddText("electron resolution");
  //pave3->AddText(Form("e: %.2f #pm %.2f (%.2f%%)",electronmeanMax,electronsigmaMax, electronsigmaMax/electronmeanMax*100));
  pave3->AddText(Form("#pi: %.2f #pm %.2f (%.2f%%)",pionmeanMax,pionsigmaMax,pionsigmaMax/pionmeanMax*100));
  pave3->AddText("separation power");
  //pave3->AddText(Form("Separation: %.2f#sigma", TMath::Abs(electronmeanMax-pionmeanMax)/((electronsigmaMax+pionsigmaMax)/2.)));
  pave3->Draw("same");

/**/

  TCanvas *dEdxMaxNoCut = new TCanvas();
  hdEdxMaxNOCut->Draw();

  TCanvas *dEdxNOCut = new TCanvas();
  hdEdxTotNOCut->Draw();

  TCanvas *dEdxNTrackCut = new TCanvas();
  hdEdxTotAfterNTrackCut->Draw();

  TCanvas *dEdxROCCut = new TCanvas();
  hdEdxTotAfterROCCut->Draw();

  TCanvas *dEdxNClCut = new TCanvas();
  hdEdxTotAfterNClCut->Draw();

  TCanvas *dEdxNClEdgeCut = new TCanvas();
  hdEdxTotAfterNClEdgeCut->Draw();

  TCanvas *dEdxTimeCut = new TCanvas();
  hdEdxTotAfterTimeCut->Draw();

  TProfile *tpf = hQmaxRelPadPos->ProfileX("",0,100);

  TCanvas *QmaxRelPadPos = new TCanvas();
  tpf->SetLineColor(kBlack);
  hQmaxRelPadPos->Draw("colz");
  tpf->Draw("same");

  if (showSec){
    TCanvas *SecClusters = new TCanvas("Secotor Clusters", "Individual Sectors", 1920, 1080);
    SecClusters->SetTitle("All Sectors");
    SecClusters->Divide(6,6);
    for (int ipad = 1; ipad<=nsectors; ++ipad){
      SecClusters->cd(ipad);
      hOcc[ipad-1]->Draw("colz");
    }

    TCanvas *AFirstHalf = new TCanvas("Sector Clusters A first half", "Individual Sectors A first half",1920,1080);
    AFirstHalf->SetTitle("A side upper half");
    AFirstHalf->Divide(3,3);
    for (int ipad = 1; ipad<=9; ++ipad){
      AFirstHalf->cd(ipad);
      hOcc[ipad-1]->Draw("colz");
    }
    TCanvas *ASecondHalf = new TCanvas("Sector Clusters A second half", "Individual Sectors A second half",1920,1080);
    ASecondHalf->SetTitle("A side lower half");
    ASecondHalf->Divide(3,3);
    for (int ipad = 1; ipad<=9; ++ipad){
      ASecondHalf->cd(ipad);
      hOcc[ipad+8]->Draw("colz");
    }
    TCanvas *CFirstHalf = new TCanvas("Sector Clusters C first half", "Individual Sectors C first half",1920,1080);
    CFirstHalf->SetTitle("C side upper half");
    CFirstHalf->Divide(3,3);
    for (int ipad = 1; ipad<=9; ++ipad){
      CFirstHalf->cd(ipad);
      hOcc[ipad+17]->Draw("colz");
    }
    TCanvas *CSecondHalf = new TCanvas("Sector Clusters C second half", "Individual Sectors C second half",1920,1080);
    CSecondHalf->SetTitle("C side lower half");
    CSecondHalf->Divide(3,3);
    for (int ipad = 1; ipad<=9; ++ipad){
      CSecondHalf->cd(ipad);
      hOcc[ipad+26]->Draw("colz");
    }

    TCanvas *QmaxTime = new TCanvas();
    hQmaxTimeMean->Draw("colz");



//    TCanvas *Sector17Pad = new TCanvas();
//    PadSec17->Draw("colz");
  }
  TrackAna->Fill();
  OutFile->Write();
  OutFile->Close();

  return 0;
}
