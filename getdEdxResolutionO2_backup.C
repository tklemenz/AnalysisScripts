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
#include "TStyle.h"
#include "TPaveText.h"

//#pragma link C++ class std::vector<o2::TPC::TrackTPC>+;

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



int getdEdxResolutionO2(TString argv)
{
  gStyle->SetOptStat(0);

  EventHeader Header;

  TString filename(argv);
  ///int dEdxExcludeRows(exRows);

  TFile *TreeFile = TFile::Open(filename.Data());
  cout<<endl<<endl<<filename<<endl<<endl;
  TTree *tree = (TTree*)gDirectory->Get("events");
  
  ///TH2D *fGainMap = 0x0;
  
  ///TFile *fGain = new TFile("~/GEMtest/code/myana/GainMap_FilesForGainMap.txt_2017-04-24_1759_0_07_trunc/GainMap_Trunc.root");

  ///fGainMap = (TH2D*)fGain->Get("GainMap");

  //cout<<"DEBUG POINT 1"<<endl;
  
  
  
  
  /* ============================================================
   * ====================== Cut parameter =======================
   * ============================================================ */
   
   
  int TrPerEv = 1;													
  int nclCut = 32;
  float CherCutLow = 0.009;
  float CherCutHigh = 0.011;
  int timeMaxLow = 25;
  int timeMaxHigh = 90;
  int cpadLow = 12;
  int cpadHigh = 16;
  float nclFracoutofCPad = .3;
  //int dEdxExcludeRows = 0;
  int excludeEdge = 1;//(dEdxExcludeRows/10)%10;
  //dEdxExcludeRows%=10;
  
  
  
  /* ============================================================
   * ============= Define histograms and graphs =================
   * ============================================================*/
   
  
  //TH1D *hNclusters   = new TH1D("hNclusters", ";Number of clusters; Counts", 65,0,65);
  TH1D *hNclustersUsed     = new TH1D("hNclustersUsed", ";Number of clusters per track; Counts", 300,0,300);

  
  TH1D *hNIROCtracks       = new TH1D("hNtracksIROC","; Number of tracks; # counts",20,0,20);
  TH1D *hNIROCtracksUsed   = new TH1D("hNtracksIROCUsed","; Number of tracks per event; # counts",50,0,50);

  TH1F *hdEdxEleTot        = new TH1F("hdEdxEleTot", "; d#it{E}/d#it{x} (a.u.); # counts", 200, 0, 400);
  TH1F *hdEdxPionTot       = new TH1F("hdEdxTot", "; d#it{E}/d#it{x} (a.u.); # counts", 200, 0, 400);
  
  TH1F *hdEdxEleMax        = new TH1F("hdEdxEleMax", "; d#it{E}/d#it{x} Q_{max} (a.u.); # counts", 100, 0, 120);
  TH1F *hdEdxPionMax       = new TH1F("hdEdxMax", "; d#it{E}/d#it{x} Q_{max} (a.u.); # counts", 100, 0, 120);
  
  TH1D *hQ              = new TH1D("hQ", "; Total cluster charge Q [ADC counts]; # counts", 600,0,600);
  TH1D *hQUsed          = new TH1D("hQUsed", "; Total cluster charge Q_{tot} [ADC counts]; # counts", 600,0,600);
  
  TH1D *hCherenkov         = new TH1D("hCherenkov", "; ADC signal; # counts", 100,0,0.06);
  TH1D *hCherenkovUsedPions     = new TH1D("hCherenkovUsed", "; ADC signal; # counts", 100,0,0.06);
  TH1D *hCherenkovUsedEle     = new TH1D("hCherenkovUsedEle", "; ADC signal; # counts", 100,0,0.06);
  
  TH2D *hPadOccupancy      = new TH2D("hOcc", "; Row; Pad", 250,-10,240,300,-100,200);
  TH2D *hPadOccupancyUsed  = new TH2D("hOccUsed", "; Row; Pad", 63,0,63,35,-5,30);
  TH2D *hPadOccupancyPad      = new TH2D("hOcc", "; Row; Pad", 250,-10,240,300,-100,200);

  
  TH2D *hQTimeMax       = new TH2D("hQTime","; TimeBin; Q_{tot}",20,0,20,300,0,1000);
  TH2D *hQTimeMaxUsed   = new TH2D("hQTimeUsed","; TimeBin; Q_{tot}",20,0,20,300,0,1000);
  


/* =============================================================================
 * ================== Loop over events and apply all cuts ======================
 * ======================= meanwhile fill histograms ===========================
 * ============================================================================= */
 
 
  std::vector<TrackTPC> *vecEvent = 0;
  //std::vector<TrackTPC> *pEvent = &Event;
  //std::vector<TrackTPC> *pEvent = 0;
  //cout<<"DEBUG POINT 2"<<endl;

  tree->SetBranchAddress("Tracks", &vecEvent);
  tree->SetBranchAddress("header", &Header);

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
  float CherenkovValue = 0;
  int runNr = 0;
  
  //cout<<"DEBUG POINT 3"<<endl;

  
  for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
    tree->GetEntry(iev);
    
    int NTracks = vecEvent->size();
    CherenkovValue = Header.cherenkovValue;
    runNr = Header.run;
    
    Tracks += NTracks;
    
    hNIROCtracks->Fill(NTracks);
    
    
    
/*========================================== CUT ================================================*/
    if (NTracks != TrPerEv) continue;																// only one-track events
    
        
    ++OneTrackEvents;	 
    //cout<<"DEBUG POINT 4"<<endl;

  for (auto trackObject : *vecEvent) {
    std::vector<Cluster> clCont;
    trackObject.getClusterVector(clCont);


/*========================================== CUT ================================================*/
      ///if (track->GetROC() != 0) continue;															// only ROC 0
      

      int ncl = clCont.size();
      bool isok = true;     
      int ncledge = 0;
	  
      onetrcl += ncl;

      for (auto clusterObject : clCont) {                     // loop over the clusters of all one-track-events
        DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
        float row = pos.getPadSecPos().getPadPos().getRow();
        float pad = pos.getPadSecPos().getPadPos().getPad();
        const int cpad = GetCPad(clusterObject);
        ushort timeMean = clusterObject.getTimeMean();
        double Q = clusterObject.getQ();
        hQTimeMax->Fill(timeMean,Q);
        hQ->Fill(Q);
        hPadOccupancy->Fill(row,cpad);
      }
      
      
      //cout<<"DEBUG POINT 5"<<endl;

	  
/*========================================== CUT ================================================*/     
     if (ncl < nclCut) continue;																	// cut on number of clusters per track
	      
	  									

      hCherenkov->Fill(CherenkovValue);
/*========================================== CUT ================================================*/
    if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;					// PID via Cherenkov
    if (CherenkovValue < CherCutLow) {
      hCherenkovUsedPions->Fill(CherenkovValue);
    }
    else if (CherenkovValue > CherCutHigh){
      hCherenkovUsedEle->Fill(CherenkovValue);
    }
     
      
	  
    for (auto clusterObject : clCont) {															// make cuts on cluster properties
      DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
      float row = pos.getPadSecPos().getPadPos().getRow();
      float pad = pos.getPadSecPos().getPadPos().getPad();
      const int cpad = GetCPad(clusterObject);
      ushort timeMean = clusterObject.getTimeMean();
      double Q = clusterObject.getQ();
	    double gainMapCorr = 1.;
	    
	    //cout<<endl<<"Row: "<<double(row)<<"\t"<<"Pad: "<<double(pad)<<"\t"<<"CPad: "<<double(cpad)<<endl;
	    		

/*========================================== CUT ================================================*/
    ///if (timeMax < timeMaxLow || timeMax > timeMaxHigh){isok = false; break;}					// cut on time max, only clusters from a certain range on z-axis can come from a particle
		
		
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
    if (((float(ncledge)/float(ncl)) > nclFracoutofCPad) && excludeEdge) continue;					// cut on detector geometry
	  

	  if (isok == true){																			// track accepted
      usedcl += ncl;
      float dEdxTot = trackObject.getTruncatedMean(.05,.7,1);
      float dEdxMax = trackObject.getTruncatedMean(.05,.7,0);
		  																			
      if (CherenkovValue < CherCutLow){
            hdEdxPionTot->Fill(dEdxTot);
            hdEdxPionMax->Fill(dEdxMax);
            ++pions;
            pioncl += ncl;
      }
      if (CherenkovValue > CherCutHigh){
            hdEdxEleTot->Fill(dEdxTot);
            hdEdxEleMax->Fill(dEdxMax);
            ++electrons;
            elecl += ncl;
      }
    //hdEdxPionTot->Fill(dEdxTot);
		++usedTracks;
		hNIROCtracksUsed->Fill(NTracks);
		hNclustersUsed->Fill(ncl);
    for (auto clusterObject : clCont) {															// loop over clusters
      DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
      float row = pos.getPadSecPos().getPadPos().getRow();
      float pad = pos.getPadSecPos().getPadPos().getPad();
      const int cpad = GetCPad(clusterObject);
      ushort timeMean = clusterObject.getTimeMean();
      double Q = clusterObject.getQ();

      hQTimeMaxUsed->Fill(timeMean,Q);
      hQUsed->Fill(Q);
      hPadOccupancyUsed->Fill(row,cpad);
      hPadOccupancyPad->Fill(row,pad);
	      
	    }	  
	  }  
  } ///break; ///if (isok==true){break;} //use only one track
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
  ///hQTimeMax->GetZaxis()->SetRangeUser(0,820);
  hQTimeMax->Draw("colz");
  
  TCanvas *c5 = new TCanvas();
  ///hQTimeMaxUsed->GetZaxis()->SetRangeUser(0,820);
  hQTimeMaxUsed->Draw("colz");
  
  TCanvas *c6 = new TCanvas();
  hQ->Draw();
  
  TCanvas *c7 = new TCanvas();
  hQUsed->Draw();
  
  TCanvas *c8 = new TCanvas();
  ///hPadOccupancy->GetZaxis()->SetRangeUser(0,3700);
  hPadOccupancy->Draw("colz");
  
  TCanvas *c9 = new TCanvas();
  ///hPadOccupancyUsed->GetZaxis()->SetRangeUser(0,3700);
  hPadOccupancyUsed->Draw("colz");
  
  TCanvas *c10 = new TCanvas();
  c10->SetLogy();
  hCherenkov->SetLineColor(kBlue);
  hCherenkov->SetLineWidth(2);
  hCherenkov->Draw();
  
  TCanvas *c11 = new TCanvas();
  c11->SetLogy();
  hCherenkovUsedPions->GetYaxis()->SetRangeUser(0.5,21000);
  hCherenkovUsedPions->SetLineColor(kBlue);
  hCherenkovUsedEle->SetLineColor(kRed);
  hCherenkovUsedPions->SetLineWidth(2);
  hCherenkovUsedEle->SetLineWidth(2);
  
  hCherenkovUsedPions->Draw();
  hCherenkovUsedEle->Draw("same");

  TCanvas *OccPad = new TCanvas();
  hPadOccupancyPad->Draw("colz");
  
  TCanvas *c15 = new TCanvas();
  TF1 *pionfit = new TF1("pionfit","gaus",hdEdxPionTot->GetXaxis()->GetXmin(),hdEdxPionTot->GetXaxis()->GetXmax());
  TF1 *electronfit = new TF1("electronfit","gaus",hdEdxEleTot->GetXaxis()->GetXmin(),hdEdxEleTot->GetXaxis()->GetXmax());
  hdEdxEleTot->SetLineColor(kRed);
  hdEdxPionTot->SetLineColor(kBlue);
  
  const Float_t frac=0.2;
  Int_t bin1=0,bin2=0;
  
  GetBinMinMax(hdEdxPionTot,frac,bin1,bin2); 
  hdEdxPionTot->Fit("pionfit","","",hdEdxPionTot->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionTot->GetXaxis()->GetBinUpEdge(bin2));
  GetBinMinMax(hdEdxEleTot,frac,bin1,bin2);
  hdEdxEleTot->Fit("electronfit","","",hdEdxEleTot->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleTot->GetXaxis()->GetBinUpEdge(bin2));
  
  //alternative fit
  //hdEdxPionTot->Fit("pionfit");
  //hdEdxEleTot->Fit("electronfit");
  
  
  hdEdxEleTot->GetFunction("electronfit")->SetLineColor(kRed);
  hdEdxPionTot->GetFunction("pionfit")->SetLineColor(kBlue);
  hdEdxPionTot->Draw();
  hdEdxEleTot->Draw("same");
  float pionmeanTot = pionfit->GetParameter(1);
  float pionsigmaTot = pionfit->GetParameter(2);
  float electronmeanTot = electronfit->GetParameter(1);
  float electronsigmaTot = electronfit->GetParameter(2);
  
  float pionres = pionsigmaTot/pionmeanTot;
  float electronres = electronsigmaTot/electronmeanTot;
  float pionsigmaerr = pionfit->GetParError(2);
  float pionmeanerr = pionfit->GetParError(1);
  float electronsigmaerr = electronfit->GetParError(2);
  float electronmeanerr = electronfit->GetParError(1);
  
  
  float separationpower = 2*(electronmeanTot-pionmeanTot)/(pionsigmaTot+electronsigmaTot);
 
  
  float pionreserror = sqrt(pow((pionsigmaerr/pionmeanTot),2)+pow(((pionsigmaTot*pionmeanerr)/pow(pionmeanTot,2)),2));
  float electronreserror = sqrt(pow((electronsigmaerr/electronmeanTot),2)+pow(((electronsigmaTot*electronmeanerr)/pow(electronmeanTot,2)),2));
 
  float separationpowererr = sqrt(pow((2*electronmeanerr/(electronsigmaTot+pionsigmaTot)),2)+pow((2*pionmeanerr/(electronsigmaTot+pionsigmaTot)),2)+
        pow(((2*electronsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2)+pow(((2*pionsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2));
  
  float elechisquareTot = electronfit->GetChisquare();
  float eleNDFTot = electronfit->GetNDF();
  float pionchisquareTot = pionfit->GetChisquare();
  float pionNDFTot = pionfit->GetNDF();

  TPaveText *pave1=new TPaveText(0.6,.7,.9,.9,"NDC");
  pave1->SetBorderSize(1);
  pave1->SetFillColor(10);
  pave1->AddText(Form("e: %.2f #pm %.2f (%.2f%%)",electronmeanTot,electronsigmaTot, electronsigmaTot/electronmeanTot*100));
  pave1->AddText(Form("#pi: %.2f #pm %.2f (%.2f%%)",pionmeanTot,pionsigmaTot,pionsigmaTot/pionmeanTot*100));
  pave1->AddText(Form("Separation: %.2f#sigma", TMath::Abs(electronmeanTot-pionmeanTot)/((electronsigmaTot+pionsigmaTot)/2.)));
  pave1->Draw("same");

  
  
  TCanvas *c16 = new TCanvas();
  TF1 *pionfitMax = new TF1("pionfitMax","gaus",hdEdxPionMax->GetXaxis()->GetXmin(),hdEdxPionMax->GetXaxis()->GetXmax());
  TF1 *electronfitMax = new TF1("electronfitMax","gaus",hdEdxEleMax->GetXaxis()->GetXmin(),hdEdxEleMax->GetXaxis()->GetXmax());
  hdEdxEleMax->SetLineColor(kRed);
  hdEdxPionMax->SetLineColor(kBlue);
  
  GetBinMinMax(hdEdxPionMax,frac,bin1,bin2);
  hdEdxPionMax->Fit("pionfitMax","","",hdEdxPionMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionMax->GetXaxis()->GetBinUpEdge(bin2));
  GetBinMinMax(hdEdxEleMax,frac,bin1,bin2);
  hdEdxEleMax->Fit("electronfitMax","","",hdEdxEleMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleMax->GetXaxis()->GetBinUpEdge(bin2));
  
  
  
  hdEdxEleMax->GetFunction("electronfitMax")->SetLineColor(kRed);
  hdEdxPionMax->GetFunction("pionfitMax")->SetLineColor(kBlue);
  hdEdxPionMax->Draw();
  hdEdxEleMax->Draw("same");
  float pionmeanMax = pionfitMax->GetParameter(1);
  float pionsigmaMax = pionfitMax->GetParameter(2);
  float electronmeanMax = electronfitMax->GetParameter(1);
  float electronsigmaMax = electronfitMax->GetParameter(2);
  
  float pionresMax = pionsigmaMax/pionmeanMax;
  float electronresMax = electronsigmaMax/electronmeanMax;
  float pionsigmaerrMax = pionfitMax->GetParError(2);
  float pionmeanerrMax = pionfitMax->GetParError(1);
  float electronsigmaerrMax = electronfitMax->GetParError(2);
  float electronmeanerrMax = electronfitMax->GetParError(1);
  
  
  float separationpowerMax = 2*(electronmeanMax-pionmeanMax)/(pionsigmaMax+electronsigmaMax);
 
  
  float pionreserrorMax = sqrt(pow((pionsigmaerrMax/pionmeanMax),2)+pow(((pionsigmaMax*pionmeanerrMax)/pow(pionmeanMax,2)),2));
  float electronreserrorMax = sqrt(pow((electronsigmaerrMax/electronmeanMax),2)+pow(((electronsigmaMax*electronmeanerrMax)/(electronmeanMax,2)),2));
 
  float separationpowererrMax = sqrt(pow((2*electronmeanerrMax/(electronsigmaMax+pionsigmaMax)),2)+pow((2*pionmeanerrMax/(electronsigmaMax+pionsigmaMax)),2)+
        pow(((2*electronsigmaerrMax*(pionmeanMax-electronmeanMax))/(pow((electronsigmaMax+pionsigmaMax),2))),2)+pow(((2*pionsigmaerrMax*(pionmeanMax-electronmeanMax))/(pow((electronsigmaMax+pionsigmaMax),2))),2));
  
  float elechisquareMax = electronfitMax->GetChisquare();
  float eleNDFMax = electronfitMax->GetNDF();
  float pionchisquareMax = pionfitMax->GetChisquare();
  float pionNDFMax = pionfitMax->GetNDF();

  TPaveText *pave2=new TPaveText(0.6,.7,.9,.9,"NDC");
  pave2->SetBorderSize(1);
  pave2->SetFillColor(10);
  pave2->AddText(Form("e: %.2f #pm %.2f (%.2f%%)",electronmeanMax,electronsigmaMax, electronsigmaMax/electronmeanMax*100));
  pave2->AddText(Form("#pi: %.2f #pm %.2f (%.2f%%)",pionmeanMax,pionsigmaMax,pionsigmaMax/pionmeanMax*100));
  pave2->AddText(Form("Separation: %.2f#sigma", TMath::Abs(electronmeanMax-pionmeanMax)/((electronsigmaMax+pionsigmaMax)/2.)));
  pave2->Draw("same");
  
  cout<<endl;cout<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<endl;
  cout<<"============================= Total Charge ======================================="<<endl;
  cout<<"Pion dE/dx resolution:          ("<<pionres*100<<" +- "<<pionreserror*100<<") %"<<endl;
  cout<<"Electron dE/dx resolution:      ("<<electronres*100<<" +- "<<electronreserror*100<<") %"<<endl;
  cout<<"Separation power:               ("<<separationpower<<" +- "<<separationpowererr<<") sigma"<<endl;
  cout<<"Pion Chisquare:                  "<<pionchisquareTot/pionNDFTot<<endl;
  cout<<"Electron Chisquare:              "<<elechisquareTot/eleNDFTot<<endl;
  cout<<"Pion mean:                       "<<pionmeanTot<<" +- "<<pionsigmaTot<<endl;
  cout<<"Electron mean:                   "<<electronmeanTot<<" +- "<<electronsigmaTot<<endl;
  cout<<endl;
  cout<<"============================== Max Charge ========================================"<<endl;
  cout<<"Pion dE/dx resolution:          ("<<pionresMax*100<<" +- "<<pionreserrorMax*100<<") %"<<endl;
  cout<<"Electron dE/dx resolution:      ("<<electronresMax*100<<" +- "<<electronreserrorMax*100<<") %"<<endl;
  cout<<"Separation power:               ("<<separationpowerMax<<" +- "<<separationpowererrMax<<") sigma"<<endl;
  cout<<"Pion Chisquare:                  "<<pionchisquareMax/pionNDFMax<<endl;
  cout<<"Electron Chisquare:              "<<elechisquareMax/eleNDFMax<<endl;
  cout<<"Pion mean:                       "<<pionmeanMax<<" +- "<<pionsigmaMax<<endl;
  cout<<"Electron mean:                   "<<electronmeanMax<<" +- "<<electronsigmaMax<<endl;
  cout<<endl;
  cout<<"============================== Statistics ========================================"<<endl;
  cout<<"Number of tracks:                "<<Tracks<<endl;
  cout<<"Number of one-track-events:      "<<OneTrackEvents<<endl;
  cout<<"Number of used tracks:           "<<usedTracks<<endl;
  cout<<"Number of pions:                 "<<pions<<endl;
  cout<<"Number of electrons:             "<<electrons<<endl;
  cout<<endl;
  /*cout<<"Number of one-track-clusters:    "<<onetrcl<<endl;
  cout<<"Number of used clusters:         "<<usedcl<<endl;
  cout<<"Number of electron clusters:     "<<elecl<<endl;
  cout<<"Number of pion clusters:         "<<pioncl<<endl;*/
  cout<<"Average clusters per pion:       "<<pioncl/pions<<endl;
  cout<<"Average clusters per electron:   "<<elecl/electrons<<endl;
  cout<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<endl;cout<<endl;
  
  //delete TreeFile;
  //if (fGain1) delete fGain1;
  //if (fGain2) delete fGain2;
  /*delete hNIROCtracks;
  delete hNIROCtracksUsed;
  delete hdEdxEleTot;
  delete hdEdxPionTot;
  delete hdEdxEleMax;
  delete hdEdxPionMax;
  delete hQ;
  delete hQUsed;
  delete hCherenkov;
  delete hCherenkovUsedPions;
  delete hCherenkovUsedEle;
  delete hPadOccupancy;
  delete hPadOccupancyUsed;
  delete hChisquare;
  delete hGainMap;
  delete hQTimeMax;
  delete hQTimeMaxUsed;
  delete c1;
  delete c2;
  delete c3;
  delete c4;
  delete c5;
  delete c6;
  delete c7;
  delete c8;
  delete c9;
  delete c10;
  delete c11;
  delete c15;
  delete c16;
  delete pionfit;
  delete electronfit; */
  return 0;
}


