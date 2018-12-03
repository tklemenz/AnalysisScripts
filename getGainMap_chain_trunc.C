#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMath.h"

#include "GEMEvent.h"
#include "GEMTrack.h"
#include "GEMCluster.h"

using namespace std;

/*
int main(int argc,char** argv){

  const char *FileList = argv[1];
  return 0;

}
*/

void GetBinMinMax(const TH1 *hist, const Float_t frac, Int_t &bin1, Int_t &bin2)
{

  const Int_t binMax=hist->GetMaximumBin();
  const Double_t contMax=hist->GetBinContent(binMax);
  bin1=binMax;
  bin2=binMax;
  while ( (bin1--)>binMax/3. ) if (hist->GetBinContent(bin1)<frac*contMax) break;
  while ( (bin2++)<binMax*3. ) if (hist->GetBinContent(bin2)<frac*contMax) break;
}

void TruncatedMean(TH1F *his, TVectorD *param, Float_t down, Float_t up)
{
  Int_t nbins = his->GetNbinsX();
  Float_t nentries = his-> GetEntries();
  Float_t sum = 0;
  Float_t mean = 0;
  Float_t sigma2 = 0;
  Float_t ncumul = 0;
  for (Int_t ibin = 1; ibin<nbins; ibin++){
    ncumul += his->GetBinContent(ibin);
    Float_t fraction = Float_t(ncumul)/Float_t(nentries);
    if (fraction>down && fraction<up){
      sum+=his->GetBinContent(ibin);
      mean+=his->GetBinCenter(ibin)*his->GetBinContent(ibin);
      sigma2+=his->GetBinCenter(ibin)*his->GetBinCenter(ibin)*his->GetBinContent(ibin);
    }
  }
  mean /= sum;
  sigma2 = TMath::Sqrt(TMath::Abs(sigma2/sum-mean*mean));
  if (param){
    (*param)[0] = his->GetMaximum();
    (*param)[1] = mean;
    (*param)[2] = sigma2; 
  }
}



int main(int argc,char** argv){
  
  const char *sarg = argv[1];
  const char *spath = argv[2];
  
  TString FileList(sarg);
  
  GEMEvent *fChainEvent=0x0;
 
  //fChain = new TChain("GEM");
  TChain fChain("GEM");
  
  fChain.SetBranchAddress("Events",&fChainEvent);
	
  TString allFiles;
  if (FileList.EndsWith(".txt")){
    allFiles=gSystem->GetFromPipe(Form("cat %s",FileList.Data()));
  }
  else {
	allFiles=gSystem->GetFromPipe(Form("ls %s",FileList.Data()));
  }
  TObjArray *arr = allFiles.Tokenize("\n");
	
  for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
	TString file=arr->At(ifile)->GetName();
	//cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
	fChain.Add(file);
  }
  cout<<endl<<endl<<endl<<"Chain ready!!"<<endl<<endl;
  TObjArray *chainEntries = fChain.GetListOfFiles();
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
  cout<<chainEntries->At(ifile)->GetTitle()<<endl;}
  
  
  /* ============================================================
   * ====================== Cut parameter =======================
   * ============================================================ */
   
   
  int TrPerEv = 1;											
  int nclCut = 30;
  int CherCutLow = 40;
  int CherCutHigh = 50;
  int timebinLow = 25;
  int timebinHigh = 90;
  
  
  
  
  /* ============================================================
   * ============= Define histograms and graphs =================
   * ============================================================*/
   
  
  //TH1D *hNclusters   = new TH1D("hNclusters", ";Number of clusters; Counts", 65,0,65);
  //TH1D *hNclustersUsed   = new TH1D("hNclustersUsed", ";Number of clusters; Counts", 65,0,65);

  
  //TH1D *hNIROCtracks = new TH1D("hNtracksIROC","; Number of tracks; # counts",10,0,10);
  //TH1D *hNIROCtracksUsed = new TH1D("hNtracksIROCUsed","; Number of tracks; # counts",10,0,10);

  
  //TH1D *hQtot       = new TH1D("hQtot", "; Total cluster charge Q_{tot} [ADC counts]; # counts", 600,0,600);
  //TH1D *hQtotUsed       = new TH1D("hQtotUsed", "; Total cluster charge Q_{tot} [ADC counts]; # counts", 600,0,600);
  
  //TH1D *hCherenkov   = new TH1D("hCherenkov", "; ADC signal (corresponds to velocity); # counts", 90,20,196);
  
  //TH2D *hPadOccupancy  = new TH2D("hOcc", "; Row; Pad", 63,0,63,35,-5,30);
  //TH2D *hPadOccupancyUsed  = new TH2D("hOccUsed", "; Row; Pad", 63,0,63,35,-5,30);
  
  //TH1D *hQtotSinglePad     = new TH1D("QtotSinglePad", "; Q_{tot}; # counts", 100,0,600);
  //TGraph *Chisquare  = new TGraph();
  
  
  TH1D *hChisquare   = new TH1D("hChisquare", "; Chisquare Landau; # counts", 400,0,200);
  TH2D *hGainMap     = new TH2D("hGainMap", "; Row; Pad", 74,-5,68,119,-59,59);
  
  TH2D *hGainMapTrunc     = new TH2D("hGainMapTrunc", "; Row; Pad", 74,-5,68,119,-59,59);
  
  TH1D *hEntries     = new TH1D("Entries", "; Entries; # counts",300,1,300000);
  
  
  //TH1D *hResidualPad = new TH1D("hResPad", "; Residual; #Counts",100,0,100);
  
  //TH2D *hQtotTimeBin = new TH2D("hQTime","; TimeBin; Q_{tot}",100,0,100,200,0,600);
  //TH2D *hQtotTimeBinUsed = new TH2D("hQTimeUsed","; TimeBin; Q_{tot}",100,0,100,200,0,600);
  
  
  
  
  /* initialize array of histograms used to get the gainmap */
  const int nrows = 63;
  const int npads = 30;
  TH1F *hQtotPads[nrows][npads];
  TH1F *hQtotPadsTrunc[nrows][npads];

  for(int irow=0; irow<nrows; ++irow) {
	  for (int ipad=0; ipad<npads; ++ipad) {
	  hQtotPads[irow][ipad] = new TH1F(Form("hQtotPad_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; #counts", 500, 0, 50);
          hQtotPadsTrunc[irow][ipad] = new TH1F(Form("hQtotPad_trunc_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; #counts", 500, 0, 50);
	  }
  }

/* =============================================================================
 * ================== Loop over events and apply all cuts ======================
 * ======================= meanwhile fill histograms ===========================
 * ============================================================================= */
 
 
  
  
  int usedTracks = 0;
  int Tracks = 0;
  int OneTrackEvents = 0;
  int pions = 0;
  int electrons = 0;
  int usedcl = 0;
  int onetrcl = 0;
  GEMEvent *fEvent=0x0;
  //TObjArray *chainEntries = fChain->GetListOfFiles();
  GEMTrack::SetRemoveRowsdEdx("9,17,25,32,39,45,51,57,58,59,60,61,62");

  cout<<endl<<endl<<endl<<"Number of files to process: "<<chainEntries->GetEntriesFast()<<endl<<endl<<endl;
  
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
	TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));	
	cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;//<<endl;
    TTree *tree = (TTree*)TreeFile->Get("GEM");
    //TTree *tree = (TTree*)chainEntries->At(ifile)->Get("GEM");
    
    
    fEvent=0x0;
    tree->SetBranchAddress("Events",&fEvent);
    //cout<<"KOMME BIS HIER!"<<endl;

  
    for (int iev=0; iev<tree->GetEntries(); ++iev){
      tree->GetEntry(iev);
      int NTracks = fEvent->GetNumberOfTracks(0);
      int CherenkovValue = fEvent->GetCherenkovValue();
      Tracks += NTracks;
    
/*===========================================CUT=================================================*/    
      if (NTracks != TrPerEv){continue;}																// only one-track events
    
    
      ++OneTrackEvents;	 
	  for (int itr=0; itr<NTracks; ++itr) {
        GEMTrack *track = const_cast <GEMTrack*>(fEvent->GetTrack(itr));
      
/*===========================================CUT=================================================*/      
        if (track->GetROC() != 0) continue;															// only ROC 0
      
      
        int ncl = track->GetNumberOfClusters();
        onetrcl += ncl;
      
        bool isok = true;
      
	    int ncledge = 0;

/*===========================================CUT=================================================*/    
        if (ncl < nclCut){isok = false; continue;}													// cut on number of clusters per track
	      
	  									

/*===========================================CUT=================================================*/
	    if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh){isok = false; continue;}	// PID via Cherenkov
      
	  
	    for (int icl = 0; icl < ncl; ++icl){															// make cuts on cluster properties
		  GEMCluster *cl = const_cast <GEMCluster*>(track->GetCluster(icl));
	      double timebin = cl->GetTimeBinW();
	    
/*===========================================CUT=================================================*/
	      if (timebin < timebinLow || timebin > timebinHigh){isok = false; break;}					// cut on time bins, only clusters from a certain range on z-axis can come from a particle		
		
									
	    }
	    if (isok == true){																			// track accepted
		  usedcl += ncl; 
		
		  float dEdx = track->GetTruncatedMean(0.,.7,1,2);
		  																			
		  if (CherenkovValue < CherCutLow){
            ++pions;
		  }
		  if (CherenkovValue > CherCutHigh){
            ++electrons;
		  }
		  ++usedTracks;

		  for (int icl = 0; icl < ncl; ++icl){														// loop over the clusters of an accepted track
		    GEMCluster *cl = const_cast <GEMCluster*>(track->GetCluster(icl));
	        const int row = cl->GetRow();
	        const int cpad = int(cl->GetCPad());
	        double QTot = cl->GetQTot();
	        hQtotPads[row][cpad]->Fill(QTot/dEdx);
	        hQtotPadsTrunc[row][cpad]->Fill(QTot/dEdx);
	      
	      
	      }	  
	    }  
      }
    }
    delete TreeFile;
    cout<<"Tree deleted"<<endl<<endl;
  }
  
  
  /* ======================================================================
   * ========== Fit every histogram in array for MPV -> GainMap ===========
   * ====================================================================== */
   
  TFile *f = TFile::Open(Form("%s/GainMap_Landau.root",spath),"RECREATE");
  TFile *chi = TFile::Open(Form("%s/ChiSquare.root",spath),"RECREATE");
  TFile *ent = TFile::Open(Form("%s/Entries.root",spath),"RECREATE");
  TFile *gTr = TFile::Open(Form("%s/GainMap_Trunc.root",spath),"RECREATE");
  //hGainMap->SetName("GainMap");
  //TCanvas *c1 = new TCanvas(); 
  
  TF1 *MPVfit = new TF1("MPV","landau"); 
  int means = 0;
  int Landau = 0;
  const Float_t frac=0.2;
  Int_t bin1=0,bin2=0;
  TVectorD param(4);
  Float_t lowTr = 0.;
  Float_t upTr = .75;
  for (int irow=0; irow<nrows; ++irow){
    for (int ipad=0; ipad<npads; ++ipad){
	  /*if (irow == 20 && ipad == 10){
		hQtotPads[irow][ipad]->SetLineWidth(2);
	    hQtotPads[irow][ipad]->SetLineColor(kBlue);
	    hQtotPads[irow][ipad]->Draw();
	  }*/
	  int entries = hQtotPads[irow][ipad]->GetEntries();
	  hEntries->Fill(entries);
	  if (entries > 40){	
		GetBinMinMax(hQtotPads[irow][ipad],frac,bin1,bin2);
        hQtotPads[irow][ipad]->Fit("MPV","NQ0","", hQtotPads[irow][ipad]->GetXaxis()->GetBinLowEdge(bin1),hQtotPads[irow][ipad]->GetXaxis()->GetBinUpEdge(bin2));
        float mpv = MPVfit->GetParameter(1);        
        float chisquare = MPVfit->GetChisquare()/MPVfit->GetNDF();
        hChisquare->Fill(chisquare);
        hGainMap->Fill(irow,ipad,mpv);
        ++Landau;


            param[0]=0;
            param[1]=0;
            param[2]=0;
            param[3]=0;
            TruncatedMean(hQtotPadsTrunc[irow][ipad],&param,lowTr,upTr);
            hGainMapTrunc->Fill(irow,ipad,param[1]);

	  }  
	  else if (entries <= 40 && entries != 0){
	    hGainMap->Fill(irow,ipad,hQtotPads[irow][ipad]->GetMean());
        ++means;
    
            param[0]=0;
            param[1]=0;
            param[2]=0;
            param[3]=0;
            TruncatedMean(hQtotPadsTrunc[irow][ipad],&param,lowTr,upTr);
            hGainMapTrunc->Fill(irow,ipad,param[1]);
	  }
         /* else if (entries > 0) {
            param[0]=0;
            param[1]=0;
            param[2]=0;
            param[3]=0;
            TruncatedMean(hQtotPadsTrunc[irow][ipad],&param,lowTr,upTr);
            hGainMapTrunc->Fill(irow,ipad,param[1]);
          }*/

    }
  }
    
  
  hGainMap->GetZaxis()->SetRangeUser(.3,1.3);
  hGainMapTrunc->GetZaxis()->SetRangeUser(.3,1.3);
  f->cd();
  hGainMap->Write("GainMap");
  delete f;
  chi->cd();
  hChisquare->Write("ChiSquare");
  delete chi;
  ent->cd();
  hEntries->Write("Entries");
  delete ent;
  gTr->cd();
  hGainMapTrunc->Write("GainMap");
  delete gTr;
  //hGainMap->Write();
  //GainMap->Write();
  //GainMap->Close();


  /* ======================================================================
   * ===================== Plot everything ================================
   * ====================================================================== */
  //TCanvas *c1 = new TCanvas();
  //hNIROCtracksUsed->Draw();
  
  //TCanvas *c2 = new TCanvas();
  //hNIROCtracks->Draw();
  
  //TCanvas *c3 = new TCanvas();
  //hNclustersUsed->Draw();
  
  //TCanvas *c4 = new TCanvas();
  //hQtotTimeBin->Draw("colz");
  
  //TCanvas *c5 = new TCanvas();
  //hQtotTimeBinUsed->Draw("colz");
  
  //TCanvas *c6 = new TCanvas();
  //hQtot->Draw();
  
  //TCanvas *c7 = new TCanvas();
  //hQtotUsed->Draw();
  
  //TCanvas *c8 = new TCanvas();
  //hPadOccupancy->Draw("colz");
  
  //TCanvas *c9 = new TCanvas();
  //hPadOccupancyUsed->Draw("colz");
  
  /*TCanvas *c10 = new TCanvas("GainMap");
  hGainMap->Draw("colz");
  
  TCanvas *c11 = new TCanvas();
  hChisquare->Draw();
  
  TCanvas *c12 = new TCanvas();
  hEntries->Draw(); */
  
  //c10->Update();
  //c10->Print("GainMap.root","root");
  
  /* ======================================================================
   * ======================= Terminal output ==============================
   * ====================================================================== */
  
  /*cout<<"Number of tracks: "<<Tracks<<endl;
  cout<<"Number of one-track-events: "<<OneTrackEvents<<endl;
  cout<<"Number of used tracks: "<<usedTracks<<endl;
  //cout<<"Number of one-track-clusters: "<<onetrcl<<endl;
  //cout<<"Number of used clusters: "<<usedcl<<endl;
  cout<<"Number of pions: "<<pions<<endl;
  cout<<"Number of electrons: "<<electrons<<endl;
  cout<<"Number of Landau fits: "<<Landau<<endl;
  cout<<"Number of means: "<<means<<endl; */
  cout<<endl<<endl<<endl<<"*** =============   done   ============= ***"<<endl<<endl<<endl;
  return 0;
}
