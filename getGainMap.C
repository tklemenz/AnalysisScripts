void GetBinMinMax(const TH1 *hist, const Float_t frac, Int_t &bin1, Int_t &bin2)
{

  const Int_t binMax=hist->GetMaximumBin();
  const Double_t contMax=hist->GetBinContent(binMax);
  bin1=binMax;
  bin2=binMax;
  while ( (bin1--)>binMax/3. ) if (hist->GetBinContent(bin1)<frac*contMax) break;
  while ( (bin2++)<binMax*3. ) if (hist->GetBinContent(bin2)<frac*contMax) break;
}



void getGainMap(const char *filename){
  
  gSystem->Load("~/GEMtest/code/libGEMEvent.so");

  TFile *TreeFile = new TFile(Form("%s", filename));
  TTree *tree = (TTree*)TreeFile->Get("GEM");
  
  
  
  /* ============================================================
   * ====================== Cut parameter =======================
   * ============================================================ */
   
   
  int TrPerEv = 1;													
  int nclCut = 32;
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
  
  
  TH1D *hChisquare   = new TH1D("hChisquare", "; Chisquare Landau; # counts", 24,0,12);
  TH2D *hGainMap     = new TH2D("hGainMap", "; Row; Pad", 63,0,63,35,-5,30);
  
  TH1D *hEntries     = new TH1D("Entries", "; Entries; # counts",6,1,3000);
  
  
  //TH1D *hResidualPad = new TH1D("hResPad", "; Residual; #Counts",100,0,100);
  
  //TH2D *hQtotTimeBin = new TH2D("hQTime","; TimeBin; Q_{tot}",100,0,100,200,0,600);
  //TH2D *hQtotTimeBinUsed = new TH2D("hQTimeUsed","; TimeBin; Q_{tot}",100,0,100,200,0,600);
  
  
  
  
  /* initialize array of histograms used to get the gainmap */
  const int nrows = 63;
  const int npads = 30;
  TH1D *hQtotPads[nrows][npads];

  for(int irow=0; irow<nrows; ++irow) {
	  for (int ipad=0; ipad<npads; ++ipad) {
	  hQtotPads[irow][ipad] = new TH1D(Form("hQtotPad_%i_%i", irow, ipad), "", 50, 0, 50);
	  }
  }

/* =============================================================================
 * ================== Loop over events and apply all cuts ======================
 * ======================= meanwhile fill histograms ===========================
 * ============================================================================= */
 
 
  GEMEvent *fEvent=0x0;
  fEvent=0x0;
  tree->SetBranchAddress("Events",&fEvent);

  GEMTrack::SetRemoveRowsdEdx("9,17,25,32,39,45,51,57,58,59,60,61,62");
  
  int usedTracks = 0;
  int Tracks = 0;
  int OneTrackEvents = 0;
  int pions = 0;
  int electrons = 0;
  int usedcl = 0;
  int onetrcl = 0;
  
  
  
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
	      
	      
	      
	    }	  
	  }  
    }
  }
  
  
  
  /* ======================================================================
   * ========== Fit every histogram in array for MPV -> GainMap ===========
   * ====================================================================== */
   
  //TFile *GainMap = new TFile("GainMap.root","RECREATE");
  //hGainMap->SetName("GainMap");
   
  TF1 *MPVfit = new TF1("MPV","landau"); 
  int means = 0;
  int Landau = 0;
  const Float_t frac=0.2;
  Int_t bin1=0,bin2=0;
  for (int irow=0; irow<nrows; ++irow){
    for (int ipad=0; ipad<npads; ++ipad){
	  if (irow == 15 && ipad == 15){
	    //TCanvas *c1 = new TCanvas();
	    hQtotPads[irow][ipad]->Draw();
	  }
	  int entries = hQtotPads[irow][ipad]->GetEntries();
	  hEntries->Fill(entries);
	  if (entries > 100){	
		GetBinMinMax(hQtotPads[irow][ipad],frac,bin1,bin2);
        hQtotPads[irow][ipad]->Fit("MPV","","", hQtotPads[irow][ipad]->GetXaxis()->GetBinLowEdge(bin1),hQtotPads[irow][ipad]->GetXaxis()->GetBinUpEdge(bin2));
        float mpv = MPV->GetParameter(1);        
        float chisquare = MPV->GetChisquare()/MPV->GetNDF();
        hChisquare->Fill(chisquare);
        hGainMap->Fill(irow,ipad,mpv);
        ++Landau;

	  }  
	  else if (entries <= 100 && entries != 0){
	    hGainMap->Fill(irow,ipad,hQtotPads[irow][ipad]->GetMean());
        ++means;

	  }            
    }
  }
    
  
  hGainMap->GetZaxis()->SetRangeUser(.5,1.2);
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
  
  TCanvas *c10 = new TCanvas("GainMap");
  hGainMap->Draw("colz");
  
  TCanvas *c11 = new TCanvas();
  hChisquare->Draw();
  
  TCanvas *c12 = new TCanvas();
  hEntries->Draw();
  
  //c10->Update();
  //c10->Print("GainMap.root","root");
  
  /* ======================================================================
   * ======================= Terminal output ==============================
   * ====================================================================== */
  
  cout<<"Number of tracks: "<<Tracks<<endl;
  cout<<"Number of one-track-events: "<<OneTrackEvents<<endl;
  cout<<"Number of used tracks: "<<usedTracks<<endl;
  //cout<<"Number of one-track-clusters: "<<onetrcl<<endl;
  //cout<<"Number of used clusters: "<<usedcl<<endl;
  cout<<"Number of pions: "<<pions<<endl;
  cout<<"Number of electrons: "<<electrons<<endl;
  cout<<"Number of Landau fits: "<<Landau<<endl;
  cout<<"Number of means: "<<means<<endl;
}
