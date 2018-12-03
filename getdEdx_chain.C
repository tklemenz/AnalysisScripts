
void GetBinMinMax(const TH1 *hist, const Float_t frac, Int_t &bin1, Int_t &bin2)
{

  const Int_t binMax=hist->GetMaximumBin();
  const Double_t contMax=hist->GetBinContent(binMax);
  bin1=binMax;
  bin2=binMax;
  while ( (bin1--)>binMax/3. ) if (hist->GetBinContent(bin1)<frac*contMax) break;
  while ( (bin2++)<binMax*3. ) if (hist->GetBinContent(bin2)<frac*contMax) break;
}




void getClusters(TString FileList, Int_t dEdxExcludeRows=2){
  
  gSystem->Load("~/GEMtest/code/libGEMEvent.so");
  //gROOT->ProcessLine(".x rootlogon.C");

  GEMEvent *fEvent=0x0;
  fEvent=0x0;

  fChain = new TChain("GEM");
  //fChain = 0x0;
  fChain->SetBranchAddress("Events",&fEvent);
	
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
	fChain->Add(file);
  }
  cout<<endl<<endl<<endl<<"Chain ready!!"<<endl<<endl;
  TObjArray *chainEntries = fChain->GetListOfFiles();
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
  cout<<chainEntries->At(ifile)->GetTitle()<<endl;}
  
  
  
  
  TFile *fGain = new TFile("~/ROC_Raw_data/GainMap2014.root");
  //TFile *fGain = new TFile("~/c1_n2.root");
  fGainMap = (AliTPCCalROC*)fGain->Get("gainMap");
  
  
  
  
  /* ============================================================
   * ====================== Cut parameter =======================
   * ============================================================ */
   
   
  int TrPerEv = 1;													
  int nclCut = 32;
  int CherCutLow = 40;
  int CherCutHigh = 50;
  int timeMaxLow = 25;
  int timeMaxHigh = 90;
  int cpadLow = 3;
  int cpadHigh = 17;
  float nclFracoutofCPad = .3;
  //int dEdxExcludeRows = 0;
  Int_t excludeEdge = (dEdxExcludeRows/10)%10;
  dEdxExcludeRows%=10;
  
  
  
  /* ============================================================
   * ============= Define histograms and graphs =================
   * ============================================================*/
   
  
  //TH1D *hNclusters   = new TH1D("hNclusters", ";Number of clusters; Counts", 65,0,65);
  TH1D *hNclustersUsed     = new TH1D("hNclustersUsed", ";Number of clusters; Counts", 65,0,65);

  
  TH1D *hNIROCtracks       = new TH1D("hNtracksIROC","; Number of tracks; # counts",10,0,10);
  TH1D *hNIROCtracksUsed   = new TH1D("hNtracksIROCUsed","; Number of tracks; # counts",10,0,10);

  TH1F *hdEdxEleTot        = new TH1F("hdEdxEleTot", "; d#it{E}/d#it{x} Q_{tot} [a.u.]; # counts", 250, 0, 500);
  TH1F *hdEdxPionTot       = new TH1F("hdEdxTot", "; d#it{E}/d#it{x} Q_{tot} [a.u.]; # counts", 250, 0, 500);
  
  TH1F *hdEdxEleMax        = new TH1F("hdEdxEleMax", "; d#it{E}/d#it{x} Q_{tot} [a.u.]; # counts", 100, 0, 120);
  TH1F *hdEdxPionMax       = new TH1F("hdEdxMax", "; d#it{E}/d#it{x} Q_{tot} [a.u.]; # counts", 100, 0, 120);
  
  TH1D *hQtot              = new TH1D("hQtot", "; Total cluster charge Q_{tot} [ADC counts]; # counts", 600,0,600);
  TH1D *hQtotUsed          = new TH1D("hQtotUsed", "; Total cluster charge Q_{tot} [ADC counts]; # counts", 600,0,600);
  
  TH1D *hCherenkov         = new TH1D("hCherenkov", "; ADC signal; # counts", 90,20,196);
  TH1D *hCherenkovUsedPions     = new TH1D("hCherenkovUsedPions", "; ADC signal; # counts", 90,20,196);
  TH1D *hCherenkovUsedEle     = new TH1D("hCherenkovUsedEle", "; ADC signal; # counts", 90,20,196);
  
  TH2D *hPadOccupancy      = new TH2D("hOcc", "; Row; Pad", 63,0,63,35,-5,30);
  TH2D *hPadOccupancyUsed  = new TH2D("hOccUsed", "; Row; Pad", 63,0,63,35,-5,30);
  
  TH1D *hChisquare         = new TH1D("hChisquare", "; Chisquare Landau; # counts", 24,0,12);
  TH2D *hGainMap           = new TH2D("hGainMap", "; Row; Pad", 63,0,63,35,-5,30);
  
  TH2D *hQtotTimeMax       = new TH2D("hQTime","; TimeBin; Q_{tot}",105,0,105,200,0,600);
  TH2D *hQtotTimeMaxUsed   = new TH2D("hQTimeUsed","; TimeBin; Q_{tot}",105,0,105,200,0,600);
  


/* =============================================================================
 * ================== Loop over events and apply all cuts ======================
 * ======================= meanwhile fill histograms ===========================
 * ============================================================================= */
  int usedTracks = 0;
  int Tracks = 0;
  int OneTrackEvents = 0;
  float pions = 0;
  float electrons = 0;
  int usedcl = 0;
  int onetrcl = 0;
  float pioncl = 0;
  float elecl = 0;
 
  GEMEvent *fEvent=0x0;
  GEMTrack::SetRemoveRowsdEdx("9,17,25,32,39,45,51,57,58,59,60,61,62");
  
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
	TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));	
	cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;//<<endl;
    TTree *tree = (TTree*)TreeFile->Get("GEM");
    
    fEvent=0x0;
    tree->SetBranchAddress("Events",&fEvent); 
  
    for (int iev=0; iev<tree->GetEntries(); ++iev){
      tree->GetEntry(iev);
    
      int NTracks = fEvent->GetNumberOfTracks(0);
      int CherenkovValue = fEvent->GetCherenkovValue();
    
      Tracks += NTracks;
    
      hNIROCtracks->Fill(NTracks);
    
    
    
/*========================================== CUT ================================================*/
      if (NTracks != TrPerEv) continue;																// only one-track events
    
        
      ++OneTrackEvents;	 
	
   	  for (int itr=0; itr<NTracks; ++itr) {
        GEMTrack *track = const_cast <GEMTrack*>(fEvent->GetTrack(itr));
      
      
/*========================================== CUT ================================================*/
        if (track->GetROC() != 0) continue;															// only ROC 0
      
      
        int ncl = track->GetNumberOfClusters();
        bool isok = true;     
	    int ncledge = 0;
	  
        onetrcl += ncl;
      
        for (int icl = 0; icl < ncl; ++icl){															// loop over the clusters of all one-track-events
		  GEMCluster *cl = const_cast <GEMCluster*>(track->GetCluster(icl));
	      const int row = cl->GetRow();
	      const double cpad = cl->GetCPad();
	      UShort_t timeMax = cl->GetTimeMax();
	      double QTot = cl->GetQTot();
	      hQtotTimeMax->Fill(timeMax,QTot);
	      hQtot->Fill(QTot);
	      hPadOccupancy->Fill(row,cpad);
	    }
      
      
	  
	  
/*========================================== CUT ================================================*/     
        if (ncl < nclCut) continue;																	// cut on number of clusters per track
	      
	  									

        hCherenkov->Fill(CherenkovValue);
/*========================================== CUT ================================================*/
	    if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;					// PID via Cherenkov
	      if (CherenkovValue < CherCutLow) {
		    hCherenkovUsedPions->Fill(CherenkovValue);
		  }
		  else {
		    hCherenkovUsedEle->Fill(CherenkovValue);
		  }
      
      
	  
	    for (int icl = 0; icl < ncl; ++icl){															// make cuts on cluster properties
		  GEMCluster *cl = const_cast <GEMCluster*>(track->GetCluster(icl));
	      const int row = cl->GetRow();
	      const Int_t pad = TMath::Nint(cl->GetPad());
	      const double cpad = cl->GetCPad();
	      UShort_t timeMax = cl->GetTimeMax();
	      double gainMapCorr = 1.;
	    		

/*========================================== CUT ================================================*/
		  if (timeMax < timeMaxLow || timeMax > timeMaxHigh){isok = false; break;}					// cut on time max, only clusters from a certain range on z-axis can come from a particle
		
		
		  if (cpad < cpadLow || cpad > cpadHigh){++ncledge;}	
		
		  if (fGainMap){
	        gainMapCorr = fGainMap->GetValue(row,pad);
		  }
		  cl->CorrectCharge(gainMapCorr);																// gain correction from gainmap
	    }

/*========================================== CUT ================================================*/ 
        if (((float(ncledge)/float(ncl)) > nclFracoutofCPad) && excludeEdge) continue;					// cut on detector geometry
	  
	  
	    if (isok == true){																			// track accepted
		  usedcl += ncl; 
		  float dEdxTot = track->GetTruncatedMean(0.,.7,1,dEdxExcludeRows);
		  float dEdxMax = track->GetTruncatedMean(0.,.7,0,dEdxExcludeRows);
		  																			
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
		  ++usedTracks;
		  hNIROCtracksUsed->Fill(NTracks);
		  hNclustersUsed->Fill(ncl);
		  for (int icl = 0; icl < ncl; ++icl){														// loop over the clusters of an accepted track
		    GEMCluster *cl = const_cast <GEMCluster*>(track->GetCluster(icl));
	        const int row = cl->GetRow();
	        const double cpad = cl->GetCPad();
	        UShort_t timeMax = cl->GetTimeMax();
	        double QTot = cl->GetQTot();
	      
	        hQtotTimeMaxUsed->Fill(timeMax,QTot);
	        hQtotUsed->Fill(QTot);
	        hPadOccupancyUsed->Fill(row,cpad);
	      
	      }	  
	    }  
      } //if (isok==true){break;} //use only one track
    }
    delete TreeFile;
    cout<<"Tree deleted"<<endl<<endl;
 } 


  /* ======================================================================
   * ===================== Plot everything ================================
   * ====================================================================== */
   
  hQtot->SetLineColor(kBlue+2);
   
  TCanvas *c1 = new TCanvas();
  hNIROCtracksUsed->Draw();
  
  TCanvas *c2 = new TCanvas();
  hNIROCtracks->Draw();
  
  TCanvas *c3 = new TCanvas();
  hNclustersUsed->Draw();
  
  TCanvas *c4 = new TCanvas();
  //hQtotTimeMax->GetZaxis()->SetRangeUser(0,820);
  hQtotTimeMax->Draw("colz");
  
  TCanvas *c5 = new TCanvas();
  //hQtotTimeMaxUsed->GetZaxis()->SetRangeUser(0,820);
  hQtotTimeMaxUsed->Draw("colz");
  
  TCanvas *c6 = new TCanvas();
  hQtot->Draw();
  
  TCanvas *c7 = new TCanvas();
  hQtotUsed->Draw();
  
  TCanvas *c8 = new TCanvas();
  hPadOccupancy->Draw("colz");
  
  TCanvas *c9 = new TCanvas();
  hPadOccupancyUsed->Draw("colz");
  
  TCanvas *c10 = new TCanvas();
  c10->SetLogy();
  hCherenkov->SetLineColor(kBlue);
  hCherenkov->Draw();
  
  TCanvas *c11 = new TCanvas();
  c11->SetLogy();
  hCherenkovUsedPions->SetLineColor(kBlue+2);
  hCherenkovUsedEle->SetLineColor(kGreen+2);
  hCherenkovUsedPions->Draw();
  hCherenkovUsedEle->Draw("same");
  
  TCanvas *c15 = new TCanvas();
  TF1 *pionfit = new TF1("pionfit","gaus",hdEdxPionTot->GetXaxis()->GetXmin(),hdEdxPionTot->GetXaxis()->GetXmax());
  TF1 *electronfit = new TF1("electronfit","gaus",hdEdxEleTot->GetXaxis()->GetXmin(),hdEdxEleTot->GetXaxis()->GetXmax());
  hdEdxEleTot->SetLineColor(kGreen+2);
  hdEdxPionTot->SetLineColor(kBlue+2);
  
  const Float_t frac=0.2;
  Int_t bin1=0,bin2=0;
  
  GetBinMinMax(hdEdxPionTot,frac,bin1,bin2); 
  hdEdxPionTot->Fit("pionfit","","",hdEdxPionTot->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionTot->GetXaxis()->GetBinUpEdge(bin2));
  GetBinMinMax(hdEdxEleTot,frac,bin1,bin2);
  hdEdxEleTot->Fit("electronfit","","",hdEdxEleTot->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleTot->GetXaxis()->GetBinUpEdge(bin2));
  
  //alternative fit
  //hdEdxPionTot->Fit("pionfit");
  //hdEdxEleTot->Fit("electronfit");
  
  
  hdEdxEleTot->GetFunction("electronfit")->SetLineColor(kGreen+2);
  hdEdxPionTot->GetFunction("pionfit")->SetLineColor(kBlue+2);
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
 
  
  float pionreserror = sqrt((pionsigmaerr/pionmeanTot)**2+((pionsigmaTot*pionmeanerr)/pionmeanTot**2)**2);
  float electronreserror = sqrt((electronsigmaerr/electronmeanTot)**2+((electronsigmaTot*electronmeanerr)/electronmeanTot**2)**2);
 
  float separationpowererr = sqrt((2*electronmeanerr/(electronsigmaTot+pionsigmaTot))**2+(2*pionmeanerr/(electronsigmaTot+pionsigmaTot))**2+
        ((2*electronsigmaerr*(pionmeanTot-electronmeanTot))/((electronsigmaTot+pionsigmaTot)**2))**2+((2*pionsigmaerr*(pionmeanTot-electronmeanTot))/((electronsigmaTot+pionsigmaTot)**2))**2);
  
  float elechisquareTot = electronfit->GetChisquare();
  float eleNDFTot = electronfit->GetNDF();
  float pionchisquareTot = pionfit->GetChisquare();
  float pionNDFTot = pionfit->GetNDF();

  
  
  TCanvas *c16 = new TCanvas();
  TF1 *pionfit = new TF1("pionfit","gaus",hdEdxPionMax->GetXaxis()->GetXmin(),hdEdxPionMax->GetXaxis()->GetXmax());
  TF1 *electronfit = new TF1("electronfit","gaus",hdEdxEleMax->GetXaxis()->GetXmin(),hdEdxEleMax->GetXaxis()->GetXmax());
  hdEdxEleMax->SetLineColor(kGreen+2);
  hdEdxPionMax->SetLineColor(kBlue+2);
  
  GetBinMinMax(hdEdxPionMax,frac,bin1,bin2); 
  hdEdxPionMax->Fit("pionfit","","",hdEdxPionMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionMax->GetXaxis()->GetBinUpEdge(bin2));
  GetBinMinMax(hdEdxEleMax,frac,bin1,bin2);
  hdEdxEleMax->Fit("electronfit","","",hdEdxEleMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleMax->GetXaxis()->GetBinUpEdge(bin2));
  
  
  
  hdEdxEleMax->GetFunction("electronfit")->SetLineColor(kGreen+2);
  hdEdxPionMax->GetFunction("pionfit")->SetLineColor(kBlue+2);
  hdEdxPionMax->Draw();
  hdEdxEleMax->Draw("same");
  float pionmeanMax = pionfit->GetParameter(1);
  float pionsigmaMax = pionfit->GetParameter(2);
  float electronmeanMax = electronfit->GetParameter(1);
  float electronsigmaMax = electronfit->GetParameter(2);
  
  float pionresMax = pionsigmaMax/pionmeanMax;
  float electronresMax = electronsigmaMax/electronmeanMax;
  float pionsigmaerr = pionfit->GetParError(2);
  float pionmeanerr = pionfit->GetParError(1);
  float electronsigmaerr = electronfit->GetParError(2);
  float electronmeanerr = electronfit->GetParError(1);
  
  
  float separationpowerMax = 2*(electronmeanMax-pionmeanMax)/(pionsigmaMax+electronsigmaMax);
 
  
  float pionreserrorMax = sqrt((pionsigmaerr/pionmeanMax)**2+((pionsigmaMax*pionmeanerr)/pionmeanMax**2)**2);
  float electronreserrorMax = sqrt((electronsigmaerr/electronmeanMax)**2+((electronsigmaMax*electronmeanerr)/electronmeanMax**2)**2);
 
  float separationpowererrMax = sqrt((2*electronmeanerr/(electronsigmaMax+pionsigmaMax))**2+(2*pionmeanerr/(electronsigmaMax+pionsigmaMax))**2+
        ((2*electronsigmaerr*(pionmeanMax-electronmeanMax))/((electronsigmaMax+pionsigmaMax)**2))**2+((2*pionsigmaerr*(pionmeanMax-electronmeanMax))/((electronsigmaMax+pionsigmaMax)**2))**2);
  
  float elechisquareMax = electronfit->GetChisquare();
  float eleNDFMax = electronfit->GetNDF();
  float pionchisquareMax = pionfit->GetChisquare();
  float pionNDFMax = pionfit->GetNDF();
  
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
}


