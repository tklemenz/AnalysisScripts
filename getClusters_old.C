
//void getClusters(const char *filename);
  
/*void ApplyGainCorrection(GEMTrack *track, GEMEvent *ev);


__________________________________________________________________________

void ApplyGainCorrection(GEMTrack *track, GEMEvent *ev){
  
  
  
  if (!fGainMap) return;
  
  for (Int_t icl=0; icl<track->GetNumberOfClusters(); ++icl){
    GEMCluster *cl = const_cast<GEMCluster*>(track->GetCluster(icl));
    const Int_t row=cl->GetRow();
    const Int_t pad=TMath::Nint(cl->GetPad());
    Double_t gainMapCorr=1.;
    if (fGainMap){
      gainMapCorr=fGainMap->GetValue(row,pad);
    }
    
    cl->CorrectCharge(gainMapCorr);
  }
  
}
*/

//__________________________________________________________________________

void getClusters(const char *filename){
  gSystem->Load("~/GEMtest/code/libGEMEvent.so");

 

  TFile *TreeFile = new TFile(Form("%s", filename));
  TTree *tree = (TTree*)TreeFile->Get("GEM");
  
  TFile *fGain = new TFile("~/ROC_Raw_data/GainMap2014.root");
  fGainMap = (AliTPCCalROC*)fGain->Get("gainMap");
  
  



  TH1D *hNclusters   = new TH1D("hNclusters", ";Number of clusters; Counts", 65,0,65);
  TH1D *hNIROCtracks = new TH1D("hNtracksIROC","; Number of tracks; # counts",10,0,10);
  TH1D *hdEdxEle     = new TH1D("hdEdxEle", "; d#it{E}/d#it{x} Q_{tot} [a.u.]; # counts", 125, 0, 250);
  TH1D *hdEdxPion    = new TH1D("hdEdx", "; d#it{E}/d#it{x} Q_{tot} [a.u.]; # counts", 125, 0, 250);
  TH1D *hQtot        = new TH1D("hQtot", "; Total cluster charge Q_{tot} [ADC counts]; # counts", 600,0,600);
  TH1D *hCherenkov   = new TH1D("hCherenkov", "; ADC signal (corresponds to velocity); # counts", 90,20,196);
  TH2D *PadResponse  = new TH2D("hResp", "; Row; Pad", 63,0,63,35,-5,30);
  //TH1D *hQtotSinglePad     = new TH1D("QtotSinglePad", "; Q_{tot}; # counts", 100,0,600);
  TGraph *Chisquare  = new TGraph();
  TH1D *hChisquare   = new TH1D("hChisquare", "; Chisquare Landau; # counts", 24,0,12);
  TH2D *hGainmap     = new TH2D("hGainmap", "; Row; Pad", 63,0,63,35,-5,30);
  //TH1D *hResidualPad = new TH1D("hResPad", "; Residual; #Counts",100,0,100);
  TH2D *hQtotTimeBin = new TH2D("hQTime","; TimeBin; Q_{tot}",100,0,100,200,0,600);
  
  
  
  const int nrows = 63;
  const int npads = 30;
  TH1D *hQtotPads[nrows][npads];

for(int irow=0; irow<nrows; ++irow) {
	for (int ipad=0; ipad<npads; ++ipad) {
	hQtotPads[irow][ipad] = new TH1D(Form("hQtotPad_%i_%i", irow, ipad), "", 600, 0, 600);
	}
}

  GEMEvent *fEvent=0x0;
  fEvent=0x0;
  tree->SetBranchAddress("Events",&fEvent);

GEMTrack::SetRemoveRowsdEdx("9,17,25,32,39,45,51,57,58,59,60,61,62");
  
  
  Double_t gainMapCorr=1.;
  for (int iev=0; iev<tree->GetEntries(); ++iev){										// loop over all events
    tree->GetEntry(iev);
    /*if (iev == 0){
	  GEMTrack *track=const_cast<GEMTrack*>(fEvent->GetTrack(0));
	  int ncl=track->GetNumberOfClusters();
	  for (int icl=0; icl<ncl; ++icl){
	    GEMCluster *cl = const_cast<GEMCluster*>(track->GetCluster(icl));
	    cout<<double(cl->GetRow())<<", "<<cl->GetPad()<<", "<<cl->GetTimeBinW()<<", "<<cl->GetPadRMS()<<", "<<cl->GetTimeBinWRMS()<<", "<<cl->GetTimeMax()<<", "<<cl->GetCPad()<<", "<<cl->GetMaxPad()<<", "
	    <<double(cl->GetNPads())<<endl;
	  } 
	  cout<<"Number of clusters: "<<ncl<<endl;
	}*/ 
    int NTracks = fEvent->GetNumberOfTracks(0);			// the testbeam data is mapped to the TPC sector 0 
    
	if (NTracks == 1) {																	// cut to 1 track events
	  GEMTrack *track=const_cast<GEMTrack*>(fEvent->GetTrack(0));
	  if (track->GetROC()!= 0) continue;												// take tracks only from ROC 0
	    
	  int ncl=track->GetNumberOfClusters();
	  
	  
	  if (ncl>30){																		// cut away tracks with less than 31 clusters
	  int CherenkovValue = fEvent->GetCherenkovValue();
        
        
        
        if (CherenkovValue < 40){														// pions
     	  hCherenkov->Fill(CherenkovValue);
	      hNIROCtracks->Fill(NTracks);
	      hNclusters->Fill(ncl);
	      for (int icl=0; icl<ncl; ++icl){
            GEMCluster *cl = const_cast<GEMCluster*>(track->GetCluster(icl));
            const Int_t row=cl->GetRow();
            const Int_t pad=TMath::Nint(cl->GetPad());
            Double_t gainMapCorr=1.;
            if (fGainMap){
              gainMapCorr=fGainMap->GetValue(row,pad);
              
            }
            cl->CorrectCharge(gainMapCorr);
            
		  }
	      float dEdxPion = track->GetTruncatedMean(0.,.7,1,0);
          hdEdxPion->Fill(dEdxPion);
          for (int icl=0; icl<ncl; ++icl){
            GEMCluster *cl = const_cast<GEMCluster*>(track->GetCluster(icl));
            hQtot->Fill(cl->GetQTot());
            PadResponse->Fill(double(cl->GetRow()),int(cl->GetCPad()));
            //hQtotTimeBin->Fill(cl->GetTimeBinW(),cl->GetQTot());
            if (int(cl->GetCPad())>4 && int(cl->GetCPad())<18){
              if (cl->GetTimeBinW() > 24 && cl->GetTimeBinW() < 91){
                hQtotPads[double(cl->GetRow())][int(cl->GetCPad())]->Fill(cl->GetQTot()/dEdxPion);
                hQtotTimeBin->Fill(cl->GetTimeBinW(),cl->GetQTot());
              //hResidualPad->Fill(cl->GetResidualPad());
		      }
            }
          }
	      
          
          
        }
        else if (CherenkovValue > 50){													// electrons
		  hCherenkov->Fill(CherenkovValue);
	      hNIROCtracks->Fill(NTracks);
	      hNclusters->Fill(ncl);
	      for (int icl=0; icl<ncl; ++icl){
            GEMCluster *cl = const_cast<GEMCluster*>(track->GetCluster(icl));
            const Int_t row=cl->GetRow();
            const Int_t pad=TMath::Nint(cl->GetPad());
            Double_t gainMapCorr=1.;
            if (fGainMap){
              gainMapCorr=fGainMap->GetValue(row,pad);
              
            }
            cl->CorrectCharge(gainMapCorr);
		  }
	      float dEdxEle = track->GetTruncatedMean(0.,.7,1,0);
          hdEdxEle->Fill(dEdxEle);
          for (int icl=0; icl<ncl; ++icl){
            GEMCluster *cl = const_cast<GEMCluster*>(track->GetCluster(icl));
            hQtot->Fill(cl->GetQTot());
            PadResponse->Fill(double(cl->GetRow()),int(cl->GetCPad()));
            //hQtotTimeBin->Fill(cl->GetTimeBinW(),cl->GetQTot());
            if (int(cl->GetCPad())>4 && int(cl->GetCPad())<18){
              if (cl->GetTimeBinW() > 24 && cl->GetTimeBinW() < 91){
              hQtotPads[double(cl->GetRow())][int(cl->GetCPad())]->Fill(cl->GetQTot()/dEdxEle);
		      hQtotTimeBin->Fill(cl->GetTimeBinW(),cl->GetQTot());
              //hResidualPad->Fill(cl->GetResidualPad());
		      }
            }
            /*double row = cl->GetRow();
            int pad = cl->GetCPad();
            if (row == 1 && pad == 0){
			  hQtotSinglePad->Fill(cl->GetQTot());	
			}*/	
            //cout<<cl->GetRow()<<", "<<cl->GetCPad()<<endl;
          }
	      
          
        } 
	      
	  }
	    
    }     
          
  } 
 /* TF1 *MPVfit = new TF1("MPV","landau",0,600); 
  int counter = 0;
  for (int irow=0; irow<nrows; ++irow){
    for (int ipad=0; ipad<npads; ++ipad){
	  if (hQtotPads[irow][ipad]->GetMean() != 0){	
        hQtotPads[irow][ipad]->Fit("MPV");
        //float chisquare = MPV->GetChisquare()/MPV->GetNDF();
        //Chisquare->SetPoint(Chisquare->GetN(),Chisquare->GetN(),chisquare);
        //hChisquare->Fill(chisquare);
        hGainmap->Fill(irow,ipad,MPV->GetParameter(1));
        //cout<<"Row: "<<irow<<", Pad: "<<ipad<<", MPV: "<<MPV->GetParameter(1)<<", ChiSquare: "<<chisquare<<endl;
        counter = counter + 1;
      }
      
    }
  }    
  
  hGainmap->GetZaxis()->SetRangeUser(0,1.15);

*/


  TCanvas *c1 = new TCanvas();
  hNIROCtracks->Draw();

  TCanvas *c2 = new TCanvas();
  hNclusters->Draw();

  TCanvas *c3 = new TCanvas();
  TF1 *pionfit = new TF1("pionfit","gaus",0,250);
  TF1 *electronfit = new TF1("electronfit","gaus",0,250);
  hdEdxEle->SetLineColor(kGreen);
  hdEdxPion->Fit("pionfit","mrg");
  hdEdxEle->Fit("electronfit","mrg");
  hdEdxEle->GetFunction("electronfit")->SetLineColor(kGreen+2);
  hdEdxPion->GetFunction("pionfit")->SetLineColor(kBlue+2);
  hdEdxPion->Draw();
  hdEdxEle->Draw("same");
  float pionmean = pionfit->GetParameter(1);
  float pionsigma = pionfit->GetParameter(2);
  float electronmean = electronfit->GetParameter(1);
  float electronsigma = electronfit->GetParameter(2);
  
  float pionres = pionsigma/pionmean;
  float electronres = electronsigma/electronmean;
  float pionsigmaerr = pionfit->GetParError(2);
  float pionmeanerr = pionfit->GetParError(1);
  float electronsigmaerr = electronfit->GetParError(2);
  float electronmeanerr = electronfit->GetParError(1);
  
  
  float separationpower = 2*(electronmean-pionmean)/(pionsigma+electronsigma);
 
  
  float pionreserror = sqrt((pionsigmaerr/pionmean)**2+((pionsigma*pionmeanerr)/pionmean**2)**2);
  float electronreserror = sqrt((electronsigmaerr/electronmean)**2+((electronsigma*electronmeanerr)/electronmean**2)**2);
 
  float separationpowererr = sqrt((2*electronmeanerr/(electronsigma+pionsigma))**2+(2*pionmeanerr/(electronsigma+pionsigma))**2+
        ((2*electronsigmaerr*(pionmean-electronmean))/((electronsigma+pionsigma)**2))**2+((2*pionsigmaerr*(pionmean-electronmean))/((electronsigma+pionsigma)**2))**2);
  
  cout<<"Pion dE/dx resolution: ("<<pionres*100<<" +- "<<pionreserror*100<<") %"<<endl;
  cout<<"Electron dE/dx resolution: ("<<electronres*100<<" +- "<<electronreserror*100<<") %"<<endl;
  cout<<"Separation power: ("<<separationpower<<" +- "<<separationpowererr<<") sigma"<<endl;
  

  TCanvas *c5 = new TCanvas();
  hQtot->Draw();

  TCanvas *c6 = new TCanvas();
  hCherenkov->Draw();
  
  TCanvas *c7 = new TCanvas();
  PadResponse->Draw("colz");
  
  //TCanvas *c8 = new TCanvas();
  //Chisquare->Draw("ap");
  
  //TCanvas *c9 = new TCanvas();
  //hChisquare->Draw();
  
  //TCanvas *c10 = new TCanvas();
  //hGainmap->Draw("colz");
  
  //TCanvas *c11 = new TCanvas();
  //hQtotSinglePad->Draw();
  
  //TCanvas *c12 = new TCanvas();
  //hResidualPad->Draw();
  
  TCanvas *c13 = new TCanvas();
  hQtotTimeBin->Draw("colz");
  
  //cout<<counter<<endl;
} 


//if(timeBin < 25 || timeBin > 90) isOK = false;
//if(cl_Pad < 3 || cl_Pad > 17) ++nclEdge; 

