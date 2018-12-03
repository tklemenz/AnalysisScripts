void getClusters(const char *filename){
	
  gSystem->Load("~/GEMtest/code/libGEMEvent.so");

  TFile *TreeFile = new TFile(Form("%s", filename));
  TTree *tree = (TTree*)TreeFile->Get("GEM");

  TH1D *hNclusters   = new TH1D("hNclusters", ";Number of clusters; Counts", 65,0,65);
  TH1D *hNIROCtracks = new TH1D("hNtracksIROC","; Number of tracks; # counts",10,0,10);
  TH1D *hdEdx        = new TH1D("hdEdx", "; d#it{E}/d#it{x} Q_{tot} [a.u.]; # counts", 125, 0, 250);
  TH1D *hQtot        = new TH1D("hQtot", "; Total cluster charge Q_{tot} [ADC counts]; # counts", 600,0,600);
  TH1D *hCherenkov   = new TH1D("hCherenkov", "; ADC signal (corresponds to velocity); # counts", 90,20,196);

  GEMEvent *fEvent=0x0;
  fEvent=0x0;
  GEM->SetBranchAddress("Events",&fEvent);

GEMTrack::SetRemoveRowsdEdx("9,17,25,32,39,45,51,57,58,59,60,61,62");

  for (int iev=0; iev<tree->GetEntries(); ++iev){
    tree->GetEntry(iev);
    hNIROCtracks->Fill(fEvent->GetNumberOfTracks(0));   // the testbeam data is mapped to the TPC sector 0
    for (int itrack=0; itrack<fEvent->GetNumberOfTracks(); ++itrack) {
	  GEMTrack *track=const_cast<GEMTrack*>(fEvent->GetTrack(itrack));
      if (track->GetROC()!= 0) continue;
      int ncl=track->GetNumberOfClusters();
      hNclusters->Fill(ncl);
      for (int icl=0; icl<ncl; ++icl){
        GEMCluster *cl = const_cast<GEMCluster*>(track->GetCluster(icl));
        hQtot->Fill(cl->GetQTot());
      }
      float dEdx = track->GetTruncatedMean(0.,.7,1,0);
      hdEdx->Fill(dEdx);
      
    }
    hCherenkov->Fill(fEvent->GetCherenkovValue());
    
  }


  TCanvas *c1 = new TCanvas();
  hNIROCtracks->Draw();

  TCanvas *c2 = new TCanvas();
  hNclusters->Draw();

  TCanvas *c3 = new TCanvas();
  hdEdx->Draw();

  TCanvas *c4 = new TCanvas();
  hQtot->Draw();

  TCanvas *c5 = new TCanvas();
  hCherenkov->Draw();
} 
