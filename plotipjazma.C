//=================================================================================
// author:  J.Nagle (jamie.nagle@colorado.edu)
// date:    July 22, 2018
// outline: Read in TTree of IP-Jazma calculation output and make some useful plots
//=================================================================================

void plotipjazma(
		 int  nuc1 = 2,
		 int  nuc2 = 197,
		 char filein[100] = "orange.root",
		 bool dilutedense = true,
		 bool qsfluc      = true,
		 double figuremaxx = 10.0,
		 bool figureonepanel = false
		 ) {

  int system = 0;
  if (nuc1 == 1 && nuc2 == 197) system = 1; // p+Au
  if (nuc1 == 2 && nuc2 == 197) system = 2; // d+Au

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // read in ROOT file with results of full IP-JAZMA calculation output
  cout << "Opening file = " << filein << endl;
  TFile *fin = new TFile(filein);
  if (!fin) {
    cout << "Error:: Input file " << filein << " not found!" << endl;
    return;
  }
  TTree *nt = static_cast <TTree *> (fin->Get("ntout"));
  if (!nt) {
    cout << "Error:  Input TTree (ntout) not found!" << endl;
    return;
  }

  // declare all leaves ?? of TTree 
  
  //=============================================================
  // DETAIL EACH CANVAS -- WHAT ARE WE PLOTTING
  //=============================================================  
  TCanvas *c = new TCanvas("c","c",10,10,900,450);
  if (!figureonepanel) c->Divide(2,1);

  // PLOT THE ENERGY DENSITY DISTRIBUTION WHICH IS THE SAME AS N_GLUON / <N_GLUON> HERE (SEE PAPER DETAILS)
  c->cd(1);
  TH1F *henergydensity = new TH1F("henergydensity","henergydensity",1000,0.0,1.0e5);
  nt->Project("henergydensity","out_energydensity","out_npart>0");
  double mean = henergydensity->GetMean();
  cout << "Initial mean energydensity = " << mean << endl;
  TH1F *henergydensityzoom = new TH1F("henergydensityzoom","henergydensityzoom",1000,0.0,10.0*mean);
  nt->Project("henergydensityzoom","out_energydensity","out_npart>0");
  mean = henergydensityzoom->GetMean();
  cout << "Zoomed mean energydensity = " << mean << endl;
  
  TH1F *henergydensityscaled = new TH1F("henergydensityscaled","henergydensityscaled",100,0.0,10.0);
  char footext[100];
  sprintf(footext,"out_energydensity/%f",mean);
  nt->Project("henergydensityscaled",footext,"out_npart>0");

  cout << "Fraction of events in lowest bin = " << henergydensityscaled->Integral(1,1)/henergydensityscaled->Integral() << endl;
  // remove the lowest bin as part of definition for minimum bias
  // note that since there is no hard cutoff for a single N-N interaction, this is always an open issue
  henergydensityscaled->SetBinContent(1,0.0);
  
  henergydensityscaled->Scale(1.0/henergydensityscaled->Integral());
  henergydensityscaled->SetXTitle("N_{gluon} / <N_{gluon}>  or dN_{ch}/d#eta / <dN_{ch}/d#eta>");
  henergydensityscaled->SetYTitle("Event Prob");
  henergydensityscaled->SetLineColor(kRed);
  henergydensityscaled->GetYaxis()->SetTitleOffset(1.2);
  henergydensityscaled->GetXaxis()->SetRangeUser(0.0,figuremaxx);
  henergydensityscaled->SetLineWidth(3);
  henergydensityscaled->DrawCopy();

  // where is the 5% highest-multiplicity line?
  double fivecut;
  TLine *tline;
  for (int i=henergydensityscaled->GetNbinsX();i>=1;i--) {
    if (henergydensityscaled->Integral(i,henergydensityscaled->GetNbinsX())/henergydensityscaled->Integral(1,henergydensityscaled->GetNbinsX())> 0.05) {
      fivecut = henergydensityscaled->GetBinCenter(i);
      cout << "5 percent cut at = " << fivecut << endl;
      tline = new TLine(fivecut,henergydensityscaled->GetMinimum(),fivecut,henergydensityscaled->GetMaximum());
      tline->SetLineStyle(2);
      tline->SetLineColor(kBlue);
      tline->SetLineWidth(2);      
      tline->Draw("same");
      break;
    }
  }

  char fooleg[100];
  if (system==1) sprintf(fooleg,"p+Au@200 GeV");
  if (system==2) sprintf(fooleg,"d+Au@200 GeV");
  if (system==3) sprintf(fooleg,"p+Pb@5 TeV");    
  TLegend *tleg = new TLegend(0.27,0.65,0.9,0.9,fooleg,"brNDC");
  if (dilutedense && !qsfluc) sprintf(fooleg,"IP-Jazma, dilute-dense, no Q_{s,0}^{2} fluc.");
  if (dilutedense && qsfluc) sprintf(fooleg,"IP-Jazma, dilute-dense, Q_{s,0}^{2} fluc.");
  if (!dilutedense  && !qsfluc) sprintf(fooleg,"IP-Jazma, dense-dense, no Q_{s,0}^{2} fluc.");    
  tleg->AddEntry(henergydensityscaled,fooleg,"l");

  // In the d+Au case, we can overlay with the MSTV calculation and the published STAR data
  // Requires input file!
  if (system == 2) {
    // overlay MSTV calculation and STAR data
    char filemstv[100] = "./fout_daumult_mstv.root";
    TFile *fmstv = new TFile(filemstv);
    if (!fmstv) {
      cout << "Error::  File " << filemstv << " not found!" << endl;
    } else {
      TH1F *hcgc = static_cast <TH1F *> (fmstv->Get("hfoo"));
      TH1F *hstar = static_cast <TH1F *> (fmstv->Get("hdata2"));
      hcgc->SetLineColor(kBlack);
      hcgc->SetLineWidth(1);
      if (!figureonepanel || !dilutedense) hcgc->Draw("same");
      hstar->SetMarkerStyle(20);
      if (!figureonepanel || !dilutedense) hstar->Draw("p,same");
      if (!figureonepanel || !dilutedense) {
	tleg->AddEntry(hstar,"STAR Data","p");
	tleg->AddEntry(hcgc,"MSTV Calc","l");
      }
    }
  }
  tleg->Draw("same");

  char foocut[100];

  if (!figureonepanel) {

    // IF MORE THAN ONE PANEL OPTION, DRAW THE MEAN QS0^2 FOR THE PROJECTILE NUCLEONS
    c->cd(2);

    TH1F *h = new TH1F("h","h",200,0.0,10.0);
    sprintf(footext,"out_qs2proj1");                                       // in p+A just the one nucleon
    if (system == 2) sprintf(footext,"(out_qs2proj1+out_qs2proj2)/2.0");   // in the d+A case, two projectile nucleons
    sprintf(foocut,"out_npart>0&&(out_energydensity/%3.2f)>%3.2f",mean,fivecut);
    nt->Project("h",footext,foocut);
    double meanqs21_withcut = h->GetMean();
    cout << "Mean Qs2 (with 0-5 cut) = " << meanqs21_withcut << endl;
    
    h->Reset();
    sprintf(footext,"out_qs2proj1");
    if (system == 2) sprintf(footext,"(out_qs2proj1+out_qs2proj2)/2.0");  
    sprintf(foocut,"out_npart>0&&(out_energydensity/%f)>%f",mean,0.0);
    nt->Project("h",footext,foocut);
    double meanqs21_nocut = h->GetMean();
    cout << "Mean Qs2 (with 0-100 cut) = " << meanqs21_nocut << endl;  
    
    TProfile *tqs21 = new TProfile("tqs21","tqs21",50,0.0,10.0,0.0,100.0);
    sprintf(footext,"out_qs2proj1:out_energydensity/%f",mean);
    if (system==2) sprintf(footext,"(out_qs2proj1+out_qs2proj2)/2.0:out_energydensity/%f",mean);  
    nt->Project("tqs21",footext,"out_npart>0");
    tqs21->SetMarkerStyle(20);
    tqs21->SetMarkerColor(kRed);
    tqs21->SetXTitle("N_{gluon} / <N_{gluon}> or dN_{ch}/d#eta / <dN_{ch}/d#eta>");
    tqs21->SetYTitle("<Q_{s,0}^{2}> for each projectile nucleon");
    tqs21->GetYaxis()->SetTitleOffset(1.2);
    tqs21->GetYaxis()->SetRangeUser(0.0,3.0);
    tqs21->GetXaxis()->SetRangeUser(0.0,figuremaxx);
    tqs21->Draw("p");
    
    if (system==1) sprintf(fooleg,"IP-Jazma p+Au@200 GeV");
    if (system==2) sprintf(fooleg,"IP-Jazma d+Au@200 GeV");
    TLegend *tleg2 = new TLegend(0.34,0.11,0.9,0.39,fooleg,"brNDC");
    sprintf(fooleg,"<Q_{s,0}^{2}> [0-100] = %3.2f",meanqs21_nocut);
    tleg2->AddEntry(tqs21,fooleg,"l");
    sprintf(fooleg,"<Q_{s,0}^{2}> [0-5] = %3.2f",meanqs21_withcut);
    tleg2->AddEntry(tqs21,fooleg,"l");
    tleg2->Draw("same");

  }
  //=========================================================================
  // END OF THE CANVAS ABOVE!
  //=========================================================================  

  //=========================================================================  
  // NOW PLOT THE AREA, THE QS2 AVERAGED OVER THE AREA
  //=========================================================================  
  
  TCanvas *crt = new TCanvas("crt","crt",10,10,900,450);
  crt->Divide(2,1);

  int tbins = 25;
  double maxrange = 10.0;
  // for comparing p+Au and d+Au, including a maxrange version with 
  bool fullrange = false;
  if (fullrange) {
    maxrange = 35.0;
    mean     = (125.0*15./12.)*2.0;
  }
  
  crt->cd(2);
  TProfile *tqs2 = new TProfile("tqs2","tqs2",tbins,0.0,maxrange,0.0,100.0);
  sprintf(foocut,"out_npart>0&&(out_energydensity/%3.2f)>%3.2f",mean,0.001);
  sprintf(footext,"out_qs2overarea:out_energydensity/%f",mean);    
  nt->Project("tqs2",footext,foocut);
  tqs2->GetYaxis()->SetTitleOffset(1.2);
  tqs2->SetMinimum(0.0);
  tqs2->SetXTitle("N_{gluon}/<N_{gluon}>");
  if (fullrange)     tqs2->SetXTitle("N_{gluon} (a.u.)");  
  tqs2->SetYTitle("< Q_{S}^{2} > over the Interaction Area");
  tqs2->SetMarkerStyle(20);
  tqs2->SetMarkerColor(kRed);
  tqs2->SetLineColor(kRed);  
  tqs2->GetXaxis()->SetRangeUser(0.0,figuremaxx);
  tqs2->Draw("p,l");

  TH1F *hqs2ave = new TH1F("hqs2ave","hqs2avearea",500,0.0,20.0);
  sprintf(foocut,"out_npart>0&&(out_energydensity/%3.2f)>%3.2f",mean,fivecut);  
  nt->Project("hqs2ave","out_qs2overarea",foocut);
  cout << "Average qs2proj within area (0-5) == " << hqs2ave->GetMean() << endl;

  if (system==1) sprintf(fooleg,"IP-Jazma p+Au@200 GeV");
  if (system==2) sprintf(fooleg,"IP-Jazma d+Au@200 GeV");    
  TLegend *tlegqs2ave = new TLegend(0.3,0.1,0.9,0.35,fooleg,"brNDC");
  sprintf(fooleg,"<Q_{S}^{2}(proj)> (0-5) = %3.2f",hqs2ave->GetMean());
  tlegqs2ave->AddEntry(tqs2,fooleg,"p");
  tlegqs2ave->Draw("same");
  
  crt->cd(1);
  TProfile *tarea = new TProfile("tarea","tarea",tbins,0.0,maxrange,0.0,100.0);
  sprintf(foocut,"out_npart>0&&(out_energydensity/%3.2f)>%3.2f",mean,0.001);
  sprintf(footext,"out_area:out_energydensity/%f",mean);    
  nt->Project("tarea",footext,foocut);
  tarea->SetMinimum(0.0);
  tarea->GetYaxis()->SetTitleOffset(1.2);
  tarea->SetXTitle("N_{gluon}/<N_{gluon}>");
  if (fullrange)     tarea->SetXTitle("N_{gluon} (a.u.)");  
  tarea->SetYTitle("Interaction Area [fm^{2}]");
  tarea->SetMarkerStyle(20);
  tarea->SetMarkerColor(kGreen+2);
  tarea->SetLineColor(kGreen+2);  
  tarea->GetXaxis()->SetRangeUser(0.0,figuremaxx);
  tarea->Draw("p,l");

  // area and average 
  TH1F *harea = new TH1F("harea","harea",500,0.0,20.0);
  sprintf(foocut,"out_npart>0&&(out_energydensity/%3.2f)>%3.2f",mean,fivecut);  
  nt->Project("harea","out_area",foocut);
  cout << "Average Area (0-5) == " << harea->GetMean() << endl;

  if (system==1) sprintf(fooleg,"IP-Jazma p+Au@200 GeV");
  if (system==2) sprintf(fooleg,"IP-Jazma d+Au@200 GeV");    
  TLegend *tlegarea = new TLegend(0.3,0.1,0.9,0.35,fooleg,"brNDC");
  sprintf(fooleg,"<Area> (0-5) = %3.2f",harea->GetMean());
  tlegarea->AddEntry(tarea,fooleg,"p");
  tlegarea->Draw("same");

  //=========================================================================
  // SEPARATE FIGURE OF DEUTERON RT DISTRIBUTION
  //=========================================================================  

  // only in the d+Au case
  if (system==2) {
    TCanvas *cdeuteronrt = new TCanvas("cdeuteronrt","cdeuteronrt",10,10,600,600);
    TProfile *trt = new TProfile("trt","trt",tbins,0.0,maxrange,0.0,100.0);
    sprintf(foocut,"out_npart>0&&(out_energydensity/%3.2f)>%3.2f",mean,0.001);
    sprintf(footext,"out_deuteronrt:out_energydensity/%f",mean);    
    nt->Project("trt",footext,foocut);
    trt->SetMinimum(0.0);
    trt->GetYaxis()->SetTitleOffset(1.2);
    trt->SetXTitle("N_{gluon}/<N_{gluon}>");
    trt->SetYTitle("< deuteron r_{T} >");
    trt->SetMarkerStyle(20);
    trt->SetMarkerColor(kBlue);
    trt->SetLineColor(kBlue);  
    trt->GetXaxis()->SetRangeUser(0.0,figuremaxx);
    trt->Draw("p,l");
  }
  
  //=========================================================================
  // SEPARATE PANEL WITH ECCENTRICITY COMPARISONS - FOR OTHER STUDIES
  //=========================================================================  
  

			       
}
