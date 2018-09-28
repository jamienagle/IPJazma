//=================================================================================
// author:  J.Nagle (jamie.nagle@colorado.edu)
// date:    July 19, 2018
// outline: Read in TTree of standard format with MC Glauber nucleon information
//          Do full IP-Jazma calculation and output summary information in TNtuple
//=================================================================================

void ipjazma(int    nuc1 = 2,                       
	     int    nuc2 = 197,
	     char   filein[100]  = "lemon.root",
	     char   fileout[100] = "orange.root",
	     int    selectevent = 9999,
	     int    nbins = 800,
	     bool   qsfluc = true,
	     bool   alphaslocal = false,
	     bool   dilutedense = false,
	     double IPSAT_RWIDTH_OVERRIDE = 0    // =0 means use the default below
	     ) {
		 
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //==============================================================================
  // Set up key IP-Jazma calculation parameters and functions 
  //==============================================================================  

  // 2-D LATTICE IN X,Y SPACE
  //-------------------------
  double binwidth = 20.0/((double)nbins);   // default running with nbins=800 --> dx,dy=0.025 fm
  TH2F *hp = new TH2F("hp","hp",nbins,-10.,10.,nbins,-10.,10.);    // stores Q2 values summing all projectile nucleon contributions
  TH2F *ht = new TH2F("ht","ht",nbins,-10.,10.,nbins,-10.,10.);    // stores Q2 values summing all target nucleon contributions
  TH2F *he = new TH2F("he","he",nbins,-10.,10.,nbins,-10.,10.);    // stores IP-Jazma resulting energy density values
  
  // IP-SAT PARAMETERS
  //-------------------------
  double IPSAT_RWIDTH = 0.32;  // see Romatschke+Romatschke figure 4.5 (result for 200 GeV)
  // ACCORDING TO MSTV EMAIL #3, RMS OF COLOR CHARGE DIST IS 0.56, AND SO IPSAT_RWIDTH = 0.40 (TRY THAT AS ALTERNATIVE)
  IPSAT_RWIDTH = 0.40;
  if (nuc1==1 && nuc2 == 208) IPSAT_RWIDTH = 0.29;  // see Romatschke+Romatschke figure 4.5 (result for 5 TeV)
  
  if (IPSAT_RWIDTH_OVERRIDE > 0) IPSAT_RIDWTH = IPSAT_RWIDTH_OVERRIDE;

  double IPSAT_Q2CENTER = 0.67;  // arbitrarily set for now...
  if (nuc1==1 && nuc2 == 208) IPSAT_Q2CENTER = 1.54;
  
  TF1 *ipsat_rdist = new TF1("ipsat_rdist","[0]*TMath::Exp(-x*x/(2.0*[1]*[1]))",0.0,10.0);  
  ipsat_rdist->SetParameters(1.0,IPSAT_RWIDTH); // this is just the shape, normalization is done below with Q2WEIGHT
  
  // number of bins for going 3.0 sigma out
  // only follow IP-Sat gaussian out to 3 sigma in this case, can change for different tests
  double rmax = 3.0 * IPSAT_RWIDTH;
  int nbinrange = 1 + (int) (rmax/binwidth);
  
  // for the option of including running alpha_s, read in parameterization data here
  TGraph *alphagraph = new TGraph("alphas.dat","%lg, %lg");

  // Qs FLUCTUATIONS - PARAMETERIZATION FROM MSTV PAPER AND THEIR VARIANCE
  //------------------------------------------------------------------------------
  TF1 *fqsfluc = new TF1("fqsfluc","(2.0/x)*TMath::Exp(-TMath::Power(TMath::Log(x*x),2)/(2.*0.5*0.5))",0.0,6.0);
  // a GetRandom returns a value for QS, and then one needs to square it for QS^2 weighting!
  
  // DILUTE-DENSE LIMIT FUNCTION F(x)
  //----------------------------------
  TH1D *dilutedensefunction;
  if (dilutedense) {
    TFile *fdilutedense = new TFile("dilutedensefunction.root");
    dilutedensefunction = static_cast <TH1D *> (fdilutedense->Get("F"));
  }

  //==============================================================================  
  // ROOT FILE INPUT with MC Glauber information
  //==============================================================================    
  TFile *fin = new TFile(filein);
  if (!fin) {
    cout << "Error:: Input file " << filein << " not found!" << endl;
    return;
  }
  TTree *lemon = static_cast <TTree *> (fin->Get("lemon"));
  if (!lemon) {
    cout << "Error:  Input TTree (lemon) not found!" << endl;
    return;
  }

  //Declaration of leaves types on the lemon tree
  Int_t           npart;
  Int_t           ncoll;
  int             nparta;
  int             npartb;
  float           b;                    // collision impact parameter
  float           eccgaus[10];
  float           eccpoint[10];  
  Int_t           nproj;
  Int_t           ntarg;
  Float_t         xproj[400];
  Float_t         yproj[400];
  Float_t         xtarg[400];
  Float_t         ytarg[400];
  
  // Set branch addresses.
  lemon->SetBranchAddress("npart",&npart);
  lemon->SetBranchAddress("nparta",&nparta);
  lemon->SetBranchAddress("npartb",&npartb);  
  lemon->SetBranchAddress("ncoll",&ncoll);
  lemon->SetBranchAddress("b",&b);
  lemon->SetBranchAddress("eccgaus",eccgaus);
  lemon->SetBranchAddress("eccpoint",eccpoint);
  lemon->SetBranchAddress("nproj",&nproj);
  lemon->SetBranchAddress("ntarg",&ntarg);
  lemon->SetBranchAddress("xproj",xproj);
  lemon->SetBranchAddress("yproj",yproj);
  lemon->SetBranchAddress("xtarg",xtarg);
  lemon->SetBranchAddress("ytarg",ytarg);

  //==============================================================================
  // output summary file for making figures
  //==============================================================================

  TFile *fout = new TFile(fileout,"RECREATE");

  TTree *ntout = new TTree("ntout","ntout");
  // variable list for output
  // global event information passed from MCGlauber input TTree (npart, ncoll, nparta, npartb, ecc)
  int       out_npart;
  int       out_ncoll;
  int       out_nparta;
  int       out_npartb;
  float     out_b;                    // collision impact parameter
  float     out_eccgaus[10];
  float     out_eccpoint[10];  
  int       out_nproj;
  int       out_ntarg;
  ntout->Branch("out_npart",&out_npart,"out_npart/I");
  ntout->Branch("out_nparta",&out_nparta,"out_nparta/I");
  ntout->Branch("out_npartb",&out_npartb,"out_npartb/I");  
  ntout->Branch("out_ncoll",&out_ncoll,"out_ncoll/I");
  ntout->Branch("out_b",&out_b,"out_b/F");
  ntout->Branch("out_eccgaus",out_eccgaus,"out_eccgaus[10]/F");
  ntout->Branch("out_eccpoint",out_eccpoint,"out_eccpoint[10]/F");  
  ntout->Branch("out_nproj",&out_nproj,"out_nproj/I");
  ntout->Branch("out_ntarg",&out_ntarg,"out_ntarg/I");
  // then add new information from IP-Jazma characterizations
  float out_energydensity;
  float out_area;
  float out_qs2overarea;
  float out_deuteronrt;
  float out_qs2proj1;
  float out_qs2proj2;
  float out_ratio_qs2projovertarg;
  float out_eccjazma[10];
  ntout->Branch("out_energydensity",&out_energydensity,"out_energydensity/F");
  ntout->Branch("out_area",&out_area,"out_area/F");
  ntout->Branch("out_qs2overarea",&out_qs2overarea,"out_qs2overarea/F");
  ntout->Branch("out_deuteronrt",&out_deuteronrt,"out_deuteronrt/F");
  ntout->Branch("out_qs2proj1",&out_qs2proj1,"out_qs2proj1/F");
  ntout->Branch("out_qs2proj2",&out_qs2proj2,"out_qs2proj2/F");
  ntout->Branch("out_ratio_qs2projovertarg",&out_ratio_qs2projovertarg,"out_qs2projovertarg/F"); 
  ntout->Branch("out_eccjazma",out_eccjazma,"out_eccjazma[10]/F");  

  //==============================================================================
  // LOOP OVER EVENTS AND CALCULATE!
  //==============================================================================  
  
  int    eventindex =  0;
  
  double QSsqr1;
  double QSsqr2;
  double ecc2gaus;
  double ecc2point;
  
  // LOOP OVER EVENTS FROM MC GLAUBER INPUT TTREE
  Long64_t nentries = lemon->GetEntries();
  Long64_t nbytes = 0;
  for (Long64_t i=0; i<nentries;i++) {

    nbytes += lemon->GetEntry(i);

    // fill the lattice with Q2 values for the projectile and target in the hp and ht histograms
    for (int inuc=0;inuc<nuc1;inuc++) {
      double Q2weight = IPSAT_Q2CENTER; 
      if (qsfluc) Q2weight = IPSAT_Q2CENTER * TMath::Power(fqsfluc->GetRandom(),2.0);  // function has mean near one, so okay...
      if (inuc == 0) QSsqr1 = Q2weight; // special storage for the deuteron case !!!!
      if (inuc == 1) QSsqr2 = Q2weight;
      
      // loop over just a window around the center, not the entire grid
      int thisbinx = hp->GetXaxis()->FindBin(xproj[inuc]);
      int thisbiny = hp->GetYaxis()->FindBin(yproj[inuc]);
      int lowbinx  = thisbinx - nbinrange;
      if (lowbinx < 1) lowbinx = 1;
      int lowbiny  = thisbiny - nbinrange;
      if (lowbiny < 1) lowbiny = 1;
      int highbinx = thisbinx + nbinrange;
      if (highbinx > nbins) highbinx = nbins;
      int highbiny = thisbiny + nbinrange;
      if (highbiny > nbins) highbiny = nbins;
      
      for (int ix=lowbinx;ix<=highbinx;ix++) {
	for (int iy=lowbiny;iy<=highbiny;iy++) {
	  double radialdist = TMath::Sqrt(TMath::Power(xproj[inuc]-hp->GetXaxis()->GetBinCenter(ix),2.0)+
					  TMath::Power(yproj[inuc]-hp->GetYaxis()->GetBinCenter(iy),2.0));
	  if (radialdist > rmax) continue;
	  double GaussianWeight = ipsat_rdist->Eval(radialdist);
	  hp->Fill(hp->GetXaxis()->GetBinCenter(ix),hp->GetYaxis()->GetBinCenter(iy), Q2weight * GaussianWeight);
	}
      }
    }
    for (int inuc=0;inuc<nuc2;inuc++) {
      double Q2weight = IPSAT_Q2CENTER;
      if (qsfluc) Q2weight = IPSAT_Q2CENTER * TMath::Power(fqsfluc->GetRandom(),2.0);   
      
      // loop over just a window around the center, not the entire grid
      int thisbinx = ht->GetXaxis()->FindBin(xtarg[inuc]);
      int thisbiny = ht->GetYaxis()->FindBin(ytarg[inuc]);
      int lowbinx  = thisbinx - nbinrange;
      if (lowbinx < 1) lowbinx = 1;
      int lowbiny  = thisbiny - nbinrange;
      if (lowbiny < 1) lowbiny = 1;
      int highbinx = thisbinx + nbinrange;
      if (highbinx > nbins) highbinx = nbins;
      int highbiny = thisbiny + nbinrange;
      if (highbiny > nbins) highbiny = nbins;
      
      for (int ix=lowbinx;ix<=highbinx;ix++) {
	for (int iy=lowbiny;iy<=highbiny;iy++) {
	  double radialdist = TMath::Sqrt(TMath::Power(xtarg[inuc]-ht->GetXaxis()->GetBinCenter(ix),2.0)+
					  TMath::Power(ytarg[inuc]-ht->GetYaxis()->GetBinCenter(iy),2.0));
	  if (radialdist > rmax) continue;	      
	  double GaussianWeight = ipsat_rdist->Eval(radialdist);
	  ht->Fill(ht->GetXaxis()->GetBinCenter(ix),ht->GetYaxis()->GetBinCenter(iy), Q2weight * GaussianWeight);
	}
      }
    }

    // calculate the energy density in each cell
    double energydensitysum = 0.0;
    double meanx            = 0.0;
    double meany            = 0.0;
    
    double ratio = 0.0;
    double rationorm = 0.0;
      
    for (int ix=1;ix<=hp->GetNbinsX();ix++) {
      for (int iy=1;iy<=hp->GetNbinsY();iy++) {
	
	double proj = hp->GetBinContent(ix,iy);
	double targ = ht->GetBinContent(ix,iy);
	double energy = proj*targ;              // dense-dense limit calculation
	
	if (dilutedense) {
	  double projectilefactor = proj;
	  masscut = 0.3;
	  // note the sqrt (Qs^2) because the function arg. is Qs/m
	  // use Cyrille's function (dilutedensefunction)
	  double targetfactor = 0.0;
	  if (targ > 0) 
	    targetfactor = dilutedensefunction->GetBinContent( dilutedensefunction->FindBin( TMath::Sqrt(targ)/masscut ) );
	  energy = projectilefactor * targetfactor;
	}
	  
	if (alphaslocal) {
	  double QSproj = TMath::Sqrt(proj);
	  double QStarg = TMath::Sqrt(targ);
	  double QSscale = QSproj;
	  if (QStarg > QSproj) QSscale = QStarg;
	  double alphasvalue = alphagraph->Eval(QSscale);
	  energy = alphasvalue * energy;
	}
	
	he->SetBinContent(ix,iy,energy);
	  
	energydensitysum += energy;

	// calculate the center-of-mass as we go along to avoid another loop
	meanx += he->GetXaxis()->GetBinCenter(ix) * energy;
	meany += he->GetYaxis()->GetBinCenter(iy) * energy;

	// Qs(proj)/Qs(targ) weighted by the energy ...
	if (targ > 0) {
	  ratio     += (proj / targ) * energy;
	  rationorm += energy;
	}
	
      }
    }   // end loop over full grid

    ratio = ratio / rationorm;
    
    meanx = meanx / energydensitysum;
    meany = meany / energydensitysum;
      
    // calculate eccentricity and save output
    // find the meanx, meany
    // then calculate mean r^2 cos(2phi), mean r^2 sin(2phi) and divide out r^2 normalization
    // phi = TMath::ATan2(y,x)
    // eccN[n] = TMath::Sqrt(sinphi[n]*sinphi[n]+cosphi[n]*cosphi[n])/rn[n];
    double sinsum[10];
    double cossum[10];
    double normsum[10];
    for (int j=0;j<10;j++) {
      sinsum[j]  = 0.0;
      cossum[j]  = 0.0;
      normsum[j] = 0.0;
    }
    
    for (int ix=1;ix<=hp->GetNbinsX();ix++) {
      for (int iy=1;iy<=hp->GetNbinsY();iy++) {
	double xtemp = he->GetXaxis()->GetBinCenter(ix) - meanx;
	double ytemp = he->GetYaxis()->GetBinCenter(iy) - meany;
	double phitemp = TMath::ATan2(ytemp,xtemp);
	double weight = he->GetBinContent(ix,iy);

	// only calculate epsilon moments 2-6 and put them in those index values
	for (int j=2;j<=6;j++) {
	  // epsilon_n weighted by r^n
	  double rtemp = TMath::Power( TMath::Power(xtemp,2) + TMath::Power(ytemp,2) , ((float)j)/2.0 );
	  sinsum[j] += weight * rtemp * TMath::Sin(((float)j)*phitemp);
	  cossum[j] += weight * rtemp * TMath::Cos(((float)j)*phitemp);
	  normsum[j] += weight * rtemp;
	}
      }
    }
    // clean up this code to use arrays!
    float eccjazma[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int j=2;j<=6;j++) {
      eccjazma[j] = TMath::Sqrt(sinsum[j]*sinsum[j] + cossum[j]*cossum[j])/normsum[j];
    }
    
    // also calculate the <Qs^2> projectile and area with some cutoff criteria
    double area = 0.0;
    double averageqs2 = 0.0;
    for (int ix=1;ix<=hp->GetNbinsX();ix++) {
      for (int iy=1;iy<=hp->GetNbinsY();iy++) {
	double qs2proj = hp->GetBinContent(ix,iy);
	double gluons  = he->GetBinContent(ix,iy);
	if (gluons > 0.0001) { // some overlap here
	  double AREAQ2CUT = 0.1;
	  if (qs2proj > AREAQ2CUT) { // some arbitrary cutoff
	    area += 1.0; // normalize at the end
	    averageqs2 += qs2proj;
	  }
	}
      }
    }
    
    averageqs2 = averageqs2 / area;
    area = area * binwidth * binwidth;  // convert from unit count to fm^2 for each little box
    
    double deuteronrt = TMath::Sqrt(TMath::Power(xproj[0]-xproj[1],2.0)+TMath::Power(yproj[0]-yproj[1],2.0));

    //============================================================================================================
    // save information out to TTree
    // first copy over information from input MCGlauber TTree for later analysis (except all nucleon positions)
    //============================================================================================================    

    out_npart = npart;
    out_ncoll = ncoll;
    out_nparta = nparta;
    out_npartb = npartb;
    out_b = b;    
    for (int j=0;j<10;j++) out_eccgaus[j] = eccgaus[j];
    for (int j=0;j<10;j++) out_eccpoint[j] = eccpoint[j];    
    out_nproj = nproj;
    out_ntarg = ntarg;

    // now fill new ipjazma calculation output variables
    
    out_energydensity = energydensitysum;
    out_area = area;
    out_qs2overarea = averageqs2;
    out_qs2proj1 = QSsqr1;
    out_qs2proj2 = QSsqr2;
    out_deuteronrt = deuteronrt;
    out_ratio_qs2projovertarg = ratio;
    for (int j=0;j<10;j++) out_eccjazma[j] = eccjazma[j];

    ntout->Fill();
    
    //=========================================================================
    // THIS CODE IS JUST FOR DISPLAYING A SINGLE EVENT AND THEN BREAKING...
    //=========================================================================
    if (eventindex == selectevent && !(selectevent > 999)) {
      // draw things and then break!
      TCanvas *cdisplay = new TCanvas("cdisplay","cdisplay",10,10,1200,450);
      cdisplay->Divide(3,1);
      cdisplay->cd(1);
      hp->SetXTitle("x coordinate [fm]");
      hp->SetYTitle("y coordinate [fm]");
      hp->Draw("colz");
      TLatex *t1 = new TLatex(-9.5,8.5,"Q_{s}^{2}(x,y) Projectile");
      t1->Draw("same");
      
      cdisplay->cd(2);
      ht->SetXTitle("x coordinate [fm]");
      ht->SetYTitle("y coordinate [fm]");
      ht->Draw("colz");
      TLatex *t2 = new TLatex(-9.5,8.5,"Q_{s}^{2}(x,y) Target");
      t2->Draw("same");	
      
      cdisplay->cd(3);
      he->SetXTitle("x coordinate [fm]");
      he->SetYTitle("y coordinate [fm]");
      
      he->Draw("colz");
      char foolatex[100];
      //      sprintf(foolatex,"#varepsilon(x,y) #propto Q_{s}^{2}(x,y)^{proj} #times Q_{s}^{2}(x,y)^{targ}");
      //      if (dilutedense) sprintf(foolatex,"#varepsilon(x,y) #propto Q_{s}^{2}(x,y)^{proj} #times F(Q_{s}(x,y)^{targ}/m)");
      sprintf(foolatex,"IP-Jazma #varepsilon(x,y)");
      TLatex *t3 = new TLatex(-9.5,8.5,foolatex);
      t3->Draw("same");
      
      TGraph *tmean = new TGraph();
      tmean->SetPoint(0,meanx,meany);
      tmean->SetMarkerStyle(29);

      cout << "Eccentricity (IP-Jazma)  = " << eccjazma[2] << endl;
      cout << "Eccentricity (MCG point) = " << eccpoint[2] << endl;
      cout << "Eccentricity (MCG gauss) = " << eccgaus[2] << endl;
      break;
      
    } // end of single event display section
    
    ht->Reset();
    hp->Reset();
    he->Reset();

    eventindex++;  // increment the event counter
    if (eventindex % 100 == 0) cout << "processed event = " << eventindex << endl;

    if (eventindex > selectevent) break; 
    
  } // END LOOP OVER EVENTS!

  //=======================================
  // WRITE OUT THE IPJAZMA SUMMARY FILE
  //=======================================
 
  ntout->Write("ntout");
  fout->Close();
  
}

// END OF PROGRAM - THAT'S ALL FOLKS!
