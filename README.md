##  IP-Jazma

Hello and welcome to the IP-Jazma code documentation provided by Jamie Nagle and Bill Zajc.   This code was used for
the calculations presented in the IP-Jazma first paper https://arxiv.org/abs/1808.01276

The steps are relatively straightforward, but of course feel free to contact me (jamie.nagle@colorado.edu) 
with questions and suggestions.

The code is divided so that one can run any Monte Carlo Glauber (MCG) code one wants.    That MCG needs to output
a ROOT TTree called "lemon" into a standard ROOT file we call "lemon.root".   The format of this TTree lemon is defined below.    Note that the more important inputs are the x,y coordinates of all the projectile and target nucleons -- and it is very important that they are in the same coordinate frame (i.e. not in the center-of-mass of each nucleus, and 
so including the impact parameter offset).    It is interesting to compare IP-Jazma results event-by-event with
MCG eccentricities and global MCG characteristics like number of participating nucleons (npart) and
 number of binary collisions (ncoll) and impact parameter (b) and so we save these out as well.

```c

  const int maxNucleons = 400;
  int       npart;
  int       ncoll;
  int       nparta;
  int       npartb;
  float     b;                    // collision impact parameter
  float     eccgaus[10];
  float     eccpoint[10];  
  int       nproj;
  int       ntarg;
  float     xproj[maxNucleons];  // x,y coordinates for all nucleons
  float     yproj[maxNucleons];  // note these must be in the global coordinate frame
  float     xtarg[maxNucleons];  // x,y coordinates for all nucleons
  float     ytarg[maxNucleons];  // note these must be in the global coordinate frame

  // create new TTree with MC Glauber information for input to IP-Jazma
  TTree *lemon = new TTree("lemon","lemon");
  lemon->Branch("npart",&npart,"npart/I");
  lemon->Branch("nparta",&nparta,"nparta/I");
  lemon->Branch("npartb",&npartb,"npartb/I");  
  lemon->Branch("ncoll",&ncoll,"ncoll/I");
  lemon->Branch("b",&b,"b/F");
  lemon->Branch("eccgaus",eccgaus,"eccgaus[10]/F");
  lemon->Branch("eccpoint",eccpoint,"eccpoint[10]/F");  
  lemon->Branch("nproj",&nproj,"nproj/I");
  lemon->Branch("ntarg",&ntarg,"ntarg/I");  
  lemon->Branch("xproj",xproj,"xproj[400]/F");
  lemon->Branch("yproj",yproj,"yproj[400]/F");
  lemon->Branch("xtarg",xtarg,"xtarg[400]/F");
  lemon->Branch("ytarg",ytarg,"ytarg[400]/F");
```

As an aside, I use a modified version of the publicly available PHOBOS MC Glauber (http://tglaubermc.hepforge.org) 
to generate my TTree lemon.  I am working on having this standard in the next release of the MCG code - should be there by
the end of October 2018 (soon).   In the meanwhile, you can use some of my lemon.root files available here:

http://www.phenix.bnl.gov/WWW/publish/nagle/IPJAZMA/LemonFiles/

Also, feel free to ask me to generate others as useful.

Then the real work is done in "ipjazma.C".    Simply open a ROOT session, compile the code, and run.   Note that
some versions of ROOT do not like the conversion of a char[] to a string and give a harmless warning.   When you run ipjazma()
you can set the input parameters and flags:

```c

    int nuc1 = 2   // number of projectile nucleons
    int nuc2 =197  // number of target nucleons
    char filein[100] = "lemon.root"  // input MCG TTree file with TTree "lemon"
    char fileout[100] = "orange.root" // output file for making nice plots and inline analysis
    selectevent = 9999   // select the number of events from lemon.root to analyze.
                         // If you select a number < 9999, it will display that event
    int nbins = 800      // number of bins in x,y lattice grid (800 x 800 corresponds to 0.025 fm x 0.025 fm)
    bool qsfluc = true   // boolean to select whether to include Qs,0^2 fluctuations for each nucleon or not
    bool alphaslocal = false // select whether the g^2 (i.e. alpha_s encoded) is fixed in the energy density formula
    bool dilutedense = false // this selects whether one uses dense-dense or dilute-dense equation for the energy density
    double IPSAT_RWIDTH_OVERRIDE = 0 // 0 means the default, otherwise pick you favorite value

```

Output of running this code is a ROOT file "orange.root".

The code can be improved in terms of speed, but for this version the priority was clarity.  Running over 10k 
d+Au events with a grid of 800 x 800 takes about 2 hours on a single core.

You can make some standard plots by running the macro "plotipjazma.C".

The code is quite reasonably documented so it should be straightforward to make other plots as you wish.   Again, I am
happy to help answer any questions.



