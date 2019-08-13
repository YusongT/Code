#include "RooNovosibirsk.h"
#include "TTree.h"
#include "TFile.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooExtendPdf.h"
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include <iostream>
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TMath.h"
#include "TAxis.h"
#include <cmath>
#include <math.h> 

//https://root-forum.cern.ch/t/error-invalid-use-of-incomplete-type/11778 //to solve the error "allocation of incomplete type 'TH1F'"
#include "TH1F.h"
#ifndef ROOT_TH1F
#include "TH1F.h"
#endif

#include <TCanvas.h>

//like a class
typedef struct {
  double value;
  double error;
} ValueWithError;//ValueWithError is the class name

typedef struct {
  ValueWithError V;
} Result;



void a201902_generate (){
  using namespace RooFit;

  double minmass=5250; double maxmass=5500;
  RooRealVar* m = new RooRealVar("m","m", 5366,minmass,maxmass);//5366?

  //all parameters are the result of fitting
  // 1st Gauss
  RooRealVar mean("mean","Mean of Gaussian",5367,minmass,maxmass) ; //2015, 2016: 5368,5367
  RooRealVar sigma("sigma","Width of Gaussian",14.62,0,100) ; //14.72,13.05//13.00, 14.62
  RooGaussian gauss("gauss","gauss(x,mean,sigma)",*m,mean,sigma) ; //Construct a gaussian p.d.f gauss. gauss and gauss(x,mean,sigma) are names.

  // 2nd Gauss
  RooRealVar sigma2("sigma2","Width of Gaussian",37.4,0,200) ; 
  sigma2.setConstant(true);
  RooGaussian gauss2("gauss2","gauss2(x,mean,sigma)",*m,mean,sigma2) ; 


  // 3rd Exp
  RooRealVar lambda("lambda", "slope", 0., -5., 5);//-0.01, -0.0008//-0.01,0
  RooExponential expo("expo", "exponential PDF", *m, lambda);


  // 4rd Novosibirsk //not shared
  RooRealVar width("width","width of ?",79.98,0,200);
  // Don't have to set constant anymore.
  //width.setConstant(true);
  RooRealVar peak("peak","peak of signal?",5441,minmass,maxmass) ; 
  //peak.setConstant(true);
  RooRealVar tail("tail","tail of ?",0.94,-1,1) ;
  //tail.setConstant(true);
  RooNovosibirsk novo("novo", "Novosibirsk", *m, peak,width,tail);




  //Composite the two gaussians
  RooRealVar f("f","f",0.85, 0.,1);//0.74, 0.85//0.75, 0.85
  RooAddPdf gausgaus("gausgaus","gaus+gaus",RooArgList(gauss,gauss2),f);

  int nevent=6723;//1062, 6723

  RooRealVar nexp("nexp","nbackground",1531,0.,nevent);//97, 1511//102,1531
  //nexp.setVal(0);  nbkg1_exp.setConstant(true);
  RooExtendPdf eexp("eexp","ebackground",expo,nexp) ;  

  RooRealVar nNovo("nNovo","nbackground",134,0.,nevent);//127,126//137,134
  //nnovo.setVal(0);  nnovo.setConstant(true);
  RooExtendPdf eNovo("eNovo","ebackground",novo,nNovo) ;


  RooRealVar nsig("nsig","nsignal",5058,0.,nevent);//2 because this is building the extended pdf of gausgaus.//838,5086//823,5058
  RooExtendPdf sigpdf("sigpdf","sigpdf",gausgaus,nsig);

  RooAddPdf model("model","gaus+gaus+exp+novop",RooArgList(sigpdf,eexp,eNovo));//use same model for V> and <0



  TRandom ran; //used later on: int nevent=ran.Poisson(1000); Is this defining a new name "ran" for "TRandom"?
  //ran.SetSeed(99999);//to really be random.

  //creat files
  for (int fileindex = 0; fileindex < 100; ++fileindex){
    
    std::string filename = "toy_" + std::to_string(fileindex) + ".root"; 
    TFile* file = new TFile(filename.c_str(),"RECREATE" );
    TTree* tree = new TTree("DecayTree", "RECREATE"); 

    //int nevent = ran.Poisson(2000);//generate a number "nevent" around 1000 that fits Poisson distribution?
    //std::cout<<nevent<<"\n"<<std::endl;//output: 1023 983 1013 973 1023 when ifile=5//these numbers are the same every time...//still the same...
    //2033, 1931, 2030, 1956, 1998
    
    //RooDataSet* dset = gausgaus.generate(*m, nevent); //generate the m for n events that fits gausgaus pdf?


    RooDataSet* dset=model.generate(*m,Extended(kTRUE));//p.23 in manual, the number of events is predicted by the p.d.f. Extended(kTRUE) introduce a Poisson fluctuation in the number of generated events.

    float mass; TBranch* branch_m = tree->Branch("scaledmass",&mass,"scaledmass/F"); 
    double V; TBranch* branch_v = tree->Branch("V",&V,"V/D");


   for (int i = 0 ; i < dset->numEntries();++i) {
      const RooArgSet* row = dset->get(i);
      mass = (float)row->getRealValue("m");  
      if (ran.Uniform() > 0.5){//.Uniform() runs from 0 to 1. 
	V = -1;//set V.
        }
      else {
        V = 1;
        } 
      tree->Fill();
    }


    tree->Write();
    file->Close();
 }


}

  
