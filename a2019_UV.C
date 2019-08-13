//change binning to 100
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
  ValueWithError A;
} Result;



Result* a2019_UV (std::string filename = "small16.root", std::string name = "data_2016",  double sval = 37.4, bool isU = true){
  using namespace RooFit;

  gROOT->ProcessLine(".x mattStyle.C");

  //new canvas for the fits. So the graph isn't covered after running it for different data files altogether.
  TCanvas* can = new TCanvas(name.c_str(),name.c_str(), 800 ,600);
  Result* myres = new Result();



  // make some strings for the categories
  std::vector<std::string> ucats = {"uplus","uminus"};
  std::vector<std::string> vcats = {"vplus","vminus"};


  // make the categories for the simultaneous fit
  RooCategory sample("sample","sample") ;
  if (isU == true){
    for (auto is: ucats  ){
      sample.defineType(is.c_str()) ;
    }
  }
  else {
    for (auto is: vcats  ){
      sample.defineType(is.c_str()) ;
    }
  }



  // open the file and get a new dataset object
  TFile* theFile = new TFile(filename.c_str());
  TTree* theTree = (TTree*)theFile->Get("DecayTree");


  // set branches 
  float mass; theTree->SetBranchAddress("scaledmass",&mass);
  double U; theTree->SetBranchAddress("U",&U);
  double V; theTree->SetBranchAddress("V",&V);


  // set up the dataset -- >empty
  double minmass = 5250; double maxmass = 5500;
  RooRealVar* m = new RooRealVar("m","m", 5366, minmass,maxmass);
  RooDataSet data("DS", "DS", RooArgSet(*m,sample));


  // now loop and file the dataset
  int nentries = theTree->GetEntries();

  for (auto i = 0; i < nentries; ++i){
    theTree->GetEntry(i);
    if (mass > minmass && mass < maxmass){
      m->setVal(mass);
      if (isU == true){
        if (U >0) {
          sample.setLabel(ucats[0].c_str());  }
        else {
          sample.setLabel(ucats[1].c_str());  }
      }
      else{
        if (V >0) {
          sample.setLabel(vcats[0].c_str());  }
        else {
          sample.setLabel(vcats[1].c_str());  }
      }
    }
      data.add(RooArgSet(*m,sample));
  }



  //Trt to see what data looks like.
  //std::cout<<"data: "<<std::endl;
  //data.Print();//RooDataSet::DS[m,sample] = 13446 entries
  
 //5399.19 2 vplus //5413.34 3 vminus
  /*
  for (int i = 0; i< data.numEntries(); ++i){
    const RooArgSet* row = data.get(i);
    //    row->Print();
    double mval = row->getRealValue("m");
    int sval = row->getCatIndex("sample");
    const char* cval = row->getCatLabel("sample");
    std::cout << mval << " " << sval  << " " << cval << std::endl; 
  }
  */

  // make the simultaneous pdf
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample);
  simPdf.Print();//RooSimultaneous::simPdf[ indexCat=sample ] = 0

  




//############################################################################
//############################################################################


  //Construct two models for positive u and negative u.

  // 1st Gauss
  RooRealVar mean("mean","Mean of Gaussian",5366,minmass,maxmass) ; 
  RooRealVar sigma("sigma","Width of Gaussian",20,0,100) ; 
  RooGaussian gauss("gauss","gauss(x,mean,sigma)",*m,mean,sigma) ; //Construct a gaussian p.d.f gauss


  // 2nd Gauss
  RooRealVar sigma2("sigma2","Width of Gaussian",sval,0,200) ; 
  sigma2.setConstant(true);
  RooGaussian gauss2("gauss2","gauss2(x,mean,sigma)",*m,mean,sigma2) ; 


  // 3rd Exp
  RooRealVar lambda("lambda", "slope", 0, -5., 5);
  RooExponential expo("expo", "exponential PDF", *m, lambda);


  // 4rd Novosibirsk //not shared
  RooRealVar width("width","width of ?",78.98,0,200);//78.98//78.97
  width.setConstant(true);
  RooRealVar peak("peak","peak of signal?",5441,minmass,maxmass);//5441//5427
  peak.setConstant(true);
  RooRealVar tail("tail","tail of ?",0.94,-1,1);//0.94//0.57
  tail.setConstant(true);
  RooNovosibirsk Novo("Novo", "Novosibirsk", *m, peak,width,tail);


  //Composite them. //modelplus
  RooRealVar f("f","f",0.7, 0.,1);
  RooAddPdf gausgaus("gausgaus","gaus+gaus",RooArgList(gauss,gauss2),f);
  //adding two gaussian fits.
  

  RooRealVar n_exp("n_exp","nbackground",500,0. ,1.1*data.numEntries());
  RooExtendPdf e_exp("e_exp","ebackground",expo,n_exp) ;  

  RooRealVar n_Novo("n_Novo","nbackground",500, 0.,1.1*data.numEntries());
  RooExtendPdf e_Novo("e_Novo","ebackground",Novo,n_Novo) ;


  RooRealVar nsig("nsig","nsignal",1000,0.,1.1*data.numEntries());
  RooExtendPdf sigpdf("sigpdf","sigpdf",gausgaus,nsig);

  RooAddPdf model_plus("model_plus","gaus+gaus+exp+novo",RooArgList(sigpdf,e_exp,e_Novo));


  //model_minus

  RooRealVar n_exp_n("n_exp_n","nbackground",500,0. ,1.1*data.numEntries());
  RooExtendPdf e_exp_n("e_exp_n","ebackground",expo,n_exp_n) ;  

  RooRealVar n_Novo_n("n_Novo_n","nbackground",500, 0., 1.1*data.numEntries());//The name in the first "" is what shows in ternimal. Regardless of the variable name.
  RooExtendPdf e_Novo_n("e_Novo_n","ebackground",Novo,n_Novo_n) ;

  RooRealVar nsig_n("nsig_n","nsignal",1000,0.,1.1*data.numEntries());
  RooExtendPdf sigpdf_n("sigpdf_n","sigpdf",gausgaus,nsig_n);

  RooAddPdf model_minus("model_minus","gaus+gaus+exp+novo_n",RooArgList(sigpdf_n,e_exp_n,e_Novo_n));


  if (isU == true){
  // make the 2 pdfs and add then
  simPdf.addPdf(model_plus,"uplus");
  simPdf.addPdf(model_minus,"uminus");
  }
  else{
  simPdf.addPdf(model_plus,"vplus");
  simPdf.addPdf(model_minus,"vminus");
  }


  RooFitResult* result = simPdf.fitTo(data,Save());

  std::cout <<"\n"<< "############The likelihood and chi2 on the overall fitting###########" <<std::endl;
  std::cout << "Likelihood: " << result->minNll() << std::endl;  
  

  //plot model_plus
  RooPlot* frame_fit = m->frame() ; //-> instead of =
  
  frame_fit->GetXaxis()->SetTitle("m [MeV/c^{2}]");
  frame_fit->GetYaxis()->SetTitle("Candidates/(1.25 MeV/c^{2})");//1.5 5250-5550, 1.25 5500
  frame_fit->SetTitle("Fitting Scaledmass");
  frame_fit->GetXaxis()->SetTitleOffset(1.2);
  frame_fit->GetYaxis()->SetTitleOffset(1.2);
  //frame_fit->GetYaxis()->SetRangeUser(0, 250);//Set range of the y axis

  data.plotOn(frame_fit,Binning(100)) ;

  //1 black, 2 red, 3 green, 4 blue (seems default), 5 yellow
  //6 pink  //7 cyan //8 less shiney green //9 purple-ish
  //line and colour: https://root.cern.ch/doc/master/classTAttLine.html
  model_plus.plotOn(frame_fit,Name("Double Gaussian"),Components("gausgaus"),LineColor(2),LineStyle(2));
  model_plus.plotOn(frame_fit,Name("Exponential"),Components("e_exp"), LineColor(3),LineStyle(10));
  model_plus.plotOn(frame_fit,Name("Novosibirsk"),Components("e_Novo"), LineColor(5),LineStyle(9));

  model_plus.plotOn(frame_fit,Name("Overall fit")) ;
  //model_minus.plotOn(frame_fit,LineColor(2)) ;


  //chi2
  double chi2 = frame_fit->chiSquare();// This part has to be after model.//chi2 need it to be binned.
  double ndf = 100 - 10; // number of bins minus number of parameters?
  double probChi2 = TMath::Prob(ndf*chi2,ndf);
  std::cout << "prob " << probChi2 << " chi2/ndf " << chi2 << std::endl; //ndf: number of degree of freedom.
  std::cout << "\n" <<std::endl;


  frame_fit->Draw() ;


  //Legend has to be after xframe->Draw();
  auto legend = new TLegend(0.95,0.95,0.6,0.75);// In example 0.1,0.7,0.48,0.9//
  // Each legend entry is made of a reference to a ROOT object, a text label and an option specifying which graphical attributes (marker/line/fill) should be displayed.
  legend->SetBorderSize(1);
  //legend->AddEntry("data","Data","P");//"P" (point) or "LP" (point and line)
  if (isU == true){
  legend->AddEntry("Overall fit","Overall fit (for U>0)","l");
  }
  else{
  legend->AddEntry("Overall fit","Overall fit (for V>0)","l");
  }
  legend->AddEntry("Double Gaussian","Double Gaussian","l");
  if (isU == true){
  legend->AddEntry("Exponential","Exponential (for U>0)","l");
  legend->AddEntry("Novosibirsk","Novosibirsk (for U>0)","l");
  }
  else{
  legend->AddEntry("Exponential","Exponential (for V>0)","l");
  legend->AddEntry("Novosibirsk","Novosibirsk (for V>0)","l");
  }
  legend->Draw();
  


  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid = frame_fit->residHist() ;

  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull = frame_fit->pullHist() ;//Style: https://root.cern.ch/doc/master/classTAttFill.html
  hpull->SetLineColor(4);//4
  hpull->SetLineStyle(1);
  hpull->SetLineWidth(2);
  hpull->SetFillColor(92);//92
  hpull->SetFillStyle(3001);//3144  
  // On the old frame draw the residual distribution and add the distribution to the frame
  //RooPlot* frame2 = m->frame(Title("Residual Distribution")) ;
  //frame2->addPlotable(hresid,"P") ;
  //data.plotOn(frame2);

  //TCanvas* can2 = new TCanvas("Residual Distribution","Residual Destribution", 800 ,600); 
  //frame2->Draw("HISTO");


  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* frame3 = m->frame(Title("Pull Distribution")) ;
  frame3->addPlotable(hpull,"L3") ;
  frame3->GetXaxis()->SetTitle("m [MeV/c^{2}]");
  frame3->GetYaxis()->SetTitle("pull");

  std::string name2 = name+"pull";
  TCanvas* can3 = new TCanvas(name2.c_str(),name2.c_str(), 800 ,600); 
  frame3->Draw();


  std::cout << "###################################################################################" <<std::endl;

  // Combining errors by derivatives.  
  double dasym;
  double Np = nsig.getVal(); 
  double dNp = nsig.getError();
  double Nn=nsig_n.getVal();
  double dNn=nsig_n.getError();
  
  dasym=sqrt(pow(2*Nn/pow((Np+Nn),2),2)*pow(dNp,2)+pow(-2*Np/pow((Np+Nn),2),2)*pow(dNn,2));

  double asym_num=(Np-Nn)/(Np+Nn);




  std::cout << "\n" <<"###################### Error with Random Parameters Method############################" << "\n" <<std::endl;


//random parameters method

  RooFitResult* res = simPdf.fitTo(data,Save());

  TH1F* histo = new TH1F("histo","histo",1000, -0.2, 0.2); 
  histo->Sumw2();

  for (int i=0; i < 10000; i++){
    const RooArgList& sample = res->randomizePars();
    RooRealVar* fnsig1 = (RooRealVar*)sample.find("nsig");
    RooRealVar* fnsig2 = (RooRealVar*)sample.find("nsig_n");
    double uasym = (fnsig1->getVal() - 
fnsig2->getVal())/(fnsig2->getVal() + fnsig1->getVal());

    histo->Fill(uasym);
  }

  double val = histo->GetMean();
  double eval = histo->GetRMS();


  myres->A.value = val;
  myres->A.error = eval;


  std::cout <<"\n"<<"A = "<< val <<std::endl;
  std::cout <<"Error of A = "<< eval <<std::endl;

  //TCanvas* can2 = new TCanvas("ramdom parameters method","RPM", 800 ,600); 
  //histo->Draw("HISTO");

  std::cout <<"A (from propagating errors) = "<< asym_num <<std::endl;
  std::cout << "dA (from propagating errors) = "<< dasym <<std::endl;


  return myres;  


}
  
