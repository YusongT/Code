//Draw a histogram instead of TGraph, in order to label the x axis with years.
void a2019_runfitsHisto(){

  gROOT->ProcessLine(".x mattStyle.C");// to modify the style of graph 


  bool isU = true;
  std::vector<std::string> runs = {"small11.root","small12.root","small15.root","small16.root", "small17.root"};//add more data

  std::vector<std::string> names = {"2011","2012","2015", "2016", "2017"};//name of canvas and label for x axis.
  //std::vector<double> sval = {38.2,37.6,40.2,40.2,40.2};// and add more here as well.
  std::vector<double> sval = {34.4,35.3,37.4,37.4,37.4};


  TH1D* graph = new TH1D("h","h",5, 1, 6 );//histogram instead of TGraph. TH1D : histograms with one double per channel. Maximum precision 14 digits

  //std::vector<std::string> Asym_results;

  for (int i = 0; i < runs.size(); ++i){
    Result* myres = a2019_UV(runs[i],names[i],sval[i],isU);
    graph->SetBinContent(i+1,myres->A.value);
    graph->SetBinError(i+1, myres->A.error);
    graph->GetXaxis()->SetBinLabel(i+1,names[i].c_str());
    //Asym_results.push_back(myres->U.value);
    //Asym_results.push_back(myres->U.error);
  }
  

  
  TCanvas* can = new TCanvas("can","can", 800., 600);

  //graph->SetMarkerStyle(20);
  graph->GetXaxis()->SetTitle("Year");
  graph->GetXaxis()->SetLabelSize(.07);
  graph->GetYaxis()->SetTitle("Asymmetry");
  graph->SetTitle("Asymmetry");
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->GetYaxis()->SetTitleOffset(1.2);
  //graph->GetYaxis()->SetLabelFont(132);
  //graph->SetTitleFont(132); //error: no member named 'SetTitleFont' in 'TGraphErrors'
  
  

  //TCanvas* can = new TCanvas("can","can", 800., 600);

  //print the fitted constant
  TFitResultPtr r= graph->Fit("pol0","S");//"S" is save
  Double_t par0   = r->Value(0); // retrieve the value for the parameter 0
  Double_t err0   = r->ParError(0); // retrieve the error for the parameter 0
  std::cout<<"par0: "<<par0<<"+/-"<<err0<<std::endl;

  //Plot the fitted constant.
  //graph->Fit("pol0");
  TF1* fun = graph->GetFunction("pol0");
  fun->SetLineColor(2);
  //graph->Draw("AP");// draw the graph.
  //fun->Draw("SAME");


  //graph->Fit("pol0");
  //TF1* fun = graph->GetFunction("pol0");
  //fun->SetLineColor(2);//make the fitted line colour red
  //graph->Draw();// draw the graph.
  //  graph->Print();
  //fun->Draw("SAME");

  //In this case can't do Afits->Print(). So have to print manually:
  //TH1.Print Name  = Afits, Entries= 5, Total sum= 0.040345
  for (int i = 1; i <= 5; ++i){
    std::cout << names[i-1] <<  " " << graph->GetBinContent(i) <<  " +/- " << graph->GetBinError(i) << std::endl; 
  }




  /*
  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid = graph->residHist() ;

  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull = graph->pullHist() ;
  hpull->SetLineColor(4);//4
  hpull->SetLineStyle(1);
  hpull->SetLineWidth(2);
  hpull->SetFillColor(92);//92
  hpull->SetFillStyle(3001);//3144


  // On the old frame draw the residual distribution and add the distribution to the frame
  RooPlot* frame2 = m->frame(Title("Residual Distribution")) ;
  frame2->addPlotable(hresid,"P") ;
  //data.plotOn(frame2);

  TCanvas* can2 = new TCanvas("Residual Distribution","Residual Destribution", 800 ,600); 
  frame2->Draw("HISTO");


  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* frame3 = m->frame(Title("Pull Distribution")) ;
  frame3->addPlotable(hpull,"L3") ;
   
  TCanvas* can3 = new TCanvas("Pull Distribution","Pull Distribution", 800 ,600); 
  frame3->Draw();
  */


  //Save the graph of the asymmetry in ofile.root.
  TFile* ofile_trial = new TFile("Fit_Result.root","RECREATE");
  graph->SetName("Afits");
  graph->Write();
  ofile_trial->Close();

  //print final result
  
  //Ufits->Print()
  
  //for (int i=0;i<Asym_results.numEntries();i++){
  //  std::cout<<Asym_results.at(i)<<std::endl;
  //}
  
  
}
