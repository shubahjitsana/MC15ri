//Make the splot
//This is a simple example showing how to make an sPlot using RooStats and RooFit.
//The method is described in the paper paper by Pvik and De Liberder, 
// Nucl. Instrum. Meth. A555 356 (2005) [arXiv:physics/0402083]
//Here the data consists of signal events and background events. The method describes how to separate signal events from background.
//So two variables are used , one is discriminating variable and another one is the control varible. Here discriminating variable is jpsi mass and the control variable is energy/momentum. At first, We generate a two component data set from two varibles stored in the ntuple. Then fit the jpsi mass with two pdfs. Then calculate the signal and background weights to plot E/p splots. Two splots are found with signal contributions and backgrounds.

//Author: Debashis Sahoo & Karim Trabelsi


//Make the splot

void splot(TString fname="allmc_stream0.root") {
using namespace RooFit;
using namespace RooStats;
gSystem->Load("libRooFit");
gSystem->Load("libRooStats.so");
RooRealVar mass("mass", "GeV",2.9,3.3);
RooRealVar rep("rep", "E/p",0,1.5);
RooDataSet* data = new RooDataSet("data", "data",RooArgSet(mass,rep));

// fit variable
Float_t jpsimass,clmx_en,mum_ptot,ratio;

// Read data file
TFile* input=new TFile(fname);
TTree* t1 = (TTree*) input->Get("h1");

// no of entries in the datafile

Int_t n_tot = (int)t1->GetEntries();
cout << " n_tot = " << n_tot << endl;
t1->SetBranchAddress("jpsimass",&jpsimass);
t1->SetBranchAddress("clmx_en",&clmx_en);
t1->SetBranchAddress("mum_ptot",&mum_ptot);
// Import the data points
  for(int i=0; i<n_tot; i++) {
      t1->GetEntry(i); 
     ratio = clmx_en/mum_ptot;
  if(jpsimass >=2.9 && jpsimass <=3.3) {   
     mass.setVal(jpsimass); 
     rep.setVal(ratio);
     data->add(RooArgSet(mass,rep));       
     }                                                     
  }

  
 cout << "make mass model" << endl;
 //make mass model for signal(jpsi) .............add jp infront of each var
 RooRealVar jpMmean("jpM#mu_{x}","jpMmean",3.09,3.05,3.125);
 RooRealVar jpMsigma1("jpM#sigma_{1}","jpmssigma1",0.006,0.0,0.01);
 RooGaussian jpMsignal1("jpMsignal1","jpsignal1",mass,jpMmean,jpMsigma1);
 RooRealVar jpMsigma2("jpM#sigma_{2}","jpsigma2",0.03,0.0,0.08);
 RooGaussian jpMsignal2("jpMsignal2","jpsignal2",mass,jpMmean,jpMsigma2);
 RooRealVar jpMsfrac("jpMsfrac","Area fraction",0.5,0.0 ,1.0);
 RooAddPdf jpMsignald("jpMsignald","jpMsignald",RooArgList(jpMsignal1,jpMsignal2),jpMsfrac);

//make mass model for background ............  
 RooRealVar bkMa("bkMa","bkMa",-2.0,-4.0,0);
 RooExponential bkMbkg("bkMbkg","bkMexp",mass,bkMa);

 RooRealVar nsig("nsig", "nsig",100000,50000,900000);
 RooRealVar nbkg("nbkg", "nbkg",100000,50000,1300000);
// RooRealVar nsig("nsig", "nsig",700000,100000,10000000);
// RooRealVar nbkg("nbkg", "nbkg",9000000,500000,10300000);


 RooAddPdf massPdf("massPdf", "", RooArgList(jpMsignald, bkMbkg), RooArgList(nsig, nbkg));

  // Fit the data for mass in order to make splots of muid

  RooFitResult * r = (RooFitResult*)massPdf.fitTo(*data," ", RooFit::Save());
  r->Print("v");

  
  // Create the sPlot data from the fit to mass.
  RooStats::SPlot* sDataX = new RooStats::SPlot("sData","An SPlot", *data, &massPdf, RooArgList(nsig, nbkg) );
  cout << "Yield of signal is " << nsig.getVal() << " From sWeights it is " << sDataX->GetYieldFromSWeight("nsig") <<endl;
  cout << "Yield of background is " << nbkg.getVal() << "  From sWeights it is " << sDataX->GetYieldFromSWeight("nbkg") <<endl;

 // cout << "signal Weight   " << sDataX->GetSWeight("nsig") 
 //      << " bkg Weight   " <<  sDataX->GetSWeight("nbkg") <<endl;
   
  for(Int_t i=0; i < 20; i++)
    {
      cout << "signal Weight   " << sDataX->GetSWeight(i,"nsig") 
		<< " bkg Weight   " << sDataX->GetSWeight(i,"nbkg") 
		<< "  Total Weight   " << sDataX->GetSumOfEventSWeight(i) 
		<< endl;
    }
  
  // create weighted data sets
  RooDataSet * dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"nsig_sw") ;
  RooDataSet * dataw_bg  = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"nbkg_sw") ;

  cout << "Making splots of the signal and background" <<endl;
  RooPlot* frame_sig_rep = rep.frame(100) ;
  frame_sig_rep->SetTitle("sPlot for the E/p distribution(data)");
  dataw_sig->plotOn(frame_sig_rep, RooFit::DataError(RooAbsData::SumW2) ) ;

  RooPlot* frame_bg_rep = rep.frame(100) ;
  frame_bg_rep->SetTitle("sPlot for the background E/p distribtuion(data)");
  dataw_bg->plotOn(frame_bg_rep, RooFit::DataError(RooAbsData::SumW2) ) ;

  // plot PDFs onto the sPlots; if the sPlot correctly reweights the data then the signal / background only shapes should
  // provide an adequate description of the data.


  TCanvas *c1 = new TCanvas("c1" ,"c1");
  c1->cd();
  frame_sig_rep->Draw();
 
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  frame_bg_rep->Draw();
  
 // c1.Print("splot_pt.pdf");
 c1->SaveAs("mc_splot_rdedx_sig.root");
 c2->SaveAs("mc_splot_rdedx_bkg.root");

 //  RooDataHist *dh1 = new RooDataHist("dh1"," dh1",RooArgSet(muid),*dataw_sig);
   TH1F *h1 = (TH1F*)dataw_sig->createHistogram("signal" ,rep);
   double scale = 1/h1->Integral();
   h1->Scale(scale);
   h1->SetLineColor(kRed);
   h1->SetLineWidth(5);

 //  RooDataHist *dh1_b = new RooDataHist("dh1_b"," dh1_b",RooArgSet(muid),*dataw_bg);
   TH1F *h1_b = (TH1F*)dataw_bg->createHistogram("background" ,rep);
   double scaleb = 1/h1_b->Integral();
   h1_b->Scale(scaleb);
   h1_b->SetLineColor(kBlue);
   h1_b->SetLineWidth(5);

  TCanvas *c3 = new TCanvas("c3" ,"c3");
  c3->cd();
  h1->Draw();
  h1_b->Draw("same");
  auto legend = new TLegend( 0.6,0.7,0.9,0.85);
  legend->AddEntry(h1,"Signal","l");
  legend->AddEntry(h1_b,"Background","l");
  legend->Draw();
   c3->SaveAs("splot_sig_bkg.root");


  }
















