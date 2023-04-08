//#ifndef __CINT__
#include "RooGlobalFunc.h"
//#endif
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooBifurGauss.h"
#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooBinning.h"
#include "TH1.h"
#include "TH2.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooLandau.h"
#include "TChain.h"
#include<cmath>
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooSuperCategory.h"
#include "RooSimultaneous.h"
#include "RooNLLVar.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <string>
using namespace RooFit ;
using namespace std;
//int main(){

void mbcfit(){
  /*******************Fit Variables***********************************/
  RooRealVar mbc("mbc","M_{bc} (GeV)",5.23,5.29);
  /**defining DATAFRAME{unbinned histogram}(FROM FIT VARIABLE) TO FIT AND PLOT**/
  RooDataSet* data=new RooDataSet("data","data",RooArgSet(mbc));
  /*******************Input root file**********************************/
  TChain* chain=new TChain();
  chain->Add("/home/belle2/ssana/MC15ri/cs/test_10_8_thrustBm_cosTBTO/signal_scaled/test.root/tree");

  Double_t  de3, md03, mbc3, r23, kid3,pid3,sig,cont_prob;
  Int_t run;
  Int_t nevt3=(int)chain->GetEntries();

  // chain->SetBranchAddress("isSignal",&sig);
  chain->SetBranchAddress("deltaE",&de3);
  chain->SetBranchAddress("Mbc",&mbc3);
  chain->SetBranchAddress("D0_bar_InvM",&md03);
  chain->SetBranchAddress("SigProb",&cont_prob);
  chain->SetBranchAddress("Kp_PID_bin_kaon",&kid3);
  // chain->SetBranchAddress("R2",&r23);
  // chain->SetBranchAddress("pi_PID_bin_pion",&pid3);
  // chain->SetBranchAddress("__run__",&run);
  

  //Loading data 
  Double_t counter =0;
  for(int l=0;l<nevt3;l++) {
    chain->GetEntry(l);
    mbc.setVal(mbc3);
    if(md03>1.85 && md03<1.88 && mbc3>5.23 && mbc3 < 5.29 && de3 < 0.1 && de3 > -0.1 && kid3 > 0.6 && cont_prob > 0.3){
      data->add(RooArgSet(mbc));
      counter++;
    }
  }

  /*****************************Mbc Fit***********************/
  // --- Build Signal PDF ---
  RooRealVar mean1("mean1","mean of Gaussian-1",5.279145,5.27,5.29);
  RooRealVar mean2("mean2","mean of Gaussian-2",5.279260,5.27,5.29);
  RooRealVar sigma1("sigma1","sigma of Gaussian-1",0., 0., 1);	
  RooRealVar sigma2("sigma2","sigma of Gaussian-2",0.01191,0.00000001,1);

  RooGaussian sig1("sig1","Gaussian-1",mbc,mean1,sigma1);  
  RooGaussian sig2("sig2","Gaussian-2",mbc,mean2,sigma2);

  RooRealVar fsig_1("frac_gaussians", "signal fraction", 0.5,0.,1.);
  RooAddPdf twoGaussians("twoGaussians", "sum of two Gaussians ",RooArgList(sig1, sig2), RooArgList(fsig_1));

  // --- Build Argus background PDF ---
  RooRealVar argpar("argpar","argus shape parameter",-34.70,-100.,-1.) ;
  RooArgusBG bkg("argus","Argus PDF",mbc,RooConst(5.29),argpar) ; //mbc background

  //Initialization of parameter before adding two pdf
  // cout<<"Total number of events which are used to fitting are : "<<counter<<endl;
  Double_t event_count = counter; 
  Double_t signal_count = counter*0.2;
  Double_t back_count = counter*0.8;
  RooRealVar n_sig("n_sig", "n_sig", signal_count, 0., event_count);//52000
  RooRealVar n_bkg("n_bkg", "n_bkg", back_count, 0., event_count);//95000
  RooAddPdf sum("sum","sum",RooArgList(twoGaussians,bkg),RooArgList(n_sig, n_bkg));//adding two pdf
  sum.fitTo(*data);
  /****************************FIT COMPLETE*************************************/
  // // Integrate sig pdf(following some hypothesis) to get yield
  // mbc.setRange("twoGaussians",5.27, 5.29);     //twoGaussians is my signal pdf    //5.27, 5.29 is the range we want to integrate
  // RooAbsReal *integral_sig = twoGaussians.createIntegral(mbc,NormSet(mbc),Range("twoGaussians"));
  
  // double  Nsig = integral_sig->getVal();
  // cout<<"Signal :"<<Nsig*n_sig.getVal()<<endl;   // nsig is the signal Yield
  // double Nsigerr = n_sig.getError()*integral_sig->getVal();  
  // cout<<"Signal error = "<<Nsigerr<<endl;
  // cout<<"Signal Area "<< setprecision(4)<<Nsig<<endl;
  /*********************Start Plotting and showing outpouts*****************/
  //Plotting fitted result
  RooPlot* deframe = mbc.frame(Title("Fitting M_{bc} of B^{#pm}"), Bins(300)) ;                          
  data->plotOn(deframe, Binning(300), DataError(RooAbsData::SumW2)) ;
  // sum.plotOn(deframe, LineColor(kBlue)	, LineStyle(kSolid)) ;
  // sum.paramOn(deframe,data,"", 2, Format("NEU"),AutoPrecision(1)), 0.7, 0.9, 0.9); //"NELU",  Prints the fitted parameter on the canvas
  sum.plotOn(deframe,Components(sig1),LineColor(kGreen),LineStyle(kDashed)) ;
  sum.plotOn(deframe,Components(sig2),LineColor(kBlack),LineStyle(kDashed)) ;
  sum.plotOn(deframe,Components(twoGaussians),LineColor(kRed),LineStyle(kDashed));
  sum.plotOn(deframe,Components(bkg),LineColor(kMagenta),LineStyle(kDashed)) ;
  sum.plotOn(deframe, LineColor(kBlue), LineStyle(kSolid)) ; // we need to write it last.
                        // otherwise pull distribution will calculate error wrt last mentioned pdf inside plotOn function
  sum.paramOn(deframe,data,"", 2, "NEU", 0.1, 0.35, 0.9); //"NELU",  Prints the fitted parameter on the canvas

  //Extract info. from fitting frame and showing
  cout<<"chisq of the fit is : "<<deframe->chiSquare()<<endl;//chi-square of the fit
  cout<<"chi-square/ndof : "<<deframe->chiSquare(7)<<endl;// Chi^2/(the number of degrees of freedom)
  RooHist* hpull = deframe->pullHist() ;
  RooPlot* frame3 = mbc.frame(Title("Pull Distribution")) ;
  hpull->SetFillColor(1);
  frame3->addPlotable(hpull,"X0B"); // "X0" is for errorbar; and "B" is for histograms
  // frame3->addPlotable(hpull,"P");


  TCanvas* c1 = new TCanvas("c1", "c1", 2550, 1500) ;
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();  
  deframe->Draw() ;
  // Adding legend
  TLegend *legend1 = new TLegend(0.35,0.7,0.5,0.9);
  TLegendEntry *entry = legend1->AddEntry("sig1","1st Gaussian pdf","l");
  entry->SetLineColor(kGreen);
  entry->SetLineStyle(kDashed);
  entry = legend1->AddEntry("sig2","2nd Gaussian pdf","l");
  entry->SetLineColor(kBlack);
  entry->SetLineStyle(kDashed);
  entry = legend1->AddEntry("twoGaussians","Combined signal pdf","l");
  entry->SetLineColor(kRed);
  entry->SetLineStyle(kDashed);
  entry = legend1->AddEntry("bkg","bkg-ARGUS","l");
  entry->SetLineColor(kMagenta);
  entry->SetLineStyle(kDashed);
  entry = legend1->AddEntry("sum","Fitted pdf","l");
  entry->SetLineColor(kBlue);
  entry->SetLineStyle(kSolid);
  legend1->Draw();

  // Printing error in "Signal yield calculation" // Problem is resolved using extra "E" inside paramOn function 
  // std::string sig_err_in_str = std::to_string(n_sig.getError());
  // TLatex* sig_err = new TLatex();
  // sig_err->SetTextSize(0.034);
  // sig_err->SetTextAlign(12);  //centered aligned
  // sig_err->DrawLatex(5.2392, 1450, "#pm");
  // sig_err->DrawLatex(5.24, 1450, sig_err_in_str); // here we need to give "TEXT" ONLY

  // // Adding arrow at (+-)3*sigma of signal pdf
  // double weightedMean = mean1.getVal()*fsig_1.getVal() + mean2.getVal()*(1.0-fsig_1.getVal());
  // double weightedSigma = std::sqrt( pow(sigma1.getVal(),2)*pow(fsig_1.getVal(),2)  +  pow(sigma2.getVal(),2)*pow((1-fsig_1.getVal()),2));
  // double mbc_min_fit = weightedMean -3*weightedSigma; 
  // double mbc_max_fit = weightedMean +3*weightedSigma;
  // double arrowHeight = n_sig.getVal()*0.03;
  // cout<<"Mbc min :"<< mbc_min_fit <<endl;
  // cout<<"Mbc max :"<< mbc_max_fit <<endl;
  // cout<<"Arrow Height :"<< arrowHeight <<endl;

	// TArrow *arr1 = new TArrow(-0.01, 100, -0.01, 400, 0.01,"|>");
	// arr1->SetLineWidth(4);
  // arr1->SetLineColor(2);
  // arr1->SetFillStyle(3008);
  // arr1->DrawArrow(mbc_min_fit, arrowHeight*0.1, mbc_min_fit, arrowHeight, 0.02,"<"); 
  // arr1->DrawArrow(mbc_max_fit, arrowHeight*0.1, mbc_max_fit, arrowHeight, 0.02,"<");

  c1->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->Draw();
  pad2->cd(); 
  frame3->Draw() ;
  frame3->SetLineStyle(9);
  frame3->GetYaxis()->SetNdivisions(505);
  frame3->GetYaxis()->SetTitle("#sqrt{#chi^{2}}"); 
  frame3->GetXaxis()->SetTitle("M_{bc} (GeV)"); 
  frame3->GetXaxis()->SetTitleSize(0.13);
  frame3->GetYaxis()->SetTitleSize(0.15);
  frame3->GetXaxis()->SetLabelSize(0.120);
  frame3->GetYaxis()->SetLabelSize(0.120); 
  frame3->GetXaxis()->SetTitleOffset(0.90);      
  frame3->GetYaxis()->SetTitleOffset(0.22);       
  frame3->GetYaxis()->SetRangeUser(-10.0, 10.0);       
  frame3->GetYaxis()->SetLimits(-10.0, 10.0);       
  frame3->GetXaxis()->SetNdivisions(505);
  frame3->GetYaxis()->CenterTitle(true);
  frame3->GetXaxis()->CenterTitle(true);
  frame3->Draw("AXISSAME");

  cout<<"Total number of events which are used to fitting are : "<<counter<<endl;
  cout<<"chisq of the fit is :"<<deframe->chiSquare()<<endl;//chi-square of the fit
  cout<<"chi-square/ndof :"<<deframe->chiSquare(7)<<endl;// Chi^2/(the number of degrees of freedom)
  cout<<"Error in calculation of signal yield : "<<n_sig.getError()<<endl;
  
  c1->Print("mbc_plot/mbc_fit_all_cs_cut6.png");
}