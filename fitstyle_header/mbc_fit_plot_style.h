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
using namespace RooFit;
using namespace std;
void mbc_fit_plot_style(RooRealVar* var_name, RooDataSet* dataset_name, 
  RooRealVar* mean1, RooRealVar* mean2, RooRealVar* sigma1, RooRealVar* sigma2,
  RooGaussian signal_pdf_1, RooGaussian signal_pdf_2, RooRealVar* fsig_1,  RooAddPdf total_signal_pdf,
  RooRealVar* argpar, RooArgusBG backround_pdf, RooRealVar* n_sig, RooRealVar* n_bkg,
  RooAddPdf* added_pdf, const char* output_name)
{
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
  RooPlot* deframe = var_name->frame(Title("Fitting M_{bc} of B^{#pm}"), Bins(300)) ;                          
  dataset_name->plotOn(deframe, Binning(300), DataError(RooAbsData::SumW2)) ;
  // sum.plotOn(deframe, LineColor(kBlue)	, LineStyle(kSolid)) ;
  // sum.paramOn(deframe,data,"", 2, Format("NEU"),AutoPrecision(1)), 0.7, 0.9, 0.9); //"NELU",  Prints the fitted parameter on the canvas
  added_pdf->plotOn(deframe,Components(signal_pdf_1),LineColor(kGreen),LineStyle(kDashed)) ;
  added_pdf->plotOn(deframe,Components(signal_pdf_2),LineColor(kBlack),LineStyle(kDashed)) ;
  added_pdf->plotOn(deframe,Components(total_signal_pdf),LineColor(kRed),LineStyle(kDashed));
  added_pdf->plotOn(deframe,Components(backround_pdf),LineColor(kMagenta),LineStyle(kDashed)) ;
  added_pdf->plotOn(deframe, LineColor(kBlue), LineStyle(kSolid)) ; // we need to write it last.
                        // otherwise pull distribution will calculate error wrt last mentioned pdf inside plotOn function
  added_pdf->paramOn(deframe,dataset_name,"", 2, "NEU", 0.1, 0.35, 0.9); //"NELU",  Prints the fitted parameter on the canvas

  //Extract info. from fitting frame and showing
  RooHist* hpull = deframe->pullHist() ;
  RooPlot* frame3 = var_name->frame(Title("Pull Distribution")) ;
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
  TLegendEntry *entry = legend1->AddEntry("signal_pdf_1","1st Gaussian pdf","l");
  entry->SetLineColor(kGreen);
  entry->SetLineStyle(kDashed);
  entry = legend1->AddEntry("signal_pdf_2","2nd Gaussian pdf","l");
  entry->SetLineColor(kBlack);
  entry->SetLineStyle(kDashed);
  entry = legend1->AddEntry("total_signal_pdf","Combined signal pdf","l");
  entry->SetLineColor(kRed);
  entry->SetLineStyle(kDashed);
  entry = legend1->AddEntry("backround_pdf","bkg-ARGUS","l");
  entry->SetLineColor(kMagenta);
  entry->SetLineStyle(kDashed);
  entry = legend1->AddEntry("added_pdf","Fitted pdf","l");
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

  // cout<<"Total number of events which are used to fitting are : "<<counter<<endl;
  cout<<"chisq of the fit is :"<<deframe->chiSquare()<<endl;//chi-square of the fit
  cout<<"chi-square/ndof :"<<deframe->chiSquare(7)<<endl;// Chi^2/(the number of degrees of freedom)
  cout<<"Error in calculation of signal yield : "<<n_sig->getError()<<endl;
  
  const char* output_filename = output_name;
  c1->Print(output_filename);
  return 0;
}