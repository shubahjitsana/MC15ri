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
using namespace std;
using namespace RooFit;
using namespace RooStats;

void splot_define_n_plot(RooRealVar* discriminating_var, RooRealVar* control_var, RooDataSet* dataset_name,
  TH1F* control_var_histogram_for_signal, TH1F* control_var_histogram_for_bkg, // TChain* chain,
  RooRealVar* n_sig, RooRealVar* n_bkg,
  RooAddPdf* added_pdf)
  {
  // Create the sPlot data from the DATA of DeltaE fitting
  RooStats::SPlot* sDataX = new RooStats::SPlot("sData","An SPlot", *dataset_name, added_pdf, RooArgList(n_sig, n_bkg) );
  cout << "Yield of signal is " << n_sig->getVal() << " From sWeights it is " << sDataX->GetYieldFromSWeight("n_sig") <<endl;
  cout << "Yield of background is " << n_bkg->getVal() << "  From sWeights it is " << sDataX->GetYieldFromSWeight("n_bkg") <<endl;
  
  // Printing different sWeight
  for(Int_t i=0; i < 20; i++){
      cout << "signal Weight   " << sDataX->GetSWeight(i,"n_sig") 
      << " bkg Weight   " << sDataX->GetSWeight(i,"n_bkg") 
      << "  Total Weight   " << sDataX->GetSumOfEventSWeight(i) 
      << endl;
  }

  // create weighted data sets
  //RooDataSet(const char* name, const char* title, RooDataSet* data, const RooArgSet& vars, const RooFormulaVar& cutVar, const char* wgtVarName = 0)
  RooDataSet * dataw_sig = new RooDataSet(dataset_name->GetName(),dataset_name->GetTitle(),dataset_name,*dataset_name->get(),0,"n_sig_sw") ;
  RooDataSet * dataw_bg  = new RooDataSet(dataset_name->GetName(),dataset_name->GetTitle(),dataset_name,*dataset_name->get(),0,"n_bkg_sw") ;
  /////////////////////////////////////////////////////////////////////////////???????????????////////////////////////////////////
  // How signal and bkg are seperately defined here???????????????????????????here how r we getting different variables from Roodataset in two different line?


  //Defining frame to plot sWeighted pdf
  cout << "Making splots of the signal and background" <<endl;
  RooPlot* frame_sig_rep = control_var->frame(50) ;  // Defining frame(from plot of bdt output) to plot signal distribution
  frame_sig_rep->SetTitle("sPlot for the Signal Probability distribution(signal)");
  dataw_sig->plotOn(frame_sig_rep, RooFit::DataError(RooAbsData::SumW2) ) ;//???????????????????????????

  RooPlot* frame_bg_rep = control_var->frame(50) ; // Defining frame(from plot of bdt output) to plot bkg distribution
  frame_bg_rep->SetTitle("sPlot for the Signal Probability distribtuion(background)");
  dataw_bg->plotOn(frame_bg_rep, RooFit::DataError(RooAbsData::SumW2) ) ;//???????????????????????????

  // plot PDFs onto the sPlots; if the sPlot correctly reweights the data then the signal / background only shapes should
  // provide an adequate description of the data.
  TCanvas *c2 = new TCanvas("c2" ,"c2");
  c2->cd();
  frame_sig_rep->Draw();

  TCanvas *c3 = new TCanvas("c3","c3");
  c3->cd();
  frame_bg_rep->Draw();

  c2->SaveAs("mc_splot_sig.png");
  c3->SaveAs("mc_splot_bkg.png");


  // Plot sWeighted distriution of control variable(MVA output)
  // in histogram along with original MVA output
  // to compare sPlot output distrubution and MVA output distribution
  control_var->setBins(50);
  TH1F *control_var_sWeighted_histogram_for_signal = (TH1F*)dataw_sig->createHistogram("signal" ,*control_var);//???????????????????????????
  double scale = 1./control_var_sWeighted_histogram_for_signal->Integral();
  control_var_sWeighted_histogram_for_signal->Scale(scale);
  control_var_sWeighted_histogram_for_signal->SetLineColor(kBlack);
  control_var_sWeighted_histogram_for_signal->SetLineWidth(5);

  TH1F *control_var_sWeighted_histogram_for_bkg = (TH1F*)dataw_bg->createHistogram("background" ,*control_var);//???????????????????????????
  double scaleb = 1./control_var_sWeighted_histogram_for_bkg->Integral();
  control_var_sWeighted_histogram_for_bkg->Scale(scaleb);
  control_var_sWeighted_histogram_for_bkg->SetLineColor(kRed);
  control_var_sWeighted_histogram_for_bkg->SetLineWidth(2);


  ///////////////////////////plotting MVA output in histogram///////////////////////////////////////////////////////////
  // To do this we need just TChain as output. But this will take nearly double time as for loop will run
  // double time
  // TH1F *control_var_histogram_for_signal = new TH1F("control_var_histogram_for_signal","",50,0,1);
  // TH1F *control_var_histogram_for_bkg = new TH1F("control_var_histogram_for_bkg","",50,0,1);
  // for(int l=0;l<(int)chain->GetEntries();l++){
  //     chain->GetEntry(l);
  //     if(sig==1) control_var_histogram_for_signal->Fill(sig_prob);//SIG
  //     if(sig==0) control_var_histogram_for_bkg->Fill(sig_prob);//bkg
  // }
  control_var_histogram_for_signal->SetLineColor(kMagenta);
  control_var_histogram_for_signal->SetLineWidth(2);
  Double_t scale2 = 1./control_var_histogram_for_signal->Integral();
  control_var_histogram_for_signal->Scale(scale2);

  control_var_histogram_for_bkg->SetLineColor(kBlue);
  control_var_histogram_for_bkg->SetLineWidth(2);
  Double_t scale4 = 1./control_var_histogram_for_bkg->Integral();
  control_var_histogram_for_bkg->Scale(scale4);
  ////////////////////////////////plotting MVA output in histogram******end//////////////////////////////////////////////


  TCanvas *c4 = new TCanvas("c4" ,"c4");
  c4->cd();

  control_var_sWeighted_histogram_for_signal->SetStats(kFALSE);
  control_var_sWeighted_histogram_for_bkg->SetStats(kFALSE);
  control_var_histogram_for_signal->SetStats(kFALSE);
  control_var_histogram_for_bkg->SetStats(kFALSE);

  control_var_sWeighted_histogram_for_signal->Draw();
  control_var_sWeighted_histogram_for_bkg->Draw("SAME");
  control_var_histogram_for_signal->Draw("SAME");//SIG
  control_var_histogram_for_bkg->Draw("SAME");//bkg
  
  TLegend* legend = new TLegend( 0.25,0.7,0.74,0.9);
  legend->AddEntry(control_var_histogram_for_signal,"Signal Probability predicted by MVA for signal","l");
  legend->AddEntry(control_var_histogram_for_bkg,"Signal Probability predicted by MVA for background","l");
  legend->AddEntry(control_var_sWeighted_histogram_for_signal,"Signal Probability predicted by sPlot for signal","l");
  legend->AddEntry(control_var_sWeighted_histogram_for_bkg,"Signal Probability predicted by sPlot for background","l");
  legend->Draw();
  c4->SaveAs("splot_sig_bkg.png");

  return 0;
  }