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

using namespace RooFit ;
using namespace std;
//int main(){

void DE_fit(){


/*******************Fit Variables***********************************/

RooRealVar deltae("deltae","#DeltaE (GeV)", -0.15, 0.15);
RooRealVar bdt("bdt","FBDT output", 0., 1.);
RooDataSet* data=new RooDataSet("data","data",RooArgSet(deltae,bdt));

/*******************Input root file**********************************/
TChain* chain=new TChain();
chain->Add("/home/niharika/belle_II/analysis/B2BII/NewRelease/CS/FBDT_test_variables/data_MC_comparison_variables/data_exp55.root/bp3");  // to take DeltaE valriable from it
chain->AddFriend("variables","/home/niharika/belle_II/analysis/B2BII/NewRelease/CS/FBDT_test_variables/data_MC_comparison_FBDT_output/expert_dropped_data.root"); // to take MVA output from it
Float_t  cs1;
Double_t o_de, o_md0, o_mbc, o_r2, o_kid,o_pid,sig,o_costbto,k_TOPAcc,pi_TOPAcc;
Int_t nevt=(int)chain->GetEntries();

chain->SetBranchAddress("deltaE",&o_de);  // linking with DeltaE distribution
chain->SetBranchAddress("variables.MVADatabaseIdentifier2", &cs1);  // linking with MVA output

//Loading data 
for(int i=0;i<nevt;i++) {
  chain->GetEntry(i);
  deltae.setVal(o_de);
  bdt.setVal(cs1);
  
    data->add(RooArgSet(deltae,bdt));
}


/*****************************Fit***********************/

RooRealVar mean_deltae("mean_{#DeltaE}","", 0.0, -0.05, 0.05);
RooRealVar f_de("f_de", "f", 1., 0., 10.);
 RooRealVar sigma0("#sigma0", "sigma", 0.00978, 0., 0.05);

RooFormulaVar sigma("#sigma", "@0*@1", RooArgList(sigma0,f_de));
RooGaussian gauss1("gauss1", "gauss1", deltae, mean_deltae, sigma0);

 RooRealVar b0("b0", "b0", -0.5910, -10., 10.);
RooRealVar b1("b1", "b1", -0.062, -10., 10.);

RooChebychev bkg("bkg","Background",deltae,RooArgSet(b0)) ;

RooRealVar nsig("nsig", "nsig", 500, -100., 10000.0);//52000
RooRealVar ncomb("ncomb", "ncombD0", 500, -100., 100000.0);//95000

RooAddPdf sum("sum","sum",RooArgList(gauss1,bkg),RooArgList(nsig, ncomb));

sum.fitTo(*data);
/*****************************Fit***********************/
//signal region

RooPlot* deframe = deltae.frame() ;
data->plotOn(deframe) ;
sum.plotOn(deframe) ;
sum.paramOn(deframe,data);
Double_t chisq = deframe->chiSquare();
RooHist* hpull = deframe->pullHist() ;
RooPlot* frame3 = deltae.frame(Title("Pull Distribution")) ;
  frame3->addPlotable(hpull,"P") ;
  
  //sum.plotOn(deframe,Components(gauss1),LineColor(kRed),LineStyle(kDashed)) ;
//sum.plotOn(deframe,Components(signal_de),LineColor(kGreen),LineStyle(kDashed)) ;
sum.plotOn(deframe,Components(gauss1),LineColor(kRed),LineStyle(kDashed)) ;
sum.plotOn(deframe,Components(bkg),LineColor(kMagenta),LineStyle(kDashed)) ;

TCanvas* c1 = new TCanvas() ;
 TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();  
  deframe->Draw() ;
  

  c1->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
pad2->Draw();
   pad2->cd(); 
 frame3->Draw() ;

cout << "chi2 de=" << chisq << endl;






/*****************************Splot***********************/
/*****************************Splot***********************/
// Create the sPlot data from the fit to mass.
  RooStats::SPlot* sDataX = new RooStats::SPlot("sData","An SPlot", *data, &sum, RooArgList(nsig, ncomb) );
  cout << "Yield of signal is " << nsig.getVal() << " From sWeights it is " << sDataX->GetYieldFromSWeight("nsig") <<endl;
  cout << "Yield of background is " << ncomb.getVal() << "  From sWeights it is " << sDataX->GetYieldFromSWeight("ncomb") <<endl;

 // cout << "signal Weight   " << sDataX->GetSWeight("nsig") 
 //      << " bkg Weight   " <<  sDataX->GetSWeight("nbkg") <<endl;
   
  for(Int_t i=0; i < 20; i++)
    {
      cout << "signal Weight   " << sDataX->GetSWeight(i,"nsig") 
		<< " bkg Weight   " << sDataX->GetSWeight(i,"ncomb") 
		<< "  Total Weight   " << sDataX->GetSumOfEventSWeight(i) 
		<< endl;
    }
  
  // create weighted data sets
  RooDataSet * dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"nsig_sw") ;
  RooDataSet * dataw_bg  = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"ncomb_sw") ;
  // How signal and bkg are seperately defined here???????????????????????????

  cout << "Making splots of the signal and background" <<endl;
  RooPlot* frame_sig_rep = bdt.frame(50) ;  // Defining frame(from plot of bdt output) to plot signal distribution
  frame_sig_rep->SetTitle("sPlot for the M_{#pi^{0}} distribution(signal)");
  dataw_sig->plotOn(frame_sig_rep, RooFit::DataError(RooAbsData::SumW2) ) ;//???????????????????????????

  RooPlot* frame_bg_rep = bdt.frame(50) ; // Defining frame(from plot of bdt output) to plot bkg distribution
  frame_bg_rep->SetTitle("sPlot for the M_{#pi^{0}} distribtuion(background)");
  dataw_bg->plotOn(frame_bg_rep, RooFit::DataError(RooAbsData::SumW2) ) ;//???????????????????????????

  // plot PDFs onto the sPlots; if the sPlot correctly reweights the data then the signal / background only shapes should
  // provide an adequate description of the data.


  TCanvas *c2 = new TCanvas("c2" ,"c2");
  c2->cd();
  frame_sig_rep->Draw();
 
  TCanvas *c3 = new TCanvas("c3","c3");
  c3->cd();
  frame_bg_rep->Draw();
  
 // c1.Print("splot_pt.pdf");
 c2->SaveAs("mc_splot_rdedx_sig.root");
 c3->SaveAs("mc_splot_rdedx_bkg.root");

 //  RooDataHist *dh1 = new RooDataHist("dh1"," dh1",RooArgSet(muid),*dataw_sig);
   bdt.setBins(50);
   TH1F *h1 = (TH1F*)dataw_sig->createHistogram("signal" ,bdt);//???????????????????????????
   double scale = 1./h1->Integral();
   h1->Scale(scale);
   h1->SetLineColor(kBlack);
   h1->SetLineWidth(3);

 //  RooDataHist *dh1_b = new RooDataHist("dh1_b"," dh1_b",RooArgSet(muid),*dataw_bg);
   TH1F *h1_b = (TH1F*)dataw_bg->createHistogram("background" ,bdt);//???????????????????????????
   double scaleb = 1./h1_b->Integral();
   h1_b->Scale(scaleb);
   h1_b->SetLineColor(kBlack);
   h1_b->SetLineWidth(3);

///////////////////////////plotting MVA output in histogram///////////////////////////////////////////////////////////
TFile *f2 = new TFile("/home/niharika/belle_II/analysis/B2BII/NewRelease/CS/FBDT_test_variables/data_MC_comparison_FBDT_output/expert_dropped_generic_Dpi.root");
TNtuple *tree2 = (TNtuple*)f2->Get("variables");
Float_t out2;
Float_t sig2;
tree2 ->SetBranchAddress("MVADatabaseIdentifier2", &out2);//bkg
tree2 ->SetBranchAddress("MVADatabaseIdentifier2_isSignal", &sig2);//signal
TH1F *hint1 = new TH1F("hint1","",50,0,1);
TH1F *hint2 = new TH1F("hint2","",50,0,1);
Long64_t nentries2 = tree2->GetEntries();
Long64_t nbytes2 = 0;//?????????????????????????????
for (Long64_t j=0; j<nentries2;j++) {
  nbytes2 += tree2->GetEntry(j);//????????????????????
  if(sig2==1) hint1 ->Fill(out2);//signal /* out2 should be sig2 BUT IF STATEMENT DOSE THE JOB*/
  if(sig2==0) hint2 ->Fill(out2);//bkg
}
hint1->SetLineColor(kRed);
hint1->SetLineWidth(3);
Double_t scale2 = 1./hint1->Integral();
hint1->Scale(scale2);

hint2->SetLineColor(kBlue);
hint2->SetLineWidth(3);
Double_t scale4 = 1./hint2->Integral();
hint2->Scale(scale4);
////////////////////////////////plotting MVA output in histogram******end//////////////////////////////////////////////



  TCanvas *c4 = new TCanvas("c4" ,"c4");
  c4->cd();
  h1->Draw("");
  h1_b->Draw("sames");
  hint1->Draw("sames");
  hint2->Draw("sames");
  auto legend = new TLegend( 0.6,0.7,0.9,0.85);
  legend->AddEntry(hint1," MC signal","l");
  legend->AddEntry(hint2," MC background","l");
  legend->AddEntry(h1,"D#pi data","l");
  //legend->AddEntry(h1_b," Data background","l");
  legend->Draw();
   c4->SaveAs("splot_sig_bkg.root");

}

