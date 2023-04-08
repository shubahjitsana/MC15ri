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

void defit_signalfolder(){
    /*******************Fit Variables***********************************/
    RooRealVar deltae_signalfolder("deltae","#DeltaE (GeV)", -0.1, 0.1);
    /**defining DATAFRAME{unbinned histogram}(FROM FIT VARIABLE) TO FIT AND PLOT**/
    RooDataSet* data_signalfolder=new RooDataSet("data","data",RooArgSet(deltae_signalfolder));
    /*******************Input root file**********************************/
    TChain* chain_signalfolder=new TChain();
    chain_signalfolder->Add("/home/belle2/ssana/MC15ri/cs/test/prescale_combine/test_signal_TM_scaled.root/tree");


    Double_t  de3_signalfolder, md03_signalfolder, mbc3_signalfolder, kid3_signalfolder,pid3_signalfolder,sig_signalfolder,sig_prob_signalfolder;
    Double_t D_s_InvM_signalfolder, D_10D_md_signalfolder, InvM1stand2ndpi_signalfolder, InvM2ndand3rdpi_signalfolder;
    Double_t InvMD0and1stpi_signalfolder, InvMD0and2ndpi_signalfolder, InvMD0and3rdpi_signalfolder, DstrminusroeD_md_signalfolder;
    Int_t run_signalfolder;
    Int_t nevt3_signalfolder=(int)chain_signalfolder->GetEntries();

    // chain->SetBranchAddress("isSignal",&sig);
    chain_signalfolder->SetBranchAddress("deltaE",&de3_signalfolder);
    chain_signalfolder->SetBranchAddress("Mbc",&mbc3_signalfolder);
    chain_signalfolder->SetBranchAddress("D0_bar_InvM",&md03_signalfolder);
    chain_signalfolder->SetBranchAddress("Kp_PID_bin_kaon",&kid3_signalfolder);
    chain_signalfolder->SetBranchAddress("SigProb",&sig_prob_signalfolder);
    chain_signalfolder->SetBranchAddress("D_s_InvM",&D_s_InvM_signalfolder);
    chain_signalfolder->SetBranchAddress("D_10D_md",&D_10D_md_signalfolder);
    chain_signalfolder->SetBranchAddress("InvM1stand2ndpi",&InvM1stand2ndpi_signalfolder);
    chain_signalfolder->SetBranchAddress("InvM2ndand3rdpi",&InvM2ndand3rdpi_signalfolder);
    chain_signalfolder->SetBranchAddress("InvMD0and1stpi",&InvMD0and1stpi_signalfolder);
    chain_signalfolder->SetBranchAddress("InvMD0and2ndpi",&InvMD0and2ndpi_signalfolder);
    chain_signalfolder->SetBranchAddress("InvMD0and3rdpi",&InvMD0and3rdpi_signalfolder);
    chain_signalfolder->SetBranchAddress("DstrminusroeD_md",&DstrminusroeD_md_signalfolder);
    // chain_signalfolder->SetBranchAddress("",&);

    
    // D0_bar_InvM >1.84 && D0_bar_InvM <1.89 && Mbc>5.27 && Mbc < 5.29 && deltaE < 0.1 && deltaE > -0.1 && Kp_PID_bin_kaon > 0.6 && SigProb < 0.86
    //Loading data 
    Double_t counter_signalfolder =0;
    for(int l=0;l<nevt3_signalfolder;l++) {
        chain_signalfolder->GetEntry(l);
        deltae_signalfolder.setVal(de3_signalfolder);
        if(mbc3_signalfolder>5.23 && mbc3_signalfolder < 5.29 && 
        de3_signalfolder < 0.1 && de3_signalfolder > -0.1 &&
        md03_signalfolder>1.85 && md03_signalfolder<1.88 &&  
        kid3_signalfolder > 0.6 && 
        (D_s_InvM_signalfolder < 1.958 || D_s_InvM_signalfolder > 1.978) &&
        (D_10D_md_signalfolder < 0.496 || D_10D_md_signalfolder > 0.628) &&
        (InvM1stand2ndpi_signalfolder < 0.489 || InvM1stand2ndpi_signalfolder > 0.506) &&
        (InvM2ndand3rdpi_signalfolder < 0.49 || InvM2ndand3rdpi_signalfolder > 0.505) &&
        (InvMD0and1stpi_signalfolder < 2.0085 || InvMD0and1stpi_signalfolder > 2.0122) &&
        (InvMD0and2ndpi_signalfolder < 2.0087 || InvMD0and2ndpi_signalfolder > 2.0118) &&
        (InvMD0and3rdpi_signalfolder < 2.009 || InvMD0and3rdpi_signalfolder > 2.0116) &&
        (DstrminusroeD_md_signalfolder < 0.14402 || DstrminusroeD_md_signalfolder > 0.14682) &&
        sig_prob_signalfolder > 0.21){
            data_signalfolder->add(RooArgSet(deltae_signalfolder));
            counter_signalfolder++;
        }
    }


    /*****************************Delta E Fit***********************/
    // --- Build Signal PDF ---
    RooRealVar mean1_signalfolder("mean1","mean of Gaussian-1",-0.001,-0.02,0.02);
    RooRealVar mean2_signalfolder("mean2","mean of Gaussian-2",0.002,-0.02,0.02);
    RooRealVar sigma1_signalfolder("sigma1","sigma of Gaussian-1",0.01, 0., 1);	
    RooRealVar sigma2_signalfolder("sigma2","sigma of Gaussian-2",0.01191,0.00000001,1);

    RooGaussian sig1_signalfolder("sig1","Gaussian-1",deltae_signalfolder,mean1_signalfolder,sigma1_signalfolder);  
    RooGaussian sig2_signalfolder("sig2","Gaussian-2",deltae_signalfolder,mean2_signalfolder,sigma2_signalfolder);

    RooRealVar fsig_1_signalfolder("frac_gaussians", "signal fraction", 0.7,0.,1.);
    RooAddPdf sum_signalfolder("twoGaussians", "sum of two Gaussians ",RooArgList(sig1_signalfolder, sig2_signalfolder), RooArgList(fsig_1_signalfolder));
    sum_signalfolder.fitTo(*data_signalfolder);
    /****************************FIT COMPLETE*************************************/
    //Plotting fitted result
    RooPlot* deframe_signalfolder = deltae_signalfolder.frame(Title("Fitting #DeltaE for signal MC"), Bins(300)) ;                          
    data_signalfolder->plotOn(deframe_signalfolder, Binning(300), DataError(RooAbsData::SumW2)) ;
    sum_signalfolder.plotOn(deframe_signalfolder,Components(sig1_signalfolder),LineColor(kGreen),LineStyle(kDashed)) ;
    sum_signalfolder.plotOn(deframe_signalfolder,Components(sig2_signalfolder),LineColor(kBlack),LineStyle(kDashed)) ;
    sum_signalfolder.plotOn(deframe_signalfolder, LineColor(kBlue), LineStyle(kSolid)) ; // we need to write it last.
                        // otherwise pull distribution will calculate error wrt last mentioned pdf inside plotOn function
    sum_signalfolder.paramOn(deframe_signalfolder,data_signalfolder,"", 2, "NEU", 0.63, 0.9, 0.9); //"NELU",  Prints the fitted parameter on the canvas

    //Extract info. from fitting frame and showing
    RooHist* hpull_signalfolder = deframe_signalfolder->pullHist() ;
    RooPlot* frame3_signalfolder = deltae_signalfolder.frame(Title("Pull Distribution")) ;
    hpull_signalfolder->SetFillColor(1);
    frame3_signalfolder->addPlotable(hpull_signalfolder,"X0B"); // "X0" is for errorbar; and "B" is for histograms


    TCanvas* c1 = new TCanvas("c1", "c1", 2550, 1500) ;
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();  
    deframe_signalfolder->Draw() ;
    // Adding legend
    TLegend *legend1_signalfolder = new TLegend(0.1,0.7,0.25,0.9);
    TLegendEntry *entry_signalfolder = legend1_signalfolder->AddEntry("sig1_signalfolder","1st Gaussian pdf","l");
    entry_signalfolder->SetLineColor(kGreen);
    entry_signalfolder->SetLineStyle(kDashed);
    entry_signalfolder = legend1_signalfolder->AddEntry("sig2_signalfolder","2nd Gaussian pdf","l");
    entry_signalfolder->SetLineColor(kBlack);
    entry_signalfolder->SetLineStyle(kDashed);
    entry_signalfolder = legend1_signalfolder->AddEntry("sum_signalfolder","Fitted pdf","l");
    entry_signalfolder->SetLineColor(kBlue);
    entry_signalfolder->SetLineStyle(kSolid);
    legend1_signalfolder->Draw();

    // // Adding arrow at (+-)3*sigma of signal pdf
    // double weightedMean = mean1_signalfolder.getVal()*fsig_1_signalfolder.getVal() + mean2_signalfolder.getVal()*(1.0-fsig_1_signalfolder.getVal());
    // double weightedSigma = std::sqrt( pow(sigma1_signalfolder.getVal(),2)*pow(fsig_1_signalfolder.getVal(),2)  +  pow(sigma2_signalfolder.getVal(),2)*pow((1-fsig_1_signalfolder.getVal()),2));
    // double deltae_min_fit = weightedMean -3*weightedSigma; 
    // double deltae_max_fit = weightedMean +3*weightedSigma;
    // double arrowHeight = n_sig_signalfolder.getVal()*0.03;
    // cout<<"DelataE min :"<< deltae_min_fit <<endl;
    // cout<<"DelataE max :"<< deltae_max_fit <<endl;
    // cout<<"Arrow Height :"<< arrowHeight <<endl;

    // TArrow *arr1 = new TArrow(-0.01, 100, -0.01, 400, 0.01,"|>");
    // arr1->SetLineWidth(4);
    // arr1->SetLineColor(2);
    // arr1->SetFillStyle(3008);
    // arr1->DrawArrow(deltae_min_fit, arrowHeight*0.1, deltae_min_fit, arrowHeight, 0.02,"<"); 
    // arr1->DrawArrow(deltae_max_fit, arrowHeight*0.1, deltae_max_fit, arrowHeight, 0.02,"<");

    c1->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->Draw();
    pad2->cd(); 
    frame3_signalfolder->Draw() ;
    frame3_signalfolder->SetLineStyle(9);
    frame3_signalfolder->GetYaxis()->SetNdivisions(505);
    frame3_signalfolder->GetYaxis()->SetTitle("#sqrt{#chi^{2}}"); 
    frame3_signalfolder->GetXaxis()->SetTitle("#DeltaE (GeV)"); 
    frame3_signalfolder->GetXaxis()->SetTitleSize(0.13);
    frame3_signalfolder->GetYaxis()->SetTitleSize(0.15);
    frame3_signalfolder->GetXaxis()->SetLabelSize(0.120);
    frame3_signalfolder->GetYaxis()->SetLabelSize(0.120); 
    frame3_signalfolder->GetXaxis()->SetTitleOffset(0.90);      
    frame3_signalfolder->GetYaxis()->SetTitleOffset(0.22);       
    frame3_signalfolder->GetYaxis()->SetRangeUser(-10.0, 10.0);       
    frame3_signalfolder->GetYaxis()->SetLimits(-10.0, 10.0);       
    frame3_signalfolder->GetXaxis()->SetNdivisions(505);
    frame3_signalfolder->GetYaxis()->CenterTitle(true);
    frame3_signalfolder->GetXaxis()->CenterTitle(true);
    frame3_signalfolder->Draw("AXISSAME");

    cout<<"Total number of events which are used to fitting are : "<<counter_signalfolder<<endl;
    cout<<"chisq of the fit is :"<<deframe_signalfolder->chiSquare()<<endl;//chi-square of the fit
    cout<<"chi-square/ndof :"<<deframe_signalfolder->chiSquare(7)<<endl;// Chi^2/(the number of degrees of freedom)

    c1->Print("de_plot/de_fit_signalfolder.png");
}