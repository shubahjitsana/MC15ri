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
using namespace RooStats;

void define_var_dataset(RooWorkspace* ws){
    /*******************Fit Variables***********************************/
    RooRealVar deltae("deltae","#DeltaE (GeV)", -0.1, 0.1);
    /**defining DATAFRAME{unbinned histogram}(FROM FIT VARIABLE) TO FIT AND PLOT**/
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(deltae));

    // import data into workspace
    ws->import(*data, Rename("data_just_defined"));    
}
void import_root_to_dataset(RooWorkspace* ws){
    // get what we need out of the workspace for import_root_to_dataset
    RooRealVar* deltae = ws->var("deltae");
    RooDataSet* data = (RooDataSet*) ws->data("data_just_defined");  
    /*******************Input root file**********************************/
    TChain* chain=new TChain();
    chain->Add("/home/belle2/ssana/MC15ri/cs/test/prescale_combine/all.root/tree");

    Double_t  de3, md03, mbc3, r23, kid3,pid3,sig,sig_prob;
    Int_t run;
    Int_t nevt3=(int)chain->GetEntries();

    chain->SetBranchAddress("deltaE",&de3);
    chain->SetBranchAddress("Mbc",&mbc3);
    chain->SetBranchAddress("D0_bar_InvM",&md03);
    chain->SetBranchAddress("Kp_PID_bin_kaon",&kid3);
    chain->SetBranchAddress("SigProb",&sig_prob);
    chain->SetBranchAddress("isSignal",&sig);
    
    // D0_bar_InvM >1.84 && D0_bar_InvM <1.89 && Mbc>5.27 && Mbc < 5.29 && deltaE < 0.1 && deltaE > -0.1 && Kp_PID_bin_kaon > 0.6 && ContProb < 0.86
    //Loading data 
    Double_t counter =0;
    for(int l=0;l<nevt3;l++) {
        chain->GetEntry(l);
        deltae->setVal(de3);
        if(mbc3>5.23 && mbc3 < 5.29 &&
         de3 < 0.1 && de3 > -0.1 &&
         md03>1.85 && md03<1.88 &&
         kid3 > 0.6 && 
         sig_prob > 0.22
         ){
            data->add(RooArgSet(*deltae));
            counter++;
        }
    } 
    std::cout<<"data load successfull"<<endl;
    // import data into workspace
    ws->import(*data, Rename("data_having_info"));
    RooRealVar get_counter("counter","counter",counter);
    ws->import(get_counter);
}
void prefit_plot(RooWorkspace* ws){
    // get what we need out of the workspace for prefit_plot
    RooRealVar* var_name = ws->var("deltae");
    RooDataSet* dataset_name = (RooDataSet*) ws->data("data_having_info");

    RooDataHist* binDataSet = new RooDataHist("binDataSet", "binDataSet", *var_name, *dataset_name);
    RooPlot* xframe1 = var_name->frame(	Title("prefit")	, Bins(200));
    binDataSet->plotOn(	xframe1	, Binning(200)		, DataError(RooAbsData::SumW2));
    xframe1->Draw();

    TCanvas *prefit_canvas = new TCanvas("prefit_canvas", "", 1500, 1500);  
    prefit_canvas -> Print("prefit_histo.png");
}
void define_pdf(RooWorkspace* ws){
    // get what we need out of the workspace for import_root_to_dataset
    RooRealVar* deltae = ws->var("deltae");
    RooRealVar* get_counter = ws->var("get_counter");
    // Double_t* counter = get_counter->getVal();
    // int counter = (int)ws->var("get_counter")->getVal();
    RooDataSet* data = (RooDataSet*) ws->data("data_having_info");   

    /*****************************Delta E Fit***********************/
    // --- Build Signal PDF ---
    RooRealVar mean1("mean1","mean of Gaussian-1",-0.001,-0.02,0.02);
    RooRealVar mean2("mean2","mean of Gaussian-2",0.002,-0.02,0.02);
    RooRealVar sigma1("sigma1","sigma of Gaussian-1",0., 0., 1);	
    RooRealVar sigma2("sigma2","sigma of Gaussian-2",0.01191,0.00000001,1);

    RooGaussian sig1("sig1","Gaussian-1",*deltae,mean1,sigma1);  
    RooGaussian sig2("sig2","Gaussian-2",*deltae,mean2,sigma2);

    RooRealVar fsig_1("frac_gaussians", "signal fraction", 0.5,0.,1.);
    RooAddPdf twoGaussians("twoGaussians", "sum of two Gaussians ",RooArgList(sig1, sig2), RooArgList(fsig_1));

    // --- Build Argus background PDF ---
    RooRealVar b1("Chbyshv-prm", "Chbyshv-prm", -0.062, -10., 10.);
    RooChebychev bkg("bkg","Background",*deltae,RooArgSet(b1)) ;

    //Initialization of parameter before adding two pdf
    // cout<<"Total number of events which are used to fitting are : "<<counter<<endl;
    Double_t event_count = get_counter->getVal(); 
    Double_t signal_count = get_counter->getVal()*0.2;
    Double_t back_count = get_counter->getVal()*0.8;
    RooRealVar n_sig("n_sig", "n_sig", signal_count, 0., event_count);
    RooRealVar n_bkg("n_bkg", "n_bkg", back_count, 0., event_count);
    // RooRealVar n_sig("n_sig", "n_sig", 145795, 14500, 146600);
    // RooRealVar n_bkg("n_bkg", "n_bkg", 157295, 156400, 158100);
    RooAddPdf sum("sum","sum",RooArgList(twoGaussians,bkg),RooArgList(n_sig, n_bkg));//adding two pdf

    ws->import(sum);
}

void de_fit(RooWorkspace* ws){
    // get what we need out of the workspace for import_root_to_dataset
    RooAbsPdf* sum = ws->pdf("sum");
    RooRealVar* get_counter = ws->var("get_counter");
    // RooRealVar* n_sig = ws->var("n_sig");
    // RooRealVar* n_bkg = ws->var("n_bkg");
    // int counter = (int)ws->var("get_counter")->getVal();
    RooDataSet* data = (RooDataSet*) ws->data("data_having_info");  


    RooFitResult * my_fit_result = (RooFitResult*)sum->fitTo(*data," ", RooFit::Save());   //this will fit only DELTA_E as we mentioned "the fitting variable(deltae)" during defining each elementary pdf
    // RooMsgService::instance().setSilentMode(true); //to turn off all informational 
         //and warning messages that would normally be printed to the console 
         //during a RooFit calculation, while still allowing error messages to be
         //printed to the console. This can be useful in cases where the user 
         //wants to suppress output messages during a batch job or automated 
         //script, or when the user is only interested in the final result of a calculation.    
    cout<<"Total number of events which are used to fitting are : "<<get_counter->getVal()<<endl;
    cout<<"Printing fitted result.........................."<<endl;
    my_fit_result->Print("v");  
}


void de_fit_and_get_sWeights(RooWorkspace* ws){
    // get what we need out of the workspace for import_root_to_dataset
    RooAbsPdf* sum = ws->pdf("sum");
    RooRealVar* get_counter = ws->var("get_counter");
    RooRealVar* n_sig = ws->var("n_sig");
    RooRealVar* n_bkg = ws->var("n_bkg");
    // int counter = (int)ws->var("get_counter")->getVal();
    RooDataSet* data = (RooDataSet*) ws->data("data_having_info");  


    RooFitResult * my_fit_result = (RooFitResult*)sum->fitTo(*data," ", RooFit::Save());   //this will fit only DELTA_E as we mentioned "the fitting variable(deltae)" during defining each elementary pdf
    // RooMsgService::instance().setSilentMode(true); //to turn off all informational 
         //and warning messages that would normally be printed to the console 
         //during a RooFit calculation, while still allowing error messages to be
         //printed to the console. This can be useful in cases where the user 
         //wants to suppress output messages during a batch job or automated 
         //script, or when the user is only interested in the final result of a calculation.    
    cout<<"Total number of events which are used to fitting are : "<<get_counter->getVal()<<endl;
    cout<<"Printing fitted result.........................."<<endl;
    my_fit_result->Print("v");  

    /*****************************Splot***********************/

    // Create the sPlot data from the DATA of DeltaE fitting
    RooStats::SPlot* sDataX = new RooStats::SPlot("sData","An SPlot", *data, sum, RooArgList(n_sig, n_bkg) );

    // Check that our weights have the desired properties, which is vlaue of Yield and sWeights are nearly same.
    cout << "Check that our weights have the desired properties, which is vlaue of Yield and sWeights are nearly same." <<endl;
    cout << "Yield of signal is " << n_sig->getVal() << " From sWeights it is " << sDataX->GetYieldFromSWeight("n_sig") <<endl;
    cout << "Yield of background is " << n_bkg->getVal() << "  From sWeights it is " << sDataX->GetYieldFromSWeight("n_bkg") <<endl;
   
    // Printing different sWeight
    for(Int_t i=0; i < 20; i++){
        cout << "signal Weight   " << sDataX->GetSWeight(i,"n_sig") 
        << " bkg Weight   " << sDataX->GetSWeight(i,"n_bkg") 
        << "  Total Weight   " << sDataX->GetSumOfEventSWeight(i) 
        << endl;
    }
    // import this new dataset with sWeights
    ws->import(*data, Rename("data_having_info_and_sWeights"));    
}


void de_plot_style(RooWorkspace* ws){
    // get what we need out of the workspace for import_root_to_dataset
    RooRealVar* var_name = ws->var("deltae");
    RooDataSet* dataset_name = (RooDataSet*) ws->data("data_having_info");
    RooRealVar* mean1 = ws->var("mean1");
    RooRealVar* mean2 = ws->var("mean2");
    RooRealVar* sigma1 = ws->var("sigma1");
    RooRealVar* sigma2 = ws->var("sigma2");
    RooAbsPdf * signal_pdf_1 = ws->pdf("sig1");
    RooAbsPdf * signal_pdf_2 = ws->pdf("sig2");
    RooRealVar* fsig_1 = ws->var("fsig_1");
    RooAbsPdf * total_signal_pdf = ws->pdf("twoGaussians");
    RooRealVar* chyb_prm = ws->var("chyb_prm");
    RooAbsPdf * backround_pdf = ws->pdf("bkg");
    RooRealVar* n_sig = ws->var("n_sig");
    RooRealVar* n_bkg = ws->var("n_bkg");
    RooAbsPdf * added_pdf = ws->pdf("sum");
    const char* output_name = "defit_using_ws.png";

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
    RooPlot* deframe = var_name->frame(Title("Fitting  #DeltaE of B^{#pm}"), Bins(300)) ;                          
    dataset_name->plotOn(deframe, Binning(300), DataError(RooAbsData::SumW2)) ;
    // sum.plotOn(deframe, LineColor(kBlue)	, LineStyle(kSolid)) ;
    // sum.paramOn(deframe,data,"", 2, Format("NEU"),AutoPrecision(1)), 0.7, 0.9, 0.9); //"NELU",  Prints the fitted parameter on the canvas
    added_pdf->plotOn(deframe,Components(*signal_pdf_1),LineColor(kGreen),LineStyle(kDashed)) ;
    added_pdf->plotOn(deframe,Components(*signal_pdf_2),LineColor(kBlack),LineStyle(kDashed)) ;
    added_pdf->plotOn(deframe,Components(*total_signal_pdf),LineColor(kRed),LineStyle(kDashed));
    added_pdf->plotOn(deframe,Components(*backround_pdf),LineColor(kMagenta),LineStyle(kDashed)) ;
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
    TLegend *legend1 = new TLegend(0.65,0.7,0.9,0.9);
    TLegendEntry *entry = legend1->AddEntry("signal_pdf_1","1st Gaussian pdf","l");
    entry->SetLineColor(kGreen);
    entry->SetLineStyle(kDashed);
    entry = legend1->AddEntry("signal_pdf_2","2nd Gaussian pdf","l");
    entry->SetLineColor(kBlack);
    entry->SetLineStyle(kDashed);
    entry = legend1->AddEntry("total_signal_pdf","Combined signal pdf","l");
    entry->SetLineColor(kRed);
    entry->SetLineStyle(kDashed);
    entry = legend1->AddEntry("backround_pdf","bkg-Chebyshev","l");
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
    frame3->GetXaxis()->SetTitle("#DeltaE (GeV)"); 
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

// void de_plot_style_after_calculating_sWeights(RooWorkspace* ws){
//     // get what we need out of the workspace for import_root_to_dataset
//     RooAbsPdf* sum = ws->pdf("sum");
//     RooRealVar* get_counter = ws->var("get_counter");
//     RooRealVar* n_sig = ws->var("n_sig");
//     RooRealVar* n_bkg = ws->var("n_bkg");
//     int counter = (int)get_counter->getVal();
//     RooDataSet* data = (RooDataSet*) ws->data("data_having_info_and_sWeights");


//     RooRealVar* var_name = ws->var("deltae");
//     RooDataSet* dataset_name = (RooDataSet*) ws->data("data_having_info");
//     RooRealVar* mean1 = ws->var("mean1");
//     RooRealVar* mean2 = ws->var("mean2");
//     RooRealVar* sigma1 = ws->var("sigma1");
//     RooRealVar* sigma2 = ws->var("sigma2");
//     RooGaussian signal_pdf_1 = ws->pdf("sig1");
//     RooGaussian signal_pdf_2 = ws->pdf("sig2");
//     RooRealVar* fsig_1 = ws->var("fsig_1");
//     RooAddPdf total_signal_pdf = ws->pdf("twoGaussians");
//     RooRealVar* chyb_prm = ws->var("chyb_prm");
//     RooChebychev backround_pdf = ws->pdf("bkg");
//     RooRealVar* n_sig = ws->var("n_sig");
//     RooRealVar* n_bkg = ws->var("n_bkg");
//     RooAddPdf* added_pdf = ws->pdf("sum");
//     const char* output_name = "defit_using_ws.png";

//     // // Integrate sig pdf(following some hypothesis) to get yield
//     // mbc.setRange("twoGaussians",5.27, 5.29);     //twoGaussians is my signal pdf    //5.27, 5.29 is the range we want to integrate
//     // RooAbsReal *integral_sig = twoGaussians.createIntegral(mbc,NormSet(mbc),Range("twoGaussians"));

//     // double  Nsig = integral_sig->getVal();
//     // cout<<"Signal :"<<Nsig*n_sig.getVal()<<endl;   // nsig is the signal Yield
//     // double Nsigerr = n_sig.getError()*integral_sig->getVal();  
//     // cout<<"Signal error = "<<Nsigerr<<endl;
//     // cout<<"Signal Area "<< setprecision(4)<<Nsig<<endl;
//     /*********************Start Plotting and showing outpouts*****************/

//     //Plotting fitted result
//     RooPlot* deframe = var_name->frame(Title("Fitting  #DeltaE of B^{#pm}"), Bins(300)) ;                          
//     dataset_name->plotOn(deframe, Binning(300), DataError(RooAbsData::SumW2)) ;
//     // sum.plotOn(deframe, LineColor(kBlue)	, LineStyle(kSolid)) ;
//     // sum.paramOn(deframe,data,"", 2, Format("NEU"),AutoPrecision(1)), 0.7, 0.9, 0.9); //"NELU",  Prints the fitted parameter on the canvas
//     added_pdf->plotOn(deframe,Components(signal_pdf_1),LineColor(kGreen),LineStyle(kDashed)) ;
//     added_pdf->plotOn(deframe,Components(signal_pdf_2),LineColor(kBlack),LineStyle(kDashed)) ;
//     added_pdf->plotOn(deframe,Components(total_signal_pdf),LineColor(kRed),LineStyle(kDashed));
//     added_pdf->plotOn(deframe,Components(backround_pdf),LineColor(kMagenta),LineStyle(kDashed)) ;
//     added_pdf->plotOn(deframe, LineColor(kBlue), LineStyle(kSolid)) ; // we need to write it last.
//                     // otherwise pull distribution will calculate error wrt last mentioned pdf inside plotOn function
//     added_pdf->paramOn(deframe,dataset_name,"", 2, "NEU", 0.1, 0.35, 0.9); //"NELU",  Prints the fitted parameter on the canvas

//     //Extract info. from fitting frame and showing
//     RooHist* hpull = deframe->pullHist() ;
//     RooPlot* frame3 = var_name->frame(Title("Pull Distribution")) ;
//     hpull->SetFillColor(1);
//     frame3->addPlotable(hpull,"X0B"); // "X0" is for errorbar; and "B" is for histograms
//     // frame3->addPlotable(hpull,"P");


//     TCanvas* c1 = new TCanvas("c1", "c1", 2550, 1500) ;
//     TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
//     pad1->Draw();             // Draw the upper pad: pad1
//     pad1->cd();  
//     deframe->Draw() ;
//     // Adding legend
//     TLegend *legend1 = new TLegend(0.65,0.7,0.9,0.9);
//     TLegendEntry *entry = legend1->AddEntry("signal_pdf_1","1st Gaussian pdf","l");
//     entry->SetLineColor(kGreen);
//     entry->SetLineStyle(kDashed);
//     entry = legend1->AddEntry("signal_pdf_2","2nd Gaussian pdf","l");
//     entry->SetLineColor(kBlack);
//     entry->SetLineStyle(kDashed);
//     entry = legend1->AddEntry("total_signal_pdf","Combined signal pdf","l");
//     entry->SetLineColor(kRed);
//     entry->SetLineStyle(kDashed);
//     entry = legend1->AddEntry("backround_pdf","bkg-Chebyshev","l");
//     entry->SetLineColor(kMagenta);
//     entry->SetLineStyle(kDashed);
//     entry = legend1->AddEntry("added_pdf","Fitted pdf","l");
//     entry->SetLineColor(kBlue);
//     entry->SetLineStyle(kSolid);
//     legend1->Draw();

//     // Printing error in "Signal yield calculation" // Problem is resolved using extra "E" inside paramOn function 
//     // std::string sig_err_in_str = std::to_string(n_sig.getError());
//     // TLatex* sig_err = new TLatex();
//     // sig_err->SetTextSize(0.034);
//     // sig_err->SetTextAlign(12);  //centered aligned
//     // sig_err->DrawLatex(5.2392, 1450, "#pm");
//     // sig_err->DrawLatex(5.24, 1450, sig_err_in_str); // here we need to give "TEXT" ONLY

//     // // Adding arrow at (+-)3*sigma of signal pdf
//     // double weightedMean = mean1.getVal()*fsig_1.getVal() + mean2.getVal()*(1.0-fsig_1.getVal());
//     // double weightedSigma = std::sqrt( pow(sigma1.getVal(),2)*pow(fsig_1.getVal(),2)  +  pow(sigma2.getVal(),2)*pow((1-fsig_1.getVal()),2));
//     // double mbc_min_fit = weightedMean -3*weightedSigma; 
//     // double mbc_max_fit = weightedMean +3*weightedSigma;
//     // double arrowHeight = n_sig.getVal()*0.03;
//     // cout<<"Mbc min :"<< mbc_min_fit <<endl;
//     // cout<<"Mbc max :"<< mbc_max_fit <<endl;
//     // cout<<"Arrow Height :"<< arrowHeight <<endl;

//     // TArrow *arr1 = new TArrow(-0.01, 100, -0.01, 400, 0.01,"|>");
//     // arr1->SetLineWidth(4);
//     // arr1->SetLineColor(2);
//     // arr1->SetFillStyle(3008);
//     // arr1->DrawArrow(mbc_min_fit, arrowHeight*0.1, mbc_min_fit, arrowHeight, 0.02,"<"); 
//     // arr1->DrawArrow(mbc_max_fit, arrowHeight*0.1, mbc_max_fit, arrowHeight, 0.02,"<");

//     c1->cd();          // Go back to the main canvas before defining pad2
//     TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
//     pad2->Draw();
//     pad2->cd(); 
//     frame3->Draw() ;
//     frame3->SetLineStyle(9);
//     frame3->GetYaxis()->SetNdivisions(505);
//     frame3->GetYaxis()->SetTitle("#sqrt{#chi^{2}}"); 
//     frame3->GetXaxis()->SetTitle("#DeltaE (GeV)"); 
//     frame3->GetXaxis()->SetTitleSize(0.13);
//     frame3->GetYaxis()->SetTitleSize(0.15);
//     frame3->GetXaxis()->SetLabelSize(0.120);
//     frame3->GetYaxis()->SetLabelSize(0.120); 
//     frame3->GetXaxis()->SetTitleOffset(0.90);      
//     frame3->GetYaxis()->SetTitleOffset(0.22);       
//     frame3->GetYaxis()->SetRangeUser(-10.0, 10.0);       
//     frame3->GetYaxis()->SetLimits(-10.0, 10.0);       
//     frame3->GetXaxis()->SetNdivisions(505);
//     frame3->GetYaxis()->CenterTitle(true);
//     frame3->GetXaxis()->CenterTitle(true);
//     frame3->Draw("AXISSAME");

//     // cout<<"Total number of events which are used to fitting are : "<<counter<<endl;
//     cout<<"chisq of the fit is :"<<deframe->chiSquare()<<endl;//chi-square of the fit
//     cout<<"chi-square/ndof :"<<deframe->chiSquare(7)<<endl;// Chi^2/(the number of degrees of freedom)
//     cout<<"Error in calculation of signal yield : "<<n_sig->getError()<<endl;

//     const char* output_filename = output_name;
//     c1->Print(output_filename);
//     return 0;
// }



// void plot_control_var_using_sWeights(RooWorkspace* ws){












//     // create weighted data sets
//     //RooDataSet(const char* name, const char* title, RooDataSet* data, const RooArgSet& vars, const RooFormulaVar& cutVar, const char* wgtVarName = 0)
//     RooDataSet * dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"n_sig_sw") ;
//     RooDataSet * dataw_bg  = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"n_bkg_sw") ;
//     /////////////////////////////////////////////////////////////////////////////???????????????////////////////////////////////////
//     // How signal and bkg are seperately defined here???????????????????????????here how r we getting different variables from Roodataset in two different line?
//     // This is very confusing but the punch line, here, is:
//     // /*The SPlot class adds (has already added) a new variable that has the name of the corresponding "{yield + "_sw"}".**/
    

//     //Defining frame to plot sWeighted pdf
//     cout << "Making splots of the signal and background" <<endl;
//     RooPlot* frame_sig_rep = bdt.frame(50) ;  // Defining frame(from plot of bdt output) to plot signal distribution
//     frame_sig_rep->SetTitle("sPlot for the Signal Probability distribution(signal)");
//     dataw_sig->plotOn(frame_sig_rep, RooFit::DataError(RooAbsData::SumW2) ) ;//???????????????????????????

//     RooPlot* frame_bg_rep = bdt.frame(50) ; // Defining frame(from plot of bdt output) to plot bkg distribution
//     frame_bg_rep->SetTitle("sPlot for the Signal Probability distribtuion(background)");
//     dataw_bg->plotOn(frame_bg_rep, RooFit::DataError(RooAbsData::SumW2) ) ;//???????????????????????????

//     // plot PDFs onto the sPlots; if the sPlot correctly reweights the data then the signal / background only shapes should
//     // provide an adequate description of the data.
//     TCanvas *c2 = new TCanvas("c2" ,"c2");
//     c2->cd();
//     frame_sig_rep->Draw();

//     TCanvas *c3 = new TCanvas("c3","c3");
//     c3->cd();
//     frame_bg_rep->Draw();

//     c2->SaveAs("mc_splot_sig.png");
//     c3->SaveAs("mc_splot_bkg.png");


//     // Plot sWeighted distriution of control variable(MVA output)
//     // in histogram along with original MVA output
//     // to compare sPlot output distrubution and MVA output distribution
//     bdt.setBins(50);
//     TH1F *h1 = (TH1F*)dataw_sig->createHistogram("signal" ,bdt);//???????????????????????????
//     double scale = 1./h1->Integral();
//     h1->Scale(scale);
//     h1->SetLineColor(kBlack);
//     h1->SetLineWidth(5);

//     TH1F *h1_b = (TH1F*)dataw_bg->createHistogram("background" ,bdt);//???????????????????????????
//     double scaleb = 1./h1_b->Integral();
//     h1_b->Scale(scaleb);
//     h1_b->SetLineColor(kRed);
//     h1_b->SetLineWidth(2);


//     ///////////////////////////plotting MVA output in histogram///////////////////////////////////////////////////////////
//     // TH1F *control_var_histogram_for_signal = new TH1F("control_var_histogram_for_signal","",50,0,1);
//     // TH1F *control_var_histogram_for_bkg = new TH1F("control_var_histogram_for_bkg","",50,0,1);
//     // for(int l=0;l<nevt3;l++){
//     //     chain->GetEntry(l);
//     //     if(sig==1) control_var_histogram_for_signal->Fill(sig_prob);//SIG
//     //     if(sig==0) control_var_histogram_for_bkg->Fill(sig_prob);//bkg
//     // }
//     control_var_histogram_for_signal->SetLineColor(kMagenta);
//     control_var_histogram_for_signal->SetLineWidth(2);
//     Double_t scale2 = 1./control_var_histogram_for_signal->Integral();
//     control_var_histogram_for_signal->Scale(scale2);

//     control_var_histogram_for_bkg->SetLineColor(kBlue);
//     control_var_histogram_for_bkg->SetLineWidth(2);
//     Double_t scale4 = 1./control_var_histogram_for_bkg->Integral();
//     control_var_histogram_for_bkg->Scale(scale4);
//     ////////////////////////////////plotting MVA output in histogram******end//////////////////////////////////////////////


//     TCanvas *c4 = new TCanvas("c4" ,"c4");
//     c4->cd();

//     h1->SetStats(kFALSE);
//     h1_b->SetStats(kFALSE);
//     control_var_histogram_for_signal->SetStats(kFALSE);
//     control_var_histogram_for_bkg->SetStats(kFALSE);

//     h1->Draw();
//     h1_b->Draw("SAME");
//     control_var_histogram_for_signal->Draw("SAME");//SIG
//     control_var_histogram_for_bkg->Draw("SAME");//bkg
    
//     auto legend = new TLegend( 0.25,0.7,0.74,0.9);
//     legend->AddEntry(control_var_histogram_for_signal,"Signal Probability predicted by MVA for signal","l");
//     legend->AddEntry(control_var_histogram_for_bkg,"Signal Probability predicted by MVA for background","l");
//     legend->AddEntry(h1,"Signal Probability predicted by sPlot for signal","l");
//     legend->AddEntry(h1_b,"Signal Probability predicted by sPlot for background","l");
//     legend->Draw();
//     c4->SaveAs("splot_sig_bkg.png");
// }





void defit_using_workspace(){
    // Create a new workspace to manage the project.
    RooWorkspace* wspace = new RooWorkspace("myWS");

    define_var_dataset(wspace);
    import_root_to_dataset(wspace);
    // prefit_plot(wspace);
    define_pdf(wspace);
    de_fit(wspace);
    de_plot_style(wspace);
    // de_fit_and_get_sWeights(wspace);
    // de_plot_style_after_calculating_sWeights(wspace);
    // plot_control_var_using_sWeights(wspace);
    // cleanup
    delete wspace;
}