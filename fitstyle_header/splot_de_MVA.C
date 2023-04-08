#include "de_fit_plot_style.h"
// #include "prefit_hist_draw.h"
// #include "my_splot.h"

void splot_de_MVA(){
    /*******************Fit Variables***********************************/
    RooRealVar deltae("deltae","#DeltaE (GeV)", -0.1, 0.1);
    RooRealVar bdt("bdt","FBDT output", 0., 1.);

    /**defining DATAFRAME{unbinned histogram}(FROM FIT VARIABLE) TO FIT AND PLOT**/
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(deltae,bdt));
    /*******************Input root file**********************************/
    TChain* chain=new TChain();
    chain->Add("/home/belle2/ssana/MC15ri/cs/test/signal_scaled/test.root/tree");

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
    TH1F *control_var_histogram_for_signal = new TH1F("control_var_histogram_for_signal","",50,0,1);
    TH1F *control_var_histogram_for_bkg = new TH1F("control_var_histogram_for_bkg","",50,0,1);
    for(int l=0;l<nevt3;l++) {
        chain->GetEntry(l);
        if(sig==1) control_var_histogram_for_signal->Fill(sig_prob);//SIG
        if(sig==0) control_var_histogram_for_bkg->Fill(sig_prob);//bkg
        deltae.setVal(de3);
        bdt.setVal(sig_prob);
        if(mbc3>5.23 && mbc3 < 5.29 &&
        de3 < 0.1 && de3 > -0.1 &&
        md03>1.85 && md03<1.88 && 
        kid3 > 0.6
        ){
            data->add(RooArgSet(deltae,bdt));
            counter++;
        }
    }

    /*****************************Delta E Fit***********************/
    // --- Build Signal PDF ---
    RooRealVar mean1("mean1","mean of Gaussian-1",-0.001,-0.02,0.02);
    RooRealVar mean2("mean2","mean of Gaussian-2",0.002,-0.02,0.02);
    RooRealVar sigma1("sigma1","sigma of Gaussian-1",0., 0., 1);	
    RooRealVar sigma2("sigma2","sigma of Gaussian-2",0.01191,0.00000001,1);

    RooGaussian sig1("sig1","Gaussian-1",deltae,mean1,sigma1);  
    RooGaussian sig2("sig2","Gaussian-2",deltae,mean2,sigma2);

    RooRealVar fsig_1("frac_gaussians", "signal fraction", 0.5,0.,1.);
    RooAddPdf twoGaussians("twoGaussians", "sum of two Gaussians ",RooArgList(sig1, sig2), RooArgList(fsig_1));

    // --- Build Argus background PDF ---
    RooRealVar b1("Chbyshv-prm", "Chbyshv-prm", -0.062, -10., 10.);
    RooChebychev bkg("bkg","Background",deltae,RooArgSet(b1)) ;

    //Initialization of parameter before adding two pdf
    // cout<<"Total number of events which are used to fitting are : "<<counter<<endl;
    Double_t event_count = counter; 
    Double_t signal_count = counter*0.2;
    Double_t back_count = counter*0.8;
    RooRealVar n_sig("n_sig", "n_sig", signal_count, 0., event_count);
    RooRealVar n_bkg("n_bkg", "n_bkg", back_count, 0., event_count);
    // RooRealVar n_sig("n_sig", "n_sig", 145795, 14500, 146600);
    // RooRealVar n_bkg("n_bkg", "n_bkg", 157295, 156400, 158100);
    RooAddPdf sum("sum","sum",RooArgList(twoGaussians,bkg),RooArgList(n_sig, n_bkg));//adding two pdf
    RooFitResult * my_fit_result = (RooFitResult*)sum.fitTo(*data," ", RooFit::Save());   //this will fit only DELTA_E as we mentioned "the fitting variable(deltae)" during defining each elementary pdf
    // RooMsgService::instance().setSilentMode(true); //to turn off all informational 
         //and warning messages that would normally be printed to the console 
         //during a RooFit calculation, while still allowing error messages to be
         //printed to the console. This can be useful in cases where the user 
         //wants to suppress output messages during a batch job or automated 
         //script, or when the user is only interested in the final result of a calculation.    
    cout<<"Total number of events which are used to fitting are : "<<counter<<endl;
    cout<<"Printing fitted result.........................."<<endl;
    my_fit_result->Print("v");
    /****************************FIT COMPLETE*************************************/

    //Plotting fitted result
    const char* output_filename = "defit_using_headerfile.png";
    de_fit_plot_style(&deltae, data, &mean1, &mean2, &sigma1, &sigma2,
    sig1, sig2, &fsig_1, twoGaussians,
    &b1, bkg, &n_sig, &n_bkg,
    &sum, output_filename);


    /*****************************Splot***********************/
    /*****************************Splot***********************/
    // splot_define_n_plot(&deltae, &bdt, data,
    // control_var_histogram_for_signal, control_var_histogram_for_bkg, // chain,
    // &n_sig, &n_bkg,
    // &sum);



    // Create the sPlot data from the DATA of DeltaE fitting
    RooStats::SPlot* sDataX = new RooStats::SPlot("sData","An SPlot", *data, &sum, RooArgList(n_sig, n_bkg) );

    // Check that our weights have the desired properties, which is vlaue of Yield and sWeights are nearly same.
    cout << "Check that our weights have the desired properties, which is vlaue of Yield and sWeights are nearly same." <<endl;
    cout << "Yield of signal is " << n_sig.getVal() << " From sWeights it is " << sDataX->GetYieldFromSWeight("n_sig") <<endl;
    cout << "Yield of background is " << n_bkg.getVal() << "  From sWeights it is " << sDataX->GetYieldFromSWeight("n_bkg") <<endl;
   
    // Printing different sWeight
    for(Int_t i=0; i < 20; i++){
        cout << "signal Weight   " << sDataX->GetSWeight(i,"n_sig") 
        << " bkg Weight   " << sDataX->GetSWeight(i,"n_bkg") 
        << "  Total Weight   " << sDataX->GetSumOfEventSWeight(i) 
        << endl;
    }

    // create weighted data sets
    //RooDataSet(const char* name, const char* title, RooDataSet* data, const RooArgSet& vars, const RooFormulaVar& cutVar, const char* wgtVarName = 0)
    RooDataSet * dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"n_sig_sw") ;
    RooDataSet * dataw_bg  = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"n_bkg_sw") ;
    /////////////////////////////////////////////////////////////////////////////???????????????////////////////////////////////////
    // How signal and bkg are seperately defined here???????????????????????????here how r we getting different variables from Roodataset in two different line?
    // This is very confusing but the punch line, here, is:
    // /*The SPlot class adds (has already added) a new variable that has the name of the corresponding "{yield + "_sw"}".**/
    

    //Defining frame to plot sWeighted pdf
    cout << "Making splots of the signal and background" <<endl;
    RooPlot* frame_sig_rep = bdt.frame(50) ;  // Defining frame(from plot of bdt output) to plot signal distribution
    frame_sig_rep->SetTitle("sPlot for the Signal Probability distribution(signal)");
    dataw_sig->plotOn(frame_sig_rep, RooFit::DataError(RooAbsData::SumW2) ) ;//???????????????????????????

    RooPlot* frame_bg_rep = bdt.frame(50) ; // Defining frame(from plot of bdt output) to plot bkg distribution
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
    bdt.setBins(50);
    TH1F *h1 = (TH1F*)dataw_sig->createHistogram("signal" ,bdt);//???????????????????????????
    double scale = 1./h1->Integral();
    h1->Scale(scale);
    h1->SetLineColor(kBlack);
    h1->SetLineWidth(5);

    TH1F *h1_b = (TH1F*)dataw_bg->createHistogram("background" ,bdt);//???????????????????????????
    double scaleb = 1./h1_b->Integral();
    h1_b->Scale(scaleb);
    h1_b->SetLineColor(kRed);
    h1_b->SetLineWidth(2);


    ///////////////////////////plotting MVA output in histogram///////////////////////////////////////////////////////////
    // TH1F *control_var_histogram_for_signal = new TH1F("control_var_histogram_for_signal","",50,0,1);
    // TH1F *control_var_histogram_for_bkg = new TH1F("control_var_histogram_for_bkg","",50,0,1);
    // for(int l=0;l<nevt3;l++){
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

    h1->SetStats(kFALSE);
    h1_b->SetStats(kFALSE);
    control_var_histogram_for_signal->SetStats(kFALSE);
    control_var_histogram_for_bkg->SetStats(kFALSE);

    h1->Draw();
    h1_b->Draw("SAME");
    control_var_histogram_for_signal->Draw("SAME");//SIG
    control_var_histogram_for_bkg->Draw("SAME");//bkg
    
    auto legend = new TLegend( 0.25,0.7,0.74,0.9);
    legend->AddEntry(control_var_histogram_for_signal,"Signal Probability predicted by MVA for signal","l");
    legend->AddEntry(control_var_histogram_for_bkg,"Signal Probability predicted by MVA for background","l");
    legend->AddEntry(h1,"Signal Probability predicted by sPlot for signal","l");
    legend->AddEntry(h1_b,"Signal Probability predicted by sPlot for background","l");
    legend->Draw();
    c4->SaveAs("splot_sig_bkg.png");
}