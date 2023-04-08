#include "de_fit_plot_style.h"
// #include "prefit_hist_draw.h"

void defit_using_fn(){
    /*******************Fit Variables***********************************/
    RooRealVar deltae("deltae","#DeltaE (GeV)", -0.1, 0.1);
    /**defining DATAFRAME{unbinned histogram}(FROM FIT VARIABLE) TO FIT AND PLOT**/
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(deltae));
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
    for(int l=0;l<nevt3;l++) {
        chain->GetEntry(l);
        deltae.setVal(de3);
        if(mbc3>5.23 && mbc3 < 5.29 &&
         de3 < 0.1 && de3 > -0.1 &&
         md03>1.85 && md03<1.88 &&
         kid3 > 0.6 && 
         sig_prob > 0.22
         ){
            data->add(RooArgSet(deltae));
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
    RooRealVar n_sig("n_sig", "n_sig", signal_count, 0., event_count);//52000
    RooRealVar n_bkg("n_bkg", "n_bkg", back_count, 0., event_count);//95000
    RooAddPdf sum("sum","sum",RooArgList(twoGaussians,bkg),RooArgList(n_sig, n_bkg));//adding two pdf
    sum.fitTo(*data);
    /****************************FIT COMPLETE*************************************/

    cout<<"Total number of events which are used to fitting are : "<<counter<<endl;
    //Plotting fitted result
    const char* output_filename = "defit_using_headerfile.png";
    de_fit_plot_style(&deltae, data, &mean1, &mean2, &sigma1, &sigma2,
    sig1, sig2, &fsig_1, twoGaussians,
    &b1, bkg, &n_sig, &n_bkg,
    &sum, output_filename);
}