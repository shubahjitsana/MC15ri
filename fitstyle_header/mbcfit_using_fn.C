#include "mbc_fit_plot_style.h"
// #include "prefit_hist_draw.h"

void mbcfit_using_fn(){
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

  // prefit_hist_draw(&mbc,data);

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
  const char* output_filename = "mbcfit_using_headerfile_2.png";
  mbc_fit_plot_style(&mbc, data, &mean1, &mean2, &sigma1, &sigma2,
  sig1, sig2, &fsig_1, twoGaussians,
  &argpar, bkg, &n_sig, &n_bkg,
  &sum, output_filename);
}