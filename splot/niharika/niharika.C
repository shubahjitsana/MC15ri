using namespace RooFit;

void splot_q()
{
  gStyle->SetOptFit(kTRUE);
  Color_t white=10;
  gStyle->SetCanvasColor(white);
  gSystem->Load("libRooFit");
  
  RooRealVar *mbc = new RooRealVar("mbc","M_{bc}",5.22,5.289);
  RooRealVar *mjpsi = new RooRealVar("mjpsi","q^{2}",0,15);

  TChain chain("h3");
  chain.Add("../../neural_network/stream1/nnout_charged.root");
  chain.Add("../../neural_network/stream1/nnout_mixed.root");
  
  Float_t x_mbc,x_delte,x_bchide8,x_nb,x_k1ma1,x_pi1ma1,x_pi2ma1,x_mu1ma1,x_mu2ma1,x_k1ma2, x_mbp,x_mjpsi,x_mbm,x_bchide4;
  chain.SetBranchAddress("mbc",&x_mbc);
  chain.SetBranchAddress("qsqr",&x_mjpsi);
  chain.SetBranchAddress("bchide8",&x_bchide8);
  chain.SetBranchAddress("bchide4",&x_bchide4);
  chain.SetBranchAddress("delte",&x_delte);
  //  chain.SetBranchAddress("NB",&x_nb);
  chain.SetBranchAddress("mode_bp",&x_mbp);
  chain.SetBranchAddress("mode_bm",&x_mbm);
  RooDataSet *data=new RooDataSet("data","", RooArgSet(*mbc,*mjpsi));
  for(int i=0;i<chain.GetEntries();i++)
    {
      chain.GetEntry(i);
      if(x_delte>-0.025&&x_delte<0.025&&x_bchide4==100&&x_mbc>5.22&&x_mbc<5.289)
	{
	  mbc->setVal(x_mbc);
	  mjpsi->setVal(x_mjpsi);
	  data->add(RooArgSet(*mbc,*mjpsi));
	  }
      }

  //mbc pdf
  RooRealVar *meanm=new RooRealVar("meanm","Mean of ist gaussian",5.278999,5.275,5.285);
  RooRealVar *sigmam=new RooRealVar("sigmam","sigma",0.0025982,0.001,0.04);
  RooGaussian *mgau=new RooGaussian("mgau"," ",*mbc,*meanm,*sigmam);
  
  RooRealVar *Margpar=new RooRealVar("Margpar","argus shape",-14.73,-100,-1);
  RooRealVar *par=new RooRealVar("par","argus shape",5.29,5.28,5.3);
  RooArgusBG *argus=new RooArgusBG("argus","ArgusPDF",*mbc, *par,*Margpar);
  //  RooArgusBG *argus=new RooArgusBG("argus","ArgusPDF",*mbc, RooConst(5.289),*Margpar);

  RooRealVar *mbsig=new RooRealVar("mbsig","",66632,0,10000000);
  RooRealVar *mBKG=new RooRealVar("mBKG","",9081,0,1000000);
  
  RooAddPdf *total=new RooAddPdf("total","total",RooArgList(*mgau,*argus),RooArgList(*mbsig,*mBKG));
  //RooAddPdf *total=new RooAddPdf("total","total",RooArgList(*argus),RooArgList(*mBKG));

  total->fitTo(*data);

  TCanvas* c2 = new TCanvas("c2","c2",1) ;
  RooPlot *xframe = data->plotOn(mbc->frame(200));
  total->plotOn(xframe);
  total->paramOn(xframe,data);
  //total->plotOn(xframe, Components(RooArgList(*mgau)),LineStyle(kDashed),LineColor(kRed));
  total->plotOn(xframe, Components(RooArgList(*argus)),LineStyle(kDashed),LineColor(kGreen));
  xframe->Draw();

 //splot
  
   std::cout << " Delta E fit is done here.. Now Splot business start "
            << std::endl;
  meanm->setConstant(true);
  sigmam->setConstant(true);
  mbsig->setConstant(true);
  mBKG->setConstant(true);
 
  RooStats::SPlot* sData = new RooStats::SPlot("sWeightedData","sWeightedData", *data, total, RooArgList(*mbsig,*mBKG) );
  std::cout << " Check Weights  \n";

  std::cout << "\n Yield of signal is "
            << mbsig->getVal() << ". From sWeight it is "
            << sData->GetYieldFromSWeight("mbsig") << endl;

  std::cout << "\n Yield of BKG is "
            << mBKG->getVal() << ". From sWeight it is "
            << sData->GetYieldFromSWeight("mBKG") << endl;

  RooDataSet * dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"mbsig_sw") ;
  RooDataSet * dataw_bg  = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"mBKG_sw") ;
  std::cout << "Making splots of the signal and background" << std::endl;

  RooPlot* frame_sig_y = mjpsi->frame(80);
  frame_sig_y->SetTitle("sPlot for the signal y distribution");
  dataw_sig->plotOn(frame_sig_y, RooFit::DataError(RooAbsData::SumW2) ) ;

  RooPlot* frame_bg_y = mjpsi->frame(80);
  frame_bg_y->SetTitle("sPlot for the background y distribtuion");
  dataw_bg->plotOn(frame_bg_y, RooFit::DataError(RooAbsData::SumW2) ) ;

   
  TCanvas *can =new TCanvas("can","canvas",1);
  frame_sig_y->Draw();

  TCanvas *can_bkg =new TCanvas("can_bkg","canvas",1);
  frame_bg_y->Draw();
  
  
}
