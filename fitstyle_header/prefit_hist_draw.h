#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
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
#include <cmath>
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooSuperCategory.h"
#include "RooSimultaneous.h"
#include "RooNLLVar.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "RooArgusBG.h"
using namespace RooFit ;
using namespace std;
void prefit_hist_draw(RooRealVar* var_name, RooDataSet* dataset_name){
  RooDataHist* binDataSet = new RooDataHist("binDataSet", "binDataSet", *var_name, *dataset_name);
  RooPlot* xframe1 = var_name->frame(	Title("prefit")	, Bins(200));
  binDataSet->plotOn(	xframe1	, Binning(200)		, DataError(RooAbsData::SumW2));
  xframe1->Draw();

  TCanvas *prefit_canvas = new TCanvas("prefit_canvas", "", 1500, 1500);  
  prefit_canvas -> Print("prefit_histo.png");
}