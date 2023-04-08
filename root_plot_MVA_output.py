#!/usr/bin/env python
# -*- coding: utf-8 -*-
# takes commadnd as
# basf2 root_plot_MVA_output.py file1 file2

import ROOT
import sys
import numpy as np
ROOT.TH1.AddDirectory(ROOT.kFALSE)

##############taking input from command line####################
chain = ROOT.TChain("tree")
for i in sys.argv[1:]:
    chain.Add(f"{i}")
    print(f"Given input files are {i}")


deltae = ROOT.vector('double')()
mbc = ROOT.vector('double')()
md0 = ROOT.vector('double')()
# kid = ROOT.vector('double')()
signal_prob = ROOT.vector('double')()
issig = ROOT.vector('double')()
iscontinuum = ROOT.vector('double')()

D_s_InvM = ROOT.vector('double')()
D_10D_md = ROOT.vector('double')()
InvM1stand2ndpi = ROOT.vector('double')()
InvM2ndand3rdpi = ROOT.vector('double')()
InvMD0and1stpi = ROOT.vector('double')()
InvMD0and2ndpi = ROOT.vector('double')()
InvMD0and3rdpi = ROOT.vector('double')()
DstrminusroeD_md = ROOT.vector('double')()


chain.SetBranchAddress("deltaE", deltae)
chain.SetBranchAddress("Mbc", mbc)
chain.SetBranchAddress("D0_bar_InvM", md0)
# chain.SetBranchAddress("Kp_PID_bin_kaon", kid)
chain.SetBranchAddress("SigProb", signal_prob)
chain.SetBranchAddress("isSignal", issig)
chain.SetBranchAddress("isContinuumEvent", iscontinuum)

chain.SetBranchAddress("D_s_InvM", D_s_InvM)
chain.SetBranchAddress("D_10D_md", D_10D_md)
chain.SetBranchAddress("InvM1stand2ndpi", InvM1stand2ndpi)
chain.SetBranchAddress("InvM2ndand3rdpi", InvM2ndand3rdpi)
chain.SetBranchAddress("InvMD0and1stpi", InvMD0and1stpi)
chain.SetBranchAddress("InvMD0and2ndpi", InvMD0and2ndpi)
chain.SetBranchAddress("InvMD0and3rdpi", InvMD0and3rdpi)
chain.SetBranchAddress("DstrminusroeD_md", DstrminusroeD_md)

histogram_continuum = ROOT.TH1F("Sig_prob", "Distribution of MVA Output",100, 0, 100)
histogram_non_continuum = ROOT.TH1F("Sig_prob", "Distribution of MVA Output",100, 0, 100)
for iEvent in range(chain.GetEntries()):
    chain.GetEntry(iEvent)
    if (
        (deltae[0] > -0.1 and deltae[0] < 0.1) and 
        (mbc[0] > 5.23 and mbc[0] < 5.29) and 
        (md0[0] > 1.85 and md0[0] < 1.88)
        # (kid[0] > 0.6) and 
        # (D_s_InvM[0] < 1.958 or D_s_InvM[0] > 1.978) and 
        # (D_10D_md[0] < 0.496 or D_10D_md[0] > 0.628) and 
        # (InvM1stand2ndpi[0] < 0.489 or InvM1stand2ndpi[0] > 0.506) and 
        # (InvM2ndand3rdpi[0] < 0.49 or InvM2ndand3rdpi[0] > 0.505) and 
        # (InvMD0and1stpi[0] < 2.0085 or InvMD0and1stpi[0] > 2.0122) and 
        # (InvMD0and2ndpi[0] < 2.0087 or InvMD0and2ndpi[0] > 2.0118) and 
        # (InvMD0and3rdpi[0] < 2.009 or InvMD0and3rdpi[0] > 2.0116) and 
        # (DstrminusroeD_md[0] < 0.14402 or DstrminusroeD_md[0] > 0.14682)
    ):
        if iscontinuum[0] == 1: histogram_continuum.Fill(signal_prob[0])
        else: histogram_non_continuum.Fill(signal_prob[0])

#################### Start Drawing ##############################
colors = [ROOT.kRed, ROOT.kViolet+10, ROOT.kGreen, ROOT.kBlack]
c1 = ROOT.TCanvas("c", "c", 3600, 2800)
legend = ROOT.TLegend(0.35,0.7,0.55,0.9)
c1.cd()

if (histogram_continuum.GetMaximum() > histogram_non_continuum.GetMaximum()):
    histogram_continuum.GetYaxis().SetRangeUser(0,(200 + histogram_continuum.GetMaximum()))
else: histogram_continuum.GetYaxis().SetRangeUser(0,(200 + histogram_non_continuum.GetMaximum()))

histogram_continuum.GetXaxis().SetTitle("Signal Probability")
histogram_continuum.GetXaxis().SetTitleSize(0.035)
histogram_continuum.GetXaxis().SetLabelSize(0.02)
histogram_continuum.GetXaxis().CenterTitle(ROOT.kTRUE)
histogram_continuum.GetYaxis().SetTitle("Events")
histogram_continuum.GetYaxis().SetTitleSize(0.035)
histogram_continuum.GetYaxis().SetTitleOffset(0.95)
histogram_continuum.GetYaxis().SetLabelSize(0.02)
histogram_continuum.GetYaxis().CenterTitle(ROOT.kTRUE)

histogram_continuum.SetStats(ROOT.kFALSE)
histogram_non_continuum.SetStats(ROOT.kFALSE)

histogram_continuum.SetLineColor(ROOT.kRed)
histogram_non_continuum.SetLineColor(ROOT.kViolet+10)

histogram_continuum.Draw()
histogram_non_continuum.Draw("SAME")

legend.AddEntry(histogram[0], "Continuum Event", "l")
legend.AddEntry(histogram[1], "Noncontinuum Event", "l")
legend.Draw()

c1.Update()
c1.SaveAs("plot_from_my_anyfile_anyvariable_bin_range/mva_output.png")