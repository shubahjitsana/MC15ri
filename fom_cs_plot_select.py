#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import sys
import numpy as np

print(f"Number of input files are {len(sys.argv)-1}")
print(f"{sys.argv[0]}")
for i in sys.argv[1:]:
    input_path = i.split("/all.root")
    fom_txt_output_path = input_path[0]
    fom_txt_output_filename = "fom_mbcde_kid_md.txt"
    fom_file = open(f"{fom_txt_output_path}/{fom_txt_output_filename}", 'w')    # defining exporting file

    fom_detail_report_filename = "fom_detail_report_mbcde_kid_md.txt"
    fom_detail_report = open(f"{fom_txt_output_path}/{fom_detail_report_filename}", 'w')    # defining exporting file

    input_file = ROOT.TFile.Open(f"{i}")    # Loading root file
    input_tree = input_file.Get('tree')
    total_event_in_charged = input_tree.GetEntries()

    # getting total number of continuum and non_continuum events
    # we need this to calculate cont_bkg rejection ratio and sig_efficiency
    total_continuum_bkg_event = 0
    total_non_continuum_sig_event = 0
    for iEvent in range(total_event_in_charged):
        input_tree.GetEntry(iEvent)
        de3 = getattr(input_tree, 'deltaE')
        mbc3 = getattr(input_tree, 'Mbc')
        md03 = getattr(input_tree, 'D0_bar_InvM')
        kid3 = getattr(input_tree, 'Kp_PID_bin_kaon')
        iscontinuum = getattr(input_tree, 'isContinuumEvent')
        if ((md03>1.85) and (md03<1.88) and (mbc3>5.27) and (mbc3 < 5.29) and (de3 < 0.02) and (de3 > -0.02) and (kid3 > 0.6)): 
            if iscontinuum == 1:total_continuum_bkg_event +=1
            else: total_non_continuum_sig_event += 1

    # fom analysis
    cs_prob_cut_array = np.array([0])
    cs_prob_fom_array = np.array([0])
    total_continuum_bkg_event_passing_cut = np.array([0])
    total_non_continuum_sig_event_passing_cut = np.array([0])
    # calculating FOM
    for cs_prob_cut in np.arange(0.01, 1, 0.01):
        sig_passing_cut = 0
        bkg_passing_cut = 0
        for iEvent in range(total_event_in_charged):
            input_tree.GetEntry(iEvent)
            de3 = getattr(input_tree, 'deltaE')
            mbc3 = getattr(input_tree, 'Mbc')
            md03 = getattr(input_tree, 'D0_bar_InvM')
            kid3 = getattr(input_tree, 'Kp_PID_bin_kaon')
            cs_prob = getattr(input_tree, 'ContProb')
            issig = getattr(input_tree, 'isSignal')
            iscontinuum = getattr(input_tree, 'isContinuumEvent')
            # if cs_prob < cs_prob_cut:
            if ((md03>1.85) and (md03<1.88) and (mbc3>5.27) and (mbc3 < 5.29) and (de3 < 0.02) and (de3 > -0.02) and (kid3 > 0.6) and (cs_prob < cs_prob_cut)):  
                # bkg += 1 if iscontinuum == 1 else sig += 1  # can't use arithmetic operator
                if iscontinuum == 1: bkg_passing_cut += 1
                else: sig_passing_cut +=1
        cs_prob_fom = sig_passing_cut/np.sqrt(sig_passing_cut + bkg_passing_cut)
        cs_prob_cut_array = np.append(cs_prob_cut_array, cs_prob_cut)
        cs_prob_fom_array = np.append(cs_prob_fom_array, cs_prob_fom)

        total_continuum_bkg_event_passing_cut = np.append(total_continuum_bkg_event_passing_cut, bkg_passing_cut)
        total_non_continuum_sig_event_passing_cut = np.append(total_non_continuum_sig_event_passing_cut, sig_passing_cut)

        fom_file.write(repr(cs_prob_cut) + '\t' + repr(cs_prob_fom) + '\n')
        fom_detail_report.write(repr(cs_prob_cut) + '\t' + repr(sig_passing_cut) + '\t' + repr(bkg_passing_cut) + '\t' + repr(cs_prob_fom) + '\n')
    fom_file.close()
    fom_detail_report.close()

    # store value in np array and plot those in Tgraph(need to flatten before plot)
    # Get max value from TGraph(also can get from np but need to access coordinate in CANVAS)
    # calculate bkg rejection and sig efficiency for MAX FOM value

    # Calculate cont_bkg rejection ratio and sig_efficiency
    fom_max = cs_prob_fom_array.max()
    fom_max_index = cs_prob_fom_array.argmax()
    # cs_prob_cut_for_fom_max = 0.01 + 0.01 * fom_max_index
    cs_prob_cut_for_fom_max = 0.01 * fom_max_index  # As cut 0.1 has been saved in 1st index
    cont_rejection_percentage = 100 * (total_continuum_bkg_event - total_continuum_bkg_event_passing_cut[fom_max_index])/total_continuum_bkg_event
    sig_efficiency = 100 * total_non_continuum_sig_event_passing_cut[fom_max_index]/total_non_continuum_sig_event

    # To use numpy array in any pyroot fuction we need to flatten it and change its datatype
    cs_prob_cut_TGraph_x = cs_prob_cut_array.flatten()
    cs_prob_fom_TGraph_y = cs_prob_fom_array.flatten()

    c1 = ROOT.TCanvas("fom", "fom", 900, 700)
    c1.SetFillColor(0)
    c1.SetGrid()

    GR2 = ROOT.TGraph(cs_prob_cut_array.size, cs_prob_cut_TGraph_x, cs_prob_fom_TGraph_y)
    GR2.SetLineColor(2)
    GR2.SetLineWidth(1)
    GR2.SetMarkerColor(9)
    GR2.SetMarkerStyle(69)
    GR2.SetMarkerSize(0.5)
    GR2.SetTitle("FOM analysis on FBDT output")
    GR2.GetXaxis().SetTitle("Cut")
    GR2.GetXaxis().SetTitleSize(0.05)
    GR2.GetYaxis().SetTitle("S/ #sqrt{S+B}")
    GR2.GetYaxis().SetTitleSize(0.05)
    GR2.Draw("ACP")

    # Drawing arrow at max fom
    arrowHeight = fom_max
    arrowHeight_y_min = arrowHeight*0.81
    arrowHeight_y_max = arrowHeight*0.97
    arr1 = ROOT.TArrow(-0.01, 100, -0.01, 400, 0.01,"|>")
    arr1.SetLineWidth(4)
    arr1.SetLineColor(2)
    arr1.SetFillStyle(3008)
    arr1.DrawArrow(cs_prob_cut_for_fom_max, arrowHeight_y_min, cs_prob_cut_for_fom_max, arrowHeight_y_max, 0.02,">") 

    # Write cont_rejection_percentage and sig_efficiency on the canvas and ALSO THE CUT VALUE
    my_info = ROOT.TLatex()
    my_info.SetTextSize(0.035)
    my_info.SetTextAlign(12)  #centered aligned
    my_info.DrawLatex(0.3,fom_max * 0.8,f"cut on FBDT output: {cs_prob_cut_for_fom_max}")
    my_info.DrawLatex(0.3,fom_max * 0.75,f"{cont_rejection_percentage:.2f}% bkg rejection")
    my_info.DrawLatex(0.3,fom_max * 0.7,f"{sig_efficiency:.2f}% sig efficiency")

    fom_plot_output_path = fom_txt_output_path
    fom_plot_output_filename = "fom_plot_mbcde_kid_md.pdf"
    c1.SaveAs(f"{fom_plot_output_path}/{fom_plot_output_filename}")