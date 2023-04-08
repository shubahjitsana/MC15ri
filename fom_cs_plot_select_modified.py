#!/usr/bin/env python
# -*- coding: utf-8 -*-
# basf2 steering.py file1.root

import ROOT
import sys
import numpy as np

print(f"Number of input files are {len(sys.argv)-1}")
print(f"{sys.argv[0]}")
for i in sys.argv[1:]:
    input_file = ROOT.TFile.Open(f"{i}")    # Loading root file
    input_tree = input_file.Get('tree')
    total_event_in_charged = input_tree.GetEntries()

    # getting total number of continuum and non_continuum events
    # we need this to calculate cont_bkg rejection ratio and sig_efficiency
    # fom analysis
    total_bkg_event_just_counting = 0
    total_sig_event_just_counting = 0
    number_of_bin = 100
    sig_passing_cut = np.zeros(number_of_bin)
    bkg_passing_cut = np.zeros(number_of_bin)
    signal_prob_cut = np.zeros(number_of_bin)
    total_sig_event = np.zeros(number_of_bin)
    total_bkg_event = np.zeros(number_of_bin)

    for index in range(number_of_bin):
        signal_prob_cut[index] = 0 + index * 0.01
    # calculating FOM
    for iEvent in range(total_event_in_charged):
        input_tree.GetEntry(iEvent)
        deltae = getattr(input_tree, 'deltaE')
        mbc = getattr(input_tree, 'Mbc')
        md0 = getattr(input_tree, 'D0_bar_InvM')
        iscontinuum = getattr(input_tree, 'isContinuumEvent')
        signal_prob = getattr(input_tree, 'SigProb')

        # kid = getattr(input_tree, 'Kp_PID_bin_kaon')
        # issig = getattr(input_tree, 'isSignal')

        # D_s_InvM = getattr(input_tree, 'D_s_InvM')
        # D_10D_md = getattr(input_tree, 'D_10D_md')
        # InvM1stand2ndpi = getattr(input_tree, 'InvM1stand2ndpi')
        # InvM2ndand3rdpi = getattr(input_tree, 'InvM2ndand3rdpi')
        # InvMD0and1stpi = getattr(input_tree, 'InvMD0and1stpi')
        # InvMD0and2ndpi = getattr(input_tree, 'InvMD0and2ndpi')
        # InvMD0and3rdpi = getattr(input_tree, 'InvMD0and3rdpi')
        # DstrminusroeD_md = getattr(input_tree, 'DstrminusroeD_md')
        if (
            (deltae > -0.02 and deltae < 0.02) and 
            (mbc > 5.27 and mbc < 5.29) and 
            (md0 > 1.85 and md0 < 1.88) #and 
            # (kid > 0.6) and 
            # (D_s_InvM < 1.958 or D_s_InvM > 1.978) and 
            # (D_10D_md < 0.496 or D_10D_md > 0.628) and 
            # (InvM1stand2ndpi < 0.489 or InvM1stand2ndpi > 0.506) and 
            # (InvM2ndand3rdpi < 0.49 or InvM2ndand3rdpi > 0.505) and 
            # (InvMD0and1stpi < 2.0085 or InvMD0and1stpi > 2.0122) and 
            # (InvMD0and2ndpi < 2.0087 or InvMD0and2ndpi > 2.0118) and 
            # (InvMD0and3rdpi < 2.009 or InvMD0and3rdpi > 2.0116) and 
            # (DstrminusroeD_md < 0.14402 or DstrminusroeD_md > 0.14682)
        ):
            if iscontinuum == 1: total_bkg_event_just_counting +=1   # getting total number of continuum and non_continuum events
            else: total_sig_event_just_counting += 1    # within signal region & passing veto cut

            for index in range(number_of_bin):
                if ((signal_prob > signal_prob_cut[index])):  
                    # bkg += 1 if iscontinuum == 1 else sig += 1  # can't use arithmetic operator
                    if iscontinuum == 1: bkg_passing_cut[index] += 1
                    else: sig_passing_cut[index] +=1
    extension = "_mbcde_kin"

    # Calculating all fom value using numpy array   
    total_sig_event = total_sig_event + total_sig_event_just_counting
    total_bkg_event = total_bkg_event + total_bkg_event_just_counting
    bkg_rejection_percentage = 100 * (total_bkg_event - bkg_passing_cut)/total_bkg_event
    sig_efficiency = 100 * sig_passing_cut/total_sig_event

    cs_prob_fom = sig_passing_cut/np.sqrt(sig_passing_cut + bkg_passing_cut)

    # Defining and exporting report file
    input_path = i.split("/all.root")
    fom_txt_output_path = input_path[0]
    fom_txt_output_filename = f"fom{extension}.txt"
    fom_file = open(f"{fom_txt_output_path}/{fom_txt_output_filename}", 'w')    # defining exporting file

    fom_detail_report_filename = f"fom_detail_report{extension}.txt"
    fom_detail_report = open(f"{fom_txt_output_path}/{fom_detail_report_filename}", 'w')    # defining exporting file

    #writing in export file
    for index in range(number_of_bin):
        fom_file.write(repr(signal_prob_cut[index]) + '\t' + repr(cs_prob_fom[index]) + '\n')
        fom_detail_report.write(repr(signal_prob_cut[index]) + '\t' + repr(sig_passing_cut[index]) +
         '\t' + repr(bkg_passing_cut[index]) + '\t' + repr(cs_prob_fom[index]) + '\t' + 
         repr(bkg_rejection_percentage[index]) + '\t' +  repr(sig_efficiency[index]) + '\n')
    fom_file.close()
    fom_detail_report.close()

    # store value in np array and plot those in Tgraph(need to flatten before plot)
    # Get max value from TGraph(also can get from np but need to access coordinate in CANVAS)
    # calculate bkg rejection and sig efficiency for MAX FOM value

    # Calculate cont_bkg rejection ratio and sig_efficiency
    fom_max = cs_prob_fom.max()
    fom_max_index = cs_prob_fom.argmax()
    # signal_prob_cut[fom_max_index] = 0.01 + 0.01 * fom_max_index
    # signal_prob_cut[fom_max_index] = 0.1 + 0.01 * fom_max_index  # As cut 0.1 has been saved in 1st index(not zeroth index)

    print(f"Maximum output is corresponding to FBDT output: {signal_prob_cut[fom_max_index]:.2f} with \
    {bkg_rejection_percentage[fom_max_index]:.2f}% bkg rejection and {sig_efficiency[fom_max_index]:.2f}% sig efficiency")

    # To use numpy array in any pyroot fuction we need to flatten it and change its datatype
    signal_prob_cut_TGraph_x = signal_prob_cut.flatten()
    cs_prob_fom_TGraph_y = cs_prob_fom.flatten()

    c1 = ROOT.TCanvas("fom", "fom", 900, 700)
    c1.SetFillColor(0)
    c1.SetGrid()

    GR2 = ROOT.TGraph(signal_prob_cut.size, signal_prob_cut_TGraph_x, cs_prob_fom_TGraph_y)
    GR2.SetLineColor(2)
    GR2.SetLineWidth(1)
    GR2.SetMarkerColor(9)
    GR2.SetMarkerStyle(69)
    GR2.SetMarkerSize(0.5)
    GR2.SetTitle(f"FOM analysis on FBDT output using cuts on {extension}")
    GR2.GetXaxis().SetTitle("Cut")
    GR2.GetXaxis().SetTitleSize(0.05)
    GR2.GetYaxis().SetTitle("S/ #sqrt{S+B}")
    GR2.GetYaxis().SetTitleSize(0.05)
    GR2.Draw("ACP")

    # Drawing arrow at max fom
    arrowHeight = fom_max
    arrowHeight_y_min = arrowHeight*0.95
    arrowHeight_y_max = arrowHeight*0.99
    arr1 = ROOT.TArrow(-0.01, 100, -0.01, 400, 0.01,"|>")
    arr1.SetLineWidth(4)
    arr1.SetLineColor(2)
    arr1.SetFillStyle(3008)
    arr1.DrawArrow(signal_prob_cut[fom_max_index], arrowHeight_y_min, signal_prob_cut[fom_max_index], arrowHeight_y_max, 0.02,">") 

    # Write bkg_rejection_percentage and sig_efficiency on the canvas and ALSO THE CUT VALUE
    my_info = ROOT.TLatex()
    my_info.SetTextSize(0.035)
    my_info.SetTextAlign(12)  #centered aligned
    my_info.DrawLatex(0.1,fom_max * 0.9,f"cut on FBDT output: {signal_prob_cut[fom_max_index]:.2f}")
    my_info.DrawLatex(0.1,fom_max * 0.86,f"{bkg_rejection_percentage[fom_max_index]:.2f}% bkg rejection")
    my_info.DrawLatex(0.1,fom_max * 0.82,f"{sig_efficiency[fom_max_index]:.2f}% sig efficiency")

    fom_plot_output_path = fom_txt_output_path
    fom_plot_output_filename = f"fom_plot{extension}.pdf"
    c1.SaveAs(f"{fom_plot_output_path}/{fom_plot_output_filename}")
