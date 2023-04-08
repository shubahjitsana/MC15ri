#!/usr/bin/env python
# -*- coding: utf-8 -*-
# takes commadnd as
# basf2 root_plot_any_variable_of_4_files.py cs/test/prescale_combine/test_signal.root cs/test/prescale_combine/test_charged.root 
# cs/test/prescale_combine/test_mixed.root cs/test/prescale_combine/test_qqbar.root 
# deltaE 100 -0.1 0.1 \#DeltaE/of/reconstructed/events

# python3 root_plot_any_variable_of_4_files.py cs/test/prescale_combine/test_qqbar.root cs/test/prescale_combine/test_signal.root cs/test/prescale_combine/test_charged.root cs/test/prescale_combine/test_mixed.root deltaE 100 -0.1 0.1 \#DeltaE/of/reconstructed/events
# python3 root_plot_any_variable_of_4_files.py cs/test/prescale_combine/test_qqbar.root cs/test/prescale_combine/test_signal.root cs/test/prescale_combine/test_charged.root cs/test/prescale_combine/test_mixed.root Mbc 100 5.23 5.29 M_{bc}/of/reconstructed/events

import ROOT
import sys
import numpy as np
ROOT.TH1.AddDirectory(ROOT.kFALSE)

##############taking input from command line####################
input_filename = []
for i in range(1,5):
    input_filename.append(sys.argv[i])
print(f"4 Input files are {input_filename}")
variable_name = sys.argv[5]
bin_number = int(sys.argv[6])
lower_range = float(sys.argv[7])
upper_range = float(sys.argv[8])
histogram_title_string = sys.argv[9].split("/")
histogram_title = ' '.join(histogram_title_string)

print(f"{variable_name} is being plotted.")
print(f"Number of bins are: {bin_number}")
print(f"Lower range of histogram is: {lower_range}")
print(f"Upper range of histogram is: {upper_range}")
print(f"Title of histogram is: {histogram_title}")    

output_filename = variable_name

######################reading data###############
input_root_file = []
for i in range(len(input_filename)):
    input_root_file.append(ROOT.TFile.Open(f"{input_filename[i]}"))

input_tree = []
for i in range(len(input_root_file)):
    input_tree.append(input_root_file[i].Get('tree'))

total_event_number = []
for i in range(len(input_tree)):
    total_event_number.append(input_tree[i].GetEntries())


histogram = []      # Definig list of 4 histograms for ploting 4 profiles
for root_file_no in range(len(input_filename)):
    histogram.append(ROOT.TH1F(f"{histogram_title}", f"Distribution of {histogram_title}",bin_number, lower_range, upper_range))

# Filling all histograms
for root_file_no in range(len(input_tree)):
    for iEvent in range(total_event_number[root_file_no]):
        input_tree[root_file_no].GetEntry(iEvent)
        
        value_of_var = getattr(input_tree[root_file_no], f"{variable_name}")
        if (variable_name == "Mbc"): 
            deltae = getattr(input_tree[root_file_no], 'deltaE')
            signal_prob = getattr(input_tree[root_file_no], 'SigProb')
        elif (variable_name == "deltaE"):
            mbc = getattr(input_tree[root_file_no], 'Mbc')
            signal_prob = getattr(input_tree[root_file_no], 'SigProb')
        elif (variable_name == "SigProb"):
            deltae = getattr(input_tree[root_file_no], 'deltaE')
            mbc = getattr(input_tree[root_file_no], 'Mbc')

        md0 = getattr(input_tree[root_file_no], 'D0_bar_InvM')
        kid = getattr(input_tree[root_file_no], 'Kp_PID_bin_kaon')
        
        issig = getattr(input_tree[root_file_no], 'isSignal')
        iscontinuum = getattr(input_tree[root_file_no], 'isContinuumEvent')

        D_s_InvM = getattr(input_tree[root_file_no], 'D_s_InvM')
        D_10D_md = getattr(input_tree[root_file_no], 'D_10D_md')
        InvM1stand2ndpi = getattr(input_tree[root_file_no], 'InvM1stand2ndpi')
        InvM2ndand3rdpi = getattr(input_tree[root_file_no], 'InvM2ndand3rdpi')
        InvMD0and1stpi = getattr(input_tree[root_file_no], 'InvMD0and1stpi')
        InvMD0and2ndpi = getattr(input_tree[root_file_no], 'InvMD0and2ndpi')
        InvMD0and3rdpi = getattr(input_tree[root_file_no], 'InvMD0and3rdpi')
        DstrminusroeD_md = getattr(input_tree[root_file_no], 'DstrminusroeD_md')

        if (variable_name == "Mbc"):
            if (
                (value_of_var > lower_range and value_of_var < upper_range) and
                (deltae > -0.1 and deltae < 0.1) and  
                (md0 > 1.85 and md0 < 1.88) and 
                (kid > 0.6) and 
                (D_s_InvM < 1.958 or D_s_InvM > 1.978) and 
                (D_10D_md < 0.496 or D_10D_md > 0.628) and 
                (InvM1stand2ndpi < 0.489 or InvM1stand2ndpi > 0.506) and 
                (InvM2ndand3rdpi < 0.49 or InvM2ndand3rdpi > 0.505) and 
                (InvMD0and1stpi < 2.0085 or InvMD0and1stpi > 2.0122) and 
                (InvMD0and2ndpi < 2.0087 or InvMD0and2ndpi > 2.0118) and 
                (InvMD0and3rdpi < 2.009 or InvMD0and3rdpi > 2.0116) and 
                (DstrminusroeD_md < 0.14402 or DstrminusroeD_md > 0.14682) and
                (signal_prob > 0.21)
            ):
                if root_file_no == 1:
                    if issig==1: histogram[root_file_no].Fill(value_of_var)
                else: histogram[root_file_no].Fill(value_of_var)
        elif (variable_name == "deltaE"):
            if (
                (value_of_var > lower_range and value_of_var < upper_range) and
                (mbc > 5.23 and mbc < 5.29) and 
                (md0 > 1.85 and md0 < 1.88) and 
                (kid > 0.6) and 
                (D_s_InvM < 1.958 or D_s_InvM > 1.978) and 
                (D_10D_md < 0.496 or D_10D_md > 0.628) and 
                (InvM1stand2ndpi < 0.489 or InvM1stand2ndpi > 0.506) and 
                (InvM2ndand3rdpi < 0.49 or InvM2ndand3rdpi > 0.505) and 
                (InvMD0and1stpi < 2.0085 or InvMD0and1stpi > 2.0122) and 
                (InvMD0and2ndpi < 2.0087 or InvMD0and2ndpi > 2.0118) and 
                (InvMD0and3rdpi < 2.009 or InvMD0and3rdpi > 2.0116) and 
                (DstrminusroeD_md < 0.14402 or DstrminusroeD_md > 0.14682) and
                (signal_prob > 0.21)
            ):
                if root_file_no == 1:
                    if issig==1: histogram[root_file_no].Fill(value_of_var)
                else: histogram[root_file_no].Fill(value_of_var)
        elif (variable_name == "SigProb"):
            if (
                (value_of_var > lower_range and value_of_var < upper_range) and
                (deltae > -0.1 and deltae < 0.1) and
                (mbc > 5.23 and mbc < 5.29) and 
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
                if root_file_no == 1:
                    if issig==1: histogram[root_file_no].Fill(value_of_var)
                else: histogram[root_file_no].Fill(value_of_var)
                
bin_size_MeV =  ((upper_range - lower_range) * 1000)/bin_number      # convert GeV to MeV


#################### Start Drawing ##############################
colors = [ROOT.kRed, ROOT.kViolet+10, ROOT.kGreen, ROOT.kBlack]
c1 = ROOT.TCanvas("c", "c", 3600, 2800)
if (variable_name == "SigProb"): legend = ROOT.TLegend(0.45,0.7,0.65,0.9)
else: legend = ROOT.TLegend(0.1,0.7,0.3,0.9)
c1.cd()

y_axis = np.zeros(4)
for root_file_no in range(len(input_filename)):
    y_axis[root_file_no] = histogram[root_file_no].GetMaximum()
# y_axis_max_index = y_axis.argmax()
y_axis_max = y_axis.max()

for root_file_no in range(len(input_filename)):
    if root_file_no == 0:
        histogram[root_file_no].GetXaxis().SetTitle(f"{histogram_title}  (GeV/c^{{2}})")
        histogram[root_file_no].GetXaxis().SetTitleSize(0.035)
        histogram[root_file_no].GetXaxis().SetLabelSize(0.02)
        histogram[root_file_no].GetXaxis().CenterTitle(ROOT.kTRUE)
        histogram[root_file_no].GetYaxis().SetTitle(f"Events/{bin_size_MeV:.2f} (MeV/c^{{2}})")
        histogram[root_file_no].GetYaxis().SetTitleSize(0.035)
        histogram[root_file_no].GetYaxis().SetTitleOffset(0.95)
        histogram[root_file_no].GetYaxis().SetLabelSize(0.02)
        histogram[root_file_no].GetYaxis().CenterTitle(ROOT.kTRUE)

        histogram[root_file_no].SetStats(ROOT.kFALSE)
        histogram[root_file_no].SetLineWidth(4)
        histogram[root_file_no].SetLineColor(colors[root_file_no])

        histogram[root_file_no].GetYaxis().SetRangeUser(0,(200 + y_axis_max))
        # if (variable_name == "SigProb"):
        #     # histogram[root_file_no].GetYaxis().SetRangeUser(0,(200 + histogram[y_axis_max_index].GetMaximum()))
        #     histogram[root_file_no].GetYaxis().SetRangeUser(0,(200 + y_axis_max))
        # else: histogram[root_file_no].GetYaxis().SetRangeUser(0,60000)
        histogram[root_file_no].Draw()
    # elif root_file_no == 1:
    #     histogram[root_file_no].Scale(0.2/8.7796)   # scale signalMC from 8.7796ab^(-1) to 200fb^(-1)
    #     histogram[root_file_no].SetStats(ROOT.kFALSE)
    #     histogram[root_file_no].SetLineWidth(4)
    #     histogram[root_file_no].SetLineColor(colors[root_file_no])
    #     histogram[root_file_no].Draw("SAME")
    else:
        histogram[root_file_no].SetStats(ROOT.kFALSE)
        histogram[root_file_no].SetLineWidth(4)
        histogram[root_file_no].SetLineColor(colors[root_file_no])
        histogram[root_file_no].Draw("SAME")

legend.AddEntry(histogram[0], "qqbar", "l")
legend.AddEntry(histogram[1], "signal", "l")
legend.AddEntry(histogram[2], "Charged", "l")
legend.AddEntry(histogram[3], "mixed", "l")
legend.Draw()


c1.Update()
c1.SaveAs(f"plot_from_my_anyfile_anyvariable_bin_range/{output_filename}_kin_veto_sigprob21_signal_scaled.png")