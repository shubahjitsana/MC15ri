#!/usr/bin/env python
# -*- coding: utf-8 -*-
# takes commadnd as
# basf2 root_plot_4_cs_variables_of_2_bkg_files_train_test.py cs/train/signal_scaled/train_deltaz_filtered_qqbar.root 
# cs/test_1_35/signal_scaled/test_deltaz_filtered_qqbar.root 
# R2 500 1.9 5.1 Invariant/mass/of/#bar{D}^{0}/and/1^{st}/pion 
# DeltaZ 500 1.9 5.1 Invariant/mass/of/#bar{D}^{0}/and/2^{nd}/pion 
# abs_qr 500 1.9 5.1 Invariant/mass/of/#bar{D}^{0}/and/3^{rd}/pion 
# thrustBm 500 0 3.7 Invariant/mass/of/1^{st}/2^{nd}/and/3^{rd}/pion

import ROOT
import sys
import numpy as np
ROOT.TH1.AddDirectory(ROOT.kFALSE)

normalization = int(sys.argv[1]) # "0" means not normalized; "1" means normalized

##############taking input from command line####################
input_filename = []
input_filename.append(sys.argv[2])
input_filename.append(sys.argv[3])
# input_filename.append("cs/train/signal_scaled/train_deltaz_filtered_qqbar.root")
# input_filename.append("cs/test_1_35/signal_scaled/test_deltaz_filtered_qqbar.root")
print(f"2 Input files are {input_filename}")

variable_name = []
variable_name.append("cosTBTO")
variable_name.append("KSFWVariables__bohoo0__bc")
variable_name.append("KSFWVariables__bohso12__bc")
variable_name.append("cosTBz")

bin_number = []
bin_number.append(200)
bin_number.append(200)
bin_number.append(200)
bin_number.append(200)

lower_range = []
lower_range.append(0)
lower_range.append(0)
lower_range.append(-0.4)
lower_range.append(0)

upper_range = []
upper_range.append(1.1)
upper_range.append(0.35)
upper_range.append(0.65)
upper_range.append(1)

histogram_title = []
histogram_title.append("#left|cos#theta_{TBTO}#right|")
histogram_title.append("KSFW H^{oo}_{0}")
histogram_title.append("KSFW H^{so}_{12}")
histogram_title.append("#left|cos#theta_{TBz}#right|")

for k in range(len(variable_name)):
    print(f"{variable_name[k]} is being plotted.")
    print(f"Number of bins are: {bin_number[k]}")
    print(f"Lower range of histogram is: {lower_range[k]}")
    print(f"Upper range of histogram is: {upper_range[k]}")
    print(f"Title of histogram is: {histogram_title[k]}")    

if normalization == 1:
    output_filename = "cs_" + sys.argv[6] + "2_" + '_'.join(variable_name) + "_normalized"
else:
    output_filename = "cs_" + sys.argv[6] + "2_" + '_'.join(variable_name)

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

histogram = []      ## will be define as 4*7
for var in range(len(variable_name)):
    histogram_row = []      ## 4-D list to create new ROW  # Definig list of 4 histograms for ploting 4 variables
    for root_file_no in range(len(input_filename)):
        if normalization == 1:
            histogram_row.append(ROOT.TH1F(f"{histogram_title[var]}", f"Normalized distribution of {histogram_title[var]}",bin_number[var], lower_range[var], upper_range[var]))
        else:
            histogram_row.append(ROOT.TH1F(f"{histogram_title[var]}", f"Distribution of {histogram_title[var]}",bin_number[var], lower_range[var], upper_range[var]))
                #Adding elements(2) through row-wise #ploting same variable for given all input files
    histogram.append(histogram_row)     # 4*2-D ###  appending ROWs(4) along column

# Filling all histograms
for root_file_no in range(len(input_tree)):
    for iEvent in range(total_event_number[root_file_no]):
        input_tree[root_file_no].GetEntry(iEvent)
        for var in range(len(variable_name)):
            histogram[var][root_file_no].Fill(getattr(input_tree[root_file_no], f"{variable_name[var]}"))

upper_range_np = np.array(upper_range)
lower_range_np = np.array(lower_range)
bin_number_np = np.array(bin_number)
bin_size_MeV =  ((upper_range_np - lower_range_np))/bin_number_np      # convert GeV to MeV CANCELLED BY ==> # * 1000


#################### Start Drawing ##############################
colors = [ROOT.kRed, ROOT.kViolet+10]
c1 = ROOT.TCanvas("c", "c", 3600, 2800)
# c1.Divide(2,2)
pad = []
pad.append(ROOT.TPad("pad1", "pad1", 0, 0.495, 0.505, 1.0))
pad.append(ROOT.TPad("pad2", "pad2", 0.495, 0.495, 1.0, 1.0))
pad.append(ROOT.TPad("pad3", "pad3", 0, 0, 0.505, 0.505))
pad.append(ROOT.TPad("pad4", "pad4", 0.495, 0, 1.0, 0.505))

legend = []
legend.append(ROOT.TLegend(0.1,0.7,0.3,0.9))
legend.append(ROOT.TLegend(0.7,0.7,0.9,0.9))
legend.append(ROOT.TLegend(0.7,0.7,0.9,0.9))
legend.append(ROOT.TLegend(0.7,0.7,0.9,0.9))

for var in range(len(variable_name)):
    # c1.cd(var+1)

    c1.cd()
    pad[var].Draw()
    pad[var].cd()

    for root_file_no in range(len(input_filename)):
        if root_file_no == 0:
            histogram[var][root_file_no].GetXaxis().SetTitle(f"{histogram_title[var]}")
            histogram[var][root_file_no].GetXaxis().SetTitleSize(0.035)
            histogram[var][root_file_no].GetXaxis().SetLabelSize(0.02)
            histogram[var][root_file_no].GetXaxis().CenterTitle(ROOT.kTRUE)
            histogram[var][root_file_no].GetYaxis().SetTitleSize(0.035)
            histogram[var][root_file_no].GetYaxis().SetTitleOffset(0.95)
            histogram[var][root_file_no].GetYaxis().SetLabelSize(0.02)
            histogram[var][root_file_no].GetYaxis().CenterTitle(ROOT.kTRUE)

            histogram[var][root_file_no].SetStats(ROOT.kFALSE)
            histogram[var][root_file_no].SetLineColor(colors[root_file_no])

            # Normalization
            if normalization == 1:
                histogram[var][root_file_no].Scale(1/total_event_number[root_file_no])
                histogram[var][root_file_no].GetYaxis().SetTitle(f"# Events per bin/# Events")
            else:
                histogram[var][root_file_no].GetYaxis().SetTitle(f"Events/{bin_size_MeV[var]:.2f}")
                histogram[var][root_file_no].GetYaxis().SetRangeUser(0,(200 + histogram[var][root_file_no].GetMaximum()))
            histogram[var][root_file_no].Draw()
        else:
            histogram[var][root_file_no].SetStats(ROOT.kFALSE)
            histogram[var][root_file_no].SetLineColor(colors[root_file_no])
            # Normalization
            if normalization == 1:
                histogram[var][root_file_no].Scale(1/total_event_number[root_file_no])
                histogram[var][root_file_no].Draw("SAME")
            else:
                histogram[var][root_file_no].Draw("SAME")

    legend[var].AddEntry(histogram[var][0], sys.argv[4], "l")
    legend[var].AddEntry(histogram[var][1], sys.argv[5], "l")
    # legend[var].AddEntry(histogram[var][0], "Traing data", "l")
    # legend[var].AddEntry(histogram[var][1], "Testing data", "l")
    legend[var].Draw()
        
    pad[var].Update()


c1.Update()
c1.SaveAs(f"plot_from_my_anyfile_anyvariable_bin_range/{output_filename}.png")