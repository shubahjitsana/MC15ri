#!/usr/bin/env python
# -*- coding: utf-8 -*-
# takes commadnd as
# basf2 root_plot_4_variables_of_7_files.py data/combined/all.root data/combined/charged.root data/combined/mixed.root 
# data/combined/uubar.root data/combined/ddbar.root data/combined/ssbar.root data/combined/ccbar.root 
# InvMD0and1stpi 500 1.9 5.1 Invariant/mass/of/#bar{D}^{0}/and/1^{st}/pion 
# InvMD0and2ndpi 500 1.9 5.1 Invariant/mass/of/#bar{D}^{0}/and/2^{nd}/pion 
# InvMD0and3rdpi 500 1.9 5.1 Invariant/mass/of/#bar{D}^{0}/and/3^{rd}/pion 
# InvM1st2nd3rdpi 500 0 3.7 Invariant/mass/of/1^{st}/2^{nd}/and/3^{rd}/pion

import ROOT
import sys
import numpy as np
ROOT.TH1.AddDirectory(ROOT.kFALSE)

##############taking input from command line####################
input_filename = []
for i in range(1,8):
    input_filename.append(sys.argv[i])
print(f"7 Input files are {input_filename}")
variable_name = []
for i in [8, 13, 18, 23]:
    variable_name.append(sys.argv[i])

bin_number = []
for i in [9, 14, 19, 24]:
    bin_number.append(int(sys.argv[i]))

lower_range = []
for i in [10, 15, 20, 25]:
    lower_range.append(float(sys.argv[i]))

upper_range = []
for i in [11, 16, 21, 26]:
    upper_range.append(float(sys.argv[i]))

histogram_title_string = []
for i in [12, 17, 22, 27]:
    histogram_title_string.append(sys.argv[i].split("/"))  # saving list in an list eliment
histogram_title = []
for i in range(len(histogram_title_string)):
    histogram_title.append(' '.join(histogram_title_string[i]))

for k in range(len(variable_name)):
    print(f"{variable_name[k]} is being plotted.")
    print(f"Number of bins are: {bin_number[k]}")
    print(f"Lower range of histogram is: {lower_range[k]}")
    print(f"Upper range of histogram is: {upper_range[k]}")
    print(f"Title of histogram is: {histogram_title[k]}")    

output_filename = '_'.join(variable_name)

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
        histogram_row.append(ROOT.TH1F(f"{histogram_title[var]}", f"Distribution of {histogram_title[var]}",bin_number[var], lower_range[var], upper_range[var]))
                #Adding elements(7) through row-wise #ploting same variable for given all input files
    histogram.append(histogram_row)     # 4*7-D ###  appending ROWs(4) along column

# Filling all histograms
for root_file_no in range(len(input_tree)):
    for iEvent in range(total_event_number[root_file_no]):
        input_tree[root_file_no].GetEntry(iEvent)
        for var in range(len(variable_name)):
            histogram[var][root_file_no].Fill(getattr(input_tree[root_file_no], f"{variable_name[var]}"))

upper_range_np = np.array(upper_range)
lower_range_np = np.array(lower_range)
bin_number_np = np.array(bin_number)
bin_size_MeV =  ((upper_range_np - lower_range_np) * 1000)/bin_number_np      # convert GeV to MeV


#################### Start Drawing ##############################
colors = [ROOT.kRed, ROOT.kViolet+10, ROOT.kGreen, ROOT.kYellow, ROOT.kMagenta, ROOT.kCyan, ROOT.kBlack]
c1 = ROOT.TCanvas("c", "c", 3600, 2800)
pad = []
pad.append(ROOT.TPad("pad1", "pad1", 0, 0.51, 0.49, 1.0))
pad.append(ROOT.TPad("pad2", "pad2", 0.5, 0.51, 1.0, 1.0))
pad.append(ROOT.TPad("pad3", "pad3", 0, 0, 0.49, 0.49))
pad.append(ROOT.TPad("pad4", "pad4", 0.5, 0, 1.0, 0.49))

legend = []
legend.append(ROOT.TLegend(0.7,0.7,0.9,0.9))
legend.append(ROOT.TLegend(0.7,0.7,0.9,0.9))
legend.append(ROOT.TLegend(0.7,0.7,0.9,0.9))
legend.append(ROOT.TLegend(0.7,0.7,0.9,0.9))

for var in range(len(variable_name)):
    c1.cd()
    pad[var].Draw()
    pad[var].cd()

    for root_file_no in range(len(input_filename)):
        if root_file_no == 0:
            histogram[var][root_file_no].GetXaxis().SetTitle(f"{histogram_title[var]}  (GeV/c^{{2}})")
            histogram[var][root_file_no].GetXaxis().SetTitleSize(0.035)
            histogram[var][root_file_no].GetXaxis().SetLabelSize(0.02)
            histogram[var][root_file_no].GetXaxis().CenterTitle(ROOT.kTRUE)
            histogram[var][root_file_no].GetYaxis().SetTitle(f"Events/{bin_size_MeV[var]:.2f} (MeV/c^{{2}})")
            histogram[var][root_file_no].GetYaxis().SetTitleSize(0.035)
            histogram[var][root_file_no].GetYaxis().SetTitleOffset(0.95)
            histogram[var][root_file_no].GetYaxis().SetLabelSize(0.02)
            histogram[var][root_file_no].GetYaxis().CenterTitle(ROOT.kTRUE)

            histogram[var][root_file_no].SetStats(ROOT.kFALSE)
            histogram[var][root_file_no].SetLineColor(colors[root_file_no])

            histogram[var][root_file_no].GetYaxis().SetRangeUser(0,(200 + histogram[var][root_file_no].GetMaximum()))
            histogram[var][root_file_no].Draw()
        else:
            histogram[var][root_file_no].SetStats(ROOT.kFALSE)
            histogram[var][root_file_no].SetLineColor(colors[root_file_no])
            histogram[var][root_file_no].Draw("SAME")

    legend[var].AddEntry(histogram[var][0], "For all data", "l")
    legend[var].AddEntry(histogram[var][1], "Only for charged", "l")
    legend[var].AddEntry(histogram[var][2], "Only for mixed", "l")
    legend[var].AddEntry(histogram[var][3], "Only for uubar", "l")
    legend[var].AddEntry(histogram[var][4], "Only for ddbar", "l")
    legend[var].AddEntry(histogram[var][5], "Only for ssbar", "l")
    legend[var].AddEntry(histogram[var][6], "Only for ccbar", "l")
    legend[var].Draw()
        
    pad[var].Update()


c1.Update()
c1.SaveAs(f"plot_from_my_anyfile_anyvariable_bin_range/{output_filename}.png")