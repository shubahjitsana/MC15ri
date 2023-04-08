#!/usr/bin/env python
# -*- coding: utf-8 -*-
# takes commadnd as "pyroot my_anyfile_anyvariable_bin_range.py file_path variable_name bin_number lower_range upper_range x_title output_filename"

import ROOT
import sys
ROOT.TH1.AddDirectory(ROOT.kFALSE)

input_filename = sys.argv[1]
variable_name = sys.argv[2]
bin_number = int(sys.argv[3])
lower_range = float(sys.argv[4])
upper_range = float(sys.argv[5])
# histogram_title = variable_name
input_from_6th_string = sys.argv[6].split("/")
# histogram_title = ""
# for i in input_from_6th_string:
#     histogram_title += i
#     histogram_title += " "
histogram_title = ' '.join(input_from_6th_string)

# Different output filename option for different files in one folder but ploting same variable
# input_filename_string = input_filename.split("data/combined/")    # ONLY need to change here depending on the pathname you are giving
# input_filename_last_string = input_filename_string[1]
# output_filename_option = input_filename_last_string.split(".root")
# output_filename = output_filename_option[0] + "_" + variable_name
output_filename = sys.argv[7]



print(f"Input file name is: {input_filename}")
print(f"{variable_name} is being plotted.")
print(f"Number of bins are: {bin_number}")
print(f"Lower range of histogram is: {lower_range}")
print(f"Upper range of histogram is: {upper_range}")
print(f"Title of histogram is: {histogram_title}")

input_root_file = ROOT.TFile.Open(f"{input_filename}")
input_tree = input_root_file.Get('tree')
total_event_number = input_tree.GetEntries()

histogram = []
histogram.append(ROOT.TH1F(f"{histogram_title}", f"Distribution of {histogram_title}", bin_number, lower_range, upper_range))
for iEvent in range(total_event_number):
    input_tree.GetEntry(iEvent)
    # variable_value = getattr(input_tree, f"{variable_name}")
    # histogram[0].Fill(variable_value)
    value_of_var = getattr(input_tree, f"{variable_name}")

    deltae = getattr(input_tree, 'deltaE')
    mbc = getattr(input_tree, 'Mbc')
    kid = getattr(input_tree, 'Kp_PID_bin_kaon')
    D_s_InvM = getattr(input_tree, 'D_s_InvM')
    D_10D_md = getattr(input_tree, 'D_10D_md')
    InvM1stand2ndpi = getattr(input_tree, 'InvM1stand2ndpi')
    InvM2ndand3rdpi = getattr(input_tree, 'InvM2ndand3rdpi')
    InvMD0and1stpi = getattr(input_tree, 'InvMD0and1stpi')
    InvMD0and2ndpi = getattr(input_tree, 'InvMD0and2ndpi')
    InvMD0and3rdpi = getattr(input_tree, 'InvMD0and3rdpi')
    DstrminusroeD_md = getattr(input_tree, 'DstrminusroeD_md')
    issignal = getattr(input_tree, 'isSignal')

    if (
        (issignal == 1) and
        (value_of_var > lower_range and value_of_var < upper_range) and
        (deltae > -0.15 and deltae < 0.15) and 
        (mbc > 5.23 and mbc < 5.29) and 
        # (md0 > 1.85 and md0 < 1.88) # and 
        # (kid > 0.6) and 
        (D_s_InvM < 1.958 or D_s_InvM > 1.978) and 
        (D_10D_md < 0.496 or D_10D_md > 0.628) and 
        (InvM1stand2ndpi < 0.489 or InvM1stand2ndpi > 0.506) and 
        (InvM2ndand3rdpi < 0.49 or InvM2ndand3rdpi > 0.505) and 
        (InvMD0and1stpi < 2.0085 or InvMD0and1stpi > 2.0122) and 
        (InvMD0and2ndpi < 2.0087 or InvMD0and2ndpi > 2.0118) and 
        (InvMD0and3rdpi < 2.009 or InvMD0and3rdpi > 2.0116) and 
        (DstrminusroeD_md < 0.14402 or DstrminusroeD_md > 0.14682)
    ):
        histogram[0].Fill(value_of_var)

colors = [ROOT.kViolet+10]
bin_size_MeV = ((upper_range - lower_range) * 1000)/bin_number      # convert GeV to MeV
c1 = ROOT.TCanvas(f"{histogram_title}", f"{histogram_title}", 900, 700)
histogram[0].SetLineColor(colors[0])
# histogram[0].SetLineColor(ROOT.kViolet+10)
histogram[0].GetXaxis().SetTitle(f"{histogram_title}  (GeV/c^{{2}})")
histogram[0].GetXaxis().SetTitleSize(0.035)
histogram[0].GetXaxis().SetLabelSize(0.02)
histogram[0].GetXaxis().CenterTitle(ROOT.kTRUE)
histogram[0].GetYaxis().SetTitle(f"Events/{bin_size_MeV:.2f} (MeV/c^{{2}})")
histogram[0].GetYaxis().SetTitleSize(0.035)
histogram[0].GetYaxis().SetTitleOffset(0.95)
histogram[0].GetYaxis().SetLabelSize(0.02)
histogram[0].GetYaxis().CenterTitle(ROOT.kTRUE)
# histogram.SetStats(ROOT.kFALSE)
histogram[0].Draw()
c1.Update()
c1.SaveAs(f"plot_from_my_anyfile_anyvariable_bin_range/{output_filename}.png")