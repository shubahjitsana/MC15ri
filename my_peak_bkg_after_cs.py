#!/usr/bin/env python
# -*- coding: utf-8 -*-
import ROOT
import sys
# import numpy as np

# get input/output path/file from the command line
# 1st 5 are for input and next 5 are for output
input_filename = sys.argv[1]
output_filename = sys.argv[2]

# Loading input root file and creating new root file
inFile = ROOT.TFile.Open(f"{input_filename}")
inTree = inFile.Get('tree')

outFile = ROOT.TFile.Open(f"{output_filename}","RECREATE"); 
outTree = inTree.CloneTree(0)
print(f"Taking Peak bkg events from {input_filename} and saving to {output_filename}.")

# Scaling the file
Number_of_selected_events = 0       # To count and print the number of selected events per file
for iEvent in range(inTree.GetEntries()):
    inTree.GetEntry(iEvent)
    deltae = getattr(inTree, 'deltaE')
    mbc = getattr(inTree, 'Mbc')
    md0 = getattr(inTree, 'D0_bar_InvM')

    kid = getattr(inTree, 'Kp_PID_bin_kaon')
    signal_prob = getattr(inTree, 'SigProb')

    issig = getattr(inTree, 'isSignal')
    D_s_InvM = getattr(inTree, 'D_s_InvM')
    D_10D_md = getattr(inTree, 'D_10D_md')
    InvM1stand2ndpi = getattr(inTree, 'InvM1stand2ndpi')
    InvM2ndand3rdpi = getattr(inTree, 'InvM2ndand3rdpi')
    InvMD0and1stpi = getattr(inTree, 'InvMD0and1stpi')
    InvMD0and2ndpi = getattr(inTree, 'InvMD0and2ndpi')
    InvMD0and3rdpi = getattr(inTree, 'InvMD0and3rdpi')
    DstrminusroeD_md = getattr(inTree, 'DstrminusroeD_md')
    if (
        (mbc > 5.23 and mbc < 5.29) and 
        (deltae > -0.1 and deltae < 0.1) and  
        (md0 > 1.85 and md0 < 1.88) and 

        (D_s_InvM < 1.958 or D_s_InvM > 1.978) and 
        (D_10D_md < 0.496 or D_10D_md > 0.628) and 
        (InvM1stand2ndpi < 0.489 or InvM1stand2ndpi > 0.506) and 
        (InvM2ndand3rdpi < 0.49 or InvM2ndand3rdpi > 0.505) and 
        (InvMD0and1stpi < 2.0085 or InvMD0and1stpi > 2.0122) and 
        (InvMD0and2ndpi < 2.0087 or InvMD0and2ndpi > 2.0118) and 
        (InvMD0and3rdpi < 2.009 or InvMD0and3rdpi > 2.0116) and 
        (DstrminusroeD_md < 0.14402 or DstrminusroeD_md > 0.14682) and

        (kid > 0.6) and 
        (signal_prob > 0.21) and
        (issig !=1)
    ):
        outTree.Fill()
        Number_of_selected_events += 1

print(f"Afer scaling {Number_of_selected_events} events out of {inTree.GetEntries()} events have been stored in {output_filename}")
outTree.AutoSave()
outFile.Close()