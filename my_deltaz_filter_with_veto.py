#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import sys
import numpy as np

# get input/output path/file from the command line
# 1st 6 are for input and next 6 are for output
input_filename = []
output_filename = []
# for i in range(1, 7):
for i in range(1, 6):
    input_filename.append(sys.argv[i])
print(f"5 Input files are {input_filename}")
# for j in range(7,13):
for j in range(6,11):
    output_filename.append(sys.argv[j])
print(f"5 Output files are {output_filename}")

for i in range(5):
    # Loading input root file and creating new root file
    inFile = ROOT.TFile.Open(f"{input_filename[i]}")
    inTree = inFile.Get('tree')

    outFile = ROOT.TFile.Open(f"{output_filename[i]}","RECREATE"); 
    outTree = inTree.CloneTree(0)
    print(f"Filtering deltaZ and applying veto at {input_filename[i]} and saving to {output_filename[i]}")

    # Selecting events
    deltaZ_is_not_NAN_veto_passed_event_number = 0
    for iEvent in range(inTree.GetEntries()):   
        inTree.GetEntry(iEvent)
        deltae = getattr(inTree, 'deltaE')
        mbc = getattr(inTree, 'Mbc')
        md0 = getattr(inTree, 'D0_bar_InvM')
        deltaz = getattr(inTree, 'DeltaZ')
        # kid = getattr(inTree, 'Kp_PID_bin_kaon')
        # D_s_InvM = getattr(inTree, 'D_s_InvM')
        # D_10D_md = getattr(inTree, 'D_10D_md')
        # InvM1stand2ndpi = getattr(inTree, 'InvM1stand2ndpi')
        # InvM2ndand3rdpi = getattr(inTree, 'InvM2ndand3rdpi')
        # InvMD0and1stpi = getattr(inTree, 'InvMD0and1stpi')
        # InvMD0and2ndpi = getattr(inTree, 'InvMD0and2ndpi')
        # InvMD0and3rdpi = getattr(inTree, 'InvMD0and3rdpi')
        # DstrminusroeD_md = getattr(inTree, 'DstrminusroeD_md')
        if not np.isnan(deltaz):
            if (
                (deltae > -0.1 and deltae < 0.1) and 
                (mbc > 5.23 and mbc < 5.29) and 
                (md0 > 1.85 and md0 < 1.88) # and 
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
                outTree.Fill()
                deltaZ_is_not_NAN_veto_passed_event_number += 1
    print(f"Afer deltaz filtering and using veto from {input_filename[i]} {deltaZ_is_not_NAN_veto_passed_event_number} events out of {inTree.GetEntries()} have been stored in {output_filename[i]}.")
    outTree.AutoSave()
    outFile.Close()