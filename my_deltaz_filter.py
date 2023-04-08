#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import sys
import numpy as np

# get input/output path/file from the command line
# 1st 6 are for input and next 6 are for output
input_filename = []
output_filename = []
for i in range(1, 7):
    input_filename.append(sys.argv[i])
print(f"6 Input files are {input_filename}")
for j in range(7,13):
    output_filename.append(sys.argv[j])
print(f"6 Output files are {output_filename}")

for i in range(6):
    # Loading input root file and creating new root file
    inFile = ROOT.TFile.Open(f"{input_filename[i]}")
    inTree = inFile.Get('tree')

    outFile = ROOT.TFile.Open(f"{output_filename[i]}","RECREATE"); 
    outTree = inTree.CloneTree(0)
    print(f"Filtering deltaZ at {input_filename[i]} and saving to {output_filename[i]}")

    # Selecting events
    deltaZ_is_not_NAN_event_number = 0
    for iEvent in range(inTree.GetEntries()):   
        inTree.GetEntry(iEvent)
        deltaz = getattr(inTree, 'DeltaZ')
        if not np.isnan(deltaz):
            outTree.Fill()
            deltaZ_is_not_NAN_event_number += 1
    print(f"Afer deltaz filtering from {input_filename[i]} {deltaZ_is_not_NAN_event_number} events out of {inTree.GetEntries()} have been stored in {output_filename[i]}.")
    outTree.AutoSave()
    outFile.Close()