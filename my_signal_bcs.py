#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import sys

# get input/output path/file from the command line
# 1st 2 are for input and next 2 are for output
input_filename = []
output_filename = []
# for i in range(1, 3):
for i in range(1, 2):
    input_filename.append(sys.argv[i])
print(f"1 Input file(s) is/are {input_filename}")
# for i in range(3, 5):
for i in range(2, 3):
    output_filename.append(sys.argv[i])
print(f"1 Ouput file(s) is/are {output_filename}")

total_signal_event_number = 0
for i in range(1):
    # Loading input root file and creating new root file
    inFile = ROOT.TFile.Open(f"{input_filename[i]}")
    inTree = inFile.Get('tree')
    print(f"Filtering signal events from {input_filename[i]} and saving to {output_filename[i]}")

    outFile = ROOT.TFile.Open(f"{output_filename[i]}","RECREATE"); 
    outTree = inTree.CloneTree(0)
    # Filtering only Signal Events
    signal_event_number = 0     # To count and print the number of selected signal
    total_event_in_file = inTree.GetEntries()    # To count and print the number of reconstructed event
    for iEvent in range(total_event_in_file):
        inTree.GetEntry(iEvent)
        signal = getattr(inTree, 'isSignal')
        if signal==1:
            outTree.Fill()
            signal_event_number += 1
    print(f"There were {total_event_in_file} reconstructed events in {input_filename[i]}")
    print(f"There were {signal_event_number} signal events which have been stored in {output_filename[i]}")
    total_signal_event_number += signal_event_number
    outTree.AutoSave()
    outFile.Close()
print(f"Total {total_signal_event_number} signal events have been stored in {output_filename}")