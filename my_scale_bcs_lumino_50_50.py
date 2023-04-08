#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file scale such a way that 50% will be signal and 50% will be continnum background
import ROOT
import sys
import numpy as np

# get input/output path/file from the command line
# 1st 5 are for input and next 5 are for output
input_filename = []
output_filename = []
for i in range(1, 6):
    input_filename.append(sys.argv[i])
print(f"5 Input files are {input_filename}")
for j in range(6,11):
    output_filename.append(sys.argv[j])
print(f"5 Output files are {output_filename}")

# Getting scales
infile_uubar = ROOT.TFile.Open(f"{input_filename[0]}")
infile_ddbar = ROOT.TFile.Open(f"{input_filename[1]}")
infile_ssbar = ROOT.TFile.Open(f"{input_filename[2]}")
infile_ccbar = ROOT.TFile.Open(f"{input_filename[3]}")
infile_bbbar = ROOT.TFile.Open(f"{input_filename[4]}")

intree_uubar = infile_uubar.Get('tree')
intree_ddbar = infile_ddbar.Get('tree')
intree_ssbar = infile_ssbar.Get('tree')
intree_ccbar = infile_ccbar.Get('tree')
intree_bbbar = infile_bbbar.Get('tree')

# Counting number of events of different folder before scaling
bbbar_Event = np.zeros(4)   # vector to multiply with cross section percentage
bbbar_Event[0] = intree_bbbar.GetEntries()
bbbar_Event[1] = intree_bbbar.GetEntries()
bbbar_Event[2] = intree_bbbar.GetEntries()
bbbar_Event[3] = intree_bbbar.GetEntries()


NEvent = np.zeros(4)
NEvent[0] = intree_uubar.GetEntries()
NEvent[1] = intree_ddbar.GetEntries()
NEvent[2] = intree_ssbar.GetEntries()
NEvent[3] = intree_ccbar.GetEntries()
print(f"The event number of u,d,s,c(bar) are : {NEvent}")

NEvent_total = NEvent.sum()
cross_section_ratio_in_percentage = np.array([0.43168, 0.10785, 0.10301, 0.35745])
print(f"The cross sections ratio of u,d,s,c(bar) are : {cross_section_ratio_in_percentage}")

New_NEvent = bbbar_Event * cross_section_ratio_in_percentage

scale = NEvent/New_NEvent                                       # Setting Scale
print(f"Previous total events in u,d,s,c,b(bar) folder were {NEvent}")
print(f"New total event according to cross section should be {New_NEvent}")
print(f"So scaling factor is {scale}")

for i in range(4):
    # Loading input root file and creating new root file
    inFile = ROOT.TFile.Open(f"{input_filename[i]}")
    inTree = inFile.Get('tree')

    outFile = ROOT.TFile.Open(f"{output_filename[i]}","RECREATE"); 
    outTree = inTree.CloneTree(0)
    print(f"Scaling from {input_filename[i]} and saving to {output_filename[i]} \
    with scale factor {scale[i]}")

    # Scaling the file
    Number_of_selected_events = 0       # To count and print the number of selected events per file
    for iEvent in range(inTree.GetEntries()):
        rand = ROOT.gRandom.Rndm()
        num = 1/scale[i]
        if rand>num:    #shouldn't be "<" symbol bcz "continue" BREAKS one iteration
            continue     
        inTree.GetEntry(iEvent)
        outTree.Fill()
        Number_of_selected_events += 1

    print(f"Afer scaling {Number_of_selected_events} events has been stored in {output_filename[i]}")
    outTree.AutoSave()
    outFile.Close()