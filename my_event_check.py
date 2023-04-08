#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import sys
import numpy as np

print(f"Number of input files are {len(sys.argv)-1}")
print(f"{sys.argv[0]}")
for i in sys.argv[1:]:
    input_file = ROOT.TFile.Open(f"{i}")
    input_tree = input_file.Get('tree')
    total_event_in_charged = input_tree.GetEntries()
    number_of_variable_to_be_printed_in_one_line = 5
    count = 0
    print(f"There are total {total_event_in_charged} reconstructed events in {i}")
#     for iEvent in range(total_event_in_charged):
#         input_tree.GetEntry(iEvent)
#         # run_1 = getattr(input_tree, '__run__')
#         # event_1 = getattr(input_tree, '__event__')
#         # mbc_1 = getattr(input_tree, 'Mbc')
#         deltae_1 = getattr(input_tree, 'deltaE')
#         deltaz_1 = getattr(input_tree, 'DeltaZ')
#         # signal_1 = getattr(input_tree, 'isSignal')


#         # print(f"{iEvent} : run number : {run_1} ; event number: {event_1} ; issignal {signal_1}; Mbc : {mbc_1} ; deltaE : {deltae_1}")

#         # Printing "DeltaZ" value
#         # if count < number_of_variable_to_be_printed_in_one_line :
#         #     print(f"deltaz : {deltaz_1};", end=" ")
#         #     count += 1
#         # else:
#         #     print(f"deltaz : {deltaz_1}")
#         #     count = 0



#         # cheching nan value of deltaZ
#         if np.isnan(deltaz_1):
#             count += 1
#     print(f"There are total {count} events out of \
# {total_event_in_charged} reconstructed events whose DetlaZ vlaue is \"nan\"")