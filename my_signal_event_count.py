#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import sys


print(f"Number of input files are {len(sys.argv)-1}")
print(f"Steering file name is : {sys.argv[0]}")
whole_signal_event_number = 0 
whole_total_event_in_charged = 0 
for i in sys.argv[1:]:
    input_file = ROOT.TFile.Open(f"{i}")
    input_tree = input_file.Get('tree')
    signal_event_number = 0 
    total_event_in_charged = input_tree.GetEntries()
    for iEvent in range(total_event_in_charged):
        input_tree.GetEntry(iEvent)

        # To find N_cand
        # signal = getattr(input_tree, 'isSignal')
        # bcs_rank = getattr(input_tree, 'BCS3_rank')
        # if ((signal==1) and ((bcs_rank==1) or (bcs_rank==2) or (bcs_rank==3))):
        #     signal_event_number += 1

        # To find number of TM events passing cuts
        deltae = getattr(input_tree, 'deltaE')
        mbc = getattr(input_tree, 'Mbc')
        md0 = getattr(input_tree, 'D0_bar_InvM')

        kid = getattr(input_tree, 'Kp_PID_bin_kaon')

        issig = getattr(input_tree, 'isSignal')
        D_s_InvM = getattr(input_tree, 'D_s_InvM')
        D_10D_md = getattr(input_tree, 'D_10D_md')
        InvM1stand2ndpi = getattr(input_tree, 'InvM1stand2ndpi')
        InvM2ndand3rdpi = getattr(input_tree, 'InvM2ndand3rdpi')
        InvMD0and1stpi = getattr(input_tree, 'InvMD0and1stpi')
        InvMD0and2ndpi = getattr(input_tree, 'InvMD0and2ndpi')
        InvMD0and3rdpi = getattr(input_tree, 'InvMD0and3rdpi')
        DstrminusroeD_md = getattr(input_tree, 'DstrminusroeD_md')
        if (
            (mbc > 5.27 and mbc < 5.29) and 
            (deltae > -0.1 and deltae < 0.1) and  
            (md0 > 1.8536 and md0 < 1.8748) and 

            (kid > 0.6) and 

            (D_s_InvM < 1.958 or D_s_InvM > 1.978) and 
            (D_10D_md < 0.496 or D_10D_md > 0.628) and 
            (InvM1stand2ndpi < 0.489 or InvM1stand2ndpi > 0.506) and 
            (InvM2ndand3rdpi < 0.49 or InvM2ndand3rdpi > 0.505) and 
            (InvMD0and1stpi < 2.0085 or InvMD0and1stpi > 2.0122) and 
            (InvMD0and2ndpi < 2.0087 or InvMD0and2ndpi > 2.0118) and 
            (InvMD0and3rdpi < 2.009 or InvMD0and3rdpi > 2.0116) and 
            (DstrminusroeD_md < 0.14402 or DstrminusroeD_md > 0.14682) #and
            # (issig ==1)
        ):
            signal_event_number += 1

    print(f"There are {total_event_in_charged} reconstructed events in {i}")
    print(f"There are {signal_event_number} signal events in {i}")
    whole_signal_event_number += signal_event_number 
    whole_total_event_in_charged += total_event_in_charged
print(f"There are total {whole_signal_event_number} signal events out of \
{whole_total_event_in_charged} reconstructed events in given list of files.")