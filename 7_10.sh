#!/bin/bash
# For testing using events like 50% as signal and 50% as continuum
sh 7_combine_for_test.sh
sh 8_bcs_deltaz_filter.sh
sh 8_bcs_a_TM.sh
sh 8_bcs_b_scale.sh
sh 9_combine_for_final_test_lumino.sh
sh 10_test.sh