#!/bin/bash

# For training using events according to overall cross section
# sh 2_combine_for_train_5pm.sh
# sh 3a_signal_only.sh
# sh 3b_scale.sh
# sh 3c_signalMC.sh
# sh 4_combine_for_final_train.sh
# sh 5_train.sh
# sh 6_rec_for_test_5pm.sh

# For training using events like 50% as signal and 50% as continuum
sh 2_combine_for_train_6pm.sh
sh 3_bcs_deltaz_filter.sh
sh 3_bcs_a_TM.sh
sh 3_bcs_b_scale.sh
sh 4_combine_for_final_train_lumino.sh
sh 5_train.sh
sh 6_rec_for_test_6pm.sh