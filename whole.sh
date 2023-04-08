#!/bin/bash
# this script compiles the whole code of continuum suppresion
# But I'm fetching problem to pass output of "bjobs" command in "if" statement


sh 1_rec_all.sh

# Starting training when all data is reconstructed
echo "Please wait at least 2 hours for completion of all data reconstruction."
echo "Meanwhile please avoid submitting another job to bsub."
sleep 2h
for i in {1..40}
do
    if [[ $(bjobs) == "No unfinished job found" ]]
    then
        echo "All data have been reconstruted."
        # For training using events like 50% of signal and 50% of continuum
        sh 2_combine_for_train_6pm.sh
        sh 3_bcs_deltaz_filter.sh
        sh 3_bcs_a_signal_only.sh
        sh 3_bcs_b_scale.sh
        sh 4_combine_for_final_train_lumino.sh
        sh 5_train.sh
        sh 6_rec_for_test_6pm.sh
        break
    else
        echo "All data haven't been reconstruted yet."
        echo "Please wait another 20 minutes for completion of all data reconstruction."
        echo "Meanwhile please avoid submitting another job to bsub."
        sleep 20m
    fi
done

# Starting testing when testing data is reconstructed
echo "Please wait at least 15 minutes for completion of testing data reconstruction."
echo "Meanwhile please avoid submitting another job to bsub."
sleep 15m
for i in {1..40}
do
    if [[ $(bjobs) == "No unfinished job found" ]]
    then
        echo "All testing data have been reconstruted."
        # For testing using events like 50% of signal and 50% of continuum
        sh 7_combine_for_test.sh
        sh 8_bcs_deltaz_filter.sh
        sh 8_bcs_a_signal_only.sh
        sh 8_bcs_b_scale.sh
        sh 9_combine_for_final_test_lumino.sh
        sh 10_test.sh
        break
    else
        echo "All testing data haven't been reconstruted yet."
        echo "Please wait another 5 minutes for completion of testing data reconstruction."
        echo "Meanwhile please avoid submitting another job to bsub."
        sleep 5m
    fi
done