#!/bin/bash
# This file selectes the signal event from charged file
# To run this we need PYROOT file named "my_signal.py" in same directory
path=$(pwd)

# path of steering file code
pyrt_file="$path/my_signal_bcs.py"                      
echo "Steering file is $pyrt_file"

# Input Directory
input_path="$path/cs/test/prescale_combine"
output_path=$input_path
# Preparing string for piping in pyroot code
in2out2=""
in2out2+="$input_path/test_deltaz_filtered_charged.root "
in2out2+="$input_path/test_deltaz_filtered_signal.root "
in2out2+="$output_path/test_deltaz_filtered_charged_TM.root "
in2out2+="$output_path/test_deltaz_filtered_signal_TM.root "

# Command for scaling
# echo "$(pyroot my_signal.py $in2out1)"       #For own system
echo "$(basf2 $pyrt_file $in2out2)"        #For BELLE II server

output_file="$output_path/test_deltaz_filtered_bbbar.root"
if [ -f "$output_file" ]
then
    echo "$output_file already exists."
    $(rm -f ${output_file})
    echo "So $output_file has been deleted and a new file will be generated."
fi
input_files="$output_path/test_deltaz_filtered_charged_TM.root "
input_files+="$output_path/test_deltaz_filtered_signal_TM.root "
echo "$(hadd ${output_file} ${input_files})"