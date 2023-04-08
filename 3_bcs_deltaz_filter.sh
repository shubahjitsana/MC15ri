#!/bin/bash
# This file scales the event from uubar,ddbar,ssbar,ccbar file
# To run this we need PYROOT file named "my_scale.py" in same directory
path=$(pwd)

# path of steering file code
pyrt_file="$path/my_deltaz_filter_with_veto.py"                      
# pyrt_file="$path/my_deltaz_filter.py"                      
echo "Steering file is $pyrt_file"

# Input Directory
input_path="$path/cs/train"
output_path=$input_path

# Preparing string for piping in pyroot code
declare -a options=("uubar" "ddbar" "ssbar" "ccbar" "charged" "signal")
in6out6=""
for opt in "${options[@]}"
do
    in6out6+="$input_path/train_$opt.root "
done
for opt in "${options[@]}"
do
    in6out6+="$output_path/train_deltaz_filtered_$opt.root "
done

# Command for scaling
# echo "$(pyroot my_scale.py $in6out6)"         #For own system
echo "$(basf2 $pyrt_file $in6out6)"          #For BELLE II server

# Deleting files containing deltaz value as "nan"
for opt in "${options[@]}"
do
    file_to_be_deleted="$input_path/train_$opt.root"
    $(rm -f ${file_to_be_deleted})
    echo "So $file_to_be_deleted has been deleted."

done