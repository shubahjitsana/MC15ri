#!/bin/bash
# This file combines files for training
# WHEN WE PUT "*" IN "ls" COMMAND, IT RETURNS FULL PATHNAME ALSO.
# SO WE DON'T NEED TO RUN "FOR" LOOP TO ADD PATH NAME
path=$(pwd)

# Output-Input Directory to save combined root files
input_path="$path/cs/train/signal_scaled"
output_path=$input_path

output_file="$output_path/train.root"
if [ -f "$output_file" ]
then
    echo "$output_file already exists."
    $(rm -f ${output_file})
    echo "So $output_file has been deleted."
fi

# Reading input files
input_files=$(ls ${input_path}/*.root)
input_file_count=$(ls ${input_path}/*.root | wc -l)
echo "Number of input files in $input_path is $input_file_count"

# Combinig files
echo "Please wait for Combinded output file from $input_path folder....."
# echo "$(hadd ${output_path}/train.root ${path_in_string})" 
# WHEN WE PUT "*" IN "ls" COMMAND, IT RETURNS FULL PATHNAME ALSO.
echo "$(hadd ${output_file} ${input_files})" 
                                # the 2nd arfument of "hadd" should be "STRING" only