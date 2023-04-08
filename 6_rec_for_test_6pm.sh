#!/bin/bash
# takes command as: ./6_rec_for_test.sh 481 801 256 244 801
# this file reconstruct from all the folder
# To run this we need steering file named "b2d_test.py" in same directory
# WHEN WE PUT "*" IN "ls" COMMAND, IT RETURNS FULL PATHNAME ALSO.
# But we need "for" loop to start reconstructing from perticular number

path=$(pwd)

# path of steering file code
steering_file="$path/b2d_test_bcs_signal_1pi.py"                      
# steering_file="$path/b2d_test_bcs.py"                      
# steering_file="$path/b2d_test_fbdt_output_with_other_variable.py"                      
echo "Steering file is $steering_file"

# path of  weight file
weightfile="$path/cs/test/MVAFastBDT.root"
echo "Path of weightfile is $weightfile"


let rec_start_charged=${1:-1}-1        #as for loop starts from next numbering
let rec_start_uubar=${2:-1}-1
let rec_start_ddbar=${3:-1}-1
let rec_start_ssbar=${4:-1}-1
let rec_start_ccbar=${5:-1}-1
let rec_start_signal=${6:-1}-1
let rec_start_mixed=${7:-1}-1
let rec_start=0

# Strating reconstruction
declare -a options=("signal" "charged" "uubar" "ddbar" "ssbar" "ccbar" "mixed") # "mixed"

for opt in "${options[@]}"
do
    # path of output files
    output_folder_name="$path/cs/test/$opt"        
    if [ -d "$output_folder_name" ]
    then
        echo "$output_folder_name already exists."
    else 
        $(mkdir -p ${output_folder_name})
        echo "$output_folder_name has been created"
    fi
    echo "Name of the folder of output file is $output_folder_name"

    # Setting different input folder for 
    # Separately reconstructing from signal MC and generic MC
    # path of input files and     # Reading file name in given folder
    if [[ $opt == "signal" ]]
    then
        input_datapath='/home/belle2/ssana/Official_signalMC'
        input_folder_name="$input_datapath/sub00"
        echo "Name of input folder is $input_folder_name"

        input_filelist=$(ls ${input_folder_name}/*.root)
        input_file_count=$(ls ${input_folder_name}/*.root | wc -l)
        echo "Number of input file in $input_folder_name is $input_file_count"
    elif [[ $opt == "uubar" || $opt == "ccbar" ]]
    then
        input_datapath='/group/belle2/dataprod/MC/MC15ri'
        input_folder_name_1="$input_datapath/$opt/sub00"
        input_folder_name_2="$input_datapath/$opt/sub01"
        echo "Name of input folder are $input_folder_name_1 and $input_folder_name_2"

        input_filelist=$(ls ${input_folder_name_1}/*.root)
        input_filelist+=" "
        input_filelist+=$(ls ${input_folder_name_2}/*.root)

        input_file_count_1=$(ls ${input_folder_name_1}/*.root | wc -l)
        input_file_count_2=$(ls ${input_folder_name_2}/*.root | wc -l)
        input_file_count=$(($input_file_count_1 + $input_file_count_2))
        echo "Number of input file in $input_folder_name_1 and $input_folder_name_2 is $input_file_count"

    else
        input_datapath='/group/belle2/dataprod/MC/MC15ri'
        input_folder_name="$input_datapath/$opt/sub00"
        echo "Name of input folder is $input_folder_name"

        input_filelist=$(ls ${input_folder_name}/*.root)
        input_file_count=$(ls ${input_folder_name}/*.root | wc -l)
        echo "Number of input file in $input_folder_name is $input_file_count"
    fi
    
    ###########################################################
    # Reading file name in given folder
    # input_filelist=$(ls ${input_folder_name}/*.root)
    # input_file_count=$(ls ${input_folder_name}/*.root | wc -l)
    # echo "Number of input file in $input_folder_name is $input_file_count"

    #Storing filename in array from long single string
    declare -a input_file_array=()
    let count=0
    for l in ${input_filelist}
    do
        input_file_array[$count]="$l"
        ((count++))
    done
    length_input_file_array=${#input_file_array[@]}

    #checking whether array dimension is ok
    if [[ $input_file_count == $count && $input_file_count == $length_input_file_array ]]
    then
        echo "Input file count and input file array dimension are same."
    else
        echo "There is some discrepancy between Input file count and input file array dimension"
        echo "Because Input file count is $input_file_count but array dimension is $length_input_file_array"
    fi

    # Finalizing the starting number of array for reconstruction
    if [[ $opt == "charged" ]]
    then
        rec_start=$rec_start_charged
    elif [[ $opt == "mixed" ]]
    then
        rec_start=$rec_start_mixed
    elif [[ $opt == "uubar" ]]
    then
        rec_start=$rec_start_uubar
    elif [[ $opt == "ddbar" ]]
    then
        rec_start=$rec_start_ddbar
    elif [[ $opt == "ssbar" ]]
    then
        rec_start=$rec_start_ssbar
    elif [[ $opt == "ccbar" ]]
    then
        rec_start=$rec_start_ccbar
    elif [[ $opt == "signal" ]]
    then
        rec_start=$rec_start_signal
    fi

    # Submiting Jobs
    let i=0
    for((i=$rec_start; i<$length_input_file_array; i++))
    do
        # if [[ $opt == "signal" && $i = "5" ]]
        # then
        #     continue    # we will keep 5th file for test and {(1-4)&6} files are for train
        # fi
        echo "$(bsub -q s "basf2 $steering_file ${input_file_array[$i]} $output_folder_name $weightfile")"
    done
done