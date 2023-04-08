#!/bin/bash
# This file combines files for training(mailny used for scaling)
# WHEN WE PUT "*" IN "ls" COMMAND, IT RETURNS FULL PATHNAME ALSO.
# SO WE DON'T NEED TO RUN "FOR" LOOP TO ADD PATH NAME
path=$(pwd)

# Create Output Directory to save combined root files
output_path="${path}/cs/test/prescale_combine"
if [ -d "$output_path" ]
then
    echo "$output_path already exists."
else 
    $(mkdir -p ${output_path})
    echo "$output_path has been created"
fi

##the maximum number of root file to be combined at once
let nRootFile=300

#combining for different folder
declare -a options=("charged" "uubar" "ddbar" "ssbar" "ccbar" "signal" "mixed") # 

for opt in "${options[@]}"
do
    # Ouput file name
    output_file="$output_path/test_$opt.root"
    if [ -f "$output_file" ]
    then
        echo "$output_file already exists."
        $(rm -f ${output_file})
        echo "So $output_file has been deleted."
    fi

    # Input Directory
    input_path="$path/cs/test/$opt"
    input_files=$(ls ${input_path}/*.root)
    input_file_count=$(ls ${input_path}/*.root | wc -l)
    echo "Number of input files in $input_path is $input_file_count"

    # Combining files
    # if [[ $opt == "signal" ]]
    # then
    #     echo "$(cp $input_files $output_file)" 
    #             #it's just copying a file from one folder to another folder with different name
    #             # For multiple input file it faces problem. So "hadd" is safe option
    #     echo "Signal MC file $input_files  has been copied to $output_file"
    # else
    echo "Please wait for combination of $nRootFile input file from $input_path folder....."

    if [ $input_file_count -lt $nRootFile ]
    then
        echo "In $input_path folder the total number($input_file_count) of files to be combine is less than $nRootFile."
        echo "So hadd is being applying just once."
        echo "$(hadd ${output_file} ${input_files})"
    else
        combine=""  #Empty array to store files
        intermediate_files=""   # To store the name of intermediate combinded files
        let nRootFilesFlag=0    # To take counts and combine when it takes value n*$nRootFile
        let newFileNo=0;
        # Combining files in a bunch containing 300 files
        for l in ${input_files}
        do 
            combine+="$l "
            ((nRootFilesFlag++))

            if [ $nRootFilesFlag -eq $nRootFile ]
            then
                ((newFileNo++))
                temp="${output_path}/test_${opt}_${newFileNo}.root"
                echo "$(hadd ${temp} ${combine})"
                intermediate_files+="$temp "
                nRootFilesFlag=0
                combine=""
            fi
        done

        ##Leftover root file
        if [ -z "$combine" ]        # Checking whether the string is empty
        then
            echo "No root file is leftover after running for loop $newFileNo*$nRootFile times"
        else
            ((newFileNo++))
            echo "Please wait for combination of last bunch of input file from $input_path folder....."
            temp="${output_path}/test_${opt}_${newFileNo}.root"
            echo "$(hadd ${temp} ${combine})"
            intermediate_files+="$temp "
        fi


        # Combining all files
        echo "$(hadd ${output_file} ${intermediate_files})" 
                                        # the 2nd arfument of "hadd" should be "STRING" only
        echo "Deleting $intermediate_files ......"
        $(rm -rf ${intermediate_files})
    fi
    # fi
done