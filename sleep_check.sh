#!/bin/bash
# if [[ $(bjobs) == "No unfinished job found" ]]
# then
#     echo "Job complete!!"
# else
#     echo "I'll sleep for 15 minutes"
#     sleep 2
#     echo "I slept 15 minutes!"
# fi
statusis=$(bjobs)
echo "this is taking a extra line"
echo "$statusis is successful"