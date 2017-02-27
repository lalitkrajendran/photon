#!/usr/bin/env bash
# This program calls the ray tracing code to generate images for
# 10 parameters at a time. This is because I am unable to run the ray
# tracing code for more than 12 cases without getting a memory erorr.

# this is the parameter for error analysis
error_analysis_parameter=$1

# this is the starting case number
start_case=$2

# set directory where the parameter files are stored
#parameter_file_path=/home/barracuda/a/lrajendr/Projects/camera_simulation/data/bos_parameters/$error_analysis_parameter
parameter_file_path=/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/parameters/$error_analysis_parameter

# set directory where ray tracing code is located
src_filepath=/home/barracuda/a/lrajendr/Projects/camera_simulation/src

# read total number of parameter files
total_num_parameter_files=$(ls $parameter_file_path/*.mat | wc -l)

# display the total number of parameter files to the user
echo "total number of files:" $total_num_parameter_files

# set number of parameter files to read in call to the ray tracing code
batch_num_parameter_files=10

# display the number of parameter files that will be read at a time to the user
echo "number of files to be read in one call:" $batch_num_parameter_files

# call ray tracing code sequentially
for i in $(seq $start_case $batch_num_parameter_files $total_num_parameter_files)
do
    echo $i
    # call ray tracing code
    cd $src_filepath
    python batch_run_bos_simulation.py $parameter_file_path/ $i $batch_num_parameter_files
done
