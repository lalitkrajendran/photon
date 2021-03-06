#!/usr/bin/env python2.7
# coding: utf-8
import glob
import numpy as np
import scipy.io as sio
from run_simulation_02 import run_simulation_02
import helper_functions
import os
import pickle
import time
from sys import argv
import sys
import os
import platform

if platform.system() == 'Linux':
    mount_directory = '/scratch/shannon/c/aether/'
else:
    mount_directory = '/Volumes/aether_c/'

# os.environ['LD_LIBRARY_PATH'] = '~/usr/lib/python2.7'
print(os.environ['LD_LIBRARY_PATH'])

# exit(0)
# get the starting index, and the number of parameter files to read at a time
script, filepath, starting_index, number_of_parameter_files = argv

# convert string variables to integers
starting_index = int(starting_index)
number_of_parameter_files = int(number_of_parameter_files)
print('number of parameter files to run in one simulation', number_of_parameter_files)

# # filepath = '/scratch/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/parameters/gaussian/'
# # filepath = os.path.join(mount_directory, 'Projects/BOS/tomo-bos/analysis/data/jhu/parameters')
# dir_name = os.path.dirname(os.path.abspath(__file__))
# # top level data directory
# top_data_directory = os.path.join(os.path.abspath(os.path.join(dir_name, '..')), 'sample-data', 'piv')
# filepath = os.path.join(top_data_directory, 'parameters')
# starting_index = 1
# number_of_parameter_files = 1

# specify directory where parameters are stored
bos_simulation_parameters_read_directory = filepath

# This is the list of camera parameters to run the simulation over
bos_parameter_list = sorted(glob.glob(os.path.join(bos_simulation_parameters_read_directory, '*.mat')))
print("Number of cases: %d" % len(bos_parameter_list))

# initialize the random number generator
np.random.seed()

# start the clock
start = time.time()

# This iterates through the different cases performing the image simulation
# for i in range(1, len(bos_parameter_list) + 1):
for i in range(starting_index,starting_index+number_of_parameter_files):
    parameter_index = i

    # This displays that the current camera simulation is being ran
    print("\n\n\n\n")
    print("Running case %d simulation . . . " % parameter_index)

    # This is the current filename of the camera parameter file to load
    parameter_filename_read = os.path.join(bos_simulation_parameters_read_directory, os.path.basename(bos_parameter_list[i - 1]))
    print("parameter_filename_read: " + parameter_filename_read)
    # This loads the current parameter data
    # NOTE: the mat file is being loaded as an object array of one dimension. this helps to maintain
    # similarity in syntax from the matlab code, since since matlab structures and python objects
    # have similar syntax for accessing their contents
    mat_contents = sio.loadmat(parameter_filename_read, struct_as_record=False, squeeze_me=True)
    if('simulation_paramters' in mat_contents.keys()):
        simulation_parameters = mat_contents['simulation_parameters']
    else:
        simulation_parameters = mat_contents

    # TODO temp
    #simulation_parameters = simulation_parameters['simulation_parameters']
    # convert object to dictionary
    if(type(simulation_parameters)!=dict):
        simulation_parameters = helper_functions._todict(simulation_parameters)

    # if there are nested mat_structs inside the dictionary, then convert them to dictionaries too
    for key in simulation_parameters.keys():
        if(type(simulation_parameters[key]) == sio.matlab.mio5_params.mat_struct):
             simulation_parameters[key] = helper_functions._todict(simulation_parameters[key])

    image_directory = simulation_parameters['output_data']['image_directory']

    if not platform.system() == 'Linux':
        image_directory = image_directory.replace('/scratch/shannon/c/aether/', mount_directory)
    # This creates the top level bos pattern image directory
    if not (os.path.exists(image_directory)):
        # This creates the bos pattern image directory
        os.makedirs(image_directory)

    # create the directories to save final light ray positions and directions
    if not (os.path.exists(os.path.join(image_directory, 'light-ray-positions'))):
        os.makedirs(os.path.join(image_directory, 'light-ray-positions'))
    simulation_parameters['output_data'][
        'lightray_positions_filepath'] = os.path.join(image_directory, 'light-ray-positions')

    if not (os.path.exists(os.path.join(image_directory, 'light-ray-directions'))):
        os.makedirs(os.path.join(image_directory, 'light-ray-directions'))
    simulation_parameters['output_data'][
        'lightray_directions_filepath'] = os.path.join(image_directory, 'light-ray-directions')

    # convert int variables to float
    for i in simulation_parameters:
        # if(type(simulation_parameters[i]) is str or type(simulation_parameters[i]) is list or type(simulation_parameters[i]) is unicode):
        if(type(simulation_parameters[i]) is bytes or type(simulation_parameters[i]) is list or type(simulation_parameters[i]) is str):
            continue
        for j in simulation_parameters[i]:
            if (type(simulation_parameters[i][j]) is int):
                simulation_parameters[i][j] = float(simulation_parameters[i][j])

    # This runs the camera simulation for the current camera
    run_simulation_02(simulation_parameters)

end = time.time()
print("TOTAL time taken (minutes): ", (end - start)/60)
