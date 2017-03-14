# coding: utf-8
import glob
import numpy as np
import scipy.io as sio
from run_piv_simulation_02 import run_piv_simulation_02
import os
import pickle
import time
from sys import argv

os.environ['LD_LIBRARY_PATH'] = '~/usr/lib/python2.7'

# The following functions convert an object to a struct so that it can be saved to a mat file
def loadmat(filename):
    '''
	this function should be called instead of direct spio.loadmat
	as it cures the problem of not properly recovering python dictionaries
	from mat files. It calls the function check keys to cure all entries
	which are still mat-objects
	'''
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):
    '''
	checks if entries in dictionary are mat-objects. If yes
	todict is called to change them to nested dictionaries
	'''
    for key in dict:
        if isinstance(dict[key], sio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):
    '''
	A recursive function which constructs from matobjects nested dictionaries
	'''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, sio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

# #get the startin index, and the number of parameter files to read at a time
# script, filepath, starting_index, number_of_parameter_files = argv
#
# # convert string variables to integers
# starting_index = int(starting_index)
# number_of_parameter_files = int(number_of_parameter_files)

filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/parameters/focusing/'
starting_index = 1
number_of_parameter_files = 1

# specify directory where parameters are stored
# bos_simulation_parameters_read_directory = '../data/bos_parameters/magnification/'
bos_simulation_parameters_read_directory = filepath

# This is the list of camera parameters to run the simulation over
bos_parameter_list = sorted(glob.glob(bos_simulation_parameters_read_directory + '*.mat'))
print "Number of cases: %d" % len(bos_parameter_list)

# initialize the random number generator
np.random.seed()

# start the clock
start = time.time()

# This iterates through the different cases performing the image simulation
# for i in range(1, len(bos_parameter_list) + 1):
for i in range(starting_index,starting_index+number_of_parameter_files):
    parameter_index = i

    # This displays that the current camera simulation is being ran
    print "\n\n\n\n"
    print "Running case %d simulation . . . " % parameter_index

    # This is the current filename of the camera parameter file to load
    parameter_filename_read = bos_simulation_parameters_read_directory + os.path.basename(
        bos_parameter_list[i - 1])
    print "parameter_filename_read: " + parameter_filename_read

    # This loads the current parameter data
    # NOTE: the mat file is being loaded as an object array of one dimension. this helps to maintain
    # similarity in syntax from the matlab code, since since matlab structures and python objects
    # have similar syntax for accessing their contents
    mat_contents = sio.loadmat(parameter_filename_read, struct_as_record=False, squeeze_me=True)
    if('piv_simulation_paramters' in mat_contents.keys()):
        piv_simulation_parameters = mat_contents['piv_simulation_parameters']
    else:
        piv_simulation_parameters = mat_contents

    # TODO temp
    #piv_simulation_parameters = piv_simulation_parameters['piv_simulation_parameters']
    # convert object to dictionary
    if(type(piv_simulation_parameters)!=dict):
        piv_simulation_parameters = _todict(piv_simulation_parameters)

    # if there are nested mat_structs inside the dictionary, then convert them to dictionaries too
    for key in piv_simulation_parameters.keys():
        if(type(piv_simulation_parameters[key]) == sio.matlab.mio5_params.mat_struct):
             piv_simulation_parameters[key] = _todict(piv_simulation_parameters[key])

    bos_pattern_image_top_directory = piv_simulation_parameters['output_data']['bos_pattern_image_directory']

    # This creates the top level bos pattern image directory
    if not (os.path.exists(bos_pattern_image_top_directory)):
        # This creates the bos pattern image directory
        os.makedirs(bos_pattern_image_top_directory)

    # create the directories to save final light ray positions and directions
    if not (os.path.exists(bos_pattern_image_top_directory + 'light-ray-positions/')):
        os.makedirs(bos_pattern_image_top_directory + 'light-ray-positions/')

    if not (os.path.exists(bos_pattern_image_top_directory + 'light-ray-directions/')):
        os.makedirs(bos_pattern_image_top_directory + 'light-ray-directions/')

    # convert int variables to float
    for i in piv_simulation_parameters:
        if(type(piv_simulation_parameters[i]) is str or type(piv_simulation_parameters[i]) is list):
            continue
        for j in piv_simulation_parameters[i]:
            if (type(piv_simulation_parameters[i][j]) is int):
                piv_simulation_parameters[i][j] = float(piv_simulation_parameters[i][j])

    # This runs the camera simulation for the current camera
    run_piv_simulation_02(piv_simulation_parameters)

end = time.time()
print "TOTAL time taken (minutes): ", (end - start)/60
