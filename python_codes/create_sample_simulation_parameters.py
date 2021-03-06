# this script creates parameters for BOS simulations with varying displacements to assess light ray pdfs


import numpy as np
import scipy.io as sio
import sys
import os
from sys import argv
import glob
import platform

if platform.system() == 'Linux':
    mount_directory = '/scratch/shannon/c/aether/'
else:
    mount_directory = '/Volumes/aether_c/'

# import path containing python modules
sys.path.append(os.path.join(mount_directory, 'Projects/BOS/general-codes/python-codes'))
import modify_plt_settings
import loadmat_functions

# import path containing image generation code
sys.path.append(os.path.join(mount_directory, 'Projects/BOS/image-generation/analysis/src/photon/python_codes'))
import create_simulation_parameters as CS

def create_write_directories(top_data_directory):
    # top level image directory
    top_image_directory = os.path.join(top_data_directory, 'images')
    # directory to save parameters
    top_parameter_directory = os.path.join(top_data_directory, 'parameters')
    # if the write directory does not exist, create it
    if (not (os.path.exists(top_parameter_directory))):
        os.makedirs(top_parameter_directory)
    # if it exists, then delete the old parameter files
    else:
        files = glob.glob(os.path.join(top_parameter_directory, '*.mat'))
        for f in files:
            os.remove(f)

    return top_image_directory, top_parameter_directory

def main():
    '''
    script, expected_displacement_pixels = argv

    # convert string to number
    expected_displacement_pixels = int(expected_displacement_pixels)
    '''
    # simulation type ('piv' or 'bos' or 'calibration')
    simulation_type = 'bos'

    # ==========================
    # read/write settings
    # ==========================
    dir_name = os.path.dirname(os.path.abspath(__file__))
    # top level data directory
    top_data_directory = os.path.join(os.path.abspath(os.path.join(dir_name, '..')), 'sample-data', simulation_type)
    # create write directory
    top_image_directory, top_parameter_directory = create_write_directories(top_data_directory)

    # ==========================
    # image generation settings
    # ==========================

    # generate a sample data structure
    simulation_parameters = CS.create_simulation_parameters(simulation_type)

    # simulation_parameters['particle_field']['perform_mie_scattering'] = False
    
    # set file containing density gradients if it is a bos simulatoin
    if simulation_type == 'bos':
        simulation_parameters['density_gradients']['density_gradient_filename'] = os.path.join(top_data_directory, 'sample-density.nrrd')

    # set the final directory where the images for each case will be stored
    simulation_parameters['output_data']['image_directory'] = top_image_directory
    print("image directory", simulation_parameters['output_data']['image_directory'])

    # parameter filename
    filename_full = os.path.join(top_parameter_directory, 'sample-parameters.mat')
    print('filename', filename_full)

    # save parameters to file
    sio.savemat(filename_full, simulation_parameters, appendmat=True, format='5', long_field_names=True)

if __name__ == "__main__":
    main()

