# this script creates parameters for BOS simulations with varying seeding densities.

import create_bos_simulation_parameters as CBS
import numpy as np
import scipy.io as sio
import sys
import os
from sys import argv


sys.path.append('/home/barracuda/a/lrajendr/Projects/camera_simulation/src/')

'''
script, expected_displacement_pixels = argv

# convert string to number
expected_displacement_pixels = int(expected_displacement_pixels)
'''

# set filepath
filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/parameters/seeding/'

# if the write directory does not exist, create it
if (not (os.path.exists(filepath))):
    os.makedirs(filepath)

# set base file name
filename_base = 'bos_simulation_parameters_'

# this is the dot size in microns
dot_size_microns = 500

# array of pixel displacements
displacement_array = np.linspace(start=0.1, stop=2.0, num=20, endpoint=True)

# this is the number of image pairs that will be create for a single parameter value. this is to ensure an adequate number of
# vectors for error analysis
number_of_image_pairs = 4
offset = 0

# set the seeding density required (dots/pixels). this needs to be less than 1
seeding_density = np.linspace(start=4.0, stop=20.0, num=5, endpoint=True)

for displacement_index in range(0,displacement_array.size):
    delta_x = displacement_array[displacement_index]
    delta_y = displacement_array[displacement_index]
    
    for i in range(0,seeding_density.size): #len(dot_size_microns)):
        # call function to create params for each dot size
        piv_simulation_parameters = CBS.create_bos_simulation_parameters()

        # modify the dot size parameter
        piv_simulation_parameters['bos_pattern']['grid_point_diameter'] = dot_size_microns

        # dof settings
        # modify the object distance
        piv_simulation_parameters['lens_design']['object_distance'] = 2500e3  # 680e3

        # this is the field of view of the camera
        field_of_view = 200e3

        # this is the field of view over which the dots will actually be generated. this is smaller than the full field of
        # view to account for edge effects of the lens
        image_field_of_view = field_of_view/2

        piv_simulation_parameters['bos_pattern']['X_Min'] = - image_field_of_view / 2
        piv_simulation_parameters['bos_pattern']['X_Max'] = + image_field_of_view / 2

        piv_simulation_parameters['bos_pattern']['Y_Min'] = - image_field_of_view / 2
        piv_simulation_parameters['bos_pattern']['Y_Max'] = + image_field_of_view / 2

        # set the number of pixels on the camera (where the particles will be generated)
        x_pixel_number_dots = 512
        y_pixel_number_dots = 512

        # calculate total number of particles in the image
        num_particles_total = seeding_density[i] * x_pixel_number_dots/32 * y_pixel_number_dots/32

        nx = np.sqrt(num_particles_total).astype('int')
        piv_simulation_parameters['bos_pattern']['x_grid_point_number'] = nx
        print 'nx', nx

        ny = np.sqrt(num_particles_total).astype('int')
        piv_simulation_parameters['bos_pattern']['y_grid_point_number'] = ny
        print 'nx', ny

        # modify the f# of the lens
        f_number = 22
        piv_simulation_parameters['lens_design']['aperture_f_number'] = f_number

        # modify lens focal length
        piv_simulation_parameters['lens_design']['focal_length'] = 200e3 #105e3


        # modify the number of 'particles' (ie lightray source points) per grid point
        piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] = 100

        piv_simulation_parameters['bos_pattern']['lightray_number_per_particle'] = 1e3

        # set image noise to zero
        piv_simulation_parameters['camera_design']['image_noise'] = 0.00

        top_write_directory = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/seeding/'

        # for each dot size, several image pairs will be created. so loop through all the cases, initialize a random number to
        # ensure a different dot pattern for each image pair and also set the final directory where the images for each case
        # will be stored
        for image_pair_index in range(0,number_of_image_pairs):
            piv_simulation_parameters['output_data']['bos_pattern_image_directory'] = top_write_directory + 'delta_x=%.2f_delta_y=%.2f/%d/%d/' % \
                                                                                      (delta_x, delta_y, seeding_density[i], offset + image_pair_index+1)

            print "output directory", piv_simulation_parameters['output_data']['bos_pattern_image_directory']

            filename_full = filename_base + '%04d' % (displacement_index*seeding_density.size*number_of_image_pairs + i*number_of_image_pairs+image_pair_index+1) + '.mat'

            piv_simulation_parameters['density_gradients'][
                'density_gradient_filename'] = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_no_noise_delta_x_%.2f_delta_y_%.2f_zmin_200_zmax_400_nz_0100.nrrd' % (delta_x, delta_y)

            # save parameters to file
            sio.savemat(filepath + filename_full, piv_simulation_parameters, appendmat=True, format='5', long_field_names=True)

