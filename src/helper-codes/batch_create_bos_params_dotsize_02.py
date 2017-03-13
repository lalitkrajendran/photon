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
filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/parameters/dotsize/'
# filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/data/bos_parameters/blur/'

# if the write directory does not exist, create it
if (not (os.path.exists(filepath))):
    os.makedirs(filepath)

# set base file name
filename_base = 'bos_simulation_parameters_'

# create array of dot sizes in pixels
# dot_size_pixels = np.array([0.1, 0.5, 1.0]) * expected_displacement_pixels
#dot_size_pixels = np.linspace(.1, .5, num=9, endpoint=True) * 10
#dot_size_pixels = np.linspace(0.05, 0.4, num=8, endpoint=True) * 10
# dot_size_pixels = np.linspace(0.0001, 0.3, num=1, endpoint=False) * 10
# dot_size_pixels = 5. # 0.001
# dot_size_pixels = np.linspace()
# print dot_size_pixels
# # set scaling factor to convert from pixels to microns (from trial and error)
# pixels_to_microns = 1e2
#
# # calculate dot size in microns
# dot_size_microns = dot_size_pixels * pixels_to_microns #* 0.001
dot_size_microns_array = np.linspace(start=200,stop=800, num=7, endpoint=True)

# array of pixel displacements
displacement_array = np.linspace(start=0.1, stop=1.0, num=10, endpoint=True)

# this is the number of image pairs that will be create for a single parameter value. this is to ensure an adequate number of
# vectors for error analysis
# number_of_image_pairs = 5
number_of_image_pairs = 5
offset = 0

# set the seeding density required (dots/32 x 32pixels)
seeding_density=10.0

for displacement_index in range(0,displacement_array.size):
    delta_x = displacement_array[displacement_index]
    delta_y = displacement_array[displacement_index]
    
    for i in range(0,dot_size_microns_array.size): #len(dot_size_microns)):
        # call function to create params for each dot size
        piv_simulation_parameters = CBS.create_bos_simulation_parameters()

        # modify the dot size parameter
        piv_simulation_parameters['bos_pattern']['grid_point_diameter'] = dot_size_microns_array[i]

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

        # # modify the minimum and and maximum X co-ordinates of the bos pattern
        # piv_simulation_parameters['bos_pattern']['X_Min'] = -360e3 #-5e4 #-2.5e4 #-25e4 #
        # piv_simulation_parameters['bos_pattern']['X_Max'] = +360e3 #+5e4 #+2.5e4 #25e4
        #
        # # modify the minimum and and maximum Y co-ordinates of the bos pattern
        # piv_simulation_parameters['bos_pattern']['Y_Min'] = -360e3 #-5e4 #-2.5e4 #-25e4
        # piv_simulation_parameters['bos_pattern']['Y_Max'] = +360e3 #+5e4 #25e4 #+2.5e4

        # set the number of pixels on the camera (where the particles will be generated)
        x_pixel_number_dots = 512
        y_pixel_number_dots = 512

        # # calculate the number of dots along x
        # # nx = 50
        # nx = seeding_density[i] * x_pixel_number
        # piv_simulation_parameters['bos_pattern']['x_grid_point_number'] = nx
        # print 'nx', nx
        #
        # # modify the number of dots along y
        # # ny = 50
        # ny = seeding_density[i] * y_pixel_number
        # piv_simulation_parameters['bos_pattern']['y_grid_point_number'] = ny
        # print 'ny', ny

        # calculate total number of particles in the image
        num_particles_total = seeding_density * x_pixel_number_dots/32 * y_pixel_number_dots/32

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

        top_write_directory = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/dotsize/'

        # for each dot size, 40 image pairs will be created. so loop through all the cases, initialize a random number to
        # ensure a different dot pattern for each image pair and also set the final directory where the images for each case
        # will be stored
        for image_pair_index in range(0,number_of_image_pairs):
            # set the location where the rendered images will be saved
            piv_simulation_parameters['output_data']['bos_pattern_image_directory'] = top_write_directory + 'delta_x=%.2f_delta_y=%.2f/%dum/%d/' % \
                                                                                      (delta_x, delta_y, dot_size_microns_array[i], offset + image_pair_index+1)            
            print "output directory", piv_simulation_parameters['output_data']['bos_pattern_image_directory']

            filename_full = filename_base + '%04d' % (displacement_index*dot_size_microns_array.size*number_of_image_pairs + i*number_of_image_pairs+image_pair_index+1) + '.mat'

            piv_simulation_parameters['density_gradients'][
                'density_gradient_filename'] = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_delta_x_%.2f_delta_y_%.2f_zmin_200_zmax_400_nz_0100.nrrd' % (delta_x, delta_y)

            # save parameters to file
            sio.savemat(filepath + filename_full, piv_simulation_parameters, appendmat=True, format='5', long_field_names=True)
