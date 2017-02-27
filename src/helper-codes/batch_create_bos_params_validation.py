# this script creates parameters for BOS simulations with varying image plane perturbations.
# purpose is to find the correct image perturbation so that the dot pattern is perfectly in focus.

import create_bos_simulation_parameters as CBS
import numpy as np
import scipy.io as sio
import sys
import os
from sys import argv


sys.path.append('/home/barracuda/a/lrajendr/Projects/camera_simulation/src/')

# set filepath
# filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/data/bos_parameters/focusing/'
filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/parameters/validation/without-error/'

# if the write directory does not exist, create it
if (not (os.path.exists(filepath))):
    os.makedirs(filepath)

# set base file name
filename_base = 'bos_simulation_parameters_'

# set dot size in microns
dot_size_microns = 5e2

# this is the first derivative of density gradient
grad_x_array = np.linspace(start=0.5, stop=5.0, num=10, endpoint=True)
# grad_x = 5.0

# this is the number of image pairs that will be create for a single parameter value. this is to ensure an adequate number of
# vectors for error analysis
# number_of_image_pairs = 5
number_of_image_pairs = 1
offset = 0

# set the seeding density required (no. of dots per 32x32 window)
seeding_density = 10.0

# create an array of random numbers that will serve as the seed to the random number generator that in turn generates
# positions of the dots in the pattern
random_number_seed = np.random.randint(low=1, high=100000, size=(number_of_image_pairs, 2))

for grad_x_index in range(0,len(grad_x_array)):
    
    # set density gradient
    grad_x = grad_x_array[grad_x_index]
    
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

    # modify the minimum and and maximum X co-ordinates of the bos pattern
    piv_simulation_parameters['bos_pattern']['X_Min'] = - image_field_of_view / 2
    piv_simulation_parameters['bos_pattern']['X_Max'] = + image_field_of_view / 2

    piv_simulation_parameters['bos_pattern']['Y_Min'] = - image_field_of_view / 2
    piv_simulation_parameters['bos_pattern']['Y_Max'] = + image_field_of_view / 2

    # set the number of pixels on the camera (where the particles will be generated)
    x_pixel_number_dots = 512
    y_pixel_number_dots = 512


    # calculate total number of particles in the image
    num_particles_total = seeding_density * x_pixel_number_dots/32 * y_pixel_number_dots/32

    # calculate the number of dots along x
    nx = np.sqrt(num_particles_total).astype('int')
    piv_simulation_parameters['bos_pattern']['x_grid_point_number'] = nx
    print 'nx', nx

    # modify the number of dots along y
    ny = np.sqrt(num_particles_total).astype('int')
    piv_simulation_parameters['bos_pattern']['y_grid_point_number'] = ny
    print 'nx', ny

    # modify the f# of the lens
    f_number = 22
    piv_simulation_parameters['lens_design']['aperture_f_number'] = f_number

    # modify lens focal length
    piv_simulation_parameters['lens_design']['focal_length'] = 200e3 #105e3

    # modify the number of 'particles' (ie lightray source points) per grid point
    piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] = 10

    # add a perturbation to the image distance
    # piv_simulation_parameters['lens_design']['perturbation'] = image_plane_perturbation
    
    # top_write_directory = '/home/shannon/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/'
    top_write_directory = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/validation/without-error/'

    # for each dot size, 40 image pairs will be created. so loop through all the cases, initialize a random number to
    # ensure a different dot pattern for each image pair and also set the final directory where the images for each case
    # will be stored
    for image_pair_index in range(0,number_of_image_pairs):
        # set the location where the rendered images will be saved
        # piv_simulation_parameters['output_data']['bos_pattern_image_directory'] = top_write_directory + '%d/%d/'\
        #                                                                         % (dot_size_microns[i], image_pair_index+1)
        piv_simulation_parameters['output_data']['bos_pattern_image_directory'] = top_write_directory + 'grad_x=%.2f/%d/' % \
                                                                                  (grad_x, offset + image_pair_index+1)

        # set the random number controlling the position of dots in the dot pattern
        piv_simulation_parameters['bos_pattern']['random_number_seed'] = random_number_seed[image_pair_index,:]

        print "output directory", piv_simulation_parameters['output_data']['bos_pattern_image_directory']

        filename_full = filename_base + '%04d' % (grad_x_index*number_of_image_pairs + image_pair_index+1) + '.mat'

        piv_simulation_parameters['density_gradients'][
            'density_gradient_filename'] = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_grad_x_%.2f_zmin_200_zmax_400.nrrd' % grad_x

        print filename_full
        
        # save parameters to file
        sio.savemat(filepath + filename_full, piv_simulation_parameters, appendmat=True, format='5', long_field_names=True)

        print "output directory", piv_simulation_parameters['output_data']['bos_pattern_image_directory']

