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
filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/data/bos_parameters/spark/'
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
# dot_size_pixels = np.linspace(0.5, 5, num=10, endpoint=True)
dot_size_pixels = 2.5

print dot_size_pixels
# set scaling factor to convert from pixels to microns (from trial and error)
pixels_to_microns = 1e2

# calculate dot size in microns
dot_size_microns = dot_size_pixels * pixels_to_microns #* 0.001
print 'dot_size_microns', dot_size_microns

# # this is the dot size in microns
# dot_size_microns = 300

# # this is the first derivative of density gradient
grad_x = 5.0

# this is the number of image pairs that will be create for a single parameter value. this is to ensure an adequate number of
# vectors for error analysis
number_of_image_pairs = 1

# set the seeding density required (dots/32x32 pixels)
seeding_density = 10.0

# set the range of focal lengths (microns) to consider
focal_length = 105e3

# create an array of random numbers that will serve as the seed to the random number generator that in turn generates
# positions of the dots in the pattern
random_number_seed = np.random.randint(low=1, high=100000, size=(number_of_image_pairs, 2))

# call function to create params for each dot size
piv_simulation_parameters = CBS.create_bos_simulation_parameters()

# modify the dot size parameter
piv_simulation_parameters['bos_pattern']['grid_point_diameter'] = dot_size_microns

# modify the object distance (microns)
piv_simulation_parameters['lens_design']['object_distance'] = 25 * 25.4 * 1e3

# this is the field of view of the camera
field_of_view = 10e3

# this is the field of view over which the dots will actually be generated. this is smaller than the full field of
# view to account for edge effects of the lens
image_field_of_view = field_of_view #/2

piv_simulation_parameters['bos_pattern']['X_Min'] = - image_field_of_view / 2
piv_simulation_parameters['bos_pattern']['X_Max'] = + image_field_of_view / 2

piv_simulation_parameters['bos_pattern']['Y_Min'] = - image_field_of_view / 2
piv_simulation_parameters['bos_pattern']['Y_Max'] = + image_field_of_view / 2

# set the number of pixels on the camera (where the particles will be generated)
x_pixel_number_dots = 256 #int(512 * focal_length_array[i] / 200.0e3)
y_pixel_number_dots = 256 #int(512 * focal_length_array[i] / 200.e3)

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

nx = np.sqrt(num_particles_total * 0.75).astype('int')
piv_simulation_parameters['bos_pattern']['x_grid_point_number'] = nx
print 'nx', nx

ny = np.sqrt(num_particles_total * 0.75).astype('int')
piv_simulation_parameters['bos_pattern']['y_grid_point_number'] = ny
print 'nx', ny

# modify the f# of the lens
f_number = 32
piv_simulation_parameters['lens_design']['aperture_f_number'] = f_number

# modify lens focal length
piv_simulation_parameters['lens_design']['focal_length'] = focal_length

# modify the number of 'particles' (ie lightray source points) per grid point
piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] = 100

top_write_directory = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/spark/'

# for each dot size, 40 image pairs will be created. so loop through all the cases, initialize a random number to
# ensure a different dot pattern for each image pair and also set the final directory where the images for each case
# will be stored
for image_pair_index in range(0,number_of_image_pairs):
    # set the location where the rendered images will be saved
    piv_simulation_parameters['output_data']['bos_pattern_image_directory'] = top_write_directory + 'grad_x=%.2f/%d/%d/' % \
                                                                              (grad_x, focal_length*1e-3, image_pair_index+1)

    # set the random number controlling the position of dots in the dot pattern
    piv_simulation_parameters['bos_pattern']['random_number_seed'] = random_number_seed[image_pair_index,:]

    print "output directory", piv_simulation_parameters['output_data']['bos_pattern_image_directory']

    filename_full = filename_base + '%04d' % (image_pair_index+1) + '.mat'

    piv_simulation_parameters['density_gradients'][
        'density_gradient_filename'] = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_grad_x_%.2f_zmin_200_zmax_400.nrrd' % grad_x

    # save parameters to file
    sio.savemat(filepath + filename_full, piv_simulation_parameters, appendmat=True, format='5', long_field_names=True)
