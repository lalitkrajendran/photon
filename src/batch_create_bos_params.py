import create_bos_simulation_parameters as CBS
import numpy as np
import scipy.io as sio
import os
from sys import argv

'''
script, expected_displacement_pixels = argv

# convert string to number
expected_displacement_pixels = int(expected_displacement_pixels)
'''

# set filepath
filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/data/bos_parameters/dot-size/'
# filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/data/bos_parameters/blur/'

# set base file name
filename_base = 'bos_simulation_parameters_'


# create array of dot sizes in pixels
# dot_size_pixels = np.array([0.1, 0.5, 1.0]) * expected_displacement_pixels
dot_size_pixels = np.linspace(.1, .5, num=9, endpoint=True) * 10

print dot_size_pixels
# set scaling factor to convert from pixels to microns (from trial and error)
pixels_to_microns = 1e2

# calculate dot size in microns
dot_size_microns = dot_size_pixels * pixels_to_microns

# # this is the dot size in microns
# dot_size_microns = 300

# this is the first derivative of density gradient
grad_x = 2.0

# this is the range of second derivatives of the density gradient to be considered
# grad_xx = np.array([15.0, 30.0, 50.0])

# this is the number of image pairs that will be create for a single dot size. this is to ensure an adequate number of
# vectors for error analysis
number_of_image_pairs = 10

# create an array of random numbers that will serve as the seed to the random number generator that in turn generates
# positions of the dots in the pattern
random_number_seed = np.random.randint(low=1, high=100000, size=(number_of_image_pairs,2))

for i in range(0,dot_size_microns.size): #len(dot_size_microns)):
    # call function to create params for each dot size
    piv_simulation_parameters = CBS.create_bos_simulation_parameters()

    # modify the dot size parameter
    piv_simulation_parameters['bos_pattern']['grid_point_diameter'] = dot_size_microns[i]

    # modify the number of dots along x
    nx = 75
    piv_simulation_parameters['bos_pattern']['x_grid_point_number'] = nx
    
    # modify the number of dots along y
    ny = 75
    piv_simulation_parameters['bos_pattern']['y_grid_point_number'] = ny

    # modify the f# of the lens
    f_number = 16
    piv_simulation_parameters['lens_design']['aperture_f_number'] = f_number
    
    ## modify the object distance
    #piv_simulation_parameters['lens_design']['object_distance'] = 680e3

    # modify the minimum and and maximum X co-ordinates of the bos pattern
    piv_simulation_parameters['bos_pattern']['X_Min'] = -2.5e4
    piv_simulation_parameters['bos_pattern']['X_Max'] = +2.5e4

    # modify the minimum and and maximum Y co-ordinates of the bos pattern
    piv_simulation_parameters['bos_pattern']['Y_Min'] = -2.5e4
    piv_simulation_parameters['bos_pattern']['Y_Max'] = +2.5e4

    # modify the file containing the density gradient data
    piv_simulation_parameters['density_gradients']['density_gradient_filename'] = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_grad_x_%02d.nrrd' % grad_x
    
    '''
    # modify the number of particle that will be used for each grid point (more particles for larger grid points)
    # the multiplying constant is the number of particles for each micron and has been determined from experience. this
    # could vary from one optical setup to the other
    piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] = int(3 * dot_size_microns[i])
    
    # set a max size of 100 to avoid memory errors
    if(piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] > 100):
        piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] = 100
    '''
    piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] = 100.


    # set the top directory which contain images of all dot sizes
    top_write_directory = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/' \
                          'dot-size/%dx%d-f%d-grad_x%.1f/' % (nx,ny,f_number,grad_x)

    # for each dot size, 40 image pairs will be created. so loop through all the cases, initialize a random number to
    # ensure a different dot pattern for each image pair and also set the final directory where the images for each case
    # will be stored
    for image_pair_index in range(0,number_of_image_pairs):
        # set the location where the rendered images will be saved
        piv_simulation_parameters['output_data']['bos_pattern_image_directory'] = top_write_directory + '%d/%d/'\
                                                                                % (dot_size_microns[i], image_pair_index+1)

        # # if the write directory does not exist, create it
        # if (not (os.path.exists(piv_simulation_parameters['output_data']['bos_pattern_image_directory']))):
        #     os.makedirs(piv_simulation_parameters['output_data']['bos_pattern_image_directory'])

        # set the random number controlling the position of dots in the dot pattern
        piv_simulation_parameters['bos_pattern']['random_number_seed'] = random_number_seed[image_pair_index,:]

        print "output directory", piv_simulation_parameters['output_data']['bos_pattern_image_directory']

        filename_full = filename_base + '%02d' % (i*number_of_image_pairs+image_pair_index+1) + '.mat'

        # save parameters to file
        sio.savemat(filepath + filename_full, piv_simulation_parameters, appendmat=True, format='5', long_field_names=True)

        # piv_simulation_parameters['output_data']['bos_pattern_image_directory'] = \
    #     '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/blur/%dx%d-f%d-dotsize%d-grad_x%d/%d' % (nx,ny,f_number,dot_size_microns,grad_x,
    #         grad_xx[i]) + '/'

    # print "output directory", piv_simulation_parameters['output_data']['bos_pattern_image_directory']



    # create the full filename where the parameters will be saved
    # filename_full = filename_base + '%02d' % (i+1) + '.mat'
    #
    # # save parameters to file
    # sio.savemat(filepath + filename_full, piv_simulation_parameters, appendmat=True, format='5', long_field_names=True)
    #
    #

