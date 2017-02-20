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
filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/data/bos_parameters/seeding/'
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
dot_size_pixels = 5. # 0.001

print dot_size_pixels
# set scaling factor to convert from pixels to microns (from trial and error)
pixels_to_microns = 1e2

# calculate dot size in microns
dot_size_microns = dot_size_pixels * pixels_to_microns #* 0.001

# # this is the dot size in microns
# dot_size_microns = 300

# this is the first derivative of density gradient
grad_x_array = np.linspace(start=0.5, stop=5.0, num=10, endpoint=True)

# this is the range of second derivatives of the density gradient to be considered
# grad_xx = np.array([15.0, 30.0, 50.0])

# this is the number of image pairs that will be create for a single parameter value. this is to ensure an adequate number of
# vectors for error analysis
# number_of_image_pairs = 5
number_of_image_pairs = 1
offset = 0

# set the seeding density required (dots/pixels). this needs to be less than 1
# seeding_density = np.array([1/64., 1/32., 1/16., 1/8.])
# seeding densities required (no. of particles/32x32 pix window)
# seeding_density = np.array([1, 4, 8, 10, 12, 16])
# seeding_density = np.linspace(start=2.0, stop=20.0, num=10, endpoint=True)
seeding_density = np.linspace(start=1.0, stop=20.0, num=20, endpoint=True)
# create an array of random numbers that will serve as the seed to the random number generator that in turn generates
# positions of the dots in the pattern
random_number_seed = np.random.randint(low=1, high=100000, size=(number_of_image_pairs, 2))


for grad_x_index in range(0,grad_x_array.size):
    grad_x = grad_x_array[grad_x_index]
    
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

        # # JHU settings
        # # modify the object distance
        # piv_simulation_parameters['lens_design']['object_distance'] = 2500e3 #680e3
        #
        # # modify the object buffer
        # piv_simulation_parameters['lens_design']['object_distance_buffer'] = 500e3 #5000e3
        #
        # field_of_view = 250e3

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
        piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] = 10

        # # % This is the number of lightrays to simulate per 'particle' (ie lightray
        # # % source point) in the calibration grid
        # piv_simulation_parameters['bos_pattern']['lightray_number_per_particle'] = 5e2

        # modify the file containing the density gradient data
        #piv_simulation_parameters['density_gradients']['density_gradient_filename'] = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_grad_x_%02d.nrrd' % grad_x
        # z_1 = 0.5e3
        # z_2 = 30e3
        # piv_simulation_parameters['density_gradients']['density_gradient_filename'] = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_grad_x_%02d_zmin_%02d_zmax_%02d.nrrd' % (grad_x, z_1/1e3, z_2/1e3)

        # piv_simulation_parameters['density_gradients']['density_gradient_filename'] = '/home/barracuda/a/lrajendr/Projects/SchlierenRayVis/Data/turbulent_channel_Mb=0.85_Reb=6900/channel_isothermal_n432237_t6.200e+02.nrrd'

        '''
        # modify the number of particle that will be used for each grid point (more particles for larger grid points)
        # the multiplying constant is the number of particles for each micron and has been determined from experience. this
        # could vary from one optical setup to the other
        piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] = int(3 * dot_size_microns[i])

        # set a max size of 100 to avoid memory errors
        if(piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] > 100):
            piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] = 100
        '''
        # piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] = 100.


        # set the top directory which contain images of all dot sizes
        #top_write_directory = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/' \
        #                      'dot-size/%dx%d-f%d-grad_x%.1f/' % (nx,ny,f_number,grad_x)

        top_write_directory = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/'

        # for each dot size, 40 image pairs will be created. so loop through all the cases, initialize a random number to
        # ensure a different dot pattern for each image pair and also set the final directory where the images for each case
        # will be stored
        for image_pair_index in range(0,number_of_image_pairs):
            # set the location where the rendered images will be saved
            # piv_simulation_parameters['output_data']['bos_pattern_image_directory'] = top_write_directory + '%d/%d/'\
            #                                                                         % (dot_size_microns[i], image_pair_index+1)
            piv_simulation_parameters['output_data']['bos_pattern_image_directory'] = top_write_directory + 'grad_x=%.2f/%d/%d/' % \
                                                                                      (grad_x, seeding_density[i], offset + image_pair_index+1)

            # # if the write directory does not exist, create it
            # if (not (os.path.exists(piv_simulation_parameters['output_data']['bos_pattern_image_directory']))):
            #     os.makedirs(piv_simulation_parameters['output_data']['bos_pattern_image_directory'])

            # set the random number controlling the position of dots in the dot pattern
            piv_simulation_parameters['bos_pattern']['random_number_seed'] = random_number_seed[image_pair_index,:]

            print "output directory", piv_simulation_parameters['output_data']['bos_pattern_image_directory']

            filename_full = filename_base + '%04d' % (grad_x_index*seeding_density.size*number_of_image_pairs + i*number_of_image_pairs+image_pair_index+1) + '.mat'

            # piv_simulation_parameters['density_gradients'][
            #     'density_gradient_filename'] = '/home/shannon/a/aether/Projects/siv/analysis/data/JHU_dataset_new/256x256x128/new_h5_reshape_nz4/nrrd/jhu_mixing_%04d.nrrd' % (image_pair_index + 1)

            piv_simulation_parameters['density_gradients'][
                'density_gradient_filename'] = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_grad_x_%.2f_zmin_200_zmax_400.nrrd' % grad_x

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

