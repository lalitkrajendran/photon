from create_bos_simulation_parameters import create_bos_simulation_parameters
import numpy as np
import scipy.io as sio

# set filepath
filepath = '/home/barracuda/a/lrajendr/Projects/data/bos_parameters/dot_size/'

# set base file name
filename_base = 'bos_simulation_parameter_'

# set expected displacement
expected_displacement_pixels = 10.0

# create array of dot sizes in pixels
dot_size_pixels = np.array([0.1, 0.5, 1.0]) * expected_displacement_pixels

# set scaling factor to convert from pixels to microns (from trial and error)
pixels_to_microns = 1e2

# calculate dot size in microns
dot_size_microns = dot_size_pixels * pixels_to_microns

for i in range(0,len(dot_size_microns)):
    # call function to create params for each dot size
    piv_simulation_parameters = create_bos_simulation_parameters()

    # modify the dot size parameter
    piv_simulation_parameters['bos_pattern']['grid_point_diameter'] = dot_size_microns[i]

    # set the location where the rendered images will be save
    piv_simulation_parameters['bos_pattern']['output_data']['bos_pattern_image_directory'] = \
        '/home/barracuda/a/lrajendr/Projects/results/images/bos/error-analysis/dot-size/' + '%02d' % dot_size_microns[i] '/'

    # create the full filename where the parameters will be saved
    filename_full = filename_base + '%02d' % (i+1) + '.mat'

    # save parameters to file
    sio.savemat(filepath + filename_full, piv_simulation_parameters, appendmat=True, format='5', long_field_names=True)



