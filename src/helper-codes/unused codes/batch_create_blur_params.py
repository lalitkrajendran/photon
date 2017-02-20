import create_bos_simulation_parameters as CBS
import numpy as np
import scipy.io as sio
from sys import argv

# set filepath
filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/data/bos_parameters/blur/'

# set base file name
filename_base = 'bos_simulation_parameters_'

# this is the dot size in microns
dot_size_microns = 300

# this is the size of blur (pixels)
blur = 0.1

# call function to create params for each dot size
piv_simulation_parameters = CBS.create_bos_simulation_parameters()

# modify the dot size parameter
piv_simulation_parameters['bos_pattern']['grid_point_diameter'] = dot_size_microns

# modify the number of dots along x
nx = 50
piv_simulation_parameters['bos_pattern']['x_grid_point_number'] = nx

# modify the number of dots along y
ny = 50
piv_simulation_parameters['bos_pattern']['y_grid_point_number'] = ny

# modify the f# of the lens
f_number = 8
piv_simulation_parameters['lens_design']['aperture_f_number'] = f_number

# modify the file containing the density gradient data
piv_simulation_parameters['density_gradients']['density_gradient_filename'] = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/blur_%.2f_f%d.nrrd' % (blur, f_number)

# set a max size of 100 to avoid memory errors
piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point'] = 100.


# set the location where the rendered images will be saved
piv_simulation_parameters['output_data']['bos_pattern_image_directory'] = \
    '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/blur/%dx%d-f%d-dotsize%d-blur%.2f/' % (nx,ny,f_number,dot_size_microns,blur)

print "output directory", piv_simulation_parameters['output_data']['bos_pattern_image_directory']

i = 0
# create the full filename where the parameters will be saved
filename_full = filename_base + '%02d' % (i+1) + '.mat'

# save parameters to file
sio.savemat(filepath + filename_full, piv_simulation_parameters, appendmat=True, format='5', long_field_names=True)



