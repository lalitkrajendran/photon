# coding: utf-8
import glob
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from run_piv_simulation_02 import run_piv_simulation_02
import os
import pickle

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


# This function creates a series of PIV images that closely resemble the data that was produced by the PIV Challenge Case E data.
top_write_directory = './test_directory/'

# Create the top write directory
if not (os.path.exists(top_write_directory)):
    # This creates the camera image directory
    os.mkdir(top_write_directory)

# This is the directory containing the camera simulation parameters to run for each camera
camera_simulation_parameters_read_directory = './piv_challenge_simulation_parameters/'

# This is the list of camera parameters to run the simulation over
camera_parameter_list = sorted(glob.glob(camera_simulation_parameters_read_directory + 'camera*.mat'))
print "Number of cameras: %d" % len(camera_parameter_list)

#### Creates Data Write Directories
# This is the name of the camera image top level directory
camera_image_top_directory = top_write_directory + 'camera_images/'
if not (os.path.exists(camera_image_top_directory)):
    # This creates the camera image directory
    os.mkdir(camera_image_top_directory)

# This is the name of the particle image top level directory
particle_image_top_directory = top_write_directory + 'camera_images/particle_images/'
# This creates the top level particle image directory
if not (os.path.exists(particle_image_top_directory)):
    # This creates the particle image directory
    os.mkdir(particle_image_top_directory)

# This creates the particle image data directories
for i in range(1, len(camera_parameter_list) + 1):
    # This is the current subdirecotry name
    current_subdirectory = particle_image_top_directory + 'camera_' + "%02d" % i + '/'
    # This creates the current camera particle image directory
    if not (os.path.exists(current_subdirectory)):
        # This creates the particle image camera directory
        os.mkdir(current_subdirectory)

# This is the name of the calibration image top level directory
calibration_image_top_directory = top_write_directory + 'camera_images/calibration_images/'
# This creates the top level calibration image directory
if not (os.path.exists(calibration_image_top_directory)):
    # This creates the calibration image directory
    os.mkdir(calibration_image_top_directory)

# This creates the calibration image data directories
for i in range(1, len(camera_parameter_list) + 1):
    # This is the current subdirecotry name
    current_subdirectory = calibration_image_top_directory + 'camera_' + "%02d" % i + '/'
    # This creates the current camera calibration image directory
    if not (os.path.exists(current_subdirectory)):
        # This creates the calibration image camera directory
        os.mkdir(current_subdirectory)

### Creates Particle Position Data
# This is a uniform velocity field to simulate
X_Velocity = 01.0e3
Y_Velocity = 01.0e3
Z_Velocity = 00.1e3

# This is the domain of the particles to be simulated
X_Min = -7.5e4
X_Max = +7.5e4
Y_Min = -7.5e4
Y_Max = +7.5e4
Z_Min = -7.5e3
Z_Max = +7.5e3

# This is the number of particles to simulate
total_particle_number = 100000
# total_particle_number = 10000

# initialize random number seed
RAND_SEED = 5
np.random.seed(RAND_SEED)

# This generates the particle positions
np.random.seed(5)
X = (X_Max - X_Min) * np.random.rand(int(total_particle_number), 1) + X_Min
np.random.seed(22)
Y = (Y_Max - Y_Min) * np.random.rand(int(total_particle_number), 1) + Y_Min
np.random.seed(51)
Z = (Z_Max - Z_Min) * np.random.rand(int(total_particle_number), 1) + Z_Min

# This generates the first frame particle positions
X1 = X - X_Velocity / 2.0
Y1 = Y - Y_Velocity / 2.0
Z1 = Z - Z_Velocity / 2.0
# Z1 = Z

# This generates the second frame particle positions
X2 = X + X_Velocity / 2.0
Y2 = Y + Y_Velocity / 2.0
Z2 = Z + Z_Velocity / 2.0
# Z2 = Z

# This creates the directory to save the particle data positions in
particle_position_data_directory = top_write_directory + 'particle_positions/'
# This creates the particle position data directory
if not (os.path.exists(particle_position_data_directory)):
    # This creates the particle position directory
    os.mkdir(particle_position_data_directory)

# This renames the first frame particle position data
X = X1
Y = Y1
Z = Z1

# This writes the first frame particle position data
particle_data_frame_0001 = {'X': X, 'Y': Y, 'Z': Z}
pickle.dump(particle_data_frame_0001, open(particle_position_data_directory + 'particle_data_frame_0001.p', 'wb'))
sio.savemat(particle_position_data_directory + 'particle_data_frame_0001.mat', particle_data_frame_0001)

# This renames the second frame particle position data
X = X2
Y = Y2
Z = Z2

# This writes the second frame particle position data
particle_data_frame_0002 = {'X': X, 'Y': Y, 'Z': Z}
pickle.dump(particle_data_frame_0002, open(particle_position_data_directory + 'particle_data_frame_0002.p', 'wb'))
sio.savemat(particle_position_data_directory + 'particle_data_frame_0002.mat', particle_data_frame_0002)

# create nrrd file
# nrrd_filename = os.path.dirname(os.path.realpath(__file__)) + '/' + 'test.nrrd'
# create_nrrd(nrrd_filename)


### Performs Camera Simulation
# This iterates through the different cameras performing the image simulation
for i in range(1, len(camera_parameter_list) + 1):
    # if(i>1):
    #   break
    camera_index = i
    # This displays that the current camera simulation is being ran
    print "\n\n\n\n"
    print "Running camera %d simulation . . . " % camera_index

    # This is the current filename of the camera parameter file to load
    parameter_filename_read = camera_simulation_parameters_read_directory + os.path.basename(
        camera_parameter_list[i - 1])
    print "parameter_filename_read: " + parameter_filename_read
    # This loads the current parameter data
    # NOTE: the mat file is being loaded as an object array of one dimension. this helps to maintain
    # similarity in syntax from the matlab code, since since matlab structures and python objects
    # have similar syntax for accessing their contents
    mat_contents = sio.loadmat(parameter_filename_read, struct_as_record=False, squeeze_me=True)
    piv_simulation_parameters = mat_contents['piv_simulation_parameters']

    # This changes the directory containing the particle locations in the
    # parameters structure

    piv_simulation_parameters.particle_field.data_directory = particle_position_data_directory;

    print "x_camera_angle: %f" % piv_simulation_parameters.camera_design.x_camera_angle

    # This changes the vector giving the frames of particle positions to load in
    # the parameters structure (this indexes into the list generated by the
    # command 'dir([data_directory,data_filename_prefix,'*.mat'])')
    piv_simulation_parameters.particle_field.frame_vector = np.linspace(1, 2, 2).astype('int')
    # This changes the number of particles to simulate out of the list of possible
    # particles (if this number is larger than the number of saved particles,
    # an error will be returned)
    piv_simulation_parameters.particle_field.particle_number = 100000
    # piv_simulation_parameters.particle_field.particle_number = 10000

    # This changes the directory to save the particle images in parameters structure
    piv_simulation_parameters.output_data.particle_image_directory = particle_image_top_directory + 'camera_%02d' % camera_index + '/'

    # This changes the directory to save the calibration grid images in
    # parameters structure
    piv_simulation_parameters.output_data.calibration_grid_image_directory = calibration_image_top_directory + 'camera_%02d' % camera_index + '/'

    # This runs the camera simulation for the current camera
    piv_simulation_parameters_temp = piv_simulation_parameters
    # convert object to dictionary
    piv_simulation_parameters = _todict(piv_simulation_parameters)

    # convert int variables to float
    for i in piv_simulation_parameters:
        for j in piv_simulation_parameters[i]:
            if (type(piv_simulation_parameters[i][j]) is int):
                piv_simulation_parameters[i][j] = float(piv_simulation_parameters[i][j])

    piv_simulation_parameters['particle_field']['beam_propogation_vector'] = \
    piv_simulation_parameters['particle_field']['beam_propogation_vector'].astype('float')

    # # save dictionary to mat-file
    sio.savemat('mat_files/piv_simulation_parameters.mat', {'piv_simulation_parameters': piv_simulation_parameters},
                long_field_names=True)

    # start camera simulation
    run_piv_simulation_02(piv_simulation_parameters)
