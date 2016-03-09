import glob
import pickle
#import os
import numpy as np
import scipy.io as sio
from perform_ray_tracing_03 import perform_ray_tracing_03
#import subprocess
#import time
#from numpy import linalg as la
#import scipy.special as scispec
#import scipy.misc as scimisc
#import skimage.io as ski_io
import matplotlib.pyplot as plt
# start matlab engine
#eng = matlab.engine.start_matlab()

def create_single_lens_optical_system():
# This function is designed to create a data structure containing the
# necessary optical elements to simulate a plenoptic camera.

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #% Single Lens Design                                                      %
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #% This creates the optical data structure
    optical_system = {}

    # This string specifies the distance unit system to use for the optical 
    # design
    optical_system['distance_units'] = 'micron'

    # This string specifies the angle unit system to use for the optical 
    # design
    optical_system['angle_units'] = 'degrees'

    # This creates the top level system of elements
    optical_system['design'] = {}

    # This adds a string stating the type of the current system
    optical_system['design']['element_type'] = 'system'

    # This adds a double giving the number of elements in the current system -
    # if the current system is a single element, the element number is Null
    optical_system['design']['element_number'] = 1

    # This is a Boolean value stating whether the current elements are coplanar
    # (this value may only be true for an element type of 'system' and will be
    # Null for single optical elements)
    optical_system['design']['elements_coplanar'] = False

    # This is a double value stating the z axis distance between the current 
    # back optical axis vertex and the front optical axis vertex of the next
    # optical element
    optical_system['design']['z_inter_element_distance'] = 0.0

    # This is a 1 x 2 double precision vector giving the optical element offset
    # from the optical axis in the x and y dimensions respectively
    optical_system['design']['axial_offset_distances'] = np.array([0.0,0.0]);

    # This is a 1 x 3 vector giving the rotation angles of the current optical
    # element or system
    optical_system['design']['rotation_angles'] = np.array([0.0,0.0,0.0]);

    # This creates a substructure to describe the geometry of the curret
    # element
    optical_system['design']['element_geometry'] = {}

    # This also specifies the geometry of the current optical element - if the
    # element type is 'system', then the geometry has a Null value
#     optical_system['design']['element_geometry'] = None

    # This creates a substructure to describe the optical properties of the 
    # current element (the value of the structure will be null if the current 
    # element is a 'system')
    optical_system['design']['element_properties'] = {}

    #%% This describes the first system within the overall optical train

    # This adds a string stating the type of the current system
    optical_system['design']['optical_element'] = [dict() for x in range(optical_system['design']['element_number'])]
    optical_system['design']['optical_element'][0]['element_type'] = 'system'

    # % This adds a double giving the number of elements in the current system -
    # % if the current system is a single element, the element number is Null
    optical_system['design']['optical_element'][0]['element_number'] = 1

    # % This is a Boolean value stating whether the current elements are coplanar
    # % (this value may only be true for an element type of 'system' and will be
    # % Null for single optical elements)
    optical_system['design']['optical_element'][0]['elements_coplanar'] = False

    # % This is a double value stating the z axis distance between the current 
    # % back optical axis vertex and the front optical axis vertex of the next
    # % optical element
    optical_system['design']['optical_element'][0]['z_inter_element_distance'] = 1.0e4

    # % This is a 1 x 2 double precision vector giving the optical element offset
    # % from the optical axis in the x and y dimensions respectively
    optical_system['design']['optical_element'][0]['axial_offset_distances'] = np.array([0.0 ,0.0 ])

    # % This is a 1 x 3 vector giving the rotation angles of the current optical
    # % element or system
    optical_system['design']['optical_element'][0]['rotation_angles'] = np.array([0.0*np.pi/180.0, 0.0 * np.pi/180.0, 0.0])

    # % This creates a substructure to describe the geometry of the curret
    # % element
    optical_system['design']['optical_element'][0]['element_geometry'] = {}

    # % This specifies the geometry of the current optical element - if the
    # % element type is 'system', then the geometry has a Null value
    #optical_system['design']['optical_element'][0]['element_geometry'] = None

    # % This creates a substructure to describe the optical properties of the 
    # % current element (the value of the structure will be null if the current 
    # % element is a 'system')
    optical_system['design']['optical_element'][0]['element_properties'] = {}

    # % This specifies the optical properties of the current element - if the
    # % element type is 'system' then the geometry has a Null value
    #optical_system['design']['optical_element'][0]['element_properties'] = None

    # %% This describes the first lens in the system

    # % This creates a structure to contain the optical elements within the
    # % current system
    #optical_system['design']['optical_element'][0]['optical_element'] = {}

    # % This adds a string stating the type of the current system
    optical_system['design']['optical_element'][0]['optical_element'] = [dict() for x in range(optical_system['design']['optical_element'][0]['element_number'])]

    optical_system['design']['optical_element'][0]['optical_element'][0]['element_type'] = 'lens'

    # % This adds a double giving the number of elements in the current system -
    # % if the current system is a single element, the element number is Null
    # optical_system['design']['optical_element'][0]['optical_element'][0]['element_number'] = None
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_number'] = np.NAN

    # % This is a Boolean value stating whether the current elements are coplanar
    # % (this value may only be true for an element type of 'system' and will be
    # % Null for single optical elements)
    optical_system['design']['optical_element'][0]['optical_element'][0]['elements_coplanar'] = []

    # % This is a double value stating the z axis distance between the current 
    # % back optical axis vertex and the front optical axis vertex of the next
    # % optical element
    optical_system['design']['optical_element'][0]['optical_element'][0]['z_inter_element_distance'] = 0.0e3

    # % This is a 1 x 2 double precision vector giving the optical element offset
    # % from the optical axis in the x and y dimensions respectively
    optical_system['design']['optical_element'][0]['optical_element'][0]['axial_offset_distances'] = np.array([0.0,0.0])

    # % This is a 1 x 3 vector giving the rotation angles of the current optical
    # % element or system
    optical_system['design']['optical_element'][0]['optical_element'][0]['rotation_angles'] = np.array([0.0,0.0,0.0])

    # % This creates a substructure to describe the geometry of the current
    # % element (the value of the structure will be null if the current element
    # % is a 'system')
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry'] = {}

    # % This is a string specifying the shape of the front surface of the current
    # % optical element which must be some function of the independent variables
    # % 'x', 'y', and 'r' where r = sqrt(x^2+y^2), constant values added or
    # % subtracted from this function are ignored
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['front_surface_shape'] = '-sqrt((200e3)^2-(x.^2+y.^2))'

    # % This is a string specifying the shape of the front surface of the current
    # % optical element which must be some function of the independent variables
    # % 'x', 'y', and 'r' where r = sqrt(x^2+y^2), constant values added or
    # % subtracted from this function are ignored
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['back_surface_shape'] = '+sqrt((400e3)^2-(x.^2+y.^2))'

    # % This is a double value specifying the pitch (diameter) of the current
    # % optical element
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['pitch'] = 100.0e3

    # % This is a Boolean value stating whether the the front surface of the
    # % current optical element is spherical
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['front_surface_spherical'] = True

    # % This is a Boolean value stating whether the the back surface of the
    # % current optical element is spherical
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['back_surface_spherical'] = True

    # % This is a double value specifying the radius of curvature of the front
    # % surface of the current optical element (if the element is not spherical,
    # % then this value is Null)
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['front_surface_radius'] = +200.0e3;

    # % This is a double value specifying the radius of curvature of the back
    # % surface of the current optical element (if the element is not spherical,
    # % then this value is Null)
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['back_surface_radius'] = -400.0e3

    # % This is a double value specifying the distance between the optical axis
    # % vertices on the current optical element
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['vertex_distance'] = 10.0e3


    # % This creates a substructure to describe the optical properties of the 
    # % current element (the value of the structure will be null if the current 
    # % element is a 'system')
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties'] = {}

    # % This adds the refractive index of the current optical element - if this
    # % is a double value, the refractive index is taken to be constant, if this
    # % is a string, then the refractive index is taken to be a function of the
    # % independent variable wavelength 'lambda'
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['refractive_index'] = 1.5

    # % This adds the Abbe number dispersion constant of the current optical
    # % element (defined at lambda_D = 589.3 nm, lambda_F = 486.1 nm, and
    # % lambda_C = 656.3 nm) - this variable must be a double value.  If the
    # % refractive index is a double, then it is assummed to be measured at
    # % lambda_D = 589.3 nm and the Abbe number is used with a reduced Cauchy
    # % formula to calculate dispersion.  If the refractive index is a
    # % string, then the Abbe number is ignored.
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['abbe_number'] = []

    # % This is the equivalent thin lens focal length for the current element (if
    # % this value is defined - if not, then the value is Null)
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['thin_lens_focal_length'] = 85.0e3

    # % This is the ratio of the intensity of the transmitted light to the
    # % incident light - this can be a double to specify a constant value or a
    # % function of 'x', 'y', and 'r' to specify a varying ratio
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['transmission_ratio'] = 1.0

    # % This is the rate of absorbance of the incident light rays (and has
    # % units of radiance/length - ie the total transmission ratio is given
    # % by one minus the integral of the absorbance rate over the path of
    # % the light ray) - this can be a double to specify a constant value or
    # % a function of 'x', 'y', and 'r' to specify a varying ratio
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['absorbance_rate'] = 0.0
    
    return optical_system


def create_camera_optical_system(piv_simulation_parameters):
    # This function creates the optical system used in running the ray tracing
    # simulation.

    # This extracts the lens focal length from the design structure (in microns)
    focal_length=piv_simulation_parameters['lens_design']['focal_length']
    
    # This extracts the lens f/# from the design structure
    aperture_f_number=piv_simulation_parameters['lens_design']['aperture_f_number']

    # This creates the default single lens optical system which can be modified
    # for running the current simulation
    optical_system = create_single_lens_optical_system()

    # This is the required pitch (aperture) of the lens to match the focal 
    # length and f number
    lens_pitch = focal_length/aperture_f_number

    # This arbitrarily defines the radius of curvature of the surfaces of the
    # lens to be equal to 100 mm (any value can be used here)
    lens_radius_of_curvature = 100.0e3
    
    # This calculates the thickness of the lens based upon the defined
    # curvature and pitch
    lens_thickness = (lens_radius_of_curvature - np.sqrt(np.power(lens_radius_of_curvature,2) - np.power(lens_pitch,2)))/2.0

    # This calculates the refractive index that would be required to produce a
    # lens with the specified parameters
    refractive_index_1 = (2.0*lens_thickness*focal_length - 2.0*focal_length*lens_radius_of_curvature - np.power(lens_radius_of_curvature,2) - lens_radius_of_curvature*np.sqrt(-4.0*lens_thickness*focal_length + np.power(2.0*focal_length + lens_radius_of_curvature,2)))/(2.0*focal_length*(lens_thickness - 2.0*lens_radius_of_curvature))
    refractive_index_2 = (2.0*lens_thickness*focal_length - 2.0*focal_length*lens_radius_of_curvature - np.power(lens_radius_of_curvature,2) + lens_radius_of_curvature*np.sqrt(-4.0*lens_thickness*focal_length + np.power(2.0*focal_length + lens_radius_of_curvature,2)))/(2.0*focal_length*(lens_thickness - 2.0*lens_radius_of_curvature))

    # This is a temporary vector of both refractive index possibilities that 
    # checks whether the values are real and greater than 1
    refractive_index_temp = np.array([np.isreal(refractive_index_1),np.isreal(refractive_index_2)])*np.array([refractive_index_1>=1,refractive_index_2>=1])*np.array([refractive_index_1,refractive_index_2])

    # This sets any values equal to zero equal to infinity since they are were 
    # identified as not valid refractive indices in the last step (and will be 
    # thrown out in the next step)
    refractive_index_temp[refractive_index_temp==0] = np.infty

    # This is the smallest of the remaining possible refractive index values
    refractive_index = np.min(refractive_index_temp)

    # This is a string specifying the shape of the front surface of the current
    # optical element which must be some function of the independent variables
    # 'x', 'y', and 'r' where r = sqrt(x^2+y^2), constant values added or
    # subtracted from this function are ignored
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['front_surface_shape'] = '-np.sqrt((' + '%f' % lens_radius_of_curvature + ')**2-(x**2+y**2))'

    # This is a string specifying the shape of the front surface of the current
    # optical element which must be some function of the independent variables
    # 'x', 'y', and 'r' where r = sqrt(x^2+y^2), constant values added or
    # subtracted from this function are ignored
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['back_surface_shape'] = '+np.sqrt((' + '%f' % lens_radius_of_curvature + ')**2-(x**2+y**2))'

    # This is a double value specifying the pitch (diameter) of the current
    # optical element
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['pitch'] = lens_pitch

    # This is a double value specifying the radius of curvature of the front
    # surface of the current optical element (if the element is not spherical,
    # then this value is Null)
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['front_surface_radius'] = +lens_radius_of_curvature

    # This is a double value specifying the radius of curvature of the back
    # surface of the current optical element (if the element is not spherical,
    # then this value is Null)
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['back_surface_radius'] = -lens_radius_of_curvature

    # This is a double value specifying the distance between the optical axis
    # vertices on the current optical element
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['vertex_distance'] = lens_thickness

    # This adds the refractive index of the current optical element - if this
    # is a double value, the refractive index is taken to be constant, if this
    # is a string, then the refractive index is taken to be a function of the
    # independent variable wavelength 'lambda'
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['refractive_index'] = refractive_index

    # This is the equivalent thin lens focal length for the current element (if
    # this value is defined - if not, then the value is Null)
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['thin_lens_focal_length'] = focal_length
    
    return optical_system

def calculate_rotation_matrix(theta_x,theta_y,theta_z):
# This function calculates the rotation matrix of the angles 'theta_x',
# 'theta_y', and 'theta_z' which are measured in radians and returns the
# rotation matrix as the output argument 'rotation_matrix'.

    # This creates the rotation matrix about the x axis
    rotation_x = np.matrix([[1.0,0.0,0.0],[0.0,np.cos(theta_x),np.sin(theta_x)],[0.0,-np.sin(theta_x),np.cos(theta_x)]])

    # This creates the rotation matrix about the y axis
    rotation_y = np.matrix([[np.cos(theta_y),0.0,-np.sin(theta_y)],[0.0, 1.0, 0.0],[np.sin(theta_y),0.0,np.cos(theta_y)]])

    # This creates the rotation matrix about the z axis
    rotation_z = np.matrix([[np.cos(theta_z),np.sin(theta_z),0.0],[-np.sin(theta_z), np.cos(theta_z), 0.0],[0.0,0.0,1.0]])

    # This is the full rotation matrix
    rotation_matrix = rotation_x * rotation_y * rotation_z

    return rotation_matrix

def rotate_coordinates(X,Y,Z,Alpha,Beta,Gamma,XC,YC,ZC):
# This function takes the coordinates specified by the arrays X, Y, and Z
# and rotates them by angles Alpha, Beta, and Gamma about the x, y, and z
# axes respectively.  The scalars XC, YC, and ZC correspond to the center
# of rotation, ie the coordinates will be rotated by an angle Alpha about
# the point (YC,ZC), an angle Beta about the point (XC,ZC), and by an angle
# Gamma about the point (XC,YC).
#
# This function was extracted from a sub-function so that multiple
# functions could call this same copy.
#
# Authors: Rod La Foy
# First Created On: 11 November 2013
# Last Modified On: 11 November 2013

    # This calculates the rotation matrix for the current angles
    R = calculate_rotation_matrix(Alpha,Beta,Gamma)

    # This translates the coordinates (X,Y,Z) so that the point (XC,YC,ZC) now
    # corresponds to the system origin
    XR = X - XC
    YR = Y - YC
    ZR = Z - ZC

    # This reshapes the arrays to linear (row) arrays for applying the rotation
    # matrix
#     XR = np.reshape(XR,(1,XR.size))
#     YR = np.reshape(YR,(1,YR.size))
#     ZR = np.reshape(ZR,(1,ZR.size))
    XYZR = np.array([XR,YR,ZR])
    
    # This applies the rotation to the coordinates
    WR = np.dot(R,XYZR)

    # This extracts the XR, YR, and ZR coordinates from the rotated matrix
    XR = WR[0][:]
    YR = WR[1][:]
    ZR = WR[2][:]

    # This reshapes the vectors into the shape of the original arrays
#     XR = np.reshape(XR,(X.size,1))
#     YR = np.reshape(YR,(Y.size,1))
#     ZR = np.reshape(ZR,(Z.size,1))

    # This translates the coordinates so that the origin is shifted back to the
    # original position
    XR = XR + XC
    YR = YR + YC
    ZR = ZR + ZC

    return [XR,YR,ZR]

def load_lightfield_data(piv_simulation_parameters,optical_system,mie_scattering_data,frame_index,lightfield_source):
# This function creates the lightfield data for performing the ray tracing
# operation.

    # lightfield_source = {}
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Extracts parameters from 'piv_simulation_parameters'                    %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # This extracts the object distance of the lens (ie the distance between
    # the lens front principal plane and the center of the focal plane) to the
    # design structure (in microns)
    object_distance = piv_simulation_parameters['lens_design']['object_distance']

    # This extracts the lens focal length from the design structure (in microns)
    focal_length = piv_simulation_parameters['lens_design']['focal_length']

    # This extracts the x angle of the camera to the particle volume
    x_camera_angle = piv_simulation_parameters['camera_design']['x_camera_angle']

    # This extracts the y angle of the camera to the particle volume
    y_camera_angle = piv_simulation_parameters['camera_design']['y_camera_angle']

    # This extracts the directory containing the particle locations from
    # the parameters structure
    data_directory = piv_simulation_parameters['particle_field']['data_directory']

    # This extracts the prefix of the particle data filenames from the
    # parameters structure
    data_filename_prefix = piv_simulation_parameters['particle_field']['data_filename_prefix']

    # This extracts the number of particles to simulate out of the list of
    # possible particles (if this number is larger than the number of
    # saved particles, an error will be returned)
    particle_number = piv_simulation_parameters['particle_field']['particle_number']

    # This extracts the Full Width Half Maximum of the laser sheet
    # Gaussian function (in microns) which will produce an illuminated
    # sheet on the XY plane
    gaussian_beam_fwhm = piv_simulation_parameters['particle_field']['gaussian_beam_fwhm']

    # This extracts a Boolean value stating whether to perform Mie scattering
    # simulation
    perform_mie_scattering = piv_simulation_parameters['particle_field']['perform_mie_scattering']

    # TEMPORARY
    if perform_mie_scattering:
        mie_scattering_data = {}

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Extracts parameters from 'optical_system'                               %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # This extracts the lens refractive index
    refractive_index = optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['refractive_index']

    # This extracts the radius of curvature for the front surface of the lens
    front_surface_radius = optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['front_surface_radius']

    # This extracts the radius of curvature for the back surface of the lens
    back_surface_radius = optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['back_surface_radius']

    # This extracts the front surface to back surface vertex distance of the
    # lens
    optical_system_length = optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['vertex_distance']

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Extracts parameters from 'mie_scattering_data'                          %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # This extracts the particle diameter distribution if specified by the
    # parameters structure, otherwise a Null value is used

    #COMMENTING OUT THE MIE SCATTERING FOR NOW
    #     if perform_mie_scattering:
    #         # This extracts the particle diameter index (into 'scattering_irradiance')
    #         # from the parameters structure
    #         print 'Hello'
    #         particle_diameter_index_distribution = mie_scattering_data['particle_diameter_index_distribution']
    #         # This sets the irradiance constant for using Mie scattering
    #         irradiance_constant = 500.0/1.0e4
    #     else:
    #         # This sets the particle diameter indices to a Null value
    #         particle_diameter_index_distribution= None
    #         # This sets the irradiance constant for not using Mie scattering
    #         irradiance_constant = 500.0

    # This sets the particle diameter indices to a Null value
    particle_diameter_index_distribution = None
    # This sets the irradiance constant for not using Mie scattering
    irradiance_constant = 500.0

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Calculates optical system properties                                    %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # This calculates the image distance
    image_distance = np.power(1.0/focal_length - 1.0/object_distance,-1.0)

    # This calculates the offsets of the principal planes
    h1_principal_plane = -(focal_length*(refractive_index-1.0)*optical_system_length)/(back_surface_radius*refractive_index)
    h2_principal_plane=-(focal_length*(refractive_index-1.0)*optical_system_length)/(front_surface_radius*refractive_index)

    # This calculates the positions of the lens vertex planes
    v2_vertex_plane = image_distance + h2_principal_plane
    v1_vertex_plane = v2_vertex_plane + optical_system_length

    # This is the object distance of the lens
    z_object = v1_vertex_plane - h1_principal_plane + object_distance

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Loads the current particle field data                                   %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # This is the list of particle data files that can be loaded
    particle_data_list = glob.glob(data_directory + data_filename_prefix + '*.mat')

    # This is the filename of the first frame specified by the
    # 'frame_vector' vector
    particle_data_filename_read = particle_data_list[frame_index-1]

    # This loads the particle data file to memory
    mat_contents = sio.loadmat(particle_data_filename_read, squeeze_me = True)
    particle_data = pickle.load(open(particle_data_filename_read[:-3] + 'p','rb'))

    # reads variable from mat_contents
    # X = np.array(mat_contents['X'])
    # Y = np.array(mat_contents['Y'])
    # Z = np.array(mat_contents['Z'])
    X = np.squeeze(particle_data['X'])
    Y = np.squeeze(particle_data['Y'])
    Z = np.squeeze(particle_data['Z'])

    # This extracts only the specified number of particles from the
    # particle position vectors
    # NOTE: During array slicing in python, the last index is particle_number -1
    X = X[0:int(particle_number)]
    Y = Y[0:int(particle_number)]
    Z = Z[0:int(particle_number)]

    # This calculates the standard deviation of the Gaussian beam function
    # based upon the Full Width Half Maximum
    gaussian_sigma = gaussian_beam_fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))

    # This calculates the intensity of the particles based upon the
    # Gaussian beam distribution (the scale factor coefficient is to ensure
    # that the particle intensity gives roughly "typical" results for default
    # parameters)
    R = irradiance_constant*(1.0/(gaussian_sigma*np.sqrt(2.0*np.pi))) * np.exp(-1.0*(np.power(Z,2.0)/(2.0*np.power(gaussian_sigma,2.0))))
    # This rotates the image coordinates by the specified angles
    [X,Y,Z] = rotate_coordinates(X,Y,Z,x_camera_angle,y_camera_angle,0.0,0.0,0.0,0.0)

    # This translates the Z coordinates of the paricles to the focal plane
    Z = Z + z_object

    # This adds in the particles to the lightfield source data
    lightfield_source['x'] = X
    lightfield_source['y'] = Y
    lightfield_source['z'] = Z
    lightfield_source['radiance'] = np.array(R)
    lightfield_source['diameter_index'] = particle_diameter_index_distribution

    return lightfield_source

def calculate_sunflower_coordinates(grid_point_diameter,lightray_number_per_grid_point):
#% This function creates two coordinate vectors x and y that consist of a
#% series of points that fill a circle centered at (0,0) with a diameter
#% equal to the grid_point_diameter.  The points will lie on concentric
#% circles such that the distance between adjacent circles is constant and
#% equal to the approximate nearest neighbor distance of all points.  The
#% number of output points is approximately equal to
#% lightray_number_per_grid_point.
#
    #% This is the area of the grid point in microns
    grid_point_area = np.pi*(grid_point_diameter/2.0)**2.0
    #% This is the average distance to the nearest lightray point on the grid
    #% point based upon the number of rays and the size of the point
    lightray_point_spacing = np.sqrt(grid_point_area/lightray_number_per_grid_point)
    
    #% This is a vector of radii to place the lightray points at
    radius_lightray_vector = np.linspace(lightray_point_spacing,(grid_point_diameter/2.0),np.round((grid_point_diameter/2.0)/lightray_point_spacing))
    
    #% This is the density of the lightray points per unit distance
    rho = 1.0/lightray_point_spacing
    
    #% This initiliases the coordiates of the lightrays
    x_lightray_coordinates=[]
    y_lightray_coordinates=[]
    
    #% This iterates through the different radii generateing the points
    for n in range(0,len(radius_lightray_vector)):
        
        # This is the radius of the current circle
        radius_current = radius_lightray_vector[n]
        
        # This is the number of points to generate within the current circle
        circle_lightray_point_number = np.round(rho*(2.0 * np.pi * radius_current))
        
        # This is a vector of angles for the points
        theta_current = (2.0*np.pi/circle_lightray_point_number)*(np.arange(1.0,circle_lightray_point_number)-1.0) + 2.0*np.pi*np.random.rand(1,1)
        
        # These are the coordinates of the current radius's points
        x_temp = radius_current * np.cos(theta_current)
        y_temp = radius_current * np.sin(theta_current)
        
        # This adds the current coordinates to the full coordinate vector
        if(n==0):
            x_lightray_coordinates = np.array(x_temp)
            y_lightray_coordinates = np.array(y_temp)
        else:
            x_lightray_coordinates = np.append(x_lightray_coordinates,x_temp)
            y_lightray_coordinates = np.append(y_lightray_coordinates,y_temp)

    
    #% This adds the origin to the vector of coordinates
    x_lightray_coordinates = np.transpose(np.append(x_lightray_coordinates,0.0))
    y_lightray_coordinates = np.transpose(np.append(y_lightray_coordinates,0.0))
    
    return [x_lightray_coordinates,y_lightray_coordinates]





def generate_calibration_lightfield_data(piv_simulation_parameters,optical_system,plane_index):
# % This function generates a structure containing the location of the each of the
# % lightfield source points as well as their irradiance and color.  Only a single value is
# % returned for each source point and no angular information is stored.  The angular data 
# % of the lightrays will be generated later - this is done to avoid out of memory errors.

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Extracts parameters from 'piv_simulation_parameters'                    %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # % This extracts the object distance of the lens (ie the distance between
    # % the lens front principal plane and the center of the focal plane) to the
    # % design structure (in microns)
    object_distance = piv_simulation_parameters['lens_design']['object_distance']
    # % This extracts the lens focal length from the design structure (in microns)
    focal_length = piv_simulation_parameters['lens_design']['focal_length']

    # % This extracts the calibration plane spacing from the structure
    calibration_plane_spacing = piv_simulation_parameters['calibration_grid']['calibration_plane_spacing']
    # % This extracts the calibration plane number from the structure
    calibration_plane_number = piv_simulation_parameters['calibration_grid']['calibration_plane_number']
    # % This extracts the grid point diameter from the structure
    grid_point_diameter = piv_simulation_parameters['calibration_grid']['grid_point_diameter']
    # % This extracts the grid point spacing from the structure
    x_grid_point_spacing = piv_simulation_parameters['calibration_grid']['x_grid_point_spacing']
    y_grid_point_spacing = piv_simulation_parameters['calibration_grid']['y_grid_point_spacing']
    # % This extracts the grid point number from the calibration structure
    x_grid_point_number = piv_simulation_parameters['calibration_grid']['x_grid_point_number']
    y_grid_point_number = piv_simulation_parameters['calibration_grid']['y_grid_point_number']
    # % This extracts the number of light source partciles per grid point from 
    # % the calibration structure
    particle_number_per_grid_point = piv_simulation_parameters['calibration_grid']['particle_number_per_grid_point'];

    # % This extracts the x angle of the camera to the particle volume
    x_camera_angle = piv_simulation_parameters['camera_design']['x_camera_angle']
    # % This extracts the y angle of the camera to the particle volume
    y_camera_angle = piv_simulation_parameters['camera_design']['y_camera_angle']

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Extracts parameters from 'optical_system'                               %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # % This extracts the lens refractive index
    refractive_index = optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['refractive_index']
    # % This extracts the radius of curvature for the front surface of the lens
    front_surface_radius=optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['front_surface_radius']
    # % This extracts the radius of curvature for the back surface of the lens
    back_surface_radius=optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['back_surface_radius']
    # % This extracts the front surface to back surface vertex distance of the
    # % lens
    optical_system_length=optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry']['vertex_distance']

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Calculates optical system properties                                    %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # % This calculates the image distance
    image_distance = (1.0/focal_length - 1.0/object_distance)**-1

    # % This calculates the offsets of the principal planes
    h1_principal_plane = -(focal_length * (refractive_index - 1.0) * optical_system_length)/(back_surface_radius * refractive_index)
    h2_principal_plane = -(focal_length * (refractive_index - 1.0) * optical_system_length)/(front_surface_radius * refractive_index)

    # % This calculates the positions of the lens vertex planes
    v2_vertex_plane = image_distance + h2_principal_plane
    v1_vertex_plane = v2_vertex_plane + optical_system_length

    # % This is the object distance of the lens
    z_object = v1_vertex_plane - h1_principal_plane + object_distance

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Generates the calibration grid data                                     %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # % This is a vector of Z positions corresponding to the locations Z
    # % locations of the calibration planes
    grid_plane_z_world_coordinate = calibration_plane_spacing*np.linspace(-(calibration_plane_number-1.0)/2.0,(calibration_plane_number-1.0)/2.0,calibration_plane_number,endpoint=True)
    # % This is the Z world coordinate of the current grid to simulate
    current_z_world_coordinate = grid_plane_z_world_coordinate[plane_index]

    # % This is a vector of the x world coordinates of the grid points
    x_grid_point_coordinate_vector = x_grid_point_spacing*np.linspace(-(x_grid_point_number-1.0)/2.0,(x_grid_point_number - 1.0)/2.0,x_grid_point_number,endpoint=True)
    # % This is a vector of the y world coordinates of the grid points
    y_grid_point_coordinate_vector = y_grid_point_spacing*np.linspace(-(y_grid_point_number-1.0)/2.0,(y_grid_point_number - 1.0)/2.0,y_grid_point_number,endpoint=True)

    # % This generates a series of points that fill the circle of the grid point
    # % uniformly
    [x_lightray_coordinates,y_lightray_coordinates] = calculate_sunflower_coordinates(grid_point_diameter,particle_number_per_grid_point);

    # % These are the full coordinate vectors of all the lightrays
    x = np.zeros(int(particle_number_per_grid_point*x_grid_point_number*y_grid_point_number))
    y = np.zeros(int(particle_number_per_grid_point*x_grid_point_number*y_grid_point_number))

    # % This is a counting index
    count = 0

    # % This iterates through the grid points generating the light rays for each
    # % point
    for x_grid_index in range(0,int(x_grid_point_number)):
        for y_grid_index in range(0,int(y_grid_point_number)):
            
            # % This is the coordinate of the current grid point
            x_grid_point_coordinate = x_grid_point_coordinate_vector[x_grid_index]
            y_grid_point_coordinate = y_grid_point_coordinate_vector[y_grid_index]

            #% These are the coordinates of the current lightrays
            x_grid_point_lightray_coordinates = x_lightray_coordinates + x_grid_point_coordinate
            y_grid_point_lightray_coordinates = y_lightray_coordinates + y_grid_point_coordinate

            # This is the vector of indices to add in the new coordinates
            index_vector = (count+1) + np.linspace(1,len(x_grid_point_lightray_coordinates),endpoint=True).astype('int')

            # This increments the counting variable
            count = count + len(x_grid_point_lightray_coordinates)
            
            # This adds the current lightrays tot the total lightray vectors
            if(x_grid_index==0 and y_grid_index==0):
                x = x_grid_point_lightray_coordinates
                y = y_grid_point_lightray_coordinates
            else:
                x = np.append(x,x_grid_point_lightray_coordinates)
                y = np.append(y,y_grid_point_lightray_coordinates)
            # x[index_vector] = x_grid_point_lightray_coordinates
            # y[index_vector] = y_grid_point_lightray_coordinates
            #
     
    # % This eliminates any coordinates that were initilized but not changed
    if (count+1)<len(x):
        # This deletes the ending values from the coordinate vectors
        # HOW TO DELETE ELEMENTS IN PYTHON?
        x = x[:count+1]
        y = y[:count+1]
#        x[count+1:-1]=[];
#        y[count+1:-1]=[];

    # % This generates a series of points that fill the circle of one quarter of 
    # % the grid point uniformly
    [x_lightray_coordinates, y_lightray_coordinates] = calculate_sunflower_coordinates(grid_point_diameter/4,particle_number_per_grid_point/16);

    # % These are the coordinates of the first origin marker
    x_grid_point_coordinate = -x_grid_point_spacing/2.0
    y_grid_point_coordinate = 0

    # % These are the coordinates of the current lightrays
    x_grid_point_lightray_coordinates = x_lightray_coordinates + x_grid_point_coordinate
    y_grid_point_lightray_coordinates = y_lightray_coordinates + y_grid_point_coordinate

    # % This adds the current lightrays to the total lightray vectors
    x = np.append(x,x_grid_point_lightray_coordinates)
    y = np.append(y,y_grid_point_lightray_coordinates)

    # % These are the coordinates of the second origin marker
    x_grid_point_coordinate = 0
    y_grid_point_coordinate = y_grid_point_spacing/2.0

    # % These are the coordinates of the current lightrays
    x_grid_point_lightray_coordinates = x_lightray_coordinates + x_grid_point_coordinate
    y_grid_point_lightray_coordinates = y_lightray_coordinates + y_grid_point_coordinate

    # % This adds the current lightrays tot he total lightray vectors
    x = np.append(x,x_grid_point_lightray_coordinates)
    y = np.append(y,y_grid_point_lightray_coordinates)

    # % This initializes the Z coordinate vector
    z = current_z_world_coordinate * np.ones(x.shape)
    # % This initializes the radiance vector
    radiance = np.ones(x.shape);

    # % This rotates the image coordinates by the specified angles
    [x,y,z] = rotate_coordinates(x,y,z,x_camera_angle,y_camera_angle,0.0,0.0,0.0,0.0)

    # % % This rotates the particles by the specified angles
    # % [x_source,y_source,z_source]=rotate_coordinates(x_source,y_source,z_source,theta_x,theta_y,theta_z,0,0,0);
    # % This translates the Z coordinates of the paricles to the focal plane
    z = z + z_object

    # % This adds in the particles to the lightfield source data
    lightfield_source = {'x': x, 'y': y, 'z': z, 'radiance': radiance, 'diameter_index': np.ones(x.shape)}

    return lightfield_source


def run_piv_simulation_02(piv_simulation_parameters):
# This function runs a PIV simulation using the thick lens, non-paraxial
# camera simulation.

    # % % This creates the simulation parameters for running the simulation
    # % piv_simulation_parameters=create_piv_simulation_parameters_02;
    optical_system = {}
    # This creates the optical system parameters used to simulate the defined lens system
    optical_system = create_camera_optical_system(piv_simulation_parameters)

    # convert Nones to nans so you can save it as a MATLAB file. then reconvert them to null in MATLAB
    # convert Nones to nans so you can save it as a MATLAB file. then reconvert them to null in MATLAB
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_number'] = np.NAN
    optical_system['design']['optical_element'][0]['optical_element'][0]['elements_coplanar'] = np.NAN
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['abbe_number'] = np.NAN
    optical_system['design']['optical_element'][0]['element_geometry'] = np.NAN
    optical_system['design']['optical_element'][0]['element_properties'] = np.NAN

    # save optical system data to file
    sio.savemat('optical_system.mat', {'optical_system': optical_system},long_field_names=True)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Particle Field Simulation                                               %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # # % This extracts the Boolean value stating whether to generate the particle
    # # % field images from the structure
    generate_particle_field_images = piv_simulation_parameters['particle_field']['generate_particle_field_images']

    # % This extracts the vector giving the frames of particle positions to load to
    # % the parameters structure (this indexes into the list generated by the
    # % command 'dir([data_directory,data_filename_prefix,'*.mat'])')
    frame_vector = piv_simulation_parameters['particle_field']['frame_vector']

    # % This adds the directory to save the particle images to parameters
    # % structure
    particle_image_directory = piv_simulation_parameters['output_data']['particle_image_directory']

    # % This extracts the number of lightrays to simulate per particle (this is roughly
    # % equivalent to the power of the laser)
    lightray_number_per_particle = piv_simulation_parameters['particle_field']['lightray_number_per_particle']

    # % This extracts the number of lightrays to propogate per iteration (this is a
    # % function of the RAM available on the computer)
    lightray_process_number = piv_simulation_parameters['particle_field']['lightray_process_number']

    # % This is the gain of the sensor in decibels to be used in the particle
    # % field simulation
    pixel_gain = piv_simulation_parameters['particle_field']['pixel_gain']

    # % This extracts a Boolean value stating whether to perform Mie scattering
    # % simulation
    perform_mie_scattering = piv_simulation_parameters['particle_field']['perform_mie_scattering']

    # % This generates the particle field images if specified by the parameters
    # % structure
    if generate_particle_field_images:

        # This displays that the particle images are being simulated
        print('\n\n')
        print('Simulating particle images . . . ')
      
        # for the time being, ignore mie scattering
        perform_mie_scattering = False
        scattering_data = [] #None
        scattering_type = 'diffuse'
        
    #     # This calculates the Mie scattering data if specified in the parameters
    #     # data structure, otherwise the Mie scattering data is set to a Null value
    #     if perform_mie_scattering:
    #         # This calculates the Mie scattering parameters data used in simulating the
    #         # scattering intensities of the simulated particles
    #         scattering_data = create_mie_scattering_data(piv_simulation_parameters)
    #         # This sets the scattering type to mie for the particle simulation
    #         scattering_type = 'mie'
    #     else:
    #         # This sets the Mie scattering data to a Null value
    #         scattering_data = None
    #         # This sets the scattering type to diffuse for the particle simulation
    #         scattering_type = 'diffuse'
        
        
         # This iterates through the frame vectors performing the ray tracing operations for each frame
        field_type = 'particle'
        # for frame_index in np.array(frame_vector):
        for frame_index in range(1,3):
            # This creates the lightfield data for performing the raytracing operation
            #print " Generating lightfield source"
            lightfield_source = dict()
                        
            lightfield_source = load_lightfield_data(piv_simulation_parameters,optical_system,scattering_data,frame_index, lightfield_source)
                
            # This adds the number of lightrays per particle to the
            # 'lightfield_source' data
            lightfield_source.update({'lightray_number_per_particle' : lightray_number_per_particle})
        
            # This adds the number of lightrays to simulateously process to the
            # 'lightfield_source' data
            lightfield_source.update({'lightray_process_number' : lightray_process_number})
                    
            #convert None to NAN
            lightfield_source['diameter_index'] = np.NAN

            # convert none to NAN just for MATLAB
            scattering_data = np.NAN

            sio.savemat('particle_' + '%02d' % frame_index + '.mat',{'piv_simulation_parameters' : piv_simulation_parameters,
                                   'optical_system' : optical_system,
                                   'pixel_gain' : pixel_gain,
                                   'scattering_data' : scattering_data,
                                   'scattering_type' : scattering_type,
                                   'lightfield_source' : lightfield_source
                                   },long_field_names=True)

            # % This performs the ray tracing to generate the sensor image
            # I=perform_ray_tracing_03(piv_simulation_parameters,optical_system,pixel_gain,scattering_data,scattering_type,lightfield_source);
            I = perform_ray_tracing_03(piv_simulation_parameters,optical_system,pixel_gain,scattering_data,scattering_type,lightfield_source,field_type)

            # This is the filename to save the image data to
            image_filename_write = particle_image_directory + 'particle_image_frame_' + '%04d' % frame_index + '.tif'

            I.tofile(image_filename_write)
            # This saves the image to memory
            # ski_io.imsave(image_filename_write,I,'tif','Compression','none')

    # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # # % Calibration Grid Simulation                                             %
    # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # # % This extracts the Boolean value stating whether to generate the calibration
    # # % images from the structure
    # generate_calibration_grid_images = piv_simulation_parameters['calibration_grid']['generate_calibration_grid_images']
    # # % This extracts the calibration plane number from the structure
    # calibration_plane_number = int(piv_simulation_parameters['calibration_grid']['calibration_plane_number'])
    # # % This extracts the directory to save the calibration grid images from
    # # % parameters structure
    # calibration_grid_image_directory = piv_simulation_parameters['output_data']['calibration_grid_image_directory']
    # # % This extracts the number of lightrays to simulate per particle (this is roughly
    # # % equivalent to the power of the laser)
    # lightray_number_per_particle = piv_simulation_parameters['calibration_grid']['lightray_number_per_particle']
    # # % This extracts the number of lightrays to propogate per iteration (this is a
    # # % function of the RAM available on the computer)
    # lightray_process_number = piv_simulation_parameters['calibration_grid']['lightray_process_number']
    # # % This is the gain of the sensor in decibels to be used in the calibration
    # # % grid simulation
    # pixel_gain = piv_simulation_parameters['calibration_grid']['pixel_gain']
    #
    # # % This generates the calibration grid images if specified by the parameters
    # # % structure
    # if generate_calibration_grid_images:
    #
    #     # % This displays that the calibration images are being simulated
    #     # fprintf('\n\n');
    #     print 'Simulating calibration images . . . '
    #
    #     # % This sets the scattering type to diffuse for the calibration grid
    #     # % simulation
    #     scattering_type = 'diffuse'
    #     # % This sets the scattering data to a Null value for the calibration grid
    #     scattering_data = None
    #
    #     field_type = 'calibration'
    #     # % This iterates through the calibration grid planes performing the ray
    #     # % tracing operations for each plane
    #     # for plane_index in range(0,calibration_plane_number):
    #     for plane_index in range(0,1):
    #
    #         # % This creates the lightfield data for performing the raytracing operation
    #         lightfield_source = generate_calibration_lightfield_data(piv_simulation_parameters,optical_system,plane_index)
    #         # % This adds the number of lightrays per particle to the
    #         # % 'lightfield_source' data
    #         lightfield_source['lightray_number_per_particle'] = lightray_number_per_particle
    #         # % This adds the number of lightrays to simulateously process to the
    #         # % 'lightfield_source' data
    #         lightfield_source['lightray_process_number'] = lightray_process_number
    #
    #         # save data to mat file
    #         # convert none to NAN just for MATLAB
    #         scattering_data = np.NAN
    #         sio.savemat('calibration_' + '%02d' % (plane_index+1) + '.mat',{'piv_simulation_parameters' : piv_simulation_parameters,
    #                                        'optical_system' : optical_system,
    #                                        'pixel_gain' : pixel_gain,
    #                                        'scattering_data' : scattering_data,
    #                                        'scattering_type' : scattering_type,
    #                                        'lightfield_source' : lightfield_source
    #                                        },long_field_names=True)
    #         # % This performs the ray tracing to generate the sensor image
    #         # I=perform_ray_tracing_03(piv_simulation_parameters,optical_system,pixel_gain,scattering_data,scattering_type,lightfield_source);
    #         I = perform_ray_tracing_03(piv_simulation_parameters,optical_system,pixel_gain,scattering_data,scattering_type,lightfield_source,field_type)
    #
    #         # % This is the filename to save the image data to
    #         image_filename_write=calibration_grid_image_directory + 'calibration_image_plane_' + '%04.0f' % plane_index + '.tif'
            # % This saves the image to memory
            # imwrite(I,image_filename_write,'tif','Compression','none');

