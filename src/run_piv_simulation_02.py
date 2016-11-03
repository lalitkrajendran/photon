import glob
import pickle
import os
import numpy as np
import scipy.io as sio
from perform_ray_tracing_03 import perform_ray_tracing_03
from numpy import linalg as la
import scipy.special as scispec
from bhmie import bhmie
import tifffile as TIFF
import matplotlib.image as mpimg
import skimage.io as ski_io
import matplotlib.pyplot as plt

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
    optical_system['design']['axial_offset_distances'] = np.array([0.0,0.0])

    # This is a 1 x 3 vector giving the rotation angles of the current optical
    # element or system
    optical_system['design']['rotation_angles'] = np.array([0.0,0.0,0.0])

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
    optical_system['design']['optical_element'][0]['element_geometry'] = np.NAN

    # % This creates a substructure to describe the optical properties of the 
    # % current element (the value of the structure will be null if the current 
    # % element is a 'system')
    optical_system['design']['optical_element'][0]['element_properties'] = {}

    # % This specifies the optical properties of the current element - if the
    # % element type is 'system' then the geometry has a Null value
    optical_system['design']['optical_element'][0]['element_properties'] = np.NAN

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
    optical_system['design']['optical_element'][0]['optical_element'][0]['elements_coplanar'] = np.NAN

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
    optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['abbe_number'] = np.NAN

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

    # Added by lalit to account for out of focus effects - did not work
    # lens_thickness = 2.0*(lens_radius_of_curvature - np.sqrt(
    #     np.power(lens_radius_of_curvature, 2) - np.power(lens_pitch/2.0, 2)))

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
    XYZR = np.array([XR,YR,ZR])
    
    # This applies the rotation to the coordinates
    WR = np.dot(R,XYZR)

    # This extracts the XR, YR, and ZR coordinates from the rotated matrix
    XR = WR[0][:]
    YR = WR[1][:]
    ZR = WR[2][:]

    # This translates the coordinates so that the origin is shifted back to the
    # original position
    XR = XR + XC
    YR = YR + YC
    ZR = ZR + ZC

    return [XR,YR,ZR]

def log_normal_pdf(x,mu,sigma):
# This function calculates the log-normal probability density function of
# the argument 'x' with a mean of 'mu' and a standard deviation of 'sigma'.

    # This calculates the log-normal probility density function
    y1 = 1.0/(x * sigma * np.sqrt(2.0 * np.pi))
    y2 = (np.log(x) - mu)
    y3 = np.exp(-np.power(y2,2.0)/(2.0*np.power(sigma,2.0)))
    y = y1 * y3

    return y

def inverse_log_normal_pdf(y,mu,sigma):
# This function calculates the inverse, ie the two values of x such that
#
#  y = pdf(x1) and y = pdf(x2)
#
# for the log-normal probability density function of the argument 'y' with
# a mean of 'mu' and a standard deviation of 'sigma'.

    # This calculates the inverse log-normal probability density function
    x1 = np.exp(mu - sigma**2 - sigma * np.sqrt(sigma**2.0 - 2.0 * mu - 2.0 * np.log(y * sigma * np.sqrt(2.0*np.pi))))
    x2 = np.exp(mu - sigma**2 + sigma * np.sqrt(sigma**2.0 - 2.0 * mu - 2.0 * np.log(y * sigma * np.sqrt(2.0*np.pi))))

    return [x1,x2]

def log_normal_cdf(x,mu,sigma):
# This function calculates the log-normal cumulative density function of
# the argument 'x' with a mean of 'mu' and a standard deviation of 'sigma'.

    # This calculates the log-normal cumulative density function
    y=(1.0 + scispec.erf((np.log(x) - mu)/(sigma * np.sqrt(2.0))))/2.0

    return  y

def calculate_log_normal_pdf_extrema(mu,sigma,t):
# % This function calculates the extrema of the log-normal probability
# % density function such that
# %
# %   t = 1 - ( cdf(x_max) - cdf(x_min) )
# %
# % ie the extrema given the portion of the log-normal pdf such that the
# % percentage of the pdf that lie outside of the extrema equals t.  The
# % input arguments are the log-normal mean 'mu', the log-normal standard
# % deviation 'sigma', and the threshhold 't'.

    # % This initializes a first estimate of the 'x_max' extrema
    x_max = np.exp(mu + sigma)

    # % This iterates through the possible values of the extrema using Newton's
    # % method to numerically calculate the extrema values
    while True:

        # This calculates the value of the log-normal pdf at this location
        y = log_normal_pdf(x_max,mu,sigma)

        # This inverts the log-normal pdf to calculate the other x value that will
        # produce the current pdf value
        [x_min,x_max] = inverse_log_normal_pdf(y,mu,sigma)

        # This calculates the difference between the current 'x_max' estimate and
        # the next iteration 'x_max' estimate
        dx1 = 1.0 - (log_normal_cdf(x_max,mu,sigma) - log_normal_cdf(x_min,mu,sigma)) - t
        dx2 = -np.exp(2.0*mu - 2.0*np.power(sigma,2.0))/np.power(x_max,2.0)
        dx3 = log_normal_pdf(x_min,mu,sigma) * dx2 - log_normal_pdf(x_max,mu,sigma)
        dx = dx1/dx3

        # This breaks the iteration loop if the change in the 'x_max' value is
        # smaller than the machine precision for 'mu'
        if np.abs(dx) < (np.finfo(type(mu)).eps*1.0e2):
            # This breaks the iteration of Newton's method
            break

        # This calculates the next estimate of the upper limit using Newton's
        # method
        x_max = x_max - dx

    return [x_min,x_max]

def calculate_particle_diameter_distribution(piv_simulation_parameters):
# This function calculates the discrete values of the possible particle
# diameters to be simulated during the experiment.

    # This extracts the mean diameter of the particles being used in the
    # simulated experiment (in microns)
    particle_diameter_mean = piv_simulation_parameters['particle_field']['particle_diameter_mean']

    # % This extracts the standard deviation of the particle diameter used in
    # % the simulated experiment (in microns)
    particle_diameter_std = piv_simulation_parameters['particle_field']['particle_diameter_std']

    # % This extracts the number of different particle sizes to model since the
    # % particle diameters are taken in discrete intervals for computational
    # % reasons
    particle_diameter_number = piv_simulation_parameters['particle_field']['particle_diameter_number']

    # % This extracts the cutoff threshhold of the log-normal cumulative density
    # % function beyond which extrema particle diameters are not calculated (ie
    # % if this is set to 0.01 then 1% of the possible particle diameters both
    # % much smaller  and much larger than the mean diameter that would be found
    # % on a continuous particle diameter range will not be included in the
    # % simulation)
    particle_diameter_cdf_threshhold = piv_simulation_parameters['particle_field']['particle_diameter_cdf_threshhold']

    # % This calculates the location parameter 'mu' of the particle diameters,
    # % this is different than the mean since the distribution is log-normal
    particle_diameter_mu = np.log(particle_diameter_mean) - 0.5 * np.log(1.0 + (particle_diameter_std/particle_diameter_mean)**2)

    # % This calculates the scale parameter 'sigma' of the particle diameters,
    # % this is different than the standard deviation since the distribution is
    # % log-normal
    particle_diameter_sigma = np.sqrt(np.log(1.0 + (particle_diameter_std/particle_diameter_mean)**2.0))

    # % This calculates the minimum and maximum particle diameters such that the
    # % ratio of particle diameters generated between the two extrema compared to
    # % all possible particle diameters equals the value given by the expression
    # % 1 - particle_diameter_cdf_threshhold
    [minimum_particle_diameter,maximum_particle_diameter]=calculate_log_normal_pdf_extrema(particle_diameter_mu,particle_diameter_sigma,particle_diameter_cdf_threshhold);

    # % This calculates the spacing between the adjacent particle diameters
    particle_diameter_spacing = (maximum_particle_diameter-minimum_particle_diameter)/particle_diameter_number

    # % This is a vector of the particle diameters
    particle_diameter_vector = minimum_particle_diameter + particle_diameter_spacing * (np.arange(0,particle_diameter_number) + 0.5)

    # % This calculates the probability density function at each of the particle
    # % diameters (to calculate the ratios of different particle diameters to
    # % one-another)
    particle_diameter_pdf = log_normal_pdf(particle_diameter_vector, particle_diameter_mu, particle_diameter_sigma)
    # % This renormalizes the particle diameter pdf so that it sums to one
    particle_diameter_pdf = particle_diameter_pdf/np.sum(particle_diameter_pdf)

    return [particle_diameter_vector,particle_diameter_pdf]

def calculate_particle_diameter_indices(piv_simulation_parameters,particle_diameter_pdf,particle_diameter_vector):
# This function calculates the distribution of the particle diameter
# indices based upon the particle diameter probability density function
# and the total particle number.

    # This extracts the number of particles to simulate out of the list of
    # possible particles
    particle_number = piv_simulation_parameters['particle_field']['particle_number']

    # This calculates the cumulative sum of the particle diameter PDF function
    # to determine the particle diameter distribution
    particle_diameter_cdf = np.cumsum(particle_diameter_pdf)

    # This adds a zero to the beginning of the cumulative distribution function
    # for indexing purposed in the following for loop
    particle_diameter_cdf = np.insert(particle_diameter_cdf,0,0.0)

    # This generates a list of random numbers to randomly determine the
    # particle diameters
    np.random.seed(232)
    random_vector = np.random.rand(particle_number)

    # This initializes a vector of the particle diameters
    particle_diameter_index_distribution = np.zeros(particle_number)

    # This iterates through the possible particle diameters assigning the
    # diameters based upon the values of the 'random_vector' variable
    for particle_diameter_index in range(0,len(particle_diameter_vector)-1):

        # These are the indices of the particles to set to the current diameter
        term1 = (particle_diameter_cdf[particle_diameter_index]<=random_vector)
        term2 = (particle_diameter_cdf[particle_diameter_index + 1] > random_vector)
        diameter_indices = np.where(np.logical_and(term1,term2))
        #diameter_indices = np.where(particle_diameter_cdf[particle_diameter_index]<=random_vector and particle_diameter_cdf[particle_diameter_index>random_vector])
        #diameter_indices = (term1&term2).astype('int')

        ## THIS COULD BE A PROBLEM?############################################
        # This sets the current indices of the particle diameter distribution
        # equal to the current diameter
        particle_diameter_index_distribution[diameter_indices] = particle_diameter_index

    return particle_diameter_index_distribution

def calculate_mie_scattering_intensity(piv_simulation_parameters,particle_diameter_vector):
# % This function calculates the intensity of the Mie scattering produced by
# % the different particle diameters given by 'particle_diameter_vector' for
# % the illumination source defined in 'piv_simulation_parameters'.

    # % This extracts the refractive index of the medium in which the particles
    # % are seeded (typically either water or air)
    medium_refractive_index = piv_simulation_parameters['particle_field']['medium_refractive_index']
    # % This extracts the refractive index of the seeding particles used in the
    # % simulation
    particle_refractive_index = piv_simulation_parameters['particle_field']['particle_refractive_index']

    # % This extracts the number of angles to calculate the Mie scattering
    # % intensity over (which is later interpolated to the precise angles for
    # % each paricle)
    mie_scattering_angle_number = piv_simulation_parameters['particle_field']['mie_scattering_angle_number']

    # % This is the wavelength of the simulated laser used for illumination of
    # % the particles
    beam_wavelength = piv_simulation_parameters['particle_field']['beam_wavelength']

    # % This initializes the scattering irradiance array for the particle the
    # % range of particle diameters
    nrows = int(2*mie_scattering_angle_number - 1)
    ncols = len(particle_diameter_vector)
    scattering_irradiance = np.zeros((nrows,ncols))

    # % This iterates through the different particle diameters calculating the
    # % Mie scattering intensities
    for particle_diameter_index in range (0,len(particle_diameter_vector)):
        # This is the radius of the current particle for which the Mie scattering
        # will be calculated
        current_particle_radius = particle_diameter_vector[particle_diameter_index]

        # This calculates the Mie scattering intensities for the particles

        # calculate size parameter
        x = 2*np.pi*current_particle_radius*medium_refractive_index/beam_wavelength
        # set refractive index of medium (air)
        ref_med = medium_refractive_index 
        # set refractive index of particle (olive oil)
        ref_part = particle_refractive_index
        # calculate relative refractive index
        refrel = ref_part/ref_med
        # specify number of angles between 0 and 90 where irradiance data is reqd
        nang = mie_scattering_angle_number
        # call bhmie function to evaluate mie scattering intensities
        [s1,s2,qext,qsca,qback,gsca,theta] = bhmie(x,refrel,nang)
        dang = 0.5*np.pi/(nang-1)
        # create range of angles (in radians)
        scattering_angle = (np.arange(1,s1.size+1)-1.0)*dang
        s11 = 0.5*(abs(s1)**2 + abs(s2)**2)

        scattering_irradiance[:,particle_diameter_index] = s11
        
    return [scattering_angle,scattering_irradiance]

def create_mie_scattering_data(piv_simulation_parameters):
# This function creates various parameters used in simulation the Mie
# scattering of the particles.  Specifically this calculate the size
# distribution of the particles, the Mie scattering intensities as a
# function of the scattering angles, and several parameters necessary to
# calculate the scattering angles of the particles with respect to the
# laser beam and the simulated cameras.

    # This is the x angle of the camera to the particle volume
    x_camera_angle = piv_simulation_parameters['camera_design']['x_camera_angle']
    # This is the y angle of the camera to the particle volume
    y_camera_angle = piv_simulation_parameters['camera_design']['y_camera_angle']

    # This is a direction vector (ie the magnitude doesn't matter) that points
    # in the direction of the laser beam propogation - this vector (at least
    # for now) is defined by a 1 x 3 array and lies in the XY plane (ie the
    # last component must be zero)
    beam_propogation_vector = piv_simulation_parameters['particle_field']['beam_propogation_vector']

    # This calculates the particle diameters that will be simulated and the
    # relative frequency of each of the diameters
    [particle_diameter_vector,particle_diameter_pdf] = calculate_particle_diameter_distribution(piv_simulation_parameters)

    # % This calculates the distribution of the particle diameter indices
    # % (indexing into 'particle_diameter_vector') based upon the particle
    # % diameter probability density function and the total particle number
    particle_diameter_index_distribution = calculate_particle_diameter_indices(piv_simulation_parameters,particle_diameter_pdf,particle_diameter_vector)

    # % This calculates the Mie scattering intensity data for the set of particle
    # % diameters and illumination source
    [scattering_angle,scattering_irradiance] = calculate_mie_scattering_intensity(piv_simulation_parameters,particle_diameter_vector)

    # This calculates the rotation matrix that transforms between the the world
    # coordinate system and the camera coordinate system
    rotation_matrix = calculate_rotation_matrix(x_camera_angle,y_camera_angle,0.0)

    # This computes the inverse rotation matrix (ie the transpose of the
    # rotation matrix)
    inverse_rotation_matrix = rotation_matrix.transpose()

    # This normalizes the laser beam propogation vector
    beam_propogation_vector = beam_propogation_vector/la.norm(beam_propogation_vector)

    # This creates a structure to store the Mie scattering data parameters
    # within
    mie_scattering_data = {}

    # This saves the paticle diameter vector into the parameters structre
    mie_scattering_data['particle_diameter_vector'] = particle_diameter_vector

    # This saves the particle diameter probability density function into the
    # parameters structure
    mie_scattering_data['particle_diameter_pdf'] = particle_diameter_pdf

    # This saves the particle diameter index (into 'scattering_irradiance') to
    # the parameters structure
    mie_scattering_data['particle_diameter_index_distribution'] = particle_diameter_index_distribution

    # This saves the scattering angle data into the parameters structure
    mie_scattering_data['scattering_angle'] = scattering_angle

    # This saves the scattering irradiance values for the different particle
    # diameters into the parameters structure
    mie_scattering_data['scattering_irradiance'] = scattering_irradiance

    # This saves the inverse rotation matrix into the parameters structure
    mie_scattering_data['inverse_rotation_matrix'] = inverse_rotation_matrix

    # This saves the normalized beam propogation direction vector into the
    # parameters structure
    mie_scattering_data['beam_propogation_vector'] = beam_propogation_vector

    return mie_scattering_data


def load_lightfield_data(piv_simulation_parameters,optical_system,mie_scattering_data,frame_index,lightfield_source):
# This function creates the lightfield data for performing the ray tracing
# operation.

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

    if perform_mie_scattering:
        # This extracts the particle diameter index (into 'scattering_irradiance')
        # from the parameters structure
        particle_diameter_index_distribution = mie_scattering_data['particle_diameter_index_distribution']
        # This sets the irradiance constant for using Mie scattering
        irradiance_constant = 500.0/1.0e4  # LKR added the extra 100 multiplication
    else:
        # This sets the particle diameter indices to a Null value
        particle_diameter_index_distribution= None
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

    # This is the number of source points that will simulated in one call to the GPU. this is a function
	# of the GPU memory, and threads-blocks-grids specifications
    lightfield_source['source_point_number'] = 10000

    # This is the number of rays to be generated for each source point
    lightfield_source['lightray_number_per_particle'] = 10000

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
        np.random.seed(714)
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

    # # % This is a vector of the x world coordinates of the grid points
    # x_grid_point_coordinate_vector = x_grid_point_spacing*np.linspace(-(x_grid_point_number-1.0)/2.0,(x_grid_point_number - 1.0)/2.0,x_grid_point_number,endpoint=True)
    # # % This is a vector of the y world coordinates of the grid points
    # y_grid_point_coordinate_vector = y_grid_point_spacing*np.linspace(-(y_grid_point_number-1.0)/2.0,(y_grid_point_number - 1.0)/2.0,y_grid_point_number,endpoint=True)

    #TODO remove
    current_z_world_coordinate = 0.0
    x_grid_point_number = 100
    y_grid_point_number = 100
    X_Min = -7.5e4
    X_Max = +7.5e4
    Y_Min = -7.5e4
    Y_Max = +7.5e4
    # x_grid_point_coordinate_vector = X_Min + (X_Max - X_Min)*np.random.rand(x_grid_point_number,1)
    # y_grid_point_coordinate_vector = Y_Min + (Y_Max - Y_Min)*np.random.rand(y_grid_point_number,1)
    x_grid_point_coordinate_vector = X_Min + (X_Max - X_Min) * np.random.rand(x_grid_point_number, y_grid_point_number)
    y_grid_point_coordinate_vector = Y_Min + (Y_Max - Y_Min) * np.random.rand(x_grid_point_number, y_grid_point_number)

    grid_point_diameter/=10
    particle_number_per_grid_point/=10

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

            # # % This is the coordinate of the current grid point
            # x_grid_point_coordinate = x_grid_point_coordinate_vector[x_grid_index]
            # y_grid_point_coordinate = y_grid_point_coordinate_vector[y_grid_index]

            #TODO remove
            # % This is the coordinate of the current grid point
            x_grid_point_coordinate = x_grid_point_coordinate_vector[x_grid_index,y_grid_index]
            y_grid_point_coordinate = y_grid_point_coordinate_vector[x_grid_index,y_grid_index]


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

    # % This eliminates any coordinates that were initilized but not changed
    if (count+1)<len(x):
        # This deletes the ending values from the coordinate vectors
        # HOW TO DELETE ELEMENTS IN PYTHON?
        x = x[:count+1]
        y = y[:count+1]

    # % This generates a series of points that fill the circle of one quarter of 
    # % the grid point diameter uniformly
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
    radiance = np.ones(x.shape)

    # % This rotates the image coordinates by the specified angles
    [x,y,z] = rotate_coordinates(x,y,z,x_camera_angle,y_camera_angle,0.0,0.0,0.0,0.0)

    # % % This rotates the particles by the specified angles
    # % [x_source,y_source,z_source]=rotate_coordinates(x_source,y_source,z_source,theta_x,theta_y,theta_z,0,0,0);
    # % This translates the Z coordinates of the paricles to the focal plane
    z = z + z_object

    # % This adds in the particles to the lightfield source data
    lightfield_source = {'x': x, 'y': y, 'z': z, 'radiance': radiance, 'diameter_index': np.ones(x.shape)}

    # This is the number of source points that will simulated in one call to the GPU. this is a function
	# of the GPU memory, and threads-blocks-grids specifications
    lightfield_source['source_point_number'] = 10000

    return lightfield_source

def generate_bos_lightfield_data(piv_simulation_parameters,optical_system):
# % This function generates a structure containing the location of the each of the
# % lightfield source points as well as their irradiance and color on a BOS texture.  Only a single value is
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

    # % This extracts the grid point diameter from the structure
    grid_point_diameter = piv_simulation_parameters['bos_pattern']['grid_point_diameter']
    # % This extracts the grid point number from the calibration structure
    x_grid_point_number = piv_simulation_parameters['bos_pattern']['x_grid_point_number']
    y_grid_point_number = piv_simulation_parameters['bos_pattern']['y_grid_point_number']
    # % This extracts the number of light source partciles per grid point from
    # % the calibration structure
    particle_number_per_grid_point = piv_simulation_parameters['bos_pattern']['particle_number_per_grid_point']

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

    # this perturbs the image distance by a small amount to see where the focal plane lies
    if('perturbation' in piv_simulation_parameters['lens_design'].keys()):
        image_distance += (1 + piv_simulation_parameters['lens_design']['perturbation']) * image_distance


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Generates the bos pattern data                                     %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    current_z_world_coordinate = 0.0

    if(piv_simulation_parameters['bos_pattern']['X_Min']):
        X_Min = piv_simulation_parameters['bos_pattern']['X_Min']
    else:
        X_Min = -7.5e4

    if(piv_simulation_parameters['bos_pattern']['X_Max']):
        X_Max = piv_simulation_parameters['bos_pattern']['X_Max']
    else:
        X_Max = +7.5e4


    if (piv_simulation_parameters['bos_pattern']['Y_Min']):
        Y_Min = piv_simulation_parameters['bos_pattern']['Y_Min']
    else:
        Y_Min = -7.5e4


    if (piv_simulation_parameters['bos_pattern']['Y_Max']):
        Y_Max = piv_simulation_parameters['bos_pattern']['Y_Max']
    else:
        Y_Max = +7.5e4

    x_grid_point_coordinate_vector = X_Min + (X_Max - X_Min) * np.random.rand(x_grid_point_number, y_grid_point_number)
    y_grid_point_coordinate_vector = Y_Min + (Y_Max - Y_Min) * np.random.rand(x_grid_point_number, y_grid_point_number)

    # % This generates a series of points that fill the circle of the grid point
    # % uniformly
    [x_lightray_coordinates,y_lightray_coordinates] = calculate_sunflower_coordinates(grid_point_diameter,particle_number_per_grid_point);

    # % These are the full coordinate vectors of all the lightrays
    if(piv_simulation_parameters['bos_pattern']['random_number_seed'][0]):
        np.random.seed(piv_simulation_parameters['bos_pattern']['random_number_seed'][0])
    x = np.zeros(int(particle_number_per_grid_point*x_grid_point_number*y_grid_point_number))

    if(piv_simulation_parameters['bos_pattern']['random_number_seed'][1]):
        np.random.seed(piv_simulation_parameters['bos_pattern']['random_number_seed'][1])
    y = np.zeros(int(particle_number_per_grid_point*x_grid_point_number*y_grid_point_number))

    # % This is a counting index
    count = 0

    # % This iterates through the grid points generating the light rays for each
    # % point
    for x_grid_index in range(0,int(x_grid_point_number)):
        for y_grid_index in range(0,int(y_grid_point_number)):

            # # % This is the coordinate of the current grid point
            # x_grid_point_coordinate = x_grid_point_coordinate_vector[x_grid_index]
            # y_grid_point_coordinate = y_grid_point_coordinate_vector[y_grid_index]

            #TODO remove
            # % This is the coordinate of the current grid point
            x_grid_point_coordinate = x_grid_point_coordinate_vector[x_grid_index,y_grid_index]
            y_grid_point_coordinate = y_grid_point_coordinate_vector[x_grid_index,y_grid_index]


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

    # % This eliminates any coordinates that were initilized but not changed
    if (count+1)<len(x):
        # This deletes the ending values from the coordinate vectors
        # HOW TO DELETE ELEMENTS IN PYTHON?
        x = x[:count+1]
        y = y[:count+1]

    # % This generates a series of points that fill the circle of one quarter of
    # % the grid point diameter uniformly
    [x_lightray_coordinates, y_lightray_coordinates] = calculate_sunflower_coordinates(grid_point_diameter/4,particle_number_per_grid_point/16);

    # % These are the coordinates of the first origin marker
    x_grid_point_coordinate = 0
    y_grid_point_coordinate = 0

    # % These are the coordinates of the current lightrays
    x_grid_point_lightray_coordinates = x_lightray_coordinates + x_grid_point_coordinate
    y_grid_point_lightray_coordinates = y_lightray_coordinates + y_grid_point_coordinate

    # % This adds the current lightrays to the total lightray vectors
    x = np.append(x,x_grid_point_lightray_coordinates)
    y = np.append(y,y_grid_point_lightray_coordinates)

    # % These are the coordinates of the second origin marker
    x_grid_point_coordinate = 0
    y_grid_point_coordinate = 0

    # % These are the coordinates of the current lightrays
    x_grid_point_lightray_coordinates = x_lightray_coordinates + x_grid_point_coordinate
    y_grid_point_lightray_coordinates = y_lightray_coordinates + y_grid_point_coordinate

    # % This adds the current lightrays tot he total lightray vectors
    x = np.append(x,x_grid_point_lightray_coordinates)
    y = np.append(y,y_grid_point_lightray_coordinates)

    # % This initializes the Z coordinate vector
    z = current_z_world_coordinate * np.ones(x.shape)
    # % This initializes the radiance vector
    radiance = np.ones(x.shape)*10

    # % This rotates the image coordinates by the specified angles
    [x,y,z] = rotate_coordinates(x,y,z,x_camera_angle,y_camera_angle,0.0,0.0,0.0,0.0)

    # % % This rotates the particles by the specified angles
    # % [x_source,y_source,z_source]=rotate_coordinates(x_source,y_source,z_source,theta_x,theta_y,theta_z,0,0,0);
    # % This translates the Z coordinates of the paricles to the focal plane
    z = z + z_object

    # % This adds in the particles to the lightfield source data
    lightfield_source = {'x': x, 'y': y, 'z': z, 'radiance': radiance, 'diameter_index': np.ones(x.shape)}

    # This is the number of source points that will simulated in one call to the GPU. this is a function
	# of the GPU memory, and threads-blocks-grids specifications
    lightfield_source['source_point_number'] = 10000

    return lightfield_source


def generate_bos_image_lightfield_data(piv_simulation_parameters,optical_system):
# % This function generates a structure containing the location of the each of the
# % lightfield source points as well as their irradiance and color on a BOS texture.  Only a single value is
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
    # % Generates the bos pattern data                                     %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    current_z_world_coordinate = 0.0

    X_Min = -7.5e4
    X_Max = +7.5e4
    Y_Min = -7.5e4
    Y_Max = +7.5e4

    # specify image filename
    filename = './bos_pattern_black_dot_936.png'
    # filename = './bos_pattern_white_dot.png'

    # read image array
    img = mpimg.imread(filename)

    # the image array will read in all 4 channels RGBA. just need the first channel because it is a grayscale image
    # and all the three RGB channels will be identical
    img = img[:, :, 0]

    # get the size of the image in pixels
    height, width = img.shape

    # estimate the pixel dimension in microns
    pixel_width = (X_Max - X_Min)/width

    # initialize the arrays that will hold the ray source co-ordinates
    x = np.zeros((1,1))
    y = np.zeros((1,1))
    radiance = np.zeros((1,1))
    # % This is a counting index
    count = 0

    # find pixel locations where the image intensity is non-zero
    rc = np.argwhere(img > 0)
    x = X_Min + (width - rc[:,1]) * pixel_width + pixel_width / 2.0
    y = Y_Max - (rc[:,0] * pixel_width + pixel_width / 2.0)
    radiance = img[img > 0].astype(np.float64)
    # radiance = np.ones(x.shape)

    # # iterate through the pixels in the image and emit light rays if the pixel intensity is non-zero
    # for row in range(0, height):
    #     for col in range(0, width):
    #         if(img[row,col] >= 0):
    #             # get x,y coordinates of the pixel
    #             x_current = X_Min + col * pixel_width + pixel_width / 2.0
    #             y_current = Y_Max - (row * pixel_width + pixel_width / 2.0)
    #
    #             if(count == 0):
    #                 x = x_current
    #                 y = y_current
    #                 radiance = img[row,col]
    #             else:
    #                 # append this position to the total light ray coordinate vector
    #                 x = np.append(x, x_current)
    #                 y = np.append(y, y_current)
    #                 radiance = np.append(radiance, img[row,col])
    #
    #             count += 1

    # % This initializes the Z coordinate vector
    z = current_z_world_coordinate * np.ones(x.shape)

    # # % This initializes the radiance vector
    # radiance = np.ones(x.shape)

    # % This rotates the image coordinates by the specified angles
    [x,y,z] = rotate_coordinates(x,y,z,x_camera_angle,y_camera_angle,0.0,0.0,0.0,0.0)

    # % % This rotates the particles by the specified angles
    # % [x_source,y_source,z_source]=rotate_coordinates(x_source,y_source,z_source,theta_x,theta_y,theta_z,0,0,0);
    # % This translates the Z coordinates of the paricles to the focal plane
    z = z + z_object

    # % This adds in the particles to the lightfield source data
    lightfield_source = {'x': x, 'y': y, 'z': z, 'radiance': radiance, 'diameter_index': np.ones(x.shape)}

    # This is the number of source points that will simulated in one call to the GPU. this is a function
	# of the GPU memory, and threads-blocks-grids specifications
    lightfield_source['source_point_number'] = 10000

    # if the source point number is greater than the number of particles, then use the use the number of particles
    # instead
    if(lightfield_source['source_point_number'] > x.size):
        lightfield_source['source_point_number'] = x.size

    return lightfield_source


def run_piv_simulation_02(piv_simulation_parameters):
# This function runs a PIV simulation using the thick lens, non-paraxial
# camera simulation.

    # % % This creates the simulation parameters for running the simulation
    # % piv_simulation_parameters=create_piv_simulation_parameters_02;
    optical_system = {}
    # This creates the optical system parameters used to simulate the defined lens system
    optical_system = create_camera_optical_system(piv_simulation_parameters)

    # # convert Nones to nans so you can save it as a MATLAB file. then reconvert them to null in MATLAB
    # optical_system['design']['optical_element'][0]['optical_element'][0]['element_number'] = np.NAN
    # optical_system['design']['optical_element'][0]['optical_element'][0]['elements_coplanar'] = np.NAN
    # optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties']['abbe_number'] = np.NAN
    # optical_system['design']['optical_element'][0]['element_geometry'] = np.NAN
    # optical_system['design']['optical_element'][0]['element_properties'] = np.NAN

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

    # generate_particle_field_images = False

    # % This generates the particle field images if specified by the parameters
    # % structure
    if generate_particle_field_images:

        # This displays that the particle images are being simulated
        print('\n\n')
        print('Simulating particle images . . . ')
      
        
        # This calculates the Mie scattering data if specified in the parameters
        # data structure, otherwise the Mie scattering data is set to a Null value
        if perform_mie_scattering:
            # This calculates the Mie scattering parameters data used in simulating the
            # scattering intensities of the simulated particles
            scattering_data = create_mie_scattering_data(piv_simulation_parameters)
            # This sets the scattering type to mie for the particle simulation
            scattering_type = 'mie'
        else:
            # This sets the Mie scattering data to a Null value
            scattering_data = None
            # This sets the scattering type to diffuse for the particle simulation
            scattering_type = 'diffuse'

        
         # This iterates through the frame vectors performing the ray tracing operations for each frame
        field_type = 'particle'
        # for frame_index in np.array(frame_vector):
        for frame_index in range(1,3):
            # This creates the lightfield data for performing the raytracing operation
            print "frame_index : %d" % frame_index
            lightfield_source = dict()
                        
            lightfield_source = load_lightfield_data(piv_simulation_parameters,optical_system,scattering_data,frame_index, lightfield_source)
                
            # convert none to NAN just for MATLAB
            if(scattering_data == None):
                scattering_data = np.NAN
            os.chdir(os.path.dirname(os.path.realpath(__file__)))

            # % This performs the ray tracing to generate the sensor image
            I, I_raw = perform_ray_tracing_03(piv_simulation_parameters,optical_system,pixel_gain,scattering_data,scattering_type,lightfield_source,field_type)

            # This is the filename to save the image data to
            image_filename_write = particle_image_directory + 'particle_image_frame_' + '%04d' % frame_index + '.tif'

            # This saves the image to memory
            TIFF.imsave(image_filename_write,I)

            # this is the filename to sve the raw image data to
            raw_image_filename_write = particle_image_directory + 'particle_image_frame_' + '%04d' % frame_index + '.bin'

            # this saves the image to memory
            I_raw.tofile(raw_image_filename_write)


    # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # # % Calibration Grid Simulation                                             %
    # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # % This extracts the Boolean value stating whether to generate the calibration
    # % images from the structure
    generate_calibration_grid_images = piv_simulation_parameters['calibration_grid']['generate_calibration_grid_images']
    # % This extracts the calibration plane number from the structure
    calibration_plane_number = int(piv_simulation_parameters['calibration_grid']['calibration_plane_number'])
    # % This extracts the directory to save the calibration grid images from
    # % parameters structure
    calibration_grid_image_directory = piv_simulation_parameters['output_data']['calibration_grid_image_directory']
    # % This extracts the number of lightrays to simulate per particle (this is roughly
    # % equivalent to the power of the laser)
    lightray_number_per_particle = piv_simulation_parameters['calibration_grid']['lightray_number_per_particle']
    # % This extracts the number of lightrays to propogate per iteration (this is a
    # % function of the RAM available on the computer)
    lightray_process_number = piv_simulation_parameters['calibration_grid']['lightray_process_number']
    # % This is the gain of the sensor in decibels to be used in the calibration
    # % grid simulation
    pixel_gain = piv_simulation_parameters['calibration_grid']['pixel_gain']

    # % This generates the calibration grid images if specified by the parameters
    # % structure
    if generate_calibration_grid_images:

        # % This displays that the calibration images are being simulated
        # fprintf('\n\n');
        print 'Simulating calibration images . . . '

        # % This sets the scattering type to diffuse for the calibration grid
        # % simulation
        scattering_type = 'diffuse'
        # % This sets the scattering data to a Null value for the calibration grid
        scattering_data = None

        field_type = 'calibration'

        #TODO remove
        calibration_plane_number = 1
        # % This iterates through the calibration grid planes performing the ray
        # % tracing operations for each plane
        for plane_index in range(1,calibration_plane_number+1):

            # % This creates the lightfield data for performing the raytracing operation
            lightfield_source = generate_calibration_lightfield_data(piv_simulation_parameters,optical_system,plane_index-1)
            # % This adds the number of lightrays per particle to the
            # % 'lightfield_source' data
            lightfield_source['lightray_number_per_particle'] = lightray_number_per_particle
            # % This adds the number of lightrays to simulateously process to the
            # % 'lightfield_source' data
            lightfield_source['lightray_process_number'] = lightray_process_number
            # save data to mat file
            # convert none to NAN just for MATLAB
            scattering_data = np.NAN

            # % This performs the ray tracing to generate the sensor image
            I, I_raw = perform_ray_tracing_03(piv_simulation_parameters,optical_system,pixel_gain,scattering_data,scattering_type,lightfield_source,field_type)

            # % This is the filename to save the image data to
            image_filename_write = calibration_grid_image_directory + 'calibration_image_plane_' + '%04d' % plane_index + '.tif'

            # % This saves the image to memory
            TIFF.imsave(image_filename_write, I)

            # this is the filename to save the raw image data to
            raw_image_filename_write = calibration_grid_image_directory + 'calibration_image_plane_' + '%04d' % plane_index + '.bin'

            # % This saves the image to memory
            I_raw.tofile(raw_image_filename_write)


    # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # # % BOS pattern Simulation                                             %
    # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # % This extracts the Boolean value stating whether to generate the calibration
    # % images from the structure
    generate_bos_pattern_images = piv_simulation_parameters['bos_pattern']['generate_bos_pattern_images']
    # % This extracts the directory to save the calibration grid images from
    # % parameters structure
    bos_pattern_image_directory = piv_simulation_parameters['output_data']['bos_pattern_image_directory']
    # % This extracts the number of lightrays to simulate per particle (this is roughly
    # % equivalent to the power of the laser)
    lightray_number_per_particle = piv_simulation_parameters['bos_pattern']['lightray_number_per_particle']
    # % This extracts the number of lightrays to propogate per iteration (this is a
    # % function of the RAM available on the computer)
    lightray_process_number = piv_simulation_parameters['bos_pattern']['lightray_process_number']
    # % This is the gain of the sensor in decibels to be used in the calibration
    # % grid simulation
    pixel_gain = piv_simulation_parameters['bos_pattern']['pixel_gain']

    # % This generates the calibration grid images if specified by the parameters
    # % structure
    if generate_bos_pattern_images:

        # % This displays that the calibration images are being simulated
        # fprintf('\n\n');
        print 'Simulating bos_pattern images . . . '

        # % This sets the scattering type to diffuse for the calibration grid
        # % simulation
        scattering_type = 'diffuse'
        # % This sets the scattering data to a Null value for the calibration grid
        scattering_data = None

        field_type = 'calibration'

        #TODO remove
        calibration_plane_number = 1

        # % This creates the lightfield data for performing the raytracing operation
        # lightfield_source = generate_bos_image_lightfield_data(piv_simulation_parameters,optical_system)
        lightfield_source = generate_bos_lightfield_data(piv_simulation_parameters, optical_system)

        # % This adds the number of lightrays per particle to the
        # % 'lightfield_source' data
        lightfield_source['lightray_number_per_particle'] = lightray_number_per_particle
        # % This adds the number of lightrays to simulateously process to the
        # % 'lightfield_source' data
        lightfield_source['lightray_process_number'] = lightray_process_number
        # save data to mat file
        # convert none to NAN just for MATLAB
        scattering_data = np.NAN

        ################################################################################################################
        # render the reference image without density gradients
        ################################################################################################################

        piv_simulation_parameters['density_gradients']['simulate_density_gradients'] = False
        # % This performs the ray tracing to generate the sensor image
        I, I_raw = perform_ray_tracing_03(piv_simulation_parameters,optical_system,pixel_gain,scattering_data,scattering_type,lightfield_source,field_type)

        # % This is the filename to save the image data to
        image_filename_write = bos_pattern_image_directory + 'bos_pattern_image_1.tif'

        # % This saves the image to memory
        TIFF.imsave(image_filename_write, I)

        # this is the filename to save the raw image data to
        raw_image_filename_write = bos_pattern_image_directory + 'bos_pattern_image_1.bin'

        # % This saves the image to memory
        I_raw.tofile(raw_image_filename_write)

        ################################################################################################################
        # render the image with density gradients
        ################################################################################################################

        piv_simulation_parameters['density_gradients']['simulate_density_gradients'] = True
        # % This performs the ray tracing to generate the sensor image
        I, I_raw = perform_ray_tracing_03(piv_simulation_parameters, optical_system, pixel_gain, scattering_data,
                                          scattering_type, lightfield_source, field_type)

        # % This is the filename to save the image data to
        image_filename_write = bos_pattern_image_directory + 'bos_pattern_image_2.tif'

        # % This saves the image to memory
        TIFF.imsave(image_filename_write, I)

        # this is the filename to save the raw image data to
        raw_image_filename_write = bos_pattern_image_directory + 'bos_pattern_image_2.bin'

        # % This saves the image to memory
        I_raw.tofile(raw_image_filename_write)

