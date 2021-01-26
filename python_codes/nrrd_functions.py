import nrrd
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import platform
import sys
import os

if platform.system() == 'Linux':
    mount_directory = '/scratch/shannon/c/aether/'
else:
    mount_directory = '/Volumes/aether_c/'

def save_nrrd(data, nrrd_filename):
    # This function saves the density stored in data['rho'] to an nrrd file
    # specified by nrrd_filename

    # specify orientation of space
    space = '3D-right-handed'

    # --------------------------
    # generate arrays to store co-ordinates along the three axes
    # --------------------------
    x = np.array(data['x'])
    y = np.array(data['y'])
    z = np.array(data['z'])

    # --------------------------
    # set origin
    # --------------------------
    x0 = x.min()
    y0 = y.min()
    z0 = z.min()

    space_orig = np.array([x0, y0, z0]).astype('float32')

    # --------------------------
    # set grid spacing
    # --------------------------
    del_x = np.diff(x)[0]
    del_y = np.diff(y)[0]
    del_z = np.diff(z)[0]

    spacings = np.array([del_x, del_y, del_z]).astype('float32')

    # --------------------------
    # spcify other relevant options
    # --------------------------
    options = {'type' : 'f4', 'space': space, 'encoding': 'raw',
               'space origin' : space_orig, 'spacings' : spacings}

    print("saving density to %s \n" % nrrd_filename)

    # --------------------------
    # save data to nrrd file
    # --------------------------
    nrrd.write(nrrd_filename, np.array(data['rho']).astype('float32'), options)


def calculate_theoretical_deflection(rho_grad, M, Z_D, del_z, rho_0, pixel_pitch): 
    # this function calculates the theoretical displacement of a dot pattern on the camera sensor for a constant density
    # gradient field using the paraxial approximation

    # this is the gladstone-dale constant (m^3/kg)
    K = 0.225e-3

    # this is the refractive index of the undisturbed medium
    n_0 = K * rho_0 + 1

    # this is the gradient of refractive index (m^-1)
    n_grad = K * rho_grad

    # angular deflection of light ray (radians)
    epsilon = 1/n_0 * n_grad * del_z

    # this is the displacement on the sensor (pixels)
    displacement = M * Z_D * epsilon / pixel_pitch

    print('density gradient (kg/m^4): %.2G' % rho_grad)
    print('refractive index gradient (/m): %.4G' % n_grad)
    print('angular deflection of ray (radians): %.6G' % epsilon)
    print('displacement on the sensor (pixels): %.2G' % displacement)


def calculate_density_gradient_from_deflection(disp, M, Z_D, del_z, rho_0, pixel_pitch):
    '''
    this function calculates the density gradients required to produce the specified
    pixel displacement at the camera sensor given the extent of the density gradient
    volume, magnification and pixel pitch

    '''

    # this is the gladstone-dale constant (m^3/kg)
    K = 0.225e-3
    
    # this is the refractive index of the undisturbed medium
    n_0 = K * rho_0 + 1
    
    # this is the required angular deflection of the ray (radians)
    epsilon = disp * pixel_pitch/(Z_D * M)
    
    # this is the required refractive index gradient
    n_grad = epsilon * n_0/del_z
        
    # this is the required gradient of density along x (kg/m^4)
    rho_grad = 1/K * n_grad

    return rho_grad


def calculate_nrrd_extents(density_field):

    # extract minimum value of co-ordinates
    min_loc = density_field[1]['space origin']
    # extract maximum value of co-ordinates
    max_loc = density_field[1]['space origin'] + (density_field[1]['sizes'] - 1) * density_field[1]['spacings']
    
    # extract number of elements
    sizes = density_field[1]['sizes']
    
    # extract spacing
    spacing = density_field[1]['spacings']
    
    # calculate extent of the whole density field
    extents = (sizes - 1) * spacing

    return extents