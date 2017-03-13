'''
This program creates a 3D density field in NRRRD format for use in ray tracing.

'''

import nrrd
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

def save_nrrd(data,nrrd_filename):
# This function saves the density stored in data['rho'] to an nrrd file
# specified by nrrd_filename

    # specify orientation of space
    space = '3D-right-handed'

    # generate arrays to store co-ordinates along the three axes
    x = np.array(data['x'])
    y = np.array(data['y'])
    z = np.array(data['z'])

    # set origin
    x0 = x.min()
    y0 = y.min()
    z0 = z.min()

    space_orig = np.array([x0, y0, z0]).astype('float32')

    # set grid spacing
    del_x = np.diff(x)[0]
    del_y = np.diff(y)[0]
    del_z = np.diff(z)[0]

    spacings = np.array([del_x, del_y, del_z]).astype('float32')

    # spcify other relevant options
    options = {'type' : 'f4', 'space': space, 'encoding': 'raw',
               'space origin' : space_orig, 'spacings' : spacings}

    print "saving density to %s \n" % nrrd_filename

    # save data to nrrd file
    nrrd.write(nrrd_filename, np.array(data['rho']).astype('float32'), options)

def calculate_theoretical_deflection(data, grad_x, z_1, z_2):
# this function calculates the theoretical displacement of a dot pattern on the camera sensor for a constant density
# gradient field using the paraxial approximation

    # distance between center of density gradient and dot pattern (m)
    Z_D = (z_1+z_2)/2 * 1e-6

    # magnification
    M = 0.087

    # this is the density of the undisturbed medium (kg/m^3)
    rho_0 = 1.225

    # this is the gladstone-dale constant (m^3/kg)
    K = 0.225e-3

    # this is the refractive index of the undisturbed medium
    n_0 = K * rho_0 + 1

    # this is the density gradient (kg/m^4)
    drho_dx = grad_x/((data['x'].max() - data['x'].min())*1e-6)

    # this is the thickness of the density gradient volume
    del_z = (data['z'].max() - data['z'].min())*1e-6

    # this is the gradient of refractive index (m^-1)
    dn_dx = K * drho_dx

    # angular deflection of light ray (radians)
    epsilon_x = 1/n_0 * dn_dx * del_z

    # this is the displacement on the sensor (m)
    displacement = Z_D * M * epsilon_x

    # this is the width of a pixel on the camera sensor (m)
    pixel_pitch = 17e-6

    print 'angular deflection of ray (radians): %.2G' % epsilon_x
    print 'displacement on the sensor (mm): %.2G' % (displacement*1e3)
    print 'displacement on the sensor (pixels): %.2G' % (displacement/pixel_pitch)


def calculate_density_gradient(delta_x, delta_y, data, rho_0, M, pixel_pitch, z_1, z_2):
    '''
    this function calculates the density gradients required to produce the specified
    pixel displacement at the camera sensor given the extent of the density gradient
    volume, magnification and pixel pitch

    '''
    
    # distance between center of density gradient and dot pattern (m)
    Z_D = (z_1+z_2)/2 * 1e-6
    
    # this is the thickness of the density gradient volume (m)
    del_z = (data['z'].max() - data['z'].min()) * 1e-6
    
    # this is the gladstone-dale constant (m^3/kg)
    K = 0.225e-3
    
    # this is the refractive index of the undisturbed medium
    n_0 = K * rho_0 + 1
    
    # this is the required angular deflection of the ray (radians)
    epsilon_x = delta_x * pixel_pitch/(Z_D * M)
    
    # this is the required refractive index gradient
    dn_dx = epsilon_x * n_0/del_z
    
    # this is the required density gradient (kg/m^4)
    drho_dx = 1/K * dn_dx

    # # this is teh required grad_x
    # grad_x = drho_dx * ((data['x'].max() - data['x'].min())*1e-6)
    # grad_y = grad_x
    #
    # print 'grad_x', grad_x
    # this is the required gradient of density along x (kg/m^4)
    grad_x = 1/K * delta_x * pixel_pitch/(M * Z_D) * n_0/del_z * ((data['x'].max() - data['x'].min())*1e-6)

    # this is the required gradient of density along x (kg/m^4)
    grad_y = 1/K * delta_y * pixel_pitch/(M * Z_D) * n_0/del_z * ((data['y'].max() - data['y'].min())*1e-6)

    return grad_x, grad_y

# dictionary to store density field
data = {}

# number of points along the density gradient volume
nx = 10
ny = 10
nz = 10

# generate array of z co-ordinates
x = np.linspace(start=-1e3, stop=1e3, num=nx, endpoint=True).astype(np.float32)
y = np.linspace(start=-1e3, stop=1e3, num=ny, endpoint=True).astype(np.float32)
z = np.linspace(start=0, stop=5, num = nz, endpoint=True).astype(np.float32)

# this is the density of the undisturbed medium (kg/m^3)
rho_0 = 1.225

# this is the gladstone-dale constant (m^3/kg)
K = 0.225e-3

# this is the refractive index of the undisturbed medium
n_0 = K * rho_0 + 1

# gradient of the square of the refractive index
alpha = 1e-1

# array of refractive indices
refractive_index_array = np.sqrt(n_0**2 + alpha * z)

# generate array of co-ordinates
data['x'] = x
data['y'] = y
data['z'] = z

X,Y = np.meshgrid(x,y,indexing='ij')

# initialize array of densities
data['rho'] = rho_0*np.ones([nx,ny,nz])

for k in range(0,nz):
    data['rho'][:,:,k] = (refractive_index_array[k] - 1)/K

nrrd_filename = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/n2linear_alpha%.2f_nz%d.nrrd' % (alpha, nz)

save_nrrd(data,nrrd_filename)


