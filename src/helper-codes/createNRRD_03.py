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

#script, first, second, third = argv

data = {}

# number of points along the density gradient volume
nx = 200
ny = 200
nz = 100

# distance of the dot pattern from the camera lens (microns)
z_object = 2500e3 #823668.800556

# distance between dot pattern and near end of density gradient volume (microns)
z_1 = 200e3 # 0.5e3
# distance between dot pattern and far end of density gradient volume (microns)
z_2 = 400e3 # 30e3

# z co-ordinate of the density grad volume nearest to the dot pattern (microns)
Z_Max = z_object - z_1
# z co-ordinate of the density grad volume nearest to the dot pattern (microns)
Z_Min = z_object - z_2

# expected field of view for the given setting (microns)
field_of_view = 175e3

data['x'] = np.linspace(-field_of_view/2,field_of_view/2,nx).astype('float32') # - X_Velocity/2
data['y'] = np.linspace(-field_of_view/2,field_of_view/2,ny).astype('float32') # - Y_Velocity/2
data['z'] = np.linspace(Z_Min, Z_Max, nz).astype('float32')

# array of required displacements (pixels)
displacement_array = np.linspace(start=0.1, stop=2.0, num=20, endpoint=True)
print "displacement array", displacement_array

# magnification
M = 0.087

# this is the width of a pixel on the camera sensor (m)
pixel_pitch = 17e-6

# this is the density of the undisturbed medium (kg/m^3)
rho_0 = 1.225

for displacement_index in range(0,len(displacement_array)):

    # generate array of co-ordinates
    x = np.array(data['x'])
    y = np.array(data['y'])
    z = np.array(data['z'])

    X,Y = np.meshgrid(x,y,indexing='ij')

    # initialize array of densities
    data['rho'] = rho_0*np.ones([nx,ny,nz])

    # set pixel displacement values for the current volume
    delta_x = displacement_array[displacement_index]
    delta_y = displacement_array[displacement_index]
    
    # calculate the corresponding values of grad_x and grad_y
    [grad_x, grad_y] = calculate_density_gradient(delta_x, delta_y, data, rho_0, M, pixel_pitch, z_1, z_2)
    
    for k in range(0,nz):
        data['rho'][:,:,k] += grad_x * (X - x.min())/(X.max() - x.min()) + grad_y * (Y - y.min())/(Y.max() - y.min())

        # add noise to the density field
        data['rho'][:,:,k] += np.random.normal(0, scale=0.05*1.225, size=data['rho'][:,:,k].shape)

        data['rho'][:,:,k] += abs(data['rho'][:,:,k].min())

    nrrd_filename = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_delta_x_%.2f_delta_y_%.2f_zmin_%02d_zmax_%02d_nz_%04d.nrrd' % (delta_x, delta_y, z_1/1e3, z_2/1e3, nz)

    calculate_theoretical_deflection(data, grad_x, z_1, z_2)
    save_nrrd(data,nrrd_filename)


