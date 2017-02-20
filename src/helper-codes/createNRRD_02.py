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

#script, first, second, third = argv

data = {}
nx = 200
ny = 200
nz = 500
z_object = 2500e3 #823668.800556
#<<<<<<< Updated upstream

#Z_Min = -3.0e2
#Z_Max = 3.0e02

#Z_Max = z_object - 1e2
#Z_Min = z_object - 50e4
#z_1 = 50e3
#z_2 = 150e3

z_1 = 200e3 # 0.5e3
z_2 = 400e3 #30e3

Z_Max = z_object - z_1
Z_Min = z_object - z_2

# X_Velocity = 01.0e3
# Y_Velocity = 01.0e3
#Z_Velocity = 00.1e3

field_of_view = 175e3

data['x'] = np.linspace(-field_of_view/2,field_of_view/2,nx).astype('float32') # - X_Velocity/2
data['y'] = np.linspace(-field_of_view/2,field_of_view/2,ny).astype('float32') # - Y_Velocity/2
# data['z'] = np.linspace(-15e3,15e3,nz).astype('float32') + z_object - Z_Velocity/2
#data['z'] = np.linspace(8.7e5,9.5e5,nz).astype('float32')
#data['z'] = np.linspace(Z_Min,Z_Max,nz).astype('float32') + z_object
data['z'] = np.linspace(Z_Min, Z_Max, nz).astype('float32')

#grad_x = 1.   # same as pixel displacement if Z_max - Z_min = 50e4
grad_x_array = np.linspace(start=5.0, stop=5.0, num=1, endpoint=True)

for grad_x in grad_x_array:

    x = np.array(data['x'])
    y = np.array(data['y'])
    z = np.array(data['z'])

    X,Y = np.meshgrid(x,y,indexing='ij')

    data['rho'] = 1.225*np.ones([nx,ny,nz])

    for k in range(0,nz):
        data['rho'][:,:,k] += grad_x * (X - x.min())/(X.max() - x.min()) #+ grad_xx * (X - x.min())**2/(X.max() - x.min())**2 # + grad_y*(Y - y.min())/Y.max()

        # add noise to the density field
        data['rho'][:,:,k] += np.random.normal(0, scale=0.05*1.225, size=data['rho'][:,:,k].shape)

        data['rho'][:,:,k] += abs(data['rho'][:,:,k].min())

    nrrd_filename = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_grad_x_%.2f_zmin_%02d_zmax_%02d_nz_%04d.nrrd' % (grad_x, z_1/1e3, z_2/1e3, nz)

    calculate_theoretical_deflection(data, grad_x, z_1, z_2)
    save_nrrd(data,nrrd_filename)


