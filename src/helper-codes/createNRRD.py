'''
This program creates a 3D density field in NRRRD format for use in ray tracing.

The refractive index field is defined as n(x) = ax^2 + bx + c. 

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

#script, first, second, third = argv

data = {}
nx = 200
ny = 200
nz = 100
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

#print "del_x", data['x'][1] - data['x'][0]

'''
# convert co-ordinates from microns to meters
data['x']/=1.0e6
data['y']/=1.0e6
data['z']/=1.0e6
'''
# specify gradients for x and y directions
#grad_x = np.float64(first)

#grad_xx = 50.0 
#grad_x = 1.   # same as pixel displacement if Z_max - Z_min = 50e4
grad_x = 5.0
grad_y = 0.0

x = np.array(data['x'])
y = np.array(data['y'])
z = np.array(data['z'])

X,Y = np.meshgrid(x,y,indexing='ij')

data['rho'] = 1.225*np.ones([nx,ny,nz])

for k in range(0,nz):
    data['rho'][:,:,k] += grad_x * (X - x.min())/(X.max() - x.min()) #+ grad_xx * (X - x.min())**2/(X.max() - x.min())**2 # + grad_y*(Y - y.min())/Y.max()

#print "multiplying factor:", (X-x.min())/(X.max() - x.min())
#print "del_rho:", (data['rho'][1,0,0] - data['rho'][0,0,0])

#nrrd_filename = '/home/barracuda/a/lrajendr/Projects/ray_tracing_density_gradients/schlieren-0.2.0-Build/const_grad.nrrd'
#nrrd_filename = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_grad_x_08.nrrd'
nrrd_filename = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/const_grad_BOS_grad_x_%.2f_zmin_%02d_zmax_%02d.nrrd' % (grad_x, z_1/1e3, z_2/1e3)

save_nrrd(data,nrrd_filename)

# # estimate ray deflection through the volume
# del_rho = data['rho'][1,0,0] - data['rho'][0,0,0]
# del_n = 0.226*1e-6/1e-3*del_rho
# del_x = np.diff(data['x'])[0]
#
# ray_deflection = del_n/del_x * abs(Z_Max - Z_Min)
#
# print "del_x (microns) : %.2E, del_rho: %.2E, del_n: %.2E" % (del_x, del_rho, del_n)
# print "theoretical deflection (radians) : %0.2E" % ray_deflection

'''
plt.figure(3)
plt.plot(data['x'],data['rho'][:,:,1])
plt.title("density")
plt.grid()
plt.show()
'''

'''
=======

X_Velocity = 01.0e3
Y_Velocity = 01.0e3
Z_Velocity = 00.1e3

data['x'] = np.linspace(-15e4,15e4,nx).astype('float32') - X_Velocity/2
data['y'] = np.linspace(-15e4,15e4,ny).astype('float32') - Y_Velocity/2
# data['z'] = np.linspace(-15e3,15e3,nz).astype('float32') + z_object - Z_Velocity/2
data['z'] = np.linspace(7.5e5,9.0e5,nz).astype('float32')
# specify gradients for x and y directions
grad_x = 5e-4
grad_y = 0.0

x = np.array(data['x'])
y = np.array(data['y'])
z = np.array(data['z'])

X,Y = np.meshgrid(x,y,indexing='ij')

data['rho'] = 1.225*np.ones([nx,ny,nz])
for k in range(0,nz):
    data['rho'][:,:,k] += grad_x*(X-x.min())/X.max() # + grad_y*(Y - y.min())/Y.max()

# data['rho'] = data['rho'].T
# plot density profile
# plt.plot(data['rho'][:,10,0])
# plt.show()
save_nrrd(data,'const_grad.nrrd')
>>>>>>> Stashed changes
'''



