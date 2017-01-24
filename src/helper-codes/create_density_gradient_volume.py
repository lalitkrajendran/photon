'''
This program creates a 3D density field in NRRRD format for use in ray tracing.

The refractive index field is defined as n(x) = K1 - K2*(x - K3)^2 

'''

import nrrd
import numpy as np
import numpy.linalg as la
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

# define a dictionary to store the density gradient data
data = {}

# this is the number of grid points along the x, y, and z axis for the volume
nx = 200
ny = 200
nz = 100

# this is the distance between the bos pattern and the camera sensor
z_object = 823668.800556

# these are the min and max bounds for the volume in the camera fitted system
Z_Max = z_object - 1e2
Z_Min = z_object - 50e4

# this is the distance of the center of the refractive index medium from the camera (microns)
Z_D = 25e3

# this is the dimension of a pixel on the camera sensor (microns)
l_p = 17

# this is the Magnification
M = 0.1761

# this is the depth of the refractive index medium (microns)
W = 50e3

# this is the maximum allowable displacement in pixels at the camera sensor
D_max = 5.

# this is the minimum x co-ordinate of the refractive index field
x_min = -5e4

# this is the maximum x co-ordinate of the refractive index field
x_max = 5e4

# this is the number of grid points along x
x_num = nx

# this is the minimum x co-ordinate of the refractive index field that is parabolic
x_parabolic_min = -10e3

# this is the maximum x co-ordinate of the refractive index field that is parabolic
x_parabolic_max = 10e3

# this is the number of grid points along x in the parbolic region
x_parabolic_num = int(np.round(x_num * 1.0 * (x_parabolic_max - x_parabolic_min) / (x_max - x_min)))

# this is the total field of view
del_x = x_parabolic_max - x_parabolic_min

# this is the f number of the camera 
f_number = 8

# this is the focal length of the camera lens (microns)
f = 105e3

# this is the object distance (microns)
s = 700e3

# this is the required blur (pixels)
blur = .1

# this is the Gladstone-Dale constant for air (m^3/kg)
K = 0.226*1e-3

# this is the density of air at STP (kg/m^3)
rho_air = 1.225

# this is the refractive index of air at STP
n_air = 1 + K*rho_air

# calculate the viewing angle in radians
theta = np.arctan(f/(2*f_number*s)) 

# calculate K3
K3 = (x_min + x_max)/2

# calculate max possible value of K2
K2_max = 1/2. * D_max * l_p /(Z_D * M * W) * 2 / del_x
print "K2_max: ", K2_max

# calculate max possible blur
blur_max = 4 / 3. * theta * W**3 / l_p * K2_max
print "max achievable blur (pixels): %.2f" % blur_max
print "actual blur (pixels): %.2f" % blur

# calculate actual value of K2
K2 = 3/4. * blur * l_p / (theta * W**3)
print "K2: ", K2

# calculate K1
K1 = n_air + K2 * del_x**2 / 4.
print "Coefficients: %.2e, %.2e, %.2e" % (K1, K2, K3)

# generate the co-ordinate along the x direction
x_1 = np.linspace(start=x_min, stop=x_parabolic_min, num=np.round((x_num - x_parabolic_num)/2.), endpoint=False)
x_parabolic = np.linspace(start=x_parabolic_min, stop=x_parabolic_max, num=x_parabolic_num, endpoint=False)
x_2 = np.linspace(start=x_parabolic_max, stop=x_max, num=np.round((x_num - x_parabolic_num)/2.), endpoint=False)

x = np.append(np.append(x_1, x_parabolic), x_2).astype(np.float32)

# generate the co-ordinate along the y direction
y = x

# generate the co-ordinate along the z direction
z = np.linspace(Z_Min, Z_Max, nz).astype('float32')

n_1 = n_air * np.ones(np.shape(x_1))
n_parabolic = K1 - K2 * np.power(x_parabolic - K3, 2)
n_2 = n_air * np.ones(np.shape(x_2))

n = np.append(np.append(n_1, n_parabolic), n_2).astype(np.float32)

# generate the displacement field
d = Z_D*M/l_p * W * np.diff(n)/np.diff(x)
d = np.append(0,d)

# generate the density field (kg/m^3)
rho = (n - 1)/K

X,Y = np.meshgrid(x,y,indexing='ij')

data['rho'] = 1.225*np.ones([nx,ny,nz]).astype(np.float32)

for k in range(0,nz):
    for j in range(0,ny):
    #data['rho'][:,:,k] += grad_x * (X - x.min())/(X.max() - x.min()) + grad_xx * (X - x.min())**2/(X.max() - x.min())**2 # + grad_y*(Y - y.min())/Y.max()
        data['rho'][:,j,k] = rho.astype(np.float32)
        #data['rho'][:,j,k] = 0.00001 * (x - x.min()) + 1.225*np.ones(x.shape)
        
nrrd_filename = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/blur_%.2f_f%d.nrrd' % (blur, f_number)

data['x'] = x.astype('float32')
data['y'] = y.astype('float32')
data['z'] = z.astype('float32')

save_nrrd(data,nrrd_filename)

plt.figure(3)
plt.plot(data['x'],data['rho'][:,:,1], label='volume')
#plt.plot(x,rho, label='line')
plt.legend()
plt.title("density")
plt.grid()
plt.draw()

# verify that the constraints are met
print "n_min: %.5f" % np.min(n)
print "n_max: %.5f" % np.max(n)
print "d_min: %.2f" % np.min(d)
print "d_max: %.2f" % np.max(d)
print "rho_min: %.2f" % np.min(rho)
print "rho_max: %.2f" % np.max(rho)

plt.show()


# estimate ray deflection through the volume
del_rho = data['rho'][1,0,0] - data['rho'][0,0,0]
del_n = 0.226*1e-6/1e-3*del_rho
del_x = np.diff(data['x'])[0]

ray_deflection = del_n/del_x * abs(Z_Max - Z_Min)

print "del_x (microns) : %.2E, del_rho: %.2E, del_n: %.2E" % (del_x, del_rho, del_n)
print "theoretical deflection (radians) : %0.2E" % ray_deflection

