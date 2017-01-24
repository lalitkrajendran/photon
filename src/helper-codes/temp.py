import nrrd
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

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


#####################################################
# set object parameters
#####################################################

# this is the distance of the center of the refractive index medium from the camera lens (microns)
Z_D = 25e3

# this is the minimum x co-ordinate of the refractive index field
x_min = -5e4

# this is the maximum x co-ordinate of the refractive index field
x_max = 5e4

# this is the total field of view
del_x = x_max - x_min

# these are the number of points along the x, y, and z axis of the 3D volume
nx = 200
ny = 200
nz = 100

# this is the dimension of a pixel on the camera sensor (microns)
l_p = 17

# this is the Magnification
M = 0.1761

# this is the depth of the refractive index medium (microns)
W = 50e3

# this is the maximum allowable displacement in pixels at the camera sensor
D_max = 5.

# this is the distance of the object from the sensor
z_object = 823668.800556    

# this is the f number of the camera 
f_number = 2.8

# this is the focal length of the camera lens (microns)
f = 105e3

# this is the object distance (microns)
s = 700e3

# this is the required blur (pixels)
blur = 0.5

# calculate the viewing angle in radians
theta = np.arctan(f/(2*f_number*s)) 

# this is the glad stone dale constant (cm^3/g)
gladstone_dale = 0.226

# set the z bounds of the volume containing refractive index gradients
Z_Max = z_object - 1e2
Z_Min = z_object - 50e4


# calculate K3
K3 = (x_min + x_max)/2

# calculate max possible value of K2
K2_max = 1/2. * D_max * l_p /(Z_D * M * W) * 2 / del_x
print "K2_max: ", K2_max

# calculate max possible blur
blur_max = 4 / 3 * theta * W**3 / l_p * K2_max
print "max achievable blur (pixels): %.2f" % blur_max
print "actual blur (pixels): %.2f" % blur

# calculate actual value of K2
K2 = 3/4. * blur * l_p / (theta * W**3)
print "K2: ", K2
# calculate K1
K1 = 1 + K2 * del_x**2 / 4.

print "Coefficients: %.2e, %.2e, %.2e" % (K1, K2, K3)

#####################################################
# generate the refractive index field
#####################################################

# define a dictionary to hold the refractive index data
data = dict()

# create arrays to store coordiantes of the volume
data['x'] = np.linspace(start=x_min, stop=x_max, num=nx, endpoint=True).astype('float32')
data['y'] = np.linspace(start=x_min, stop=x_max, num=nx, endpoint=True).astype('float32')
data['z'] = np.linspace(start=Z_Min, stop=Z_Max, num=nz, endpoint=True).astype('float32')

x = np.array(data['x'])
y = np.array(data['y'])
z = np.array(data['z'])

X,Y = np.meshgrid(x,y,indexing='ij')

data['rho'] = np.ones([nx,ny,nz]).astype(np.float32)

for k in range(0,nz):
#    data['rho'][:,:,k] += grad_x * (X - x.min())/(X.max() - x.min()) + grad_xx * (X - x.min())**2/(X.max() - x.min())**2
     data['rho'][:,:,k] =  (K1 - K2 * np.power(x_loc - K3, 2)) / (gladstone_dale * 1e-3)
     

'''
plt.figure(1)
plt.plot(n)
plt.title("n")
plt.grid()
plt.draw()

# generate the displacement field
d = Z_D*M/l_p * W * np.diff(n)/np.diff(x_loc)

plt.figure(2)
plt.plot(d)
plt.title("d")
plt.grid()
plt.draw()

# verify that the constraints are met
print "n_min: %.5f" % np.min(n)
print "n_max: %.5f" % np.max(n)
print "d_min: %.2f" % np.min(d)
print "d_max: %.2f" % np.max(d)
'''
#plt.show()
