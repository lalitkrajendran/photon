import nrrd
import numpy as np
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

data = {}
nx = 200
ny = 200
nz = 100
z_object = 823668.800556

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
#plt.plot(data['rho'][:,10,0])
#plt.show()

nrrd_filename = '/home/barracuda/a/lrajendr/Projects/ray_tracing_density_gradients/schlieren-0.2.0-Build/const_grad.nrrd'
save_nrrd(data,nrrd_filename)




