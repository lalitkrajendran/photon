import nrrd
import numpy as np

def save_nrrd(data,nrrd_filename):

    space = '3D-right-handed'

    x = np.array(data['x'])
    y = np.array(data['y'])
    z = np.array(data['z'])

    # set origin
    x0 = x.min()
    y0 = y.min()
    z0 = z.min()

    space_orig = np.array([x0, y0, z0]).astype('float32')

    # set spacing
    del_x = np.diff(x)[0]
    del_y = np.diff(y)[0]
    del_z = np.diff(z)[0]

    spacings = np.array([del_x, del_y, del_z]).astype('float32')

    options = {'type' : 'f4', 'space': space, 'encoding': 'raw',
               'space origin' : space_orig, 'spacings' : spacings}

    print "saving density to %s \n" % nrrd_filename
    nrrd.write(nrrd_filename, np.array(data['rho']).astype('float32'), options)

data = {}
nx = 200
ny = 200
nz = 50
data['x'] = np.linspace(-15e4,15e4,nx).astype('float32')
data['y'] = np.linspace(-15e4,15e4,ny).astype('float32')
data['z'] = np.linspace(-15e3,15e3,nz).astype('float32')

data['rho'] = 1.225*np.ones([nx,ny,nz])

save_nrrd(data,'test.nrrd')




