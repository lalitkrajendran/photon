#!/usr/bin/env python
"""
Created on Fri Feb 12 17:50:59 2016

This program reads in an h5 file and does the following : 

1.) remaps the data onto a uniform grid

2.) converts data to dimensional values

3.) saves the density information in NRRD format for use in SchlierenRay.

4.) saves all the variables at the corresponding datapoints in a MAT file

5.) generates a shell script file for batch processing 

Dependencies : pynrrd (https://pypi.python.org/pypi/pynrrd)

@author: lalit rajendran (lrajendr@purdue.edu)
"""
import numpy as np
import h5py
import glob
from scipy.interpolate import RegularGridInterpolator as rgi
import nrrd
import scipy.io as sio
import sys


def get_vars_from_hdf5(filepath):
# this function opens the h5 file specified by the filepath variable and 
# extracts its contents in the form of a dictionary. created by zhang.
    
    print "Reading data from file\n"
    # filename to be read
    f = h5py.File(filepath, 'r')
    
    d={}
    
    # read and store data in a dictionary
    d['rho'] = f.get('rho')[...]    # density
    d['u'] = f.get('u')[...]        # x component of velocity
    d['v'] = f.get('v')[...]        # y component of velocity
    d['w'] = f.get('w')[...]        # z component of velocity
    d['p'] = f.get('p')[...]        # pressure
    d['T'] = f.get('T')[...]        # temperature
    
    # read co-ordinate data
    d['x'] = f.get('x')[...]        
    d['y'] = f.get('y')[...]
    d['z'] = f.get('z')[...]

    # read the number of grid points along each co-ordinate axis
    d['Ni'] = f.get('Ni')[...]
    d['Nj'] = f.get('Nj')[...]
    d['Nk'] = f.get('Nk')[...]

    # read the time step and time information
    d['time'] = np.float(f.get('time')[...])
    d['time_step'] = np.int(f.get('time_step')[...])
    
    # close the file    
    f.close()
    
    return d


def remap_data_on_query_points(data1,xyz2):
# this function queries the dataset stored in the data1 dictionary at the 
# locations specified by the grid points xyz2 and stores the data in the
# variable dataq

    print "interpolating data to a uniform grid\n"
    x1 = np.array(data1['x'])
    y1 = np.array(data1['y'])
    z1 = np.array(data1['z'])

    x2 = np.array(xyz2['x'])
    y2 = np.array(xyz2['y'])
    z2 = np.array(xyz2['z'])

    dataq = {}

    # construct a new structured volume with uniform and equal spacing along
    # all the three axes
    xq,yq,zq = np.meshgrid(x2,y2,z2,indexing = 'ij')

    # create array of string to query the data dictionary
    var = ['rho', 'u', 'v', 'w', 'p', 'T']

    # query data for various flow parameters and store in dataq
    for current_var in var:        
        my_interpolating_function = rgi((x1,y1,z1), np.array(data1[current_var]))
        dataq[current_var] = my_interpolating_function(np.array([xq,yq,zq]).T)

    # store the grid points
    dataq['x'] = x2
    dataq['y'] = y2
    dataq['z'] = z2
    
    return dataq

def scale_variables(data_ndim):
# this function converts the non-dimensional values in the input data set to dimensional 
# values using user defined reference values

    print "scaling data to dimensional values\n"
    data_dim = {}

    # define reference variables
    ref_var = {'rho' : 1.225, 'T' : 300.0, 'L' : 0.1}
    const_var = {'gamma' : 1.4, 'R' : 287}
    ref_var['V'] = np.sqrt(const_var['gamma'] * const_var['R'] * ref_var['T'])

    # calculate scaling values for the variables
    scale_var = {'rho' : np.float(ref_var['rho']), 'u' : np.float(ref_var['V']), 
             'v' : np.float(ref_var['V']), 'w' : np.float(ref_var['V']), 
            'p' : np.float(ref_var['rho']) * np.float(ref_var['V'])**2,
            'T' : np.float(ref_var['T'])}    
  
    var = ['rho', 'u', 'v', 'w', 'p', 'T']
    # convert variables to dimensional values
    for var_i in var:
        data_dim[var_i] = np.float(scale_var[var_i]) * np.array(data_ndim[var_i])

    data_dim['x'] = data_ndim['x']
    data_dim['y'] = data_ndim['y']
    data_dim['z'] = data_ndim['z']
          
    return data_dim    

def generate_uniform_grid(xyz):
# this function converts the non-uniform co-ordinate grid co-ordinates to 
# have uniform spacing fo mapping

    x = np.array(xyz['x'])
    y = np.array(xyz['y'])
    z = np.array(xyz['z'])
    
    # obtain grid spacing along the three axes
    x_diff = np.diff(x)[0]
    y_diff = np.diff(y)[0]
    z_diff = np.diff(z)[0]
    print "x_diff: %f, y_diff: %f, z_diff: %f\n" % (x_diff, y_diff, z_diff)

    # set the new grid spacing to be the mean grid spacing 
    x_diff_2 = np.mean(np.diff(x))
    y_diff_2 = np.mean(np.diff(y))
    z_diff_2 = np.mean(np.diff(z))
    
    # find number of points to be generated with the new grid spacing
    N_x_2 = (np.max(x) - np.min(x))/x_diff_2 + 1
    N_y_2 = (np.max(y) - np.min(y))/y_diff_2 + 1
    N_z_2 = (np.max(z) - np.min(z))/z_diff_2 + 1
    
    # generate new uniformly spaced grid points along the three axes 
    x2 = np.linspace(np.min(x),np.max(x),N_x_2)
    y2 = np.linspace(np.min(y),np.max(y),N_y_2)
    z2 = np.linspace(np.min(z),np.max(z),N_z_2)

    xyz_2 = {'x' : np.array(x2), 'y' : np.array(y2), 'z' : np.array(z2)}

    return xyz_2

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
    nrrd.write(nrrd_filename, np.array(data_scaled['rho']).astype('float32'), options)

# specify filepath and filenames
filepath = '/home/barracuda/a/lrajendr/Projects/SchlierenRayVis/Data/turbulent_channel_Mb=0.85_Reb=6900/'
filenames = glob.glob(filepath + '*.h5')

# calculate number of files to be read
Nt = len(filenames)
print "filepath: %s\n" % filepath
print "No. of files : %d\n" % Nt

# specify number of iterations to perform during debugging
Nt = 0
#print 'Debug Mode, Nt = %d' % Nt


# specify name of batchfile for writing bash commands for batch processing
batch_file = open(filepath + 'batchfile.sh','w')


# Loop through all files, reading them in, converting them to nrrd, saving the data as mat and
# adding an execution command to the shell batch file.
for count,file in enumerate(filenames):
    if(count>Nt):
        break
    
    # obtain variables stored in the file
    data = get_vars_from_hdf5(file)
    print ('Acquired data from: %d %s' % (count+1, file))

    # for the first iteration, establish the mapping routine
    if count == 0:
        xyz = {'x' : np.array(data['x']), 'y' : np.array(data['y']), 'z' : np.array(data['z'])}
        xyz_2 = generate_uniform_grid(xyz)
    
    # query for the values of the flow quantities at locations on the new grid        
    data_remapped = remap_data_on_query_points(data,xyz_2)

    # convert variables to dimensional values
    data_scaled = scale_variables(data_remapped)
    
    # save density as nrrd file
    nrrd_filename = file[0:-2] + 'nrrd' # 0:-2 goes till first to last but TWO elements. last index (-2) is ignored
    save_nrrd(data_scaled,nrrd_filename)
    
    # save scaled, interpolated data as a mat file
    mat_filename = file[0:-2] + 'mat'   
    sio.savemat(mat_filename,data_scaled)

    # append command in shell file
    data_scalar = 1000.0
    cutoff_scalar = 10.0
    tiff_filename = file[0:-2] + 'tiff'
    batch_command = 'cuda-memcheck ./schlieren' + ' ' + nrrd_filename + ' ' + str(data_scalar) + ' ' + str(cutoff_scalar) + ' ' + tiff_filename + '\n'
    batch_file.write(batch_command)


batch_file.close()
                    

                    
    
    







    
