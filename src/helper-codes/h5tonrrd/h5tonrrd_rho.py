# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 19:01:43 2016

This program reads in an h5 file and does the following : 

1.) remaps the data onto a uniform grid

2.) converts data to dimensional values

3.) saves the density information in NRRD format for use in SchlierenRay.

this is a modified versino of h5tonrrd.py to only read and write the density only

Dependencies : pynrrd (https://pypi.python.org/pypi/pynrrd)

@author: lalit rajendran (lrajendr@purdue.edu)
"""
import numpy as np
import h5py
import glob
from scipy.interpolate import RegularGridInterpolator as rgi
import nrrd
import os
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

    # read co-ordinate data
    d['x'] = f.get('x')[...]        
    d['y'] = f.get('y')[...]
    d['z'] = f.get('z')[...]

    # transpose the density array if required
    if(d['rho'].shape[0] != d['x'].size):
        d['rho'] = np.transpose(d['rho'])

    # read the number of grid points along each co-ordinate axis
    if('Ni' in f.items()):
        d['Ni'] = f.get('Ni')[...]
    if ('Nj' in f.items()):
        d['Nj'] = f.get('Nj')[...]
    if('Nk' in f.items()):
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
    var = ['rho']

    # query data for various flow parameters and store in dataq
    for current_var in var:        
        # my_interpolating_function = rgi((x1,y1,z1), np.array(data1[current_var]))
        my_interpolating_function = rgi((x1, y1, z1), data1[current_var])
        # dataq[current_var] = my_interpolating_function(np.array([xq,yq,zq]).T) #.T)
        dataq[current_var] = my_interpolating_function(np.stack([xq, yq, zq], axis=-1))

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

    # calculate scaling that will make the dataset fit the field of view
    # shift data to 0,0,0
    data_ndim['x'] -= data_ndim['x'].mean()
    data_ndim['y'] -= data_ndim['y'].mean()
    data_ndim['z'] -= data_ndim['z'].mean()

    # this is the field of view of the camera in microns
    field_of_view = 100e3 #450e3 #1600e3 #100e3

    # scale the coordinates to fit the field of view
    scaling_distance = field_of_view / (data_ndim['x'].max() - data_ndim['x'].min())


    # define reference variables
    # density - kg/m^3, L - microns
    ref_var = {'rho' : 4 * 1.225, 'L' : scaling_distance}



    # calculate scaling values for the variables
    scale_var = {'rho' : np.float(ref_var['rho'])}
  
    var = ['rho']
    # convert variables to dimensional values
    for var_i in var:
        data_dim[var_i] = np.float(scale_var[var_i]) * np.array(data_ndim[var_i])

    # convert lengths to dimensional values
    data_dim['x'] = data_ndim['x'] * ref_var['L']
    data_dim['y'] = data_ndim['y'] * ref_var['L']
    data_dim['z'] = data_ndim['z'] * ref_var['L'] * 8 * (data_ndim['x'].max() - data_ndim['x'].min())/(data_ndim['z'].max() - data_ndim['z'].min())
          
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
    
    # space = '3D-right-handed'
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
    del_x = np.diff(np.squeeze(x))[0]
    del_y = np.diff(np.squeeze(y))[0]
    del_z = np.diff(np.squeeze(z))[0]

    spacings = np.array([del_x, del_y, del_z]).astype('float32')
    
    options = {'type' : 'f4', 'space': space, 'encoding': 'raw', 
               'space origin' : space_orig, 'spacings' : spacings}

    print "saving density to %s \n" % nrrd_filename    
    nrrd.write(nrrd_filename, np.array(data['rho']).astype('float32'), options)

# this is teh file path to read files from
# filepath = '/home/barracuda/a/lrajendr/Projects/SchlierenRayVis/Data/turbulent_channel_Mb=0.85_Reb=6900/'
read_filepath = '/home/shannon/a/aether/Projects/siv/analysis/data/JHU_dataset_new/256x256x128/new_h5_reshape_nz4/'

# this is the filepath to write files to
write_filepath = '/home/shannon/a/aether/Projects/siv/analysis/data/JHU_dataset_new/256x256x128/new_h5_reshape_nz4/nrrd/'

# populate a  list of filenames that end with .h5
filenames = glob.glob(read_filepath + '*.h5')

filenames = filenames[1:]
# calculate number of files to be read
Nt = len(filenames)
print "filepath: %s\n" % read_filepath
print "No. of files : %d\n" % Nt

# specify number of iterations to perform during debugging
# Nt = 0
#print 'Debug Mode, Nt = %d' % Nt


# Loop through all files, reading them in, converting them to nrrd, saving the data as mat and
# adding an execution command to the shell batch file.
for count,file in enumerate(filenames):
    if(count>Nt):
        break
    
    # obtain variables stored in the file
    data = get_vars_from_hdf5(file)
    print ('Acquired data from: %d %s' % (count+1, file))
    # data['rho'] = (np.flipud(data['rho']))
    data['rho'] = (np.flipud(np.transpose(data['rho'], axes=(1,0,2))))
    # reshape density data
    # data['rho'] = np.reshape(data['rho'], (256,192,128))
    
    # # for the first iteration, establish the mapping routine
    # if count == 0:
    #     xyz = {'x' : np.array(data['x']), 'y' : np.array(data['y']), 'z' : np.array(data['z'])}
    #     xyz_2 = generate_uniform_grid(xyz)
    #
    # # query for the values of the flow quantities at locations on the new grid
    # data_remapped = remap_data_on_query_points(data,xyz_2)

    # convert variables to dimensional values
    # data_scaled = scale_variables(data_remapped)
    data_scaled = scale_variables(data)

    # shift origin to center of volume
    data_scaled['x'] -= np.mean(data_scaled['x'])
    data_scaled['y'] -= np.mean(data_scaled['y'])
    # data_scaled['z'] -= np.mean(data_scaled['z'])

    # account for object distance
    object_distance = 2000e3 # 2500e3 #823668.800556 #2500e3
    # data_scaled['z'] = object_distance - np.flipud(data_scaled['z']) - 100e3
    data_scaled['z'] = object_distance - np.flipud(data_scaled['z'] - data_scaled['z'].min()) #- 100e3

    # save density as nrrd file
    nrrd_filename = write_filepath + os.path.basename(file)[0:-2] + 'nrrd' # 0:-2 goes till first to last but TWO elements. last index (-2) is ignored
    save_nrrd(data_scaled,nrrd_filename)
    

