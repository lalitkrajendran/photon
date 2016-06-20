'''
this program compares the particle positions used for ray tracing by the matlab
and the python versions of the camera simulation code
'''

import numpy as np
import matplotlib.pyplot as plt
import pickle 
import scipy.io as sio
import glob
from numpy import linalg as la

# this is the path to the file containing the particle positions used by the matlab code
matlab_path = '/home/barracuda/a/lrajendr/Projects/camera_simulation_package_02/test_directory_250000/particle_positions/'
# this is the path to the file containing the particle positions used by the python code
python_path = '/home/barracuda/a/lrajendr/Projects/camera_simulation/test_directory/particle_positions/'

# This is the name of the file containing the particle positions
particle_filename = 'particle_data_frame_0001'

# this is the complete filename
python_filename = python_path + particle_filename + '.p'
matlab_filename = matlab_path + particle_filename + '.mat'

# display the filenames to the user
print 'python filename: ' + python_filename
print 'matlab filename: ' + matlab_filename

# load particle positions from python code
pos_python = pickle.load(open(python_filename, 'rb'))
# load particle positions from matlab code
pos_matlab = sio.loadmat(matlab_filename, squeeze_me = True)

# this is the dictionary containing the differences in positions for all three co-ordinates
diff = {'X': None, 'Y': None, 'Z': None } 
diff['X'] = np.squeeze(pos_python['X']) - np.asarray(pos_matlab['X'])
diff['Y'] = np.squeeze(pos_python['Y']) - np.asarray(pos_matlab['Y'])
diff['Z'] = np.squeeze(pos_python['Z']) - np.asarray(pos_matlab['Z'])

# this is the dictionary containing the l2 norm of the differences between the 
# two sets of vectors
diff_norm = {'X': None, 'Y': None, 'Z': None }
diff_norm['X'] = la.norm(diff['X'])
diff_norm['Y'] = la.norm(diff['Y'])
diff_norm['Z'] = la.norm(diff['Z'])


'''
for key in pos_python:
  print 'key: ', key
  # calculate difference between the two sets of position vectors
  diff[key] = pos_python[key] - pos_matlab[key]
  # calculate the l2 norm of the differences
  diff_norm[key] = la.norm(diff[key])
'''

# display the results
print 'l2 norm: del_X: %.2G, del_Y: %.2G, del_Z: %.2G' % (diff_norm['X'], diff_norm['Y'], diff_norm['Z'])

