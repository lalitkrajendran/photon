'''
this program reads in a set of positions from file and makes a scatter plot
'''

import numpy as np
import matplotlib.pyplot as plt
import pickle
import scipy.io as sio


# this is the file containing the ray positions from the matlab code
matlab_path = '/home/barracuda/a/lrajendr/Projects/photon/analysis/src/matlab_camera_simulation/from_kolmogorov/matlab_camera_simulation/final_pos_matlab.bin'

# this is the file containing the ray positions from the python cod
python_path = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/results/final_pos.bin'

# load positions from matlab file
pos_matlab = np.fromfile(matlab_path, dtype=np.float32)

# load positions from python file
pos_python = np.fromfile(python_path, dtype=np.float32)

# reshape arrays
pos_matlab = np.reshape(pos_matlab, (pos_matlab.size/2, 2))
pos_python = np.reshape(pos_python, (pos_python.size/2, 2))

# print number of points in each array
print "number of points: matlab - %d, python - %d" % (pos_matlab.shape[0], pos_python.shape[0])

# print arrays
print pos_matlab
print pos_python


# -----------------------------------------------------------------------------
# load particle positions from matlab and python
# -----------------------------------------------------------------------------
# this is the path to the file containing the particle positions used by the matlab code
matlab_path = '/home/barracuda/a/lrajendr/Projects/photon/analysis/src/matlab_camera_simulation/from_kolmogorov/matlab_camera_simulation/test_directory/particle_positions/'

#matlab_path = '/home/barracuda/a/lrajendr/Projects/camera_simulation_package_02/test_directory_250000/particle_positions/'
# this is the path to the file containing the particle positions used by the python code
python_path = '/home/barracuda/a/lrajendr/Projects/camera_simulation/test_directory_100/particle_positions/'

# This is the name of the file containing the particle positions
particle_filename = 'particle_data_frame_0001'

# this is the complete filename
python_filename = python_path + particle_filename + '.p'
matlab_filename = matlab_path + particle_filename + '.mat'

# display the filenames to the user
print 'python filename: ' + python_filename
print 'matlab filename: ' + matlab_filename

# load particle positions from python code
pos_python_particles = pickle.load(open(python_filename, 'rb'))
# load particle positions from matlab code
pos_matlab_particles = sio.loadmat(matlab_filename, squeeze_me = True)

# this is the dictionary containing the differences in positions for all three co-ordinates
diff = {'X': None, 'Y': None, 'Z': None } 
diff['X'] = np.squeeze(pos_python_particles['X']) - np.asarray(pos_matlab_particles['X'])
diff['Y'] = np.squeeze(pos_python_particles['Y']) - np.asarray(pos_matlab_particles['Y'])
diff['Z'] = np.squeeze(pos_python_particles['Z']) - np.asarray(pos_matlab_particles['Z'])

# make scatter plot
plt.figure(1)
plt.scatter(pos_matlab[:,0],pos_matlab[:,1],marker = '*', color = 'r')
plt.scatter(pos_matlab_particles['X'],pos_matlab_particles['Y'], facecolors = 'none', edgecolors = 'r')
plt.title('matlab')

plt.figure(2)
plt.scatter(pos_python[:,0],pos_python[:,1],marker = '*', color = 'b')
plt.scatter(pos_python_particles['X'],pos_python_particles['Y'], facecolors = 'none', edgecolors = 'b')
plt.title('python')

plt.figure(3)
plt.scatter(pos_matlab[:,0],pos_matlab[:,1],marker = '*', color = 'r')
plt.scatter(pos_python[:,0],pos_python[:,1],marker = '*', color = 'b')
plt.title('red - matlab, blue - python')


# make a line plot of ray positions (x,y and z) with their indices
plt.figure(4)
index_array = range(0,pos_matlab.shape[0])
particle_index_array = range(0, pos_matlab_particles['X'].size)
for i in range(0,pos_matlab_particles['X'].size):
    particle_index_array[i]*=10
plt.plot(index_array, pos_matlab[:,0], color='r', label='x')
#plt.plot(index_array, pos_matlab[:,1], color='g', label='y')
plt.plot(particle_index_array, pos_matlab_particles['X'], color='k', marker = '*', label = 'x_p')
#plt.plot(particle_index_array, pos_matlab_particles['Y'], color='g', marker = '*', label = 'y_p')
plt.legend()
plt.title('matlab')

plt.figure(5)
plt.plot(index_array, pos_python[:,0], color='r', label='x')
#plt.plot(index_array, pos_python[:,1], color='g', label='y')
plt.plot(particle_index_array, pos_python_particles['X'], color='k', marker = '*', label = 'x_p')
#plt.plot(particle_index_array, pos_python_particles['Y'], color='g', marker = '*', label = 'y_p')
plt.legend()
plt.title('python')

plt.show()
'''

# make scatter plot
plt.figure(1)
plt.scatter(pos_matlab[:,0],pos_matlab[:,1],marker = '*', color = 'r')
plt.title('matlab')

plt.figure(2)
plt.scatter(pos_python[:,0],pos_python[:,1],marker = '*', color = 'b')
plt.title('python')

plt.figure(3)
plt.scatter(pos_matlab[:,0],pos_matlab[:,1],marker = '*', color = 'r')
plt.scatter(pos_python[:,0],pos_python[:,1],marker = '*', color = 'b')
plt.title('red - matlab, blue - python')

plt.show()
'''
