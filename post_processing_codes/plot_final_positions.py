'''
this program reads in a set of positions from file and makes a scatter plot
'''

import numpy as np
import matplotlib.pyplot as plt


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

# print arrays
print pos_matlab
print pos_python




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
