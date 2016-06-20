'''
this program reads in the images from the matlab and python versions of the code,
saves the superposed images to a file, and computes the l2 norm of the intensity
difference between the two images
'''

import numpy as np
import sys

sys.path.append('/home/barracuda/a/lrajendr/Projects/camera_simulation/')

import tifffile as TIFF
from numpy import linalg as la
import os
import glob

# this is the path to the file containing the particle positions used by the matlab code
matlab_path = '/home/barracuda/a/lrajendr/Projects/camera_simulation_package_02/test_directory_250000/camera_images/'
# this is the path to the file containing the particle positions used by the python code
python_path = '/home/barracuda/a/lrajendr/Projects/camera_simulation/test_directory/camera_images/'

# this is a flag to specify if superposed images should be saved or not
save_images = True

# types of images to read
field_types = ['particle' , 'calibration']

# this is the number of cameras
num_cameras = 5

for field in field_types:
  for camera_index in range(1,num_cameras+1):
    # this is the folder name
    m_folder_name = matlab_path + field + '_images/' + 'camera_%02d' % camera_index + '/'
    py_folder_name = python_path + field + '_images/' + 'camera_%02d' % camera_index + '/'
        
    print 'camera: ', camera_index+1

    # load list of all bin files
    m_filename_list = sorted(glob.glob(m_folder_name + '*.tif'))
    #py_filename_list = sorted(glob.glob(py_folder_name + ))
    
    for file_index in range(0,len(m_filename_list)):
      # read matlab image
      I_m = TIFF.imread(m_filename_list[file_index])
      print 'matlab file: ', os.path.basename(m_filename_list[file_index])

      # read python image
      #I_py = TIFF.imread(py_filename_list[file_index])
      if(field == 'calibration'):
        I_py = TIFF.imread(py_folder_name + os.path.basename(m_filename_list[file_index])[:-5] + str(file_index) + '.tif')
        print 'python file: ', py_folder_name + os.path.basename(m_filename_list[file_index])[:-5] + str(file_index) + '.tif'

      else:
        I_py = TIFF.imread(py_folder_name + os.path.basename(m_filename_list[file_index]))
        print 'python file: ', py_folder_name + os.path.basename(m_filename_list[file_index])


      # compute difference between the two images
      diff = (I_m - I_py).astype(np.float32)
      # compute l2 norm of differences
      diff_norm = la.norm(diff)
      
      print "norm: %0.2G" %  diff_norm, 
      
      if(save_images):
        # define array to save superposed image
        I_superposed = np.zeros((1024,1024,3)).astype(np.uint16)
        # store matlab image in red channel
        I_superposed[:,:,0] = I_m
        # store python image in green channel
        I_superposed[:,:,1] = I_py
        
        # this is the name of the file where the image will be saved
        image_filename_write = py_folder_name + os.path.basename(m_filename_list[file_index])[:-4] + '_superposed.tif' 
        
        print 'writing to', image_filename_write
        
        # save image to file
        TIFF.imsave(image_filename_write, I_superposed)










