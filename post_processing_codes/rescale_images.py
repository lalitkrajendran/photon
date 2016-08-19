'''
this program rescales the raw image files containing floating point values
to the range (0, 2^bit_depth-1) and saves them as TIFF files
'''

import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('/home/barracuda/a/lrajendr/Projects/camera_simulation/')
import tifffile as tiff

filepath = '/home/shannon/b/aether/Projects/piv-development/photon/analysis/data/test_sayantan_parameters/images/25e4/camera_images/particle_images/'

num_cameras = 1
num_frames = 1

pixel_gain = 25.0
pixel_bit_depth = 10

for camera_index in range(1, num_cameras+1):
  for frame_index in range(1, num_frames+1):
    # this displays the camera number and the frame number of the file that will be read
    print 'camera : %d, frame: %d' % (camera_index, frame_index)

    # this is the name of the binary file containing the raw intensity values
    filename = 'camera_0%d' % camera_index + '/particle_image_frame_000%d.bin' % frame_index
    
    # this reads the image from the file into an array
    I = np.fromfile(filepath+filename, dtype=np.float32)
   
    # this reshapes the array to a 2D matrix
    I = np.reshape(I, (1024, 1024))

    # This rescales the image intensity to account for the pixel gain
    I *= 10**(pixel_gain/20.0)

    # This rescales the image to values ranging from 0 to 2^bit_depth-1
    I = (2**pixel_bit_depth - 1) * I / (np.nanmax(I) - 1.0)          #(2**16 - 1.0)

    # This rounds the values of the image to the nearest integer
    I = np.round(I)
    
    # This rescales the image to the full 16 bit range (but only bit_depth bits
    # of information are now used)
    I *= (2**16 - 1.0) / (2**pixel_bit_depth - 1.0)

    # This converts the image from double precision to 16 bit precision
    I2 = I.astype(np.uint16)

    # this is the filename of the tiff file where the image will be saved
    tiff_filename = 'camera_0%d' % camera_index + '/particle_image_rescaled_frame_000%d.tif' % frame_index

    # this saves the image to file
    tiff.imsave(filepath + tiff_filename, I2)


