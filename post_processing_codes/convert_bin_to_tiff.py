'''
this program reads in a set of .bin files containing the image array,
and saves them as a 16 bit tiff file after contrast stretching
'''

import numpy as np
import glob
import sys
import matplotlib.pyplot as plt

sys.path.append('/home/barracuda/a/lrajendr/Projects/camera_simulation/')
import tifffile as TIFF

from skimage import exposure

# specify directory contatining images
read_directory = '/home/barracuda/a/lrajendr/Projects/camera_simulation/test_directory_100000_bream_vector_z/camera_images/'
#read_directory = '/home/barracuda/a/lrajendr/Projects/photon/analysis/src/matlab_camera_simulation/from_kolmogorov/matlab_camera_simulation/test_directory_250000/camera_images/'

# types of images to read
#field_types = ['particle' , 'calibration']
field_types = ['particle']
# this is the number of cameras
num_cameras = 1

# this is the pixel_gain
pixel_gain = 25.0

# this is the pixel_bit_depth
pixel_bit_depth = 16.0

for field in field_types:
  for camera_index in range(1,num_cameras+1):
    # this is the folder name
    folder_name = read_directory + field + '_images/' + 'camera_%02d' % camera_index + '/'
    
    print 'reading from ' + folder_name

    # load list of all bin files
    filename_list = sorted(glob.glob(folder_name + '*.bin'))
    print filename_list
    for filename in filename_list:
      #print 'read_file : ' + filename
      # read file contents
      I = np.fromfile(filename, dtype = np.float32)
      
      # reshape array
      I = np.reshape(I, (1024,1024))

      # transpose array since it was saved in matlab
      I = np.transpose(I)

      # % This rescales the image intensity to account for the pixel gain
      I *= 10 ** (pixel_gain/ 20.0)

      # % This rescales the image to values ranging from 0 to 2^bit_depth-1
      I = (2 ** pixel_bit_depth - 1) * I / (np.nanmax(I) - 1)

      # % This rounds the values of the image to the nearest integer
      I = np.round(I)

      # % This rescales the image to the full 16 bit range (but only bit_depth bits
      # % of information are now used)
      I *= (2 ** 16 - 1.0) / (2 ** pixel_bit_depth - 1.0)

      # % This converts the image from double precision to 16 bit precision
      I = np.uint16(I)

      plt.imshow(I,cmap='gray')
      plt.show()
      '''
      # Contrast stretching (http://homepages.inf.ed.ac.uk/rbf/HIPR2/stretch.htm)

      # These are the percentile thresholds of the intensity values to rescale the image
      threshold_lo = 0
      threshold_hi = 99.8

      # These are the intensities correspoding to the high and low thresholds
      p_lo, p_hi = np.percentile(I, (threshold_lo, threshold_hi))

      # This rescales the image
      I = exposure.rescale_intensity(I, in_range=(p_lo, p_hi))

      # this is the name of the file to save the image
      image_filename_write = filename[:-3] + 'tif'

      print 'writing to ' + image_filename_write
      
      #print 'I.shape: ', I.shape
      # this saves the image to memory
      TIFF.imsave(image_filename_write, I)
      '''
