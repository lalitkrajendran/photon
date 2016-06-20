'''
This program reads a set of images from the specified directory, and replaces
them with a transposed version
'''


import numpy as np
import tifffile as TIFF

# type of image
image_type = 'calibration'

# top level image directory
image_directory = './test_directory/camera_images/' + image_type + '_images/'
print 'directory: ' + image_directory

# number of cameras
num_cameras = 5

# number of planes
num_planes = 7

# number of frames
num_frames = 2

print 'transposing images'

if(image_type == 'calibration'):
  for camera_index in range(1, num_cameras+1):
    print 'camera_index', camera_index
    for plane_index in range(0, num_planes):
      print 'plane_index', plane_index
      # file name to read
      image_name = image_directory + 'camera_' + '%02d' % camera_index + '/calibration_image_plane_' + '%04d' % plane_index + '.tif'
      # read image into an array
      img = TIFF.imread(image_name)
      # tranpose the array
      img = np.transpose(img)
      # save image to file
      TIFF.imsave(image_name,img)

else:
  for camera_index in range(1, num_cameras+1):
    print 'camera_index', camera_index
    for frame_index in range(1, num_frames+1):
      print 'frame_index', frame_index,
      # filename
      image_name = image_directory + 'camera_' + '%02d' % camera_index + '/particle_image_frame_' + '%04d' % frame_index + '.tif'
      # read image into an array
      img = TIFF.imread(image_name)
      # tranpose the array
      img = np.transpose(img)
      # save image to file
      TIFF.imsave(image_name,img)
    print '\n'

