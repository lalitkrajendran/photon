# this program compares two images to show the effect of density gradients
from skimage import io,color
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('/home/barracuda/a/lrajendr/Projects/camera_simulation/src/')

import tifffile as TIFF

# load images
# img1 = TIFF.imread('/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dns/0/2/bos_pattern_image_1.tif')
# img2 = TIFF.imread('/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dns/0/2/bos_pattern_image_2.tif')

# filepath containing images
img_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dof/f16/'

# load images
img1 = TIFF.imread(img_filepath + 'bos_pattern_image_1.tif')
img2 = TIFF.imread(img_filepath + 'bos_pattern_image_2.tif')

# convert rgb images to grayscale
img1 = color.rgb2gray(img1)
img2 = color.rgb2gray(img2)

# create multi color image to see the difference
t = np.zeros((1024,1024,3))

# store first image in red channel
t[:,:,0] = img1

# store second image in green channel
t[:,:,1] = img2

# display image
#plt.imshow(t)
plt.imshow(t)
# plt.imsave('bos_comparison.png',t)
plt.show()


