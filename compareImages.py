# this program compares two images to show the effect of density gradients
from skimage import io,color
import numpy as np
import matplotlib.pyplot as plt

# load images
img1 = io.imread('no_grad.png')
img2 = io.imread('const_grad.png')

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
plt.imshow(t)
plt.show()


