'''
This program creates a BOS target consisting of white or black dots on a background, saves the pattern to file and also
calculates and saves the co-ordinates of the white points in the image based on a local co-ordinate system that can be
specified.
'''

import numpy as np
import matplotlib.pyplot as plt
import tifffile as tiff

# this is the path to the folder where the image will be stored
img_filepath = './images/'
# this is the name of the file that the image will be saved as
img_filename = 'bos_img.tif'

# this specifies the color of the dots ('b' - black, 'w' - white)
dot_color = 'w'

# this is the background color (opposite of the dot pattern color)
bgd_color = 'b' if (dot_color == 'w') else 'w'

# this is the number of dots to generate
num_dots = 10

# this is teh dimensions of the dot pattern
N = 1024

# create a set of random numbers for the location of the dots
bos_pattern = np.random.rand((N, N))

