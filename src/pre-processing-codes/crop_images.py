# this program reads in a set of images and crops them to remove the
# blurring at the edges due to real lens effects

import numpy as np
import sys
sys.path.append('/home/barracuda/a/lrajendr/Projects/camera_simulation/src/')

import tifffile as TIFF
import glob
import os.path as ospath
from progressbar import ProgressBar
from sys import argv

script, grad_x, density = argv

# this is the filepath where the images will be read
# read_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/processing/' + case_name + '/reordered-images/'
read_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/grad_x=' + grad_x + '/' + density + '/processing/reordered-images/'

# this is the filepath where the cropped images will be saved
# save_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/processing/' + case_name + '/cropped-images/'
save_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/grad_x=' + grad_x + '/' + density + '/processing/cropped-images/'

# these are the boundaries of the cropped region
xmin = 256
xmax = 768
ymin = 256
ymax = 768

# this loads a list of all the images in the directory
files = sorted(glob.glob(read_filepath + '*.tif'))

# this is the number of images in the directory
N = len(files);
print 'number of images: %d\n' % N

print 'beginning to crop images'

# this initializes the progress bar
pbar = ProgressBar(maxval=len(files)).start()

# this loops through all the files in the directory and crops them
for i in range(0,N):  
    # read the image
    img_ref = TIFF.imread(files[i])

    # this crops the image
    img_crop = img_ref[xmin:xmax, ymin:ymax]

    # this is the name of the image file without the path 
    [path, filename] = ospath.split(files[i])
    
    # this is the name of the file under which the image will be saved
    save_filename = filename
    
    # this saves the cropped image to file
    TIFF.imsave(save_filepath+save_filename, img_crop);
    
    # this displays the progress to the user
    pbar.update(i)

pbar.finish()

print 'finished cropping images\n'

'''
% % this displays the images before and after cropping
% T = zeros(1024,1024,3);
% T(:,:,1) = img_ref;
% T(xmin:xmax,ymin:ymax,2) = img_crop;
% figure(1)
% imshow(T)
'''
