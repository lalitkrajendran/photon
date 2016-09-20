import numpy as np
import matplotlib.pyplot as plt

filename = "image_array.bin"
filename2 = "image_array_const_grad_x_0_2.bin"

# read image from file
img = np.fromfile(filename,dtype=np.float64)
img2 = np.fromfile(filename2,dtype=np.float64)

# specify x and y pixel numbers
x_pixel_number = 1024
y_pixel_number = 1024

# reshape image
img = np.reshape(img,(x_pixel_number, y_pixel_number))
img2 = np.reshape(img2,(x_pixel_number, y_pixel_number))

# rescale images
img *= img/np.nanmax(img)*0.95*255
img2 *= img2/np.nanmax(img2)*0.95*255

# create variable to store each image in a single channel
t = np.zeros((x_pixel_number,y_pixel_number,3))
t[:,:,0] = img
t[:,:,1] = img2

#plt.imsave('image_comparison.png',t,vmin=0,vmax = 0.25*np.nanmax(img))
plt.imsave('image_comparison.png',t)


'''
# plot image
plt.imshow(img,cmap='gray',vmin=0,vmax=0.25*np.nanmax(img))
plt.show()
'''
