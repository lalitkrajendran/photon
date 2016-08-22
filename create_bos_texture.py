# This program creates a BOS pattern consisting of black or white dots

import numpy as np
import matplotlib.pyplot as plt
# these are the minimum and maximum X,Y co-ordinates (microns) in the object space for the points on the pattern
X_min = -7.5e4
X_max = 7.5e4
Y_min = -7.5e4
Y_max = 7.5e4

# this is the Z co-ordinate of the pattern (in microns)
Z = 0.0

# this is the diameter of a dot in the bos pattern in microns
dot_diameter = 1.5e2

# this is the spacing for the grid over which the dots will be randomly distributed
grid_spacing = 5.0 * dot_diameter

# these are the number of grid points over X and Y
num_grid_points = (X_max - X_min)/grid_spacing

# these are the total number of dots to generate along each axis
num_dots = 500

# this randomly generates integers corresponding to the X and Y indices in the grid
xy_loc_int = np.random.randint(0, high=num_grid_points, size=(num_dots,2))

# this converts integer locations to positions in real space
xy_loc_float = np.zeros(xy_loc_int.shape, dtype=np.float32)
xy_loc_float[:,0] = X_min + xy_loc_int[:,0]*grid_spacing
xy_loc_float[:,1] = X_min + xy_loc_int[:,1]*grid_spacing


# this makes a scatter plot
fig = plt.figure(1)
ax = fig.add_subplot(111, axisbg='black')
ax.scatter(xy_loc_float[:,0], xy_loc_float[:,1], s=20, c='w', cmap='gray')

ax.set_xlim([X_min, X_max])
ax.set_ylim([Y_min, Y_max])
ax.set_autoscale_on(False)
ax.set_xticks([])
ax.set_yticks([])

plt.savefig('bos_scatter_test.png', bbox_inches='tight')
plt.show()
#plt.imsave('bos_scatter_test.png')