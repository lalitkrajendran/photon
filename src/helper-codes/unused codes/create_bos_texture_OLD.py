# This program creates a BOS pattern consisting of black or white dots

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

# these are the minimum and maximum X,Y co-ordinates (microns) in the object space for the points on the pattern
X_min = -7.5e4
X_max = 7.5e4
Y_min = -7.5e4
Y_max = 7.5e4

# this is the Z co-ordinate of the pattern (in microns)
Z = 0.0

# this is the diameter of a dot in the bos pattern in microns
dot_diameter = 3.0e3

# this is the spacing for the grid over which the dots will be randomly distributed
grid_spacing = dot_diameter

# these are the number of grid points over X and Y
num_grid_points = (X_max - X_min)/grid_spacing

# this is the ratio of the number of dots to the number of grid points. this ratio should be less than 1
density_dots = 0.9

# this is the number of dots to be generated along each dimension
num_dots = density_dots * num_grid_points

# this sets the DPI (no.of pixels per inch). this is same as the PPI of the monitor (pixels per inch)
dpi = 72

# this sets the background color of the image ('w' for white, 'k' for black)
bg_color = 'w'

# this sets the dot pattern color. this should be the opposite of the bg_color for the dot to be visible.
# ('w' for white, 'k' for black)
dot_color = 'k'

# this calculates the size of the figure in inches
figure_size_inches = ((X_max - X_min)/25.4e3, (Y_max - Y_min)/25.4e3)

# this calculates the size of the figure in pixels
figure_size_pixels = tuple([np.round(dpi*x) for x in figure_size_inches])

# this calculates the diameter of the dot in pixels
dot_diameter_pixels = dot_diameter * 1.0/25.4e3 * dpi
dot_diameter_points = dot_diameter * 1.0/25.4e3 * 72

# this calculates the area of the circle in the scatter plot in pixels
dot_area_pixels = np.pi * (dot_diameter_pixels/2.0)**2.0
dot_area_points = np.pi * (dot_diameter_points/2)**2.0


# this displays the parameters of the texture to the user
print "number of grid points:", num_grid_points
print "number of dots:", num_dots
print "dot diameter in pixels:", dot_diameter_pixels
print "area of the dot in pixels^2:", dot_area_pixels
print "size of the figure in pixels:", figure_size_pixels

# this randomly generates integers corresponding to the X and Y indices in the grid
xy_loc_int = np.random.randint(0, high=num_grid_points, size=(num_dots,2))

# this converts integer locations to positions in real space
xy_loc_float = np.zeros(xy_loc_int.shape, dtype=np.float32)
xy_loc_float[:,0] = X_min + xy_loc_int[:,0]*grid_spacing
xy_loc_float[:,1] = Y_min + xy_loc_int[:,1]*grid_spacing

# this creates a figure where the circles will be drawn

fig = plt.figure(1) #, figsize=(6,6))
fig.set_size_inches(figure_size_inches)
fig.set_dpi(dpi)
ax = plt.Axes(fig, [0., 0., 1., 1.], axisbg=bg_color)
# ax.set_axis_off()
fig.add_axes(ax)

# this draws circles whose centers are the co-ordinates in xy_loc_float and the area is the dot area in pixels
ax.scatter(xy_loc_float[:,0], xy_loc_float[:,1], s=dot_area_points, marker = 'o', c=dot_color, edgecolors=dot_color)
ax.set_xlim([X_min, X_max])
ax.set_ylim([Y_min, Y_max])
# ax.set_autoscale_on(False)
ax.set_xticks([])
ax.set_yticks([])

# plt.axis('equal')
# plt.show()

fig.savefig('bos_scatter_test_dot_black.png', bbox_inches='tight', pad_inches=0, dpi=dpi)