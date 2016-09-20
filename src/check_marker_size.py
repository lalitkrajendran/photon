'''
this program is used to study the relation between the marker size specified in the plot command and the actual size
in the image
'''


import numpy as np
import matplotlib.pyplot as plt

# define arrays
xy_loc_float = np.zeros((5,2))
xy_loc_float[:,0] = np.arange(start=0, stop=5, step=1, dtype=np.float32)
xy_loc_float[:,1] = np.arange(start=0, stop=5, step=1, dtype=np.float32)

# figure size in inches
figure_size_inches = 6

# dpi of the monitor
dpi = 96

# make scatter plot
fig = plt.figure(1) #, figsize=(6,6))
fig.set_size_inches(figure_size_inches)
fig.set_dpi(dpi)
ax = plt.Axes(fig, [0., 0., 1., 1.], axisbg='black')
# ax.set_axis_off()
fig.add_axes(ax)

ax.scatter(xy_loc_float[:,0], xy_loc_float[:,1], s=dot_area_pixels, c='w', edgecolors='none')
ax.set_xlim([X_min, X_max])
ax.set_ylim([Y_min, Y_max])
# ax.set_autoscale_on(False)
ax.set_xticks([])
ax.set_yticks([])

plt.axis('equal')
# plt.show()

fig.savefig('bos_scatter_test_circles.png', bbox_inches='tight', pad_inches=0, dpi=dpi)