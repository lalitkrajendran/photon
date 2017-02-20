import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, )
        Radius of circle in data unit.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)
        `c` can be a 2-D array in which the rows are RGB or RGBA, however.
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls),
        norm, cmap, transform, etc.

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
    plt.colorbar()

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """


    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None
    if 'fc' in kwargs: kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs: kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs: kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs: kwargs.setdefault('linewidth', kwargs.pop('lw'))

    patches = [Circle((x_, y_), s_) for x_, y_, s_ in np.broadcast(x, y, s)]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        collection.set_array(np.asarray(c))
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    if c is not None:
        plt.sci(collection)
    return collection

# these are the minimum and maximum X,Y co-ordinates (microns) in the object space for the points on the pattern
X_min = -7.5e4
X_max = 7.5e4
Y_min = -7.5e4
Y_max = 7.5e4

# this is the Z co-ordinate of the pattern (in microns)
Z = 0.0

# this is the diameter of a dot in the bos pattern in microns
dot_diameter = 2e2

# this is the spacing for the grid over which the dots will be randomly distributed
grid_spacing = 2 * dot_diameter

# these are the number of grid points over X and Y
num_grid_points = np.round((X_max - X_min)/grid_spacing)

# this is the ratio of the number of dots to the number of grid points. this ratio should be less than 1
density_dots = 0.9

# this is the number of dots to be generated along each dimension
num_dots = np.round(density_dots * num_grid_points)

# this sets the DPI (no.of pixels per inch). this is an integer multiple of PPI of the monitor (pixels per inch) if DPI is in
# multiples of 72.
dpi = 936

# this scales the total size of the figure
scale_factor = 1

# this sets the background color of the image ('w' for white, 'k' for black)
bg_color = 'w'

# this sets the dot pattern color. this should be the opposite of the bg_color for the dot to be visible.
# ('w' for white, 'k' for black)
dot_color = 'k' if bg_color=='w' else 'w'

# this calculates the size of the figure in inches
figure_size_inches = (scale_factor*(X_max - X_min)/25.4e3, scale_factor*(Y_max - Y_min)/25.4e3)

# this calculates the size of the figure in pixels
figure_size_pixels = tuple([np.round(dpi*x) for x in figure_size_inches])

# this calculates the diameter of the dot in pixels
dot_diameter_image = dot_diameter * 1.0/25.4e3 * dpi

# this calculates the area of the circle in the scatter plot in pixels
dot_area_image = np.pi * (dot_diameter_image/2)**2.0

# this displays the parameters of the texture to the user
print "number of grid points:", num_grid_points
print "number of dots:", num_dots
print "dot diameter in pixels:", dot_diameter_image
print "area of the dot in pixels^2:", dot_area_image
print "size of the figure in pixels:", figure_size_pixels


fig = plt.figure(1)
fig.set_size_inches(figure_size_inches)
# ax = plt.subplot(aspect='equal')
ax = plt.Axes(fig, [0., 0., 1., 1.], axisbg=bg_color)
# ax.set_axis_off()
fig.add_axes(ax)

# this randomly generates integers corresponding to the X and Y indices in the grid
xy_loc_int = np.random.randint(0, high=num_grid_points, size=(num_dots,2))

# this converts integer locations to positions in real space
xy_loc_float = np.zeros(xy_loc_int.shape, dtype=np.float32)
xy_loc_float[:,0] = X_min + xy_loc_int[:,0]*grid_spacing
xy_loc_float[:,1] = Y_min + xy_loc_int[:,1]*grid_spacing

out = circles(xy_loc_float[:,0], xy_loc_float[:,1], 0.5*dot_diameter*np.ones((num_dots,1)), c=dot_color, alpha=0.5, ec='none')

ax.set_xlim([X_min, X_max])
ax.set_ylim([Y_min, Y_max])
ax.set_xticks([])
ax.set_yticks([])

fig.savefig('bos_pattern_black_dot_936.png', bbox_inches='tight', pad_inches=0, dpi = dpi)

# NOTE : The plot will not have the axis equal, but that is ok. the saved image has the correct resolution along both
# axes
# plt.show()