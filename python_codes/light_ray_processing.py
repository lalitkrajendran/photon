import numpy as np
import scipy.io as sio
import sys
import os.path
import matplotlib.pyplot as plt
import platform
from scipy.interpolate import griddata

if platform.system() == 'Linux':
    mount_directory = '/scratch/shannon/c/aether/'
else:
    mount_directory = '/Volumes/aether_c/'

# import path containing python modules
sys.path.append(os.path.join(mount_directory, 'Projects/BOS/general-codes/python-codes'))
import modify_plt_settings
import loadmat_functions

# modify plot settings
import matplotlib
matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt = modify_plt_settings.modify_plt_settings(plt)

# get color cycle
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
# convert to rgba array
colors = matplotlib.colors.to_rgba_array(colors)


def load_image_from_bin(filename, shape):
    # Function to load raw image from ray tracing
    #
    # INPUTS:
    # filename (str): name of the .bin file
    # shape (tuple): row, col size of the image
    #
    # OUTPUTS:
    # I (float): 2D image array
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    I_temp = np.fromfile(filename, dtype='float32')

    I = np.reshape(I_temp, newshape=shape)

    return I


def load_intermediate_light_ray_positions(filename, num_intermediate_positions):    
    # Function to load final light ray positions
    #
    # INPUTS:
    # filename: path to .bin file containing light ray positions
    # num_intermediate_positions: number of positions saved for each ray
    #
    # OUTPUTS:
    # pos_dict (dict): ray positions (um) along x,y,z (image co-ord) and number of rays
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    pos_temp = np.fromfile(file=filename, dtype='float32')

    num_rays = int(pos_temp.size/(3 * num_intermediate_positions))

    pos = np.reshape(a=pos_temp, newshape=(num_rays, num_intermediate_positions, 3))

    pos_dict = {'x': pos[:, :, 0], 'y': pos[:, :, 1], 'z': pos[:, :, 2], 'num_rays': num_rays}

    return pos_dict


def load_light_ray_positions(filename):    
    # Function to load final light ray positions
    #
    # INPUTS:
    # filename: path to .bin file containing light ray positions
    #
    # OUTPUTS:
    # pos_dict (dict): ray positions (um) along x,y,z (image co-ord) and number of rays
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    pos_temp = np.fromfile(file=filename, dtype='float32')

    num_rays = int(pos_temp.size/3)

    pos = np.reshape(a=pos_temp, newshape=(num_rays, 3))

    pos_dict = {'x': pos[:, 0], 'y': pos[:, 1], 'z': pos[:, 2], 'num_rays': num_rays}

    return pos_dict


def load_intermediate_light_ray_directions(filename, num_intermediate_positions):    
    # Function to load light ray directions after density gradients
    #
    # INPUTS:
    # filename: path to .bin file containing light ray directions
    #
    # OUTPUTS:
    # dir_dict (dict): ray angles (rad.) along x,y,z (image co-ord) and number of rays
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)

    dir_temp = np.fromfile(file=filename, dtype='float32')

    num_rays = int(dir_temp.size/(3 * num_intermediate_positions))

    dirs = np.reshape(a=dir_temp, newshape=(num_rays, num_intermediate_positions, 3))

    dir_dict = {'x': np.arccos(dirs[:, :, 0]), 'y': np.arccos(dirs[:, :, 1]), 'z': np.arccos(dirs[:, :, 2]), 'num_rays': num_rays}

    return dir_dict


def load_light_ray_directions(filename):    
    # Function to load light ray directions after density gradients
    #
    # INPUTS:
    # filename: path to .bin file containing light ray directions
    #
    # OUTPUTS:
    # dir_dict (dict): ray angles (rad.) along x,y,z (image co-ord) and number of rays
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)

    dir_temp = np.fromfile(file=filename, dtype='float32')

    num_rays = int(dir_temp.size/3)

    dirs = np.reshape(a=dir_temp, newshape=(num_rays, 3))

    dir_dict = {'x': np.arccos(dirs[:, 0]), 'y': np.arccos(dirs[:, 1]), 'z': np.arccos(dirs[:, 2]), 'num_rays': num_rays}

    return dir_dict


def load_light_ray_data(folder):
    # Function to light ray positions and directions
    #
    # INPUTS:
    # folder: directory containing light ray data for reference and gradient images
    #
    # OUTPUTS:
    # pos1, pos2: ray positions for reference and gradient images
    # dir1, dir2: ray directions for reference and gradient images
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    # -----------------
    # load positions
    # -----------------    
    # reference image
    filename = os.path.join(folder, 'light-ray-positions', 'im1', 'pos_0000.bin')
    pos1 = load_light_ray_positions(filename)

    # gradient image
    filename = os.path.join(folder, 'light-ray-positions', 'im2', 'pos_0000.bin')
    pos2 = load_light_ray_positions(filename)

    # -----------------
    # load directions
    # -----------------    
    # reference image
    filename = os.path.join(folder, 'light-ray-directions', 'im1', 'dir_0000.bin')
    dir1 = load_light_ray_directions(filename)

    # gradient image
    filename = os.path.join(folder, 'light-ray-directions', 'im2', 'dir_0000.bin')
    dir2 = load_light_ray_directions(filename)

    return pos1, pos2, dir1, dir2


def load_intermediate_light_ray_data(folder, num_intermediate_positions):
    # Function to light ray positions and directions
    #
    # INPUTS:
    # folder: directory containing light ray data for reference and gradient images
    # num_intermediate_positions: number of intermediate positions saved inside the density gradients
    #
    # OUTPUTS:
    # ipos, idir: intermediate ray positions and directions
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    # -----------------
    # load positions
    # -----------------    
    # gradient image
    filename = os.path.join(folder, 'light-ray-positions', 'im2', 'intermediate_pos_0000.bin')
    ipos = load_intermediate_light_ray_positions(filename, num_intermediate_positions)

    # -----------------
    # load directions
    # -----------------    
    # gradient image
    filename = os.path.join(folder, 'light-ray-directions', 'im2', 'intermediate_dir_0000.bin')
    idir = load_intermediate_light_ray_directions(filename, num_intermediate_positions)

    return ipos, idir


def calculate_lightray_deflections(pos1, pos2, dir1, dir2):
    # Function to calculate lightray deflections
    #
    # INPUTS:
    # pos1, pos2: ray positions for reference and gradient images
    # dir1, dir2: ray directions for reference and gradient images
    #
    # OUTPUTS:
    # d_pos: ray displacements
    # d_dir: ray angular deflections
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
        
    # co-ordinates
    coords = pos1.keys()

    # initialize dictionaries to hold results
    d_pos = dict()
    d_dir = dict()
    
    # loop through co-ordinates and calculate displacements and deflections
    for i, coord in enumerate(coords):
        # displacements
        # d_pos[coord] = pos2[coord] - pos1[coord]
        d_pos[coord] = pos1[coord] - pos2[coord]
        # angular deflections
        d_dir[coord] = dir2[coord] - dir1[coord]

    return d_pos, d_dir


def calculate_dot_average(a, num_rays_per_dot):
    # Function to calculate dot averaged properties for the light rays
    #
    # INPUTS:
    # a: any property for which the average is to be calculated
    # num_rays_per_dot: number of light rays launched from each dot
    #
    # OUTPUTS:
    # a_dot: dot averaged property
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)

    # ensure num rays per dot is int
    num_rays_per_dot = int(num_rays_per_dot)        
    
    # co-ordinates
    coords = list(a.keys())

    # number of rays
    num_rays = a[coords[0]].size

    # intialize dictionary
    a_dot = dict()
    
    # calculate average
    for i, coord in enumerate(coords):
        if coord == 'num_rays':
            continue
        a_dot[coord] = np.add.reduceat(a[coord], range(0, num_rays, num_rays_per_dot)) / num_rays_per_dot

    return a_dot


def convert_pos_to_pix(pos, pos0, pixel_pitch):
    # Function to convert physical positions (um) to pixels
    #
    # INPUTS:
    # pos (dict): ray positions (um)
    # pos0 (dict): sensor origin (um)
    # pixel_pitch (float): size of a pixel (um)
    #
    # OUTPUTS:
    # pos_i (dict): positions in pixels
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    # co-ordinates
    coords = ['x', 'y']

    # initialize dictionary
    pos_i = dict()

    # convert to pixels
    for coord in coords:
        pos_i[coord] = (pos[coord] - pos0[coord]) / pixel_pitch

    return pos_i


def plot_dot_displacements_quiver(pos, d_pos, x_lim, y_lim, skip=1, scale=None):
    # Function to plot dot displacements
    #
    # INPUTS:
    # pos: ray positions on reference image
    # d_pos: ray displacements
    # x_lim: axis limits along x
    # y_lim: axis limits along y
    # skip: number of vectors to skip
    # scale: scaling factor for quiver plot (smaller scale to get larger arrows)
    #
    # OUTPUTS:
    # fig, ax: figure and axes for the plotted figure
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
        
    fig = plt.figure()
    plt.quiver(pos['x'][::skip], pos['y'][::skip], d_pos['x'][::skip], d_pos['y'][::skip], scale=scale,
                    color=colors[0, :], headwidth=5)
    ax = fig.axes[0]
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.set_xlabel('X (pix.)')
    ax.set_ylabel('Y (pix.)')
    
    return fig, ax


def plot_dot_displacements_contour(pos, d_pos, x_lim, y_lim, grid_spacing, skip=1):
    # Function to plot dot displacements as contour
    #
    # INPUTS:
    # pos: ray positions on reference image
    # d_pos: ray displacements
    # x_lim: axis limits along x
    # y_lim: axis limits along y
    # grid_spacing: grid spacing to interpolate deflections
    # skip: number of grid points to skip in the contour plot
    #
    # OUTPUTS:
    # fig, ax: figure and axes for the plotted figure
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)

    # create 2D co-ordinate grid
    X, Y = create_coordinate_grid(x_lim, y_lim, grid_spacing)
    
    # interpolate displacements onto grid
    pos_grid, d_pos_grid = interpolate_deflections_to_grid(pos, d_pos, X, Y)

    # calculate displacement magnitude
    d_pos_mag = np.sqrt(d_pos_grid['x']**2 + d_pos_grid['y']**2)
    
    # plot contours 
    # fig = plt.figure()
    plt.pcolormesh(pos_grid['x'][::skip], pos_grid['y'][::skip], d_pos_mag[::skip], cmap='plasma')
    plt.colorbar()
    # ax = fig.axes[0]
    ax = plt.gca()
    ax.set_aspect('equal')#, adjustable='box')
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.set_xlabel('X (pix.)')
    ax.set_ylabel('Y (pix.)')

    # return fig, ax

def create_coordinate_grid(x_lim, y_lim, grid_spacing): 
    # Function to create a 2D co-ordinate grid
    #
    # INPUTS:
    # x_lim, y_lim: axis limits
    # grid_spacing: spacing of the grid to interpolate the displacements    
    #
    # OUTPUTS:
    # X, Y: x,y co-ordinates on the grid (2D arrays)
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
        
    # calculate number of grid points
    num_grid_points_x = int((x_lim[1] - x_lim[0]) / grid_spacing)
    num_grid_points_y = int((y_lim[1] - y_lim[0]) / grid_spacing)
    
    # create 2D co-ordinate grid
    x = x_lim[0] + np.linspace(start=0, stop=num_grid_points_x-1, num=num_grid_points_x, endpoint=True) * grid_spacing
    y = y_lim[0] + np.linspace(start=0, stop=num_grid_points_y-1, num=num_grid_points_y, endpoint=True) * grid_spacing
    X, Y = np.meshgrid(x, y)

    return X, Y

def interpolate_deflections_to_grid(pos, d_pos, X, Y):
    # Function to interpolate dot deflections to a grid
    #
    # INPUTS:
    # pos: dot positions
    # d_pos: dot displacements
    # X, Y: co-ordinate grid on which interpolation is to be performed
    #
    # OUTPUTS:
    # pos_grid: grid co-ordinates
    # d_pos_grid: interpolated displacements
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    # collect points with non-nan indices
    valid_idx_pos = np.logical_and(np.isfinite(pos['x']), np.isfinite(pos['y']))
    valid_idx_dpos = np.logical_and(np.isfinite(d_pos['x']), np.isfinite(d_pos['y']))
    valid_idx = np.logical_and(valid_idx_pos, valid_idx_dpos)
    x_valid = pos['x'][valid_idx]
    y_valid = pos['y'][valid_idx]
    dx_valid = d_pos['x'][valid_idx]
    dy_valid = d_pos['y'][valid_idx]

    # create array of points
    points = np.transpose(np.vstack((x_valid, y_valid)))

    # interpolate onto grid
    dx_grid = griddata(points, dx_valid, (X, Y), method='linear', fill_value=0)
    dy_grid = griddata(points, dy_valid, (X, Y), method='linear', fill_value=0)

    # dx_grid = griddata(points, dx_valid, (X, Y), method='cubic', fill_value=0)
    # dy_grid = griddata(points, dy_valid, (X, Y), method='cubic', fill_value=0)

    # prepare variables to return
    pos_grid = dict()
    pos_grid['x'] = X
    pos_grid['y'] = Y

    d_pos_grid = dict()
    d_pos_grid['x'] = dx_grid
    d_pos_grid['y'] = dy_grid

    return pos_grid, d_pos_grid


def remove_edge_dots(pos, dir, buffer, num_pixels):
    # Function to remove ray information for dots on the boundary
    # 
    # INPUTS:
    # pos: ray positions
    # buffer: margin buffer (pix.)
    # num_pixels: number of pixels on the camera sensor
    #
    # OUTPUTS:
    # pos2: ray positions after removing edge dots
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
        
    coords = pos.keys()
    ids = dict()

    # identify indices lying outside the buffer
    for coord in coords:
        # ids[coord] = np.argwhere(np.logical_or(pos[coord] < buffer, pos[coord] > (num_pixels-buffer)))
        ids[coord] = np.argwhere(np.logical_or(pos[coord] < buffer, pos[coord] > (num_pixels-buffer)))
            
    # find all indices            
    if ids['x'].size == 0:
        ids_del = ids['y']
    elif ids['y'].size == 0:
        ids_del = ids['x']
    else:
        ids_del = np.unique(np.concatenate((ids['x'], ids['y'])))
        # ids_del = np.argwhere(np.logical_or(ids['x'], ids['y']))
        
    # delete indices
    pos2 = pos
    dir2 = dir
    for coord in coords:
        pos2[coord][ids_del] = np.nan
        dir2[coord][ids_del] = np.nan

    return pos2, dir2


def load_and_display_image(folder):
    # Function to load and display the dot image
    #
    # INPUTS:
    # folder: folder containing the tif images
    #
    # OUTPUTS:
    # fig, ax: figure and axes from plot
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    # load image
    im_ref = plt.imread(os.path.join(folder, 'bos_pattern_image_1.tif'))
    im_grad = plt.imread(os.path.join(folder, 'bos_pattern_image_2.tif'))

    # display image
    fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(8, 8))
    plt.axes(ax[0])
    plt.imshow(im_ref, cmap='Greys_r')
    plt.title('Ref')
    plt.axes(ax[1])
    plt.imshow(im_grad, cmap='Greys_r')
    plt.title('Grad')

    return fig, ax


def calculate_sensor_origin(num_pixels, pixel_pitch):
    # Function to calculate physical location of the sensor origin
    #
    # INPUTS:
    # num_pixels: number of pixels on the camera sensor
    # pixel_pitch: size of a sensor
    #
    # OUTPUTS:
    # pos0: physical co-ordinates of the origin pixel
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    # calculate physical location of 0,0 pixel 
    pos0 = {'x': -(num_pixels/2 - 1) * pixel_pitch, 'y': -(num_pixels/2 - 1) * pixel_pitch, 'z': 0.0}
    
    return pos0

    
def process_lightray_data(folder, save_results=False, display_progress=False):
    # Function to load and process lightray data
    #
    # INPUTS:
    # folder: directory containing light ray data
    # save_results: save dot position and deflections to a mat file in the same folder
    # display_progress: display intermediate progress to user
    #
    # OUTPUTS:
    # pos_dot: dot averaged ray positions 
    # dir_dot: dot averaged ray directions 
    # d_pos_dot: dot averaged ray displacements
    # d_dir_dot: dot averaged ray deflections
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
            
    # # load and display image
    # fig, ax = load_and_display_image(folder)

    # --------------------------
    # load image generation parameters
    # --------------------------
    if display_progress:
        print('Loading Parameters')
    parameters = loadmat_functions.loadmat(os.path.join(folder, 'parameters.mat'))

    # --------------------------
    # calculate camera sensor origin
    # --------------------------
    num_pixels = parameters['camera_design']['x_pixel_number']
    pixel_pitch = parameters['camera_design']['pixel_pitch']
    pos0 = calculate_sensor_origin(num_pixels, pixel_pitch)

    # --------------------------
    # load light ray position and directions
    # --------------------------
    if display_progress:
        print('Loading Light Ray Data')    
    pos1, pos2, dir1, dir2 = load_light_ray_data(folder)

    # --------------------------
    # convert light ray positions to pixels
    # --------------------------
    pos1 = convert_pos_to_pix(pos1, pos0, pixel_pitch)
    pos2 = convert_pos_to_pix(pos2, pos0, pixel_pitch)

    # --------------------------
    # calculate dot averages
    # --------------------------
    if display_progress:
        print('Calculating Dot Averages')
    num_rays_per_dot = int(parameters['bos_pattern']['lightray_number_per_particle'])
    pos1_dot = calculate_dot_average(a=pos1, num_rays_per_dot=num_rays_per_dot)
    pos2_dot = calculate_dot_average(a=pos2, num_rays_per_dot=num_rays_per_dot)

    dir1_dot = calculate_dot_average(a=dir1, num_rays_per_dot=num_rays_per_dot)
    dir2_dot = calculate_dot_average(a=dir2, num_rays_per_dot=num_rays_per_dot)

    # --------------------------
    # calculate light ray deflections
    # --------------------------
    print('Calculating Deflections')
    d_pos_dot, d_dir_dot = calculate_lightray_deflections(pos1_dot, pos2_dot, dir1_dot, dir2_dot)

    # --------------------------
    # remove edge dots
    # --------------------------
    dot_diameter = parameters['camera_design']['diffraction_diameter']
    # pos1_dot, dir1_dot = remove_edge_dots(pos1_dot, dir1_dot, buffer=2*dot_diameter, num_pixels=num_pixels)
    # pos2_dot, dir2_dot = remove_edge_dots(pos2_dot, dir2_dot, buffer=2*dot_diameter, num_pixels=num_pixels)

    # # --------------------------
    # # remove very large displacements
    # # --------------------------
    # print('Removing Large Deflections')
    # indices = np.argwhere(np.logical_or(abs(d_pos_dot['x']) > 5, abs(d_pos_dot['y']) > 5))
    # # d_pos_dot['x'][indices] = np.nan
    # # d_pos_dot['y'][indices] = np.nan
    # # d_dir_dot['x'][indices] = np.nan
    # # d_dir_dot['y'][indices] = np.nan
    # d_pos_dot['x'] = np.delete(d_pos_dot['x'], indices)
    # d_pos_dot['y'] = np.delete(d_pos_dot['y'], indices)
    # d_dir_dot['x'] = np.delete(d_dir_dot['x'], indices)
    # d_dir_dot['y'] = np.delete(d_dir_dot['y'], indices)

    # display min, max displacements
    print('x: %.2f to %.2f pix.' % (np.nanmin(d_pos_dot['x']), np.nanmax(d_pos_dot['x'])))
    print('y: %.2f to %.2f pix.' % (np.nanmin(d_pos_dot['y']), np.nanmax(d_pos_dot['y'])))

    # display displacements
    # fig, ax = plot_dot_displacements_quiver(pos1_dot, d_pos_dot, [1, num_pixels], [1, num_pixels], skip=2, scale=0.5e2)

    # --------------------------
    # save results to file
    # --------------------------
    # prepare results for saving
    pos_dot = [pos1_dot, pos2_dot]
    dir_dot = [dir1_dot, dir2_dot]
    
    if save_results:
        # mat_file = {'x': pos1_dot['x'], 'y': pos1_dot['y'], 'dx': d_pos_dot['x'], 'dy': d_pos_dot['y']}
        # mat_file = {'pos': pos_dot, 'dir': dir_dot, 'd_pos': d_pos_dot, 'd_dir': d_dir_dot}
        mat_file = {'pos1': pos1_dot, 'pos2': pos2_dot, 'dir1': dir1_dot, 'dir2': dir2_dot, 'd_pos': d_pos_dot, 'd_dir': d_dir_dot}
        sio.savemat(os.path.join(folder, 'dot-deflections.mat'), mat_file, long_field_names=True)

    return pos_dot, dir_dot, d_pos_dot, d_dir_dot

