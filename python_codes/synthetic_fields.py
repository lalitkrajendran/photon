import numpy as np
import nrrd

def create_coordinate_grid(n=101, x_range=[-0.5, 0.5], y_range=[-0.5, 0.5]):
    # Function to create a 2D co-ordinate grid
    #
    # INPUTS:
    # n: number of grid points along 1 direction
    # x_range, y_range: min and max values of x, y co-ordinate grid
    #
    # OUTPUTS:
    # X, Y: 2D co-ordinate grid
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)    
    
    # create 1D co-ordinate arrays
    x = np.linspace(start=x_range[0], stop=x_range[1], num=n, endpoint=True)
    y = np.linspace(start=y_range[0], stop=y_range[1], num=n, endpoint=True)

    # create 2D grid
    X, Y = np.meshgrid(x, y, indexing='xy')

    return X, Y


def create_coordinate_grid_3d(n=101, x_range=[-0.5, 0.5], y_range=[-0.5, 0.5], z_range=[-0.5, 0.5]):
    # Function to create a 3D co-ordinate grid
    #
    # INPUTS:
    # n: number of grid points along 1 direction
    # x_range, y_range, z_range: min and max values of x, y, z co-ordinate grid
    #
    # OUTPUTS:
    # X, Y, Z: 3D co-ordinate grid
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)    
    
    # create 1D co-ordinate arrays
    x = np.linspace(start=x_range[0], stop=x_range[1], num=n, endpoint=True)
    y = np.linspace(start=y_range[0], stop=y_range[1], num=n, endpoint=True)
    z = np.linspace(start=z_range[0], stop=z_range[1], num=n, endpoint=True)

    # create grid
    X, Y, Z = np.meshgrid(x, y, z, indexing='xy')

    return X, Y, Z


def create_sine_field(n=101, peak=1.0, l=10, x_range=[-0.5, 0.5], y_range=[-0.5, 0.5]):
    # Function to create a sinusoidal field
    #
    # INPUTS:
    # n: number of grid points
    # peak: max value of the field
    # l: wavelength
    # x_range, y_range: min and max values of x, y co-ordinate grid
    #
    # OUTPUTS:
    # X, Y: 2D co-ordinate grid
    # f: scalar field
    # f_x, f_y: gradient field
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    # create co-ordinate grid
    X, Y = create_coordinate_grid(n, x_range, y_range)

    # calculate wavenumber
    k = 2 * np.pi / l
    
    # create scalar field
    # f = peak * np.sin(k * X) * np.sin(k * Y)
    f = peak * np.cos(k * X) * np.cos(k * Y)
    
    # create gradient field
    # f_x = peak * k * np.cos(k * X) * np.sin(k * Y)
    # f_y = peak * k * np.sin(k * X) * np.cos(k * Y)  
    f_x = - peak * k * np.sin(k * X) * np.cos(k * Y)  
    f_y = - peak * k * np.cos(k * X) * np.sin(k * Y)  

    return X, Y, f, f_x, f_y


def create_sine_field_3d(n=101, peak=1.0, l=10, x_range=[-0.5, 0.5], y_range=[-0.5, 0.5], z_range=[-0.5, 0.5]):
    # Function to create a 3D sinusoidal field
    #
    # INPUTS:
    # n: number of grid points
    # peak: max value of the field
    # l: wavelength
    # x_range, y_range: min and max values of x, y co-ordinate grid
    #
    # OUTPUTS:
    # X, Y: 3D co-ordinate grid
    # f: scalar field
    # f_x, f_y, f_z: gradient field
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    # create co-ordinate grid
    X, Y, Z = create_coordinate_grid_3d(n, x_range, y_range, z_range)

    # calculate origin
    x0, y0, z0 = np.mean(x_range), np.mean(y_range), np.mean(z_range)

    # calculate wavenumber
    k = 2 * np.pi / l
    
    # create scalar field
    # f = peak * np.sin(k * X) * np.sin(k * Y)
    f = peak * np.cos(k * X) * np.cos(k * Y) * np.cos(k * Z)
    
    # create gradient field
    # f_x = peak * k * np.cos(k * X) * np.sin(k * Y)
    # f_y = peak * k * np.sin(k * X) * np.cos(k * Y)  
    f_x = - peak * k * np.sin(k * (X - x0)) * np.cos(k * (Y - y0)) * np.cos(k * (Z - z0))
    f_y = - peak * k * np.cos(k * (X - x0)) * np.sin(k * (Y - y0)) * np.cos(k * (Z - z0))
    f_z = - peak * k * np.cos(k * (X - x0)) * np.cos(k * (Y - y0)) * np.sin(k * (Z - z0))
    
    return X, Y, Z, f, f_x, f_y, f_z


def create_gaussian_field(n=101, peak=1.0, peak_loc=[0.0, 0.0], std=0.1, x_range=[-0.5, 0.5], y_range=[-0.5, 0.5]):
    # Function to create a sinusoidal field
    #
    # INPUTS:
    # n: number of grid points
    # peak: max value of the field
    # peak_loc [1x2 list]: x, y peak location
    # x_range, y_range: min and max values of x, y co-ordinate grid
    #
    # OUTPUTS:
    # X, Y: 2D co-ordinate grid
    # f: scalar field
    # f_x, f_y: gradient field
    #
    # AUTHOR:
    # Lalit Rajendran (lrajendr@purdue.edu)
    
    # create co-ordinate grid
    X, Y = create_coordinate_grid(n, x_range, y_range)

    # create scalar field
    f = peak * np.exp( - 1.0/(2 * std**2) * ( (X - peak_loc[0])**2 + (Y - peak_loc[1])**2))
    
    # create gradient field
    f_x = peak * - 1.0/(2 * std**2) * 2 * (X - peak_loc[0]) * np.exp( - 1.0/(2 * std**2) * ((X - peak_loc[0])**2 + (Y - peak_loc[1])**2))
    f_y = peak * - 1.0/(2 * std**2) * 2 * (Y - peak_loc[1]) * np.exp( - 1.0/(2 * std**2) * ((X - peak_loc[0])**2 + (Y - peak_loc[1])**2))

    return X, Y, f, f_x, f_y


def save_nrrd(data, nrrd_filename):
    # This function saves the density stored in data['rho'] to an nrrd file
    # specified by nrrd_filename

    # specify orientation of space
    space = '3D-right-handed'

    # generate arrays to store co-ordinates along the three axes
    x = np.array(data['x'])
    y = np.array(data['y'])
    z = np.array(data['z'])

    # set origin
    x0 = x.min()
    y0 = y.min()
    z0 = z.min()

    space_orig = np.array([x0, y0, z0]).astype('float32')

    # set grid spacing
    del_x = np.diff(x)[0]
    del_y = np.diff(y)[0]
    del_z = np.diff(z)[0]

    spacings = np.array([del_x, del_y, del_z]).astype('float32')

    # spcify other relevant options
    options = {'type' : 'f4', 'space': space, 'encoding': 'raw',
               'space origin' : space_orig, 'spacings' : spacings}

    print("saving density to %s \n" % nrrd_filename)

    # save data to nrrd file
    nrrd.write(nrrd_filename, np.array(data['rho']).astype('float32'), options, compression_level=1)


def calculate_theoretical_deflection(rho_grad, M, Z_D, del_z, rho_0, pixel_pitch): #data, grad_x, z_1, z_2, M, pixel_pitch):
    # this function calculates the theoretical displacement of a dot pattern on the camera sensor for a constant density
    # gradient field using the paraxial approximation

    # this is the gladstone-dale constant (m^3/kg)
    K = 0.225e-3

    # this is the refractive index of the undisturbed medium
    n_0 = K * rho_0 + 1

    # this is the gradient of refractive index (m^-1)
    n_grad = K * rho_grad

    # angular deflection of light ray (radians)
    epsilon = 1/n_0 * n_grad * del_z

    # this is the displacement on the sensor (pixels)
    displacement = M * Z_D * epsilon / pixel_pitch

    print('density gradient (kg/m^4): %.2G' % rho_grad)
    print('refractive index gradient (/m): %.4G' % n_grad)
    print('angular deflection of ray (radians): %.6G' % epsilon)
    print('displacement on the sensor (pixels): %.2G' % displacement)


def calculate_density_gradient(disp, M, Z_D, del_z, rho_0, pixel_pitch): #delta_x, delta_y, data, rho_0, M, Z_D, pixel_pitch):
    '''
    this function calculates the density gradients required to produce the specified
    pixel displacement at the camera sensor given the extent of the density gradient
    volume, magnification and pixel pitch

    '''

    # this is the gladstone-dale constant (m^3/kg)
    K = 0.225e-3
    
    # this is the refractive index of the undisturbed medium
    n_0 = K * rho_0 + 1
    
    # this is the required angular deflection of the ray (radians)
    epsilon = disp * pixel_pitch/(Z_D * M)
    
    # this is the required refractive index gradient
    n_grad = epsilon * n_0/del_z
        
    # this is the required gradient of density along x (kg/m^4)
    rho_grad = 1/K * n_grad

    return rho_grad


def calculate_density_noise(displacement_noise_std, M, Z_D, delta_x, delta_z, rho_0, pixel_pitch):
    '''
    This function calculates the standard deviation of the noise to be added to the
    density field to obtain a specified standard deviation of the noise in the displacement
    field.
    
    This is calculated based on the propogation of errors in the BOS image generation methodology.
    
    It is assumed that the noise is a random variable drawn from a zero-mean Gaussian 
    distribution.
    
    INPUT:
    displacement_noise_std: standard deviation of required noise in the displacement field | scalar | float | pix.
    M: magnification of the optical system | scalar | float | unitless
    Z_D: distance between the target and the mid-point of the density gradient field | scalar | float | m
    delta_x: grid spacing along x for points in the density gradient field | scalar | float | m
    delta_z: extent of the density gradient field along z | scalar | float | m
    rho_0: ambient density | scalar | float | kg/m^3
    pixel_pitch: dimension of a single pixel on the camera sensor | scalar | float | m
    
    OUTPUT:
    rho_noise_std: standard deviation of noise to be added to the density field | scalar | float | kg/m^3
    '''
    
    # gladstone dale constant (m^3/kg)
    K = 0.225e-3
    
    # calculate ambient refractive index
    n_0 = K * rho_0 + 1
    
    # calculate standard deviation of noise to be added to the density field (added factor of 2 to account for central differencing in calculation of the refractive index gradient)
    rho_noise_std = 2 * displacement_noise_std * pixel_pitch * delta_x / (M * Z_D * K/n_0 * np.sqrt(2.0) * delta_z)
    
    return rho_noise_std