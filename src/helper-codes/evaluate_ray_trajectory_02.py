'''
The purpose of this code is to read teh final positions and directions of the light rays from file and compare them to
theoretical value to compute the error
'''

import numpy as np
import matplotlib.pyplot as plt

######### calculate trajectory from theory #########

# this is the density of the undisturbed medium (kg/m^3)
rho_0 = 1.225

# this is the gladstone-dale constant (m^3/kg)
K = 0.225e-3

# this is the refractive index of the undisturbed medium
n_0 = K * rho_0 + 1

# gradient of the square of the refractive index (m^-1)
alpha = 1e-1

# this is the angle between the initial ray direction and the y axis (radians)
# theta_0 = np.pi/2 - np.radians(80.0)
theta_0 = np.radians(45.0)
print "theta (rad)", theta_0

# this is a temporary variable that depends on n0 and theta0
xi_0 = n_0 * np.cos(theta_0)

################ generate an array of z co-ordinates ###################
nz = 100

z = np.linspace(start=0, stop=10, num = nz, endpoint=True)

################## calculate light ray trajectory #####################
r = 2 * xi_0/alpha * (np.sqrt(-xi_0**2 + n_0**2 + alpha * z) - np.sqrt(-xi_0**2 + n_0**2))

print "curvature", abs(np.diff(np.diff(r))).max()
# print "curvature", np.diff(np.diff(r))
################## plot results #######################################
# title_string = r'$\alpha = %d, \theta_0 = %.0f^\circ$' % (alpha, np.degrees(theta_0))
# figure_save_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/results/validation/plots/'

# plt.plot(z,r)
# plt.xlabel('z')
# plt.ylabel('r')
# plt.title(title_string)
# plt.savefig(figure_save_filepath + 'trajectory.png')
# plt.show()

#################### load light ray positions from file #####################

# light ray filepath
lightray_filepath = '/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/'

nz_array = np.array([10, 100, 1000])


for nz_index in range(0, nz_array.size):

    nz = nz_array[nz_index]

    filename_pos = 'lightrayPos_f_%d.bin' % nz
    filename_dir = 'lightrayDir_f_%d.bin' % nz

    # load files
    pos_f = np.fromfile(lightray_filepath + filename_pos, dtype=np.float32)
    dir_f = np.fromfile(lightray_filepath + filename_dir, dtype=np.float32)

    # reshape arrays
    pos_f = np.reshape(pos_f, (pos_f.size/3, 3))
    dir_f = np.reshape(dir_f, (dir_f.size/3, 3))

    if(nz_index == 0):
        theta_0 = np.radians(np.linspace(start=55.0, stop=10.0, num=pos_f.shape[0], endpoint=True)).T

        # z_code = np.zeros((theta_0.size, nz_array.size))
        r_code = np.zeros((theta_0.size, nz_array.size))
        r_theory = np.zeros(r_code.shape)
        error_r = np.zeros(r_code.shape)

    ############### compare theoretical positions to actual positions and calculate error ################

    xi_0 = n_0 * np.cos(theta_0).T

    z_code = pos_f[:,2]
    r_code[:,nz_index] = pos_f[:,1]
    r_theory[:,nz_index] = 2 * xi_0/alpha * (np.sqrt(-xi_0**2 + n_0**2 + alpha * z_code) - np.sqrt(-xi_0**2 + n_0**2))

    error_r[:,nz_index] = pos_f[:,1] - r_theory[:,nz_index]

# plot results
figure_save_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/results/validation/plots/'

plt.figure()
plt.hold(True)
plt.plot(np.degrees(theta_0), r_theory[:,2], 'k', label='theory')
c = ['r', 'g', 'b']
for nz_index in range(0, nz_array.size):
    nz = nz_array[nz_index]

    plt.plot(np.degrees(theta_0), r_code[:,nz_index], c[nz_index]+'*', label='nz=%d' % nz)

plt.xlabel(r'$\theta \;\;(deg.)$')
plt.ylabel(r'$y$')
plt.legend()
plt.grid(True)
plt.savefig(figure_save_filepath + 'n2linear-effect-of-nz.png')
# plt.show()

