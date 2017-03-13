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
theta_0 = np.radians(55.0)
print "theta (rad)", theta_0

# this is a temporary variable that depends on n0 and theta0
xi_0 = n_0 * np.cos(theta_0)

################ generate an array of z co-ordinates ###################
nz = 100

z = np.linspace(start=0, stop=5.25, num = nz, endpoint=True)
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

# position file name
filename_pos = 'lightrayPos_f_10.bin'
filename_dir = 'lightrayDir_f_10.bin'

# load files
pos_f = np.fromfile(lightray_filepath + filename_pos, dtype=np.float32)
dir_f = np.fromfile(lightray_filepath + filename_dir, dtype=np.float32)

# reshape arrays
pos_f = np.reshape(pos_f, (pos_f.size/3, 3))
dir_f = np.reshape(dir_f, (dir_f.size/3, 3))

################ plot the trajectory and the final light ray positions #####
title_string = r'$\alpha = %.1f, \theta_0 = %.0f^\circ$' % (alpha, np.degrees(theta_0))
figure_save_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/results/validation/plots/'

plt.plot(z,r, label='theory')
plt.hold(True)
plt.plot(pos_f[0,2], pos_f[0,1], 'r*', label='code')
plt.xlim([0, 1.1*pos_f[0,2]])
plt.legend(bbox_to_anchor=(0.25,1))
plt.grid(True)
plt.xlabel('z')
plt.ylabel('r')
plt.title(title_string)

plt.savefig(figure_save_filepath + 'trajectory.png')
plt.show()

# ############### compare theoretical positions to actual positions and calculate error ################
#
# theta_0 = np.radians(np.linspace(start=55.0, stop=10.0, num=pos_f.shape[0], endpoint=True)).T
# xi_0 = n_0 * np.cos(theta_0).T
#
# z_code = pos_f[:,2]
# r_theory = 2 * xi_0/alpha * (np.sqrt(-xi_0**2 + n_0**2 + alpha * z_code) - np.sqrt(-xi_0**2 + n_0**2))
#
# error_r = pos_f[:,1] - r_theory
#
# # plot results
# figure_save_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/results/validation/plots/'
#
# plt.figure()
# plt.plot(np.degrees(theta_0), pos_f[:,1], 'b*', label='code')
# plt.plot(np.degrees(theta_0), r_theory, 'r', label='theory')
# plt.plot(np.degrees(theta_0), abs(error_r), 'g*', label='error')
# plt.xlabel(r'$\theta \;\;(deg.)$')
# plt.ylabel(r'$y$')
# plt.legend()
# plt.grid(True)
# plt.savefig(figure_save_filepath + 'n2linear-error.png')
# # plt.show()

