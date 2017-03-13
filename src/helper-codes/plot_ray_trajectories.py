'''
The purpose of this program is to plot the analytical solution to a light ray trajectory in a mdium with an n^2 linear
refractive index variation
'''

import numpy as np
import matplotlib.pyplot as plt

# this is the density of the undisturbed medium (kg/m^3)
rho_0 = 1.225

# this is the gladstone-dale constant (m^3/kg)
K = 0.225e-3

# this is the refractive index of the undisturbed medium
n_0 = K * rho_0 + 1

# gradient of the square of the refractive index (m^-1)
alpha = 1e-1

# this is the angle between the initial ray direction and the y axis (radians)
theta_0 = np.pi/2 - np.radians(10.0)
print "theta (rad)", theta_0

# this is a temporary variable that depends on n0 and theta0
xi_0 = n_0 * np.cos(theta_0)

################ generate an array of z co-ordinates ###################
nz = 100

z = np.linspace(start=0, stop=1e3, num = nz, endpoint=True)
################## calculate light ray trajectory #####################
r = 2 * xi_0/alpha * (np.sqrt(-xi_0**2 + n_0**2 + alpha * z) - np.sqrt(-xi_0**2 + n_0**2))

print "curvature", abs(np.diff(np.diff(r))).max()
# print "curvature", np.diff(np.diff(r))
################## plot results #######################################
title_string = r'$\alpha = %d, \theta_0 = %.0f^\circ $' % (alpha, np.degrees(theta_0))
figure_save_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/results/validation/plots/'

plt.plot(z,r)
plt.xlabel('z')
plt.ylabel('r')
plt.title(title_string)
plt.savefig(figure_save_filepath + 'trajectory.png')
plt.show()
