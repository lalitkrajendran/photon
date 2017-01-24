
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

# this is the distance of the center of the refractive index medium from the camera (microns)
Z_D = 25e3

# this is the dimension of a pixel on the camera sensor (microns)
l_p = 17

# this is the Magnification
M = 0.1761

# this is the depth of the refractive index medium (microns)
W = 50e3

# this is the maximum allowable displacement in pixels at the camera sensor
D_max = 5.

# this is the minimum x co-ordinate of the refractive index field
x_min = -5e4

# this is the maximum x co-ordinate of the refractive index field
x_max = 5e4

# this is the number of grid points along x
x_num = 200

# this is the minimum x co-ordinate of the refractive index field that is parabolic
x_parabolic_min = -5e3

# this is the maximum x co-ordinate of the refractive index field that is parabolic
x_parabolic_max = 5e3

# this is the number of grid points along x in the parbolic region
x_parabolic_num = int(np.round(x_num * 1.0 * (x_parabolic_max - x_parabolic_min) / (x_max - x_min)))

# this is the total field of view
del_x = x_parabolic_max - x_parabolic_min

# this is the f number of the camera 
f_number = 2.8

# this is the focal length of the camera lens (microns)
f = 105e3

# this is the object distance (microns)
s = 700e3

# this is the required blur (pixels)
blur = 5

# this is the Gladstone-Dale constant for air (m^3/kg)
K = 0.226*1e-3

# this is the density of air at STP (kg/m^3)
rho_air = 1.225

# this is the refractive index of air at STP
n_air = 1 + K*rho_air

# calculate the viewing angle in radians
theta = np.arctan(f/(2*f_number*s)) 

'''
# this is the value of alpha
alpha = D_max * l_p /(Z_D * M * W) * 4/3. * theta * W**3 / (blur * l_p)

print "alpha: ", alpha
print "del_x: ", del_x
print "2*alpha > del_x :" + str(2*alpha > del_x)
'''

# calculate K3
K3 = (x_min + x_max)/2

# calculate max possible value of K2
K2_max = 1/2. * D_max * l_p /(Z_D * M * W) * 2 / del_x
print "K2_max: ", K2_max

# calculate max possible blur
blur_max = 4 / 3. * theta * W**3 / l_p * K2_max
print "max achievable blur (pixels): %.2f" % blur_max
print "actual blur (pixels): %.2f" % blur

# calculate actual value of K2
K2 = 3/4. * blur * l_p / (theta * W**3)
print "K2: ", K2

# calculate K1
K1 = n_air + K2 * del_x**2 / 4.
print "Coefficients: %.2e, %.2e, %.2e" % (K1, K2, K3)

'''
# generate the refractive index field
x_loc = np.linspace(start=x_min, stop=x_max, num=x_num, endpoint=True)
n = K1 - K2 * np.power(x_loc - K3, 2)

# generate the displacement field
d = Z_D*M/l_p * W * np.diff(n)/np.diff(x_loc)

# generate the density field (kg/m^3)
rho = (n - 1)/K
'''

# generate the refractive index field
x_1 = np.linspace(start=x_min, stop=x_parabolic_min, num=np.round((x_num - x_parabolic_num)/2.), endpoint=False)
x_parabolic = np.linspace(start=x_parabolic_min, stop=x_parabolic_max, num=x_parabolic_num, endpoint=False)
x_2 = np.linspace(start=x_parabolic_max, stop=x_max, num=np.round((x_num - x_parabolic_num)/2.), endpoint=False)
print np.shape(x_1)
print np.shape(x_parabolic)
print np.shape(x_2)
#x = np.array([x_1, x_parabolic, x_2])
x = np.append(np.append(x_1, x_parabolic), x_2)

n_1 = n_air * np.ones(np.shape(x_1))
n_parabolic = K1 - K2 * np.power(x_parabolic - K3, 2)
n_2 = n_air * np.ones(np.shape(x_2))
n = np.append(np.append(n_1, n_parabolic), n_2)

'''
print "repeated elements"
print "x"
print np.repeat(x,2)
print "n"
print np.repeat(n,2)
'''

# generate the displacement field
d = Z_D*M/l_p * W * np.diff(n)/np.diff(x)
d = np.append(0,d)

# generate the density field (kg/m^3)
rho = (n - 1)/K

plt.figure(1)
plt.plot(x,n, '.')
plt.title("n")
plt.grid()
plt.draw()

plt.figure(2)
plt.plot(x,d)
plt.title("d")
plt.grid()
plt.draw()

plt.figure(3)
plt.plot(x,rho)
plt.title("density")
plt.grid()
plt.draw()


# verify that the constraints are met
print "n_min: %.5f" % np.min(n)
print "n_max: %.5f" % np.max(n)
print "d_min: %.2f" % np.min(d)
print "d_max: %.2f" % np.max(d)
print "rho_min: %.2f" % np.min(rho)
print "rho_max: %.2f" % np.max(rho)

plt.show()
