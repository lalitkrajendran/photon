'''
The purpose of this script is to generate a refractive index field that will provide

'''

import numpy as np
import numpy.linalg as la

# this is the distance of the center of the refractive index medium from the camera (microns)
Z_D = 25e3

# this is the dimension of a pixel on the camera sensor (microns)
l_p = 17

# this is the Magnification
M = 0.1761

# this is the depth of the refractive index medium (microns)
W = 50e3

# this is the minimum allowable displacement in pixels at the camera sensor
D_min = 0.

# this is the maximum allowable displacement in pixels at the camera sensor
D_max = 8.

# this is the minimum x co-ordinate of the refractive index field
x_min = -5e4

# this is the maximum x co-ordinate of the refractive index field
x_max = 5e4

# this is the number of grid points along x
x_num = 200

# this is the f number of the camera 
f_number = 8

# this is the focal length of the camera lens (microns)
f = 105e3

# this is the object distance (microns)
s = 700e3

# this is the matrix containing the coefficients of the unknowns
A = np.array([[x_min**2/2, x_min, 1], [x_min, 1, 0], [x_max, 1, 0]])
#print "A"
#print A

# this is the vector containing the values on the right hand side
b = np.array([1, D_min*l_p/(Z_D*M*W), D_max*l_p/(Z_D*M*W)])
#print "b"
#print b

# this solves the system of equations
x = la.solve(A,b)
#print "x"
#print x

# generate the refractive index field
x_loc = np.linspace(start=x_min, stop=x_max, num=x_num, endpoint=True)
#print "x_loc"
#print x_loc

n = x[0]*np.power(x_loc,2)/2 + x[1]*x_loc + x[2]
#print "n"
#print n

# generate the displacement field
d = Z_D*M/l_p * W * np.diff(n)/np.diff(x_loc)
#print "d"
#print d

# verify that the constraints are met
print "n(x_min): ", n[0]
print "d(x_min): ", d[0]
print "d(x_max): ", d[-1]

# calculate the viewing angle in radians
theta = np.arctan(f/(2*f_number*s)) 
print "viewing angle (degrees):", np.degrees(theta)

# calculate the size of blur spot
blur = 2/3. * x[0] * theta * W**3 / l_p
print "blur (pixels):", blur
