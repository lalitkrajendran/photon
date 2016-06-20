#cython: boundscheck=False
import numpy as np
cimport numpy as np

cdef double[:,:] increment_pixels(double[:,:] I, double[:,:] pixel_weights, int[:,:] ii_indices, int[:,:] jj_indices, \
double[:] radiance,double[:] cos_4_alpha, int N):
    # % These nested for loops complete the same calculation as the above
    # % block of code, but require about 2.5 times longer to complete
    # %
    # % This iterates through the pixel locations upon which a ray intersected
    # % incrementing the image intensities
    cdef int n
    cdef int p

    for n in range(0, N):
        # % This iterates through the four adjacent pixels to the
        # % intersection point of the light ray
        for p in range(0, 4):
            # % This interpolates the ray intensities onto the pixels
            # % surrounding the pixel intersection point
            # print "n: %d, p: %d, ii_indices[n][p]: %d, jj_indices[n][p]: %d, pixel_weights[n][p]: %f " % (n,p,ii_indices[n][p],jj_indices[n][p],pixel_weights[n][p])
            I[ii_indices[n][p]-1][jj_indices[n][p]-1] += pixel_weights[n][p] * radiance[n] * \
                                              cos_4_alpha[n]

    return I


def increment_pixel_radiance(np.ndarray I, np.ndarray pixel_weights, np.ndarray ii_indices, np.ndarray jj_indices, \
np.ndarray radiance,np.ndarray cos_4_alpha, int N):
    # % These nested for loops complete the same calculation as the above
    # % block of code, but require about 2.5 times longer to complete
    # %
    # % This iterates through the pixel locations upon which a ray intersected
    # % incrementing the image intensities

    cdef double[:,:] I2 = I
    cdef double[:,:] pixel_weights_2 = pixel_weights
    cdef int[:,:] ii_indices_2 = ii_indices
    cdef int[:,:] jj_indices_2 = jj_indices
    cdef double[:] radiance_2 = radiance
    cdef double[:] cos_4_alpha_2 = cos_4_alpha

    I2 = increment_pixels(I2, pixel_weights_2, ii_indices_2, jj_indices_2, radiance_2,cos_4_alpha_2, N)

    return np.asarray(I2)
