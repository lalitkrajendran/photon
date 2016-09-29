import numpy as np
import ctypes

# ctypes pointer to 2D array
    _floatpp = npct.ndpointer(dtype=np.uintp, flags='C')
    _floatp = npct.ndpointer(dtype=np.float32,flags='C')
    _doublep = npct.ndpointer(dtype=np.float64, flags='C')
    _intp = npct.ndpointer(dtype=np.int32, flags='C')  # define a class to represent the scattering_data dictionary
class scattering_data_struct(ctypes.Structure):
    _fields_ = [
        ('inverse_rotation_matrix', ctypes.c_float * 9),  # _floatpp),
        ('beam_propogation_vector', ctypes.c_float * 3),  # ctypes.POINTER(ctypes.c_float)),
        ('scattering_angle', _floatp),
        ('scattering_irradiance', _floatp),
        ('num_angles', ctypes.c_int),
        ('num_diameters', ctypes.c_int)
    ]


class lightfield_source_struct(ctypes.Structure):
    _fields_ = [
        ('lightray_number_per_particle', ctypes.c_int),
        ('source_point_number', ctypes.c_int),
        ('diameter_index', _intp),
        ('radiance', _doublep),
        ('x', _floatp),
        ('y', _floatp),
        ('z', _floatp),
        ('num_rays', ctypes.c_int),
        ('num_particles', ctypes.c_int)

    ]

    # define a class to hold element_geometry


class element_geometry_struct(ctypes.Structure):
    _fields_ = [
        ('front_surface_radius', ctypes.c_float),
        # ('front_surface_shape', ctypes.c_char_p),
        ('front_surface_spherical', ctypes.c_bool),
        ('back_surface_radius', ctypes.c_float),
        # ('back_surface_shape', ctypes.c_char_p),
        ('back_surface_spherical', ctypes.c_bool),
        ('pitch', ctypes.c_float),
        ('vertex_distance', ctypes.c_double)
    ]

    # define a class to hold element properties


class element_properties_struct(ctypes.Structure):
    _fields_ = [
        ('abbe_number', ctypes.c_float),
        ('absorbance_rate', ctypes.c_float),
        ('refractive_index', ctypes.c_double),
        ('thin_lens_focal_length', ctypes.c_float),
        ('transmission_ratio', ctypes.c_float)
    ]


    # define a class to hold element_data


class element_data_struct(ctypes.Structure):
    _fields_ = [
        ('axial_offset_distances', ctypes.c_double * 2),
        ('element_geometry', element_geometry_struct),
        ('element_number', ctypes.c_float),
        ('element_properties', element_properties_struct),
        # ('element_type', ctypes.c_char_p),
        ('element_type', ctypes.c_char),
        ('elements_coplanar', ctypes.c_float),
        ('rotation_angles', ctypes.c_double * 3),
        ('z_inter_element_distance', ctypes.c_float)
    ]


# create structure to hold camera_design information
class camera_design_struct(ctypes.Structure):
    _fields_ = [
        ('pixel_bit_depth', ctypes.c_int),
        ('pixel_gain', ctypes.c_float),
        ('pixel_pitch', ctypes.c_float),
        ('x_camera_angle', ctypes.c_float),
        ('y_camera_angle', ctypes.c_float),
        ('x_pixel_number', ctypes.c_int),
        ('y_pixel_number', ctypes.c_int),
        ('z_sensor', ctypes.c_float)
    ]
