import numpy as np
from numpy import linalg as la
import run_piv_simulation_02
import scipy.io as sio
from progressbar import ProgressBar
import subprocess
import os
import time
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import scipy.integrate as siint
from scipy.interpolate import interp1d
import string

def create_element_coordinate_data(optical_element):
    # % This function creates data to describe the location of the current
    # % optical element within the entire optical system.
    #
    # % This extracts the optical element type (which must be 'lens', 'aperture',
    # % or 'mirror')
    element_type = optical_element['element_type']
    # % This checks whether the optical element is of the appropriate type
    if element_type != 'lens' and element_type != 'aperture' and element_type != 'mirror':
        # % This returns an error stating that the optical element type is
        # % incorrect (probably due to the element in fact being a system)
        raise ValueError('The optical element type is not ''lens'', ''aperture'', or ''mirror''.')

    # % This extracts the pitch of the optical element
    pitch = optical_element['element_geometry']['pitch']

    # % These are the coordinates of the center of the current element
    x_element_center = 0.0
    y_element_center = 0.0
    z_element_center = 0.0

    # % This is a vector of the element center coordiantes
    element_center = np.array([x_element_center, y_element_center, z_element_center])

    # % This is the pitch of the lens about the center coordinate
    element_pitch = pitch / 2.0

    # % These are the parameters define the plane upon which the current
    # % element is centered
    a_element_plane_parameter = 0.0
    b_element_plane_parameter = 0.0
    c_element_plane_parameter = 1.0
    d_element_plane_parameter = 0.0

    # % This is a vector of the parameters defining the plane of the current
    # % element
    element_plane_parameters = np.array(
        [a_element_plane_parameter, b_element_plane_parameter, c_element_plane_parameter, d_element_plane_parameter])

    return [element_center, element_pitch, element_plane_parameters]


def create_single_element_parameters(optical_element):
    # % This function creates the patch objects for a single optical element that
    # % must have a type of 'lens', 'aperture', or 'mirror' and returns the patch
    # % arrays in the cell arrays 'x_patch_array', 'y_patch_array', and
    # % 'element_type_array'.
    #
    # % This is the offset of the current element from the optical axis
    axial_offset_distances = optical_element['axial_offset_distances']

    # % This is the rotation angles for the current optical element
    rotation_angles = optical_element['rotation_angles']
    # % These are the rotation angles in each dimension
    x_rotation_angle = rotation_angles[0]
    y_rotation_angle = rotation_angles[1]
    z_rotation_angle = rotation_angles[2]

    # % This calculates the patch coordiantes for the current optical element
    [element_center, element_pitch, element_plane_parameters] = create_element_coordinate_data(optical_element)

    # This calculates the rotation matrix
    rotation_matrix = run_piv_simulation_02.calculate_rotation_matrix(x_rotation_angle, y_rotation_angle,
                                                                      z_rotation_angle)

    # % This rotates the vector defining the plane of the current optical element
    rotated_plane_parameters = rotation_matrix * element_plane_parameters[0:3].reshape([3, 1])
    # % This extracts the rotated optical element plane vector
    element_plane_parameters[0:3] = np.transpose(rotated_plane_parameters)

    # % This translates the patch coordinates by the axis offset
    element_center[0] = element_center[0] + axial_offset_distances[0]
    element_center[1] = element_center[1] + axial_offset_distances[1]

    # % This ensures that the element plane parameters first three elements are a
    # % unit vecotr (because unit vectors are cool)
    element_plane_parameters /= la.norm(element_plane_parameters[0:4])

    return [element_center, element_pitch, element_plane_parameters]


def create_element_coordinate_arrays(optical_element, elements_coplanar, element_type, element_count, element_data):
    # % This function recursively scans the optical elements structure given by
    # % 'optical_element' to produce patch cell arrays given by 'x_patch_array',
    # % 'y_patch_array', and 'element_type_array' to produce a series of patch
    # % objects to plot.

    # % This is the number of elements within the current optical element
    # % structure
    element_number = optical_element['element_number']
    #
    # % This initializes the parameters describing the optical elements positions
    element_center = []
    element_pitch = []
    element_plane_parameters = []
    element_thickness = []
    element_system_index = []

    # % This initializes the optical element system index
    element_system_index_temp = 0

    # % This is the total distance of all of the element plotted so far
    total_element_distance = 0.0

    # % This is the distance across the sub-systems of the optical elements
    system_distance = 0.0
    #
    # % This iterates through the sub optical elements
    for element_index in range(0, element_number):

        # % This extracts the current optical element of the current system
        current_optical_element = optical_element['optical_element'][element_index]

        # % This is the type of the current optical element
        current_element_type = current_optical_element['element_type']

        # % If the current element type is a system, then this recursively call
        # % this function
        if (current_element_type == 'system'):

            # % This is a Boolean value stating whether the current elements are coplanar
            # % (this value may only be true for an element type of 'system' and will be
            # % Null for single optical elements)
            elements_coplanar = current_optical_element['elements_coplanar']

            # % This is the distance between the current element and the next element
            # % along the z axis
            z_inter_element_distance = current_optical_element['z_inter_element_distance']

            # % This is the offset of the current element from the optical axis
            axial_offset_distances = current_optical_element['axial_offset_distances']

            # % This is the rotation angles for the current optical element
            rotation_angles = current_optical_element['rotation_angles']
            # % These are the rotation angles in each dimension
            x_rotation_angle = rotation_angles[0]
            y_rotation_angle = rotation_angles[1]
            z_rotation_angle = rotation_angles[2]

            theta_x = x_rotation_angle
            theta_y = y_rotation_angle
            theta_z = z_rotation_angle

            # % This calculates the X axis rotation matrix
            x_rotation_matrix = np.matrix(
                [[1.0, 0.0, 0.0], [0.0, np.cos(theta_x), np.sin(theta_x)], [0.0, -np.sin(theta_x), np.cos(theta_x)]])

            # % This calculates the Y axis rotation matrix
            y_rotation_matrix = np.matrix(
                [[np.cos(theta_y), 0.0, -np.sin(theta_y)], [0.0, 1.0, 0.0], [np.sin(theta_y), 0.0, np.cos(theta_y)]])

            # % This calculates the Z axis rotation matrix
            z_rotation_matrix = np.matrix(
                [[np.cos(theta_z), np.sin(theta_z), 0.0], [-np.sin(theta_z), np.cos(theta_z), 0.0], [0.0, 0.0, 1.0]])

            # % This calculates the full rotation matrix
            rotation_matrix = x_rotation_matrix * y_rotation_matrix * z_rotation_matrix;

            # % This recursively calls this function to produce the patch objects
            # % of the sub-system
            [element_center_temp, element_pitch_temp, element_plane_parameters_temp, element_type,
             element_thickness_temp,
             total_element_distance_temp, element_system_index_temp, element_count,
             element_data] = create_element_coordinate_arrays(current_optical_element, elements_coplanar, element_type,
                                                              element_count, element_data)

            # % These are the coordinates about which the rotation occurs
            x0 = 0.0
            y0 = 0.0
            z0 = total_element_distance_temp / 2.0

            # % This is the vector of point about which the rotation occurs
            rotation_origin = np.array([x0, y0, z0])

            # % This iterates through optical element parameters for the current
            # % system, rotating them and offsetting them by the specified
            # % quantitites for the current element
            for element_temp_index in range(0, np.shape(element_center_temp)[0]):
                # % This extracts the circle center coordinates
                x0 = element_center_temp[element_temp_index][0]
                y0 = element_center_temp[element_temp_index][1]
                z0 = element_center_temp[element_temp_index][2]

                # % This extracts the circle pitch
                R = element_pitch_temp[element_temp_index]
                # % This extracts the plane parameters
                a = element_plane_parameters_temp[element_temp_index][0]
                b = element_plane_parameters_temp[element_temp_index][1]
                c = element_plane_parameters_temp[element_temp_index][2]
                d = element_plane_parameters_temp[element_temp_index][3]

                # % This is a point that is known to lie on the plane
                xp = 0.0
                yp = 0.0
                zp = -d / c

                # % This is the normal vector of the plane
                plane_normal = np.array([a, b, c])
                # % This calculates the rotated normal vector of the plane
                plane_normal = rotation_matrix * plane_normal.reshape([3, 1])
                # % This extracts the new plane normal components
                a_temp = plane_normal.item(0)
                b_temp = plane_normal.item(1)
                c_temp = plane_normal.item(2)

                # % This is the vector of the point known to lie on the plane
                plane_point = np.array([xp, yp, zp])
                # % This calculates the transoformed point on the plane
                plane_point = np.squeeze(np.asarray(
                    rotation_matrix * (plane_point - rotation_origin).reshape([3, 1]) + rotation_origin.reshape(
                        [3, 1])))
                # % This extracts the transformed coordinates
                xp = plane_point[0]
                yp = plane_point[1]
                zp = plane_point[2]

                # % This is the new constant parameter term
                d_temp = -(a_temp * xp + b_temp * yp + c_temp * zp)

                # % This rotates the points lieing on the center of the circle
                element_center_temp[element_temp_index][:] = np.transpose(rotation_matrix * (
                    (element_center_temp[element_temp_index][:] - rotation_origin).reshape(
                        [3, 1])) + rotation_origin.reshape([3, 1]))

                # % This translates the patch coordinates by the axis offset
                element_center_temp[element_temp_index][0] = element_center_temp[element_temp_index][0] + \
                                                             axial_offset_distances[0]
                element_center_temp[element_temp_index][1] = element_center_temp[element_temp_index][1] + \
                                                             axial_offset_distances[1]

                # % This translates the current element by the total element distance
                element_center_temp[element_temp_index][2] = element_center_temp[element_temp_index][
                                                                 2] + system_distance

                # % This translates the plane parameters constant offset to
                # % account for the axial offsets and the system distances
                d_temp = d_temp - (
                    a_temp * axial_offset_distances[0] + b_temp * axial_offset_distances[1] + c_temp * system_distance);

                # % This saves the new plane parameters to the full plane
                # % parameters array
                element_plane_parameters_temp[element_temp_index][:] = np.array([a_temp, b_temp, c_temp, d_temp])

            # % This increments the total system distance by the current system
            # % distance
            system_distance = system_distance + total_element_distance_temp + z_inter_element_distance

        else:

            # % This is the thickness of the optical element from the front optical axis
            # % vertex to the back optical axis vertex
            vertex_distance = current_optical_element['element_geometry']['vertex_distance']

            # % This is the distance between the current element and the next element
            # % along the z axis
            z_inter_element_distance = current_optical_element['z_inter_element_distance']

            # % This creates the optical element parameters for the current
            # % element
            [element_center_temp, element_pitch_temp, element_plane_parameters_temp] = create_single_element_parameters(
                current_optical_element)

            # % This adds the current total element distance to the element
            # % center z axis
            element_center_temp[2] = element_center_temp[2] + total_element_distance

            # % This translates the plane offset parameter
            element_plane_parameters_temp[3] = -element_plane_parameters_temp[2] * element_center_temp[2]

            # % This tests whether the current element is the last element in the
            # % system and if the elements are co-planar, this adds the vertex
            # % distance of the the last element
            if element_index < element_number:
                # % This increments the total element distance by the thickness of the
                # % current element plus the distance to the next element
                total_element_distance = total_element_distance + int(
                    not (elements_coplanar)) * vertex_distance + z_inter_element_distance
            else:
                # % This increments the total element distance by the thickness of the
                # % current element plus the distance to the next element
                total_element_distance = total_element_distance + vertex_distance + z_inter_element_distance

            # % This adds the current element type to the full element type cell array
            if not (element_type):
                # % If the full element type cell array is empty, this adds the first element
                # % as a cell
                element_type.append(current_element_type)
            else:
                # % If the full element type array is non-empty, this adds the current element
                # % type to the end of the cell array
                element_type.append(current_element_type)

            # % This sets the temporary optical element thickness value equal to the vertex
            # % distance
            element_thickness_temp = vertex_distance

            # % This increments the optical element system index if the system is
            # % not coplanar
            if not elements_coplanar:
                # % This increments the optical element system index
                element_system_index_temp = element_system_index_temp + 1

            # % This increments the counting variable giving the index of the
            # % current optical element
            element_count = element_count + 1

            # % This saves the current optical elements data parameters into the
            # % output data cell array
            element_data = np.append(element_data, current_optical_element)

        # % This adds the new optical element centers onto the list of previous
        # % optical centers
        element_center.append(element_center_temp)
        # % This adds the new optical element radii onto the list of the previous
        # % optical radii
        element_pitch.append(element_pitch_temp)
        # % This adds the new optical element plane parameters onto the list of
        # % the previous optical element plane parameters
        element_plane_parameters.append(element_plane_parameters_temp)
        # % This adds the new optical element thickness value onto the list of the previous
        # % optical element thickness values
        element_thickness.append(element_thickness_temp)
        # % This adds the new optical element system index onto the list of the
        # % previous system indices
        element_system_index.append(element_system_index_temp)

    return [element_center, element_pitch, element_plane_parameters, element_type, element_thickness,
            total_element_distance, element_system_index, element_count, element_data]


def generate_lightfield_angular_data(lens_pitch, image_distance, scattering_data, scattering_type, lightfield_source,
                                     lightray_number_per_particle, n_min, n_max):
    # % This function generates the lightfield data for the source points specified by the
    # % structure lightfield_source.  The data is only generated for the source points from
    # % n_min to n_max.  The parameter lightray_number_per_particle is the number of rays to generate for each
    # % source point.

    # % If the scattering type is 'mie' then the Mie scattering data is loaded,
    # % otherwise nothing is loaded from the scattering data
    if scattering_type == 'mie':
        #% This saves the scattering angle data into the parameters structure
        mie_scattering_angle = scattering_data['scattering_angle']
        #% This saves the scattering irradiance values for the different particle
        #% diameters into the parameters structure
        mie_scattering_irradiance=scattering_data['scattering_irradiance']
        #% This extracts the inverse rotation matrix from the parameters structure
        inverse_rotation_matrix=scattering_data['inverse_rotation_matrix']
        #% This extracts the normalized beam propogation direction vector from the
        #% parameters structure
        beam_propogation_vector=scattering_data['beam_propogation_vector']

    # % This is the number of points to calculate the lightfield structure for
    source_point_number = n_max - n_min + 1

    # % This initializes the data structure
    lightfield_data = {}
    lightfield_data['x'] = np.zeros(source_point_number * lightray_number_per_particle)
    lightfield_data['y'] = np.zeros(source_point_number * lightray_number_per_particle)
    lightfield_data['z'] = np.zeros(source_point_number * lightray_number_per_particle)
    lightfield_data['theta'] = np.zeros(source_point_number * lightray_number_per_particle)
    lightfield_data['phi'] = np.zeros(source_point_number * lightray_number_per_particle)
    lightfield_data['radiance'] = np.zeros(source_point_number * lightray_number_per_particle)

    # % This is the radius of the current lens for which the light rays are being
    # % generated
    R = lens_pitch / 2.0

    # % This iterates through the "particle" locations generating the light rays
    for n in range(n_min,n_max+1):
        # % This is the current source point for which the lightfield is being generated
        x_current = np.squeeze(np.asarray(lightfield_source['x']))[n]
        y_current = np.squeeze(np.asarray(lightfield_source['y']))[n]
        z_current = np.squeeze(np.asarray(lightfield_source['z']))[n]

        # % This creates random radial coordinates for the lightrays to intersect
        # % on the lens
        np.random.seed(523)
        r = R * np.sqrt(np.random.rand(1, lightray_number_per_particle))

        # % This creates random angular coordinates for the lightrays to
        # % intersect on the lens
        np.random.seed(105)
        psi = 2.0 * np.pi * np.random.rand(1, lightray_number_per_particle)

        # % This calculates the random cartesian coordinate of the points the
        # % rays will intersect on the lens
        x_lens = r * np.cos(psi)
        y_lens = r * np.sin(psi)

        # % This calculates the angular data using Mie scattering if specified in the
        # % parameters structure, otherwise the particles are assummed to have
        # % uniform irradiance
        if scattering_type=='mie':

            # % This extracts the current particle diameter index
            diameter_index=lightfield_source['diameter_index'][n]

            # % This calculates the lightrays direction vectors (in the camera
            # % coordiante system)
            ray_direction_vector = np.array([np.squeeze(x_lens-x_current),np.squeeze(y_lens-y_current),image_distance*np.ones(lightray_number_per_particle)-np.squeeze(z_current)])

            # % This normalizes the ray direction vectors
            ray_direction_vector /=np.sqrt(ray_direction_vector[0,:]**2+ray_direction_vector[1,:]**2+ray_direction_vector[2,:]**2)[np.newaxis,:]

            # % This rotates the lightrays direction vectors by the inverse of the
            # % camera rotation array so that the ray is now in the world coordinate
            # % system
            ray_direction_vector=inverse_rotation_matrix*ray_direction_vector
            # % This calculates the angle that the light ray direction vectors make
            # % with the laser propogation direction vector
            ray_scattering_angles=np.squeeze(np.arccos(beam_propogation_vector*ray_direction_vector))

            # % This calculates the Mie scattering irradiances at the currently
            # % scattered angles and with the current particle diameter
            set_interp = interp1d(mie_scattering_angle,mie_scattering_irradiance[:,diameter_index])
            ray_scattering_irradiance = np.squeeze(set_interp(ray_scattering_angles))
            # ray_scattering_irradiance=interp1d(mie_scattering_angle,mie_scattering_irradiance(:,diameter_index),ray_scattering_angles,'linear')

            # % This calculates the total irradiance for the current particle's rays
            irradiance_current=ray_scattering_irradiance*lightfield_source['radiance'][n]

        elif scattering_type =='diffuse':

            # % This specifies the total irradiance for the current particle's
            # % rays to be uniform
            irradiance_current = lightfield_source['radiance'][n]

        # % This calculates the x angles for the light rays
        theta_temp = -(np.squeeze(x_lens)- x_current) / (image_distance - z_current)
        # % This calculates the y angles for the light rays
        phi_temp = -(np.squeeze(y_lens) - y_current) / (image_distance - z_current)

        # % This is the vector of indicies to save the light rays to
        index_vector = np.array((n - n_min) * lightray_number_per_particle + np.arange(0, lightray_number_per_particle)).astype('int')
        # index_vector -= 1

        # % This saves the light rays to the data structure
        lightfield_data['x'][index_vector] = x_current
        lightfield_data['y'][index_vector] = y_current
        lightfield_data['z'][index_vector] = z_current
        lightfield_data['theta'][index_vector] = theta_temp
        lightfield_data['phi'][index_vector] = phi_temp
        lightfield_data['radiance'][index_vector] = irradiance_current

    return lightfield_data


def ray_sphere_intersection(xc, yc, zc, R, vx, vy, vz, x0, y0, z0, surface):
    # % This function returns the point at which a ray first intersects a sphere
    # % in 3-space.  The sphere is defined by it's center coordinates given by
    # % the argumnets 'xc', 'yc', and 'zc' and it's radius given by the argument
    # % 'R'.  Thus the equation for the surface of the sphere is
    # %
    # %   (x-xc)^2+(y-yc)^2+(z-zc)^2=R^2
    # %
    # % The ray is defined by the direction vector given by the arguments 'vx',
    # % 'vy', and 'vz' and it's origin coordinate given by the arguments 'x0',
    # % 'y0', and 'z0'.  Thus the parametric equations defining the ray are
    # %
    # %   x(t)=x0+vx*t
    # %   y(t)=y0+vy*t
    # %   z(t)=z0+vz*t
    # %
    # % The string 'surface' may equal either 'front' or 'back' and states
    # % whether the light ray is entering the front of the lens or exiting the
    # % back of the lens.
    # %
    # % The the output arguments 'xi', 'yi', and 'zi' give the coordinate of the
    # % first intersection or the ray with the sphere (ie there will be two
    # % intersections and the one output is the one that is closest to the origin
    # % of the ray that is also positive).
    # %
    # % The input arguments can be of any size, but they must all have the same
    # % dimensions.  The output arguments will be the same size as the input
    # % argument.  If a particular ray does not intersection a particular sphere,
    # % then NaN values will be returned in the output arguments.

    # % This defines the coeffients of the quadratic polynomial
    alpha = vx ** 2 + vy ** 2 + vz ** 2
    beta = 2 * (vx * (x0 - xc) + vy * (y0 - yc) + vz * (z0 - zc))
    gamma = (x0 - xc) ** 2 + (y0 - yc) ** 2 + (z0 - zc) ** 2 - R ** 2

    # % This calculates the square root argument of the quadratic formula
    square_root_arguments = beta ** 2 - 4.0 * alpha * gamma

    # % These are the indices of the real solutions (ie the cases where the rays
    # % actually intersect the sphere surface)
    intersection_indices = np.argwhere(square_root_arguments >= 0.0)

    # % This initializes the two possible independent time variables
    t1 = np.empty(alpha.size)
    t1.fill(np.NAN)
    t2 = np.empty(alpha.size)
    t2.fill(np.NAN)

    # % This yields the solutions for the independent time variable for the cases
    # % that do intersect the sphere
    t1[intersection_indices] = (-beta[intersection_indices] + np.sqrt(square_root_arguments[intersection_indices])) / (
    2.0 * alpha[intersection_indices])
    t2[intersection_indices] = (-beta[intersection_indices] - np.sqrt(square_root_arguments[intersection_indices])) / (
    2.0 * alpha[intersection_indices])

    # % This chooses the intersection time based upon whether the front surface
    # % of the lens or the back surface of the lens is being intersected
    if surface == 'front':
        # % If the front surface has positive curvature than the first
        # % intersection point should be taken, otherwise if the surface has
        # % negative curvature, then the second intersection point should be
        # % taken
        if R > 0:
            # % This gives the first intersection time variable (ie the minimum
            # % of 't1' and 't2' that is non-NaN, in the case where both 't1'
            # % and 't2' have NaN values, the output of 't' is also NaN)
            t = np.fmin(t1, t2)
        elif R <= 0:
            # % This gives the second intersection time variable (ie the maximum
            # % of 't1' and 't2' that is non-NaN, in the case where both 't1'
            # % and 't2' have NaN values, the output of 't' is also NaN)
            t = np.fmax(t1, t2)

    elif surface == 'back':
        # % If the back surface has positive curvature than the second
        # % intersection point should be taken, otherwise if the surface has
        # % negative curvature, then the first intersection point should be
        # % taken
        if R > 0:
            # % This gives the second intersection time variable (ie the maximum
            # % of 't1' and 't2' that is non-NaN, in the case where both 't1'
            # % and 't2' have NaN values, the output of 't' is also NaN)
            # %t=nanmax(t1,t2); TODO
            t = np.fmin(t1, t2)

        elif R <= 0:
            # % This gives the first intersection time variable (ie the minimum
            # % of 't1' and 't2' that is non-NaN, in the case where both 't1'
            # % and 't2' have NaN values, the output of 't' is also NaN)
            # %t=nanmin(t1,t2); TODO
            t = np.fmax(t1, t2)

    # % This calculates the 3D coordinate of the intersection point
    xi = x0 + vx * t
    yi = y0 + vy * t
    zi = z0 + vz * t

    return [xi, yi, zi]


def measure_distance_to_optical_axis(xi, yi, zi, lens_center, plane_parameters):
    # % This function measures the distance from the point given by the M x 1
    # % vectors 'xi', 'yi', and 'zi' to the optical axis of the optical element
    # % with a center located at the 3 x N array given by 'lens_center' which
    # % has the lens plane defined by the 4 x N array 'plane_parameters'.  The
    # % coordinates defined by 'lens_center' must satisfy the parameters defined
    # % by 'plane_parameters', ie
    # %
    # %   plane_parameters(ii,1) * lens_center(ii,1) +
    # %       plane_parameters(ii,2) * lens_center(ii,2) +
    # %       plane_parameters(ii,3) * lens_center(ii,3) + plane_parameters(ii,4) = 0
    # %
    # % must be satisfied.  The output argument 'optical_axis_distance' is the
    # % distance from the points given by 'xi, 'yi', and 'zi' to the optical
    # % axis passing through the point 'lens_center' that is normal to the plane
    # % defined by 'plane_parameters' and will be M x N in size.
    if(len(plane_parameters.shape)==1):
        plane_parameters = np.reshape(plane_parameters,(plane_parameters.shape[0],1))
        lens_center = np.reshape(lens_center,(lens_center.shape[0],1))

    # % This extracts the normal vector components of the plane parameters
    a = plane_parameters[0,:]
    b = plane_parameters[1,:]
    c = plane_parameters[2,:]

    # % This extracts the lens center (that must lie within the plane)
    x0 = lens_center[0,:]
    y0 = lens_center[1,:]
    z0 = lens_center[2,:]

    # % This calculates the minimum time (of the parametric equations describing
    # % the optical axis)
    # t_minimum = bsxfun(@rdivide,(bsxfun(@times,a,bsxfun(@minus,xi,x0))+bsxfun(@times,b,bsxfun(@minus,yi,y0))+bsxfun(@times,c,bsxfun(@minus,zi,z0))),(a.^2+b.^2+c.^2));
    t_minimum = (a * (xi - x0) + b * (yi - y0) + c * (zi - z0)) / (a ** 2 + b ** 2 + c ** 2)

    # % This is the nearest point on the optical axis to the points defined by
    # % 'xi', 'yi', and 'zi'
    x_optical_axis = x0 + a * t_minimum
    y_optical_axis = y0 + b * t_minimum
    z_optical_axis = z0 + c * t_minimum

    # % This is the distance to the optical axis
    optical_axis_distance = np.sqrt(
        (xi - x_optical_axis) ** 2 + (yi - y_optical_axis) ** 2 + (zi - z_optical_axis) ** 2)

    return optical_axis_distance


def convert_world_coordinates_to_lens_coordinates(x, y, z, xc, yc, zc, a, b, c):
    # % This function converts the world coordinate arrays given by 'x', 'y', and
    # % 'z' to the lens coordinates output as 'xL', 'yL', and 'zL'.  The input
    # % arguments 'xc', 'yc', and 'zc' give the center coordinate of the lens
    # % while the arguments 'a', 'b', and 'c' define the plane parallel to the
    # % lens, ie the plane defined by
    # %
    # %   a * x + b * y + c * z + d = 0
    # %
    # % will contain the point (xc,yc,zc) and is normal to the optical axis of
    # % the lens (which passes through the point (xc,yc,zc)).
    # %
    # % The input arguments 'x', 'y', and 'z' may be of any size, but must all
    # % have the same dimensions.  The output arguments 'xL', 'yL', and 'zL'
    # % will have the same size as the input arguments 'x', 'y', and 'z'.  The
    # % input arguments 'xc', 'yc', 'zc', 'a', 'b', and 'c' must all be scalars.

    # % This is defines a temporary first normalization quantity (which is used
    # % to define the second normalization quantity)
    bc_norm_squared = b ** 2 + c ** 2
    # % This is defines the second normalization quantity (which is used several
    # % times)
    abc_norm = np.sqrt(a ** 2 + bc_norm_squared)
    # % This is defines the first normalization quantity (which is used several
    # % times)
    bc_norm = np.sqrt(bc_norm_squared)
    #
    # % This calculates the transformed x axis lens coordinates
    xL = (x - xc) * (bc_norm / abc_norm) - (y - yc) * ((a * b) / (bc_norm * abc_norm)) - (z - zc) * (
    (a * c) / (bc_norm * abc_norm))
    # % This calculates the transformed y axis lens coordinates
    yL = (y - yc) * (c / bc_norm) - (z - zc) * (b / bc_norm)
    # % This calculates the transformed z axis lens coordinates
    zL = (x - xc) * (a / abc_norm) + (y - yc) * (b / abc_norm) + (z - zc) * (c / abc_norm)

    return [xL, yL, zL]


def propogate_rays_through_single_element(optical_element, element_center, element_plane_parameters, light_ray_data):
    # % This function calculates the propogation of a set of light rays through a
    # % single optical element.

    # % This extracts the optical element type
    element_type = optical_element['element_type']

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Extraction of the light ray propogation data                            %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # % This extracts the propogation direction of the light rays
    ray_propogation_direction = light_ray_data['ray_propogation_direction']

    # % This extracts the light ray source coordinates
    ray_source_coordinates = light_ray_data['ray_source_coordinates']

    # % This extracts the wavelength of the light rays
    ray_wavelength = light_ray_data['ray_wavelength']

    # % This extracts the light ray radiance
    ray_radiance = light_ray_data['ray_radiance']

    # % If the element type is 'lens', this extracts the optical properties
    # % required to perform the ray propogation
    if element_type == 'lens':

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # % Extraction of the lens optical properties                           %
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        # % This extracts the pitch of the lens
        element_pitch = optical_element['element_geometry']['pitch']

        # % This is the thickness of the optical element from the front optical axis
        # % vertex to the back optical axis vertex
        element_vertex_distance = optical_element['element_geometry']['vertex_distance']

        # % This extracts the front surface radius of curvature
        element_front_surface_curvature = optical_element['element_geometry']['front_surface_radius']

        # % This extracts the back surface radius of curvature
        element_back_surface_curvature = optical_element['element_geometry']['back_surface_radius']

        # % This extracts the refractive index of the current optical element
        element_refractive_index = optical_element['element_properties']['refractive_index']

        # % This extracts the Abbe number of the current optical element
        element_abbe_number = optical_element['element_properties']['abbe_number']

        # % This extracts the transmission ratio of the current optical
        # % element
        element_transmission_ratio = optical_element['element_properties']['transmission_ratio']

        # % This extracts the absorbance rate of the current optical element
        element_absorbance_rate = optical_element['element_properties']['absorbance_rate']

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # % Propogation of the light rays through the lens front surface        %
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        # % This extracts the individual plane parameters
        a = element_plane_parameters[0]
        b = element_plane_parameters[1]
        c = element_plane_parameters[2]
        d = element_plane_parameters[3]

        # % This extracts the center point coordinates
        xc = element_center[0]
        yc = element_center[1]
        zc = element_center[2]

        # % This calculates the square of the plane normal vector magnitude
        norm_vector_magnitude = np.sqrt(a ** 2 + b ** 2 + c ** 2)

        # % This is the current offset to the center of the front spherical
        # % surface
        ds = +element_vertex_distance / 2.0 - element_front_surface_curvature

        # % This calculates the center point coordinates of the front surface
        # % spherical shell
        xc_front_surface = xc + a * ds / norm_vector_magnitude
        yc_front_surface = yc + b * ds / norm_vector_magnitude
        zc_front_surface = zc + c * ds / norm_vector_magnitude

        # % This calculates the points at which the light rays intersect the
        # % front surface of the lens
        [x_intersect, y_intersect, z_intersect] = ray_sphere_intersection(xc_front_surface, yc_front_surface,
                                                                          zc_front_surface,
                                                                          element_front_surface_curvature,
                                                                          ray_propogation_direction[:,0],
                                                                          ray_propogation_direction[:,1],
                                                                          ray_propogation_direction[:,2],
                                                                          ray_source_coordinates[:,0],
                                                                          ray_source_coordinates[:,1],
                                                                          ray_source_coordinates[:,2], 'front')

        # % This calculates how far the intersection point is from the optical
        # % axis of the lens (for determining if the rays hit the lens outside of
        # % the pitch and thus would be destroyed)
        optical_axis_distance = measure_distance_to_optical_axis(x_intersect, y_intersect, z_intersect,
                                                                 element_center.T, element_plane_parameters.T)

        # % This gives the indices of the light rays that actually intersect the
        # % lens (ie the rays that are within the lens pitch)
        intersect_lens_indices = optical_axis_distance <= (element_pitch / 2.0)

        # % This sets any of the light ray directions outside of the domain of
        # % the lens to NaN values
        ray_propogation_direction[np.logical_not(intersect_lens_indices)][:] = np.NAN

        # % This sets any of the light ray wavelengths outside of the  domain of
        # % the lens to NaN values
        ray_wavelength[np.logical_not(intersect_lens_indices)] = np.NAN

        # % This sets any of the light ray radiances outside of the  domain of
        # % the lens to NaN values
        ray_radiance[np.logical_not(intersect_lens_indices)] = np.NAN

        # % This sets the intersection points of any of the light rays outside of
        # % the domain of the lens to NaN values
        x_intersect[np.logical_not(intersect_lens_indices)] = np.NAN
        y_intersect[np.logical_not(intersect_lens_indices)] = np.NAN
        z_intersect[np.logical_not(intersect_lens_indices)] = np.NAN

        # % This calculates the normal vectors of the lens at the intersection
        # % points
        lens_normal_vectors = np.array([x_intersect - xc_front_surface, y_intersect - yc_front_surface,
                               z_intersect - zc_front_surface])

        # % This normalizes the lens normal vectors to have unit magnitudes
        # %lens_normal_vectors=lens_normal_vectors/norm(lens_normal_vectors,2);
        lens_normal_vectors = lens_normal_vectors.T
        lens_normal_vectors /= np.sqrt(
            lens_normal_vectors[:,0] ** 2 + lens_normal_vectors[:,1] ** 2 + lens_normal_vectors[:,2] ** 2)[:,np.newaxis]

        # % If the refractive index is a constant double value, this directly
        # % calculates the refractive index ratio, otherwise the ratio is
        # % calculated as a function of the wavelength
        if isinstance(element_refractive_index,(int,long,float,complex)):
            # % If the Abbe number is defined, this calculates the Cauchy formula
            # % approxiation to the refractive index, otherwise, the refractive
            # % index is defined to be constant
            if not(np.isnan(element_abbe_number)):
                # % This defines the three optical wavelengths used in defining
                # % the Abbe number
                lambda_D = 589.3
                lambda_F = 486.1
                lambda_C = 656.3

                # % This is the ratio of the incident refractive index (1 since the
                # % inter-lens media is assummed to be air) to the transmitted
                # % refractive index
                refractive_index_ratio = 1.0 / (
                element_refractive_index + (1. / (ray_wavelength ** 2) - 1 / (lambda_D ** 2)) *
                ((element_refractive_index - 1) / (element_abbe_number * (1 / (lambda_F ** 2) - 1 / (lambda_C ** 2)))))

            else:
                # % This is the ratio of the incident refractive index (1 since the
                # % inter-lens media is assummed to be air) to the transmitted
                # % refractive index
                refractive_index_ratio = 1.0 / element_refractive_index*np.ones(ray_propogation_direction.shape[0])
        else:
            # % This defines the wavelength of the light as the variable 'lambda'
            # % for evaluation of the refractive index function
            ray_lambda = ray_wavelength
            # % This evaluates the string defining the refractive index in terms
            # % of the independent variable lambda TODO
            element_refractive_index_double = eval('element_refractive_index_double=' + element_refractive_index)
            # % This is the ratio of the incident refractive index (1 since the
            # % inter-lens media is assummed to be air) to the transmitted
            # % refractive index
            refractive_index_ratio = 1.0 / element_refractive_index_double*np.ones(ray_propogation_direction.shape[0])

        # % This is the scaled cosine of the angle of the incident light ray
        # % vectors and the normal vectors of the lens (ie the dot product of the
        # % vectors)
        ray_dot_product = -np.diag(np.dot(ray_propogation_direction, lens_normal_vectors.T))
        # % This calculates the radicand in the refraction ray propogation
        # % direction equation
        refraction_radicand = 1.0 - (refractive_index_ratio ** 2) * (1.0 - ray_dot_product ** 2)
        # % This calculates the new light ray direction vectors (this is a
        # % standard equation in optics relating incident and transmitted light
        # % ray vectors)
        # ray_propogation_direction = bsxfun(@times,refractive_index_ratio,ray_propogation_direction)+bsxfun(@times,(refractive_index_ratio.*ray_dot_product-sqrt(refraction_radicand)),lens_normal_vectors);
        ray_propogation_direction = ray_propogation_direction*refractive_index_ratio[:,np.newaxis]  + \
                                    (refractive_index_ratio * ray_dot_product - np.sqrt(refraction_radicand))[:,np.newaxis] * lens_normal_vectors

        # % This normalizes the ray propogation direction so that it's magnitude
        # % equals one
        # ray_propogation_direction=bsxfun(@rdivide,ray_propogation_direction,sqrt(ray_propogation_direction(:,1).^2+ray_propogation_direction(:,2).^2+ray_propogation_direction(:,3).^2));
        ray_propogation_direction /= np.sqrt(
            ray_propogation_direction[:,0] ** 2 + ray_propogation_direction[:,1] ** 2 + ray_propogation_direction[:,
                2] ** 2)[:,np.newaxis]

        # % This sets the new light ray origin to the intersection point with the
        # % front surface of the lens
        ray_source_coordinates = np.array([x_intersect, y_intersect, z_intersect])
        ray_source_coordinates = ray_source_coordinates.T
        # % Any rays that have complex values due to the radicand in the above
        # % equation being negative experience total internal reflection.  The
        # % values of these rays are set to NaN for now.
        tir_indices = np.less(refraction_radicand, 0)

        # % This sets any of the light ray directions experiencing total
        # % internal reflection to NaN values
        ray_propogation_direction[tir_indices,:] = np.NAN

        # % This sets any of the light ray origin coordinates experiencing total
        # % internal reflection to NaN values
        ray_source_coordinates[tir_indices,:] = np.NAN

        # % This sets any of the light ray wavelengths experiencing total
        # % internal reflection to NaN values
        ray_wavelength[tir_indices] = np.NAN

        # % This sets any of the light ray radiances experiencing total internal
        # % reflection to NaN values
        ray_radiance[tir_indices] = np.NAN

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # % Propogation of the light rays through the lens back surface         %
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        # % This is the current offset to the center of the back spherical
        # % surface
        ds = -element_vertex_distance / 2 - element_back_surface_curvature

        # % This calculates the center point coordinates of the back surface
        # % spherical shell
        xc_back_surface = xc + a * ds / norm_vector_magnitude
        yc_back_surface = yc + b * ds / norm_vector_magnitude
        zc_back_surface = zc + c * ds / norm_vector_magnitude

        # % This calculates the points at which the light rays intersect the
        # % back surface of the lens
        [x_intersect, y_intersect, z_intersect] = ray_sphere_intersection(xc_back_surface, yc_back_surface,
                                                                          zc_back_surface,
                                                                          element_back_surface_curvature,
                                                                          ray_propogation_direction[:,0],
                                                                          ray_propogation_direction[:,1],
                                                                          ray_propogation_direction[:,2],
                                                                          ray_source_coordinates[:,0],
                                                                          ray_source_coordinates[:,1],
                                                                          ray_source_coordinates[:,2], 'back')

        # % This calculates how far the intersection point is from the optical
        # % axis of the lens (for determining if the rays hit the lens outside of
        # % the pitch and thus would be destroyed)
        optical_axis_distance = measure_distance_to_optical_axis(x_intersect, y_intersect, z_intersect,
                                                                 element_center.T, element_plane_parameters.T)

        # % This gives the indices of the light rays that actually intersect the
        # % lens (ie the rays that are within the lens pitch)
        intersect_lens_indices = np.less_equal(optical_axis_distance,(element_pitch / 2))

        # % This sets any of the light ray directions outside of the domain of
        # % the lens to NaN values
        ray_propogation_direction[np.logical_not(intersect_lens_indices),:] = np.NAN

        # % This sets any of the light ray origin coordinates outside of the
        # % domain of the lens to NaN values
        ray_source_coordinates[np.logical_not(intersect_lens_indices),:] = np.NAN

        # % This sets any of the light ray wavelengths outside of the  domain of
        # % the lens to NaN values
        ray_wavelength[np.logical_not(intersect_lens_indices)] = np.NAN

        # % This sets any of the light ray radiances outside of the  domain of
        # % the lens to NaN values
        ray_radiance[np.logical_not(intersect_lens_indices)] = np.NAN

        # % This sets the intersection points of any of the light rays outside of
        # % the domain of the lens to NaN values
        x_intersect[np.logical_not(intersect_lens_indices)] = np.NAN
        y_intersect[np.logical_not(intersect_lens_indices)] = np.NAN
        z_intersect[np.logical_not(intersect_lens_indices)] = np.NAN

        # % This calculates the normal vectors of the lens at the intersection
        # % points
        lens_normal_vectors = -1 * np.array(
            [x_intersect - xc_back_surface, y_intersect - yc_back_surface, z_intersect - zc_back_surface])
        lens_normal_vectors = lens_normal_vectors.T

        # % This normalizes the lens normal vectors to have unit magnitudes
        lens_normal_vectors /= np.sqrt(
            lens_normal_vectors[:,0] ** 2 + lens_normal_vectors[:,1] ** 2 + lens_normal_vectors[:,2] ** 2)[:,np.newaxis]

        # % If the refractive index is a constant double value, this directly
        # % calculates the refractive index ratio, otherwise the ratio is
        # % calculated as a function of the wavelength
        if isinstance(element_refractive_index,(int,long,float,complex)):
            # % If the Abbe number is defined, this calculates the Cauchy formula
            # % approxiation to the refractive index, otherwise, the refractive
            # % index is defined to be constant
            if not(np.isnan(element_abbe_number)):
                # % This defines the three optical wavelengths used in defining
                # % the Abbe number
                lambda_D = 589.3
                lambda_F = 486.1
                lambda_C = 656.3
                # % This is the ratio of the incident refractive index to the transmitted
                # % refractive index (1 since the  inter-lens media is assummed to be
                # % air)
                refractive_index_ratio = element_refractive_index + (1. / (ray_wavelength ** 2) - 1 / (
                lambda_D ** 2)) * ((element_refractive_index - 1) / (
                element_abbe_number * (1 / (lambda_F ** 2) - 1 / (lambda_C ** 2))))
            else:
                # % This is the ratio of the incident refractive index to the transmitted
                # % refractive index (1 since the  inter-lens media is assummed to be
                # % air)
                refractive_index_ratio = element_refractive_index*np.ones(ray_propogation_direction.shape[0])
        else:
            # % This defines the wavelength of the light as the variable 'lambda'
            # % for evaluation of the refractive index function
            ray_lambda = ray_wavelength
            # % This evaluates the string defining the refractive index in terms
            # % of the independent variable lambda TODO
            element_refractive_index_double = eval('element_refractive_index_double=' + element_refractive_index)
                # % This is the ratio of the incident refractive index (1 since the
            # % inter-lens media is assummed to be air) to the transmitted
            # % refractive index
            refractive_index_ratio = element_refractive_index_double*np.ones(ray_propogation_direction.shape[0])

        # % This is the scaled cosine of the angle of the incident light ray
        # % vectors and the normal vectors of the lens (ie the dot product of the
        # % vectors)
        ray_dot_product = -np.diag(np.dot(ray_propogation_direction, lens_normal_vectors.T))

        # % This calculates the radicand in the refraction ray propogation
        # % direction equation
        refraction_radicand = 1.0 - (refractive_index_ratio ** 2) * (1.0 - ray_dot_product ** 2)
        # % This calculates the new light ray direction vectors (this is a
        # % standard equation in optics relating incident and transmitted light
        # % ray vectors)
        ray_propogation_direction = refractive_index_ratio[:,np.newaxis] * ray_propogation_direction + (
                                                                                         refractive_index_ratio * ray_dot_product - np.sqrt(
                                                                                             refraction_radicand))[:,np.newaxis] * lens_normal_vectors
        # % This normalizes the ray propogation direction so that it's magnitude
        # % equals one
        ray_propogation_direction /= np.sqrt(
            ray_propogation_direction[:,0] ** 2 + ray_propogation_direction[:,1] ** 2 + ray_propogation_direction[:,
                2] ** 2)[:,np.newaxis]

        # % If the absorbance rate is non-zero, this calculates how much of the
        # % radiance is absorbed by the lens, otherwise the output radiance is
        # % just scaled by the transmission ratio
        if element_absorbance_rate:
            # % This calculates the distance that the rays traveled through
            # % the lens
            propogation_distance = np.sqrt((x_intersect - ray_source_coordinates[:,0]) ** 2 + (
            y_intersect - ray_source_coordinates[:,1]) ** 2 + (z_intersect - ray_source_coordinates[:,2]) ** 2)

            # % If the absorbance rate is a simple constant, this
            if isinstance(element_absorbance_rate,(int,long,float,complex)):
                # % This calculates the new light ray radiance values
                ray_radiance = (1.0 - element_absorbance_rate) * ray_radiance * propogation_distance

            else:
                # % These are the initial coordinates of the light ray passing
                # % through the lens in world coordinates
                x1 = ray_source_coordinates[:,0]
                y1 = ray_source_coordinates[:,1]
                z1 = ray_source_coordinates[:,2]

                # % These are the final coordinates of the light ray passing
                # % through the lens in world coordinates
                x2 = x_intersect
                y2 = y_intersect
                z2 = z_intersect

                # % This converts the initial light ray coordinates from the
                # % world coordinate system to the lens coordinate system
                [xL1, yL1, zL1] = convert_world_coordinates_to_lens_coordinates(x1, y1, z1, xc, yc, zc, a, b, c)

                # % This converts the final light ray coordinates from the
                # % world coordinate system to the lens coordinate system
                [xL2, yL2, zL2] = convert_world_coordinates_to_lens_coordinates(x2, y2, z2, xc, yc, zc, a, b, c)

                # % This calculates the radial lens coordinates for the initial
                # % light ray coordinates
                rL1 = np.sqrt(xL1 ** 2 + yL1 ** 2)
                # % This calculates the radial lens coordinates for the final
                # % light ray coordinates
                rL2 = np.sqrt(xL2 ** 2 + yL2 ** 2)

                # % This defines the independent variables that may be used in
                # % the definition of the absorbance function
                x_substitution = '(xL1+(xL2-xL1)*t)'
                y_substitution = '(yL1+(yL2-yL1)*t)'
                z_substitution = '(zL1+(zL2-zL1)*t)'
                r_substitution = '(rL1+(rL2-rL1)*t)'

                # % This replaces any instances of the string 'x' in the
                # % absorbance function string with the string 'x_substitution'
                element_absorbance_rate = element_absorbance_rate.replace('x', x_substitution)
                # % This replaces any instances of the string 'y' in the
                # % absorbance function string with the string 'y_substitution'
                element_absorbance_rate = element_absorbance_rate.replace('y', y_substitution)
                # % This replaces any instances of the string 'z' in the
                # % absorbance function string with the string 'z_substitution'
                element_absorbance_rate = element_absorbance_rate.replace('z', z_substitution)
                # % This replaces any instances of the string 'r' in the
                # % absorbance function string with the string 'r_substitution'
                element_absorbance_rate = element_absorbance_rate.replace('r', r_substitution)

                # % This defines the function handle to the absorbance function
                absorbance_function_handle = eval(
                    "lambda " + x_substitution + "," + y_substitution + "," + z_substitution + "," +
                    r_substitution + ": " + element_absorbance_rate)

                # % This calculates the integral of the absorbance function
                # % across the lens
                element_transmission_ratio = 1 - propogation_distance * siint.quad(absorbance_function_handle, 0.0, 1.0)

                # % This rescales the output radiance by the transmission ratio
                ray_radiance = element_transmission_ratio * ray_radiance

        else:
            # % This rescales the output radiance by the transmission ratio
            ray_radiance = element_transmission_ratio * ray_radiance

        # % This sets the new light ray origin to the intersection point with the
        # % front surface of the lens
        ray_source_coordinates = np.array([x_intersect, y_intersect, z_intersect])
        ray_source_coordinates = ray_source_coordinates.T
        # % Any rays that have complex values due to the radicand in the above
        # % equation being negative experience total internal reflection.  The
        # % values of these rays are set to NaN for now.
        tir_indices = np.less(refraction_radicand,0)
        #
        # % This sets any of the light ray directions experiencing total
        # % internal reflection to NaN values
        ray_propogation_direction[tir_indices,:] = np.NAN
        # % This sets any of the light ray origin coordinates experiencing total
        # % internal reflection to NaN values
        ray_source_coordinates[tir_indices,:] = np.NAN
        # % This sets any of the light ray wavelengths experiencing total
        # % internal reflection to NaN values
        ray_wavelength[tir_indices] = np.NAN
        # % This sets any of the light ray radiances experiencing total internal
        # % reflection to NaN values
        ray_radiance[tir_indices] = np.NAN

    # % If the element type is 'aperture', this extracts the optical properties
    # % required to perform the ray propogation
    elif element_type == 'aperture':

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # % Extraction of the aperture optical properties                       %
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        # % This extracts the pitch of the aperture stop
        element_pitch = optical_element['element_geometry']['pitch']

        # % This is the thickness of the optical element from the front optical axis
        # % vertex to the back optical axis vertex
        element_vertex_distance = optical_element['element_geometry']['vertex_distance']

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # % Propogation of the light rays through the aperture front surface    %
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        # % This extracts the individual plane parameters
        a = element_plane_parameters[0]
        b = element_plane_parameters[1]
        c = element_plane_parameters[2]
        d = element_plane_parameters[3]

        # % This calculates the square of the plane normal vector magnitude
        norm_vector_magnitude = np.sqrt(a ** 2 + b ** 2 + c ** 2)

        # % This is the current offset to the center of the front spherical
        # % surface
        ds = -element_vertex_distance / 2.0

        # % This calculates the transformed plane parameters (only d changes)
        d_temp = d - ds * norm_vector_magnitude

        # % This is the independent intersection time between the light rays and
        # % the first plane of the aperture stop
        intersection_time = -(
        a * ray_source_coordinates[:,0] + b * ray_source_coordinates[:,1] + c * ray_source_coordinates[:,
            2] + d_temp) / (a * ray_propogation_direction[:,0] + b * ray_propogation_direction[:,1] + c *
                            ray_propogation_direction[:,2])

        # % This calculates the intersection points
        x_intersect = ray_source_coordinates[:,0] + ray_propogation_direction[:,0] * intersection_time
        y_intersect = ray_source_coordinates[:,1] + ray_propogation_direction[:,1] * intersection_time
        z_intersect = ray_source_coordinates[:,2] + ray_propogation_direction[:,2] * intersection_time

        # % This calculates how far the intersection point is from the optical
        # % axis of the lens (for determining if the rays hit the lens outside of
        # % the pitch and thus would be destroyed)
        optical_axis_distance = measure_distance_to_optical_axis(x_intersect, y_intersect, z_intersect,
                                                                 element_center.T, element_plane_parameters.T)

        # % This gives the indices of the light rays that actually intersect the
        # % aperture stop (ie the rays that are within the aperture pitch)
        intersect_aperture_indices = (optical_axis_distance <= (element_pitch / 2))

        # % This sets any of the light ray directions outside of the domain of
        # % the lens to NaN values
        ray_propogation_direction[np.logical_not(intersect_aperture_indices),:] = np.NAN
        # % This sets any of the light ray origin coordinates outside of the
        # % domain of the lens to NaN values
        ray_source_coordinates[np.logical_not(intersect_aperture_indices),:] = np.NAN
        # % This sets any of the light ray wavelengths outside of the  domain of
        # % the lens to NaN values
        ray_wavelength[np.logical_not(intersect_aperture_indices)] = np.NAN
        # % This sets any of the light ray radiances outside of the  domain of
        # % the lens to NaN values
        ray_radiance[np.logical_not(intersect_aperture_indices)] = np.NAN

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # % Propogation of the light rays through the aperture back surface     %
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        # % This is the current offset to the center of the front spherical
        # % surface
        ds = +element_vertex_distance / 2

        # % This calculates the transformed plane parameters (only d changes)
        d_temp = d - ds * norm_vector_magnitude

        # % This is the independent intersection time between the light rays and
        # % the first plane of the aperture stop
        intersection_time = -(
        a * ray_source_coordinates[:,0] + b * ray_source_coordinates[:,1] + c * ray_source_coordinates[:,
            2] + d_temp) / (a * ray_propogation_direction[:,0] + b * ray_propogation_direction[:,1] + c *
                            ray_propogation_direction[:,2])

        # % This calculates the intersection points
        x_intersect = ray_source_coordinates[:,0] + ray_propogation_direction[:,0] * intersection_time
        y_intersect = ray_source_coordinates[:,1] + ray_propogation_direction[:,1] * intersection_time
        z_intersect = ray_source_coordinates[:,2] + ray_propogation_direction[:,2] * intersection_time

        # % This calculates how far the intersection point is from the optical
        # % axis of the lens (for determining if the rays hit the lens outside of
        # % the pitch and thus would be destroyed)
        optical_axis_distance = measure_distance_to_optical_axis(x_intersect, y_intersect, z_intersect,
                                                                 element_center.T, element_plane_parameters.T)

        # % This gives the indices of the light rays that actually intersect the
        # % aperture stop (ie the rays that are within the aperture pitch)
        intersect_aperture_indices = (optical_axis_distance <= (element_pitch / 2))

        # % This sets any of the light ray directions outside of the domain of
        # % the lens to NaN values
        ray_propogation_direction[np.logical_not(intersect_aperture_indices),:] = np.NAN
        # % This sets any of the light ray origin coordinates outside of the
        # % domain of the lens to NaN values
        ray_source_coordinates[np.logical_not(intersect_aperture_indices),:] = np.NAN
        # % This sets any of the light ray wavelengths outside of the  domain of
        # % the lens to NaN values
        ray_wavelength[np.logical_not(intersect_aperture_indices)] = np.NAN
        # % This sets any of the light ray radiances outside of the  domain of
        # % the lens to NaN values
        ray_radiance[np.logical_not(intersect_aperture_indices)] = np.NAN

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Saving the light ray propogation data                                   %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # % This extracts the propogation direction of the light rays
    light_ray_data['ray_propogation_direction'] = ray_propogation_direction
    #
    # % This extracts the light ray source coordinates
    light_ray_data['ray_source_coordinates'] = ray_source_coordinates
    #
    # % This extracts the wavelength of the light rays
    light_ray_data['ray_wavelength'] = ray_wavelength
    #
    # % This extracts the light ray radiance
    light_ray_data['ray_radiance'] = ray_radiance

    return light_ray_data


def propogate_rays_through_multiple_elements(optical_element, element_center, element_plane_parameters, light_ray_data):
    # % This function calculates the propogation of a set of light rays through a
    # % single optical element.

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Extraction of the light ray propogation data                            %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # % This extracts the propogation direction of the light rays
    ray_propogation_direction = light_ray_data['ray_propogation_direction']
    # % This extracts the light ray source coordinates
    ray_source_coordinates = light_ray_data['ray_source_coordinates']
    # % This extracts the wavelength of the light rays
    ray_wavelength = light_ray_data['ray_wavelength']
    # % This extracts the light ray radiance
    ray_radiance = light_ray_data['ray_radiance']

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Determination of optical element light ray intersection matching        %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # % This finds the unique set of planes defining the lens elements to
    # % intersect the light rays with
    unique_plane_element_parameters = np.vstack({tuple(row) for row in element_plane_parameters})
    # unique_plane_element_parameters = unique(element_plane_parameters,'rows') TODO: CHECK

    # % This is the number of unique plans that the light rays intersect
    unique_plane_number = unique_plane_element_parameters.shape[0]

    # % This initializes a vector of intersection times to calculate the order in
    # % which the light rays intersect the various optical elements
    unique_plane_intersection_time = np.zeros([ray_source_coordinates.shape[0], unique_plane_number])

    # % This initializes the light ray plane intersection point arrays
    x_intersect_approximate = np.zeros([ray_source_coordinates.shape[0], unique_plane_number])
    y_intersect_approximate = np.zeros([ray_source_coordinates.shape[0], unique_plane_number])
    z_intersect_approximate = np.zeros([ray_source_coordinates.shape[0], unique_plane_number])

    # % This iterates through the unique planes and calculates when the light
    # % rays hit each plane since the light rays must be sequentially propogated
    # % through the optical elements, but the ordering of the elements may be a
    # % bit more poorly defined for multiple elements (although in practice, this
    # % probably won't happen very often)
    for plane_parameters_index in range(0, unique_plane_number):
        # % This extracts the individual plane parameters
        a = unique_plane_element_parameters[plane_parameters_index][0]
        b = unique_plane_element_parameters[plane_parameters_index][1]
        c = unique_plane_element_parameters[plane_parameters_index][2]
        d = unique_plane_element_parameters[plane_parameters_index][3]

        # % This is the independent intersection time between the light rays and
        # % the current plane of optical elements
        unique_plane_intersection_time[:,plane_parameters_index] = -(
        a * ray_source_coordinates[:,0] + b * ray_source_coordinates[:,1] + c * ray_source_coordinates[:,2] + d) / (
                                                                    a * ray_propogation_direction[:,0] + b *
                                                                    ray_propogation_direction[:,1] + c *
                                                                    ray_propogation_direction[:,2])

        # % This calculates the intersection points of the light rays with the
        # % optical element planes
        x_intersect_approximate[:,plane_parameters_index] = ray_source_coordinates[:,0] + \
                                                             ray_propogation_direction[:,0] * \
                                                             unique_plane_intersection_time[plane_parameters_index]
        y_intersect_approximate[:,plane_parameters_index] = ray_source_coordinates[:,1] + \
                                                             ray_propogation_direction[:,1] * \
                                                             unique_plane_intersection_time[plane_parameters_index]
        z_intersect_approximate[:,plane_parameters_index] = ray_source_coordinates[:,2] + \
                                                             ray_propogation_direction[:,2] * \
                                                             unique_plane_intersection_time[plane_parameters_index]

    # % This sorts the plane intersection times to determine which planes should
    # % be tested first
    unique_plane_index_order = np.argsort(unique_plane_intersection_time, 0)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Light ray propogation through multiple elements                         %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # % This initializes the structure to contain the light ray data for
    # % the current optical element
    light_ray_data_temp = {}

    # % This iterates through the different planes, calculating the location of
    # % the optical elements that the light rays intersect
    for plane_parameters_index in range(0, unique_plane_number):

        # % These are the current intersection points
        x_intersect_approxiate_current = x_intersect_approximate[:,unique_plane_index_order[plane_parameters_index]]
        y_intersect_approxiate_current = y_intersect_approximate[:,unique_plane_index_order[plane_parameters_index]]
        z_intersect_approxiate_current = z_intersect_approximate[:,unique_plane_index_order[plane_parameters_index]]

        # % This is a vector giving the indices of the current optical elements
        optical_element_indices = np.argwhere(
            unique_plane_element_parameters[plane_parameters_index][:] == element_plane_parameters)

        # % These are the centers of the optical elements on the current plane
        xc_current = element_center[optical_element_indices][0]
        yc_current = element_center[optical_element_indices][1]
        zc_current = element_center[optical_element_indices][2]

        # % This finds the closest lens center to each approximate intersection
        # % point (there are a small number of cases where this may give an
        # % incorrect match since the intersection points used here are based
        # % upon the center points of the lenses and not the front surfaces, but
        # % this should work well for the large majority of possible systems)
        nbrs = NearestNeighbors(n_neighbors=1, metric='euclidean').fit([xc_current, yc_current, zc_current])
        distances, indices = nbrs.kneighbors(
            [x_intersect_approxiate_current, y_intersect_approxiate_current, z_intersect_approxiate_current])

        nearest_neighbor_index = np.squeeze(indices)
        # nearest_neighbor_index=knnsearch([xc_current,yc_current,zc_current],[x_intersect_approxiate_current,y_intersect_approxiate_current,z_intersect_approxiate_current],'K',1,'Distance','euclidean');

        # % This iterates through the individual optical elements extracting the
        # % set of light rays that likely intersect the element and propogating
        # % them through the element
        for element_index in range(0, optical_element_indices.size):
            # % This is the set of indices into the light rays to propogate
            # % through the current optical element
            light_ray_indices = (nearest_neighbor_index == element_index)

            # % This extracts the propogation direction of the light rays
            light_ray_data_temp['ray_propogation_direction'] = ray_propogation_direction[light_ray_indices][:]
            # % This extracts the light ray source coordinates
            light_ray_data_temp['ray_source_coordinates'] = ray_source_coordinates[light_ray_indices][:]
            # % This extracts the wavelength of the light rays
            light_ray_data_temp['ray_wavelength'] = ray_wavelength[light_ray_indices]
            # % This extracts the light ray radiance
            light_ray_data_temp['ray_radiance'] = ray_radiance[light_ray_indices]

            # % This extracts the current optical element data
            current_optical_element = optical_element[optical_element_indices[element_index]]
            # % This extracts the current optical element plane parameters
            current_plane_parameters = element_plane_parameters[optical_element_indices[element_index]][:]
            # % This extracts the current center of the optical element
            current_element_center = element_center[optical_element_indices[element_index]][:]

            # % This propogates the light rays through the single optical element
            light_ray_data_temp = propogate_rays_through_single_element(current_optical_element, current_element_center,
                                                                        current_plane_parameters, light_ray_data_temp)

            # % This extracts the propogation direction of the light rays
            ray_propogation_direction[light_ray_indices][:] = light_ray_data_temp['ray_propogation_direction']
            # % This extracts the light ray source coordinates
            ray_source_coordinates[light_ray_indices][:] = light_ray_data_temp['ray_source_coordinates']
            # % This extracts the wavelength of the light rays
            ray_wavelength[light_ray_indices] = light_ray_data_temp['ray_wavelength']
            # % This extracts the light ray radiance
            ray_radiance[light_ray_indices] = light_ray_data_temp['ray_radiance']

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Saving the light ray propogation data                                   %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # % This extracts the propogation direction of the light rays
    light_ray_data['ray_propogation_direction'] = ray_propogation_direction
    # % This extracts the light ray source coordinates
    light_ray_data['ray_source_coordinates'] = ray_source_coordinates
    # % This extracts the wavelength of the light rays
    light_ray_data['ray_wavelength'] = ray_wavelength
    # % This extracts the light ray radiance
    light_ray_data['ray_radiance'] = ray_radiance

    return light_ray_data


def propogate_rays_through_optical_system(element_data, element_center, element_plane_parameters, element_system_index,
                                          light_ray_data):
    # % This function propogates the light ray data defined by the structure
    # % 'light_ray_data' through the optical system defined by the input
    # % arguments.

    # % This is the number of sequential optical elements within the total
    # % optical system that the light rays must be interatively passed through
    sequential_element_number = np.max(element_system_index)

    # % Since the the optical is defined from the sensor moving outward (ie the
    # % opposite direction in which the light will enter a camera system), this
    # % reverses the indexin of the optical elements so that the first element
    # % index corresponds to the first optical element that the light will hit
    element_system_index = np.squeeze(sequential_element_number - element_system_index) # subtracted 1

    # % This iterates through the sequential optical elements propogating the
    # % light rays through each successive element (or system of coplanar
    # % elements)
    for element_index in range(0, sequential_element_number):

        # % These are the indices of the current element or elements to propogate
        # % the light rays through
        current_element_indices = np.squeeze(np.argwhere(element_system_index == element_index))

        # % This is the number of elements that the light rays need to be
        # % simultaneously propogated through
        simultaneous_element_number = current_element_indices.size

        # % If there is only a single element that the rays are to be propogated
        # % through, this propogates the rays through the single element;
        # % otherwise the light rays are simultaneously propogated through the
        # % multiple elements
        if simultaneous_element_number == 1:

            # % This extracts the current optical element data
            current_optical_element = element_data[current_element_indices]
            # % This extracts the current optical element plane parameters
            current_plane_parameters = element_plane_parameters[current_element_indices,:]
            # % This extracts the current center of the optical element
            current_element_center = element_center[current_element_indices,:]

            # % This propogates the light rays through the single optical element
            light_ray_data = propogate_rays_through_single_element(current_optical_element, current_element_center,
                                                                   current_plane_parameters, light_ray_data);

        else:

            # % This initializes a cell array to contain the optical element data
            current_optical_element = np.zeros((1, 1), dtype=np.object)
            # % This iterates through the individual optical elements extracting
            # % the optical element data
            for simultaneous_element_index in range(1, simultaneous_element_number + 1):
                # % This extracts the current optical element data
                current_optical_element[simultaneous_element_index - 1][:] = element_data[
                    current_element_indices[simultaneous_element_index - 1]]
            # % This extracts the current optical element plane parameters
            current_plane_parameters = element_plane_parameters[current_element_indices][:]
            # % This extracts the current center of the optical element
            current_element_center = element_center[current_element_indices][:]

            # % This propogates the light rays through the multiple optical
            # % elements
            light_ray_data = propogate_rays_through_multiple_elements(current_optical_element, current_element_center,
                                                                      current_plane_parameters, light_ray_data)

    return light_ray_data


def intersect_sensor_better(camera_design, x, y):
    # % This function determines the indicies of the pixel on the sensor at which
    # % the ray intersects.

    # % This is the pixel pitch [micron]
    pixel_pitch = camera_design['pixel_pitch']
    # % This is the number of pixels in the x-direction
    x_pixel_number = camera_design['x_pixel_number']
    # % This is the number of pixels in the y-direction
    y_pixel_number = camera_design['y_pixel_number']
    #
    # % This is the coordinate of pixel (1,1) [0][0]
    pixel_1_x = -pixel_pitch * (x_pixel_number - 1) / 2.0
    pixel_1_y = -pixel_pitch * (y_pixel_number - 1) / 2.0

    # % This is the number of pixel diameters the point (x,y) is from the center
    # % of the (0,0) pixel
    d_x = (x - pixel_1_x) / pixel_pitch + 1.5
    d_y = (y - pixel_1_y) / pixel_pitch + 1.5

    # % These are the coordinates of the point that is 1/2 pixel less than the
    # % center coordinate
    d_y_lower = d_y - 0.5
    d_x_lower = d_x - 0.5

    # % These are the percentages of overlap for the upper right corner of the
    # % pixel (actually this is the distance of overlap - but if the pixels have
    # % dimensions of 1 x 1 then this is the same) in each direction
    d_ii_ul = np.ceil(d_y_lower) - d_y_lower
    d_jj_ul = np.ceil(d_x_lower) - d_x_lower
    # % This is the area of overlap with the upper right pixel
    w_ul = (d_ii_ul) * (d_jj_ul)

    # % These are the percentages of overlap for the upper left corner of the
    # % pixel (actually this is the distance of overlap - but if the pixels have
    # % dimensions of 1 x 1 then this is the same) in each direction
    d_ii_ur = np.ceil(d_y_lower) - d_y_lower
    d_jj_ur = 1 - d_jj_ul
    # % This is the area of overlap with the upper left pixel
    w_ur = (d_ii_ur) * (d_jj_ur)

    # % These are the percentages of overlap for the lower left corner of the
    # % pixel (actually this is the distance of overlap - but if the pixels have
    # % dimensions of 1 x 1 then this is the same) in each direction
    d_ii_ll = 1 - d_ii_ul
    d_jj_ll = np.ceil(d_x_lower) - d_x_lower
    # % This is the area of overlap with the lower right pixel
    w_ll = (d_ii_ll) * (d_jj_ll)
    #
    # % These are the percentages of overlap for the lower right corner of the
    # % pixel (actually this is the distance of overlap - but if the pixels have
    # % dimensions of 1 x 1 then this is the same) in each direction
    d_ii_lr = 1 - d_ii_ul
    d_jj_lr = 1 - d_jj_ul
    # % This is the area of overlap with the lower left pixel
    w_lr = (d_ii_lr) * (d_jj_lr)

    # % These are the integral values of the upper left pixel
    ii_ul = np.ceil(d_y_lower) - 1
    jj_ul = np.ceil(d_x_lower) - 1
    # % These are the integral values of the upper right pixel
    ii_ur = np.ceil(d_y_lower) - 1
    jj_ur = jj_ul + 1
    # % These are the integral values of the lower left pixel
    ii_ll = ii_ul + 1
    jj_ll = np.ceil(d_x_lower) - 1
    # % These are the integral values of the lower right pixel
    ii_lr = ii_ul + 1
    jj_lr = jj_ul + 1

    # % This is a vector of the ii coordinates
    ii_indices = np.array([ii_ul, ii_ur, ii_ll, ii_lr]).T
    # % This is a vector of the jj coordinates
    jj_indices = np.array([jj_ul, jj_ur, jj_ll, jj_lr]).T
    # % This is a vector of pixel weights
    pixel_weights = np.array([w_ul, w_ur, w_ll, w_lr]).T

    # % This sets the weights of the jj indices that are out of bounds to NaN
    # % values
    for r in range(0,pixel_weights.shape[0]):
        for c in range(0,pixel_weights.shape[1]):
            if(jj_indices[r,c]<0) or (x_pixel_number<=jj_indices[r,c]):
                pixel_weights[r,c] = np.NAN
            if(ii_indices[r,c]<0) or  (y_pixel_number<=ii_indices[r,c]):
                pixel_weights[r,c] = np.NAN
    # pixel_weights[(jj_indices < 1) or (x_pixel_number < jj_indices)] = np.NAN
    # % This sets the weights of the ii indices that are out of bounds to NaN
    # % values
    # pixel_weights[(ii_indices < 1) or (y_pixel_number < ii_indices)] = np.NAN

    return [ii_indices, jj_indices, pixel_weights]

def trace_rays_through_density_gradients(light_ray_data):
    # This function calls the SchlierenRay code to trace rays through a medium with
    # density gradients

    # specify file names to save lightray position and direction
    light_ray_pos_filename = 'lightrayPos_i.bin'
    light_ray_dir_filename = 'lightrayDir_i.bin'

    # change working directory to where schlieren is located
    #subprocess.call('cd ../ray_tracing_density_gradients/schlieren-0.2.0-Build/', shell=True)
    os.chdir('/home/barracuda/a/lrajendr/Projects/ray_tracing_density_gradients/schlieren-0.2.0-Build/')    

    # save light ray data to file
    light_ray_data['ray_source_coordinates'].astype('float32').tofile(light_ray_pos_filename)
    light_ray_data['ray_propogation_direction'].astype('float32').tofile(light_ray_dir_filename)

    # call CUDA program to trace rays through density gradients
    subprocess.call('cuda-memcheck ./schlieren ' + 'const_grad.nrrd' + ' ' + light_ray_pos_filename + ' ' + light_ray_dir_filename
                    + ' '+ str(light_ray_data['ray_source_coordinates'].shape[0]), shell=True)

    # read updated position and direction data from file
    temp_pos_array = np.fromfile('lightrayPos_f.bin',dtype='float32')
    temp_dir_array = np.fromfile('lightrayDir_f.bin',dtype='float32')

    # reshape arrays
    light_ray_data['ray_source_coordinates'] = np.reshape(temp_pos_array,(temp_pos_array.size/3,3))
    light_ray_data['ray_propogation_direction'] = np.reshape(temp_dir_array,(temp_dir_array.size/3,3))

    # move back to directory where the python codes are located
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    print "current directory : " + os.getcwd()

    return light_ray_data

def perform_ray_tracing_03(piv_simulation_parameters, optical_system, pixel_gain, scattering_data, scattering_type,
                           lightfield_source, field_type):
    # % This function calculates the ray tracing of the input lightrays
    # % generating the output image 'I'.
    # %
    # % This version is similar to version 01 except that it generates the rays
    # % in a more efficient manner so that fewer rays do not intersect the lens.
    # % This version does not contain the parallelization features offered by
    # % version 02.

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Extracts parameters from 'piv_simulation_parameters'                    %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # % This extracts the object distance of the lens (ie the distance between
    # % the lens front principal plane and the center of the focal plane) to the
    # % design structure (in microns)
    object_distance = piv_simulation_parameters['lens_design']['object_distance']
    # % This extracts the lens focal length from the design structure (in microns)
    focal_length = piv_simulation_parameters['lens_design']['focal_length']
    # % This extracts the lens f/# from the design structure
    aperture_f_number = piv_simulation_parameters['lens_design']['aperture_f_number']
    # % This is the number of pixels in the x-direction
    x_pixel_number = piv_simulation_parameters['camera_design']['x_pixel_number']
    # % This is the number of pixels in the y-direction
    y_pixel_number = piv_simulation_parameters['camera_design']['y_pixel_number']
    # % This is the bit depth of the camera sensor (which must be an integer
    # % less then or equal to 16)
    pixel_bit_depth = piv_simulation_parameters['camera_design']['pixel_bit_depth']
    # % This is the wavelength of the laser used for illumination of the
    # % particles (in microns)
    beam_wavelength = piv_simulation_parameters['particle_field']['beam_wavelength']
    #
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Extracts parameters from 'optical_system'                               %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # % This extracts the lens refractive index
    refractive_index = optical_system['design']['optical_element'][0]['optical_element'][0]['element_properties'][
        'refractive_index']
    # % This extracts the radius of curvature for the front surface of the lens
    front_surface_radius = optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry'][
        'front_surface_radius']
    # % This extracts the front surface to back surface vertex distance of the
    # % lens
    optical_system_length = optical_system['design']['optical_element'][0]['optical_element'][0]['element_geometry'][
        'vertex_distance']
    #
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Extracts parameter from 'lightfield_source'                             %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # % This extracts the number of lightrays per particle to the
    # % 'lightfield_source' data
    lightray_number_per_particle = lightfield_source['lightray_number_per_particle']
    # % This extracts the number of lightrays to simulateously process to the
    # % 'lightfield_source' data
    lightray_process_number = lightfield_source['lightray_process_number']
    #
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Calculates optical system properties                                    %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # % This calculates the image distance
    image_distance = (1.0 / focal_length - 1.0 / object_distance) ** (-1.0)
    #
    # % This calculates the offsets of the principal planes
    h2_principal_plane = -(focal_length * (refractive_index - 1.0) * optical_system_length) / (
        front_surface_radius * refractive_index)
    #
    # % This calculates the positions of the lens vertex planes
    v2_vertex_plane = image_distance + h2_principal_plane
    v1_vertex_plane = v2_vertex_plane + optical_system_length
    #
    # % This is the Z position of the sensor
    z_sensor = 0.0
    #
    # % This is the Z position of the lens center
    z_lens = (v1_vertex_plane + v2_vertex_plane) / 2.0
    #
    # % This calculates the pitch of the lens
    lens_pitch = focal_length / aperture_f_number
    #
    # % This initializes the element type array to Null so that the recursive function can do
    # % the magics
    elements_coplanar = []
    # % This initializes the elements coplanar variable to Null so that the recursive function
    # % can do some other magics
    element_type = []
    # % This initializes a counting variable to index the optical element data
    # % cell array
    element_count = 0
    # % This initializes the optical element data cell array
    # element_data = cell(1,1);
    element_data = np.zeros((1, 1), dtype=np.object)
    #
    # % This calculates the coordinates of the circles and planes that define the
    # % optical elements in the current system
    # [element_center,~,element_plane_parameters,~,~,~,element_system_index,~,element_data]=create_element_coordinate_arrays(optical_system.design,elements_coplanar,element_type,element_count,element_data);
    # %   [element_center,~,element_plane_parameters,~,~,~,element_system_index,~,element_data]=create_element_coordinate_arrays(optical_system.design,elements_coplanar,element_type,element_count,element_data);
    [element_center, _, element_plane_parameters, _, _, _, element_system_index, _,
     element_data] = create_element_coordinate_arrays(optical_system['design'],
                                                      elements_coplanar,
                                                      element_type, element_count,
                                                      element_data)

    # remove first element because in the cell for element_data. because i don't understand python syntax well enough
    # to know why it's added in the first place.
    element_data = element_data[1:]
    # [element_center, element_pitch, element_plane_parameters, element_type, element_thickness, total_element_distance,
    #  element_system_index, element_count, element_data] = create_element_coordinate_arrays(optical_system['design'],
    #                                                                                        elements_coplanar,
    #                                                                                        element_type, element_count,
    #                                                                                        element_data)

    element_center = np.squeeze(element_center)
    # element_pitch = np.squeeze(element_pitch)
    element_plane_parameters = np.squeeze(element_plane_parameters)
    if(len(element_plane_parameters.shape)==1):
        element_plane_parameters = np.reshape(element_plane_parameters,(1,element_plane_parameters.shape[0]))
    if(len(element_center.shape)==1):
        element_center = np.reshape(element_center,(1,element_center.shape[0]))

    # element_type = np.squeeze(element_type)
    # element_thickness = np.squeeze(element_thickness)

    # % This offsets the lens to account for the sensor position
    element_plane_parameters[:,3] = element_plane_parameters[:,3] - element_plane_parameters[:,2] * z_lens
    element_center[:,2] = element_center[:,2] + z_lens

    if (field_type == 'particle'):
        savefile = 'mat_files/start_particle.mat'
    else:
        savefile = 'mat_files/start_calibration.mat'

    sio.savemat(savefile, {'image_distance': image_distance,
                           'h2_principal_plane': h2_principal_plane,
                           'v2_vertex_plane': v2_vertex_plane,
                           'v1_vertex_plane': v1_vertex_plane,
                           'z_sensor': z_sensor,
                           'z_lens': z_lens,
                           'lens_pitch': lens_pitch,
                           'element_center': element_center,
                           'element_plane_parameters': element_plane_parameters,
                           'element_type': element_type,
                           'element_system_index': element_system_index,
                           'element_count': element_count,
                           'element_data': element_data,
                           'elements_coplanar': elements_coplanar}, long_field_names=True)

    # I = 1
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Begins raytracing operations                                            %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # % This initializes the sensor image
    I = np.zeros([x_pixel_number, y_pixel_number])

    # # TODO change this
    lightray_process_number = 10
    lightray_number_per_particle = 10

    # % This generates an array of indices into the source points to calculate the lightfield
    lightfield_N = lightfield_source['x'].size / np.ceil(lightray_process_number*1.0 / lightray_number_per_particle)
    if(lightfield_N<1.0):
        # lightfield_vector = np.linspace(0,lightfield_source['x'].size,endpoint=False)
        lightfield_vector = np.r_[0:lightfield_source['x'].size-1]
    else:
        # lightfield_vector = np.linspace(0, lightfield_source['x'].size, lightfield_N,endpoint=False)
        lightfield_vector = np.r_[0:lightfield_source['x'].size-1:lightfield_N]
    # % This checks whether the last segment of indices is at the end of the
    # % vector and if not, adds them
    if lightfield_vector[-1] != lightfield_source['x'].size-1:
        lightfield_vector = np.append(lightfield_vector, lightfield_source['x'].size-1)

    light_ray_data = {}
    # % This iterates through the rays in blocks of less then or equal to
    # % lightray_process_number in size
    pbar = ProgressBar(maxval=len(lightfield_vector)-2).start()
    for m in range(0, len(lightfield_vector) - 1):
        # % This displays the progress of the sensor rendering
        pbar.update(m)
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # % Generate lightfield and propogate it through the optical system     %
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        # % This generates the lightfield data for the current set of source points
        lightfield_data = generate_lightfield_angular_data(lens_pitch, image_distance, scattering_data, scattering_type,
                                                           lightfield_source, lightray_number_per_particle,
                                                           int(lightfield_vector[m]), int(lightfield_vector[m + 1]))
        # % This extracts the light ray source coordinates
        light_ray_data['ray_source_coordinates'] = np.array(
            [lightfield_data['x'], lightfield_data['y'], lightfield_data['z']]).T
        # % This extracts the propogation direction of the light rays
        light_ray_data['ray_propogation_direction'] = np.array(
            [lightfield_data['theta'], lightfield_data['phi'], -np.ones(lightfield_data['theta'].size)]).T
        # % This ensures that the light ray propogation direction vector has a
        # % unit magnitude
        light_ray_data['ray_propogation_direction'] /= la.norm(light_ray_data['ray_propogation_direction'], axis=1)[:,
                                                       None]
        # % This extracts the wavelength of the light rays
        light_ray_data['ray_wavelength'] = beam_wavelength * np.ones(lightfield_data['theta'].size)
        # % This extracts the light ray radiance
        # %light_ray_data.ray_radiance=ones(size(lightfield_data.theta))';
        light_ray_data['ray_radiance'] = (1.0 / aperture_f_number ** 2) * np.transpose(lightfield_data['radiance'])

        light_ray_data = trace_rays_through_density_gradients(light_ray_data)


        # % This propogates some imaginary light rays through the optical system
        light_ray_data = propogate_rays_through_optical_system(element_data, element_center, element_plane_parameters,
                                                               element_system_index, light_ray_data)

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # % Propogation of the light rays to the sensor                         %
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        # % This extracts the propogation direction of the light rays
        ray_propogation_direction = light_ray_data['ray_propogation_direction']
        # % This extracts the light ray source coordinates
        ray_source_coordinates = light_ray_data['ray_source_coordinates']

        # % This extracts the individual plane parameters for the sensor
        a = 0.0
        b = 0.0
        c = 1.0
        d = -z_sensor

        # % This is the independent intersection time between the light rays and
        # % the first plane of the aperture stop
        intersection_time = -(
        a * ray_source_coordinates[:,0] + b * ray_source_coordinates[:,1] + c * ray_source_coordinates[:,2] + d) / (
                            a * ray_propogation_direction[:,0] + b * ray_propogation_direction[:,1] + c *
                            ray_propogation_direction[:,2])

        # % This calculates the intersection points
        x_intersect = ray_source_coordinates[:,0] + ray_propogation_direction[:,0] * intersection_time
        y_intersect = ray_source_coordinates[:,1] + ray_propogation_direction[:,1] * intersection_time
        z_intersect = ray_source_coordinates[:,2] + ray_propogation_direction[:,2] * intersection_time

        # % This sets the new light ray origin to the intersection point with the
        # % front surface of the lens
        ray_source_coordinates = np.array([x_intersect, y_intersect, z_intersect])
        ray_source_coordinates = ray_source_coordinates.T
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # % Add lightray radiance to sensor integration                         %
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        # % This is the angle between the lightray and the sensor (ie a lightray
        # % normal to the sensor would yield an angle of zero)
        alpha = np.arctan(np.sqrt((ray_propogation_direction[:,0] / ray_propogation_direction[:,2]) ** 2 + (
        ray_propogation_direction[:,1] / ray_propogation_direction[:,2]) ** 2))
        # % This calculates the cos^4(alpha) term which controls the contribution
        # % of the incident light rays onto the measured energy in the sensor
        cos_4_alpha = np.cos(alpha) ** 4

        # % This calculates the indicies of the pixel on the sensor that the ray
        # % intersects and the relative weighting between the pixels
        [ii_indices, jj_indices, pixel_weights] = intersect_sensor_better(piv_simulation_parameters['camera_design'],
                                                                          ray_source_coordinates[:,0],
                                                                          ray_source_coordinates[:,1])

        # % These nested for loops complete the same calculation as the above
        # % block of code, but require about 2.5 times longer to complete
        # %
        # % This iterates through the pixel locations upon which a ray intersected
        # % incrementing the image intensities
        for n in range(0, ray_source_coordinates.shape[0]):
            # % This iterates through the four adjacent pixels to the
            # % intersection point of the light ray
            for p in range(0, 4):
                # % This adds the ray to the upper left intersection pixel (if the
                # % weighting is not a NaN value)
                if not(np.isnan(pixel_weights[n][p])):
                    # % This interpolates the ray intensities onto the pixels
                    # % surrounding the pixel intersection point
                    I[ii_indices[n][p]][jj_indices[n][p]] += pixel_weights[n][p] * light_ray_data['ray_radiance'][n] * \
                                                             cos_4_alpha[n]

    pbar.finish()
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # % Rescales and resamples image for export                                 %
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    '''
    # % This rescales the image intensity to account for the pixel gain
    I *= 10 ** (pixel_gain/ 20.0)
    #
    # % This rescales the image to values ranging from 0 to 2^bit_depth-1
    I = (2 ** pixel_bit_depth - 1) * I / (2 ** 16 - 1)
    # % This rounds the values of the image to the nearest integer
    I = np.round(I)
    # % This rescales the image to the full 16 bit range (but only bit_depth bits
    # % of information are now used)
    I *= (2 ** 16 - 1.0) / (2 ** pixel_bit_depth - 1.0)

    # % This converts the image from double precision to 16 bit precision
    I = np.uint16(I)

    #plt.imshow(I,cmap = plt.get_cmap('gray'),vmin=0,vmax=2**16-1)
    #plt.show()
    '''
    return I
