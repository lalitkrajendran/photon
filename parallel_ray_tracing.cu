/*
 * parallel_ray_tracing.cu
 *
 *  Created on: Apr 20, 2016
 *      Author: lrajendr
 */
#include <stdio.h>
#include <fstream>
#include <string>
#include "parallel_ray_tracing.h"
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <string.h>
#include <vector_types.h>
#include <vector_functions.h>
#include "float3_operators.h"
#include <iostream>
#include <stdlib.h>
#include <numeric>
using namespace std;
//#incldue <cutil_math.h>

cudaArray *data_array = 0;
texture<float, 2> mie_scattering_irradiance;

__device__ float random_single(unsigned int seed)
{

  /* CUDA's random number library uses curandState_t to keep track of the seed value
     we will store a random state for every thread  */
  curandState_t state;

  /* the seed can be the same for each core, here we pass the time in from the CPU */
  /* the sequence number should be different for each core (unless you want all
                             cores to get the same sequence of numbers for some reason - use thread id! */
  /* the offset is how much extra we advance in the sequence for each call, can be 0 */

  /* we have to initialize the state */
  curand_init(seed, blockIdx.x, 0, &state);

  float rand_num = curand_uniform(&state);
  return rand_num;
}

__global__ void generate_lightfield_angular_data(float lens_pitch, float image_distance,
		scattering_data_t scattering_data, int scattering_type, lightfield_source_t lightfield_source,
                                     int lightray_number_per_particle, int n_min, int n_max, float beam_wavelength,
                                     float aperture_f_number, light_ray_data_t* light_ray_data)
{
	/*
		This function generates the light field data for the source points specified by the
		structure lightfield_source.  The data is only generated for the source points from
		n_min to n_max.  The parameter lightray_number_per_particle is the number of rays to generate for each
		source point.
	*/

	//--------------------------------------------------------------------------------------
	// compute indices to access in lightfield_source and lightfield_data
	//--------------------------------------------------------------------------------------

	// find global thread ID
	int local_thread_id = threadIdx.x;

	float del_scattering_angle = (scattering_data.scattering_angle[1] - scattering_data.scattering_angle[0])*180.0/M_PI;
	float min_scattering_angle = scattering_data.scattering_angle[0];

	// get id of particle which is the source of light rays
	int particle_id = blockIdx.x*blockDim.x + local_thread_id;

	// get id of ray emitted by the particle
	int local_ray_id = blockIdx.z;
	int global_ray_id = local_ray_id + particle_id*lightray_number_per_particle;

	if(global_ray_id >= lightfield_source.num_particles*lightray_number_per_particle)
		return;

	// get source coordinates of the light ray
	float x_current = lightfield_source.x[particle_id];;
	float y_current = lightfield_source.y[particle_id];;
	float z_current = lightfield_source.z[particle_id];;

	//--------------------------------------------------------------------------------------
	// compute direction of propagation of light ray
	//--------------------------------------------------------------------------------------

	// generate random points on the lens
	float random_number_1 = random_single(particle_id * local_ray_id * global_ray_id);
	float random_number_2 = random_single(particle_id + local_ray_id + global_ray_id);
	float x_lens = 0.5*lens_pitch*random_number_1*cos(2*M_PI*random_number_2);
	float y_lens = 0.5*lens_pitch*random_number_1*sin(2*M_PI*random_number_2);

	// calculate the x angles for the light rays
	float theta_temp = -(x_lens - x_current) / (image_distance - z_current);
	// calculate the y angles for the light rays
	float phi_temp = -(y_lens - y_current) / (image_distance - z_current);



	//--------------------------------------------------------------------------------------
	// compute irradiance of the light ray
	//--------------------------------------------------------------------------------------

	int diameter_index;
	float3 ray_direction_vector, temp_vector;
	float dot_vector[3];
	float ray_scattering_angles, ray_scattering_irradiance;
	double irradiance_current;
	float3 beam_propagation_vector;
	// if scattering_type is mie, then use mie scattering data
	if(scattering_type)
	{
		//% This extracts the normalized beam propagation direction vector from the
		//% parameters structure
		beam_propagation_vector.x=scattering_data.beam_propagation_vector[0];
		beam_propagation_vector.y=scattering_data.beam_propagation_vector[1];
		beam_propagation_vector.z=scattering_data.beam_propagation_vector[2];

		// % This extracts the current particle diameter index
		diameter_index = lightfield_source.diameter_index[particle_id];
		// % This calculates the light ray's direction vectors (in the camera
		// % coordinate system)
		ray_direction_vector = make_float3(x_lens-x_current,y_lens-y_current,image_distance-z_current);

		// % This normalizes the ray direction vectors
		ray_direction_vector = normalize(ray_direction_vector);

		// % This rotates the light rays direction vectors by the inverse of the
        // % camera rotation array so that the ray is now in the world coordinate
        // % system
		for(int i = 0; i < 3; i++)
		{
			temp_vector.x = scattering_data.inverse_rotation_matrix[i*3 + 0];
			temp_vector.y = scattering_data.inverse_rotation_matrix[i*3 + 1];
			temp_vector.z = scattering_data.inverse_rotation_matrix[i*3 + 2];

			dot_vector[i] = dot(temp_vector,ray_direction_vector);
		}
		ray_direction_vector = make_float3(dot_vector[0],dot_vector[1],dot_vector[2]);

		// % This calculates the angle that the light ray direction vectors make
		// % with the laser propagation direction vector in radians
		ray_scattering_angles = angleBetween(beam_propagation_vector,ray_direction_vector);
		// % This calculates the Mie scattering irradiance at the current
		// % scattered angle and with the current particle diameter
		int lookup_angle = (int)(ray_scattering_angles - scattering_data.scattering_angle[0])/del_scattering_angle;
		ray_scattering_irradiance = scattering_data.scattering_irradiance[lookup_angle*27 + diameter_index];

//		ray_scattering_irradiance = tex2D(mie_scattering_irradiance,diameter_index,lookup_angle);
		// % This calculates the total irradiance for the current particle's rays
		irradiance_current=ray_scattering_irradiance*lightfield_source.radiance[particle_id];
	}
	// if not mie scattering, then set irradiance to be uniform
	else
	{
		// % This specifies the total irradiance for the current particle's
		// % rays to be uniform
		irradiance_current = lightfield_source.radiance[particle_id];
	}

	// save the light rays to the light field data structure
	light_ray_data[global_ray_id].ray_source_coordinates = make_float3(x_current,y_current,z_current);
	light_ray_data[global_ray_id].ray_propagation_direction = normalize(make_float3(theta_temp,phi_temp,1.0));
	light_ray_data[global_ray_id].ray_wavelength = beam_wavelength;
	light_ray_data[global_ray_id].ray_radiance = 1/(aperture_f_number*aperture_f_number)*irradiance_current;

}

__device__ float3 ray_sphere_intersection(float3 pos_c, float R, float3 dir_i, float3 pos_i, char surface)
{
	// % This function returns the point at which a ray first intersects a sphere
    // % in 3-space.  The sphere is defined by it's center coordinates given by
    // % the argumnets 'xc', 'yc', and 'zc' and it's radius given by the argument
    // % 'R'.  Thus the equation for the surface of the sphere is
    // %
    // %   (x-xc)^2+(y-yc)^2+(z-zc)^2=R^2
    // %
    // % The ray is defined by the direction vector given by the arguments 'vx',
    // % 'vy', and 'vz' and it's origin coordinate given by the arguments 'x0',
    // % 'y0', and 'z0'.  Thus the parametric equations defining the ray are
    // %
    // %   x(t)=x0+vx*t
    // %   y(t)=y0+vy*t
    // %   z(t)=z0+vz*t
    // %
    // % The string 'surface' may equal either 'front' or 'back' and states
    // % whether the light ray is entering the front of the lens or exiting the
    // % back of the lens.
    // %
    // % The the output arguments 'xi', 'yi', and 'zi' give the coordinate of the
    // % first intersection or the ray with the sphere (ie there will be two
    // % intersections and the one output is the one that is closest to the origin
    // % of the ray that is also positive).
    // %
    // % The input arguments can be of any size, but they must all have the same
    // % dimensions.  The output arguments will be the same size as the input
    // % argument.  If a particular ray does not intersection a particular sphere,
    // % then NaN values will be returned in the output arguments.

    // % This defines the coefficients of the quadratic polynomial
    float alpha = dot(dir_i,dir_i);
    float beta = 2 * dot(dir_i,(pos_i-pos_c));
    float gamma = dot(pos_i-pos_c,pos_i-pos_c) - R*R;

    // % This calculates the square root argument of the quadratic formula
    float square_root_arguments = beta * beta - 4.0 * alpha * gamma;

    float3 pos_f;
    // If solution is not real, then exit function
    if(square_root_arguments<0.0)
    {
    	pos_f = make_float3(nan,nan,nan);
    	return pos_f;
    }

    // % This yields the solutions for the independent time variable for the cases
    // % that do intersect the sphere
    float t1 = (-beta + sqrt(square_root_arguments)) / (2.0 * alpha);
    float t2 = (-beta - sqrt(square_root_arguments)) / (2.0 * alpha);
    float t;
    // % This chooses the intersection time based upon whether the front surface
    // % of the lens or the back surface of the lens is being intersected
    if (surface == 'f'){
        // % If the front surface has positive curvature than the first
        // % intersection point should be taken, otherwise if the surface has
        // % negative curvature, then the second intersection point should be
        // % taken
        if (R > 0){
        	// % This gives the first intersection time variable (ie the minimum
			// % of 't1' and 't2' that is non-NaN, in the case where both 't1'
			// % and 't2' have NaN values, the output of 't' is also NaN)
			t = (t1 <= t2 ? t1 : t2);
        }

        else{
            // % This gives the second intersection time variable (ie the maximum
            // % of 't1' and 't2' that is non-NaN, in the case where both 't1'
            // % and 't2' have NaN values, the output of 't' is also NaN)
            t = (t1 >= t2 ? t1 : t2);
        }
    }

    else{
        // % If the back surface has positive curvature than the second
        // % intersection point should be taken, otherwise if the surface has
        // % negative curvature, then the first intersection point should be
        // % taken
        if (R > 0){
            // % This gives the second intersection time variable (ie the maximum
            // % of 't1' and 't2' that is non-NaN, in the case where both 't1'
            // % and 't2' have NaN values, the output of 't' is also NaN)
            // %t=nanmax(t1,t2);
            // NOTE: using nanmax gives wrong results (no particles make it through) even though it is consistent with
            // the comments
            t = t1 <= t2 ? t1 : t2;
        }
        else {
            // % This gives the first intersection time variable (ie the minimum
            // % of 't1' and 't2' that is non-NaN, in the case where both 't1'
            // % and 't2' have NaN values, the output of 't' is also NaN)
            // %t=nanmin(t1,t2);
            // NOTE: using nanmin gives wrong results (no particles make it through) even though it is consistent with
            // the comments (it is because, the sign of the radius of curvature changes from the front to the back surface.
            // so, for example, a convex lens has a positive front curvature and a negative back curvature)
            t = t1 >= t2 ? t1 : t2;
        }
    }
    // % This calculates the 3D coordinate of the intersection point
    pos_f = pos_i + dir_i*t;

    return pos_f;
}

__device__ float measure_distance_to_optical_axis(float3 pos_i, float3 pos_0, float plane_parameters)
{
    // % This function measures the distance from the point given by the M x 1
    // % vectors 'xi', 'yi', and 'zi' to the optical axis of the optical element
    // % with a center located at the 3 x N array given by 'lens_center' which
    // % has the lens plane defined by the 4 x N array 'plane_parameters'.  The
    // % coordinates defined by 'lens_center' must satisfy the parameters defined
    // % by 'plane_parameters', ie
    // %
    // %   plane_parameters(ii,1) * lens_center(ii,1) +
    // %       plane_parameters(ii,2) * lens_center(ii,2) +
    // %       plane_parameters(ii,3) * lens_center(ii,3) + plane_parameters(ii,4) = 0
    // %
    // % must be satisfied.  The output argument 'optical_axis_distance' is the
    // % distance from the points given by 'xi, 'yi', and 'zi' to the optical
    // % axis passing through the point 'lens_center' that is normal to the plane
    // % defined by 'plane_parameters' and will be M x N in size.

	// % This extracts the normal vector components of the plane parameters
    float a = plane_parameters.x;
    float b = plane_parameters.y;
    float c = plane_parameters.z;

    // % This calculates the minimum time (of the parametric equations describing
    // % the optical axis)
    float t_minimum = dot(make_float3(a,b,c),pos_i-pos_0)/(a*a + b*b + c*c);

    // % This is the nearest point on the optical axis to the points defined by
    // % 'xi', 'yi', and 'zi'
    float3 pos_optical_axis = pos_0 + make_float3(a,b,c)*t_minimum;

    // % This is the distance to the optical axis
    float optical_axis_distance = sqrt(dot(pos_i-pos_optical_axis,pos_i-pos_optical_axis));

    return optical_axis_distance;
}


__device__ light_ray_data_t propagate_rays_through_single_element(element_data_t optical_element, float3 element_center,
		   float4 element_plane_parameters, light_ray_data_t light_ray_data)
{

	//# % This function calculates the propagation of a set of light rays through a
	//# % single optical element.

	//# % This extracts the optical element type
	char* element_type = optical_element.element_type;

	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//# % Extraction of the light ray propagation data                            %
	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//# % This extracts the propagation direction of the light rays
	float3 ray_propagation_direction = light_ray_data.ray_propagation_direction;

	//# % This extracts the light ray source coordinates
	float3 ray_source_coordinates = light_ray_data.ray_source_coordinates;

	//# % This extracts the wavelength of the light rays
	ray_wavelength = light_ray_data.ray_wavelength;

	//# % This extracts the light ray radiance
	ray_radiance = light_ray_data.ray_radiance;

	//# % If the element type is 'lens', this extracts the optical properties
	//# % required to perform the ray propagation
	if (strcmp(element_type,'lens') == 0)
	{
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//# % Extraction of the lens optical properties                           %
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		//# % This extracts the pitch of the lens
		float element_pitch = optical_element.element_geometry.pitch;

		//# % This is the thickness of the optical element from the front optical axis
		//# % vertex to the back optical axis vertex
		double element_vertex_distance = optical_element.element_geometry.vertex_distance;

		//# % This extracts the front surface radius of curvature
		float element_front_surface_curvature = optical_element.element_geometry.front_surface_radius;

		//# % This extracts the back surface radius of curvature
		float element_back_surface_curvature = optical_element.element_geometry.back_surface_radius;

		//# % This extracts the refractive index of the current optical element
		double element_refractive_index = optical_element.element_properties.refractive_index;

		//# % This extracts the Abbe number of the current optical element
		float element_abbe_number = optical_element.element_properties.abbe_number;

		//# % This extracts the transmission ratio of the current optical
		//# % element
		float element_transmission_ratio = optical_element.element_properties.transmission_ratio;

		//# % This extracts the absorbance rate of the current optical element
		float element_absorbance_rate = optical_element.element_properties.absorbance_rate;

		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//# % propagation of the light rays through the lens front surface        %
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		//# % This extracts the individual plane parameters
		float a = element_plane_parameters.x;
		float b = element_plane_parameters.y;
		float c = element_plane_parameters.z;
		float d = element_plane_parameters.w;

		//# % This extracts the center point coordinates
		float3 pos_c = element_center;

		//# % This calculates the square of the plane normal vector magnitude
		float norm_vector_magnitude = sqrt(a*a + b*b + c*c);

		//# % This is the current offset to the center of the front spherical
		//# % surface
		float ds = +element_vertex_distance / 2.0 - element_front_surface_curvature;

		//# % This calculates the center point coordinates of the front surface
		//# % spherical shell
		float3 pos_front_surface = pos_c + make_float3(a,b,c)*ds/norm_vector_magnitude;

		//# % This calculates the points at which the light rays intersect the
		//# % front surface of the lens
		float pos_intersect = ray_sphere_intersection(pos_front_surface,element_front_surface_curvature,
							ray_propagation_direction,ray_source_coordinates, 'f');

		//# % This calculates how far the intersection point is from the optical
		//# % axis of the lens (for determining if the rays hit the lens outside of
		//# % the pitch and thus would be destroyed)
		float optical_axis_distance = measure_distance_to_optical_axis(pos_intersect,element_center,element_plane_parameters);

		//# % This gives the indices of the light rays that actually intersect the
		//# % lens (ie the rays that are within the lens pitch)
		//intersect_lens_indices = optical_axis_distance <= (element_pitch / 2.0)
		//
		//# % This sets any of the light ray directions outside of the domain of
		//# % the lens to NaN values
		//ray_propagation_direction[np.logical_not(intersect_lens_indices)][:] = np.NAN
		//
		//# % This sets any of the light ray wavelengths outside of the  domain of
		//# % the lens to NaN values
		//ray_wavelength[np.logical_not(intersect_lens_indices)] = np.NAN
		//
		//# % This sets any of the light ray radiances outside of the  domain of
		//# % the lens to NaN values
		//ray_radiance[np.logical_not(intersect_lens_indices)] = np.NAN
		//
		//# % This sets the intersection points of any of the light rays outside of
		//# % the domain of the lens to NaN values
		//x_intersect[np.logical_not(intersect_lens_indices)] = np.NAN
		//y_intersect[np.logical_not(intersect_lens_indices)] = np.NAN
		//z_intersect[np.logical_not(intersect_lens_indices)] = np.NAN
		//
		//# % This calculates the normal vectors of the lens at the intersection
		//# % points
		//lens_normal_vectors = np.array([x_intersect - xc_front_surface, y_intersect - yc_front_surface,
		//					   z_intersect - zc_front_surface])
		//
		//# % This normalizes the lens normal vectors to have unit magnitudes
		//# %lens_normal_vectors=lens_normal_vectors/norm(lens_normal_vectors,2);
		//lens_normal_vectors = lens_normal_vectors.T
		//lens_normal_vectors /= np.sqrt(
		//	lens_normal_vectors[:,0] ** 2 + lens_normal_vectors[:,1] ** 2 + lens_normal_vectors[:,2] ** 2)[:,np.newaxis]
		//
		//# % If the refractive index is a constant double value, this directly
		//# % calculates the refractive index ratio, otherwise the ratio is
		//# % calculated as a function of the wavelength
		//if isinstance(element_refractive_index,(int,long,float,complex)):
		//	# % If the Abbe number is defined, this calculates the Cauchy formula
		//	# % approxiation to the refractive index, otherwise, the refractive
		//	# % index is defined to be constant
		//	if not(np.isnan(element_abbe_number)):
		//		# % This defines the three optical wavelengths used in defining
		//		# % the Abbe number
		//		lambda_D = 589.3
		//		lambda_F = 486.1
		//		lambda_C = 656.3
		//
		//		# % This is the ratio of the incident refractive index (1 since the
		//		# % inter-lens media is assummed to be air) to the transmitted
		//		# % refractive index
		//		refractive_index_ratio = 1.0 / (
		//		element_refractive_index + (1. / (ray_wavelength ** 2) - 1 / (lambda_D ** 2)) *
		//		((element_refractive_index - 1) / (element_abbe_number * (1 / (lambda_F ** 2) - 1 / (lambda_C ** 2)))))
		//
		//	else:
		//		# % This is the ratio of the incident refractive index (1 since the
		//		# % inter-lens media is assummed to be air) to the transmitted
		//		# % refractive index
		//		refractive_index_ratio = 1.0 / element_refractive_index*np.ones(ray_propagation_direction.shape[0])
		//else:
		//	# % This defines the wavelength of the light as the variable 'lambda'
		//	# % for evaluation of the refractive index function
		//	ray_lambda = ray_wavelength
		//	# % This evaluates the string defining the refractive index in terms
		//	# % of the independent variable lambda TODO
		//	element_refractive_index_double = eval('element_refractive_index_double=' + element_refractive_index)
		//	# % This is the ratio of the incident refractive index (1 since the
		//	# % inter-lens media is assummed to be air) to the transmitted
		//	# % refractive index
		//	refractive_index_ratio = 1.0 / element_refractive_index_double*np.ones(ray_propagation_direction.shape[0])
		//
		//# % This is the scaled cosine of the angle of the incident light ray
		//# % vectors and the normal vectors of the lens (ie the dot product of the
		//# % vectors)
		//#ray_dot_product = -np.diag(np.dot(ray_propagation_direction, lens_normal_vectors.T))
		//ray_dot_product = -np.einsum('ij,ij->i',ray_propagation_direction,lens_normal_vectors)
		//
		//# % This calculates the radicand in the refraction ray propagation
		//# % direction equation
		//refraction_radicand = 1.0 - (refractive_index_ratio ** 2) * (1.0 - ray_dot_product ** 2)
		//# % This calculates the new light ray direction vectors (this is a
		//# % standard equation in optics relating incident and transmitted light
		//# % ray vectors)
		//# ray_propagation_direction = bsxfun(@times,refractive_index_ratio,ray_propagation_direction)+bsxfun(@times,(refractive_index_ratio.*ray_dot_product-sqrt(refraction_radicand)),lens_normal_vectors);
		//ray_propagation_direction = ray_propagation_direction*refractive_index_ratio[:,np.newaxis]  + \
		//							(refractive_index_ratio * ray_dot_product - np.sqrt(refraction_radicand))[:,np.newaxis] * lens_normal_vectors
		//
		//# % This normalizes the ray propagation direction so that it's magnitude
		//# % equals one
		//# ray_propagation_direction=bsxfun(@rdivide,ray_propagation_direction,sqrt(ray_propagation_direction(:,1).^2+ray_propagation_direction(:,2).^2+ray_propagation_direction(:,3).^2));
		//ray_propagation_direction /= np.sqrt(
		//	ray_propagation_direction[:,0] ** 2 + ray_propagation_direction[:,1] ** 2 + ray_propagation_direction[:,
		//		2] ** 2)[:,np.newaxis]
		//
		//# % This sets the new light ray origin to the intersection point with the
		//# % front surface of the lens
		//ray_source_coordinates = np.array([x_intersect, y_intersect, z_intersect])
		//ray_source_coordinates = ray_source_coordinates.T
		//# % Any rays that have complex values due to the radicand in the above
		//# % equation being negative experience total internal reflection.  The
		//# % values of these rays are set to NaN for now.
		//tir_indices = np.less(refraction_radicand, 0)
		//
		//# % This sets any of the light ray directions experiencing total
		//# % internal reflection to NaN values
		//ray_propagation_direction[tir_indices,:] = np.NAN
		//
		//# % This sets any of the light ray origin coordinates experiencing total
		//# % internal reflection to NaN values
		//ray_source_coordinates[tir_indices,:] = np.NAN
		//
		//# % This sets any of the light ray wavelengths experiencing total
		//# % internal reflection to NaN values
		//ray_wavelength[tir_indices] = np.NAN
		//
		//# % This sets any of the light ray radiances experiencing total internal
		//# % reflection to NaN values
		//ray_radiance[tir_indices] = np.NAN
		//
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//# % propagation of the light rays through the lens back surface         %
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//#
		//# % This is the current offset to the center of the back spherical
		//# % surface
		//ds = -element_vertex_distance / 2 - element_back_surface_curvature
		//
		//# % This calculates the center point coordinates of the back surface
		//# % spherical shell
		//xc_back_surface = xc + a * ds / norm_vector_magnitude
		//yc_back_surface = yc + b * ds / norm_vector_magnitude
		//zc_back_surface = zc + c * ds / norm_vector_magnitude
		//
		//# % This calculates the points at which the light rays intersect the
		//# % back surface of the lens
		//[x_intersect, y_intersect, z_intersect] = ray_sphere_intersection(xc_back_surface, yc_back_surface,
		//																  zc_back_surface,
		//																  element_back_surface_curvature,
		//																  ray_propagation_direction[:,0],
		//																  ray_propagation_direction[:,1],
		//																  ray_propagation_direction[:,2],
		//																  ray_source_coordinates[:,0],
		//																  ray_source_coordinates[:,1],
		//																  ray_source_coordinates[:,2], 'back')
		//
		//
		//# % This calculates how far the intersection point is from the optical
		//# % axis of the lens (for determining if the rays hit the lens outside of
		//# % the pitch and thus would be destroyed)
		//optical_axis_distance = measure_distance_to_optical_axis(x_intersect, y_intersect, z_intersect,
		//														 element_center.T, element_plane_parameters.T)
		//
		//# % This gives the indices of the light rays that actually intersect the
		//# % lens (ie the rays that are within the lens pitch)
		//intersect_lens_indices = np.less_equal(optical_axis_distance,(element_pitch / 2))
		//
		//# % This sets any of the light ray directions outside of the domain of
		//# % the lens to NaN values
		//ray_propagation_direction[np.logical_not(intersect_lens_indices),:] = np.NAN
		//
		//# % This sets any of the light ray origin coordinates outside of the
		//# % domain of the lens to NaN values
		//ray_source_coordinates[np.logical_not(intersect_lens_indices),:] = np.NAN
		//
		//# % This sets any of the light ray wavelengths outside of the  domain of
		//# % the lens to NaN values
		//ray_wavelength[np.logical_not(intersect_lens_indices)] = np.NAN
		//
		//# % This sets any of the light ray radiances outside of the  domain of
		//# % the lens to NaN values
		//ray_radiance[np.logical_not(intersect_lens_indices)] = np.NAN
		//
		//# % This sets the intersection points of any of the light rays outside of
		//# % the domain of the lens to NaN values
		//x_intersect[np.logical_not(intersect_lens_indices)] = np.NAN
		//y_intersect[np.logical_not(intersect_lens_indices)] = np.NAN
		//z_intersect[np.logical_not(intersect_lens_indices)] = np.NAN
		//
		//# % This calculates the normal vectors of the lens at the intersection
		//# % points
		//lens_normal_vectors = -1 * np.array(
		//	[x_intersect - xc_back_surface, y_intersect - yc_back_surface, z_intersect - zc_back_surface])
		//lens_normal_vectors = lens_normal_vectors.T
		//
		//# % This normalizes the lens normal vectors to have unit magnitudes
		//lens_normal_vectors /= np.sqrt(
		//	lens_normal_vectors[:,0] ** 2 + lens_normal_vectors[:,1] ** 2 + lens_normal_vectors[:,2] ** 2)[:,np.newaxis]
		//
		//# % If the refractive index is a constant double value, this directly
		//# % calculates the refractive index ratio, otherwise the ratio is
		//# % calculated as a function of the wavelength
		//if isinstance(element_refractive_index,(int,long,float,complex)):
		//	# % If the Abbe number is defined, this calculates the Cauchy formula
		//	# % approxiation to the refractive index, otherwise, the refractive
		//	# % index is defined to be constant
		//	if not(np.isnan(element_abbe_number)):
		//		# % This defines the three optical wavelengths used in defining
		//		# % the Abbe number
		//		lambda_D = 589.3
		//		lambda_F = 486.1
		//		lambda_C = 656.3
		//		# % This is the ratio of the incident refractive index to the transmitted
		//		# % refractive index (1 since the  inter-lens media is assummed to be
		//		# % air)
		//		refractive_index_ratio = element_refractive_index + (1. / (ray_wavelength ** 2) - 1 / (
		//		lambda_D ** 2)) * ((element_refractive_index - 1) / (
		//		element_abbe_number * (1 / (lambda_F ** 2) - 1 / (lambda_C ** 2))))
		//	else:
		//		# % This is the ratio of the incident refractive index to the transmitted
		//		# % refractive index (1 since the  inter-lens media is assummed to be
		//		# % air)
		//		refractive_index_ratio = element_refractive_index*np.ones(ray_propagation_direction.shape[0])
		//else:
		//	# % This defines the wavelength of the light as the variable 'lambda'
		//	# % for evaluation of the refractive index function
		//	ray_lambda = ray_wavelength
		//	# % This evaluates the string defining the refractive index in terms
		//	# % of the independent variable lambda TODO
		//	element_refractive_index_double = eval('element_refractive_index_double=' + element_refractive_index)
		//		# % This is the ratio of the incident refractive index (1 since the
		//	# % inter-lens media is assummed to be air) to the transmitted
		//	# % refractive index
		//	refractive_index_ratio = element_refractive_index_double*np.ones(ray_propagation_direction.shape[0])
		//
		//# % This is the scaled cosine of the angle of the incident light ray
		//# % vectors and the normal vectors of the lens (ie the dot product of the
		//# % vectors)
		//#ray_dot_product = -np.diag(np.dot(ray_propagation_direction, lens_normal_vectors.T))
		//ray_dot_product = -np.einsum('ij,ij->i',ray_propagation_direction,lens_normal_vectors)
		//
		//# % This calculates the radicand in the refraction ray propagation
		//# % direction equation
		//refraction_radicand = 1.0 - (refractive_index_ratio ** 2) * (1.0 - ray_dot_product ** 2)
		//# % This calculates the new light ray direction vectors (this is a
		//# % standard equation in optics relating incident and transmitted light
		//# % ray vectors)
		//ray_propagation_direction = refractive_index_ratio[:,np.newaxis] * ray_propagation_direction + (
		//																				 refractive_index_ratio * ray_dot_product - np.sqrt(
		//																					 refraction_radicand))[:,np.newaxis] * lens_normal_vectors
		//# % This normalizes the ray propagation direction so that it's magnitude
		//# % equals one
		//ray_propagation_direction /= np.sqrt(
		//	ray_propagation_direction[:,0] ** 2 + ray_propagation_direction[:,1] ** 2 + ray_propagation_direction[:,
		//		2] ** 2)[:,np.newaxis]
		//
		//# % If the absorbance rate is non-zero, this calculates how much of the
		//# % radiance is absorbed by the lens, otherwise the output radiance is
		//# % just scaled by the transmission ratio
		//if element_absorbance_rate:
		//	# % This calculates the distance that the rays traveled through
		//	# % the lens
		//	propagation_distance = np.sqrt((x_intersect - ray_source_coordinates[:,0]) ** 2 + (
		//	y_intersect - ray_source_coordinates[:,1]) ** 2 + (z_intersect - ray_source_coordinates[:,2]) ** 2)
		//
		//	# % If the absorbance rate is a simple constant, this
		//	if isinstance(element_absorbance_rate,(int,long,float,complex)):
		//		# % This calculates the new light ray radiance values
		//		ray_radiance = (1.0 - element_absorbance_rate) * ray_radiance * propagation_distance
		//
		//	else:
		//		# % These are the initial coordinates of the light ray passing
		//		# % through the lens in world coordinates
		//		x1 = ray_source_coordinates[:,0]
		//		y1 = ray_source_coordinates[:,1]
		//		z1 = ray_source_coordinates[:,2]
		//
		//		# % These are the final coordinates of the light ray passing
		//		# % through the lens in world coordinates
		//		x2 = x_intersect
		//		y2 = y_intersect
		//		z2 = z_intersect
		//
		//		# % This converts the initial light ray coordinates from the
		//		# % world coordinate system to the lens coordinate system
		//		[xL1, yL1, zL1] = convert_world_coordinates_to_lens_coordinates(x1, y1, z1, xc, yc, zc, a, b, c)
		//
		//		# % This converts the final light ray coordinates from the
		//		# % world coordinate system to the lens coordinate system
		//		[xL2, yL2, zL2] = convert_world_coordinates_to_lens_coordinates(x2, y2, z2, xc, yc, zc, a, b, c)
		//
		//		# % This calculates the radial lens coordinates for the initial
		//		# % light ray coordinates
		//		rL1 = np.sqrt(xL1 ** 2 + yL1 ** 2)
		//		# % This calculates the radial lens coordinates for the final
		//		# % light ray coordinates
		//		rL2 = np.sqrt(xL2 ** 2 + yL2 ** 2)
		//
		//		# % This defines the independent variables that may be used in
		//		# % the definition of the absorbance function
		//		x_substitution = '(xL1+(xL2-xL1)*t)'
		//		y_substitution = '(yL1+(yL2-yL1)*t)'
		//		z_substitution = '(zL1+(zL2-zL1)*t)'
		//		r_substitution = '(rL1+(rL2-rL1)*t)'
		//
		//		# % This replaces any instances of the string 'x' in the
		//		# % absorbance function string with the string 'x_substitution'
		//		element_absorbance_rate = element_absorbance_rate.replace('x', x_substitution)
		//		# % This replaces any instances of the string 'y' in the
		//		# % absorbance function string with the string 'y_substitution'
		//		element_absorbance_rate = element_absorbance_rate.replace('y', y_substitution)
		//		# % This replaces any instances of the string 'z' in the
		//		# % absorbance function string with the string 'z_substitution'
		//		element_absorbance_rate = element_absorbance_rate.replace('z', z_substitution)
		//		# % This replaces any instances of the string 'r' in the
		//		# % absorbance function string with the string 'r_substitution'
		//		element_absorbance_rate = element_absorbance_rate.replace('r', r_substitution)
		//
		//		# % This defines the function handle to the absorbance function
		//		absorbance_function_handle = eval(
		//			"lambda " + x_substitution + "," + y_substitution + "," + z_substitution + "," +
		//			r_substitution + ": " + element_absorbance_rate)
		//
		//		# % This calculates the integral of the absorbance function
		//		# % across the lens
		//		element_transmission_ratio = 1 - propagation_distance * siint.quad(absorbance_function_handle, 0.0, 1.0)
		//
		//		# % This rescales the output radiance by the transmission ratio
		//		ray_radiance = element_transmission_ratio * ray_radiance
		//
		//else:
		//	# % This rescales the output radiance by the transmission ratio
		//	ray_radiance = element_transmission_ratio * ray_radiance
		//
		//# % This sets the new light ray origin to the intersection point with the
		//# % front surface of the lens
		//ray_source_coordinates = np.array([x_intersect, y_intersect, z_intersect])
		//ray_source_coordinates = ray_source_coordinates.T
		//# % Any rays that have complex values due to the radicand in the above
		//# % equation being negative experience total internal reflection.  The
		//# % values of these rays are set to NaN for now.
		//tir_indices = np.less(refraction_radicand,0)
		//#
		//# % This sets any of the light ray directions experiencing total
		//# % internal reflection to NaN values
		//ray_propagation_direction[tir_indices,:] = np.NAN
		//# % This sets any of the light ray origin coordinates experiencing total
		//# % internal reflection to NaN values
		//ray_source_coordinates[tir_indices,:] = np.NAN
		//# % This sets any of the light ray wavelengths experiencing total
		//# % internal reflection to NaN values
		//ray_wavelength[tir_indices] = np.NAN
		//# % This sets any of the light ray radiances experiencing total internal
		//# % reflection to NaN values
		//ray_radiance[tir_indices] = np.NAN
	}
	//# % If the element type is 'aperture', this extracts the optical properties
	//# % required to perform the ray propagation
	//elif element_type == 'aperture':
	//
	//    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//    # % Extraction of the aperture optical properties                       %
	//    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//    #
	//    # % This extracts the pitch of the aperture stop
	//    element_pitch = optical_element['element_geometry']['pitch']
	//
	//    # % This is the thickness of the optical element from the front optical axis
	//    # % vertex to the back optical axis vertex
	//    element_vertex_distance = optical_element['element_geometry']['vertex_distance']
	//
	//    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//    # % propagation of the light rays through the aperture front surface    %
	//    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//    #
	//    # % This extracts the individual plane parameters
	//    a = element_plane_parameters[0]
	//    b = element_plane_parameters[1]
	//    c = element_plane_parameters[2]
	//    d = element_plane_parameters[3]
	//
	//    # % This calculates the square of the plane normal vector magnitude
	//    norm_vector_magnitude = np.sqrt(a ** 2 + b ** 2 + c ** 2)
	//
	//    # % This is the current offset to the center of the front spherical
	//    # % surface
	//    ds = -element_vertex_distance / 2.0
	//
	//    # % This calculates the transformed plane parameters (only d changes)
	//    d_temp = d - ds * norm_vector_magnitude
	//
	//    # % This is the independent intersection time between the light rays and
	//    # % the first plane of the aperture stop
	//    intersection_time = -(
	//    a * ray_source_coordinates[:,0] + b * ray_source_coordinates[:,1] + c * ray_source_coordinates[:,
	//        2] + d_temp) / (a * ray_propagation_direction[:,0] + b * ray_propagation_direction[:,1] + c *
	//                        ray_propagation_direction[:,2])
	//
	//    # % This calculates the intersection points
	//    x_intersect = ray_source_coordinates[:,0] + ray_propagation_direction[:,0] * intersection_time
	//    y_intersect = ray_source_coordinates[:,1] + ray_propagation_direction[:,1] * intersection_time
	//    z_intersect = ray_source_coordinates[:,2] + ray_propagation_direction[:,2] * intersection_time
	//
	//    # % This calculates how far the intersection point is from the optical
	//    # % axis of the lens (for determining if the rays hit the lens outside of
	//    # % the pitch and thus would be destroyed)
	//    optical_axis_distance = measure_distance_to_optical_axis(x_intersect, y_intersect, z_intersect,
	//                                                             element_center.T, element_plane_parameters.T)
	//
	//    # % This gives the indices of the light rays that actually intersect the
	//    # % aperture stop (ie the rays that are within the aperture pitch)
	//    intersect_aperture_indices = (optical_axis_distance <= (element_pitch / 2))
	//
	//    # % This sets any of the light ray directions outside of the domain of
	//    # % the lens to NaN values
	//    ray_propagation_direction[np.logical_not(intersect_aperture_indices),:] = np.NAN
	//    # % This sets any of the light ray origin coordinates outside of the
	//    # % domain of the lens to NaN values
	//    ray_source_coordinates[np.logical_not(intersect_aperture_indices),:] = np.NAN
	//    # % This sets any of the light ray wavelengths outside of the  domain of
	//    # % the lens to NaN values
	//    ray_wavelength[np.logical_not(intersect_aperture_indices)] = np.NAN
	//    # % This sets any of the light ray radiances outside of the  domain of
	//    # % the lens to NaN values
	//    ray_radiance[np.logical_not(intersect_aperture_indices)] = np.NAN
	//
	//    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//    # % propagation of the light rays through the aperture back surface     %
	//    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//
	//    # % This is the current offset to the center of the front spherical
	//    # % surface
	//    ds = +element_vertex_distance / 2
	//
	//    # % This calculates the transformed plane parameters (only d changes)
	//    d_temp = d - ds * norm_vector_magnitude
	//
	//    # % This is the independent intersection time between the light rays and
	//    # % the first plane of the aperture stop
	//    intersection_time = -(
	//    a * ray_source_coordinates[:,0] + b * ray_source_coordinates[:,1] + c * ray_source_coordinates[:,
	//        2] + d_temp) / (a * ray_propagation_direction[:,0] + b * ray_propagation_direction[:,1] + c *
	//                        ray_propagation_direction[:,2])
	//
	//    # % This calculates the intersection points
	//    x_intersect = ray_source_coordinates[:,0] + ray_propagation_direction[:,0] * intersection_time
	//    y_intersect = ray_source_coordinates[:,1] + ray_propagation_direction[:,1] * intersection_time
	//    z_intersect = ray_source_coordinates[:,2] + ray_propagation_direction[:,2] * intersection_time
	//
	//    # % This calculates how far the intersection point is from the optical
	//    # % axis of the lens (for determining if the rays hit the lens outside of
	//    # % the pitch and thus would be destroyed)
	//    optical_axis_distance = measure_distance_to_optical_axis(x_intersect, y_intersect, z_intersect,
	//                                                             element_center.T, element_plane_parameters.T)
	//
	//    # % This gives the indices of the light rays that actually intersect the
	//    # % aperture stop (ie the rays that are within the aperture pitch)
	//    intersect_aperture_indices = (optical_axis_distance <= (element_pitch / 2))
	//
	//    # % This sets any of the light ray directions outside of the domain of
	//    # % the lens to NaN values
	//    ray_propagation_direction[np.logical_not(intersect_aperture_indices),:] = np.NAN
	//    # % This sets any of the light ray origin coordinates outside of the
	//    # % domain of the lens to NaN values
	//    ray_source_coordinates[np.logical_not(intersect_aperture_indices),:] = np.NAN
	//    # % This sets any of the light ray wavelengths outside of the  domain of
	//    # % the lens to NaN values
	//    ray_wavelength[np.logical_not(intersect_aperture_indices)] = np.NAN
	//    # % This sets any of the light ray radiances outside of the  domain of
	//    # % the lens to NaN values
	//    ray_radiance[np.logical_not(intersect_aperture_indices)] = np.NAN
	//
	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//# % Saving the light ray propagation data                                   %
	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//
	//# % This extracts the propagation direction of the light rays
	//light_ray_data['ray_propagation_direction'] = ray_propagation_direction
	//#
	//# % This extracts the light ray source coordinates
	//light_ray_data['ray_source_coordinates'] = ray_source_coordinates
	//#
	//# % This extracts the wavelength of the light rays
	//light_ray_data['ray_wavelength'] = ray_wavelength
	//#
	//# % This extracts the light ray radiance
	//light_ray_data['ray_radiance'] = ray_radiance


}

__global__ void propagate_rays_through_optical_system(element_data_t element_data, float3* element_center, float4* element_plane_parameters,
		int* element_system_index,int num_elements, int num_rays, light_ray_data_t* light_ray_data)
{
	// % This function propagates the light ray data defined by the structure
	// % 'light_ray_data' through the optical system defined by the input
	// % arguments.

	//--------------------------------------------------------------------------------------
	// compute indices to access in light_ray_data
	//--------------------------------------------------------------------------------------

	// find global thread ID
	int local_thread_id = threadIdx.x;

	float del_scattering_angle = (scattering_data.scattering_angle[1] - scattering_data.scattering_angle[0])*180.0/M_PI;
	float min_scattering_angle = scattering_data.scattering_angle[0];

	// get id of particle which is the source of light rays
	int particle_id = blockIdx.x*blockDim.x + local_thread_id;

	// get id of ray emitted by the particle
	int local_ray_id = blockIdx.z;
	int global_ray_id = local_ray_id + particle_id*lightray_number_per_particle;

	// if the current ray ID is greater than the total number of rays, then exit
	if(global_ray_id >= num_rays)
		return;

	int k;
	// % This is the number of sequential optical elements within the total
	// % optical system that the light rays must be interatively passed through
	int sequential_element_number = 0;
	for(k = 0; k < num_elements; k++)
	{
		if(sequential_element_number<=element_system_index[k])
			sequential_element_number = element_system_index[k];
	}

	// % Since the the optical is defined from the sensor moving outward (i.e. the
	// % opposite direction in which the light will enter a camera system), this
	// % reverses the indexin of the optical elements so that the first element
	// % index corresponds to the first optical element that the light will hit
	for(k = 0; k < num_elements; k++)
		element_system_index[k] = sequential_element_number - element_system_index[k]; // subtracted 1

	// % This iterates through the sequential optical elements propagating the
	// % light rays through each successive element (or system of coplanar
	// % elements)
	int element_index;
	for(element_index = 0; element_index < sequential_element_number; element_index++)
	{
		// These are the indices of the current element or elements to propagate
		// the light rays through
//		current_element_indices = np.squeeze(np.argwhere(element_system_index == element_index))
		int current_element_indices[num_elements] = {0};
		int element_ctr = 0;
		for(k = 0; k < num_elements; k++)
		{
			if(element_system_index[k]==element_index)
			{
				current_element_indices[k] = k;
				element_ctr++;
			}
		}

		// % This is the number of elements that the light rays need to be
		// % simultaneously propagated through
		simultaneous_element_number = element_ctr;

		// % If there is only a single element that the rays are to be propagated
		// % through, this propagates the rays through the single element;
		// % otherwise the light rays are simultaneously propagated through the
		// % multiple elements
		if(simultaneous_element_number == 1)

			// % This extracts the current optical element data
			element_data_t current_optical_element = element_data[current_element_indices[0]];
			// % This extracts the current optical element plane parameters
			float4 current_plane_parameters = element_plane_parameters[current_element_indices[0]];
			// % This extracts the current center of the optical element
			float3 current_element_center = element_center[current_element_indices[0]];

			// % This propagates the light rays through the single optical element
			light_ray_data[global_ray_id] = propagate_rays_through_single_element(current_optical_element, current_element_center,
																   current_plane_parameters, light_ray_data[global_ray_id]);

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

			# % This propagates the light rays through the multiple optical
			# % elements
			light_ray_data = propagate_rays_through_multiple_elements(current_optical_element, current_element_center,
																	  current_plane_parameters, light_ray_data)

	}
}


extern "C"{

void read_from_file()
{
	float lens_pitch, image_distance, beam_wavelength, aperture_f_number;
	scattering_data_t scattering_data_p;
	int scattering_type;
	lightfield_source_t lightfield_source_p;
	int lightray_number_per_particle;
	int n_min; int n_max;

	/*
	 * 	This function saves the data passed from python to a file.
	 * 	This is done to enable debugging this program within eclipse.
	 */

	int k,l;

	//--------------------------------------------------------------------------------------
	// read scalars from file
	//--------------------------------------------------------------------------------------

	// open file
	string filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/scalars.bin";
	std::ifstream file_scalars(filename.c_str(), std::ios::in |
				std::ios::binary);
	// lens_pitch
	file_scalars.read((char*)&lens_pitch, sizeof(float));

	// image_distance
	file_scalars.read((char*)&image_distance, sizeof(float));

	// scattering_type
	file_scalars.read((char*)&scattering_type, sizeof(int));

	// n_min
	file_scalars.read((char*)&n_min, sizeof(int));

	// n_max
	file_scalars.read((char*)&n_max, sizeof(int));

	// lightray_number_per_particle
	file_scalars.read((char*)&lightray_number_per_particle, sizeof(int));

	// beam_wavelength
	file_scalars.read((char*)&beam_wavelength, sizeof(float));

	// aperture_f_number
	file_scalars.read((char*)&aperture_f_number, sizeof(float));

	file_scalars.close();

	//--------------------------------------------------------------------------------------
	// read scattering data
	//--------------------------------------------------------------------------------------

	// open file
	filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/scattering_data.bin";
	std::ifstream file_scattering(filename.c_str(), std::ios::in |
			std::ios::binary);
	// inverse rotation matrix
	for(k = 0; k < 9; k++)
		file_scattering.read ((char*)&scattering_data_p.inverse_rotation_matrix[k], sizeof(float));


	// beam_propagation_vector
	for(k = 0; k < 3; k++)
		file_scattering.read ((char*)&scattering_data_p.beam_propagation_vector[k], sizeof(float));

	// num_angles
	file_scattering.read ((char*)&scattering_data_p.num_angles, sizeof(int));

	// num_diameters
	file_scattering.read ((char*)&scattering_data_p.num_diameters, sizeof(int));

	// scattering_angle
	scattering_data_p.scattering_angle = (float *) malloc(scattering_data_p.num_angles*sizeof(float));
	for(k = 0; k < scattering_data_p.num_angles; k++)
			file_scattering.read ((char*)&scattering_data_p.scattering_angle[k], sizeof(float));

	// scattering_irradiance
	scattering_data_p.scattering_irradiance = (float *) malloc(scattering_data_p.num_angles * scattering_data_p.num_diameters*sizeof(float));
	for(k = 0; k < scattering_data_p.num_angles * scattering_data_p.num_diameters; k++)
			file_scattering.read ((char*)&scattering_data_p.scattering_irradiance[k], sizeof(float));

	file_scattering.close();

	//--------------------------------------------------------------------------------------
	// read lightfield_source data
	//--------------------------------------------------------------------------------------

	// open file
	filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/lightfield_source.bin";
	std::ifstream file_lightfield_source(filename.c_str(), std::ios::in |
			std::ios::binary);

	// lightray_number_per_particle
	file_lightfield_source.read ((char*)&lightfield_source_p.lightray_number_per_particle, sizeof(int));

	// lightray_process_number
	file_lightfield_source.read ((char*)&lightfield_source_p.lightray_process_number, sizeof(int));

	// num_particles
	file_lightfield_source.read ((char*)&lightfield_source_p.num_particles, sizeof(int));

	// num_rays
	file_lightfield_source.read ((char*)&lightfield_source_p.num_rays, sizeof(int));

	// diameter_index
	lightfield_source_p.diameter_index = (int *) malloc(lightfield_source_p.num_particles * sizeof(int));
	for(k = 0; k < lightfield_source_p.num_particles; k++)
		file_lightfield_source.read ((char*)&lightfield_source_p.diameter_index[k], sizeof(int));

	// radiance
	lightfield_source_p.radiance = (double *) malloc(lightfield_source_p.num_particles * sizeof(double));
	for(k = 0; k < lightfield_source_p.num_particles; k++)
		file_lightfield_source.read ((char*)&lightfield_source_p.radiance[k], sizeof(double));

	// x
	lightfield_source_p.x = (float *) malloc(lightfield_source_p.num_particles*sizeof(float));
	for(k = 0; k < lightfield_source_p.num_particles; k++)
		file_lightfield_source.read ((char*)&lightfield_source_p.x[k], sizeof(float));

	// y
	lightfield_source_p.y = (float *) malloc(lightfield_source_p.num_particles*sizeof(float));
	for(k = 0; k < lightfield_source_p.num_particles; k++)
		file_lightfield_source.read ((char*)&lightfield_source_p.y[k], sizeof(float));

	// z
	lightfield_source_p.z = (float *) malloc(lightfield_source_p.num_particles*sizeof(float));
	for(k = 0; k < lightfield_source_p.num_particles; k++)
		file_lightfield_source.read ((char*)&lightfield_source_p.z[k], sizeof(float));

	file_lightfield_source.close();


//	char* scattering_type_str;
//	if(scattering_type)
//		strcpy(scattering_type_str,"mie");
//	else
//		strcpy(scattering_type_str,"diffuse");
	char scattering_type_str[] = "mie";

	//--------------------------------------------------------------------------------------
	// save optical element data
	//--------------------------------------------------------------------------------------
	// open file
	filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/optical_elements.bin";
	std::ifstream file_optical_elements(filename.c_str(), std::ios::in |
			std::ios::binary);

	// number of elements
	int num_elements;
	file_optical_elements.read((char*)&num_elements,sizeof(int));

	// element_center
	//double (*element_center)[3] = (double *) malloc(sizeof(double*)*num_elements);
	double element_center[num_elements][3];
	for(k = 0; k < num_elements; k++)
		for(l = 0; l < 3; l++)
			file_optical_elements.read((char*)&element_center[k][l],sizeof(double));

	// element_data
	element_data_t element_data[num_elements];
	std::string temp_string;
	for(k = 0; k < num_elements; k++){
		// axial offset distance
		file_optical_elements.read((char*)&element_data[k].axial_offset_distances[0],sizeof(double));
		file_optical_elements.read((char*)&element_data[k].axial_offset_distances[1],sizeof(double));

		// element_geometry
		file_optical_elements.read((char*)&element_data[k].element_geometry.front_surface_radius,sizeof(float));
//		file_optical_elements.read((char*)&element_data[k].element_geometry.front_surface_shape,40*sizeof(char));

		file_optical_elements.read((char*)&element_data[k].element_geometry.front_surface_spherical,sizeof(bool));
		file_optical_elements.read((char*)&element_data[k].element_geometry.back_surface_radius,sizeof(float));
//		file_optical_elements.read((char*)&element_data[k].element_geometry.back_surface_shape,strlen(element_data[k].element_geometry.back_surface_shape)*sizeof(char));


		file_optical_elements.read((char*)&element_data[k].element_geometry.back_surface_spherical,sizeof(bool));
		file_optical_elements.read((char*)&element_data[k].element_geometry.pitch,sizeof(float));
		file_optical_elements.read((char*)&element_data[k].element_geometry.vertex_distance,sizeof(double));

//		file_optical_elements.read((char*)&element_data[k].element_geometry,sizeof(element_geometry_t));
		// element_number
		file_optical_elements.read((char*)&element_data[k].element_number,sizeof(float));

		// element_properties
		file_optical_elements.read((char*)&element_data[k].element_properties.abbe_number,sizeof(float));
		file_optical_elements.read((char*)&element_data[k].element_properties.absorbance_rate,sizeof(float));
		file_optical_elements.read((char*)&element_data[k].element_properties.refractive_index,sizeof(double));
		file_optical_elements.read((char*)&element_data[k].element_properties.thin_lens_focal_length,sizeof(float));
		file_optical_elements.read((char*)&element_data[k].element_properties.transmission_ratio,sizeof(float));

//		// element_type
//		file_optical_elements.read((char*)&element_data[k].element_type,strlen(element_data[k].element_type)*sizeof(char));

		// elements_coplanar
		file_optical_elements.read((char*)&element_data[k].elements_coplanar,sizeof(float));

		// rotation_angles
		for(l = 0; l < 3; l++)
			file_optical_elements.read((char*)&element_data[k].rotation_angles[l],sizeof(double));

		// z_inter_element_distance
		file_optical_elements.read((char*)&element_data[k].z_inter_element_distance,sizeof(float));

	}

	std::string temp_string_1 = "-np.sqrt((100000.000000)**2-(x**2+y**2))";

////	int len = strlen(temp_string);
////	element_data[0].element_geometry.front_surface_shape = "-np.sqrt((100000.000000)**2-(x**2+y**2))";
//	//	element_data[0].element_geometry.front_surface_shape = (char *) malloc(strlen(temp_string)*sizeof(char));
//	strcpy(element_data[0].element_geometry.front_surface_shape,temp_string.c_str());//"-np.sqrt((100000.000000)**2-(x**2+y**2))");
	element_data[0].element_geometry.front_surface_shape = &temp_string_1[0];
	printf("front_surface_shape: %s\n",element_data[0].element_geometry.front_surface_shape);
////	element_data[0].element_geometry.front_surface_shape[strlen(element_data[0].element_geometry.front_surface_shape)] = '\0';
//
	std::string temp_string_2 = "+np.sqrt((100000.000000)**2-(x**2+y**2))";
//	strcpy(element_data[0].element_geometry.back_surface_shape,temp_string.c_str());//"+np.sqrt((100000.000000)**2-(x**2+y**2))");
////	element_data[0].element_geometry.back_surface_shape = (char *) malloc(strlen(temp_string_2)*sizeof(char));
////	strcpy(element_data[0].element_geometry.back_surface_shape,temp_string_2.c_str());
	element_data[0].element_geometry.back_surface_shape = &temp_string_2[0];
	printf("back_surface_shape: %s\n",element_data[0].element_geometry.back_surface_shape);

	// element_type
	temp_string = "lens";
//	strcpy(element_data[0].element_type,"lens");//"lens");
	element_data[0].element_type = &temp_string[0];
	printf("element_type: %s\n", element_data[0].element_type);

	// element plane parameters
	double element_plane_parameters[num_elements][4];
	for(k = 0; k < num_elements; k++)
		for(l = 0; l < 4; l++)
			file_optical_elements.read((char*)&element_plane_parameters[k][l],sizeof(double));

	// element system index
	int element_system_index[num_elements];
	for(k = 0; k < num_elements; k++)
		file_optical_elements.read((char*)&element_system_index[k],sizeof(int));

	file_optical_elements.close();

	// convert 2d array to array of pointers
	double (*element_center_p)[3] = element_center;
	double (*element_plane_parameters_p)[4] = element_plane_parameters;
	// call the ray tracing function
	start_ray_tracing(lens_pitch, image_distance,&scattering_data_p, scattering_type_str,&lightfield_source_p,
			lightray_number_per_particle,n_min, n_max,beam_wavelength,aperture_f_number,
			num_elements, element_center_p,element_data,element_plane_parameters_p,element_system_index);

}

void save_to_file(float lens_pitch, float image_distance,
		scattering_data_t* scattering_data_p, char* scattering_type_str,
		lightfield_source_t* lightfield_source_p, int lightray_number_per_particle,
		int n_min, int n_max,float beam_wavelength, float aperture_f_number,
		int num_elements, double (*element_center)[3],element_data_t* element_data_p,
		double (*element_plane_parameters)[4], int *element_system_index)
{
	/*
	 * 	This function saves the data passed from python to a file.
	 * 	This is done to enable debugging this program within eclipse.
	 */

	int k,l;
	int scattering_type = strcmp(scattering_type_str,"mie")==0 ? 1 : 0;
	//--------------------------------------------------------------------------------------
	// save scalars to file
	//--------------------------------------------------------------------------------------
	printf("saving scalars to file\n");
	// open file
	string filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/scalars.bin";
	std::ofstream file_scalars(filename.c_str(), std::ios::out |
				std::ios::binary);
	// lens_pitch
	file_scalars.write((char*)&lens_pitch, sizeof(float));

	// image_distance
	file_scalars.write((char*)&image_distance, sizeof(float));

	// scattering_type
	file_scalars.write((char*)&scattering_type, sizeof(int));

	// n_min
	file_scalars.write((char*)&n_min, sizeof(int));

	// n_max
	file_scalars.write((char*)&n_max, sizeof(int));

	// lightray_number_per_particle
	file_scalars.write((char*)&lightray_number_per_particle, sizeof(int));

	// beam_wavelength
	file_scalars.write((char*)&beam_wavelength, sizeof(float));

	// aperture_f_number
	file_scalars.write((char*)&aperture_f_number, sizeof(float));

	file_scalars.close();

	//--------------------------------------------------------------------------------------
	// save scattering data to file
	//--------------------------------------------------------------------------------------
	if(strcmp(scattering_type_str,"mie")==0){

	// open file
	filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/scattering_data.bin";
	std::ofstream file_scattering(filename.c_str(), std::ios::out |
			std::ios::binary);
	// inverse rotation matrix
	printf("inverse rotation matrix: \n");
	for(k = 0; k < 9; k++){
		if(k%3==0)
			printf("\n");
		printf("%f ",scattering_data_p->inverse_rotation_matrix[k]);
		scattering_data_p->inverse_rotation_matrix[k] = (float) scattering_data_p->inverse_rotation_matrix[k];
			file_scattering.write ((char*)&scattering_data_p->inverse_rotation_matrix[k], sizeof(float));

//		}
		printf("\n");
	}
	printf("\n");
	// beam_propagation_vector
	printf("beam propagation vector\n");
	for(k = 0; k < 3; k++){
		printf("%f ", scattering_data_p->beam_propagation_vector[k]);
		file_scattering.write ((char*)&scattering_data_p->beam_propagation_vector[k], sizeof(float));
	}
	// num_angles
	file_scattering.write ((char*)&scattering_data_p->num_angles, sizeof(int));
	// num_diameters
	file_scattering.write ((char*)&scattering_data_p->num_diameters, sizeof(int));
	// scattering_angle
	for(k = 0; k < scattering_data_p->num_angles; k++)
			file_scattering.write ((char*)&scattering_data_p->scattering_angle[k], sizeof(float));
	// scattering_irradiance
	for(k = 0; k < scattering_data_p->num_angles * scattering_data_p->num_diameters; k++)
			file_scattering.write ((char*)&scattering_data_p->scattering_irradiance[k], sizeof(float));

	file_scattering.close();
	}
	//--------------------------------------------------------------------------------------
	// save lightfield_source data to file
	//--------------------------------------------------------------------------------------

	// open file
	filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/lightfield_source.bin";
	std::ofstream file_lightfield_source(filename.c_str(), std::ios::out |
			std::ios::binary);

	// lightray_number_per_particle
	file_lightfield_source.write ((char*)&lightfield_source_p->lightray_number_per_particle, sizeof(int));
	// lightray_process_number
	file_lightfield_source.write ((char*)&lightfield_source_p->lightray_process_number, sizeof(int));
	// num_particles
	file_lightfield_source.write ((char*)&lightfield_source_p->num_particles, sizeof(int));
	// num_rays
	file_lightfield_source.write ((char*)&lightfield_source_p->num_rays, sizeof(int));
	// diameter_index
	for(k = 0; k < lightfield_source_p->num_particles; k++)
		file_lightfield_source.write ((char*)&lightfield_source_p->diameter_index[k], sizeof(int));
	// radiance
	for(k = 0; k < lightfield_source_p->num_particles; k++)
		file_lightfield_source.write ((char*)&lightfield_source_p->radiance[k], sizeof(double));
	// x
	for(k = 0; k < lightfield_source_p->num_particles; k++)
		file_lightfield_source.write ((char*)&lightfield_source_p->x[k], sizeof(float));
	// y
		for(k = 0; k < lightfield_source_p->num_particles; k++)
			file_lightfield_source.write ((char*)&lightfield_source_p->y[k], sizeof(float));
	// z
	for(k = 0; k < lightfield_source_p->num_particles; k++)
		file_lightfield_source.write ((char*)&lightfield_source_p->z[k], sizeof(float));

	file_lightfield_source.close();

	//--------------------------------------------------------------------------------------
	// save optical element data
	//--------------------------------------------------------------------------------------
	// open file
	filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/optical_elements.bin";
	std::ofstream file_optical_elements(filename.c_str(), std::ios::out |
			std::ios::binary);
	printf("\n");
	// number of elements
	file_optical_elements.write((char*)&num_elements,sizeof(int));
	printf("num_elements: %d\n",num_elements);

	// element_center
	printf("element_center: \n");
	for(k = 0; k < num_elements; k++)
	{
		for(l = 0; l < 3; l++)
		{
			file_optical_elements.write((char*)&element_center[k][l],sizeof(double));
			printf("%f ",element_center[k][l]);
		}
		printf("\n");
	}

	// element_data
	printf("element_data\n");
	for(k = 0; k < num_elements; k++){
		printf("element_number: %d\n", k+1);
		// axial offset distance
		file_optical_elements.write((char*)&element_data_p[k].axial_offset_distances[0],sizeof(double));
		file_optical_elements.write((char*)&element_data_p[k].axial_offset_distances[1],sizeof(double));
		printf("axial_offset_distances: %f, %f\n", element_data_p[k].axial_offset_distances[0],element_data_p[k].axial_offset_distances[1]);

		// element_geometry
		printf("element_geometry:\n");
		file_optical_elements.write((char*)&element_data_p[k].element_geometry.front_surface_radius,sizeof(float));
		printf("front_surface_radius: %f\n", element_data_p[k].element_geometry.front_surface_radius);
//		file_optical_elements.write((char*)&element_data_p[k].element_geometry.front_surface_shape,strlen(element_data_p->element_geometry.front_surface_shape)*sizeof(char));
		printf("front_surface_shape: %s\n", element_data_p[k].element_geometry.front_surface_shape);
		file_optical_elements.write((char*)&element_data_p[k].element_geometry.front_surface_spherical,sizeof(bool));
		printf("front_surface_spherical: %s\n",element_data_p[k].element_geometry.front_surface_spherical ? "true":"false");
		file_optical_elements.write((char*)&element_data_p[k].element_geometry.back_surface_radius,sizeof(float));
		printf("back_surface_radius: %f\n", element_data_p[k].element_geometry.back_surface_radius);
//		file_optical_elements.write((char*)&element_data_p[k].element_geometry.back_surface_shape,strlen(element_data_p->element_geometry.back_surface_shape)*sizeof(char));
		printf("back_surface_shape: %s\n", element_data_p[k].element_geometry.back_surface_shape);
		file_optical_elements.write((char*)&element_data_p[k].element_geometry.back_surface_spherical,sizeof(bool));
		printf("back_surface_spherical: %s\n",element_data_p[k].element_geometry.back_surface_spherical ? "true":"false");
		file_optical_elements.write((char*)&element_data_p[k].element_geometry.pitch,sizeof(float));
		printf("pitch: %f\n",element_data_p[k].element_geometry.pitch);
		file_optical_elements.write((char*)&element_data_p[k].element_geometry.vertex_distance,sizeof(double));
		printf("vertex_distance: %f\n",element_data_p[k].element_geometry.vertex_distance);

		// element_number
		file_optical_elements.write((char*)&element_data_p[k].element_number,sizeof(float));
		printf("element_number: %f\n",element_data_p[k].element_number);

		// element_properties
		printf("element_properties:\n");
		file_optical_elements.write((char*)&element_data_p[k].element_properties.abbe_number,sizeof(float));
		printf("abbe_number: %f\n",element_data_p[k].element_properties.abbe_number);
		file_optical_elements.write((char*)&element_data_p[k].element_properties.absorbance_rate,sizeof(float));
		printf("absorbance_rage: %f\n",element_data_p[k].element_properties.absorbance_rate);
		file_optical_elements.write((char*)&element_data_p[k].element_properties.refractive_index,sizeof(double));
		printf("refractive_index: %f\n", element_data_p[k].element_properties.refractive_index);
		file_optical_elements.write((char*)&element_data_p[k].element_properties.thin_lens_focal_length,sizeof(float));
		printf("thin_lens_focal_length: %f\n",element_data_p[k].element_properties.thin_lens_focal_length);
		file_optical_elements.write((char*)&element_data_p[k].element_properties.transmission_ratio,sizeof(float));
		printf("transmission_ratio : %f\n",element_data_p[k].element_properties.transmission_ratio);

		// element_type
//		file_optical_elements.write((char*)&element_data_p[k].element_type,strlen(element_data_p->element_type)*sizeof(char));
		printf("element_type: %s\n",element_data_p[k].element_type);

		// elements_coplanar
		file_optical_elements.write((char*)&element_data_p[k].elements_coplanar,sizeof(float));
		printf("element_coplanar: %f\n",element_data_p[k].elements_coplanar);

		// rotation_angles
		printf("rotation_angles: ");
		for(l = 0; l < 3; l++)
		{
			file_optical_elements.write((char*)&element_data_p[k].rotation_angles[l],sizeof(double));
			printf("%f ",element_data_p[k].rotation_angles[l]);
		}
		printf("\n");

		// z_inter_element_distance
		file_optical_elements.write((char*)&element_data_p[k].z_inter_element_distance,sizeof(float));
		printf("z_inter_element_distance: %f\n",element_data_p[k].z_inter_element_distance);

	}

	// element plane parameters
	printf("element_plane_parameters: \n");
	for(k = 0; k < num_elements; k++)
	{
		for(l = 0; l < 4; l++)
		{
			file_optical_elements.write((char*)&element_plane_parameters[k][l],sizeof(double));
			printf("%f ", element_plane_parameters[k][l]);
		}

		printf("\n");
	}

	// element system index
	printf("element_system_index:\n");
	for(k = 0; k < num_elements; k++)
	{
		file_optical_elements.write((char*)&element_system_index[k],sizeof(int));
		printf("%d ",element_system_index[k]);
	}

	file_optical_elements.close();

}

int add(int a, int b)
{
	return a+b;
}

void start_ray_tracing(float lens_pitch, float image_distance,
		scattering_data_t* scattering_data_p, char* scattering_type_str,
		lightfield_source_t* lightfield_source_p, int lightray_number_per_particle,
		int n_min, int n_max,float beam_wavelength, float aperture_f_number,
		int num_elements, double (*element_center)[3],element_data_t* element_data_p,
				double (*element_plane_parameters)[4], int *element_system_index)
{
	// create instance of structure using the pointers
	scattering_data_t scattering_data = *scattering_data_p;
	lightfield_source_t lightfield_source = *lightfield_source_p;

	int source_point_number = n_max - n_min + 1;

	// allocate space for the light field variables on the CPU

	int N = lightray_number_per_particle*source_point_number;

	//--------------------------------------------------------------------------------------
	// allocate space on GPU for lightfield_source
	//--------------------------------------------------------------------------------------

	// declare pointers to device arrays
	float* d_source_x;
	float* d_source_y;
	float* d_source_z;
	double *d_source_radiance;
	int *d_source_diameter_index;
	int num_particles = lightfield_source.num_particles;

	// allocate space for device arrays on GPU
	//cudaMalloc((void **)&gpuData, sizeof(float)*size);
	float *gpuData;
	int size = 10;
	cudaThreadSynchronize();
	cudaMalloc((void **)&gpuData, sizeof(float)*size);
	cudaMalloc((void **)&d_source_x,sizeof(float)*num_particles);
	cudaMalloc((void **)&d_source_y,num_particles*sizeof(float));
	cudaMalloc((void **)&d_source_z,num_particles*sizeof(float));
	cudaMalloc((void **)&d_source_radiance,num_particles*sizeof(double));
	cudaMalloc((void **)&d_source_diameter_index,num_particles*sizeof(int));

	// copy data to GPU
	cudaMemcpy(d_source_x,lightfield_source.x,num_particles*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(d_source_y,lightfield_source.y,num_particles*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(d_source_z,lightfield_source.z,num_particles*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(d_source_radiance,lightfield_source.radiance,num_particles*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(d_source_diameter_index,lightfield_source.diameter_index,num_particles*sizeof(int),cudaMemcpyHostToDevice);

	// make copy of host structure
	lightfield_source_t  lightfield_source_copy = lightfield_source;

	// point host structure to device array
	lightfield_source.x = d_source_x;
	lightfield_source.y = d_source_y;
	lightfield_source.z = d_source_z;
	lightfield_source.radiance = d_source_radiance;
	lightfield_source.diameter_index = d_source_diameter_index;

	//--------------------------------------------------------------------------------------
	// allocate space on GPU for scattering_data
	//--------------------------------------------------------------------------------------

	// declare pointers to device arrays
	float *d_scattering_angle;
	float* d_scattering_irradiance;

	// allocate space for device arrays on GPU
	cudaMalloc((void**)&d_scattering_angle,scattering_data.num_angles*sizeof(float));
	cudaMalloc((void**)&d_scattering_irradiance,scattering_data.num_angles*scattering_data.num_diameters*sizeof(float));

	// copy data to GPU
	cudaMemcpy(d_scattering_angle,scattering_data.scattering_angle,scattering_data.num_angles*sizeof(float)
	,cudaMemcpyHostToDevice);
	cudaMemcpy(d_scattering_irradiance,scattering_data.scattering_irradiance,scattering_data.num_angles*scattering_data.num_diameters*sizeof(float)
		,cudaMemcpyHostToDevice);

	// make copy of host structure
	scattering_data_t scattering_data_copy = scattering_data;

	// point host structure to device array
	scattering_data.scattering_angle = d_scattering_angle;

	cudaMalloc((void**)&data_array,scattering_data.num_angles*scattering_data.num_diameters*sizeof(float));
	cudaMemcpy(data_array,scattering_data.scattering_irradiance,
				scattering_data.num_angles*scattering_data.num_diameters*sizeof(float),cudaMemcpyHostToDevice);
	cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();
	cudaBindTexture2D( NULL, mie_scattering_irradiance,
	                               data_array,
	                               desc, scattering_data.num_diameters, scattering_data.num_angles,
	                               sizeof(float) * scattering_data.num_diameters );
	scattering_data.scattering_irradiance = d_scattering_irradiance;

	//--------------------------------------------------------------------------------------
	// allocate space on GPU for light_ray_data
	//--------------------------------------------------------------------------------------

//	// allocate space for light_ray_data on CPU
//	light_ray_data_t light_ray_data;
//	light_ray_data.ray_source_coordinates = (float3 *) malloc(N*sizeof(float3));
//	light_ray_data.ray_propagation_direction = (float3 *) malloc(N*sizeof(float3));
//	light_ray_data.ray_wavelength = (float *) malloc(N*sizeof(float));
//	light_ray_data.ray_radiance = (double *) malloc(N*sizeof(double));
//	light_ray_data.num_lightrays = N;
//	// declare pointers to GPU arrays
//	float3 *d_ray_source_coordinates, *d_ray_propagation_direction;
//	float *d_ray_wavelength;
//	double *d_ray_radiance;
//
//	// allocate memory on GPU
//	cudaMalloc((void**)&d_ray_source_coordinates, N*sizeof(float3));
//	cudaMalloc((void**)&d_ray_propagation_direction, N*sizeof(float3));
//	cudaMalloc((void**)&d_ray_wavelength, N*sizeof(float));
//	cudaMalloc((void**)&d_ray_radiance, N*sizeof(double));
//
//	// initialize arrays to zero
//	cudaMemset(d_ray_source_coordinates,0.0,N*sizeof(float3));
//	cudaMemset(d_ray_propagation_direction,0.0,N*sizeof(float3));
//	cudaMemset(d_ray_wavelength,0.0,N*sizeof(float));
//	cudaMemset(d_ray_radiance,0.0,N*sizeof(double));
//
//	// copy contents of light_ray_data structure
//	light_ray_data_t light_ray_data_copy = light_ray_data;
//	// point structure to device arrays
//	light_ray_data.ray_source_coordinates = d_ray_source_coordinates;
//	light_ray_data.ray_propagation_direction = d_ray_propagation_direction;
//	light_ray_data.ray_wavelength = d_ray_wavelength;
//	light_ray_data.ray_radiance = d_ray_radiance;

	int scattering_type = 0;
	if(strcmp(scattering_type_str,"mie")==0)
		scattering_type = 1;
	int num_rays = N;

//	light_ray_data_t light_ray_data[num_rays];
	light_ray_data_t* d_light_ray_data;

	cudaMalloc((void**)&d_light_ray_data,num_rays*sizeof(light_ray_data_t));
	// allocate threads per block
	dim3 block(10,1,1);
	// allocate blocks per grid
	dim3 grid(source_point_number/block.x,1,lightray_number_per_particle);

	// call kernel
	generate_lightfield_angular_data<<<grid,block>>>(lens_pitch, image_distance,scattering_data,
			scattering_type, lightfield_source,lightray_number_per_particle, n_min, n_max,
			beam_wavelength,aperture_f_number,d_light_ray_data);

	cudaThreadSynchronize();

	light_ray_data_t light_ray_data[num_rays];
	// copy light_ray_data back to host
	cudaMemcpy(light_ray_data,d_light_ray_data,num_rays*sizeof(light_ray_data_t),cudaMemcpyDeviceToHost);

	cudaThreadSynchronize();

	// display first and last few elements of lightfield_data

//	printf("lightfield_data contents\n");
//	printf("ray_source_coordinates (1st): %f, %f, %f\n",light_ray_data[0].ray_source_coordinates.x,light_ray_data[0].ray_source_coordinates.y,light_ray_data[0].ray_source_coordinates.z);
//	printf("ray_source_coordinates (last): %f, %f, %f\n",light_ray_data[num_rays-1].ray_source_coordinates.x,light_ray_data[num_rays-1].ray_source_coordinates.y,light_ray_data[N-1].ray_source_coordinates.z);
//	printf("ray_propagation_direction (1st): %f, %f, %f\n",light_ray_data[0].ray_propagation_direction.x,light_ray_data[0].ray_propagation_direction[0].y,light_ray_data.ray_propagation_direction[0].z);
//	printf("ray_propagation_direction (last): %f, %f, %f\n",light_ray_data.ray_propagation_direction[N-1].x,light_ray_data.ray_propagation_direction[N-1].y,light_ray_data.ray_propagation_direction[N-1].z);
//	printf("ray_wavelength (1st, last): %f, %f\n",light_ray_data.ray_wavelength[0],light_ray_data.ray_wavelength[N-1]);
//	printf("ray_radiance (1st, last): %f, %f\n",light_ray_data.ray_radiance[0],light_ray_data.ray_radiance[N-1]);

	// free pointers
	cudaFree(gpuData);
	cudaFree(d_source_x);
	cudaFree(d_source_y);
	cudaFree(d_source_z);
	cudaFree(d_source_radiance);
	cudaFree(d_source_diameter_index);

	cudaFree(d_scattering_angle);
	cudaFree(d_scattering_irradiance);
	cudaFree(data_array);

	/*
	 * propagate Rays through the optical system
	 */

	/*
	 *  convert co-ordinate arrays to float3 and float4 arrays so that they can
	 *  be represented by a 1D array.
	 */

	float3* element_center_2 = (float3 *) malloc(num_elements*sizeof(float3));
	float4* element_plane_parameters_2 = (float4 *) malloc(num_elements*sizeof(float4));

	for(k = 0; k < num_elements; k++)
	{
		element_center_2[k] = make_float3(element_center[k][0],element_center[k][1],element_center[k][2]);
		element_plane_parameters_2 = make_float4(element_plane_parameters[k][0],element_plane_parameters[k][1],element_plane_parameters[k][2],element_plane_parameters[k][3],element_plane_parameters[k][4]);
	}

	// allocate space on GPU for coordinate arrays
	float3* d_element_center;
	float4* d_element_plane_parameters;
	int* d_element_system_index;

	cudaMalloc((void**)&d_element_center,num_elements*sizeof(float3));
	cudaMalloc((void**)&d_element_plane_parameters,num_elements*sizeof(float4));
	cudaMalloc((void**)&d_element_system_index,num_elements*sizeof(int));

	// copy data from GPU to CPU
	cudaMemcpy(d_element_center,element_center_2,num_elements*sizeof(float3),cudaMemcpyHostToDevice);
	cudaMemcpy(d_element_plane_parameters,element_plane_parameters_2,num_elements*sizeof(float4),cudaMemcpyHostToDevice);
	cudaMemcpy(d_element_system_index,element_system_index,num_elements*sizeof(int),cudaMemcpyHostToDevice);

	element_data_t element_data = *element_data_p;

	propagate_rays_through_optical_system<<<grid,block>>>>(element_data,d_element_center,d_element_plane_parameters,d_element_system_index,num_elements,num_rays,d_light_ray_data);

	cudaFree(d_light_ray_data);
	udaFree(d_element_center);
	cudaFree(d_element_plane_parameters);
	cudaFree(d_element_system_index);


}



}


