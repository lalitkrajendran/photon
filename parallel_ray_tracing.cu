/*
 * parallel_ray_tracing.cu
 *
 *  Created on: April 20, 2016
 *      Author: lrajendr
 *
 *	This program
 *
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
#include <time.h>


using namespace std;
//#incldue <cutil_math.h>

#define CUDART_NAN_F            __int_as_float(0x7fffffff)
#define CUDART_NAN              __longlong_as_double(0xfff8000000000000ULL)

cudaArray *data_array = 0;
texture<float, 2> mie_scattering_irradiance;

__shared__ lightfield_source_single_t lightfield_source_shared;
__constant__ camera_design_t camera_design_const;


__device__ float random_single(unsigned int seed, int id)
{
	/*
	 * this function generates a single random float for a uniform distribution using
	 * the cuRAND library
	 *
	 */

  /* CUDA's random number library uses curandState_t to keep track of the seed value
     we will store a random state for every thread  */
  curandState_t state;
//  int id = threadIdx.x + blockIdx.x * blockDim.x;

  /* the seed can be the same for each core, here we pass the time in from the CPU */
  /* the sequence number should be different for each core (unless you want all
                             cores to get the same sequence of numbers for some reason - use thread id! */
  /* the offset is how much extra we advance in the sequence for each call, can be 0 */

  /* we have to initialize the state */
//  curand_init(seed, blockIdx.x, 0, &state);
  curand_init(seed,0,id,&state);

  float rand_num = curand_uniform(&state);
  return rand_num;
}

__device__ light_ray_data_t generate_lightfield_angular_data(float lens_pitch, float image_distance,
		scattering_data_t scattering_data, int scattering_type, lightfield_source_single_t lightfield_source,
                                     int lightray_number_per_particle, int n_min, int n_max, float beam_wavelength,
                                     float aperture_f_number, int num_rays,
                                     float rand_num_1,float rand_num_2)

{
	/*
		This function generates the light field data for the source points specified by the
		structure lightfield_source.  The data is only generated for the source points from
		n_min to n_max.  The parameter lightray_number_per_particle is the number of rays to generate for each
		source point.
	*/

	// get difference in scattering angle
	float del_scattering_angle = (scattering_data.scattering_angle[1] - scattering_data.scattering_angle[0])*180.0/M_PI;

	// get source coordinates of the light ray
	float x_current = lightfield_source.x;
	float y_current = lightfield_source.y;
	float z_current = lightfield_source.z;

	//--------------------------------------------------------------------------------------
	// compute direction of propagation of light ray
	//--------------------------------------------------------------------------------------

//	// generate random numbers corresponding to radial and angular coordinates on the lens
//	float random_number_1 = random_single(particle_id * local_ray_id * global_ray_id);
//	float random_number_2 = random_single(particle_id + local_ray_id + global_ray_id);
//	unsigned int seed1 = 1105;
//	unsigned int seed2 = 4092;
//	float random_number_1 = random_single(seed1,local_ray_id); //global_ray_id);
//	float random_number_2 = random_single(seed2,local_ray_id); //global_ray_id);

//	// generate random points on the lens where the rays should intersect
//	float x_lens = 0.5*lens_pitch*sqrt(random_number_1)*cos(2*M_PI*random_number_2);
//	float y_lens = 0.5*lens_pitch*sqrt(random_number_1)*sin(2*M_PI*random_number_2);

	// generate random points on the lens where the rays should intersect
	float x_lens = 0.5*lens_pitch*sqrt(rand_num_1)*cos(2*M_PI*rand_num_2);
	float y_lens = 0.5*lens_pitch*sqrt(rand_num_1)*sin(2*M_PI*rand_num_2);

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
		diameter_index = lightfield_source.diameter_index; //[current_source_point_number];
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
		irradiance_current=ray_scattering_irradiance*lightfield_source.radiance; //[current_source_point_number];
	}
	// if not mie scattering, then set irradiance to be uniform (diffuse)
	else
	{
		// % This specifies the total irradiance for the current particle's
		// % rays to be uniform
		irradiance_current = lightfield_source.radiance; //[current_source_point_number];
	}



	// save the light rays to the light ray data structure
	light_ray_data_t light_ray_data;

	light_ray_data.ray_source_coordinates = make_float3(x_current,y_current,z_current);
	light_ray_data.ray_propagation_direction = normalize(make_float3(theta_temp,phi_temp,-1.0));
	light_ray_data.ray_wavelength = beam_wavelength;
	light_ray_data.ray_radiance = 1/(aperture_f_number*aperture_f_number)*irradiance_current;

	return light_ray_data;
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
//    printf("alpha: %f, beta: %f, gamma: %f, square_root_arguments: %f\n",alpha,beta,gamma,square_root_arguments);
    float3 pos_f;
    // If solution is not real, then exit function
    if(square_root_arguments<0.0)
    {
    	pos_f = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);
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

__device__ float measure_distance_to_optical_axis(float3 pos_i, float3 pos_0, float4 plane_parameters)
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
	char element_type = optical_element.element_type;

	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//# % Extraction of the light ray propagation data                            %
	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//# % This extracts the propagation direction of the light rays
	float3 ray_propagation_direction = light_ray_data.ray_propagation_direction;

	//# % This extracts the light ray source coordinates
	float3 ray_source_coordinates = light_ray_data.ray_source_coordinates;

	//# % This extracts the wavelength of the light rays
	float ray_wavelength = light_ray_data.ray_wavelength;

	//# % This extracts the light ray radiance
	double ray_radiance = light_ray_data.ray_radiance;

	float element_pitch;
	double element_vertex_distance;
	float a,b,c,d,norm_vector_magnitude, optical_axis_distance,ds;
	float3 pos_c,pos_intersect;
	//# % If the element type is 'lens', this extracts the optical properties
	//# % required to perform the ray propagation
	if (element_type == 'l')
	{
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//# % Extraction of the lens optical properties                           %
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		//# % This extracts the pitch of the lens
		element_pitch = optical_element.element_geometry.pitch;

		//# % This is the thickness of the optical element from the front optical axis
		//# % vertex to the back optical axis vertex
		element_vertex_distance = optical_element.element_geometry.vertex_distance;

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
		a = element_plane_parameters.x;
		b = element_plane_parameters.y;
		c = element_plane_parameters.z;
		d = element_plane_parameters.w;

		//# % This extracts the center point coordinates
		pos_c = element_center;

		//# % This calculates the square of the plane normal vector magnitude
		norm_vector_magnitude = sqrt(a*a + b*b + c*c);

		//# % This is the current offset to the center of the front spherical
		//# % surface
		ds = +element_vertex_distance / 2.0 - element_front_surface_curvature;

		//# % This calculates the center point coordinates of the front surface
		//# % spherical shell
		float3 pos_front_surface = pos_c + make_float3(a,b,c)*ds/norm_vector_magnitude;

		//# % This calculates the points at which the light rays intersect the
		//# % front surface of the lens
		pos_intersect = ray_sphere_intersection(pos_front_surface,element_front_surface_curvature,
							ray_propagation_direction,ray_source_coordinates, 'f');

		//# % This calculates how far the intersection point is from the optical
		//# % axis of the lens (for determining if the rays hit the lens outside of
		//# % the pitch and thus would be destroyed)
		optical_axis_distance = measure_distance_to_optical_axis(pos_intersect,element_center,element_plane_parameters);

		//# % This checks if the given light ray actually intersects the lens
		//# % (ie the rays that are within the lens pitch)
		if(optical_axis_distance > element_pitch/2.0)
		{

			light_ray_data.ray_source_coordinates = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);
			//# % This sets any of the light ray directions outside of the domain of
			//# % the lens to NaN values
			light_ray_data.ray_propagation_direction = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			//# % This sets any of the light ray wavelengths outside of the  domain of
			//# % the lens to NaN values
			light_ray_data.ray_wavelength = CUDART_NAN_F;

			//# % This sets any of the light ray radiances outside of the  domain of
			//# % the lens to NaN values
			light_ray_data.ray_radiance = CUDART_NAN;

			//# % This sets the intersection points of any of the light rays outside of
			//# % the domain of the lens to NaN values
			pos_intersect = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			return light_ray_data;
		}

		//# % This calculates the normal vectors of the lens at the intersection
		//# % points
		float3 lens_normal_vectors = pos_intersect - pos_front_surface;

		//# % This normalizes the lens normal vectors to have unit magnitudes
		//# %lens_normal_vectors=lens_normal_vectors/norm(lens_normal_vectors,2);
		lens_normal_vectors = normalize(lens_normal_vectors);

		//# % If the refractive index is a constant double value, this directly
		//# % calculates the refractive index ratio, otherwise the ratio is
		//# % calculated as a function of the wavelength
		// (I DON'T KNOW HOW TO DO STRING EVALUATION IN CUDA. SO I AM ONLY GOING TO
		// CONSIDER THE FIRST CASE - LKR.)
		// TODO : implement eval

		//# % If the Abbe number is defined, this calculates the Cauchy formula
		//# % approximation to the refractive index, otherwise, the refractive
		//# % index is defined to be constant
		float refractive_index_ratio, lambda_D, lambda_F, lambda_C;
		if(!isnan(element_abbe_number)){

			//# % This defines the three optical wavelengths used in defining
			//# % the Abbe number
			lambda_D = 589.3;
			lambda_F = 486.1;
			lambda_C = 656.3;

			//# % This is the ratio of the incident refractive index (1 since the
			//# % inter-lens media is assummed to be air) to the transmitted
			//# % refractive index
			refractive_index_ratio = 1.0 /(element_refractive_index + (1. / (ray_wavelength*ray_wavelength)
					- 1 / (lambda_D*lambda_D)) *((element_refractive_index - 1) / (element_abbe_number * (1 / (lambda_F*lambda_F) - 1 / (lambda_C*lambda_C)))));

		}
		else{
		//# % This is the ratio of the incident refractive index (1 since the
		//# % inter-lens media is assumed to be air) to the transmitted
		//# % refractive index
		refractive_index_ratio = 1.0 / element_refractive_index;

		}

		//# % This is the scaled cosine of the angle of the incident light ray
		//# % vectors and the normal vectors of the lens (ie the dot product of the
		//# % vectors)
		float ray_dot_product = -dot(ray_propagation_direction,lens_normal_vectors);

		//# % This calculates the radicand in the refraction ray propagation
		//# % direction equation
		float refraction_radicand = 1.0 - (refractive_index_ratio*refractive_index_ratio)
				* (1.0 - ray_dot_product*ray_dot_product);

		//# % Any rays that have complex values due to the radicand in the above
		//# % equation being negative experience total internal reflection.  The
		//# % values of these rays are set to NaN for now.
		if(refraction_radicand < 0.0){

			//# % This sets any of the light ray directions experiencing total
			//# % internal reflection to NaN values
			light_ray_data.ray_propagation_direction = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			//# % This sets any of the light ray origin coordinates experiencing total
			//# % internal reflection to NaN values
			light_ray_data.ray_source_coordinates = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			//# % This sets any of the light ray wavelengths experiencing total
			//# % internal reflection to NaN values
			light_ray_data.ray_wavelength = CUDART_NAN_F;

			//# % This sets any of the light ray radiance experiencing total internal
			//# % reflection to NaN values
			light_ray_data.ray_radiance = CUDART_NAN;

			return light_ray_data;
		}

		//# % This calculates the new light ray direction vectors (this is a
		//# % standard equation in optics relating incident and transmitted light
		//# % ray vectors)
		ray_propagation_direction = ray_propagation_direction*refractive_index_ratio +
				(refractive_index_ratio*ray_dot_product - sqrt(refraction_radicand))*lens_normal_vectors;

		//# % This normalizes the ray propagation direction so that it's magnitude
		//# % equals one
		ray_propagation_direction = normalize(ray_propagation_direction);

		//# % This sets the new light ray origin to the intersection point with the
		//# % front surface of the lens
		ray_source_coordinates = pos_intersect;

		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//# % propagation of the light rays through the lens back surface         %
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		//# % This is the current offset to the center of the back spherical
		//# % surface
		ds = -element_vertex_distance / 2 - element_back_surface_curvature;

		//# % This calculates the center point coordinates of the back surface
		//# % spherical shell
		float3 pos_back_surface = pos_c + make_float3(a,b,c)*ds/norm_vector_magnitude;

		//# % This calculates the points at which the light rays intersect the
		//# % back surface of the lens

		pos_intersect = ray_sphere_intersection(pos_back_surface,element_back_surface_curvature,
							  ray_propagation_direction,ray_source_coordinates, 'b');

		//# % This calculates how far the intersection point is from the optical
		//# % axis of the lens (for determining if the rays hit the lens outside of
		//# % the pitch and thus would be destroyed)
		optical_axis_distance = measure_distance_to_optical_axis(pos_intersect,element_center,element_plane_parameters);

		//# % This gives the indices of the light rays that actually intersect the
		//# % lens (ie the rays that are within the lens pitch)
		//intersect_lens_indices = np.less_equal(optical_axis_distance,(element_pitch / 2))

		//# % This checks if the given light ray actually intersects the lens
		//# % (ie the rays that are within the lens pitch)
		if(optical_axis_distance > element_pitch/2.0){

			light_ray_data.ray_source_coordinates = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);
			//# % This sets any of the light ray directions outside of the domain of
			//# % the lens to NaN values
			light_ray_data.ray_propagation_direction = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			//# % This sets any of the light ray wavelengths outside of the  domain of
			//# % the lens to NaN values
			light_ray_data.ray_wavelength = CUDART_NAN_F;

			//# % This sets any of the light ray radiances outside of the  domain of
			//# % the lens to NaN values
			light_ray_data.ray_radiance = CUDART_NAN;

			//# % This sets the intersection points of any of the light rays outside of
			//# % the domain of the lens to NaN values
			pos_intersect = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			return light_ray_data;
		}

		//# % This calculates the normal vectors of the lens at the intersection
		//# % points
		lens_normal_vectors = -(pos_intersect - pos_back_surface);

		//# % This normalizes the lens normal vectors to have unit magnitudes
		lens_normal_vectors = normalize(lens_normal_vectors);

		//# % If the refractive index is a constant double value, this directly
		//# % calculates the refractive index ratio, otherwise the ratio is
		//# % calculated as a function of the wavelength
		// (I DON'T KNOW HOW TO DO STRING EVALUATION IN CUDA. SO I AM ONLY GOING TO
		// CONSIDER THE FIRST CASE - LKR.)
		// TODO : implement eval
		//# % If the Abbe number is defined, this calculates the Cauchy formula
		//# % approxiation to the refractive index, otherwise, the refractive
		//# % index is defined to be constant
		if (!isnan(element_abbe_number)){
			//	# % This defines the three optical wavelengths used in defining
			//	# % the Abbe number
			lambda_D = 589.3;
			lambda_F = 486.1;
			lambda_C = 656.3;

			//	# % This is the ratio of the incident refractive index to the transmitted
			//	# % refractive index (1 since the  inter-lens media is assummed to be
			//	# % air)
			refractive_index_ratio = element_refractive_index + (1.0 / (ray_wavelength*ray_wavelength)
					- 1.0/(lambda_D*lambda_D)) * ((element_refractive_index - 1)/(element_abbe_number * (1 / (lambda_F*lambda_F) - 1 / (lambda_C*lambda_C))));
		}
		else{
			//# % This is the ratio of the incident refractive index to the transmitted
			//# % refractive index (1 since the  inter-lens media is assummed to be
			//# % air)
			refractive_index_ratio = element_refractive_index;
		}

		//# % This is the scaled cosine of the angle of the incident light ray
		//# % vectors and the normal vectors of the lens (ie the dot product of the
		//# % vectors)
		ray_dot_product = -dot(ray_propagation_direction,lens_normal_vectors);

		//# % This calculates the radicand in the refraction ray propagation
		//# % direction equation
		refraction_radicand = 1.0-(refractive_index_ratio*refractive_index_ratio)*(1.0-ray_dot_product*ray_dot_product);

		//# % Any rays that have complex values due to the radicand in the above
		//# % equation being negative experience total internal reflection.  The
		//# % values of these rays are set to NaN for now.
		if(refraction_radicand < 0.0){

			//# % This sets any of the light ray directions experiencing total
			//# % internal reflection to NaN values
			light_ray_data.ray_propagation_direction = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			//# % This sets any of the light ray origin coordinates experiencing total
			//# % internal reflection to NaN values
			light_ray_data.ray_source_coordinates = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			//# % This sets any of the light ray wavelengths experiencing total
			//# % internal reflection to NaN values
			light_ray_data.ray_wavelength = CUDART_NAN_F;

			//# % This sets any of the light ray radiances experiencing total internal
			//# % reflection to NaN values
			light_ray_data.ray_radiance = CUDART_NAN;

			return light_ray_data;
		}


		//# % This calculates the new light ray direction vectors (this is a
		//# % standard equation in optics relating incident and transmitted light
		//# % ray vectors)
		ray_propagation_direction = refractive_index_ratio*ray_propagation_direction +
				(refractive_index_ratio*ray_dot_product - sqrt(refraction_radicand))
				*lens_normal_vectors;

		//# % This normalizes the ray propagation direction so that it's magnitude
		//# % equals one
		ray_propagation_direction = normalize(ray_propagation_direction);

		//# % If the absorbance rate is non-zero, this calculates how much of the
		//# % radiance is absorbed by the lens, otherwise the output radiance is
		//# % just scaled by the transmission ratio
		if (element_absorbance_rate!=0){

			//	# % This calculates the distance that the rays traveled through
			//	# % the lens
			float propagation_distance = sqrt(dot(pos_intersect - ray_source_coordinates,pos_intersect - ray_source_coordinates));
			//# % If the absorbance rate is a simple constant, this
			//# % This calculates the new light ray radiance values
			ray_radiance = (1.0 - element_absorbance_rate) * ray_radiance * propagation_distance;

			// TODO: implement eval
		}

		else{
			//# % This rescales the output radiance by the transmission ratio
			ray_radiance = element_transmission_ratio * ray_radiance;
		}

		//# % This sets the new light ray origin to the intersection point with the
		//# % front surface of the lens
		ray_source_coordinates = pos_intersect;

	}

	//# % If the element type is 'aperture', this extracts the optical properties
	//# % required to perform the ray propagation
	else
	{

		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//# % Extraction of the aperture optical properties                       %
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		//# % This extracts the pitch of the aperture stop
		element_pitch = optical_element.element_geometry.pitch;

		//# % This is the thickness of the optical element from the front optical axis
		//# % vertex to the back optical axis vertex
		element_vertex_distance = optical_element.element_geometry.vertex_distance;
		//
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//# % propagation of the light rays through the aperture front surface    %
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//#
		//# % This extracts the individual plane parameters
		a = element_plane_parameters.x;
		b = element_plane_parameters.y;
		c = element_plane_parameters.z;
		d = element_plane_parameters.w;

		//# % This calculates the square of the plane normal vector magnitude
		norm_vector_magnitude = sqrt(a*a + b*b + c*c);

		//# % This is the current offset to the center of the front spherical
		//# % surface
		ds = -element_vertex_distance / 2.0;

		//# % This calculates the transformed plane parameters (only d changes)
		float d_temp = d - ds * norm_vector_magnitude;

		//# % This is the independent intersection time between the light rays and
		//# % the first plane of the aperture stop
		float intersection_time = -(dot(make_float3(a,b,c),ray_source_coordinates) + d_temp)
										/ dot(make_float3(a,b,c),ray_propagation_direction);

		//# % This calculates the intersection points
		pos_intersect = ray_source_coordinates + ray_propagation_direction*intersection_time;

		//# % This calculates how far the intersection point is from the optical
		//# % axis of the lens (for determining if the rays hit the lens outside of
		//# % the pitch and thus would be destroyed)
		optical_axis_distance = measure_distance_to_optical_axis(pos_intersect,element_center,element_plane_parameters);

		//# % This checks if the given light ray actually intersects the lens
		//# % (ie the rays that are within the lens pitch)
		if(optical_axis_distance > element_pitch/2.0){

			light_ray_data.ray_source_coordinates = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);
			//# % This sets any of the light ray directions outside of the domain of
			//# % the lens to NaN values
			light_ray_data.ray_propagation_direction = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			//# % This sets any of the light ray wavelengths outside of the  domain of
			//# % the lens to NaN values
			light_ray_data.ray_wavelength = CUDART_NAN_F;

			//# % This sets any of the light ray radiances outside of the  domain of
			//# % the lens to NaN values
			light_ray_data.ray_radiance = CUDART_NAN;

			//# % This sets the intersection points of any of the light rays outside of
			//# % the domain of the lens to NaN values
			pos_intersect = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			return light_ray_data;
		}

		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//# % propagation of the light rays through the aperture back surface     %
		//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		//# % This is the current offset to the center of the front spherical
		//# % surface
		ds = +element_vertex_distance / 2;

		//# % This calculates the transformed plane parameters (only d changes)
		d_temp = d - ds * norm_vector_magnitude;

		//# % This is the independent intersection time between the light rays and
		//# % the first plane of the aperture stop
		//intersection_time = -(
		//a * ray_source_coordinates[:,0] + b * ray_source_coordinates[:,1] + c * ray_source_coordinates[:,
		//	2] + d_temp) / (a * ray_propagation_direction[:,0] + b * ray_propagation_direction[:,1] + c *
		//					ray_propagation_direction[:,2])
		intersection_time = -(dot(make_float3(a,b,c),ray_source_coordinates) + d_temp) /
				dot(make_float3(a,b,c),ray_propagation_direction);

		//# % This calculates the intersection points
		pos_intersect = ray_source_coordinates + ray_propagation_direction * intersection_time;

		//# % This calculates how far the intersection point is from the optical
		//# % axis of the lens (for determining if the rays hit the lens outside of
		//# % the pitch and thus would be destroyed)
		optical_axis_distance = measure_distance_to_optical_axis(pos_intersect,element_center,element_plane_parameters);

		//# % This checks if the given light ray actually intersects the lens
		//# % (ie the rays that are within the lens pitch)
		if(optical_axis_distance > element_pitch/2.0){

			light_ray_data.ray_source_coordinates = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);
			//# % This sets any of the light ray directions outside of the domain of
			//# % the lens to NaN values
			light_ray_data.ray_propagation_direction = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			//# % This sets any of the light ray wavelengths outside of the  domain of
			//# % the lens to NaN values
			light_ray_data.ray_wavelength = CUDART_NAN_F;

			//# % This sets any of the light ray radiances outside of the  domain of
			//# % the lens to NaN values
			light_ray_data.ray_radiance = CUDART_NAN;

			//# % This sets the intersection points of any of the light rays outside of
			//# % the domain of the lens to NaN values
			pos_intersect = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

			return light_ray_data;
		}
	}

	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//# % Saving the light ray propagation data                                   %
	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//# % This extracts the propagation direction of the light rays
	light_ray_data.ray_propagation_direction = ray_propagation_direction;

	//# % This extracts the light ray source coordinates
	light_ray_data.ray_source_coordinates = ray_source_coordinates;

	//# % This extracts the wavelength of the light rays
	light_ray_data.ray_wavelength = ray_wavelength;

	//# % This extracts the light ray radiance
	light_ray_data.ray_radiance = ray_radiance;

	return light_ray_data;
}

__device__ void argsort(float* array, int num_elements, int* arg_array)
{
	/*
	 * 	This function implements the bubble sort routine, and returns an array containing
	 * 	the indices that would sort the original array
	 * 	(taken from http://www.programmingsimplified.com/c/source-code/c-program-bubble-sort)
	 */

	int c, d;
	int n = num_elements;
	float swap; int arg_swap;
	// initialize index array
	for(c = 0; c < n; c++)
		arg_array[c] = c;

	// perform sorting
	for (c = 0 ; c < ( n - 1 ); c++)
	{
		for (d = 0 ; d < n - c - 1; d++)
		{
		  if (array[d] > array[d+1]) /* For decreasing order use < */
		  {
			swap       = array[d];
			array[d]   = array[d+1];
			array[d+1] = swap;

			arg_swap       = arg_array[d];
			arg_array[d]   = arg_array[d+1];
			arg_array[d+1] = arg_swap;

		  }
		}
	}

}

__device__ light_ray_data_t propagate_rays_through_multiple_elements(element_data_t* optical_element, float3* element_center,
		float4* element_plane_parameters, int simultaneous_element_number, light_ray_data_t light_ray_data)
{
	//# % This function calculates the propagation of a set of light rays through
	//# % multiple optical elements.

	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//# % Extraction of the light ray propagation data                            %
	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//# % This extracts the propagation direction of the light rays
	float3 ray_propagation_direction = light_ray_data.ray_propagation_direction;
	//# % This extracts the light ray source coordinates
	float3 ray_source_coordinates = light_ray_data.ray_source_coordinates;
	//# % This extracts the wavelength of the light rays
	float ray_wavelength = light_ray_data.ray_wavelength;
	//# % This extracts the light ray radiance
	double ray_radiance = light_ray_data.ray_radiance;

	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//# % Determination of optical element light ray intersection matching        %
	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//# % This finds the unique set of planes defining the lens elements to
	//# % intersect the light rays with
	int num_unique_planes = 0;
	float4* unique_plane_element_parameters = (float4*) malloc(simultaneous_element_number*sizeof(float4));
	int k,l;

	// take each plane parameter set and scan through the whole array. if there is no match,
	// then add it to unique_plane array. else, continue
	int flag = 0;
	for(k = 0; k < simultaneous_element_number; k++)
	{
		for(l = 0; l < simultaneous_element_number; l++)
		{
			// if the indices are the same, then go to the next iteration
			if(k==l)
				continue;
			// if a duplicate is found, then update the flag and quit this loop
			if(element_plane_parameters[k]==element_plane_parameters[l])
			{
				flag++;
				break;
			}

		}

		// if flag is still zero, then this element is unique. add the plane parameters
		// to the unique array, and update the number of unique planes.
		if(flag==0)
		{
			unique_plane_element_parameters[num_unique_planes] = element_plane_parameters[k];
			num_unique_planes++;
		}

		// re-initialize flag to zero
		flag = 0;

	}

	//# % This is the number of unique planes that the light rays intersect
	int unique_plane_number = num_unique_planes;

	//# % This initializes a vector of intersection times to calculate the order in
	//# % which the light rays intersect the various optical elements
	float* unique_plane_intersection_time = (float *) malloc(unique_plane_number*sizeof(float));
	memset(unique_plane_intersection_time,0.0,unique_plane_number*sizeof(float));

	//# % This initializes the light ray plane intersection point arrays
	float3* pos_intersect_approximate = (float3 *) malloc(unique_plane_number * sizeof(float3));
	memset(pos_intersect_approximate,0,unique_plane_number * sizeof(float3));

	//# % This iterates through the unique planes and calculates when the light
	//# % rays hit each plane since the light rays must be sequentially propagated
	//# % through the optical elements, but the ordering of the elements may be a
	//# % bit more poorly defined for multiple elements (although in practice, this
	//# % probably won't happen very often)
	int plane_parameters_index;
	for(plane_parameters_index = 0; plane_parameters_index < unique_plane_number; plane_parameters_index++)
	{
		//# % This extracts the individual plane parameters
		float a = unique_plane_element_parameters[plane_parameters_index].x;
		float b = unique_plane_element_parameters[plane_parameters_index].y;
		float c = unique_plane_element_parameters[plane_parameters_index].z;
		float d = unique_plane_element_parameters[plane_parameters_index].w;

		//# % This is the independent intersection time between the light rays and
		//# % the current plane of optical elements
		unique_plane_intersection_time[plane_parameters_index] = -(dot(make_float3(a,b,c),ray_source_coordinates) + d)/
				dot(make_float3(a,b,c),ray_propagation_direction);

		//# % This calculates the intersection points of the light rays with the
		//# % optical element planes
		pos_intersect_approximate[plane_parameters_index] = ray_source_coordinates +
				ray_propagation_direction*unique_plane_intersection_time[plane_parameters_index];

	}

	//# % This sorts the plane intersection times to determine which planes should
	//# % be tested first
	int* unique_plane_index_order = (int *) malloc(unique_plane_number*sizeof(int));
	argsort(unique_plane_intersection_time,unique_plane_number,unique_plane_index_order);

	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//# % Light ray propagation through multiple elements                         %
	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//# % This initializes the structure to contain the light ray data for
	//# % the current optical element
	light_ray_data_t light_ray_data_temp;
	float3 pos_intersect_approximate_current;

	//# % This iterates through the different planes, calculating the location of
	//# % the optical elements that the light rays intersect

	// initialize array to hold optical_element_indices
	int* optical_element_indices = (int *) malloc(sizeof(int)*simultaneous_element_number);
	int num_optical_element_indices = 0;

	for(plane_parameters_index = 0; plane_parameters_index < unique_plane_number; plane_parameters_index++)
	{
		//# % These are the current intersection points
		pos_intersect_approximate_current = pos_intersect_approximate[unique_plane_index_order[plane_parameters_index]];

		// find the indices of elements that have the same plane parameters as the current
		// unique plane element
		for(k = 0; k < simultaneous_element_number; k++)
		{
			// if the plane parameters of the current unique plane element matches that of the
			// plane with index k, then add the index to the optical_element_indices array
			if(unique_plane_element_parameters[plane_parameters_index] == element_plane_parameters[k])
			{
//				// expand the array to store additional elements
//				if(num_optical_element_indices!=0)
//				{
//					optical_element_indices = (int *) realloc(optical_element_indices,num_optical_element_indices*sizeof(int));
//				}

				optical_element_indices[num_optical_element_indices] = k;
				num_optical_element_indices++;
			}
		}

		//# % These are the centers of the optical elements on the current plane
//		float3* pos_c_current = (float3 *) malloc(num_optical_element_indices * sizeof(float3));
		float3* pos_c_current = new float3[num_optical_element_indices];
		for(k = 0; k < num_optical_element_indices; k++)
			pos_c_current[k] = element_center[optical_element_indices[k]];

		//TODO: Implement knnsearch
		//# % This finds the closest lens center to each approximate intersection
		//# % point (there are a small number of cases where this may give an
		//# % incorrect match since the intersection points used here are based
		//# % upon the center points of the lenses and not the front surfaces, but
		//# % this should work well for the large majority of possible systems)
		//nbrs = NearestNeighbors(n_neighbors=1, metric='euclidean').fit([xc_current, yc_current, zc_current])
		//distances, indices = nbrs.kneighbors(
		//	[x_intersect_approxiate_current, y_intersect_approxiate_current, z_intersect_approxiate_current])
		//
		//nearest_neighbor_index = np.squeeze(indices)
		//# nearest_neighbor_index=knnsearch([xc_current,yc_current,zc_current],[x_intersect_approxiate_current,y_intersect_approxiate_current,z_intersect_approxiate_current],'K',1,'Distance','euclidean');
		//
		//# % This iterates through the individual optical elements extracting the
		//# % set of light rays that likely intersect the element and propagating
		//# % them through the element
		//for element_index in range(0, optical_element_indices.size):
		//	# % This is the set of indices into the light rays to propagate
		//	# % through the current optical element
		//	light_ray_indices = (nearest_neighbor_index == element_index)
		//
		//	# % This extracts the propagation direction of the light rays
		//	light_ray_data_temp['ray_propagation_direction'] = ray_propagation_direction[light_ray_indices][:]
		//	# % This extracts the light ray source coordinates
		//	light_ray_data_temp['ray_source_coordinates'] = ray_source_coordinates[light_ray_indices][:]
		//	# % This extracts the wavelength of the light rays
		//	light_ray_data_temp['ray_wavelength'] = ray_wavelength[light_ray_indices]
		//	# % This extracts the light ray radiance
		//	light_ray_data_temp['ray_radiance'] = ray_radiance[light_ray_indices]
		//
		//	# % This extracts the current optical element data
		//	current_optical_element = optical_element[optical_element_indices[element_index]]
		//	# % This extracts the current optical element plane parameters
		//	current_plane_parameters = element_plane_parameters[optical_element_indices[element_index]][:]
		//	# % This extracts the current center of the optical element
		//	current_element_center = element_center[optical_element_indices[element_index]][:]
		//
		//	# % This propagates the light rays through the single optical element
		//	light_ray_data_temp = propagate_rays_through_single_element(current_optical_element, current_element_center,
		//																current_plane_parameters, light_ray_data_temp)
		//
		//	# % This extracts the propagation direction of the light rays
		//	ray_propagation_direction[light_ray_indices][:] = light_ray_data_temp['ray_propagation_direction']
		//	# % This extracts the light ray source coordinates
		//	ray_source_coordinates[light_ray_indices][:] = light_ray_data_temp['ray_source_coordinates']
		//	# % This extracts the wavelength of the light rays
		//	ray_wavelength[light_ray_indices] = light_ray_data_temp['ray_wavelength']
		//	# % This extracts the light ray radiance
		//	ray_radiance[light_ray_indices] = light_ray_data_temp['ray_radiance']

		free(pos_c_current);
	}
	// free pointers
	free(unique_plane_intersection_time);
	free(unique_plane_element_parameters);
	free(pos_intersect_approximate);
	free(unique_plane_index_order);


	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//# % Saving the light ray propagation data                                   %
	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//# % This extracts the propagation direction of the light rays
	light_ray_data.ray_propagation_direction = ray_propagation_direction;
	//# % This extracts the light ray source coordinates
	light_ray_data.ray_source_coordinates = ray_source_coordinates;
	//# % This extracts the wavelength of the light rays
	light_ray_data.ray_wavelength = ray_wavelength;
	//# % This extracts the light ray radiance
	light_ray_data.ray_radiance = ray_radiance;

	return light_ray_data;
}

__device__ light_ray_data_t propagate_rays_through_optical_system(element_data_t* element_data, float3* element_center, float4* element_plane_parameters,
		int* element_system_index,int num_elements, int num_rays, int lightray_number_per_particle, light_ray_data_t light_ray_data)

{
	// % This function propagates the light ray data defined by the structure
	// % 'light_ray_data' through the optical system defined by the input
	// % arguments.

	int k;
	// % This is the number of sequential optical elements within the total
	// % optical system that the light rays must be iteratively passed through
	int sequential_element_number = 0;
	for(k = 0; k < num_elements; k++)
	{
		if(sequential_element_number<=element_system_index[k])
			sequential_element_number = element_system_index[k];
	}

	// % Since the the optical is defined from the sensor moving outward (i.e. the
	// % opposite direction in which the light will enter a camera system), this
	// % reverses the indexing of the optical elements so that the first element
	// % index corresponds to the first optical element that the light will hit

	// % This iterates through the sequential optical elements propagating the
	// % light rays through each successive element (or system of coplanar
	// % elements)
	int element_index;
	int* current_element_indices = (int *) malloc(num_elements*sizeof(int));
//	printf("sequential_element_number:%d\n",sequential_element_number);
	for(element_index = 0; element_index < sequential_element_number; element_index++)
	{
		// These are the indices of the current element or elements to propagate
		// the light rays through
//		current_element_indices = np.squeeze(np.argwhere(element_system_index == element_index))

		int element_ctr = 0;
		for(k = 0; k < num_elements; k++)
		{
//			if(element_system_index_local[k]==element_index)
			if(sequential_element_number - element_system_index[k]==element_index)
			{
				current_element_indices[k] = k;
				element_ctr++;
			}
		}

		// % This is the number of elements that the light rays need to be
		// % simultaneously propagated through
		int simultaneous_element_number = element_ctr;
//		printf("simultaneous_element_number: %d\n",simultaneous_element_number);


		// % If there is only a single element that the rays are to be propagated
		// % through, this propagates the rays through the single element;
		// % otherwise the light rays are simultaneously propagated through the
		// % multiple elements

		if(simultaneous_element_number == 1)
		{

//			// % This extracts the current optical element data
//			element_data_t current_optical_element_single = element_data[current_element_indices[0]];
//			// % This extracts the current optical element plane parameters
//			float4 current_plane_parameters_single = element_plane_parameters[current_element_indices[0]];
//			// % This extracts the current center of the optical element
//			float3 current_element_center_single = element_center[current_element_indices[0]];

			light_ray_data = propagate_rays_through_single_element(element_data[current_element_indices[0]], element_center[current_element_indices[0]],
																			   element_plane_parameters[current_element_indices[0]],light_ray_data);

		}
		else
		{

			element_data_t* current_optical_element = (element_data_t *) malloc(simultaneous_element_number*sizeof(element_data_t));
			float4* current_plane_parameters = (float4 *) malloc(simultaneous_element_number*sizeof(float4));
			float3* current_element_center = (float3 *) malloc(simultaneous_element_number*sizeof(float3));

			//# % This initializes a cell array to contain the optical element data
//			element_data_t current_optical_element[simultaneous_element_number];
			//# % This iterates through the individual optical elements extracting
			//# % the optical element data
			int simultaneous_element_index;
			for(simultaneous_element_index = 0; simultaneous_element_index <= simultaneous_element_number;simultaneous_element_index++)
			{
				// % This extracts the current optical element data
				current_optical_element[simultaneous_element_index] = element_data[
					current_element_indices[simultaneous_element_index]];
				//# % This extracts the current optical element plane parameters
				current_plane_parameters[simultaneous_element_index] = element_plane_parameters[current_element_indices[simultaneous_element_index]];
				//# % This extracts the current center of the optical element
				current_element_center[simultaneous_element_index] = element_center[current_element_indices[simultaneous_element_index]];

			}

			//# % This propagates the light rays through the multiple optical
			//# % elements
			light_ray_data = propagate_rays_through_multiple_elements(current_optical_element, current_element_center,
																	current_plane_parameters, simultaneous_element_number,light_ray_data);

			// free allocated memory
			free(current_optical_element);
			free(current_plane_parameters);
			free(current_element_center);

		}



	}

	// free allocated memory
	free(current_element_indices);
//	free(element_system_index_local);

	return light_ray_data;


}

__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                                          (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                        __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

//__device__ pixel_data_t intersect_sensor(light_ray_data_t light_ray_data,camera_design_t camera_design, int lightray_number_per_particle, int num_rays)
__device__ pixel_data_t intersect_sensor(light_ray_data_t light_ray_data,int lightray_number_per_particle, int num_rays)
{

	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//# % Propagation of the light rays to the sensor                         %
	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//# % This extracts the propagation direction of the light rays
	float3 ray_propagation_direction = light_ray_data.ray_propagation_direction;
	//# % This extracts the light ray source coordinates
	float3 ray_source_coordinates = light_ray_data.ray_source_coordinates;

	//# % This extracts the individual plane parameters for the sensor
	float a = 0.0;
	float b = 0.0;
	float c = 1.0;
//	float d = -camera_design.z_sensor;
	float d = -camera_design_const.z_sensor;

	//# % This is the independent intersection time between the light rays and
	//# % the first plane of the aperture stop
	float intersection_time = -(dot(make_float3(a,b,c),ray_source_coordinates) + d)/
			dot(make_float3(a,b,c),ray_propagation_direction);

	//# % This calculates the intersection points
	float3 pos_intersect = ray_source_coordinates + ray_propagation_direction * intersection_time;

	//# % This sets the new light ray origin to the intersection point with the
	//# % front surface of the lens
	ray_source_coordinates = pos_intersect;

	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//# % Add lightray radiance to sensor integration                         %
	//# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//# % This is the angle between the lightray and the sensor (ie a lightray
	//# % normal to the sensor would yield an angle of zero)
	//alpha = np.arctan(np.sqrt((ray_propagation_direction[:,0] / ray_propagation_direction[:,2]) ** 2 + (
	//    ray_propagation_direction[:,1] / ray_propagation_direction[:,2]) ** 2))
	float alpha = atan(sqrt((ray_propagation_direction.x/ray_propagation_direction.z)*(ray_propagation_direction.x/ray_propagation_direction.z)
			+ (ray_propagation_direction.y/ray_propagation_direction.z)*(ray_propagation_direction.y/ray_propagation_direction.z)));

	//# % This calculates the cos^4(alpha) term which controls the contribution
	//# % of the incident light rays onto the measured energy in the sensor
	double cos_4_alpha = cos(alpha)*cos(alpha)*cos(alpha)*cos(alpha);

	//# % This calculates the indices of the pixel on the sensor that the ray
	//# % intersects and the relative weighting between the pixels

	//# % This is the pixel pitch [micron]
	float pixel_pitch = camera_design_const.pixel_pitch;
	//# % This is the number of pixels in the x-direction
	int x_pixel_number = camera_design_const.x_pixel_number;
	//# % This is the number of pixels in the y-direction
	int y_pixel_number = camera_design_const.y_pixel_number;

	//# % This is the coordinate of pixel (1,1) [0][0]
	float pixel_1_x = -pixel_pitch * (x_pixel_number - 1) / 2.0;
	float pixel_1_y = -pixel_pitch * (y_pixel_number - 1) / 2.0;

	//# % This is the number of pixel diameters the point (x,y) is from the center
	//# % of the (0,0) pixel
	float x = ray_source_coordinates.x;
	float y = ray_source_coordinates.y;
	float d_x = (x - pixel_1_x) / pixel_pitch + 1.5;
	float d_y = (y - pixel_1_y) / pixel_pitch + 1.5;

	//# % These are the coordinates of the point that is 1/2 pixel less than the
	//# % center coordinate
	float d_y_lower = d_y - 0.5;
	float d_x_lower = d_x - 0.5;

	//# % These are the percentages of overlap for the upper right corner of the
	//# % pixel (actually this is the distance of overlap - but if the pixels have
	//# % dimensions of 1 x 1 then this is the same) in each direction
	double d_ii_ul = ceil(d_y_lower) - d_y_lower;
	double d_jj_ul = ceil(d_x_lower) - d_x_lower;
	//# % This is the area of overlap with the upper right pixel
	double w_ul = (d_ii_ul) * (d_jj_ul);

	//# % These are the percentages of overlap for the upper left corner of the
	//# % pixel (actually this is the distance of overlap - but if the pixels have
	//# % dimensions of 1 x 1 then this is the same) in each direction
	double d_ii_ur = ceil(d_y_lower) - d_y_lower;
	double d_jj_ur = 1 - d_jj_ul;

	//# % This is the area of overlap with the upper left pixel
	double w_ur = (d_ii_ur) * (d_jj_ur);

	//# % These are the percentages of overlap for the lower left corner of the
	//# % pixel (actually this is the distance of overlap - but if the pixels have
	//# % dimensions of 1 x 1 then this is the same) in each direction
	double d_ii_ll = 1 - d_ii_ul;
	double d_jj_ll = ceil(d_x_lower) - d_x_lower;

	//# % This is the area of overlap with the lower right pixel
	double w_ll = (d_ii_ll) * (d_jj_ll);

	//# % These are the percentages of overlap for the lower right corner of the
	//# % pixel (actually this is the distance of overlap - but if the pixels have
	//# % dimensions of 1 x 1 then this is the same) in each direction
	double d_ii_lr = 1 - d_ii_ul;
	double d_jj_lr = 1 - d_jj_ul;

	//# % This is the area of overlap with the lower left pixel
	double w_lr = (d_ii_lr) * (d_jj_lr);

	//# % These are the integral values of the upper left pixel
	int ii_ul = int(ceil(d_y_lower) - 1);
	int jj_ul = int(ceil(d_x_lower) - 1);

	//# % These are the integral values of the upper right pixel
	int ii_ur = int(ceil(d_y_lower) - 1);
	int jj_ur = int(jj_ul + 1);

	//# % These are the integral values of the lower left pixel
	int ii_ll = int(ii_ul + 1);
	int jj_ll = int(ceil(d_x_lower) - 1);

	//# % These are the integral values of the lower right pixel
	int ii_lr = int(ii_ul + 1);
	int jj_lr = int(jj_ul + 1);

	pixel_data_t pixel_data;
	//	//# % This is a vector of the ii coordinates
	pixel_data.ii_indices = make_int4(ii_ul,ii_ur,ii_ll,ii_lr);
	//	//# % This is a vector of the jj coordinates
	pixel_data.jj_indices = make_int4(jj_ul,jj_ur,jj_ll,jj_lr);
	//	//# % This is a vector of pixel weights
	pixel_data.pixel_weights = make_double4(w_ul, w_ur, w_ll, w_lr);
	pixel_data.cos_4_alpha = cos_4_alpha;

	return pixel_data;

}

__global__ void parallel_ray_tracing(float lens_pitch, float image_distance,
		scattering_data_t scattering_data, int scattering_type, lightfield_source_t lightfield_source,
		 int lightray_number_per_particle, int n_min, int n_max, float beam_wavelength,
		 float aperture_f_number, int num_rays,
		 float* rand_array_1,float* rand_array_2,element_data_t* element_data,
		 float3* element_center, float4* element_plane_parameters,
		int* element_system_index,int num_elements,
		double* image_array)
		//		camera_design_t* camera_design_p, double* image_array)

{
	//--------------------------------------------------------------------------------------
	// compute indices to access in lightfield_source and lightfield_data
	//--------------------------------------------------------------------------------------
	// find global thread ID
	int local_thread_id = threadIdx.x;

	// get id of particle which is the source of light rays
	int local_particle_id = blockIdx.y;

	// get id of ray emitted by the particle
	int local_ray_id = blockIdx.x*blockDim.x + local_thread_id;
	int global_ray_id = local_ray_id + local_particle_id*lightray_number_per_particle;

	// if the ray id is greater than the total number of rays to be simulated, exit.
	if(global_ray_id >= num_rays)
		return;
	// particle id in the light field source array
	int current_source_point_number = n_min + local_particle_id;

	// populate shared memory
	if(local_thread_id == 0)
	{
		/* since all the threads in a given block correspond to rays from the same source
		 * point, store the source point information in shared memory
		 */
		lightfield_source_shared.x = lightfield_source.x[current_source_point_number];
		lightfield_source_shared.y = lightfield_source.y[current_source_point_number];
		lightfield_source_shared.z = lightfield_source.z[current_source_point_number];
		lightfield_source_shared.radiance = lightfield_source.radiance[current_source_point_number];
		lightfield_source_shared.diameter_index = lightfield_source.diameter_index[current_source_point_number];

	}

	//--------------------------------------------------------------------------------------
	// generate light rays
	//--------------------------------------------------------------------------------------

	light_ray_data_t light_ray_data = generate_lightfield_angular_data(lens_pitch, image_distance,scattering_data,
				scattering_type, lightfield_source_shared,lightray_number_per_particle, n_min, n_max,
				beam_wavelength,aperture_f_number,num_rays,rand_array_1[local_ray_id],rand_array_2[local_ray_id]);
	//--------------------------------------------------------------------------------------
	// propagate rays through the optical system
	//--------------------------------------------------------------------------------------

	light_ray_data = propagate_rays_through_optical_system(element_data, element_center,
			element_plane_parameters,element_system_index,num_elements,num_rays,
			lightray_number_per_particle,light_ray_data);

	//--------------------------------------------------------------------------------------
	// intersect rays with sensor
	//--------------------------------------------------------------------------------------

//	camera_design_t camera_design = *camera_design_p;
//
//	// perform ray intersection with the sensor and the radiance integration on the gpu
//	pixel_data_t pixel_data = intersect_sensor(light_ray_data,camera_design,
//			lightray_number_per_particle,num_rays);

	// perform ray intersection with the sensor and the radiance integration on the gpu
	pixel_data_t pixel_data = intersect_sensor(light_ray_data,
			lightray_number_per_particle,num_rays);


//	//# % This is the number of pixels in the x-direction
//	int x_pixel_number = camera_design.x_pixel_number;
//	//# % This is the number of pixels in the y-direction
//	int y_pixel_number = camera_design.y_pixel_number;

	//# % This is the number of pixels in the x-direction
	int x_pixel_number = camera_design_const.x_pixel_number;
	//# % This is the number of pixels in the y-direction
	int y_pixel_number = camera_design_const.y_pixel_number;


	int image_index, k;
	double pixel_increment;
	int ii_indices[4] = {pixel_data.ii_indices.x,pixel_data.ii_indices.y,pixel_data.ii_indices.z,pixel_data.ii_indices.w};
	int jj_indices[4] = {pixel_data.jj_indices.x,pixel_data.jj_indices.y,pixel_data.jj_indices.z,pixel_data.jj_indices.w};
	double pixel_weights[4] = {pixel_data.pixel_weights.x,pixel_data.pixel_weights.y,pixel_data.pixel_weights.z,pixel_data.pixel_weights.w};
	double cos_4_alpha = pixel_data.cos_4_alpha;

	for(k = 0; k < 4; k++)
	{
		if(ii_indices[k]<0 || ii_indices[k]>=y_pixel_number || jj_indices[k]<0 || jj_indices[k]>=x_pixel_number)
			continue;

		image_index = (ii_indices[k]-1)*x_pixel_number + jj_indices[k]-1;
		pixel_increment = pixel_weights[k]*light_ray_data.ray_radiance*cos_4_alpha;

//		// (possible race condition)
//		image_array[image_index] += pixel_increment;

		// dumb way to do it
		atomicAdd(&image_array[image_index],pixel_increment);
	}

}

extern "C"{

void read_from_file()
{
	float lens_pitch, image_distance, beam_wavelength, aperture_f_number;
	scattering_data_t scattering_data;
	int scattering_type;
	lightfield_source_t lightfield_source;
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
		file_scattering.read ((char*)&scattering_data.inverse_rotation_matrix[k], sizeof(float));


	// beam_propagation_vector
	for(k = 0; k < 3; k++)
		file_scattering.read ((char*)&scattering_data.beam_propagation_vector[k], sizeof(float));

	// num_angles
	file_scattering.read ((char*)&scattering_data.num_angles, sizeof(int));

	// num_diameters
	file_scattering.read ((char*)&scattering_data.num_diameters, sizeof(int));

	// scattering_angle
	scattering_data.scattering_angle = (float *) malloc(scattering_data.num_angles*sizeof(float));
	for(k = 0; k < scattering_data.num_angles; k++)
			file_scattering.read ((char*)&scattering_data.scattering_angle[k], sizeof(float));

	// scattering_irradiance
	scattering_data.scattering_irradiance = (float *) malloc(scattering_data.num_angles * scattering_data.num_diameters*sizeof(float));
	for(k = 0; k < scattering_data.num_angles * scattering_data.num_diameters; k++)
			file_scattering.read ((char*)&scattering_data.scattering_irradiance[k], sizeof(float));

	file_scattering.close();

	//--------------------------------------------------------------------------------------
	// read lightfield_source data
	//--------------------------------------------------------------------------------------

	// open file
	filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/lightfield_source.bin";
	std::ifstream file_lightfield_source(filename.c_str(), std::ios::in |
			std::ios::binary);

	// lightray_number_per_particle
	file_lightfield_source.read ((char*)&lightfield_source.lightray_number_per_particle, sizeof(int));

	// lightray_process_number
	file_lightfield_source.read ((char*)&lightfield_source.lightray_process_number, sizeof(int));

	// num_particles
	file_lightfield_source.read ((char*)&lightfield_source.num_particles, sizeof(int));

	// num_rays
	file_lightfield_source.read ((char*)&lightfield_source.num_rays, sizeof(int));

	// diameter_index
	lightfield_source.diameter_index = (int *) malloc(lightfield_source.num_particles * sizeof(int));
	for(k = 0; k < lightfield_source.num_particles; k++)
		file_lightfield_source.read ((char*)&lightfield_source.diameter_index[k], sizeof(int));

	// radiance
	lightfield_source.radiance = (double *) malloc(lightfield_source.num_particles * sizeof(double));
	for(k = 0; k < lightfield_source.num_particles; k++)
		file_lightfield_source.read ((char*)&lightfield_source.radiance[k], sizeof(double));

	// x
	lightfield_source.x = (float *) malloc(lightfield_source.num_particles*sizeof(float));
	for(k = 0; k < lightfield_source.num_particles; k++)
		file_lightfield_source.read ((char*)&lightfield_source.x[k], sizeof(float));

	// y
	lightfield_source.y = (float *) malloc(lightfield_source.num_particles*sizeof(float));
	for(k = 0; k < lightfield_source.num_particles; k++)
		file_lightfield_source.read ((char*)&lightfield_source.y[k], sizeof(float));

	// z
	lightfield_source.z = (float *) malloc(lightfield_source.num_particles*sizeof(float));
	for(k = 0; k < lightfield_source.num_particles; k++)
		file_lightfield_source.read ((char*)&lightfield_source.z[k], sizeof(float));

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
		file_optical_elements.read((char*)&element_data[k].element_type,sizeof(char));

		// elements_coplanar
		file_optical_elements.read((char*)&element_data[k].elements_coplanar,sizeof(float));

		// rotation_angles
		for(l = 0; l < 3; l++)
			file_optical_elements.read((char*)&element_data[k].rotation_angles[l],sizeof(double));

		// z_inter_element_distance
		file_optical_elements.read((char*)&element_data[k].z_inter_element_distance,sizeof(float));

	}

//	std::string temp_string_1 = "-np.sqrt((100000.000000)**2-(x**2+y**2))";
//
//////	int len = strlen(temp_string);
//////	element_data[0].element_geometry.front_surface_shape = "-np.sqrt((100000.000000)**2-(x**2+y**2))";
////	//	element_data[0].element_geometry.front_surface_shape = (char *) malloc(strlen(temp_string)*sizeof(char));
////	strcpy(element_data[0].element_geometry.front_surface_shape,temp_string.c_str());//"-np.sqrt((100000.000000)**2-(x**2+y**2))");
//	element_data[0].element_geometry.front_surface_shape = &temp_string_1[0];
//	printf("front_surface_shape: %s\n",element_data[0].element_geometry.front_surface_shape);
//////	element_data[0].element_geometry.front_surface_shape[strlen(element_data[0].element_geometry.front_surface_shape)] = '\0';
////
//	std::string temp_string_2 = "+np.sqrt((100000.000000)**2-(x**2+y**2))";
////	strcpy(element_data[0].element_geometry.back_surface_shape,temp_string.c_str());//"+np.sqrt((100000.000000)**2-(x**2+y**2))");
//////	element_data[0].element_geometry.back_surface_shape = (char *) malloc(strlen(temp_string_2)*sizeof(char));
//////	strcpy(element_data[0].element_geometry.back_surface_shape,temp_string_2.c_str());
//	element_data[0].element_geometry.back_surface_shape = &temp_string_2[0];
//	printf("back_surface_shape: %s\n",element_data[0].element_geometry.back_surface_shape);
//
//	// element_type
//	temp_string = "lens";
////	strcpy(element_data[0].element_type,"lens");//"lens");
//	element_data[0].element_type = &temp_string[0];
//	printf("element_type: %s\n", element_data[0].element_type);

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

	//--------------------------------------------------------------------------------------
	// read camera_design data
	//--------------------------------------------------------------------------------------

	// declare structure to hold the data
	camera_design_t camera_design;

	// open file
	filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/camera_design.bin";
	std::ifstream file_camera_design(filename.c_str(), std::ios::in |
			std::ios::binary);
	printf("\n");

	// pixel bit depth
	file_camera_design.read((char*)&camera_design.pixel_bit_depth,sizeof(int));
//	printf("pixel bit depth: %d\n",camera_design.pixel_bit_depth);

	// pixel gain
	file_camera_design.read((char*)&camera_design.pixel_gain,sizeof(float));
//	printf("pixel gain: %f\n",camera_design.pixel_gain);

	// pixel pitch
	file_camera_design.read((char*)&camera_design.pixel_pitch,sizeof(float));
//	printf("pixel pitch: %f\n",camera_design.pixel_pitch);

	// x_camera_angle
	file_camera_design.read((char*)&camera_design.x_camera_angle,sizeof(float));
//	printf("x camera angle: %f\n",camera_design.x_camera_angle);

	// y_camera_angle
	file_camera_design.read((char*)&camera_design.y_camera_angle,sizeof(float));
//	printf("y camera angle: %f\n",camera_design.y_camera_angle);

	// x_pixel_number
	file_camera_design.read((char*)&camera_design.x_pixel_number,sizeof(int));
//	printf("x pixel number: %d\n",camera_design.x_pixel_number);

	// y_pixel_number
	file_camera_design.read((char*)&camera_design.y_pixel_number,sizeof(int));
//	printf("y pixel number: %d\n",camera_design.y_pixel_number);

	// z sensor
	file_camera_design.read((char*)&camera_design.z_sensor,sizeof(float));
//	printf("z sensor: %f\n",camera_design.z_sensor);

	file_camera_design.close();

	//--------------------------------------------------------------------------------------
	// read image array
	//--------------------------------------------------------------------------------------
	// open file
	filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/image_array.bin";
	std::ifstream file_image_array(filename.c_str(), std::ios::in |
			std::ios::binary);
	printf("\n");

	double* image_array = (double *) malloc(camera_design.x_pixel_number*camera_design.y_pixel_number*sizeof(double));

	// write image array to file
	for(k = 0; k < camera_design.x_pixel_number*camera_design.y_pixel_number; k++)
		file_image_array.read((char*)&image_array[k],sizeof(double));

	file_image_array.close();



	// call the ray tracing function
	start_ray_tracing(lens_pitch, image_distance,&scattering_data, scattering_type_str,&lightfield_source,
			lightray_number_per_particle,n_min, n_max,beam_wavelength,aperture_f_number,
			num_elements, element_center_p,element_data,element_plane_parameters_p,element_system_index,&camera_design,image_array);

}

void save_to_file(float lens_pitch, float image_distance,
		scattering_data_t* scattering_data_p, char* scattering_type_str,
		lightfield_source_t* lightfield_source_p, int lightray_number_per_particle,
		int n_min, int n_max,float beam_wavelength, float aperture_f_number,
		int num_elements, double (*element_center)[3],element_data_t* element_data_p,
		double (*element_plane_parameters)[4], int *element_system_index, camera_design_t* camera_design_p,
		double* image_array)
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
//		printf("front_surface_shape: %s\n", element_data_p[k].element_geometry.front_surface_shape);
		file_optical_elements.write((char*)&element_data_p[k].element_geometry.front_surface_spherical,sizeof(bool));
		printf("front_surface_spherical: %s\n",element_data_p[k].element_geometry.front_surface_spherical ? "true":"false");
		file_optical_elements.write((char*)&element_data_p[k].element_geometry.back_surface_radius,sizeof(float));
		printf("back_surface_radius: %f\n", element_data_p[k].element_geometry.back_surface_radius);
//		file_optical_elements.write((char*)&element_data_p[k].element_geometry.back_surface_shape,strlen(element_data_p->element_geometry.back_surface_shape)*sizeof(char));
//		printf("back_surface_shape: %s\n", element_data_p[k].element_geometry.back_surface_shape);
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
		printf("element_type: %c\n",element_data_p[k].element_type);
		file_optical_elements.write((char*)&element_data_p[k].element_type,sizeof(char));

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

	//--------------------------------------------------------------------------------------
	// save camera_design data
	//--------------------------------------------------------------------------------------
	// open file
	filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/camera_design.bin";
	std::ofstream file_camera_design(filename.c_str(), std::ios::out |
			std::ios::binary);
	printf("\n");

	// pixel bit depth
	file_camera_design.write((char*)&camera_design_p->pixel_bit_depth,sizeof(int));
	printf("pixel bit depth: %d\n",camera_design_p->pixel_bit_depth);

	// pixel gain
	file_camera_design.write((char*)&camera_design_p->pixel_gain,sizeof(float));
	printf("pixel gain: %f\n",camera_design_p->pixel_gain);

	// pixel pitch
	file_camera_design.write((char*)&camera_design_p->pixel_pitch,sizeof(float));
	printf("pixel pitch: %f\n",camera_design_p->pixel_pitch);

	// x_camera_angle
	file_camera_design.write((char*)&camera_design_p->x_camera_angle,sizeof(float));
	printf("x camera angle: %f\n",camera_design_p->x_camera_angle);

	// y_camera_angle
	file_camera_design.write((char*)&camera_design_p->y_camera_angle,sizeof(float));
	printf("y camera angle: %f\n",camera_design_p->y_camera_angle);

	// x_pixel_number
	file_camera_design.write((char*)&camera_design_p->x_pixel_number,sizeof(int));
	printf("x pixel number: %d\n",camera_design_p->x_pixel_number);

	// y_pixel_number
	file_camera_design.write((char*)&camera_design_p->y_pixel_number,sizeof(int));
	printf("y pixel number: %d\n",camera_design_p->y_pixel_number);

	// z sensor
	file_camera_design.write((char*)&camera_design_p->z_sensor,sizeof(float));
	printf("z sensor: %f\n",camera_design_p->z_sensor);

	file_camera_design.close();

	//--------------------------------------------------------------------------------------
	// save image array
	//--------------------------------------------------------------------------------------
	// open file
	filename = "/home/barracuda/a/lrajendr/Projects/parallel_ray_tracing/data/image_array.bin";
	std::ofstream file_image_array(filename.c_str(), std::ios::out |
			std::ios::binary);
	printf("\n");

	int num_pixels = camera_design_p->x_pixel_number * camera_design_p->y_pixel_number;

	printf("before intersecting rays with sensor\n");
	printf("image array\n");
	printf("first three: %f, %f, %f\n",image_array[0],image_array[1],image_array[2]);
	printf("last three: %f, %f, %f\n",image_array[num_pixels-3],image_array[num_pixels-2],image_array[num_pixels-1]);

	// write image array to file
	for(k = 0; k < camera_design_p->x_pixel_number*camera_design_p->y_pixel_number; k++)
		file_image_array.write((char*)&image_array[k],sizeof(double));

	file_image_array.close();
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
				double (*element_plane_parameters)[4], int *element_system_index,
				camera_design_t* camera_design_p, double* image_array)
{
	/*
	 * This function receives the ray tracing parameters from the python code, allocates
	 * memory on the gpu, specifies the thread and block information and calls the gpu
	 * functions that perform the various steps in the ray tracing process./
	 *
	 * INPUTS:
	 * lens_pitch - diameter of the lens
	 * image_distance - distance between the lens and the image location of the particle
	 * 					volume. calculated using the thin lens equation in the python code.
	 * scattering_data_p - pointer to structure containing mie scattering data
	 * scattering_type_str - "mie" for mie scattering or "diffuse" for uniform light radiance
	 * lightfield_source_p -pointer to structure containing info about light field generation
	 * lightray_number_per_particle - number of light rays to be generated for each particle
	 * n_min - id of the first particle in the set of particles that will be simulated in this call
	 * n_max - id of the last particle in the set of particles that will be simulated in this call
	 * beam_wavelength - wavelength of the laser beam in microns
	 * aperture_f_number - f# of the lens/aperture (= focal length / pitch)
	 *
	 * num_elements - number of optical elements that the light ray has to be propagated through
	 * element_center - array containing the 3D positions of all optical elements
	 * element_data_p - pointer to structure containing various properties of optical elements
	 * element_plane_parameters - parameters define the plane upon which the element is centered
	 * 							  (first three elements are a unit vector
	 * element_system_index -
	 *
	 *
	 *
	 *
	 */



	// create instance of structure using the pointers
	scattering_data_t scattering_data = *scattering_data_p;
	lightfield_source_t lightfield_source = *lightfield_source_p;

	// number of particles that will simulated in a single call
//	int source_point_number = n_max - n_min + 1;
	int source_point_number = 1000;
	// number of the light rays to be generated and traced in a single call
	int num_rays = source_point_number*lightray_number_per_particle;

	// counter variable for all the for loops in this function
	int k;
	//--------------------------------------------------------------------------------------
	// allocate space on GPU for lightfield_source
	//--------------------------------------------------------------------------------------

	// declare pointers to device arrays
	float* d_source_x = 0;
	float* d_source_y = 0;
	float* d_source_z = 0;
	double *d_source_radiance = 0;
	int *d_source_diameter_index = 0;
	int num_particles = lightfield_source.num_particles;
	printf("num_particles: %d\n",num_particles);

	// allocate space for device arrays on GPU
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
	float *d_scattering_angle = 0;
	float* d_scattering_irradiance = 0;

	// allocate space for device arrays on GPU
	cudaMalloc((void**)&d_scattering_angle,scattering_data.num_angles*sizeof(float));
	cudaMalloc((void**)&d_scattering_irradiance,scattering_data.num_angles*scattering_data.num_diameters*sizeof(float));

	// copy data to GPU
	cudaMemcpy(d_scattering_angle,scattering_data.scattering_angle,scattering_data.num_angles*sizeof(float)
	,cudaMemcpyHostToDevice);
	cudaMemcpy(d_scattering_irradiance,scattering_data.scattering_irradiance,scattering_data.num_angles*scattering_data.num_diameters*sizeof(float)
		,cudaMemcpyHostToDevice);

	// point host structure to device array
	scattering_data.scattering_angle = d_scattering_angle;

//	// setup texture array to hold scattering irradiance data
//	cudaMalloc((void**)&data_array,scattering_data.num_angles*scattering_data.num_diameters*sizeof(float));
//	cudaMemcpy(data_array,scattering_data.scattering_irradiance,
//				scattering_data.num_angles*scattering_data.num_diameters*sizeof(float),cudaMemcpyHostToDevice);
//	cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();
//	cudaBindTexture2D( NULL, mie_scattering_irradiance,
//	                               data_array,
//	                               desc, scattering_data.num_diameters, scattering_data.num_angles,
//	                               sizeof(float) * scattering_data.num_diameters );
	scattering_data.scattering_irradiance = d_scattering_irradiance;


	int scattering_type = 0;
	if(strcmp(scattering_type_str,"mie")==0)
		scattering_type = 1;

	//--------------------------------------------------------------------------------------
	// read random numbers from file
	//--------------------------------------------------------------------------------------

	// allocate space for cpu arrays to hold the random numbers
	float* h_rand1 = (float*) malloc(lightray_number_per_particle*sizeof(float));
	float* h_rand2 = (float*) malloc(lightray_number_per_particle*sizeof(float));

	// declare gpu arrays to store the random numbers
	float* d_rand1 = 0;
	float* d_rand2 = 0;

	// allocate space for GPU arrays to hold the random numbers
	cudaMalloc((void**)&d_rand1,sizeof(float)*lightray_number_per_particle);
	cudaMalloc((void**)&d_rand2,sizeof(float)*lightray_number_per_particle);

	// open files
	string filename = "/home/barracuda/a/lrajendr/Projects/camera_simulation/random1.bin";
	std::ifstream file_rand1(filename.c_str(),ios::in|ios::binary);

	filename = "/home/barracuda/a/lrajendr/Projects/camera_simulation/random2.bin";
	std::ifstream file_rand2(filename.c_str(),ios::in|ios::binary);

	// read random numbers from file
	for(k = 0; k < lightray_number_per_particle; k++)
	{
		file_rand1.read((char*)&h_rand1[k],sizeof(float));
		file_rand2.read((char*)&h_rand2[k],sizeof(float));
	}

	// check if data has been received properly
//	printf("random arrays\n");
//	printf("rand_1 1st, last : %f, %f\n",h_rand1[0],h_rand1[lightray_number_per_particle-1]);
//	printf("rand_2 1st, last : %f, %f\n",h_rand2[0],h_rand2[lightray_number_per_particle-1]);

	// close files
	file_rand1.close();
	file_rand2.close();

	// copy data to device
	cudaMemcpy(d_rand1,h_rand1,sizeof(float)*lightray_number_per_particle,cudaMemcpyHostToDevice);
	cudaMemcpy(d_rand2,h_rand2,sizeof(float)*lightray_number_per_particle,cudaMemcpyHostToDevice);

	//--------------------------------------------------------------------------------------
	// setup blocks and threads for ray tracing on the GPU
	//--------------------------------------------------------------------------------------

	// allocate number of threads per block
	dim3 block(512,1);

	// calculate number of blocks required
	int grid_x;
	// if number of light rays to be traced is less than or equal to the number of threads
	// in one block then just one block is enough
	if(lightray_number_per_particle<=block.x)
		grid_x = 1;
	// if number of light rays to be traced is greater than the number of threads in one block
	// and it is an exact multiple, then calculate the number of blocks required
	else if(lightray_number_per_particle>block.x && lightray_number_per_particle%block.x ==0)
		grid_x = lightray_number_per_particle/block.x;
	// if number of light rays to be traced is greater than the number of threads in one block
	// and it is NOT an exact multiple, then use an extra block
	else
		grid_x = lightray_number_per_particle/block.x + 1;

	/* allocate number of blocks per grid	 */
	dim3 grid(grid_x,source_point_number);
	printf("grid: %d, %d\n",grid_x,source_point_number);

	// print the number of rays
	printf("num_rays: %d\n",num_rays);

	//--------------------------------------------------------------------------------------
	// setup memory to propagate rays through the optical system
	//--------------------------------------------------------------------------------------

	//  convert coordinate arrays to float3 and float4 arrays
	float3* element_center_2 = (float3 *) malloc(num_elements*sizeof(float3));
	float4* element_plane_parameters_2 = (float4 *) malloc(num_elements*sizeof(float4));

	for(k = 0; k < num_elements; k++)
	{
		element_center_2[k] = make_float3(element_center[k][0],element_center[k][1],element_center[k][2]);
		element_plane_parameters_2[k] = make_float4(element_plane_parameters[k][0],element_plane_parameters[k][1],element_plane_parameters[k][2],element_plane_parameters[k][3]);
	}


	// declare arrays to hold element coordinate information the gpu
	float3* d_element_center = 0;
	float4* d_element_plane_parameters = 0;
	int* d_element_system_index = 0;

	// allocate space on GPU for coordinate arrays
	cudaMalloc((void**)&d_element_center,num_elements*sizeof(float3));
	cudaMalloc((void**)&d_element_plane_parameters,num_elements*sizeof(float4));
	cudaMalloc((void**)&d_element_system_index,num_elements*sizeof(int));

	// copy data from GPU to CPU
	cudaMemcpy(d_element_center,element_center_2,num_elements*sizeof(float3),cudaMemcpyHostToDevice);
	cudaMemcpy(d_element_plane_parameters,element_plane_parameters_2,num_elements*sizeof(float4),cudaMemcpyHostToDevice);
	cudaMemcpy(d_element_system_index,element_system_index,num_elements*sizeof(int),cudaMemcpyHostToDevice);

	element_data_t* d_element_data = 0;
	cudaMalloc((void**)&d_element_data,num_elements*sizeof(element_data_t));
	cudaMemcpy(d_element_data,element_data_p,num_elements*sizeof(element_data_t),cudaMemcpyHostToDevice);

	for(k = 0; k < num_elements; k++)
	{
//		cudaMalloc((void**)&d_element_data[k].axial_offset_distances,2*sizeof(double));
		cudaMemcpy(d_element_data[k].axial_offset_distances,element_data_p[k].axial_offset_distances,2*sizeof(double),cudaMemcpyHostToDevice);

//		cudaMalloc((void**)&d_element_data[k].rotation_angles,3*sizeof(double));
		cudaMemcpy(d_element_data[k].rotation_angles,element_data_p[k].rotation_angles,3*sizeof(double),cudaMemcpyHostToDevice);

	}

	//--------------------------------------------------------------------------------------
	// allocate memory for intersecting rays with sensor
	//--------------------------------------------------------------------------------------

	// declare pointer to the structure that will hold the camera design information
	camera_design_t* d_camera_design = 0;

	// allocate space on the GPU for the structure that will hold the camera design information
	cudaMalloc((void **)&d_camera_design, sizeof(camera_design_t));

	// copy camera design information to the gpu
	cudaMemcpy(d_camera_design,camera_design_p,sizeof(camera_design_t),cudaMemcpyHostToDevice);

	// declare pointer to the array that will hold the image
	double* d_image_array = 0;

	// allocate space on the gpu for the array that will hold the image
	cudaMalloc((void **)&d_image_array, sizeof(double)*camera_design_p->x_pixel_number*camera_design_p->y_pixel_number);

	// copy image array contents from the previous iteration to the gpu
	cudaMemcpy(d_image_array,image_array,sizeof(double)*camera_design_p->x_pixel_number*camera_design_p->y_pixel_number,cudaMemcpyHostToDevice);

	// number of pixels in the image
	int num_pixels = camera_design_p->x_pixel_number * camera_design_p->y_pixel_number;

	// copy camera design to constant memory
	cudaMemcpyToSymbol(camera_design_const, camera_design_p,sizeof(camera_design_t));

	clock_t begin, end;
	double time_spent;

	begin = clock();

	for(k = 0; k < num_particles/source_point_number; k++)
	{
		n_min = k*source_point_number;
		printf("k = %d\n",k);
		parallel_ray_tracing<<<grid,block>>>(lens_pitch, image_distance,scattering_data,
							scattering_type, lightfield_source,lightray_number_per_particle, n_min, n_max,
							beam_wavelength,aperture_f_number,num_rays,d_rand1,d_rand2,
							d_element_data, d_element_center,d_element_plane_parameters,
							d_element_system_index,num_elements,
							d_image_array);
		cudaDeviceSynchronize();
//		break;
	}

	end = clock();

	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

	printf("time spent in GPU: %f\n",time_spent);

	// copy image data to CPU
	cudaMemcpy(image_array,d_image_array,sizeof(double)*camera_design_p->x_pixel_number*camera_design_p->y_pixel_number,cudaMemcpyDeviceToHost);

	// display image array contents after intersecting the rays with the sensor
	printf("finished intersecting rays with sensor\n");
	printf("image array\n");
	printf("first three: %G, %G, %G\n",image_array[0],image_array[1],image_array[2]);
	printf("last three: %G, %G, %G\n",image_array[num_pixels-3],image_array[num_pixels-2],image_array[num_pixels-1]);

	printf("saving image array to file");
	filename = "/home/barracuda/a/lrajendr/Projects/camera_simulation/image_array.bin";
	std::ofstream file_image_array(filename.c_str(),ios::out | ios::binary);
	for(k = 0; k < num_pixels; k++)
		file_image_array.write((char*)&image_array[k],sizeof(double));

	file_image_array.close();

	// Destroy pointers to arrays
	free(h_rand1);
	free(h_rand2);
	cudaFree(d_rand1);
	cudaFree(d_rand2);

	// free GPU arrays
	// free pointers
	cudaFree(d_source_x);
	cudaFree(d_source_y);
	cudaFree(d_source_z);
	cudaFree(d_source_radiance);
	cudaFree(d_source_diameter_index);

	cudaFree(d_scattering_angle);
	cudaFree(d_scattering_irradiance);
	cudaFree(data_array);

	// free arrays
	cudaFree(d_element_center);
	cudaFree(d_element_plane_parameters);
	cudaFree(d_element_system_index);
	cudaFree(d_element_data);

	free(element_center_2);
	free(element_plane_parameters_2);

	cudaFree(d_camera_design);
	cudaFree(d_image_array);


}



}


