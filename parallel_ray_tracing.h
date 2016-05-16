/*
 * parallel_ray_tracing.h
 *
 *  Created on: Apr 20, 2016
 *      Author: lrajendr
 */

#ifndef PARALLEL_RAY_TRACING_H_
#define PARALLEL_RAY_TRACING_H_

// this structure holds scattering data information
typedef struct scattering_data_t
{
	// rotation matrix
	float inverse_rotation_matrix[9];
	// direction of propagation of the laser beam
	float beam_propagation_vector[3];
	// scattering angle at which Mie scattering info is available
	float* scattering_angle;
	// scattering irradiance values for each scattering angle
	float* scattering_irradiance;
	// number of angles between 0 and 180 degrees for which Mie scattering
	// info is available
	int num_angles;
	// number of particle diameters over which scattering data is available
	int num_diameters;

}scattering_data_t;

// this structure holds information about the light field generated by the source
typedef struct lightfield_source_t
{
	// number of rays to be generated for each light ray
	int lightray_number_per_particle;
	// number of light rays to be processed at a time (limited by available memory)
	int lightray_process_number;
	// this array contains the index corresponding to the particles diameter. this is used
	// to access the mie scattering irradiance for the given particle from the
	// scattering_irradiance array stored in the scattering_data structure
	int* diameter_index;
	// this is the radiance for a light ray generated from each particle
	double* radiance;
	// this is the x location of each particle
	float* x;
	// this is the y location of each particle
	float* y;
	// this is the z location of each particle
	float* z;
	// total number of rays (<=lightray_process_number)
	int num_rays;
	// number of particles to simulate
	int num_particles;
}lightfield_source_t;

// this structure contains the light field data for all the particles
typedef struct lightfield_data_t
{
	// x location of light ray
	float* x;
	// y location of light ray
	float* y;
	// z location of light ray
	float* z;
	// radiance of light ray
	double* radiance;
	// angle of propagation
	float* theta;
	// angle of propagation
	float* phi;
	// number of light rays
	int num_lightrays;
}lightfield_data_t;

// this structure contains the light field generated by a single particle
typedef struct lightfield_source_single_t
{
	// this contains the index corresponding to the particles diameter. this is used
	// to access the mie scattering irradiance for the given particle from the
	// scattering_irradiance array stored in the scattering_data structure
	int diameter_index;
	// this is the radiance for a light ray generated from each particle
	double radiance;
	// this is the x location of each particle
	float x;
	// this is the y location of each particle
	float y;
	// this is the z location of each particle
	float z;
}lightfield_source_single_t;


// this structure represents a single light ray and contains its properties
typedef struct light_ray_data_t
{
	// location of the light ray vector)
	float3 ray_source_coordinates;
	// direction of the light ray vector
	float3 ray_propagation_direction;
	// wavelength of the light ray
	float ray_wavelength;
	// radiance of the light ray
	double ray_radiance;
}light_ray_data_t;

// this structure contains the geometric properties of an optical element
typedef struct element_geometry_t
{
	// radius of the front surface of the lens/aperture
	float front_surface_radius;
	//char* front_surface_shape;
	// truth value that indicates whether the surface is spherical
	bool front_surface_spherical;
	// radius of the back surface of the lens/aperture
	float back_surface_radius;
	//char* back_surface_shape;
	// truth value that indicates whether the surface is spherical
	bool back_surface_spherical;
	// diameter of the element
	float pitch;
	// distance between the
	double vertex_distance;
}element_geometry_t;

// this structure contains the optical properties of an optical element
typedef struct element_properties_t
{
	float abbe_number;
	float absorbance_rate;
	double refractive_index;
	float thin_lens_focal_length;
	float transmission_ratio;
}element_properties_t;

// this structure contains all the information about a given optical element
typedef struct element_data_t
{
	// offset of this element from the optical axis
	double axial_offset_distances[2];

	element_geometry_t element_geometry;
	float element_number;
	element_properties_t element_properties;
	//char* element_type;
	char element_type;
	// Boolean value stating whether the current elements are co-planar
    // (this value may only be true for an element type of 'system' and will be
	// Null for single optical elements)
	float elements_coplanar;
	// the rotation angles for the current optical element
	double rotation_angles[3];
	// distance between the current element and the next element along the z axis
	float z_inter_element_distance;
}element_data_t;

// this structure contains the data that characterizes the camera
typedef struct camera_design_t
{
	// number of bits used to quantize and store the intensity
	int pixel_bit_depth;
	// conversion factor between the photon count of a pixel and its intensity
	float pixel_gain;
	// length of one side of a square pixel (microns)
	float pixel_pitch;
	// angle that the camera z makes with the global x axis
	float x_camera_angle;
	// angle that the camera z makes with the global y axis
	float y_camera_angle;
	// number of pixels along the camera x axis
	int x_pixel_number;
	// number of pixels along the camera y axis
	int y_pixel_number;
	// distance between the lens and the sensor in the camera coordinate system
	float z_sensor;
}camera_design_t;

struct pixel_data_t
{
	int4 ii_indices;
	int4 jj_indices;
	double4 pixel_weights;
	double cos_4_alpha;
};
__device__ void parallel_ray_tracing(float , float , scattering_data_t* , int ,
		lightfield_source_t* ,int , int , int , float , float , light_ray_data_t* , int ,
		 float* ,float* ,element_data_t* , float3* , float4* , int* ,int , camera_design_t* ,
		 double* );

__device__ light_ray_data_t generate_lightfield_angular_data(float ,float,scattering_data_t* ,
		int , lightfield_source_single_t , int , int , int, float, float,light_ray_data_t,int,
		float*, float*);

__device__ light_ray_data_t propagate_rays_through_optical_system(element_data_t*, float3* , float4* ,
		int* , int , int , int , light_ray_data_t);

__device__ pixel_data_t intersect_sensor(light_ray_data_t ,camera_design_t , int, int);

//__global__ void generate_lightfield_angular_data(float ,float,scattering_data_t* ,
//		int , lightfield_source_t* , int , int , int, float, float,light_ray_data_t*,int,
//		float*, float*);
//
//__global__ void propagate_rays_through_optical_system(element_data_t*, float3* , float4* ,
//		int* , int , int , int , light_ray_data_t*);
//
//__global__ void intersect_sensor(light_ray_data_t* ,camera_design_t* , double* , int, int);

__device__ float random_single(unsigned int );

__device__ light_ray_data_t propogate_rays_through_multiple_elements(element_data_t* , float3*,
		float4* , int, light_ray_data_t );
__device__ light_ray_data_t propagate_rays_through_single_element(element_data_t , float3,
		   float4, light_ray_data_t );
__device__ float3 ray_sphere_intersection(float3 , float , float3 , float3 , char );

__device__ float measure_distance_to_optical_axis(float3 , float3 , float );
__device__ void argsort(float* , int , int*);
__device__ double atomicAdd(double* , double);





extern "C"
{

void save_to_file(float , float ,scattering_data_t* , char* ,lightfield_source_t* ,
		int ,int , int, float, float,int , double (*element_center)[3],element_data_t* ,
		double (*element_plane_parameters)[4], int* ,camera_design_t* , double*);
void read_from_file();
int add(int a, int b);
void start_ray_tracing(float , float ,scattering_data_t* , char* ,lightfield_source_t* ,
		int ,int , int, float, float, int , double (*element_center)[3],element_data_t*,
		double (*element_plane_parameters)[4], int* ,camera_design_t* , double*);

}



#endif /* PARALLEL_RAY_TRACING_H_ */
