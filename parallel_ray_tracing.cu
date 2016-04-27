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
                                     float aperture_f_number, light_ray_data_t light_ray_data)
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
//	int block_id = blockIdx.x;
	int local_thread_id = threadIdx.x;
//	int global_thread_id = block_id + local_thread_id;

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
	float3 beam_propogation_vector;
	// if scattering_type is mie, then use mie scattering data
	if(scattering_type)
	{
		//% This extracts the normalized beam propagation direction vector from the
		//% parameters structure
		beam_propogation_vector.x=scattering_data.beam_propogation_vector[0];
		beam_propogation_vector.y=scattering_data.beam_propogation_vector[1];
		beam_propogation_vector.z=scattering_data.beam_propogation_vector[2];

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
		ray_scattering_angles = angleBetween(beam_propogation_vector,ray_direction_vector);
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
	light_ray_data.ray_source_coordinates[global_ray_id] = make_float3(x_current,y_current,z_current);
	light_ray_data.ray_propagation_direction[global_ray_id] = normalize(make_float3(theta_temp,phi_temp,1.0));
	light_ray_data.ray_wavelength[global_ray_id] = beam_wavelength;
	light_ray_data.ray_radiance[global_ray_id] = 1/(aperture_f_number*aperture_f_number)*irradiance_current;
//	d_lightfield_data.x[global_ray_id] = x_current;
//	d_lightfield_data.y[global_ray_id] = y_current;
//	d_lightfield_data.z[global_ray_id] = z_current;
//	d_lightfield_data.theta[global_ray_id] = theta_temp;
//	d_lightfield_data.phi[global_ray_id] = phi_temp;
//	d_lightfield_data.radiance[global_ray_id] = irradiance_current;

}


int main(int argc, char** argv)
{

	printf("Hello world\n");
//	float lens_pitch,image_distance;
//	scattering_data_t* scattering_data_p;
//	char* scattering_type_str;
//	lightfield_source_t* lightfield_source_p;
//	int lightray_number_per_particle, n_min, n_max;
//	lightfield_data_t* lightfield_data_p;
//
//	start_ray_tracing(lens_pitch,image_distance,scattering_data_p,scattering_type_str,
//			lightfield_source_p, lightray_number_per_particle,n_min,n_max,lightfield_data_p);

	return 0;
}

//int add(int a, int b)
//{
//	return a+b;
//}


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
	// save scalars to file
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


	// beam_propogation_vector
	for(k = 0; k < 3; k++)
		file_scattering.read ((char*)&scattering_data_p.beam_propogation_vector[k], sizeof(float));

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
	// save lightfield_source data to file
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

	start_ray_tracing(lens_pitch, image_distance,&scattering_data_p, scattering_type_str,&lightfield_source_p,lightray_number_per_particle,n_min, n_max,beam_wavelength,aperture_f_number);


}

void save_to_file(float lens_pitch, float image_distance,
		scattering_data_t* scattering_data_p, char* scattering_type_str,
		lightfield_source_t* lightfield_source_p, int lightray_number_per_particle,
		int n_min, int n_max,float beam_wavelength, float aperture_f_number)
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
	// beam_propogation_vector
	printf("beam propagation vector\n");
	for(k = 0; k < 3; k++){
		printf("%f ", scattering_data_p->beam_propogation_vector[k]);
		file_scattering.write ((char*)&scattering_data_p->beam_propogation_vector[k], sizeof(float));
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


}

int add(int a, int b)
{
	return a+b;
}

void start_ray_tracing(float lens_pitch, float image_distance,
		scattering_data_t* scattering_data_p, char* scattering_type_str,
		lightfield_source_t* lightfield_source_p, int lightray_number_per_particle,
		int n_min, int n_max,float beam_wavelength, float aperture_f_number)
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

	// allocate space for light_ray_data on CPU
	light_ray_data_t light_ray_data;
	light_ray_data.ray_source_coordinates = (float3 *) malloc(N*sizeof(float3));
	light_ray_data.ray_propagation_direction = (float3 *) malloc(N*sizeof(float3));
	light_ray_data.ray_wavelength = (float *) malloc(N*sizeof(float));
	light_ray_data.ray_radiance = (double *) malloc(N*sizeof(double));
	light_ray_data.num_lightrays = N;
	// declare pointers to GPU arrays
	float3 *d_ray_source_coordinates, *d_ray_propagation_direction;
	float *d_ray_wavelength;
	double *d_ray_radiance;

	// allocate memory on GPU
	cudaMalloc((void**)&d_ray_source_coordinates, N*sizeof(float3));
	cudaMalloc((void**)&d_ray_propagation_direction, N*sizeof(float3));
	cudaMalloc((void**)&d_ray_wavelength, N*sizeof(float));
	cudaMalloc((void**)&d_ray_radiance, N*sizeof(double));

	// initialize arrays to zero
	cudaMemset(d_ray_source_coordinates,0.0,N*sizeof(float3));
	cudaMemset(d_ray_propagation_direction,0.0,N*sizeof(float3));
	cudaMemset(d_ray_wavelength,0.0,N*sizeof(float));
	cudaMemset(d_ray_radiance,0.0,N*sizeof(double));

	// copy contents of light_ray_data structure
	light_ray_data_t light_ray_data_copy = light_ray_data;
	// point structure to device arrays
	light_ray_data.ray_source_coordinates = d_ray_source_coordinates;
	light_ray_data.ray_propagation_direction = d_ray_propagation_direction;
	light_ray_data.ray_wavelength = d_ray_wavelength;
	light_ray_data.ray_radiance = d_ray_radiance;

	int scattering_type = 0;
	if(strcmp(scattering_type_str,"mie")==0)
		scattering_type = 1;


	// allocate threads per block
	dim3 block(10,1,1);
	// allocate blocks per grid
	dim3 grid(source_point_number/block.x,1,lightray_number_per_particle);

	// call kernel
	generate_lightfield_angular_data<<<grid,block>>>(lens_pitch, image_distance,scattering_data,
			scattering_type, lightfield_source,lightray_number_per_particle, n_min, n_max,
			beam_wavelength,aperture_f_number,light_ray_data);

	cudaThreadSynchronize();

	// copy light_ray_data back to host
	cudaMemcpy(light_ray_data_copy.ray_source_coordinates,light_ray_data.ray_source_coordinates,N*sizeof(float3),cudaMemcpyDeviceToHost);
	cudaMemcpy(light_ray_data_copy.ray_propagation_direction,light_ray_data.ray_propagation_direction,N*sizeof(float3),cudaMemcpyDeviceToHost);
	cudaMemcpy(light_ray_data_copy.ray_wavelength,light_ray_data.ray_wavelength,N*sizeof(float),cudaMemcpyDeviceToHost);
	cudaMemcpy(light_ray_data_copy.ray_radiance,light_ray_data.ray_radiance,N*sizeof(double),cudaMemcpyDeviceToHost);

	cudaThreadSynchronize();

	// point original light_ray_data to updated host arrays
	light_ray_data.ray_source_coordinates = light_ray_data_copy.ray_source_coordinates;
	light_ray_data.ray_propagation_direction = light_ray_data_copy.ray_propagation_direction;
	light_ray_data.ray_wavelength = light_ray_data_copy.ray_wavelength;
	light_ray_data.ray_radiance = light_ray_data_copy.ray_radiance;

	// display first and last few elements of lightfield_data
	N = light_ray_data.num_lightrays;
	printf("lightfield_data contents\n");
	printf("ray_source_coordinates (1st): %f, %f, %f\n",light_ray_data.ray_source_coordinates[0].x,light_ray_data.ray_source_coordinates[0].y,light_ray_data.ray_source_coordinates[0].z);
	printf("ray_source_coordinates (last): %f, %f, %f\n",light_ray_data.ray_source_coordinates[N-1].x,light_ray_data.ray_source_coordinates[N-1].y,light_ray_data.ray_source_coordinates[N-1].z);
	printf("ray_propagation_direction (1st): %f, %f, %f\n",light_ray_data.ray_propagation_direction[0].x,light_ray_data.ray_propagation_direction[0].y,light_ray_data.ray_propagation_direction[0].z);
	printf("ray_propagation_direction (last): %f, %f, %f\n",light_ray_data.ray_propagation_direction[N-1].x,light_ray_data.ray_propagation_direction[N-1].y,light_ray_data.ray_propagation_direction[N-1].z);
	printf("ray_wavelength (1st, last): %f, %f\n",light_ray_data.ray_wavelength[0],light_ray_data.ray_wavelength[N-1]);
	printf("ray_radiance (1st, last): %f, %f\n",light_ray_data.ray_radiance[0],light_ray_data.ray_radiance[N-1]);


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

	cudaFree(d_ray_source_coordinates);
	cudaFree(d_ray_propagation_direction);
	cudaFree(d_ray_wavelength);
	cudaFree(d_ray_radiance);

}



}


