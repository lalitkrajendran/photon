/*
 * trace_rays_through_density_gradients.h
 *
 * This file contains the functions that are involved in simulating ray deflection in
 * a variable density environment
 *
 * These functions are taken directly from SchlierenRay (https://github.com/TACC/SchlierenRay/tree/master)
 * with some modifications to suit the current program
 *
 *  Created on: May 20, 2016
 *      Author: lrajendr
 */

#ifndef TRACE_RAYS_THROUGH_DENSITY_GRADIENTS_H_
#define TRACE_RAYS_THROUGH_DENSITY_GRADIENTS_H_

/*
  For more information, please see: http://software.sci.utah.edu

  The MIT License

  Copyright (c) 2012-2013
  Scientific Computing and Imaging Institute, University of Utah

  License for the specific language governing rights and limitations under
  Permission is hereby granted, free of charge, to any person obtaining a
  copy of this software and associated documentation files (the "Software"),
  to deal in the Software without restriction, including without limitation
  the rights to use, copy, modify, merge, publish, distribute, sublicense,
  and/or sell copies of the Software, and to permit persons to whom the
  Software is furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included
  in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
  DEALINGS IN THE SOFTWARE.
*/

#include <vector>
#include <teem/nrrd.h>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <float.h>
#include <assert.h>
#include <cstring>
#include <stdint.h>
#include <cuda_runtime.h>
#include <cuda.h>

#include "float3_operators.h"
#include "parallel_ray_tracing.h"
#include "helper_cuda.h"

#define CUDART_NAN_F            __int_as_float(0x7fffffff)
#define CUDART_NAN              __longlong_as_double(0xfff8000000000000ULL)

cudaArray* data_array = 0;

using namespace std;
// this is a texture that will contain the refractive index gradient data
texture<float4, 3> tex_data;

typedef struct {

  // this array will be used to store the density data
  float* data;
  // this is the name of the file containing the density data
  char* filename;
  // these are the number of data points along x, y and z
  int sizex, sizey, sizez;
  // these are the minimum and maximum coordinates of the volume containing the density
  // data
  double xmin, xmax, ymin, ymax, zmin, zmax;
  // these are the grid spacings for the volume containing the density data
  double del_x, del_y, del_z;

} DataFile;

__device__ bool IntersectWithVolume(float3* ray_pos, float3 ray_dir, float3 p1, float3 p2)
{
	/*
	 * This function checks whether a light ray given its initial position and direction
	 * will strike the volume containing the density gradients
	 *
	 * INPUTS:
	 * ray_pos - pointer to vector containing the position of the light ray
	 * ray_dir - direction of propagation of the light ray
	 * p1 - corner of the volume containing the minimum value of all three coordinates
	 * p2 - corner of the volume containing the maximum value of all three coordinates
	 *
	 * OUTPUT:
	 * truth value - True/False on whether the ray intersects the volume or not
	 * ray_pos - updated ray position if the ray does strike the volume
	 *
	 */


	float d1 = p1.x;
	float d2 = p2.x;
	float tnear =  - (FLT_MAX - 1);
	float tfar = FLT_MAX;
	float t1= (d1 - ray_pos->x)/(ray_dir.x);
	float t2 = (d2 - ray_pos->x)/(ray_dir.x);
	if (t1 >t2)
	{
		float temp = t1;
		t1 = t2;
		t2 = temp;
	}

	if (t1 > tnear)
		tnear = t1;
	if (t2 < tfar)
		tfar = t2;
	if (tnear > tfar) //miss
		return false;
	if (tfar < 0.0) // box behind ray
		return false;

	t1 = (p1.y - ray_pos->y)/(ray_dir.y);
	t2 = (p2.y - ray_pos->y)/(ray_dir.y);
	if(t1 >t2)
	{
		float temp = t1;
		t1 = t2;
		t2 = temp;
	}
	if (t1 > tnear)
		tnear = t1;
	if (t2 < tfar)
		tfar = t2;
	if (tnear > tfar) //miss
		return false;
	if (tfar < 0.0) // box behind ray
		return false;

	t1 = (p1.z - ray_pos->z)/(ray_dir.z);
	t2 = (p2.z - ray_pos->z)/(ray_dir.z);
	if (t1 >t2)
	{
		float temp = t1;
		t1 = t2;
		t2 = temp;
	}

	float t;
	if (t1 >= 0 && t1 > tnear)
		tnear = t1;
	if (t2 < tfar)
		tfar = t2;
	if (tnear > tfar) //miss
		return false;
	if (tfar < 0.0) // box behind ray
		return false;
	else if (tnear < 0)  //I dunnno if this is right... put this in for rays starting in box
		t = tfar;
	else
		t = tnear;

	ray_pos->x += ray_dir.x*t;
	ray_pos->y += ray_dir.y*t;
	ray_pos->z += ray_dir.z*t;

	return true;
}

__device__ float3 calculate_lookup_index(float3 pos, density_grad_params_t params, float3 lookup_scale)
{
	// calculate distance between the ray location and the minimum corner of the volume
	float3 offset = pos - params.min_bound;

	float3 lookupfn;

	// calculate the normalized lookup index corresponding to the current ray position
	lookupfn.x = lookup_scale.x*offset.x;
	lookupfn.y = lookup_scale.y*offset.y;
	lookupfn.z = lookup_scale.z*offset.z;

	// calculate the lookup index
	float3 lookup = {static_cast<float>(lookupfn.x*params.data_width), static_cast<float>(lookupfn.y*params.data_height), static_cast<float>(lookupfn.z*params.data_depth)};

	return lookup;
}

__device__ bool ray_inside_box(float3 pos, density_grad_params_t params,
		float3 lookup)
{
	/*
	 * check if light ray is inside the density gradient volume.
	 * returns true if it is or false otherwise
	 */


	// if the lookup index lies outside the volume, exit the loop
	if(pos.x <= params.min_bound.x || pos.y <= params.min_bound.y || pos.z <= params.min_bound.z ||
			pos.x >= params.max_bound.x || pos.y >= params.max_bound.y || pos.z >= params.max_bound.z )
		return false;

	if(lookup.x <= 0 || lookup.y <= 0 || lookup.z <= 0 ||
			lookup.x >= params.data_width-1 || lookup.y >= params.data_height-1 || lookup.z >= params.data_depth-1 )
		return false;

	return true;

}

__device__ light_ray_data_t trace_rays_through_density_gradients(light_ray_data_t light_ray_data, density_grad_params_t params)
{

	// this is the corner of the volume containing the minimum coordinates
	float3 min_bound = params.min_bound;

	// this is the corner of the volume containing the maximum coordinates
	float3 max_bound = params.max_bound;


	// this is the scaling factor to covert the coordinates to integer index locations
	float3 lookup_scale = {1.0f/(max_bound.x-min_bound.x), 1.0f/(max_bound.y - min_bound.y), 1.0f/(max_bound.z-min_bound.z)};

 	float max_scale = max(max(float(params.data_width), float(params.data_height)), float(params.data_depth));

 	// this is the initial ray position
  	float3 pos = light_ray_data.ray_source_coordinates;
 	// this is the initial ray direction
  	float3 dir = light_ray_data.ray_propagation_direction;

 	// if initial position is outside the volume, intersect ray with volume
 	if(pos.x <= min_bound.x || pos.y <= min_bound.y || pos.z <= min_bound.z ||
 	pos.x >= max_bound.x || pos.y >= max_bound.y || pos.z >= max_bound.z )
 	{
 		// if the ray did not intersect the volume, then set all the properties to NAN
 		// and exit the function
 		if(!IntersectWithVolume(&pos, dir, params.min_bound, params.max_bound))
 		{
 			//# % This sets any of the light ray positions outside of the domain
 			//# % to NaN values
 			light_ray_data.ray_source_coordinates = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

 			//# % This sets any of the light ray directions outside of the domain
 			//# % to NaN values
 			light_ray_data.ray_propagation_direction = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

 			return light_ray_data;
 		}
 	}

 	int i = 0;
 	int insideBox = 1;
 	float refractive_index = 1.000277;
 	float3 lookup;
 	float4 val;

 	float3 A, B, C, D;
 	float3 R_n, T_n;
 	float delta_t;
 	//--------------------------------------------------------------------------------------
 	// trace light ray through the variable density medium
 	//--------------------------------------------------------------------------------------
 	while(insideBox==1)
 	{
 		// update loop index
 		i = i+1;

 		/************ Update ray position using RK4 method *********/
		/*
 		 * This function updates the position and direction of a light ray in a refractive index gradient
 		 * field by using a discretized version of fermat's equation using an RK4 method.
 		 * Taken from: Sharma et. al., Applied Optics (1982). Variables are identical to those in the paper.
 		 *
 		 * light_ray_data - structure containing the light ray properties
 		 * params - structure containing the density gradient properties
 		 * lookup_scale - used to calculate the array index in the density gradient texture from the
 		 * 			 absolute position of the light ray
 		 *
 		 *
 		 */

 		pos = light_ray_data.ray_source_coordinates;
 		dir = light_ray_data.ray_propagation_direction;
 		if(i==1)
 			pos = pos + dir * params.step_size/refractive_index;

 		// calculate lookup index to access refractive index gradient value
		lookup = calculate_lookup_index(pos, params, lookup_scale);
		// check if ray is inside volume
		if(!ray_inside_box(pos, params, lookup))
			return light_ray_data;
		// extract refractive index gradients as well as the local refractive index
		val = tex3D(tex_data, lookup.x, lookup.y, lookup.z);

		// get light ray position
 		R_n = pos;
 		// step size used in the RK4 integration (val.w is the refractive index)
 		delta_t = params.step_size/val.w;
 		// calculate optical ray direction vector
 		T_n = val.w * light_ray_data.ray_propagation_direction;

 		// calculate coefficients
 		D = make_float3(val.w * val.x, val.w * val.y, val.w * val.z);
 		A = delta_t * D;

 		pos = R_n + delta_t/2.0 * T_n + 1/8.0 * delta_t * A;
 		lookup = calculate_lookup_index(pos, params, lookup_scale);
 		// check if ray is inside volume
		if(!ray_inside_box(pos, params, lookup))
			return light_ray_data;

		val = tex3D(tex_data, lookup.x, lookup.y, lookup.z);
 		D = make_float3(val.w * val.x, val.w * val.y, val.w * val.z);
 		B = delta_t * D;

 		pos = R_n + delta_t * T_n + 1/2.0 * delta_t * B;
 		lookup = calculate_lookup_index(pos, params, lookup_scale);
 		// check if ray is inside volume
		if(!ray_inside_box(pos, params, lookup))
			return light_ray_data;
 		val = tex3D(tex_data, lookup.x, lookup.y, lookup.z);
 		D = make_float3(val.w * val.x, val.w * val.y, val.w * val.z);
 		C = delta_t * D;

 		// calculate new positions and directions
 		R_n = R_n + delta_t * (T_n + 1/6.0 * (A + 2*B));
// 		float3 T_n_increment = 1/6.
 		T_n = T_n + 1/6.0 * (A + 4*B + C);

 		// store the new position and direction in the light ray vector
 		light_ray_data.ray_source_coordinates = R_n;
 		light_ray_data.ray_propagation_direction = normalize(T_n/val.w);


// 		// retrieve the refractive index gradient at the given location
// 		val = tex3D(tex_data, round(lookup.x), round(lookup.y), round(lookup.z)); //*params.dataScalar;

// 		// calculate the change in ray direction
// 		normal = make_float3(-val.x,-val.y,val.z);
// 		// update the ray direction
// 		dir = dir + params.step_size*normal;

 		// find the change in the ray trajectory
// 		val = eulerIntegrate(pos, params, lookup_scale);
// 		val = rk4Integrate(pos, params, lookup_scale);
//
// 		// update the ray direction
// 		dir = dir + make_float3(-val.x, -val.y, val.z);

// 		// get the refractive index at the current location
// 		refractive_index = val.w;
//
// 		// normalize the direction to ensure that it is a unit vector
// 		dir = normalize(dir);
 		// update the ray position
//	    pos = pos + dir*params.step_size/refractive_index;
//	    pos = pos + dir*params.step_size;

   }
//
// 	// this is the final location of the ray
// 	light_ray_data.ray_source_coordinates = pos;
// 	// this is the final direction of the ray
// 	light_ray_data.ray_propagation_direction = normalize(dir);

// 	printf("refractive_index: %f\n", refractive_index);

 	return light_ray_data;
}

extern "C"{

void Host_Init(density_grad_params_t* paramsp, density_grad_params_t* dparams)
{
	/*
	 * This function sets up a texture lookup to the GPU array containing the density
	 * gradient data
	 */

	printf("Setting up data texture\n");

	//setup data texture
	tex_data.addressMode[0] = cudaAddressModeClamp;
	tex_data.addressMode[1] = cudaAddressModeClamp;
	tex_data.filterMode = cudaFilterModeLinear;
	tex_data.normalized = false;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
	cudaExtent extent = make_cudaExtent(paramsp->data_width, paramsp->data_height, paramsp->data_depth);
	checkCudaErrors( cudaMalloc3DArray(&data_array, &channelDesc, extent) );
	cudaMemcpy3DParms copyParams = {0};
	copyParams.srcPtr = make_cudaPitchedPtr((void*)paramsp->data, extent.width*sizeof(float4), extent.width, extent.height);
	copyParams.dstArray = data_array;
	copyParams.extent = extent;
	copyParams.kind = cudaMemcpyHostToDevice;
	checkCudaErrors(  cudaMemcpy3D(&copyParams) );

	cudaBindTextureToArray(tex_data, data_array, channelDesc);

}

void loadNRRD(DataFile* datafile, int data_min, int data_max)
{

	/*
	 * This function reads the density data from the NRRD file in addition to the
	 * variables containing information about the extent of the volume, and the grid
	 * spacing. It also converts the density to refractive index.
	 * Information about the nrrd file format is available at :
	 * http://teem.sourceforge.net/nrrd/lib.html
	 */

	printf("loading file %s : ", datafile->filename);
	// create pointer to an object representing the nrrd file
	Nrrd* nrrd = nrrdNew();

	// if the file does not exist, print error and exit
	if(nrrdLoad(nrrd, datafile->filename, 0))
	{
		char* err=biffGetDone(NRRD);
		cerr << "Failed to open \"" + string(datafile->filename) + "\":  " + string(err) << endl;
		exit(__LINE__);
	}

	// obtain number of grid points along each axis
	int sizex, sizey, sizez;
	sizex = nrrd->axis[0].size;
	sizey = nrrd->axis[1].size;
	sizez = nrrd->axis[2].size;

	// get min, max and grid spacing for each axis
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double del_x, del_y, del_z;

	xmin = nrrd->spaceOrigin[0];
	del_x = nrrd->axis[0].spacing;
	xmax = xmin + (sizex-1)*del_x;

	ymin = nrrd->spaceOrigin[1];
	del_y = nrrd->axis[1].spacing;
	ymax = ymin + (sizey-1)*del_y;

	zmin = nrrd->spaceOrigin[2];
	del_z = nrrd->axis[2].spacing;
	zmax = zmin + (sizez-1)*del_z;

	printf("\n******************** Co-ordinate System Info ******************************\n");
	printf("xmin: %g, xmax: %g, N_x: %d, del_x: %f\n", xmin,xmax,sizex,del_x);
	printf("ymin: %g, ymax: %g, N_y: %d, del_y: %f\n", ymin,ymax,sizey,del_y);
	printf("zmin: %g, zmax: %g, N_z: %d, del_z: %f\n", zmin,zmax,sizez,del_z);


	// not sure what these statements do
	if (data_max > sizez)
	data_max = sizez;
	if (sizez > (data_max-data_min))
	sizez = data_max-data_min;

	printf(" size: %f %f %f ", float(sizex), float(sizey), float(sizez));

	// initialize data array and max and min variables
	float* data = new float[sizex*sizey*sizez];
	float min = FLT_MAX;
	float max = -FLT_MAX;
	float* dataNrrd = (float*)nrrd->data;
	float* datai = data;
	float data_fudge = 1.0;
	// set GladStone Dale constant (cm^3/g) for refractive index calculation
	float K = 0.226;

	for(int i = 0; i < sizex; i++)
	{
		for(int j = 0; j < sizey; j++)
		{
			for( int k = 0; k < sizez; k++)
			{
				// read density data from file
				*datai = (*dataNrrd)*data_fudge;

				// convert density to refractive index
				*datai = 1.0 + K/1000.0*(*datai);

				// update max and min values
				if (*datai > max)
				  max = *datai;
				if (*datai < min)
				  min = *datai;

				datai++;
				dataNrrd++;

			}
		}
	}

	// transfer data to pointer
	datafile->data = data;
	datafile->sizex = sizex;
	datafile->sizey = sizey;
	datafile->sizez = sizez;

	datafile->xmin = xmin;
	datafile->xmax = xmax;
	datafile->ymin = ymin;
	datafile->ymax = ymax;
	datafile->zmin = zmin;
	datafile->zmax = zmax;
	datafile->del_x = del_x;
	datafile->del_y = del_y;
	datafile->del_z = del_z;

	// close file
	nrrdNuke(nrrd);
	printf("  ...done\n");
	printf("Min: %f, Max: %f \n",min, max);
}
}


density_grad_params_t setData(float* data, density_grad_params_t _params)
{

	/*
	 * This function calculates the refractive index gradient for all the grid points where
	 * density data is available
	 */

	int data_width = _params.data_width;
	int data_height = _params.data_height;
	int data_depth = _params.data_depth;

	printf("setData(%d, %d, %d, %d)\n", (u_int64_t)data, data_width, data_height, data_depth);
	_params.data_min = FLT_MAX;

	// this is the array containing the refractive index gradients
	_params.data = new float4[data_width*data_height*data_depth];



	// calculate grid spacings for gradient calculation
	float grid_x = (_params.max_bound.x - _params.min_bound.x)/data_width;
	float grid_y = (_params.max_bound.y - _params.min_bound.y)/data_height;
	float grid_z = (_params.max_bound.z - _params.min_bound.z)/data_depth;

	printf("grid_x: %f, grid_y: %f, grid_z: %f\n", grid_x, grid_y, grid_z);

	// loop over all the grid points and compute the refractive index gradient
	for(size_t z = 0; z < data_depth; z++)
	{
		for(size_t y = 0; y < data_height; y++)
		{
		  for(size_t x = 0; x < data_width; x++) {
			size_t DELTA = 1;
			float3 lookup = {x,y,z};
			if (lookup.x < DELTA || lookup.y < DELTA || lookup.z < DELTA ||
				lookup.x >= data_width-DELTA || lookup.y >= data_height -DELTA || lookup.z >=data_depth-DELTA)
			  continue;

			// obtain refractive index values from points lying on either side of the current grid
			// point along x, y and z

			float3 sample1, sample2;
			lookup = make_float3(x-1,y,z);
			sample1.x = data[size_t(lookup.z*data_width*data_height + lookup.y*data_width + lookup.x)];
			lookup = make_float3(x+1,y,z);
			sample2.x = data[size_t(lookup.z*data_width*data_height + lookup.y*data_width + lookup.x)];

			lookup = make_float3(x,y-1,z);
			sample1.y = data[size_t(lookup.z*data_width*data_height + lookup.y*data_width + lookup.x)];
			lookup = make_float3(x,y+1,z);
			sample2.y = data[size_t(lookup.z*data_width*data_height + lookup.y*data_width + lookup.x)];

			lookup = make_float3(x,y,z-1);
			sample1.z = data[size_t(lookup.z*data_width*data_height + lookup.y*data_width + lookup.x)];
			lookup = make_float3(x,y,z+1);
			sample2.z = data[size_t(lookup.z*data_width*data_height + lookup.y*data_width + lookup.x)];

			// calculate the refractive index gradient
			float3 normal;
			normal.x = (sample2.x - sample1.x)/(2*grid_x);
			normal.y = (sample2.y - sample1.y)/(2*grid_y);
			normal.z = (sample2.z - sample1.z)/(2*grid_z);

			// this is the array index where the refractive index gradient value will be stored
			int data_loc = (z*data_width*data_height + y*data_width + x);

			// assign refractive index gradient values to the array
			_params.data[data_loc].x = normal.x;
			_params.data[data_loc].y = normal.y;
			_params.data[data_loc].z = normal.z;
			_params.data[data_loc].w = data[size_t(z*data_width*data_height + y*data_width + x)];

			if(_params.data[data_loc].w == 0)
				printf("_params.data[data_loc].w == 0 at data_loc: %d\n", data_loc);

			if (_params.data[(z*data_width*data_height + y*data_width + x)].w < _params.data_min)
			  _params.data_min = _params.data[(z*data_width*data_height + y*data_width + x)].w;

		  }
		}
	}

	return _params;
}

density_grad_params_t readDatafromFile(char* filename)
{
    /*
     * This function reads density data from the file and sets up the simulation parameters
     * for the density gradients.
     *
     * INPUT:
     * filename - name of the nrrd file containing the density data
     * 			  (NOTE: Data should be in single precision!)
     *
     * OUTPUT:
     * _params - structure containing the parameters required to simulate ray deflection in
     * 			a variable density environment
     *
     */

	    // this is the structure that will hold the simulation parameters
    density_grad_params_t _params;

    // file object to read data
    vector<string> files;

    // data structure that will contain the density and grid information from the file
    DataFile* dataFiles[200];

    // add the current filename to the list
    files.push_back(string(filename));

    // create a datafile element
    for(int file = 0; file < files.size(); file++)
    {
        dataFiles[file] = new DataFile();
        dataFiles[file]->filename = new char[256];
    }

    char** input_files;
    input_files = new char*[files.size()];

    // these are the minimum and maximum number of data points to read from the file
    // (not sure what these are for)
	int data_min = 0;
	int data_max = 1024;

    for( int i = 0; i < files.size(); i++)
    {
        cout << "file: " << files[i] << endl;
        input_files[i] = new char[files[i].length()];
        strcpy(input_files[i], files[i].c_str());
        strcpy(dataFiles[i]->filename, input_files[i]);
        loadNRRD(dataFiles[i],data_min,data_max);
    }

    // this is the address to the location where the density data are stored
    float* data = dataFiles[0]->data;

    // get minimum and maximum coordinates along all three axes
    _params.min_bound = make_float3(dataFiles[0]->xmin, dataFiles[0]->ymin, dataFiles[0]->zmin);
    _params.max_bound = make_float3(dataFiles[0]->xmax, dataFiles[0]->ymax, dataFiles[0]->zmax);

    // get number of grid points along all three axes where density data are available
    _params.data_width = dataFiles[0]->sizex;
    _params.data_height = dataFiles[0]->sizey;
    _params.data_depth = dataFiles[0]->sizez;

    cout << "setting up renderer\n";
    _params = setData(data, _params);

    // set step size to be minimum of the three grid spacings
    float step_size = fmin(dataFiles[0]->del_x,dataFiles[0]->del_y);
    step_size = step_size < dataFiles[0]->del_z ? step_size : dataFiles[0]->del_z;
    printf("step_size: %f\n", step_size);

    _params.step_size = dataFiles[0]->del_z;

    return _params;

}

__global__ void check_trace_rays_through_density_gradients(light_ray_data_t* light_ray_data, density_grad_params_t params)
{
	/*
	 * This function is used to check the ray deflection produced by the
	 * trace_rays_through_density_gradients routine for a set of rays that are parallel
	 * to the z axis
	 */

	int id = threadIdx.x;

	light_ray_data[id] = trace_rays_through_density_gradients(light_ray_data[id],params);

}

#endif /* TRACE_RAYS_THROUGH_DENSITY_GRADIENTS_H_ */
