/*
 * trace_rays_through_density_gradients.h
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

cudaArray* data_array = 0, *texture_array = 0;

using namespace std;
texture<float4, 3> tex_data;

typedef struct {
  float* data;
  char* filename;
  int sizex, sizey, sizez;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double del_x, del_y, del_z;

} DataFile;



__device__ bool IntersectWithVolume(float3* ray_pos, float3 ray_dir, float3 p1, float3 p2)
{
//	float3 ray_pos = *ray_pos_p;

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

//	ray_pos = ray_pos+ray_dir*t;
	ray_pos->x += ray_dir.x*t;
	ray_pos->y += ray_dir.y*t;
	ray_pos->z += ray_dir.z*t;


	return true;
}


__device__ light_ray_data_t trace_rays_through_density_gradients(light_ray_data_t light_ray_data, density_grad_params_t params)
{

 	float3 min_bound = params.min_bound, max_bound = params.max_bound;
 	float3 lookup_scale = {1.0f/(max_bound.x-min_bound.x), 1.0f/(max_bound.y - min_bound.y), 1.0f/(max_bound.z-min_bound.z)};
 	int data_width = params.data_width, data_height = params.data_height, data_depth = params.data_depth;

 	float max_scale = max(max(float(params.data_width), float(params.data_height)), float(params.data_depth));

  	// calculate grid spacings for gradient calculation
  	float grid_x = (max_bound.x - min_bound.x)/data_width;
  	float grid_y = (max_bound.y - min_bound.y)/data_height;
  	float grid_z = (max_bound.z - min_bound.z)/data_depth;

  	//printf("grid_x: %f, grid_y: %f, grid_z: %f\n", grid_x, grid_y, grid_z);

 	float3 pos = light_ray_data.ray_source_coordinates;
 	float3 dir = light_ray_data.ray_propagation_direction;

 	// if initial position is outside the volume, intersect ray with volume
 	if(pos.x <= min_bound.x || pos.y <= min_bound.y || pos.z <= min_bound.z ||
 	pos.x >= max_bound.x || pos.y >= max_bound.y || pos.z >= max_bound.z )
 	{
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

 	//printf("Intersected Volume\n");

 	float3 normal;

 	int i = 0;
 	int insideBox = 1;
 	float3 lookupfn;
 	// Trace Ray through volume
 	while(insideBox==1)
 	{
 		i = i+1;

 		float3 offset = pos-min_bound;
// 		float3 lookupfn = lookup_scale*offset; // normalized lookup
 		lookupfn.x = lookup_scale.x*offset.x;
 		lookupfn.y = lookup_scale.y*offset.y;
 		lookupfn.z = lookup_scale.z*offset.z;


 		float3 lookup = {static_cast<float>(lookupfn.x*params.data_width), static_cast<float>(lookupfn.y*params.data_height), static_cast<float>(lookupfn.z*params.data_depth)};

 		if(pos.x < min_bound.x || pos.y < min_bound.y || pos.z < min_bound.z ||
 		pos.x > max_bound.x || pos.y > max_bound.y || pos.z > max_bound.z )
 		{
 		 //printf("pos: %f, %f, %f\n", pos.x,pos.y,pos.z);
 		 break;
 		}

 		float4 val = tex3D(tex_data, round(lookup.x), round(lookup.y), round(lookup.z)); //*params.dataScalar;

 		normal = make_float3(val.x/(2*grid_x),val.y/(2*grid_y),val.z/(2*grid_z));
 		//normal = make_float3(val.x/grid_x,val.y/grid_y,val.z/grid_z);
 		if(normal.x!=0 || normal.y!=0 || normal.z!=0)
// 			printf("normal: %f, %f, %f\n",normal.x,normal.y,normal.z);
 		//#if !LINE_OF_SIGHT
 		dir = dir + params.step_size*normal;
 		dir = normalize(dir);
	    pos = pos + dir*params.step_size; ///old_index;

       }

 	light_ray_data.ray_source_coordinates = pos;
 	light_ray_data.ray_propagation_direction = dir;

 	return light_ray_data;
}

extern "C"{

void Host_Init(density_grad_params_t* paramsp, density_grad_params_t* dparams)
{
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
	/* This function reads the density data from the NRRD file in addition to
	variables containing information about the extent of the volume and the grid
	spacing. It also converts the density to refractive index.
	Information about the nrrd file format is available at :
	http://teem.sourceforge.net/nrrd/lib.html
	*/

	printf("loading file %s : ", datafile->filename);
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
	xmax = xmin + sizex*del_x;

	ymin = nrrd->spaceOrigin[1];
	del_y = nrrd->axis[1].spacing;
	ymax = ymin + sizey*del_y;

	zmin = nrrd->spaceOrigin[2];
	del_z = nrrd->axis[2].spacing;
	zmax = zmin + sizez*del_z;

	printf("\n******************** Co-ordinate System Info ******************************\n");
	printf("xmin: %f, xmax: %f, N_x: %d, del_x: %f\n", xmin,xmax,sizex,del_x);
	printf("ymin: %f, ymax: %f, N_y: %d, del_y: %f\n", ymin,ymax,sizey,del_y);
	printf("zmin: %f, zmax: %f, N_z: %d, del_z: %f\n", zmin,zmax,sizez,del_z);


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


density_grad_params_t setData(float* data, int data_width, int data_height, int data_depth, density_grad_params_t _params)
{

	printf("setData(%d, %d, %d, %d)\n", (u_int64_t)data, data_width, data_height, data_depth);
	_params.data_min = FLT_MAX;

	//compute gradient
	_params.data = new float4[data_width*data_height*data_depth];
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
			float3 normal;
			normal.x = sample1.x - sample2.x;
			normal.y = sample1.y - sample2.y;
			normal.z = sample1.z - sample2.z;

			float4& datap = _params.data[size_t(z*data_width*data_height + y*data_width + x)];
			datap.x = -normal.x;
			datap.y = -normal.y;
			datap.z = -normal.z;
			datap.w = data[size_t(z*data_width*data_height + y*data_width + x)];

			/*
			if(normal.x !=0 || normal.y !=0 || normal.z!=0)
				printf("lookup: %f,%f,%f normal: %f, %f, %f \n", lookup.x, lookup.y, lookup.z, normal.x, normal.y, normal.z);
			*/

			if (datap.w < _params.data_min)
			  _params.data_min = datap.w;
		  }
		}
	}
	float m = fmax(float(data_width), float(data_height));
	m = fmax(m, float(data_depth));
	float3 dataScale = make_float3(float(data_width)/m, float(data_height)/m, float(data_depth)/m);
	//_params.min_bound = -dataScale/2.0f;
	//_params.max_bound = dataScale/2.0f;

	_params.data_width = data_width;
	_params.data_height = data_height;
	_params.data_depth = data_depth;

	return _params;
}


density_grad_params_t readDatafromFile(char* filename)
{
    float step_size = 0.1;
    int data_min = 0;
    int data_max = 1024;

    density_grad_params_t _params;

    vector<string> files; //, tempFiles, pressFiles;
    DataFile* dataFiles[200];

    files.push_back(string(filename));

    for(int file = 0; file < files.size(); file++)
    {
        dataFiles[file] = new DataFile();
        dataFiles[file]->filename = new char[256];
    }

    char** input_files;
    input_files = new char*[files.size()];


    for( int i = 0; i < files.size(); i++)
    {
        cout << "file: " << files[i] << endl;
        input_files[i] = new char[files[i].length()];
        strcpy(input_files[i], files[i].c_str());
        strcpy(dataFiles[i]->filename, input_files[i]);
        loadNRRD(dataFiles[i],data_min,data_max);
    }

    float* data = dataFiles[0]->data;
    int zrange = dataFiles[0]->sizez;


    // get min and max bound along all three axes
    _params.min_bound = make_float3(dataFiles[0]->xmin, dataFiles[0]->ymin, dataFiles[0]->zmin);
    _params.max_bound = make_float3(dataFiles[0]->xmax, dataFiles[0]->ymax, dataFiles[0]->zmax);

    _params.data_width = dataFiles[0]->sizex;
    _params.data_height = dataFiles[0]->sizey;
    _params.data_depth = dataFiles[0]->sizez;

    cout << "setting up renderer\n";

    _params = setData(data, dataFiles[0]->sizex,dataFiles[0]->sizey,zrange,_params);

    // set step size to be minimum of the three grid spacings
    step_size = fmin(dataFiles[0]->del_x,dataFiles[0]->del_y);
    step_size = step_size < dataFiles[0]->del_z ? step_size : dataFiles[0]->del_z;
    printf("step_size: %f\n", step_size);
    _params.step_size = step_size;

    return _params;

}

__global__ void check_trace_rays_through_density_gradients(light_ray_data_t* light_ray_data, density_grad_params_t params)
{
	int id = threadIdx.x;

 	float3 min_bound = params.min_bound, max_bound = params.max_bound;
 	float3 lookup_scale = {1.0f/(max_bound.x-min_bound.x), 1.0f/(max_bound.y - min_bound.y), 1.0f/(max_bound.z-min_bound.z)};
 	int data_width = params.data_width, data_height = params.data_height, data_depth = params.data_depth;

 	float max_scale = max(max(float(params.data_width), float(params.data_height)), float(params.data_depth));

  	// calculate grid spacings for gradient calculation
  	float grid_x = (max_bound.x - min_bound.x)/data_width;
  	float grid_y = (max_bound.y - min_bound.y)/data_height;
  	float grid_z = (max_bound.z - min_bound.z)/data_depth;

  	//printf("grid_x: %f, grid_y: %f, grid_z: %f\n", grid_x, grid_y, grid_z);

 	float3 pos = light_ray_data[id].ray_source_coordinates;
 	float3 dir = light_ray_data[id].ray_propagation_direction;

 	// if initial position is outside the volume, intersect ray with volume
 	if(pos.x <= min_bound.x || pos.y <= min_bound.y || pos.z <= min_bound.z ||
 	pos.x >= max_bound.x || pos.y >= max_bound.y || pos.z >= max_bound.z )
 	{
 		if(!IntersectWithVolume(&pos, dir, params.min_bound, params.max_bound))
 		{
 			//# % This sets any of the light ray positions outside of the domain
 			//# % to NaN values
 			light_ray_data[id].ray_source_coordinates = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);

 			//# % This sets any of the light ray directions outside of the domain
 			//# % to NaN values
 			light_ray_data[id].ray_propagation_direction = make_float3(CUDART_NAN_F,CUDART_NAN_F,CUDART_NAN_F);


 		}
 	}

 	//printf("Intersected Volume\n");

 	float3 normal;

 	int i = 0;
 	int insideBox = 1;
 	float3 lookupfn;
 	// Trace Ray through volume
 	while(insideBox==1)
 	{
 		i = i+1;

 		float3 offset = pos-min_bound;
// 		float3 lookupfn = lookup_scale*offset; // normalized lookup
 		lookupfn.x = lookup_scale.x*offset.x;
 		lookupfn.y = lookup_scale.y*offset.y;
 		lookupfn.z = lookup_scale.z*offset.z;


 		float3 lookup = {static_cast<float>(lookupfn.x*params.data_width), static_cast<float>(lookupfn.y*params.data_height), static_cast<float>(lookupfn.z*params.data_depth)};

 		if(pos.x < min_bound.x || pos.y < min_bound.y || pos.z < min_bound.z ||
 		pos.x > max_bound.x || pos.y > max_bound.y || pos.z > max_bound.z )
 		{
 		 //printf("pos: %f, %f, %f\n", pos.x,pos.y,pos.z);
 		 break;
 		}

 		float4 val = tex3D(tex_data, round(lookup.x), round(lookup.y), round(lookup.z)); //*params.dataScalar;

 		normal = make_float3(val.x/(2*grid_x),val.y/(2*grid_y),val.z/(2*grid_z));
 		//normal = make_float3(val.x/grid_x,val.y/grid_y,val.z/grid_z);
 		if(normal.x!=0 || normal.y!=0 || normal.z!=0)
// 			printf("normal: %f, %f, %f\n",normal.x,normal.y,normal.z);
 		//#if !LINE_OF_SIGHT
 		dir = dir + params.step_size*normal;
 		dir = normalize(dir);
	    pos = pos + dir*params.step_size; ///old_index;

       }

 	light_ray_data[id].ray_source_coordinates = pos;
 	light_ray_data[id].ray_propagation_direction = (dir);

}

#endif /* TRACE_RAYS_THROUGH_DENSITY_GRADIENTS_H_ */
