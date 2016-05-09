/*
 * float3_operators.h
 *
 *	This file contains functions that define arithmetic operations on float3 and float4
 *	variables using operator overloading.
 *
 *  Created on: Apr 22, 2016
 *      Author: lrajendr
 */

#ifndef FLOAT3_OPERATORS_H_
#define FLOAT3_OPERATORS_H_

#include <string>

__device__ float3 operator/(const float3 &a, const float &b) {

return make_float3(a.x/b, a.y/b, a.z/b);

}

__device__ float3 operator*(const float3 &a, const float &b) {

return make_float3(a.x*b, a.y*b, a.z*b);

}

__device__ float3 operator*(const float &b,const float3 &a) {

return make_float3(a.x*b, a.y*b, a.z*b);

}

__device__ float3 operator+(const float3 &a, const float3 &b) {

return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);

}

__device__ float3 operator-(const float3 &a, const float3 &b) {

return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);

}

__device__ float3 operator-(const float3 &a) {

return make_float3(-a.x, -a.y, -a.z);

}

__device__ bool operator==(const float3&a, const float3&b) {
	if(a.x==b.x && a.y==b.y && a.z==b.z)
		return true;

	return false;
}

__device__ bool operator==(const float4&a, const float4&b) {
	if(a.x==b.x && a.y==b.y && a.z==b.z && a.w==b.w)
		return true;

	return false;
}


__device__ float3 normalize(const float3&a) {
    float m = sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
    return a/m;
}

__device__ float3 cross(const float3&a, const float3&b) {
    return make_float3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}


__device__ float dot(const float3&a, const float3&b) {
    return (a.x*b.x + a.y*b.y + a.z*b.z);
}

__device__ float angleBetween(const float3& v1, const float3& v2)
{
  // returns angle between two vectors in degrees
  float dp = dot(v1, v2);
  return acosf(dp)*180.0/M_PI;
}

void print3(char* msg,float3 v)
{
//#if __DEVICE_EMULATION__
  printf("%s",msg);
  printf(": %f %f %f\n", v.x, v.y, v.z);
//#endif
}
#endif /* FLOAT3_OPERATORS_H_ */
