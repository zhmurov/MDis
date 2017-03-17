#pragma once

#include <math.h>
#include <vector_types.h>
#include "Cuda.h"
#include "LinearAlgebra.h"

namespace quaternion {
	CUDA_FUNC inline float4 sum(float4 a, float4 b) { return a + b; }
	CUDA_FUNC inline float4 mul(float4 a, float4 b) {
		return make_float4(
				a.y*b.z - a.z*b.y + a.w*b.x + b.w*a.x,
				- a.x*b.z + a.z*b.x + a.w*b.y + b.w*a.y,
				a.x*b.y - a.y*b.x + a.w*b.z + b.w*a.z,
				a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
				);
	}
	CUDA_FUNC inline float4 mul_vect(float4 a, float4 b) {
		return make_float4(
				a.y*b.z - a.z*b.y,
				- a.x*b.z + a.z*b.x,
				a.x*b.y - a.y*b.x, 0.0f
				);
	}
	CUDA_FUNC inline float4 conj(float4 a) {
		return make_float4(-a.x, -a.y, -a.z, a.w);
	}
	CUDA_FUNC inline float magn(float4 a) {
		return sqrtf(a.x*a.x + a.y*a.y + a.z*a.z + a.w*a.w);
	}
	CUDA_FUNC inline float norm(float4 a) {
		return a.x*a.x + a.y*a.y + a.z*a.z + a.w*a.w;
	}
	CUDA_FUNC inline float4 inverse(float4 a) {
		return conj(a) * (1.0f / norm(a));
	}
	// abs("a") = angle of rotation
	// direction "a" = axis of rotation
	CUDA_FUNC inline float4 makeRotationAround(float4 a) {
		float angle = sqrtf(a.x*a.x + a.y*a.y + a.z*a.z);
		if (angle < 1e-5f)
			return make_float4(0, 0, 0, 1);
		float sin_a = sinf(angle*0.5f);
		return make_float4(sin_a*a.x/angle, sin_a*a.y/angle, sin_a*a.z/angle, (float)cosf(angle*0.5f));
	}
	CUDA_FUNC inline float4 toRotationVector(float4 q) {
		float4 ret = q * (2.0f * acosf(q.w) / sqrtf(1 - q.w * q.w));
		ret.w = 0.0f;
		return ret;
	}
	CUDA_FUNC inline void toMatrix(float4 q, Matrix<float, 3, 3>& m) {
		m(0,0) = 1.0f - (q.y * q.y * 2 + q.z * q.z * 2);
		m(1,0) = q.x * q.y * 2 - q.w * q.z * 2;
		m(2,0) = q.x * q.z * 2 + q.w * q.y * 2;
		m(0,1) = q.x * q.y * 2 + q.w * q.z * 2;
		m(1,1) = 1.0f - (q.x * q.x * 2 + q.z * q.z * 2);
		m(2,1) = q.y * q.z * 2 - q.w * q.x * 2;
		m(0,2) = q.x * q.z * 2 - q.w * q.y * 2;
		m(1,2) = q.y * q.z * 2 + q.w * q.x * 2;
		m(2,2) = 1.0f - (q.x * q.x * 2 + q.y * q.y * 2);
	}
}

