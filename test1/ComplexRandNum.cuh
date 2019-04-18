#ifndef _COMPLEXRANDNUM_CUH_
#define _COMPLEXRANDNUM_CUH_

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <curand.h>
#include <curand_kernel.h>
#include <time.h>
#include <math.h>
#include <CUDACOMPLEX.cuh>

//both have to malloc for the result
void getRandomNumber(CudaComplex *result, int num, int state_num);				//get random number in host
__device__ void getRand(CudaComplex *N, curandState *state, int num_rand);		//get random number in device
__device__ CudaComplex rand_xorshift(unsigned int* rng_state);
__device__ CudaComplex rand_linear(unsigned int *rng_state);
#endif