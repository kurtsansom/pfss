#ifndef CUDA_INTERFACE_H
#define CUDA_INTERFACE_H

#include <cuda.h>
#include <cuda_runtime.h>
#include <string>
#include <iostream>

using namespace std;


// this assumes CUDA compute capability >= 1.2
#ifdef CUDACC
	#if CUDACC==35 or CUDACC==30
		#define MAX_TPB 1024
		#define MAX_BPSM 16
		#define NUM_BLOCKS_WITH_MAX_THREADS 2
	#elif CUDACC==21 or CUDACC==20
		#define MAX_TPB 1024
		#define MAX_BPSM 8
		#define NUM_BLOCKS_WITH_MAX_THREADS 1
	#else
		#define MAX_TPB 512
		#define MAX_BPSM 8
		#define NUM_BLOCKS_WITH_MAX_THREADS 2
	#endif
#else
	#define MAX_TPB 256
	#define MAX_BPSM 1
	#define NUM_BLOCKS_WITH_MAX_THREADS 1
#endif


#define MAX_NUM_BLOCK 10000;

void checkCUDAerror(const char *msg);

/*! \brief checks for any errors related to CUDA and exits, if there are any (you need to fix this or stop using CUDA)	*/
void checkCUDAerror(const string &file, const int &line);

#ifdef __NVCC__
__host__ __device__
#endif
bool isInside_elliptic(uint index, double a2, double b2, double c2, double *d_r, double *d_theta, double *d_phi, \
                                  double *d_x, double *d_y, double *d_z);

#endif // CUDA_INTERFACE_H
