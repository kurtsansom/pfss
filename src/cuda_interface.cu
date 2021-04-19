#include <cuda.h>
#include <cuda_runtime.h>

#include "src/cuda_interface.h"
#include "stdio.h"

#define MAX_NUM_BLOCK 10000;


void checkCUDAerror(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if(cudaSuccess != err)
    {
        fprintf(stderr, "CUDA error: %s: %s!\n", msg, cudaGetErrorString(err));
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
}

void checkCUDAerror(const string &file, const int &line)
{
    cudaError_t err = cudaGetLastError();
    if(cudaSuccess != err)
    {
    	cerr << file << line << ": CUDA error: " <<  cudaGetErrorString(err);
        exit(EXIT_FAILURE);
    }
}

bool isInside_elliptic(uint index, double a2, double b2, double c2, double *d_r, double *d_theta, double *d_phi, \
                                  double *d_x, double *d_y, double *d_z){

    if(d_x[index]*d_x[index] / a2 + d_y[index] * d_y[index] / b2 + d_z[index] * d_z[index] / c2 <= 1.0)
        return true;
    else
        return false;
}


