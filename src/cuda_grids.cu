#include <cuda.h>
#include <cuda_runtime.h>

#include "src/grids.h"
#include "src/ellipticalGrid.h"
#include "src/cuda_interface.h"


__global__ void
__launch_bounds__(512, 1)
kernel_print(SphericalGrid *d_grid)
{
	printf("lowerR: %E\n", d_grid->lowerR);
}


void SphericalGrid::initNULL_GPU()
{
    d_memStruct	= NULL;
    onDevice    = false;

    checkCUDAerror("SphericalGrid::initNULL_GPU");
}

void SphericalGrid::clear_GPU()
{
    if(d_memStruct != NULL)
    	cudaFree(d_memStruct);

    checkCUDAerror("SphericalGrid::clear_GPU");

    onDevice = false;
    initNULL_GPU();
}

/*! this is only responsible for initializing the pointers of the object in GPU memory
 *
 */
void SphericalGrid::initGPUmemStruct()
{
	long unsigned index 			= sizeof(SphericalGrid);
	long unsigned newCudaMemAddr 	= (long unsigned)d_memStruct + index;

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->pos), &newCudaMemAddr, sizeof(Vec3D *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Vec3D);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->B), &newCudaMemAddr, sizeof(Vec3D *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Vec3D);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->psi), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->relError), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->temp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_imjk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ipjk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijmk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijpk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijkm), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijkp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_imjmk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_imjpk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ipjmk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ipjpk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_imjkm), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_imjkp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ipjkm), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ipjkp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijmkm), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijmkp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijpkm), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijpkp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->g), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->h), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->p), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->p_r), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->p_t), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->p_p), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	checkCUDAerror("SphericalGrid::initGPUmemStruct");
}

/*! Allocates GPU memory for object. Amount of memory is (runtime) constant. All pointers to arrays / objects are
 *  initialized to point to their reserved memory addresses
 */
bool SphericalGrid::allocateGPUmem()
{
	clear_GPU();

	//				= Memory Structure		+ 					pos / B 			+ psi / relError / temp + g, h, p, p_r, p_t, p_p
	//																				+ 19 scale factors
	size_t memsize 	= sizeof(SphericalGrid) + numGridPoints * ( 2 * sizeof(Vec3D) 	+ 22 * sizeof(hcFloat) 	+ 6 * sizeof(Matrix3x3));
	SphericalGrid empty;
	cudaMalloc((void **) &d_memStruct, memsize);
	cudaMemcpy((char *)d_memStruct, &empty, sizeof(SphericalGrid), cudaMemcpyHostToDevice);
	//initGPUmemStruct();						// TODO: unnecessary due to pushToGPU?

	checkCUDAerror("SphericalGrid::allocateGPUmem");

	return true;
}

bool SphericalGrid::pushToGPU(){

	if(numGridPoints == 0)
	{
		printf("ERROR! SphericalGrid::pushToDevice: numGridPoints=0 -> nothing to do!\n");
		return false;
	}

    allocateGPUmem();

	long unsigned index 			= 0;
	long unsigned newCudaMemAddr 	= (long unsigned)d_memStruct + index;

	cudaMemcpy((char *)newCudaMemAddr, this, sizeof(SphericalGrid), cudaMemcpyHostToDevice);
	initGPUmemStruct();

	index = sizeof(SphericalGrid);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, pos, numGridPoints * sizeof(Vec3D), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Vec3D);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, B, numGridPoints * sizeof(Vec3D), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Vec3D);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, psi, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, relError, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, temp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_imjk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ipjk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijmk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijpk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijkm, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijkp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_imjmk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_imjpk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ipjmk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ipjpk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_imjkm, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_imjkp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ipjkm, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ipjkp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijmkm, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijmkp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijpkm, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijpkp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, g, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, h, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, p, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, p_r, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, p_t, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, p_p, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	onDevice = true;

	checkCUDAerror("SphericalGrid::pushToGPU");

	return true;
}

bool SphericalGrid::pullFromGPU()
{
	if(!onDevice)
		return false;

	clear_CPU();

	SphericalGrid *tAlloc = new SphericalGrid;
	cudaMemcpy(tAlloc, d_memStruct, sizeof(SphericalGrid), cudaMemcpyDeviceToHost);

	init(tAlloc->sinLatGrid, tAlloc->maxSinLat, tAlloc->minSinLat, tAlloc->lowerR, tAlloc->upperR, tAlloc->numR, false, 1.0);

	cudaMemcpy(this->pos, 		tAlloc->pos, 		numGridPoints * sizeof(Vec3D), 	cudaMemcpyDeviceToHost);
	cudaMemcpy(this->B, 		tAlloc->B, 			numGridPoints * sizeof(Vec3D), 	cudaMemcpyDeviceToHost);
	cudaMemcpy(this->psi, 		tAlloc->psi, 		numGridPoints * sizeof(hcFloat), cudaMemcpyDeviceToHost);
	cudaMemcpy(this->relError, 	tAlloc->relError, 	numGridPoints * sizeof(hcFloat), cudaMemcpyDeviceToHost);
	cudaMemcpy(this->temp, 		tAlloc->temp, 		numGridPoints * sizeof(hcFloat), cudaMemcpyDeviceToHost);
	// the following should not be necessary to fetch from GPU
	/*
	cudaMemcpy(this->g, 		tAlloc->g,	 		numGridPoints * sizeof(Matrix3x3), cudaMemcpyDeviceToHost);
	cudaMemcpy(this->h, 		tAlloc->h,	 		numGridPoints * sizeof(Matrix3x3), cudaMemcpyDeviceToHost);
	cudaMemcpy(this->p, 		tAlloc->p,	 		numGridPoints * sizeof(Matrix3x3), cudaMemcpyDeviceToHost);
	cudaMemcpy(this->p_r, 		tAlloc->p_r, 		numGridPoints * sizeof(Matrix3x3), cudaMemcpyDeviceToHost);
	cudaMemcpy(this->p_t, 		tAlloc->p_t, 		numGridPoints * sizeof(Matrix3x3), cudaMemcpyDeviceToHost);
	cudaMemcpy(this->p_p, 		tAlloc->p_p, 		numGridPoints * sizeof(Matrix3x3), cudaMemcpyDeviceToHost);//*/

	tAlloc->initNULL();
	delete tAlloc;

	checkCUDAerror("CudaAllocator::pullFromGPU");

	return true;
}

void SphericalGrid::extract_relError(){

	SphericalGrid *tAlloc = new SphericalGrid;;
	cudaMemcpy(tAlloc, d_memStruct, sizeof(SphericalGrid), cudaMemcpyDeviceToHost);

	cudaMemcpy(relError, tAlloc->relError, numGridPoints * sizeof(hcFloat), cudaMemcpyDeviceToHost);

	tAlloc->initNULL();
	delete tAlloc;

	checkCUDAerror("SphericalGrid::extractRelError");
}


/*! this is only responsible for initializing the pointers of the object in GPU memory
 *	TODO: so much code duplication
 */
void EllipticalGrid::initGPUmemStruct(){

	long unsigned index 			= sizeof(EllipticalGrid);
	long unsigned newCudaMemAddr 	= (long unsigned)d_memStruct + index;

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->pos), &newCudaMemAddr, sizeof(Vec3D *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Vec3D); // TODO Vec3D* ?

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->B), &newCudaMemAddr, sizeof(Vec3D *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Vec3D);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->psi), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->relError), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->temp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_imjk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ipjk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijmk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijpk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijkm), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijkp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_imjmk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_imjpk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ipjmk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ipjpk), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_imjkm), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_imjkp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ipjkm), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ipjkp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijmkm), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijmkp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijpkm), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(d_memStruct->s_ijpkp), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->g), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->h), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->p), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->p_r), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->p_t), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->p_p), &newCudaMemAddr, sizeof(Matrix3x3 *), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy(&(((EllipticalGrid*)d_memStruct)->a), &newCudaMemAddr, sizeof(hcFloat *), cudaMemcpyHostToDevice);
	index += numR * sizeof(hcFloat);

	checkCUDAerror("EllipticalGrid::initGPUmemStruct");
}

/*! Allocates GPU memory for object. Amount of memory is (runtime) constant. All pointers to arrays / objects are
 *  initialized to point to their reserved memory addresses
 *
 */
bool EllipticalGrid::allocateGPUmem(){

	clear_GPU();

	//				= Memory Structure		 + 					pos / B 			+ psi / relError / temp	+ g,h,p,p_r,p_t,p_p			+ a
	//																				+ 19 scale factors
	size_t memsize 	= sizeof(EllipticalGrid) + numGridPoints * ( 2 * sizeof(Vec3D) 	+ 22 * sizeof(hcFloat) 	+ 6 * sizeof(Matrix3x3))	+ numR * sizeof(hcFloat);
	EllipticalGrid empty;
	cudaMalloc((void **) &d_memStruct, memsize);
	cudaMemcpy((char *)d_memStruct, &empty, sizeof(EllipticalGrid), cudaMemcpyHostToDevice);

	checkCUDAerror("EllipticalGrid::allocateGPUmem");

	return true;
}

bool EllipticalGrid::pushToGPU(){

	if(numGridPoints == 0)
	{
		printf("ERROR! EllipticalGrid::pushToDevice: numGridPoints=0 -> nothing to do!\n");
		return false;
	}

    allocateGPUmem();

	long unsigned index 			= 0;
	long unsigned newCudaMemAddr 	= (long unsigned)d_memStruct + index;

	cudaMemcpy((char *)newCudaMemAddr, this, sizeof(EllipticalGrid), cudaMemcpyHostToDevice);
	initGPUmemStruct();

	index = sizeof(EllipticalGrid);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, pos, numGridPoints * sizeof(Vec3D), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Vec3D);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, B, numGridPoints * sizeof(Vec3D), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Vec3D);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, psi, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, relError, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, temp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_imjk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ipjk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijmk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijpk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijkm, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijkp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_imjmk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_imjpk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ipjmk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ipjpk, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_imjkm, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_imjkp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ipjkm, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ipjkp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijmkm, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijmkp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijpkm, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, s_ijpkp, numGridPoints * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(hcFloat);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, g, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, h, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, p, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, p_r, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, p_t, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, p_p, numGridPoints * sizeof(Matrix3x3), cudaMemcpyHostToDevice);
	index += numGridPoints * sizeof(Matrix3x3);

	newCudaMemAddr = (long unsigned)d_memStruct + index;
	cudaMemcpy((char *)newCudaMemAddr, a, numR * sizeof(hcFloat), cudaMemcpyHostToDevice);
	index += numR * sizeof(hcFloat);

	onDevice = true;

	checkCUDAerror("EllipticalGrid::pushToGPU");

	return true;
}

bool EllipticalGrid::pullFromGPU()
{
	if(!onDevice)
		return false;

	clear_CPU();

	EllipticalGrid *tAlloc = new EllipticalGrid;
	cudaMemcpy(tAlloc, d_memStruct, sizeof(EllipticalGrid), cudaMemcpyDeviceToHost);

	init(tAlloc->sinLatGrid, tAlloc->maxSinLat, tAlloc->minSinLat, tAlloc->lowerR, tAlloc->upperR,	tAlloc->numR, false, 1.0);

	cudaMemcpy(this->pos, 		tAlloc->pos, 		numGridPoints * sizeof(Vec3D), 		cudaMemcpyDeviceToHost);
	cudaMemcpy(this->B, 		tAlloc->B, 			numGridPoints * sizeof(Vec3D), 		cudaMemcpyDeviceToHost);
	cudaMemcpy(this->psi, 		tAlloc->psi, 		numGridPoints * sizeof(hcFloat), 	cudaMemcpyDeviceToHost);
	cudaMemcpy(this->relError, 	tAlloc->relError, 	numGridPoints * sizeof(hcFloat), 	cudaMemcpyDeviceToHost);
	cudaMemcpy(this->temp, 		tAlloc->temp, 		numGridPoints * sizeof(hcFloat), 	cudaMemcpyDeviceToHost);
	// the following should not be necessary to pull from GPU
	/*
	cudaMemcpy(this->g, 		tAlloc->g, 			numGridPoints * sizeof(Matrix3x3), 	cudaMemcpyDeviceToHost);
	cudaMemcpy(this->h, 		tAlloc->h, 			numGridPoints * sizeof(Matrix3x3), 	cudaMemcpyDeviceToHost);
	cudaMemcpy(this->p, 		tAlloc->p, 			numGridPoints * sizeof(Matrix3x3), 	cudaMemcpyDeviceToHost);
	cudaMemcpy(this->p_r, 		tAlloc->p_r,		numGridPoints * sizeof(Matrix3x3), 	cudaMemcpyDeviceToHost);
	cudaMemcpy(this->p_t, 		tAlloc->p_t,		numGridPoints * sizeof(Matrix3x3), 	cudaMemcpyDeviceToHost);
	cudaMemcpy(this->p_p, 		tAlloc->p_p,		numGridPoints * sizeof(Matrix3x3), 	cudaMemcpyDeviceToHost);//*/
	cudaMemcpy(this->a, 		tAlloc->a, 			numR * sizeof(hcFloat), 			cudaMemcpyDeviceToHost);

	tAlloc->initNULL();
	delete tAlloc;
	checkCUDAerror("EllipticalGrid::pullFromGPU");

	return true;
}

void EllipticalGrid::extract_relError()
{
	EllipticalGrid *tAlloc = new EllipticalGrid;;
	cudaMemcpy(tAlloc, d_memStruct, sizeof(EllipticalGrid), cudaMemcpyDeviceToHost);

	cudaMemcpy(relError, tAlloc->relError, numGridPoints * sizeof(hcFloat), cudaMemcpyDeviceToHost);

	tAlloc->initNULL();
	delete tAlloc;

	checkCUDAerror("EllipticalGrid::extract_relError");
}
