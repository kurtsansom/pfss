#include "src/cuda_interface.h"
#include "src/laplaceSolver.h"
#include "src/ellipticalGrid.h"

#include "engine/hcImage.h"

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          CUDA Kernels
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

/**********************************************************************************************************************
 *  compute lower boundary
 *
 *  this is correct for the general curvilinear case and the spheric case
 **********************************************************************************************************************/
__global__ void
__launch_bounds__(1024, 1)
kernel_computeLowerBoundary(LaplaceSolver &solver)
{
    uint idx 			= blockIdx.x * blockDim.x + threadIdx.x;
    SphericalGrid &gr	= *solver.grid;
    uint numR       	= gr.numR;
    uint numTheta   	= gr.numTheta;
    uint numPhi     	= gr.numPhi;

    // compute position in grid
    uint k      = idx % numPhi;
    uint j      = (idx - k) / numPhi;
    uint i      = 0;
    uint ind    = k * numR * numTheta + j * numR + i;
    uint texInd	= (numTheta - 1 - j) * numPhi + k;

    if(texInd < numTheta * numPhi)
    {
    	Vec3D pos		= gr.pos[ind];
		Vec3D pos_p		= gr.pos[ind+1];
		Vec3D pos_pp	= gr.pos[ind+2];
		hcFloat psi_p	= gr.psi[ind+1];
		hcFloat psi_pp	= gr.psi[ind+2];

        hcFloat B_l     = solver.d_image[texInd];
        hcFloat B_r     = B_l / __sinf(pos[1]);

        hcFloat dp		= pos_p[0]  - pos[0];
		hcFloat dpp		= pos_pp[0] - pos[0];

#ifdef RSCALE
        dp	/= r_sol;
		dpp	/= r_sol;
#endif

        hcFloat a		= 1 / (dpp*dpp - dp*dp);
        hcFloat psi   	= a * (B_r * (dp*dpp*dpp - dpp*dp*dp) + psi_p*dpp*dpp - psi_pp*dp*dp);

        gr.setPsi(ind, psi);
    }
}

/**********************************************************************************************************************
 *  iterate spheric grid
 *
 *  this version is specifically taylored to the spheric case. It is faster than the general purpose curvilinear
 *  version.
 **********************************************************************************************************************/

__global__ void
__launch_bounds__(1024, 1)
kernel_iterate_spheric(LaplaceSolver &solver, uint prevBlocks)
{
	const hcFloat PI  = 3.1415926535897;

	uint idx 		= (blockIdx.x + prevBlocks) * blockDim.x + threadIdx.x;

	SphericalGrid &gr = *solver.grid;

	uint numR       = gr.numR;
	uint numT	   	= gr.numTheta;
	uint numP   	= gr.numPhi;

	Vec3D *pos		= gr.getPosArray();
	hcFloat *psi	= gr.getPsiArray();

	// compute position in grid
	uint i = idx % numR;
	uint j = ((idx   - i) / numR) % numT;
	uint k = ((idx-i - j  * numR) / numR) / numT;

	uint ind = k * numR * numT + j * numR + i;

	if(idx < numT * numP * numR && i > 0)
	{
		if(i==numR-1)
		{
			gr.setTemp(ind, 0.0);
			gr.setRelError(ind, 0.0);
		}
		else
		{
			uint ind_r_m	= i==0				? gr.getIndex(i,j,k)					: gr.getIndex(i-1,j,k);
			uint ind_r_p    = i==numR-1			? gr.getIndex(i,j,k) 					: gr.getIndex(i+1,j,k);

			uint ind_t_m	= j==0				? 0										: gr.getIndex(i,j-1,k);
			uint ind_t_p	= j==numT-1			? 0										: gr.getIndex(i,j+1,k);

			uint ind_p_m  	= k==0		   		? gr.getIndex(i,j,numP-1)	 			: gr.getIndex(i,j,k-1);
			uint ind_p_p  	= k==numP-1 		? gr.getIndex(i,j,0) 					: gr.getIndex(i,j,k+1);


			hcFloat dr_p   	= pos[ind_r_p][0] - pos[ind][0];
			hcFloat dr_m   	= pos[ind][0]     - pos[ind_r_m][0];

			hcFloat dt_p	= (j==numT-1		? PI - pos[ind][1] 		  				: pos[ind_t_p][1]  	- pos[ind][1]);
			hcFloat dt_m	= (j==0				? pos[ind][1]			  				: pos[ind][1]       - pos[ind_t_m][1]);

			hcFloat dp_p	= (k==numP-1 		? pos[ind_p_p][2]  + 2*PI - pos[ind][2] : pos[ind_p_p][2]   - pos[ind][2]);
			hcFloat dp_m	= (k==0				? 2*PI - pos[ind_p_m][2]  + pos[ind][2] : pos[ind][2]		- pos[ind_p_m][2]);

			hcFloat r		= pos[ind][0];
			hcFloat theta	= pos[ind][1];
			hcFloat sin2	= __sinf(theta) * __sinf(theta);

#ifdef RSCALE
			dr_p    		/= r_sol;
			dr_m    		/= r_sol;
			r				/= r_sol;
#endif

			hcFloat psi_pole= 0.0;

			for(uint l=0; l<numP; ++l)
			{
				uint ind_pole 	= gr.getIndex(i,j,l);
				psi_pole 		+= psi[ind_pole];
			}
			psi_pole 			/= numP;

			hcFloat psi_ipjk	= psi[ind_r_p];
			hcFloat psi_imjk	= psi[ind_r_m];

			hcFloat psi_ijpk	= ( j==numT-1 		? psi_pole 		: psi[ind_t_p]);
			hcFloat psi_ijmk	= ( j==0			? psi_pole 		: psi[ind_t_m]);

			hcFloat psi_ijkp	= psi[ind_p_p];
			hcFloat psi_ijkm	= psi[ind_p_m];

			/*-----
			//----- first version, working but not as accurate as below
			//-----

			hcFloat h_r	= dr_p + dr_m;
			hcFloat h_t	= dt_p + dt_m;
			hcFloat h_p	= dp_p + dp_m;

			hcFloat A_r      = 1 / (dr_p * dr_m * h_r);
			hcFloat A_theta  = 1 / (dt_p * dt_m * h_t);
			hcFloat A_phi    = 1 / (dp_p * dp_m * h_p);

			hcFloat d 	= 2*A_r 	/ r 		* (dr_m*dr_m - dr_p*dr_p)
						+ 2*A_r 				* (dr_m + dr_p)
						+ A_theta 	/ (r*r) 	* cos(theta)/sin(theta) * (dt_m*dt_m - dt_p*dt_p)
						+ 2*A_theta / (r*r) 	* (dt_m + dt_p)
						+ 2*A_phi 	/ (r*r*sin2)* (dp_m + dp_p);

			hcFloat rhs 	= 2 * A_r 	/ r 		* (dr_m*dr_m*psi_ipjk - dr_p*dr_p*psi_imjk)
						+ 2*A_r 				* (dr_m * psi_ipjk + dr_p * psi_imjk)
						+ A_theta 	/ (r*r) 	* cos(theta) / sin(theta) * (dt_m*dt_m * psi_ijpk - dt_p*dt_p * psi_ijmk)
						+ 2*A_theta / (r*r) 	* (dt_m * psi_ijpk + dt_p * psi_ijmk)
						+ 2*A_phi 	/ (r*r*sin2)* (dp_m * psi_ijkp + dp_p * psi_ijkm);

			//-----
			////-- higher accuracy here (but numerically unstable due to two "interlocked" solutions)
			//-----

			hcFloat f_r		= - dr_m*dr_pp 	/ (dr_p*dr_p * (dr_pp + dr_p))
							+	dr_m*dr_m 	/ (dr_p*dr_p * h_r)
							-   2			/ h_r
							+	dr_p*dr_p	/ (dr_m*dr_m * h_r)
							-	dr_p*dr_mm	/ (dr_m*dr_m * (dr_mm + dr_m));

			hcFloat f_t		= - dt_m*dt_pp 	/ (dt_p*dt_p * (dt_pp + dt_p))
							+	dt_m*dt_m 	/ (dt_p*dt_p * h_t)
							-   2			/ h_t
							+	dt_p*dt_p	/ (dt_m*dt_m * h_t)
							-	dt_p*dt_mm	/ (dt_m*dt_m * (dt_mm + dt_m));

			hcFloat f_p		= - dp_m*dp_pp 	/ (dp_p*dp_p * (dp_pp + dp_p))
							+	dp_m*dp_m 	/ (dp_p*dp_p * h_p)
							-   2			/ h_p
							+	dp_p*dp_p	/ (dp_m*dp_m * h_p)
							-	dp_p*dp_mm	/ (dp_m*dp_m * (dp_mm + dp_m));


			hcFloat d		= 2 		/ (r*h_r)				* (dr_p/dr_m - dr_m/dr_p)
							+ cos(theta)/ (r*r*sin(theta)*h_t)	* (dt_p/dt_m - dt_m/dt_p)
							+ f_r		/  h_r
							+ f_t		/ (r*r*h_t)
							+ f_p		/ (r*r*sin2*h_p);

			d *= -1;

			hcFloat rhs		= 2	/ (r*h_r)		* (   psi_ipjk  * dr_m/dr_p  - psi_imjk * dr_p/dr_m)
							+ 1	/  h_r			* (  (psi_ippjk * dr_m/dr_pp + psi_ipjk * (dr_m*dr_pp/(dr_p*dr_p) - dr_m/dr_pp)) 	   	/ (dr_pp + dr_p)
												   + (psi_immjk * dr_p/dr_mm + psi_imjk * (dr_p*dr_mm/(dr_m*dr_m) - dr_p/dr_mm)) 		/ (dr_mm + dr_m)
												   - (psi_ipjk  * (dr_m*dr_m/(dr_p*dr_p) - 1) + psi_imjk * (dr_p*dr_p/(dr_m*dr_m) - 1)) /  h_r)
							+ 1 / (r*r*h_t)		* (  cos(theta)/sin(theta) * (dt_m/dt_p * psi_ijpk - dt_p/dt_m * psi_ijmk)
												   + (psi_ijppk * dt_m/dt_pp + psi_ijpk * (dt_m*dt_pp/(dt_p*dt_p) - dt_m/dt_pp)) 	   	/ (dt_pp + dt_p)
												   + (psi_ijmmk * dt_p/dt_mm + psi_ijmk * (dt_p*dt_mm/(dt_m*dt_m) - dt_p/dt_mm))		/ (dt_mm + dt_m)
												   - (psi_ijpk  * (dt_m*dt_m/(dt_p*dt_p) - 1) + psi_ijmk * (dt_p*dt_p/(dt_m*dt_m) - 1)) /  h_t)
							+ 1	/ (r*r*sin2*h_p)* (  (psi_ijkpp * dp_m/dp_pp + psi_ijkp * (dp_m*dp_pp/(dp_p*dp_p) - dp_m/dp_pp)) 	   	/ (dp_pp + dp_p)
												   + (psi_ijkmm * dp_p/dp_mm + psi_ijkm * (dp_p*dp_mm/(dp_m*dp_m) - dp_p/dp_mm)) 		/ (dp_mm + dp_m)
												   - (psi_ijkp  * (dp_m*dp_m/(dp_p*dp_p) - 1) + psi_ijkm * (dp_p*dp_p/(dp_m*dp_m) - 1)) /  h_p);
			//-----
			/*////-- higher accuracy here
			//-----

			hcFloat ar	= 1 / (dr_p*dr_m*dr_m + dr_p*dr_p*dr_m);
			hcFloat at	= 1 / (dt_p*dt_m*dt_m + dt_p*dt_p*dt_m);
			hcFloat ap	= 1 / (dp_p*dp_m*dp_m + dp_p*dp_p*dp_m);

			hcFloat fr	= 2 * ar * (dr_p*dr_p - dr_m*dr_m) / r - 2 * ar * (dr_p + dr_m);
			hcFloat ft	= __cosf(theta) * at * (dt_p*dt_p - dt_m*dt_m) / (__sinf(theta) * r*r) - 2 * at * (dt_p + dt_m) / (r*r);
			hcFloat fp	= - 2 * ap * (dp_p + dp_m) / (r*r*sin2);

			hcFloat d	= -(fr + ft + fp);

			hcFloat rhs	= 2 * ar   			* (psi_imjk * dr_p * (1 - dr_p / r) 					+ psi_ipjk * dr_m * (1 + dr_m / r)						)
						+ at/(r*r) 			* (psi_ijmk * dt_p * (2 - __cosf(theta)*dt_p / __sinf(theta))	+ psi_ijpk * dt_m * (2 + __cosf(theta)*dt_m / __sinf(theta))	)
						+ 2*ap/(r*r*sin2)	* (psi_ijkm * dp_p 										+ psi_ijkp * dp_m										);//*/

			//*/

			gr.setTemp(ind, rhs / d);
			gr.setRelError(ind, fabsf((gr.getTemp(ind) - gr.getPsi(ind)) / gr.getTemp(ind)));
		}
   }
}

__global__ void
__launch_bounds__(1024, 1)
kernel_setPsi(LaplaceSolver &solver, uint prevBlocks)
{
    uint idx 		= (blockIdx.x + prevBlocks) * blockDim.x + threadIdx.x;
    uint numR       = solver.grid->numR;
    uint numTheta   = solver.grid->numTheta;
    uint numPhi     = solver.grid->numPhi;
    // compute position in grid
    uint i 			= idx%numR;
    uint j 			= ((idx-i)/numR)%numTheta;
    uint k 			= ((idx-i - j * numR) / numR) / numTheta;

    if(idx<numTheta*numPhi*numR && i<numR-1)
    {
    	uint ind 	= solver.grid->getIndex(i,j,k);
    	hcFloat tmp	= solver.grid->getTemp(ind);
        solver.grid->setPsi(ind, tmp);
    }
}

/**********************************************************************************************************************
 *  iterate elliptic
 **********************************************************************************************************************/
__global__ void
__launch_bounds__(1024, 1)
kernel_iterate_elliptic(LaplaceSolver &solver, uint prevBlocks)
{
    uint idx 		= (blockIdx.x + prevBlocks) * blockDim.x + threadIdx.x;

	SphericalGrid &gr = *solver.grid;

    uint numR       = gr.numR;
    uint numT	   	= gr.numTheta;
    uint numP   	= gr.numPhi;

    Vec3D *pos		= gr.getPosArray();
    hcFloat *psi	= gr.getPsiArray();

    // compute position in grid
    uint i = idx % numR;
    uint j = ((idx   - i) / numR) % numT;
    uint k = ((idx-i - j  * numR) / numR) / numT;

    uint ind = k * numR * numT + j * numR + i;


    if(idx < numT * numP * numR && i > 0)
    {
    	if(i==numR-1)
        {
    		gr.setTemp(ind, 0.0);
    		gr.setRelError(ind, 0.0);
        }
        else
        {
        	uint km				= k==0		? numP-1: k-1;
        	uint kp				= k==numP-1	? 0		: k+1;

        	uint ind_imjk		= 				  	  gr.getIndex(i-1,	j,	 k );
			uint ind_ipjk   	= 				  	  gr.getIndex(i+1,	j,	 k );

			uint ind_ijmk		= j==0		? 0		: gr.getIndex(i,	j-1, k );
			uint ind_ijpk		= j==numT-1	? 0		: gr.getIndex(i,	j+1, k );

			uint ind_ijkm  		= 				  	  gr.getIndex(i,	j,	 km);
			uint ind_ijkp  		= 				  	  gr.getIndex(i,	j,	 kp);

			uint ind_imjmk		= j==0		? 0		: gr.getIndex(i-1,	j-1, k );
			uint ind_imjpk		= j==numT-1	? 0		: gr.getIndex(i-1,	j+1, k );
			uint ind_ipjmk		= j==0		? 0		: gr.getIndex(i+1,	j-1, k );
			uint ind_ipjpk		= j==numT-1	? 0		: gr.getIndex(i+1,	j+1, k );

			uint ind_imjkm		= 				  	  gr.getIndex(i-1,	j,	 km);
			uint ind_imjkp		= 				  	  gr.getIndex(i-1,	j,	 kp);
			uint ind_ipjkm		= 				  	  gr.getIndex(i+1,	j,	 km);
			uint ind_ipjkp		= 				  	  gr.getIndex(i+1,	j,	 kp);

			uint ind_ijmkm		= j==0		? 0		: gr.getIndex(i,	j-1, km);
			uint ind_ijmkp		= j==0		? 0		: gr.getIndex(i,	j-1, kp);
			uint ind_ijpkm		= j==numT-1	? 0		: gr.getIndex(i,	j+1, km);
			uint ind_ijpkp		= j==numT-1	? 0		: gr.getIndex(i,	j+1, kp);

			hcFloat psi_imjk	= psi[ind_imjk];
			hcFloat psi_ipjk	= psi[ind_ipjk];

			hcFloat psi_pole	= 0.0;
			hcFloat psi_pole_m	= 0.0;
			hcFloat psi_pole_p	= 0.0;
			for(uint l=0; l<numP; ++l)
			{
				uint ind_pole 	= l * numR * numT + j * numR + i;
				uint ind_pole_m	= l * numR * numT + j * numR + i-1;
				uint ind_pole_p	= l * numR * numT + j * numR + i+1;

				psi_pole 		+= psi[ind_pole];
				psi_pole_m 		+= psi[ind_pole_m];
				psi_pole_p 		+= psi[ind_pole_p];
			}
			psi_pole 			/= numP;
			psi_pole_m 			/= numP;
			psi_pole_p 			/= numP;

			hcFloat psi_ijmk	= j==0			? psi_pole 		: psi[ind_ijmk];
			hcFloat psi_ijpk	= j==numT-1 	? psi_pole 		: psi[ind_ijpk];

			hcFloat psi_ijkm	= psi[ind_ijkm];
			hcFloat psi_ijkp	= psi[ind_ijkp];

			hcFloat psi_imjmk	= j==0			? psi_pole_m	: psi[ind_imjmk];
			hcFloat psi_imjpk	= j==numT-1		? psi_pole_m	: psi[ind_imjpk];
			hcFloat psi_ipjmk	= j==0			? psi_pole_p	: psi[ind_ipjmk];
			hcFloat psi_ipjpk	= j==numT-1		? psi_pole_p	: psi[ind_ipjpk];

			hcFloat psi_imjkm	= psi[ind_imjkm];
			hcFloat psi_imjkp	= psi[ind_imjkp];
			hcFloat psi_ipjkm	= psi[ind_ipjkm];
			hcFloat psi_ipjkp	= psi[ind_ipjkp];

			hcFloat psi_ijmkm	= j==0			? psi_pole		: psi[ind_ijmkm];
			hcFloat psi_ijmkp	= j==0			? psi_pole		: psi[ind_ijmkp];
			hcFloat psi_ijpkm	= j==numT-1		? psi_pole		: psi[ind_ijpkm];
			hcFloat psi_ijpkp	= j==numT-1		? psi_pole		: psi[ind_ijpkp];

			hcFloat s_ijk		= gr.s_ijk[ind];
			hcFloat s_imjk 		= gr.s_imjk[ind];
			hcFloat s_ipjk 		= gr.s_ipjk[ind];

			hcFloat s_ijmk 		= gr.s_ijmk[ind];
			hcFloat s_ijpk 		= gr.s_ijpk[ind];

			hcFloat s_ijkm 		= gr.s_ijkm[ind];
			hcFloat s_ijkp 		= gr.s_ijkp[ind];

			hcFloat s_imjmk 	= gr.s_imjmk[ind];
			hcFloat s_imjpk 	= gr.s_imjpk[ind];
			hcFloat s_ipjmk 	= gr.s_ipjmk[ind];
			hcFloat s_ipjpk 	= gr.s_ipjpk[ind];

			hcFloat s_imjkm 	= gr.s_imjkm[ind];
			hcFloat s_imjkp 	= gr.s_imjkp[ind];
			hcFloat s_ipjkm 	= gr.s_ipjkm[ind];
			hcFloat s_ipjkp 	= gr.s_ipjkp[ind];

			hcFloat s_ijmkm 	= gr.s_ijmkm[ind];
			hcFloat s_ijmkp 	= gr.s_ijmkp[ind];
			hcFloat s_ijpkm 	= gr.s_ijpkm[ind];
			hcFloat s_ijpkp 	= gr.s_ijpkp[ind];

            hcFloat rhs		= psi_imjmk * s_imjmk	+ psi_imjpk * s_imjpk
            				+ psi_ipjmk * s_ipjmk	+ psi_ipjpk * s_ipjpk
            				+ psi_imjkm * s_imjkm	+ psi_imjkp * s_imjkp
            				+ psi_ipjkm * s_ipjkm	+ psi_ipjkp * s_ipjkp
            				+ psi_ijmkm * s_ijmkm 	+ psi_ijmkp * s_ijmkp
            				+ psi_ijpkm * s_ijpkm 	+ psi_ijpkp * s_ijpkp
            				+ psi_imjk  * s_imjk 	+ psi_ipjk  * s_ipjk
            				+ psi_ijmk  * s_ijmk	+ psi_ijpk  * s_ijpk
            				+ psi_ijkm  * s_ijkm	+ psi_ijkp  * s_ijkp;

            /*											// relaxation
			hcFloat newVal	= -1/s_ijk * rhs;
			hcFloat psi		= gr.getPsi(ind);
			hcFloat relax	= 0.2;						// smaller than 1 -> under relaxation, otherwise -> overrelaxation
			hcFloat retval	= (1.0-relax)*psi + relax * newVal;
			/*/
            hcFloat retval	= -1/s_ijk * rhs;//*/

            gr.setTemp(ind, retval);
            gr.setRelError(ind, fabsf((gr.getTemp(ind) - gr.getPsi(ind)) / gr.getTemp(ind)));

            if (isnan(gr.getTemp(ind)) || isinf(gr.getTemp(ind)))
				printf("\nNaN/INF idx: %u, i/j/k: %u(%u)/%u(%u)/%u(%u)",	idx, i, numR, j, numT, k, numP);
        }
   }
}


//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          LaplaceSolver member functions
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void LaplaceSolver::initNULL_GPU(){

    d_image                 = NULL;
    d_surfaceShapeLookup    = NULL;
}

void LaplaceSolver::clear_GPU(){

    if(d_image != NULL)			        cudaFree(d_image);
    if(d_surfaceShapeLookup != NULL)	cudaFree(d_surfaceShapeLookup);
    if(grid != NULL)			    	grid->clear_GPU();

    checkCUDAerror("LaplaceSolver::clear_GPU");
    onDevice = false;
    initNULL_GPU();
}

bool LaplaceSolver::pushToDevice(hcImageFITS &photBoundary)
{
    // TODO test if there is enough memory available
    clear_GPU();

    cudaMalloc((void **)&d_solver, sizeof(LaplaceSolver));
    cudaMemcpy(d_solver, this, sizeof(LaplaceSolver), cudaMemcpyHostToDevice);

    hcFloat *imgData = new hcFloat[photBoundary.width * photBoundary.height];
    for(uint x=0; x<photBoundary.width; ++x)
    	for(uint y=0; y<photBoundary.height; ++y)
    		imgData[y*photBoundary.width +x] = photBoundary(x,y);

    cudaMalloc((void **)&d_image,   photBoundary.width * photBoundary.height * sizeof(hcFloat  ));
    cudaMemcpy(d_image, imgData, 	photBoundary.width * photBoundary.height * sizeof(hcFloat  ), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_solver->d_image), &d_image,                                        sizeof(hcFloat *), cudaMemcpyHostToDevice);

    grid->pushToGPU();

    if(!grid->isElliptical())	cudaMemcpy(&(d_solver->grid), &(grid->d_memStruct), sizeof(SphericalGrid*), cudaMemcpyHostToDevice);
	else						cudaMemcpy(&(d_solver->grid), &(((EllipticalGrid*)grid)->d_memStruct), sizeof(EllipticalGrid*), cudaMemcpyHostToDevice);

    uint null 	= 0;
    cudaMemcpy(&(d_solver->d_solver), &null, sizeof(uint *), cudaMemcpyHostToDevice);

    onDevice = true;
    checkCUDAerror("LaplaceSolver::pushToDevice");
    return true;
}

void LaplaceSolver::computeLowerBoundaryPsi_GPU()
{
    uint TPB                = MAX_TPB;
    uint numBlocks_part1    = grid->numTheta   * grid->numPhi / TPB + ((grid->numTheta * grid->numPhi)%TPB == 0 ? 0 : 1);

   	kernel_computeLowerBoundary<<<numBlocks_part1, TPB, 0>>>(*this->d_solver);
    cudaDeviceSynchronize();
    checkCUDAerror("Kernel_computeLowerBoundary!\n");
}

void LaplaceSolver::iterate_GPU()
{
    uint TPB 			= MAX_TPB;
    uint BPI 			= MAX_NUM_BLOCK;
    uint numGridPoints  = grid->numTheta * grid->numPhi * grid->numR;
    uint numBlocks      = numGridPoints / TPB + (numGridPoints%TPB == 0 ? 0 : 1);
    uint numSubdiv      = numBlocks     / BPI + (numBlocks%BPI     == 0 ? 0 : 1);

	for(uint num=0;num<numSubdiv;++num)
	{
		checkCUDAerror(__FILE__, __LINE__);
#ifdef SPHERICUNITVEC
		if(!grid->isElliptical())	kernel_iterate_spheric<<< BPI, TPB, 0>>>(*this->d_solver, num * BPI);
		else
#endif
									kernel_iterate_elliptic<<<BPI, TPB, 0>>>(*this->d_solver, num * BPI);
		cudaDeviceSynchronize();
		checkCUDAerror(__FILE__, __LINE__);

		kernel_setPsi<<<BPI, TPB, 0>>>(*this->d_solver, num * BPI);
		cudaDeviceSynchronize();
		checkCUDAerror(__FILE__, __LINE__);
	}
}

int LaplaceSolver::computeSolution_GPU(uint maxNumIterations, float errorThreshold, hcImageFITS &photBoundary, bool verbose)
{
	printStdOutMess(__FILE__, __LINE__, "CUDA-version of PFSS solver started");
    grid->clearValues();

    // first, compute the lower boundary potential via B_r = - dPsi / dr -> Psi_0 = Phsi_1 + dr * B_r
    // => Psi_0(t=0) = dr * B_l * csc(theta)
    pushToDevice(photBoundary);
    computeLowerBoundaryPsi_GPU();
  	grid->extract_relError();

    // now perform the main loop. Compute next value for Psi and abort
    // if the difference to the step before is below threshold
    hcFloat maxError		= errorThreshold;
    uint counter            = 0;
    uint threshCheckSteps   = 100;

    uint loopnum 			= 0;
    while(loopnum < maxNumIterations)
    {
        ++counter;
        iterate_GPU();
        computeLowerBoundaryPsi_GPU();

        if(counter == threshCheckSteps)
        {
			grid->extract_relError();

            counter 	= 0;
            maxError 	= 0.0;

            for(uint i=1;i<grid->numR;++i)
                for(uint j=0;j<grid->numTheta;++j)
                    for(uint k=0;k<grid->numPhi;++k)
                    {
                        uint ind = grid->getIndex(i,j,k);
                        if(grid->getRelError(ind) > maxError)	maxError = grid->getRelError(ind);
                    }
        }

#ifdef VERBOSE
			cout << "\r                                                                    \r";
			cout << "\tstep " << loopnum << " / " << maxNumIterations << ", error / threshold: " << toStr(maxError) << " / " << toStr(errorThreshold);
#endif

        if(maxError < errorThreshold)
		{
#ifdef VERBOSE
        	cout << "\n";
#endif
			break;
		}
        ++loopnum;
    }

   	grid->pullFromGPU();

   	for(uint i=0;i<grid->numR;++i)
        for(uint j=0;j<grid->numTheta;++j)
            for(uint k=0;k<grid->numPhi;++k)
            {
            	uint ind 		= grid->getIndex(i,j,k);
            	Vec3D pos		= grid->pos[ind];
				grid->B[ind]	= grid->getBFromPsi(i,j,k);
            }

    clear_GPU();

    solutionComputed = true;

    return loopnum;
}

