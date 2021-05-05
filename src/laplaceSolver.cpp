#include "src/laplaceSolver.h"
#include "src/ellipticalGrid.h"

#include "engine/hcImageFITS.h"
#include "engine/hcImage.h"
#include "engine/math/hcFunction.h"

#ifdef GUI
#include "src/fluxtube.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <fstream>


//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          LaplaceSolver
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

LaplaceSolver::LaplaceSolver()
{
	initNULL();
}

LaplaceSolver::LaplaceSolver(const LaplaceSolver &solver)
{
    initNULL_CPU();
    printf("ERROR! LaplaceSolver cpy constructor not implemented!\n");
    exit(1);
}


LaplaceSolver::~LaplaceSolver()
{
	clear();
}

LaplaceSolver &LaplaceSolver::operator=(const LaplaceSolver &other)
{
    if(this == &other)
        return *this;

    printf("ERROR! LaplaceSolver::assignment operator TODO!");
    exit(1);

    return *this;
}

void LaplaceSolver::initNULL()
{
    initNULL_CPU();
#ifdef CUDA
    initNULL_GPU();
#endif
}

void LaplaceSolver::clear()
{
#ifdef CUDA
    clear_GPU();
#endif
    clear_CPU();
}

void LaplaceSolver::initNULL_CPU()
{
    solutionComputed    = false;
    this->grid 			= NULL;
}

void LaplaceSolver::clear_CPU()
{
	delete this->grid;
    initNULL_CPU();
}

bool LaplaceSolver::init(bool sinLatGrid, hcFloat maxSinLat, hcFloat r_ss, uint numR, hcFloat ell)
{
    clear_CPU();

    if(ell==0.0)
    {
    	printErrMess(__FILE__, __LINE__, "ellipticity 0.0 not allowed");
    	return false;
    }

    this->grid = ell == 1.0 ? new SphericalGrid : new EllipticalGrid;
    this->grid->init(sinLatGrid, maxSinLat, -maxSinLat, r_sol, r_ss, numR, true, ell);

    return true;
}

void SphericalGrid::iterateElliptic_gridPoint(uint i, uint j, uint k)
{
	SphericalGrid &gr 	= *this;
	uint numR       	= gr.numR;
	uint numT	   		= gr.numTheta;
	uint numP   		= gr.numPhi;
	uint ind 			= k * numR * numT + j * numR + i;
	hcFloat *psi		= gr.getPsiArray();

	if(i==numR-1)
	{
		gr.setTemp(ind, 0.0);
		gr.setRelError(ind, 0.0);
	}
	else if(i > 0)
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

		hcFloat retval	= -1/s_ijk * rhs;

		gr.setTemp(ind, retval);
		gr.setRelError(ind, fabsf((gr.getTemp(ind) - gr.getPsi(ind)) / gr.getTemp(ind)));

		if (std::isnan(gr.getTemp(ind)) || std::isinf(gr.getTemp(ind)))
		{
			printErrMess(__FILE__, __LINE__, "solver got nan/inf value");
			exit(1);
		}
	}
}

/*! entry function for multicore PFSS iteration function
 */
void *LaplaceSolver::iterateElliptic_threadEntry(void *parameter)
{
	threadParamLaplace *param 	= (threadParamLaplace*)parameter;
	SphericalGrid *grid			= param->grid;
	uint id 					= param->idx;
	bool debug 					= false;

	if(debug)
	{
		cout << "------------------------------------------------------\n";
		cout << "Thread:     " << param->threadID << "\n";
		cout << "NumRunning: " << *param->numRunningThreads << "\n";
		fflush(stdout);
	}

	SphericalGrid &gr 	= *grid;
	uint numR       	= gr.numR;
	uint numT	   		= gr.numTheta;
	uint numP   		= gr.numPhi;

	uint idx = id;
	while( (idx < id + param->gppt) && (idx < numR*numT*numP) )
	{
		// compute position in grid
		uint i = idx % numR;
		uint j = ((idx   - i) / numR) % numT;
		uint k = ((idx-i - j  * numR) / numR) / numT;

		grid->iterateElliptic_gridPoint(i, j, k);

		++idx;
	}

	pthread_mutex_lock(param->runningMutex);
	--(*param->numRunningThreads);
	*param->threadRunning = 0;
	pthread_mutex_unlock(param->runningMutex);
	pthread_exit(NULL);
}

void LaplaceSolver::iterateElliptic_MT()
{
	SphericalGrid &gr 	= *grid;
	uint numR      		= gr.numR;
	uint numT	   		= gr.numTheta;
	uint numP   		= gr.numPhi;
	uint numGridPoints	= numR*numT*numP;
	const uint gppt		= 1000;				// grid points to be computed per thread in each step

	uint numSubdiv      = numGridPoints / gppt + (numGridPoints % gppt == 0 ? 0 : 1);

	bool *workedUpon 	= new bool[numSubdiv];

	for(uint i=0; i<numSubdiv; ++i)
		workedUpon[i] = false;

	pthread_t 				threads[		NUMTHREADS];
	volatile _Atomic_word 	threadRunning[	NUMTHREADS];
	threadParamLaplace		tParams[		NUMTHREADS];

	pthread_mutex_t runningMutex 	= PTHREAD_MUTEX_INITIALIZER;
	volatile _Atomic_word numRunningThreads = 0;

	for(uint i=0; i<NUMTHREADS; ++i)
	{
		threadRunning[i] = 0;
		tParams[i].init(i, &numRunningThreads, &runningMutex, &threadRunning[i], grid);
	}

	bool workLeft = true;

	while(workLeft)
	{
		workLeft = false;
		if(numRunningThreads < NUMTHREADS-1) // -1 should not be here, but then "if(j==numMaxThreads-1)" is triggered sometimes
		{
			//for(uint i=0; i<numR*numT*numP; ++i)
			for(uint i=0; i<numSubdiv; ++i)
			{
				bool breaker = false;

				if(!workedUpon[i])
				{
					workedUpon[i] 		= true;

					for(uint j=0; j<NUMTHREADS; ++j)
					{
						if(threadRunning[j] == 0)
						{
							pthread_mutex_lock(&runningMutex);
							++numRunningThreads;
							threadRunning[j] = 1;
							pthread_mutex_unlock(&runningMutex);
							tParams[j].set(i*gppt, gppt);

							int rc = pthread_create(&threads[j], NULL, LaplaceSolver::iterateElliptic_threadEntry, (void*)(&tParams[j]));

							if(rc)
							{
								printErrMess(__FILE__, __LINE__, "return code from pthread_create() is " + to_string(rc));
								exit(1);
							}
							pthread_detach(threads[j]);
							breaker = true;
							break;
						}

						if(j==NUMTHREADS-1)
						{
							printErrMess(__FILE__, __LINE__, "No free thread found! (you should not be able to see this. If you do, I fucked up.... sorry)");
							cerr << "numRunningThreads: " << numRunningThreads 	<< "\n";
							cerr << "numMaxThreads:     " << NUMTHREADS 		<< "\n";
							for(uint k=0; k<NUMTHREADS; ++k)
								cerr << "Thread " << k << " running: " << threadRunning[k] << "\n";
							exit(1);
						}
					}
				}

				if(breaker)	break;
			}
		}

		for(uint i=0;i<numSubdiv;++i)
			if(!workedUpon[i])
			{
				workLeft = true;
				break;
			}
	}

	while(numRunningThreads > 0)	usleep(100);

	delete [] workedUpon;
}

void LaplaceSolver::iterateElliptic_ST()
{
	SphericalGrid &gr 	= *grid;
	uint numR      		 = gr.numR;
	uint numT	   		= gr.numTheta;
	uint numP   		= gr.numPhi;

    for(uint i=1; i<numR; ++i)
    	for(uint j=0; j<numT; ++j)
    		for(uint k=0; k<numP; ++k)
    			gr.iterateElliptic_gridPoint(i, j, k);
}

void LaplaceSolver::iterateElliptic()
{
#if NUMTHREADS > 1
	iterateElliptic_MT();
#else
	iterateElliptic_ST();
#endif

	for(uint i=1; i<grid->numR-1; ++i)
		for(uint j=0; j<grid->numTheta; ++j)
			for(uint k=0; k<grid->numPhi; ++k)
			{
				uint ind = grid->getIndex(i,j,k);
				grid->setPsi(ind, grid->getTemp(ind));
			}
}

void LaplaceSolver::iterateSpheric()
{
	uint numR       = grid->numR;
	uint numTheta   = grid->numTheta;
	uint numPhi     = grid->numPhi;

	Vec3D *pos		= grid->getPosArray();
	hcFloat *psi	= grid->getPsiArray();

    for(uint i=1; i < numR; ++i)
        for(uint j=0; j < numTheta; ++j)
            for(uint k=0; k < numPhi; ++k)
            {
            	uint ind 			= k * numR * numTheta + j * numR + i;

				uint ind_r_m    	= (i==0			? k * numR * numTheta + j * numR + i 		: k * numR * numTheta + j * numR + i - 1);
				uint ind_r_p    	= (i==numR-1	? k * numR * numTheta + j * numR + i 		: k * numR * numTheta + j * numR + i + 1);

				uint ind_theta_m	= k * numR * numTheta + (j-1) * numR + i;
				uint ind_theta_p	= k * numR * numTheta + (j+1) * numR + i;

				uint ind_phi_m  	= (k==0		   		? (numPhi-1) * numR * numTheta + j * numR + i : (k-1) * numR * numTheta + j * numR + i);
				uint ind_phi_p  	= (k==numPhi-1 		?                                j * numR + i : (k+1) * numR * numTheta + j * numR + i);


				hcFloat dr_p    	= 					  pos[ind_r_p][0] - pos[ind][0];
				hcFloat dr_m    	= 					  pos[ind][0]     - pos[ind_r_m][0];

				hcFloat dt_p		= (j==numTheta-1	? PI - pos[ind][1] 		  	: pos[ind_theta_p][1]  - pos[ind][1]);
				hcFloat dt_m		= (j==0				? pos[ind][1]			  	: pos[ind][1]          - pos[ind_theta_m][1]);

				hcFloat dp_p		= (k==numPhi-1 		? pos[ind_phi_p][2]  + 2*PI - pos[ind][2] 		: pos[ind_phi_p][2]   - pos[ind][2]);
				hcFloat dp_m		= (k==0				? 2*PI - pos[ind_phi_m][2]  + pos[ind][2] 		: pos[ind][2]         - pos[ind_phi_m][2]);

				hcFloat r			= pos[ind][0];
				hcFloat theta		= pos[ind][1];
				hcFloat sin2		= sin(theta) * sin(theta);

				hcFloat psi_pole	= 0.0;

				for(uint l=0; l<numPhi; ++l)
				{
					uint ind_pole 	= l * numR * numTheta + j 	  * numR + i;
					psi_pole 		+= psi[ind_pole];
				}
				psi_pole 			/= numPhi;

				hcFloat psi_ipjk	= psi[ind_r_p];
				hcFloat psi_imjk	= psi[ind_r_m];

				hcFloat psi_ijpk	= ( j==numTheta-1 	? psi_pole 		: psi[ind_theta_p]);
				hcFloat psi_ijmk	= ( j==0			? psi_pole 		: psi[ind_theta_m]);

				hcFloat psi_ijkp	= psi[ind_phi_p];
				hcFloat psi_ijkm	= psi[ind_phi_m];

				hcFloat ar	= 1 / (dr_p*dr_m*dr_m + dr_p*dr_p*dr_m);
				hcFloat at	= 1 / (dt_p*dt_m*dt_m + dt_p*dt_p*dt_m);

				hcFloat br	= 2 / (dr_p*dr_p*dr_m + dr_p*dr_m*dr_m);
				hcFloat bt	= 2 / (dt_p*dt_p*dt_m + dt_p*dt_m*dt_m);
				hcFloat bp	= 2 / (dp_p*dp_p*dp_m + dp_p*dp_m*dp_m);

				hcFloat fr	= 2 * ar * (dr_p*dr_p - dr_m*dr_m) / r - br * (dr_p + dr_m);
				hcFloat ft	= cos(theta) * at * (dt_p*dt_p - dt_m*dt_m) / (sin(theta) * r*r) - bt * (dt_p + dt_m) / (r*r);
				hcFloat fp	= -bp * (dp_p + dp_m) / (r*r*sin2);

				hcFloat d	= -(fr + ft + fp);

				hcFloat rhs	= psi_imjk * (-2 * ar * dr_p*dr_p / r + br * dr_p) 							+ psi_ipjk * (2 * ar * dr_p * dr_p / r + br * dr_m)
							+ psi_ijmk * (-cos(theta)*at*dt_p*dt_p / (sin(theta)*r*r) + bt*dt_p/(r*r))	+ psi_ijpk * (cos(theta)*at*dt_m*dt_m / (sin(theta)*r*r) + bt*dt_m/(r*r))
							+ psi_ijkm * (bp * dp_p / (r*r*sin2))										+ psi_ijkp * (bp * dp_m / (r*r*sin2));

				grid->setTemp(ind, rhs / d);
                grid->setRelError(ind, fabs((grid->getTemp(ind) - grid->getPsi(ind)) / grid->getTemp(ind)));
           }

    for(uint i=1; i<grid->numR-1; ++i)
        for(uint j=0; j<grid->numTheta; ++j)
            for(uint k=0; k<grid->numPhi; ++k)
            {
                uint ind = grid->getIndex(i,j,k);
                grid->setPsi(ind, grid->getTemp(ind));
            }
}

void LaplaceSolver::iterate_CPU()
{
#ifdef SPHERICUNITVEC
	if(!grid->isElliptical())
		iterateSpheric();
	else
#endif
		iterateElliptic();
}

void LaplaceSolver::computeLowerBoundary_CPU(hcImageFITS &boundary)
{
	EllipticalGrid &gr	= *(EllipticalGrid*)grid;
	uint numR			= gr.numR;
	uint numT			= gr.numTheta;
	uint numP			= gr.numPhi;

    for(uint j=0; j<numT; ++j)
        for(uint k=0; k<numP; ++k)
        {
            uint ind    	= k * numR * numT + j * numR;

            Vec3D pos		= gr.pos[ind];
			Vec3D pos_p		= gr.pos[ind+1];
			Vec3D pos_pp	= gr.pos[ind+2];
			hcFloat psi_p	= gr.psi[ind+1];
			hcFloat psi_pp	= gr.psi[ind+2];

			hcFloat B_l     = boundary(k, numT-1-j);
			hcFloat B_r     = B_l / sin(pos[1]); // radial approach

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

hcFloat LaplaceSolver::getMaxIncrement()
{
	hcFloat maxIncrement = 0.0;
	for(uint i=1;i<grid->numR;++i)
		for(uint j=0;j<grid->numTheta;++j)
			for(uint k=0;k<grid->numPhi;++k)
			{
				uint ind = grid->getIndex(i,j,k);
				if(grid->getRelError(ind) > maxIncrement)
					maxIncrement = grid->getRelError(ind);
			}
	return maxIncrement;
}

int LaplaceSolver::computeSolution(hcImageFITS &photBoundary)
{
	uint numMaxIterations	= 1000000;
	hcFloat threshold		= 1E-2;

	hcDate start, stop;
	start.setFromSystemTime();
	int retval;

#ifdef CUDA
    retval = computeSolution_GPU(numMaxIterations, threshold, photBoundary);
#else
    retval = computeSolution_CPU(numMaxIterations, threshold, photBoundary);
#endif
    stop.setFromSystemTime();
    uint seconds 	= (stop-start)/hcDate::facSec;
    printStdOutMess(__FILE__, __LINE__, "Computation time: " + to_string(seconds) + " s, solution steps: " + to_string(retval));
    return retval;
}

int LaplaceSolver::computeSolution_CPU(uint maxNumIterations, hcFloat threshold, hcImageFITS &photBoundary)
{
#if NUMTHREADS > 1
	printStdOutMess(__FILE__, __LINE__, "CPU-version (NUMTHREADS=" + to_string(NUMTHREADS) + ") of PFSS solver started");
#else
	printStdOutMess(__FILE__, __LINE__, "CPU-version (single-trheaded) of PFSS solver started");
#endif

    grid->clearValues();

    // first, compute the lower boundary potential via B_r = - dPsi / dr -> Psi_0 = Phsi_1 + dr * B_r
    // => Psi_0(t=0) = dr * B_l * csc(theta)
    computeLowerBoundary_CPU(photBoundary);

    // now perform the main loop. Compute next value for Psi and abort
    // if the difference to the step before is below threshold
    hcFloat maxError		= threshold;
    uint counter            = 0;
    uint threshCheckSteps   = 100;

    uint loopnum 			= 0;
    while(loopnum < maxNumIterations)
    {
        ++counter;
        iterate_CPU();
        computeLowerBoundary_CPU(photBoundary);

        if(counter == threshCheckSteps)
        {
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
			cout << "\r                                                                                                                           \r";
			cout << "\tstep " << loopnum << "/" << maxNumIterations << ", error/threshold: " << toStr(maxError) << "/" << toStr(threshold);fflush(stdout);
#endif

        if(maxError < threshold)
		{
#ifdef VERBOSE
        	cout << "\n";
#endif
			break;
		}
        ++loopnum;
    }

   	for(uint i=0;i<grid->numR;++i)
        for(uint j=0;j<grid->numTheta;++j)
            for(uint k=0;k<grid->numPhi;++k)
            {
            	uint ind 		= grid->getIndex(i,j,k);
            	Vec3D pos		= grid->pos[ind];
				grid->B[ind]	= grid->getBFromPsi(i,j,k);
            }

    solutionComputed = true;

    return loopnum;

	/*
	printStdOutMess(__FILE__, __LINE__, "CUDA-version of PFSS solver started");
    uint loopnum = 0;
    grid->clearValues();

    // first, compute the lower boundary potential via B_r = - dPsi / dr -> Psi_0 = Phsi_1 + dr * B_r
    // => Psi_0(t=0) = dr * B_l * csc(theta)
    computeLowerBoundary_CPU(photBoundary);

    // now perform the main loop. Compute next value for Psi and abort
    // if the difference to the step before is low
    double maxIncrement     = threshold;
    uint counter            = 0;
    uint errorCheckSteps    = 100;

    while(loopnum < maxNumIterations)
    {
        ++counter;
        iterate();
        computeLowerBoundary_CPU(photBoundary);

        if(counter >= errorCheckSteps)
        {
            counter 		= 0;
            maxIncrement	= getMaxIncrement();
        }

        if(maxIncrement < threshold)	break;
        ++loopnum;
    }

    for(uint i=0;i<grid->numR;++i)
		for(uint j=0;j<grid->numTheta;++j)
			for(uint k=0;k<grid->numPhi;++k)
			{
				uint ind 		= grid->getIndex(i,j,k);
				Vec3D pos		= grid->pos[ind];
				grid->B[ind]	= grid->getBFromPsi(i,j,k);
			}

    solutionComputed = true;

    return loopnum;//*/
}

void LaplaceSolver::dump() const{

    printf("Dumping LaplaceSolver:\n");
    printf("Solution computed: %u\n", solutionComputed);
    grid->dump();
}
