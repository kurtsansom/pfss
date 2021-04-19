#ifndef LAPLACESOLVER_H
#define LAPLACESOLVER_H

#include "engine/hcSet.h"
#include "engine/hcImageFITS.h"

#include "src/grids.h"

#ifdef GUI
#include "src/sourceSurfShape.h"
#endif

struct threadParameterLaplace
{
    uint threadID;
    pthread_mutex_t *runningMutex;
    SphericalGrid *grid;
    volatile _Atomic_word *numRunningThreads;
    volatile _Atomic_word *threadRunning;

    uint idx;	/*!< \brief id of grid point to be computed, as in CUDA version								*/
    uint gppt;	/*!< \brief grid points per thread															*/

    threadParameterLaplace(){init();}

    void init(	uint threadID=0, volatile _Atomic_word *numRunningThreads=NULL, pthread_mutex_t *runningMutex=NULL,
        		volatile _Atomic_word *threadRunning=NULL, SphericalGrid *grid=NULL)
	{
		this->threadID      	= threadID;
		this->runningMutex		= runningMutex;
		this->grid     			= grid;
		this->numRunningThreads	= numRunningThreads;
		this->threadRunning 	= threadRunning;
		this->idx				= 0;
		this->gppt				= 0;
	}

    void set(const uint &idx, const uint &gppt)
    {
    	this->idx  = idx;
    	this->gppt = gppt;
    }
};
typedef struct threadParameterLaplace threadParamLaplace;

class LaplaceSolver {
public:
	SphericalGrid *grid;
    bool solutionComputed;        				/*!< \brief has some computation been done on this?     						*/

    LaplaceSolver();                            /*!< \brief std constructor                             						*/
    LaplaceSolver(const LaplaceSolver &solver); /*!< \brief cpy constructor                             						*/
    ~LaplaceSolver();

    LaplaceSolver &operator=(const LaplaceSolver &other);

    bool init(bool sinLatGrid, hcFloat maxSinLat, hcFloat r_ss, uint numR, hcFloat ell);

    void iterateSpheric();
    	/*!< iteration procedure only for spherical grid with unit base vectors													*/

    static void *iterateElliptic_threadEntry(void *parameter);
    	/*!< \brief thread entry point forelliptic iteration																	*/

    void iterateElliptic_MT();
    	/*!< iteration procedure for general spheric/elliptic grid with general (non-unit) base vectors, multi-threaded			*/

    void iterateElliptic_ST();
    	/*!< iteration procedure for general spheric/elliptic grid with general (non-unit) base vectors, single-threaded		*/

    void iterateElliptic();
    	/*!< iteration procedure for general spheric/elliptic grid with general (non-unit) base vectors							*/

    void iterate_CPU();
    	/*!< iteration switcher for spheric/elliptic grid with unit/non-unit base vectors										*/

    int computeSolution(hcImageFITS &photImage);
        /*!< \brief computes solution of Laplace's equation via finite differences                      						*/

#ifdef CUDA
    LaplaceSolver *d_solver;
    SphericalGrid *d_grid;
    hcFloat *d_image;           /*!< \brief phtotospheric boundary on GPU                               						*/
    int *d_surfaceShapeLookup;  /*!< \brief lookup for arbitrary source surface shape on GPU            						*/
    bool onDevice;              /*!< \brief arrays already allocated and ready to compute?              						*/

    void initNULL_GPU();
    void clear_GPU();

    bool pushToDevice(hcImageFITS &photBoundary);
        /*!< \brief allocates memory on device and copies data from host into it                        						*/

    int computeSolution_GPU(uint maxNumIterations, float errorThreshold, hcImageFITS &photBoundary, bool verbose=false);
        /*!< \brief computes finite difference laplace equation on GPU                                  						*/

    void computeLowerBoundaryPsi_GPU();
        /*!< \brief computes the lower boundary via photospheric image data stored in device RAM        						*/

    void iterate_GPU();
        /*!< \brief iteration step for solving Laplace's equation on GPU                                						*/
#endif

    void dump() const;

private:

    void initNULL();
    void clear();

    void initNULL_CPU();
    void clear_CPU();

    int computeSolution_CPU(uint maxNumIterations, hcFloat threshold, hcImageFITS &photBoundary);
        /*!< \brief computes finite difference laplace equation on CPU                                  						*/

    void computeLowerBoundary_CPU(hcImageFITS &boundary);
    	/*!< \brief computes the lower boundary via photospheric image data stored in host ram          						*/

    hcFloat getMaxIncrement();
};

#endif // LAPLACESOLVER_H
