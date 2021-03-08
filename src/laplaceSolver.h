#ifndef LAPLACESOLVER_H
#define LAPLACESOLVER_H

#include "engine/hcSet.h"
#include "engine/hcImageFITS.h"

#include "src/grids.h"

#ifdef GUI
#include "src/sourceSurfShape.h"
#endif

class LaplaceSolver {
public:
	SphericalGrid *grid;
    bool solutionComputed;        				/*!< \brief has some computation been done on this?     		*/

    LaplaceSolver();                            /*!< \brief std constructor                             		*/
    LaplaceSolver(const LaplaceSolver &solver); /*!< \brief cpy constructor                             		*/
    ~LaplaceSolver();

    LaplaceSolver &operator=(const LaplaceSolver &other);

    bool init(bool sinLatGrid, hcFloat maxSinLat, hcFloat r_ss, uint numR, hcFloat ell);

    void iterateSpheric();
    	/*!< iteration procedure only for spherical grid with unit base vectors									*/

    void iterateElliptic();
    	/*!< iteration procedure for general spheric/elliptic grid with general (non-unit) base vectors			*/

    void iterate();
    	/*!< iteration switcher for spheric/elliptic grid with unit/non-unit base vectors						*/

    int computeSolution(hcImageFITS &photImage);
        /*!< \brief computes solution of Laplace's equation via finite differences                      		*/

#ifdef CUDA
    LaplaceSolver *d_solver;
    SphericalGrid *d_grid;
    hcFloat *d_image;           /*!< \brief phtotospheric boundary on GPU                               		*/
    int *d_surfaceShapeLookup;  /*!< \brief lookup for arbitrary source surface shape on GPU            		*/
    bool onDevice;              /*!< \brief arrays already allocated and ready to compute?              		*/

    void initNULL_GPU();
    void clear_GPU();

    bool pushToDevice(hcImageFITS &photBoundary);
        /*!< \brief allocates memory on device and copies data from host into it                        		*/

    int computeSolution_GPU(uint maxNumIterations, float errorThreshold, hcImageFITS &photBoundary, bool verbose=false);
        /*!< \brief computes finite difference laplace equation on GPU                                  		*/

    void computeLowerBoundaryPsi_GPU();
        /*!< \brief computes the lower boundary via photospheric image data stored in device RAM        		*/

    void iterate_GPU();
        /*!< \brief iteration step for solving Laplace's equation on GPU                                		*/
#endif

    void dump() const;

private:

    void initNULL();
    void clear();

    void initNULL_CPU();
    void clear_CPU();

    int computeSolution_CPU(uint maxNumIterations, hcFloat threshold, hcImageFITS &photBoundary);
        /*!< \brief computes finite difference laplace equation on CPU                                  		*/

    void computeLowerBoundaryPsi_CPU(hcFloat *image);
        /*!< \brief computes the lower boundary via photospheric image data stored in host ram          		*/

    void computeLowerBoundary_CPU(hcImageFITS &boundary);

    hcFloat getMaxIncrement();
};

#endif // LAPLACESOLVER_H
