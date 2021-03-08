#ifndef MAGMAPPING_H
#define MAGMAPPING_H

#include "src/magline.h"
#include "src/laplaceSolver.h"

#ifdef GUI
#include "src/instrument.h"
#include "engine/hcOrbit.h"
#endif

typedef int _Atomic_word;

bool convertMagMap2bmp(const char *filename);

void drawEcliptic(hcImageRGBA &img, hcFloat maxSinLat, bool sinLatFormat, uint crNum);

class SphericalGrid;
class MagMapping;
struct threadParameterMagMapping
{
    uint threadID;
    pthread_mutex_t *runningMutex;
    SphericalGrid *grid;
    volatile _Atomic_word *numRunningThreads;
    volatile _Atomic_word *threadRunning;

    Magline *magline;			/*!< \brief pointer to magline to be worked upon														*/
	Vec3D *posStart;			/*!< \brief pointer to position from which to start tracking field line, in computational coordinates	*/

    threadParameterMagMapping(){}

    void init(	uint threadID, volatile _Atomic_word *numRunningThreads, pthread_mutex_t * runningMutex,
        		volatile _Atomic_word *threadRunning, SphericalGrid *grid)
	{
		this->threadID      	= threadID;
		this->runningMutex		= runningMutex;
		this->grid     			= grid;
		this->numRunningThreads	= numRunningThreads;
		this->threadRunning 	= threadRunning;
	}

    void set(Magline *magline, Vec3D *posStart)
    {
    	this->magline	= magline;
        this->posStart	= posStart;
    }
};
typedef struct threadParameterMagMapping threadParamMag;

class MagMapping{
public:

	bool sinLatFormat;	/*!< \brief azimuthal spacing linear in latitude (false) or sine(latitude) (true) 	*/
	bool compCoords;	/*!< \brief height level in computationalCoords										*/
    uint numTheta;		/*!< \brief number of "pixels" in meridional direction								*/
	uint numPhi;		/*!< \brief number of "pixels" in zonal direction									*/

    MagMapping();									/*!< \brief std constructor								*/
    MagMapping(const MagMapping &other);			/*!< \brief cpy constructor								*/
    MagMapping(bool sinLatFormat,
    		   bool compCoords,
    		   uint numTheta,
    		   uint numPhi,
    		   hcFloat height);						/*!< \brief constructor									*/
    virtual ~MagMapping();							/*!< \brief destructor									*/

    MagMapping &operator=(const MagMapping &other);	/*!< \brief assignment operator							*/
    Magline &operator()(uint indTheta, uint indPhi);/*!< \brief returns magline at position					*/

    void init(bool sinLatFormat,
    		  bool compCoords,
    		  uint numTheta,
    		  uint numPhi,
    		  hcFloat height);						/*!< \brief initializer									*/

    uint getNumTheta();
    uint getNumPhi();

    hcFloat getHeight() const;
        	/*!< returns base height of this map															*/

	Magline **maglines;

    bool exportBinary(const char *filename);
        /*!< \brief exports mapping to binary file                                                      	*/

    bool importBinary(const string &filename);
        /*!< \brief imports mapping from binary file                                                    	*/

    hcFloat diffFootpoints(MagMapping &other);
    	/*!< \brief computes averaged squared distance of photospheric footpoints							*/

    hcFloat getOpenFlux();
    	/*!< \brief computes the open magnetic flux through mapping											*/

    virtual bool createAtHeight(	LaplaceSolver &solver, hcFloat height, uint numThetaIn, uint numPhiIn,
                                    hcFloat maxSinLatIn, bool sinLatFormat, bool compCoords);
    	/*!< \brief create Mapping at desired height level													*/

    bool createAtHeight_SP(LaplaceSolver &solver, hcFloat height, uint numThetaIn, uint numPhiIn,
               	   	   	   hcFloat maxSinLatIn, bool sinLatFormat, bool compCoords);
    	/*!< \brief single core version of magnetic field line creation										*/

    bool createAtHeight_MP(LaplaceSolver &solver, hcFloat height, uint numThetaIn, uint numPhiIn,
                           hcFloat maxSinLatIn, bool sinLatFormat, bool compCoords);
    	/*!< \brief multi core version of magnetic field line creation										*/

    static void *createAtHeight_threadEntryPoint(void *parameter);
    	/*!< \brief thread entry point for createAtHeight_MP												*/

    void createMappingHeader(char *header, hcFloat height, hcFloat *otherHeights, uint numHeights, hcFloat lowerR, hcFloat upperR);
    	/*!< \brief gives human readable explanation of what is in exported ASCII file						*/

    bool exportImage(string filename);
    	/*!< \brief creates polarity image of this mapping													*/

    bool exportFootpointImage(string filename);
    	/*!< \brief creates polarity image of photospheric footpoints for each magline in this mapping		*/

    bool exportASCII(string filename, hcFloat *heights, uint numHeights, hcFloat lowerR, hcFloat sourceSurfaceR);
		/*!< \brief exports mag mapping with connection to several height levels to txt- and image-file 	*/

    bool exportExpansionFactorImage(SphericalGrid *grid, const string &fn_fits, const string &fn_bitmap);
        /*!< \brief export expansion factor of magnetic field lines in this map								*/

    bool exportCurvatureAtHeight(const char *filename, hcFloat height);
    	/*!< \brief exports the curvature of the magnetic field lines at height	(TODO not working)			*/

    bool exportMultiCurvature(const char *filename, hcFloat height, hcFloat *heights, uint numHeights, hcFloat lowerR, hcFloat sourceSurfaceR);
		/*!< \brief exports the curvature of the magnetic field lines at height								*/

    void dump();

#ifdef GUI
    bool createProjectedView(const string fn, Imager &imager, Imager &photImg, const hcObliquity &solarRotAxis, const hcDate &date);
#endif
private:

    Vec2D *coords;      /*!< \brief theta / phi coordinates of pixel                                    						*/
    hcFloat height;     /*!< \brief height position of map                                              						*/

    void initNULL();
    void clear();

    uint index(uint indTheta, uint indPhi);
};

#endif // MAGMAPPING_H
