#ifndef MAGLINE_H
#define MAGLINE_H

#include "engine/math/hcVec.h"
#include "src/grids.h"

#define maxNumPointsPerMagline 1000
#define MAGLINE_NUM_POSITIONS 10

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Magline - declaration (Warning: already tried CUDA implementation to speed things up, doesnt work)
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

class Magline{
public:

	static uint colorInvalid;
	static uint colorClosed;
	static uint colorPositive;
	static uint colorNegative;

    uint numSamples;    /*!< \brief number of points where the magline has been sampled             				*/
    Vec3D *magdata;     /*!< \brief array of magnetic field strengths where magline has been sampled (spherical) 	*/
    Vec3D *posdata;     /*!< \brief array of positions where magline has been sampled (spherical)   				*/
    hcFloat *curvature;	/*!< \brief unsigned curvature of the magline at all sample points							*/
    char flag;			/*!< \brief 0=invalid, 1=closed, 2=positive, 3=negative, 4=man.+, 5=man.-					*/

    bool closed;        /*!< \brief does the magline connect two points on the phtotosphere?        				*/
    bool valid;         /*!< \brief are there valid values in magdata and posdata?                  				*/
    bool polarity;		/*!< \brief true if coming out of surface, false if going in, don't care if valid=false		*/

    Magline();                      							/*!< \brief std constructor         				*/
    Magline(const Magline &other);  							/*!< \brief cpy constructor         				*/
    virtual ~Magline();											/*!< \brief destructor								*/

    Magline &operator=(const Magline &other);					/*!< \brief assignment operator						*/

    void init();
    void initNULL();
    void clear();

    bool computeCurvature();

    Vec3D getMagVec(uint num);
    	/*!< \brief returns magnetic vector at given numerical position in array										*/

    Vec3D getPosVec(uint num);
    	/*!< \brief returns positional vector (spherical coordinates) at given numerical position in array				*/

    bool getValuesAtHeight(hcFloat height, Vec3D &posData, Vec3D &magData);
        /*!< \brief (linearily) interpolates values at sample points to obtain intermediate values at a specific height */

    int getAllValuesAtHeight(hcFloat height, Vec3D *posData, Vec3D *magData);
    	/*!< \brief returns number of positions where this maglines pierces through a specific height					*/

    bool getCurvatureAtHeight(hcFloat height, hcFloat &curvature);
    	/*!< \brief (linearily) interpolates curvature values at height													*/

    unsigned char createMaglineThroughPos(	SphericalGrid &grid, const Vec3D &posParam,
    										bool debug, hcFloat epsMin=1E1, hcFloat epsMax=1E2);
        /*!< \brief tracks a magnetic field line through a specific position                                            */

    bool lowerDistLTupperDist(const Magline &other);
        /*!< \brief checks if the distance between the lower endings of two lines is less than their upper distance     */
    
    bool isInSameFluxtubeAs(const Magline &other);
        /*!< \brief calls lowerDistLTupperDist, this criterion is just sufficient                                       */

    bool exportBinary(char **array);
        /*!< \brief dumps instance into *array                                                                          */

    bool importBinary(char *array);
        /*!< \brief imports instance from array                                                                         */

    static void initStaticMembers();

    void dump() const;

};

/*! \brief finds (unfortunately not the nearest) boundary intersection of some magline at pos
 *
 */

bool createMaglineOnTheFly(Vec3D &pos, const SphericalGrid &grid,
		unsigned char &error, bool createMagline,
		Vec3D *poslines, Vec3D *maglines, uint *numSamples, bool direction,
		bool debug, hcFloat epsMin, hcFloat epsMax);

bool createMaglineOnTheFly_RK4(Vec3D &pos, const SphericalGrid &grid,
		unsigned char &error, bool createMagline,
		Vec3D *poslines, Vec3D *maglines, uint *numSamples, bool direction,
		bool debug);

bool createMaglineOnTheFly_RK5(Vec3D &pos, const SphericalGrid &grid,
		unsigned char &error, bool createMagline,
		Vec3D *poslines, Vec3D *maglines, uint *numSamples, bool direction,
		bool debug);

bool createMaglineOnTheFly_RKF45(Vec3D &pos, const SphericalGrid &grid,
		unsigned char &error, bool createMagline,
		Vec3D *poslines, Vec3D *maglines, uint *numSamples, bool direction,
		bool debug, hcFloat epsMin, hcFloat epsMax);

bool createMaglineOnTheFly_RKF12(Vec3D &pos, const SphericalGrid &grid,
		unsigned char &error, bool createMagline,
		Vec3D *poslines, Vec3D *maglines, uint *numSamples, bool direction,
		bool debug);

bool createMaglineOnTheFly_RKF23(Vec3D &pos, const SphericalGrid &grid,
		unsigned char &error, bool createMagline,
		Vec3D *poslines, Vec3D *maglines, uint *numSamples, bool direction,
		bool debug);

/*! \brief requests the magnetic field at next position and handles the different errors that may occur doing so
 */
bool maglineTracingHandler(Vec3D &B, Vec3D &pos, Vec3D &pos_cart, Vec3D pos_old, unsigned char &error, const SphericalGrid &grid, bool debug=false);

#endif
