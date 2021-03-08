#ifndef ELLIPTICALGRID_H
#define ELLIPTICALGRID_H

#include "engine/math/hcVec.h"

#include "src/grids.h"
#include "src/magMapping.h"

/*! \brief grid structure for numerical computation in 3D space
 */
class EllipticalGrid : public SphericalGrid {
public:

	bool prolate;

	virtual bool isElliptical() const;    /*!< \brief determine whether the grid is spherical or elliptical							*/

    virtual Vec3D getPos(uint index, bool ellipticCoords = false) const;

	virtual Vec3D getB(uint index, bool ellipticCoords = false) const;

	virtual void setB(uint index, Vec3D value, bool ellipticCoords = false);

#ifdef __NVCC__
    __host__ __device__
#endif
    Vec3D getSphericalPos(uint index);						/*!< \brief returns position in spherical coordinates, same as getPos, but not virtual				*/

#ifdef __NVCC__
    __host__ __device__
#endif
	hcFloat* getEllAArray() const;							/*!< \brief parameter a of ellipsis           								*/

#ifdef __NVCC__
    __host__ __device__
#endif
	hcFloat getEllA(uint index) const;

    EllipticalGrid();										/*!< \brief std constructor												*/
    EllipticalGrid(const EllipticalGrid &grid);				/*!< \brief cpy constructor												*/
    virtual ~EllipticalGrid();								/*!< \brief destructor													*/

   EllipticalGrid &operator=(const EllipticalGrid &other);	/*!< \brief assignment operator											*/
   EllipticalGrid &operator=(const SphericalGrid &other);	/*!< \brief assignment operator											*/

    void initNULL_CPU();

    void initNULL();

    void clear_CPU();

    void clear();

    void dump() const;

#ifdef __NVCC__
    __host__ __device__
#endif
    hcFloat getEllParamsFromPos(Vec3D pos, bool posElliptic) const;
    	/*!< \brief computes elliptic parameters a for arbitrary position pos													*/

    void convertMagMapping(MagMapping &map);

    virtual void init(bool sinLatGrid, hcFloat maxSinLat, hcFloat minSinLat, hcFloat lowerR,
    				hcFloat upperR, uint numR, bool clearGPU=true, hcFloat a = 1.0);

#ifdef __NVCC__
    __host__ __device__
#endif
    bool getGradientVectors(Vec3D cartPos, Vec3D &er, Vec3D &et, Vec3D &ep, bool prolate) const;
    	/*!< \brief returns cartesian gradient (contravariant) basis vectors in physical domain									*/

#ifdef __NVCC__
    __host__ __device__
#endif
    bool getTangentVectors(Vec3D cartPos, Vec3D &er, Vec3D &et, Vec3D &ep, bool prolate) const;
    	/*!< \brief returns cartesian tangent (covariant) basis vectors in physical domain										*/

#ifdef CUDA

    virtual bool allocateGPUmem();			/*!< \brief allocate memory on device															*/

	virtual void initGPUmemStruct();		/*!< \brief initialize device pointers in data structure										*/

	virtual bool pushToGPU();				/*!< \brief allocates memory on device and copies arrays from host                              */

	virtual bool pullFromGPU();		        /*!< \brief copies data from device to host                                                     */

	virtual void extract_relError();		/*!< \brief extracts rel error array from device                                                */
#endif

protected:
    hcFloat *a;						/*!< \brief parameter a of ellipsis																*/

};

#endif
