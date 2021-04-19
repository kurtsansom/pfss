#ifndef GRIDS_H
#define GRIDS_H

#include "engine/math/hcVec.h"
#include "engine/math/hcMatrix.h"
#include "engine/hcTools.h"

#include "src/pfss.h"

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Spherical Grid - declaration
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

/*! \brief grid structure for numerical computation in 3D space
 */
class SphericalGrid {
public:

	uint sizeofFloat;				/*!< \brief 4-single, 8-double, 16-long double(not supported in CUDA yet)						*/

    uint numR;              		/*!< \brief number of r-steps                                                  					*/
    uint numTheta;          		/*!< \brief number of theta-steps                                              					*/
    uint numPhi;            		/*!< \brief number of phi-steps                                                					*/

    uint numGridPoints;             /*!< \brief overall grid points                                                         		*/

    bool sinLatGrid;                /*!< \brief is the grid equally spaced in latitude (=false) or sin(lat) (=true)         		*/
    hcFloat maxSinLat;				/*!< \brief sin(latitude) of northernmost pixel border											*/
    hcFloat minSinLat;				/*!< \brief sin(latitude) of southernmost pixel border											*/

    hcFloat lowerR;					/*!< \brief lowest radial coordinate (photosphere)												*/
    hcFloat upperR;					/*!< \brief uppermost radial coordinate (source surface)										*/
    hcFloat geometricFactor;		/*!< \brief factor by which radial spacing increases from shell to shell						*/

    SphericalGrid();										/*!< \brief std constructor												*/
    SphericalGrid(const SphericalGrid &grid);				/*!< \brief cpy constructor												*/
    virtual ~SphericalGrid();								/*!< \brief destructor													*/

    virtual SphericalGrid &operator=(const SphericalGrid &other);	/*!< \brief assignment operator									*/

    virtual void initNULL();
    virtual void initNULL_CPU();

    virtual void clear();
    virtual void clear_CPU();

    virtual void init(bool sinLatGrid, hcFloat maxSinLat, hcFloat minSinLat,
            		  hcFloat lowerR, hcFloat upperR, uint numR, bool clearGPU=true, hcFloat a = 1.0);
    	/*!< \brief initializes the grid																							*/

    static void getOptimalGridParameters(uint numR, hcFloat lowerR, hcFloat upperR, hcFloat &geometricFactor, uint &numT, uint &numP, bool high=false);
    	/*!< \brief computes parameters so that distances in all directions between grid points are more or less equal				*/

    void initMembers(uint numR, uint numT, uint numP,
    		bool sinLatGrid, hcFloat maxSinLat, hcFloat minSinLat, hcFloat lowerR, hcFloat upperR);
    	/*!< \brief initializes member variables, usable by elliptical and spherical grid											*/

    void getScalingFactors();
    	/*!< \brief computes scaling factors to be used by the Laplace solver														*/

    void clearValues();
    	/*!< \brief clears the values (B_r, relError, temp, etc.) not the positions (pos)											*/

#ifdef CUDA

    SphericalGrid *d_memStruct;

    bool onDevice;					/*!< \brief data pushed to device?																*/

    void initNULL_GPU();			/*!< \brief initialize device variables to well defined state									*/

    virtual bool allocateGPUmem();	/*!< \brief allocate memory on device															*/

    virtual void initGPUmemStruct();/*!< \brief initialize device pointers in data structure										*/

    void clear_GPU();				/*!< \brief free device memory																	*/

    virtual bool pushToGPU();		/*!< \brief allocates memory on device and copies arrays from host                              */

    virtual bool pullFromGPU();		/*!< \brief copies data from device to host                                                     */

    virtual void extract_relError();/*!< \brief extracts rel error array from device                                                */

#endif

    void diff(const SphericalGrid &other);	/*!< \brief computes difference in psi array											*/

    hcFloat maxSpacingR(uint &index);		/*!< \brief returns index of grid point where spacing in r-direction is tallest			*/

    hcFloat maxSpacingTheta(uint &index);	/*!< \brief returns index of grid point where spacing in theta-direction is tallest		*/

    hcFloat maxSpacingPhi(uint &index);		/*!< \brief returns index of grid point where spacing in phi-direction is tallest		*/

    virtual bool isElliptical() const;    	/*!< \brief determine whether the grid is spherical or elliptical						*/

    virtual Vec3D getPos(uint index, bool ellipticCoords) const;

	virtual Vec3D getB(uint index, bool ellipticCoords = false) const;

	virtual void setB(uint index, Vec3D value, bool ellipticCoords = false);

	virtual void dumpCoords(uint fixed = 0, bool ellipticCoords = false) const;

	void iterateElliptic_gridPoint(uint i, uint j, uint k);
		/*!< \brief solve Laplace equation on a single grid point according to the elliptic computation scheme						*/

#ifdef __NVCC__
    __host__ __device__
#endif
	hcFloat getPsi(uint index) const;

#ifdef __NVCC__
    __host__ __device__
#endif
	void setPsi(uint index, hcFloat value);

#ifdef __NVCC__
    __host__ __device__
#endif
	void setRelError(uint index, hcFloat value);

#ifdef __NVCC__
    __host__ __device__
#endif
	hcFloat getTemp(uint index) const;

#ifdef __NVCC__
    __host__ __device__
#endif
	void setTemp(uint index, hcFloat value);

    hcFloat getRelError(uint index) const;

#ifdef __NVCC__
    __host__ __device__
#endif
	hcFloat* getPsiArray() const;
    	/*!< \brief scalar potential      																							*/

	hcFloat* getRelErrorArray() const;
		/*!< \brief relative error compered to preciding step           															*/

#ifdef __NVCC__
    __host__ __device__
#endif
	Vec3D* getPosArray() const;
    	/*!< \brief position in world coordinates (spherical)           															*/

    Vec3D* getBArray() const;
    	/*!< \brief magnetic field components (spherical)           																*/

    hcFloat* getTempArray() const;
    	/*!< \brief temp value for computations on spec. grid point           														*/

    PFSSsolution_SHC *harmSolution;
    	/*!< \brief PFSS solution via spherical harmonic function approach, temporary testing variable								*/

#ifdef __NVCC__
    __host__ __device__ __forceinline__
#endif
    float getStepSize() const;
        /*!< \brief returns the optimal step size for following values (e.g., Magnetic field lines) throughout the grid     		*/

#ifdef __NVCC__
    __host__ __device__
#endif
    inline uint getIndex(uint i, uint j, uint k) const;
    uint getIndexPhiPlus(uint ind);

#ifdef __NVCC__
    __host__ __device__ __forceinline__
#endif
    void getIndicesFromIndex(uint ind, uint &i, uint  &j, uint &k) const;
    	/*!< \brief if ind = k * numR * numTheta + j * numR + i is given, computes i, j, k									*/

    void printGridIndicesFromIndex(uint ind);
    void computeLowerBoundaryPsi(float *image);

    virtual unsigned char getNearestNeighbors(Vec3D &pos, Vec<8, hcFloat> &indices, bool debug = false) const;

    virtual bool getInterpolatedB(Vec3D &B, Vec3D &pos, unsigned char &error, bool debug=false) const;

    Vec3D getBFromPsi(uint r, uint t, uint p);

    void compDerivR(Matrix3x3 *arr, Matrix3x3 *deriv);	/*!< \brief computes derivative w.r.t. r of quantity stored at arr			*/
    void compDerivT(Matrix3x3 *arr, Matrix3x3 *deriv);	/*!< \brief computes derivative w.r.t. theta of quantity stored at arr		*/
    void compDerivP(Matrix3x3 *arr, Matrix3x3 *deriv);	/*!< \brief computes derivative w.r.t. phi of quantity stored at arr		*/

    bool exportEntireGrid(const char *filename);		/*!< \brief export grid to binary file										*/
    bool importEntireGrid(const char *filename);		/*!< \brief import grid from binary file									*/

    bool evaluateCSSS(const char *filename);

    virtual void dump() const;

    Vec3D *pos;						/*!< \brief position in computational coordinates (spherical)      								*/
	Vec3D *B;						/*!< \brief magnetic field components (spherical)               								*/
	hcFloat *psi;                   /*!< \brief scalar potential                                    								*/

	hcFloat *relError;              /*!< \brief relative error compered to preciding step           								*/
	hcFloat *temp;                  /*!< \brief temp value for computations on spec. grid point										*/

	Matrix3x3 *g;					/*!< \brief metric coefficients																	*/
	Matrix3x3 *h;					/*!< \brief dual metric coefficients															*/
	Matrix3x3 *p;					/*!< \brief metric helper variable (sqrt(g)*g)													*/
	Matrix3x3 *p_r;					/*!< \brief derivative of p w.r.t. r															*/
	Matrix3x3 *p_t;					/*!< \brief derivative of p w.r.t. theta														*/
	Matrix3x3 *p_p;					/*!< \brief derivative of p w.r.t. phi															*/

	hcFloat *s_ijk;
	hcFloat *s_imjk;
	hcFloat *s_ipjk;

	hcFloat *s_ijmk;
	hcFloat *s_ijpk;

	hcFloat *s_ijkm;
	hcFloat *s_ijkp;

	hcFloat *s_imjmk;
	hcFloat *s_imjpk;
	hcFloat *s_ipjmk;
	hcFloat *s_ipjpk;

	hcFloat *s_imjkm;
	hcFloat *s_imjkp;
	hcFloat *s_ipjkm;
	hcFloat *s_ipjkp;

	hcFloat *s_ijmkm;
	hcFloat *s_ijmkp;
	hcFloat *s_ijpkm;
	hcFloat *s_ijpkp;

	hcFloat *s_ijmmk;
	hcFloat *s_imjmmk;
	hcFloat *s_ipjmmk;
	hcFloat *s_ijmmkm;
	hcFloat *s_ijmmkp;

	hcFloat *s_ijppk;
	hcFloat *s_imjppk;
	hcFloat *s_ipjppk;
	hcFloat *s_ijppkm;
	hcFloat *s_ijppkp;
};



//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          SphericalGrid - inline implementation
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

/*! for now assumes equidistant grid! this should probably be amended sometime
 */
inline float SphericalGrid::getStepSize() const
{
    uint ind0 = getIndex(0,0,0);
    uint ind1 = getIndex(1,0,0);
    return((pos[ind1][0]-pos[ind0][0]) / 10);
}

inline uint SphericalGrid::getIndex(uint i, uint j, uint k) const
{
	if(i >= numR || j >= numTheta || k >=numPhi)
		printf("SphericalGrid::getIndex out of bounds! i:%u/%u j: %u/%u k: %u/%u\n", i,numR,j,numTheta,k,numPhi);

    return k * numR * numTheta + j * numR + i;
}

inline void SphericalGrid::getIndicesFromIndex(uint ind, uint &i, uint  &j, uint &k) const
{
    i=ind%numR;
    uint ind1 = (ind-i)/numR;

    j = ind1%numTheta;
    k = (ind1-j)/numTheta;
}

/*! \brief computes the nearest grid cells in the grid
 *
 * 	@param position for which nearest neighbors shall be found (spherical/computational coordinates)
 * 	@param indices	vector containing grid indices for nearest neighbors
 *
 *
 *  returns bit-coded clipping id if pos is outside grid:
 *
 *  (1 << 0) = 1   	pos.r below lower boundary
 *  (1 << 1) = 2   	pos.r above upper boundary
 *
 *  (1 << 2) = 4   	pos.theta above highest latitude (north pole)	in these two cases, nearest neighbors are being computed in accordance with the
 *  (1 << 3) = 8   	pos.theta below lowest latitude (south pole) 	polar boundary conditions
 *
 *  (1 << 7) = 128  different error, should not be possible to reach
 *
 *  indices has to be 8-dim vector, numbering is done just like with the grid:
 *
 *  indices[0]  -> rm, tm, pm
 *  indices[1]  -> rp, tm, pm
 *  indices[2]  -> rm, tp, pm
 *  indices[3]  -> rp, tp, pm
 *
 *  indices[4]  -> rm, tm, pp
 *  indices[5]  -> rp, tm, pp
 *  indices[6]  -> rm, tp, pp
 *  indices[7]  -> rp, tp, pp
 *
 */
inline unsigned char SphericalGrid::getNearestNeighbors(Vec3D &position, Vec<8, hcFloat> &indices, bool debug) const
{
	const hcFloat eps = 1E-4;

    while(position[2] < 0.0) 		position(2) += 2 * PI;
    while(position(2) >= 2 * PI)	position(2) -= 2 * PI;

    if(position[1] < 0.0)
    {
        if(fabs(position[1]) > PI)	return (1 << 7);

        position(1) = -position[1];
        position(2)-= PI;

        if(position[2] < 0.0)		position(2) += 2 * PI;
    }

    if(position[1] > PI)
    {
        if(position[1] >= 2 * PI)	return (1 << 7);

        position(1) = PI - (position[1] - PI);
        position(2) -= PI;

        if(position[2] < 0)			position(2) += 2 * PI;
    }

    indices.content[0] = 0;
    indices.content[1] = 0;
    indices.content[2] = 0;
    indices.content[3] = 0;
    indices.content[4] = 0;
    indices.content[5] = 0;
    indices.content[6] = 0;
    indices.content[7] = 0;

    uint i = 0;
    uint j = 0;
    uint k = 0;

    bool elliptical = isElliptical(); // TODO: should always be true, regardless of ellipticgrid or sphericalgrid

	while(getPos(getIndex(i, 0, 0), elliptical)[0]  <= position[0] && ++i < numR);
	while(getPos(getIndex(0, j, 0), elliptical)[1]  <= position[1] && ++j < numTheta);
	while(getPos(getIndex(0, 0, k), elliptical)[2]  <= position[2] && ++k < numPhi);

    if(k==numPhi)
        k = 0;

    //*
    if(i==numR && (        position[0] <= getPos(getIndex(numR-1, 0, 0), elliptical)[0]
                   || fabs(position[0] -  getPos(getIndex(numR-1, 0, 0), elliptical)[0])/position[0] < eps))//*/
    /*/
    if(i==numR && (        position[0]/r_sol <= getPos(getIndex(numR-1, 0, 0), elliptical)[0] / r_sol
                       || fabs(position[0] / r_sol -  getPos(getIndex(numR-1, 0, 0), elliptical)[0] / r_sol)/(position[0] / r_sol) < 1E-6))//*/
    	i -= 1;

    unsigned char retval = 0;
    retval += (i==0      	? 1 : 0) << 0;
    retval += (i==numR      ? 1 : 0) << 1;
    retval += (j==0         ? 1 : 0) << 2;
    retval += (j==numTheta  ? 1 : 0) << 3;

    if(i>0 && i<numR)
	{
		uint kk 		= k==0 ? numPhi-1 : k-1;
		uint phi_mm 	= kk;
		uint phi_mp 	= j==0 || j==numTheta ? (kk + 1) 			  % numPhi	: k;
		uint phi_pm 	= j==0 || j==numTheta ? (kk + numPhi / 2)     % numPhi	: kk;
		uint phi_pp 	= j==0 || j==numTheta ? (kk + numPhi / 2 + 1) % numPhi	: k;
		uint theta_m 	= j==0 ? 0 	: (j==numTheta ? numTheta-1 : j-1);
		uint theta_p	= j==0 ? 0 	: (j==numTheta ? numTheta-1 : j  );

		indices.content[0] = getIndex(i-1, theta_m, phi_mm);
		indices.content[1] = getIndex(i  , theta_m, phi_mm);
		indices.content[2] = getIndex(i-1, theta_p, phi_pm);
		indices.content[3] = getIndex(i  , theta_p, phi_pm);
		indices.content[4] = getIndex(i-1, theta_m, phi_mp);
		indices.content[5] = getIndex(i  , theta_m, phi_mp);
		indices.content[6] = getIndex(i-1, theta_p, phi_pp);
		indices.content[7] = getIndex(i  , theta_p, phi_pp);

		//if(debug)
		if(false)
		{
			Vec3D n0 = getPos(indices[0], elliptical);
			Vec3D n1 = getPos(indices[1], elliptical);
			Vec3D n2 = getPos(indices[2], elliptical);
			Vec3D n3 = getPos(indices[3], elliptical);
			Vec3D n4 = getPos(indices[4], elliptical);
			Vec3D n5 = getPos(indices[5], elliptical);
			Vec3D n6 = getPos(indices[6], elliptical);
			Vec3D n7 = getPos(indices[7], elliptical);

			Vec3D m0 = getB(indices[0], elliptical);
			Vec3D m1 = getB(indices[1], elliptical);
			Vec3D m2 = getB(indices[2], elliptical);
			Vec3D m3 = getB(indices[3], elliptical);
			Vec3D m4 = getB(indices[4], elliptical);
			Vec3D m5 = getB(indices[5], elliptical);
			Vec3D m6 = getB(indices[6], elliptical);
			Vec3D m7 = getB(indices[7], elliptical);

			uint i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4,i5,j5,k5,i6,j6,k6,i7,j7,k7;
			getIndicesFromIndex(indices[0], i0, j0, k0);
			getIndicesFromIndex(indices[1], i1, j1, k1);
			getIndicesFromIndex(indices[2], i2, j2, k2);
			getIndicesFromIndex(indices[3], i3, j3, k3);
			getIndicesFromIndex(indices[4], i4, j4, k4);
			getIndicesFromIndex(indices[5], i5, j5, k5);
			getIndicesFromIndex(indices[6], i6, j6, k6);
			getIndicesFromIndex(indices[7], i7, j7, k7);

			printf("getNN: position: %E/%E/%E\n", position[0], position[1], position[2]);
			printf("getNN: 0: %u/%u/%u %E/%E/%E %E/%E/%E\n", i0, j0, k0, n0[0], n0[1], n0[2], m0[0], m0[1], m0[2]);
			printf("getNN: 1: %u/%u/%u %E/%E/%E %E/%E/%E\n", i1, j1, k1, n1[0], n1[1], n1[2], m1[0], m1[1], m1[2]);
			printf("getNN: 2: %u/%u/%u %E/%E/%E %E/%E/%E\n", i2, j2, k2, n2[0], n2[1], n2[2], m2[0], m2[1], m2[2]);
			printf("getNN: 3: %u/%u/%u %E/%E/%E %E/%E/%E\n", i3, j3, k3, n3[0], n3[1], n3[2], m3[0], m3[1], m3[2]);
			printf("getNN: 4: %u/%u/%u %E/%E/%E %E/%E/%E\n", i4, j4, k4, n4[0], n4[1], n4[2], m4[0], m4[1], m4[2]);
			printf("getNN: 5: %u/%u/%u %E/%E/%E %E/%E/%E\n", i5, j5, k5, n5[0], n5[1], n5[2], m5[0], m5[1], m5[2]);
			printf("getNN: 6: %u/%u/%u %E/%E/%E %E/%E/%E\n", i6, j6, k6, n6[0], n6[1], n6[2], m6[0], m6[1], m6[2]);
			printf("getNN: 7: %u/%u/%u %E/%E/%E %E/%E/%E\n", i7, j7, k7, n7[0], n7[1], n7[2], m7[0], m7[1], m7[2]);

			fflush(stdout);
		}
	}

    return retval;
}

/*! \brief  computes an interpolated value for B from the surrounding grid points via
 *          trilinear interpolation (r->theta->phi)
 *
 *  \param B 		(out) 		the magnetic field at pos (spherical coordinate system)
 *  \param position	(in/out)	the position (spherical/computational coordinates) where the magnetic field is requested
 *  \param error	(out)		error code (see below)
 *
 *  returns true only if a valid magnetic field vector has been found
 *
 *  The error codes returned in error are:
 *
 *  (1 << 0) = 1        pos.r below lower boundary
 *  (1 << 1) = 2        pos.r above upper boundary
 *
 *  (1 << 7) = 128      different error, should not be possible to reach
 *
 */
inline bool SphericalGrid::getInterpolatedB(Vec3D &B_out, Vec3D &position, unsigned char &error, bool debug) const
{
	if(debug)
	{
		stringstream sstream;
		sstream << __FILE__ << ":" << __LINE__ << " pos: " << position[0] << "/" << position[1] << "/" << position[2] << "\n";
		cout << sstream.str();
	}

    Vec<8, hcFloat> ind;
    unsigned char retval;

    if(!position.isValid() || fabs(position[1]) > 20 || fabs(position[2]) > 20)
    {
        printf("ERROR! SphericalGrid::getInterpolatedB: You supplied an invalid position!\npos: %E / %E / %E\n", position[0], position[1], position[2]);
        return false;
    }

    retval 	= getNearestNeighbors(position, ind, debug);
    error 	= retval;

    if(error & (1<<7))
	{
		printf("ERROR! SphericalGrid:getInterpolatedB: This should not happen!\nerror: %u\npos: %E / %E / %E", error, position[0], position[1], position[2]);
		return false;
	}

    if(error & ((1<<1) + (1<<0)))  	// position above/below upper/lower boundary
    {
    	if(debug) cout << "getInterpolated: pos above/below grid, return false\n";
		return false;
    }

    bool elliptical		= isElliptical();
    EllipticalGrid &egr = *(EllipticalGrid*)this;

    uint ind_000 	= ind[0];
	uint ind_100 	= ind[1];
	uint ind_010 	= ind[2];
	uint ind_110 	= ind[3];
	uint ind_001 	= ind[4];
	uint ind_101 	= ind[5];
	uint ind_011 	= ind[6];
	uint ind_111 	= ind[7];

	Vec3D pos_000	= getPos(ind_000, elliptical);
	Vec3D pos_100	= getPos(ind_100, elliptical);
	Vec3D pos_010	= getPos(ind_010, elliptical);
	Vec3D pos_110	= getPos(ind_110, elliptical);
	Vec3D pos_001	= getPos(ind_001, elliptical);
	Vec3D pos_101	= getPos(ind_101, elliptical);
	Vec3D pos_011	= getPos(ind_011, elliptical);
	Vec3D pos_111	= getPos(ind_111, elliptical);

	Vec3D B_000		= getB(ind_000, elliptical);
	Vec3D B_100		= getB(ind_100, elliptical);
	Vec3D B_010		= getB(ind_010, elliptical);
	Vec3D B_110		= getB(ind_110, elliptical);
	Vec3D B_001		= getB(ind_001, elliptical);
	Vec3D B_101		= getB(ind_101, elliptical);
	Vec3D B_011		= getB(ind_011, elliptical);
	Vec3D B_111		= getB(ind_111, elliptical);

    if(error)
    {
    	if(debug)
    	{
    		stringstream sstream;
    		sstream << __FILE__ << ":" << __LINE__ << " Error " << error << "encountered\n";
    		cout << sstream.str();
    	}
        if(error & ((1<<2) + (1<<3)))   // position lies outside of theta boundaries, this can be remedied by utilizing
        {                               // polar boundaries

			hcFloat dr    	=  pos_100[0]	- pos_000[0];
			hcFloat drm 	= (pos_100[0]	- position[0]) / dr;

            hcFloat dt 		= error & (1<<2) 				? pos_000[1]  						:  PI - pos_000[1];
            hcFloat dtm		= error & (1<<2) 				? position[1] / dt 					: (PI - position[1]) / dt;

            hcFloat dp		= pos_000[2] 	> pos_001[2] 	? 2*PI - pos_000[2]  + pos_001[2] 	: pos_001[2] - pos_000[2];
            hcFloat dpm		= (position[2] 	> pos_001[2] 	? 2*PI - position[2] + pos_001[2] 	: pos_001[2] - position[2]) / dp;

            Vec3D pos_plow((double)pos_000[0],  (error & (1<<2) ? 0.0 : PI), 0.0);
            Vec3D pos_phigh((double)pos_100[0], (error & (1<<2) ? 0.0 : PI), 0.0);

            pos_000			= pos_000.convCoordSpher2Cart();
            pos_100			= pos_100.convCoordSpher2Cart();
            pos_001			= pos_001.convCoordSpher2Cart();
            pos_101			= pos_101.convCoordSpher2Cart();

            pos_plow		= pos_plow.convCoordSpher2Cart();
            pos_phigh		= pos_phigh.convCoordSpher2Cart();

            B_000			= B_000.convVecSpher2Cart(pos_000);
            B_100			= B_100.convVecSpher2Cart(pos_100);
            B_001			= B_001.convVecSpher2Cart(pos_001);
            B_101			= B_101.convVecSpher2Cart(pos_101);

			Vec3D B_plow;
			Vec3D B_phigh;
			uint theta_ind = error & (1<<2) ? 0 : numTheta - 1;

			uint x, y, z;
			getIndicesFromIndex(ind_000, x, y, z);

            for(uint l=0;l<numPhi;++l)
            {
                uint index_low 	= getIndex(x, theta_ind, l);
                uint index_high = getIndex(x+1, theta_ind, l);
                Vec3D poslow	= getPos(index_low, elliptical);
				Vec3D poshigh	= getPos(index_high, elliptical);
                Vec3D Blow		= getB(index_low, elliptical);
                Vec3D Bhigh		= getB(index_high, elliptical);

                poslow			= poslow.convCoordSpher2Cart();
                poshigh			= poshigh.convCoordSpher2Cart();

                Blow			= Blow.convVecSpher2Cart(poslow);
                Bhigh			= Bhigh.convVecSpher2Cart(poshigh);

                B_plow  += Blow;
                B_phigh += Bhigh;
            }

            B_plow  /= numPhi;
            B_phigh /= numPhi;

            Vec3D B_mid_r0 		= B_000  	* drm + B_100    	* (1 - drm);
            Vec3D B_mid_r2 		= B_001  	* drm + B_101    	* (1 - drm);
            Vec3D B_mid_pole 	= B_plow 	* drm + B_phigh  	* (1 - drm);

            Vec3D B_mid_t0 		= B_mid_r0 	* dtm + B_mid_pole 	* (1 - dtm);
            Vec3D B_mid_t1 		= B_mid_r2 	* dtm + B_mid_pole 	* (1 - dtm);

            B_out  				= B_mid_t0 	* dpm + B_mid_t1 	* (1 - dpm);

            B_out				= B_out.convVecCart2Spher(position.convCoordSpher2Cart());

            if(harmSolution != NULL)  	harmSolution->eval(position, B_out); // TODO: ugly,away with this
            if(error & (1<<2))			error -= (1<<2);
            if(error & (1<<3))			error -= (1<<3);

            return true;
        }

        printf("SphericalGrid::getInterpolatedB: This should not be reachable!\n");
        return false;
    }

    B_000	= elliptical ? B_000.convVecEll2Cart(pos_000.convCoordEll2Cart(egr), egr) : B_000.convVecSpher2Cart(pos_000.convCoordSpher2Cart());
	B_100	= elliptical ? B_100.convVecEll2Cart(pos_100.convCoordEll2Cart(egr), egr) : B_100.convVecSpher2Cart(pos_100.convCoordSpher2Cart());
	B_010	= elliptical ? B_010.convVecEll2Cart(pos_010.convCoordEll2Cart(egr), egr) : B_010.convVecSpher2Cart(pos_010.convCoordSpher2Cart());
	B_110	= elliptical ? B_110.convVecEll2Cart(pos_110.convCoordEll2Cart(egr), egr) : B_110.convVecSpher2Cart(pos_110.convCoordSpher2Cart());
	B_001	= elliptical ? B_001.convVecEll2Cart(pos_001.convCoordEll2Cart(egr), egr) : B_001.convVecSpher2Cart(pos_001.convCoordSpher2Cart());
	B_101	= elliptical ? B_101.convVecEll2Cart(pos_101.convCoordEll2Cart(egr), egr) : B_101.convVecSpher2Cart(pos_101.convCoordSpher2Cart());
	B_011	= elliptical ? B_011.convVecEll2Cart(pos_011.convCoordEll2Cart(egr), egr) : B_011.convVecSpher2Cart(pos_011.convCoordSpher2Cart());
	B_111	= elliptical ? B_111.convVecEll2Cart(pos_111.convCoordEll2Cart(egr), egr) : B_111.convVecSpher2Cart(pos_111.convCoordSpher2Cart());

	hcFloat drm = (pos_100[0]	- position[0]) / (pos_100[0]	- pos_000[0]);
	hcFloat dtm = (pos_010[1] 	- position[1]) / (pos_010[1]	- pos_000[1]);
	hcFloat dpm	= (position[2] 	> pos_001[2] 	? 2*PI - position[2] + pos_001[2] 	: pos_001[2] - position[2])
			    / (pos_000[2] 	> pos_001[2] 	? 2*PI - pos_000[2]  + pos_001[2] 	: pos_001[2] - pos_000[2]);

    Vec3D B_mid_r0 	= B_000    * drm + B_100    * (1 - drm);
    Vec3D B_mid_r1 	= B_010    * drm + B_110    * (1 - drm);
    Vec3D B_mid_r2 	= B_001    * drm + B_101    * (1 - drm);
    Vec3D B_mid_r3 	= B_011    * drm + B_111    * (1 - drm);

    Vec3D B_mid_t0 	= B_mid_r0 * dtm + B_mid_r1 * (1 - dtm);
    Vec3D B_mid_t1 	= B_mid_r2 * dtm + B_mid_r3 * (1 - dtm);

    cout.precision(4);
	if(false)
	{
		cout << "elliptical: " << elliptical << "\n";
		cout << "drm: " << drm << " dtm: " << dtm << " dpm: " << dpm << "\n";
		cout << "position: " << position[0] << "/" << position[1] << "/" << position[2] << "\n\n";
		cout << "pos_000: " << pos_000[0] << "/" << pos_000[1] << "/" << pos_000[2] << "\n";
		cout << "pos_100: " << pos_100[0] << "/" << pos_100[1] << "/" << pos_100[2] << "\n";
		cout << "pos_010: " << pos_010[0] << "/" << pos_010[1] << "/" << pos_010[2] << "\n";
		cout << "pos_110: " << pos_110[0] << "/" << pos_110[1] << "/" << pos_110[2] << "\n";
		cout << "pos_001: " << pos_001[0] << "/" << pos_001[1] << "/" << pos_001[2] << "\n";
		cout << "pos_101: " << pos_101[0] << "/" << pos_101[1] << "/" << pos_101[2] << "\n";
		cout << "pos_011: " << pos_011[0] << "/" << pos_011[1] << "/" << pos_011[2]	<< "\n";
		cout << "pos_111: " << pos_111[0] << "/" << pos_111[1] << "/" << pos_111[2]	<< "\n\n";
		cout << "B_000: " << B_000[0] << "/" << B_000[1] << "/" << B_000[2] << "\n";
		cout << "B_100: " << B_100[0] << "/" << B_100[1] << "/" << B_100[2] << "\n";
		cout << "B_010: " << B_010[0] << "/" << B_010[1] << "/" << B_010[2] << "\n";
		cout << "B_110: " << B_110[0] << "/" << B_110[1] << "/" << B_110[2] << "\n";
		cout << "B_001: " << B_001[0] << "/" << B_001[1] << "/" << B_001[2] << "\n";
		cout << "B_101: " << B_101[0] << "/" << B_101[1] << "/" << B_101[2] << "\n";
		cout << "B_011: " << B_011[0] << "/" << B_011[1] << "/" << B_011[2]	<< "\n";
		cout << "B_111: " << B_111[0] << "/" << B_111[1] << "/" << B_111[2]	<< "\n\n";
		cout << "B_mid_r0: " << B_mid_r0[0] << "/" << B_mid_r0[1] << "/" << B_mid_r0[2] << "\n";
		cout << "B_mid_r1: " << B_mid_r1[0] << "/" << B_mid_r1[1] << "/" << B_mid_r1[2] << "\n";
		cout << "B_mid_r2: " << B_mid_r2[0] << "/" << B_mid_r2[1] << "/" << B_mid_r2[2] << "\n";
		cout << "B_mid_r3: " << B_mid_r3[0] << "/" << B_mid_r3[1] << "/" << B_mid_r3[2] << "\n";
		cout << "B_mid_t0: " << B_mid_t0[0] << "/" << B_mid_t0[1] << "/" << B_mid_t0[2] << "\n";
		cout << "B_mid_t1: " << B_mid_t1[0] << "/" << B_mid_t1[1] << "/" << B_mid_t1[2] << "\n";
		cout << "B_out:    " << B_out[0]   << "/" << B_out[1] 	<< "/" << B_out[2] 		<< "\n\n";
	}

    B_out  			= B_mid_t0 * dpm + B_mid_t1 * (1 - dpm);

    B_out			= elliptical ? B_out.convVecCart2Ell(position.convCoordEll2Cart(egr), egr) : B_out.convVecCart2Spher(position.convCoordSpher2Cart());

    if(harmSolution != NULL)  	harmSolution->eval(position, B_out); // TODO: ugly,away with this

    error = retval;
    if(debug)
    {
    	stringstream sstream;
    	sstream << __FILE__ << ":" << __LINE__ << " returning true, error: " << (uint)error << ", mag: " << B_out[0] << "/" << B_out[0] <<  "/" << B_out[0] << "\n";
    	cout << sstream.str();
    }
    return true;
}

#endif // GRIDS_H
