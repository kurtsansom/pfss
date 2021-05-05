#include "src/magline.h"
#include "src/ellipticalGrid.h"

#include "engine/math/hcFunction.h"
#include "engine/math/hcLine.h"
#include "engine/math/hcSphere.h"

#include <stdio.h>

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Magline
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

Magline::Magline()
{
	initNULL();
}

Magline::Magline(const Magline &other)
{
	initNULL();
	init();
	if(other.numSamples > maxNumPointsPerMagline)
	{
		cerr << __FILE__ << ":" << __LINE__ << ": Other magline has more data samples that is allowed for this object (" << other.numSamples << ")\n";
		this->numSamples = maxNumPointsPerMagline;
	}
	else this->numSamples 	= other.numSamples;
	this->closed 			= other.closed;
	this->valid	 			= other.valid;
	this->polarity			= other.polarity;

	for (uint i=0; i<numSamples; ++i)
	{
		this->posdata[i] 	= other.posdata[i];
		this->magdata[i] 	= other.magdata[i];
	}

}

Magline::~Magline() {

	clear();
}

Magline &Magline::operator=(const Magline &other)
{
	if (this == &other)	return *this;
	init();
	numSamples 	= other.numSamples;
	closed 		= other.closed;
	valid 		= other.valid;
	polarity	= other.polarity;

	for (uint i = 0; i < numSamples; ++i)
	{
		this->magdata[i] 	= other.magdata[i];
		this->posdata[i] 	= other.posdata[i];
	}
	return *this;
}

void Magline::init()
{
	clear();

	magdata 	= new Vec3D[maxNumPointsPerMagline];
	posdata		= new Vec3D[maxNumPointsPerMagline];
}

void Magline::initNULL()
{
	magdata 	= NULL;
	posdata 	= NULL;
	closed 		= false;
	valid 		= false;
	polarity	= false;
	numSamples 	= 0;
}

void Magline::clear()
{
	delete[] magdata;
	delete[] posdata;

	initNULL();
}
/*! follows the magnetic field line at supplied position pos to some boundary and stores positional and magnetic
 * 	information in the arrays poslines and maglines respectively
 *
 * \param pos           Position through which the field line traverses (spherical/computational coordinates)
 * \param grid          Spherical grid which holds the magnetic field strength
 * \param error         return value stating the problem if something goes awry (see getNearestNeighbor)
 * \param stepSize      stepping size between two array entries for tracing the magnetic field line
 * \param createMagline should the array entries be created or is only the final position of interest?
 * \param poslines      array containing positional data
 * \param maglines      array containing magnetic field strength data
 * \param numSamples    number of entries created in these arrays
 * \param debug         flag for debugging output
 *
 * \return validity of magnetic field line
 *
 * error code:
 *
 * (1 << 0)		magnetic field line contains only starting point (might be valid if tracking starts at source surface)
 * (1 << 1)		magnetic field line contains too many samples
 * (1 << 2)		invalid starting position
 * (1 << 3)		loopcounter reached, maybe a circular line?
 *
 */
bool createMaglineOnTheFly(Vec3D &pos, const SphericalGrid &grid,
		unsigned char &error, bool createMagline,
		Vec3D *poslines, Vec3D *maglines, uint *numSamples, bool direction,
		bool debug, hcFloat epsMin, hcFloat epsMax)
{
	//return createMaglineOnTheFly_RKF12(pos, grid, error, createMagline, poslines, maglines, numSamples, direction, debug);
	//return createMaglineOnTheFly_RKF23(pos, grid, error, createMagline, poslines, maglines, numSamples, direction, debug);
	return createMaglineOnTheFly_RKF45(pos, grid, error, createMagline, poslines, maglines, numSamples, direction, debug, epsMin, epsMax);
	//return createMaglineOnTheFly_RK5(pos, grid, error, createMagline, poslines, maglines, numSamples, direction, debug);
	//return createMaglineOnTheFly_RK4(pos, grid, error, createMagline, poslines, maglines, numSamples, direction, debug);
}


/*! TODO. s.o.
 * 	return last, not in if-constructs?
 */
bool createMaglineOnTheFly_RKF45(Vec3D &pos, const SphericalGrid &grid,	unsigned char &error, bool createMagline,
		Vec3D *poslines, Vec3D *maglines, uint *numSamples, bool direction, bool debug, hcFloat epsMin, hcFloat epsMax)
{
	Vec3D B;
	hcFloat upperR 			= grid.upperR;
	hcFloat lowerR			= grid.lowerR;
	unsigned char errorInt 	= 0;
	error					= 0;
	bool elliptic			= grid.isElliptical();
	EllipticalGrid &egr		= *(EllipticalGrid*)&(const_cast<SphericalGrid&>(grid));
	bool valid	 			= grid.getInterpolatedB(B, pos, errorInt, debug);

	if(!valid)
	{
		if(errorInt & 240)	// non-recoverable error
		{
			cerr << __FILE__ << ":" << __LINE__ << ": Position through which to compute magline is not valid. " << "Error: " << errorInt << ", pos(r/theta/phi): " << pos[0] << "/" << pos[1] << "/" << pos[2] << "\n";
			return false;
		}

		if(errorInt & 3)	// position above/below outer/inner boundary
		{
			bool upperBoundaryViolated 	= errorInt & (1<<1);
			hcFloat dist 				= upperBoundaryViolated ? fabs(pos[0] - upperR) / upperR : fabs(pos[0] - lowerR) / lowerR;

			if(dist < 1E-9)
			{
				hcFloat r_bu 	= pos[0];
				pos.content[0]	= upperBoundaryViolated ? upperR * 0.99999 : lowerR * 1.0000001;
				valid			= grid.getInterpolatedB(B, pos, errorInt, debug); 	// Try slightly different position
				pos.content[0]	= r_bu;
			}
		}
	}

	if (!valid||errorInt)
	{
		cerr << __FILE__ << ":" << __LINE__ << ": The magline could not be followed to the boundary correctly. valid: " << valid << " , error: " << errorInt << ", pos(r/theta/phi): " << pos[0] << "/" << pos[1] << "/" << pos[2] << "\n";
		error = (1<<2);
		return false;
	}

	if( createMagline )
	{
		poslines[0] 	= pos;
		maglines[0]		= B;
	}

	float step 			= grid.getStepSize();
	uint i 				= 1;
	bool startFromLower = fabs(pos[0] - lowerR) / lowerR < 1E-5 ? true 	: false;
	bool startFromUpper	= fabs(pos[0] - upperR) / upperR < 1E-5	? true 	: false;
	int sign			=    (createMagline && ((direction && !(B[0] > 0.0)) || (!direction && (B[0] > 0.0))))
						  || (!createMagline && !(B[0] > 0.0)) ? -1 : 1;

	if((startFromLower && !direction) || (startFromUpper && direction))
	{
		*numSamples = 1;
		return true;
	}

	Vec3D pos_old 		= pos;
	Vec3D pos_old_cart 	= elliptic ? pos.convCoordEll2Cart(egr) : pos.convCoordSpher2Cart();
	bool validNextPoint = maglineTracingHandler(B, pos, pos_old_cart, pos_old, errorInt, grid);
	pos_old				= pos;
	hcFloat hMin		= step / 10;		// the smaller this value, the less invalid maglines will be produced as the interpolation
	hcFloat h			= hMin;				// will be more accurate at the equatorial region where the magnetic polarity changes
	uint counter 		= 0;
	uint maxLoop		= 10*maxNumPointsPerMagline;

	while (++counter<=maxLoop)
	{
		validNextPoint 	= maglineTracingHandler(B, pos, pos_old_cart, pos_old, errorInt, grid);
		pos_old			= pos;

		while(validNextPoint)
		{
			Vec3D pos_cart, k1, k2, k3, k4, k5, k6;
			// --------------------------------------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------------------  RKF45 - step 1
			// --------------------------------------------------------------------------------------------------------------------
			pos 			= pos_old;
			validNextPoint	= maglineTracingHandler(B, pos, pos_cart, pos_old, errorInt, grid);
			pos_old_cart	= elliptic ? pos.convCoordEll2Cart(egr) : pos.convCoordSpher2Cart();
			k1 				= sign*(elliptic ? B.convVecEll2Cart(pos_cart, egr) : B.convVecSpher2Cart(pos_cart))*h/B.length();
			if(!validNextPoint)	break;
			// --------------------------------------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------------------  RKF45 - step 2
			// --------------------------------------------------------------------------------------------------------------------
			pos 			= pos_old_cart + 1.0 / 4.0 * k1;
			pos 			= elliptic ? pos.convCoordCart2Ell(egr) : pos.convCoordCart2Spher();
			validNextPoint	= maglineTracingHandler(B, pos, pos_cart, pos_old, errorInt, grid);
			k2 				= sign*(elliptic ? B.convVecEll2Cart(pos_cart, egr) : B.convVecSpher2Cart(pos_cart))*h/B.length();
			if(!validNextPoint)	break;
			// --------------------------------------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------------------  RKF45 - step 3
			// --------------------------------------------------------------------------------------------------------------------
			pos 			= pos_old_cart + 3.0/32.0 * k1 + 9.0/32.0 * k2;
			pos 			= elliptic ? pos.convCoordCart2Ell(egr) : pos.convCoordCart2Spher();
			validNextPoint	= maglineTracingHandler(B, pos, pos_cart, pos_old, errorInt, grid);
			k3 				= sign*(elliptic ? B.convVecEll2Cart(pos_cart, egr) : B.convVecSpher2Cart(pos_cart))*h/B.length();
			if(!validNextPoint)	break;
			// --------------------------------------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------------------  RKF45 - step 4
			// --------------------------------------------------------------------------------------------------------------------
			pos 			= pos_old_cart + 1932.0/2197.0 * k1 - 7200.0/2197.0 * k2 + 7296.0/2197.0 * k3;
			pos 			= elliptic ? pos.convCoordCart2Ell(egr) : pos.convCoordCart2Spher();
			validNextPoint	= maglineTracingHandler(B, pos, pos_cart, pos_old, errorInt, grid);
			k4 				= sign*(elliptic ? B.convVecEll2Cart(pos_cart, egr) : B.convVecSpher2Cart(pos_cart))*h/B.length();
			if(!validNextPoint)	break;
			// --------------------------------------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------------------  RKF45 - step 5
			// --------------------------------------------------------------------------------------------------------------------
			pos 			= pos_old_cart + 439.0/216 * k1 - 8.0 * k2 + 3680.0/513.0 * k3 - 845.0/4104.0 * k4;
			pos 			= elliptic ? pos.convCoordCart2Ell(egr) : pos.convCoordCart2Spher();
			validNextPoint	= maglineTracingHandler(B, pos, pos_cart, pos_old, errorInt, grid);
			k5 				= sign*(elliptic ? B.convVecEll2Cart(pos_cart, egr) : B.convVecSpher2Cart(pos_cart))*h/B.length();
			if(!validNextPoint)	break;
			// --------------------------------------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------------------  RKF45 - step 6
			// --------------------------------------------------------------------------------------------------------------------
			pos 			= pos_old_cart - 8.0/27.0 * k1 + 2.0 * k2 - 3544.0/2565.0 * k3 + 1859.0/4104.0 * k4 - 11.0/40.0 * k5;
			pos 			= elliptic ? pos.convCoordCart2Ell(egr) : pos.convCoordCart2Spher();
			validNextPoint	= maglineTracingHandler(B, pos, pos_cart, pos_old, errorInt, grid);
			k6 				= sign*(elliptic ? B.convVecEll2Cart(pos_cart, egr) : B.convVecSpher2Cart(pos_cart))*h/B.length();
			if(!validNextPoint)	break;
			// --------------------------------------------------------------------------------------------------------------------
			// -------------------------------------------------------------------------------------------------  RKF45 - iteration
			// --------------------------------------------------------------------------------------------------------------------
			pos 			= pos_old_cart;
			Vec3D dirY 		= 16.0/135.0 * k1 + 6656.0/12825.0 * k3 + 28561.0/56430.0 * k4 - 9.0/50.0 * k5 + 2.0/55.0 * k6;
			Vec3D dirZ		= 25.0/216.0 * k1 + 1408.0/2565.0  * k3 + 2197.0/4104.0   * k4 - 1.0/5.0  * k5;
			Vec3D posZ		= pos_old_cart;
			posZ			+= dirZ;
			pos 			+= dirY;
			Vec3D distVec	= posZ - pos;
			hcFloat dist	= distVec.length();

			if(dist < epsMax || h < hMin)							// accuracy demand met
			{
				pos 			= elliptic ? pos.convCoordCart2Ell(egr) : pos.convCoordCart2Spher();
				validNextPoint	= maglineTracingHandler(B, pos, pos_cart, pos_old, errorInt, grid);

				if((dist < epsMin || h < hMin) && validNextPoint)	// accuracy is so high that we can increase stepsize
					h = h * 2.0;

				break;
			}
			else
				h = h / 2.0;
		}

		if(!validNextPoint && h > hMin)	// next point lies outside of grid, try reducing step size and recomputing
		{
			h 	= h/2.0;
			pos = pos_old;
			continue;
		}

		if (createMagline)
		{
			if (i < maxNumPointsPerMagline)
			{
				maglines[i] = B;
				poslines[i] = pos;
				++i;
				*numSamples = i;
			}
			else
			{
				Vec3D posC	= elliptic ? pos.convCoordEll2Cart(egr) : pos.convCoordSpher2Cart();
				error 		+= (1 << 1);
				*numSamples = maxNumPointsPerMagline;
				valid 		= false;
				break;
			}
		}

		if(!validNextPoint)
		{
			valid = true;
			break;
		}
	}

	if(i == 1)
	{
		printErrMess(__FILE__, __LINE__, "magnetic field line empty");
		error += (1 << 0);
	}

	if(counter >= maxLoop)
	{
		printErrMess(__FILE__, __LINE__, "loopcounter reached, circular magline?");
		error += (1 << 3);
		valid = false;
	}

	return valid;
}

/*!
 *  \param B 		(in/out) 	magnetic vector as returned by getInterpolatedB (spherical/computational coordinate system)
 *  \param pos		(in/out) 	position as returned by getInterpolatedB	    (spherical/computational coordinate system)
 *  \param pos_cart	(out)		position										(cartesian/physical coordinate system)
 *  \param pos_old  (in)		last valid position of magline					(spherical/computational coordinate system)
 *  \param error	(in/out) 	error returned by getInterpolatedB
 *
 *  \return true -> everything is fine, nothing to be done, false -> magnetic field requested outside of grid
 */
bool maglineTracingHandler(Vec3D &B, Vec3D &pos, Vec3D &pos_cart, Vec3D pos_old,
						unsigned char &error, const SphericalGrid &grid)
{
	bool valid			= grid.getInterpolatedB(B, pos, error);
	bool elliptic		= grid.isElliptical();
	EllipticalGrid &egr	= *(EllipticalGrid*)&(const_cast<SphericalGrid&>(grid));
	pos_cart			= elliptic ? pos.convCoordEll2Cart(egr) : pos.convCoordSpher2Cart();	// cartesian/physical coordinates TODO s.u.

	if (error)
	{
		if (error & 240) //  pos is outside of grid, not just below lower / above upper boundary
		{
			cerr << __FILE__ << ":" << __LINE__ << ": Magfield requested at invalid position. " << "valid: " << valid << ", error: " << error << ", pos(r/theta/phi): " << pos[0] << "/" << pos[1] << "/" << pos[2] << "\n";
			return false;
		}

		if(error & 3)	// if the instruction flow reaches this it means that the magline has reached some boundary
		{				// compute intersection of magline with boundary sphere and return this position as last
						// point for magline.
			pos_cart				= pos.convCoordSpher2Cart();	// now it is in cartesian/computational coordinates TODO s.o.
			bool hitUpperBoundary 	= error & (1 << 1);
			pos_old 				= pos_old.convCoordSpher2Cart();
			hcLine3D line(pos_old, pos_cart);	// line connecting the last two positions in magline
			hcSphere<3> sphere;
			sphere.init(Vec3D(0.0, 0.0, 0.0), (hitUpperBoundary ? grid.upperR : grid.lowerR));

			Vec3D res0, res1;					// find the intersection of magline with boundary sphere
			if (line.getIntersectionsWithSphere(sphere, res0, res1, hitUpperBoundary) != 2)
			{
				printErrMess(__FILE__, __LINE__, "magline intersects photosphere, though no intersection point could be found (might be low numerical accuracy, ignore this message if it is rare)");
				return false;
			}
			// get the nearest of two (line-sphere) intersection points
			pos 	= pos_old.dist(res0) < pos_old.dist(res1) ? res0 : res1;
			pos 	= pos.convCoordCart2Spher();
			valid 	= grid.getInterpolatedB(B, pos, error);

			if (!valid || error)
			{
				printErrMess(__FILE__, __LINE__, "Last step searching for inside point failed, valid: " + to_string(valid) + ", error: " + to_string(error) + ", upperR: " + toStr(grid.upperR) + ", lowerR: " + toStr(grid.lowerR) + ", pos(r/theta/phi): " + toStr(pos[0]) + "/" + toStr(pos[1]) + "/" + toStr(pos[2]));
				return false;
			}
			return false;
		}
	}

	return true;
}

/*!	@param height heliocentric height where this magline is to be sampled
 * 	@param posData position of sampled magline (spherical coordinates)
 * 	@param magData magnetic direction vector at height (spherical coordinates)
 */
bool Magline::getValuesAtHeight(hcFloat height, Vec3D &posData, Vec3D &magData)
{
	hcFloat eps = 1E-2;	// TODO: numerical accuracy per DEFINE ?

	if (height < this->posdata[0][0])
	{
		if (fabs(height - this->posdata[0][0]) / this->posdata[0][0] > eps)
		{
			cerr << __FILE__ << ":" << __LINE__ << ": requested height(" << height << " is below lowest entry in Magline(" << posdata[0][0] << "). relError=" << fabs(height - this->posdata[0][0]) / this->posdata[0][0] << " \n";
			return false;
		}
		else
		{
			posData = this->posdata[0];
			magData = this->magdata[0];
			return true;
		}
	}

	for (uint i = 0; i < numSamples; ++i)
	{
		if (this->posdata[i][0] == height)
		{
			hcFloat r 		= this->posdata[i][0];
			hcFloat theta 	= this->posdata[i][1];
			hcFloat phi 	= this->posdata[i][2];

			hcFloat br 		= this->magdata[i][0];
			hcFloat btheta 	= this->magdata[i][1];
			hcFloat bphi 	= this->magdata[i][2];

			posData			= Vec3D(r, theta, phi);
			magData			= Vec3D(br, btheta, bphi);

			return true;
		}

		if (this->posdata[i][0] > height)
		{
			hcFloat frac 	= (this->posdata[i][0] - height)	/ (this->posdata[i][0] - this->posdata[i - 1][0]);
			hcFloat phi_i 	= this->posdata[i][2];
			hcFloat phi_im 	= this->posdata[i - 1][2];

			if (fabs(phi_i - phi_im) > PI)
			{
				if (phi_i > phi_im)	phi_i  -= 2 * PI;
				else				phi_im -= 2 * PI;
			}

			hcFloat r 		= frac * this->posdata[i - 1][0] 	+ (1 - frac) * this->posdata[i][0];
			hcFloat theta 	= frac * this->posdata[i - 1][1] 	+ (1 - frac) * this->posdata[i][1];
			hcFloat phi 	= frac * phi_im 					+ (1 - frac) * phi_i;

			if (phi < 0)phi += 2 * PI;

			hcFloat br 		= frac * this->magdata[i - 1][0] + (1 - frac) * this->magdata[i][0];
			hcFloat btheta 	= frac * this->magdata[i - 1][1] + (1 - frac) * this->magdata[i][1];
			hcFloat bphi 	= frac * this->magdata[i - 1][2] + (1 - frac) * this->magdata[i][2];

			posData			= Vec3D(r, theta, phi);
			magData			= Vec3D(br, btheta, bphi);

			return true;
		}
	}

	return false;
}

/*!	@param height heliocentric height where this magline is to be sampled
 * 	@param posData positions of sampled magline (spherical coordinates)
 * 	@param magData magnetic direction vectors at height (spherical coordinates)
 * 	@retval number of positions where magline pierces through height
 * 			if more than MAGLINE_NUM_POSITIONS positions exist, only the first
 * 			MAGLINE_NUM_POSITIONS will be reported
 *
 * 	posData and magData are pre-allocated arrays of size MAGLINE_NUM_POSITIONS
 *
 */
int Magline::getAllValuesAtHeight(hcFloat height, Vec3D *posData, Vec3D *magData)
{
	hcFloat eps = 1E-2;
	int retval	= 0;

	if (height < this->posdata[0][0])
	{
		if (fabs(height - this->posdata[0][0]) / this->posdata[0][0] > eps)
		{
			printf("ERROR! Magline::getAllValuesAtHeight(): requested height is below lowest entry in Magline!\n");
			printf("Height: %E, posdata[0][0]: %E\n", height, posdata[0][0]);
			printf("relError: %E\n\n", fabs(height - this->posdata[0][0]) / this->posdata[0][0]);
			return 0;
		}
		else
		{
			posData[0] = this->posdata[0];
			magData[0] = this->magdata[0];
			++retval;
		}
	}

	if (this->posdata[0][0] == height)
	{
		hcFloat r 		= this->posdata[0][0];
		hcFloat theta 	= this->posdata[0][1];
		hcFloat phi 	= this->posdata[0][2];

		hcFloat br 		= this->magdata[0][0];
		hcFloat btheta 	= this->magdata[0][1];
		hcFloat bphi 	= this->magdata[0][2];

		posData[0]		= Vec3D(r, theta, phi);
		magData[0]		= Vec3D(br, btheta, bphi);

		++retval;
	}

	for (uint i = 1; i < numSamples; ++i)
	{
		Vec3D &pos 		= this->posdata[i];
		Vec3D &pos_m 	= this->posdata[i-1];
		Vec3D &mag 		= this->magdata[i];
		Vec3D &mag_m	= this->magdata[i-1];

		bool pierced = ((pos_m[0] < height && pos[0] >= height) || (pos_m[0] > height && pos[0] <= height));
		if (pierced)
		{
			if(retval >= MAGLINE_NUM_POSITIONS)
				break;

			hcFloat frac 	= (pos[0] - height)	/ (pos[0] - pos_m[0]);
			hcFloat phi_i 	= pos[2];
			hcFloat phi_im 	= pos_m[2];

			if (fabs(phi_i - phi_im) > PI)
			{
				if (phi_i > phi_im)
					phi_i -= 2 * PI;
				else
					phi_im -= 2 * PI;
			}

			hcFloat r 		= frac * pos_m[0] 	+ (1 - frac) * pos[0];
			hcFloat theta 	= frac * pos_m[1] 	+ (1 - frac) * pos[1];
			hcFloat phi 	= frac * phi_im 					+ (1 - frac) * phi_i;

			if (phi < 0)
				phi += 2 * PI;

			hcFloat br 		= frac * mag_m[0] + (1 - frac) * mag[0];
			hcFloat btheta 	= frac * mag_m[1] + (1 - frac) * mag[1];
			hcFloat bphi 	= frac * mag_m[2] + (1 - frac) * mag[2];

			posData[retval]	= Vec3D(r, theta, phi);
			magData[retval]	= Vec3D(br, btheta, bphi);

			++retval;
		}
	}

	return retval;
}

/*! produces the magnetic and positional values of a magnetic field line.
 *
 *  @param grid utilized for tracking field line
 *  @param posParam position through which to draw the magline (spherical/elliptical coordinates)
 *  @param debug gives debugging information throughout computation
 *  @param stepSize distance between two data samples in the visualization
 *
 *  \return error value (see below)
 *
 *  (1 << 0) = 1    maglineDown empty
 *  (1 << 1) = 2    maglineDown has too many samples
 *
 *  (1 << 2) = 4	maglineUp empty
 *  (1 << 3) = 8	maglineUp has too many samples
 *
 *  (1 << 4) = 16	both directions have no samples (except first point)
 *  (1 << 5) = 32	both directions combined have too many samples
 *
 *  (1 << 6) = 64   magline connects two points on the source surface
 *  (1 << 7) = 128	start position not valid
 *
 */
unsigned char Magline::createMaglineThroughPos(	SphericalGrid &grid, const Vec3D &posParam,
												bool debug, hcFloat epsMin, hcFloat epsMax)
{
	init();
	hcFloat eps 		= 1E-5;

	if (!posParam.isValid())
	{
		cerr << __FILE__ << "/" << __LINE__ << " posParam is not valid: ";
		cerr << posParam[0] << "/" << posParam[1] << "/" << posParam[2] << "\n";
		return (1 << 7);
	}

	Vec3D pos = posParam;

	if (pos[0] < grid.lowerR)
	{
		if (fabs(pos[0] - grid.lowerR) / grid.lowerR < eps)	pos(0) = grid.lowerR;
		else
		{
			printf("ERROR! Magline::createMaglineThroughPos: Creation of magline outside of grid requested (below lower boundary)!\npos: %E, Lower boundary: %E, r-lowerR: %E\n", pos[0], grid.lowerR, pos[0] - grid.lowerR);
			return (1 << 7);
		}
	}

	if (pos[0] > grid.upperR)
	{
		if (fabs(pos[0] - grid.upperR) / grid.upperR < eps)	pos(0) = grid.upperR;
		else
		{
			printf("ERROR! Magline::createMaglineThroughPos: Creation of magline outside of grid requested (above upper boundary)!\npos: %E, Upper boundary: %E, r-upperR: %E\n", pos[0], grid.upperR, pos[0] - grid.upperR);
			return (1 << 7);
		}
	}

	bool validator 			= true;
	Vec3D startPos_down 	= pos;
	Vec3D startPos_up 		= pos;
	Vec3D *positions_down 	= new Vec3D[maxNumPointsPerMagline];
	Vec3D *mags_down		= new Vec3D[maxNumPointsPerMagline];
	Vec3D *positions_up		= new Vec3D[maxNumPointsPerMagline];
	Vec3D *mags_up			= new Vec3D[maxNumPointsPerMagline];
	uint num_down			= 0;
	uint num_up				= 0;
	unsigned char retval	= 0;
	unsigned char error 	= 0;

	if(debug)
		cout << "createMagline through pos " << pos[0] << "/" << pos[1] << "/" << pos[2] << "\n";

	validator &= createMaglineOnTheFly(startPos_down, grid, error, true, positions_down, mags_down, &num_down, false, false, epsMin, epsMax);
	if (error)
	{
		retval += ((error & (1 << 0) ? 1 : 0) << 0);
		retval += ((error & (1 << 1) ? 1 : 0) << 1);

		if(debug)	cout << "Magline::createMaglineThroughPos: Tracking fieldline DOWN yielded error " << error << "\n";
	}

	error = 0;
	validator &= createMaglineOnTheFly(startPos_up,   grid, error, true, positions_up,   mags_up,   &num_up,   true,  false, epsMin, epsMax);
	if (error)
	{
		retval += ((error & (1 << 0) ? 1 : 0) << 2);
		retval += ((error & (1 << 1) ? 1 : 0) << 3);

		if(debug)	cout << "Magline::createMaglineThroughPos: Tracking fieldline UP yielded error " << error << "\n";
	}

	if((int)num_down + (int)num_up - 1 > maxNumPointsPerMagline )	{
		printf("ERROR! Magline::createMaglineThroughPos: combined number of samples (%u) too high!\n", num_down+num_up );
		retval 		+= (1 << 5);
		this->valid = false;
	}
	else if((int)num_down + (int)num_up - 1 < 0)
	{
		printf("ERROR! Magline::createMaglineThroughPos: magline empty!\n" );
		retval		+= (1 << 4);
		exit(1);
		this->valid = false;
	}
	else
	{
		for(uint i=0; i<num_down;++i)
		{
			this->posdata[i] 	= positions_down[num_down-1-i];
			this->magdata[i]	= mags_down[num_down-1-i];
		}
		for(uint i=1; i<num_up; ++i)
		{
			this->posdata[num_down+i-1]	= positions_up[i];
			this->magdata[num_down+i-1]	= mags_up[i];
		}
	}

	delete [] positions_down;
	delete [] positions_up;
	delete [] mags_down;
	delete [] mags_up;

	numSamples = (num_down + num_up) < 1 ? 1 : num_down + num_up - 1;

	if(!validator)
	{
		this->valid = false;
		return retval;
	}

	// sort entries of magline from lowest to highest entry
	if(posdata[0][0] > posdata[numSamples-1][0])
	{
		Vec3D *temp_posdata;
		Vec3D *temp_magdata;

		temp_posdata = posdata;
		temp_magdata = magdata;

		posdata = new Vec3D[numSamples];
		magdata = new Vec3D[numSamples];

		for(uint i=0; i<numSamples; ++i)
		{
			posdata[i] = temp_posdata[numSamples-1-i];
			magdata[i] = temp_magdata[numSamples-1-i];
		}

		delete [] temp_posdata;
		delete [] temp_magdata;
	}

	this->valid 		= true;
	bool startFromLower = (fabs(posdata[0][0] 			 		- grid.lowerR) / grid.lowerR < eps ? 1 : 0);
	bool endAtLower 	= (fabs(posdata[this->numSamples-1][0] 	- grid.lowerR) / grid.lowerR < eps ? 1 : 0);
	closed				= (startFromLower ? (endAtLower ? true : false) : (endAtLower ? false : true));

	if(!startFromLower && !endAtLower)
	{
		this->valid	= false;
		error 		+= (1<<6);
		return error;
	}

	polarity = magdata[0][0] > 0.0 ? true : false;
	return error;		// here data is still stored in computational coordinates
}

Vec3D Magline::getMagVec(uint num)
{
	if (num < numSamples) return magdata[num];

	printf("ERROR! Magline::getMagVec: requested index (%u) out of bounds (numSamples: %u)!\n",	num, numSamples);
	return Vec3D(0.0, 0.0, 0.0);
}

Vec3D Magline::getPosVec(uint num)
{
	if (num < numSamples) return posdata[num];

	printf("ERROR! Magline::getPosVec: requested index (%u) out of bounds (numSamples: %u)!\n",	num, numSamples);
	return Vec3D(0.0, 0.0, 0.0);
}

/*! lower distance between maglines is less than (LT) the upper distance of the same maglines
 *
 *  \param other Magline to be compared to
 *
 *  \return true, if distance at the lower boundary is less than distance at upper boundary, false otherwise
 *
 */
bool Magline::lowerDistLTupperDist(const Magline &other)
{
	float eps = 1E-5;

	if (!this->valid || !other.valid) {
		printf("ERROR! Magline::lowerDistLTupperDist: At least one magline is invalid! (this: %u, other: %u)\n",
				this->valid, other.valid);
		return false;
	}

	if (this->closed || other.closed) {
		printf(
				"ERROR! Magline::lowerDistLTupperDist: At least one magline is closed! (this: %u, other: %u)\n",
				this->closed, other.closed);
		return false;
	}

	if (   fabs(posdata[0][0] - other.posdata[0][0]) / posdata[0][0] > eps
		|| fabs(posdata[numSamples - 1][0]- other.posdata[other.numSamples - 1][0])	/ posdata[numSamples - 1][0] > eps)
	{
		printf("ERROR! Magline::lowerDistLTupperDist: Maglines not comparable!\n");
		printf("this->posdata[0]:\n");
		posdata[0].dump();
		printf("other->posdata[0]:\n");
		other.posdata[0].dump();
		printf("this->posdata[numSamples-1]:\n");
		posdata[numSamples - 1].dump();
		printf("other->posdata[other.numSamples-1]:\n");
		other.posdata[other.numSamples - 1].dump();
		return false;
	}

	if (distOnUnitSphere(posdata[0], other.posdata[0])
			< distOnUnitSphere(posdata[numSamples - 1],
					other.posdata[other.numSamples - 1]))
		return true;
	else
		return false;
}

bool Magline::isInSameFluxtubeAs(const Magline &other) {

	return lowerDistLTupperDist(other);
}


bool Magline::exportBinary(std::ofstream &stream)
{
	bool retval 		= true;
	uint sizeofFloat	= sizeof(hcFloat);
	stream.write(reinterpret_cast<char*>(&sizeofFloat),	sizeof(uint));
	stream.write(reinterpret_cast<char*>(&numSamples),	sizeof(uint));
	stream.write(reinterpret_cast<char*>(&valid),		sizeof(bool));
	stream.write(reinterpret_cast<char*>(&closed),		sizeof(bool));
	stream.write(reinterpret_cast<char*>(&polarity),	sizeof(bool));

	for (uint i=0; i<numSamples; ++i)
		for (uint k=0; k<3; ++k)
		{
			stream.write(reinterpret_cast<char*>(&posdata[i](k)), sizeofFloat);
			stream.write(reinterpret_cast<char*>(&magdata[i](k)), sizeofFloat);
		}

	return retval;
}

bool Magline::importBinary(std::ifstream &stream)
{
	init();

	bool retval 		= true;
	uint sizeoffloat	= 0;
	stream.read(reinterpret_cast<char*>(&sizeoffloat),	sizeof(uint));
	if(sizeoffloat != 4 && sizeoffloat != 8) return false;
	char *tempFloat		= sizeoffloat == 4 ? reinterpret_cast<char*>(new float()): reinterpret_cast<char*>(new double());
	stream.read(reinterpret_cast<char*>(&numSamples),	sizeof(uint));
	stream.read(reinterpret_cast<char*>(&valid),		sizeof(bool));
	stream.read(reinterpret_cast<char*>(&closed),		sizeof(bool));
	stream.read(reinterpret_cast<char*>(&polarity),		sizeof(bool));

	for (uint i=0; i<numSamples; ++i)
		for (uint k=0; k<3; ++k)
		{
			stream.read(reinterpret_cast<char*>(tempFloat),	sizeoffloat); posdata[i](k) = *(reinterpret_cast<hcFloat*>(tempFloat));
			stream.read(reinterpret_cast<char*>(tempFloat),	sizeoffloat); magdata[i](k) = *(reinterpret_cast<hcFloat*>(tempFloat));
		}

	delete tempFloat;
	return retval;
}

void Magline::dump() const
{
	cout << "Magline::dump():\n";
	//for (uint i = 0; i < numSamples; ++i)
	//	cout << "i: " << i << ", pos: " << posdata[i][0] << "/" << posdata[i][1] << "/" << posdata[i][2] << "\n";
	cout << "Number of samples: " << numSamples 				<< "\n";
	cout << "Magline valid:     " << valid 						<< "\n";
	cout << "Magline closed:    " << closed 					<< "\n";
	cout << "Position at lower boundary:\n";
	if (numSamples > 0)	posdata[0].dump();
	else				cout << "No data\n";
	cout << "Magnetic field at lower boundary:\n";
	if (numSamples > 0) magdata[0].dump();
	else				cout << "No data\n";

	for(uint i=0; i<numSamples; ++i)
	{
		cout << i << " " << posdata[i][0] << "/" <<  posdata[i][1] << "/" << posdata[i][2] << " " << magdata[i][0] << "/" <<  magdata[i][1] << "/" << magdata[i][2] << "\n";
	}
}

void Magline::initStaticMembers()
{
#ifdef PRESENTATION
	colorInvalid	= char2RGBA8(0,0,0,255);
#else
	colorInvalid	= char2RGBA8(255,255,255,255);
#endif
	colorClosed 	= char2RGBA8(0,255,255,255);
	colorPositive 	= char2RGBA8(0,0,255,255);
	colorNegative 	= char2RGBA8(255,0,0,255);
}
