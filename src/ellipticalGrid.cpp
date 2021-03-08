#include "src/ellipticalGrid.h"

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          EllipticalGrid
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

EllipticalGrid::EllipticalGrid() : SphericalGrid()
{
    initNULL();
}

EllipticalGrid::EllipticalGrid(const EllipticalGrid &other) : SphericalGrid()
{
    initNULL();
    *this = other;
}

EllipticalGrid::~EllipticalGrid()
{
	clear();
}

EllipticalGrid &EllipticalGrid::operator=(const EllipticalGrid &other)
{
	if(this == &other)
		return *this;

	clear();
	if(other.numGridPoints == 0)
		return *this;

	numR            = other.numR;
	numTheta        = other.numTheta;
	numPhi          = other.numPhi;

	a				= new hcFloat[numR];

	numGridPoints   = other.numGridPoints;

	sinLatGrid      = other.sinLatGrid;
	maxSinLat       = other.maxSinLat;
	minSinLat       = other.minSinLat;

	lowerR          = other.lowerR;
	upperR          = other.upperR;
	geometricFactor = other.geometricFactor;

	this->pos		= new Vec3D[numGridPoints];
	this->B			= new Vec3D[numGridPoints];
	this->psi       = new hcFloat[numGridPoints];
	this->temp      = new hcFloat[numGridPoints];
	this->relError  = new hcFloat[numGridPoints];

	memcpy(pos,     other.getPosArray(),    numGridPoints * sizeof(Vec3D));
	memcpy(B,   	other.getBArray(),    	numGridPoints * sizeof(Vec3D));
	memcpy(psi,     other.getPsiArray(),   	numGridPoints * sizeof(hcFloat));
	memcpy(relError,other.getRelErrorArray(),numGridPoints * sizeof(hcFloat));
	memcpy(temp,    other.getTempArray(),   numGridPoints * sizeof(hcFloat));
	memcpy(a,     	other.getEllAArray(),   numR * sizeof(hcFloat));

	return *this;
}

EllipticalGrid &EllipticalGrid::operator=(const SphericalGrid &other)
{
	if(this == &other)
		return *this;

	clear();
	if(other.numGridPoints == 0)
		return *this;

	numR            = other.numR;
	numTheta        = other.numTheta;
	numPhi          = other.numPhi;

	a				= NULL;

	numGridPoints   = other.numGridPoints;

	sinLatGrid      = other.sinLatGrid;
	maxSinLat       = other.maxSinLat;
	minSinLat       = other.minSinLat;

	lowerR          = other.lowerR;
	upperR          = other.upperR;
	geometricFactor = other.geometricFactor;

	pos				= new Vec3D[numGridPoints];
	B				= new Vec3D[numGridPoints];
	psi         	= new hcFloat[numGridPoints];
	temp        	= new hcFloat[numGridPoints];
	relError    	= new hcFloat[numGridPoints];

	memcpy(pos,     other.getPosArray(),      	numGridPoints * sizeof(Vec3D));
	memcpy(B,   	other.getBArray(),    		numGridPoints * sizeof(Vec3D));
	memcpy(psi,     other.getPsiArray(),      	numGridPoints * sizeof(hcFloat));
	memcpy(relError,other.getRelErrorArray(), 	numGridPoints * sizeof(hcFloat));
	memcpy(temp,    other.getTempArray(),     	numGridPoints * sizeof(hcFloat));

	this->a			= new hcFloat[this->numR];

	for(uint i=0; i < numR; i++)
		this->a[i] 	= 1.0;

	return *this;
}

void EllipticalGrid::initNULL()
{
    initNULL_CPU();
#ifdef CUDA
    initNULL_GPU();
#endif
}

void EllipticalGrid::initNULL_CPU()
{
	SphericalGrid::initNULL_CPU();
	prolate	= false;
	a 		= NULL;
}

void EllipticalGrid::clear_CPU()
{
	SphericalGrid::clear_CPU();
	delete [] a;
	initNULL_CPU();
}

void EllipticalGrid::clear()
{
    clear_CPU();
#ifdef CUDA
    clear_GPU();
#endif
}

/*! additional information like different domains (block regular) and distribution of
 *  grid points have to be supplied
 *
 *  @param rDistribution (0 - equally spaced, 1 - geometric series)
 */
void EllipticalGrid::init(bool sinLatGrid, hcFloat maxSinLat, hcFloat minSinLat, hcFloat lowerR, hcFloat upperR,
						  uint numR, bool clearGPU, hcFloat ell)
{
	//clear();
	SphericalGrid::init(sinLatGrid, maxSinLat, minSinLat, lowerR, upperR, numR, clearGPU, ell);

	this->a				= new hcFloat[numR];
	hcFloat ru			= upperR;
	hcFloat rl			= lowerR;
#ifdef RSCALE
	ru					/= r_sol;
	rl					/= r_sol;
#endif
	if(ell < 0.0)
	{
		cerr << __FILE__ << ":" << __LINE__ << " ellipticity parameter < 0.0 (" << ell << ")\n";
		return;
	}
	else if(ell < 1.0)
	{
		prolate = true;
		ell 	= 1.0/ell;
	}
	else prolate = false;

	hcFloat as			= (ell-1.0)/(ru*ru-rl*rl);
	hcFloat alpha		= (1.0-as*rl*rl) / 2.0;

	for(uint i=0;i<numR;++i)
	{
		hcFloat r 	= pos[getIndex(i, 0, 0)][0];
		hcFloat ra	= r;
#ifdef RSCALE
		ra			/= r_sol;
#endif
		hcFloat ad	= 2*alpha + as*ra*ra;
		this->a[i]	= ad;
	}

	for(uint i=0;i<numR;++i)
		for(uint j=0;j<numTheta;++j)
			for(uint k=0;k<numPhi;++k)
			{
				uint ind 		= getIndex(i, j, k);
				Vec3D posE		= pos[ind];
				Vec3D posC		= posE.convCoordEll2Cart(*this);
				Vec3D erU, etU, epU, erL, etL, epL;
				getGradientVectors(posC, erU, etU, epU, prolate);
				getTangentVectors( posC, erL, etL, epL, prolate);

				g[ind].loadZeroes();
				p[ind].loadZeroes();
				p_r[ind].loadZeroes();
				p_t[ind].loadZeroes();
				p_p[ind].loadZeroes();

				g[ind](0,0)	= erL * erL;
				g[ind](0,1)	= erL * etL;
				g[ind](0,2)	= erL * epL;

				g[ind](1,0)	= g[ind](0,1);
				g[ind](1,1)	= etL * etL;
				g[ind](1,2)	= etL * epL;

				g[ind](2,0)	= g[ind](0,2);
				g[ind](2,1)	= g[ind](1,2);
				g[ind](2,2)	= epL * epL;

				h[ind](0,0)	= erU * erU;
				h[ind](0,1)	= erU * etU;
				h[ind](0,2)	= erU * epU;

				h[ind](1,0)	= h[ind](0,1);
				h[ind](1,1)	= etU * etU;
				h[ind](1,2)	= etU * epU;

				h[ind](2,0)	= h[ind](0,2);
				h[ind](2,1)	= h[ind](1,2);
				h[ind](2,2)	= epU * epU;

				*(Matrix<3,3,hcFloat>*)(&p[ind]) = h[ind] * sqrt(g[ind].det());
			}

	compDerivR(p, p_r);
	compDerivT(p, p_t);
	compDerivP(p, p_p);

	// make sure we don't have an invalid grid
	hcFloat rpos 	= 0.0;
	hcFloat s		= (ell-1.0)/(numR-1);

	for(uint i=0; i<numR; ++i)
	{
		hcFloat di	= i*s + 1.0;
		Vec3D posE 	= getPos(getIndex(i, 0, 0),true);
		Vec3D posS 	= posE;
		posS.convertEllipticalCoordsToCartesian(di, prolate);
		posS = posS.convCoordCart2Spher();

		if(posS[0] < rpos)
		{
			cerr << __FILE__ << "/" << __LINE__ << " The a-parameter (" << ell << ") produces an invalid grid! Cannot work with this!\n";
			fflush(stdout);
			exit(1);
		}
		rpos = posS[0];
	}

	getScalingFactors();
}

bool EllipticalGrid::isElliptical() const
{
	return true;
}

Vec3D EllipticalGrid::getPos(uint index, bool ellipticCoords) const
{
	if (index < numGridPoints)
	{
		if(ellipticCoords)
			return pos[index];
		else
		{
			uint i, j, k;
			getIndicesFromIndex(index, i, j, k);
			Vec3D sPos = pos[index];
			sPos.convertEllipticalCoordsToCartesian(a[i], prolate);
			sPos = sPos.convCoordCart2Spher();
			return sPos;
		}
	}
	else
	{
		printf("ERROR! EllipticalGrid::getPos: index (%u) out of bounds (%u)\n", index, numGridPoints);
		return Vec3D(0.0, 0.0, 0.0);
	}
}

Vec3D EllipticalGrid::getB(uint index, bool computationalCoords) const
{
	if (index < numGridPoints)
	{
		if(computationalCoords)
			return B[index];
		else
		{
			Vec3D sB 	= B[index];
			sB 			= sB.convVecEll2Spher(pos[index], *this);
			return sB;
		}
	}
	else
	{
		printf("ERROR! EllipticalGrid::getB: index (%u) out of bounds (%u)\n", index, numGridPoints);
		return Vec3D(0.0, 0.0, 0.0);
	}
}

void EllipticalGrid::setB(uint index, Vec3D value, bool ellipticCoords)
{
	if (index < numGridPoints)
	{
		if(ellipticCoords)
			this->B[index] = value;
		else
		{
			Vec3D cPos 	= pos[index].convCoordEll2Cart(*this);
			Vec3D sB 	= value.convVecSpher2Cart(cPos);
			sB 			= sB.convVecCart2Ell(cPos, *this);
			this->B[index] = sB;
		}
	}
	else
		printf("ERROR! EllipticalGrid::setB: index (%u) out of bounds (%u)\n", index, numGridPoints);
}

Vec3D EllipticalGrid::getSphericalPos(uint index)
{
	Vec3D sPos  = pos[index].convCoordEll2Cart(*this);
	sPos		= sPos.convCoordCart2Spher();
	return sPos;
}

hcFloat* EllipticalGrid::getEllAArray() const
{
	return this->a;
}

hcFloat EllipticalGrid::getEllA(uint index) const
{
	if (index < numR)
		return this->a[index];
	else
	{
		printf("ERROR! EllipticalGrid::getEllA: index (%u) out of bounds (%u)\n", index, numR);
		return 0.0;
	}
}

void EllipticalGrid::dump() const
{
	cout << "Dumping EllipticalGrid:--------------------------------\n";
	cout << setfill(' ') << left;
	cout << setw(20) << "numR:" 			<< numR 			<< "\n";
	cout << setw(20) << "numTheta:" 		<< numTheta 		<< "\n";
	cout << setw(20) << "numPhi:" 			<< numPhi 			<< "\n";
	cout << setw(20) << "numGridPoints:" 	<< numGridPoints 	<< "\n";
	cout << setw(20) << "prolate:" 			<< prolate 			<< "\n";
	cout << setw(20) << "a:" 				<< this->getEllA(this->numR - 1) << "\n";
	cout << setw(20) << "sinLatGrid:" 		<< sinLatGrid 		<< "\n";
	cout << setw(20) << "maxSinLat:" 		<< maxSinLat 		<< "\n";
	cout << setw(20) << "lowerR:" 			<< lowerR/r_sol 	<< "r_sol\n";
	cout << setw(20) << "upperR:" 			<< upperR/r_sol 	<< "r_sol\n";
	cout << setw(20) << "geomFactor:" 		<< geometricFactor 	<< "r_sol\n";
}

/*	@param pos 				position vector for which a to be computed (in spherical/elliptic coordinates, see below)
 * 	@param posElliptical	pos given in elliptic(true) or spherical(false) coordinates
 * 	@param a				elliptic parameter in x-direction
 *
 */
hcFloat EllipticalGrid::getEllParamsFromPos(Vec3D pos, bool posElliptical) const
{
	uint i;

	if(pos[0] < lowerR) return getEllA(0);

	if(posElliptical) 	// pos in elliptical coordinates
	{
		if(pos[0] > upperR) return getEllA(numR-1);

		i = 0;
		while(++i < numR && pos[0] > this->pos[getIndex(i, 0, 0)][0]);
	}
	else				// pos in spherical coordinates
	{
		Vec3D pos_cart 	= pos.convCoordSpher2Cart();
		i				= 1;
		bool condition	= false;
		do
		{
			hcFloat a 		= getEllA(i);
			Vec3D posCopy	= pos_cart;
			posCopy.convertCartesianCoordsToElliptical(a, prolate);
			condition		= posCopy[0] <= this->pos[getIndex(i, 0, 0)][0];
		}
		while(!condition && ++i < numR);
	}

	if(i==numR) return getEllA(numR-1);

	hcFloat a0		= getEllA(i-1);
	hcFloat a1		= getEllA(i);

	hcFloat r_low	= this->pos[getIndex(i-1, 0, 0)][0];
	hcFloat r_high 	= this->pos[getIndex(i  , 0, 0)][0];

	if(!posElliptical)
	{
		Vec3D pos_c	= pos.convCoordSpher2Cart();
		hcFloat x 	= pos_c[0];
		hcFloat y 	= pos_c[1];
		hcFloat z 	= pos_c[2];
		hcFloat s0	= prolate ? sqrt(r_low *r_low  / ((x*x) + (y*y) + (z*z)/(a0*a0))) : sqrt(r_low *r_low  / ((x*x)/(a0*a0) + (y*y)/(a0*a0) + (z*z)));
		hcFloat s1	= prolate ? sqrt(r_high*r_high / ((x*x) + (y*y) + (z*z)/(a1*a1))) : sqrt(r_high*r_high / ((x*x)/(a1*a1) + (y*y)/(a1*a1) + (z*z)));
		Vec3D pos0	= pos_c * s0;
		Vec3D pos1 	= pos_c * s1;
		pos0		= pos0.convCoordCart2Spher();
		pos1		= pos1.convCoordCart2Spher();
		r_low		= pos0[0];
		r_high		= pos1[0];
	}

	hcFloat drm 	= (pos[0] - r_low) / (r_high - r_low);
	hcFloat a 		= a0 + (a1 - a0) * drm;
	return a;
}

void SphericalGrid::compDerivR(Matrix3x3 *arr, Matrix3x3 *deriv)
{
	Matrix3x3 zeros;

	for(uint i=0; i<numR; ++i)
		for(uint j=0; j<numTheta; ++j)
			for(uint k=0; k<numPhi; ++k)
			{
				uint ind 		= getIndex(i,j,k);
				uint ind_p 		= i>=numR-1		? 0					: getIndex(i+1,j,k);
				uint ind_pp 	= i>=numR-2		? 0 				: getIndex(i+2,j,k);
				uint ind_m 		= i<=0			? 0 				: getIndex(i-1,j,k);
				uint ind_mm 	= i<=1			? 0					: getIndex(i-2,j,k);

				Matrix3x3 f		= arr[ind];
				Matrix3x3 fp	= i < numR-1	? arr[ind_p]		: zeros;
				Matrix3x3 fpp	= i < numR-2	? arr[ind_pp]		: zeros;
				Matrix3x3 fm	= i > 0			? arr[ind_m]		: zeros;
				Matrix3x3 fmm	= i > 1 		? arr[ind_mm]		: zeros;

				hcFloat x		= pos[ind][0];
				hcFloat xp		= i < numR-1	? pos[ind_p][0]		: 0.0;
				hcFloat xpp		= i < numR-2	? pos[ind_pp][0]	: 0.0;
				hcFloat xm		= i > 0			? pos[ind_m][0]		: 0.0;
				hcFloat xmm		= i > 1			? pos[ind_mm][0]	: 0.0;

				hcFloat dxp		= xp  - x;
				hcFloat dxm		= x   - xm;
				hcFloat dxpp 	= xpp - x;
				hcFloat dxmm	= x   - xmm;
#ifdef RSCALE
				dxp		/= r_sol;
				dxm		/= r_sol;
				dxpp 	/= r_sol;
				dxmm	/= r_sol;
#endif

				hcFloat dxp2	= dxp*dxp;
				hcFloat dxpp2	= dxpp*dxpp;
				hcFloat dxm2	= dxm*dxm;
				hcFloat dxmm2	= dxmm*dxmm;

				deriv[ind].loadZeroes();

				if(i==0)			*(Matrix<3,3,hcFloat>*)(&deriv[ind]) =  1/(dxp*dxpp2 - dxp2*dxpp) * (f*(-dxpp2+dxp2) + fp*dxpp2 - fpp*dxp2);
				else if(i==numR-1)	*(Matrix<3,3,hcFloat>*)(&deriv[ind]) = -1/(dxm*dxmm2 - dxm2*dxmm) * (f*(-dxmm2+dxm2) + fm*dxmm2 - fmm*dxm2);
				else				*(Matrix<3,3,hcFloat>*)(&deriv[ind]) =  1/(dxp*dxm2  + dxp2*dxm ) * (fm*(-dxp2) - f*(dxm2-dxp2) + fp*dxm2);
			}
}

void SphericalGrid::compDerivT(Matrix3x3 *arr, Matrix3x3 *deriv)
{
	for(uint i=0; i<numR; ++i)
		for(uint j=0; j<numTheta; ++j)
		{
			Matrix3x3 poleN;
			Matrix3x3 poleS;

			if(j==0)
			{
				for(uint k2=0; k2<numPhi;++k2)
				{
					uint indexN = getIndex(i, 0, k2);
					poleN		+= arr[indexN];
				}
				poleN /= numPhi;
			}
			else if(j==numTheta-1)
			{
				for(uint k2=0; k2<numPhi;++k2)
				{
					uint indexS = getIndex(i, numTheta-1, k2);
					poleS		+= arr[indexS];
				}
				poleS /= numPhi;
			}

			for(uint k=0; k<numPhi; ++k)
			{
				uint ind 		= getIndex(i,j,k);
				uint ind_p 		= j==numTheta-1	? 0	: getIndex(i,j+1,k);
				uint ind_m 		= j==0			? 0 : getIndex(i,j-1,k);

				Matrix3x3 f		= arr[ind];
				Matrix3x3 fp	= j < numTheta-1	? arr[ind_p]	: poleS;
				Matrix3x3 fm	= j > 0				? arr[ind_m]	: poleN;

				hcFloat x		= pos[ind][1];
				hcFloat xp		= j < numTheta-1	? pos[ind_p][1]	: PI;
				hcFloat xm		= j > 0				? pos[ind_m][1]	: 0.0;

				hcFloat dxp		= xp  - x;
				hcFloat dxm		= x   - xm;

				hcFloat dxm2	= dxm*dxm;
				hcFloat dxp2	= dxp*dxp;

				*(Matrix<3,3,hcFloat>*)(&deriv[ind]) =  1/(dxp*dxm2  + dxp2*dxm ) * (fm*(-dxp2) - f*(dxm2-dxp2) + fp*dxm2);
			}
		}
}

void SphericalGrid::compDerivP(Matrix3x3 *arr, Matrix3x3 *deriv)
{
	for(uint i=0; i<numR; ++i)
		for(uint j=0; j<numTheta; ++j)
			for(uint k=0; k<numPhi; ++k)
			{
				uint ind 		= getIndex(i,j,k);
				uint ind_p 		= getIndex(i,j,(k==numPhi-1 ? 0 		: k+1));
				uint ind_m 		= getIndex(i,j,(k==0		? numPhi-1	: k-1));

				Matrix3x3 f		= arr[ind];
				Matrix3x3 fp	= arr[ind_p];
				Matrix3x3 fm	= arr[ind_m];

				hcFloat x		= pos[ind][2];
				hcFloat xp		= pos[ind_p][2];
				hcFloat xm		= pos[ind_m][2];

				hcFloat dxp		= k==numPhi-1 	? 2*PI - x + xp : xp -x;
				hcFloat dxm		= k==0			? x + 2*PI - xm : x  - xm;

				hcFloat a		= 1 / (dxp*dxm*dxm + dxp*dxp*dxm);

				*(Matrix<3,3,hcFloat>*)(&deriv[ind]) =
					fm * a * (-dxp*dxp) - f * a * (dxm*dxm - dxp*dxp) + fp * a * dxm*dxm;
			}
}

Vec3D SphericalGrid::getBFromPsi(uint i, uint j, uint k)
{
	uint ind_ijk 	= getIndex(i,j,k);
	uint ind_imjk 	= i > 0				? getIndex(i-1,j,k)	: 0;
	uint ind_immjk 	= i > 1				? getIndex(i-2,j,k)	: 0;
	uint ind_ipjk 	= i < numR-1		? getIndex(i+1,j,k)	: 0;
	uint ind_ippjk 	= i < numR-2		? getIndex(i+2,j,k)	: 0;

	uint ind_ijmk 	= j > 0				? getIndex(i,j-1,k)	: 0;
	uint ind_ijpk 	= j < numTheta-1	? getIndex(i,j+1,k)	: 0;

	uint ind_ijkm	= getIndex(i,j,(k==0		? numPhi-1	: k-1));
	uint ind_ijkp 	= getIndex(i,j,(k==numPhi-1 ? 0 		: k+1));

	hcFloat f		= psi[ind_ijk];
	hcFloat fr_m	= i > 0				? psi[ind_imjk]		: 0.0;
	hcFloat fr_mm	= i > 1 			? psi[ind_immjk]	: 0.0;
	hcFloat fr_p	= i < numR-1		? psi[ind_ipjk]		: 0.0;
	hcFloat fr_pp	= i < numR-2		? psi[ind_ippjk]	: 0.0;

	hcFloat poleN	= 0.0;
	hcFloat poleS	= 0.0;
	for(uint k2=0; k2<numPhi;++k2)
	{
		poleN		+= psi[getIndex(i, 0,          k2)];
		poleS		+= psi[getIndex(i, numTheta-1, k2)];
	}
	poleN /= numPhi;
	poleS /= numPhi;

	hcFloat ft_m	= j > 0				? psi[ind_ijmk]		: poleN;
	hcFloat ft_p	= j < numTheta-1	? psi[ind_ijpk]		: poleS;

	hcFloat fp_p	= psi[ind_ijkp];
	hcFloat fp_m	= psi[ind_ijkm];

	Vec3D position	= pos[ind_ijk];
	hcFloat r		= position[0];
	hcFloat t		= position[1];
	hcFloat p		= position[2];

	hcFloat rm		= i > 0			? pos[ind_imjk][0]	: 0.0;
	hcFloat rmm		= i > 1			? pos[ind_immjk][0]	: 0.0;
	hcFloat rp		= i < numR-1	? pos[ind_ipjk][0]	: 0.0;
	hcFloat rpp		= i < numR-2	? pos[ind_ippjk][0]	: 0.0;

	hcFloat tm		= j > 0			? pos[ind_ijmk][1]	: 0.0;
	hcFloat tp		= j < numTheta-1? pos[ind_ijpk][1]	: PI;

	hcFloat pm		= pos[ind_ijkm][2];
	hcFloat pp		= pos[ind_ijkp][2];

#ifdef RSCALE
	r				/= r_sol;
	rm 				/= r_sol;
	rmm 			/= r_sol;
	rp 				/= r_sol;
	rpp 			/= r_sol;
#endif

	hcFloat drm		= r    - rm;
	hcFloat drmm	= r    - rmm;
	hcFloat drp		= rp   - r;
	hcFloat drpp	= rpp  - r;
	hcFloat drp2	= drp  * drp;
	hcFloat drm2	= drm  * drm;
	hcFloat drpp2	= drpp * drpp;
	hcFloat drmm2	= drmm * drmm;

	hcFloat dtm		= t    - tm;
	hcFloat dtp		= tp   - t;
	hcFloat dtm2	= dtm  * dtm;
	hcFloat dtp2	= dtp  * dtp;

	hcFloat dpm		= k==0			? p + 2*PI - pm : p  - pm;
	hcFloat dpp		= k==numPhi-1 	? 2*PI - p + pp : pp - p;
	hcFloat dpm2	= dpm*dpm;
	hcFloat dpp2	= dpp*dpp;

	hcFloat del_r = (i==0 		?  	1.0/(drp*drpp2 - drp2*drpp) * (f*(-drpp2+drp2) + fr_p*drpp2    - fr_pp*drp2) :
					(i==numR-1 	?  -1.0/(drm*drmm2 - drm2*drmm) * (f*(-drmm2+drm2) + fr_m*drmm2    - fr_mm*drm2) :
									1.0/(drp*drm2  + drp2*drm ) * (fr_m*(-drp2)    - f*(drm2-drp2) + fr_p*drm2)));
	hcFloat del_t =  				1.0/(dtp*dtm2  + dtp2*dtm ) * (ft_m*(-dtp2)    - f*(dtm2-dtp2) + ft_p*dtm2);
	hcFloat del_p =  				1.0/(dpp*dpm2  + dpp2*dpm ) * (fp_m*(-dpp2)    - f*(dpm2-dpp2) + fp_p*dpm2);

	/*
	if(i==0)			del_r	=  1.0/(drp*drpp2 - drp2*drpp) * (f*(-drpp2+drp2) + fr_p*drpp2    - fr_pp*drp2);
	else if(i==numR-1)	del_r	= -1.0/(drm*drmm2 - drm2*drmm) * (f*(-drmm2+drm2) + fr_m*drmm2    - fr_mm*drm2);
	else				del_r	=  1.0/(drp*drm2  + drp2*drm ) * (fr_m*(-drp2)    - f*(drm2-drp2) + fr_p*drm2);//*/

	Vec3D retval;
#ifdef SPHERICUNITVEC
	if(!isElliptical())
		retval			= Vec3D((hcFloat)-del_r, (hcFloat)-1.0/r * del_t, (hcFloat)-1.0/(r*sin(t)) * del_p);
	else
#endif
		retval			= Vec3D(-del_r, -del_t, -del_p);

	return retval;
}

bool EllipticalGrid::getGradientVectors(Vec3D cartPos, Vec3D &er, Vec3D &et, Vec3D &ep, bool prolate) const
{
	Vec3D posE		= cartPos.convCoordCart2Ell(*this);
	hcFloat rd		= posE[0];
	hcFloat thetad	= posE[1];
	hcFloat phid	= posE[2];

	hcFloat rud		= upperR;
	hcFloat rld		= lowerR;
#ifdef RSCALE
	rd 				/= r_sol;
	rud				/= r_sol;
	rld				/= r_sol;
#endif

	hcFloat std		= sin(thetad);
	hcFloat spd		= sin(phid);
	hcFloat ctd		= cos(thetad);
	hcFloat cpd		= cos(phid);

	hcFloat as		= (this->a[numR-1]-1.0)/(rud*rud-rld*rld);
	hcFloat alpha	= (1.0-as*rld*rld) / 2.0;
	hcFloat ad		= 2*alpha + as*rd*rd;
	hcFloat dad_drd = 2*as*rd;

	hcFloat drd_dx	= prolate ? 1.0/(ad + ctd*ctd*dad_drd*rd)	* (ad * std * cpd) 	: 1.0/(ad + std*std*dad_drd*rd)	* ( std * cpd);
	hcFloat drd_dy	= prolate ? 1.0/(ad + ctd*ctd*dad_drd*rd)	* (ad * std * spd)	: 1.0/(ad + std*std*dad_drd*rd)	* ( std * spd);
	hcFloat drd_dz	= prolate ? 1.0/(ad + ctd*ctd*dad_drd*rd)	*  ctd				: 1.0/(ad + std*std*dad_drd*rd)	* ( ctd * ad );

	hcFloat dtd_dx	= prolate ? 1.0/(rd * (ad+rd*ctd*ctd*dad_drd)) 	* ( ctd * cpd * (ad+rd*dad_drd)) : 1.0/(rd*(ad + std*std*dad_drd*rd)) * ( ctd * cpd);
	hcFloat dtd_dy	= prolate ? 1.0/(rd * (ad+rd*ctd*ctd*dad_drd)) 	* ( ctd * spd * (ad+rd*dad_drd)) : 1.0/(rd*(ad + std*std*dad_drd*rd)) * ( ctd * spd);
	hcFloat dtd_dz	= prolate ? 1.0/(rd * (ad+rd*ctd*ctd*dad_drd)) 	* (-std)						 : 1.0/(rd*(ad + std*std*dad_drd*rd)) * (-std * (ad+rd*dad_drd));

	hcFloat dpd_dx	= prolate ? 1.0/(rd*std) * (-spd)			: 1.0/(ad*rd*std) * (-spd);
	hcFloat dpd_dy	= prolate ? 1.0/(rd*std) * (cpd)			: 1.0/(ad*rd*std) * ( cpd);
	hcFloat dpd_dz	= 0.0;

	er				= 													 Vec3D(drd_dx, drd_dy, drd_dz);
	et				= rd<num_eps || std<num_eps	? Vec3D(1.0, 0.0, 0.0) : Vec3D(dtd_dx, dtd_dy, dtd_dz);
	ep				= rd<num_eps || std<num_eps ? Vec3D(0.0, 1.0, 0.0) : Vec3D(dpd_dx, dpd_dy, dpd_dz);

	return true;
}

bool EllipticalGrid::getTangentVectors(Vec3D cartPos, Vec3D &er, Vec3D &et, Vec3D &ep, bool prolate) const
{
	Vec3D posE		= cartPos.convCoordCart2Ell(*this);

	hcFloat rd		= posE[0];
	hcFloat thetad	= posE[1];
	hcFloat phid	= posE[2];

	hcFloat rud		= upperR;
	hcFloat rld		= lowerR;
#ifdef RSCALE
	rd 				/= r_sol;
	rud				/= r_sol;
	rld				/= r_sol;
#endif

	hcFloat std		= sin(thetad);
	hcFloat spd		= sin(phid);
	hcFloat ctd		= cos(thetad);
	hcFloat cpd		= cos(phid);

	hcFloat as		= (this->a[numR-1]-1.0)/(rud*rud-rld*rld);
	hcFloat alpha	= (1.0-as*rld*rld) / 2.0;
	hcFloat ad		= 2*alpha + as*rd*rd;
	hcFloat dad_drd = 2*as*rd;

	hcFloat dx_drd	= prolate ? std*cpd					: std*cpd * (ad + dad_drd*rd);
	hcFloat dy_drd	= prolate ? std*spd 				: std*spd * (ad + dad_drd*rd);
	hcFloat dz_drd	= prolate ? ctd * (ad+rd*dad_drd)	: ctd;

	hcFloat dx_dtd	= prolate ? rd * ctd*cpd			: rd *  ad*ctd*cpd;
	hcFloat dy_dtd	= prolate ? rd * ctd*spd			: rd *  ad*ctd*spd;
	hcFloat dz_dtd	= prolate ? rd * ad*(-std)			: rd * (-std);

	hcFloat dx_dpd	= prolate ? rd * std * (-spd)		: ad*rd*std * (-spd);
	hcFloat dy_dpd	= prolate ? rd * std * (cpd)		: ad*rd*std *   cpd;
	hcFloat dz_dpd	= 0.0;

	er				= Vec3D(dx_drd, dy_drd, dz_drd);
	et				= Vec3D(dx_dtd, dy_dtd, dz_dtd);
	ep				= Vec3D(dx_dpd, dy_dpd, dz_dpd);

	return true;
}


//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Vec3D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


Vec3D Vec3D::convCoordCart2Ell(const EllipticalGrid &grid) const
{
	Vec3D retval 	= *this;
	Vec3D pos		= this->convCoordCart2Spher();
	hcFloat a		= grid.getEllParamsFromPos(pos, false);

	if(grid.prolate)
	{
		retval(2)	/= a;
	}
	else
	{
		retval(0) 	/= a;
		retval(1) 	/= a;
	}
	retval 			= retval.convCoordCart2Spher();
	return retval;
}

Vec3D Vec3D::convCoordEll2Cart(const EllipticalGrid &grid) const
{
	Vec3D retval 	= *this;
	hcFloat a		= grid.getEllParamsFromPos(retval, true);
	retval			= retval.convCoordSpher2Cart();
	//printf("prolate: %u, a: %E\n", grid.prolate, a);
	//printf("%E %E %E\n", retval[0], retval[1], retval[2]);
	if(grid.prolate)
	{
		retval(2)	*= a;
	}
	else
	{
		retval(0)	*= a;
		retval(1)	*= a;
	}
	//printf("%E %E %E\n", retval[0], retval[1], retval[2]);
	return retval;
}

Vec3D Vec3D::convVecCart2Ell(const Vec3D &cartPos, const EllipticalGrid &grid) const
{
	// TODO: this method is inefficient because I didn't have time to find
	// the appropriate transformations and worked with SLE instead

	Vec3D er, et, ep;
	grid.getGradientVectors(cartPos, er, et, ep, grid.prolate);

	Matrix<3,4,double> M;
	M(0,0) = er[0];	M(0,1) = et[0];	M(0,2) = ep[0]; M(0,3) = content[0];
	M(1,0) = er[1];	M(1,1) = et[1];	M(1,2) = ep[1]; M(1,3) = content[1];
	M(2,0) = er[2];	M(2,1) = et[2];	M(2,2) = ep[2]; M(2,3) = content[2];

	Vec3D retval;
	uint sol = M.solveSLE(retval);

	/*
	Vec3D inv = retval.convVecEll2Cart(cartPos, grid);
	if(fabs((*this)[0]-inv[0])>1E-4 || fabs((*this)[1]-inv[1])>1E-4 || fabs((*this)[2]-inv[2])>1E-4)
		printf("VecC2E:\n%E %E %E\n%E %E %E\n", (*this)[0], (*this)[1], (*this)[2], inv[0], inv[1], inv[2]);//*/

	return retval;
}

Vec3D Vec3D::convVecEll2Cart(const Vec3D &cartPos, const EllipticalGrid &grid) const
{
	Vec3D er, et, ep;
	grid.getGradientVectors(cartPos, er, et, ep, grid.prolate);

	Vec3D retval 	= er * content[0] + et * content[1] + ep * content[2];

	return retval;
}

Vec3D Vec3D::convVecSpher2Ell(const Vec3D &posSpher, const EllipticalGrid &grid) const
{
	Vec3D posC = posSpher.convCoordSpher2Cart();
	return this->convVecSpher2Cart(posC).convVecCart2Ell(posC, grid);
}

Vec3D Vec3D::convVecEll2Spher(const Vec3D &posEll, const EllipticalGrid &grid) const
{
	Vec3D posC	= posEll.convCoordEll2Cart(grid);
	return this->convVecEll2Cart(posC, grid).convVecCart2Spher(posC);
}

Vec3D Vec3D::convCoordSpher2Ell(EllipticalGrid &grid) const
{
	return this->convCoordSpher2Cart().convCoordCart2Ell(grid);
}

Vec3D Vec3D::convCoordEll2Spher(EllipticalGrid &grid) const
{
	return this->convCoordEll2Cart(grid).convCoordCart2Spher();
}

void EllipticalGrid::convertMagMapping(MagMapping &map)
{
	for(uint j=0; j<map.numTheta; ++j)
		for(uint k=0; k<map.numPhi; ++k)
		{
			Magline &magl 	= map(j,k);

			for(uint i=0; i<magl.numSamples; ++i)
			{
				Vec3D posE		= magl.posdata[i];
				Vec3D magE		= magl.magdata[i];

				Vec3D posC		= posE.convCoordEll2Cart(*this);
				Vec3D magC 		= magE.convVecEll2Cart(posC, *this);

				Vec3D posS		= posC.convCoordCart2Spher();
				Vec3D magS		= magC.convVecCart2Spher(posC);

				Vec3D magS2		= magE.convVecEll2Spher(posE, *this);

				magl.posdata[i]	= posS;
				magl.magdata[i]	= magS;

				/*
				if(j==2&&k==10)
				{
					if(i==0)
					{
						Vec3D rE, tE, pE, rS, tS, pS;
						getSphericGradientVectors(posC, rS, tS, pS);
						getGradientVectors(posC, rE, tE, pE);
						printf("%E %E %E - %E %E %E - %E %E %E\n%E %E %E - %E %E %E - %E %E %E\n\n",
								rE[0], rE[1], rE[2], tE[0], tE[1], tE[2], pE[0], pE[1], pE[2],
								rS[0], rS[1], rS[2], tS[0], tS[1], tS[2], pS[0], pS[1], pS[2]);
					}
					printf("%u E:%E %E %E S1:%E %E %E\n", i, magE[0], magE[1], magE[2], magS[0], magS[1], magS[2]);
				}//*/


				/*// TODO: this was written by Lasse and not tested, so there are probably still errors in here
				Vec3D* Spline, *dSpline;
				uint numParamSamples = 0;

				if(i == magl.numSamples - 1 && fabs(posE[0] - solver.grid->upperR) / solver.grid->upperR < 1)
				{
					if(debug)
						printf("MagMapping::createAtHeight: Samples before: %u\n", magl.numSamples);
					// Cubic Spline Interpolation to a sphere
					// Determine Radius of the smallest sphere containing the whole ellipsoid
					hcFloat high_a 	= eGrid->getEllA(solver.grid->numR - 1);
					hcFloat high_b 	= eGrid->getEllB(solver.grid->numR - 1);
					hcFloat high_c 	= eGrid->getEllC(solver.grid->numR - 1);
					hcFloat rSphere = high_a >= high_b ? (high_a >= high_c ? high_a : high_c) : (high_b >= high_c ? high_b : high_c);
					if(debug)
					{
						printf("MagMapping::createAtHeight: Sphere Radius: %E\n", rSphere);
						printf("MagMapping::createAtHeight: Scale Parameter: %E\n", (rSphere/(posS[0]/posE[0])));
					}
					// Determine start and end positions and the directions of the magline there
					Vec3D startPos 	= posC;
					Vec3D startDir 	= -magC;
					startDir 		/=  startDir.length();
					// Mit normalen Vektor der Ellipse vergleichen nach: https://math.stackexchange.com/questions/1927334/how-to-compute-the-normal-to-the-ellipsoid-at-the-point-on-the-surface-of-ellips
					if(debug)
					{
						Vec3D normal 	= Vec3D(2*high_a*startPos[0], 2*high_b*startPos[1], 2*high_c*startPos[2]);
						normal 			/= normal.length();
						printf("\nMagMapping::createAtHeight: startDir(cart,normed): %E %E %E\n", startDir[0], startDir[1], startDir[2]);
						printf("MagMapping::createAtHeight: normalEl(cart,normed): %E %E %E\n\n", normal[0], normal[1], normal[2]);
					}
					Vec3D endPos 	= posS;
					endPos(0) 		= endPos[0] * (rSphere/(posS[0]/posE[0])); // TODO ???
					Vec3D endDir 	= endPos;
					endPos.convertSphericalCoordsToCartesian();
					endDir.convertSphericalToCartesian(endPos);
					endDir 			/= endDir.length();
					startDir 		*= startPos.dist(endPos)/2.0;
					endDir  		*= startPos.dist(endPos)/2.0;
					if(debug)
					{
						printf("MagMapping::createAtHeight: startPos(cart): %E %E %E\n", startPos[0], startPos[1],startPos[2]);
						printf("MagMapping::createAtHeight: startDir(cart): %E %E %E\n", startDir[0], startDir[1],startDir[2]);
						printf("MagMapping::createAtHeight: endPos(cart): %E %E %E\n", endPos[0], endPos[1],endPos[2]);
						printf("MagMapping::createAtHeight: endDir(cart): %E %E %E\n", endDir[0], endDir[1],endDir[2]);
					}
					// Using cubic hermite splines like http://run.usc.edu/cs420-s14/lec08-splines/08-splines-6up.pdf, Page 3
					// Parametric spline model p(t) = At³ + Bt² + Ct + D, where p, A, B, C, D are vectors and t is scalar parameter between 0 and 1
					// Calculate the parameters
					if(!(startPos.dist(endPos)/1E7 < 1E-2))
					{
						Vec3D D = startPos;
						Vec3D C = startDir;
						Vec3D B = (3 * endPos) - endDir - (2 * C) - (3 * D);
						Vec3D A = endPos - B - C - D;
						// Determine the number of samples
						numParamSamples = (startPos.dist(endPos)/1E7) + 1;
						// Check if there are enough samples left to store the data to the magline
						if(magl.numSamples + numParamSamples > maxNumPointsPerMagline)
							numParamSamples = maxNumPointsPerMagline - magl.numSamples;

						if(debug)
							printf("MagMapping::createAtHeight: Additional Samples: %u\n\n", numParamSamples);
						// Numeric interpolation of the spline and its direction
						Spline 	= new Vec3D[numParamSamples];
						dSpline = new Vec3D[numParamSamples];

						for(uint m = 0; m < numParamSamples; ++m)
						{
							hcFloat t 		= m * (1.0 / (numParamSamples - 1));
							Vec3D splineC 	= (A * pow(t, 3)) + (B * pow(t, 2)) + (C * t) + D;
							Vec3D splineS 	= splineC;
							splineS.convertCartesianCoordsToSpherical();
							Spline[m] 		= splineS;
							Vec3D dsplineC 	= (3 * A * pow(t, 2)) + (2 * B * t) + C;
							Vec3D dsplineS 	= dsplineC;
							dsplineS.convertCartesianToSpherical(splineC);
							dSpline[m] 		= dsplineS;

							if(debug)
								printf("%u %E %E %E %E %E %E\n", i + m + 1, splineC[0], splineC[1], splineC[2], splineS[0], splineS[1], splineS[2]);
						}
					}
					else if(debug)
						printf("MagMapping::createAtHeight: Additional Samples: %u\n", numParamSamples);
				}

				// Store Spline data to magline
				if (numParamSamples > 0)
				{
					Vec3D *temp_posdata;
					Vec3D *temp_magdata;

					temp_posdata = magl.posdata;
					temp_magdata = magl.magdata;

					magl.posdata = new Vec3D[magl.numSamples + numParamSamples];
					magl.magdata = new Vec3D[magl.numSamples + numParamSamples];

					uint m;
					for(m = 0; m < magl.numSamples; ++m)
					{
						magl.posdata[m] = temp_posdata[m];
						magl.magdata[m] = temp_magdata[m];
					}
					for(m = magl.numSamples; m < magl.numSamples + numParamSamples; ++m)
					{
						magl.posdata[m] = Spline[m - magl.numSamples];
						magl.magdata[m] = dSpline[m - magl.numSamples];
					}

					magl.numSamples += numParamSamples;

					if(debug)
						printf("MagMapping::createAtHeight: Samples after: %u\n", magl.numSamples);

					delete [] temp_posdata;
					delete [] temp_magdata;

					delete [] Spline;
					delete [] dSpline;

					break;
				}//*/
			}
		}
}

