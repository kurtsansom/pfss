#include "src/grids.h"
#include "src/pfss.h"

#include <fstream>
#include <stdio.h>

SphericalGrid::SphericalGrid()
{
    initNULL();
}

SphericalGrid::SphericalGrid(const SphericalGrid &other)
{
    initNULL();
    *this = other;
}

SphericalGrid::~SphericalGrid()
{
	clear();
}

SphericalGrid &SphericalGrid::operator=(const SphericalGrid &other)
{
    if(this == &other)
        return *this;

    clear();
    if(other.numGridPoints == 0)
        return *this;

    numR            = other.numR;
    numTheta        = other.numTheta;
    numPhi          = other.numPhi;

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

    return *this;
}

void SphericalGrid::initNULL()
{
    initNULL_CPU();
#ifdef CUDA
    initNULL_GPU();
#endif
}

void SphericalGrid::initNULL_CPU()
{

	pos				= NULL;
	B				= NULL;
    psi             = NULL;
    temp            = NULL;
    relError        = NULL;
    g				= NULL;
    h				= NULL;
    p				= NULL;
    p_r				= NULL;
    p_t				= NULL;
    p_p				= NULL;
    s_ijk 			= NULL;
	s_imjk 			= NULL;
	s_ipjk 			= NULL;
	s_ijmk 			= NULL;
	s_ijpk 			= NULL;
	s_ijkm 			= NULL;
	s_ijkp 			= NULL;
	s_imjmk 		= NULL;
	s_imjpk 		= NULL;
	s_ipjmk 		= NULL;
	s_ipjpk 		= NULL;
	s_imjkm 		= NULL;
	s_imjkp 		= NULL;
	s_ipjkm 		= NULL;
	s_ipjkp 		= NULL;
	s_ijmkm 		= NULL;
	s_ijmkp 		= NULL;
	s_ijpkm 		= NULL;
	s_ijpkp 		= NULL;

	sizeofFloat		= 0;

    numR            = 0;
    numTheta        = 0;
    numPhi          = 0;

    sinLatGrid      = false;
    maxSinLat       = 0.0;
    minSinLat       = 0.0;

    numGridPoints   = 0;

    lowerR          = 0.0;
    upperR          = 0.0;
    geometricFactor = 1.0;

    harmSolution	= NULL;
}


void SphericalGrid::clear()
{
    clear_CPU();
#ifdef CUDA
    clear_GPU();
#endif
}

void SphericalGrid::clear_CPU()
{
	delete [] pos;
	delete [] B;
    delete [] psi;
    delete [] temp;
    delete [] relError;
    delete [] g;
    delete [] h;
    delete [] p;
    delete [] p_r;
    delete [] p_t;
    delete [] p_p;
    delete [] s_ijk;
    delete [] s_imjk;
    delete [] s_ipjk;
    delete [] s_ijmk;
    delete [] s_ijpk;
    delete [] s_ijkm;
    delete [] s_ijkp;
    delete [] s_imjmk;
    delete [] s_imjpk;
    delete [] s_ipjmk;
    delete [] s_ipjpk;
    delete [] s_imjkm;
    delete [] s_imjkp;
    delete [] s_ipjkm;
    delete [] s_ipjkp;
    delete [] s_ijmkm;
    delete [] s_ijmkp;
    delete [] s_ijpkm;
    delete [] s_ijpkp;

    initNULL_CPU();
}

/*! additional information like different domains (block regular) and distribution of
 *  grid points have to be supplied
 *
 *  @param rDistribution (false - equally spaced, true - geometric series)
 */
void SphericalGrid::init(bool sinLatGrid, hcFloat maxSinLat, hcFloat minSinLat, hcFloat lowerR, hcFloat upperR,
						 uint numR, bool clearGPU, hcFloat a)
{
	// TODO: ifdef
	if(clearGPU)	clear();
	else			clear_CPU();

	uint numT, numP;
	hcFloat geometricFactor;
	getOptimalGridParameters(numR, lowerR, upperR, geometricFactor, numT, numP);
	initMembers(numR, numT, numP, sinLatGrid, maxSinLat, minSinLat, lowerR, upperR);
	this->geometricFactor 	= geometricFactor;

	hcFloat dr 			= (1-this->geometricFactor) / (1-pow(this->geometricFactor, numR-1)) * (upperR - lowerR);
	hcFloat dSinTheta   = (maxSinLat-minSinLat) / (numT-1);

	// this structures take care of interpolation in theta and phi direction. Interpolation is linear in phi AND theta
	// (instead of linear in phi and sin(theta))
	for(uint i=0;i<numR;++i)
	{
		hcFloat r 	= i==0 ? lowerR : lowerR + dr * (1-pow(this->geometricFactor, i)) / (1-this->geometricFactor);
		hcFloat ra	= r;
#ifdef RSCALE
		ra			/= r_sol;
#endif
		for(uint j=0;j<numT;++j)
		{
			hcFloat theta = PI/2.0 - asin(maxSinLat - j * dSinTheta);

			for(uint k=0;k<numP;++k)
			{
				hcFloat phi		= k             * 2 * PI / numP;
#ifdef PHINOTCENTERED
				hcFloat phi 	= (k + 1.0/2.0) * 2 * PI / numP;
#endif
				uint ind 		= k * numR * numT + j * numR + i;

				psi[ind]        = 0.0;
				temp[ind]       = 0.0;
				relError[ind]   = 0.0;
				B[ind]			= Vec3D(0.0, 0.0, 0.0);
				pos[ind]		= Vec3D(r, theta, phi);

				g[ind].loadZeroes();
				p[ind].loadZeroes();
				p_r[ind].loadZeroes();
				p_t[ind].loadZeroes();
				p_p[ind].loadZeroes();
				g[ind](0,0) = 1;
				g[ind](1,1)	= ra*ra;
				g[ind](2,2) = ra*ra*sin(theta)*sin(theta);

				h[ind](0,0) = 1;
				h[ind](1,1)	= 1/(ra*ra);
				h[ind](2,2) = 1/(ra*ra*sin(theta)*sin(theta));

				*(Matrix<3,3,hcFloat>*)(&p[ind])		= h[ind] * sqrt(g[ind].det());
			}
		}
	}

	compDerivR(p, p_r);
	compDerivT(p, p_t);
	compDerivP(p, p_p);

	getScalingFactors();
}

void SphericalGrid::initMembers(uint numR, uint numT, uint numP,
		bool sinLatGrid, hcFloat maxSinLat, hcFloat minSinLat, hcFloat lowerR, hcFloat upperR)
{
	if(numR==0 || numT == 0 || numP==0)
	{
		printf("ERROR! SphericalGrid::init Number of grid points in each dimension must be >0!\nnumR: %u, numTheta: %u, numPhi: %u\n", numR, numT, numP);
		return;
	}

	sizeofFloat			= sizeof(hcFloat);

	numGridPoints   	= numR * numT * numP;

	pos					= new Vec3D[numGridPoints];
	B					= new Vec3D[numGridPoints];
	psi         		= new hcFloat[numGridPoints];
	temp        		= new hcFloat[numGridPoints];
	relError    		= new hcFloat[numGridPoints];
	g					= new Matrix3x3[numGridPoints];
	h					= new Matrix3x3[numGridPoints];
	p					= new Matrix3x3[numGridPoints];
	p_r					= new Matrix3x3[numGridPoints];
	p_t					= new Matrix3x3[numGridPoints];
	p_p					= new Matrix3x3[numGridPoints];
	s_ijk				= new hcFloat[numGridPoints];
	s_imjk				= new hcFloat[numGridPoints];
	s_ipjk				= new hcFloat[numGridPoints];
	s_ijmk				= new hcFloat[numGridPoints];
	s_ijpk				= new hcFloat[numGridPoints];
	s_ijkm				= new hcFloat[numGridPoints];
	s_ijkp				= new hcFloat[numGridPoints];
	s_imjmk				= new hcFloat[numGridPoints];
	s_imjpk				= new hcFloat[numGridPoints];
	s_ipjmk				= new hcFloat[numGridPoints];
	s_ipjpk				= new hcFloat[numGridPoints];
	s_imjkm				= new hcFloat[numGridPoints];
	s_imjkp				= new hcFloat[numGridPoints];
	s_ipjkm				= new hcFloat[numGridPoints];
	s_ipjkp				= new hcFloat[numGridPoints];
	s_ijmkm				= new hcFloat[numGridPoints];
	s_ijmkp				= new hcFloat[numGridPoints];
	s_ijpkm				= new hcFloat[numGridPoints];
	s_ijpkp				= new hcFloat[numGridPoints];

	this->numR          = numR;
	this->numTheta      = numT;
	this->numPhi        = numP;
	this->sinLatGrid    = sinLatGrid;
	this->maxSinLat     = maxSinLat;
	this->minSinLat     = minSinLat;

	this->lowerR        = lowerR;
	this->upperR        = upperR;
}

void SphericalGrid::getScalingFactors()
{
	uint numR	= this->numR;
	uint numT	= this->numTheta;
	uint numP	= this->numPhi;

	for(uint i=1; i<numR-1; ++i)		// we dont need scaling factors for spherical boundaries
		for(uint j=0; j<numT; ++j)
			for(uint k=0; k<numP; ++k)
			{
				uint km			= k==0		? numP-1: k-1;
				uint kp			= k==numP-1	? 0		: k+1;

				uint ind		= getIndex(i,j,k);

				uint ind_imjk	= 				  	  getIndex(i-1,	j,	 k );
				uint ind_ipjk   = 				  	  getIndex(i+1,	j,	 k );

				uint ind_ijmk	= j==0		? 0		: getIndex(i,	j-1, k );
				uint ind_ijpk	= j==numT-1	? 0		: getIndex(i,	j+1, k );

				uint ind_ijkm  	= 				  	  getIndex(i,	j,	 km);
				uint ind_ijkp  	= 				  	  getIndex(i,	j,	 kp);

				uint ind_imjmk	= j==0		? 0		: getIndex(i-1,	j-1, k );
				uint ind_imjpk	= j==numT-1	? 0		: getIndex(i-1,	j+1, k );
				uint ind_ipjmk	= j==0		? 0		: getIndex(i+1,	j-1, k );
				uint ind_ipjpk	= j==numT-1	? 0		: getIndex(i+1,	j+1, k );

				uint ind_imjkm	= 				  	  getIndex(i-1,	j,	 km);
				uint ind_imjkp	= 				  	  getIndex(i-1,	j,	 kp);
				uint ind_ipjkm	= 				  	  getIndex(i+1,	j,	 km);
				uint ind_ipjkp	= 				  	  getIndex(i+1,	j,	 kp);

				uint ind_ijmkm	= j==0		? 0		: getIndex(i,	j-1, km);
				uint ind_ijmkp	= j==0		? 0		: getIndex(i,	j-1, kp);
				uint ind_ijpkm	= j==numT-1	? 0		: getIndex(i,	j+1, km);
				uint ind_ijpkp	= j==numT-1	? 0		: getIndex(i,	j+1, kp);

				Vec3D pos_ijk	= pos[ind];

				Vec3D pos_imjk	= pos[ind_imjk];
				Vec3D pos_ipjk	= pos[ind_ipjk];

				Vec3D pos_ijmk	= pos[ind_ijmk];
				Vec3D pos_ijpk	= pos[ind_ijpk];

				Vec3D pos_ijkm	= pos[ind_ijkm];
				Vec3D pos_ijkp	= pos[ind_ijkp];

				hcFloat dr_p    = pos_ipjk[0] - pos_ijk[0] ;
				hcFloat dr_m    = pos_ijk[0]  - pos_imjk[0];

			#ifdef RSCALE
				dr_p    /=  r_sol;
				dr_m    /=  r_sol;
			#endif

				hcFloat dt_p	= j==numT-1	? PI - pos_ijk[1] 	  				: pos_ijpk[1]  	- pos_ijk[1];
				hcFloat dt_m	= j==0		? pos_ijk[1]		  				: pos_ijk[1] 	- pos_ijmk[1];

				hcFloat dp_p	= k==numP-1 ? pos_ijkp[2]  + 2*PI - pos_ijk[2] 	: pos_ijkp[2]  	- pos_ijk[2];
				hcFloat dp_m	= k==0		? 2*PI - pos_ijkm[2]  + pos_ijk[2] 	: pos_ijk[2]    - pos_ijkm[2];

				hcFloat dr_p2	= dr_p*dr_p;
				hcFloat dr_m2	= dr_m*dr_m;
				hcFloat dt_p2	= dt_p*dt_p;
				hcFloat dt_m2	= dt_m*dt_m;
				hcFloat dp_p2	= dp_p*dp_p;
				hcFloat dp_m2	= dp_m*dp_m;

				hcFloat Cr		= 1 / (dr_p*dr_m2 + dr_p2*dr_m);
				hcFloat Ct		= 1 / (dt_p*dt_m2 + dt_p2*dt_m);
				hcFloat Cp		= 1 / (dp_p*dp_m2 + dp_p2*dp_m);

				hcFloat t_r			= Cr*(p_r[ind](0,0) + p_t[ind](1,0) + p_p[ind](2,0));
				hcFloat t_t			= Ct*(p_r[ind](0,1) + p_t[ind](1,1) + p_p[ind](2,1));
				hcFloat t_p			= Cp*(p_r[ind](0,2) + p_t[ind](1,2) + p_p[ind](2,2));
				hcFloat t_rt		= 2*p[ind](0,1)*Cr*Ct;
				hcFloat t_rp		= 2*p[ind](0,2)*Cr*Cp;
				hcFloat t_tp		= 2*p[ind](1,2)*Ct*Cp;

				hcFloat s_imjk		= -t_r*dr_p2 + t_rt*dr_p2*(dt_m2-dt_p2) + t_rp*dr_p2*(dp_m2-dp_p2) + 2*p[ind](0,0)*Cr*dr_p;
				hcFloat s_ipjk		=  t_r*dr_m2 - t_rt*dr_m2*(dt_m2-dt_p2) - t_rp*dr_m2*(dp_m2-dp_p2) + 2*p[ind](0,0)*Cr*dr_m;

				hcFloat s_ijmk		= -t_t*dt_p2 + t_rt*(dr_m2-dr_p2)*dt_p2 + t_tp*dt_p2*(dp_m2-dp_p2) + 2*p[ind](1,1)*Ct*dt_p;
				hcFloat s_ijpk		=  t_t*dt_m2 - t_rt*(dr_m2-dr_p2)*dt_m2 - t_tp*dt_m2*(dp_m2-dp_p2) + 2*p[ind](1,1)*Ct*dt_m;

				hcFloat s_ijkm		= -t_p*dp_p2 + t_rp*(dr_m2-dr_p2)*dp_p2 + t_tp*(dt_m2-dt_p2)*dp_p2 + 2*p[ind](2,2)*Cp*dp_p;
				hcFloat s_ijkp		=  t_p*dp_m2 - t_rp*(dr_m2-dr_p2)*dp_m2 - t_tp*(dt_m2-dt_p2)*dp_m2 + 2*p[ind](2,2)*Cp*dp_m;

				hcFloat s_imjmk		=  t_rt*dr_p2*dt_p2;
				hcFloat s_imjpk		= -t_rt*dr_p2*dt_m2;
				hcFloat s_ipjmk		= -t_rt*dr_m2*dt_p2;
				hcFloat s_ipjpk		=  t_rt*dr_m2*dt_m2;

				hcFloat s_imjkm		=  t_rp*dr_p2*dp_p2;
				hcFloat s_imjkp		= -t_rp*dr_p2*dp_m2;
				hcFloat s_ipjkm		= -t_rp*dr_m2*dp_p2;
				hcFloat s_ipjkp		=  t_rp*dr_m2*dp_m2;

				hcFloat s_ijmkm		=  t_tp*dt_p2*dp_p2;
				hcFloat s_ijmkp		= -t_tp*dt_p2*dp_m2;
				hcFloat s_ijpkm		= -t_tp*dt_m2*dp_p2;
				hcFloat s_ijpkp		=  t_tp*dt_m2*dp_m2;

				hcFloat s_ijk		= -t_r *(dr_m2-dr_p2) - t_t*(dt_m2-dt_p2) - t_p*(dp_m2-dp_p2)
									+  t_rt*(dr_m2-dr_p2)*(dt_m2-dt_p2)
									+  t_rp*(dr_m2-dr_p2)*(dp_m2-dp_p2)
									+  t_tp*(dt_m2-dt_p2)*(dp_m2-dp_p2)
									- 2*p[ind](0,0)*Cr*(dr_p+dr_m) - 2*p[ind](1,1)*Ct*(dt_p+dt_m) - 2*p[ind](2,2)*Cp*(dp_p+dp_m);

				//printf("%u %u %u %E %E\n", i,j,k, dr_p, pos_ipjk[0]);

				this->s_ijk[ind]	= s_ijk;
				this->s_imjk[ind] 	= s_imjk;
				this->s_ipjk[ind] 	= s_ipjk;

				this->s_ijmk[ind] 	= s_ijmk;
				this->s_ijpk[ind] 	= s_ijpk;

				this->s_ijkm[ind] 	= s_ijkm;
				this->s_ijkp[ind] 	= s_ijkp;

				this->s_imjmk[ind] 	= s_imjmk;
				this->s_imjpk[ind] 	= s_imjpk;
				this->s_ipjmk[ind] 	= s_ipjmk;
				this->s_ipjpk[ind] 	= s_ipjpk;

				this->s_imjkm[ind] 	= s_imjkm;
				this->s_imjkp[ind] 	= s_imjkp;
				this->s_ipjkm[ind] 	= s_ipjkm;
				this->s_ipjkp[ind] 	= s_ipjkp;

				this->s_ijmkm[ind] 	= s_ijmkm;
				this->s_ijmkp[ind] 	= s_ijmkp;
				this->s_ijpkm[ind] 	= s_ijpkm;
				this->s_ijpkp[ind] 	= s_ijpkp;
			}
}


void SphericalGrid::clearValues()
{
    if(pos == NULL)
    	return;

    for(uint i=0; i<numGridPoints; i++)
    {
        psi[i]      = 0.0;
        temp[i]     = 0.0;
        relError[i] = 0.0;
    }
}

bool SphericalGrid::isElliptical() const
{
	return false;
}

Vec3D SphericalGrid::getPos(uint index, bool computationalCoords) const
{
	if (index < numGridPoints)
		return pos[index];
	else
	{
		printf("ERROR! SphericalGrid::getPos: index (%u) out of bounds (%u)\n", index, numGridPoints);
		return Vec3D(0.0, 0.0, 0.0);
	}
}

Vec3D* SphericalGrid::getPosArray() const
{
	return pos;
}

Vec3D* SphericalGrid::getBArray() const
{
	return B;
}

Vec3D SphericalGrid::getB(uint index, bool ellipticCoords) const
{
	if (index < numGridPoints)
		return B[index];
	else
	{
		printf("ERROR! SphericalGrid::getB: index (%u) out of bounds (%u)\n", index, numGridPoints);
		return Vec3D(0.0, 0.0, 0.0);
	}
}

void SphericalGrid::setB(uint index, Vec3D value, bool ellipticCoords)
{
	if (index < numGridPoints)
		B[index] = value;
	else
		printf("ERROR! SphericalGrid::setB: index (%u) out of bounds (%u)\n", index, numGridPoints);
}

hcFloat* SphericalGrid::getPsiArray() const
{
	return this->psi;
}

hcFloat SphericalGrid::getPsi(uint index) const
{
	if (index < numGridPoints)
		return psi[index];
	else
	{
		printf("ERROR! SphericalGrid::getPsi: index (%u) out of bounds (%u)\n", index, numGridPoints);
		return 0.0;
	}
}

void SphericalGrid::setPsi(uint index, hcFloat value)
{
	if (index < numGridPoints)
		psi[index] = value;
	else
		printf("ERROR! SphericalGrid::setPsi: index (%u) out of bounds (%u)\n", index, numGridPoints);
}

hcFloat* SphericalGrid::getRelErrorArray() const
{
	return this->relError;
}

hcFloat SphericalGrid::getRelError(uint index) const
{
	if (index < numGridPoints)
		return relError[index];
	else
	{
		printf("ERROR! SphericalGrid::getRelError: index (%u) out of bounds (%u)\n", index, numGridPoints);
		return 0.0;
	}
}

void SphericalGrid::setRelError(uint index, hcFloat value)
{
	if (index < numGridPoints)
		relError[index] = value;
	else
		printf("ERROR! SphericalGrid::setRelError: index (%u) out of bounds (%u)\n", index, numGridPoints);
}

hcFloat* SphericalGrid::getTempArray() const
{
	return this->temp;
}

hcFloat SphericalGrid::getTemp(uint index) const
{
	if (index < numGridPoints)
	    return temp[index];
	else
	{
		printf("ERROR! SphericalGrid::getTemp: index (%u) out of bounds (%u)\n", index, numGridPoints);
		return 0.0;
	}
}

void SphericalGrid::setTemp(uint index, hcFloat value)
{
	if (index < numGridPoints)
		temp[index] = value;
	else
		printf("ERROR! SphericalGrid::setTemp: index (%u) out of bounds (%u)\n", index, numGridPoints);
}

uint SphericalGrid::getIndexPhiPlus(uint ind)
{
    uint i,j,k;
    getIndicesFromIndex(ind, i, j, k);
    k = (k + 1) % numPhi;
    return getIndex(i,j,k);
}

void SphericalGrid::printGridIndicesFromIndex(uint ind)
{
    uint i,j,k;
    uint ind1;

    i		= ind%numR;
    ind1 	= (ind-i)/numR;

    j 		= ind1%numTheta;
    k 		= (ind1-j)/numTheta;

    printf("ind: %u, i: %u, j: %u, k: %u\n", ind, i,j,k);
}


bool SphericalGrid::exportEntireGrid(const char *filename)
{
    createFolderStructureTo(filename);

    std::ofstream solutionFile(filename, std::ios::out | std::ios::binary);

    if(sizeofFloat != sizeof(hcFloat))
    {
    	cerr << __FILE__ << ":" << __LINE__ << " Grid was loaded with different floating point precision than this program was compiled for.\n";
    	cerr << "Storing it again would destroy accuracy claim of this class, so it is forbidden\n";
    	return false;
    }

    solutionFile.write(reinterpret_cast<char*>(&sizeofFloat),sizeof(uint));
    solutionFile.write(reinterpret_cast<char*>(pos),      	numGridPoints*sizeof(Vec3D));
    solutionFile.write(reinterpret_cast<char*>(B),      	numGridPoints*sizeof(Vec3D));
    solutionFile.close();

    return true;
}

bool SphericalGrid::importEntireGrid(const char *filename)
{
    if(!checkFileEx(filename, "SphericalGrid::importEntireGrid"))
        return false;

    std::ifstream solutionFile;
    solutionFile.open(filename, std::ios::in | std::ios::binary);

    solutionFile.read(reinterpret_cast<char*>(&sizeofFloat), sizeof(uint));
    if(sizeofFloat != 4 && sizeofFloat != 8 && sizeofFloat != 16)			// just in case file has been stored using old scheme
    {
    	solutionFile.close();
        solutionFile.open(filename, std::ios::in | std::ios::binary);
        sizeofFloat = 4;
    }

    if(sizeof(hcFloat)==sizeofFloat)
    {
		solutionFile.read(reinterpret_cast<char*>(pos),       	numGridPoints*sizeof(Vec3D));
		solutionFile.read(reinterpret_cast<char*>(B),       	numGridPoints*sizeof(Vec3D));
    }
    else
    {
    	if(sizeofFloat==4)
    	{
    		Vec<3,float> *posT 	= new Vec<3,float> [numGridPoints];
    		Vec<3,float> *BT	= new Vec<3,float> [numGridPoints];

    		solutionFile.read(reinterpret_cast<char*>(posT),       	numGridPoints*sizeof(Vec<3,float>));
			solutionFile.read(reinterpret_cast<char*>(BT),       	numGridPoints*sizeof(Vec<3,float>));

			for(uint i=0; i<numR; ++i)
				for(uint j=0; j<numTheta; ++j)
					for(uint k=0; k<numPhi; ++k)
					{
						uint index		= getIndex(i,j,k);
						pos[index]		= posT[index];
						B[index]		= BT[index];
					}

			delete [] posT;
			delete [] BT;
    	}
    	else if(sizeofFloat==8)
    	{
    		Vec<3,double> *posT = new Vec<3,double> [numGridPoints];
			Vec<3,double> *BT	= new Vec<3,double> [numGridPoints];

			solutionFile.read(reinterpret_cast<char*>(posT),       	numGridPoints*sizeof(Vec<3,double>));
			solutionFile.read(reinterpret_cast<char*>(BT),       	numGridPoints*sizeof(Vec<3,double>));

			for(uint i=0; i<numR; ++i)
				for(uint j=0; j<numTheta; ++j)
					for(uint k=0; k<numPhi; ++k)
					{
						uint index		= getIndex(i,j,k);
						pos[index]		= posT[index];
						B[index]		= BT[index];
					}

			delete [] posT;
			delete [] BT;
    	}
    	else if(sizeofFloat==16)
    	{
    		Vec<3,long double> *posT= new Vec<3,long double> [numGridPoints];
			Vec<3,long double> *BT	= new Vec<3,long double> [numGridPoints];

			solutionFile.read(reinterpret_cast<char*>(posT),       	numGridPoints*sizeof(Vec<3,long double>));
			solutionFile.read(reinterpret_cast<char*>(BT),       	numGridPoints*sizeof(Vec<3,long double>));

			for(uint i=0; i<numR; ++i)
				for(uint j=0; j<numTheta; ++j)
					for(uint k=0; k<numPhi; ++k)
					{
						uint index		= getIndex(i,j,k);
						pos[index]		= posT[index];
						B[index]		= BT[index];
					}

			delete [] posT;
			delete [] BT;
    	}
    	else
    	{
    		cerr << __FILE__ << ":" << __LINE__ << " importing grid with sizeofFloat=" << sizeofFloat << " is not supported\n";
    	}
    }
    solutionFile.close();

    return true;
}

void SphericalGrid::diff(const SphericalGrid &other)
{
	if(this->numR != other.numR || this->numTheta != other.numTheta || this->numPhi != other.numPhi)
	{
		printf("ERROR! SphericalGrid::diff: dimenstion do not match!\n\t\tthis\tother\nnumR\t%u\t%u\nnumTheta:\t%u\t%u\nnumPhi:\t%u\t%u\n", this->numR, other.numR, this->numTheta, other.numTheta, this->numPhi, other.numPhi);
		return;
	}

	for(uint r=0; r<numR; ++r)
		for(uint t=0; t<numTheta; ++t)
			for(uint p=0; p<numPhi; ++p)
			{
				uint ind	= getIndex(r, t, p);
				psi[ind]	= fabs(psi[ind] - other.psi[ind]);
			}
}

hcFloat SphericalGrid::maxSpacingR(uint &index)
{
	uint maxDist 	= 0.0;
	bool ell		= isElliptical();

	for(uint r=0; r<numR; ++r)
		for(uint theta=0; theta<numTheta; ++theta)
			for(uint phi=0; phi<numPhi; ++phi)
			{
				uint ind 	= getIndex(r, theta, phi);
				uint ind_p	= r < numR-1 	? getIndex(r+1, theta, phi) : getIndex(r, theta, phi);
				uint ind_m	= r > 0			? getIndex(r-1, theta, phi) : getIndex(r, theta, phi);

				Vec3D pos_h	= getPos(ind, ell).convCoordSpher2Cart();
				Vec3D pos_p = getPos(ind_p, ell).convCoordSpher2Cart();
				Vec3D pos_m = getPos(ind_m, ell).convCoordSpher2Cart();

				hcFloat dist_p	= pos_p.dist(pos_h);
				hcFloat dist_m	= pos_m.dist(pos_h);

				if(dist_p > maxDist)
				{
					maxDist = dist_p;
					index	= ind;
				}

				if(dist_m > maxDist)
				{
					maxDist	= dist_m;
					index	= ind;
				}
			}

	return maxDist;
}

hcFloat SphericalGrid::maxSpacingTheta(uint &index)
{
	uint maxDist 	= 0.0;
	bool ell		= isElliptical();

	for(uint r=0; r<numR; ++r)
		for(uint theta=0; theta<numTheta; ++theta)
			for(uint phi=0; phi<numPhi; ++phi)
			{
				uint ind 	= getIndex(r, theta, phi);
				uint ind_p	= theta < numTheta - 1	? getIndex(r, theta+1, phi) : getIndex(r, theta, phi);
				uint ind_m	= theta > 0				? getIndex(r, theta-1, phi) : getIndex(r, theta, phi);

				Vec3D pos_h	= getPos(ind, ell).convCoordSpher2Cart();
				Vec3D pos_p = getPos(ind_p, ell).convCoordSpher2Cart();
				Vec3D pos_m = getPos(ind_m, ell).convCoordSpher2Cart();

				hcFloat dist_p	= pos_p.dist(pos_h);
				hcFloat dist_m	= pos_m.dist(pos_h);

				if(dist_p > maxDist)
				{
					maxDist = dist_p;
					index	= ind;
				}

				if(dist_m > maxDist)
				{
					maxDist	= dist_m;
					index	= ind;
				}
			}

	return maxDist;
}

hcFloat SphericalGrid::maxSpacingPhi(uint &index)
{
	uint maxDist 	= 0.0;
	bool ell		= isElliptical();

	for(uint r=0; r<numR; ++r)
		for(uint theta=0; theta<numTheta; ++theta)
			for(uint phi=0; phi<numPhi; ++phi)
			{
				uint ind 	= getIndex(r, theta, phi);
				uint ind_p	= phi < numPhi - 1	? getIndex(r, theta, phi+1) : getIndex(r, theta, 0);
				uint ind_m	= phi > 0			? getIndex(r, theta, phi-1) : getIndex(r, theta, numPhi-1);

				Vec3D pos_h	= getPos(ind, ell).convCoordSpher2Cart();
				Vec3D pos_p = getPos(ind_p, ell).convCoordSpher2Cart();
				Vec3D pos_m = getPos(ind_m, ell).convCoordSpher2Cart();

				hcFloat dist_p	= pos_p.dist(pos_h);
				hcFloat dist_m	= pos_m.dist(pos_h);

				if(dist_p > maxDist)
				{
					maxDist = dist_p;
					index	= ind;
				}

				if(dist_m > maxDist)
				{
					maxDist	= dist_m;
					index	= ind;
				}
			}

	return maxDist;
}

bool SphericalGrid::evaluateCSSS(const char *filename)
{
	CSSS_magfield csss(filename);

	for(uint i=0; i<numR; ++i)
		for(uint j=0; j<numTheta; ++j)
			for(uint k=0; k<numPhi; ++k)
			{
				uint ind = getIndex(i, j, k);
				csss.eval(pos[ind], B[ind]);
			}

	return true;
}

void SphericalGrid::dump() const
{
	printf("Dumping SphericalGrid:--------------------------------\n");
    printf("numR:\t\t%u\n", numR);
    printf("numTheta:\t%u\n", numTheta);
    printf("numPhi:\t\t%u\n", numPhi);

    printf("numGridPoints:\t%u\n", numGridPoints);

    printf("sinLatGrid:\t%u\n", sinLatGrid);
    printf("maxSinLat:\t%E\n", maxSinLat);
    printf("minSinLat:\t%E\n", minSinLat);

    printf("lowerR:\t%E\n", lowerR);
    printf("upperR:\t%E\n", upperR);
    printf("geometricFactor:\t%E\n", geometricFactor);
#ifdef CUDA
    printf("OnDevice: %u\n", onDevice);
#endif
}
void SphericalGrid::dumpCoords(uint fixed, bool ellipticCoords) const
{
	printf("----- Dumping coords of Spherical Grid in spherical coords and converted to cartesian\n");
	printf("ind\ti\tj\tk\tr\ttheta\tphi\tx\ty\tz\n");
	for(uint i = (fixed != 4 ? 0 : numR - 1); i < (fixed != 1 ? numR : 1); ++i)
	{
		for(uint j = 0; j < (fixed != 2 ? numTheta : 1); ++j)
		{
			for(uint k = 0; k < (fixed != 3 ? numPhi : 1); ++k)
			{
				uint ind 	= k * numR * numTheta + j * numR + i;
				Vec3D Pos 	= getPos(ind, false);
				Vec3D cPos 	= Pos.convCoordSpher2Cart();
				printf("%u\t%u\t%u\t%u\t%E\t%E\t%E\t%E\t%E\t%E\n", ind, i, j, k, Pos[0], Pos[1], Pos[2], cPos[0], cPos[1], cPos[2]);
			}
		}
	}
}


static void getOptimalGridParameters_hi(uint numR, hcFloat lowerR, hcFloat upperR, hcFloat &geometricFactor, uint &numT, uint &numP)
{
	hcFloat q		= 1.0000;
	hcFloat delta	= 0.000001;
	bool eval		= false;

	do																	// this is the condition for the same ratio dr/dtheta at photosphere
	{																	// and source surface
		q 			+= delta;
		hcFloat lhs	= (q-1.0)/(pow(q, numR-1)-pow(q, numR-2));
		hcFloat rhs	= lowerR/upperR;
		eval		= lhs > rhs;
	}
	while(eval);

	geometricFactor = q;
	numT			= PI*lowerR/(upperR-lowerR)*(1-pow(q, numR-1))/(1-q)+1;		// this is the condition for (roughly) the same spacing in theta- and r-direction at the photosphere
	numP 			= 2*numT;
}

static void getOptimalGridParameters_lo(uint numR, hcFloat lowerR, hcFloat upperR, hcFloat &geometricFactor, uint &numT, uint &numP)
{
	hcFloat q		= 1.0000;
	hcFloat delta	= 0.000001;
	bool eval		= false;

	do																	// this is the condition for the same ratio dr/dtheta at photosphere
	{																	// and source surface
		q 			+= delta;
		hcFloat lhs	= (q-1.0)/(pow(q, numR-1)-pow(q, numR-2));
		hcFloat rhs	= lowerR/upperR;
		eval		= lhs > rhs;
	}
	while(eval);

	geometricFactor = q;
	numT	 		= 2.5*numR;
	numP			= 5*numR;
}

void SphericalGrid::getOptimalGridParameters(uint numR, hcFloat lowerR, hcFloat upperR, hcFloat &geometricFactor,
		uint &numT, uint &numP, bool hi)
{
	if(hi)	getOptimalGridParameters_hi(numR, lowerR, upperR, geometricFactor, numT, numP);
	else	getOptimalGridParameters_lo(numR, lowerR, upperR, geometricFactor, numT, numP);
}
