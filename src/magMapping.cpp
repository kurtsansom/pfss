#include "src/magMapping.h"
#include "src/ellipticalGrid.h"
#include "src/filenames.h"

#include "engine/hcImage.h"
#include "engine/hcTime.h"

#include <fstream>
#include <iostream>
#include <fstream>
#include <unistd.h>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace std;

#ifdef GUI
#include "engine/math/hcCoordTransform.h"
#endif

bool convertMagMap2bmp(const char *filename)
{
	if(!doesFileExist(filename))
	{
		printf("readMagneticHeightMap: file '%s' does not exist!\n", filename);
		return false;
	}

	char line[10000];

	regex re2("[0-9]\\.[0-9]{6}E[-+][0-9]{2} ([0-9]\\.[0-9]{6}E[-+][0-9]{2}) [NPCI].*", regex::icase);
	regex re3("[0-9]\\.[0-9]{6}E[-+][0-9]{2} [0-9]\\.[0-9]{6}E[-+][0-9]{2} ([NPCI]).*", regex::icase);

	ifstream file(filename, ios_base::in);

	uint n 			= 0;
	uint numTheta	= 0;
	uint numPhi		= 0;

	while(file.getline(line, 10000))
	{
		if(line[0] == '#' || strlen(line) == 0)
			continue;

		cmatch what;
		regex_search(line, what, re2);
		float second = lexical_cast<float>(what[1]);

		if(second == 0.0 && n>0 && numPhi == 0)
			numPhi = n;

		++n;
	}

	numTheta = n / numPhi;

	file.close();
	file.open(filename, ios_base::in);

	hcImageRGBA image(numPhi, numTheta);

	uint x = 0;
	uint y = 0;

	while(file.getline(line, 10000))
	{
		if(line[0] == '#' || line[0] =='\0')
			continue;

		cmatch what;
		char id;

		regex_search(line, what, re3);
		id = lexical_cast<char>(what[1]);

		if(id=='I')			image(x, y) = char2RGBA8(0,  0,  0,  255);
		else if(id=='C')	image(x, y) = char2RGBA8(0,  0,  255,255);
		else if(id=='N')	image(x, y) = char2RGBA8(255,0,  0,  255);
		else if(id=='P')	image(x, y) = char2RGBA8(0,  255,0,  255);

		++x;

		if(x==numPhi)
		{
			x=0;
			++y;
		}
	}

	char bmpFilename[1000];
	sprintf(bmpFilename, "%s.bmp", filename);
	image.save(bmpFilename);

	return true;
}

#ifdef GUI
void drawEcliptic(hcImageRGBA &img, hcFloat maxSinLat, bool sinLatFormat, uint crNum)
{
	//uint color						= char2RGBA8(0, 0, 0, 255);
	uint cr							= 0;
	uint cg							= 0;
	uint cb							= 0;
	hcFloat d						= 1 * deg2rad;
	hcFloat dSinLat					= 2*maxSinLat 		/ (img.height - 1);
	hcFloat dLat					= 2*asin(maxSinLat)	/ (img.height - 1);

	hcDate time;

	for(uint x=0; x<img.width; ++x)
	{
		hcFloat phi			= x * 2 * PI / img.width;
		hcFloat carrTime	= 360 * crNum - phi*rad2deg;

		time.setFromCarringtonTime(carrTime);

		for(uint y=0; y<img.height; ++y)
		{
			hcFloat theta		= sinLatFormat ? PI/2 - asin(maxSinLat - y * dSinLat) : PI/2-asin(maxSinLat) + y*dLat;
			Vec2D posHGC		= Vec2D(phi, (hcFloat)(PI/2.0 - theta));
			Vec2D posECL		= hcTFcelHGC2ECL(posHGC, time);

			if(fabs(posECL[1]) < d)
			{
				hcFloat factor 	= fabs(fabs(posECL[1]) - d)/d;
				uint mask_r		= 255 << 0;
				uint mask_g		= 255 << 8;
				uint mask_b		= 255 << 16;
				uint color_img	= img(x, img.height-y-1);
				uint img_r		= (mask_r & color_img) >> 0;
				uint img_g		= (mask_g & color_img) >> 8;
				uint img_b		= (mask_b & color_img) >> 16;
				uint r			= (factor * cr + (1-factor) * img_r);
				uint g			= (factor * cg + (1-factor) * img_g);
				uint b			= (factor * cb + (1-factor) * img_b);
				img(x, img.height-y-1) = char2RGBA8(r,g,b,255);
			}
		}
	}
}
#endif

MagMapping::MagMapping()
{
    initNULL();
}

MagMapping::MagMapping(const MagMapping &other)
{
    initNULL();
    *this = other;
}

MagMapping::MagMapping(const PFSSsolutionInfo &info, bool sinLatFormat, bool compCoords, hcFloat maxSinLat, uint numTheta, uint numPhi, hcFloat height)
{
    initNULL();
    init(info, sinLatFormat, compCoords, maxSinLat, numTheta, numPhi, height);
}

MagMapping::~MagMapping()
{
    clear();
}

MagMapping &MagMapping::operator=(const MagMapping &other)
{
    if(this == &other)
        return *this;

    clear();
    init(other.info, other.sinLatFormat, other.compCoords, other.maxSinLat, other.numTheta, other.numPhi, other.getHeight());

    for(uint j=0; j<numTheta; ++j)
        for(uint k=0; k<numPhi; ++k)
            *maglines[index(j,k)] = *other.maglines[index(j,k)];

    return *this;
}

Magline &MagMapping::operator()(uint indTheta, uint indPhi)
{
    if(maglines == NULL)
    {
        cerr << __FILE__ << "/" << __LINE__ << " maglines not initialized!\n";
        exit(1);
    }

    if(indTheta >= numTheta || indPhi >= numPhi)
    {
        printf("ERROR! MagMapping::operator(): index out of bounds! indTheta: %u / %u, indPhi: %u / %u\n", indTheta, numTheta, indPhi, numPhi);
        return *maglines[0];
    }

    return *maglines[index(indTheta, indPhi)];
}


void MagMapping::initNULL()
{
	info			= PFSSsolutionInfo();
    maglines    	= NULL;
    coords      	= NULL;
    height			= 0.0;
    compCoords		= true;
    sinLatFormat 	= true;
    numTheta    	= 0;
    numPhi      	= 0;
}

void MagMapping::clear()
{
    if(maglines != NULL)
    {
        for(uint j=0; j<numTheta; ++j)
            for(uint k=0; k<numPhi; ++k)
                delete maglines[index(j,k)];
    }
    delete [] maglines;
    delete [] coords;

    initNULL();
}

void MagMapping::init(const PFSSsolutionInfo &info, bool sinLatFormat, bool compCoords, hcFloat maxSinLat, uint numTheta, uint numPhi, hcFloat height)
{
    clear();
    this->info			= info;
    this->sinLatFormat	= sinLatFormat;
    this->maxSinLat		= maxSinLat;
    //this->dLat			= ?				// TODO: non sine-latitude?
    hcFloat dLat		= 0.0;
    this->numTheta  	= numTheta;
    this->numPhi    	= numPhi;
    this->height		= height;
    this->compCoords	= compCoords;

    maglines    		= new Magline*[numTheta * numPhi];
    coords      		= new Vec3D[   numTheta * numPhi];

    hcFloat dSinLat		= 2*maxSinLat/(numTheta-1);

    for(uint j=0; j<numTheta; ++j)
        for(uint k=0; k<numPhi; ++k)
        {
        	hcFloat theta		= sinLatFormat ? PI/2 - asin(maxSinLat - j * dSinLat) : PI/2-asin(maxSinLat) + j*dLat;
			hcFloat phi         = k * 2 * PI / numPhi;
			coords[index(j, k)] = Vec3D(height, theta, phi);
            maglines[index(j,k)]= new Magline();
        }
}

uint MagMapping::index(uint indTheta, uint indPhi)
{
	return indTheta * numPhi + indPhi;
}

bool MagMapping::isComputed() const
{
	return maglines != NULL;
}

bool MagMapping::checkComputed() const
{
	if(!isComputed())
	{
		printErrMess(__FILE__, __LINE__, "magmapping not computed");
		return false;
	}
	return true;
}

// TODO: insert in hcTools?
bool MagMapping::checkExportedFile(const string &filename) const
{
	if(doesFileExist(filename))
	{
		printStdOutMess(__FILE__, __LINE__, "file '" + filename + "' does already exist, do not export");
		return false;
	}

	if(!createFolderStructureTo(filename.c_str()))
	{
		printErrMess(__FILE__, __LINE__, "file '" + filename + "' cannot be created");
		return false;
	}

	return true;
}

uint MagMapping::getNumTheta()
{
    return numTheta;
}

uint MagMapping::getNumPhi()
{
    return numPhi;
}

Vec3D MagMapping::getCoords(uint j, uint k)
{
	return coords[index(j,k)];
}

bool MagMapping::setCoords(uint j, uint k, const Vec3D &coord)
{
	if(j>=numTheta || k>=numPhi)
	{
		printErrMess(__FILE__, __LINE__, "pixel index out of bounds, theta: " + to_string(j) + "/" + to_string(numTheta) + " phi: " + to_string(k) + "/" + to_string(numPhi));
		return false;
	}

	coords[index(j,k)] = coord;
	return true;
}

hcFloat MagMapping::getHeight() const
{
	return height;
}

bool MagMapping::exportBinary(const string &filename)
{
    if(!createFolderStructureTo(filename.data()) || boost::filesystem::is_directory(filename.data()))
    {
        printErrMess(__FILE__, __LINE__, "file '" + filename + "' cannot not be created");
        return false;
    }

    bool retval 		= true;
    uint sizeofFloat 	= sizeof(hcFloat);
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::binary);

    file.write(reinterpret_cast<char*>(&sizeofFloat),	sizeof(uint));
    retval &= info.exportBinary(file);
    file.write(reinterpret_cast<char*>(&sinLatFormat),	sizeof(bool));
    file.write(reinterpret_cast<char*>(&maxSinLat),		sizeofFloat);
    file.write(reinterpret_cast<char*>(&compCoords),	sizeof(bool));
    file.write(reinterpret_cast<char*>(&numTheta),		sizeof(uint));
    file.write(reinterpret_cast<char*>(&numPhi),		sizeof(uint));
    file.write(reinterpret_cast<char*>(&height),		sizeofFloat);

    for(uint i=0; i<numTheta*numPhi; ++i)
    {
    	file.write(reinterpret_cast<char*>(&coords[i](0)),sizeofFloat);
        file.write(reinterpret_cast<char*>(&coords[i](1)),sizeofFloat);
        file.write(reinterpret_cast<char*>(&coords[i](2)),sizeofFloat);
    }

    for(uint j=0; j<numTheta; ++j)
        for(uint k=0; k<numPhi; ++k)
        	retval &= operator ()(j, k).exportBinary(file);
    file.close();

    return retval;
}


bool MagMapping::importBinary(const string &filename)
{
	if(!doesFileExist(filename)) return false;

    PFSSsolutionInfo inf;
    ifstream stream;
	stream.open(filename, std::ios::out | std::ios::binary);
    bool retval 		= true;
    uint sizeoffloat	= 0;
	stream.read(reinterpret_cast<char*>(&sizeoffloat),	sizeof(uint));
	if(sizeoffloat != 4 && sizeoffloat != 8) return false;
	char *tempFloat		= sizeoffloat == 4 ? reinterpret_cast<char*>(new float()): reinterpret_cast<char*>(new double());
    retval &= inf.importBinary(stream);
    stream.read(reinterpret_cast<char*>(&sinLatFormat),	sizeof(bool));
    stream.read(reinterpret_cast<char*>(tempFloat),		sizeoffloat); maxSinLat = sizeoffloat==4 ? *(reinterpret_cast<float*>(tempFloat)) : *(reinterpret_cast<double*>(tempFloat));
    stream.read(reinterpret_cast<char*>(&compCoords),	sizeof(bool));
    stream.read(reinterpret_cast<char*>(&numTheta),		sizeof(uint));
    stream.read(reinterpret_cast<char*>(&numPhi),		sizeof(uint));
    stream.read(reinterpret_cast<char*>(tempFloat),		sizeoffloat); height 	= sizeoffloat==4 ? *(reinterpret_cast<float*>(tempFloat)) : *(reinterpret_cast<double*>(tempFloat));

    init(inf, sinLatFormat, compCoords, maxSinLat, numTheta, numPhi, height);

    for(uint i=0; i<numTheta*numPhi; ++i)
    {
        stream.read(reinterpret_cast<char*>(tempFloat),	sizeoffloat); coords[i](0) = sizeoffloat==4 ? *(reinterpret_cast<float*>(tempFloat)) : *(reinterpret_cast<double*>(tempFloat));
        stream.read(reinterpret_cast<char*>(tempFloat),	sizeoffloat); coords[i](1) = sizeoffloat==4 ? *(reinterpret_cast<float*>(tempFloat)) : *(reinterpret_cast<double*>(tempFloat));
        stream.read(reinterpret_cast<char*>(tempFloat),	sizeoffloat); coords[i](2) = sizeoffloat==4 ? *(reinterpret_cast<float*>(tempFloat)) : *(reinterpret_cast<double*>(tempFloat));
    }

    for(uint j=0; j<numTheta; ++j)
        for(uint k=0; k<numPhi; ++k)
        	retval &= operator ()(j, k).importBinary(stream);
    stream.close();

    delete tempFloat;
    return retval;
}

/*! compares photospheric footpoints of back mapped magnetic field lines
 *
 * 	This function only makes sense if it is a backmapping from source surface down to the photosphere
 *
 */
hcFloat MagMapping::diffFootpoints(MagMapping &other)
{
	hcFloat retval = 0.0;

	if(this->numPhi != other.numPhi || this->numTheta != other.numTheta)
	{
		printf("ERROR! MagMapping::diffFootpoints: Dimensions do not match!\nthis.numTheta: %u, other.numTheta: %u\nthis.numPhi: %u, other.numPhi:%u\n", this->numTheta, other.numTheta, this->numPhi, other.numPhi);
		return retval;
	}

	uint numPolChange	= 0;

	for(uint t=0; t<numTheta; ++t)
		for(uint p=0; p<numPhi; ++p)
		{
			uint ind 		= index(t, p);
			Magline &magl0 	= this->operator()(t, p);
			Magline &magl1 	= other.operator()(t, p);

			if(!magl0.valid || !magl1.valid)		continue;
			if(magl0.polarity != magl1.polarity)	++numPolChange;

			hcFloat theta0 	= this->coords[ind][0];
			hcFloat theta1 	= other.coords[ind][0];
			hcFloat phi0 	= this->coords[ind][1];
			hcFloat phi1 	= other.coords[ind][1];

			if(theta0 != theta1 || phi0 != phi1)
			{
				printf("ERROR MagMapping::diffFootpoints:\nt: %u, p: %u, theta0: %E, theta1: %E, phi0: %E, phi1: %E\n", t, p, theta0, theta1, phi0, phi1);
				return 0.0;
			}

			hcFloat t0		= magl0.posdata[0][1];
			hcFloat p0		= magl0.posdata[0][2];

			hcFloat t1 		= magl1.posdata[0][1];
			hcFloat p1 		= magl1.posdata[0][2];

			hcFloat dt		= fabs(t0-t1);
			hcFloat dp		= min(fabs(p0-p1), fabs(p0-p1-2*PI));

			hcFloat dist 	= dt*dt + dp*dp;

			retval += dist;
		}

	retval /= (numTheta * numPhi);

	return retval;
}

hcFloat MagMapping::getOpenFlux()
{
	hcFloat retval 		= 0.0;
	hcFloat maxSinLat	= sin(PI/2.0 - coords[index(0,0)][0]);

	for(uint t=0; t<numTheta; ++t)
		for(uint p=0; p<numPhi; ++p)
		{
			Magline &magl 	= this->operator()(t, p);
			if(magl.closed || !magl.valid) continue;										// consider only maglines that are open

			Vec3D &pos				= coords[index(t,p)];
			hcFloat dSinTheta		= 2*maxSinLat/(numTheta-1);
			hcFloat upperPixBound	= PI/2.0 - asin(maxSinLat - (t-1.0/2.0) * dSinTheta);
			hcFloat lowerPixBound	= PI/2.0 - asin(maxSinLat - (t+1.0/2.0) * dSinTheta);

			hcFloat dTheta			= lowerPixBound - upperPixBound;
			hcFloat dPhi			= 2 * PI / numPhi;

			hcFloat area			= height*height * dTheta * dPhi;

			Vec3D *posArr			= new Vec3D[MAGLINE_NUM_POSITIONS];
			Vec3D *magArr			= new Vec3D[MAGLINE_NUM_POSITIONS];
			uint num				= magl.getAllValuesAtHeight(height, posArr, magArr);

			if(num == 0)
			{
				printStdOutMess(__FILE__, __LINE__, "no intersection points found, his message should not be possible");
				continue;
			}

			uint k=0;
			hcFloat dT,dP;
			if(num > 1)
			{
				while(k<num)
				{
					dT		= (pos[0] - posArr[k][1]);
					dP		= min(fabs(pos[1]-posArr[k][2]),fabs(2*PI-fabs(pos[1]-posArr[k][2])));
					if(dT < 1E-4 && dP < 1E-4)	break;
					++k;
				}
			}
			else
			{
				dT		= (pos[0] - posArr[k][1]);
				dP		= min(fabs(pos[1]-posArr[k][2]),fabs(2*PI-fabs(pos[1]-posArr[k][2])));
			}

			if(k==num-1 && (dT > 1E-4 || dP > 1E-4))
			{
				cout.precision(4);cout << fixed;
				cout << "position does not coincide with array:\n";
				cout << "pos: " << pos[0] 		<< "/" << pos[1] 		<< "\n";
				cout << "arr: " << posArr[0][1] << "/" << posArr[0][2]	<< "\n\n";
				fflush(stdout);
				continue;
			}

			Vec3D &B 		= magArr[0];
			hcFloat flux	= B[0] * area;
			retval 			+= fabs(flux);

			delete [] posArr;
			delete [] magArr;
		}

	return retval;
}

/*! computed and stored entirely in computatinal coordinates(vectors, transformation to physical coordinates/vectors is
 *  done afterwards when re-importing data
 *
 * \param solver        gives the environment in which maglines shall be created
 * \param height        the height above photosphere where equidistant "pixels" are evaluated
 * \param numTheta      the number of pixels in theta direction (if 0, the entire mapping is computed in the native grid of solver)
 * \param numPhi        the number of pixels in phi direction
 * \param maxSinLat     the highest point in mapping
 * \param minSinLat     the lowest point in mapping
 * \param sinLatGrid    evenly spaced in sinLat or in lat?
 * \param compCoords	height given in computational coordinates (true) or world coordinates (false)
 * \return              success
 */
bool MagMapping::createAtHeight(const PFSSsolutionInfo &info, LaplaceSolver &solver, hcFloat height, uint numThetaIn, uint numPhiIn,
                                hcFloat maxSinLatIn, bool sinLatFormatIn, bool compCoords)
{
	printStdOutMess(__FILE__, __LINE__, "started mapping with NUMTHREADS=" + to_string(NUMTHREADS));

	if(!solver.solutionComputed)
	{
		printErrMess(__FILE__, __LINE__, "solver has not been used for computation");
		return false;
	}

	hcDate timeStart, timeStop;

	uint numTheta		= numThetaIn == 0 ? solver.grid->numTheta 	: numThetaIn;
	uint numPhi			= numThetaIn == 0 ? solver.grid->numPhi		: numPhiIn;
	hcFloat maxSinLat	= numThetaIn == 0 ? solver.grid->maxSinLat	: maxSinLatIn;
	bool sinLatGrid		= numThetaIn == 0 ? solver.grid->sinLatGrid	: sinLatFormatIn;
	bool retval 		= false;

#if NUMTHREADS > 1
		retval = createAtHeight_MP(info, solver, height, numTheta, numPhi, maxSinLat, sinLatGrid, compCoords);
#else
		retval = createAtHeight_SP(info, solver, height, numTheta, numPhi, maxSinLat, sinLatGrid, compCoords);
#endif

	if(solver.grid->isElliptical())
	{
		EllipticalGrid* eGrid = (EllipticalGrid*)solver.grid;
		eGrid->convertMagMapping(*this);
	}

	timeStop.setFromSystemTime();
	string message 	= "took " + toStr((timeStop-timeStart)/hcDate::facSec) + " s to map with NUMTHREADS=" + to_string(NUMTHREADS);
	printStdOutMess(__FILE__, __LINE__ , message);
	return retval;
}

bool MagMapping::createAtHeight_SP(const PFSSsolutionInfo &info, LaplaceSolver &solver, hcFloat height, uint numTheta, uint numPhi,
                                hcFloat maxSinLat, bool sinLatFormat, bool compCoords)
{
	hcFloat dSinLat = 2*maxSinLat		/ (numTheta - 1);
	hcFloat dLat	= 2*asin(maxSinLat)	/ (numTheta - 1);

    init(info, sinLatFormat, compCoords, maxSinLat, numTheta, numPhi, height);

    for(uint j=0; j<numTheta; ++j)
        for(uint k=0; k<numPhi; ++k)
        {
            Magline &magline	= operator()(j, k);
            hcFloat theta		= sinLatFormat ? PI/2 - asin(maxSinLat - j * dSinLat) : PI/2-asin(maxSinLat) + j*dLat;
            hcFloat phi         =  k * 2 * PI / numPhi;
            Vec3D pos			= Vec3D(height, theta, phi);
            Vec3D posSpher		= solver.grid->isElliptical() ? pos.convCoordSpher2Ell(*(dynamic_cast<EllipticalGrid*>(solver.grid))) : pos;
            Vec3D posComp		= compCoords ? pos : posSpher;
            //coords[index(j, k)] = Vec2D(posComp[1], posComp[2]);
            coords[index(j, k)] = posComp;
            magline.createMaglineThroughPos(*solver.grid, posComp, false);
        }

    return true;
}

/*! multi-threaded version of magnetic field line creation
 */
bool MagMapping::createAtHeight_MP(const PFSSsolutionInfo &info, LaplaceSolver &solver, hcFloat height, uint numTheta, uint numPhi,
                                hcFloat maxSinLat, bool sinLatFormat, bool compCoords)
{
	init(info, sinLatFormat, compCoords, maxSinLat, numTheta, numPhi, height);

	uint numMaglines				= numTheta * numPhi;
	bool *workedUpon 				= new bool[numMaglines];
	Vec3D *posComp					= new Vec3D[numMaglines];

	SphericalGrid *grid				= solver.grid;
	EllipticalGrid *egrid			= dynamic_cast<EllipticalGrid*>(grid);
	hcFloat dSinLat					= 2*maxSinLat 		/ (numTheta - 1);
	hcFloat dLat					= 2*asin(maxSinLat)	/ (numTheta - 1);

	for(uint i=0;i<numMaglines;++i)
	{
		uint nPhi			= i 		 	% numPhi;
		uint nTheta			= (i - nPhi) 	/ numPhi;
		hcFloat theta		= sinLatFormat ? PI/2 - asin(maxSinLat - nTheta * dSinLat) : PI/2-asin(maxSinLat) + nTheta*dLat;
		hcFloat phi			= nPhi    		* 2 * PI / numPhi;
		Vec3D pos			= Vec3D(height, theta, phi);
		Vec3D posSpher		= grid->isElliptical() ? pos.convCoordSpher2Ell(*egrid) : pos;
		Vec3D posC			= compCoords ? pos : posSpher;

		//coords[i]			= Vec2D(posC[1], posC[2]);
		coords[i]			= posC;
		posComp[i]			= posC;
		workedUpon[i] 		= false;
	}

	pthread_t 				threads[		NUMTHREADS];
	volatile _Atomic_word 	threadRunning[	NUMTHREADS];
	threadParamMag 			tParams[		NUMTHREADS];

	pthread_mutex_t runningMutex 	= PTHREAD_MUTEX_INITIALIZER;
	volatile _Atomic_word numRunningThreads = 0;

	for(uint i=0; i<NUMTHREADS; ++i)
	{
		threadRunning[i] = 0;
		tParams[i].init(i, &numRunningThreads, &runningMutex, &threadRunning[i], solver.grid);
	}

	bool workLeft = true;

	while(workLeft)
	{
		workLeft = false;
		if(numRunningThreads < NUMTHREADS-1) // -1 should not be here, but then "if(j==numMaxThreads-1)" is triggered sometimes
		{
			for(uint i=0; i<numMaglines; ++i)
			{
				bool breaker = false;

				if(!workedUpon[i])
				{
					workedUpon[i] 		= true;
					uint nPhi			= i 		 	% numPhi;
					uint nTheta			= (i - nPhi) 	/ numPhi;
					Magline *magline	= &(operator()(nTheta, nPhi));
					Vec3D *posC			= &posComp[i];

					for(uint j=0; j<NUMTHREADS; ++j)
					{
						if(threadRunning[j] == 0)
						{
							pthread_mutex_lock(&runningMutex);
							++numRunningThreads;
							threadRunning[j] = 1;
							pthread_mutex_unlock(&runningMutex);
							tParams[j].set(magline, posC);

							int rc = pthread_create(&threads[j], NULL, MagMapping::createAtHeight_threadEntryPoint, (void*)(&tParams[j]));

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
							cerr << __FILE__ << "/" << __LINE__ << " No free thread found! (you should not be able to see this. If you do, I fucked up.... sorry)\n";
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

		for(uint i=0;i<numMaglines;++i)
			if(!workedUpon[i])
			{
				workLeft = true;
				break;
			}
	}

	while(numRunningThreads > 0)
		sleep(1);

	delete [] workedUpon;
	delete [] posComp;

	return true;
}

/*! entry function for multicore creation of magnetic field lines
 */
void *MagMapping::createAtHeight_threadEntryPoint(void *parameter)
{
	threadParamMag *param 	= (threadParamMag*)parameter;
	SphericalGrid *grid		= param->grid;
	bool debug 				= false;

	if(debug)
	{
		cout << "------------------------------------------------------\n";
		cout << "Thread:     " << param->threadID << "\n";
		cout << "NumRunning: " << *param->numRunningThreads << "\n";
		//cout << "Pos:        " << *(param->posStart) << "/" << *(param->posStart)[1] << "/" << *(param->posStart)[2] << "\n";

		fflush(stdout);
	}
	param->magline->createMaglineThroughPos(*grid, *param->posStart, debug);

	pthread_mutex_lock(param->runningMutex);
	--(*param->numRunningThreads);
	*param->threadRunning = 0;
	pthread_mutex_unlock(param->runningMutex);
	pthread_exit(NULL);
}

void MagMapping::createMappingHeader(char *header, hcFloat height, hcFloat *otherHeights, uint numHeights, hcFloat lowerR, hcFloat upperR)
{
	char tempstr[1000];
	char tempstr2[1000];
	hcDate now;
	now.setFromSystemTime();
	string timestamp = now.toString();

	header[0] = '\0';
	sprintf(header, "# pfss version %u.%u\n", PFSS_VER_MAJ, PFSS_VER_MIN);
	strcat(header, "# author:\t\tMartin A. Kruse (kruse@physik.uni-kiel.de)\n#\n");


	tempstr2[0]      = '\0';
	sprintf(tempstr2, "# time stamp: %s\n#\n", timestamp.data());
	strcat(header, tempstr2);

	strcat(header, "# magnetic configuration map\n");
	strcat(header, "# format:\n#\n");
	strcat(header, "#\ttheta(height) phi(height) X theta(phot) phi(phot) theta(ss) phi(ss) theta(height1) phi(height1) ....\n#\n");
	strcat(header, "# where X is the polarity flag (P positive, N negative, C closed, I invalid)\n");
	strcat(header, "# in case of closed field lines, the coordinates theta(ss) and phi(ss) give the coordinates of the second\n");
	strcat(header, "# footpoint on the photosphere. A '-' at some coordinate position means that the magnetic field line\n");
	strcat(header, "# does not extend to that specific height\n#\n");

	sprintf(tempstr, "# sinLatFormat: %u (0 - theta spacing linear in latitude, 1 - linear in sine(latitude))\n", sinLatFormat);
	strcat(header, tempstr);
	sprintf(tempstr, "# computationalCoordinates: %u (0 - height in world coordinates, 1 - height in computational coordinates)\n", sinLatFormat);
	strcat(header, tempstr);

	tempstr[0]      = '\0';
	sprintf(tempstr, "# height of map (m):\t\t%E\n# number of intermed. height levels: %u\n# height levels (m):\n#\n", height, numHeights);
	strcat(header, tempstr);

	tempstr[0]      = '\0';
	sprintf(tempstr, "# 0\t%E (base height of this map - m)\n", height);
	strcat(header, tempstr);

	tempstr[0]      = '\0';
	sprintf(tempstr, "# 1\t%E (height of lower boundary - m)\n", lowerR);
	strcat(header, tempstr);

	tempstr[0]      = '\0';
	sprintf(tempstr, "# 2\t%E (height of source sourface - m)\n", upperR);
	strcat(header, tempstr);

	for(uint i=0; i<numHeights; ++i)
	{
		tempstr[0]      = '\0';
		sprintf(tempstr, "# %u\t%E\n", i+3, otherHeights[i]);
		strcat(header, tempstr);
	}

	tempstr[0] = '\0';
	sprintf(tempstr, "#\n# each line thus has the following format:\n#\n");
	strcat(tempstr, "# theta(height) phi(height) X theta(phot) phi(phot) theta(ss) phi(ss) ");
	strcat(header, tempstr);
	for(uint i=0; i<numHeights; ++i)
	{
		tempstr[0] = '\0';
		sprintf(tempstr, "theta(%E) phi(%E) ", otherHeights[i], otherHeights[i]);
		strcat(header, tempstr);
	}
	strcat(header, "\n#\n");
}

bool MagMapping::exportImagePolarity()
{
	string fn = getFilename_magMappingImg(info, height, sinLatFormat, compCoords, numTheta, numPhi);
	if(!checkComputed())		return false;
	if(!checkExportedFile(fn))	return false;

	hcImageRGBA img(numPhi, numTheta);

	for(uint j=0; j<numTheta; ++j)
		for(uint k=0; k<numPhi; ++k)
		{
			uint xi			= k;
			uint yi			= numTheta - j - 1;
			Magline &magL   = *maglines[index(j, k)];                       	// magline considered
			//Vec2D &coord    = coords[index(j, k)];							// pivot coordinates of magline

			if(!magL.valid)				img(xi, yi) = Magline::colorInvalid;	// determine pixel colors
			else
			{
				if(magL.closed)			img(xi, yi)	= Magline::colorClosed;
				else
				{
					if(magL.polarity)	img(xi, yi) = Magline::colorPositive;
					else				img(xi, yi)	= Magline::colorNegative;
				}
			}
		}

	return img.save(fn);
}

bool MagMapping::exportImageMagfield()
{
	string fn_r = getFilename_magMappingMagfield_r(info, height, sinLatFormat, compCoords, numTheta, numPhi);
	string fn_t = getFilename_magMappingMagfield_t(info, height, sinLatFormat, compCoords, numTheta, numPhi);
	string fn_p = getFilename_magMappingMagfield_p(info, height, sinLatFormat, compCoords, numTheta, numPhi);
	if(!checkComputed())																	return false;
	if(!checkExportedFile(fn_r) && !checkExportedFile(fn_t) && !checkExportedFile(fn_p))	return false;

	hcFloat dSinLat		= 2*maxSinLat 		/ (numTheta - 1);
	//hcFloat dLat		= 2*asin(maxSinLat)	/ (numTheta - 1);
	hcFloat dLong		= 2 * PI / numPhi * 360/ 2.0 / PI;
	hcFloat crpix1		= (numPhi) 		/ 2.0;
	hcFloat crpix2		= (numTheta+1) 	/ 2.0;
	hcFloat crval1		= 180;
	hcFloat crval2		= 0;
	hcFloat nullval		= 0.0;

	hcImageFITS img_r(numPhi, numTheta);
	hcImageFITS img_t(numPhi, numTheta);
	hcImageFITS img_p(numPhi, numTheta);

	PFSSsolutionInfo info;
	hcFloat height;
	bool sinLatFormat, compCoords;
	uint resTheta, resPhi;

	getParamFromFN_magMapping(fn_r, info, height, sinLatFormat, compCoords, resTheta, resPhi);

	EllipticalGrid egr;
	if(info.method == METH_ELLIPTICAL)
		egr.init(info.sinLatFormat, info.maxSinLat, -info.maxSinLat, r_sol, info.rss, info.numR, false, info.ell);

	for(uint j=0; j<numTheta; ++j)
		for(uint k=0; k<numPhi; ++k)
		{
			Magline &magl 	= operator()(j, k);
			if(!magl.valid)
			{
				img_r(k, numTheta - j - 1) = 0.0;
				img_t(k, numTheta - j - 1) = 0.0;
				img_p(k, numTheta - j - 1) = 0.0;
				continue;
			}

			hcFloat val_r = 0.0;
			hcFloat val_t = 0.0;
			hcFloat val_p = 0.0;

			for(uint m=0; m<magl.numSamples; ++m)
				if((magl.posdata[m] - coords[index(j,k)]).length() < 1E-3)
				{
					val_r	= magl.magdata[m][0];
					val_t	= magl.magdata[m][1];
					val_p	= magl.magdata[m][2];
					break;
				}

			img_r(k, numTheta - j - 1) = val_r;
			img_t(k, numTheta - j - 1) = val_t;
			img_p(k, numTheta - j - 1) = val_p;
		}
	bool retval 	= true;

	retval &= img_r.save(fn_r);
	retval &= img_t.save(fn_t);
	retval &= img_p.save(fn_p);

	hcDate now;

	for(uint i=0; i<3; ++i)
	{
		hcImageFITS &img = (i==0 ? img_r : (i==1 ? img_t : img_p));
		stringstream comment;
		comment << "This file was generated by the PFSS program (version: " << PFSS_VER_MAJ << ":" << PFSS_VER_MIN << ") developed at IEAP, Kiel university, Author: Martin A. Kruse (kruse@physik.uni-kiel.de), File generated: " << now.toString() << " ";
		img.writeKeyComment(comment.str()); comment.str("");

		if(info.method == METH_ELLIPTICAL)
		{
			if(compCoords)
			{
				comment << "WARNING: This map uses computational coordinates in an elliptical grid. ";
				comment << "World coordinates presented in this file have to be transformed using the coordinate transformations from elliptic to cartesian/spheric coordinates.";
			}
			comment << "The elliptical PFSS solver and its data products have not been thoroughly reviewed.";
		}
		comment << "Parameters of the PFSS solution: Model: " << getStringFromModelID(info.model) << ", group: " << getStringFromGroupID(info.group) << ", method: " << getStringFromMethodID(info.method);
		comment << "rss: " << info.rss / r_sol << " r_sol, ellipticity: " << info.ell << ", SHC order: " << info.orderSHC << ", numR: " << info.numR << ", numTeta: " << info.numTheta << ", numPhi: " << info.numPhi;
		comment << ", dailyID: " << info.dailyID << ", computation time: " << info.computationTime << ", date computed: " << info.dateComputed.toString();
		img.writeKeyComment(comment.str()); comment.str("");

	//	expansion.writeKeyFloat("DSINLAT", "Spacing in Sine-Latitude direction", dSinLat);
		img.writeKeyFloat("MAXSL", "max(sin(latitude)), northern sin(latitude) of image boundary", maxSinLat);
		img.writeKeyString("CTYPE1", "longitude", "CRLN-CEA");
		img.writeKeyString("CTYPE2", "sine latitude", "CRLT-CEA");
		//expansion.writeKeyString("CTYPE1", "", "HPLN-TAN");
		//expansion.writeKeyString("CTYPE2", "", "HPLT-TAN");
		img.writeKeyFloat("CDELT1", "coordinate increment along axis 1", dLong);
		img.writeKeyFloat("CDELT2", "coordinate increment along axis 2", dSinLat);
		img.writeKeyFloat("CRPIX1", "coord. system reference pixel 1", crpix1);
		img.writeKeyFloat("CRPIX2", "coord. system reference pixel 2", crpix2);
		img.writeKeyFloat("CRVAL1", "coordinate at ref. pixel 1", crval1);
		img.writeKeyFloat("CRVAL2", "coordinate at ref. pixel 2", crval2);
		img.writeKeyFloat("CROTA1", "coord. system rotation angle 1", nullval);
		img.writeKeyFloat("CROTA2", "coord. system rotation angle 2", nullval);
		img.writeKeyString("BUNIT", "unit of pixel values", "GAUSS");
		img.writeKeyString("CUNIT1", "coordinate unit 1", "degree");
		img.writeKeyString("CUNIT2", "coordinate unit 2", "sinlat");
		img.writeKeyString("WCSNAME", "World Coordinate system name", "Heliographic");
	}
	return retval;
}

bool MagMapping::exportImageFootpoint(string fn)
{
	if(!checkComputed())	return false;

	if(!createFolderStructureTo(fn.c_str()))
	{
		printErrMess(__FILE__, __LINE__, "requested file '" + fn + "' cannot be created");
		return false;
	}

	hcImageRGBA imgFoot(numPhi, numTheta);									// debug footpoint image at photosphere

	for(uint j=0; j<numTheta; ++j)											// initialize footpoint img black
		for(uint k=0; k<numPhi; ++k)
			imgFoot(k,j) = char2RGBA8(0,0,0,255);

	hcFloat maxSinLat	= sin(PI/2.0 - coords[index(0, 0)][0]);
	hcFloat minSinLat	= sin(PI/2.0 - coords[index(numTheta-1, 0)][0]);	// TODO symmetric
	hcFloat tfX			= (numPhi-1)	/ (2*PI);
	hcFloat tfY			= (numTheta-1) 	/ (maxSinLat - minSinLat);

	for(uint j=0; j<numTheta; ++j)
		for(uint k=0; k<numPhi; ++k)
		{
			Magline &magL   = *maglines[index(j, k)];                       // magline considered
			Vec2D foot		= Vec2D(magL.posdata[0][1], magL.posdata[0][2]);// footpoint on photosphere (one of them if magline closed)
			int x			= (uint)round(foot[1]*tfX);						// img position of footpoint
			uint y			= (uint)round((sin(PI/2.0-foot[0])-minSinLat)*tfY);

			if(magL.valid && y<numTheta)
			{
				if(magL.closed)			imgFoot(x, y) 	= Magline::colorClosed;
				else
				{
					if(magL.polarity)	imgFoot(x,y) 	= Magline::colorPositive;
					else				imgFoot(x,y)	= Magline::colorNegative;
				}
			}
		}

	return imgFoot.save(fn.c_str());
}

/*! outputs textfiles where every line corresponds to one pixel position of the map
 *
 *  the output format is described in the header of the output file
 */
bool MagMapping::exportASCII(string filename, hcFloat *heights, uint numHeights, hcFloat lowerR, hcFloat upperR)
{
	if(!checkComputed())	return false;

	if(!createFolderStructureTo(filename.c_str()))
	{
		printErrMess(__FILE__, __LINE__, "requested file '" + filename + "' cannot be created");
		return false;
	}

    FILE *output = fopen(filename.c_str(), "w");

    char outputstring[100000];
    char tempstr[10000];

    createMappingHeader(outputstring, height, heights, numHeights, lowerR, upperR);
    fprintf(output, "%s", outputstring);

    for(uint j=0; j<numTheta; ++j)
        for(uint k=0; k<numPhi; ++k)
        {
            Magline &magL   = *maglines[index(j, k)];                       // magline considered
            //Vec2D &coord    = coords[index(j, k)];							// pivot coordinates of magline
            Vec3D &coord    = coords[index(j, k)];							// pivot coordinates of magline
            uint indL		= 0;
            uint indU		= magL.numSamples - 1;

            outputstring[0] = '\0';
            sprintf(outputstring, "%E %E %E ", coord[0], coord[1], coord[2]);

            if(!magL.valid)                	strcat(outputstring, "I - - / - -");
            else
            {
                tempstr[0] = '\0';

                if(magL.closed)				sprintf(tempstr, "C ");
                else
                {
                    if(magL.polarity)		sprintf(tempstr, "P ");
                    else					sprintf(tempstr, "N ");
                }
                strcat(outputstring, tempstr);
                tempstr[0] = '\0';
                sprintf(tempstr, "%E %E / %E %E ", magL.posdata[indL][1], magL.posdata[indL][2], magL.posdata[indU][1], magL.posdata[indU][2]);
                strcat(outputstring, tempstr);

                for(uint i=0; i<numHeights; ++i)
                {
                    tempstr[0] = '\0';
                    Vec3D pos[MAGLINE_NUM_POSITIONS], mag[MAGLINE_NUM_POSITIONS];
                    int numPos = magL.getAllValuesAtHeight(heights[i], pos, mag);
                    if(numPos == 0)
                        sprintf(tempstr, "/ 0 ; - - ");
                    else
                    {
                    	sprintf(tempstr, "/ %u ", numPos);
                    	char tstr[30];
                    	for(int l=0; l<numPos; ++l)
                    	{
                    		sprintf(tstr, "; %E %E ", pos[l][1], pos[l][2]);
                    		strcat(tempstr, tstr);
                    	}
                    }

                    strcat(outputstring, tempstr);
                }
            }
            fprintf(output, "%s\n", outputstring);
        }

    fclose(output);
    return true;
}


bool MagMapping::exportExpansionFactorImage()
{
	string fn = getFilename_magMappingExpansion(info, height, sinLatFormat, compCoords, numTheta, numPhi);
	if(!checkComputed())		return false;
	if(!checkExportedFile(fn))	return false;

	hcFloat dSinLat		= 2*maxSinLat 		/ (numTheta - 1);
	//hcFloat dLat		= 2*asin(maxSinLat)	/ (numTheta - 1);
	hcFloat dLong		= 2 * PI / numPhi * 360/ 2.0 / PI;
	hcFloat crpix1		= (numPhi) / 2.0;
	hcFloat crpix2		= (numTheta+1) / 2.0;
	hcFloat crval1		= 180;
	hcFloat crval2		= 0;
	hcFloat nullval		= 0.0;

	hcImageFITS expansion(numPhi, numTheta);

	for(uint j=0; j<numTheta; ++j)
		for(uint k=0; k<numPhi; ++k)
		{
			Magline &magl 	= operator()(j, k);
			if(magl.closed || !magl.valid)
			{
				expansion(k, numTheta - j - 1) = 0.0;
				continue;
			}

			hcFloat rl		= magl.posdata[0][0];
			hcFloat ru		= magl.posdata[magl.numSamples-1][0];
			hcFloat bl		= magl.magdata[0].length();
			hcFloat bu		= magl.magdata[magl.numSamples-1].length();
			hcFloat factor 	= bl * rl * rl / (bu * ru * ru);
			expansion(k, numTheta - j - 1)	= factor;
		}
	bool retval 	= true;

	retval &= expansion.save(fn);

	/*//	This creates a bitmap image of the expansion factor, please use this only for debugging
	hcFloat valMax	= 10000;
	hcFloat valMin	= 1;
	hcFloat conv	= 255 / log10(valMax/valMin);
	hcImageRGBA bitmap(numPhi, numTheta);

	for(uint j=0; j<numTheta; ++j)
		for(uint k=0; k<numPhi; ++k)
		{
			uint x 			= k;
			uint y 			= numTheta - j - 1;
			hcFloat val_in 	= expansion(x,y) < valMin ? valMin : (expansion(x,y) > valMax ? valMax : expansion(x,y));
			hcFloat val_out	= conv * (log10(val_in)-log10(valMin));
			bitmap(x,y)		= char2RGBA8(val_out, val_out, val_out, 255);
		}
	retval &= bitmap.save(fn_bitmap.data());//*/

	hcDate now; now.setFromSystemTime();
	stringstream comment;
	comment << "This file was generated by the pfss program (version: " << PFSS_VER_MAJ << ":" << PFSS_VER_MIN << ") developed at IEAP, Kiel university, Author: Martin A. Kruse (kruse@physik.uni-kiel.de), File generated: " << now.toString() << " ";
	expansion.writeKeyComment(comment.str()); comment.str("");

	PFSSsolutionInfo info;
	hcFloat height;
	bool sinLatFormat, compCoords;
	uint resTheta, resPhi;

	getParamFromFN_magMapping(fn, info, height, sinLatFormat, compCoords, resTheta, resPhi);

	if(info.method == METH_ELLIPTICAL)
	{
		if(compCoords)
		{
			comment << "WARNING: This map uses computational coordinates in an elliptical grid. ";
			comment << "World coordinates presented in this file have to be transformed using the coordinate transformations from elliptic to cartesian/spheric coordinates.";
		}
		comment << "The elliptical PFSS solver and its data products have not been thoroughly revied.";
	}
	comment << "Parameters of the PFSS solution: Model: " << getStringFromModelID(info.model) << ", group: " << getStringFromGroupID(info.group) << ", method: " << getStringFromMethodID(info.method);
	comment << "rss: " << info.rss / r_sol << " r_sol, ellipticity: " << info.ell << ", SHC order: " << info.orderSHC << ", numR: " << info.numR << ", numTeta: " << info.numTheta << ", numPhi: " << info.numPhi;
	comment << ", dailyID: " << info.dailyID << ", computation time: " << info.computationTime << ", date computed: " << info.dateComputed.toString();
	expansion.writeKeyComment(comment.str()); comment.str("");
	expansion.writeKeyFloat("MAXSL", "max(sin(latitude)), northern sin(latitude) of image boundary", maxSinLat);
	expansion.writeKeyString("CTYPE1", "longitude", "CRLN-CEA");
	expansion.writeKeyString("CTYPE2", "sine latitude", "CRLT-CEA");
	//expansion.writeKeyString("CTYPE1", "", "HPLN-TAN");
	//expansion.writeKeyString("CTYPE2", "", "HPLT-TAN");
	expansion.writeKeyFloat("CDELT1", "coordinate increment along axis 1", dLong);
	expansion.writeKeyFloat("CDELT2", "coordinate increment along axis 2", dSinLat);
	expansion.writeKeyFloat("CRPIX1", "coord. system reference pixel 1", crpix1);
	expansion.writeKeyFloat("CRPIX2", "coord. system reference pixel 2", crpix2);
	expansion.writeKeyFloat("CRVAL1", "coordinate at ref. pixel 1", crval1);
	expansion.writeKeyFloat("CRVAL2", "coordinate at ref. pixel 2", crval2);
	expansion.writeKeyFloat("CROTA1", "coord. system rotation angle 1", nullval);
	expansion.writeKeyFloat("CROTA2", "coord. system rotation angle 2", nullval);
	expansion.writeKeyString("CUNIT1", "coordinate unit 1", "degree");
	expansion.writeKeyString("CUNIT2", "coordinate unit 2", "sinlat");
	expansion.writeKeyString("WCSNAME", "World Coordinate system name", "Heliographic");
	return retval;
}

void MagMapping::dump(uint indent)
{
	uint numValid = 0;
	for(uint i=0; i<numTheta; ++i)
		for(uint j=0; j<numPhi; ++j)
			if(operator()(i,j).valid) ++numValid;

	stringstream ind;
	if(indent > 0) ind << setw(indent) << setfill(' ') << " ";
	cout << ind.str() << "Dumping MagMapping:\n";
	info.dump(indent+1);
	cout << ind.str() << setw(20) << setfill(' ' ) << std::left << "numTheta:" 	<< numTheta << "\n";
	cout << ind.str() << setw(20) << setfill(' ' ) << "numPhi:" 				<< numPhi 	<< "\n";
	cout << ind.str() << setw(20) << setfill(' ' ) << "numValid:"				<< numValid << "\n";
}

/*
void plotLine(hcImageRGBA &img, hcImageFloat &zbuf, hcImageBool &occ, const uint &color, const Vec3D &posStart, const Vec3D &posEnd)
{
	uint width	= img.width;
	uint height	= img.height;
	hcFloat tfx	= (width-1)/(2.0);
	hcFloat tfy	= (height-1)/(2.0);
	Vec3D slope	= posEnd-posStart;

	int x0		= round((posStart[0] + 1.0)*tfx);
	int x1		= round((posEnd[0]   + 1.0)*tfx);
	int y0		= round((posStart[1] + 1.0)*tfy);
	int y1		= round((posEnd[1]   + 1.0)*tfy);

	int dx		= x1-x0 >= 0 ? 1 : -1;
	int dy		= y1-y0 >= 0 ? 1 : -1;

	bool b		= fabs(x1-x0)>=fabs(y1-y0);

	for(int m=b?x0:y0;b?(dx==1?m<=x1:m>=x1):(dy==1?m<=y1:m>=y1); m+=b?dx:dy)
	{
		hcFloat posX	= m/tfx - 1.0;
		hcFloat posY	= m/tfy - 1.0;
		hcFloat tx		= (posX - posStart[0])/slope[0];
		hcFloat ty		= (posY - posStart[1])/slope[1];
		Vec3D pos		= posStart + (b?tx:ty)*slope;
		if(b) posY		= pos[1];
		else  posX		= pos[0];
		hcFloat posZ	= pos[2];
		int x			= b ? m							: round((posX + 1.0)*tfx);
		int y			= b ? round((posY + 1.0)*tfy) 	: m;

		if(x>=0 && x<width && y>=0 && y<height)
			if(posZ < zbuf(x,y) && (occ(x,y) != false || posZ < 0.0))
			{
				zbuf(x,y) 	= posZ;
				img(x,y)	= color;
			}
	}
}

void createOcculterMap(hcImageBool &occ, Imager &imager, hcFloat clipRadius, hcFloat sizeHor, hcFloat sizeVert,
		const hcDate &date, hcFloat occRadius, bool invert)
{
	Vec3D lookDir, up, right;
	imager.getBaseVectors(lookDir, up, right, date, clipRadius);
	Vec3D pos			= imager.getPos(date);
	Vec3D target		= imager.getTarget(date);
	hcCamera3D cam		= imager.getView(date, clipRadius);
	hcFloat distCam		= (target - pos).length();
	hcFloat near		=  distCam - clipRadius;
	hcFloat far			=  distCam + clipRadius;

	hcPerspectiveProjection3D projSun(near, far, -r_sol, r_sol, -r_sol, r_sol);

	hcFloat tfX			= sizeHor/(occ.width-1);
	hcFloat tfY			= sizeVert/(occ.height-1);

	for(uint x=0; x<occ.width; ++x)
		for(uint y=0; y<occ.height; ++y)
		{
			hcFloat posX	= x*tfX-sizeHor/2.0;
			hcFloat posY	= y*tfY-sizeVert/2.0;
			hcFloat posZ	= near;
			Vec3D posPix	= pos + posZ*lookDir + posX*right + posY*up;

			Vec4D posPixHom	= Vec4D(posPix[0], posPix[1], posPix[2], 1.0f);
			Vec4D projC2	= projSun * cam * posPixHom;
			projC2			/= projC2[3];

			hcFloat projX	= projC2[0];
			hcFloat projY	= projC2[1];

			//cout << x << "/" << y << " " << projX << "/" << projY << "\n";

			occ(x,y)	= (sqrt(projX*projX + projY*projY) < occRadius ? (invert ? true : false) : (invert ? false : true));
		}
}///

bool MagMapping::createProjectedView(const string fn, Imager &corImg, Imager &photImg, const hcObliquity &solarRotAxis, const hcDate &date)
{
	cout << "MagMapping::createProjectedView started\n";
	fflush(stdout);

	if(maglines == NULL)
	{
		cerr << "MagMapping::createProjectedView: magneticLines not initialized!\n";
		return false;
	}

	uint width			= 2000;
	uint height			= 1000;
	hcFloat aspRatio	= (hcFloat)(width)/height;

	hcFloat clipRadius	= 5*r_sol;	// horizontal clipradius
	Vec3D target		= corImg.getTarget(date);
	Vec3D pos			= corImg.getPos(date);
	hcFloat distCam		= (target - pos).length();
	hcFloat fov_2		= asin(clipRadius/distCam);
	hcFloat d			= (distCam - clipRadius) * tan(fov_2);
	hcFloat dh			= d;
	hcFloat dv			= d/aspRatio;
	hcFloat near		=  distCam - clipRadius;
	hcFloat far			=  distCam + clipRadius;

	Matrix4x4 solarOrientation = solarRotAxis.getOrientation(date);
	hcPerspectiveProjection3D project(near, far, -dh, dh, -dv, dv);
	hcCamera3D cam(project, pos);
	cam.setPos(pos);
	cam.lookat(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 0.0, 1.0));
	hcImageRGBA img(width, height);
	hcImageBool occ(width, height);
	hcImageFloat zbuf(width, height);

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			img(x,y) 	= char2RGBA8(0,0,0,255);
			zbuf(x,y)	= 1.0;
		}

	// ---------------------------------------------------------------------------------------
	// 																			Overlay Coronagraph image
	// ---------------------------------------------------------------------------------------

	hcFloat sizeHor		= 2*dh;
	hcFloat sizeVert	= 2*dv;
	string obsFN;

	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 2*clipRadius/r_sol, true);
	corImg.drawImage(img, occ, sizeHor, sizeVert, date, clipRadius, obsFN);

	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 1.2, true);
	photImg.drawImage(img, occ, sizeHor, sizeVert, date, clipRadius, obsFN);

	// ---------------------------------------------------------------------------------------
	// 																			print projected magmap
	// ---------------------------------------------------------------------------------------

	MagMapping projMap = *this;

	for(uint x=0; x<numPhi; ++x)
		for(uint y=0; y<numTheta; ++y)
		{
			Magline &magl 		= operator()(y,x);
			Magline &projMagl	= projMap(y,x);

			//printf("%u/%u\t%u %u\n\t\t%u %u\n", x,y,magl.closed, magl.polarity, projMagl.closed, projMagl.polarity);

			for(uint i=0; i<magl.numSamples; ++i)
			{
				Vec3D pos 		= magl.posdata[i].convCoordSpher2Cart();
				Vec4D homPos	= Vec4D(pos[0], pos[1], pos[2], (float)1.0);
				Vec4D projPos	= cam.projection * cam * solarOrientation * homPos;
				projPos			/= projPos[3];
				Vec3D projPos3D = Vec3D(projPos[0], projPos[1], projPos[2]);

				projMagl.posdata[i] = projPos3D;
			}
		}

	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 1.0, false);



	for(uint x=0; x<numPhi; ++x)
		for(uint y=0; y<numTheta; ++y)
		{
			Magline &magl	= projMap(y,x);
			uint color		= 	 !magl.valid 	? Magline::colorInvalid 	:
								(magl.closed 	? Magline::colorClosed 		:
								(!magl.polarity ? Magline::colorNegative	:
												  Magline::colorPositive));

			for(uint i=1; i<magl.numSamples; ++i)
			{
				Vec3D &posStart	= magl.posdata[i-1];
				Vec3D &posEnd	= magl.posdata[i];

				plotLine(img, zbuf, occ, color, posStart, posEnd);
			}
		}

	// ---------------------------------------------------------------------------------------
	// 																			Draw the sun
	// ---------------------------------------------------------------------------------------

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			if(y==0){ cout << x << "\n";fflush(stdout);}
			hcFloat posX	= x*tfX-d;
			hcFloat posY	= y*tfY-d;
			hcFloat posZ	= near;
			Vec3D posPix	= pos + posZ*lookDir + posX*right + posY*up;

			hcFloat res0, res1;
			hcLine3D line(pos, posPix);
			hcDSphere3D sun(Vec3D(0.0, 0.0, 0.0), r_sol);
			if(line.intersectsSphere3D(sun, res0, res1, false)==2)
			{
				Vec3D pos0, pos1;
				line.getPos(res0, pos0);
				line.getPos(res1, pos1);


				//pos0.convCoordCart2Spher().dump();
				//pos1.convCoordCart2Spher().dump();

				Vec4D hom0(pos0[0], pos0[1], pos0[2], 1.0f);
				Vec4D hom1(pos1[0], pos1[1], pos1[2], 1.0f);

				Vec4D proj0	= cam.projection * cam * hom0;
				Vec4D proj1	= cam.projection * cam * hom1;
				proj0		/= proj0[3];
				proj1		/= proj1[3];

				if(x>=0 && x<width && y>=0 && y<height)
					if(min(proj0[2], proj1[2]) < zbuf(x,y))
					{
						zbuf(x,y) 	= min(proj0[2], proj1[2]);
						img(x,y)	= char2RGBA8(150,150,150,255);
					}
			}
		}///


	img.save(fn.data());

	return true;
}//*/
