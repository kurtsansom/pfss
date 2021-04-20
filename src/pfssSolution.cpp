#include "src/pfssSolution.h"
#include "src/carRotInfo.h"
#include "src/pfss.h"
#include "src/ellipticalGrid.h"
#include "src/filenames.h"

#include <fstream>
#include <string>
#include <random>
#include <iomanip>

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/format.hpp"

#ifdef GUI
#include "src/euv.h"
#endif

namespace fs = boost::filesystem;
using namespace boost;
using namespace std;

extern string dirData;

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          PFSSsolution
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

PFSSsolution::PFSSsolution()
{
	initNULL();
	init();
}

PFSSsolution::PFSSsolution(const PFSSsolution &other)
{
	initNULL();
    *this = other;
}

PFSSsolution &PFSSsolution::operator=(const PFSSsolution &other)
{
    if(this == &other)
        return *this;

    info         = other.info;
    solver       = other.solver;
    photBoundary = other.photBoundary;

    return *this;
}

void PFSSsolution::initNULL()
{}

void PFSSsolution::clear()
{
	initNULL();
}

void PFSSsolution::init()
{
	clear();
}

bool PFSSsolution::isDirSet()
{
	if(dirData.length() == 0)
	{
		cerr << __FILE__ << ":" << __LINE__ << ": output directory has not been set.\n";
		return false;
	}
	return true;
}

string PFSSsolution::getFilenamePhotMagfield()
{
	return dirData + "/" + to_string(info.CRnum) + "/" + getFilename_photMagfield(info);
}

bool PFSSsolution::loadPhotBoundary(const string &filename, uint seed, uint scaleMethod)
{
    if(!checkFileEx(filename, "PFSSsolution::loadPhotBoundary"))
    	return false;
    if(!photBoundary.load(filename))
    	return false;

    if(seed!=0)	photBoundary.addSignedFractionNoise(0.1, seed); 	// testing noise impact on computation

    if(!photBoundary.remeshImage(info.numPhi, info.numTheta, scaleMethod))
		return false;
    if(!solver.init(photBoundary.synInfo.sinLatFormat, photBoundary.synInfo.maxSinLat, info.rss, info.numR, info.ell))
    	return false;

    ((SynopticInfo*)(&info))->operator =(photBoundary.synInfo);

	string photMagfieldFN = getFilenamePhotMagfield();
	createFolderStructureTo(photMagfieldFN.data());
	photBoundary.save(photMagfieldFN.data());

	if(!photBoundary.load(photMagfieldFN.data()))
		return false;

    return true;
}

/*! stores computed PFSS solution to disk and stores meta-information into cfg-files
 *  (cfg file for information on this specific solution and one entry in super-cfg
 *  so that the program knows that there is a solution for this specific Carrington rotation
 *  and instrument)
 */
bool PFSSsolution::save()
{
    if(!doesFileExist(photBoundary.filename))
    {
    	cerr << __FILE__ << ":" << __LINE__ << ": File \n\t'" << photBoundary.filename << "'\n\tfor photospheric boundary does not exist.\n";fflush(stderr);
        return false;
    }

    if(!isDirSet()) return false;

    if(info.CRnum < 1600)
    {
    	cerr << __FILE__ << ":" << __LINE__ << ": Carrington rotation number (" << info.CRnum <<") not valid.\n";fflush(stderr);
        return false;
    }

    string gridFilename = getFilenameGrid();
    string cfgFilename	= getFilenameConfig();

    solver.grid->exportEntireGrid(gridFilename.data());

    hcFloat a = 1.0;
    if(solver.grid->isElliptical())
    {
    	EllipticalGrid &eGrid = *static_cast<EllipticalGrid*>(solver.grid);
    	uint numR	= eGrid.numR;
    	a			= eGrid.getEllA(numR-1);
    }

    string timeStamp 	= info.dateComputed.toString();

    ofstream config(cfgFilename);

    config << "rotationNum " 		<< info.CRnum 								<< "\n";
	config << "instrument " 		<< getStringFromOriginID(info.instrument) 	<< "\n";
	config << "model " 				<< getStringFromModelID(info.model) 		<< "\n";
	config << "method " 			<< getStringFromMethodID(info.method) 		<< "\n";
	config << "computed_by " 		<< getStringFromGroupID(info.group) 		<< "\n";
	config << "produced_at " 		<< timeStamp.data() 						<< "\n";
	config << "comptime " 			<< info.computationTime 					<< "\n";
	config << "solutionSteps "		<< info.solutionSteps	 					<< "\n";
	config << "orderSHC " 			<< info.orderSHC 							<< "\n";
	config << "numR " 				<< solver.grid->numR 						<< "\n";
	config << "numTheta " 			<< solver.grid->numTheta 					<< "\n";
	config << "numPhi " 			<< solver.grid->numPhi 						<< "\n";
	config << "sinLatGrid " 		<< solver.grid->sinLatGrid 					<< "\n";
	config << "maxSinLat " 			<< solver.grid->maxSinLat 					<< "\n";
	config << "upperR " 			<< solver.grid->upperR 						<< "\n";
	config << "geometricFactor " 	<< solver.grid->geometricFactor 			<< "\n";
	config << "ellipticity " 		<< a 										<< "\n";
	config << "sizeofFloat " 		<< info.sizeofFloat 						<< "\n";

    config.close();
    insertSolutionInSuperConfig();

    return true;
}

/*! given a rotation/instrument-cfg filename and a path where the solution is stored,
 *  this function imports a previously computed solution from disk *
 */
bool PFSSsolution::load(const string &cfgFN)
{
    if(!checkFileEx(cfgFN, "PFSSsolution::loadPrecomputed"))
    {
        cerr << __FILE__ << ":" << __LINE__ << ": File '" << cfgFN <<"' does not exist.\n";
        return false;
    }

    std::ifstream config(cfgFN);
    uint sizeofFloat		= 4;
    uint computationTime	= 0;
    hcFloat maxSinLat     	= 0.0;
    string boundaryFilename	= "";
    string gridFilename		= "";
    string timestamp		= "";
    string line				= "";

    getParamFromFN_pfssSolutionInfo(string(cfgFN), info);

    while(getline(config, line))
    {
    	string option 		= "";
    	string argument 	= "";
    	stringstream stream(line);
    	stream >> option >> argument;
    	if(option == "produced_at")			timestamp			= argument;
        if(option == "comptime")        	computationTime 	= (uint)atoi(argument.data());
		if(option == "maxSinLat")           maxSinLat 			= atof(argument.data());
        if(option == "sizeofFloat")			sizeofFloat 		= (uint)atoi(argument.data());
    }
    config.close();
    info.sizeofFloat		= sizeofFloat;
    info.sinLatFormat		= true;
    info.maxSinLat			= maxSinLat;
    info.dateComputed		= hcDate(timestamp);
    info.computationTime 	= computationTime;

    solver.init(info.sinLatFormat, info.maxSinLat, info.rss, info.numR, info.ell);
    SphericalGrid &gr 		= *solver.grid;
    EllipticalGrid &egr 	= *(dynamic_cast<EllipticalGrid*>(&gr));
    bool isEll				= gr.isElliptical();
    boundaryFilename 		= getFilenamePhotMagfield();
    gridFilename			= getFilenameGrid();

    if(!photBoundary.load(boundaryFilename.data()))			return false;
    //TODO: this does not load maxsinlat and doesnt init members like scaling factors
    if(!solver.grid->importEntireGrid(gridFilename.data())) return false;

    solver.solutionComputed = true;
    return true;
}

bool PFSSsolution::computeKielSHC(const string &filename, uint order, hcFloat r_ss, uint compResR)
{
	if (!doesFileExist(filename))
	{
		cerr << __FILE__ << ":" << __LINE__ << " input file \n'" << filename << "'\ndoes not exist\n";
		return false;
	}

	if(!isDirSet()) return false;

	SynopticInfo synInf;
	uint numT, numP;
	hcFloat geomFac;
	SphericalGrid::getOptimalGridParameters(compResR, r_sol, r_ss, geomFac, numT, numP);
	info.init(synInf, sizeof(hcFloat), MODEL_PFSS, METH_SHC, GROUP_KIEL, r_ss, 1.0, order, compResR, numT, numP);

	if(!loadPhotBoundary(filename))
	{
		cerr << __FILE__ << ":" << __LINE__ << " File\n'" << filename << "'\ncannot not be loaded. Continue with next file.\n";
		return false;
	}

	string coeffFN = dirData + "/" + to_string(info.CRnum) + "/" + getFilename_harmonicCoeff(info);

	if(doesFileExist(coeffFN))
	{
		printStdOutMess(__FILE__, __LINE__, "Solution '" + coeffFN + "' already exists, skip recomputation.");
		return true;
	}

	PFSSsolution_SHC_sun coeffMagfield; // TODO: here you need to add r_ss to constructor
	coeffMagfield.determineCoefficientsFromPhotMagfield(order, r_ss, photBoundary);
	solver.solutionComputed = true;
	coeffMagfield.exportCoefficients(coeffFN.data());

	for(uint r=0;r<solver.grid->numR;++r)
		for(uint t=0;t<solver.grid->numTheta;++t)
			for(uint p=0;p<solver.grid->numPhi;++p)
			{
				uint ind 		= solver.grid->getIndex(r, t, p);
				Vec3D result;
				Vec3D pos = solver.grid->getPos(ind, false);
				coeffMagfield.eval(pos, result);
				solver.grid->setB(ind, result);
			}
	save();

	return true;
}

bool PFSSsolution::computeKielGrid(const string &filename, const string &optionalID, hcFloat r_ss, uint compResR, hcFloat ell, uint scaleMethod)
{
	if (!doesFileExist(filename))
	{
		printErrMess(__FILE__, __LINE__, "file '" + filename + "' does not exist.");
		return false;
	}

	if(!isDirSet()) return false;

	SynopticInfo synInf;
	uint numT, numP;
	hcFloat geomFac;
	SphericalGrid::getOptimalGridParameters(compResR, r_sol, r_ss, geomFac, numT, numP);
	info.init(	synInf, sizeof(hcFloat), MODEL_PFSS, ell==1.0?METH_NUMERIC:METH_ELLIPTICAL, GROUP_KIEL,
				r_ss, ell, 0, compResR, numT, numP);

	if(!loadPhotBoundary(filename, 0, scaleMethod))
		return false;

	string gridFilename = getFilenameGrid();
	string cfgFilename	= getFilenameConfig();

	if(doesFileExist(gridFilename.data()))
	{
		printStdOutMess(__FILE__, __LINE__, "Solution '" + gridFilename + "' already exists, skip recomputation.");
		return (load(cfgFilename.data()) ? true : false);
	}
	printStdOutMess(__FILE__, __LINE__,  "Solution '" + gridFilename + "'does not exist, start computation.");

	hcDate start, end;
	start.setFromSystemTime();

	int retval = solver.computeSolution(photBoundary);

	if(retval<0)	return false;

	end.setFromSystemTime();
	info.dateComputed.setFromSystemTime();
	info.computationTime 	= (uint)((end-start)/hcDate::facSec);
	info.solutionSteps		= retval;

	return save();
}

// TODO not working so far
void PFSSsolution::loadAndMapCSSS(	const char *boundaryFN, const char *coeffFN,
									uint compRadialRes, uint imgThetaRes, uint imgPhiRes, bool mapIntermediateHeights)
{
	if (!doesFileExist(boundaryFN) || !doesFileExist(coeffFN))
	{
		printf("ERROR!\tPFSSsoluion::evalCSSS\n\t input file '%s' or '%s' does not exist\n", boundaryFN, coeffFN);
		return;
	}

	if(!isDirSet()) return;

	SynopticInfo synInf;
	uint numT, numP;
	hcFloat geomFac;
	SphericalGrid::getOptimalGridParameters(compRadialRes, r_sol, 2.5*r_sol, geomFac, numT, numP);
	info.init(	synInf, sizeof(hcFloat), MODEL_CSSS, METH_NUMERIC, GROUP_BALA, 2.5*r_sol, 0.0, 0, compRadialRes, numT, numP);

	if(!loadPhotBoundary(boundaryFN))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": File '" << boundaryFN << "' cannot not be loaded. Continue with next file.\n";
		return;
	}
	solver.grid->evaluateCSSS(coeffFN);
	save();
	multiMapSolution(imgThetaRes, imgPhiRes, mapIntermediateHeights, true, true);

	printStdOutMess(__FILE__, __LINE__, "evalCSSS concluded");
}

string PFSSsolution::getFilenameConfig()
{
	return dirData + "/" + to_string(info.CRnum) + "/" + getFilename_pfssSolutionConfig(info);
}

string PFSSsolution::getFilenameGrid()
{
	return dirData + "/" + to_string(info.CRnum) + "/" + getFilename_pfssSolutionBin(info);
}

bool PFSSsolution::insertSolutionInSuperConfig()
{
	string cfgFN	= getFilenameConfig();
	string superFN	= getFilename_superconfig();
    char cfgFNin[1000];

    if(!createFolderStructureTo(superFN.data())) return false;
    std::ifstream config(superFN);
    stringstream temp;

    bool entryFound     = false;
    while(config >> cfgFNin)
    {
    	if(cfgFN != string(cfgFNin))
    		temp << cfgFNin	<< "\n";
    	else
    	{
    		temp << cfgFN	<< "\n";
    		entryFound = true;
    	}
    }

    if(!entryFound)	temp << cfgFN << "\n";

    config.close();
    std::ofstream supercfg(superFN);
    supercfg << temp.str();
    supercfg.close();

    return true;
}

LaplaceSolver &PFSSsolution::getSolver()
{
    return solver;
}

/*!	Magnetic field lines are invoked at a given height in a quasi equidistant manner (sine-latitude grid).
 *  The field lines are then traced through the heliosphere and sampled at other given heights.
 *  This way we can track the "evolution" of the magnetic field from different starting points.
 *
 *  @param path 		to output folder
 *  @param height 		at which quasi-equidistant field lines are to be invoked
 *  @param heights		other height levels to which the lines will be traced
 *  @param numHeights 	number of entries in heights
 *  @param maxSinLat	highest sine(latitude)
 *  @param minSinLat	lowest sine(latitude)
 *  @param numTheta		number of grid points in latitudinal direction
 *  @param numPhi		number of grid points in azimuthal direction
 *  @param filename		filename for multimapping to be exported to
 */
bool PFSSsolution::mapHeightLevel(
		hcFloat height, hcFloat *heights, uint numHeights,
		hcFloat maxSinLat, uint numTheta, uint numPhi,
		bool sinLatFormat, bool compCoords, bool exportASCII)
{
	if(!isDirSet()) return false;

	LaplaceSolver &solver = getSolver();

	uint numTheta_m			= numTheta==0 ? solver.grid->numTheta 	: numTheta;
	uint numPhi_m			= numTheta==0 ? solver.grid->numPhi		: numPhi;
	hcFloat maxSinLat_m		= numTheta==0 ? solver.grid->maxSinLat	: maxSinLat;
	bool sinLatGrid			= numTheta==0 ? solver.grid->sinLatGrid	: sinLatFormat;

	stringstream filenameAscii, filenameMap, filenameImg, filenameFootImg, filenameExpansion, filenameExpBitmap;
	string oDir	= dirData + "/" + to_string(info.CRnum) + "/";
	filenameAscii		<< oDir << getFilename_magMappingASCII(				info, height, sinLatGrid, compCoords, numTheta_m, numPhi_m);
	filenameMap			<< oDir << getFilename_magMappingBin(				info, height, sinLatGrid, compCoords, numTheta_m, numPhi_m);
	filenameImg			<< oDir << getFilename_magMappingImg(				info, height, sinLatGrid, compCoords, numTheta_m, numPhi_m);
	filenameFootImg		<< oDir << getFilename_magMappingFootImg(			info, height, sinLatGrid, compCoords, numTheta_m, numPhi_m);
	filenameExpansion	<< oDir << getFilename_magMappingExpansion(			info, height, sinLatGrid, compCoords, numTheta_m, numPhi_m);
	filenameExpBitmap	<< oDir << getFilename_magMappingExpansionBitmap(	info, height, sinLatGrid, compCoords, numTheta_m, numPhi_m);

	MagMapping map;
	if(doesFileExist(filenameMap.str()))
	{
		printStdOutMess(__FILE__, __LINE__, "map at height " + toStr(height) + " m does already exist, skip computation.");
		if(!map.importBinary(filenameMap.str().data()))
		{
			cerr << __LINE__ << "/" << __FILE__ << ": file " << filenameMap.str() << " cannot be loaded.\n";
			return false;
		}
	}
	else
	{
		printStdOutMess(__FILE__, __LINE__, "map at height " + toStr(height) + " m does not exist, start computation.");
		map.createAtHeight(info, solver, height, numTheta_m, numPhi_m, maxSinLat_m, sinLatGrid, compCoords);
		map.exportBinary(filenameMap.str().data());
	}

	if(exportASCII)	map.exportASCII(filenameAscii.str(), heights, numHeights, solver.grid->lowerR, solver.grid->upperR);

	map.exportImage(filenameImg.str());
	//map.exportFootpointImage(filenameFootImg);
	map.exportExpansionFactorImage(solver.grid, filenameExpansion.str(), filenameExpBitmap.str());

    return true;
}

string PFSSsolution::getMagneticMappingFilename(const MagMapping &magmap)
{
	if(!isDirSet()) return "";
	string retval = dirData + "/" + to_string(info.CRnum) + "/" + getFilename_magMappingBin(	info, magmap.getHeight(), magmap.sinLatFormat, magmap.compCoords, magmap.numTheta, magmap.numPhi);
	return retval;
}

/*! creates magnetic field maps with footpoints at several heights and
 *  mapped to all other heights in this list.
 */
bool PFSSsolution::multiMapSolution(uint numTheta, uint numPhi, bool computeIntermediateHeightLevels, bool sinLatFormat, bool compCoords)
{
	uint numLevels		= 2;
	hcFloat *heights 	= new hcFloat[numLevels];
	hcFloat upperBound	= 2.5*r_sol;
	hcFloat lowerBound	= r_sol;
	hcFloat dr			= (upperBound - lowerBound) / (numLevels+1);

    for(uint i=0; i<numLevels; ++i)
    {
    	/*// ----------------------------------------------------------- first implementation
    	uint ind	= solver.grid->getIndex(i, 0, 0);
    	heights[i]	= solver.grid->getPos(ind, isEll)[0];//*/

    	//*// ----------------------------------------------------------- equal distances from lower to upper (exclusive)
    	heights[i] 	= lowerBound + (i+1)*dr;//*/

    	/*// ----------------------------------------------------------- only lower/upper levels
    	dr			= 0.05*r_sol;
    	heights[i] 	= lowerBound + (i+1)*dr;//*/
    }

    hcDate start, stop;
    start.setFromSystemTime();

    hcFloat lowerR	= solver.grid->lowerR;
	hcFloat upperR	= solver.grid->upperR;

    if(!mapHeightLevel(upperR, heights, numLevels, 14.5/15.0, numTheta, numPhi, sinLatFormat, true,  false)) return false;
    if(!mapHeightLevel(lowerR, heights, numLevels, 14.5/15.0, numTheta, numPhi, sinLatFormat, true,  false)) return false;

    if(!compCoords)
	if(!mapHeightLevel(upperR, heights, numLevels, 14.5/15.0, numTheta, numPhi, sinLatFormat, false, false)) return false;

	if(computeIntermediateHeightLevels)
		for(uint i=0; i<numLevels; ++i)
			if(!mapHeightLevel(heights[i], heights, numLevels, 14.5/15.0, numTheta, numPhi, sinLatFormat, compCoords, false)) return false;

    delete [] heights;

    stop.setFromSystemTime();
    uint seconds = (stop-start)/hcDate::facSec;
    printStdOutMess(__FILE__, __LINE__, "took " + to_string(seconds) + " s to compute all heights.");
    return true;
}

#ifdef GUI
bool PFSSsolution::insertFootpointData(const char *fn, FootpointData &data)
{
	if(!doesFileExist(fn))
		createFolderStructureTo(fn);

	hcSortedListStorage<FootpointData> list;

	std::ifstream stream(fn);
	string line;
	while(getline(stream, line))
		list.insertElement(*(new FootpointData(line)));
	stream.close();

	for(uint i=0; i<list.numElements; ++i)
		if(data == *(list.elements[i]))	list.removeElement(*list.elements[i]);

	list.insertElement(*(new FootpointData(data)));

	boost::filesystem::remove(fn);
	ofstream file(fn);
	for(uint i=0; i<list.numElements; ++i)
		file << list.elements[i]->toString() << "\n";
	file.close();

	return true;
}

bool PFSSsolution::EUVanalysis(	euvID id, hcImageFITS &euv, bool EIT, hcFloat latThresh,
								MagMapping &mapBack, MagMapping &mapForw, FootpointData &retval)
{
	if(!solver.solutionComputed)
	{
		cerr << __FILE__ << ":" << __LINE__ << ": Solution has not been computed. No analysis possible!\n";
		return false;
	}

	uint numX			= euv.width;
	uint numY			= euv.height;

	std::mt19937 rng(0);
	//LaplaceSolver &sol 	= solver;

	hcFloat maxLatitude	= (EIT ?  83 :  90) * 2*PI / 360;		// EIT images are linear in latitude as opposed to sine latitude
	hcFloat minLatitude	= (EIT ? -83 : -90) * 2*PI / 360;
	hcFloat tfX			= (numX-1)	/ (2*PI);
	hcFloat tfY			= (numY-1) 	/ (maxLatitude - minLatitude);

	string obs 			= EIT ? "EIT" : "AIA";
	string fnIMG		= getFilename_EUVimg(outDir, obs, id, info, latThresh);
	string fnFoot 		= getFilename_EUVfootpoints(outDir, obs, id, info, latThresh);
	string fnForwOpen	= getFilename_EUVforwOpen(outDir, obs, id, info, latThresh);
	string fnForwClose	= getFilename_EUVforwClose(outDir, obs, id, info, latThresh);

	std::uniform_int_distribution<int> genX(0, numX-1);
	std::uniform_int_distribution<int> genY(0, numY-1);

	hcImageFITS imgEUVCopy(numX, numY);
	hcImageFITS imgFoot(numX, numY);
	hcImageFITS imgOpen(numX, numY);
	hcImageFITS imgClosed(numX, numY);

	// collect all pixels from EUV map that are identified to be footpoints
	// of magnetic field lines tracked down from source surface
	for(uint j=0; j<mapBack.numTheta; ++j)
		for(uint k=0; k<mapBack.numPhi; ++k)
		{
			Magline &magl 	= mapBack(j, k);
			hcFloat lat		= PI/2.0-magl.posdata[0][1];				// latitude of footpoint
			int x			= (uint)round(magl.posdata[0][2]*tfX);		// img position of footpoint
			int y			= (int)round((lat-minLatitude)*tfY);

			if(y>=0 && y<(int)numY && lat >= -latThresh && lat <= latThresh && fabs(euv(x,y) - 1E-20)>1E-6)
			{
				//foot(x,y) 	= euv(x,y);
				imgFoot(x,y)= euv(x,y);
			}
		}

	uint numNAN = 0;

	// 1 - create copy of EUV map within specified latitude threshold (euvCopy)
	// 2 - create EUV copies of forward mapped maglines from photosphere discerning between open and closed maglines (forwOpen, forwClosed)
	for(uint x=0; x<numX; ++x)
		for(uint y=0; y<numY; ++y)
		{
			Magline &magl 	= mapForw(numY-y-1, x);
			hcFloat lat		= PI/2.0-magl.posdata[0][1];				// latitude of footpoint

			if(lat >= -latThresh && lat <= latThresh)	imgEUVCopy(x,y) = euv(x,y);

			if(fabs(euv(x,y) - 1E-20)<1E-6)				++numNAN;	// number of invalid pixels independent of latitude threshold
			else if(lat >= -latThresh && lat <= latThresh && magl.valid)
			{
				if(magl.closed)	imgClosed(x,y)	=  euv(x,y);
				else			imgOpen(x,y)	=  euv(x,y);
			}
		}

	ImageStatistics statsEUVCopy	= imgEUVCopy.getImageStatistics();
	ImageStatistics statsImgFoot	= imgFoot.getImageStatistics();
	ImageStatistics statsImgOpen	= imgOpen.getImageStatistics();
	ImageStatistics statsImgClosed	= imgClosed.getImageStatistics();

	uint numRandRuns = 50;
	ImageStatistics stats[numRandRuns];
	ImageStatistics statsRandMean;

	for(uint n=0; n<numRandRuns; ++n)
		{
			hcImageFITS imgFootRand(numX, numY);

			for(uint k=0; k<statsImgFoot.numPixels; ++k)	// take the same number of random pixels as number of footpoints in the backmapping
			{
				int xr,yr;
				hcFloat latr;

				do
				{
					xr		= genX(rng);
					yr		= genY(rng);
					latr	= yr/tfY + minLatitude;
				}while(fabs(euv(xr,yr)-1E-20)<1E-6 && latr >= -latThresh && latr <= latThresh);
				imgFootRand(xr,yr)	= euv(xr,yr);
			}
			stats[n] = imgFootRand.getImageStatistics();

			if(n==numRandRuns-1)
			{
				for(uint i=0; i<numRandRuns; ++i)
					statsRandMean += stats[i];

				statsRandMean /= numRandRuns;
			}
		}

	hcFloat ratioOC		= statsImgOpen.mean/statsImgClosed.mean;
	hcFloat ratioRand	= statsImgFoot.mean/statsRandMean.mean;

	imgEUVCopy.save(fnIMG.data());
	imgFoot.save(fnFoot.data());
	imgOpen.save(fnForwOpen.data());
	imgClosed.save(fnForwClose.data());

	SphericalGrid &gr = *solver.grid;
	/*
	hcFloat ellipticity = 1.0;
	if(gr.isElliptical())
	{
		EllipticalGrid &egr = *(EllipticalGrid *)(&gr);
		ellipticity			= egr.getEllA(egr.numR-1);
	}//*/

	retval 		= FootpointData(info, statsEUVCopy, statsImgFoot, statsImgOpen, statsImgClosed, statsRandMean,
								mapBack.numTheta, mapBack.numPhi, gr.maxSinLat, latThresh,
								id, euv.height, euv.width, numNAN, ratioRand, ratioOC);

	return true;
}


void rescaleAnalysis(const string &fn)
{
	EUVdata *euvOrig	= new EUVdata(fn);
	EUVdata *euvRescale	= new EUVdata(fn);

	if(euvOrig == NULL)
		return;

	euvRescale->data.rescale(360,167);

	uint numTo			= euvOrig->data.height;
	uint numTr			= euvRescale->data.height;
	uint numPo			= euvOrig->data.width;
	uint numPr			= euvRescale->data.width;

	Vec2D *coordsO		= new Vec2D[numTo*numPo];
	Vec2D *coordsR		= new Vec2D[numTr*numPr];

	hcFloat maxSinLat	= 1.0;
	hcFloat maxLatitude	= PI/2.0;
	hcFloat minLatitude	= -PI/2.0;
	hcFloat dLatO		= (asin(maxSinLat)-asin(-maxSinLat)) / (numTo-1);
	hcFloat dLatR		= (asin(maxSinLat)-asin(-maxSinLat)) / (numTr-1);

	for(uint j=0; j<numTo; ++j)
		for(uint k=0; k<numPo; ++k)
		{
			hcFloat theta		= PI/2-asin(maxSinLat) + j*dLatO;
			hcFloat phi			= (k + 1.0/2.0) * 2 * PI / numPo;
			coordsO[j*numPo+k]	= Vec2D(theta, phi);
		}

	for(uint j=0; j<numTr; ++j)
		for(uint k=0; k<numPr; ++k)
		{
			hcFloat theta		= PI/2-asin(maxSinLat) + j*dLatR;
			hcFloat phi			= (k + 1.0/2.0) * 2 * PI / numPr;
			coordsR[j*numPr+k]	= Vec2D(theta, phi);
		}

	//hcFloat tfXo		= (numPo-1)	/ (2*PI);
	//hcFloat tfYo		= (numTo-1)	/ (maxLatitude - minLatitude);
	hcFloat tfXr		= (numPr-1)	/ (2*PI);
	hcFloat tfYr		= (numTr-1)	/ (maxLatitude - minLatitude);

	hcFloat **pixels	= new hcFloat*[numTr*numPr];
	uint *numPixels		= new uint[numTr*numPr];

	for(uint x=0; x<numPr; ++x)
		for(uint y=0; y<numTr; ++y)
		{
			uint ind		= y*numPr+x;
			pixels[ind] 	= new hcFloat[200];
			numPixels[ind]	= 0;
		}

	for(uint jo=0; jo<numTo; ++jo)
		for(uint ko=0; ko<numPo; ++ko)
		{
			hcFloat lat	= PI/2.0-coordsO[jo*numPo+ko][0];				// latitude of footpoint
			hcFloat lon	= coordsO[jo*numPo+ko][1];
			int x		= (uint)round(lon*tfXr);						// img position of footpoint
			int y		= (int)round((lat-minLatitude)*tfYr);

			uint ind	= y*numPr+x;
			pixels[ind][numPixels[ind]++] = euvOrig->data(ko,jo);
		}

	hcImageFITS numImg, diff, onePixel;
	onePixel.init(10,10);
	numImg.init(numPr,numTr);
	diff.init(numPr,numTr);

	for(uint x=0; x<10; ++x)
		for(uint y=0; y<10; ++y)
			onePixel(x,y) = 0.0/0.0;

	for(uint x=0; x<numPr; ++x)
		for(uint y=0; y<numTr; ++y)
		{
			uint ind		= y*numPr+x;
			numImg(x,y) 	= numPixels[ind];
			hcFloat mean, std;
			getMean(pixels[ind], numPixels[ind], mean, std);
			diff(x,y)		= euvRescale->data(x,y) - mean;

			if(x==3 && y==4)
			{
				printf("numPixels: %u\n", numPixels[ind]);
				for(uint t=0; t<numPixels[ind]; ++t)
				{
					uint yi	= t/10;
					uint xi = t%10;
					onePixel(xi,yi) = pixels[ind][t];
				}
			}
		}

	onePixel.save("onePixel.fits");
	euvRescale->data.save("rescaled.fits");
	numImg.save("numImg.fits");
	diff.save("diffMean-Rescaled.fits");
}


bool PFSSsolution::footpointAnalysis(hcFloat latThresh)
{
	cout << "-----------------------------------------------------------------------------------------------------\n";
	cout << "PFSSsolution::footpoint analysis:\n";
	cout << "latThresh:  " << latThresh 	<< "\n";
	cout << "-----------------------------------------------------------------------------------------------------\n";
	fflush(stdout);

	if(!solver.solutionComputed)
	{
		cerr << __FILE__ << ":" << __LINE__ << ": Solution has not been computed. No analysis possible\n";
		return false;
	}

	uint cr				= info.CRnum;
	bool EIT 			= cr <= 2055;
	uint numX			= 360;
	uint numY			= 167;
	uint numTheta		= 200;
	uint numPhi			= 400;

	stringstream fnSummary, fnMapBack, fnMapForw;
	fnMapBack	<< outDir << "/" << cr << "/" << getFilename_magMappingBin(info, info.rss, 	true,	true, 	numTheta, 	numPhi);
	fnMapForw	<< outDir << "/" << cr << "/" << getFilename_magMappingBin(info, r_sol, 	false, 	true,	numY, 		numX);
	fnSummary 	<< outDir << "/" << cr << "/" << getFilename_EUVfootpointSummary(info, latThresh);

	if(doesFileExist(fnSummary.str()))
	{
		cout << "Footpoint analysis file does already exist:\n";
		cout << fnSummary.str() << "\n";
		return false;
	}

	if(cr < 1916 || cr > 2186 || (cr > 2055 && cr < 2097))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": No EUV data for CR=" << cr << "\n";
		return false;
	}

	string fn171 	= getFilename_EUVdata(cr, W171);
	string fn193 	= getFilename_EUVdata(cr, W193);
	string fn195 	= getFilename_EUVdata(cr, W195);
	string fn284 	= getFilename_EUVdata(cr, W284);
	string fn304 	= getFilename_EUVdata(cr, W304);

	EUVdata *w171	= new EUVdata(fn171);
	EUVdata *w193	= EIT ? NULL 				: new EUVdata(fn193);
	EUVdata *w195	= EIT ? new EUVdata(fn195) 	: NULL;
	EUVdata *w284	= EIT ? new EUVdata(fn284) 	: NULL;
	EUVdata *w304	= new EUVdata(fn304);

	if(w171->crNum == 0 || (EIT && w195->crNum == 0) || (!EIT && w193->crNum == 0) || (EIT && w284->crNum == 0) || w304->crNum == 0)
	{
		cerr << __FILE__ << ":" << __LINE__ << ": EUV data for rotation " << cr << " does not exist. Abort footpoint analysis.\n";
		return false;
	}

	w171->data.rescale(numX, numY);
	w304->data.rescale(numX, numY);
	if(EIT)
	{
		w284->data.rescale(numX, numY);
		w195->data.rescale(numX, numY);
	}
	else
		w193->data.rescale(numX, numY);

	for(uint x=0;x<numX; ++x)
		for(uint y=0;y<numY; ++y)
		{
			if(			w171->data(x,y) < 0.0)	w171->data(x,y) = 1E-20;
			if(!EIT && 	w193->data(x,y) < 0.0) 	w193->data(x,y) = 1E-20;
			if( EIT && 	w195->data(x,y) < 0.0) 	w195->data(x,y) = 1E-20;
			if( EIT && 	w284->data(x,y) < 0.0) 	w284->data(x,y) = 1E-20;
			if(			w304->data(x,y) < 0.0)	w304->data(x,y) = 1E-20;
		}

	hcFloat max171, max193, max195, max284, max304;
				max171	= w171->data.getHistogramMaximum(40);
	if(!EIT) 	max193	= w193->data.getHistogramMaximum(40);
	if(EIT)  	max195	= w195->data.getHistogramMaximum(40);
	if(EIT)  	max284	= w284->data.getHistogramMaximum(40);
				max304	= w304->data.getHistogramMaximum(40);

	for(uint x=0; x<numX; ++x)
		for(uint y=0; y<numY; ++y)
		{
						w171->data(x,y)	/= max171;
			if(!EIT) 	w193->data(x,y)	/= max193;
			if(EIT)	 	w195->data(x,y)	/= max195;
			if(EIT)  	w284->data(x,y)	/= max284;
						w304->data(x,y)	/= max304;
		}

	std::mt19937 rng(0);
	LaplaceSolver &sol = solver;

	hcFloat maxLatitude	= (EIT ?  83 :  90) * 2*PI / 360;				// EIT/AIA images are linear in latitude as opposed to sine latitude
	//hcFloat minLatitude	= (EIT ? -83 : -90) * 2*PI / 360;

	MagMapping map, mapForw;

	// create backmapping from source surface down to photosphere
	if(!doesFileExist(fnMapBack.str()))
	{
		map.createAtHeight(info, sol, sol.grid->upperR, numTheta, numPhi, sol.grid->maxSinLat, true, true);
		if(sol.grid->isElliptical())
		{
			EllipticalGrid* eGrid = (EllipticalGrid*)sol.grid;
			eGrid->convertMagMapping(map);
		}
		map.exportBinary(fnMapBack.str().data());
	}
	else
		map.importBinary(fnMapBack.str().data());

	// create forward mapping from photosphere to source surface
	if(!doesFileExist(fnMapForw.str()))
	{
		mapForw.createAtHeight(info, sol, sol.grid->lowerR, numY, numX, sin(maxLatitude), false, true);
		if(sol.grid->isElliptical())
		{
			EllipticalGrid* eGrid = (EllipticalGrid*)sol.grid;
			eGrid->convertMagMapping(mapForw);
		}
		mapForw.exportBinary(fnMapForw.str().data());
	}
	else
		mapForw.importBinary(fnMapForw.str().data());

	cout << "PFSSsolution::footpointAnalysis: Mapping done, starting analysis...\n";
	fflush(stdout);

	FootpointData ft171, ft193, ft195, ft284, ft304;
	EUVanalysis(W171, w171->data, EIT, latThresh, map, mapForw, ft171);
	EUVanalysis(W304, w304->data, EIT, latThresh, map, mapForw, ft304);
	if(EIT)
	{
		EUVanalysis(W195, w195->data, EIT, latThresh, map, mapForw, ft195);
		EUVanalysis(W284, w284->data, EIT, latThresh, map, mapForw, ft284);
	}
	else
	{
		EUVanalysis(W193, w193->data, EIT, latThresh, map, mapForw, ft193);
	}

	hcDate now;
	now.setFromSystemTime();
	string instr		= getStringFromOriginID(info.instrument);
	string timestamp 	= now.toString();

	FILE *summary = fopen(fnSummary.str().data(), "w");
	fprintf(summary, "%s footpoint summary file written at %s\n\n", EIT ? "EIT" : "AIA", timestamp.data());
	fprintf(summary, "Instrument: %s\n", instr);
	fprintf(summary, "Dimensions of solver:\t%u x %u x %u\n", sol.grid->numR, sol.grid->numTheta, sol.grid->numPhi);
	fprintf(summary, "MaxSinLat and MinSinLat of solver: %E / %E\n", sol.grid->maxSinLat, sol.grid->minSinLat);
	fprintf(summary, "Dimensions of magnetic mapping: %u x %u\n", map.numPhi, map.numTheta);
	fprintf(summary, "\tmaxSinLat: %E\n\tminSinLat:: %E\n", sol.grid->maxSinLat, sol.grid->minSinLat);
	fprintf(summary, "Image\twidth x height\t# inv. pixels\tmean value\t\t\tmedian value\n");
	fprintf(summary, "171\t\t%u x %u\t\t%u\t\t\t%E\t\t\t%E\n", w171->data.width, w171->data.height, ft171.euvNumInvalidPixels, ft171.statsEUVcopy.mean, ft171.statsEUVcopy.perc.perc50);
	if(EIT)
	{
		fprintf(summary, "193\t\tN/A\n");
		fprintf(summary, "195\t\t%u x %u\t\t%u\t\t\t%E\t\t\t%E\n", w195->data.width, w195->data.height, ft195.euvNumInvalidPixels, ft195.statsEUVcopy.mean, ft195.statsEUVcopy.perc.perc50);
		fprintf(summary, "284\t\t%u x %u\t\t%u\t\t\t%E\t\t\t%E\n", w284->data.width, w284->data.height, ft284.euvNumInvalidPixels, ft284.statsEUVcopy.mean, ft284.statsEUVcopy.perc.perc50);
	}
	else
	{
		fprintf(summary, "193\t\t%u x %u\t\t%u\t\t\t%E\t\t\t%E\n", w193->data.width, w193->data.height, ft193.euvNumInvalidPixels, ft193.statsEUVcopy.mean, ft193.statsEUVcopy.perc.perc50);
		fprintf(summary, "195\t\tN/A\n");
		fprintf(summary, "284\t\tN/A\n");
	}
	fprintf(summary, "304\t\t%u x %u\t\t%u\t\t\t%E\t\t\t%E\n", w304->data.width, w304->data.height, ft304.euvNumInvalidPixels, ft304.statsEUVcopy.mean, ft304.statsEUVcopy.perc.perc50);
	fprintf(summary, "Consider footpoints between latitudes %E and %E\n", -latThresh, latThresh);
	fprintf(summary, "wavelength\t# footp.\tmean bright. footp.\tmean bright. rand\tratio(back/rand)\t# open magl. (forw)\t# closed magl. (forw.)\tave. bright. open magl.\tave. bright. closed magl.\tratio(open/closed)\n");
	fprintf(summary, "171\t\t\t%u\t\t%E\t\t%E\t\t%E\t\t\t%u\t\t\t\t%u\t\t\t\t%E\t\t\t\t%E\t\t\t\t%E\n",
			ft171.statsFoot.numPixels, ft171.statsFoot.mean, ft171.statsRandMean.mean, ft171.ratioRand, ft171.statsOpen.numPixels,
			ft171.statsClosed.numPixels, ft171.statsOpen.mean, ft171.statsClosed.mean, ft171.ratioOC);
	if(EIT)
	{
		fprintf(summary, "193\t\t\tN/A\n");
		fprintf(summary, "195\t\t\t%u\t\t%E\t\t%E\t\t%E\t\t\t%u\t\t\t\t%u\t\t\t\t%E\t\t\t\t%E\t\t\t\t%E\n",
				ft195.statsFoot.numPixels, ft195.statsFoot.mean, ft195.statsRandMean.mean, ft195.ratioRand, ft195.statsOpen.numPixels,
				ft195.statsClosed.numPixels, ft195.statsOpen.mean, ft195.statsClosed.mean, ft195.ratioOC);
		fprintf(summary, "284\t\t\t%u\t\t%E\t\t%E\t\t%E\t\t\t%u\t\t\t\t%u\t\t\t\t%E\t\t\t\t%E\t\t\t\t%E\n",
				ft284.statsFoot.numPixels, ft284.statsFoot.mean, ft284.statsRandMean.mean, ft284.ratioRand, ft284.statsOpen.numPixels,
				ft284.statsClosed.numPixels, ft284.statsOpen.mean, ft284.statsClosed.mean, ft284.ratioOC);
	}
	else
	{
		fprintf(summary, "193\t\t\t%u\t\t%E\t\t%E\t\t%E\t\t\t%u\t\t\t\t%u\t\t\t\t%E\t\t\t\t%E\t\t\t\t%E\n",
				ft193.statsFoot.numPixels, ft193.statsFoot.mean, ft193.statsRandMean.mean, ft193.ratioRand, ft193.statsOpen.numPixels,
				ft193.statsClosed.numPixels, ft193.statsOpen.mean, ft193.statsClosed.mean, ft193.ratioOC);
		fprintf(summary, "195\t\t\tN/A\n");
		fprintf(summary, "284\t\t\tN/A\n");
	}
	fprintf(summary, "304\t\t\t%u\t\t%E\t\t%E\t\t%E\t\t\t%u\t\t\t\t%u\t\t\t\t%E\t\t\t\t%E\t\t\t\t%E\n",
			ft304.statsFoot.numPixels, ft304.statsFoot.mean, ft304.statsRandMean.mean, ft304.ratioRand, ft304.statsOpen.numPixels,
			ft304.statsClosed.numPixels, ft304.statsOpen.mean, ft304.statsClosed.mean, ft304.ratioOC);
	fclose(summary);

	char fnData[1000];
	sprintf(fnData, "%s/footpointData", outDir);
	insertFootpointData(fnData, ft171);
	insertFootpointData(fnData, ft304);
	if(EIT)
	{
		insertFootpointData(fnData, ft195);
		insertFootpointData(fnData, ft284);
	}
	else
		insertFootpointData(fnData, ft193);

	delete w171;
	delete w193;
	delete w195;
	delete w284;
	delete w304;
	return true;
}
#endif

/*!	Even though the magnetic field is to be computed analytically employing the spherical harmonic coefficient approach,
 * 	further processing (like mapping of magnetic field lines from source surface down to the photosphere) is done on a
 * 	3D grid.
 *
 * 	@param 	filename		synoptic photospheric magnetogram to be loaded
 * 	@param	order			maximum order of coefficients to be computed
 * 	@param	r_ss			heliocentric distance of source surface (in m)
 * 	@param	compResR		radial resolution of computational grid
 * 	@param	compResTheta	meridional resolutions of computational grid
 * 	@param	compResPhi		zonal resolution of computational grid
 * 	@param	rDist			distribution of radial grid shells (use 0 for equidistant spacing in r-direction)
 * 	@param	geomIncFactor	geometric increment factor for non-equidistant r-spacing
 * 	@param	mapResTheta		meridional resolution of magnetic mapping
 * 	@param	mapResPhi		zonal resolution of magnetic mapping
 * 	@param	mapIntermediateHeights	true, if for each r-shell of the comp grid a magnetic mapping is to be computed
 *
 * 	@return true, if computation of PFSS and mapping successful
 *
 */
bool PFSSsolution::computeAndMapKielSHC(const string &filename, uint order, hcFloat r_ss,
		uint compResR, uint mapResTheta, uint mapResPhi, bool mapIntermediateHeights)
{
	if(!computeKielSHC(filename, order, r_ss, compResR)) 						return false;
	if(!multiMapSolution(mapResTheta, mapResPhi, mapIntermediateHeights, true, true)) return false;
	return true;
}

/*!
 *
 *	@param	filename		filename of the photospheric magnetogram to be analyzed
 *	@param	optionalID		optional identifier to be used in output filenames
 *	@param	r_ss			heliocentric position of surce surface
 *	@param	compResR		radial resolution of computational grid
 * 	@param	compResTheta	meridional resolutions of computational grid
 * 	@param	compResPhi		zonal resolution of computational grid
 * 	@param	rDist			distribution of radial grid shells (use 0 for equidistant spacing in r-direction)
 * 	@param	geomIncFactor	geometric increment factor for non-equidistant r-spacing
 * 	@param	mapResTheta		meridional resolution of magnetic mapping
 * 	@param	mapResPhi		zonal resolution of magnetic mapping
 * 	@param	mapIntermediateHeights	true, if for each r-shell of the comp grid a magnetic mapping is to be computed
 * 	@param	surfaceShape	0 - spherical grid, 1 - elliptical grid
 * 	@param	surfShapeA		stretching parameter for elliptical grid
 * 	@param	scaleMethod		rescaling algorithm for photospheric data
 */
bool PFSSsolution::computeAndMapKielGrid(
		const string &filename, const char *optionalID, hcFloat r_ss,
		uint compResR, 		uint mapResTheta, uint mapResPhi, 	bool mapIntermediateHeights,
		hcFloat ellipticity, uint scaleMethod)
{
	if(!computeKielGrid(filename, optionalID, r_ss, compResR, ellipticity, scaleMethod))	return false;
	if(!multiMapSolution(mapResTheta, mapResPhi, mapIntermediateHeights, true, true))		return false;
	return true;
}

/*!	@param	filename		configuration file to be loaded
 * 	@param	mapResTheta		meridional resolution of magnetic mapping
 * 	@param	mapResPhi		zonal resolution of magnetic mapping
 * 	@param	mapIntermediateHeights	true, if for each r-shell of the comp grid a magnetic mapping is to be computed *
 */
bool PFSSsolution::loadAndMapKielGrid(const string &filename, uint mapResTheta, uint mapResPhi, bool mapIntermediateHeights)
{
	if(!load(filename)) 																return false;
	if(!multiMapSolution(mapResTheta, mapResPhi, mapIntermediateHeights, true, true)) 	return false;
	return true;
}

/*! A computational grid is necessary for mapping magnetic field lines
 *
 * 	@param	photFilename			Photospheric synoptic magnetogram
 * 	@param	StanfordCoeffFilename	Spherical harmonic coefficients to be evaluated
 * 	@param	compResR				radial resolution of computational grid
 * 	@param	compResTheta			meridional resolutions of computational grid
 * 	@param	compResPhi				zonal resolution of computational grid
 * 	@param	rDist					distribution of radial grid shells (use 0 for equidistant spacing in r-direction)
 * 	@param	geomIncFactor			geometric increment factor for non-equidistant r-spacing
 * 	@param	mapResTheta				meridional resolution of magnetic mapping
 * 	@param	mapResPhi				zonal resolution of magnetic mapping
 * 	@param	mapIntermediateHeights	true, if for each r-shell of the comp grid a magnetic mapping is to be computed
 */
bool PFSSsolution::loadAndMapStanfordSHC(const char *photFilename, const char *StanfordCoeffFilename,
		uint compResR, uint mapResTheta, uint mapResPhi, bool mapIntermediateHeights)
{
	hcFloat r_ss = 2.5*r_sol;
	if (!doesFileExist(photFilename) || !doesFileExist(StanfordCoeffFilename))
	{
		printf("ERROR!\tPFSSsoluion::loadAndMapStanfordSHC\n\t input file '%s' does not exist or is directory\n", photFilename);
		return false;
	}

	if(!isDirSet()) return false;

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::loadAndMapStanfordSHC:\n-- photFile:\t'%s'\n-- coeffFile:\t '%s'\n", photFilename, StanfordCoeffFilename);
	printf("--------------------------------------------------------------------\n\n");

	SynopticInfo synInf;
	uint numT, numP;
	hcFloat geomFac;
	SphericalGrid::getOptimalGridParameters(compResR, r_sol, r_ss, geomFac, numT, numP);
	info.init(synInf, sizeof(hcFloat), MODEL_PFSS, METH_SHC, GROUP_STANFORD, r_ss, 0.0, 9, compResR, numT, numP);

	if(!loadPhotBoundary(photFilename))	return false;

	PFSSsolution_SHC_sun coeffMagfield;
	coeffMagfield.importCoefficients(StanfordCoeffFilename, true, 2.5*r_sol);

	string coeffFN = dirData + "/" + to_string(info.CRnum) + "/" + getFilename_harmonicCoeff(info);
	coeffMagfield.exportCoefficients(coeffFN.data());

	for(uint r=0;r<solver.grid->numR;++r)
		for(uint t=0;t<solver.grid->numTheta;++t)
			for(uint p=0;p<solver.grid->numPhi;++p)
			{
				bool debug = false;
				uint ind = solver.grid->getIndex(r, t, p);
				Vec3D result = solver.grid->getB(ind);
				coeffMagfield.eval(solver.grid->getPos(ind, false), result);
				solver.grid->setB(ind, result);
			}

	solver.solutionComputed = true;
	save();
	multiMapSolution(mapResTheta, mapResPhi, mapIntermediateHeights, true, true);

	printf("--------------------------------------------------------------------------\n");
	printf("--- INFO: PFSSsolution::loadAndMapStanfordSHC concluded. \n");
	printf("--------------------------------------------------------------------------\n\n");

	return true;
}


void plotLine(hcImageRGBA &img, hcImageFloat &zbuf, hcImageBool &occ, const uint &color, const Vec3D &posStart, const Vec3D &posEnd)
{
	int width	= img.width;
	int height	= img.height;
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
			if(posZ < zbuf(x,y) && (occ(x,y) == true || posZ < 0.0))
			{
				zbuf(x,y) 	= posZ;
				img(x,y)	= color;
			}
	}
}

#ifdef GUI
//*
void createOcculterMap(hcImageBool &occ, hcImageRGBA &imager, hcFloat clipRadius, hcFloat sizeHor, hcFloat sizeVert,
		const hcDate &date, hcFloat occRadius, bool invert)
{
	Vec3D lookDir, up, right;
	//mager.getBaseVectors(lookDir, up, right, date, clipRadius);
	//Vec3D pos			= imager.getPos(date);
	//Vec3D target		= imager.getTarget(date);
	//hcCamera3D cam		= imager.getView(date, clipRadius);

	Vec3D pos		= Vec3D(30.0, 30.0, 30.0);
	Vec3D target 	= Vec3D(0.0, 0.0, 0.0);
	lookDir			= target - pos;
	up 				= Vec3D(0.0, 0.0, 1.0);
	right 			= Vec3D(1.0, 0.0, 0.0);
	hcCamera3D cam;
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

			Vec4D posPixHom	= Vec4D(posPix[0], posPix[1], posPix[2], (hcFloat)1.0);
			Vec4D projC2	= projSun * cam * posPixHom;
			projC2			/= projC2[3];

			hcFloat projX	= projC2[0];
			hcFloat projY	= projC2[1];

			occ(x,y)	= (sqrt(projX*projX + projY*projY) < occRadius ? (invert ? true : false) : (invert ? false : true));
		}
}//*/

bool drawImageThesis(hcImageRGBA &canvas, hcImageBool &occ, hcFloat sizeHor, hcFloat sizeVert,
		hcCamera3D &camIn, hcFloat radiusAroundTarget, string &inImage)
{
	//ObservationImage *obs = getObservation(date);


	//fnObservation	= inImage;
	//hcImageFITS obsImg(fnObservation.data());
	hcImageRGBA obsImg;
	cout << inImage << "\n";
	obsImg.load(inImage);
	obsImg.dump();fflush(stdout);

	//hcCamera3D cam							= getView(date, radiusAroundTarget);
	hcCamera3D cam							= camIn;
	//Vec3D pos								= getPos(date);
	Vec3D pos								= cam.pos;

	Vec3D lookDir 							= Vec3D(0.0, 0.0, 0.0) - pos;	lookDir.normalize();
	Vec3D up								= cam.up; 						up.normalize();
	Vec3D right								= lookDir.cp(up); 				right.normalize();

	//percentileDataStruct percentile 		= obsImg.getPercentiles();
	//hcFloat maxV							= percentile.perc90;
	//hcFloat minV							= percentile.perc05;

	//hcFloat colorScale 						= 255 / (maxV-minV);
	//hcFloat colorScale 						= 255;
	hcPerspectiveProjection3D &projImager	= cam.projection;

	hcFloat distCam							= (Vec3D(0.0, 0.0, 0.0) - pos).length();
	hcFloat near							= distCam - radiusAroundTarget;
	hcFloat tfX								= sizeHor /(canvas.width - 1);		// pixel -> world
	hcFloat tfY								= sizeVert/(canvas.height- 1);
	hcFloat tfx								= (obsImg.width - 1)/(2.0);			// normalized (-1,1) -> pixel
	hcFloat tfy								= (obsImg.height- 1)/(2.0);

	for(uint x=0; x<canvas.width; ++x)
		for(uint y=0; y<canvas.height; ++y)
		{
			hcFloat posX	= x*tfX-sizeHor/2.0;
			hcFloat posY	= y*tfY-sizeVert/2.0;
			hcFloat posZ	= near;
			Vec3D posPix	= pos + posZ*lookDir + posX*right + posY*up;	// cartesian world coordinates

			Vec4D posPixHom	= Vec4D(posPix[0], posPix[1], posPix[2], (hcFloat)1.0);
			Vec4D projC2	= projImager * cam * posPixHom;
			projC2			/= projC2[3];

			hcFloat projX	= projC2[0];
			hcFloat projY	= projC2[1];

			if(fabs(projX) < 1.0 && fabs(projY) < 1.0)
			{
				uint pixX			= (projX + 1.0) * tfx;
				uint pixY			= (projY + 1.0) * tfy;

				cout << "projX: " << projX 	<< "\n";
				cout << "projY: " << projY 	<< "\n";
				cout << "tfx:   " << tfx	<< "\n";
				cout << "tfy:   " << tfx	<< "\n";
				fflush(stdout);

				uint val			= obsImg(pixX, pixY);
				//uint col			= min((hcFloat)255.0, (max((hcFloat)(val-minV), (hcFloat)0.0))*colorScale);
				//uint color			= char2RGBA8(col, col, col, 255);
				uint color			= val;

				if(occ(x,y))
				{
					//if(sizeof(T) == sizeof(hcImageRGBA))
						canvas(x,y)	= color;
					//else
					//	canvas(x,y)	= val;
				}
 			}
		}

	return true;
}

bool PFSSsolution::createProjectedView(hcImageRGBA &corImg, uint width, uint height, uint numMaglX, uint numMaglY,
		hcImageRGBA &photImg, const hcObliquity &solarRotAxis, hcFloat clipRadius, const hcDate &date, bool forwardMap)
{
	cout << "PFSSsolution::createProjectedView started\n";
	fflush(stdout);

	hcFloat aspRatio	= (hcFloat)(width)/height;
	//Vec3D target		= corImg.getTarget(date);
	//Vec3D pos			= corImg.getPos(date);
	Vec3D pos;
	pos(0) = 15.0*r_sol;
	pos(1) = 0.0;
	pos(2) = 0.0;
	Vec3D target		= Vec3D(0.0, 0.0, 0.0);

	hcFloat distCam		= (target - pos).length();
	hcFloat fov_2		= asin(clipRadius/distCam);
	hcFloat d			= (distCam - clipRadius) * tan(fov_2);
	hcFloat dh			= d;
	hcFloat dv			= d/aspRatio;
	hcFloat near		= distCam - clipRadius;
	hcFloat far			= distCam + clipRadius;
	hcFloat sizeHor		= 2*dh;
	hcFloat sizeVert	= 2*dv;

	Matrix4x4 solarOrientation 	= solarRotAxis.getOrientation(date);
	hcPerspectiveProjection3D project(near, far, -dh, dh, -dv, dv);
	hcCamera3D cam(project, pos);
	//cam.setPos(pos);
	//cam.lookat(Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 0.0, 1.0));
	cam.setPos(pos);
	cam.lookat(Vec3D(0.0, 0.0, 0.0));
	//cam.projection.changeFOV(0.23);
	cam.update();

	hcImageRGBA img(width, height);
	hcImageBool occ(width, height);
	hcImageFloat zbuf(width, height);

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			img(x,y) 	= char2RGBA8(0,0,0,255);
			occ(x,y)	= false;
			zbuf(x,y)	= 1.0;
		}

	bool isElliptic 		= solver.grid->isElliptical();
	EllipticalGrid* eGrid 	= (EllipticalGrid*)solver.grid;
	hcFloat ellipticity		= isElliptic ? eGrid->getEllA(eGrid->getIndex(eGrid->numR-1, 0, 0)) : 1.0;

	//bool EIT 				= false;
	//hcFloat maxLatitude		= (EIT ?  83 :  90) * 2*PI / 360;
	hcFloat maxLatitude		= 75 * 2*PI / 360;


	hcFloat heightMap 		= forwardMap ? solver.grid->lowerR : solver.grid->upperR;

	MagMapping map;
	/*
	char fnDisk[1000];
	getMagneticMappingFilename(true, numMaglY, numMaglX, heightMap, fnDisk);
	string mapFN	= str(format("%1%.map") % fnDisk);//*/
	string oDir 	= outDir;
	stringstream mapFN;
	mapFN << outDir << "/" << info.CRnum << "/" << getFilename_magMappingBin(info, heightMap, true, true, numMaglY, numMaglX);

	if(!doesFileExist(mapFN.str().data()))
	{
		map.createAtHeight(info, solver, heightMap, numMaglY,	numMaglX, sin(maxLatitude), false, true);
		if(isElliptic)		eGrid->convertMagMapping(map);

		map.exportASCII(mapFN.str().data(), NULL, 0, solver.grid->lowerR, solver.grid->upperR);
		map.exportBinary(mapFN.str().data());
	}
	else
		if(!map.importBinary(mapFN.str().data()))
		{
			cerr << "PFSSsolution::createProjectedView: Magnetic map '" << mapFN.str() << "' does exist but could not be loaded!\n";
			return false;
		}

	cout << "Mapping created\n";
	map.dump();
	map.exportImage("test.bmp");


	hcDate now;	now.setFromSystemTime();
	stringstream dir, fn, fnImg, fnInfo, fnLOS, fnImgCor;
	dir		<< "../data/images/" << date.getCarringtonRotationNum() << "/";
	fn 		<< dir.str() << date.getCarringtonRotationNum() << "_" << now.toString() << "_" << "corrImg" << "-" << "photImg" << "_" << (forwardMap ? "forward" : "backward");
	fnImg	<< fn.str() << ".bmp";
	fnInfo	<< fn.str() << ".info";
	fnLOS	<< dir.str() << date.getCarringtonRotationNum() << "_" << now.toString() << "_magneticLOS.fits";
	fnImgCor<< dir.str() << date.getCarringtonRotationNum() << "_" << now.toString() << "_coronagraph.fits";

	// ---------------------------------------------------------------------------------------
	// 											Overlay Coronagraph and photospheric image
	// ---------------------------------------------------------------------------------------

	string eclipseFN = "../data/eclipse.png";

	//*
	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 2*clipRadius/r_sol, true);
	cout << "Draw eclipse\n";fflush(stdout);
	drawImageThesis(img, occ, sizeHor, sizeVert, cam, clipRadius, eclipseFN);//*/

	//*
	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 1.2, true);
	//photImg.drawImage(img, occ, sizeHor, sizeVert, date, clipRadius, photObsFN);//*/

	// ---------------------------------------------------------------------------------------
	// 																print projected magmap
	// ---------------------------------------------------------------------------------------

	//*
	MagMapping projMap = map;
	projMap.dump();

	for(uint x=0; x<map.numPhi; ++x)
		for(uint y=0; y<map.numTheta; ++y)
		{
			Magline &magl 		= map(y,x);
			Magline &projMagl	= projMap(y,x);

			for(uint i=0; i<magl.numSamples; ++i)
			{
				Vec3D pos 		= magl.posdata[i].convCoordSpher2Cart();
				Vec4D homPos	= Vec4D(pos[0], pos[1], pos[2], (hcFloat)1.0);
				Vec4D projPos	= cam.projection * cam * solarOrientation * homPos;
				projPos			/= projPos[3];
				Vec3D projPos3D = Vec3D(projPos[0], projPos[1], projPos[2]);

				projMagl.posdata[i] = projPos3D;

				/*
				cam.projection.dump();
				cam.dump();
				solarOrientation.dump();
				homPos.dump();//*/
			}
			//cout << x << "/" << y << " " << magl.posdata[0][0] << "/" << magl.posdata[0][1] << "/" << magl.posdata[0][2] << "\n";
			//cout << x << "/" << y << " " << projMagl.posdata[0][0] << "/" << projMagl.posdata[0][1] << "/" << projMagl.posdata[0][2] << "\n\n";
			//exit(1);
		}


	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 1.0, false);

	for(uint x=0; x<map.numPhi; ++x)
		for(uint y=0; y<map.numTheta; ++y)
		{
			Magline &magl	= projMap(y,x);
			uint color		= 	 !magl.valid 	? Magline::colorInvalid 	:
								(magl.closed 	? Magline::colorClosed 		:
								(!magl.polarity ? Magline::colorNegative	:
												  Magline::colorPositive));
			if(!magl.valid) continue;
			for(uint i=1; i<magl.numSamples; ++i)
			{
				Vec3D &posStart	= magl.posdata[i-1];
				Vec3D &posEnd	= magl.posdata[i];

				plotLine(img, zbuf, occ, color, posStart, posEnd);
			}
		}
	//*/


	// ---------------------------------------------------------------------------------------
	// 																			Draw the sun
	// ---------------------------------------------------------------------------------------
	/*

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
		}//*/

	createFolderStructureTo(fnImg.str().data());
	img.save(fnImg.str().data());

	//hcImageFITS corIMG(corObsFN.data());
	//hcImageFITS photIMG(photObsFN.data());

	ofstream info(fnInfo.str());
	info << "Filename:\t\t\t\t"				<< fnImg.str()		<< "\n";
	info << "\tHorizontal res.:\t"			<< img.width 		<< "\n";
	info << "\tVertical res.:\t\t"			<< img.height 		<< "\n";
	info << "\tcreated at:\t\t\t"			<< now.toString() 	<< "\n";
	info << "PFSS solution:\n";
	info << "\tEllipticity:\t\t"			<< ellipticity 		<< "\n";
	info << "\tnumR:\t\t\t\t"				<< solver.grid->numR<< "\n";
	info << "\tnumT:\t\t\t\t"				<< solver.grid->numTheta<< "\n";
	info << "\tnumP:\t\t\t\t"				<< solver.grid->numPhi	<< "\n";
	info << "Magnetic mapping:\t\t"			<< (forwardMap?"forward":"backward") << "\n";
	info << "\tZonal res.:\t\t\t" 			<< map.numPhi 		<< "\n";
	info << "\tMeridional res.:\t"			<< map.numTheta		<< "\n";
	//info << "Solar surface image:\t"		<< corObsFN 		<< "\n";
	//info << "\tHorizontal res.:\t"			<< corIMG.width		<< "\n";
	//info << "\tVertical res.:\t\t"			<< corIMG.height	<< "\n";
	//info << "Photospheric image:\t\t"		<< photObsFN		<< "\n";
	//info << "\tHorizontal res.:\t"			<< photIMG.width	<< "\n";
	//info << "\tVertical res.:\t\t"			<< photIMG.height	<< "\n";

	return true;
}

/*
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

			Vec4D posPixHom	= Vec4D(posPix[0], posPix[1], posPix[2], (hcFloat)1.0);
			Vec4D projC2	= projSun * cam * posPixHom;
			projC2			/= projC2[3];

			hcFloat projX	= projC2[0];
			hcFloat projY	= projC2[1];

			occ(x,y)	= (sqrt(projX*projX + projY*projY) < occRadius ? (invert ? true : false) : (invert ? false : true));
		}
}

bool PFSSsolution::createWhitelightCorona(Imager &corImg, uint width, uint height, hcFloat clipRadius,
		hcDate date, string fnLOS, string fnImgCor)
{
	//cout << "PFSSsolution::createWLcorona:\nwidth/height:\t" << width << "/" << height << "\n";
	//cout << "clipRad:\t" << clipRadius << "\ndate:\t" << date.toSpiceString() << "\n";

	hcFloat aspRatio	= (hcFloat)(width)/height;
	Vec3D target		= corImg.getTarget(date);
	Vec3D pos			= corImg.getPos(date);
	hcFloat distCam		= (target - pos).length();
	hcFloat fov_2		= asin(clipRadius/distCam);
	hcFloat d			= (distCam - clipRadius) * tan(fov_2);
	hcFloat dh			= d;
	hcFloat dv			= d/aspRatio;
	hcFloat near		= distCam - clipRadius;
	hcFloat far			= distCam + clipRadius;
	hcFloat sizeHor		= 2*dh;
	hcFloat sizeVert	= 2*dv;

	hcImageBool occ(width, height);
	hcImageFITS imgLOS(width,height), imgCor(width, height);

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			occ(x,y)	= false;
			imgLOS(x,y)	= 0.0;
			imgCor(x,y)	= 0.0;
		}

	string corObsFN;
	hcDate start, end;
	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 1.0, true);
	start.setFromSystemTime();
#ifdef NUMTHREADS
	corImg.drawMagneticLOS_MP(imgLOS, occ, sizeHor, sizeVert, date, clipRadius, *solver.grid);
#else
	corImg.drawMagneticLOS(imgLOS, occ, sizeHor, sizeVert, date, clipRadius, *solver.grid);
#endif
	end.setFromSystemTime();
	cout << "drawMagLOS took " << (uint)((end-start)/hcDate::facSec) << " s\n";
	fflush(stdout);
	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 2*clipRadius/r_sol, true);

	corImg.drawImage(imgCor, occ, sizeHor, sizeVert, date, clipRadius, corObsFN);

	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 1.0, false);
	undrawImage(imgCor, occ);

	imgLOS.normalize();
	imgCor.normalize();
	imgLOS.save(fnLOS.data());
	imgCor.save(fnImgCor.data());

	return true;
}


bool PFSSsolution::createProjectedView(Imager &corImg, uint width, uint height, uint numMaglX, uint numMaglY,
		Imager &photImg, const hcObliquity &solarRotAxis, hcFloat clipRadius, const hcDate &date, bool forwardMap)
{
	cout << "PFSSsolution::createProjectedView started\n";
	fflush(stdout);

	hcFloat aspRatio	= (hcFloat)(width)/height;
	Vec3D target		= corImg.getTarget(date);
	Vec3D pos			= corImg.getPos(date);
	hcFloat distCam		= (target - pos).length();
	hcFloat fov_2		= asin(clipRadius/distCam);
	hcFloat d			= (distCam - clipRadius) * tan(fov_2);
	hcFloat dh			= d;
	hcFloat dv			= d/aspRatio;
	hcFloat near		= distCam - clipRadius;
	hcFloat far			= distCam + clipRadius;
	hcFloat sizeHor		= 2*dh;
	hcFloat sizeVert	= 2*dv;

	Matrix4x4 solarOrientation 	= solarRotAxis.getOrientation(date);
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
			occ(x,y)	= false;
			zbuf(x,y)	= 1.0;
		}

	bool isElliptic 		= solver.grid->isElliptical();
	EllipticalGrid* eGrid 	= (EllipticalGrid*)solver.grid;
	hcFloat ellipticity		= isElliptic ? eGrid->getEllA(eGrid->getIndex(eGrid->numR-1, 0, 0)) : 1.0;

	bool EIT 				= false;
	hcFloat maxLatitude		= (EIT ?  83 :  90) * 2*PI / 360;
	hcFloat minLatitude		= (EIT ? -83 : -90) * 2*PI / 360;

	hcFloat heightMap 		= forwardMap ? solver.grid->lowerR : solver.grid->upperR;

	MagMapping map;
	/*
	char fnDisk[1000];
	getMagneticMappingFilename(true, numMaglY, numMaglX, heightMap, fnDisk);
	string mapFN	= str(format("%1%.map") % fnDisk);///
	string oDir 	= outDir;
	stringstream mapFN;
	mapFN << outDir << "/" << info.CRnum << "/" << getFilename_magMappingBin(info, heightMap, true, true, numMaglY, numMaglX);

	if(!doesFileExist(mapFN.str().data()))
	{
		map.createAtHeight(solver, heightMap, numMaglY,	numMaglX, sin(maxLatitude), false, true);
		if(isElliptic)
			eGrid->convertMagMapping(map);

		map.exportASCII(mapFN.str().data(), NULL, 0, solver.grid->lowerR, solver.grid->upperR);
		map.exportBinary(mapFN.str().data());
	}
	else
		if(!map.importBinary(mapFN.str().data()))
		{
			cerr << "PFSSsolution::createProjectedView: Magnetic map '" << mapFN.str() << "' does exist but could not be loaded!\n";
			return false;
		}

	hcDate now;	now.setFromSystemTime();
	stringstream dir, fn, fnImg, fnInfo, fnLOS, fnImgCor;
	dir		<< "../data/images/" << date.getCarringtonRotationNum() << "/";
	fn 		<< dir.str() << date.getCarringtonRotationNum() << "_" << now.toString() << "_" << getStringFromInstrumentID(corImg.id) << "-" << getStringFromInstrumentID(photImg.id) << "_" << (forwardMap ? "forward" : "backward");
	fnImg	<< fn.str() << ".bmp";
	fnInfo	<< fn.str() << ".info";
	fnLOS	<< dir.str() << date.getCarringtonRotationNum() << "_" << now.toString() << "_magneticLOS.fits";
	fnImgCor<< dir.str() << date.getCarringtonRotationNum() << "_" << now.toString() << "_coronagraph.fits";

	// ---------------------------------------------------------------------------------------
	// 											Overlay Coronagraph and photospheric image
	// ---------------------------------------------------------------------------------------

	string corObsFN;
	string photObsFN;

	//*
	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 2*clipRadius/r_sol, true);
	corImg.drawImage(img, occ, sizeHor, sizeVert, date, clipRadius, corObsFN);///

	//*
	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 1.2, true);
	photImg.drawImage(img, occ, sizeHor, sizeVert, date, clipRadius, photObsFN);//*/

	// ---------------------------------------------------------------------------------------
	// 																print projected magmap
	// ---------------------------------------------------------------------------------------

	/*
	MagMapping projMap = map;

	for(uint x=0; x<map.numPhi; ++x)
		for(uint y=0; y<map.numTheta; ++y)
		{
			Magline &magl 		= map(y,x);
			Magline &projMagl	= projMap(y,x);

			for(uint i=0; i<magl.numSamples; ++i)
			{
				Vec3D pos 		= magl.posdata[i].convCoordSpher2Cart();
				Vec4D homPos	= Vec4D(pos[0], pos[1], pos[2], (hcFloat)1.0);
				Vec4D projPos	= cam.projection * cam * solarOrientation * homPos;
				projPos			/= projPos[3];
				Vec3D projPos3D = Vec3D(projPos[0], projPos[1], projPos[2]);

				projMagl.posdata[i] = projPos3D;
			}
		}

	createOcculterMap(occ, corImg, clipRadius, sizeHor, sizeVert, date, 1.0, false);

	for(uint x=0; x<map.numPhi; ++x)
		for(uint y=0; y<map.numTheta; ++y)
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
	//*/


	// ---------------------------------------------------------------------------------------
	// 																			Draw the sun
	// ---------------------------------------------------------------------------------------
	/*

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

	createFolderStructureTo(fnImg.str().data());
	img.save(fnImg.str().data());

	hcImageFITS corIMG(corObsFN.data());
	hcImageFITS photIMG(photObsFN.data());

	ofstream info(fnInfo.str());
	info << "Filename:\t\t\t\t"				<< fnImg.str()		<< "\n";
	info << "\tHorizontal res.:\t"			<< img.width 		<< "\n";
	info << "\tVertical res.:\t\t"			<< img.height 		<< "\n";
	info << "\tcreated at:\t\t\t"			<< now.toString() 	<< "\n";
	info << "PFSS solution:\n";
	info << "\tEllipticity:\t\t"			<< ellipticity 		<< "\n";
	info << "\tnumR:\t\t\t\t"				<< solver.grid->numR<< "\n";
	info << "\tnumT:\t\t\t\t"				<< solver.grid->numTheta<< "\n";
	info << "\tnumP:\t\t\t\t"				<< solver.grid->numPhi	<< "\n";
	info << "Magnetic mapping:\t\t"			<< (forwardMap?"forward":"backward") << "\n";
	info << "\tZonal res.:\t\t\t" 			<< map.numPhi 		<< "\n";
	info << "\tMeridional res.:\t"			<< map.numTheta		<< "\n";
	info << "Solar surface image:\t"		<< corObsFN 		<< "\n";
	info << "\tHorizontal res.:\t"			<< corIMG.width		<< "\n";
	info << "\tVertical res.:\t\t"			<< corIMG.height	<< "\n";
	info << "Photospheric image:\t\t"		<< photObsFN		<< "\n";
	info << "\tHorizontal res.:\t"			<< photIMG.width	<< "\n";
	info << "\tVertical res.:\t\t"			<< photIMG.height	<< "\n";

	return true;
}//*/
#endif //GUI

void PFSSsolution::dump() const
{
	printStdOutMess(__FILE__, __LINE__, "dumping PFSSsolution:");
    info.dump();
}
