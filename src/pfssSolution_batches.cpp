#include "src/pfssSolution.h"
#include "src/filenames.h"

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/format.hpp"

#include <iomanip>

using namespace boost;
using namespace filesystem;

/*! @param	inDir					direcory to be worked upon
 * 	@param	r_ss					heliocenric position of source surface
 *	@param	compResR				radial resolution of computational grid
 * 	@param	compResTheta			meridional resolutions of computational grid
 * 	@param	compResPhi				zonal resolution of computational grid
 * 	@param	rDist					distribution of radial grid shells (use 0 for equidistant spacing in r-direction)
 * 	@param	geomIncFactor			geometric increment factor for non-equidistant r-spacing
 * 	@param	mapResTheta				meridional resolution of magnetic mapping
 * 	@param	mapResPhi				zonal resolution of magnetic mapping
 * 	@param	mapIntermediateHeights	true, if for each r-shell of the comp grid a magnetic mapping is to be computed
 */
bool PFSSsolution::batchKielSHC(const char *inDir, hcFloat r_ss,
		uint compResR, uint mapResTheta, uint mapResPhi, bool mapIntermediateHeights)
{
    if (!exists(inDir) || !is_directory(inDir))
    {
    	printf("ERROR!\tPFSSsoluion::batchKielSHC\n\t input directory '%s' does not exist\n", inDir);
        return false;
    }

    if (outDir[0] == '\0')
	{
    	printf("ERROR!\tPFSSsolution::batchKielSHC\n\t output directory has not been set!\n");
    	return false;
	}

	cout << "--------------------------------------------------------------------\n";
	cout << "-- PFSS::batchKielSHC on dir: '" << inDir << "'\n";
	cout << "-- r_ss:        " << r_ss /r_sol	<< " r_sol\n";
	cout << "-- compResR:    " << compResR 		<< "\n";
	cout << "-- mapResTheta: " << mapResTheta 	<< "\n";
	cout << "-- mapResPhi:   " << mapResPhi 	<< "\n";
	cout << "--------------------------------------------------------------------\n\n";
	fflush(stdout);

    path directory(inDir);
    directory_iterator end_iter;

	for (directory_iterator dir_iter(inDir) ; dir_iter != end_iter ; ++dir_iter)
        if (is_regular_file(dir_iter->status()))
        {
        	char filename[1000];
        	sprintf(filename, dir_iter->path().c_str());
            printf("Working on file\n'%s'\n\n", filename);
            computeAndMapKielSHC(filename, 9, r_ss, compResR, mapResTheta, mapResPhi, mapIntermediateHeights);
        }
        else
        {
            printf("ERROR! HelioMagfield::batchKielSHC: file %s is not a regular file!\n", dir_iter->path().c_str());
            continue;
        }

    printf("--------------------------------------------------------------------------\n");
    printf("--- INFO: PFSSsolution::batchKielSHC concluded. \n");
    printf("--------------------------------------------------------------------------\n\n");
    return true;
}

/*! opens inDir and computes PFSS model for each photospheric map as well as
 *  magnetic field line maps. This function only does the
 *  computational work, for the GUI version please see HelioMagfield::batchRotationComputing
 *
 *	@param	inDir					direcory to be worked upon
 *	@param	r_ss					heliocenric position of source surface
 *	@param	compResR				radial resolution of computational grid
 * 	@param	compResTheta			meridional resolutions of computational grid
 * 	@param	compResPhi				zonal resolution of computational grid
 * 	@param	rDist					distribution of radial grid shells (use 0 for equidistant spacing in r-direction)
 * 	@param	geomIncFactor			geometric increment factor for non-equidistant r-spacing
 * 	@param	mapResTheta				meridional resolution of magnetic mapping
 * 	@param	mapResPhi				zonal resolution of magnetic mapping
 * 	@param	mapIntermediateHeights	true, if for each r-shell of the comp grid a magnetic mapping is to be computed
 *
 */
bool PFSSsolution::batchKielGrid(const string &inDir, hcFloat r_ss,
		uint compResR, uint mapResTheta, uint mapResPhi, bool mapIntermediateHeights, hcFloat ellipticity)
{
    if (!exists(inDir) || !is_directory(inDir))
    {
    	printf("ERROR!\tPFSSsoluion::batchKielGrid\n\t input directory '%s' does not exist\n", inDir);
        return false;
    }

    if (outDir[0] == '\0')
	{
    	printf("ERROR!\tPFSSsolution::batchKielGrid\n\t output directory has not been set!\n");
    	return false;
	}

	cout << "--------------------------------------------------------------------\n";
	cout << "-- PFSS::batchKielGrid on dir: '" << inDir << "'\n";
	cout << "-- r_ss:        " << r_ss /r_sol	<< " r_sol\n";
	cout << "-- ellipticity: " << ellipticity	<< "\n";
	cout << "-- compResR:    " << compResR 		<< "\n";
	cout << "-- mapResTheta: " << mapResTheta 	<< "\n";
	cout << "-- mapResPhi:   " << mapResPhi 	<< "\n";
	cout << "--------------------------------------------------------------------\n\n";
	fflush(stdout);

	path directory(inDir);
	directory_iterator end_iter;

	for (directory_iterator dir_iter(inDir) ; dir_iter != end_iter ; ++dir_iter)
        if (is_regular_file(dir_iter->status()))
        {
        	string filename = dir_iter->path().string();
        	cout << "Working on file '" << filename << "'\n";fflush(stdout);
            //computeAndMapKielGrid(filename, "", r_ss, compResR, mapResTheta, mapResPhi, mapIntermediateHeights, ellipticity);

            computeKielGrid(filename, "", r_ss, compResR, ellipticity);
        }
        else
        {
            printf("ERROR! HelioMagfield::batchKielGrid: file %s is not a regular file!\n", dir_iter->path().c_str());
            continue;
        }

    printf("--------------------------------------------------------------------------\n");
    printf("--- PFSSsolution::batchKielGrid concluded. \n");
    printf("--------------------------------------------------------------------------\n\n");

    return true;
}

#ifdef GUI
bool PFSSsolution::batchFootpointAnalysis(	const string &inDir, hcFloat r_ss, uint compResR, methodID method,
											hcFloat ellipticity, uint orderSHC, hcFloat latThresh)
{
    if (!exists(inDir) || !is_directory(inDir))
    {
    	printf("ERROR!\tPFSSsoluion::batchFootpointAnalysis\n\t input directory '%s' does not exist\n", inDir);
        return false;
    }

    if (outDir[0] == '\0')
	{
    	printf("ERROR!\tPFSSsolution::batchFootpointAnalysis\n\t output directory has not been set!\n");
    	return false;
	}

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::batchFootpointAnalysis_KielGrid on dir: '%s'\n", inDir.data());
	printf("--------------------------------------------------------------------\n\n");

	path directory(inDir);
	directory_iterator end_iter;

	uint numberOfComputations 	= 0;
	uint count 					= 0;
	hcDate startTime;
	startTime.setFromSystemTime();

	for (directory_iterator dir_iter(inDir) ; dir_iter != end_iter ; ++dir_iter)
		if (is_regular_file(dir_iter->status()))
			++numberOfComputations;

	for (directory_iterator dir_iter(inDir) ; dir_iter != end_iter ; ++dir_iter)
	{
		string filename = dir_iter->path().c_str();

        if (is_regular_file(dir_iter->status()))
        {
            cout << "Working on file " << filename << " with parameters\n";
            cout << "r_ss:       " << r_ss/r_sol 	<< "\n";
            cout << "compResR:   " << compResR 		<< "\n";
            cout << "ell:        " << ellipticity 	<< "\n";
            cout << "orderSHC:   " << orderSHC 		<< "\n\n";

            if(!photBoundary.load(filename.data()))
            {
            	cerr << __FILE__ << ":" << __LINE__ << " Synoptic map " << filename << " cannot be loaded.\n";
            	return false;
            }
            modelID mod = MODEL_PFSS;
            uint numT, numP;
            hcFloat rfac;
            SphericalGrid::getOptimalGridParameters(compResR, r_sol, r_ss, rfac, numT, numP, false);
            PFSSsolutionInfo solInf(photBoundary.synInfo, sizeof(hcFloat), MODEL_PFSS, method, GROUP_KIEL, r_ss, ellipticity, 9, compResR, numT, numP);
            stringstream fpsummary;
            fpsummary << "../data/PFSS/output/" << photBoundary.synInfo.CRnum << "/" << getFilename_EUVfootpointSummary(solInf,latThresh);

            if(doesFileExist(fpsummary.str()))
            {
            	cout << "Footpoint analysis file " << fpsummary.str() << " does already exist. Skip recomputation.\n";
            	continue;
            }
            else
            	cout << "File " << fpsummary.str() << " does not exist.\n";

            uint cr = photBoundary.synInfo.CRnum;
            if(cr < 1916 || cr > 2186 || (cr > 2055 && cr < 2097))
            {
            	cout << "There is no EUV data for CR " << cr << "\n";
            	continue;
            }

            if(method==METH_NUMERIC || method==METH_ELLIPTICAL)
            	computeKielGrid(filename.data(), "", r_ss, compResR, ellipticity);
            else if(method==METH_SHC)
            	computeKielSHC(filename.data(), orderSHC, r_ss, 35);
            else
            {
            	cerr << __FILE__ << ":" << __LINE__ << ": method " << getStringFromMethodID(method) << " not supported for footpoint analysis\n";
            	return false;
            }
            cout << "PFSSsolution::batchFootpointAnalysis: " << getStringFromMethodID(info.method) << " " << info.maxSinLat << "\n";
			footpointAnalysis(latThresh);
			++count;

			computationTimeEstimate(numberOfComputations, count, startTime);
        }
        else
        {
            cerr << __FILE__ << ":" << __LINE__ << ": file " << dir_iter->path().c_str() << " is not a regular file.\n";
            continue;
        }
	}

    printf("--------------------------------------------------------------------------\n");
    printf("--- INFO: PFSSsolution::batchFootpointAnalysis_KielGrid concluded. \n");
    printf("--------------------------------------------------------------------------\n\n");

    return true;
}

bool PFSSsolution::batchFootpointAnalysisAllComputed(hcFloat latThresh)
{
	if(!exists(outDir))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": outDir not set or does not exist!\n";
		return false;
	}

	path directory(outDir);
	directory_iterator end_iter, endIter2;

	uint solutionCount = 0;

	for (directory_iterator dir_iter(outDir) ; dir_iter != end_iter ; ++dir_iter)
		if(is_directory(dir_iter->status()))
			for(directory_iterator innerIter(dir_iter->path().c_str()); innerIter != endIter2; ++innerIter)
			{
				string filename = innerIter->path().c_str();
				if(isFileType_pfssSolutionConfig(filename))
						++solutionCount;
			}

	hcDate startTime;
	startTime.setFromSystemTime();

	uint solutionsCompleted = 0;

	for (directory_iterator dir_iter(outDir) ; dir_iter != end_iter ; ++dir_iter)
		if(is_directory(dir_iter->status()))
			for(directory_iterator innerIter(dir_iter->path().c_str()); innerIter != endIter2; ++innerIter)
			{
				string filename = innerIter->path().c_str();
				if(isFileType_pfssSolutionConfig(filename))
				{
					++solutionsCompleted;
					cout << "Working on file " << filename << "\n"; fflush(stdout);
					PFSSsolutionInfo solInf;
					getParamFromFN_pfssSolutionInfo(filename, solInf);
					stringstream fnsummary;
					fnsummary << outDir << "/" << solInf.CRnum << "/" << getFilename_EUVfootpointSummary(solInf, INFINITY);

					uint cr = solInf.CRnum;
					if(cr < 1916 || cr > 2186 || (cr > 2055 && cr < 2097))
					{
						cout << "There is no EUV data for CR " << cr << "\n";
						continue;
					}

					if(doesFileExist(fnsummary.str()))
					{
						cout << "Footpoint analysis file " << fnsummary.str() << " already exists. Skip recomputation.\n";
						continue;
					}

					if(load(filename.data()))	footpointAnalysis(latThresh);

					computationTimeEstimate(solutionCount, solutionsCompleted, startTime);
				}
			}


	return true;
}

/*
bool PFSSsolution::batchCoronagraphAnalysis(hcSortedList<hcDate> &dates, const hcObliquity &solarRotAxis,
		uint compResR, uint mapResTheta, uint mapResPhi, Imager &coronagraph, Imager &euv, hcFloat ellipticity)
{
	cout <<"--------------------------------------------------------------------\n";
	cout << "-- PFSSsolution::batchCoronagraphAnalysis \n";
	cout << "--------------------------------------------------------------------\n\n";
	fflush(stdout);

	uint width 			= 100;
	uint height			= 50;
	uint numMaglX		= 100;
	uint numMaglY		= 50;
	hcFloat clipRadius	= 5*r_sol;

	for(uint i=0; i<dates.numElements; ++i)
	{
		hcDate date 		= *dates.elements[i];
		uint crNum			= date.getCarringtonRotationNum();
		ObservationImage *obsC1 	= coronagraph.getObservation(date);
		ObservationImage *obsEUV	= euv.getObservation(date);

		cout << "PFSSsolution::batchCoronagraphAnalysis: Analyze date " << date.toSpiceString() << "\n";

		if(obsC1 == NULL)
		{
			cout << "PFSSsolution::batchCoronagraphAnalysis: no observations for date " << date.toSpiceString() << "\n";
			cout << "obsC1:\t" << obsC1 << "\n";
			cout << "obsEUV:\t"<< obsEUV << "\n";
			continue;
		}

		stringstream magnetogramMDI;
		magnetogramMDI << "../data/photMaps/MDI/synop_Ml_0." << crNum << ".fits";

		if(doesFileExist(magnetogramMDI.str().data()))
			computeKielGrid(magnetogramMDI.str().data(), "", 2.5*r_sol, compResR, ellipticity);
		else
		{
			cout << "PFSSsolution::batchCoronagraphAnalysis: MDI map " << magnetogramMDI.str() << " does not exist\n";
			continue;
		}

		bool forwardMap = true;
		hcDate now;	now.setFromSystemTime();
		stringstream dir, fn, fnImg, fnInfo, fnLOS, fnImgCor;
		dir		<< "../data/images/" << date.getCarringtonRotationNum() << "/";
		fn 		<< dir.str() << crNum << "_" << now.toString() << "_" << getStringFromInstrumentID(coronagraph.id) << "-" << getStringFromInstrumentID(euv.id) << "_" << (forwardMap ? "forward" : "backward");
		fnImg	<< fn.str() << ".bmp";
		fnInfo	<< fn.str() << ".info";
		fnLOS	<< dir.str() << date.getCarringtonRotationNum() << "_" << now.toString() << "_magneticLOS.fits";
		fnImgCor<< dir.str() << date.getCarringtonRotationNum() << "_" << now.toString() << "_coronagraph.fits";

		createProjectedView(coronagraph, width, height, numMaglX, numMaglY, euv, solarRotAxis, clipRadius, date, forwardMap);
		createWhitelightCorona(coronagraph, width, height, clipRadius, date, fnLOS.str(), fnImgCor.str());
	}

	return true;
}//*/
#endif

/*! opens inDir and computes filters on magnetic photospheric maps
 *
 */
void PFSSsolution::filterPhotMaps(const char *inDir)
{
    path directory(inDir);
    directory_iterator end_iter;

    if (!exists(inDir) || !is_directory(inDir))
    {
    	printf("ERROR!\tPFSSsolution::filterPhotMaps\n\t input directory '%s' does not exist\n", inDir);
        return;
    }

    if (outDir[0] == '\0')
	{
    	printf("ERROR!\tPFSSsolution::filterPhotMaps\n\t output directory has not been set!\n");
    	return;
	}

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::filterPhotMaps in dir: '%s'\n", inDir);
	printf("--------------------------------------------------------------------\n\n");

	for (directory_iterator dir_iter(inDir) ; dir_iter != end_iter ; ++dir_iter)
    {
        if (!is_regular_file(dir_iter->status()))
        {
        	printf("ERROR! PFSSsolution::filterPhotMaps: file %s is not a regular file!\n", dir_iter->path().c_str());
			continue;
        }

		char filename[1000];
		sprintf(filename, dir_iter->path().c_str());

		char fn[1000];
		boost::regex reg(".*/(.*)\\.FITS", boost::regex::icase);
		boost::cmatch what;
		bool retval = boost::regex_search(filename, what, reg);
		sprintf(fn, what[1].str().c_str());

		hcImageFITS img(filename);
		hcImageFITS imgMean(img);
		hcImageFITS imgMedian(img);
		imgMean.meanFilter(5, true);
		imgMedian.medianFilter(5, true);

		char outFN[1000];
		sprintf(outFN, "%s/%s.FITS", outDir, fn);
		img.save(outFN);

		sprintf(outFN, "%s/%s_mean.FITS", outDir, fn);
		imgMean.save(outFN);

		sprintf(outFN, "%s/%s_median.FITS", outDir, fn);
		imgMedian.save(outFN);

		hcFloat meanDiff 	= imgMean.meanSquaredDiff(img);
		hcFloat medianDiff	= imgMedian.meanSquaredDiff(img);

		printf("%E\t%E\t%s\n", meanDiff, medianDiff, meanDiff < medianDiff ? "<" : ">-------");
    }

    printf("--------------------------------------------------------------------------\n");
    printf("--- INFO: PFSSsolution::filterPhotMaps concluded. \n");
    printf("--------------------------------------------------------------------------\n\n");
}

void PFSSsolution::compare2Stanford(const char *photFilename, const char *StanfordCoeffFilename,
		uint compRadialRes, uint compThetaRes, uint compPhiRes, bool rDistribution, hcFloat geomIncFactor,
		uint imgThetaRes,   uint imgPhiRes,    bool mapIntermediateHeights)
{
	if (!doesFileExist(photFilename) || !doesFileExist(StanfordCoeffFilename))
	{
		printf("ERROR!\tPFSSsoluion::compare2Stanford\n\t input file '%s' does not exist or is directory\n", photFilename);
		return;
	}

	if (outDir[0] == '\0')
	{
		printf("ERROR!\tPFSSsolution::compare2Stanford\n\t output directory has not been set!\n");
		return;
	}

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::compare2Stanford:\n-- photFile:\t'%s'\n-- coeffFile:\t '%s'\n", photFilename, StanfordCoeffFilename);
	printf("--------------------------------------------------------------------\n\n");

	computeAndMapKielSHC(photFilename, 9, 2.5*r_sol,				// Kiel SHC
			compRadialRes, imgThetaRes, imgPhiRes, mapIntermediateHeights);

	loadAndMapStanfordSHC(photFilename, StanfordCoeffFilename,		// Stanford SHC
			compRadialRes, imgThetaRes,   imgPhiRes,    mapIntermediateHeights);

	computeAndMapKielGrid(photFilename, "", 2.5*r_sol,				// Kiel grid
			compRadialRes, imgThetaRes, imgPhiRes, mapIntermediateHeights, 0.0);

	printf("--------------------------------------------------------------------------\n");
	printf("--- INFO: PFSSsolution::compare2Stanford concluded. \n");
	printf("--------------------------------------------------------------------------\n\n");
}

void PFSSsolution::compareSHCorders(const string &photFilename, hcFloat r_ss,
		uint compRadialRes, uint imgThetaRes, uint imgPhiRes, bool mapIntermediateHeights)
{
	if (!doesFileExist(photFilename.data()))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": File " << photFilename << " does not exist.\n";
		return;
	}

	if (outDir[0] == '\0')
	{
		cerr << __FILE__ << ":" << __LINE__ << ": output directory has not been set.\n";
		return;
	}

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::compareSHCorders on file: '%s'\n", photFilename.data());
	printf("--------------------------------------------------------------------\n\n");

	for(uint k=5; k<36; ++k)
		computeAndMapKielSHC(photFilename, k, r_ss, compRadialRes, imgThetaRes, imgPhiRes, mapIntermediateHeights);

	//computeAndMapKielGrid(filename, "", r_ss, compRadialRes, imgThetaRes, imgPhiRes, mapIntermediateHeights, 1.0);

	printf("--------------------------------------------------------------------------\n");
	printf("--- INFO: PFSSsolution::compareSHCorders concluded. \n");
	printf("--------------------------------------------------------------------------\n\n");
}

void PFSSsolution::compareRdist(const char *inDir, uint compRadialRes, uint compThetaRes, uint compPhiRes,
								uint imgThetaRes, uint imgPhiRes)
{
	path directory(inDir);
	directory_iterator end_iter;

	if (!exists(inDir) || !is_directory(inDir))
	{
		printf("ERROR!\tPFSSsoluion::paramStudyNoise\n\t input directory '%s' does not exist\n", inDir);
		return;
	}

	if (outDir[0] == '\0')
	{
		printf("ERROR!\tPFSSsolution::paramStudyNoise\n\t output directory has not been set!\n");
		return;
	}

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::compareRdist on dir: '%s'\n", inDir);
	printf("--------------------------------------------------------------------\n\n");

	char filename[1000];
	sprintf(filename, "%s/synop_Ml_0.2046.fits", inDir);

	printf("Working on file\n'%s'\n\n", filename);

	char idString[1000];
	hcFloat highestFac 	= 1.07;
	uint numRdists		= 10;

	for(uint l=0; l<numRdists; ++l)
	{
		hcFloat rFactor = 1.0 + (highestFac-1.0)/(numRdists-1) * l;
		printf("Compute rDist %E\n", rFactor);

		/*
		sprintf(idString, "rFactor_%E", rFactor);

		if(!loadPhotBoundary(filename, idString, MODEL_PFSS, METH_NUMERIC, GROUP_KIEL, 2.5*r_sol, compRadialRes))
		{
			printf("ERROR! PFSSsolution::computeAndMapKielSHC: File %s could not be loaded! Continue with next file!\n", filename);
			return;
		}//*/

		cout << "PFSSsolution::compareRdist needs optionalID for pfsssolutioninfo\n";
		exit(1);

		solver.computeSolution(photBoundary);
		save();
		multiMapSolution(imgThetaRes, imgPhiRes, false, true, true);

		SphericalGrid &gr	= *solver.grid;
		hcFloat upperR		= gr.upperR;
		hcFloat lowerR		= gr.lowerR;

		uint ind0			= gr.getIndex(0,0,0);
		uint ind1			= gr.getIndex(1,0,0);

		hcFloat dr 			= (1-highestFac) / (1-pow(highestFac, compRadialRes-1)) * (upperR - lowerR);
		hcFloat minLowerDist= dr * (1-pow(highestFac, 1)) / (1-highestFac);
		uint numR			= (upperR - lowerR) / minLowerDist;

		hcImageFITS	normedImage;
		normedImage.init(numR, compThetaRes);
		for(uint i=0; i<numR; ++i)
		{
			hcFloat r	= lowerR + i * (upperR - lowerR) / numR;
			uint k=0;
			while(++k<numR && gr.pos[gr.getIndex(k, 0, 0)][0] < r);

			hcFloat r_lo	= gr.pos[gr.getIndex(k-1, 0, 0)][0];
			hcFloat r_hi	= gr.pos[gr.getIndex(k  , 0, 0)][0];
			hcFloat frac	= (r_hi - r) / (r_hi - r_lo);

			for(uint j=0; j<compThetaRes; ++j)
			{
				hcFloat val_lo	= gr.psi[gr.getIndex(k-1, j, 0)];
				hcFloat val_hi	= gr.psi[gr.getIndex(k  , j, 0)];
				hcFloat val		= frac * val_lo + (1-frac) * val_hi;
				normedImage(i,j) = val;
			}
		}

		char fnRadImage[1000];
		sprintf(fnRadImage, "%s/1_radialNormedMap_%s.fits", outDir, idString);
		normedImage.save(fnRadImage);

		hcImageFITS	unnormedImage;
		unnormedImage.init(compRadialRes, compThetaRes);
		for(uint i=0; i<compRadialRes; ++i)
		{
			for(uint j=0; j<compThetaRes; ++j)
			{
				hcFloat val		= gr.psi[gr.getIndex(i, j, 0)];
				unnormedImage(i,j) = val;
			}
		}

		sprintf(fnRadImage, "%s/2_radialunNormedMap_%s.fits", outDir, idString);
		unnormedImage.save(fnRadImage);
	}

	for(uint k=9; k<10; ++k)
		computeAndMapKielSHC(filename, k, 2.5*r_sol, compRadialRes, imgThetaRes, imgPhiRes, false);


	printf("--------------------------------------------------------------------------\n");
	printf("--- INFO: PFSSsolution::compareRdist concluded. \n");
	printf("--------------------------------------------------------------------------\n\n");
}



void PFSSsolution::compareAnaInter(const char *inDir)
{
	uint compRadialRes		= 50;
	uint compThetaRes		= 100;
	uint compPhiRes			= 200;

	uint imgThetaRes		= 100;
	uint imgPhiRes			= 200;

	path directory(inDir);
	directory_iterator end_iter;

	if (!exists(inDir) || !is_directory(inDir))
	{
		printf("ERROR!\tPFSSsoluion::paramStudyNoise\n\t input directory '%s' does not exist\n", inDir);
		return;
	}

	if (outDir[0] == '\0')
	{
		printf("ERROR!\tPFSSsolution::paramStudyNoise\n\t output directory has not been set!\n");
		return;
	}

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::compareAnaInter on dir: '%s'\n", inDir);
	printf("--------------------------------------------------------------------\n\n");

	for (directory_iterator dir_iter(inDir) ; dir_iter != end_iter ; ++dir_iter)
	{
		if (is_regular_file(dir_iter->status()))
		{
			char filename[1000];
			sprintf(filename, dir_iter->path().c_str());

			printf("Working on file\n'%s'\n\n", filename);



			/*
			char idString[1000];
			char numbers[4];
			sprintf(numbers, "0%u", 9);
			sprintf(idString, "maxHarmonicOrder%s_", numbers);


			if(!loadPhotBoundary(filename, idString, MODEL_PFSS, METH_SHC, GROUP_KIEL, 2.5*r_sol,
					compRadialRes, compThetaRes, compPhiRes))
			{
				printf("ERROR! PFSSsolution::computeAndMapKielSHC: File %s could not be loaded! Continue with next file!\n", filename);
				return;
			}//*/

			cout << "PFSSsolution::compareAnaInter needs optionalID for pfsssolutioninfo\n";
			exit(1);

			printf("Computing coefficients for principal order %u\n", 9);

			stringstream coeffFN;
			coeffFN << outDir << "/" << info.CRnum << "/" << getFilename_harmonicCoeff(info);

			PFSSsolution_SHC_sun coeffMagfield;
			coeffMagfield.determineCoefficientsFromPhotMagfield(9, 2.5*r_sol, photBoundary);
			solver.solutionComputed = true;
			coeffMagfield.exportCoefficients(coeffFN.str().data());

			for(uint r=0;r<solver.grid->numR;++r)
				for(uint t=0;t<solver.grid->numTheta;++t)
					for(uint p=0;p<solver.grid->numPhi;++p)
					{
						bool debug 		= false;
						uint ind 		= solver.grid->getIndex(r, t, p);
						Vec3D result 	= solver.grid->getB(ind);
						coeffMagfield.eval(solver.grid->getPos(ind, false), result);
						solver.grid->setB(ind, result);
					}

			save();
			multiMapSolution(imgThetaRes, imgPhiRes, false, true, true);

			solver.grid->harmSolution = &coeffMagfield;
			save();
			multiMapSolution(imgThetaRes, imgPhiRes, false, true, true);
		}
		else
		{
			printf("ERROR! HelioMagfield::paramStudyNoise: file %s is not a regular file!\n", dir_iter->path().c_str());
			continue;
		}
	}

	printf("--------------------------------------------------------------------------\n");
	printf("--- INFO: PFSSsolution::compareAnaInter concluded. \n");
	printf("--------------------------------------------------------------------------\n\n");
}

/*! conducts a parameter study of noise in the magnetograms on files in directory inDir
 *
 */
void PFSSsolution::paramStudyNoise(const char *inDir, uint compRadialRes, uint compThetaRes, uint compPhiRes,
								uint imgThetaRes,  uint imgPhiRes, uint numNoiseComputations)
{
    path directory(inDir);
    directory_iterator end_iter;

    if (!exists(inDir) || !is_directory(inDir))
    {
    	printf("ERROR!\tPFSSsoluion::paramStudyNoise\n\t input directory '%s' does not exist\n", inDir);
        return;
    }

    if (outDir[0] == '\0')
	{
    	printf("ERROR!\tPFSSsolution::paramStudyNoise\n\t output directory has not been set!\n");
    	return;
	}

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::paramStudyNoise on dir: '%s'\n", inDir);
	printf("--------------------------------------------------------------------\n\n");

	for (directory_iterator dir_iter(inDir) ; dir_iter != end_iter ; ++dir_iter)
    {
        if (is_regular_file(dir_iter->status()))
        {
        	char filename[1000];
        	sprintf(filename, dir_iter->path().c_str());

            printf("Working on file\n'%s'\n\n", filename);

            for(uint i=0; i<=numNoiseComputations; ++i)
            {
            	char idString[1000];

				if(i==0)				sprintf(idString, "noNoise_");
				else
				{
					char numbers[4];
					if(		i<10)		sprintf(numbers, "000%u", i);
					else if(i<100)		sprintf(numbers, "00%u", i);
					else if(i<1000)		sprintf(numbers, "0%u", i);
					else				sprintf(numbers, "%u", i);

					sprintf(idString, "noise%s_", numbers);

					//photBoundary.addSignedFractionNoise(0.1, i);
					//photBoundary.addPixelNoise(0.9, i);
				}

				computeAndMapKielGrid(filename, idString, 2.5*r_sol, compRadialRes, false, imgThetaRes, imgPhiRes, false, 1.0);
            }
        }
        else
        {
            printf("ERROR! HelioMagfield::paramStudyNoise: file %s is not a regular file!\n", dir_iter->path().c_str());
            continue;
        }
    }

    printf("--------------------------------------------------------------------------\n");
    printf("--- INFO: PFSSsolution::paramStudyNoise concluded. \n");
    printf("--------------------------------------------------------------------------\n\n");
}

/*! conducts a parameter study of the solver threshold level on input file named filename
 *
 */
void PFSSsolution::paramStudyThresh(const char *filename, uint compRadialRes, uint compThetaRes, uint compPhiRes,
		uint imgThetaRes,  uint imgPhiRes)
{
    if (!doesFileExist(filename))
    {
    	printf("ERROR!\tPFSSsoluion::paramStudyThresh\n\t input file '%s' does not exist\n", filename);
        return;
    }

    if (outDir[0] == '\0')
	{
    	printf("ERROR!\tPFSSsolution::paramStudyNoise\n\t output directory has not been set!\n");
    	return;
	}

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::paramStudyThresh on file: '%s'\n", filename);
	printf("--------------------------------------------------------------------\n\n");

	uint numThresholds = 30;
	hcFloat threshes[numThresholds];
	uint times[numThresholds];

	for(uint i=1; i<=numThresholds; ++i)
	{
		char idString[1000];

		hcFloat thresh;
		if(		i<=10)		thresh = i*1E-4;
		else if(i<=20)		thresh = (i-10)*1E-3;
		else				thresh = (i-20)*1E-2;

		sprintf(idString, "thresh%E_", thresh);

		computeAndMapKielGrid(filename, idString, 2.5*r_sol, compRadialRes, false, imgThetaRes, imgPhiRes, false, 1.0);

		save();
		multiMapSolution(imgThetaRes, imgPhiRes, false, true, true);
	}

    printf("--------------------------------------------------------------------------\n");
    printf("--- INFO: PFSSsolution::paramStudyThresh concluded. \n");
    printf("--------------------------------------------------------------------------\n\n");
}

/*! conducts a parameter study of the height of the source surface on files in directory inDir
 *
 */
void PFSSsolution::paramStudyRss(const char *inDir, uint compRadialRes, uint compThetaRes, uint compPhiRes,
								uint imgThetaRes,  uint imgPhiRes)
{
    path directory(inDir);
    directory_iterator end_iter;

    if (!exists(inDir) || !is_directory(inDir))
    {
    	printf("ERROR!\tPFSSsolution::paramStudyRss\n\t input directory '%s' does not exist\n", inDir);
        return;
    }

    if (outDir[0] == '\0')
	{
    	printf("ERROR!\tPFSSsolution::paramStudyRss\n\t output directory has not been set!\n");
    	return;
	}

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::paramStudyRss on dir: '%s'\n", inDir);
	printf("--------------------------------------------------------------------\n\n");

	for (directory_iterator dir_iter(inDir) ; dir_iter != end_iter ; ++dir_iter)
    {
        if (is_regular_file(dir_iter->status()))
        {
        	char filename[1000];
        	sprintf(filename, dir_iter->path().c_str());

            printf("Working on file\n'%s'\n\n", filename);

            for(uint i=4; i<5; ++i)
            {
            	hcFloat rSourceSurface;

            	if(     i==1)	rSourceSurface = 1.50 * r_sol;
            	else if(i==2)	rSourceSurface = 2.00 * r_sol;
				else if(i==3)	rSourceSurface = 2.50 * r_sol;
				else if(i==4)	rSourceSurface = 3.00 * r_sol;
				else if(i==5)	rSourceSurface = 3.25 * r_sol;
				else if(i==6)	rSourceSurface = 5.00 * r_sol;
				else if(i==7)	rSourceSurface = 10.0 * r_sol;

            	computeAndMapKielGrid(filename, "", rSourceSurface, compRadialRes, false, imgThetaRes, imgPhiRes, false, 1.0);
            }
        }
        else
        {
            printf("ERROR! PFSSsolution::paramStudyRss: file %s is not a regular file!\n", dir_iter->path().c_str());
            continue;
        }
    }

    printf("--------------------------------------------------------------------------\n");
    printf("--- INFO: PFSSsolution::paramStudyRss concluded. \n");
    printf("--------------------------------------------------------------------------\n\n");
}

/*! conducts a parameter study of the resolution of the numerical grid
 *
 */
void PFSSsolution::paramStudyRes()
{
    if (outDir[0] == '\0')
	{
    	cerr << __FILE__ << ":" << __LINE__ << " output directory has not been set\n";
    	return;
	}

    stringstream filename;
    filename << "../data/PFSS/input/synop_Ml_0.2066.fits";
    uint imgThetaRes	= 200;
    uint imgPhiRes		= 400;
    uint nlow			= 10;
	uint nhigh			= 45;
	hcFloat rss			= 2.5*r_sol;

	cout << "--------------------------------------------------------------------\n";
	cout << "-- PFSS::paramStudyRes on file:\n'" << filename.str() << "'\n";
	cout << "--------------------------------------------------------------------\n\n";

	for(uint i=nlow; i<=nhigh; ++i)
	{
		cout << "Working on radialRes " << i << "\n";fflush(stdout);
		computeAndMapKielGrid(filename.str().data(), "", rss, i, imgThetaRes, imgPhiRes, false, 2.0);
	}

	/*
	for(uint i=nlow; i<=nhigh; ++i)
	{
		cout << "Working on radialRes " << i << "\n";fflush(stdout);
		computeAndMapKielGrid(filename.str().data(), "", 2*r_sol, i, imgThetaRes, imgPhiRes, false, 1.0);
	}/*///

    cout << "--------------------------------------------------------------------------\n";
    cout << "--- INFO: PFSSsolution::paramStudyRes concluded. \n";
    cout << "--------------------------------------------------------------------------\n\n";
}

void PFSSsolution::paramStudyScaleMethod(const char *filename)
{
	if (!doesFileExist(filename))
	{
		printf("ERROR!\tPFSSsoluion::paramStudyAngularRes\n\t input file '%s' does not exist\n", filename);
		return;
	}

	if (outDir[0] == '\0')
	{
		printf("ERROR!\tPFSSsolution::paramStudyAngularRes\n\t output directory has not been set!\n");
		return;
	}

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::paramStudyAngularRes on file '%s'\n", filename);
	printf("--------------------------------------------------------------------\n\n");

	for(uint k=0; k<8; ++k)
	{
		char id[1000];
		if(k==0)		sprintf(id, "scaleMethod_ownCos_");
		else if(k==1)	sprintf(id, "scaleMethod_ownLin_");
		else if(k==2)	sprintf(id, "scaleMethod_box_");
		else if(k==3)	sprintf(id, "scaleMethod_bilin_");
		else if(k==4)	sprintf(id, "scaleMethod_bicub_");
		else if(k==5)	sprintf(id, "scaleMethod_bspline_");
		else if(k==6)	sprintf(id, "scaleMethod_catmul_");
		else if(k==7)	sprintf(id, "scaleMethod_lanczos_");

		for(uint i=0; i<1; ++i)
		{
			uint compRadialRes		= 40+i;
			uint imgThetaRes		= 100;
			uint imgPhiRes			= 200;

			computeAndMapKielGrid(filename, id, 2.5*r_sol, compRadialRes, false, imgThetaRes, imgPhiRes, false, 1.0);
		}
	}

	printf("--------------------------------------------------------------------------\n");
	printf("--- INFO: PFSSsolution::paramStudyAngularRes concluded. \n");
	printf("--------------------------------------------------------------------------\n\n");
}

void PFSSsolution::paramStudyRadialRes(const char *filename)
{

	if (!doesFileExist(filename))
	{
		printf("ERROR!\tPFSSsoluion::paramStudyRadialRes\n\t input file '%s' does not exist\n", filename);
		return;
	}

	if (outDir[0] == '\0')
	{
		printf("ERROR!\tPFSSsolution::paramStudyRadialRes\n\t output directory has not been set!\n");
		return;
	}

	printf("--------------------------------------------------------------------\n");
	printf("-- PFSS::paramStudyRadialRes on file '%s'\n", filename);
	printf("--------------------------------------------------------------------\n\n");

	for(uint i=0; i<30; ++i)
	{
		uint compRadialRes		= 10 + i;
		uint imgThetaRes		= 100;
		uint imgPhiRes			= 200;

		computeAndMapKielGrid(filename, "", 2.5*r_sol, compRadialRes, false, imgThetaRes, imgPhiRes, false, 1.0);
		computeAndMapKielGrid(filename, "", 2.5*r_sol, compRadialRes, true, imgThetaRes, imgPhiRes, false, 1.0);
	}

	printf("--------------------------------------------------------------------------\n");
	printf("--- INFO: PFSSsolution::paramStudyRadialRes concluded. \n");
	printf("--------------------------------------------------------------------------\n\n");
}
