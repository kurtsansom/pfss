#ifndef PFSSSOLUTION_H
#define PFSSSOLUTION_H

#include "engine/hcTools.h"

#include "src/carRotInfo.h"
#include "src/synPhotMagfield.h"
#include "src/laplaceSolver.h"
#include "src/magMapping.h"

#include "src/enum.h"
#include "src/pfssSolutionInfo.h"

#ifdef GUI
#include "src/footpointData.h"
#endif

class PFSSsolution{
public:

	string outDir;					/*!< \brief the directory to store solution specific data into				*/
	methodID method;				/*!< \brief SHC, numeric spheric or numeric elliptic						*/

    LaplaceSolver solver;           /*!< \brief the actual solver                       						*/
    PFSSsolutionInfo info;    		/*!< \brief information about solution										*/
    SynPhotMagfield photBoundary;   /*!< \brief boundary condition for solver           						*/

    PFSSsolution();                             /*!< \brief std constructor             						*/
    PFSSsolution(const PFSSsolution &other);    /*!< \brief cpy constructor             						*/
    virtual ~PFSSsolution(){}

    void initNULL();
    void clear();
    void init();

    string getFilenamePhotMagfield();
    	/*!< \brief get filename for photospheric synoptic magfield file										*/

    string getFilenameGrid();
    	/*!< \brief get filename for binary grid data															*/

    string getFilenameConfig();
    	/*!< \brief get filename for configuration file															*/

    PFSSsolution &operator=(const PFSSsolution &other);

    LaplaceSolver &getSolver();

    bool loadPhotBoundary(const string &filename, uint seed=0, uint scaleMethod=7);
        /*!< \brief loads magnetic map of photosphere and initializes this object       						*/

    bool save();
        /*!< \brief stores solution computed by solver into file located at path        						*/

    bool load(const string &cfgFN);
        /*!< \brief loads precomputed solution into this object                         						*/

    bool computeKielSHC(const string &filename, uint order, hcFloat r_ss, uint compResR);
    	/*!< \brief compute solution according to Spherical Harmonic Coefficient approach						*/

    bool computeKielGrid(const string &filename, const string &optionalID, hcFloat r_ss,
    		uint compResR, hcFloat surfShapeA, uint scaleMethod=7);
        /*!< \brief computes finite difference laplace equation via GPU or CPU          						*/

    bool mapHeightLevel(hcFloat height, hcFloat *heights, uint numHeights,
    		hcFloat maxSinLat, uint numTheta, uint numPhi,
    		bool sinLatFormat, bool computationalCoords, bool exportASCII);
    	/*!< \brief routine for computing magnetic height maps in memory										*/

    void getMagneticMappingFilename(const MagMapping &magmap, char *filename);
    	/*!< \brief getting filename for export/import of magnetic mappings										*/

    bool setOutDir(const string &outDir);
		/*!< sets output directory for solution storage															*/

    bool multiMapSolution(	uint numTheta, uint numPhi, bool computeIntermediateHeightLevels,
    						bool sinLatFormat, bool computaionalCoords);
    	/*!< \brief computes magnetic field line maps to several height levels									*/

    void loadAndMapCSSS(const char *boundaryFN, const char *coeffFN,
    		uint compRadialRes, uint imgThetaRes, uint imgPhiRes, bool mapIntermediateHeights);
    	/*!< \brief loads Bala CSSS file and computes mapping (TODO: not working so far)						*/

    void filterPhotMaps(const char *inDir);

    bool computeAndMapKielSHC(const string &filename, uint order, hcFloat r_ss,
    		uint compResR, uint mapResTheta, uint mapResPhi, bool mapIntermediateHeights);
    	/*!< \brief computes spherical harmonic coefficients for given photospheric magfield					*/

    bool computeAndMapKielGrid(const string &filename, const char *optionalID, hcFloat r_ss,
    		uint compResR, uint mapResTheta, uint mapResPhi, bool mapIntermediateHeights,
			hcFloat ellipticity, uint scaleMethod=7);
    	/*!< \brief computes and maps Kiel grid approach														*/

    bool loadAndMapKielGrid(const char *filename, uint mapResTheta, uint mapResPhi, bool mapIntermediateHeights);
		/*!< \brief loads Kiel grid solution and computes mapping												*/

	bool loadAndMapStanfordSHC(const char *photFilename, const char *StanfordCoeffFilename,
			uint compResR, uint imgThetaRes, uint imgPhiRes, bool mapIntermediateHeights);
		/*!< \brief loads Stanford SHC file and computes mapping												*/

// ---------------------------------------------------------------------------------------------------------------
// pfssSolution_batches
// ---------------------------------------------------------------------------------------------------------------

	bool batchKielSHC(const char *inDir, hcFloat r_ss,
			uint compResR, uint mapResTheta, uint mapResPhi, bool mapIntermediateHeights);
		/*!< \brief computes and maps Kiel SHC approach for all files in inDir									*/

	bool batchKielGrid(const string &inDir, hcFloat r_ss, uint compResR,
			uint mapResTheta, uint mapResPhi, bool mapIntermediateHeights, hcFloat surfShapeA);
		/*!< \brief computes and maps Kiel grid approach for all files in inDir									*/

	void paramStudyNoise(const char *inDir,
			uint compRadialRes, uint compThetaRes, uint compPhiRes,
			uint imgThetaRes,   uint imgPhiRes, uint numNoiseComputations);
		/*!< \brief analyses impact of noise on solver															*/

	void paramStudyThresh(const char *inDir,
			uint compRadialRes, uint compThetaRes, uint compPhiRes,
			uint imgThetaRes,   uint imgPhiRes);
		/*!< \brief analyses impact of solver threshold level on results										*/

	void paramStudyRss(const char *inDir,
			uint compRadialRes, uint compThetaRes, uint compPhiRes,
			uint imgThetaRes,   uint imgPhiRes);
		/*!< \brief analyses impact of source surface height on results											*/

	void paramStudyRes();
		/*!< \brief analyzes impact of grid point numbers on results											*/

	void paramStudyScaleMethod(const char *filename);
		/*!< \brief analyses impact on image scaling algorithm on results										*/

	void paramStudyRadialRes(const char *filename);

	void compareAnaInter(const char *inDir);

	void compareSHCorders(const string &inDir, hcFloat r_ss,
			uint compRadialRes, uint imgThetaRes, uint imgPhiRes, bool mapIntermediateHeights=false);

	void compareRdist(const char *inDir,
			uint compRadialRes, uint compThetaRes, uint compPhiRes,
			uint imgThetaRes,   uint imgPhiRes);

	void compare2Stanford(const char *photFilename, const char *StanfordCoeffFilename,
			uint compRadialRes, uint compThetaRes, uint compPhiRes, bool rDistribution,
			hcFloat geometricFactor, uint imgThetaRes, uint imgPhiRes, bool mapIntermediateHeights);
		/*!< \brief computes grid method, shc method and loads Stanford file for comparison						*/

	bool footpointAnalysis(hcFloat latThresh=INFINITY);

	bool batchFootpointAnalysis(	const string &inDir, hcFloat r_ss, uint compResR, methodID method,
									hcFloat ellipticity, uint orderSHC, hcFloat latThresh=INFINITY);
		/*!< \brief computes PFSS according to SHC approach and performs footpoint analysis on EIT maps			*/

	bool batchFootpointAnalysisAllComputed(hcFloat latThresh=INFINITY);
		/*!< \brief (re-)computes footpoint analysis for all solutions stored in outDir							*/

// ---------------------------------------------------------------------------------------------------------------
// GUI functions
// ---------------------------------------------------------------------------------------------------------------
#ifdef GUI
	bool insertFootpointData(const char *fn, FootpointData &data);

    bool EUVanalysis(	euvID id, hcImageFITS &euv, bool EIT, hcFloat latThresh,
						MagMapping &mapBack, MagMapping &mapForw, FootpointData &retval);
    	/*!< \brief analysis of EUV footpoint brightness														*/

	bool batchCoronagraphAnalysis(hcSortedList<hcDate> &dates, const hcObliquity &solarRotAxis,
			uint compResR, uint mapResTheta, uint mapResPhi, Imager &coronagraph, Imager &euv, hcFloat ellipticity);

	bool createWhitelightCorona(Imager &corImg, uint width, uint height, hcFloat clipRadius,
			hcDate date, string fnLOS, string fnImgCor);

	//bool createProjectedView(Imager &corImg, uint width, uint height, uint numMaglX, uint numMaglY,
	//		Imager &photImg, const hcObliquity &solarRotAxis, hcFloat clipRadius, const hcDate &date, bool forwardMap);

	bool createProjectedView(hcImageRGBA &corImg, uint width, uint height, uint numMaglX, uint numMaglY,
			hcImageRGBA &photImg, const hcObliquity &solarRotAxis, hcFloat clipRadius, const hcDate &date, bool forwardMap);
#endif

// ---------------------------------------------------------------------------------------------------------------
// pfssSolution_batches end
// ---------------------------------------------------------------------------------------------------------------

    void dump() const;

private:    

    bool insertSolutionInSuperConfig();
        /*!< \brief adds this solution to the list of all computed solutions            */

};

void rescaleAnalysis(const string &fn);

#endif
