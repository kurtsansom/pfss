#ifndef FILENAMES_H
#define FILENAMES_H

#include "src/enum.h"
#include "src/ellipticalGrid.h"
#include "src/pfssSolutionInfo.h"

#include <string>
using namespace std;


// ---------------------------------------------------------------------------------------------------------------
// directories
// ---------------------------------------------------------------------------------------------------------------

bool setConfig(const string &filename);

string getFilename_crlist();

// ---------------------------------------------------------------------------------------------------------------
// directories
// ---------------------------------------------------------------------------------------------------------------

string getDir_EUVdata(uint crNum);
	/*!< \brief returns directory of AIA/EIT observations																			*/

// ---------------------------------------------------------------------------------------------------------------
// PFSS filenames
// ---------------------------------------------------------------------------------------------------------------

string getFilename_pfssSolutionInfo(const PFSSsolutionInfo &info);
	/*!< \brief returns part of filename that identifies the computed PFSS-/CSSS-solution											*/

string getFilename_pfssSolutionConfig(const PFSSsolutionInfo &info);

string getFilename_pfssSolutionBin(const PFSSsolutionInfo &info);

string getFilename_harmonicCoeff(const PFSSsolutionInfo &info);

string getFilename_photMagfield(const PFSSsolutionInfo &info);

string getFilename_magMappingBin(	const PFSSsolutionInfo &info, hcFloat height,
									bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi);
	/*!< \brief returns filename for binary magnetic mapping files																	*/

string getFilename_magMappingASCII(	const PFSSsolutionInfo &info, hcFloat height,
									bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi);
	/*!< \brief returns filename for ASCII magnetic mapping files																	*/

string getFilename_magMappingImg(	const PFSSsolutionInfo &info, hcFloat height,
									bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi);
	/*!< \brief returns filename for magnetic mapping polarity images																*/

string getFilename_magMappingFootImg(	const PFSSsolutionInfo &info, hcFloat height,
										bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi);
	/*!< \brief returns filename for magnetic footpoint images																		*/

string getFilename_magMappingExpansion(	const PFSSsolutionInfo &info, hcFloat height,
										bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi);
	/*!< \brief returns filename for magnetic flux tube expansion factor image														*/

string getFilename_magMappingExpansionMax(	const PFSSsolutionInfo &info, hcFloat height,
											bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi);
	/*!< \brief returns filename for magnetic flux tube expansion factor image (max value between phot and ss)						*/

string getFilename_magMappingExpansionDiff(	const PFSSsolutionInfo &info, hcFloat height,
											bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi);
	/*!< \brief returns filename for difference image betwween expansion factor and max expansion factor							*/

string getFilename_magMappingMagfield_r(const PFSSsolutionInfo &info, hcFloat height,
										bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi);
	/*!< \brief returns filename for magnetic field image, radial component															*/

string getFilename_magMappingMagfield_t(const PFSSsolutionInfo &info, hcFloat height,
										bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi);
	/*!< \brief returns filename for magnetic field image, azimuthal component														*/

string getFilename_magMappingMagfield_p(const PFSSsolutionInfo &info, hcFloat height,
										bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi);
	/*!< \brief returns filename for magnetic field image, zonal component															*/

string getFilename_magMappingExpansionBitmap(	const PFSSsolutionInfo &info, hcFloat height,
												bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi);
	/*!< \brief returns filename for magnetic flux tube expansion factor image	(bitmap)											*/

string getFilename_superconfig();
	/*!< \brief get filename for PFSS solutions collection file																		*/

// ---------------------------------------------------------------------------------------------------------------
// in-situ/model comparison analysis filenames
// ---------------------------------------------------------------------------------------------------------------

string getFilenameAnalysisXB(	string dir, spacecraftID scid, uint crStart, uint crEnd,
								instrumentID magID, instrumentID xbID, hcFloat ss,
								methodID method, hcFloat ellipticity);

string getFilenameAnalysisSWspeed(	string dir, spacecraftID scid, uint crStart, uint crEnd,
									instrumentID magID, instrumentID xbID, hcFloat ss,
									methodID method, hcFloat ellipticity);

// ---------------------------------------------------------------------------------------------------------------
// EUV analysis filenames
// ---------------------------------------------------------------------------------------------------------------

string getFilename_EUVdata(uint crNum, euvID euv);
	/*!< \brief returns filename of EIT/AIA synoptic EUV map																		*/

string getFilename_EUVfootpointSummary(const PFSSsolutionInfo &info, hcFloat latThresh);
	/*!< \brief returns filename for EUV analysis summary file																		*/

string getFilename_EUVimg(const string &obs, euvID id, const PFSSsolutionInfo &info, hcFloat latThresh);
	/*!< \brief returns filename for EUV image																						*/

string getFilename_EUVfootpoints(const string &obs, euvID id, const PFSSsolutionInfo &info, hcFloat latThresh);
	/*!< \brief returns filename for EUV footpoint image																			*/

string getFilename_EUVforwOpen(const string &obs, euvID id,	const PFSSsolutionInfo &info, hcFloat latThresh);
	/*!< \brief returns filename for image depicting open photospheric maglines														*/

string getFilename_EUVforwClose(const string &obs, euvID id, const PFSSsolutionInfo &info, hcFloat latThresh);
	/*!< \brief returns filename for image depicting closed photospheric maglines													*/

string getFilename_AIA(const string &instDir, euvID id, uint crNum);
	/*< \brief returns AIA filename for specified Carrington rotation and wavelength (id)											*/

// ---------------------------------------------------------------------------------------------------------------
// functions for extracting parameters from filenames
// ---------------------------------------------------------------------------------------------------------------

originID getOriginIDfromPhotFilename(const string &filename);
	/*!< \brief identify magnetometer from synoptic photospheric magnetogram filename												*/

bool getParamFromFN_pfssSolutionInfo(const	string &filename, PFSSsolutionInfo &info);
	/*!< \brief extracts pfssSolutionInfo parameters from filename																	*/

bool getParamFromFN_magMapping(		const string &filename, PFSSsolutionInfo &info,
									hcFloat &height, bool &sinLatFormat, bool &compCoords,
									uint &resTheta, uint &resPhi);


bool isFileType_pfssSolutionConfig(const string &filename);

bool isFileType_pfssSolutionBin(const string &filename);

bool isFileType_harmonicCoeff(const string &filename);

bool isFileType_photMagfield(const string &filename);

bool isFileType_magMappingBin(const string &filename);

bool isFileType_magMappingAscii(const string &filename);

bool isFileType_magMappingImg(const string &filename);

bool isFileType_magMappingFootImg(const string &filename);


// ---------------------------------------------------------------------------------------------------------------
// file headers
// ---------------------------------------------------------------------------------------------------------------

string getHeader_main();
	/*!< \brief program version, author and license information																		*/

//string getHeader_expansionFactorImage();
	/*!< \brief header information for expansion factor images																		*/


void renameFiles(const string &inDir, const string &outDir);

#endif
