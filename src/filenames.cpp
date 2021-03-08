#include "src/filenames.h"

#include "engine/hcConstants.h"

#include <iostream>
#include <iomanip>
#include <math.h>

#include "boost/regex.hpp"
#include "boost/filesystem.hpp"

using namespace boost;
using namespace boost::filesystem;

//				  CR        inst  group 	model      			rss        							method  					ell|orderSHC						nr       nt       np
regex pat_pfssSolutionInfo("([0-9]{4})_(.*)_(.*)_((?:PFSS)|(?:CSSS))([0-9]{1}\\.[0-9]{2})_((?:GridEll)|(?:GridSph)|(?:SHC))((?:[0-9]{1}\\.[0-9]{2})|(?:[0-9]{2}))?_([0-9]*)x([0-9]*)x([0-9]*)_", regex::icase);

regex pat_magMap("((?:LinLat)|(?:SinLat))Map([0-9]{1}\\.[0-9]{2})((?:comp)|(?:world))_res([0-9]*)x([0-9]*)", regex::icase);


// ---------------------------------------------------------------------------------------------------------------
// directories
// ---------------------------------------------------------------------------------------------------------------

string getDir_EUVdata(uint crNum)
{
	return crNum <= 2055 ? "../data/photMaps/eit/" : "../data/photMaps/aia/";
}

// ---------------------------------------------------------------------------------------------------------------
// PFSS filenames
// ---------------------------------------------------------------------------------------------------------------

string getFilename_pfssSolutionInfo(const PFSSsolutionInfo &info)
{
    if(info.dailyID >= 100000000)
    {
    	cerr << __FILE__ << ":" << __LINE__ << ": DailyID too big (" << info.dailyID << ")\n";
    	return "";
    }

    stringstream dailyIDstr, methodstr, retval, modelstr, groupstr;
    methodstr.precision(2);
    retval.precision(2);
    if(info.instrument == ORIGIN_SOHOMDI_DAILY || info.instrument == ORIGIN_VER_DAILY)
    	dailyIDstr << setw(8) << setfill('0') << info.dailyID;

    retval << setw(4) << setfill('0') << info.CRnum 	<< "_";
    retval << getStringFromOriginID(info.instrument) 	<< dailyIDstr.str() 	<< "_";
    retval << getStringFromGroupID(info.group)			<< "_";
    retval << getStringFromModelID(info.model) 			<< fixed << info.rss/r_sol 		<< "_";
    retval << getStringFromMethodID(info.method);
    //retval.precision(4);
    if(info.method==METH_ELLIPTICAL)	retval 			<< fixed 					<< info.ell;
    //retval.precision(2);
    if(info.method==METH_SHC)			retval 			<< setw(2) << setfill('0') 	<< info.orderSHC;
    retval << "_";
    retval << info.numR << "x" << info.numTheta << "x" << info.numPhi 	<< "_";

    PFSSsolutionInfo inf_m;
    if(!getParamFromFN_pfssSolutionInfo(retval.str(), inf_m))
    {
    	cerr << __FILE__ << ":" << __LINE__ << ": Regex does not match. Filename:\n";
		cout << retval.str() << "\n";
		exit(1);
    }

    return retval.str();
}

string getFilename_pfssSolutionConfig(const PFSSsolutionInfo &info)
{
	stringstream retval;
	retval << getFilename_pfssSolutionInfo(info) << "config.cfg";
	return retval.str();
}

string getFilename_pfssSolutionBin(const PFSSsolutionInfo &info)
{
	stringstream retval;
	retval << getFilename_pfssSolutionInfo(info) << "grid.bin";
	return retval.str();
}

string getFilename_harmonicCoeff(const PFSSsolutionInfo &info)
{
	stringstream retval;
	retval << getFilename_pfssSolutionInfo(info) << "harm.coeff";
	return retval.str();
}

string getFilename_photMagfield(const PFSSsolutionInfo &info)
{
	stringstream retval;
	retval << getFilename_pfssSolutionInfo(info) << "synPhotMagfield.fits";
	return retval.str();
}

string getFilename_magMapping(	const PFSSsolutionInfo &info, hcFloat height,
								bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi)
{
	stringstream retval;
	retval.precision(2);
	retval << getFilename_pfssSolutionInfo(info);
	retval << (sinLatFormat	? "SinLat":"LinLat") << "Map" << fixed << height/r_sol;
	retval << (compCoords 	? "comp"  :"world" ) << "_";
	retval << "res" << resTheta << "x" << resPhi;

	PFSSsolutionInfo inf_m;
	if(!getParamFromFN_magMapping(retval.str(), inf_m, height, sinLatFormat, compCoords, resTheta, resPhi))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": Regex does not match. Filename:\n";
		cout << retval.str() << "\n";
		exit(1);
	}

	return retval.str();
}

string getFilename_magMappingExt(	const PFSSsolutionInfo &info, hcFloat height,
									bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi, const string &extension)
{
	stringstream retval;
	retval << getFilename_magMapping(info, height, sinLatFormat, compCoords, resTheta, resPhi);
	retval << extension;
	return retval.str();
}

string getFilename_magMappingBin(	const PFSSsolutionInfo &info, hcFloat height,
									bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi)
{
	return getFilename_magMappingExt(info, height, sinLatFormat, compCoords, resTheta, resPhi, ".map");
}

string getFilename_magMappingASCII(	const PFSSsolutionInfo &info, hcFloat height,
									bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi)
{
	return getFilename_magMappingExt(info, height, sinLatFormat, compCoords, resTheta, resPhi, ".asc");
}

string getFilename_magMappingImg(	const PFSSsolutionInfo &info, hcFloat height,
									bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi)
{
	return getFilename_magMappingExt(info, height, sinLatFormat, compCoords, resTheta, resPhi, ".bmp");
}

string getFilename_magMappingFootImg(	const PFSSsolutionInfo &info, hcFloat height,
										bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi)
{
	return getFilename_magMappingExt(info, height, sinLatFormat, compCoords, resTheta, resPhi, "_foot.bmp");
}

string getFilename_magMappingExpansion(	const PFSSsolutionInfo &info, hcFloat height,
										bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi)
{
	return getFilename_magMappingExt(info, height, sinLatFormat, compCoords, resTheta, resPhi, "_expansion.fits");
}

string getFilename_magMappingExpansionBitmap(	const PFSSsolutionInfo &info, hcFloat height,
												bool sinLatFormat, bool compCoords, uint resTheta, uint resPhi)
{
	return getFilename_magMappingExt(info, height, sinLatFormat, compCoords, resTheta, resPhi, "_expansion.bmp");
}

string getFilenameAnalysisXB(	string dir, spacecraftID scid, uint crStart, uint crEnd,
								instrumentID magID, instrumentID xbID, hcFloat ss,
								methodID method, hcFloat ellipticity)
{
	string instrMagNameFull		= getStringFromInstrumentID(magID);
	string instrSpeedNameFull	= getStringFromInstrumentID(xbID);

	boost::regex pattern("(?:STEREO|ACE|ULYSSES)-(.*)", boost::regex::icase);
	match_results<string::const_iterator> what;

	string instrMagName, instrSpeedName;
	if(boost::regex_search(instrMagNameFull, what, pattern))	instrMagName = what[1];
	else
	{
		cout << __FILE__ << ":" << __LINE__ << " InstrMagName could not be found.\n";
		return "";
	}

	if(boost::regex_search(instrSpeedNameFull, what, pattern))	instrSpeedName = what[1];
	else
	{
		cout << __FILE__ << ":" << __LINE__ << " InstrSpeedName could not be found.\n";
		return "";
	}

	stringstream fn;
	fn.precision(2);
	fn << dir << "/XuBorovsky_" << getStringFromSpacecraftID(scid) << "_";
	fn << "CR" << crStart << "-CR" << crEnd << "_";
	fn << instrMagName << "_" << instrSpeedName << "_";
	fn << "rss" << fixed << ss << "_";
	fn << (method==METH_SHC?"SHC":(method==METH_NUMERIC?"Spheric":"Elliptic"));
	if(method==METH_ELLIPTICAL)	fn  << "_" << "ell" << fixed << ellipticity;

	return fn.str();
}

string getFilenameAnalysisSWspeed(	string dir, spacecraftID scid, uint crStart, uint crEnd,
									instrumentID magID, instrumentID xbID, hcFloat ss,
									methodID method, hcFloat ellipticity)
{
	string instrMagNameFull		= getStringFromInstrumentID(magID);
	string instrSpeedNameFull	= getStringFromInstrumentID(xbID);

	boost::regex pattern("(?:STEREO|ACE|ULYSSES)-(.*)", boost::regex::icase);
	match_results<string::const_iterator> what;

	string instrMagName, instrSpeedName;
	if(boost::regex_search(instrMagNameFull, what, pattern))	instrMagName = what[1];
	else
	{
		cout << __FILE__ << ":" << __LINE__ << " InstrMagName could not be found.\n";
		return "";
	}

	if(boost::regex_search(instrSpeedNameFull, what, pattern))	instrSpeedName = what[1];
	else
	{
		cout << __FILE__ << ":" << __LINE__ << " InstrSpeedName could not be found.\n";
		return "";
	}

	stringstream fn;
	fn.precision(2);
	fn << dir << "/SWspeed_" << getStringFromSpacecraftID(scid) << "_";
	fn << "CR" << crStart << "-CR" << crEnd << "_";
	fn << instrMagName << "_" << instrSpeedName << "_";
	fn << "rss" << fixed << ss << "_";
	fn << (method==METH_SHC?"SHC":(method==METH_NUMERIC?"Spheric":"Elliptic"));
	if(method==METH_ELLIPTICAL)	fn  << "_" << "ell" << fixed << ellipticity;

	return fn.str();
}

string getFilename_superconfig(const string &path)
{
	stringstream retval;
	retval << path << "/superconfig.cfg";
	return retval.str();
}

// ------------------------------------------------------------------------------------------------------------------------
// EUV filenames
// ------------------------------------------------------------------------------------------------------------------------

string getFilename_EUVdata(uint crNum, euvID euv)
{
	bool EIT 			= crNum <= 2055;
	string directory 	= getDir_EUVdata(crNum);
	string wlstr		= euv==W171?"171":(euv==W193?"193":(euv==W195?"195":(euv==W284?"284":(euv==W304?"304":"unknown"))));
	stringstream retval;
	retval << directory << "CR" << crNum << "_" << wlstr << (EIT?(crNum<=1975?"A_a":(crNum<=1989?"A_b":"A_c")):"") << ".fits";
	return retval.str();
}

string getFilename_EUVfootpointSummary(const PFSSsolutionInfo &info, hcFloat latThresh)
{
	stringstream retval;
	retval.precision(2);
	retval << getFilename_pfssSolutionInfo(info) << "_fpSummary_thresh" << fixed << latThresh << ".fs";
	return retval.str();
}

string getFilename_EUVprefix(const string &outDir, const string &obs, euvID id, const PFSSsolutionInfo &info, hcFloat latThresh)
{
	stringstream retval;
	retval.precision(6);
	retval << outDir << "/" << info.CRnum << "/" << obs << getStringFromEUVid(id) << "_lat" << fixed << latThresh << "_";
	return retval.str();
}

string getFilename_EUVimg(const string &outDir, const string &obs, euvID id, const PFSSsolutionInfo &info, hcFloat latThresh)
{
	stringstream retval;
	retval.precision(6);
	retval << getFilename_EUVprefix(outDir, obs, id, info, latThresh) << "euvmap.fits";
	return retval.str();
}

string getFilename_EUVfootpoints(const string &outDir, const string &obs, euvID id, const PFSSsolutionInfo &info, hcFloat latThresh)
{
	stringstream retval;
	retval.precision(6);
	retval << getFilename_EUVprefix(outDir, obs, id, info, latThresh) << getFilename_pfssSolutionInfo(info) << "_backmap.fits";
	return retval.str();
}

string getFilename_EUVforwOpen(const string &outDir, const string &obs, euvID id, const PFSSsolutionInfo &info, hcFloat latThresh)
{
	stringstream retval;
	retval.precision(6);
	retval << getFilename_EUVprefix(outDir, obs, id, info, latThresh) << getFilename_pfssSolutionInfo(info) << "_forwmapOpen.fits";
	return retval.str();
}

string getFilename_EUVforwClose(const string &outDir, const string &obs, euvID id, const PFSSsolutionInfo &info, hcFloat latThresh)
{
	stringstream retval;
	retval.precision(6);
	retval << getFilename_EUVprefix(outDir, obs, id, info, latThresh) << getFilename_pfssSolutionInfo(info) << "_forwmapClose.fits";
	return retval.str();
}

string getFilename_AIA(const string &instDir, euvID id, uint crNum)
{
	string wl =  id == W171 ? "171" :
				(id == W195	? "195"	:
				(id == W304 ? "304"	: ""));
	if(wl == "")
	{
		cerr << __FILE__ << ":" << __LINE__ << " wavelength " << getStringFromEUVid(id) << " not measured with AIA.\n";
		return "";
	}

	stringstream retval;
	retval << instDir << "/" << "CR" << crNum << "_" << wl << ".fits";
	return retval.str();
}

// ------------------------------------------------------------------------------------------------------------------------
// get parameters from filenames
// ------------------------------------------------------------------------------------------------------------------------

originID getOriginIDfromPhotFilename(const string &fn)
{
	boost::regex pattern(".*WSO\\.[0-9]{4}\\.F\\.txt$", boost::regex::icase);
	if (boost::regex_match(fn, pattern))	return ORIGIN_WSO;

	//pattern.assign(".*hmi\\.synoptic_ml_720s\\.[0-9]{4}\\.synopml\\.fits$", boost::regex::icase);
	pattern.assign(".*hmi\\.synoptic_ml.*\\.[0-9]{4}.*\\.fits$", boost::regex::icase);
	if (boost::regex_match(fn, pattern))	return ORIGIN_SDOHMI;

	pattern.assign(".*m[0-9]{4}f\\.fits$", boost::regex::icase);
	if (boost::regex_match(fn, pattern))	return ORIGIN_NSOKPVT;

	pattern.assign(".*synop_ml_0\\.[0-9]{4}\\.fits$", boost::regex::icase);
	if (boost::regex_match(fn, pattern))	return ORIGIN_SOHOMDI;

	pattern.assign(".*mrmqs[0-9]{6}t[0-9]{4}c[0-9]{4}_000\\.fits$", boost::regex::icase);
	if (boost::regex_match(fn, pattern))	return ORIGIN_NSOGONG;

	pattern.assign(".*mdisynop\\.[0-9]{6}\\.fits$", boost::regex::icase);
	if (boost::regex_match(fn, pattern))	return ORIGIN_SOHOMDI_DAILY;

	pattern.assign(".*verdaily[0-9]{8}\\.fits$", boost::regex::icase);
	if (boost::regex_match(fn, pattern))	return ORIGIN_VER_DAILY;

	pattern.assign(".*[0-9]{4}_.*_.*\\.fits$", boost::regex::icase);
	if (boost::regex_match(fn, pattern))	return ORIGIN_OWN;

    cerr << __FILE__ << ":" << __LINE__ << ": instrument cannot not be identified from filename\n" << fn << "\n\n";
    return ORIGIN_UNKNOWN;
}

bool getParamFromFN_pfssSolutionInfo(string filename, PFSSsolutionInfo &info)
{
	cmatch what;
	if(regex_search(filename.data(), what, pat_pfssSolutionInfo))
	{
		char *endC;
		hcFloat ell 		= 1.0;
		uint orderSHC		= 0;
		uint cr				= strtol(what[1].str().data(), &endC, 10);
		string obsstr		= what[2];
		groupID group		= getGroupIDfromString(what[3].str());
		modelID model		= getModelIDfromString(what[4].str());
		hcFloat rss			= strtof(what[5].str().data(), &endC);
		methodID method		= getMethodIDfromString(what[6].str());
		if(method == METH_ELLIPTICAL)	ell 		= strtof(what[7].str().data(), &endC);
		if(method == METH_SHC)			orderSHC	= strtol(what[7].str().data(), &endC, 10);
		uint numR			= strtol(what[8].str().data(), &endC, 10);
		uint numT			= strtol(what[9].str().data(), &endC, 10);
		uint numP			= strtol(what[10].str().data(), &endC, 10);

		if(ell<1.0)	ell = 1.0/(round(1.0/ell*100)/100.0);	// TODO: this has to be done because only three digits are being stored in filename

		//cout << "getParamFromFN_pfssSolutionInfo ell: " << ell << "\n";

		regex daily("(.*)([0-9]*)", 				regex::icase);
		cmatch wha;

		originID magnetometer;
		uint dailyID;
		if(regex_search(obsstr.data(), wha, daily))
		{
			magnetometer	= getOriginIDfromString(wha[1].str());
			dailyID 		= strtol(wha[2].str().data(), &endC, 10);
		}
		else
		{
			magnetometer	= getOriginIDfromString(obsstr);
			dailyID			= 0;
		}

		SynopticInfo synInf;
		synInf.init(magnetometer, 0.0, cr, dailyID); //maxSinLat is not part of the filename, so we cannot extract it
		info.init(synInf, sizeof(hcFloat), model, method, group, rss*r_sol, ell, orderSHC, numR, numT, numP);

		return true;
	}

	//cout << "regex no match\n";
	return false;
}

bool getParamFromFN_magMapping(	const string &filename, PFSSsolutionInfo &info,
								hcFloat &height, bool &sinLatFormat, bool &compCoords,
								uint &resTheta, uint &resPhi)
{
	if(!getParamFromFN_pfssSolutionInfo(filename, info))
		return false;

	cmatch what;
	if(regex_search(filename.data(), what, pat_magMap))
	{
		char *endC;
		sinLatFormat 	= what[1].str() == "SinLat" ? true : false;
		height			= strtof(what[2].str().data(), &endC) * r_sol;
		compCoords		= what[3].str() == "comp"	? true : false;
		resTheta		= strtol(what[4].str().data(), &endC, 10);
		resPhi			= strtol(what[5].str().data(), &endC, 10);

		return true;
	}
	return false;
}

bool isFileType_pfssSolutionConfig(const string &filename)
{
	stringstream re;
	re << pat_pfssSolutionInfo << "config\\.cfg$";

	regex pattern(re.str(), regex::icase);
	cmatch what;
	if(regex_search(filename.data(), what, pattern)) return true;
	return false;
}

bool isFileType_pfssSolutionBin(const string &filename)
{
	stringstream re;
	re << pat_pfssSolutionInfo << "grid\\.bin$";

	regex pattern(re.str(), regex::icase);
	cmatch what;
	if(regex_search(filename.data(), what, pattern)) return true;
	return false;
}

bool isFileType_harmonicCoeff(const string &filename)
{
	stringstream re;
	re << pat_pfssSolutionInfo << "harm\\.coeff$";
	regex pattern(re.str(), regex::icase);
	cmatch what;
	if(regex_search(filename.data(), what, pattern)) return true;
	return false;
}

bool isFileType_photMagfield(const string &filename)
{
	stringstream re;
	re << pat_pfssSolutionInfo << "synPhotMagfield\\.fits$";
	regex pattern(re.str(), regex::icase);
	cmatch what;
	if(regex_search(filename.data(), what, pattern)) return true;
	return false;
}

bool isFileType_magMappingBin(const string &filename)
{
	stringstream re;
	re << pat_magMap.str() << "\\.map$";
	regex pattern(re.str(), regex::icase);

	cmatch what;
	if(regex_search(filename.data(), what, pattern)) return true;

	return false;
}

bool isFileType_magMappingAscii(const string &filename)
{
	stringstream re;
	re << pat_magMap.str() << "\\.asc$";
	regex pattern(re.str(), regex::icase);

	cmatch what;
	if(regex_search(filename.data(), what, pattern)) return true;

	return false;
}

bool isFileType_magMappingImg(const string &filename)
{
	stringstream re;
	re << pat_magMap.str() << "\\.bmp$";
	regex pattern(re.str(), regex::icase);

	cmatch what;
	if(regex_search(filename.data(), what, pattern)) return true;

	return false;
}

bool isFileType_magMappingFootImg(const string &filename)
{
	stringstream re;
	re << pat_magMap.str() << "_foot\\.bmp$";
	regex pattern(re.str(), regex::icase);

	cmatch what;
	if(regex_search(filename.data(), what, pattern)) return true;

	return false;
}
