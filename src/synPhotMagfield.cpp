#include "src/synPhotMagfield.h"
#include "src/carRotInfo.h"
#include "src/filenames.h"

#include "boost/regex.hpp"
#include "boost/lexical_cast.hpp"
#include "math.h"
#include <iomanip>

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          SynopticInfo
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

SynopticInfo::SynopticInfo()
{
	initNULL();
}

SynopticInfo::SynopticInfo(const SynopticInfo &other)
{
	initNULL();
	operator=(other);
}

SynopticInfo::SynopticInfo(originID id, hcFloat maxSinLat, uint crNum, uint dailyID)
{
	initNULL();
	init(id, maxSinLat, crNum, dailyID);
}

SynopticInfo::~SynopticInfo()
{
	clear();
}

SynopticInfo &SynopticInfo::operator=(const SynopticInfo &other)
{
	if(this == &other)
		return *this;

	instrument			= other.instrument;
	sinLatFormat		= other.sinLatFormat;
	maxSinLat			= other.maxSinLat;
	dailyID				= other.dailyID;

	((crInfo*)this)->operator =(other);

	return *this;
}

bool SynopticInfo::operator==(const SynopticInfo &other)
{
	bool retval = true;
	retval &= crInfo::operator==(other);
	retval &= instrument 	== other.instrument;
	retval &= sinLatFormat	== other.sinLatFormat;
	retval &= maxSinLat		== other.maxSinLat;
	retval &= dailyID		== other.dailyID;
	return retval;
}

bool SynopticInfo::operator>(const SynopticInfo &other)
{
	if(		crInfo::operator>(other))		return true;
	else if(crInfo::operator<(other))		return false;

	if(		instrument > other.instrument)	return true;
	else if(instrument < other.instrument)	return false;

	if(		maxSinLat > other.maxSinLat)	return true;
	else if(maxSinLat < other.maxSinLat)	return false;

	if(		dailyID > other.dailyID)		return true;
	else if(dailyID < other.dailyID)		return false;

	return false;
}

bool SynopticInfo::operator<(const SynopticInfo &other)
{
	return !(operator>(other) || operator==(other));
}

bool SynopticInfo::exportBinary(ofstream &stream)
{
	bool retval 		= true;
	uint sizeoffloat	= sizeof(hcFloat);
	stream.write(reinterpret_cast<char*>(&sizeoffloat),	sizeof(uint));
	retval &= crInfo::exportBinary(stream);
	stream.write(reinterpret_cast<char*>(&instrument),	sizeof(uint));
	stream.write(reinterpret_cast<char*>(&dailyID),		sizeof(uint));
	stream.write(reinterpret_cast<char*>(&sinLatFormat),sizeof(bool));
	stream.write(reinterpret_cast<char*>(&maxSinLat),	sizeoffloat);
	return retval;
}

bool SynopticInfo::importBinary(ifstream &stream)
{
	bool retval 		= true;
	uint sizeoffloat	= 0;
	stream.read(reinterpret_cast<char*>(&sizeoffloat),	sizeof(uint));
	if(sizeoffloat != 4 && sizeoffloat != 8) return false;
	char *tempFloat		= sizeoffloat == 4 ? reinterpret_cast<char*>(new float()): reinterpret_cast<char*>(new double());
	retval &= crInfo::importBinary(stream);
	stream.read(reinterpret_cast<char*>(&instrument),	sizeof(uint));
	stream.read(reinterpret_cast<char*>(&dailyID),		sizeof(uint));
	stream.read(reinterpret_cast<char*>(&sinLatFormat),	sizeof(bool));
	stream.read(reinterpret_cast<char*>(tempFloat),	sizeoffloat); maxSinLat = sizeoffloat==4 ? *(reinterpret_cast<float*>(tempFloat)) : *(reinterpret_cast<double*>(tempFloat));
	delete tempFloat;
	return retval;
}

void SynopticInfo::initNULL()
{
	dailyID				= 0;
	sinLatFormat		= true;
	instrument			= ORIGIN_UNKNOWN;
	maxSinLat			= 1.0;
	crInfo::initNULL();
}

void SynopticInfo::clear()
{
	initNULL();
}

void SynopticInfo::init(originID id, hcFloat maxSinLat, uint CRnum, uint dailyID)
{
	clear();

	this->dailyID		= dailyID;
	this->sinLatFormat	= true;
	this->instrument 	= id;
	this->maxSinLat		= maxSinLat;
	((crInfo*)this)->init(CRnum);
}

void SynopticInfo::dump(uint indent) const
{
	stringstream ind;
	if(indent > 0) ind << setw(indent) << setfill(' ') << " ";
	cout << ind.str() << "Dumping SynopticInfo:\n";fflush(stdout);
	crInfo::dump(indent+1);
	cout << left;
	cout << ind.str() << setw(20) << setfill(' ') << "Instrument:" 		<< getStringFromOriginID(instrument) 	<< "\n";
	cout << ind.str() << setw(20) << setfill(' ') << "SinLatFormat:" 	<< (sinLatFormat ? "true" : "false")	<< "\n";
	cout << ind.str() << setw(20) << setfill(' ') << "maxSinLat:" 		<< maxSinLat							<< "\n";
	cout << ind.str() << setw(20) << setfill(' ') << "DailyID:" 		<< dailyID								<< "\n";
}


//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          SynPhotMagfield
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

SynPhotMagfield::SynPhotMagfield() : hcImageFITS()
{
	initNULL();
}

SynPhotMagfield::SynPhotMagfield(const SynPhotMagfield &other)
{
	initNULL();
	printf("ERROR! MagCRHandler cpy const not implemented!\n");
	exit(1);
}

SynPhotMagfield &SynPhotMagfield::operator=(const SynPhotMagfield &other)
{
	if (this == &other)
		return *this;

	printf("ERROR! MagCRHandler assignment op not implemented!\n");
	exit(1);

	return *this;
}

void SynPhotMagfield::initNULL()
{
	synInfo.clear();
}

bool SynPhotMagfield::load(const string &filename)
{
	bool retval = false;
	originID id = getOriginIDfromPhotFilename(filename);

	if (id == ORIGIN_WSO)
	{
		char temp[1000] = { 0 };
		char a 			= '0';
		uint i 			= 0;
		uint occ 		= 0;

		while (a != '\0')
		{
			a = filename[i];
			if (a == '.')	occ = i;
			++i;
		}
		strncpy(temp, filename.data() + occ + 1, 4);
		temp[4] = '\0';

		if (!strcmp(temp, "FITS"))		retval = openWSOconverted(	filename);
		else							retval = openWSO(			filename);
	}
	else if(id == ORIGIN_NSOGONG)		retval = openNSOGONG(		filename);
	else if(id == ORIGIN_NSOKPVT)		retval = openKPVT(			filename);
	else if(id == ORIGIN_SOHOMDI)		retval = openSOHOMDI(		filename);
	else if(id == ORIGIN_SOHOMDI_DAILY)	retval = openSOHOMDI_DAILY(	filename);
	else if(id == ORIGIN_VER_DAILY)		retval = openVerDAILY(		filename);
	else if(id == ORIGIN_SDOHMI)		retval = openSDOHMI(		filename);
	else if(id == ORIGIN_OWN)			retval = openOwnFormat(		filename);
	else	cerr << __FILE__ << ":" << __LINE__ << ": File could not be opened. Unknown instrument. Filename: '" << filename << "'\n\n";

	return retval;
}

bool SynPhotMagfield::openOwnFormat(const string &filename)
{
	std::string fn(filename);
	boost::regex re0(".*([0-9]{4})_(OWN|WSO|KPVT|GONG|HMI|MDI|MDI_DAILY|VER_DAILY)_.*\\.fits", boost::regex::icase);

	boost::smatch what;
	uint CRnum;
	std::string instrument;

	if(boost::regex_search(fn, what, re0))
	{
		CRnum 		= boost::lexical_cast<int>(what[1]);
		instrument	= std::string(what[2]);
	}
	else
	{
		cerr << __FILE__ << ":" << __LINE__ << ": Filename could not be matched!\n";
		return false;
	}

	bool retval = true;

	if(!hcImageFITS::load(filename))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": File\n\t'" << filename << "'\n\tcannot not be opened\n";
		retval = false;
	}

	if(!retval) return false;

	originID id;

	if(     !strcmp(instrument.c_str(), "WSO"))			id = ORIGIN_WSO;
	else if(!strcmp(instrument.c_str(), "MDI"))			id = ORIGIN_SOHOMDI;
	else if(!strcmp(instrument.c_str(), "HMI"))			id = ORIGIN_SDOHMI;
	else if(!strcmp(instrument.c_str(), "GONG"))		id = ORIGIN_NSOGONG;
	else if(!strcmp(instrument.c_str(), "KPVT"))		id = ORIGIN_NSOKPVT;
	else if(!strcmp(instrument.c_str(), "MDI_DAILY"))	id = ORIGIN_SOHOMDI_DAILY;
	else if(!strcmp(instrument.c_str(), "VER_DAILY"))	id = ORIGIN_VER_DAILY;
	else if(!strcmp(instrument.c_str(), "OWN"))			id = ORIGIN_OWN;
	else
	{
		cerr << __FILE__ << ":" << __LINE__ << " Instrument name '" << instrument.c_str() << "' cannot not be identified\n";
		return false;
	}

	hcFloat maxSinLat,minSinLat;

	if (!readKeyFloat(string("MAXSINLA"), maxSinLat)|| !readKeyFloat(string("MINSINLA"), minSinLat))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": MAXSINLAT or MINSINLAT key not found in file.\n";
		return false;
	}

	synInfo.init(id, maxSinLat, CRnum, 0);
	return true;
}

bool SynPhotMagfield::openKPVT(const string &filename)
{
	boost::regex pattern(".*m([0-9]{4})f\\.fits", boost::regex::icase);

	if (!hcImageFITS::load(filename) || !boost::regex_match(filename, pattern))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  File '" << filename << "' could not be opened.\n";
		return false;
	}

	hcFloat bscale = 0.0;
	if (!readKeyFloat("BSCALE", bscale))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  KPVT files do have a BSCALE-field.\nThe file you provided (" << filename << ") does not.\n";
		return false;
	}

	hcFloat bzero = 0.0;
	if (!readKeyFloat("BZERO", bzero))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  KPVT files do have a BZERO-field.\nThe file you provided (" << filename << ") does not.\n";
		return false;
	}

	boost::smatch what;
	boost::regex_search(filename, what, pattern);
	uint CRnum = boost::lexical_cast<int>(what[1]);

	// TODO maxSinLat
	synInfo.init(ORIGIN_NSOKPVT, 1.0, CRnum, 0);

	cropPolesFromImage(1);

	float dSinLat = 2*synInfo.maxSinLat / (height - 1);
	for (uint x = 0; x < width; ++x)
		for (uint y = 0; y < height; ++y)
		{
			hcFloat sinLat 		= synInfo.maxSinLat - y * dSinLat;
			data[y * width + x] = (data[y * width + x])	* cos(asin(sinLat)); // TODO: see readme ftp://vso.nso.edu/kpvt/synoptic/README
		}

	return true;
}

bool SynPhotMagfield::openWSO(const string &filename)
{
	char fitsFilename[1000];
	fitsFilename[0] = '\0';
	sprintf(fitsFilename, "%s.FITS", filename.data());
	SynPhotMagfield tempHandler;
	tempHandler.convertWSOtxtToWSOfits(filename, fitsFilename);

	bool retval = openWSOconverted((fitsFilename));
	remove(fitsFilename);
	return retval;
}

bool SynPhotMagfield::openWSOconverted(const string &filename) {

	boost::regex re(".*WSO\\.([0-9]{4})\\.F\\.txt\\.fits", boost::regex::icase);

	if (!hcImageFITS::load(filename) || !boost::regex_match(filename, re))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  File '" << filename << "' cannot be opened.\n";
		return false;
	}

	boost::smatch what;
	boost::regex_search(filename, what, re);
	uint CRnum = boost::lexical_cast<int>(what[1]);

	synInfo.init(ORIGIN_WSO, 14.5/15.0, CRnum, 0);

	for (uint x = 0; x < width; ++x)
		for (uint y = 0; y < height; ++y)
			data[y * width + x] = (data[y * width + x]);

	return true;
}

bool SynPhotMagfield::openSOHOMDI(const string &filename)
{
	boost::regex re(".*\\.([0-9]{4})\\.fits", boost::regex::icase);

	if (!hcImageFITS::load(filename) || !boost::regex_match(filename, re))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  File '" << filename << "' cannot be opened.\n";
		return false;
	}

	hcFloat dSinLat = 0.0;
	if (!readKeyFloat("CDELT2", dSinLat))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  SOHOMDI files do have a CDELT2 field.\nThe file you provided (" << filename << ") does not.\n";
		return false;
	}

	hcFloat bscale 	= 0.0;
	hcFloat bzero	= 0.0;

	if (!readKeyFloat("BSCALE", bscale))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  BSCALE field not found in file '" << filename << "'.\n";
		return false;
	}
	if (!readKeyFloat("BZERO", bzero))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  BZERO field not found in file '" << filename << "'!\n";
		return false;
	}

	for(uint j=0; j<height;++j)
		for(uint k=0; k<width;++k)
			data[j*width + k] = bzero + bscale*data[j*width + k];

	boost::smatch what;
	boost::regex_search(filename, what, re);
	uint CRnum = boost::lexical_cast<int>(what[1]);

	hcFloat maxSinLat = (height - 1) / 2.0 * dSinLat;
	synInfo.init(ORIGIN_SOHOMDI, maxSinLat, CRnum, 0);

	cropPolesFromImage(15);

	return true;
}

bool SynPhotMagfield::openSOHOMDI_DAILY(const string &filename)
{
	boost::regex re(".*mdisynop\\.([0-9]{6})\\.fits", boost::regex::icase);

	if (!hcImageFITS::load(filename) || !boost::regex_match(filename, re))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": openSOHOMDI_DAILY: File '" << filename << "' could not be opened!\n";
		return false;
	}

	hcFloat dSinLat = 0.0;
	if (!readKeyFloat("CDELT2", dSinLat))
	{
		cerr << __FILE__ << ":" << __LINE__ << " Field CDELT2 cannot be found in file\n'" << filename << "'\n";
		return false;
	}

	hcFloat lonLast	= 0.0;
	if (!readKeyFloat("LON_LAST", lonLast))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  SOHOMDI files do have a LON_LAST field.\nThe file you provided (" << filename << ") does not!\n";
		return false;
	}

	hcFloat lonFirst	= 0.0;
	if (!readKeyFloat("LON_FRST", lonFirst))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  SOHOMDI files do have a LON_FRST field.\nThe file you provided (" << filename << ") does not!\n";
		return false;
	}

	boost::smatch what;
	boost::regex_search(filename, what, re);
	uint mdiNum = boost::lexical_cast<int>(what[1]);

	hcFloat cr		= lonLast / 360.0;
	uint crNum		= (cr-floor(cr) > 0.5 ? ceil(cr) : floor(cr));
	hcFloat crTime	= crNum * 360.0;
	int shift		= (int)((crTime-lonLast)*width/360.0);
	shiftXdirection(shift);

	hcFloat maxSinLat = (height - 1) / 2.0 * dSinLat;

	synInfo.init(ORIGIN_SOHOMDI_DAILY, maxSinLat, crNum, mdiNum);

	cropPolesFromImage(10);

	return true;
}

bool SynPhotMagfield::openNSOGONG(const string &filename)
{
	boost::regex re(".*c([0-9]{4})_.*", boost::regex::icase);

	if (!hcImageFITS::load(filename) || !boost::regex_match(filename, re))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  File '" << filename << "' cannot not be opened!\n";
		return false;
	}

	hcFloat dSinLat = 0.0;
	if (!readKeyFloat("CDELT2", dSinLat))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  SOHOMDI files do have a CDELT2 field.\nThe file you provided (" << filename << ") does not!\n";
		return false;
	}

	boost::smatch what;
	boost::regex_search(filename, what, re);
	uint CRnum = boost::lexical_cast<int>(what[1]);

	hcFloat maxSinLat = (height - 1) / 2.0 * dSinLat;

	synInfo.init(ORIGIN_NSOGONG, maxSinLat, CRnum, 0);

	cropPolesFromImage(1);

	return true;
}

bool SynPhotMagfield::openSDOHMI(const string &filename)
{
	boost::regex re(".*hmi\\.synoptic_ml.*\\.([0-9]{4}).*\\.fits$", boost::regex::icase);	// TODO: put regex in filenames.cpp

	if (!hcImageFITS::load(filename) || !boost::regex_match(filename, re))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  File '" << filename << "' cannot not be opened!\n";
		return false;
	}

	boost::smatch what;
	boost::regex_search(filename, what, re);
	uint CRnum 		= boost::lexical_cast<int>(what[1]);
	hcFloat dSinLat = 0.0;

	if (!readKeyFloat("CDELT2", dSinLat))
	{
		cerr << __FILE__ << ":" << __LINE__ << ":  openSDOHMI files do have a CDELT2 field.\nThe file you provided (" << filename << ") does not!\n";
		return false;
	}

	hcFloat maxSinLat = (height - 1) / 2.0 * dSinLat;

	synInfo.init(ORIGIN_SDOHMI, maxSinLat, CRnum, 0);

	cropPolesFromImage(15);

	return true;
}

bool SynPhotMagfield::openVerDAILY(const string &filename)
{
	boost::regex re(".*verdaily([0-9]{8})\\.fits", boost::regex::icase);

	if (!hcImageFITS::load(filename) || !boost::regex_match(filename, re))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": File '" << filename << "' cannot not be opened.\n";
		return false;
	}

	boost::smatch what;
	boost::regex_search(filename, what, re);

	boost::regex re_year(".*verdaily([0-9]{4})[0-9]{4}\\.fits", boost::regex::icase);
	boost::regex_search(filename, what, re_year);
	uint year	= boost::lexical_cast<int>(what[1]);

	boost::regex re_month(".*verdaily[0-9]{4}([0-9]{2})[0-9]{2}\\.fits", boost::regex::icase);
	boost::regex_search(filename, what, re_month);
	uint month	= boost::lexical_cast<int>(what[1]);

	boost::regex re_day(".*verdaily[0-9]{6}([0-9]{2})\\.fits", boost::regex::icase);
	boost::regex_search(filename, what, re_day);
	uint day	= boost::lexical_cast<int>(what[1]);

	hcDate date(year, month-1, day-1, 12, 0, 0, HC_UTC, HC_GREGORIAN);
	uint CRnum = date.getCarringtonRotationNum();

	uint dailyNum	= year * 10000 + month * 100 + day;

	synInfo.init(ORIGIN_VER_DAILY, 1.0, CRnum, dailyNum);

	cropPolesFromImage(50);

	return true;
}

#ifdef GUI
#include "src/grids.h"
bool SynPhotMagfield::createLOSfromGrid(SphericalGrid &grid)
{
	uint numTheta 		= grid.numTheta;
	uint numPhi			= grid.numPhi;

	hcFloat maxSinLat	= grid.maxSinLat;
	hcFloat minSinLat	= grid.minSinLat;

	((hcImageFloat*)this)->init(numPhi, numTheta);

	hcCamera3D cam;
	Vec3D observer(AU, PI/2.0, 0.0);
	observer = observer.convCoordSpher2Cart();
	hcFloat dSinLat	= (maxSinLat - minSinLat) / (numTheta - 1);

	for (uint t=0; t<numTheta; ++t)
	{
		hcFloat theta 	= PI/2.0 - asin(maxSinLat - t * dSinLat);
		hcFloat phi		= 0.0;										// we only look at pixels when they are at center meridian
		Vec3D pos((hcFloat)r_sol, theta, phi);
		pos = pos.convCoordSpher2Cart();

		cam.setPos(pos);
		cam.lookat(observer);
		cam.update();

		Matrix3x3 rotMat;
		rotMat.content[0]	= cam.orientation[0];
		rotMat.content[1]	= cam.orientation[1];
		rotMat.content[2]	= cam.orientation[2];

		rotMat.content[3]	= cam.orientation[4];
		rotMat.content[4]	= cam.orientation[5];
		rotMat.content[5]	= cam.orientation[6];

		rotMat.content[6]	= cam.orientation[8];
		rotMat.content[7]	= cam.orientation[9];
		rotMat.content[8]	= cam.orientation[10];

		for(uint p=0; p<numPhi; ++p)
		{
			uint indG		= grid.getIndex(0, t, p);
			Vec3D gridPos	= grid.getPos(indG, false);
			Vec3D mag		= grid.getB(indG, false);
			//hcFloat sign	= mag[0] > 0.0 ? 1.0 : -1.0;
			mag 			= mag.convVecSpher2Cart(pos);
			Vec3D proj		= rotMat * mag;
			uint ind		= t * width + p;
			data[ind]		= -proj[2];

			//hcFloat Blos	= sin(theta) * mag.length() * sign;
		}
	}
	// TODO

	synInfo.init(ORIGIN_OWN, maxSinLat, 1600, 0);

	return true;
}

bool SynPhotMagfield::loadDipole(hcFloat dipole, uint numTheta, uint numPhi, hcFloat maxSinLat, hcFloat minSinLat)
{
	if(maxSinLat == 0.0 && minSinLat==0.0)
	{
		maxSinLat = getMaxSinLat(numTheta, 1.0);
		minSinLat = -getMaxSinLat(numTheta, 1.0);
	}

	if(maxSinLat > getMaxSinLat(numTheta, 1.0) || minSinLat < -getMaxSinLat(numTheta, 1.0))
	{
		printf("ERROR! SynPhotMagfield::loadDipole: numTheta: %u, maxSinLat: %E, getMaxSinLat(): %E\n", numTheta, maxSinLat, getMaxSinLat(numTheta, 1.0));
		return false;
	}

	((hcImageFloat*)this)->init(numPhi, numTheta);

	hcCamera3D cam;
	Vec3D observer(AU, PI/2.0, 0.0);
	observer = observer.convCoordSpher2Cart();
	hcFloat dSinLat	= (maxSinLat - minSinLat) / (numTheta - 1);

	for (uint t=0; t<numTheta; ++t)
	{
		hcFloat theta 	= PI/2.0 - asin(maxSinLat - t * dSinLat);
		hcFloat phi		= 0.0;										// we only look at pixels when they are at center meridian
		Vec3D pos((hcFloat)r_sol, theta, phi);
		pos = pos.convCoordSpher2Cart();

		cam.setPos(pos);
		cam.lookat(observer);
		cam.update();

		for(uint p=0; p<numPhi; ++p)
		{
			uint ind	= t * width + p;
			data[ind]	= 0.0;
		}
	}
	// TODO

	synInfo.init(ORIGIN_OWN, maxSinLat, 1600, 0);

	return true;
}
#endif

void SynPhotMagfield::removeMonopole()
{
#ifndef REMOVEMONOPOLE
	return;
#endif

	hcFloat numerator 	= 0.0;
	hcFloat delta 		= 0.0;
	hcFloat dSinLat		= 2*synInfo.maxSinLat / (height-1);

	for (uint j=0; j<height; ++j)
	{
		hcFloat theta 		= PI / 2.0 - asin(synInfo.maxSinLat - (height - 1 - j) * dSinLat);
		hcFloat cosecans 	= 1.0 / sin(theta);
		hcFloat tempB 		= 0.0;

		numerator 			+= cosecans;

		for (uint k = 0; k < width; ++k)
			tempB += data[j * width + k];

		delta += tempB * cosecans;
	}

	delta /= (numerator * width);

	for (uint j = 0; j < height; ++j)
		for (uint k = 0; k < width; ++k)
			data[j * width + k] -= delta;
}

/*! Note: maxSinLat and minSinLat only relevant if scaleMethod=0 or scaleMethod=1
 *
 */
bool SynPhotMagfield::remeshImage(uint newWidth, uint newHeight, uint scaleMethod)
{
	if (data == NULL)
	{
		printf("ERROR! SynPhotMagfield::remeshImage: image not initialized!\n");
		return false;
	}

	if(newWidth == 0 || newHeight == 0)
	{
		printf("ERROR! SynPhotMagfield::remeshImage: newWidth (%u) x newHeight (%u) not valid values!\n", newWidth, newHeight);
		return false;
	}

	hcImageFITS newImage(newWidth, newHeight);

	bool retval = false;

	// WARNING: the first two methods should not be used, there are some implementation issues
	if(     scaleMethod == 0)	retval = rescaleSinLatGrid(*this, synInfo.maxSinLat, newImage, newWidth, newHeight, synInfo.maxSinLat, true);
	else if(scaleMethod == 1)	retval = rescaleSinLatGrid(*this, synInfo.maxSinLat, newImage, newWidth, newHeight, synInfo.maxSinLat, false);
	else if(scaleMethod == 2)	retval = rescaleImage(*this, newImage, newWidth, newHeight, FILTER_BOX);		// TODO: maxSinLat not employed
	else if(scaleMethod == 3)	retval = rescaleImage(*this, newImage, newWidth, newHeight, FILTER_BILINEAR);	// TODO: maxSinLat not employed
	else if(scaleMethod == 4)	retval = rescaleImage(*this, newImage, newWidth, newHeight, FILTER_BICUBIC);	// TODO: maxSinLat not employed
	else if(scaleMethod == 5)	retval = rescaleImage(*this, newImage, newWidth, newHeight, FILTER_BSPLINE);	// TODO: maxSinLat not employed
	else if(scaleMethod == 6)	retval = rescaleImage(*this, newImage, newWidth, newHeight, FILTER_CATMULLROM); // TODO: maxSinLat not employed
	else if(scaleMethod == 7)	retval = rescaleImage(*this, newImage, newWidth, newHeight, FILTER_LANCZOS3);	// TODO: maxSinLat not employed

	hcImageFloat::init(newWidth, newHeight);

	for(uint x=0;x<width;++x)
		for(uint y=0;y<height;++y)
			this->operator ()(x,y) = newImage(x,y);

	removeMonopole();

	return retval;
}

void SynPhotMagfield::cropPolesFromImage(uint numPixelsToCrop)
{
	if (2 * numPixelsToCrop >= height)
	{
		cerr << __FILE__ << ":" << __LINE__ << ": Number of pixels to omit (" << numPixelsToCrop << ") bigger than image itself (" << height << ")!\n";
		return;
	}

	hcFloat dSinLat 	= 2*synInfo.maxSinLat / (height - 1);
	hcFloat *imageNew 	= new hcFloat[width * (height - 2*numPixelsToCrop)];

	for (uint x = 0; x < width; ++x)
		for (uint y = numPixelsToCrop; y < height - numPixelsToCrop; ++y)
			imageNew[(y - numPixelsToCrop) * width + x] = data[y * width + x];

	height -= 2 * numPixelsToCrop;

	((hcImageFloat*)this)->init(width, height);
	delete [] data;
	data = imageNew;

	hcFloat maxSinLat 	= (height - 1) * dSinLat / 2.0;

	synInfo.init(synInfo.instrument, maxSinLat, synInfo.CRnum, synInfo.dailyID);
}

bool SynPhotMagfield::convertWSOtxtToWSOfits(const string &infile, const string &outfile)
{
	if (!checkFileEx(infile, "MagCRHandler"))
		return false;

	ifstream in(infile);
	if (!in.is_open())
	{
		cerr << __FILE__ << "/" << __LINE__ << ": File '" <<  infile <<"' exists, but cannot be opended!\n";
		return false;
	}

	in.seekg(0, in.end);
	int length = in.tellg();
	in.seekg(0, in.beg);

	char linein[1000];
	linein[0] = '\0';
	in.getline(linein, 1000);   // we don't need the first three lines
	in.getline(linein, 1000);
	in.getline(linein, 1000);

	int j = -1;
	int i = -1;

	uint ntheta 	= 30;
	uint nphi 		= 72;
	init(nphi, ntheta);
	//hcFloat *image 	= new hcFloat[ntheta * nphi];
	bool lastColumn	= false;

	while (in.tellg() < length)
	{
		uint l = 0;
		uint k = 0;

		linein[0] = '\0';
		in.getline(linein, 1000);
		if (linein[0] == '\0')
			continue;

		char temp[100] = { 0 };
		strncpy(temp, linein, 2);

		char a = linein[l];
		bool breaker = false;

		while (a != '\0')
		{
			if (!strcmp(temp, "CT"))
			{
				++j;
				//*							ignore last column
				if (j == (int)nphi)
				{
					breaker = true;
					break;
				}
				/*/							use last column to smooth first one
				if (j == nphi)
				{
					j = 1;
					lastColumn = true;
				}//*/

				i = -1;
				while (a != ' ' && a != '\t' && a != '\0')
					a = linein[++l];

				if (a == '\0')
				{
					printf("ERROR! MagCRHandler::convertWSOtxtToWSOfits: Phi index: %i, line not valid!\n",	j);
					breaker = true;
					break;
				}
			}

			while (!isdigit(a) && a != '-' && a != '\0')
				a = linein[++l];

			if (a == '\0')
				continue;

			k = l;
			++i;
			a = linein[++l];

			while (a != ' ' && a != '\t' && a != '\0' && a != '-')
				a = linein[++l];

			temp[0] = '\0';
			strncpy(temp, linein + k, l - k);
			if(!lastColumn)
				operator()(nphi - j - 1, ntheta - i - 1) = atof(temp) * 1E-2;
			else
			{
				operator()(nphi - j - 1, ntheta - i - 1) += atof(temp) * 1E-2;
				operator()(nphi - j - 1, ntheta - i - 1) /= 2.0;
			}
		}
		if (breaker)
			break;
	}

	save(outfile);

	return true;
}

bool SynPhotMagfield::save(const string &filename)
{
	if(!hcImageFITS::save(filename))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": could not save file '" << filename << "'\n";
		return false;
	}

	char comment[80];
	float mean 		= 0.0;
	float minVal 	= 100000;
	float maxVal 	= -minVal;

	for (uint x = 0; x < width; ++x)
		for (uint y = 0; y < height; ++y)
		{
			float val = data[y * width + x];
			if (isnan(val))
			{
				printf("NaN-value found while exporting FITS data:\n");
				printf("x: %u, y: %u, width: %u, height: %u\n", x, y, width, height);
				exit(1);
			}
			mean += val;
			if (val < minVal)
				minVal = val;
			if (val > maxVal)
				maxVal = val;
		}
	mean /= (width*height);

	comment[0] = '\0';
	sprintf(comment, "Image mean");
	writeKeyFloat("IMGMN01", comment, mean);

	comment[0] = '\0';
	sprintf(comment, "Image Min");
	writeKeyFloat("IMGMIN01", comment, minVal);

	comment[0] = '\0';
	sprintf(comment, "Image Max");
	writeKeyFloat("IMGMAX01", comment, maxVal);

	comment[0] = '\0';
	sprintf(comment, "Carrington time / rad");
	writeKeyFloat("CDELT1", comment, (2*PI)/width);

	comment[0] = '\0';
	sprintf(comment, "Sine latitude / 1");
	writeKeyFloat("CDELT2", comment, 2*synInfo.maxSinLat/(height-1));

	comment[0] = '\0';
	sprintf(comment, "Sine latitude of center of north-most pixel");
	writeKeyFloat("MAXSINLA", comment, synInfo.maxSinLat);

	comment[0] = '\0';
	sprintf(comment, "Sine latitude of center of south-most pixel)");
	writeKeyFloat("MINSINLA", comment, -synInfo.maxSinLat);

	return true;
}

void SynPhotMagfield::addHomWhiteNoise(float sigma, uint seed)
{
	hcImageFloat::addHomWhiteNoise(sigma, seed);
}

void SynPhotMagfield::addPixelNoise(float fraction, uint seed)
{
	hcImageFloat::addPixelNoise(fraction, seed);
}

void SynPhotMagfield::addSignedFractionNoise(float fraction, uint seed)
{
	hcImageFloat::addSignedFractionNoise(fraction, seed);
}

void SynPhotMagfield::dump() const
{
	cout << "Dumping SynPhotMagfield:\n";
	synInfo.dump();
	cout << "width:  " << width << "\n";
	cout << "height: " << height << "\n";
}
