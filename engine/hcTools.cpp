#include "engine/hcTools.h"
#include "engine/hcTime.h"

#include <stdio.h>

#ifdef GUI
#include "engine/math/hcCoordTransform.h"
#endif

string toStr(hcFloat num)
{
	stringstream ss;
	ss << num;
	return ss.str();
}

string genStreamMessage(const string &file, const int &line, const string &message)
{
	hcDate now;
	string pre = now.getTOD() + "|" + file + ":" + to_string(line);
	stringstream ss;
	ss << left; ss << setfill(' ') << setw(41);
	ss << pre << "| " << message << "\n";
	return ss.str();
}

void printStdOutMess(const string &file, const int &line, const string &message)
{
	string mess = genStreamMessage(file, line, message);
	cout << mess;
}

void printErrMess(const string &file, const int &line, const string &message)
{
	string mess = genStreamMessage(file, line, message);
	cout << mess;
}


uint char2RGBA8(unsigned char red, unsigned char green, unsigned char blue, unsigned char alpha)
{
    return(red  + green * 0x00000100  + blue * 0x00010000 + alpha * 0x01000000);
}

void printHexCodedInt(uint in)
{
    int mask = 0xF;
    printf("%x%x%x%x%x%x%x%x\n", (in>>28)&mask, (in>>24)&mask, (in>>20)&mask, (in>>16)&mask, (in>>12)&mask, (in>>8)&mask, (in>>4)&mask, (in>>0)&mask);
}

void RGBA82char(uint color, unsigned char &red, unsigned char &green, unsigned char &blue, unsigned char &alpha)
{
    red   =  color & 0x000000FF;
    green = (color & 0x0000FF00) >> 8;
    blue  = (color & 0x00FF0000) >> 16;
    alpha = (color & 0xFF000000) >> 24;
}

float posTransform(float oldMinValue, float oldMaxValue, float newMinValue, float newMaxValue, float oldPos)
{
    float scalefactor = (newMaxValue - newMinValue) / (oldMaxValue - oldMinValue);
    return(scalefactor * (oldPos - oldMinValue) + newMinValue);
}

bool doesFileExist(const char* filename)
{
    FILE* fptr = fopen(filename, "r");

    if (fptr != NULL)
    {
        fclose(fptr);
        return true;
    }

    return false;
}

bool doesFileExist(const string &filename)
{
    FILE* fptr = fopen(filename.data(), "r");

    if (fptr != NULL)
    {
        fclose(fptr);
        return true;
    }

    return false;
}

bool directoryExists(const string &dirname)
{
    if(boost::filesystem::exists(dirname)) // && isdirectory
        return true;
    return false;
}

bool createDir(const string &dirname){

    uint length = 0;
    char dirname2[1000] = {0};
    char dir[1000] = {0};

    while(dirname.data()[length] != '\0')
        ++length;

    if(dirname.data()[length-1] != '/' && dirname.data()[length-1] != '\\')
    {
        sprintf(dirname2, "%s/", dirname.data());
        ++length;
    }
    else
        sprintf(dirname2, "%s", dirname.data());

    uint i = 0;
    while(i<length)
    {
        if((dirname2[i] == '/' || dirname2[i] == '\\') && i>0)
        {
            dir[0] = '\0';
            strncat(dir, dirname2, i);
            if(!directoryExists(dir))
            {
                if(!boost::filesystem::create_directory(dir))
                    return false;
            }
        }
        ++i;
    }

    return true;
}

bool checkFileEx(const string &filename, const string &funcName)
{
	//if(!doesFileExist(filename)) cout << "File " << filename << " does not exist!\n";

    if(!doesFileExist(filename))	return false;
    return true;
}

void sort(hcFloat *arr, uint num)
{
	hcDate start, end;
	start.setFromSystemTime();

	sort_std(arr, num);

	end.setFromSystemTime();
}

void sort_bubble(hcFloat *arr, uint num)
{
	if(num==0)
	{
		printf("Sort requested on empty array!\n");
		return;
		//exit(1);
	}

	bool switched = true;

	while(switched)
	{
		switched = false;
		for(uint i=0; i<num-1; ++i)
		{
			if(arr[i] > arr[i+1])
			{
				hcFloat temp 	= arr[i];
				arr[i] 			= arr[i+1];
				arr[i+1]		= temp;
				switched 		= true;
			}
		}
	}
}

void sort_std(hcFloat *arr, uint num)
{
	std::sort(arr, arr+num);
}

bool createFolderStructureTo(const char *filename){

    uint length = 0;
    uint i      = 0;
    uint occ    = 0;

    while(filename[length] != '\0')
        ++length;

    char dir[1000] = {0};

    while(i<length)
    {
        if(filename[i] == '/' || filename[i] == '\\')
            occ = i;
        ++i;
    }

    if(occ == 0)
        return true;

    strncat(dir, filename, occ);
    return createDir(dir);
}

string strParser(const string &input, const string &divider, uint numEntry)
{
	uint numDiv		= divider.length();
	uint len		= input.length();
	uint count 		= 0;
	uint pos		= 0;
	uint posi		= 0;
	string retval	= "";

	if(numDiv == 0)
	{
		cerr << __FILE__ << ":" << __LINE__ << " Divider string is empty.\n";
		return retval;
	}

	if(len == 0)
	{
		cerr << __FILE__ << ":" << __LINE__ << " Input string is empty!\n";
		return retval;
	}

	char div[numDiv];
	for(uint i=0; i<numDiv; ++i)
		div[i] = divider.at(i);

	while(pos < len)
	{
		while(pos < len)
		{
			bool divFound	= false;
			for(uint i=0; i<numDiv; ++i)
				if(input.at(pos) == div[i]) divFound = true;

			if(divFound) break;
			++pos;
		}

		retval	= input.substr(posi, pos-posi);

		while(pos < len)
		{
			bool divFound	= false;
			for(uint i=0; i<numDiv; ++i)
				if(input.at(pos) == div[i]) divFound = true;

			if(!divFound) break;
			++pos;
		}
		posi	= pos;

		if(count == numEntry)
			return retval;

		++count;
	}

	retval = "\n";
	return retval;
}

uint numTabEntries(const string &input, const string &divider)
{
	uint i=0;
	while(true)
	{
		//if(strParser(input, divider, i).length() == 0)
		if(strParser(input, divider, i) == "\n")
			return i;
		++i;
	}
}

hcFloat cyclicPhi(hcFloat phi)
{
	return phi	+ (phi<0 ? (fabs(trunc(phi/(2*PI)))+1)*2*PI : (phi>=2*PI ?  -trunc(phi/(2*PI))*2*PI : 0.0));
}

void computationTimeEstimate(uint numberOfComputations, uint numberComputationsDone, const hcDate &startTime)
{
	hcDate nowTime;
	nowTime.setFromSystemTime();

	uint secondsWorked 	= (nowTime-startTime)/hcDate::facSec;
	hcFloat secsPerSol 	= secondsWorked/(hcFloat)numberComputationsDone;
	uint secondsLeft	= secsPerSol * (numberOfComputations-numberComputationsDone);
	uint secondsTotal	= secondsWorked+secondsLeft;

	uint elapsedDays	= secondsWorked / 86400;
	uint elapsedHours	= (secondsWorked - elapsedDays*86400) / 3600;
	uint elapsedMinutes	= (secondsWorked - elapsedDays*86400 - elapsedHours*3600) / 60;
	uint elapsedSeconds	= (secondsWorked - elapsedDays*86400 - elapsedHours*3600 - elapsedMinutes*60);

	uint leftDays		= secondsLeft / 86400;
	uint leftHours		= (secondsLeft - leftDays*86400) / 3600;
	uint leftMinutes	= (secondsLeft - leftDays*86400 - leftHours*3600) / 60;
	uint leftSeconds	= (secondsLeft - leftDays*86400 - leftHours*3600 - leftMinutes*60);

	uint totalDays		= secondsTotal / 86400;
	uint totalHours		= (secondsTotal - totalDays*86400) / 3600;
	uint totalMinutes	= (secondsTotal - totalDays*86400 - totalHours*3600) / 60;
	uint totalSeconds	= (secondsTotal - totalDays*86400 - totalHours*3600 - totalMinutes*60);

	cout << right << setw(35) << setfill('-') << "\n";
	cout << setw(4) << setfill(' ') << numberComputationsDone << " / " << numberOfComputations << " solutions computed\n";
	cout << "Time elapsed:      " << setw(2) << setfill('0') << elapsedDays 	<< "d" << setw(2) << elapsedHours 	<< "h" << setw(2) << elapsedMinutes  	<< "m" << setw(2) << elapsedSeconds 	<< "s\n";
	cout << "Est. time left:    " << setw(2) <<					leftDays		<< "d" << setw(2) << leftHours		<< "h" << setw(2) << leftMinutes		<< "m" << setw(2) << leftSeconds 		<< "s\n";
	cout << "Est. overall time: " << setw(2) << 				totalDays		<< "d" << setw(2) << totalHours		<< "h" << setw(2) << totalMinutes		<< "m" << setw(2) << totalSeconds		<< "s\n";
	cout << setw(35) << setfill('-') << "\n\n";
	fflush(stdout);
}

/*!	@param 	theta	colatitude at which rotation speed is to be obtained 	(rad)
 * 	@retval			differential rotation speed at colatitude theta			(rad/s)
 *
 */
hcFloat solarDiffRotSpeed(hcFloat theta)
{
	hcFloat dts = 24*60*60;
	hcFloat A   = 14.713/dts;
	hcFloat B   = -2.396/dts;
	hcFloat C  	= -1.787/dts;
	hcFloat lat	= PI / 2 - theta;
    return (A + B * pow(sin(lat), 2) + C * pow(sin(lat), 4)) * deg2rad;
}

#ifdef GUI
/*!	@param posHAE_s			Heliocentriec Aries Ecliptic spherical coordinates					(m / rad / rad)
 * 	@param radiusSphere		sphere to be mapped down to from heliographicPos					(m)
 * 	@param solarWindSpeed	speed of solar wind package to be mapped down						(m/s)
 * 	@param constRotSpeed	use differential rot speed (false) or 1/24.47d						(1)
 *
 * 	@retval	Heliocentric Aries Ecliptic spherical position of backmapped solar wind package		(m / rad / rad)
 */
Vec3D parkerBackmapHAE_s(const Vec3D &posHAE_s, const hcDate &date, hcFloat radiusSphere, hcFloat solarWindSpeed, bool constRotSpeed)
{
	return parkerBackmapHAE_c(posHAE_s.convCoordSpher2Cart(), date, radiusSphere, solarWindSpeed, constRotSpeed).convCoordCart2Spher();
}

/*!	@param posHAE_c			Heliocentriec Aries Ecliptic cartesian coordinates					(m / m / m)
 * 	@param radiusSphere		sphere to be mapped down to from heliographicPos					(m)
 * 	@param solarWindSpeed	speed of solar wind package to be mapped down						(m/s)
 * 	@param constRotSpeed	use differential rot speed (false) or 1/24.47d						(1)
 *
 * 	@retval	Heliocentric Aries Ecliptic cartesian position of backmapped solar wind package		(m / m / m)
 */
Vec3D parkerBackmapHAE_c(const Vec3D &posHAE_c, const hcDate &date, hcFloat radiusSphere, hcFloat solarWindSpeed, bool constRotSpeed)
{
	Vec3D posHGC_s		= HAE2HGC_c(posHAE_c, date).convCoordCart2Spher();
	Vec3D posHGC_BM_c	= parkerBackmapHGC_s(posHGC_s, date, radiusSphere, solarWindSpeed, constRotSpeed).convCoordSpher2Cart();
	Vec3D posHAE_BM_c	= HGC2HAE_c(posHGC_BM_c, date);
	return posHAE_BM_c;
}
#endif

/*!	@param posHGC_s			Heliographic spherical coordinates					(m / rad / rad)
 * 	@param radiusSphere		sphere to be mapped down to from heliographicPos	(m)
 * 	@param solarWindSpeed	speed of solar wind package to be mapped down		(m/s)
 * 	@param constRotSpeed	use differential rot speed (false) or 1/24.47d		(1)
 *
 * 	@retval	Heliographic spherical position of backmapped solar wind package	(m / rad / rad)
 */
Vec3D parkerBackmapHGC_s(const Vec3D &posHGC_s, const hcDate &date, hcFloat radiusSphere, hcFloat solarWindSpeed, bool constRotSpeed)
{
	if(solarWindSpeed < 1)
	{
		cerr << __FILE__ << ":" << __LINE__ << " solar wind speed(" << solarWindSpeed << ") to low for valid computation\n";
		return Vec3D(0.0, 0.0, 0.0);
	}

	if(!posHGC_s.isValid())
	{
		cerr << __FILE__ << ":" << __LINE__ << " you supplied an invalid heliocentric position: " << posHGC_s[0] << "/" << posHGC_s[1] << "/"<< posHGC_s[2] << "\n";
		return Vec3D(0.0, 0.0, 0.0);
	}

	hcFloat omega   		= constRotSpeed ? 2*PI/(24.47*24*60*60) : solarDiffRotSpeed(posHGC_s[1]);
	hcFloat dr      		= (posHGC_s[0] - radiusSphere);
	hcFloat dPhi 			= omega * dr / solarWindSpeed;
	Vec3D posHGC_bm_s		= Vec3D((hcFloat)radiusSphere, posHGC_s[1], (hcFloat)cyclicPhi(posHGC_s[2]+dPhi));
	return posHGC_bm_s;
}
