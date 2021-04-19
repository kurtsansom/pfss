#ifndef HCTOOLS_H
#define HCTOOLS_H

#include "engine/hcConstants.h"
#include "engine/hcTime.h"

#include "math.h"
#include "boost/filesystem.hpp"

#include <sstream>
#include <iomanip>

#ifndef num_eps
#define num_eps 1E-7
#endif

typedef unsigned int uint;

struct percentileDataStruct
{
	hcFloat perc00, perc01, perc02, perc05, perc10;
    hcFloat perc20, perc25, perc30, perc40, perc50;
    hcFloat perc60, perc70, perc75, perc80, perc90;
    hcFloat perc95, perc98, perc99, perc100;

    percentileDataStruct(	const string &perc00, const string &perc01, const string &perc02, const string &perc05,
    						const string &perc10, const string &perc20, const string &perc25, const string &perc30,
    						const string &perc40, const string &perc50, const string &perc60, const string &perc70,
    						const string &perc75, const string &perc80, const string &perc90, const string &perc95,
    						const string &perc98, const string &perc99, const string &perc100)
    {
    	char *endC;
    	this->perc00	= strtof(perc00.data(), &endC);	this->perc01	= strtof(perc01.data(), &endC);	this->perc02	= strtof(perc02.data(), &endC);
    	this->perc05	= strtof(perc05.data(), &endC);	this->perc10	= strtof(perc10.data(), &endC);	this->perc20	= strtof(perc20.data(), &endC);
    	this->perc25	= strtof(perc25.data(), &endC);	this->perc30	= strtof(perc30.data(), &endC);	this->perc40	= strtof(perc40.data(), &endC);
    	this->perc50	= strtof(perc50.data(), &endC);	this->perc60	= strtof(perc60.data(), &endC);	this->perc70	= strtof(perc70.data(), &endC);
    	this->perc75	= strtof(perc75.data(), &endC);	this->perc80	= strtof(perc80.data(), &endC);	this->perc90	= strtof(perc90.data(), &endC);
    	this->perc95	= strtof(perc95.data(), &endC);	this->perc98	= strtof(perc98.data(), &endC);	this->perc99	= strtof(perc99.data(), &endC);
    	this->perc100	= strtof(perc100.data(), &endC);
    }

    percentileDataStruct &operator=(const percentileDataStruct &other)
    {
    	if(this == &other)	return *this;

    	perc00 = other.perc00;	perc01 = other.perc01;	perc02 = other.perc02;	perc05 = other.perc05;	perc10 = other.perc10;
    	perc20 = other.perc20;	perc25 = other.perc25;	perc30 = other.perc30;	perc40 = other.perc40;	perc50 = other.perc50;
    	perc60 = other.perc60;	perc70 = other.perc70;	perc75 = other.perc75;	perc80 = other.perc80; 	perc90 = other.perc90;
    	perc95 = other.perc95; 	perc98 = other.perc98; 	perc99 = other.perc99; 	perc100 = other.perc100;

    	return *this;
    }

    percentileDataStruct()
    {
    	perc00 = 0.0; 	perc01 = 0.0; 	perc02 = 0.0;  	perc05 = 0.0;  	perc10 = 0.0;
    	perc20 = 0.0;  	perc25 = 0.0;  	perc30 = 0.0;  	perc40 = 0.0;  	perc50 = 0.0;
    	perc60 = 0.0;  	perc70 = 0.0;  	perc75 = 0.0;  	perc80 = 0.0;  	perc90 = 0.0;
    	perc95 = 0.0;  	perc98 = 0.0;  	perc99 = 0.0;  	perc100 = 0.0;

    }

    void operator+=(const percentileDataStruct &other)	// TODO: does this even make sense?
	{
    	this->perc00 += other.perc00;
    	this->perc01 += other.perc01;
    	this->perc02 += other.perc02;
    	this->perc05 += other.perc05;
    	this->perc10 += other.perc10;
    	this->perc20 += other.perc20;
    	this->perc25 += other.perc25;
    	this->perc30 += other.perc30;
    	this->perc40 += other.perc40;
    	this->perc50 += other.perc50;
    	this->perc60 += other.perc60;
    	this->perc70 += other.perc70;
    	this->perc75 += other.perc75;
    	this->perc80 += other.perc80;
    	this->perc90 += other.perc90;
    	this->perc95 += other.perc95;
    	this->perc98 += other.perc98;
    	this->perc99 += other.perc99;
    	this->perc100 += other.perc100;
	}

    void operator/=(hcFloat factor)						// TODO: does this even make sense?
	{
    	this->perc00 /= factor;
    	this->perc01 /= factor;
    	this->perc02 /= factor;
    	this->perc05 /= factor;
    	this->perc10 /= factor;
    	this->perc20 /= factor;
    	this->perc25 /= factor;
    	this->perc30 /= factor;
    	this->perc40 /= factor;
    	this->perc50 /= factor;
    	this->perc60 /= factor;
    	this->perc70 /= factor;
    	this->perc75 /= factor;
    	this->perc80 /= factor;
    	this->perc90 /= factor;
    	this->perc95 /= factor;
    	this->perc98 /= factor;
    	this->perc99 /= factor;
    	this->perc100 /= factor;

	}

    string toString() const
    {
    	stringstream retval;
    	retval.precision(4);
    	retval << fixed << setw(8) << setfill(' ') << perc00 << " ";
    	retval << setw(8) << perc01 << " " << setw(8) << perc02 << " " << setw(8) << perc05 << " ";
    	retval << setw(8) << perc10 << " " << setw(8) << perc20 << " " << setw(8) << perc25 << " " << setw(8) << perc30 << " ";
    	retval << setw(8) << perc40 << " " << setw(8) << perc50 << " " << setw(8) << perc60 << " " << setw(8) << perc70 << " ";
    	retval << setw(8) << perc75 << " " << setw(8) << perc80 << " " << setw(8) << perc90 << " " << setw(8) << perc95 << " ";
    	retval << setw(8) << perc98 << " " << setw(8) << perc99 << " " << setw(8) << perc100;
    	return retval.str();
    }

    void dump()
    {
    	printf("Percentile 10: %E\nPercentile 20: %E\nPercentile 25: %E\nPercentile 30: %E\nPercentile 40: %E\nPercentile 50: %E\nPercentile 60: %E\nPercentile 70: %E\nPercentile 75: %E\nPercentile 80: %E\nPercentile 90: %E\nPercentile 100: %E\n",
    			perc10, perc20, perc25, perc30, perc40, perc50, perc60, perc70, perc75, perc80, perc90, perc99);
    }
};
typedef struct percentileDataStruct percentiles;

/*! \brief returns nicely formatted floating point number as string																		*/
string toStr(hcFloat num);

/*! \brief prints status information about program flow to stdout																		*/
void printStdOutMess(const string &file, const int &line, const string &message);

/*! \brief prints error message to stderr																								*/
void printErrMess(const string &file, const int &line, const string &message);

/*! \brief encodes the given color values into one uint to be transfered to the shaders 												*/
uint char2RGBA8(unsigned char red, unsigned char green, unsigned char blue, unsigned char alpha);

/*! \brief prints hexadecimal representation of given integer  																			*/
void printHexCodedInt(uint in);

/*! \brief converts 4-ch. RGBA 8-bit data to 32-bit integer		TODO: retval?															*/
void RGBA82char(uint color, unsigned char &red, unsigned char &green, unsigned char &blue, unsigned char &alpha);

/* \brief make folder structure to desirec location TODO: param string																	*/
bool createFolderStructureTo(const char *filename);

/*! \brief returns the numEntrys tabulated entry in input where entries are separated by one of the chars in divider					*/
string strParser(const string &input, const string &divider, uint numEntry);

/*! \brief returns the number of tabular entries in input string where the entries are divided by one of the chars in divider			*/
uint numTabEntries(const string &input, const string &divider);

/*! \brief coordinate transform (e.g. 0-255 -> -128-127)																				*/
float posTransform(float oldMinValue, float oldMaxValue, float newMinValue, float newMaxValue, float oldPos);

/*! \brief checks if file filename exists                                                                                          		*/
bool doesFileExist(const char* filename);

/*! \brief checks if file filename exists                                                                                          		*/
bool doesFileExist(const string &filename);

/*! \brief checks if directory dirname exists                                                                                      		*/
bool directoryExists(const string &dirname);

bool createDir(const string &dirname);

bool checkFileEx(const string &filename, const string &funcName);

/*! \brief sort floating point values in @param arr with @param num elements															*/
void sort(hcFloat *arr, uint num);

/* \brief bubble-sort on @parm arr with @param num elements																				*/
void sort_bubble(hcFloat *arr, uint num);

/*! \brief calls std::sort on @param arr with @param num elements																		*/
void sort_std(hcFloat *arr, uint num);

/*! \brief maps phi to interval [0,2PI) with cyclic boundary																			*/
hcFloat cyclicPhi(hcFloat phi);

//void computationTimeEstimate(uint numberOfComputations, uint numberComputationsDone, const hcDate &startTime);
/*! \brief determines the differential rotational speed of the sun at the given colatitude                                   			*/
hcFloat solarDiffRotSpeed(hcFloat theta);

/*! \brief computes backmapped position (HAE spheric) of solar wind package measured at posHAE_s  down to radiusSphere					*/
Vec3D parkerBackmapHAE_s(const Vec3D &posHAE_s, const hcDate &date, hcFloat radiusSphere,
		hcFloat solarWindSpeed, bool constRotSpeed=true);

/*! \brief computes backmapped position (HAE cartesian) of solar wind package measured at posHAE_s  down to radiusSphere				*/
Vec3D parkerBackmapHAE_c(const Vec3D &posHAE_c, const hcDate &date, hcFloat radiusSphere,
		hcFloat solarWindSpeed, bool constRotSpeed=true);

/*! \brief computes backmapped position (HGC spheric) of solar wind package measured at posHGC_s down to radiusSphere					*/
Vec3D parkerBackmapHGC_s(const Vec3D &posHGC_s, const hcDate &date, hcFloat radiusSphere,
		hcFloat solarWindSpeed, bool constRotSpeed=true);

#endif
