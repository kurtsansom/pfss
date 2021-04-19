#ifndef HCTIME_H
#define HCTIME_H

#include <stdio.h>
#include <iostream>
#include <string>

#include "boost/format.hpp"
#include "boost/filesystem.hpp"

using namespace std;
using namespace boost;

typedef __int128_t int128;
typedef unsigned int uint;

enum hcCalendarIdentifier {HC_GREGORIAN, HC_JULIAN};
typedef enum hcCalendarIdentifier hcCalendarID;

enum hcTimeStandard {HC_TT, HC_UTC, HC_TAI};
typedef enum hcTimeStandard hcTimeStd;

enum hcEpochIdentifier {HC_J2000_0, HC_B1950_0, HC_B1900_0};
typedef enum hcEpochIdentifier hcEpoch;

void printInt128(int128 value);

class hcDate{
public:

	static const int128 facSec = 1000000000;  /*!< \brief time is stored with nano-second accuracy 									*/

    hcCalendarID calendarID;        /*!< \brief type of calendar used (Gregorian, Julian,...)           							*/
    hcTimeStandard timeStandard;    /*!< \brief time standard used to express the date (TT, UTC,...)    							*/

    int128 absoluteTime;            /*!< \brief Absolute time in seconds from epoch J2000.0, Jan. 01, 2000 : 12:00:00 given in Terrestrial Time */
    int128 equTTtime;               /*!< \brief what the TT clock would show at this date 											*/
    int128 timeSinceJ2000;          /*!< \brief seconds * facSec from epoch in timeStandard 										*/

    long int year;  /*!< \brief christian year  	*/
    uint doy;       /*!< \brief day of year     	*/

    uint month;     /*!< \brief month of year   	*/
    uint dom;       /*!< \brief day of month    	*/
    uint hour;		/*!< \brief hour of day			*/
    uint minute;	/*!< \brief minute of hour		*/
    uint second;	/*!< \brief second of minute	*/
    uint millisec;	/*!< \brief ms of second		*/
    uint microsec;	/*!< \brief us of ms			*/
    uint nanosec;	/*!< \brief ns of us			*/

    hcDate();																			/*!< \brief std constructor					*/
    hcDate(const hcDate &other);														/*!< \brief cpy constructor					*/
    hcDate(long int year, uint month, uint day, uint hour, uint minute, uint second,
            hcTimeStandard timeID = HC_TT, hcCalendarID calendarID = HC_GREGORIAN);
    hcDate(long int year, uint doy, uint hour, uint minute, uint second,
			hcTimeStandard timeID = HC_TT, hcCalendarID calendarID = HC_GREGORIAN);
    hcDate(double year, hcTimeStandard timeID = HC_TT,
    		hcCalendarID calendarID = HC_GREGORIAN);
    hcDate(string timestamp);

    void initNULL();

    hcDate &operator=(const hcDate &other);

    hcDate &operator+=(int128 timeDiff);

    hcDate &operator-=(int128 timeDiff);

    int128 operator-(const hcDate &other) const;

    bool operator>(const hcDate &other) const;
    	/*!< \brief determines if this comes after other																*/

    bool operator>=(const hcDate &other) const;
    	/*!< \brief determines if this comes after other																*/

    bool operator<(const hcDate &other) const;
    	/*!< \brief determines if other comes after this																*/

    bool operator<=(const hcDate &other) const;
    	/*!< \brief determines if other comes after this																*/

    bool operator==(const hcDate &other) const;
    	/*!< \brief determines equality																					*/

    void setFromSystemTime();
    	/*!< \brief reads the time from the system and converts it to hcDate											*/

    bool setFromTimeStamp(const string &timestamp);
        /*!< \brief reads time from string produced by toString		                                                   */

    void set(long int year, uint doy, uint hour, uint minute, uint second,
             hcTimeStandard timeID = HC_TT, hcCalendarID calendarID = HC_GREGORIAN);
        /*!< \brief set to a specific point in time                                                                     */

    bool set(long int year, uint month, uint day, uint hour, uint minute, uint second,
             hcTimeStandard timeID = HC_TT, hcCalendarID calendarID = HC_GREGORIAN);
        /*!< \brief set to a specific point in time (Relative to me. Eat that, Einstein!!) 								*/

    void setFromTT(int128 absoluteTime);
    	/*!< \brief set date from absolute TT time since J2000															*/

    void setFromJD(long int julianDay, double frac);
    	/*!< \brief convert julienDay.frac to hcDate																	*/

    void setFromUnix(int128 unixTime);
    	/*!< \brief set date from UNIX time (seconds since 1970-01-01, 00:00:00 UTC)									*/

    void setFromCarringtonTime(const double &crTime);
    	/*!< \brief set date from Carrington time																		*/

    bool isLeapYear();
    	/*!< \brief tells if this year is a leap year																	*/

    uint monthLength(uint numMonth);
        /*!< \brief returns the length of the given month in days 														*/

    void computeInternalTT();

    void getJulianDate(long &julianDayNum, double &frac) const;

    void getModifiedJulianDate(long &mjd, double &mjd_frac) const;

    int128 getUnixTime() const;
    	/*!< \brief get UNIX time seconds * facsec since Jan. 01, 1970 00:00:00 UTC										*/

    double getCarringtonLongitude() const;
    	/*!< \brief get Carrington Longitude in degrees (sub-earth point, measured from 0° at TODO ???)					*/

    double getCarringtonTime() const;
    	/*!< \brief get Carrington Longitude in degrees (sub-earth point, measured from 0° at TODO ???)					*/

    string getTOD() const;
    	/*!< \brief returns human readable time of day																	*/

    uint getCarringtonRotationNum() const;
		/*!< brief get Carrington rotation number corresponding to this date											*/

    hcDate getCarringtonRotStartDate();
    	/*!< \brief find start time of Carrington rotation number crNum													*/

    hcDate getCarringtonRotEndDate();
    	/*!< \brief find end time of Carrington rotation number crNum													*/

    uint computeDOY();
        /*!< \brief computes day of year from values set in the fields year, month and day                              */

    void computeDayAndMonthFromDOY();
        /*!< \brief computes day and month from the fields doy and year                                                 */

    void setFromInternalTT();

    void convert2(hcTimeStandard std);

    uint numLSsinceBeginning() const;

    bool isLeapSec() const;

    uint getWeekDay() const;

    string toString() const;

    string toSpiceString() const;

    bool exportBinary(std::ofstream &stream);

    bool importBinary(std::ifstream &Sstream);

    void str(char *out);

    void dump() const;

};

hcDate operator+(hcDate lhs, int128 timeDiff);

hcDate operator-(hcDate lhs, int128 timeDiff);

int128 getTTfromJD(long int jd);

long int getJDfromTT(int128 tttime);


bool getJDfromDate(long int year, uint month, uint dom, uint hour, uint minute, uint second, hcCalendarID calendarID,
		long &jd, double &frac);
	/*!< \brief computes Julian Date and fraction from calendar values																	*/

hcDate getDateFromCarringtonTime(const double &crTime);

hcDate getCarringtonRotationStartDate(uint crNum);

hcDate getCarringtonRotationEndDate(uint crNum);


#endif // HCTIME_H
