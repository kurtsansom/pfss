#include "engine/hcTime.h"
#include "engine/hcConstants.h"

#include <stdio.h>
#include <math.h>
#include <fstream>

#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/format.hpp"
#include "boost/regex.hpp"

using namespace boost;
using namespace boost::gregorian;
using namespace boost::posix_time;
using namespace std;

const double JD_of_J2000  = 2451545;
const long int leapSeconds[]    = { -867931157, -852033556, -820497555, -788961554, -757425553, -725803152,
									-694267151, -662731150, -631195149, -583934348, -552398347, -520862346,
									-457703945, -378734344, -315575943, -284039942, -236779141, -205243140,
									-173707139, -126273538, -79012737,  -31579136,   189345665,  284040066,
									394372867};

const int128 epochOffset_TAI 	= 32184000000;
const int128 epochOffset_UTC 	= 42184000000;			// what the UTC clock shows at J2000
const int128 epochOffset_UNIX	= 946727957816000000;	// what the UNIX clock shows at J2000


void printInt128(int128 value)
{
    if(value == 0)
    {
       cout << "0";
       return;
    }

    char string[70] 	= {0};
    char tempstr[41] 	= {0};

    bool negative		= value<0 ? true : false;
    if(negative) value 	= -value;

    uint counter 		= 0;
    while (value != 0)
    {
        int lastDigit 		= value % 10;
        tempstr[counter++] 	= "0123456789"[value % 10];   // save last digit
        value 				/= 10;                        // drop it
    }
    tempstr[counter] 	= '\0';
    uint secCount 		= counter + (counter - 1) / 3;
    uint pCount 		= 0;

    string[0]			= negative ? '-' : '+';

    for(uint i=0;i<counter;++i)
    {
    	if(pCount == 3)
    	{
    		string[secCount--] 	= '.';
    		pCount 				= 0;
    	}
    	string[secCount--] = tempstr[i];
    	++pCount;
    }

    cout << string;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                         hcDate
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

hcDate::hcDate()
{
    initNULL();
    setFromSystemTime();
}

hcDate::hcDate(const hcDate &other)
{
	initNULL();
	*this = other;
}

hcDate::hcDate(long year, uint month, uint day, uint hour, uint minute, uint second,
                 hcTimeStandard timeID, hcCalendarID calendarID)
{
	initNULL();
	set(year, month, day, hour, minute, second, timeID, calendarID);
}

hcDate::hcDate(long year, uint doy, uint hour, uint minute, uint second,
                 hcTimeStandard timeID, hcCalendarID calendarID)
{
	initNULL();
	set(year, doy, hour, minute, second, timeID, calendarID);
}

hcDate::hcDate(double year, hcTimeStandard timeID, hcCalendarID calendarID)
{
	initNULL();
	long yr 		= trunc(year);
	hcDate d		= hcDate(yr,0,0,0,0,0);
	uint numDays	= d.isLeapYear() ? 366 : 365;
	uint doy		= trunc((year-yr)*numDays);
	uint s			= trunc((year-yr-(double)doy/numDays)*86400*numDays);
	uint hour		= s/3600;
	uint minute		= (s-3600*hour)/60;
	uint second		= (s-3600*hour-60*minute);

	set(yr, doy, hour, minute, second, timeID, calendarID);
}

hcDate::hcDate(string timestamp)
{
	initNULL();
	setFromTimeStamp(timestamp);
}

void hcDate::initNULL()
{
    calendarID      = HC_GREGORIAN;
    timeStandard    = HC_TT;

    absoluteTime    = 0;
    timeSinceJ2000	= 0;
    equTTtime		= 0;

    year            = 0;
    doy             = 0;

    month           = 0;
    dom             = 0;
    hour            = 0;
    minute          = 0;
    second          = 0;
    millisec        = 0;
    microsec        = 0;
    nanosec         = 0;
}

hcDate &hcDate::operator =(const hcDate &other)
{
	if(this == &other)
		return *this;

	calendarID 		= other.calendarID;
	timeStandard	= other.timeStandard;

	absoluteTime 	= other.absoluteTime;
	equTTtime		= other.equTTtime;
	timeSinceJ2000	= other.timeSinceJ2000;

	year			= other.year;
	doy				= other.doy;

	month			= other.month;
	dom				= other.dom;
	hour			= other.hour;
	minute			= other.minute;
	second			= other.second;
	millisec		= other.millisec;
	microsec		= other.microsec;
	nanosec			= other.nanosec;

	return *this;
}

hcDate &hcDate::operator+=(int128 timeDiff)
{
	absoluteTime += timeDiff;
	convert2(timeStandard);
	return *this;
}

hcDate &hcDate::operator-=(int128 timeDiff)
{

	absoluteTime -= timeDiff;
	convert2(timeStandard);
	return *this;
}

int128 hcDate::operator-(const hcDate &other) const
{
    return this->absoluteTime - other.absoluteTime;
}

hcDate operator+(hcDate lhs, int128 timeDiff)
{
	lhs += timeDiff;
	return lhs;
}

hcDate operator-(hcDate lhs, int128 timeDiff)
{
	lhs -= timeDiff;
	return lhs;
}

bool hcDate::operator>(const hcDate &other) const
{
	return this->absoluteTime > other.absoluteTime ? true : false;
}

bool hcDate::operator>=(const hcDate &other) const
{
	return this->absoluteTime >= other.absoluteTime ? true : false;
}

bool hcDate::operator<(const hcDate &other) const
{
	return this->absoluteTime < other.absoluteTime ? true : false;
}

bool hcDate::operator<=(const hcDate &other) const
{
	return this->absoluteTime <= other.absoluteTime ? true : false;
}

bool hcDate::operator==(const hcDate &other) const
{
	return this->absoluteTime == other.absoluteTime;
}

void hcDate::set(long year, uint doy, uint hour, uint minute, uint second,
                 hcTimeStandard timeID, hcCalendarID calendarID){

    if(calendarID != HC_GREGORIAN)
    {
        printf("ERROR! hcDate::set: Only gregorian calendar supported yet! (not enum %u)\n", calendarID);
        initNULL();
        return;
    }

    if(!(timeID == HC_TT || timeID == HC_UTC || timeID == HC_TAI))
    {
        printf("ERROR! hcDate::set: Time standard (enum %u) not supported!\n", timeID);
        initNULL();
        return;
    }

    this->calendarID    = calendarID;
    this->timeStandard  = timeID;
    this->year          = year;
    this->doy           = doy;
    computeDayAndMonthFromDOY();
    this->hour          = hour;
    this->minute        = minute;
    this->second        = second;
    this->millisec      = 0;
    this->microsec      = 0;
    this->nanosec       = 0;

    if(hour > 23 || minute > 59
                 || (second > 60 && timeID == HC_UTC)
                 || (second > 59 && timeID == HC_TT)
                 || (second > 59 && timeID == HC_TAI))
    {
        printf("ERROR! hcDate::set: Invalid arguments (HH:MM:SS = %u:%u:%u), TimeID: %u!\n", hour, minute, second, timeID);
        initNULL();
        return;
    }

    /*
    if(month > 11)
    {
        printf("ERROR! hcDate::set: Would you kindly tell me, what the %u'th month of the year is called?\n", month+1);
        initNULL();
        return;
    }//*/

    if(isLeapYear() ? doy > 365 : doy >= 365)
    {
        printf("ERROR! hcDate::set: DOY (%u) >= 366!\n", doy);
        dump();
        initNULL();
        return;
    }

    //computeDOY();
    computeInternalTT();
}

bool hcDate::set(long year, uint month, uint day, uint hour, uint minute, uint second,
                 hcTimeStandard timeID, hcCalendarID calendarID){

    if(calendarID != HC_GREGORIAN)
    {
    	cerr << __FILE__ << ":" << __LINE__ << " Only Gregorian calendar supported yet! (not enum " << calendarID << "\n";
        initNULL();
        return false;
    }

    if(!(timeID == HC_TT || timeID == HC_UTC || timeID == HC_TAI))
    {
        cerr << __FILE__ << ":" << __LINE__ << " Time standard (enum " << timeID << ") not supported!\n";
        initNULL();
        return false;
    }

    this->calendarID    = calendarID;
    this->timeStandard  = timeID;
    this->year          = year;
    this->month         = month;
    this->dom           = day;
    this->hour          = hour;
    this->minute        = minute;
    this->second        = second;
    this->millisec      = 0;
    this->microsec      = 0;
    this->nanosec       = 0;

    if(hour > 23 || minute > 59
                 || (second > 60 && timeID == HC_UTC)
                 || (second > 59 && timeID == HC_TT)
                 || (second > 59 && timeID == HC_TAI))
    {
        printf("ERROR! hcDate::set: Invalid arguments (HH:MM:SS = %u:%u:%u), TimeID: %u!\n", hour, minute, second, timeID);
        initNULL();
        return false;
    }

    if(month > 11)
    {
        printf("ERROR! hcDate::set: Would you kindly tell me, what the %u'th month of the year is called?\n", month+1);
        initNULL();
        return false;
    }

    if(day >= monthLength(month))
    {
        printf("ERROR! hcDate::set: Month %u of year %lu does only have %u days. You gave as parameter day = %u!\n",
               month, year, monthLength(month), day);
        initNULL();
        return false;
    }

    computeDOY();
    computeInternalTT();    

    return true;
}

void hcDate::setFromSystemTime()
{
    ptime Jan1st1970(date(1970, 1, 1));
   	ptime Now = microsec_clock::universal_time();
   	time_duration diff = Now - Jan1st1970;
   	setFromUnix(diff.total_seconds() * facSec + diff.fractional_seconds() * 1000);

}

bool hcDate::setFromTimeStamp(const string &timestamp)
{
	string regex1 = "([0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2}):([0-9]{2}):([0-9]{2})";
	string regex2 = "([0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2})([0-9]{2})([0-9]{2})Z";

	boost::regex pattern1(regex1, boost::regex::icase);
	boost::regex pattern2(regex2, boost::regex::icase);

	match_results<string::const_iterator> what;

	if(boost::regex_search(timestamp, what, pattern1))
	{
		//cout << timestamp << " is pattern 1\n";
		char *end;
		uint year 	= strtol(what[1].str().data(), &end, 10);
		uint month 	= strtol(what[2].str().data(), &end, 10);
		uint day  	= strtol(what[3].str().data(), &end, 10);
		uint hour	= strtol(what[4].str().data(), &end, 10);
		uint minute	= strtol(what[5].str().data(), &end, 10);
		uint second	= strtol(what[6].str().data(), &end, 10);

		return set(year, month-1, day-1, hour, minute, second, HC_UTC, HC_GREGORIAN);
	}

	if(boost::regex_search(timestamp, what, pattern2))
	{
		//cout << timestamp << "is pattern 2\n";
		char *end;
		uint year 	= strtol(what[1].str().data(), &end, 10);
		uint month 	= strtol(what[2].str().data(), &end, 10);
		uint day  	= strtol(what[3].str().data(), &end, 10);
		uint hour	= strtol(what[4].str().data(), &end, 10);
		uint minute	= strtol(what[5].str().data(), &end, 10);
		uint second	= strtol(what[6].str().data(), &end, 10);

		return set(year, month-1, day-1, hour, minute, second, HC_UTC, HC_GREGORIAN);
	}

	return false;
}

void hcDate::setFromTT(int128 absoluteTime)
{
	initNULL();
	this->absoluteTime = absoluteTime;
	convert2(HC_TT);
}

void hcDate::setFromUnix(int128 unixTime)
{
	initNULL();
	this->absoluteTime = unixTime - epochOffset_UNIX;

	uint counter = 0;
	uint numLSinList = sizeof(leapSeconds) / sizeof(leapSeconds[0]);

	while(counter < numLSinList && leapSeconds[counter] * facSec - 816000000 <= absoluteTime)
	{
		absoluteTime  += 1 * facSec;
		++counter;
	}

	convert2(HC_TT);
	convert2(HC_UTC);
}

void hcDate::setFromJD(long int julianDayNum, double frac)
{
    this->calendarID    = HC_GREGORIAN;
    this->timeStandard  = HC_TT;
    this->absoluteTime  = getTTfromJD(julianDayNum);
    this->absoluteTime  += (int128)(frac * 86400 * facSec);
    this->equTTtime     = absoluteTime;
    this->timeSinceJ2000= absoluteTime;
    setFromInternalTT();
}

void hcDate::setFromCarringtonTime(const double &crTime)
{
	*this = getDateFromCarringtonTime(crTime);
}

bool hcDate::isLeapYear()
{
    if(calendarID == HC_GREGORIAN)
    {
        if(year%100 == 0)
        {
            if(year%400 == 0)
                return true;
            else
                return false;
        }
        else if(year%4 == 0)
            return true;

        return false;
    }

    cerr << __FILE__ << ":" << __LINE__ << " Calendar ID (" << calendarID <<") not supported (Only HC_GEREGORIAN = "<< HC_GREGORIAN << ")\n";
    return false;
}

uint hcDate::monthLength(uint numMonth)
{
    const int monthLength[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    if(calendarID == HC_GREGORIAN)
    {
        if(isLeapYear() && numMonth == 1)
            return monthLength[numMonth] + 1;
        else
            return monthLength[numMonth];
    }

    printf("ERROR! hcDate::monthLength(): Calendar ID (%u) not supported!\n", calendarID);
    return 0;
}

hcDate hcDate::getCarringtonRotStartDate()
{
	return getCarringtonRotationStartDate(getCarringtonRotationNum());
}

hcDate hcDate::getCarringtonRotEndDate()
{
	return getCarringtonRotationEndDate(getCarringtonRotationNum());
}

uint hcDate::computeDOY(){

    doy = 0;
    for(uint i=0;i<month;++i)
        doy += monthLength(i);

    doy += dom;
    return doy;
}

void hcDate::computeDayAndMonthFromDOY()
{
    dom         = 0;
    month       = 0;
    uint doy    = this->doy;

    for(uint i=0; i<12; ++i)
    {
        if((int)(doy - monthLength(i)) >= 0)
        {
            doy -= monthLength(i);
            ++month;
        }
        else
            break;
    }

    dom = doy;
}

void hcDate::computeInternalTT()
{
	long int a = (14 - (month+1)) / 12;
	long int y = year + 4800 - a;
	long int m = (month+1) + 12 * a - 3;

	long int JD = (dom+1) + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;

	equTTtime 	= getTTfromJD(JD) - 43200 * facSec;
	equTTtime += hour 		* 3600 	* facSec;
	equTTtime += minute 	* 60 	* facSec;
	equTTtime += second     * facSec;
	equTTtime += millisec 	* (facSec / 1000);
	equTTtime += microsec 	* (facSec / 1000000);
	equTTtime += nanosec  	* (facSec / 1000000000);

	timeSinceJ2000	= equTTtime;

    if(timeStandard == HC_TT)
    {        
        absoluteTime 	= timeSinceJ2000;
    }

    if(timeStandard == HC_TAI)
    {
        timeSinceJ2000 	= equTTtime + epochOffset_TAI;
        absoluteTime 	= equTTtime + epochOffset_TAI;
    }

    if(timeStandard == HC_UTC)
    {
        timeSinceJ2000 	= equTTtime + epochOffset_UTC;
        absoluteTime	= equTTtime + epochOffset_UTC;

        uint counter = 0;
        uint numLSinList= sizeof(leapSeconds) / sizeof(leapSeconds[0]);

        while(counter < numLSinList && leapSeconds[counter] * facSec - 816000000 <= absoluteTime)
        {
            absoluteTime  += 1 * facSec;
            ++counter;
        }
        timeSinceJ2000 += 2 * counter * facSec;

        if(second == 60)
		{
			timeSinceJ2000 	-= 2 * facSec;
			absoluteTime	-= 1 * facSec;

			if(!isLeapSec())
			{
				printf("ERROR! hcDate::computeInternalTT: This is not a leap seacond. You provided 60 as parameter for second!\n");
				initNULL();
			}
		}
    }
}

void hcDate::getJulianDate(long int &julianDayNum, double &frac) const
{
	hcDate temp = *this;
	temp.convert2(HC_UTC);
	getJDfromDate(temp.year, temp.month, temp.dom, temp.hour, temp.minute, temp.second, temp.calendarID, julianDayNum, frac);
}

void hcDate::getModifiedJulianDate(long &mjd, double &mjd_frac) const
{
	long jd;
	double jdfrac;
	getJulianDate(jd, jdfrac);
	mjd_frac	= jdfrac - 0.5;
	mjd 		= jd - 2400000;

	if(mjd_frac < 0.0)
	{
		mjd_frac+= 1.0;
		mjd 	-= 1;
	}
}

hcFloat getCarFit(hcFloat x)
{
	// the following procedure has been adopted as explained at http://umtof.umd.edu/pm/crn/CARRTIME.HTML
	hcFloat a1	= 0.26109;
	hcFloat a2	= -3.61907E-5;
	hcFloat a3	= 365.26103;
	hcFloat a4	= 1.91670;
	hcFloat a5	= -0.13059;
	hcFloat a6	= -0.00501286;
	hcFloat a7	= -0.0835779;
	hcFloat a8	= -0.17530;
	hcFloat a9	= -0.00560569;
	hcFloat a10	= -0.0134916;

	hcFloat ret	= a1 + a2*x
				+ a4 * sin(2*PI*x/a3) + a5 * sin(4*PI*x/a3) + a6 * sin(6*PI*x/a3)
				+ a7 * cos(2*PI*x/a3) + a8 * cos(4*PI*x/a3) + a9 * cos(6*PI*x/a3) - (x>364 && x<713 ? a10 : 0.0);

	return ret;
}

double hcDate::getCarringtonLongitude() const
{
	//349.03 degrees at 0000 UT on 1 January 1995
	hcDate ep(1995, 0, 0, 0, 0, 0, HC_UTC, HC_GREGORIAN);
	double x 	= (*this - ep)/facSec/86400.0;					// days since that date
	hcFloat lon	= 349.03 - (360.0 * x / 27.2753) + getCarFit(x);
	int n		= floor(lon/360.0); 							// number of full CR since ep
	hcFloat ret = lon - n*360;
	return ret;
}

double hcDate::getCarringtonTime() const
{
	//Actually the canonical zero meridian used today is the one that passed through the ascending node of the solar equator
	//on the ecliptic at Greenwich mean noon on January 1, 1854 (Julian Day 239 8220.0)
	//hcDate ep(1854, 0, 0, 12, 0, 0, HC_UTC, HC_GREGORIAN);	// this is not the starting time
	//ep.setFromJD(2398167, 0.40193);							// TODO: is this the starting point or 1854-01-01 ?

	double longitude	= getCarringtonLongitude();
	uint n				= getCarringtonRotationNum();
	double retval 		= 360*((int)n-1) + (360.0 - longitude);
	return retval;
}

uint hcDate::getCarringtonRotationNum() const
{
	hcDate ep(1995, 0, 0, 0, 0, 0, HC_UTC, HC_GREGORIAN);
	double x 		= (*this - ep)/facSec/86400.0;					// days since that date
	hcFloat lon		= 349.03 - (360.0 * x / 27.2753) + getCarFit(x);
	int n			= floor(lon/360.0);
	uint retval		= 1891 - n;
	return retval;
}

int128 hcDate::getUnixTime() const
{
	int128 result = absoluteTime + epochOffset_UNIX - numLSsinceBeginning() * facSec;
	if(isLeapSec())
		result += 1 * facSec;
	return result;
}

void hcDate::setFromInternalTT()
{
	bool thisIsLeapSecond = isLeapSec();
	if(thisIsLeapSecond && timeStandard == HC_UTC)
		equTTtime -= 1 * facSec;

	long int JD = getJDfromTT(equTTtime);

	long int a, b, c, d, e, m;
	a = JD+32044;

	if(4 * a + 3 < 0 && (4 * a + 3) % 146097 != 0)
		b = (4 * a + 3) / 146097 - 1;
	else
		b = (4 * a + 3) / 146097;

	if(146097 * b < 0 && (146097 * b) % 4 != 0)
		c = a - (146097 * b) / 4 - 1;
	else
		c = a - (146097 * b) / 4;

	if(4 * c + 3 < 0 && ( 4 * c + 3 ) % 1461 != 0 )
		d = (4 * c + 3 ) / 1461 - 1;
	else
		d = (4 * c + 3 ) / 1461;

	if(1461 * d < 0 && (1461 * d) % 4 != 0)
		e = c - (1461 * d) / 4 - 1;
	else
		e = c - (1461 * d) / 4;

	if(5 * e + 2 < 0 && (5 * e + 2) % 153 != 0)
		m = (5 * e + 2) / 153 - 1;
	else
		m = (5 * e + 2) / 153;


	year 	= 100 * b + d - 4800 + m / 10;
	month 	= m + 3 - 12 * (m / 10) - 1;
	dom		= e - (153 * m + 2) / 5;

	int128 tempVal;

	tempVal = equTTtime - getTTfromJD(JD) + 43200 * facSec;

	if(tempVal >= 86400 * facSec)
	{
		tempVal -= 86400 * facSec;

		if(dom + 1 >= monthLength(month))
		{
			if(month == 11)
			{
				year 	+= 1;
				month 	=  0;
				dom 	=  0;
			}
			else
			{
				month	+= 1;
				dom 	=  0;
			}
		}
		else
		{
			dom	+= 1;
		}
	}

	hour = tempVal / 3600 / facSec;

	tempVal	= tempVal - 3600 * facSec * hour;						// seconds after hour started
	minute	= tempVal / 60 / facSec;

	tempVal = tempVal - 60   * facSec * minute;						// seconds after minute started
	second	= tempVal / facSec;

	tempVal = tempVal -        facSec * second;                     // milliseconds after second started
	millisec = (tempVal * 1000 ) / facSec;

	tempVal = tempVal - (millisec * facSec) / 1000;                  // you get the idea
	microsec = (tempVal * 1000000 ) / facSec;

	tempVal = tempVal - (microsec  * facSec ) / 1000000;
	nanosec = (tempVal * 1000000000) / facSec;

	if(timeStandard == HC_UTC && thisIsLeapSecond)
	{
		second 		+= 1;
		equTTtime 	+= 1 * facSec;
	}

	computeDOY();
}

void hcDate::convert2(hcTimeStandard std)
{
	if(std == HC_TT)
	{
		equTTtime		= absoluteTime;
		timeSinceJ2000	= absoluteTime;
		timeStandard	= std;
		setFromInternalTT();
		return;
	}

	if(std == HC_TAI)
	{
		equTTtime		= absoluteTime - epochOffset_TAI;
		timeSinceJ2000	= absoluteTime;
		timeStandard	= std;
		setFromInternalTT();
		return;
	}

	if(std == HC_UTC)
	{
		equTTtime		= absoluteTime - epochOffset_UTC;
		timeSinceJ2000	= absoluteTime;
		timeStandard	= std;

		uint counter = 0;
		uint numLSinList= sizeof(leapSeconds) / sizeof(leapSeconds[0]);

		while(counter < numLSinList && leapSeconds[counter] * facSec - 816000000 <= absoluteTime)
		{
			equTTtime  		-= 1 * facSec;
			timeSinceJ2000 	+= 1 * facSec;
			++counter;
		}

		if(isLeapSec())
		{
			equTTtime		+= 1 * facSec;
			timeSinceJ2000 	-= 1 * facSec;
		}
		setFromInternalTT();
		return;
	}

    printf("ERROR! hcDate::convert2(): Conversion from time Standard (%u) to (%u) not supported!\n", timeStandard, std);
}

uint hcDate::numLSsinceBeginning() const
{
	uint counter = 0;
	uint numLSinList= sizeof(leapSeconds) / sizeof(leapSeconds[0]);

	while(counter < numLSinList && leapSeconds[counter] * facSec - 816000000 <= absoluteTime)
		++counter;

	return counter;
}

bool hcDate::isLeapSec() const
{
	uint numLSinList= sizeof(leapSeconds) / sizeof(leapSeconds[0]);
	uint counter = 0;

	while(counter < numLSinList && leapSeconds[counter] * facSec - 816000000 <= absoluteTime)
	{
		if(absoluteTime - (leapSeconds[counter] * facSec - 816000000) < 1 * facSec)
		{
			return true;
		}
		++counter;
	}

	return false;
}

/*	0 - monday
 * 	1 - tuesday
 * 	2 - wednesday
 * 	3 - thursday
 * 	4 - friday
 * 	5 - saturday
 * 	6 - sunday
 */
uint hcDate::getWeekDay() const
{
	long jd;
	double frac;
	getJulianDate(jd, frac);
	return jd%7;
}

string hcDate::toString() const
{
	hcDate date = *this;
	date.convert2(HC_UTC);
	return(boost::str(format("%04i-%02i-%02iT%02i%02i%02iZ") % date.year % (date.month+1) % (date.dom+1) % date.hour % date.minute % date.second));
}

string hcDate::toSpiceString() const
{
	hcDate date = *this;
	date.convert2(HC_UTC);
	return(boost::str(format("%04i-%02i-%02iT%02i:%02i:%02i") % date.year % (date.month+1) % (date.dom+1) % date.hour % date.minute % date.second));
}

bool hcDate::exportBinary(std::ofstream &stream)
{
	stream.write(reinterpret_cast<char*>(&absoluteTime),sizeof(int128));
	stream.write(reinterpret_cast<char*>(&timeStandard),sizeof(int));
	stream.write(reinterpret_cast<char*>(&calendarID),sizeof(int));
	return true;
}

bool hcDate::importBinary(std::ifstream &stream)
{
	int128 absT;
	stream.read(reinterpret_cast<char*>(&absT), sizeof(int128));
	setFromTT(absT);

	hcTimeStandard timeStd;
	stream.read(reinterpret_cast<char*>(&timeStd), sizeof(int));
	convert2(timeStd);

	stream.read(reinterpret_cast<char*>(&calendarID), sizeof(int));

	return true;
}

void hcDate::dump() const{

    stringstream calendar, timestd, dayOfWeek;
    uint dow = getWeekDay();

    calendar 	<< (calendarID == HC_GREGORIAN ? "GREGORIAN" : (calendarID == HC_JULIAN ? "JULIAN" : "UNKNOWN"));
    timestd  	<< (timeStandard == HC_TT ? "TT" : (timeStandard == HC_TAI ? "TAI" : (timeStandard == HC_UTC ? "UTC" : "UNKNOWN")));
    dayOfWeek	<< (dow==0?"Monday":(dow==1?"Tuesday":(dow==2?"Wednesday":(dow==3?"Thursday":(dow==4?"Friday":(dow==5?"Saturday":"Sunday"))))));

    long jd;
    double jdFrac;
    getJulianDate(jd, jdFrac);

    cout << "--- Dumping date:---\n\n";
    cout << setw(30) << setfill(' ') << left << "CalendarID:" 		<< calendar.str() << "\n";
    cout << setw(30) << setfill(' ') << left << "Time Standard:" 	<< timestd.str() << "\n";
    cout << setw(30) << setfill(' ') << left << "Day of Week:" 		<< dayOfWeek.str() << "\n";
    cout << setw(30) << setfill(' ') << left << "YYYY-MM-DD_HH:MM:SS:MS:US:NS:";
    cout << setw(4)  << setfill('0') << right << year 	<< "-" << setw(2) << setfill('0') << right << (month + 1) 	<< "-";
    cout << setw(2)  << setfill('0') << right << (dom+1) << "_";
    cout << setw(2)  << setfill('0') << right << hour 	<< ":" << setw(2) << setfill('0') << right << minute      	<< ":";
    cout << setw(2)  << setfill('0') << right << second  << ":" << setw(3) << setfill('0') << right << millisec 		<< ":";
    cout << setw(3)  << setfill('0') << right << microsec<< ":" << setw(3) << setfill('0') << right << nanosec 		<< "\n";
    cout << setw(30) << setfill(' ') << left << "toString:" << toString() << "\n";
    cout << setw(30) << setfill(' ') << left << "Julian Day: " << jd << "\n";
    cout << setw(30) << setfill(' ') << left << "Julian Day Fraction: " << jdFrac << "\n";
    cout << setw(30) << setfill(' ') << left << "UNIX:";printInt128(getUnixTime()); cout << "\n";
    cout << setw(30) << setfill(' ') << left << "AbsoluteTime (TT - s):";printInt128(absoluteTime / facSec);cout << "\n";
    cout << setw(30) << setfill(' ') << left << "AbsoluteTime (TT - ms):";printInt128(absoluteTime * 1000 / facSec);cout << "\n";
    cout << setw(30) << setfill(' ') << left << "AbsoluteTime (TT - ns):";printInt128(absoluteTime);cout << "\n";
    cout << setw(30) << setfill(' ') << left << "Nanoseconds since epoch:";printInt128(timeSinceJ2000);cout << "\n";
    cout << setw(30) << setfill(' ') << left << "Equ TT time (ns)";printInt128(equTTtime);cout << "\n";
    cout << "--------------------------------------------------------------------\n\n";

}

/*	WARNING: Julian Day is calculated from UT, not TT as is done here.
 * 	The difference is tiny, though this needs to be fixed if high accuracy
 * 	computations are necessary. Maybe introduce equUTCtime?
 */
int128 getTTfromJD(long int jd)
{
	return (jd - 2451545) * hcDate::facSec * 86400;
}

/*	WARNING: Julian Day is calculated from UT, not TT as is done here.
 * 	The difference is tiny, though this needs to be fixed if high accuracy
 * 	computations are necessary.  Maybe introduce equUTCtime?
 */
long int getJDfromTT(int128 tttime)
{
	//long retval = 2451545 + tttime / 86400 / hcDate::facSec;
	//if(tttime < 0 && tttime % (86400 * hcDate::facSec) != 0) retval -= 1;

	long double jd 	= 2451545.0 + tttime / 86400.0 / hcDate::facSec;
	long retval		= floor(jd);

	return retval;
}

bool getJDfromDate(long int year, uint month, uint dom, uint hour, uint minute, uint second, hcCalendarID calendarID,
		long &jd, double &frac)
{
	if(calendarID != HC_JULIAN && calendarID != HC_GREGORIAN)
	{
		cerr << __FILE__ << ":" << __LINE__ << " calendarID " << calendarID << " not supported\n";
		return false;
	}

	++month;
	++dom;

	uint y 		= month > 2 ? year 		: year-1;
	uint m		= month > 2 ? month 	: month+12;
	hcFloat d	= dom + hour/24.0 + minute/1440.0 + second/86400.0;
	int b		= calendarID == HC_JULIAN ? 0 : 2 - floor(y/100.0) + floor(y/400.0);

	long double jdFrac	= floor(365.25*(y+4716)) + floor(30.6001*(m+1)) + d + b - 1524.5;
	jd 					= floor(jdFrac);
	frac				= jdFrac - jd;

	return true;
}

hcDate getDateFromCarringtonTime(const double &crTime)
{
	int crNum		= floor(crTime/360.0);
	hcFloat crLong	= 360 - (crTime - crNum*360);

	hcDate ep(1995, 0, 0, 0, 0, 0, HC_UTC, HC_GREGORIAN);
	uint epRotNum	= 1891;
	hcDate date 	= ep;

	int128 dist		= 3600*hcDate::facSec;
	hcFloat x		= 0;

	while((int)(dist/hcDate::facSec) != 0)
	{
		x				= (crLong - (int)(crNum - epRotNum)*360.0 - 349.03 - getCarFit(x))*(-27.2753/360.0);
		hcDate newDate	= ep + x * hcDate::facSec * 86400;
		dist			= date - newDate;
		date			= newDate;
	}

	return date;
}

hcDate getCarringtonRotationStartDate(uint crNum)
{
	return getDateFromCarringtonTime(360.0*crNum + (double)0.01);
}

hcDate getCarringtonRotationEndDate(uint crNum)
{
	return getDateFromCarringtonTime(360.0*(crNum+1) - (double)0.01);
}
