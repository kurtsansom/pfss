#include "src/imageStatistics.h"

#include <iomanip>

typedef unsigned int uint;

ImageStatistics::ImageStatistics()
{
	initNULL();
}

ImageStatistics::ImageStatistics(const ImageStatistics &other)
{
	initNULL();
	operator=(other);
}

ImageStatistics::ImageStatistics(percentiles perc, uint numPixels, hcFloat mean, hcFloat stddev)
{
	this->perc		= perc;
	this->numPixels	= numPixels;
	this->mean		= mean;
	this->stddev	= stddev;
}

ImageStatistics &ImageStatistics::operator=(const ImageStatistics &other)
{
	if(this == &other)	return *this;

	numPixels	= other.numPixels;
	perc		= other.perc;
	mean		= other.mean;
	stddev		= other.stddev;

	return *this;
}

void ImageStatistics::operator+=(const ImageStatistics &other)
{
	numPixels 	+= other.numPixels;
	mean		+= other.mean;
	stddev		+= other.mean;
	perc		+= other.perc;
}

void ImageStatistics::operator/=(hcFloat factor)
{
	numPixels	/= factor;
	mean		/= factor;
	stddev		/= factor;
	perc		/= factor;
}

void ImageStatistics::initNULL()
{
	numPixels	= 0;
	perc 		= percentiles();
	mean		= 0.0;
	stddev		= 0.0;
}

string ImageStatistics::toString() const
{
	stringstream retval;
	retval.precision(4);
	retval << setw(6) << setfill(' ') << numPixels 	<< " ";	// 0
	retval << fixed << mean 						<< " ";	// 1
	retval << fixed << stddev						<< " "; // 2
	retval << perc.toString();								// 3->21
	return retval.str();
}

void ImageStatistics::dump() const
{
	cout << left;
	cout << "Dumping ImageStatistics:\n";
	//perc.dump();
	cout << setw(20) << setfill(' ') << "# pixels:" << numPixels 	<< "\n";
	cout << setw(20) << setfill(' ') << "mean:" 	<< mean 		<< "\n";
	cout << setw(20) << setfill(' ') << "stddev:" 	<< stddev 		<< "\n";
}
