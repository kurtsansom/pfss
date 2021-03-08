#ifndef IMAGESTATISTICS_H
#define IMAGESTATISTICS_H

#include "engine/hcTools.h"

using namespace std;

typedef unsigned int uint;

class ImageStatistics{
public:

	uint numPixels;		/*!< \brief number of non-zero pixels																			*/
	percentiles perc;	/*!< \brief percentiles of content																				*/
	hcFloat mean;		/*!< \brief mean value of image content																			*/
	hcFloat stddev;		/*!< \brief standard deviation of mean of image content															*/

	ImageStatistics();
	ImageStatistics(const ImageStatistics &other);
	ImageStatistics(percentiles perc, uint numPixels, hcFloat mean, hcFloat stddev);

	ImageStatistics &operator=(const ImageStatistics &other);
		/*!< \brief assignment operator																									*/

	void operator+=(const ImageStatistics &other);
		/*!< \brief operator for computing "average" statistics																			*/

	void operator/=(hcFloat factor);
		/*!< \brief operator for computing "average" statistics																			*/

	void initNULL();

	string toString() const;

	void dump() const;
};

#endif
