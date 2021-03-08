#ifndef HCIMAGEFITS_H
#define HCIMAGEFITS_H

#include "engine/hcTime.h"
#include "engine/hcImage.h"

#include "src/carRotInfo.h"
#include "src/imageStatistics.h"

#include "fitsio.h"
#include "gsl/gsl_fft_real.h"

typedef unsigned int uint;

void printerror(int status);

/*! \brief handler for 2D FITS image files
 */
class hcImageFITS : public hcImageFloat{
public:

	static pthread_mutex_t mutexFits;
    string filename;
    fitsfile *filePtr;    	/*!< pointer to input FITS file  									*/

    hcImageFITS();							/*!< std constructor								*/
    hcImageFITS(const hcImageFITS &other);	/*!< cpy constructor								*/
    hcImageFITS(uint width, uint height);
    hcImageFITS(const string &filename);    	/*!< creates handler by opening a file 				*/
    virtual ~hcImageFITS();					/*!< destructor										*/

    hcImageFITS &operator=(const hcImageFITS &other);	/*!< assignment operator				*/

    void clear();
    void initNULL();

    virtual bool load(const string &filename);
    	/*!< \brief load FITS image from file													*/

    virtual bool save(const string &filename);
    	/*!< \brief write FITS image to file													*/

    bool dumpAllKeys();
    	/*!< \brief print all stored keys to screen												*/

    bool readKeyString(string keyname, string &retval);
    	/*!< \brief read key keyname and give back as string									*/

    bool readKeyFloat(string keyname, hcFloat &value);
    	/*!< \brief read key keyname and convert to floating point value before returning		*/

    bool writeKeyFloat(const string &keyname, const string &comment, hcFloat value);

    hcFloat getHistogramMaximum(uint numBins);

    ImageStatistics getImageStatistics() const;

    percentileDataStruct getPercentiles() const;

    void normalize();

    void rescale(uint newWidth, uint newHeight);

    hcFloat crosscor(const hcImageFITS &other, uint i, uint j);

    void fft(hcImageFITS &real, hcImageFITS &imag);

    bool fft_inv(hcImageFITS &real, hcImageFITS &imag);

    static void initStaticMembers();
};

bool rescaleImage(hcImageFITS &dataIn, hcImageFITS &dataOut, uint widthOut, uint heightOut,
		FREE_IMAGE_FILTER filter=FILTER_LANCZOS3, bool debug=false);

bool rescaleSinLatGrid(hcImageFITS &imgIn,  hcFloat maxSLin,
		  	  	       hcImageFITS &imgOut, uint wOut, uint hOut, hcFloat maxSLout,
		  	  	       bool linSinLat=true, bool debug=false);
	/*!< \brief simple rescaling of sine(latitude) images																			*/

void getMean(const hcImageFITS &img, hcFloat &mean, hcFloat &std);

void getMean(hcFloat *data, uint numDataPoints, hcFloat &mean, hcFloat &std);

void getPercentiles(hcImageFITS &img, percentiles &perc);

#endif
