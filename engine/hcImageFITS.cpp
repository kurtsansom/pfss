#include "engine/hcImageFITS.h"

//#include "fftw3.h"

void printerror(int status)
{
	if (status)
	{
		cout << __FILE__ << ":" << __LINE__ << ": fits_report_error " << status << "\n";
		fits_report_error(stderr, status); /* print error report */
		fflush(stdout);fflush(stderr);
	}
	return;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImageFITS
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

hcImageFITS::hcImageFITS()
{
	initNULL();
}

hcImageFITS::hcImageFITS(const hcImageFITS &other)
{
	initNULL();

	operator=(other);
}

hcImageFITS::hcImageFITS(uint width, uint height)
{
	initNULL();
	init(width, height);
}

hcImageFITS::hcImageFITS(const string &filename)
{
	initNULL();
	load(filename);
}

hcImageFITS::~hcImageFITS()
{
	clear();
}

void hcImageFITS::clear()
{
	int status = 0;

	if (filePtr != NULL)
	{
		pthread_mutex_lock(&hcImageFITS::mutexFits);
		if (fits_close_file(filePtr, &status))
		{
			printf("ERROR! FITSHandler::clear(): infile could not be closed!\n");
			printerror(status);
		}
		pthread_mutex_unlock(&hcImageFITS::mutexFits);
	}

	initNULL();
}

void hcImageFITS::initNULL()
{
	filename[0] = '\0';
	filePtr 	= NULL;
}

hcImageFITS &hcImageFITS::operator =(const hcImageFITS &other)
{
	if(this == &other)
		return *this;

	initNULL();
	hcImageFloat::operator=(other);

	return *this;
}

bool hcImageFITS::load(const string &fn)
{
	if (!checkFileEx(fn, "FITShandler::load"))	return false;

	clear();
	this->filename		= fn;
	bool result			= true;
	long firstPixel 	= 1;
	float nullval		= 1E-20;
	int bugCounter 		= 0;
	long *numAxes 		= new long[2];

	int status			= 0;
	int nfound			= 0;
	int anynull			= 0;

	pthread_mutex_lock(&hcImageFITS::mutexFits);
	while ((bugCounter++ < 10) && fits_open_file(&filePtr, fn.data(), READWRITE, &status));
	stringstream sstream;

	if (bugCounter > 10)
	{
		sstream.str("");
		sstream << __FILE__ << ":" << __LINE__ << ": bug counter reached. Filename:\n\t'" << fn << "'\n";
		cerr << sstream.str();
		printerror(status);
		result = false;
	}

	if (result && fits_read_keys_lng(filePtr, "NAXIS", 1, 2, numAxes, &nfound, &status))
	{
		sstream.str("");
		sstream << __FILE__ << ":" << __LINE__ << ": key 'NAXIS' cannot be read. Filename:\n\t'" << fn << "\n";
		printerror(status);
		result = false;
		sstream << __FILE__ << ":" << __LINE__ << ": try again...\n";
		cerr << sstream.str();
		clear();
		status 		= 0;
		bugCounter 	= 0;
		result 		= true;
		while ((bugCounter++ < 10) && fits_open_file(&filePtr, fn.data(), READWRITE, &status));

		if(bugCounter > 10)
		{
			sstream.str("");
			sstream << __FILE__ << ":" << __LINE__ << ": cannot open file in second try.\n";
			cerr << sstream.str();
			printerror(status);
			result = false;
		}

		if (result && fits_read_keys_lng(filePtr, "NAXIS", 1, 2, numAxes, &nfound, &status))
		{
			sstream.str("");
			sstream << __FILE__ << ":" << __LINE__ << ": second try failed.\n";
			cerr << sstream.str();
			printerror(status);
			result = false;
		}
	}

	if(result)
	{
		((hcImageFloat*)this)->init(numAxes[0], numAxes[1]);

		uint numPixels = width * height;
		float *inArray = new float[numPixels];		// this is necessary because data is stored in single precision,
													// hcImageFloat might be double

		if(fits_read_img(filePtr, TFLOAT, firstPixel, numPixels, &nullval, inArray, &anynull, &status))
		{
			cerr << __FILE__ << ":" << __LINE__ << ": image could not be read.\n";
			printerror(status);
			result = false;
		}

		if(result)	for(uint i=0; i<numPixels; ++i)	data[i] = inArray[i];
		delete [] inArray;
	}
	pthread_mutex_unlock(&hcImageFITS::mutexFits);

	delete [] numAxes;

	if(!result)
	{
		fflush(stdout);
		exit(1);
	}
	return result;
}

bool hcImageFITS::dumpAllKeys()
{
	if(filePtr == NULL)
	{
		cerr << __FILE__ << ":" << __LINE__ << " File not loaded\n";
		return false;
	}
	char card[FLEN_CARD];
	int status = 0,  nkeys, ii;  /* MUST initialize status */

	pthread_mutex_lock(&hcImageFITS::mutexFits);
	fits_get_hdrspace(filePtr, &nkeys, NULL, &status);

	for (ii = 1; ii <= nkeys; ii++)
	{
		fits_read_record(filePtr, ii, card, &status); /* read keyword */
		cout << card << "\n";
	}

	//fits_close_file(fptr, &status);

	if (status)	fits_report_error(stderr, status);
	pthread_mutex_unlock(&hcImageFITS::mutexFits);
	return(status==0);
}

bool hcImageFITS::readKeyString(const string &keyname, string &retval)
{
	retval = "";

	if (filePtr == NULL)
	{
		cerr << "hcImageFITS::readKey: FITShandler has not been initialized with a file!\n";
		return false;
	}

	int status 		= 0;
	bool result 	= true;
	char *comment 	= new char[100];
	char *val		= new char[1000];

	pthread_mutex_lock(&hcImageFITS::mutexFits);
	if(fits_read_key(filePtr, TSTRING, keyname.data(), val, comment, &status))
	{
		//cerr << "ERROR! hcImageFITS::readKeyFloat: Keyword '"<<keyname<<"' not found!\n";
		//printerror(status);
		result = false;
	}
	pthread_mutex_unlock(&hcImageFITS::mutexFits);

	retval = val;
	delete [] comment;
	delete [] val;
	return result;
}

bool hcImageFITS::readKeyFloat(const string &keyname, hcFloat &value)
{
	if (filePtr == NULL)
	{
		cerr << "hcImageFITS::readKey: FITShandler has not been initialized with a file!\n";
		return false;
	}

	int status 		= 0;
	bool result 	= true;
	char *comment 	= new char[100];
	float val		= -1.0;

	pthread_mutex_lock(&hcImageFITS::mutexFits);
	if (fits_read_key(filePtr, TFLOAT, keyname.data(), &val, comment, &status))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": Keyword '"<< keyname <<"' not found!\n";
		printerror(status);
		result = false;
	}
	pthread_mutex_unlock(&hcImageFITS::mutexFits);

	/*
	char card[1000];
	fits_read_card(filePtr, keyname.data(), card, &status);
	cout << __FILE__ << "/" << __LINE__ << ": card: '" << card << "'\n";//*/

	value = val;

	delete [] comment;

	return result;
}

bool hcImageFITS::writeKeyFloat(const string &keyname, const string &comment, hcFloat value)
{
	int status 	= 0;
	bool result = true;
	float val 	= value;

	if(filePtr == NULL)
	{
		cerr << __FILE__ << ":" << __LINE__ << ": filePtr not initialized. Save file first.\n";
		return false;
	}

	stringstream payload;
	payload << std::setw(7) << std::left << keyname << " = " << value << "\0";
	if(comment.length() > 0) payload << " / " << comment << "\0";

	char card[80];
	snprintf(card, 80, "%s", payload.str().data());

	pthread_mutex_lock(&hcImageFITS::mutexFits);
	if (fits_update_card(filePtr, keyname.data(), card, &status))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": key " << keyname << " cannot be written.\n";
		printerror(status);
		result = false;
	}
	pthread_mutex_unlock(&hcImageFITS::mutexFits);

	return result;
}

bool hcImageFITS::writeKeyString(const string &keyname, const string &comment, const string &value)
{
	int status 	= 0;
	string val 	= value;

	if(filePtr == NULL)
	{
		cerr << __FILE__ << ":" << __LINE__ << ": outefilePtr not initialized. Save file first\n";
		return false;
	}

	pthread_mutex_lock(&hcImageFITS::mutexFits);
	if (fits_update_key(filePtr, TSTRING, keyname.data(), const_cast<void*>((void*)val.c_str()), comment.data(), &status))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": key " << keyname << "could not be written.\n";
		printerror(status);
		return false;
	}
	pthread_mutex_unlock(&hcImageFITS::mutexFits);

	return true;
}

bool hcImageFITS::writeKeyComment(const string &comment)
{
	int status 	= 0;
	bool result = true;

	if(filePtr == NULL)
	{
		cerr << __FILE__ << ":" << __LINE__ << ": filePtr not initialized. Save file first.\n";
		return false;
	}

	pthread_mutex_lock(&hcImageFITS::mutexFits);
	if (fits_write_comment(filePtr, comment.data(), &status))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": comment '" << comment << "' cannot be written.\n";
		printerror(status);
		result = false;
	}
	pthread_mutex_unlock(&hcImageFITS::mutexFits);

	return result;
}

bool hcImageFITS::save(const string &filename)
{
	bool retval			= true;
	long numElements 	= width * height;
	long fpixel 		= 1;
	long naxis 			= 2;
	long naxes[2] 		= { width, height };
	int status 			= 0;
	int bitpix 			= FLOAT_IMG;

	pthread_mutex_lock(&hcImageFITS::mutexFits);
	if(filePtr == NULL || filename != this->filename)
	{
		remove(filename.data());
		if (fits_create_file(&filePtr, filename.data(), &status))
		{
			printf("ERROR! hcImageFITS::writeFile: Creating file failed!\n");
			printerror(status);
			retval = false;
		}
		else if (fits_create_img(filePtr, bitpix, naxis, naxes, &status))
		{
			printf("ERROR! hcImageFITS::writeFile: Creating image failed!\n");
			printerror(status);
			retval = false;
		}

		if(!retval)
			return false;
	}

	float *outArray = new float[numElements];

	for(uint x=0; x<width;++x)
		for(uint y=0; y<height; ++y)
			outArray[y * width + x] = operator()(x,y);

	if (fits_write_img(filePtr, TFLOAT, fpixel, numElements, outArray, &status))
	{
		printf("ERROR! hcImageFITS::writeFile: writing image failed\n");
		printerror(status);
		retval = false;
	}

	delete [] outArray;

	if (fits_close_file(filePtr, &status))
	{
		printf("ERROR! hcImageFITS::writeFile: closing file failed\n");
		printerror(status);
		retval = false;
	}
	else	filePtr = NULL;
	pthread_mutex_unlock(&hcImageFITS::mutexFits);

	if(!retval)	return false;

	return hcImageFITS::load(filename);
}

hcFloat hcImageFITS::getHistogramMaximum(uint numBins)
{
	percentiles perc = getPercentiles();
	//getMedian(*this, perc);
	hcFloat maxVal		= perc.perc90;
	hcFloat dVal		= maxVal/numBins;
	uint *binCount		= new uint[numBins];

	for(uint i=0; i<numBins; ++i)
		binCount[i] = 0;

	for(uint x=0; x<width; ++x)
		for(uint y=0;y<height; ++y)
		{
			hcFloat val	= operator ()(x,y);
			if(val >= 0.0 && val <= maxVal && fabs(val-1E-20) > 1E-6)
			{
				uint ind = floor(val/dVal);
				++binCount[ind];
			}
		}

	uint max 	= 0;
	uint maxI	= 0;

	for(uint i=0; i<numBins; ++i)
		if(binCount[i] > max)
		{
			max 	= binCount[i];
			maxI	= i;
		}

	hcFloat retval = dVal*(2*maxI+1)/2.0;

	delete [] binCount;

	return retval;
}

ImageStatistics hcImageFITS::getImageStatistics() const
{
	ImageStatistics retval;
	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
			if(data[y * width + x] != 0.0)	++retval.numPixels;

	retval.perc	= getPercentiles();
	getMean(*this, retval.mean, retval.stddev);

	return retval;
}

percentileDataStruct hcImageFITS::getPercentiles() const
{
	percentileDataStruct retval;

	uint count = 0;
	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			if(fabs(content(x,y) - 1E-20)>1E-6 && content(x,y) > 0.0)
				++count;
		}

	hcFloat *arr	= new hcFloat[count];
	count			= 0;

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
			if(fabs(content(x,y) - 1E-20)>1E-6 && content(x,y) > 0.0)
				arr[count++] = content(x,y);

	if(count == 0)
	{
		cerr << "hcimageFITS::getPercentiles: no data to sort\n";
		return retval;
	}

	sort(arr, count);

	retval.perc00 = arr[0];
	retval.perc01 = arr[(int)(count*0.01)];
	retval.perc02 = arr[(int)(count*0.02)];
	retval.perc05 = arr[(int)(count*0.05)];
	retval.perc10 = arr[(int)(count*0.10)];
	retval.perc20 = arr[(int)(count*0.20)];
	retval.perc25 = arr[(int)(count*0.25)];
	retval.perc30 = arr[(int)(count*0.30)];
	retval.perc40 = arr[(int)(count*0.40)];
	retval.perc50 = arr[(int)(count*0.50)];
	retval.perc60 = arr[(int)(count*0.60)];
	retval.perc70 = arr[(int)(count*0.70)];
	retval.perc75 = arr[(int)(count*0.75)];
	retval.perc80 = arr[(int)(count*0.80)];
	retval.perc90 = arr[(int)(count*0.90)];
	retval.perc95 = arr[(int)(count*0.95)];
	retval.perc98 = arr[(int)(count*0.98)];
	retval.perc99 = arr[(int)(count*0.99)];
	retval.perc100 = arr[count-1];

	return retval;
}

void hcImageFITS::normalize()
{
	hcFloat omin	= numeric_limits<hcFloat>::infinity();
	hcFloat omax	= -numeric_limits<hcFloat>::infinity();

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			hcFloat val = operator ()(x,y);
			if(val > omax)
				omax = val;
			if(val < omin)
				omin = val;
		}

	hcFloat nmax	= 1.0;
	hcFloat nmin	= 0.0;
	hcFloat tf 		= (nmax-nmin)/(omax-omin);

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			hcFloat val 	= operator ()(x,y);
			hcFloat nval	= (val-omin) * tf + nmin;
			operator()(x,y) = nval;
		}
}

void hcImageFITS::rescale(uint newWidth, uint newHeight)
{
	hcImageFITS newImage;
	rescaleImage(*this, newImage, newWidth, newHeight, FILTER_LANCZOS3, true);
	init(newWidth, newHeight);

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
			this->operator ()(x,y) = newImage(x,y);
}

hcFloat hcImageFITS::crosscor(const hcImageFITS &other, uint i, uint j)
{
	hcImageFITS otherCpy(other);

	if(other.width != width || other.height != height)
	{
		cout << __FILE__ << ":" << __LINE__ << " dimension of other image (" << other.width << "/" << other.height << ") does not match. Rescale\n";
		otherCpy.rescale(width, height);
	}

	hcFloat thisAve 	= 0.0;
	hcFloat otherAve	= 0.0;

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			hcFloat thisVal 	= operator()(x,y);
			hcFloat otherVal	= otherCpy(x,y);

			if(!isnan((float)thisVal))
				thisAve += thisVal;

			if(!isnan((float)otherVal))
				otherAve += otherVal;
		}

	thisAve 	/= width * height;
	otherAve	/= width * height;

	hcFloat fDenom	= 0.0f;
	hcFloat gDenom	= 0.0f;

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			hcFloat f	 = operator()(x,y) 	- thisAve;
			hcFloat g	 = otherCpy(x,y) 	- otherAve;
			f *= f;
			g *= g;

			fDenom	+= f;
			gDenom	+= g;
		}

	hcFloat denom	= sqrt(fDenom * gDenom);
	hcFloat correl	= 0.0;

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
			correl += (operator()(x, y) - thisAve) * (otherCpy(x,y) - otherAve);

	correl /= denom;

	return correl;
}

/*
void hcImageFITS::fft(hcImageFITS &real, hcImageFITS &imag)
{
	double *in 			= new double[width*height];
	fftw_complex *out 	= (fftw_complex*) 	fftw_malloc(sizeof(fftw_complex) * width * height);

	real.init(width, height);
	imag.init(width, height);

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
			//in[y*width+x] = data[y*width+x];
			in[y + height*x] = data[y*width+x];


	fftw_plan p			= fftw_plan_dft_r2c_2d(width, height, in, out, FFTW_ESTIMATE);

	fftw_execute(p);
	fftw_destroy_plan(p);

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			real(x,y) = out[y + height*x][0];
			imag(x,y) = out[y + height*x][1];
		}

	fftw_free(out);
	delete [] in;
}

bool hcImageFITS::fft_inv(hcImageFITS &real, hcImageFITS &imag)
{
	if(real.width != imag.width || real.height != imag.height)
	{
		cout << __FILE__ << ":" << __LINE__ << ": dimensions of real and imaginary image do not match. Width: "
				<< real.width << "/" << imag.width << " height: " << real.height << "/" << imag.height << "\n";
		return false;
	}

	init(real.width, real.height);

	double *out 		= new double[width*height];
	fftw_complex *in 	= (fftw_complex*) 	fftw_malloc(sizeof(fftw_complex) * width * height);

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			in[y + height*x][0] = real(x,y);
			in[y + height*x][1] = imag(x,y);
		}


	fftw_plan p			= fftw_plan_dft_c2r_2d(width, height, in, out, FFTW_ESTIMATE);

	fftw_execute(p);
	fftw_destroy_plan(p);

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
			operator()(x,y) = out[y + height*x] / (width*height);

	fftw_free(in);
	delete [] out;

	return true;
}//*/

void hcImageFITS::initStaticMembers()
{
	mutexFits =  PTHREAD_MUTEX_INITIALIZER;
}

bool rescaleImage(hcImageFITS &dataIn, hcImageFITS &dataOut, uint widthOut, uint heightOut, FREE_IMAGE_FILTER filter, bool debug)
{
	uint widthIn	= dataIn.width;
	uint heightIn	= dataIn.height;
	float *imgCopy1 = new float[widthIn*heightIn];
	float *imgCopy2 = new float[widthOut*heightOut];

	dataOut.init(widthOut, heightOut);

	for(uint x=0; x<widthIn; ++x)
		for(uint y=0; y<heightIn; ++y)
			imgCopy1[y*widthIn+x] = dataIn(x,y);

	FIBITMAP *imgIn		= FreeImage_AllocateT(FIT_FLOAT, widthIn, heightIn);
	uint inRowBytes 	= widthIn  * sizeof(float);

	for (int i=0; i<(int)heightIn; ++i)
	{
	    BYTE* ptr2Line = FreeImage_GetScanLine(imgIn,i);
	    memcpy(ptr2Line, reinterpret_cast<BYTE*>(imgCopy1)+(inRowBytes*i), inRowBytes);
	}

	FIBITMAP *imgOut 	= FreeImage_Rescale(imgIn, widthOut, heightOut, filter);
	uint outRowBytes	= widthOut * sizeof(float);

	for (int i=0; i<(int)heightOut; ++i)
	{
		BYTE* ptr2Line = FreeImage_GetScanLine(imgOut,i);
		memcpy(reinterpret_cast<BYTE*>(imgCopy2)+(outRowBytes*i), ptr2Line, outRowBytes);
	}

	for(uint x=0; x<widthOut; ++x)
		for(uint y=0; y<heightOut; ++y)
			 dataOut(x,y) = imgCopy2[y*widthOut+x];

	delete [] imgCopy1;
	delete [] imgCopy2;

	return true;
}

/*!	TODO:	this implementation works like a box filter which seems to be quite bad for
 * 			rescaling magneograms. Also take a look at the definition of phi-coordinates. According to
 * 			Todd Hoeksema, pixel coordinates in FITS-files denote the CENTER-positions, whereas here
 * 			we assume the leftmost pixelboundary of the first pixel to be the boundary of the magnetogram:
 *
 * 				0				2PI
 * 				|---magnetogram---|
 * 			   |-|-|-|-|-|-|-|-|-|		<-- according to Todd
 * 			   	|-|-|-|-|-|-|-|-|-|		<-- according to what I would think makes more sense
 *
 * @param linSinLat 	interpolation/integration linear in sine(latitude) (true) or linear in latitude (false)
 *
 * 	  theta_0 	-> north-most pixel (near theta=0.0)
 * 	  theta_N-1 -> south-nmost pixel(near theta=PI)
 *
 *
 *    tb_m  ---------    	(theta boundary minus - upper pixel boundary)
 *    				|
 *    				|
 *    				|
 *    t			-----		(theta - position of pixel)
 *    			    |
 *    			    |
 *    			    |
 *    tb_p  ---------	    (theta boundary plus - lower pixel boundary)
 *
 */
bool rescaleSinLatGrid(hcImageFITS &imgIn,  hcFloat maxSLin,
					   hcImageFITS &imgOut, uint wOut, uint hOut, hcFloat maxSLout, bool linSinLat, bool debug)
{
#define MAX(a, b) ((a) > (b) ? ( a) : (b))
#define MIN(a, b) ((a) > (b) ? ( b) : (a))

	uint wIn	= imgIn.width;
	uint hIn	= imgIn.height;

	imgOut.init(wOut, hOut);

	hcFloat minSLin			= -maxSLin;
	hcFloat minSLout		= -maxSLout;

	hcFloat dSinLatO		= (maxSLout - minSLout) / (hOut-1);
	hcFloat dSinLatI		= (maxSLin  - minSLin)	/ (hIn-1);
	hcFloat northBorderO 	= maxSLout + dSinLatO * 1.0/2.0;
	hcFloat southBorderO 	= minSLout - dSinLatO * 1.0/2.0;

	//hcFloat thetaMinI		= PI/2.0 - asin(maxSLin);
	//hcFloat thetaMaxI		= PI/2.0 - asin(minSLin);

	//hcFloat thetaMaxO		= PI/2.0 - asin(minSLout);
	//hcFloat thetaMinO		= PI/2.0 - asin(maxSLout);

	hcFloat maxCosThetaI	= maxSLin;
	hcFloat minCosThetaI	= minSLin;
	hcFloat maxCosThetaO	= maxSLout;
	//hcFloat minCosThetaO	= minSLout;

	hcFloat dCosThetaI		= -dSinLatI;
	hcFloat dCosThetaO		= -dSinLatO;

	if(northBorderO > 1.0 || southBorderO < -1.0)
	{
		printf("ERROR! rescaleSinLatGrid: Parameters not valid!\n");
		printf("sinLat of northern image border: %E\nsinLat of southern image border: %E\n", northBorderO, southBorderO);
		printf("maxSinLat: %E\nminSinLat: %E\ndSinLat: %E\n", maxSLout, minSLout, dSinLatO);
		printf("hIn: %u, hout: %u\n", hIn, hOut);
		return false;
	}

	bool interpolateTheta	= fabs(dCosThetaI) >= fabs(dCosThetaO);	// new image has higher resolution in (cos(theta))-direction
	bool interpolatePhi		= 2 * PI / wIn >= 2 * PI / wOut;		// new image has higher resolution in phi - direction

	for(uint nThetaO=0; nThetaO<hOut; ++nThetaO)
	{
		for(uint nPhiO=0; nPhiO<wOut; ++nPhiO)
		{
			//uint indO				= (hOut - nThetaO - 1) * wOut + nPhiO;

			hcFloat cosThetaO		= maxCosThetaO + (nThetaO          ) * dCosThetaO;
			hcFloat upperPixBoundO	= maxCosThetaO + (nThetaO - 1.0/2.0) * dCosThetaO;	//cos(theta) of upper pixel boundary
			hcFloat lowerPixBoundO	= maxCosThetaO + (nThetaO + 1.0/2.0) * dCosThetaO;	//cos(theta) of lower pixel boundary

			hcFloat phiO			= (nPhiO + 1.0/2.0) * 2 * PI / wOut;	// WARNING: see comments in the header of this function
			hcFloat leftPixBoundO	=  nPhiO        	* 2 * PI / wOut;	// WARNING: see comments in the header of this function
			hcFloat rightPixBoundO	= (nPhiO + 1.0) 	* 2 * PI / wOut;	// WARNING: see comments in the header of this function

			if(lowerPixBoundO > maxCosThetaI || upperPixBoundO < minCosThetaI)	// 		lower pixel boundary of out-img lies north of northernmost pixel of in-img
			{																	// OR   upper pixel boundary of out-img lies south of southernmost pixel in in-img
				imgOut(nPhiO, hOut - nThetaO - 1) = 0.0;
				continue;
			}

			if(interpolateTheta && interpolatePhi)	// new image has a higher resolution in both (cos(theta))- and phi-direction
			{
				uint nThetaI 	= 1;
				uint nPhiI		= 1;

				while(maxCosThetaI + (nThetaI          ) * dCosThetaI > cosThetaO && ++nThetaI < hIn-1);

				hcFloat cosTheta_m	= maxCosThetaI + (nThetaI - 1.0    ) * dCosThetaI;
				hcFloat cosTheta_p	= maxCosThetaI + (nThetaI          ) * dCosThetaI;
				hcFloat dctm		= linSinLat																		?
									  (cosThetaO - cosTheta_m) 				/ dCosThetaI							:	// relative delta cos(theta) to midpoint
									  (acos(cosThetaO) - acos(cosTheta_m))  / (acos(cosTheta_p) - acos(cosTheta_m));	// relative delta theta to midpoint

				while((nPhiI + 1.0/2.0) * 2 * PI / wIn    < phiO && ++nPhiI < wIn);
				hcFloat phi_m		= (nPhiI - 1.0/2.0) * 2 * PI / wIn;
				hcFloat dpm			= (phiO      - phi_m)  					/ (2 * PI / wIn);							// relative delta phi to midpoint
				hcFloat val_mm		= imgIn(nPhiI - 1				, hIn - 1 - (nThetaI - 1)	);
				hcFloat val_pm		= imgIn(nPhiI - 1				, hIn - 1 - (nThetaI)		);
				hcFloat val_mp		= imgIn((nPhiI==wIn ? 0 : nPhiI), hIn - 1 - (nThetaI - 1)	);
				hcFloat val_pp		= imgIn((nPhiI==wIn ? 0 : nPhiI), hIn - 1 - (nThetaI)    	);

				hcFloat valMid_m	= (1 - dctm) * val_mm    + dctm * val_pm;
				hcFloat valMid_p	= (1 - dctm) * val_mp    + dctm * val_pp;

				hcFloat val			= (1 - dpm)  * valMid_m	 + dpm  * valMid_p;

				imgOut(nPhiO, hOut - nThetaO - 1) = val;
			}
			else if(!interpolateTheta && !interpolatePhi)	// new image has lower resolution in both theta- and phi-direction
			{
				hcFloat val 	= 0.0;
				hcFloat denom	= 0.0;

				for(uint nThetaI=0; nThetaI<hIn; ++nThetaI)
				{
					hcFloat upperPixBoundI	= maxCosThetaI + (nThetaI - 1.0/2.0) * dCosThetaI;
					hcFloat lowerPixBoundI	= maxCosThetaI + (nThetaI + 1.0/2.0) * dCosThetaI;

					if(lowerPixBoundI > upperPixBoundO || upperPixBoundI < lowerPixBoundO)
						continue;

					for(uint nPhiI=0; nPhiI<wIn; ++nPhiI)
					{
						hcFloat leftPixBoundI	=  nPhiI  		    * 2 * PI / wIn;
						hcFloat rightPixBoundI	= (nPhiI + 1.0    ) * 2 * PI / wIn;
						hcFloat pixValI			= imgIn(nPhiI ,hIn - 1 - nThetaI);

						if(rightPixBoundI <= leftPixBoundO || leftPixBoundI >= rightPixBoundO)
							continue;

						hcFloat dPhi 		= 	rightPixBoundO >= rightPixBoundI 										?
												MIN(rightPixBoundI-leftPixBoundO, rightPixBoundI-leftPixBoundI) 		:
												rightPixBoundO - leftPixBoundI;

						hcFloat dCosTheta	= 	lowerPixBoundO <= lowerPixBoundI 										?
												MIN(upperPixBoundO - lowerPixBoundI, upperPixBoundI - lowerPixBoundI) 	:
												upperPixBoundI - lowerPixBoundO;

						if(!linSinLat)
							dCosTheta		= 	lowerPixBoundO <= lowerPixBoundI 										?
												MIN(acos(lowerPixBoundI)-acos(upperPixBoundO), acos(lowerPixBoundI) - acos(upperPixBoundI)):
												acos(lowerPixBoundO) - acos(upperPixBoundI);

						hcFloat area		= dCosTheta * dPhi;
						hcFloat addval		= area * pixValI;
						val					+= addval;
						denom				+= area;
					}
				}

				val 		/= denom;
				imgOut(nPhiO, hOut - nThetaO - 1) = val;
			}
			else if(interpolateTheta && !interpolatePhi) 	// new image has higher resolution in theta direction but lower res. in phi direction
			{
				printf("rescaleSinLatGrid interpolate theta - this case is not working\n");
				exit(1);
				hcFloat valUp		= 0.0;
				hcFloat valDown 	= 0.0;
				hcFloat denomUp		= 0.0;
				hcFloat denomDown	= 0.0;

				for(uint nThetaI=0; nThetaI<hIn; ++nThetaI)
				{
					hcFloat cosThetaI		= maxCosThetaI + (nThetaI          ) * dCosThetaI;

					if(fabs(cosThetaO-cosThetaI) > fabs(dCosThetaI))
						continue;

					for(uint nPhiI=0; nPhiI<wIn; ++ nPhiI)
					{
						hcFloat leftPixBoundI	=  nPhiI  			* 2 * PI / wIn;
						hcFloat rightPixBoundI	= (nPhiI + 1.0) 	* 2 * PI / wIn;
						hcFloat pixValI			= imgIn(nPhiI, hIn - 1 - nThetaI);

						if(rightPixBoundI <= leftPixBoundO || leftPixBoundI >= rightPixBoundO)
							continue;

						hcFloat dPhi 		= 	rightPixBoundO >= rightPixBoundI ?
												MIN(rightPixBoundI-leftPixBoundO, rightPixBoundI-leftPixBoundI) :
												rightPixBoundO - leftPixBoundI;

						(cosThetaO - cosThetaI < 0 ? valUp 	: valDown) 	+= dPhi*pixValI;
						(cosThetaO - cosThetaI < 0 ? denomUp: denomDown)+= dPhi;
					}
				}

				valUp 	/= denomUp;
				valDown	/= denomDown;

				uint nThetaI = 1;
				while(maxCosThetaI + nThetaI * dCosThetaI > cosThetaO && ++nThetaI < hIn-1);

				hcFloat cosTheta_m	= maxCosThetaI + (nThetaI - 1.0    ) * dCosThetaI;

				hcFloat dctm		= (cosThetaO - cosTheta_m) 	/ dCosThetaI;
				hcFloat val 		= (1 - dctm) * valUp + dctm * valDown;
				imgOut(nPhiO, hOut - nThetaO - 1) = val;
			}
			else// !interpolateTheta && interpolate Phi		// new image has higher resolution in phi direction but lower res in theta direction
			{
				printf("rescaleSinLatGrid interpolate phi - this case is not working\n");
				exit(1);
				hcFloat valLeft		= 0.0;
				hcFloat valRight 	= 0.0;
				hcFloat denomLeft	= 0.0;
				hcFloat denomRight	= 0.0;

				for(uint nThetaI=0; nThetaI<hIn; ++nThetaI)
				{
					hcFloat upperPixBoundI	= maxCosThetaI + (nThetaI - 1.0/2.0) * dCosThetaI;
					hcFloat lowerPixBoundI	= maxCosThetaI + (nThetaI + 1.0/2.0) * dCosThetaI;

					if(lowerPixBoundI > upperPixBoundO || upperPixBoundI < lowerPixBoundO)
						continue;

					hcFloat dCosTheta	= 	lowerPixBoundO <= lowerPixBoundI ?
											MIN(upperPixBoundO - lowerPixBoundI, upperPixBoundI - lowerPixBoundI) :
											upperPixBoundI - lowerPixBoundO;

					for(uint nPhiI=0; nPhiI<wIn; ++ nPhiI)
					{
						hcFloat phiI			= (nPhiI + 1.0/2.0)	* 2 * PI / wIn;
						if(fabs(phiO-phiI) > 2 * PI / wIn)
							continue;

						hcFloat pixValI			= imgIn(nPhiI, hIn - 1 - nThetaI);

						(phiO-phiI > 0 ? valLeft	: valRight)		+= dCosTheta * pixValI;
						(phiO-phiI > 0 ? denomLeft	: denomRight)	+= dCosTheta;
					}
				}

				valLeft  /= denomLeft;
				valRight /= denomRight;

				uint nPhiI = 1;
				while((nPhiI + 1.0/2.0) * 2 * PI / wIn < phiO && ++nPhiI < wIn);
				hcFloat phi_m	= (nPhiI - 1.0/2.0) * 2 * PI / wIn;		// WARNING: see comments in the header of this function
				hcFloat dpm		= (phiO  - phi_m) 	/(2 * PI / wIn);	// WARNING: see comments in the header of this function

				hcFloat val		= (1 - dpm)  * valLeft	+ dpm * valRight;
				imgOut(nPhiO, hOut - nThetaO - 1) = val;
			}
		}
	}

	// the following takes care of pixels that lie outside the boundaries of the input image
	// it is just a simple averaging of neighboring pixels in the next line
	int nThetaO 		= -1;
	hcFloat cosThetaO	= 0.0;
	do
	{
		++nThetaO;
		cosThetaO = maxCosThetaO + (nThetaO          ) * dCosThetaO;
	}
	while(cosThetaO > maxCosThetaI);

	while(--nThetaO >= 0)
	{
		for(uint nPhiO = 0; nPhiO<wOut; ++nPhiO)
		{
			hcFloat val_m	= imgOut((nPhiO==0 	 	? wOut-1 	: nPhiO-1)	, hOut - 1 - (nThetaO + 1));
			hcFloat val		= imgOut(nPhiO									, hOut - 1 - (nThetaO + 1));
			hcFloat val_p	= imgOut((nPhiO==wOut-1 ? 0 		: nPhiO+1)	, hOut - 1 - (nThetaO + 1));

			imgOut(nPhiO,hOut - nThetaO - 1 ) = 1.0/3.0 * (val_m + val + val_p);
		}
	}

	while(cosThetaO > minCosThetaI)
	{
		++nThetaO;
		cosThetaO =  maxCosThetaO + nThetaO * dCosThetaO;
	}

	--nThetaO;
	while(++nThetaO < (int)hOut)
	{
		for(uint nPhiO=0; nPhiO<wOut; ++nPhiO)
		{
			hcFloat val_m	= imgOut((nPhiO==0 	 ? wOut-1 	: nPhiO-1)	, hOut - 1 - (nThetaO-1));
			hcFloat val		= imgOut(nPhiO								, hOut - 1 - (nThetaO-1));
			hcFloat val_p	= imgOut((nPhiO==wOut-1 ? 0 	: nPhiO+1)	, hOut - 1 - (nThetaO-1));

			imgOut(nPhiO, hOut - nThetaO - 1) = 1.0/3.0 * (val_m + val + val_p);
		}
	}

	return true;
}

void getMean(const hcImageFITS &img, hcFloat &mean, hcFloat &std)
{
	uint count 	= 0;
	mean		= 0.0;
	std			= 0.0;

	for(uint x=0; x<img.width; ++x)
		for(uint y=0; y<img.height; ++y)
		{
			if(fabs(img.content(x,y) - 1E-20)>1E-6 && img.content(x,y) > 0.0)
			{
				++count;
				mean += img.content(x,y);
			}
		}
	mean /= count;

	for(uint x=0; x<img.width; ++x)
		for(uint y=0; y<img.height; ++y)
			if(fabs(img.content(x,y) - 1E-20)>1E-6 && img.content(x,y) > 0)
				std	+= (img.content(x,y) - mean) * (img.content(x,y) - mean);

	std = sqrt(std/count);
}

void getMean(hcFloat *data, uint numDataPoints, hcFloat &mean, hcFloat &std)
{
	uint count 	= 0;
	mean		= 0.0;
	std			= 0.0;

	for(uint x=0; x<numDataPoints; ++x)
		if(fabs(data[x] - 1E-20)>1E-6 && data[x] > 0.0)
		//if(!isnan(data[x]))
		{
			++count;
			mean += data[x];
		}

	mean /= count;

	for(uint x=0; x<numDataPoints; ++x)
		if(fabs(data[x] - 1E-20)>1E-6 && data[x] > 0.0)
		//if(!isnan(data[x]))
			std	+= (data[x] - mean) * (data[x] - mean);

	std = sqrt(std/count);
}

void getPercentiles(hcImageFITS &img, percentiles &perc)
{
	perc = img.getPercentiles();
}
