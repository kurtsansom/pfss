#ifndef SYNPHOTMAGFIELD_H
#define SYNPHOTMAGFIELD_H

#include "engine/hcTime.h"
#include "engine/hcImageFITS.h"

#include "src/carRotInfo.h"

typedef unsigned int uint;

class SphericalGrid;

bool createVerenaDailyMap(hcDate date, const char *folder);

/*! \brief information on photospheric magnetic field data such as instrument which was used, sin(latitude)-format, ...
 */
class SynopticInfo : public crInfo{
public:

	originID instrument;						/*!< instrument used for synoptic map					*/
	uint dailyID;								/*!< identification for daily synoptic maps				*/
	bool sinLatFormat;							/*!< y-Axis in synoptic map given in sin(latitude)?		*/
	hcFloat maxSinLat;

	SynopticInfo();								/*!< std constructor									*/
	SynopticInfo(const SynopticInfo &other);	/*!< cpy constructor									*/
	SynopticInfo(originID id, hcFloat maxSinLat, uint CRnum, uint dailyID=0);
	virtual ~SynopticInfo();					/*!< destructor											*/

	SynopticInfo &operator=(const SynopticInfo &other);		/*!< assignment operator					*/

	bool operator==(const SynopticInfo &other);	/*!< \brief comparison operator							*/

	bool operator>(const SynopticInfo &other);	/*!< \brief comparison operator							*/

	bool operator<(const SynopticInfo &other);	/*!< \brief comparison operator							*/

	void initNULL();
	void clear();
	void init(originID id, hcFloat maxSinLat, uint CRnum, uint dailyID);

	void dump() const;
};

/*! \brief handles synoptic photospheric magnetograms stored in GAUSS
 *
 * 	1G = 1E-4T = 100uT
 */
class SynPhotMagfield : public hcImageFITS{
public:

	SynopticInfo synInfo;

    SynPhotMagfield();                             	/*!< \brief std constructor                             */
    SynPhotMagfield(const SynPhotMagfield &other);  /*!< \brief cpy constructor                             */
    virtual ~SynPhotMagfield(){}

    SynPhotMagfield &operator=(const SynPhotMagfield &other);

    void initNULL();

    virtual bool load(		const string &filename);
    bool openKPVT(			const string &filename);
    bool openWSO(			const string &filename);
    bool openWSOconverted(	const string &filename);
    bool openSOHOMDI(		const string &filename);
    bool openSOHOMDI_DAILY(	const string &filename);
    bool openNSOGONG(		const string &filename);
    bool openSDOHMI(		const string &filename);
    bool openVerDAILY(		const string &filename);
    bool openOwnFormat(		const string &filename);

    bool createLOSfromGrid(SphericalGrid &grid);
    	/*!< \brief computes line-of-sight magnetogram from magnetic field values of grid										*/

    bool loadDipole(hcFloat dipole, uint numTheta, uint numPhi,
    				hcFloat maxSinLat=0.0, hcFloat minSinLat=-0.0);
    	/*!< \brief loads artificial dipole to line-of-sight magnetogram														*/

    void removeMonopole();
    	/*!< \brief removes monopole component from magnetogram																	*/

    bool remeshImage(uint newWidth, uint newHeight, uint scaleMethod);
    	/*!< \brief rescales magnetogram																						*/

    void cropPolesFromImage(uint numPixelsToCrop);
    	/*!< \brief remove north-most and south-most lines from magnetogram														*/

    bool convertWSOtxtToWSOfits(const string &infile, const string &outfile);
    	/*!< \brief converts Stanford ASCII files to FITS files																	*/

    virtual bool save(const string &filename);
    	/*!< \brief stores photospheric magfield in FITS format to be re-imported later on										*/

    virtual void addHomWhiteNoise(float sigma, uint seed);
		/*!< \brief adds white noise by adding gaussian distributed value														*/

	virtual void addPixelNoise(float fraction, uint seed);
		/*!< \brief adds noise to each pixel seperately by adding normal distributed values with sigma = fraction * pixelValue	*/

	virtual void addSignedFractionNoise(float fraction, uint seed);
		/*!< \brief adds +/- fraction * pixel value to each pixel																*/

    void dump() const;
};

#endif
