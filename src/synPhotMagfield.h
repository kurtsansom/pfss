#ifndef SYNPHOTMAGFIELD_H
#define SYNPHOTMAGFIELD_H

#include "engine/hcTime.h"
#include "engine/hcImageFITS.h"

#include "src/carRotInfo.h"

typedef unsigned int uint;

class SphericalGrid;

/*! \brief information on photospheric magnetic field data such as instrument which was used, sin(latitude)-format, ...
 */
class SynopticInfo : public crInfo{
public:

	originID instrument;						/*!< instrument used for synoptic map					*/
	uint dailyID;								/*!< identification for daily synoptic maps				*/
	bool sinLatFormat;							/*!< y-Axis in synoptic map given in sin(latitude)?		*/
	hcFloat maxSinLat;							/*!< sin(latitude) of northern grid boundary			*/

	SynopticInfo();								/*!< std constructor									*/
	SynopticInfo(const SynopticInfo &other);	/*!< cpy constructor									*/
	SynopticInfo(originID id, hcFloat maxSinLat, uint CRnum, uint dailyID=0);
	virtual ~SynopticInfo();					/*!< destructor											*/

	SynopticInfo &operator=(const SynopticInfo &other);		/*!< assignment operator					*/

	bool operator==(const SynopticInfo &other);	/*!< \brief comparison operator							*/

	bool operator>(const SynopticInfo &other);	/*!< \brief comparison operator							*/

	bool operator<(const SynopticInfo &other);	/*!< \brief comparison operator							*/

	bool exportBinary(ofstream &stream);	/*!< \brief export instance to stream						*/

	bool importBinary(ifstream &stream);	/*!< \brief import instance from stream						*/

	void initNULL();
	void clear();
	void init(originID id, hcFloat maxSinLat, uint CRnum, uint dailyID);

	void dump(uint indent=0) const;
		/*!< \brief dumps information on this instance to stdout										*/
};

/*! \brief handles synoptic photospheric magnetograms stored in GAUSS
 *
 * 	1G = 1E-4T = 100uT
 */
class SynPhotMagfield : public hcImageFITS{
public:

	SynopticInfo synInfo;							/*!< \brief information on the synoptic data		*/

    SynPhotMagfield();                             	/*!< \brief std constructor                         */
    SynPhotMagfield(const SynPhotMagfield &other);  /*!< \brief cpy constructor                         */
    virtual ~SynPhotMagfield(){}					/*!< \brief destructor								*/

    SynPhotMagfield &operator=(const SynPhotMagfield &other);	/*!< \brief assignment operator			*/

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
    	/*!< \brief computes line-of-sight magnetogram from magnetic field values of grid			*/

    bool loadDipole(hcFloat dipole, uint numTheta, uint numPhi,
    				hcFloat maxSinLat=0.0, hcFloat minSinLat=-0.0);
    	/*!< \brief loads artificial dipole to line-of-sight magnetogram							*/

    void removeMonopole();
    	/*!< \brief removes monopole component from magnetogram										*/

    bool remeshImage(uint newWidth, uint newHeight, uint scaleMethod);
    	/*!< \brief rescales magnetogram															*/

    void cropPolesFromImage(uint numPixelsToCrop);
    	/*!< \brief remove north-most and south-most lines from magnetogram							*/

    bool convertWSOtxtToWSOfits(const string &infile, const string &outfile);
    	/*!< \brief converts Stanford ASCII files to FITS files										*/

    virtual bool save(const string &filename);
    	/*!< \brief stores photospheric magfield in FITS format to be re-imported later on			*/

    void dump() const;
    	/*!< \brief dumps information on this instance to stdout									*/
};

#endif
