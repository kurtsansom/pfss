#ifndef PFSS_H
#define PFSS_H
#include "engine/math/hcFunction.h"

#ifdef GUI
#include "engine/hcScene3D.h"
#include "engine/overlay.h"
#endif

#include "src/synPhotMagfield.h"

Vec3D getDipoleField(Vec3D pos, hcFloat dipole, bool debug=false);

Vec3D getDipoleCurrentField(Vec3D pos, hcFloat dipole, bool debug=false);


/*! \brief solution of the PFSS model via spherical harmonic coefficients (SHC)
 *
 */
class PFSSsolution_SHC{
public:
    AssocLegendrePoly_sun** assLegPoly; /*!< \brief contains all the associated Legendre polynomials required for the solution 			*/
    uint coeffOrder;    				/*!< \brief highest order of assoc. Legendre Polynomial utilized 								*/
    hcFloat **g;                 		/*!< \brief g coefficients 																		*/
    hcFloat **h;                 		/*!< \brief h coefficients 																		*/
    hcFloat sourceSurfaceFactor; 		/*!< \brief location of source surface in multiples of r_sol 									*/

    PFSSsolution_SHC(const char *filename, bool WSOfile, hcFloat r_ss);	/*!< \brief constructor											*/
	PFSSsolution_SHC();													/*!< \brief std constructor										*/
	PFSSsolution_SHC(const PFSSsolution_SHC &other);					/*!< \brief cpy constructor										*/
	virtual ~PFSSsolution_SHC();										/*!< \brief destructor											*/

	PFSSsolution_SHC &operator=(const PFSSsolution_SHC &other);	/*!< \brief assignment operator											*/

	void initNULL();
	void clear();
	virtual void init(uint order, hcFloat r_ss) = 0;

    virtual void eval(const Vec3D &pos, Vec3D &result) = 0 ;
    									/*!< \brief evaluates the magnetic field at spherical position pos								*/

    bool importCoefficients(const char *filename, bool WSOfile, hcFloat r_ss);
    									/*!< \brief reads coefficients from human readable file											*/

    void exportCoefficients(const char *filename);
    									/*!< \brief writes coefficients to human readable file											*/

    void determineCoefficientsFromPhotMagfield(uint order, hcFloat r_ss, SynPhotMagfield &photMagfield);
    									/*!< \brief computes spherical harmonic coefficients from synoptic magnetogram					*/
};

class PFSSsolution_SHC_sun : public PFSSsolution_SHC{
public:

	PFSSsolution_SHC_sun();												/*!< \brief std constructor										*/
	PFSSsolution_SHC_sun(const PFSSsolution_SHC_sun &other);			/*!< \brief cpy constructor										*/
	~PFSSsolution_SHC_sun();											/*!< \brief destructor											*/

	PFSSsolution_SHC_sun &operator=(const PFSSsolution_SHC_sun &other);	/*!< \brief assignment operator									*/
	PFSSsolution_SHC_sun operator-(const PFSSsolution_SHC_sun &other);	/*!< \brief difference between two solutions (difference of coefficients)	*/

	virtual void init(uint order, hcFloat r_ss);

	virtual void eval(const Vec3D &pos, Vec3D &result);
};

class PFSSsolution_SHC_hoek : public PFSSsolution_SHC{
public:

	PFSSsolution_SHC_hoek();
	PFSSsolution_SHC_hoek(const PFSSsolution_SHC_hoek &other);
	~PFSSsolution_SHC_hoek();

	PFSSsolution_SHC_hoek &operator=(const PFSSsolution_SHC_hoek &other);
	PFSSsolution_SHC_hoek operator-(const PFSSsolution_SHC_hoek &other);	/*!< \brief difference between two solutions (difference of coefficients)	*/

	virtual void init(uint order, hcFloat r_ss);

	virtual void eval(const Vec3D &pos, Vec3D &result);
};

/*! \brief solution of the CSSS model (not fully implemented/not working!!)
 *
 *  Gets as parameter a file with coefficients for associated Legendre polynomials.
 *  Returns on request the solution (magnetic field strength B) at a specific location
 *  in spherical (preferred) or cartesian coordinates
 */
class CSSS_magfield{
public:
    AssocLegendrePoly_sun** assLegPoly; /*!< \brief contains all the associated Legendre polynomials required for the solution	*/
    uint order;    						/*!< \brief highest order of assoc. Legendre Polynomial utilized 						*/
    double **g;                 		/*!< \brief coefficients 																*/
    double **h;                 		/*!< \brief coefficients 																*/
    double sourceSurfaceFactor; 		/*!< \brief location of source surface in multiples of r_sol							*/

    CSSS_magfield();								/*!< \brief std constructor													*/
    CSSS_magfield(const CSSS_magfield &other);		/*!< \brief cpy constructor													*/
    CSSS_magfield(const char *filename);
	~CSSS_magfield();								/*!< \brief destructor														*/

	CSSS_magfield &operator=(const CSSS_magfield &other);	/*!< \brief assignment operator										*/

    void initNULL();
    void clear();
    void init(uint order);

    void eval(const Vec3D &pos, Vec3D &result);

    bool load(const char *filename);
    void exportCoefficients(const char *filename);
    //void determineCoefficientsFromFITSsineLat(uint order, const char *FITS_filename, double maxSinLat, double minSinLat, bool accountForMissingFlux = PFSS_MISSINGFLUX);  /*!< computes coefficients g, h for potential function with no source surface (source surf->inf)*/

    void dump();
};

#ifdef GUI
/*! \brief visualization of one or more PFSS solutions
 *
 *  implements features for single image (3D model) and interactive animations of one or more
 *  PFSS solutions
 */
class PFSS_visualization : public hcScene3D{
public:

    hcDObject3D sun;
    hcDObject3D sources;

    Overlay overlay;
    PFSSsolution_SHC_sun magfield;
    hcDObject3D *magLines;
    double drawScale;

    uint numFrames;     		/*!< \brief number of Frames for an animation 									*/
    double lowerSS;     		/*!< \brief lowest source surface in animation 									*/
    double upperSS;     		/*!< \brief uppermost source surface in animation 								*/
    uint lowerOrder;
    uint upperOrder;
    bool animationPlaying; 		/*!< \brief states, if the animation shall be played 							*/
    uint animPointer;   		/*!< \brief points to the actual frame of the animation (position in *maglines) */
    uint fps;           		/*!< \brief frames per second 													*/
    double accumulatedTimediff; /*!< \brief time since last animation update 									*/
    int ssAnimScreenID; 		/*!< \brief object identifier of the source surface animation screen in overlay */
    int textScreenID;  			/*!< \brief object identifier of the upper surface animation screen in overlay 	*/
    int contourProjScreenID;

    PFSS_visualization(char* filename);
    PFSS_visualization();
    virtual ~PFSS_visualization();

    void init();
    void clear();

    void initFromSynopticFITSsineLat(uint order, const char *filename, double maxSinLat, double minSinLat);  /*!< \brief computes harmonic coefficients from synoptic (SOHO) FITS photospheric magnetic data */
    void initFromSynopticWSO(uint order, const char *filename);   /*!< \brief computes harmonic coefficients from synoptic (WSO) ASCII photospheric magnetic data */
    void initFromCoeffFile(const char* filename);                 /*!< \brief loads harmonic coefficients from WSO coefficient file */
    void compareCoeffFiles(const char *infile1, const char *infile2, const char *outfile);

    //void createWSOlikeContourMap(char *filename, bool microTesla, uint mapNo = 0);
    void initScene();
    int loadAnimation(const char *path);                                    /*!< \brief loads precomputed animation (of several PFSS solutions) */
    //void exportMagfieldAtSS(const char *filenameSS); /*!< \brief exports the coronal magnetic field at the source surface to FITS file */
    void exportRadialMagfieldAtHeight(float r, const char *filename);
    void exportRadialProjectedMagfieldAtHeight(float r, const char *filename);
    void draw(GLuint programID);
};
#endif

#endif // PFSS_H
