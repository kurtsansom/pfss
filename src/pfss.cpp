#include "engine/math/hcMatrix.h"
#include "engine/hcImageFITS.h"
//#include "engine/math/coordTransform.h"

#include "src/pfss.h"

//#include "extern/soil/SOIL.h"

#include <stdio.h>
#include <string.h>
#include <fstream>

#include <boost/regex.hpp>
#include "boost/lexical_cast.hpp"
#include <fstream>

using namespace std;
using namespace boost;

/*! returns the magnetic field vector at given position (spherical coordinates) for a dipole
 *  at the origin oriented along the positive z-axis. Return value given in spherical coordinates.
 *  mathematics found at http://www.physicsinsights.org/dipole_field_1.html
 */
Vec3D getDipoleField(Vec3D pos, hcFloat dipole, bool debug)
{
	hcFloat r		= pos[0] / r_sol;
	hcFloat theta	= pos[1];
	hcFloat stheta	= sin(theta);
	hcFloat ctheta	= cos(theta);
	hcFloat fac		= dipole / (r * r * r);
	hcFloat Br		= fac * 2 * ctheta;
	hcFloat Bt		= fac * stheta;

	Vec3D retval((hcFloat)Br, (hcFloat)Bt, (hcFloat)0.0f);

	if(debug)
		printf("r: %E, t: %E, fac: %E, Br: %E, Bt: %E\n", r, theta, fac, Br, Bt);

	return retval;
}

/*! mathematics found in "The Model Magnetic Configuration of the Extended Corona
in the Solar Wind Formation Region - Igor Veselovsky and Olga Panasenco"
 */
Vec3D getDipoleCurrentField(Vec3D pos, hcFloat dipole, bool debug)
{
	hcFloat r		= pos[0] / r_sol;
	hcFloat theta	= pos[1];
	hcFloat stheta	= sin(theta);
	hcFloat ctheta	= cos(theta);
	hcFloat fac		= dipole / (r * r * r);

	hcFloat mu0		= 4*PI*1E-7;
	hcFloat mu		= mu0 * dipole;
	hcFloat r0		= 2.5;
	hcFloat B0		= -1E-8;
	hcFloat Ba		= (B0*r0*r0) * (B0*r0*r0) * (B0*r0*r0) / (4 * mu * mu);
	hcFloat a		= 2 * mu / (B0*r0*r0);
	hcFloat rho		= r / a;

	hcFloat Br		= Ba / (rho * rho) * (ctheta / rho + (theta - PI/2.0 < 0 ? -1.0 : 1.0));
	hcFloat Bt		= 1 / (2 * rho * rho * rho) * Ba * stheta;

	if(debug)
	{
		printf("r: %E, t: %E, Ba: %E, a: %E, rho: %E, Br: %E, Bt: %E\n", r, theta, Ba, a, rho, Br, Bt);
		//printf("(B0*r0*r0): %E, ^3: %E, denom: %E, ba: %E\n", (B0*r0*r0), (B0*r0*r0) * (B0*r0*r0) * (B0*r0*r0), (4 * dipole * dipole), Ba);
	}

	Vec3D retval((hcFloat)Br, (hcFloat)Bt, (hcFloat)0.0f);

	return retval;
}

PFSSsolution_SHC::PFSSsolution_SHC(const char *filename, bool WSOfile, hcFloat r_ss)
{
	initNULL();
	importCoefficients(filename, WSOfile, r_ss);
}

PFSSsolution_SHC::PFSSsolution_SHC()
{
    initNULL();
}

PFSSsolution_SHC::~PFSSsolution_SHC()
{
    clear();
}

PFSSsolution_SHC::PFSSsolution_SHC(const PFSSsolution_SHC &other)
{
	printf("PFSS_magfield: cpy constructor not implemented!\n");
	initNULL();
}

PFSSsolution_SHC &PFSSsolution_SHC::operator=(const PFSSsolution_SHC &other)
{
	if(this == &other)
		return *this;

	printf("PFSSsolution_SHC: assignment operator not implemented!\n");

	return *this;
}

void PFSSsolution_SHC::initNULL()
{
	g 					= NULL;
	h 					= NULL;
	assLegPoly 			= NULL;
	coeffOrder 			= 0;
}

void PFSSsolution_SHC::clear()
{
	if(g!=NULL)
	{
		for(uint i=0;i<coeffOrder;++i)
		{
			delete [] g[i];
			delete [] h[i];
			delete [] assLegPoly[i];
		}
		delete [] g;
		delete [] h;
		delete [] assLegPoly;
	}

	initNULL();
}

/*! computes spherical harmonic coefficients from synoptic magnetogram
 *
 * 	@param 	order			maximum order of coefficients to be computed
 * 	@param	r_ss			heliocentric position of source surface (in m)
 * 	@param	photMagfield	data structure containing synoptic photospheric (LOS) magnetic field strength
 */
void PFSSsolution_SHC::determineCoefficientsFromPhotMagfield(uint order, hcFloat r_ss, SynPhotMagfield &photMagfield)
{
    init(order, r_ss);

    uint width 			= photMagfield.width;
    uint height 		= photMagfield.height;
    hcFloat	maxSinLat	= photMagfield.synInfo.maxSinLat;
    hcFloat minSinLat	= -photMagfield.synInfo.maxSinLat;
    hcFloat numerator 	= 1.0 / (width * height);

    for(uint l=0;l<=coeffOrder;++l)
        for(uint m=0;m<=l;++m)
        {
            g[l][m] = 0.0;
            h[l][m] = 0.0;
            for(uint y=0;y<height;++y)
            {
            	hcFloat dSinTheta 	= (maxSinLat - minSinLat) / (height-1);
                hcFloat theta		= PI / 2.0 - asin(maxSinLat - (height - 1 - y) * dSinTheta);
                hcFloat legPolyEval = assLegPoly[l][m].operator ()(cos(theta));
                hcFloat cosecans	= 1.0 / sin(theta);

                for(uint x=0;x<width;++x)
                {
                    hcFloat phi		= m * (x+1.0/2.0) * 2.0 * PI / width;
					hcFloat D_xy 	= photMagfield(x,y) * cosecans;
					hcFloat gm		= D_xy * cos(phi) * legPolyEval;
					hcFloat hm		= D_xy * sin(phi) * legPolyEval;
                    g[l][m] 		+= gm;
                    h[l][m] 		+= hm;
                }
            }

            /*//according to Altschuler & Newkirk for source surf at inf
            g[l][m] *= (2.0*l + 1) / (l + 1) * numerator;
            h[l][m] *= (2.0*l + 1) / (l + 1) * numerator;

            /*/

            // according to Xudong Sun:
            g[l][m] *= (2.0*l + 1) * numerator;
            h[l][m] *= (2.0*l + 1) * numerator;//*/
        }
}

bool PFSSsolution_SHC::importCoefficients(const char *filename, bool WSOfile, hcFloat r_ss)
{
    FILE* coeffFile;
    char str[1000];
    char tempstr[100];
    uint k,z,order;
    int highestOrder = 0;

	if (!doesFileExist(filename))
	{
		printf("ERROR! PFSSsolution_SHC::load input file '%s' does not exist\n", filename);
		init(0, 2.5*r_sol);
		return false;
	}

    coeffFile = fopen(filename, "r");

    uint i = 0;
    while(strcmp(fgets(str, 1000, coeffFile), "\n"))
    {
        while(str[i++] != ' ');
        strncpy(tempstr, str, --i);
        tempstr[i] = '\0';
        if (atoi(tempstr) > highestOrder)
            highestOrder = atoi(tempstr);
    }

    init((uint)highestOrder, r_ss);

    fseek(coeffFile,0, SEEK_SET);
    uint j=0;
    while (fgets(str, 1000, coeffFile) != NULL)
    {
        if(j==0)        // load g's
        {
            if(str[0] == '\n')
                ++j;
            else
            {
                i = 0;

                tempstr[0] = '\0';
                while(str[i++] != ' ');
                z = i;
                strncpy(tempstr, str, i-1);
                tempstr[i-1] = '\0';
                order=(uint)atoi(tempstr);
                k=0;

                while(k<=order)
                {
                    while(str[i] == ' ' || str[i] == '\t')
                        ++i;
                    z=i;
                    while(str[i] != ' ' && str[i++] != '\n');
                    strncpy(tempstr, str+z, i-1);
                    tempstr[i-z] = '\0';
                    g[order][k++] = strtof(tempstr, NULL);
                    tempstr[0] = '\0';
                    z=i;
                }
                ++order;
            }
        }
        else            // load h's
        {
            if(str[0] == '\n')
                break;
            else
            {
                i = 0;

                tempstr[0] = '\0';
                while(str[i++] != ' ');
                z = i;
                strncpy(tempstr, str, i-1);
                tempstr[i-1] = '\0';
                order=(uint)atoi(tempstr);
                k=0;

                while(k<=order)
                {
                    while(str[i] == ' ' || str[i] == '\t')
                        ++i;
                    z=i;
                    while(str[i] != ' ' && str[i++] != '\n');
                    strncpy(tempstr, str+z, i-1);
                    tempstr[i-1-z] = '\0';
                    h[order][k++] = strtof(tempstr, NULL);
                    tempstr[0] = '\0';
                    z=i;
                }

                ++order;
            }
        }

    }
    fclose(coeffFile);

    if(WSOfile)				// WSO files employ microTesla as units
    {
    	for(uint l=0; l<=coeffOrder; ++l)
    		for(uint m=0; m<=l; ++m)
    		{
    			g[l][m] /= 100;
    			h[l][m] /= 100;
    		}
    }

    return true;
}

void PFSSsolution_SHC::exportCoefficients(const char *filename)
{
    char tempstr[10000];
    char tmp[100];

    FILE *outfile = fopen(filename, "w");
    for(uint l=0;l<=coeffOrder;++l)
    {
        tempstr[0] = '\0';
        sprintf(tempstr, "%u ", l);
        for(uint m=0;m<=l;++m)
        {
            tmp[0] = '\0';
            sprintf(tmp, "%f ", g[l][m]);
            strcat(tempstr, tmp);
        }
        if(l==0)
            strcat(tempstr, "\t\t g's");
        strcat(tempstr, "\n");
        fprintf(outfile, "%s", tempstr);
    }
    fprintf(outfile, "\n");
    for(uint l=0;l<=coeffOrder;++l)
    {
        tempstr[0] = '\0';
        sprintf(tempstr, "%u ", l);
        for(uint m=0;m<=l;++m)
        {
            tmp[0] = '\0';
            sprintf(tmp, "%f ", h[l][m]);
            strcat(tempstr, tmp);
        }
        if(l==0)
            strcat(tempstr, "\t\t h's");
        strcat(tempstr, "\n");
        fprintf(outfile, "%s", tempstr);
    }
    fclose(outfile);
}

//---------------------------------------------------------------------------------------------------------------
//---   PFSSsolution_SHC_sun
//---------------------------------------------------------------------------------------------------------------

PFSSsolution_SHC_sun::PFSSsolution_SHC_sun()
{
	initNULL();
}

PFSSsolution_SHC_sun::PFSSsolution_SHC_sun(const PFSSsolution_SHC_sun &other)
{
	initNULL();
	operator =(other);
}

PFSSsolution_SHC_sun::~PFSSsolution_SHC_sun()
{
	clear();
}

PFSSsolution_SHC_sun &PFSSsolution_SHC_sun::operator=(const PFSSsolution_SHC_sun &other)
{
	if(this == &other)
		return *this;

	init(other.coeffOrder, other.sourceSurfaceFactor*r_sol);

	for(uint l=0; l<=coeffOrder; ++l)
		for(uint m=0; m<=l; ++m)
		{
			g[l][m]	= other.g[l][m];
			h[l][m]	= other.h[l][m];
		}

	return *this;
}

PFSSsolution_SHC_sun PFSSsolution_SHC_sun::operator-(const PFSSsolution_SHC_sun &other)
{
	PFSSsolution_SHC_sun retval;

	if(this->coeffOrder != other.coeffOrder || this->sourceSurfaceFactor != other.sourceSurfaceFactor)
	{
		printf("PFSS_magfield::operator- order does not match (%u != %u) or r_ss is different.\n", this->coeffOrder, other.coeffOrder);
		retval.init(0, 2.5*r_sol);
	}
	else
	{
		retval.init(this->coeffOrder, this->sourceSurfaceFactor);
		for(uint l=0; l<=coeffOrder; ++l)
			for(uint m=0; m<=l; ++m)
			{
				hcFloat g1 		= fabs(this->g[l][m] - other.g[l][m]) / fabs(this->g[l][m]);
				hcFloat g2 		= fabs(this->g[l][m] - other.g[l][m]) / fabs(other.g[l][m]);
				hcFloat h1 		= fabs(this->h[l][m] - other.h[l][m]) / fabs(this->h[l][m]);
				hcFloat h2 		= fabs(this->h[l][m] - other.h[l][m]) / fabs(other.h[l][m]);

				retval.g[l][m] 	= this->g[l][m]==0 ? 0.0 : g1 < g2 ? g1 : g2;
				retval.h[l][m] 	= this->h[l][m]==0 ? 0.0 : h1 < h2 ? h1 : h2;
			}
	}

	return retval;

}

void PFSSsolution_SHC_sun::init(uint order, hcFloat r_ss)
{
	clear();
	coeffOrder 	= order;
	sourceSurfaceFactor = r_ss/r_sol;
	g 			= new hcFloat*[coeffOrder + 1];
	h 			= new hcFloat*[coeffOrder + 1];
	assLegPoly 	= new AssocLegendrePoly_sun*[coeffOrder + 1];

	for(uint i=0;i<=coeffOrder;++i)
	{
		g[i] 			= new hcFloat[i+1];
		h[i] 			= new hcFloat[i+1];
		assLegPoly[i] 	= new AssocLegendrePoly_sun[i+1];

		for(uint j=0;j<=i;++j)
		{
			AssocLegendrePoly_sun temp(i, j);
			assLegPoly[i][j] = temp;
		}
	}
}

void PFSSsolution_SHC_sun::eval(const Vec3D &pos, Vec3D &result)
{
	hcFloat r			= pos[0];
	hcFloat theta		= pos[1];
	hcFloat phi			= pos[2];

	hcFloat r_ss 		= sourceSurfaceFactor * r_sol;     // source surface
	hcFloat retval_r 	= 0.0;
	hcFloat retval_t 	= 0.0;
	hcFloat retval_p 	= 0.0;

	for(uint l=0;l<=coeffOrder;++l)
		for(uint m=0;m<=l;++m)
		{   // equation 16 from Notes on PFSS Extrapolation, Xudong Sun
			retval_r += assLegPoly[l][m](cos(theta))
					  * ( g[l][m] * cos(m*phi) + h[l][m] * sin(m*phi) )
					  * pow(r_sol / r, (hcFloat)(l+2))
					  * (l+1+l*pow(r / r_ss, (hcFloat)(2*l+1))) / (l+1+l*pow(r_sol / r_ss, (hcFloat)(2*l+1)));

			// equation 17
			retval_t += -assLegPoly[l][m].getDeriv(theta)
					  * ( g[l][m] * cos(m*phi) + h[l][m] * sin(m*phi) )
					  * pow(r_sol / r, (hcFloat)(l+2))
					  * (1 - pow(r / r_ss, (hcFloat)(2*l+1))) / (l+1+l*pow(r_sol / r_ss, (hcFloat)(2*l+1)));

			// equation 18
			retval_p += assLegPoly[l][m](cos(theta))
					  * ( g[l][m] * sin(m*phi) - h[l][m] * cos(m*phi) )
					  * pow(r_sol / r, (hcFloat)(l+2))
					  * (1 - pow(r / r_ss, (hcFloat)(2*l+1))) / (l+1+l*pow(r_sol / r_ss, (hcFloat)(2*l+1)));
		}

	result(0) = retval_r;
	result(1) = retval_t;
	result(2) = retval_p;
}


//---------------------------------------------------------------------------------------------------------------
//---   PFSSsolution_SHC_hoek
//---------------------------------------------------------------------------------------------------------------

PFSSsolution_SHC_hoek::PFSSsolution_SHC_hoek()
{
	initNULL();
}

PFSSsolution_SHC_hoek::PFSSsolution_SHC_hoek(const PFSSsolution_SHC_hoek &other)
{
	initNULL();
	operator =(other);
}

PFSSsolution_SHC_hoek::~PFSSsolution_SHC_hoek()
{
	clear();
}

PFSSsolution_SHC_hoek &PFSSsolution_SHC_hoek::operator=(const PFSSsolution_SHC_hoek &other)
{
	if(this == &other)
		return *this;

	printf("PFSSsolution_SHC: assignment operator not implemented!\n");

	return *this;
}

PFSSsolution_SHC_hoek PFSSsolution_SHC_hoek::operator-(const PFSSsolution_SHC_hoek &other)
{
	PFSSsolution_SHC_hoek retval;

	if(this->coeffOrder != other.coeffOrder)
	{
		printf("PFSS_magfield::operator- order does not match (%u != %u)\n", this->coeffOrder, other.coeffOrder);
		retval.init(0, 2.5*r_sol);
	}
	else
	{
		retval.init(this->coeffOrder, this->sourceSurfaceFactor);
		for(uint l=0; l<=coeffOrder; ++l)
			for(uint m=0; m<=l; ++m)
			{
				retval.g[l][m] = this->g[l][m]==0 ? 0.0 : (this->g[l][m] - other.g[l][m]) / this->g[l][m];
				retval.h[l][m] = this->h[l][m]==0 ? 0.0 : (this->h[l][m] - other.h[l][m]) / this->h[l][m];
			}
	}

	return retval;

}

void PFSSsolution_SHC_hoek::init(uint order, hcFloat r_ss)
{
	clear();
	coeffOrder 	= order;
	sourceSurfaceFactor = r_ss / r_sol;
	g 			= new hcFloat*[coeffOrder + 1];
	h 			= new hcFloat*[coeffOrder + 1];
	assLegPoly 	= new AssocLegendrePoly_sun*[coeffOrder + 1];

	for(uint i=0;i<=coeffOrder;++i)
	{
		g[i] 			= new hcFloat[i+1];
		h[i] 			= new hcFloat[i+1];
		assLegPoly[i] 	= new AssocLegendrePoly_sun[i+1];

		for(uint j=0;j<=i;++j)
		{
			AssocLegendrePoly_sun temp(i, j);
			assLegPoly[i][j] = temp;
		}
	}
}

void PFSSsolution_SHC_hoek::eval(const Vec3D &pos, Vec3D &result)
{
	hcFloat r			= pos[0];
	hcFloat theta		= pos[1];
	hcFloat phi			= pos[2];

	hcFloat r_ss 		= sourceSurfaceFactor * r_sol;     // source surface
	hcFloat retval_r 	= 0.0;
	hcFloat retval_t 	= 0.0;
	hcFloat retval_p 	= 0.0;

	for(uint l=0;l<=coeffOrder;++l)
		for(uint m=0;m<=l;++m)
		{
			hcFloat c	= -pow(r_sol / r_ss, (hcFloat)(l+2));

			retval_r += assLegPoly[l][m](cos(theta))
					  * ( g[l][m] * cos(m*phi) + h[l][m] * sin(m*phi) )
					  * ((l+1) * pow(r_sol / r, (hcFloat)(l+2)) - c * l * pow(r / r_ss, (hcFloat)((int)l-1)));

			retval_t += -assLegPoly[l][m].getDeriv(theta)
					  * ( g[l][m] * cos(m*phi) + h[l][m] * sin(m*phi) )
					  * (pow(r_sol / r, (hcFloat)(l+2)) + c * pow(r / r_ss, (hcFloat)((int)l-1)));

			retval_p += assLegPoly[l][m](cos(theta))
					  * ( g[l][m] * sin(m*phi) - h[l][m] * cos(m*phi) )
					  * (pow(r_sol / r, (hcFloat)(l+2)) + c * pow(r / r_ss, (hcFloat)((int)l-1))) * m / sin(theta);
		}

	result(0) = retval_r;
	result(1) = retval_t;
	result(2) = retval_p;
}

//---------------------------------------------------------------------------------------------------------------
//---   CSSS
//---------------------------------------------------------------------------------------------------------------



CSSS_magfield::CSSS_magfield(const char *filename)
{
	initNULL();
    if(!load(filename))
    	printf("ERROR! CSSS_magfield constructor could not load file '%s'\n", filename);
}

CSSS_magfield::CSSS_magfield()
{
	initNULL();
    printf("ERROR! CSSS_magfield std constructor not implemented!\n");
}

CSSS_magfield::CSSS_magfield(const CSSS_magfield &other)
{
	initNULL();
	printf("ERROR! CSSS_magfield cpy constructor not implemented!\n");
}

CSSS_magfield::~CSSS_magfield()
{
    clear();
}

void CSSS_magfield::initNULL()
{
	g 					= NULL;
	h 					= NULL;
	assLegPoly 			= NULL;
	sourceSurfaceFactor = 2.5;
	order 			= 0;
}

void CSSS_magfield::clear()
{
	if(g != NULL)
	{
		for(uint i=0;i<=order;++i)
		{
			delete [] g[i];
			delete [] h[i];
			delete [] assLegPoly[i];
		}
		delete [] g;
		delete [] h;
		delete assLegPoly;
	}

	initNULL();
}


void CSSS_magfield::init(uint order)
{
	clear();

    this->order	= order;
    g 			= new double*[order + 1];
    h 			= new double*[order + 1];
    assLegPoly 	= new AssocLegendrePoly_sun*[order + 1];

    for(uint l=0;l<=order;++l)
    {
        g[l] 			= new double[l+1];
        h[l] 			= new double[l+1];
        assLegPoly[l] 	= new AssocLegendrePoly_sun[l+1];
        for(uint m=0;m<=l;++m)
        {
            AssocLegendrePoly_sun temp(l, m);
            assLegPoly[l][m] 	= temp;
            g[l][m]				= 0.0;
            h[l][m]				= 0.0;
        }
    }
}

void CSSS_magfield::exportCoefficients(const char *filename)
{
    char tempstr[10000];
    char tmp[100];
    FILE *outfile = fopen(filename, "w");
    for(uint l=0;l<=order;++l)
    {
        tempstr[0] = '\0';
        sprintf(tempstr, "%u ", l);
        for(uint m=0;m<=l;++m)
        {
            tmp[0] = '\0';
            sprintf(tmp, "%f ", g[l][m]);
            strcat(tempstr, tmp);
        }
        if(l==0)
            strcat(tempstr, "\t\t g's");
        strcat(tempstr, "\n");
        fprintf(outfile, "%s", tempstr);
    }
    fprintf(outfile, "\n");
    for(uint l=0;l<=order;++l)
    {
        tempstr[0] = '\0';
        sprintf(tempstr, "%u ", l);
        for(uint m=0;m<=l;++m)
        {
            tmp[0] = '\0';
            sprintf(tmp, "%f ", h[l][m]);
            strcat(tempstr, tmp);
        }
        if(l==0)
            strcat(tempstr, "\t\t h's");
        strcat(tempstr, "\n");
        fprintf(outfile, "%s", tempstr);
    }
    fclose(outfile);
}

/*
void CSSS_magfield::determineCoefficientsFromFITSsineLat(uint order, const char *FITS_filename,
		double maxSinLat, double minSinLat, bool accountForMissingFlux)
{

    uint x, y, l, m;
    initArrays(order);
    hcImageFITS fitsIMG(FITS_filename);
    uint width 	= fitsIMG.width;
    uint height = fitsIMG.height;

    double dTheta = PI / (2 * height);
    double dSinTheta = (maxSinLat - minSinLat) / (height-1);
    double dPhi = 2*PI/(width);

    // compute missing-flux-term
    double delta = 0.0;
    double tempB;
    double theta;
    double cosecans;
    double numerator = 0.0;

    for(y=0;y<height;++y)
    {
        tempB = 0.0;
        // this is for equally spaced theta
        //theta = (2.0*(height - 1 - y) + 1 ) * dTheta;
        // this is for equally spaced sin(latitude) -- theta = 180° - (lat + 90°)
        theta = PI / 2 - asin(maxSinLat - (height - 1- y) * dSinTheta);

        cosecans = 1.0 / sin(theta);
        numerator += cosecans;


        for(x=0;x<width;++x)
            if(!isnan(fitsIMG(x,y)) && !(fitsIMG(x,y) == -32768.0))
            {
                tempB += fitsIMG(x,y);
            }
        delta += tempB * cosecans;
    }

    delta /= (numerator * width);
    if(accountForMissingFlux == false)
        delta = 0;

    // now compute the coefficients
    numerator = 1.0 / (width * height);
    double legPolyEval;
    double phi;
    double D_xy;

    for(l=0;l<=coeffOrder;++l)
        for(m=0;m<=l;++m)
        {
            g[l][m] = 0.0;
            h[l][m] = 0.0;
            for(y=0;y<height;++y)
            {
                // this is for equally spaced theta
                //theta = (2*(height - 1 - y) + 1 ) * dTheta;
                //this is for equally spaced sin(latitude)
                theta 		= PI / 2 - asin(maxSinLat - (height - 1- y) * dSinTheta);

                legPolyEval = assLegPoly[l][m].operator ()(cos(theta));
                cosecans 	= 1 / sin(theta);

                for(x=0;x<width;++x)
                {
                    phi = m * x * dPhi;
                    if(!isnan(fitsIMG(x,y)) && !(fitsIMG(x,y) == -32768.0))
                        D_xy = fitsIMG(x,y) - delta;
                    else
                        D_xy = -delta;
                    g[l][m] += D_xy * cos(phi) * cosecans * legPolyEval;
                    h[l][m] += D_xy * sin(phi) * cosecans * legPolyEval;
                }
            }

            //according to Altschuler & Newkirk for source surf at inf
            //g[l][m] *= (2.0*l + 1) / (l + 1) * numerator;
            //h[l][m] *= (2.0*l + 1) / (l + 1) * numerator;

            //according to Xudong Sun:
            g[l][m] *= (2.0*l + 1) * numerator * 100;
            h[l][m] *= (2.0*l + 1) * numerator * 100;
        }
}//*/

bool CSSS_magfield::load(const char *filename)
{
	if(!doesFileExist(filename))
		return false;

	char line[1000];
	ifstream file(filename);

	uint order = 0;

	while (!file.eof())
	{
		file.getline(line, 1000);

		if((line[0] == '#') || line[0] == '\n')
			continue;

		regex re(".*([0-9]{3}).*", regex::icase);
		cmatch what;
		regex_search(line, what, re);

		if(what[0].matched)
		{
			order = boost::lexical_cast<int>(what[1]);
			break;
		}

		re.assign(".*([0-9]{2}).*", regex::icase);
		regex_search(line, what, re);

		if(what[0].matched)
		{
			order = boost::lexical_cast<int>(what[1]);
			break;
		}

		re.assign(".*([0-9]{1}).*", regex::icase);
		regex_search(line, what, re);

		if(what[0].matched)
		{
			order = boost::lexical_cast<int>(what[1]);
			break;
		}

		printf("ERROR! loadCSSScoeff: Order of solution could not be found!\n");
		return false;
	}

	init(order);

	char lstr[1000], mstr[1000], gstr[1000], hstr[1000];

	while(file >> lstr >> mstr >> gstr >> hstr)
	{
		uint l 		= atoi(lstr);
		uint m 		= atoi(mstr);
		double gval = atof(gstr);
		double hval	= atof(hstr);

		if((l>order) || (m>l))
		{
			printf("Error! loadCSSScoeff: l: %u, m: %u\n", l, m);
			return false;
		}

		g[l][m]	= gval;
		h[l][m]	= hval;
	}
	return true;
}

void CSSS_magfield::eval(const Vec3D &pos, Vec3D &result)
{
	hcFloat r			= pos[0];
	hcFloat theta		= pos[1];
	hcFloat phi			= pos[2];

    double r_ss 		= sourceSurfaceFactor * r_sol;     // source surface
    double retval_r 	= 0.0;
    double retval_t 	= 0.0;
    double retval_phi 	= 0.0;

    for(uint l=0;l<=order;++l)
        for(uint m=0;m<=l;++m)
        {   // equation 16 from Notes on PFSS Extrapolation, Xudong Sun
            retval_r += assLegPoly[l][m](cos(theta))
                    	* ( g[l][m] * cos(m*phi) + h[l][m] * sin(m*phi) )
                    	* pow(r_sol / r, (double)(l+2))
                    	* (l+1+l*pow(r / r_ss, (double)(2*l+1))) / (l+1+l*pow(r_sol / r_ss, (double)(2*l+1)));

            // equation 17 (there has to be theta instead of cos(theta) below!)
            retval_t += -assLegPoly[l][m].getDeriv(theta)
                    	* ( g[l][m] * cos(m*phi) + h[l][m] * sin(m*phi) )
                    	* pow(r_sol / r, (double)(l+2))
                    	* (1 - pow(r / r_ss, (double)(2*l+1))) / (l+1+l*pow(r_sol / r_ss, (double)(2*l+1)));

            // equation 18
            retval_phi += assLegPoly[l][m](cos(theta))
                    	  * ( g[l][m] * sin(m*phi) - h[l][m] * cos(m*phi) )
                    	  * pow(r_sol / r, (double)(l+2))
                    	  * (1 - pow(r / r_ss, (double)(2*l+1))) / (l+1+l*pow(r_sol / r_ss, (double)(2*l+1)));//*/

        }
    result(0) = retval_r;
    result(1) = retval_t;
    result(2) = retval_phi;
}

void CSSS_magfield::dump()
{
	printf("CSSS harmonic coefficients (up to order %u):\n----------------------------------------------------------\n", order);
	for(uint l=1; l<=order; ++l)
		for(uint m=0; m<=l; ++m)
			printf("l:\t%u, m:\t%u, g_lm:\t%E, h_lm:\t%E\n", m, l, g[l][m], h[l][m]);

}

//---------------------------------------------------------------------------------------------------------------
//---   PFSS_visualization
//---------------------------------------------------------------------------------------------------------------

#ifdef GUI

PFSS_visualization::PFSS_visualization(){

	drawScale			= 0.0;
	upperOrder 			= 0;
	lowerOrder			= 0;
	upperSS				= 0;
	lowerSS				= 0;
	textScreenID		= 0;
	ssAnimScreenID		= 0;
	contourProjScreenID	= 0;
	accumulatedTimediff = 0;
	animPointer			= 0;
	animationPlaying	= false;
	fps					= 0;
	numFrames			= 0;
	magLines			= NULL;
}

PFSS_visualization::PFSS_visualization(char *filename) : hcScene3D(){

    initFromCoeffFile(filename);
}

PFSS_visualization::~PFSS_visualization()
{
    clear();
}

void PFSS_visualization::clear()
{
    if(magLines != NULL)
    {
        delete [] magLines;
        magLines = NULL;
    }
}

void PFSS_visualization::init()
{
    clear();
    hcScene3D::init();
    magLines 				= NULL;
    accumulatedTimediff 	= 0.0;
    drawScale 				= 1 / (1.5 * r_sol);
    double ssSize 			= magfield.sourceSurfaceFactor * r_sol * drawScale;
    double sunSize			 = r_sol * drawScale;
    sun.appendFragment(*(new hcDSphere3D(origin3D, 1.0)));
    sources.appendFragment(*(new hcDSphere3D(origin3D, 1.0)));
    sources.TFmodel.loadIdentity();
    sources.TFmodel.scaleHom(ssSize);
    sun.TFmodel.loadIdentity();
    sun.TFmodel.scaleHom(sunSize);
    animationPlaying = false;
}

void PFSS_visualization::initFromSynopticFITSsineLat(uint order, const char *filename, double maxSinLat, double minSinLat)
{
    init();
    SynPhotMagfield synPhot;
    synPhot.load(filename);
    magfield.determineCoefficientsFromPhotMagfield(order, 2.5*r_sol, synPhot);

    sun.texContainer.rmTexture();
    sun.texContainer.addFITSTexFromFile(filename);

    sources.texContainer.rmTexture();
    sources.texContainer.addMonoTexture(char2RGBA8(255,255,255,50));

    initScene();

}

void PFSS_visualization::initFromSynopticWSO(uint order, const char *filename)
{
    std::ifstream in(filename);

    char instrings[31][1000];
    in.getline(instrings[0], 1000);
    in.getline(instrings[0], 1000);

    uint ntheta 		= 30;
    uint nphi 			= 73;
    //hcFloat *imageData 	= new hcFloat[ntheta * nphi];
    hcImageFITS img(nphi, ntheta);

    uint i;
    uint j=0;
    while(in >> instrings[0] >> instrings[1] >> instrings[2] >> instrings[3] >> instrings[4] >> instrings[5] >> instrings[6] >> instrings[7] >> instrings[8] >> instrings[9] >> instrings[10] >> instrings[11] >> instrings[12] >> instrings[13] >> instrings[14] >> instrings[15] >> instrings[16] >> instrings[17] >> instrings[18] >> instrings[19] >> instrings[20] >> instrings[21] >> instrings[22] >> instrings[23] >> instrings[24] >> instrings[25] >> instrings[26] >> instrings[27] >> instrings[28] >> instrings[29] >> instrings[30])
    {
        if(strncmp(instrings[0], "CTxxxxxxxx", 2))
        {
            printf("ERROR! PFSS_visualization::initFromSynopticWSO: Data alignment mismatch! Check the WSO-file!\n");
            printf("j: %u, instrings[0]: %s\n", j, instrings[0]);
            exit(1);
        }

        for(i=0;i<ntheta;++i)
        {
            img(nphi - j - 1, ntheta - i - 1) = atof(instrings[i+1]);
        }
        j++;
    }
    char outfilename[1000];
    outfilename[0] = '\0';
    strcat(outfilename, filename);
    strcat(outfilename, ".FITS\0");

    img.save(outfilename);
    //delete [] imageData;
    initFromSynopticFITSsineLat(order, outfilename, 14.5/15, -14.5/15);
}

void PFSS_visualization::initFromCoeffFile(const char *filename)
{
    init();
    magfield.importCoefficients(filename, true, 2.5*r_sol);

    sun.texContainer.rmTexture();
    sun.texContainer.addTexFromFile("../textures/spicules_sst_big.jpg");

    sources.texContainer.rmTexture();
    sources.texContainer.addMonoTexture(char2RGBA8(255,255,255,50));

    initScene();
}

//*/

/*! \brief compares two coefficient files
 *
 *  treats the second file as an deviation from the first and prints the
 *  relative errors for each coefficient in outfile
 */
void PFSS_visualization::compareCoeffFiles(const char *infile1, const char *infile2, const char *outfile)
{
    FILE *coeffFile1, *coeffFile2, *compFile;
    char str[1000];
    char tempstr[100];
    uint i,j,k,z,order;
    int highestOrder1 	= 0;
    int highestOrder2 	= 0;
    tempstr[0] 			= '\0';
    str[0] 				= '\0';

    coeffFile1 = fopen(infile1, "r");
    if(!coeffFile1)
    {
        printf("ERROR! PFSSsolution_SHC_sun::compareCoeffFiles: File %s does not exist!\n", infile1);
        return;
    }
    i=0;
    while(strcmp(fgets(str, 1000, coeffFile1), "\n"))
    {
        while(str[i++] != ' ');
        strncpy(tempstr, str, --i);
        tempstr[i] = '\0';
        if (atoi(tempstr) > highestOrder1)
            highestOrder1 = atoi(tempstr);
    }

    coeffFile2 = fopen(infile2, "r");
    if(!coeffFile2)
    {
        printf("ERROR! PFSSsolution_SHC_sun::compareCoeffFiles: File %s does not exist!\n", infile1);
        return;
    }
    i=0;
    while(strcmp(fgets(str, 1000, coeffFile2), "\n"))
    {
        while(str[i++] != ' ');
        strncpy(tempstr, str, --i);
        tempstr[i] = '\0';
        if (atoi(tempstr) > highestOrder2)
            highestOrder2 = atoi(tempstr);
    }

    if(highestOrder1 != highestOrder2)
    {
        printf("ERROR! PFSSsolution_SHC_sun::compareCoeffFiles: Coefficient order does not match! (%i != %i)\n", highestOrder1, highestOrder2);
        return;
    }

    fseek(coeffFile1,0, SEEK_SET);
    fseek(coeffFile2,0, SEEK_SET);

    double **g1 = (double**)malloc(sizeof(double) * (highestOrder1+1));
    double **h1 = (double**)malloc(sizeof(double) * (highestOrder1+1));
    double **g2 = (double**)malloc(sizeof(double) * (highestOrder1+1));
    double **h2 = (double**)malloc(sizeof(double) * (highestOrder1+1));
    for(i=0;i<=highestOrder1;++i)
    {
        g1[i] = (double*) malloc(sizeof(double) * (i+1));
        h1[i] = (double*) malloc(sizeof(double) * (i+1));
        g2[i] = (double*) malloc(sizeof(double) * (i+1));
        h2[i] = (double*) malloc(sizeof(double) * (i+1));
    }

    j=0;
    while (fgets(str, 1000, coeffFile1) != NULL)
    {
        if(j==0)        // load g's
        {
            if(str[0] == '\n')
                ++j;
            else
            {
                i = 0;
                tempstr[0] = '\0';
                while(str[i++] != ' ');
                z = i;
                strncpy(tempstr, str, i-1);
                tempstr[i-1] = '\0';
                order=(uint)atoi(tempstr);
                k=0;

                while(k<=order)
                {
                    while(str[i] == ' ' || str[i] == '\t')
                        ++i;
                    z=i;
                    while(str[i] != ' ' && str[i++] != '\n');
                    strncpy(tempstr, str+z, i-1);
                    tempstr[i-1-z] = '\0';
                    g1[order][k++] = strtof(tempstr, NULL);
                    tempstr[0] = '\0';
                    z=i;
                }
                ++order;
            }
        }
        else            // load h's
        {
            if(str[0] == '\n')
                break;
            else
            {
                i = 0;

                tempstr[0] = '\0';
                while(str[i++] != ' ');
                z = i;
                strncpy(tempstr, str, i-1);
                tempstr[i-1] = '\0';
                order=(uint)atoi(tempstr);
                k=0;

                while(k<=order)
                {
                    while(str[i] == ' ' || str[i] == '\t')
                        ++i;
                    z=i;
                    while(str[i] != ' ' && str[i++] != '\n');
                    strncpy(tempstr, str+z, i-1);
                    tempstr[i-1-z] = '\0';
                    h1[order][k++] = strtof(tempstr, NULL);
                    tempstr[0] = '\0';
                    z=i;
                }
                ++order;
            }
        }
    }

    j=0;
    while (fgets(str, 1000, coeffFile2) != NULL)
    {
        if(j==0)        // load g's
        {
            if(str[0] == '\n')
                ++j;
            else
            {
                i = 0;

                tempstr[0] = '\0';
                while(str[i++] != ' ');
                z = i;
                strncpy(tempstr, str, i-1);
                tempstr[i-1] = '\0';
                order=(uint)atoi(tempstr);
                k=0;

                while(k<=order)
                {
                    while(str[i] == ' ' || str[i] == '\t')
                        ++i;
                    z=i;
                    while(str[i] != ' ' && str[i++] != '\n');
                    strncpy(tempstr, str+z, i-1);
                    tempstr[i-1-z] = '\0';
                    g2[order][k++] = strtof(tempstr, NULL);
                    tempstr[0] = '\0';
                    z=i;
                }
                ++order;
            }
        }
        else            // load h's
        {
            if(str[0] == '\n')
                break;
            else
            {
                i = 0;

                tempstr[0] = '\0';
                while(str[i++] != ' ');
                z = i;
                strncpy(tempstr, str, i-1);
                tempstr[i-1] = '\0';
                order=(uint)atoi(tempstr);
                k=0;

                while(k<=order)
                {
                    while(str[i] == ' ' || str[i] == '\t')
                        ++i;
                    z=i;
                    while(str[i] != ' ' && str[i++] != '\n');
                    strncpy(tempstr, str+z, i-1);
                    tempstr[i-1-z] = '\0';
                    h2[order][k++] = strtof(tempstr, NULL);
                    tempstr[0] = '\0';
                    z=i;
                }
                ++order;
            }
        }
    }

    fclose(coeffFile1);
    fclose(coeffFile2);
    compFile = fopen(outfile, "w");

    double error;
    double gAccum = 0.0;
    double hAccum = 0.0;
    for(i=0;i<=highestOrder1;++i)
    {
        fprintf(compFile, "%i ", i);
        for(j=0;j<=i;++j)
        {
            error = fabs((g2[i][j]-g1[i][j])/g1[i][j]);
            fprintf(compFile, "%f\t", error);
            if(!isnan((long double)error))
                gAccum += error;
        }
        fprintf(compFile, "\n");
    }

    fprintf(compFile, "\n");
    for(i=0;i<=highestOrder1;++i)
    {
        fprintf(compFile, "%i ", i);
        for(j=0;j<=i;++j)
        {
            error = fabs((h2[i][j]-h1[i][j])/h1[i][j]);
            fprintf(compFile, "%f\t", error);
            if(!isnan((long double)error))
                hAccum += error;
        }
        fprintf(compFile, "\n");
    }

    fprintf(compFile, "\nAccumulated error in g: %f\nAccumulated error in h: %f\n", gAccum, hAccum);
    fclose(compFile);

    for(i=0;i<=highestOrder1;++i)
    {
        free(g1[i]);
        free(h1[i]);
        free(g2[i]);
        free(h2[i]);
    }
    free(g1);
    free(h1);
    free(g2);
    free(h2);
}

void PFSS_visualization::initScene()
{
#define maxNumPointsPerMagLine 1000

    uint numTheta 	= 10;
    uint numPhi 	= 20;

    uint i,j,k;
    Vec3D temp;
    Vec3D pos;
    Vec3D temp_mag;
    Vec3D *posData;
    Vec3D gridPos;
    double theta, phi;

    appendObject(sun);

    magLines = new hcDObject3D[1];
    magLines[0].init();
    magLines[0].texContainer.init(3);
    magLines[0].texContainer.addMonoTexture(char2RGBA8(255, 0, 0, 255)); //downward maglines
    magLines[0].texContainer.addMonoTexture(char2RGBA8(0, 255, 0, 255)); //upward maglines
    magLines[0].texContainer.addMonoTexture(char2RGBA8(255, 255, 255, 30)); //source surface

    hcDLineStrip3D lineStrip;

    for(j=0;j<=numTheta;++j)
        for(k=0;k<numPhi;++k)
        {
            if((j==0 || j==numTheta) && k > 0) 			// compute poles only for k=0
                break;
            bool upwards;
            theta 				= j * PI / numTheta;	// TODO please review this
            phi 				= k * 2 * PI / numPhi;	// TODO please review this
            gridPos.content[0] 	= r_sol;
            gridPos.content[1] 	= theta;
            gridPos.content[2] 	= phi;

            posData 			= new Vec3D[maxNumPointsPerMagLine];

            for(i=0;i<maxNumPointsPerMagLine;++i)
            {
                if(i==0)
                {
                    magfield.eval(gridPos, temp_mag);

                    if(temp_mag.content[0] < 0)
                        upwards = false;
                    else
                        upwards = true;

                    pos = gridPos.convCoordSpher2Cart();

                    pos.scale(drawScale);
                    temp = temp_mag;
                    //temp_mag.convertSphericalToCartesian(temp);
                    temp_mag = temp_mag.convVecSpher2Cart(gridPos.convCoordSpher2Cart());
                }
                else
                {
                    if((pos.length() / drawScale < r_sol - 1E7) || (pos.length() / drawScale > magfield.sourceSurfaceFactor * r_sol))
                    {
                        lineStrip.init(i, posData);
                        lineStrip.drawOptions.computeLighting = false;

                        uint num = magLines[0].appendFragment(lineStrip);
                        if(!upwards)
                            magLines[0][num].texturePointer = 0;
                        else
                            magLines[0][num].texturePointer = 1;
                        break;
                    }

                    temp = pos;
                    temp.scale(1/drawScale);
                    temp = temp.convCoordCart2Spher();
                    magfield.eval(temp, temp_mag);
                }

                if(!upwards)
                {
                    //temp.loadZeroes();
                    //temp.minus(temp_mag);
                    temp = -temp_mag;
                    temp_mag = temp;
                }

                posData[i] = pos;
                temp_mag.scale(0.01 / temp_mag.length());
                //pos.plus(temp_mag);
                pos	+= temp_mag;
            }
            delete [] posData;
        }
    appendObject(magLines[0]);
    appendObject(sources);
}

int PFSS_visualization::loadAnimation(const char *path)
{
    uint i, j;
    numFrames 	= 1;
    lowerSS 	= 0;
    upperSS 	= 0;
    lowerOrder 	= 0;
    upperOrder 	= 0;

    char cfgFilename[1000];
    char inputFilename[1000];
    cfgFilename[0] = '\0';
    sprintf(cfgFilename, "%s%s", path, "config.cfg");
    std::ifstream config(cfgFilename);

    char option[1000], argument[1000];
    option[0] = '\0';
    argument[0] = '\0';

    while(config >> option >> argument)
    {
        if (!strcmp(option, "numFrames"))
            numFrames = atof(argument);
        if (!strcmp(option, "lowerSS"))
            lowerSS = atof(argument);
        if (!strcmp(option, "upperSS"))
            upperSS = atof(argument);
        if (!strcmp(option, "fps"))
            fps = atof(argument);
        if (!strcmp(option, "lowerOrder"))
            lowerOrder = (uint)atoi(argument);
        if (!strcmp(option, "upperOrder"))
            upperOrder = (uint)atoi(argument);
    }

    clear();
    hcScene3D::init();
    //appendElement(sun);
    magLines = new hcDObject3D[numFrames];

    for(i=0;i<numFrames;++i)
    {
        inputFilename[0] = '\0';
        sprintf(inputFilename, "%smaglines%u.obj", path, i);

        magLines[i].loadMesh(inputFilename);
        magLines[i].visible = false;
        appendObject(magLines[i]);

    }

    overlay.init();
    Vec2D ssMagScreenPos(0.38, 0.8);
    ssAnimScreenID = overlay.createAnimationScreen(ssMagScreenPos, 1.2, 0.3, numFrames);
    for(i=0;i<numFrames;++i)
    {
        inputFilename[0] = '\0';
        sprintf(inputFilename, "%sssmagfield%u.fits.contour.BMP", path, i);
        overlay[ssAnimScreenID].texContainer.addTexFromFile(inputFilename);
    }

    Fonthandler ftHandler;
    char fileNam[1000] = "../data/Roboto-Regular.ttf";
    ftHandler.init(fileNam, 20);

    uint reqWidth, reqHeight;
    float reqRelWidth, reqRelHeight;
    //uint *textlayer;


    float relHeight = 0.1;
    float relWidth = 0.5;
    uint pixelHeight = relHeight * heightPixelScale;
    uint pixelWidth = relWidth * widthPixelScale;

    Vec2D textScreenPos(-0.65, 0.9);
    textScreenID = overlay.createAnimationScreen(textScreenPos, relWidth, relHeight, numFrames);

    //textlayer = (uint*)malloc(sizeof(uint) * pixelHeight * pixelWidth);
    char text[100];

    hcImageRGBA textlayer(pixelWidth, pixelHeight);

    for(i=0;i<numFrames;++i)
    {
        text[0] = '\0';
        if(upperSS == 0)    // animation is concerned with order of spherical functions
            sprintf(text, "order: %u", lowerOrder+i);
        else
            sprintf(text, "R = %f", lowerSS + i * (float)(upperSS-lowerSS) / numFrames);


        j = 0;
        while(text[j++] != '\0');
        --j;

        ftHandler.getImgDim(text, 0, 0, reqWidth, reqHeight);
        reqRelHeight = reqHeight / heightPixelScale;
        reqRelWidth = reqWidth / heightPixelScale;
        if(reqRelHeight > relHeight)
            printf("ERROR: PFSS_visualization::loadAnimation(const char *path): Resolution of window too low for text rendering\n(height = %E, required: %E)!\n", relHeight, reqRelHeight);
        if(reqRelWidth > relWidth)
            printf("ERROR: PFSS_visualization::loadAnimation(const char *path): Resolution of window too low for text rendering\n(width = %e, required: %E)!\n", relWidth, reqRelWidth);

        ftHandler.printTextToImage(text, char2RGBA8(255, 255, 255, 255), 0, 0, pixelWidth, pixelHeight, textlayer);
        overlay[textScreenID].texContainer.addTexture2D(textlayer, false);

    }//*/
    //free(textlayer);

    Vec2D contourPosProj(0.38, 0.5);
    contourProjScreenID = overlay.createAnimationScreen(contourPosProj, 1.2, 0.3, numFrames);
    for(i=0;i<numFrames;++i)
    {
        inputFilename[0] = '\0';
        sprintf(inputFilename, "%smagfield25%u.fits.contour.BMP", path, i);
        overlay[contourProjScreenID].texContainer.addTexFromFile(inputFilename);
    }

    return 1;
}

#define outputwidth 800
#define outputheight 200

void PFSS_visualization::exportRadialMagfieldAtHeight(float r, const char *filename)
{
	hcFloat dPhi 		= 2 * PI / outputwidth;
    hcFloat dTheta 		= PI / (outputheight-1);
    //hcFloat *imageData 	= new hcFloat[outputwidth * outputheight];

    hcImageFITS img(outputwidth, outputheight);

    Vec3D pos, mag;
    pos.content[0] = r;

    for(uint x=0;x<outputwidth;++x)
        for(uint y=0;y<outputheight;++y)
        {
            pos.content[1] = PI - dTheta * y;	// TODO review this
            pos.content[2] = dPhi * x;			// TODO review this
            magfield.eval(pos, mag);
            //imageData[y * outputwidth + x] = mag[0];
            img(x, y) = mag[0];
        }

    img.save(filename);

    //delete [] imageData;
}

void PFSS_visualization::exportRadialProjectedMagfieldAtHeight(float r, const char *filename)
{
    uint x,y;
    hcFloat dPhi 		= 2 * PI / outputwidth;
    hcFloat dTheta 		= PI / (outputheight-1);
    //hcFloat *imageData 	= new hcFloat[outputwidth * outputheight];
    hcImageFITS img(outputwidth, outputheight);

    Vec3D pos, mag;
    pos.content[0] = magfield.sourceSurfaceFactor * r_sol;
    hcFloat distance = r - pos[0];

    for(x=0;x<outputwidth;++x)
        for(y=0;y<outputheight;++y)
        {
            pos.content[1] = PI - dTheta * y;	// TODO review this
            pos.content[2] = dPhi * x;			// TODO review this
            magfield.eval(pos, mag);
            img(x,y) = mag[0] * pow(pos[0]/r, 2);
        }

    img.save(filename);

    //delete [] imageData;
}

void PFSS_visualization::draw(GLuint programID)
{
    hcScene3D::draw(programID);
    overlay.draw();
}
#endif









