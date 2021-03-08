#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "engine/math/hcVec.h"

#include "math.h"

double max(double a, double b);

double min(double a, double b);

uint max(uint a, uint b);

uint min(uint a, uint b);

double distOnSphere(const Vec3D &pos0, const Vec3D &pos1);
    /*!< \brief computes the length of the great arc connecting pos0 and pos1               */

double distOnUnitSphere(const Vec3D &pos0, const Vec3D &pos1);
    /*!< \brief computes the (angular) distance between two points on the unit sphere       */

void getPosOnGreatArc(const Vec3D &pos0, const Vec3D pos1, float relPos, Vec3D &result);
    /*!< \brief computes the position on the great arc between pos0 and pos1                */

double haversin(double theta);
    /*!< \brief haversine function                                                          */

double ahaversin(double h);
    /*!< \brief inverse haversine function                                                  */

double factorial(double n);
    /*!< \brief n!                                                                          */

hcFloat binomial(hcFloat n, uint k);
    /*!< \brief computes n over k                                                           */

//bool loadCSSScoeff(double ***g, double ***h, const char *filename);
	/*!< \brief loads coefficients for Spherical functions computed by Bala via CSSS		*/

/*!
 *  \brief polynomial of order n
 */
class Polynomial{
public:

    double *factor;
    uint n;

    Polynomial(uint order=0);				/*!< \brief std constructor					*/
    Polynomial(const Polynomial &other);	/*!< \brief cpy constructor					*/
    Polynomial(double *factor, uint n);
    ~Polynomial();							/*!< \brief destructor						*/

    Polynomial &operator=(const Polynomial &other);		/*!< \brief assignment operator	*/

    void initNULL();
    void clear();
    void init(uint n);

    double operator()(double x);
    void derivative(uint m);
    void scale(double scale);
    void dump();
};

class LegendrePoly : public Polynomial{
public:

    LegendrePoly(uint order = 0);					/*!< \brief std constructor			*/
    LegendrePoly(const LegendrePoly &other);		/*!< \brief cpy constructor			*/
    ~LegendrePoly();								/*!< \brief destructor				*/

    LegendrePoly &operator=(const LegendrePoly &other);	/*!< \brief assignment operator	*/

    void initNULL();
    void clear();
    void init(uint n);
};


/*! \brief conventional associated Legendre polynomials (as described by, e.g., Wikipedia)
 */
class AssocLegendrePoly{
public:

	uint l;
    uint m;
    LegendrePoly dmPl;		/*!< \brief m'th derivative of the Legendre polynomial											*/
    LegendrePoly dm_1Pl;   	/*!< \brief the (m+1)st derivative of the Legendre polynomial (necessary for theta component) 	*/

    AssocLegendrePoly(uint l=0, uint m=0);				/*!< \brief std constructor											*/
    AssocLegendrePoly(const AssocLegendrePoly &other);	/*!< \brief cpy constructor											*/
    virtual ~AssocLegendrePoly();						/*!< \brief destructor												*/

    AssocLegendrePoly &operator=(const AssocLegendrePoly &other);	/*!< \brief assignment op								*/
    virtual hcFloat operator()(hcFloat x);							/*!< \brief evaluator									*/

    void initNULL();
    void clear();
    void init(uint l, uint m);

    virtual hcFloat getDeriv(hcFloat theta);	/*!< \brief returns first derivative at cos(theta)							*/

    void dump();
};

/*! \brief associated Legendre polynomial as used by Xudong Sun for WSO data
 */
class AssocLegendrePoly_sun : public AssocLegendrePoly{
public:

    AssocLegendrePoly_sun(uint l=0, uint m=0);						/*!< \brief std constructor								*/
    AssocLegendrePoly_sun(const AssocLegendrePoly_sun &other);		/*!< \brief cpy constructor								*/
    virtual ~AssocLegendrePoly_sun();								/*!< \brief destructor									*/

    AssocLegendrePoly_sun &operator=(const AssocLegendrePoly_sun &other);	/*!< \brief assignment operator					*/
    virtual hcFloat operator()(hcFloat x);

    virtual hcFloat getDeriv(hcFloat theta);
};


#endif // FUNCTIONS_H
