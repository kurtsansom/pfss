#include "engine/math/hcFunction.h"
#include "engine/hcTools.h"

#include "math.h"
#include <stdlib.h>
#include "stdio.h"
#include <string.h>
#include <fstream>

#include "boost/lexical_cast.hpp"




long double max(long double a, long double b)
{
    return a > b ? a : b;
}

long double min(long double a, long double b)
{
    return a < b ? a : b;
}

double max(double a, double b)
{
    return a > b ? a : b;
}

double min(double a, double b)
{
    return a < b ? a : b;
}

uint max(uint a, uint b)
{
    return a > b ? a : b;
}

uint min(uint a, uint b)
{
    return a < b ? a : b;
}

/*! vec1, vec2 given in 3D spherical coordinates.
 *  The first components have to match in order to be on the same sphere, otherwise
 *  this function does not make sense
 */
double distOnSphere(const Vec3D &pos0, const Vec3D &pos1)
{
    float eps = 1E-2; // TODO: make more accurate

    if(fabs(pos0[0] - pos1[0]) / pos0[0] > eps)
    {
        printf("ERROR! distOnSphere: 3D points do not lie on the same spherical surface!\nr1: %E, r2: %E\n", pos0[0], pos1[0]);
        printf("pos0:\n");
        pos0.dump();
        printf("pos1:\n");
        pos1.dump();
        printf("fabs(pos0[0] - pos1[0]) / pos0[0]: %E\n", fabs(pos0[0] - pos1[0]) / pos0[0]);
        exit(1);
        return 0.0;
    }

    double theta1   = pos0[1];
    double theta2   = pos1[1];
    double phi1     = pos0[2];
    double phi2     = pos1[2];

    //haversin formula
    double h 		= haversin(theta1-theta2) + sin(theta1) * sin(theta2) * haversin(phi2-phi1);
    double result 	= pos0[0] * ahaversin(h);

    return result;
}

/*! the same as distOnSphere but the first components of the vectors are not checked
 *  as all points are assumed to lie on the unit sphere
 */
double distOnUnitSphere(const Vec3D &pos0, const Vec3D &pos1)
{
    double theta1   = pos0[1];
    double theta2   = pos1[1];
    double phi1     = pos0[2];
    double phi2     = pos1[2];

    //haversin formula
    double h 		= haversin(theta1-theta2) + sin(theta1) * sin(theta2) * haversin(phi2-phi1);
    double result 	= ahaversin(h);

    return result;
}

/*! relPos is between 0 and 1, where 0 points to pos0, 1 to pos1 and 0.5 to the midpoint of
 *  the great arc connecting pos0 and pos1. All vectors (pos0, pos1 and result) are given
 *  in spherical coordinates
 */
void getPosOnGreatArc(const Vec3D &pos0, const Vec3D pos1, float relPos, Vec3D &result)
{
    hcFloat eps 	= 1E-5;
    hcFloat theta 	= pos0[1];
    hcFloat phi   	= pos0[2];

    if(fabs(pos0[0]-pos1[0]) > eps * pos0[0])
    {
        printf("ERROR! getPosOnGreatArc: pos0 and pos1 lie on different spheres!\n");
        printf("pos0[0]: %E\npos1[0]: %E\nDiff-r: %E\n", pos0[0], pos1[0], fabs(pos0[0]-pos1[0]));
        printf("pos0:\n");
        pos0.dump();
        printf("pos1:\n");
        pos1.dump();
        return;
    }

    Vec3D pos1transform(pos1);
    pos1transform.transformSphericalCoordSystem(theta, phi);

    hcFloat theta1 	= pos1transform[1] * relPos;
    hcFloat phi1 	= pos1transform[2];

    result			= Vec3D(pos0[0], theta1, phi1);
    result.transformSphericalCoordSystemBack(-theta, -phi);
}

double haversin(double theta)
{
    return (1 - cos(theta)) / 2.0;
}

double ahaversin(double h)
{
    double root = sqrt(h);
    return(2*asin(root));
}

double xsquared(double x)
{
    return x*x;
}

double xpown(double x, double *params)
{
    return params[0] * pow(x, params[1]);
}

double func_identity(double x)
{
    return x;
}

/*! does neither work with uint nor long uint, so just keep it double
 *
 */
double factorial(double n)
{
    if(n==0 || n==1 )
        return 1;
    else
        return(n*factorial(n-1));
}

/*
double binomial(double n, double k)
{
    return((double)factorial(n) / (factorial(k) * factorial(n-k)));
}*/

hcFloat binomial(hcFloat n, uint k)
{
	hcFloat retval = 1.0;

	for(uint i=1; i<=k; ++i)
		retval *= (n+1.0-i) / i;

	return retval;
}

/*
bool loadCSSScoeff(double ***g, double ***h, const char *filename)
{
	if(!doesFileExist(filename))
		return false;

	char line[1000];
	ifstream file(filename);

	uint order = 0;;

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

	*g = new double*[order+1];
	*h = new double*[order+1];

	for(uint l=0; l<=order;++l)
	{
		(*g)[l] = new double[l+1];
		(*h)[l] = new double[l+1];
	}

	(*g)[0][0] = order;
	(*h)[0][0] = order;

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

		(*g)[l][m]	= gval;
		(*h)[l][m]	= hval;
	}
	return true;
}//*/

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Polynomial
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

Polynomial::Polynomial(uint order)
{
    initNULL();
    init(order);
}

Polynomial::Polynomial(double *factor, uint n)
{
    initNULL();
    init(n);

    for(uint i=0; i<=n; ++i)
        this->factor[i] = factor[i];
}

Polynomial::Polynomial(const Polynomial &other)
{
    initNULL();
    *this = other;
}

Polynomial::~Polynomial()
{
    clear();
}

Polynomial &Polynomial::operator=(const Polynomial &other)
{
    if(this == &other)
        return *this;

    init(other.n);
    for(uint i=0; i<=n; ++i)
        factor[i]   = other.factor[i];

    return *this;
}

double Polynomial::operator()(double x)
{
    double result = 0.0;

    for(uint i=0; i<=n; ++i)
    {
        if(isnan((long double)factor[i]) || isnan((long double)pow(x, (double)i)))
        {
            printf("ERROR! Polynomial::operator():\n");
            printf("i: %u, factor[i]: %E, pow: %E, x: %E\n", i, factor[i], pow(x, (double)i), x);
            exit(1);
        }
        result += factor[i] * pow(x,(double)i);
    }

    return result;
}

void Polynomial::initNULL()
{
    this->n = 0;
    factor  = NULL;
}

void Polynomial::clear()
{
    delete [] factor;
    initNULL();
}

void Polynomial::init(uint n)
{
    clear();
    this->n	= n;
    factor 	= new double[n+1];

    for(uint i=0; i<n+1; ++i)
    	factor[i] = 0.0;
}

/*! \brief transforms this to its m'th derivative
 */
void Polynomial::derivative(uint m)
{
    if (m>n)
    {
        init(0);
        factor[0] = 0.0;
        return;
    }

    if(m==0)
        return;

    double *factors = new double[n];

    for(uint i=0; i<n; ++i)
        factors[i] = this->factor[i+1] * (i+1);

    init(n-1);

    for(uint i=0; i<=n; ++i)
        this->factor[i] = factors[i];

    delete [] factors;

    derivative(m-1);
}

void Polynomial::scale(double scale)
{
    for(uint i=0;i<=n;++i)
        factor[i] *= scale;
}

void Polynomial::dump()
{
    printf("Dumping polynomial of order %u:\nP(x) = ", this->n);
    for(uint i=0; i<=n; ++i)
        printf("%E * x^%u + ", factor[i], i);
    printf("\n");
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          LegendrePoly
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


LegendrePoly::LegendrePoly(uint order) :
    Polynomial(order)
{
    initNULL();
	init(order);
}

LegendrePoly::LegendrePoly(const LegendrePoly &other) :
    Polynomial(other)
{}

LegendrePoly::~LegendrePoly()
{
	clear();
}

LegendrePoly &LegendrePoly::operator=(const LegendrePoly &other){

    if(this == &other)
        return *this;

    Polynomial::operator=(other);

    return *this;
}

void LegendrePoly::initNULL()
{}

void LegendrePoly::clear()
{
	Polynomial::clear();
	initNULL();
}

void LegendrePoly::init(uint order)
{
	clear();
	Polynomial::init(order);

	for(uint i=0; i<=order/2; ++i)
	{
		hcFloat sign 		= i%2==0 ? 1.0 : -1.0;
		hcFloat retval		= sign * factorial(2*order-2*i) / (factorial(i) * factorial(order-i) * factorial(order-2*i));
		factor[order-2*i] 	= retval;

		if(std::isnan((double)(factor[order-2*i])))
			printf("Legendrepoly::init: i: %u, order: %u, factorial(2*order-2*i): %E\n", i, order, factorial(2*order-2*i));
	}

	scale(1.0/pow(2.0,n));

	/*
	for(uint k=0; k<=order; ++k)
	{
		hcFloat retval = pow(2.0, n) * binomial(n, k) * binomial((n+k-1.0)/2.0, n);
		printf("n: %u, k: %u, 1: %E, 2: %E, diff: %E\n", n, k, factor[k], retval, factor[k] - retval);
	}
	printf("\n");//*/
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Associated Legendre Polynomials (conventional)
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


AssocLegendrePoly::AssocLegendrePoly(uint l, uint m) :
	dmPl(l),
	dm_1Pl(l)
{
	initNULL();
    init(l, m);
}

AssocLegendrePoly::AssocLegendrePoly(const AssocLegendrePoly &other) :
	dmPl(other.l),
	dm_1Pl(l)
{
	initNULL();
	init(other.l, other.m);
}

AssocLegendrePoly::~AssocLegendrePoly()
{
	clear();
}

AssocLegendrePoly &AssocLegendrePoly::operator=(const AssocLegendrePoly &other)
{
	if(this == &other)
		return *this;

	init(other.l, other.m);

	return *this;
}

hcFloat AssocLegendrePoly::operator ()(hcFloat x)
{
	hcFloat retval 	= pow(1.0-x*x, m/2.0) * dmPl(x);
	return retval;
}

void AssocLegendrePoly::initNULL()
{
	this->l = 0;
	this->m	= 0;
}

void AssocLegendrePoly::clear()
{
	dmPl.clear();
	dm_1Pl.clear();
	initNULL();
}

void AssocLegendrePoly::init(uint l, uint m)
{
	clear();
	dmPl.init(l);
	dmPl.derivative(m);
	dm_1Pl = dmPl;
	dm_1Pl.derivative(1);
	this->l = l;
	this->m = m;
}

hcFloat AssocLegendrePoly::getDeriv(hcFloat theta)
{
	hcFloat retval	= (cos(theta) * m * pow(sin(theta), (hcFloat)((int)m-1)) * dmPl(cos(theta)) - pow(sin(theta), (hcFloat)(m+1)) * dm_1Pl(cos(theta)));

	return retval;
}

void AssocLegendrePoly::dump()
{
	printf("AssocLegendrePoly_conv::dump l: %u, m: %u, dmPl:\n", l, m);
	dmPl.dump();
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Associated Legendre Polynomials (Sun)
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


AssocLegendrePoly_sun::AssocLegendrePoly_sun(uint l, uint m) :
    AssocLegendrePoly(l,m)
{
	initNULL();
	init(l, m);
}

AssocLegendrePoly_sun::AssocLegendrePoly_sun(const AssocLegendrePoly_sun &other) :
	AssocLegendrePoly(other.l, other.m)
{
	initNULL();
	AssocLegendrePoly::operator =(other);
}

AssocLegendrePoly_sun::~AssocLegendrePoly_sun()
{
	clear();
}

AssocLegendrePoly_sun &AssocLegendrePoly_sun::operator=(const AssocLegendrePoly_sun &other)
{
	if(this == &other)
		return *this;

	AssocLegendrePoly::operator =(other);

	return *this;
}

hcFloat AssocLegendrePoly_sun::operator ()(hcFloat x)
{
	hcFloat q 		= m==0 		? 1.0 : 2.0;
	hcFloat clm		= sqrt(q * (hcFloat)factorial(l-m) / (hcFloat)factorial(l+m));
	hcFloat retval	= clm * AssocLegendrePoly::operator ()(x);
    return retval;
}

hcFloat AssocLegendrePoly_sun::getDeriv(hcFloat theta)
{
	hcFloat q 		= (m==0 ? 1.0 : 2.0);
	hcFloat clm		= sqrt(q * (hcFloat)factorial(l-m)/(hcFloat)factorial(l+m));
	hcFloat term1	= m * cos(theta) * pow(sin(theta), (hcFloat)((int)m-1)) * dmPl(cos(theta));
	hcFloat term2	= pow(sin(theta), (hcFloat)(m+1)) * dm_1Pl(cos(theta));
    hcFloat retval	= clm * (term1 - term2);

    return retval;
}
