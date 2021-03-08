#include "engine/math/hcVec.h"
#include "engine/math/hcMatrix.h"
#include "engine/hcTools.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef CUDA
#include "src/cuda_interface.h"
#endif

void getSphericGradientVectors(const Vec3D &cartPos, Vec3D &er, Vec3D &et, Vec3D &ep)
{
	Vec3D posS		= cartPos.convCoordCart2Spher();
	//printf("cart:\n");
	//cartPos.dump();
	//printf("spher:\n");
	//posS.dump();

	hcFloat r		= posS[0];
#ifdef RSCALE
	r /= r_sol;
#endif
	hcFloat t 		= posS[1];
	hcFloat p 		= posS[2];

	hcFloat st		= sin(t);
	hcFloat sp 		= sin(p);
	hcFloat ct		= cos(t);
	hcFloat cp		= cos(p);

#ifdef SPHERICUNITVEC
	er 		= Vec3D(												   			 st*cp, 		 st*sp, 		 ct);
	et 		= r<num_eps 				? Vec3D(1.0, 0.0, 0.0) : Vec3D(			 ct*cp,			 ct*sp, 		-st);
	ep 		= r<num_eps || st<num_eps 	? Vec3D(0.0, 1.0, 0.0) : Vec3D(			-sp, 			 cp, 	(hcFloat)0.0);
#else
	er 		= Vec3D(												   			 st*cp, 		 st*sp, 		 ct);
	et 		= r<num_eps 				? Vec3D(1.0, 0.0, 0.0) : Vec3D( 1/r		*ct*cp,	1/r		*ct*sp, 	-1/r*st);
	ep 		= r<num_eps || st<num_eps 	? Vec3D(0.0, 1.0, 0.0) : Vec3D(-1/(r*st)*sp, 	1/(r*st)*cp, 	(hcFloat)0.0);
#endif
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Vec2D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

Vec2D::Vec2D(){}

template<class S>
Vec2D::Vec2D(S x, S y){

	this->content[0] = x;
	this->content[1] = y;
}

Vec2D::Vec2D(const Vec2D &other) : Vec<2, hcFloat>(other){}

template<class S>
Vec2D::Vec2D(const Vec<2, S> &other) : Vec<2, hcFloat>(other){}

Vec2D::~Vec2D(){}

Vec2D &Vec2D::operator=(const Vec2D &other)
{
    if(this == &other)
        return *this;

    Vec<2, hcFloat>::operator=(other);

	return *this;
}

template<class S>
Vec2D &Vec2D::operator=(const Vec<2, S> &other)
{
    Vec<2, hcFloat>::operator=(other);

	return *this;
}

Vec2D operator+(Vec2D lhs, const Vec2D &rhs)
{
	lhs += rhs;
	return lhs;
}


Vec2D operator-(Vec2D lhs, const Vec2D &rhs)
{
	lhs -= rhs;
	return lhs;
}

hcFloat operator*(const Vec2D &lhs, const Vec2D &rhs)
{
	return lhs.sp(rhs);
}

template<class S>
void Vec2D::rotate(S phi)
{
    hcFloat temp 	= content[0];
    hcFloat sine 	= sin(phi);
    hcFloat cosine 	= cos(phi);

    content[0] = cosine * content[0] - sine * content[1];
    content[1] = sine * temp         + cosine * content[1];
}

template<class S>
void Vec2D::set(S x, S y)
{
	content[0] = x;
	content[1] = y;
}


//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Vec3D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

Vec3D::Vec3D() : Vec<3, hcFloat>(){}

Vec3D::Vec3D(const Vec3D &other) : Vec<3, hcFloat>(other){}

template<class S>
Vec3D::Vec3D(const Vec<3, S> &other) : Vec<3, hcFloat>(other){}

template<class S>
Vec3D::Vec3D(S x, S y, S z) : Vec<3, hcFloat>()
{
	content[0] = x;
	content[1] = y;
	content[2] = z;
}

Vec3D::~Vec3D(){}

Vec3D &Vec3D::operator=(const Vec3D &other)
{
    if(this == &other)
        return *this;

    Vec<3, hcFloat>::operator=(other);

	return *this;
}

template<class S>
Vec3D &Vec3D::operator=(const Vec<3, S> &other)
{
    Vec<3, hcFloat>::operator=(other);

	return *this;
}


Vec3D operator+(Vec3D lhs, const Vec3D &rhs)
{
    lhs += rhs;
	return lhs;
}


Vec3D operator-(Vec3D lhs, const Vec3D &rhs)
{
    lhs -= rhs;
	return lhs;
}

template<class S>
void Vec3D::rotateX(S angle)
{
    Matrix4x4 rotmat;
    rotmat.loadIdentity();
    rotmat.rotatex(angle);

    Vec4D tempvec, resultvec;
    for(uint i=0;i<3;++i)
        tempvec(i) = content[i];
    tempvec(3) = 1.0;
    resultvec = rotmat * tempvec;

    for(uint i=0;i<3;++i)
        content[i] = resultvec[i];
}

template<class S>
void Vec3D::rotateY(S angle)
{
    Matrix4x4 rotmat;
    rotmat.loadIdentity();
    rotmat.rotatey(angle);

    Vec4D tempvec, resultvec;
    for(uint i=0;i<3;++i)
        tempvec(i) = content[i];
    tempvec(3) = 1.0;
    resultvec = rotmat * tempvec;

    for(uint i=0;i<3;++i)
        content[i] = resultvec[i];
}

template<class S>
void Vec3D::rotateZ(S angle)
{
    Matrix4x4 rotmat;
    rotmat.loadIdentity();
    rotmat.rotatez(angle);

    Vec4D tempvec, resultvec;
    for(uint i=0;i<3;++i)
        tempvec(i) = content[i];
    tempvec(3) = 1.0;
    resultvec = rotmat * tempvec;

    for(uint i=0;i<3;++i)
        content[i] = resultvec[i];
}

// --------------------------------------------------------------------------------------------------------
// in-situ conversion functions
// --------------------------------------------------------------------------------------------------------

void Vec3D::convertEllipticalCoordsToCartesian(hcFloat a, bool prolate)
{
	*this = this->convCoordSpher2Cart();
	if(prolate)
	{
		content[2] *= a;
	}
	else
	{
		content[0] *= a;
		content[1] *= a;
	}
}

void Vec3D::convertCartesianCoordsToElliptical(hcFloat a, bool prolate)
{
	if(prolate)
	{
		content[2] /= a;
	}
	else
	{
		content[0] /= a;
		content[1] /= a;
	}
	*this = this->convCoordCart2Spher();
}

// --------------------------------------------------------------------------------------------------------
// in-situ conversion functions end
// --------------------------------------------------------------------------------------------------------

// --------------------------------------------------------------------------------------------------------
// conversion functions
// --------------------------------------------------------------------------------------------------------

Vec3D Vec3D::convCoordSpher2Cart() const
{
	hcFloat r = content[0];
	hcFloat t = content[1];
	hcFloat p = content[2];

	hcFloat x = fabs(r) < num_eps ? 0.0 : r * sin(t) * cos(p);
	hcFloat y = fabs(r) < num_eps ? 0.0 : r * sin(t) * sin(p);
	hcFloat z = fabs(r) < num_eps ? 0.0 : r * cos(t);

	Vec3D retval;
	retval(0) = x;
	retval(1) = y;
	retval(2) = z;

	return retval;
}

Vec3D Vec3D::convCoordCart2Spher() const
{
	hcFloat x = content[0];
	hcFloat y = content[1];
	hcFloat z = content[2];

	//printf("%E/%E/%E\n", x, y, z);

	hcFloat r = fabs(x)<num_eps && fabs(y)<num_eps && fabs(z)<num_eps ? 0.0 : sqrt(x*x + y*y + z*z);
	hcFloat t = fabs(x)<num_eps && fabs(y)<num_eps && fabs(z)<num_eps ? 0.0 : acos(z/r);
	hcFloat p = fabs(x)<num_eps && fabs(y)<num_eps && fabs(z)<num_eps ? 0.0 : atan2(y, x);

	while (p < 0)
		p += 2 * PI;

	Vec3D retval;
	retval(0) = r;
	retval(1) = t;
	retval(2) = p;

	//printf("%E/%E/%E\n", retval(0), retval(1), retval(2));

	return retval;
}

Vec3D Vec3D::convVecCart2Spher(const Vec3D &cartPos)
{
	Vec3D er, et, ep;
	getSphericGradientVectors(cartPos, er, et, ep);

	//er.dump();
	//et.dump();
	//ep.dump();
	//exit(1);

	Vec3D vecCart = *this;

	Matrix<3,4,double> M;
	M(0,0) = er[0];	M(0,1) = et[0];	M(0,2) = ep[0]; M(0,3) = vecCart[0];
	M(1,0) = er[1];	M(1,1) = et[1];	M(1,2) = ep[1]; M(1,3) = vecCart[1];
	M(2,0) = er[2];	M(2,1) = et[2];	M(2,2) = ep[2]; M(2,3) = vecCart[2];

	Vec3D retval;
	uint sol = M.solveSLE(retval);

	return retval;
}

Vec3D Vec3D::convVecSpher2Cart(const Vec3D &cartPos) const
{
	Vec3D er, et, ep;
	getSphericGradientVectors(cartPos, er, et, ep);

	Vec3D retval 	= er * content[0] + et * content[1] + ep * content[2];
	return retval;
}

/* 	@param e1x		origin unit vector 1
 * 	@param e1y		origin unit vector 2
 * 	@param e1z		origin unit vector 3
 * 	@param e2x		destination unit vector 1
 * 	@param e2y		destination unit vector 2
 * 	@param e2z		destination unit vector 3
 *
 */
Vec3D Vec3D::convCoordCart2Cart(Vec3D e1x, Vec3D e1y, Vec3D e1z, Vec3D e2x, Vec3D e2y, Vec3D e2z) const
{
	e1x.scale(1.0/e1x.length());
	e1y.scale(1.0/e1y.length());
	e1z.scale(1.0/e1z.length());

	e2x.scale(1.0/e2x.length());
	e2y.scale(1.0/e2y.length());
	e2z.scale(1.0/e2z.length());



	hcFloat cos_11	= e2x * e1x;
	hcFloat cos_12	= e2x * e1y;
	hcFloat cos_13	= e2x * e1z;

	hcFloat cos_21	= e2y * e1x;
	hcFloat cos_22	= e2y * e1y;
	hcFloat cos_23	= e2y * e1z;

	hcFloat cos_31	= e2z * e1x;
	hcFloat cos_32	= e2z * e1y;
	hcFloat cos_33	= e2z * e1z;

	Matrix3x3 rotmat;

	rotmat(0,0)		= cos_11;
	rotmat(0,1)		= cos_12;
	rotmat(0,2)		= cos_13;

	rotmat(1,0)		= cos_21;
	rotmat(1,1)		= cos_22;
	rotmat(1,2)		= cos_23;

	rotmat(2,0)		= cos_31;
	rotmat(2,1)		= cos_32;
	rotmat(2,2)		= cos_33;

	//rotmat.dump();

	return rotmat * *this;
}//*/

// --------------------------------------------------------------------------------------------------------
// conversion functions end
// --------------------------------------------------------------------------------------------------------

void Vec3D::stereoProjUnitSphere(Vec2D tangent)
{
    Vec3D cartTangentPoint;
    cartTangentPoint(0) = 1.0;
    cartTangentPoint(1) = tangent[0];
    cartTangentPoint(2) = tangent[1];

    Vec3D pPlusT;
    pPlusT = *this;
    pPlusT += cartTangentPoint;

    double s = 2.0  / cartTangentPoint.sp(pPlusT);
    cartTangentPoint.scale(s-1);
    this->scale(s);
    //this->plus(cartTangentPoint);
    *this += cartTangentPoint;
}

template<class S>
void Vec3D::transformSphericalCoordSystem(S theta, S phi)
{
    Matrix4x4 rotmat;
    rotmat.loadIdentity();
    rotmat.rotatey(-theta);

    hcFloat radius 	= content[0];
	content[0] 		= 1.0;
    content[2] 		= cyclicPhi(content[2] - phi);
    *this 			= this->convCoordSpher2Cart();

    Vec4D t1, t2;

    for(uint j=0;j<3;++j)
        t1.content[j] = content[j];

    t1.content[3] 	= 1.0;
    t2 				= rotmat * t1;

    for(uint j=0;j<3;++j)
        content[j] = t2[j];

    *this 			= this->convCoordCart2Spher();
    content[0] 		= radius;
}

template<class S>
void Vec3D::transformSphericalCoordSystemBack(S theta, S phi)
{
    Matrix4x4 rotmat;
    rotmat.loadIdentity();
    rotmat.rotatey(-theta);

    Vec4D temp1, temp2;

    hcFloat radius = content[0];
    content[0] = 1.0;

    *this = this->convCoordSpher2Cart();

    for(uint j=0;j<3;++j)
        temp1.content[j] = content[j];
    temp1.content[3] = 1.0;
    temp2 = rotmat * temp1;
    for(uint j=0;j<3;++j)
        content[j] = temp2[j];

    *this = this->convCoordCart2Spher();

    content[2] -= phi;
    if(content[2] < 0)
        content[2] += 2 * PI;
    if(content[2] >= 2 * PI)
        content[2] -= 2 * PI;

    content[0] = radius;
}


void Vec3D::cp(const Vec3D &v, Vec3D *res) const
{
    res->content[0] = this->content[1] * v.content[2] - this->content[2] * v.content[1];
    res->content[1] = this->content[2] * v.content[0] - this->content[0] * v.content[2];
    res->content[2] = this->content[0] * v.content[1] - this->content[1] * v.content[0];
}


Vec3D Vec3D::cp(const Vec3D &other) const
{
    Vec3D retval;
    retval(0) = content[1] * other[2] - content[2] * other[1];
    retval(1) = content[2] * other[0] - content[0] * other[2];
    retval(2) = content[0] * other[1] - content[1] * other[0];

    return retval;
}


hcFloat Vec3D::getAngle(const Vec3D &other) const
{
    return acos(sp(other) / (this->length() * other.length()));
}

// TEST THIS!

hcFloat Vec3D::getAngle2(const Vec3D &other) const
{
    return atan2(this->cp(other).length(), (*this) * other);
}


void Vec3D::dumpSphericalCoords()
{
	  printf("%E r_sol\n%E PI\n%E PI\n\n", content[0] / r_sol, content[1] / PI, content[2] / PI);
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Vec4D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

Vec4D::Vec4D(){}

template<class S>
Vec4D::Vec4D(S x, S y, S z, S w) : Vec<4, hcFloat>(){

	content[0] = x;
	content[1] = y;
	content[2] = z;
	content[3] = w;
}

template<class S>
void Vec4D::setPos(S x, S y){

	content[0] = x;
	content[1] = y;
}

template<class S>
void Vec4D::setTex(S u, S v){

	content[2] = u;
	content[3] = v;
}


Vec4D::Vec4D(const Vec4D &other) : Vec<4, hcFloat>(other){}

template<class S>
Vec4D::Vec4D(const Vec<4, S> &other) : Vec<4, hcFloat>(other){}

Vec4D::~Vec4D(){}

Vec4D &Vec4D::operator=(const Vec4D &other){

    if(this == &other)
        return *this;

    Vec<4, hcFloat>::operator=(other);

	return *this;
}

template<class S>
Vec4D &Vec4D::operator=(const Vec<4, S> &other){

    Vec<4, hcFloat>::operator=(other);

	return *this;
}


Vec4D operator+(Vec4D lhs, const Vec4D &rhs){

	lhs += rhs;
	return lhs;
}


Vec4D operator-(Vec4D lhs, const Vec4D &rhs){

	lhs -= rhs;
	return lhs;
}

Vec4D::Vec4D(const Vec3D &other){

    content[0] = other[0];
    content[1] = other[1];
    content[2] = other[2];
    content[3] = 1.0;
}

template<class S>
void Vec4D::set(S x, S y, S z, S w){

    content[0] = x;
    content[1] = y;
    content[2] = z;
    content[3] = w;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Vec5D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

Vec5D::Vec5D() : Vec<5, hcFloat>(){}

template<class S>
Vec5D::Vec5D(S x, S y, S z, S u, S v) : Vec<5, hcFloat>(){

	content[0] = x;
	content[1] = y;
	content[2] = z;
	content[3] = u;
	content[4] = v;
}

Vec5D::~Vec5D(){}

template<class S>
void Vec5D::set(S x, S y, S z, S u, S v){

	content[0] = x;
	content[1] = y;
	content[2] = z;
	content[3] = u;
	content[4] = v;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Explicit template instantiation
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

template Vec2D::Vec2D<float>(float x, float y);
template Vec2D::Vec2D<float>(const Vec<2, float> &other);
template Vec2D &Vec2D::operator=<float>(const Vec<2, float> &other);
template void Vec2D::rotate<float>(float phi);
template void Vec2D::set<float>(float x, float y);

template Vec2D::Vec2D<double>(double x, double y);
template Vec2D::Vec2D<double>(const Vec<2, double> &other);
template Vec2D &Vec2D::operator=<double>(const Vec<2, double> &other);
template void Vec2D::rotate<double>(double phi);
template void Vec2D::set<double>(double x, double y);

template Vec2D::Vec2D<long double>(long double x, long double y);
template Vec2D::Vec2D<long double>(const Vec<2, long double> &other);
template Vec2D &Vec2D::operator=<long double>(const Vec<2, long double> &other);
template void Vec2D::rotate<long double>(long double phi);
template void Vec2D::set<long double>(long double x, long double y);



template Vec3D::Vec3D<float>(const Vec<3, float> &other);
template Vec3D::Vec3D<float>(float x, float y, float z);
template Vec3D &Vec3D::operator=<float>(const Vec<3, float> &other);
template void Vec3D::rotateX<float>(float angle);
template void Vec3D::rotateY<float>(float angle);
template void Vec3D::rotateZ<float>(float angle);
template void Vec3D::transformSphericalCoordSystem<float>(float theta, float phi);
template void Vec3D::transformSphericalCoordSystemBack<float>(float theta, float phi);

template Vec3D::Vec3D<double>(const Vec<3, double> &other);
template Vec3D::Vec3D<double>(double x, double y, double z);
template Vec3D &Vec3D::operator=<double>(const Vec<3, double> &other);
template void Vec3D::rotateX<double>(double angle);
template void Vec3D::rotateY<double>(double angle);
template void Vec3D::rotateZ<double>(double angle);
template void Vec3D::transformSphericalCoordSystem<double>(double theta, double phi);
template void Vec3D::transformSphericalCoordSystemBack<double>(double theta, double phi);

template Vec3D::Vec3D<long double>(const Vec<3, long double> &other);
template Vec3D::Vec3D<long double>(long double x, long double y, long double z);
template Vec3D &Vec3D::operator=<long double>(const Vec<3, long double> &other);
template void Vec3D::rotateX<long double>(long double angle);
template void Vec3D::rotateY<long double>(long double angle);
template void Vec3D::rotateZ<long double>(long double angle);
template void Vec3D::transformSphericalCoordSystem<long double>(long double theta, long double phi);
template void Vec3D::transformSphericalCoordSystemBack<long double>(long double theta, long double phi);

template Vec4D::Vec4D<float>(float x, float y, float z, float w);
template Vec4D::Vec4D<float>(const Vec<4, float> &other);
template Vec4D& Vec4D::operator=<float>(const Vec<4, float> &other);
template void Vec4D::set<float>(float x, float y, float z, float w);

template Vec4D::Vec4D<double>(double x, double y, double z, double w);
template Vec4D::Vec4D<double>(const Vec<4, double> &other);
template Vec4D& Vec4D::operator=<double>(const Vec<4, double> &other);
template void Vec4D::set<double>(double x, double y, double z, double w);

template Vec4D::Vec4D<long double>(long double x, long double y, long double z, long double w);
template Vec4D::Vec4D<long double>(const Vec<4, long double> &other);
template Vec4D& Vec4D::operator=<long double>(const Vec<4, long double> &other);
template void Vec4D::set<long double>(long double x, long double y, long double z, long double w);

template Vec5D::Vec5D<float>(float x, float y, float z, float u, float v);
template void Vec5D::set<float>(float x, float y, float z, float u, float v);

template void Vec5D::set<double>(double x, double y, double z, double u, double v);
template Vec5D::Vec5D<double>(double x, double y, double z, double u, double v);
