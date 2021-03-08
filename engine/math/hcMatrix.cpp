
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "engine/math/hcVec.h"
#include "engine/math/hcMatrix.h"
#include "engine/hcConstants.h"

#ifdef GUI
#include "engine/hcQuaternion.h"
#endif

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Matrix2x2
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

Matrix2x2::Matrix2x2() : MatrixNxN<2, hcFloat>(){}

Matrix2x2::Matrix2x2(const Matrix2x2 &other) : MatrixNxN<2, hcFloat>(other){}

template<class S>
Matrix2x2::Matrix2x2(const MatrixNxN<2, S> &other) : MatrixNxN<2, hcFloat>(other){}

Matrix2x2::~Matrix2x2(){}

template<class S>
void Matrix2x2::scalex(S scale){

	scaleAxis(0,scale);
}

template<class S>
void Matrix2x2::scaley(S scale){

	scaleAxis(1,scale);
}

template<class S>
void Matrix2x2::scale(S scalex, S scaley){

	scaleAxis(0,scalex);
	scaleAxis(1,scaley);
}


Matrix2x2 &Matrix2x2::operator=(const Matrix2x2 &other){

	if(this == &other)
		return * this;

	MatrixNxN<2, hcFloat>::operator=(other);

	return *this;
}

template<class S>
Matrix2x2 &Matrix2x2::operator=(const MatrixNxN<2, S> &other){

	MatrixNxN<2, hcFloat>::operator=(other);

	return *this;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Matrix3x3
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

Matrix3x3::Matrix3x3() : MatrixNxN<3, hcFloat>(){}

Matrix3x3::Matrix3x3(const Matrix3x3 &other) : MatrixNxN<3, hcFloat>(other){}

template<class S>
Matrix3x3::Matrix3x3(const MatrixNxN<3, S> &other) : MatrixNxN<3, hcFloat>(other){}

template<class S>
Matrix3x3::Matrix3x3(const Vec<3, S> &vec1, const Vec<3, S> &vec2, const Vec<3, S> &vec3){

    for(uint i=0;i<3;++i)
    {
        content[i*3 + 0] = vec1[i];
        content[i*3 + 1] = vec2[i];
        content[i*3 + 2] = vec3[i];
    }
}

Matrix3x3::~Matrix3x3(){}

template<class S>
void Matrix3x3::scale(S scalex, S scaley, S scalez){

    scaleAxis(0,scalex);
    scaleAxis(1,scaley);
    scaleAxis(2,scalez);
}

Matrix3x3 &Matrix3x3::operator=(const Matrix3x3 &other){

	if(this == &other)
		return *this;

	MatrixNxN<3, hcFloat>::operator=(other);

	return *this;
}

template<class S>
Matrix3x3 &Matrix3x3::operator=(const MatrixNxN<3, S> &other){

	MatrixNxN<3, hcFloat>::operator=(other);

	return *this;
}

template<class S>
void Matrix3x3::scalex(S scale){

	scaleAxis(0,scale);
}

template<class S>
void Matrix3x3::scaley(S scale){

	scaleAxis(1, scale);
}

template<class S>
void Matrix3x3::scalez(S scale)
{
	scaleAxis(2, scale);
}

void Matrix3x3::convertSphericalToCartesian(const Vec3D &cartPos)
{
	Vec3D posS	= cartPos.convCoordCart2Spher();
    double r	= posS[0];
    double t	= posS[1];
    double p	= posS[2];

    this->loadZeroes();

    content[0] = sin(t) * cos(p);
    content[1] = cos(t) * cos(p);
    content[2] = -sin(p);

    content[3] = sin(t) * sin(p);
    content[4] = cos(t) * sin(p);
    content[5] = cos(p);

    content[6] = cos(t);
    content[7] = -sin(t);
    content[8] = 0.0;
}

void Matrix3x3::convertCartesianToSpherical(const Vec3D &cartPos)
{

	Vec3D posS	= cartPos.convCoordCart2Spher();
	double r	= posS[0];
	double t	= posS[1];
	double p	= posS[2];

    content[0] = sin(t) * cos(p);
    content[1] = sin(t) * sin(p);
    content[2] = cos(t);

    content[3] = cos(t) * cos(p);
    content[4] = cos(t) * sin(p);
    content[5] = -sin(t);

    content[6] = -sin(p);
    content[7] = cos(p);
    content[8] = 0.0;
}

template<class S>
void Matrix3x3::loadRotationX(S theta)
{
	loadZeroes();
	content[0]	= 1.0;
	content[4]	= cosf(theta);
	content[5]	= -sinf(theta);
	content[7]	= sinf(theta);
	content[8]	= cosf(theta);
}

template<class S>
void Matrix3x3::loadRotationY(S theta)
{
	loadZeroes();
	content[0]	= cosf(theta);
	content[2]	= sinf(theta);
	content[4]	= 1.0;
	content[6]	= -sinf(theta);
	content[8]	= cosf(theta);
}

template<class S>
void Matrix3x3::loadRotationZ(S theta)
{
	loadZeroes();
	content[0]	= cosf(theta);
	content[1]	= -sinf(theta);
	content[3]	= sinf(theta);
	content[4]	= cosf(theta);
	content[8]	= 1.0;
}

template<class S>
void Matrix3x3::rotateAroundAxis(Vec3D axis, S angle)
{
	hcFloat cost	= cos(angle);
	hcFloat sint	= sin(angle);
	axis.scale(1.0/axis.length());
	hcFloat ux		= axis[0];
	hcFloat uy		= axis[1];
	hcFloat uz		= axis[2];

	operator()(0, 0) 	= cost + ux*ux*(1-cost);
	operator()(0, 1) 	= ux*uy*(1-cost)-uz*sint;
	operator()(0, 2) 	= ux*uz*(1-cost)+uy*sint;

	operator()(1, 0) 	= uy*ux*(1-cost)+uz*sint;
	operator()(1, 1) 	= cost+uy*uy*(1-cost);
	operator()(1, 2) 	= uy*uz*(1-cost)-ux*sint;

	operator()(2, 0) 	= uz*ux*(1-cost)-uy*sint;
	operator()(2, 1) 	= uz*uy*(1-cost)+ux*sint;
	operator()(2, 2) 	= cost+uz*uz*(1-cost);
}

// see medelung 1964 Die mathematischen Hilfsmittel des Physikers
template<class S>
void Matrix3x3::loadEulerTransform(S omega, S theta, S phi)
{
	loadZeroes();
	content[0]	=  cosf(phi)*cosf(omega) - sinf(phi)*sinf(omega)*cosf(theta);
	content[1]	=  cosf(phi)*sinf(omega) + sinf(phi)*cosf(omega)*cosf(theta);
	content[2]	=  sinf(phi)*sinf(theta);

	content[3]	= -sinf(phi)*cosf(omega) - cosf(phi)*sinf(omega)*cosf(theta);
	content[4]	= -sinf(phi)*sinf(omega) + cosf(phi)*cosf(omega)*cosf(theta);
	content[5]	=  cosf(phi)*sinf(theta);

	content[6]	=  sinf(omega)*sinf(theta);
	content[7]	= -cosf(omega)*sinf(theta);
	content[8]	=  cosf(theta);
}


//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Matrix4x4
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

Matrix4x4::Matrix4x4() : MatrixNxN<4, hcFloat>() {}

Matrix4x4::Matrix4x4(const Matrix4x4 &other) : MatrixNxN<4, hcFloat>(other){}

template<class S>
Matrix4x4::Matrix4x4(const MatrixNxN<4, S> &other) : MatrixNxN<4, hcFloat>(other){}

Matrix4x4::~Matrix4x4(){}

Matrix4x4 &Matrix4x4::operator=(const Matrix4x4 &other){

	if(this == &other)
		return *this;

	MatrixNxN<4, hcFloat>::operator=(other);

	return *this;
}

template<class S>
Matrix4x4 &Matrix4x4::operator=(const MatrixNxN<4, S> &other){

	MatrixNxN<4, hcFloat>::operator=(other);

	return *this;
}

template<class S>
void Matrix4x4::scalex(S scale){

	scaleAxis(0,scale);
}

template<class S>
void Matrix4x4::scaley(S scale){

	scaleAxis(1,scale);
}

template<class S>
void Matrix4x4::scalez(S scale){

	scaleAxis(2,scale);
}

template<class S>
void Matrix4x4::scalew(S scale){

	scaleAxis(3,scale);
}

template<class S>
void Matrix4x4::scale(S scalex, S scaley, S scalez, S scalew){

    scaleAxis(0,scalex);
    scaleAxis(1,scaley);
    scaleAxis(2,scalez);
    scaleAxis(3,scalew);
}

void Matrix4x4::loadTFMatrix(const Vec3D &right, const Vec3D &up, const Vec3D &back, const Vec3D &pos){

    uint n = 4;
    for(uint i=0;i<n-1;++i)
    {
        this->content[i*n+0] = right.content[i];
        this->content[i*n+1] = up.content[i];
        this->content[i*n+2] = back.content[i];
        this->content[i*n+3] = pos.content[i];
    }
    this->content[3*4 + 0] = 0.0;
    this->content[3*4 + 1] = 0.0;
    this->content[3*4 + 2] = 0.0;
    this->content[3*4 + 3] = 1.0;

}

template<class S>
void Matrix4x4::rotatex(S phi){

    Matrix4x4 rotmat;
    rotmat.content[0] 	= 1;
    rotmat.content[5] 	= cosf(phi);
    rotmat.content[6] 	= -sinf(phi);
    rotmat.content[9] 	= sinf(phi);
    rotmat.content[10] 	= cosf(phi);
    rotmat.content[15] 	= 1;
    *this = rotmat * *this;
}

template<class S>
void Matrix4x4::rotatey(S phi){

    Matrix4x4 rotmat;
    rotmat.content[0] 	= cosf(phi);
    rotmat.content[2] 	= sinf(phi);
    rotmat.content[5] 	= 1;
    rotmat.content[8] 	= -sinf(phi);
    rotmat.content[10] 	= cosf(phi);
    rotmat.content[15] 	= 1;
    *this = rotmat * *this;
}

template<class S>
void Matrix4x4::rotatez(S phi){

    Matrix4x4 rotmat;
    rotmat.content[0] 	= cosf(phi);
    rotmat.content[1] 	= -sinf(phi);
    rotmat.content[4] 	= sinf(phi);
    rotmat.content[5] 	= cosf(phi);
    rotmat.content[10] 	= 1;
    rotmat.content[15] 	= 1;
    *this = rotmat * *this;
}

// http://math.stackexchange.com/questions/293116/rotating-one-3-vector-to-another
// Rodrigues rotation formula
void Matrix4x4::rotateAtoB(const Vec3D a, const Vec3D b)
{
    double eps  = 1E-9;

    double theta1 = a.getAngle(b);
    double theta2 = a.getAngle2(b);

    //printf("theta1: %E\ntheta2: %E\n", theta1, theta2);

    if(theta2 < eps)
    {
        loadIdentity();
        return;
    }

    Vec3D x;
    if(PI - theta2 < eps)
    {
        if(fabs(a[0]) < fabs(a[1]))
        	x = fabs(a[0]) < fabs(a[2]) ? Vec3D(1.0, 0.0, 0.0) : Vec3D(0.0, 0.0, 1.0);
        else
        	x = fabs(a[1]) < fabs(a[2]) ? Vec3D(0.0, 1.0, 0.0) : Vec3D(0.0, 0.0, 1.0);
    }
    else
    {
        x = a.cp(b);
        x.normalize();
    }

    Matrix3x3 A;
    A.content[0] =  0.0;
    A.content[1] = -x[2];
    A.content[2] =  x[1];
    A.content[3] =  x[2];
    A.content[4] =  0.0;
    A.content[5] = -x[0];
    A.content[6] = -x[1];
    A.content[7] =  x[0];
    A.content[8] =  0.0;

    Matrix3x3 rot;
    rot.loadIdentity();

    rot += sin(theta2) * A + (1 - cos(theta2)) * A * A;

    this->loadZeroes();
    this->operator()(0, 0) = rot(0, 0);
    this->operator()(0, 1) = rot(0, 1);
    this->operator()(0, 2) = rot(0, 2);

    this->operator()(1, 0) = rot(1, 0);
    this->operator()(1, 1) = rot(1, 1);
    this->operator()(1, 2) = rot(1, 2);

    this->operator()(2, 0) = rot(2, 0);
    this->operator()(2, 1) = rot(2, 1);
    this->operator()(2, 2) = rot(2, 2);

    this->operator()(3, 3) = 1.0;
}

template<class S>
void Matrix4x4::translate(S x, S y, S z){

    Vec3D temp(x,y,z);
    ((Matrix*)(this))->translate(temp);
}

#ifdef GUI
void Matrix4x4::set(const aiMatrix4x4 &mat){

    content[0] 	= mat.a1;
    content[1] 	= mat.a2;
    content[2] 	= mat.a3;
    content[3] 	= mat.a4;

    content[4]	= mat.b1;
    content[5] 	= mat.b2;
    content[6] 	= mat.b3;
    content[7] 	= mat.b4;

    content[8] 	= mat.c1;
    content[9] 	= mat.c2;
    content[10] = mat.c3;
    content[11] = mat.c4;

    content[12] = mat.d1;
    content[13] = mat.d2;
    content[14] = mat.d3;
    content[15] = mat.d4;
}
#endif

#ifdef GUI
void Matrix4x4::scaleRotTrans(const Vec3D &scale, const hcQuaternion &rot, const Vec3D &trans){

    content[0] 	= scale[0];
    content[1] 	= 0.0;
    content[2] 	= 0.0;
    content[3] 	= 0.0;

    content[4] 	= 0.0;
    content[5] 	= scale[1];
    content[6] 	= 0.0;
    content[7] 	= 0.0;

    content[8] 	= 0.0;
    content[9] 	= 0.0;
    content[10] = scale[2];
    content[11] = 0.0;

    content[12] = 0.0;
    content[13] = 0.0;
    content[14] = 0.0;
    content[15] = 1.0;

    Matrix4x4 rotation;
    rot.toMatrix(rotation);

    *this = rotation * *this;

    translate(trans[0], trans[1], trans[2]);
}
#endif

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Matrix5x5
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

Matrix5x5::Matrix5x5() : MatrixNxN<5, hcFloat>() {}

Matrix5x5::Matrix5x5(const Matrix5x5 &other) : MatrixNxN<5, hcFloat>(other){}

template<class S>
Matrix5x5::Matrix5x5(const MatrixNxN<5, S> &other) : MatrixNxN<5, hcFloat>(other){}

Matrix5x5::~Matrix5x5(){}

Matrix5x5 &Matrix5x5::operator=(const Matrix5x5 &other){

	if(this == &other)
		return *this;

	MatrixNxN<5, hcFloat>::operator=(other);

	return *this;
}

template<class S>
Matrix5x5 &Matrix5x5::operator=(const MatrixNxN<5, S> &other){

	MatrixNxN<5, hcFloat>::operator=(other);

	return *this;
}

void Matrix5x5::loadTFMatrix(const Vec4D &right, const Vec4D &up, const Vec4D &back, const Vec4D &over, const Vec4D &pos){

    uint n = 5;
    for(uint i=0;i<n-1;++i)
    {
        this->content[i*n+0] = right[i];
        this->content[i*n+1] = up[i];
        this->content[i*n+2] = back[i];
        this->content[i*n+3] = over[i];
        this->content[i*n+4] = pos[i];
    }
    this->content[3*4 + 0] = 0.0;
    this->content[3*4 + 1] = 0.0;
    this->content[3*4 + 2] = 0.0;
    this->content[3*4 + 3] = 0.0;
    this->content[3*4 + 4] = 1.0;

}

template<class S>
void Matrix5x5::rotateXY(S phi){

    Matrix5x5 rotmat;
    rotmat.content[0] 	= cosf(phi);
    rotmat.content[1] 	= sinf(phi);
    rotmat.content[5] 	= -sinf(phi);
    rotmat.content[6] 	= cosf(phi);
    rotmat.content[12] 	= 1.0;
    rotmat.content[18] 	= 1.0;
    rotmat.content[24] 	= 1.0;
    Matrix5x5 test;
    test = rotmat * *this;
    *this = rotmat * *this;
}

template<class S>
void Matrix5x5::rotateXZ(S phi){

    Matrix5x5 rotmat;
    rotmat.content[0] 	= cosf(phi);
    rotmat.content[2] 	= -sinf(phi);
    rotmat.content[6] 	= 1.0;
    rotmat.content[10] 	= sinf(phi);
    rotmat.content[12] 	= cosf(phi);
    rotmat.content[18] 	= 1.0;
    rotmat.content[24] 	= 1.0;
    *this = rotmat * *this;
}

template<class S>
void Matrix5x5::rotateYZ(S phi){

    Matrix5x5 rotmat;
    rotmat.content[0] 	= 1.0;
    rotmat.content[6] 	= cosf(phi);
    rotmat.content[7] 	= sinf(phi);
    rotmat.content[11] 	= -sinf(phi);
    rotmat.content[12] 	= cosf(phi);
    rotmat.content[18] 	= 1.0;
    rotmat.content[24] 	= 1.0;
    *this = rotmat * *this;
}

template<class S>
void Matrix5x5::rotateYW(S phi){

    Matrix5x5 rotmat;
    rotmat.content[0] 	= 1.0;
    rotmat.content[6] 	= cosf(phi);
    rotmat.content[8] 	= -sinf(phi);
    rotmat.content[12] 	= 1.0;
    rotmat.content[16] 	= sinf(phi);
    rotmat.content[18] 	= cosf(phi);
    rotmat.content[24] 	= 1.0;
    *this = rotmat * *this;
}

template<class S>
void Matrix5x5::rotateXW(S phi){

    Matrix5x5 rotmat;
    rotmat.content[0] 	= cosf(phi);
    rotmat.content[3] 	= sinf(phi);
    rotmat.content[6] 	= 1.0;
    rotmat.content[12] 	= 1.0;
    rotmat.content[15] 	= -sinf(phi);
    rotmat.content[18] 	= cosf(phi);
    rotmat.content[24] 	= 1.0;
    *this = rotmat * *this;
}

template<class S>
void Matrix5x5::rotateZW(S phi){

    Matrix5x5 rotmat;
    rotmat.content[0] 	= 1.0;
    rotmat.content[6] 	= 1.0;
    rotmat.content[12] 	= cosf(phi);
    rotmat.content[13] 	= -sinf(phi);
    rotmat.content[17] 	= sinf(phi);
    rotmat.content[18] 	= cosf(phi);
    rotmat.content[24] 	= 1.0;
    *this = rotmat * *this;
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

template MatrixNxN<2, float> operator*(const MatrixNxN<2, float> &lhs, const MatrixNxN<2, float> &rhs);
template MatrixNxN<2, float> operator*(const MatrixNxN<2, float> &lhs, const MatrixNxN<2, double> &rhs);
template MatrixNxN<2, double> operator*(const MatrixNxN<2, double> &lhs, const MatrixNxN<2, float> &rhs);
template MatrixNxN<2, double> operator*(const MatrixNxN<2, double> &lhs, const MatrixNxN<2, double> &rhs);
template MatrixNxN<2, long double> operator*(const MatrixNxN<2, long double> &lhs, const MatrixNxN<2, long double> &rhs);

template MatrixNxN<3, float> operator*(const MatrixNxN<3, float> &lhs, const MatrixNxN<3, float> &rhs);
template MatrixNxN<3, float> operator*(const MatrixNxN<3, float> &lhs, const MatrixNxN<3, double> &rhs);
template MatrixNxN<3, double> operator*(const MatrixNxN<3, double> &lhs, const MatrixNxN<3, float> &rhs);
template MatrixNxN<3, double> operator*(const MatrixNxN<3, double> &lhs, const MatrixNxN<3, double> &rhs);
template MatrixNxN<3, long double> operator*(const MatrixNxN<3, long double> &lhs, const MatrixNxN<3, long double> &rhs);

template MatrixNxN<4, float> operator*(const MatrixNxN<4, float> &lhs, const MatrixNxN<4, float> &rhs);
template MatrixNxN<4, float> operator*(const MatrixNxN<4, float> &lhs, const MatrixNxN<4, double> &rhs);
template MatrixNxN<4, double> operator*(const MatrixNxN<4, double> &lhs, const MatrixNxN<4, float> &rhs);
template MatrixNxN<4, double> operator*(const MatrixNxN<4, double> &lhs, const MatrixNxN<4, double> &rhs);

template Matrix2x2::Matrix2x2<float>(const MatrixNxN<2, float> &other);
template void Matrix2x2::scalex<float>(float scale);
template void Matrix2x2::scaley<float>(float scale);
template void Matrix2x2::scale<float>(float scalex, float scaley);
template Matrix2x2 &Matrix2x2::operator=<float>(const MatrixNxN<2, float> &other);

template Matrix2x2::Matrix2x2<double>(const MatrixNxN<2, double> &other);
template void Matrix2x2::scalex<double>(double scale);
template void Matrix2x2::scaley<double>(double scale);
template void Matrix2x2::scale<double>(double scalex, double scaley);
template Matrix2x2 &Matrix2x2::operator=<double>(const MatrixNxN<2, double> &other);



template Matrix3x3::Matrix3x3<float>(const MatrixNxN<3, float> &other);
template Matrix3x3::Matrix3x3<float>(const Vec<3, float> &vec1, const Vec<3, float> &vec2, const Vec<3, float> &vec3);
template Matrix3x3 &Matrix3x3::operator=<float>(const MatrixNxN<3, float> &other);
template void Matrix3x3::scale<float>(float scalex, float scaley, float scalez);
template void Matrix3x3::scalex<float>(float scale);
template void Matrix3x3::scaley<float>(float scale);
template void Matrix3x3::scalez<float>(float scale);
template void Matrix3x3::loadRotationX<float>(float theta);
template void Matrix3x3::loadRotationY<float>(float theta);
template void Matrix3x3::loadRotationZ<float>(float theta);
template void Matrix3x3::rotateAroundAxis<float>(Vec3D axis, float theta);
template void Matrix3x3::loadEulerTransform<float>(float omega, float theta, float phi);

template Matrix3x3::Matrix3x3<double>(const MatrixNxN<3, double> &other);
template Matrix3x3::Matrix3x3<double>(const Vec<3, double> &vec1, const Vec<3, double> &vec2, const Vec<3, double> &vec3);
template Matrix3x3 &Matrix3x3::operator=<double>(const MatrixNxN<3, double> &other);
template void Matrix3x3::scale<double>(double scalex, double scaley, double scalez);
template void Matrix3x3::scalex<double>(double scale);
template void Matrix3x3::scaley<double>(double scale);
template void Matrix3x3::scalez<double>(double scale);
template void Matrix3x3::loadRotationX<double>(double theta);
template void Matrix3x3::loadRotationY<double>(double theta);
template void Matrix3x3::loadRotationZ<double>(double theta);
template void Matrix3x3::rotateAroundAxis<double>(Vec3D axis, double theta);
template void Matrix3x3::loadEulerTransform<double>(double omega, double theta, double phi);

template Matrix3x3::Matrix3x3<long double>(const MatrixNxN<3, long double> &other);
template Matrix3x3::Matrix3x3<long double>(const Vec<3, long double> &vec1, const Vec<3, long double> &vec2, const Vec<3, long double> &vec3);
template Matrix3x3 &Matrix3x3::operator=<long double>(const MatrixNxN<3, long double> &other);
template void Matrix3x3::scale<long double>(long double scalex, long double scaley, long double scalez);
template void Matrix3x3::scalex<long double>(long double scale);
template void Matrix3x3::scaley<long double>(long double scale);
template void Matrix3x3::scalez<long double>(long double scale);
template void Matrix3x3::loadRotationX<long double>(long double theta);
template void Matrix3x3::loadRotationY<long double>(long double theta);
template void Matrix3x3::loadRotationZ<long double>(long double theta);
template void Matrix3x3::rotateAroundAxis<long double>(Vec3D axis, long double theta);
template void Matrix3x3::loadEulerTransform<long double>(long double omega, long double theta, long double phi);



template Matrix4x4::Matrix4x4<float>(const MatrixNxN<4, float> &other);
template Matrix4x4 &Matrix4x4::operator=<float>(const MatrixNxN<4, float> &other);
template void Matrix4x4::scalex<float>(float scale);
template void Matrix4x4::scaley<float>(float scale);
template void Matrix4x4::scalez<float>(float scale);
template void Matrix4x4::scalew<float>(float scale);
template void Matrix4x4::scale<float>(float scalex, float scaley, float scalez, float scalew);
template void Matrix4x4::rotatex<float>(float phi);
template void Matrix4x4::rotatey<float>(float phi);
template void Matrix4x4::rotatez<float>(float phi);
template void Matrix4x4::translate<float>(float x, float y, float z);

template Matrix4x4::Matrix4x4<double>(const MatrixNxN<4, double> &other);
template Matrix4x4 &Matrix4x4::operator=<double>(const MatrixNxN<4, double> &other);
template void Matrix4x4::scalex<double>(double scale);
template void Matrix4x4::scaley<double>(double scale);
template void Matrix4x4::scalew<double>(double scale);
template void Matrix4x4::scale<double>(double scalex, double scaley, double scalez, double scalew);
template void Matrix4x4::rotatex<double>(double phi);
template void Matrix4x4::rotatey<double>(double phi);
template void Matrix4x4::rotatez<double>(double phi);
template void Matrix4x4::translate<double>(double x, double y, double z);

template Matrix4x4::Matrix4x4<long double>(const MatrixNxN<4, long double> &other);
template Matrix4x4 &Matrix4x4::operator=<long double>(const MatrixNxN<4, long double> &other);
template void Matrix4x4::scalex<long double>(long double scale);
template void Matrix4x4::scaley<long double>(long double scale);
template void Matrix4x4::scalew<long double>(long double scale);
template void Matrix4x4::scale<long double>(long double scalex, long double scaley, long double scalez, long double scalew);
template void Matrix4x4::rotatex<long double>(long double phi);
template void Matrix4x4::rotatey<long double>(long double phi);
template void Matrix4x4::rotatez<long double>(long double phi);
template void Matrix4x4::translate<long double>(long double x, long double y, long double z);

template Matrix5x5::Matrix5x5<float>(const MatrixNxN<5, float> &other);
template Matrix5x5 &Matrix5x5::operator=<float>(const MatrixNxN<5, float> &other);
template void Matrix5x5::rotateXY<float>(float phi);
template void Matrix5x5::rotateXZ<float>(float phi);
template void Matrix5x5::rotateYZ<float>(float phi);
template void Matrix5x5::rotateYW<float>(float phi);
template void Matrix5x5::rotateXW<float>(float phi);
template void Matrix5x5::rotateZW<float>(float phi);

template Matrix5x5::Matrix5x5<double>(const MatrixNxN<5, double> &other);
template Matrix5x5 &Matrix5x5::operator=<double>(const MatrixNxN<5, double> &other);
template void Matrix5x5::rotateXY<double>(double phi);
template void Matrix5x5::rotateXZ<double>(double phi);
template void Matrix5x5::rotateYZ<double>(double phi);
template void Matrix5x5::rotateYW<double>(double phi);
template void Matrix5x5::rotateXW<double>(double phi);
template void Matrix5x5::rotateZW<double>(double phi);
