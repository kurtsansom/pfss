#ifndef vectors_h
#define vectors_h

#include "stdio.h"

#include <cmath>
#include <fstream>

#include "engine/hcTime.h"

class EllipticalGrid;

typedef unsigned int uint;

#ifndef num_eps
#define num_eps 1E-7
#endif

using namespace std;

/*! \brief implementation of the mathematical (finite dimensionality) vector concept
 *
 */
template<uint dim, class T>
class Vec{
public:

      T content[dim];						/*!< \brief elements of the vector							*/
#ifdef __NVCC__
      __host__ __device__
#endif
      Vec();								/*!< \brief std constructor									*/

template<class S>
#ifdef __NVCC__
      __host__ __device__
#endif
      Vec(const Vec<dim, S>& other);		/*!< \brief cpy constructor									*/

#ifdef __NVCC__
      __host__ __device__
#endif
      ~Vec();								/*!< \brief destructor										*/

#ifdef __NVCC__
      __host__ __device__
#endif
      T& operator()(uint n);

#ifdef __NVCC__
      __host__ __device__
#endif
      T operator[](uint n) const;

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
	Vec<dim, T> &operator=(const Vec<dim, S> &other);		/*!< \brief assignment operator												*/

#ifdef __NVCC__
      __host__ __device__
#endif
      Vec &operator*=(T scale);

#ifdef __NVCC__
      __host__ __device__
#endif
      Vec &operator/=(T scale);

#ifdef __NVCC__
      __host__ __device__
#endif
      Vec &operator+=(const Vec &other);

#ifdef __NVCC__
      __host__ __device__
#endif
      Vec &operator-=(const Vec &other);

#ifdef __NVCC__
      __host__ __device__
#endif
      Vec operator-();

#ifdef __NVCC__
      __host__ __device__
#endif
      bool operator==(const Vec &vec) const;

#ifdef __NVCC__
      __host__ __device__
#endif
      bool isAlmostEqual(const Vec<dim, T> &other, T eps=1E-6);
	/*!< \brief checks if vectors are the same within numerical uncertainty bounds														*/

#ifdef __NVCC__
      __host__ __device__
#endif
      Vec<dim, T> &normalize();

#ifdef __NVCC__
      __host__ __device__
#endif
      void loadZeroes();

#ifdef __NVCC__
      __host__ __device__
#endif
      void loadHom(const Vec<dim-1, T> &vec);

#ifdef __NVCC__
      __host__ __device__
#endif
      T sp(const Vec &v) const;		           /*!< \brief std-scalar product with other vector                        					*/

#ifdef __NVCC__
      __host__ __device__
#endif
      double sp_double(const Vec &other) const;

#ifdef __NVCC__
      __host__ __device__
#endif
      T dist(const Vec &other);         		/*!< \brief distance to point in std-norm 												*/

#ifdef __NVCC__
      __host__ __device__
#endif
      T length() const;                 		/*!< \brief distance to origin in std-norm / length of vector 							*/

#ifdef __NVCC__
      __host__ __device__
#endif
      void scale(T factor);             		/*!< \brief scales the vector by factor 												*/

#ifdef __NVCC__
      __host__ __device__
#endif
      void zero();                          	/*!< \brief sets all components to 0 													*/

      /*
#ifdef __NVCC__
      __host__ __device__
#endif
      void forceCopy(const Vec &v);//*/         	/*!< copies components from Vector and overwrites the dimensionality of this			*/

#ifdef __NVCC__
      __host__ __device__
#endif
      bool isNullVector() const;

#ifdef __NVCC__
      __host__ __device__
#endif
      bool isValid() const;                 /*!< \brief checks if some element is INF or NAN */

      bool exportBinary(std::ofstream &stream);

      bool importBinary(std::ifstream &stream, uint sizeofFloat=0);

#ifdef __NVCC__
      __host__ __device__
#endif
      void dump() const;
};

template<uint dim, class T>
#ifdef __NVCC__
      __host__ __device__
#endif
Vec<dim, T> operator+(Vec<dim, T> lhs, const Vec<dim, T> &rhs);

template<uint dim, class T>
#ifdef __NVCC__
      __host__ __device__
#endif
Vec<dim, T> operator-(Vec<dim, T> lhs, const Vec<dim, T> &rhs);

template<uint dim, class T>
#ifdef __NVCC__
      __host__ __device__
#endif
Vec<dim, T> operator*(Vec<dim, T> lhs, double scale);

template<uint dim, class T>
#ifdef __NVCC__
      __host__ __device__
#endif
Vec<dim, T> operator*(double scale, Vec<dim, T> lhs);

template<uint dim, class T>
#ifdef __NVCC__
      __host__ __device__
#endif
Vec<dim, T> operator/(Vec<dim, T> lhs, double scale);

template<uint dim, class T>
#ifdef __NVCC__
      __host__ __device__
#endif
T operator*(const Vec<dim, T> &lhs, const Vec<dim, T> &rhs);

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Vec
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

template<uint dim, class T>
Vec<dim, T>::Vec()
{
	for (uint i=0; i<dim; ++i)
		content[i] = 0.0f;
}

template<uint dim, class T>
template<class S>
Vec<dim, T>::Vec(const Vec<dim, S>& other)
{
  for(uint i=0;i<dim;++i)
	  content[i] = other.content[i];
}

template<uint dim, class T>
Vec<dim, T>::~Vec(){}

template<uint dim, class T>
T &Vec<dim, T>::operator()(uint i)
{
	if(i < dim)
		return(content[i]);
	else
	{
		printf("ERROR: Vec::operator(): Index for vector access out of bounds (i: %u, dim: %u)\n", i, dim);
		return(content[0]);
	}
}

template<uint dim, class T>
T Vec<dim, T>::operator[](uint i) const
{
	if(i < dim)
		return(content[i]);
	else
	{
		printf("ERROR: Vec::operator[]: Index for vector access out of bounds (i: %u, dim: %u)\n", i, dim);
		return(content[0]);
	}
}

template<uint dim, class T> template<class S>
Vec<dim, T> &Vec<dim, T>::operator=(const Vec<dim, S> &other)
{

	for(uint i=0; i<dim; ++i)
		content[i] = other.content[i];

	return *this;
}

template<uint dim, class T>
Vec<dim, T> &Vec<dim, T>::operator+=(const Vec<dim, T> &other)
{
	for(uint i=0; i<dim; ++i)
		content[i] 	+= other.content[i];

	return *this;
}

template<uint dim, class T>
Vec<dim, T> &Vec<dim, T>::operator-=(const Vec<dim, T> &other)
{
	for(uint i=0; i<dim; ++i)
		content[i] 	-= other.content[i];

	return *this;
}

template<uint dim, class T>
Vec<dim, T> Vec<dim, T>::operator-()
{
  Vec<dim, T> retval(*this);

  for(uint i=0; i<dim; ++i)
	  retval(i) = -retval[i];

  return retval;
}

template<uint dim, class T>
Vec<dim, T> &Vec<dim, T>::operator*=(T scale)
{
	for(uint i=0; i<dim; ++i)
		content[i] 	*= scale;

	return *this;
}

template<uint dim, class T>
Vec<dim, T> &Vec<dim, T>::operator/=(T scale)
{
	for(uint i=0; i<dim; ++i)
		content[i] 	/= scale;

	return *this;
}

template<uint dim, class T>
Vec<dim, T> operator+(Vec<dim, T> lhs, const Vec<dim, T> &rhs)
{
	lhs += rhs;
	return lhs;
}

template<uint dim, class T>
Vec<dim, T> operator-(Vec<dim, T> lhs, const Vec<dim, T> &rhs)
{
	lhs -= rhs;
	return lhs;
}

template<uint dim, class T>
Vec<dim, T> operator*(Vec<dim, T> lhs, double scale)
{
	lhs *= scale;
	return lhs;
}

template<uint dim, class T>
Vec<dim, T> operator*(double scale, Vec<dim, T> lhs)
{
	lhs *= scale;
	return lhs;
}

template<uint dim, class T>
Vec<dim, T> operator/(Vec<dim, T> lhs, double scale)
{
	lhs /= scale;
	return lhs;
}

template<uint dim, class T>
T operator*(const Vec<dim, T> &lhs, const Vec<dim, T> &rhs)
{
	return lhs.sp(rhs);
}

template<uint dim, class T>
bool Vec<dim, T>::operator==(const Vec<dim, T> &other) const
{
  for(uint i=0; i<dim; ++i)
	  if(content[i] != other[i])
		  return false;

  return true;
}

template<uint dim, class T>
void Vec<dim, T>::loadHom(const Vec<dim-1, T> &vec)
{
	for(uint i=0;i<dim-1; ++i)
		content[i] = vec[i];

	this->content[dim-1] = 1.0;
}

template<uint dim, class T>
bool Vec<dim, T>::isAlmostEqual(const Vec<dim, T> &other, T eps)
{
	for(uint i=0; i<dim; ++i)
		if(fabs(content[i] - other[i]) > eps)
			return false;

	return true;
}

template<uint dim, class T>
void Vec<dim, T>::dump() const
{
	printf("Dumping hcVec (%u-dimensional):\n", dim);

	for(uint i = 0; i<dim; ++i)
		printf("%u: %E\n", i, content[i]);
	printf("\n");
}

template<uint dim, class T>
Vec<dim, T> &Vec<dim, T>::normalize()
{
	T length = this->length();

	//if(length < num_eps)
	if(length < 1E-10)
	{
		printf("ERROR! Vec::normalize: length of vector is 0!\n");
		return *this;
	}

	for(uint i=0;i<dim;++i)
		this->content[i] = this->content[i] / length;

	return *this;
}

template<uint dim, class T>
void Vec<dim, T>::loadZeroes()
{
  for(uint i=0; i<dim; ++i)
	  content[i] = 0.0;
}

template<uint dim, class T>
T Vec<dim, T>::sp(const Vec<dim, T> &v) const
{
	T result = 0.0;
	for(uint i=0; i<dim; ++i)
		result += content[i] * v.content[i];
	return result;

}

template<uint dim, class T>
double Vec<dim, T>::sp_double(const Vec<dim, T> &other) const
{
	double result = 0.0;

	for(uint i=0; i<dim; ++i)
		result += (double)content[i] * (double)other.content[i];

	return result;
}

template<uint dim, class T>
T Vec<dim, T>::dist(const Vec<dim, T> &other)
{
	T result = 0.0f;

	for(uint i=0; i<dim; ++i)
		result += (other[i] - content[i]) * (other[i] - content[i]);

	return sqrt(result);
}

template<uint dim, class T>
T Vec<dim, T>::length() const
{
	T temp = 0.0f;
	for(uint i=0; i<dim; ++i)
		temp += content[i] * content[i];
	return sqrt(temp);
}

template<uint dim, class T>
void Vec<dim, T>::scale(T factor)
{
	for(uint i=0; i<dim; ++i)
		content[i] *= factor;
	return;
}

template<uint dim, class T>
void Vec<dim, T>::zero()
{
	for(uint i=0; i<dim; ++i)
		content[i] = 0;
}

/*
template<uint dim, class T>
void Vec<dim, T>::forceCopy(const Vec<dim, T> &v){

	for(uint i=0;i<dim;++i)
		content[i] = v.content[i];
}//*/

template<uint dim, class T>
bool Vec<dim, T>::isNullVector() const
{
	for(uint i=0; i<dim; ++i)
		if(fabs(content[i]) > num_eps)
			return 0;
	return 1;
}

template<uint dim, class T>
bool Vec<dim, T>::isValid() const
{
	for(uint i=0; i<dim; ++i)
		if(isnan((float)content[i]) || INFINITY == content[i] || -INFINITY == content[i])
			return false;

	return true;
}

template<uint dim, class T>
bool Vec<dim, T>::exportBinary(std::ofstream &stream)
{
	for(uint i=0; i<dim; ++i)
		stream.write(reinterpret_cast<char*>(&content[i]),sizeof(hcFloat));

	return true;
}

template<uint dim, class T>
bool Vec<dim, T>::importBinary(std::ifstream &stream, uint sizeofFloat)
{
	if(sizeofFloat==sizeof(hcFloat))
	{
		for(uint i=0; i<dim; ++i)
			stream.read(reinterpret_cast<char*>(&content[i]),sizeof(hcFloat));
	}
	else
	{
		cerr << __FILE__ << ":" << __LINE__ << " sizeofFloat (" << sizeofFloat << ") != sizeof(hcFloat) (" << sizeof(hcFloat) << ") not supported.\n";
		return false;
	}
	return true;
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

/*! \brief 3D Vectors
 *
 */
class Vec2D : public Vec<2, hcFloat>{
public:
  
#ifdef __NVCC__
	__host__ __device__
#endif
    Vec2D();									/*!< \brief std constructor							*/

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    Vec2D(S x, S y);

#ifdef __NVCC__
	__host__ __device__
#endif
    Vec2D(const Vec2D &other);					/*!< \brief cpy constructor							*/

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    Vec2D(const Vec<2, S> &other);		/*!< \brief pseudo-cpy constructor					*/

#ifdef __NVCC__
	__host__ __device__
#endif
    ~Vec2D();									/*!< \brief destructor								*/

#ifdef __NVCC__
	__host__ __device__
#endif
    Vec2D &operator=(const Vec2D &other);

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    Vec2D &operator=(const Vec<2, S> &other);

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    void rotate(S phi);

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    void set(S x, S y);

//virtual bool importBinary(std::ofstream &stream);
};


//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Vec3D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

/*! \brief 3D Vectors
 *
 */
class Vec3D : public Vec<3, hcFloat>{
public:

#ifdef __NVCC__
	__host__ __device__
#endif
    Vec3D();		            			/*!< \brief std constructor                                                     	*/

#ifdef __NVCC__
	__host__ __device__
#endif
    Vec3D(const Vec3D &other);  			/*!< \brief cpy constructor                                                     	*/

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    Vec3D(const Vec<3, S> &other);  		/*!< \brief pseudo-cpy constructor                                                  */

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    Vec3D(S x, S y, S z);

#ifdef __NVCC__
	__host__ __device__
#endif
    ~Vec3D();		            			/*!< \brief destructor			                                                   	*/

#ifdef __NVCC__
	__host__ __device__
#endif
    Vec3D &operator=(const Vec3D &other);

	template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    Vec3D &operator=(const Vec<3, S> &other);

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    void rotateX(S angle);
        /*!< \brief  rotate around x-axis in global cartesian coordinate system */

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    void rotateY(S angle);
        /*!< \brief  rotate around y-axis in global cartesian coordinate system */

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    void rotateZ(S angle);
        /*!< \brief  rotate around z-axis in global cartesian coordinate system */

#ifdef __NVCC__
	__host__ __device__
#endif
    void stereoProjUnitSphere(Vec2D tangent);
        /*!< \brief stereographic projection of this (in cartCoords ) onto plane tangent at unit sphere in spherical coordinates (Vec2D tangent) */

#ifdef __NVCC__
	__host__ __device__
#endif
	void convertEllipticalCoordsToCartesian(hcFloat a, bool prolate);
	/*!< \brief transforms elliptical coordinates to cartesian coordinates                                   						*/

#ifdef __NVCC__
	__host__ __device__
#endif
	void convertCartesianCoordsToElliptical(hcFloat a, bool prolate);
	/*!< \brief transforms cartesian coordinates to elliptical coordinates                                   						*/

// --------------------------------------------------------------------------------------------------------
// conversion functions
// --------------------------------------------------------------------------------------------------------
#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convCoordSpher2Cart() const;

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convCoordCart2Spher() const;

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convCoordCart2Ell(const EllipticalGrid &grid) const;

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convCoordEll2Cart(const EllipticalGrid &grid) const;

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convCoordSpher2Ell(EllipticalGrid &grid) const;

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convCoordEll2Spher(EllipticalGrid &grid) const;

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convVecSpher2Cart(const Vec3D &cartPos) const;

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convVecCart2Spher(const Vec3D &cartPos);

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convCoordCart2Cart(Vec3D e1x, Vec3D e1y, Vec3D e1z, Vec3D e2x, Vec3D e2y, Vec3D e2z) const;
	/*!< \brief convert position vector from cartesian coordinate system 1 to cartesian coordinate system 2								*/

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convVecSpher2Ell(const Vec3D &posSpher, const EllipticalGrid &grid) const;

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convVecEll2Spher(const Vec3D &posEll, const EllipticalGrid &grid) const;

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convVecEll2Cart(const Vec3D &cartPos, const EllipticalGrid &grid) const;

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convVecCart2Ell(const Vec3D &cartPos, const EllipticalGrid &grid) const;

#ifdef __NVCC__
	__host__ __device__
#endif
Vec3D convVecCart2Ell2(const Vec3D &cartPos, EllipticalGrid &grid) const;

Vec3D convVecHAE2GSE(const hcDate &date);
	/*!< \brief convert vecter in Heliocentric Aries Ecliptic cartesian coordinates to Geocentric Solar Ecliptic cartesian coordsinates	*/

Vec3D convVecGSE2HAE(const hcDate &date);
	/*!< \brief convert vecter in Geocentric Solar Ecliptic cartesian coordsinates to Heliocentric Aries Ecliptic cartesian coordinates	*/

Vec3D convVecHAE2SW(const Vec3D &pos_HAE_c, hcFloat swSpeed, hcFloat lowerR, const hcDate &date) const;
	/*!< \brief convert vector from Heliocentric Aries Ecliptic cartesian coordinates to solar wind frame								*/


// --------------------------------------------------------------------------------------------------------
// conversion functions end
// --------------------------------------------------------------------------------------------------------

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    void transformSphericalCoordSystem(S theta, S phi);
        /*!< \brief transforms point given in spherical coordinates into another spherical system               */

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    void transformSphericalCoordSystemBack(S theta, S phi);
        /*!< \brief transform back point given in spherical coordinates (see transformSphericalCoordSystem)     */

#ifdef __NVCC__
	__host__ __device__
#endif
    void cp(const Vec3D &v, Vec3D *res) const;
        /*!< \brief cross product of two Vec3D                                                                  */

#ifdef __NVCC__
	__host__ __device__
#endif
    Vec3D cp(const Vec3D &other) const;
        /*!< \brief cross product of two Vec3D                                                                  */

#ifdef __NVCC__
	__host__ __device__
#endif
    hcFloat getAngle(const Vec3D &other) const;

#ifdef __NVCC__
	__host__ __device__
#endif
    hcFloat getAngle2(const Vec3D &other) const;

#ifdef __NVCC__
	__host__ __device__
#endif
    void dumpSphericalCoords();

	//virtual bool importBinary(std::ofstream &stream);
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Vec4D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
/*! \brief 4D Vectors
 *
 */
class Vec4D : public Vec<4, hcFloat>{
public:

#ifdef __NVCC__
	__host__ __device__
#endif
    Vec4D();				                    /*!< \brief std constructor     			*/

#ifdef __NVCC__
	__host__ __device__
#endif
    Vec4D(const Vec4D &other);                  /*!< \brief cpy constructor    				*/

#ifdef __NVCC__
	__host__ __device__
#endif
    Vec4D(const Vec3D &other);                  /*!< \brief homogen-cpy constructor    		*/

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    Vec4D(const Vec<4, S> &other);        		/*!< \brief pseudo-cpy constructor    		*/

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    Vec4D(S x, S y, S z, S w);

#ifdef __NVCC__
	__host__ __device__
#endif
    ~Vec4D();				                    /*!< \brief destructor     					*/

#ifdef __NVCC__
	__host__ __device__
#endif
    Vec4D &operator=(const Vec4D &other);

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    Vec4D &operator=(const Vec<4, S> &other);

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    void set(S x, S y, S z, S w);

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    void setPos(S x, S y);

template<class S>
#ifdef __NVCC__
	__host__ __device__
#endif
    void setTex(S u, S v);

	//virtual bool importBinary(std::ofstream &stream);
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Vec5D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
/*! \brief 5D Vectors
 *
 */
class Vec5D : public Vec<5, hcFloat>{
public:

    Vec5D();

    template<class S>
    Vec5D(S x, S y, S z, S u, S v);
    ~Vec5D();

    template<class S>
    void set(S x, S y, S z, S u, S v);

    //virtual bool importBinary(std::ofstream &stream);
};


template class Vec<2, float>;

template Vec<2, float> operator+<2, float>(Vec<2, float> lhs, const Vec<2, float> &rhs);
template Vec<2, float> operator-<2, float>(Vec<2, float> lhs, const Vec<2, float> &rhs);
template Vec<2, float> operator*<2, float>(Vec<2, float> lhs, double scale);
template Vec<2, float> operator*<2,float>(double scale, Vec<2, float> lhs);
template Vec<2, float> operator/<2, float>(Vec<2, float> lhs, double scale);
template float operator*<2, float>(const Vec<2, float> &lhs, const Vec<2, float> &rhs);

template Vec<2, double> operator+<2, double>(Vec<2, double> lhs, const Vec<2, double> &rhs);
template Vec<2, double> operator-<2, double>(Vec<2, double> lhs, const Vec<2, double> &rhs);
template Vec<2, double> operator*<2, double>(Vec<2, double> lhs, double scale);
template Vec<2, double> operator*<2,double>(double scale, Vec<2, double> lhs);
template Vec<2, double> operator/<2, double>(Vec<2, double> lhs, double scale);
template double operator*<2, double>(const Vec<2, double> &lhs, const Vec<2, double> &rhs);


template Vec<3, float> operator+<3, float>(Vec<3, float> lhs, const Vec<3, float> &rhs);
template Vec<3, float> operator-<3, float>(Vec<3, float> lhs, const Vec<3, float> &rhs);
template Vec<3, float> operator*<3, float>(Vec<3, float> lhs, double scale);
template Vec<3, float> operator*<3,float>(double scale, Vec<3, float> lhs);
template Vec<3, float> operator/<3, float>(Vec<3, float> lhs, double scale);
template float operator*<3, float>(const Vec<3, float> &lhs, const Vec<3, float> &rhs);

template Vec<3, double> operator+<3, double>(Vec<3, double> lhs, const Vec<3, double> &rhs);
template Vec<3, double> operator-<3, double>(Vec<3, double> lhs, const Vec<3, double> &rhs);
template Vec<3, double> operator*<3, double>(Vec<3, double> lhs, double scale);
template Vec<3, double> operator*<3,double>(double scale, Vec<3, double> lhs);
template Vec<3, double> operator/<3, double>(Vec<3, double> lhs, double scale);
template double operator*<3, double>(const Vec<3, double> &lhs, const Vec<3, double> &rhs);


template Vec<4, float> operator+<4, float>(Vec<4, float> lhs, const Vec<4, float> &rhs);
template Vec<4, float> operator-<4, float>(Vec<4, float> lhs, const Vec<4, float> &rhs);
template Vec<4, float> operator*<4, float>(Vec<4, float> lhs, double scale);
template Vec<4, float> operator*<4,float>(double scale, Vec<4, float> lhs);
template Vec<4, float> operator/<4, float>(Vec<4, float> lhs, double scale);
template float operator*<4, float>(const Vec<4, float> &lhs, const Vec<4, float> &rhs);

template Vec<4, double> operator+<4, double>(Vec<4, double> lhs, const Vec<4, double> &rhs);
template Vec<4, double> operator-<4, double>(Vec<4, double> lhs, const Vec<4, double> &rhs);
template Vec<4, double> operator*<4, double>(Vec<4, double> lhs, double scale);
template Vec<4, double> operator*<4,double>(double scale, Vec<4, double> lhs);
template Vec<4, double> operator/<4, double>(Vec<4, double> lhs, double scale);
template double operator*<4, double>(const Vec<4, double> &lhs, const Vec<4, double> &rhs);


template Vec<5, float> operator+<5, float>(Vec<5, float> lhs, const Vec<5, float> &rhs);
template Vec<5, float> operator-<5, float>(Vec<5, float> lhs, const Vec<5, float> &rhs);
template Vec<5, float> operator*<5, float>(Vec<5, float> lhs, double scale);
template Vec<5, float> operator*<5,float>(double scale, Vec<5, float> lhs);
template Vec<5, float> operator/<5, float>(Vec<5, float> lhs, double scale);
template float operator*<5, float>(const Vec<5, float> &lhs, const Vec<5, float> &rhs);

template Vec<5, double> operator+<5, double>(Vec<5, double> lhs, const Vec<5, double> &rhs);
template Vec<5, double> operator-<5, double>(Vec<5, double> lhs, const Vec<5, double> &rhs);
template Vec<5, double> operator*<5, double>(Vec<5, double> lhs, double scale);
template Vec<5, double> operator*<5,double>(double scale, Vec<5, double> lhs);
template Vec<5, double> operator/<5, double>(Vec<5, double> lhs, double scale);
template double operator*<5, double>(const Vec<5, double> &lhs, const Vec<5, double> &rhs);

#ifdef __NVCC__
	__host__ __device__
#endif
	void getSphericGradientVectors(const Vec3D &cartPos, Vec3D &er, Vec3D &et, Vec3D &ep);
		/*!< \brief returns contravariant spherical basis vectors												*/


#endif
