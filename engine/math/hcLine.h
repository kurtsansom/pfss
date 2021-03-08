#ifndef lines_h
#define lines_h

#include "engine/math/hcVec.h"
#include "engine/math/hcMatrix.h"
#include "engine/math/hcPlane.h"

#include <stdlib.h>
#include <stdio.h>

/*! \brief contains mathematical model and functions for nD lines
 *
 *  The nD-Line is parametrized by by two nD-vectors
 *  (position and direction vectors): line = pos + lambda * dir
 *  lambda is element of [0,1] if line is initialized by giving two points
 *  to be connected by the line
 */
template<uint dim>
class hcLine{
public:

    Vec<dim, hcFloat> pos;
    Vec<dim, hcFloat> dir;

#ifdef __NVCC__
    __host__ __device__
#endif
    hcLine();													/*!< \brief std-constructor						*/

#ifdef __NVCC__
    __host__ __device__
#endif
    hcLine(const hcLine<dim> &other);							/*!< \brief cpy-constructor						*/

#ifdef __NVCC__
    __host__ __device__
#endif
    hcLine(const Vec<dim, hcFloat> &pos, const Vec<dim, hcFloat> &dir);			/*!< \brief direction constructor				*/

#ifdef __NVCC__
      __host__ __device__
#endif
    ~hcLine();													/*!< \brief destructor							*/

#ifdef __NVCC__
    __host__ __device__
#endif
    hcLine<dim> &operator=(const hcLine<dim> &other);			/*!< \brief assignment operator					*/
 
#ifdef __NVCC__
    __host__ __device__
#endif
    int intersectsL(const hcLine<dim> &l, Vec<dim, hcFloat> &result) const;
    	/*!< \brief tests the intersectio of this line with another line										*/

#ifdef __NVCC__
    __host__ __device__
#endif
    int intersectsP(const Vec<dim, hcFloat> &vec, float &result) const;
    	/*!< \brief tests if a point lies on this line															*/

#ifdef __NVCC__
    __host__ __device__
#endif
    int createLineThroughPoints(const Vec<dim, hcFloat> &point1, const Vec<dim, hcFloat> &point2);
    	/*!< \brief creates a line through point1 and point2													*/

#ifdef __NVCC__
    __host__ __device__
#endif
    float pointProjection(const Vec<dim, hcFloat> &point) const;

#ifdef __NVCC__
    __host__ __device__
#endif
    void set(const Vec<dim, hcFloat> &vec1, const Vec<dim, hcFloat> &vec2);

#ifdef __NVCC__
    __host__ __device__
#endif
    int setV1(const Vec<dim, hcFloat> &vec);

#ifdef __NVCC__
    __host__ __device__
#endif
    int setV2(const Vec<dim, hcFloat> &vec);

#ifdef __NVCC__
    __host__ __device__
#endif
    Vec<dim, hcFloat> getPos(hcFloat lambda) const;

#ifdef __NVCC__
    __host__ __device__
#endif
    void dump() const;
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Line
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

template<uint dim>
hcLine<dim>::hcLine()
{
	printf("hcLine<uint dim> std constructor should not be used!\n");
}

template<uint dim>
hcLine<dim>::hcLine(const Vec<dim, hcFloat> &v1, const Vec<dim, hcFloat> &v2) :
	pos(v1), dir(v2)
{

	/*
	for(uint i=0; i<dim; ++i)
	{
		this->pos.content[i] = v1.content[i];
		this->dir.content[i] = v2.content[i];
	}//*/
}

template<uint dim>
hcLine<dim>::hcLine(const hcLine<dim> &other)
{
    pos = other.pos;
    dir = other.dir;
}

template<uint dim>
hcLine<dim>::~hcLine(){}

template<uint dim>
void hcLine<dim>::set(const Vec<dim, hcFloat> &vec1, const Vec<dim, hcFloat> &vec2)
{
	pos = vec1;
	dir = vec2;
}

/*! \brief checks if the nD-point denoted by vec is on this Line
 *
 *  \retval 0 point not on Line
 *  \retval 1 point on Line, line parameter is returned in result
 *
 */
template<uint dim>
int hcLine<dim>::intersectsP(const Vec<dim, hcFloat> &vec, float &result) const
{
    float eps = 1E-7;

	result = 0.0;
	for(uint i=0; i<dim; ++i)
		if(fabs(dir.content[i]) > eps)
		{
			result = (vec.content[i] - pos.content[i]) / dir.content[i];
			break;
		}

	for(uint i=0; i<dim; ++i)
	{
		if(fabs(dir.content[i]) < eps)
		{
			if(fabs(vec.content[i] - pos.content[i]) > eps)
				return 0;
		}
		else
			if (fabs(result-(vec.content[i] - pos.content[i]) / dir.content[i]) > eps) // there is no solution in this case
				return 0;
	}

    return 1;
}

template<uint dim>
hcLine<dim> &hcLine<dim>::operator=(const hcLine<dim> &other)
{
	if(this == &other)
		return *this;

	pos = other.pos;
	dir = other.dir;

	return *this;
}

/*! \brief computes the intersection with another nD-line
 *
 * \retval 0 no intersection found
 * \retval 1 exactly one intersection found
 * \retval 2 lines are identical
 *
 *  the Vector result contains the intersection point in case 1
 *  and the position vector of this Line in the case 2
 */
template<uint dim>
int hcLine<dim>::intersectsL(const hcLine<dim> &l, Vec<dim, hcFloat> &result) const
{
	Vec<dim, hcFloat> tv1(l.pos);
	Vec<dim, hcFloat> tv2;
	Vec<dim, hcFloat> res;
	tv2	-= l.dir;
	tv1	-= pos;
	Matrix<dim, 3, hcFloat> m;
	for(uint i=0; i<dim; ++i)
	{
		m(i, 0) = dir[i];
		m(i, 1)	= tv2[i];
		m(i, 2) = tv1[i];
	}
	uint erg = m.solveSLE(res); // TODO: actually, we don't need the entire result-vector, the first entry would suffice -> make more efficient?
	if (erg == 1)
	{
		result 	= getPos(res(0));
		return 1;
	}
	else if (erg == 2)
	{
		for(uint i=0; i<dim; ++i)
			result(i) = pos[i];
		return 2;
	}
	else
		return 0;
}

template<uint dim>
float hcLine<dim>::pointProjection(const Vec<dim, hcFloat> &point) const
{
    Vec<dim, hcFloat> temp 	= point;
    temp			-= pos;
    return(dir.sp(temp) / dir.sp(dir));
}

/*! \brief creates ND-line through two ND-points, if they are not too close to each other
 *
 *  \retval 0 points too close or dimension mismatch
 *  \retval 1 line successfully created
 */
template<uint dim>
int hcLine<dim>::createLineThroughPoints(const Vec<dim, hcFloat> &point1, const Vec<dim, hcFloat> &point2)
{
    Vec<dim, hcFloat> direction(point2);
    direction.minus(point1);

    if(direction.length() < num_eps)
    {
        printf("ERROR! Line::createLineThroughPoints: Points are (nearly) identical!\n");
        point1.dump();
        point2.dump();
        return 0;
    }
    pos	= point1;
    dir = direction;
    return 1;
}

/*! \brief computes the nD-Vetor obtained by inserting lambda in this Line equation
 */
template<uint dim>
Vec<dim, hcFloat> hcLine<dim>::getPos(hcFloat lambda) const
{
	Vec<dim, hcFloat> retval;

    for(uint i=0; i<dim; ++i)
        retval(i) = pos[i] + lambda * dir[i];

    return retval;
}

template<uint dim>
int hcLine<dim>::setV1(const Vec<dim, hcFloat> &vec){

	pos = vec;
	return 1;
}

template<uint dim>
int hcLine<dim>::setV2(const Vec<dim, hcFloat> &vec)
{
	dir = vec;
	return 1;
}

template<uint dim>
void hcLine<dim>::dump() const
{
    printf("Dumping hcLine (dim=%u):\n", dim);
    printf("pos: ");
    pos.dump();
    printf("dir: ");
    dir.dump();
}


template<uint dim> class hcSphere;//*/

class hcLine2D : public hcLine<2>{
public:

    hcLine2D(const Vec2D &pos = Vec2D(0.0, 0.0),
             const Vec2D &dir = Vec2D(1.0, 1.0));           /*!< \brief std constructor         */
    hcLine2D(const hcLine2D &other);                        /*!< \brief cpy constructor         */
    hcLine2D(float x1, float y1, float x2, float y2);

    hcLine2D &operator=(const hcLine2D &other);
};

class hcLine3D : public hcLine<3>{
public:

#ifdef __NVCC__
    __host__ __device__
#endif
    hcLine3D(const Vec3D &pos0 = Vec3D(0.0, 0.0, 0.0),      /*!< \brief std constructor         */
             const Vec3D &pos1 = Vec3D(1.0, 1.0, 1.0));

#ifdef __NVCC__
    __host__ __device__
#endif
    hcLine3D(const hcLine3D &other);						/*!< \brief cpy constructor			*/

#ifdef __NVCC__
    __host__ __device__
#endif
    ~hcLine3D();											/*!< \brief destructor				*/

#ifdef __NVCC__
    __host__ __device__
#endif
    hcLine3D &operator=(const hcLine3D &other);				/*!< \brief assignment operator		*/

#ifdef __NVCC__
    __host__ __device__
#endif
    void init(const Vec3D &pos0, const Vec3D &pos1);

#ifdef __NVCC__
    __host__ __device__
#endif
    int intersectsSphere3D(const hcSphere<3> &sphere, hcFloat &result0, hcFloat &result1, bool getInsidePos);
        /*!< \brief returns the parameter of the line where it intersects the sphere */

#ifdef __NVCC__
    __host__ __device__
#endif
    Vec3D intersectsPlane3D(const hcPlane3D &plane, uint &result);
        /*!< \brief returns the point of intersection */

#ifdef __NVCC__
    __host__ __device__
#endif
    int getIntersectionsWithSphere(const hcSphere<3> &sphere, Vec3D &result0, Vec3D &result1, bool getInsidePos=false);
        /*!< \brief gives the two points where the line pireces the sphere (or the point where it touches the sphere) */
};


#endif
