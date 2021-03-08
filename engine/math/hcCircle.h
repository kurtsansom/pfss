#ifndef circle_h
#define circle_h

#include <stdlib.h>
#include <stdio.h>

#include "engine/math/hcVec.h"
#include "engine/math/hcLine.h"

/*! \brief mathematical model and functions for circles in nD-space
 */
template<uint dim>
class hcCircle{
public:

    Vec<dim, hcFloat> pos;
    float radius;

#ifdef __NVCC__
    __host__ __device__
#endif
    hcCircle();													/*!< \brief std constructor					*/

#ifdef __NVCC__
    __host__ __device__
#endif
    hcCircle(const hcCircle &other);							/*!< \brief cpy constructor					*/

#ifdef __NVCC__
    __host__ __device__
#endif
    hcCircle(const Vec<dim, hcFloat> &pos, hcFloat radius);				/*!< \brief constructor						*/

#ifdef __NVCC__
    __host__ __device__
#endif
    ~hcCircle();												/*!< \brief destructor						*/

#ifdef __NVCC__
    __host__ __device__
#endif
    hcCircle<dim> &operator=(const hcCircle<dim> &other);		/*!< \brief assignment operator				*/

#ifdef __NVCC__
    __host__ __device__
#endif
    void init(const Vec<dim, hcFloat> &pos, hcFloat radius);

#ifdef __NVCC__
    __host__ __device__
#endif
    void dump() const;
};

/***************************************************************************************
 *  Circle
 **************************************************************************************/

template<uint dim>
hcCircle<dim>::hcCircle() : pos(){

    radius = 0.0;
}

template<uint dim>
hcCircle<dim>::hcCircle(const hcCircle &other) : pos(other.pos){

    radius = other.radius;
}

template<uint dim>
hcCircle<dim>::hcCircle(const Vec<dim, hcFloat> &pos, hcFloat radius){

    init(pos, radius);
}

template<uint dim>
hcCircle<dim>::~hcCircle(){}

template<uint dim>
hcCircle<dim> &hcCircle<dim>::operator=(const hcCircle<dim> &other){

	if(this == &other)
		return *this;

	this->radius 	= other.radius;
	this->pos		= other.pos;

	return *this;
}

template<uint dim>
void hcCircle<dim>::init(const Vec<dim, hcFloat> &pos, hcFloat radius){

    this->radius 	= radius;
    this->pos 		= pos;
}

template<uint dim>
void hcCircle<dim>::dump() const{

    printf("Dumping Circle:\n");
    printf("Center:\n");
    pos.dump();
    printf("Radius: %E\n\n", radius);
}


/*! \brief mathematical model and functions for circles in 2D-space
 *
 *  this class is not fit for circles oriented in nD-space
 */
template<uint dim>
class hcCircleFlat{
public:

    Vec<dim, hcFloat> pos;
    float radius;

    hcCircleFlat(	const Vec2D &pos 	= Vec2D(0.0, 0.0),  	/*!< \brief std constructor */
    				float radius     	= 1.0);

    hcCircleFlat(const hcCircleFlat &other);            		/*!< \brief cpy constructor */


    hcCircleFlat &operator=(const hcCircleFlat &other);

    void init(const Vec2D& pos, float radius);
    void set(const Vec2D& pos, float radius);
    int intersectsL2D(const hcLine2D &line, Vec2D &result0, Vec2D &result1) const;
    int intersectsC2D(const hcCircleFlat &circle, Vec2D &result1, Vec2D &result2, hcLine2D &resultline) const;
    int getTangentVecForced(const Vec2D &in, Vec2D &result) const;
};


/***************************************************************************************
 *  CircleFlat
 **************************************************************************************/

template<uint dim>
hcCircleFlat<dim>::hcCircleFlat(const Vec2D &pos, float radius)
{
    init(pos, radius);
}

template<uint dim>
hcCircleFlat<dim>::hcCircleFlat(const hcCircleFlat<dim> &other) :
    pos(other.pos)
{
    this->radius    = other.radius;
}

template<uint dim>
void hcCircleFlat<dim>::init(const Vec2D &pos, float radius){

    this->pos       = pos;
    this->radius    = radius;
}

template<uint dim>
hcCircleFlat<dim> &hcCircleFlat<dim>::operator=(const hcCircleFlat<dim> &other){

    if(this == &other)
        return *this;

    init(other.pos, other.radius);

    return *this;
}

template<uint dim>
void hcCircleFlat<dim>::set(const Vec2D& pos, float radius)
{
    this->pos 		= pos;
    this->radius 	= radius;
}

/*! \brief computes the intersections point(s) with a line for the std-norm
 *
 * \retval 0 no intersection of line with circle
 * \retval 1 two intersections, stored in result1 and result2
 * \retval 2 line tangential, tangential point stored in result1
 */
template<uint dim>
int hcCircleFlat<dim>::intersectsL2D(const hcLine2D &line, Vec2D &result0, Vec2D &result1) const
{
	// the following stuff is obtained by inserting equation for line p = (v1 + t * v2)
	// into cirlce equation (p - pos_circle)^2 = radius^2
	// after fiddeling around a bit, the following variables offer themselves
	Vec<dim, hcFloat> diffvec(line.pos);
	//diffvec.minus(pos);
	diffvec			-= pos;
	float s2 		= line.dir.sp(line.dir);
	float diff2 	= diffvec.sp(diffvec);
	float numerator = 1 / s2;
	float quadterm 	= - diffvec.sp(line.dir) * numerator;
	float test 		= (radius * radius - diff2) * numerator + quadterm * quadterm;

	if(test < 0)                         // no intersection
		return 0;
	else if (test == 0)
	{
		result0 = line.getPos(quadterm);   // tangential point
		return 2;
	}
	else
	{
		float root 		= sqrt(test);      // two intersection points
		float lambda0 	= quadterm + root;
		float lambda1 	= quadterm - root;
		result0			= line.getPos(lambda0);
		result1			= line.getPos(lambda1);
		return 1;
	}
}

/*! computes the intersection of 2 circles in 2D-space (std-norm)
 * \retval 0 no intersection / areas disjunct or one circle in the other
 * \retval 1 two intersection points, stored in result1 and result2, the line produced by these two
 *      points is stored in resultline
 * \retval 2 tangential point, stored in result1
 */
template<uint dim>
int hcCircleFlat<dim>::intersectsC2D(const hcCircleFlat &circle, Vec2D &result0, Vec2D &result1, hcLine2D &resultline) const{

	Vec<dim, hcFloat> diffvec(circle.pos);
	diffvec		-= pos;
	float dist 	= diffvec.length();

	if ((dist < circle.radius + radius) && (dist > abs(radius - circle.radius)))       // main case, intersection with two points
	{
		float inv_dist 	= 1 / dist;
		// position of midpoint of chord connecting the intersection-points
		// relative to this->pos
		float chordpos 	= 0.5 * (dist + inv_dist * (radius * radius - circle.radius * circle.radius));
		diffvec.scale(chordpos / dist);

		resultline.pos 	= pos;
		resultline.pos	+= diffvec;			// v1 now holds the pos.vector to the midpoint of chord

		if (diffvec.content[0] == 0)      	// now get a vector perp. to diffvec as v2
			resultline.dir = Vec2D(1.0, 0.0);
		else
			resultline.dir = Vec2D(- diffvec.content[1] / diffvec.content[0], (hcFloat)1.0);

		this->intersectsL2D(resultline, result0, result1);
		return 1;
	}
	else
	{
		if ((dist == circle.radius + radius) || (dist == abs(radius - circle.radius))) // circles touch each other
		{
			diffvec.scale(radius / dist);
			result0 = pos;
			result0 += diffvec;
			return 2;
		}
		else                                                                             // no intersection
			return 0;
	}
}

/*! \brief computes the tangent vector to a point on the circle
 *
 *  No checks are made during computation. This speeds up the process
 *  and gets rid of numerical artifacts, though it has to be assured
 *  beforehand that the point really lies on the circle. Otherwise
 *  following computations might yield rubbish
 */
template<uint dim>
int hcCircleFlat<dim>::getTangentVecForced(const Vec2D &in, Vec2D &result) const{

    //Vec2D tempvec(in);
    //tempvec.minus(pos);
	Vec2D tempvec 	= in;
	tempvec 		-= pos;

    result(0) 		= -tempvec[1];
    result(1) 		=  tempvec[0];
    return 1;
}

#endif
