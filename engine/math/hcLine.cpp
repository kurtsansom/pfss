#include "engine/math/hcLine.h"
#include "engine/math/hcVec.h"
#include "engine/math/hcMatrix.h"
#include "engine/math/hcSphere.h"

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Line2D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

hcLine2D::hcLine2D(float x1, float y1, float x2, float y2) : hcLine<2>()
{
    pos.content[0] = x1;
    pos.content[1] = y1;
    dir.content[0] = x2;
    dir.content[1] = y2;
}


hcLine2D::hcLine2D(const Vec2D &v1, const Vec2D &v2) : hcLine<2>()
{
    this->pos.content[0] = v1.content[0];
    this->pos.content[1] = v1.content[1];
    this->dir.content[0] = v2.content[0];
    this->dir.content[1] = v2.content[1];
}

hcLine2D::hcLine2D(const hcLine2D &other) : hcLine<2>(other){

    *this = other;
}

hcLine2D &hcLine2D::operator=(const hcLine2D &other){

    if(this == &other)
        return *this;

    hcLine::operator=(other);

    return *this;
}


//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          Line3D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

hcLine3D::hcLine3D(const Vec3D &pos0, const Vec3D &pos1) : hcLine<3>(pos0, pos1-pos0)
{
    init(pos0, pos1);
}

hcLine3D::hcLine3D(const hcLine3D &other) : hcLine<3>(other)
{}

hcLine3D::~hcLine3D(){}

hcLine3D &hcLine3D::operator=(const hcLine3D &other){

    if(this == &other)
        return *this;

    hcLine::operator=(other);

    return *this;
}

void hcLine3D::init(const Vec3D &pos0, const Vec3D &pos1)
{
	this->pos = pos0;
	this->dir = pos1;
	this->dir -= pos0;
}

/*! returns the line parameter lambda0 and lambda1 where the line pierces the sphere
*
* @param sphere the sphere to be pierced
* @param result0 first intersection or tangent parameter of sphere and line
* @param result1 second intersection parameter of sphere and line
* @param getInsidePos numerical parameter, determines if result0/1 shall be slightly above or below sphere surface
*        this parameter affects only the pirecing (retval = 1) and not not the tangent (retval = 2) case
*
* \retval 0 no intersection of line with sphere
* \retval 1 line tangential, tangential point stored in result0
* \retval 2 two intersections, stored in result0 and result1
*/
int hcLine3D::intersectsSphere3D(const hcSphere<3> &sphere, hcFloat &result0, hcFloat &result1, bool getInsidePos){

    const float eps = 1E-6;
    Vec3D d = this->dir;

    float dx    = d[0];
    float dy    = d[1];
    float dz    = d[2];

    Vec3D p0 	= this->pos;
    p0 			-= sphere.pos;

    float px   = p0[0];
    float py   = p0[1];
    float pz   = p0[2];

    float d2    = d.sp(d);
    float r2    = sphere.radius * sphere.radius;
    float p2    = p0.sp(p0);

    // quadratic extension
    float first = - (px * dx + py * dy + pz * dz) / d2;
    float square= (px * dx + py * dy + pz * dz) / d2;
    square      = square * square;
    float root  = sqrt((r2 - p2) / d2 + square);

    result0     = first - root;
    result1     = first + root;

    if(isnan(result0))
    {
        if(isnan(result1))
            return 0;       // there is no real solution to the quadratic equation
        else
        {
            printf("ERROR! Line3D::intersectsSphere3D: this case should not be possible");
            result0 = result1;
            return 1;
        }
    }
    else
    {
        if(fabs(result0 - result1) < eps)   // line tangential to sphere
            return 1;
        else                                // default case with two intersection points
            return 2;
    }
}

/*! returns the point of intersection (if existing, zero-vec otherwise)
*
* @param plane the plane to be pierced
* @param result 0 -> no intersection, 1 -> one intersection, 2 -> line lies in plane
*
*/
Vec3D hcLine3D::intersectsPlane3D(const hcPlane3D &plane, uint &result)
{
	hcFloat sp	= this->dir * plane.normal;

	if(fabs(sp) < 1E-6)		// plane and line are parallel
		if(plane.normal * (this->pos - plane.pos) < 1E-6)	// line lies in plane
		{
			result	= 2;
			return this->pos;
		}
		else
		{
			result	= 0;
			return Vec3D(0.0f, 0.0f, 0.0f);
		}

	Vec3D p1		= plane.pos;
	Vec3D n			= plane.normal;
	Vec3D p2		= this->pos;
	Vec3D d			= this->dir;

	hcFloat x		= (p1*n - p2*n) / (d*n);
	result			= 1;
	Vec3D retval	= p2 + x*d;
	return retval;
}

/*! returns the positions where the line pierces the sphere
*
* @param sphere the sphere to be pierced
* @param result0 first intersection or tangent point of sphere and line
* @param result1 second intersection point of sphere and line
* @param getInsidePos numerical parameter, determines if result0/1 shall be slightly above or below sphere surface
*
* \retval 0 no intersection of line with sphere
* \retval 1 line tangential, tangential point stored in result0
* \retval 2 two intersections, stored in result0 and result1
*/
int hcLine3D::getIntersectionsWithSphere(const hcSphere<3> &sphere, Vec3D &result0, Vec3D &result1, bool getInsidePos)
{
    float eps = 1E-6;

    hcFloat lambda0, lambda1;
    int result = intersectsSphere3D(sphere, lambda0, lambda1, getInsidePos);

    if(result == 0)
    {
        result0.zero();
        result1.zero();
        return 0;
    }

    if(result == 2)
    {
        result0	= getPos(lambda0);
        result1	= getPos(lambda1);

        result0	= result0.convCoordCart2Spher();
        result1 = result1.convCoordCart2Spher();

        if(getInsidePos)
        {
            result0(0) = sphere.radius - sphere.radius * eps;
            result1(0) = sphere.radius - sphere.radius * eps;
        }
        else
        {
            result0(0) = sphere.radius + sphere.radius * eps;
            result1(0) = sphere.radius + sphere.radius * eps;
        }

        result0 = result0.convCoordSpher2Cart();
        result1 = result1.convCoordSpher2Cart();

        // --------------------------------------------------------------------------------
        // this is just a test and might be omitted for better program performance later on
        // --------------------------------------------------------------------------------
        Vec3D pos0 = result0.convCoordCart2Spher();
        Vec3D pos1 = result1.convCoordCart2Spher();

        if(getInsidePos)
        {
            if(pos0[0] >= sphere.radius)
            {
                printf("ERROR! Line3D::getIntersectionsWithSphere: Test failed! Position inside sphere requested, but pos0 is outside!\n");
                pos0.dump();
                return 0;
            }
            if(pos1[0] >= sphere.radius)
            {
                printf("ERROR! Line3D::getIntersectionsWithSphere: Test failed! Position inside sphere requested, but pos1 is outside!\n");
                pos1.dump();
                return 0;
            }
        }
        else
        {
            if(pos0[0] <= sphere.radius)
            {
                printf("ERROR! Line3D::getIntersectionsWithSphere: Test failed! Position outside sphere requested, but pos0 is inside!\n");
                pos0.dump();
                return 0;
            }
            if(pos1[0] <= sphere.radius)
            {
                printf("ERROR! Line3D::getIntersectionsWithSphere: Test failed! Position outside sphere requested, but pos1 is inside!\n");
                pos1.dump();
                return 0;
            }
        }
        // --------------------------------------------------------------------------------
        // this is just a test and might be omitted for better program performance later on
        // --------------------------------------------------------------------------------
        return 2;
    }

    if(result == 1)     // TODO: has to be amended in accordance with getInsidePos parameter
    {
        result0 = getPos(lambda0);
        result1	= Vec3D(0.0, 0.0, 0.0);
        return 1;
    }

    printf("ERROR! Line3D::getIntersectionsWithSphere(): This case (result == %i) thould not be possible!\n", result);
    return -1;
}




