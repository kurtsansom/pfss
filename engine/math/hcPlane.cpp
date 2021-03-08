#include "engine/math/hcPlane.h"
#include "engine/math/hcLine.h"

hcPlane3D::hcPlane3D(const Vec3D &pos, const Vec3D &normal){

    initNULL();
    init(pos, normal);
}

hcPlane3D::hcPlane3D(const hcPlane3D &other){

    *this = other;
}

hcPlane3D::hcPlane3D(const Vec3D &pos0, const Vec3D &pos1, const Vec3D &pos2){

    initNULL();
    init(pos0, pos1, pos2);
}

hcPlane3D &hcPlane3D::operator=(const hcPlane3D &other){

    if(this == &other)
        return *this;

    this->pos   = other.pos;
    this->normal= other.normal;

    return *this;
}

void hcPlane3D::clear(){

    initNULL();
}

void hcPlane3D::initNULL(){

    pos.zero();
    normal.zero();
}

bool hcPlane3D::operator==(const hcPlane3D &other)
{
    printf("ERROR! hcPlane3D::operator== not implemented!\n");
    return false;
}

bool hcPlane3D::init(const Vec3D &pos, const Vec3D &normal)
{
	this->normal 	= normal;
	this->pos 		= pos;

    if(normal.isNullVector())
    {
        printf("ERROR! hcPlane3D::set(const Vec3D &normal): Normal is Null-Vector!\n");
        return false;
    }

    this->normal.normalize();
    return true;
}

bool hcPlane3D::init(const Vec3D &pos0, const Vec3D &pos1, const Vec3D &pos2)
{
    Vec3D distVec0 = pos1 - pos0;
    Vec3D distVec1 = pos2 - pos0;

    hcLine<3> line(pos0, distVec0);

    float result;
    int isOnLine = line.intersectsP(pos2, result);
    if(isOnLine == 1)
    {
        printf("Line:\n");
        line.dump();
        printf("Pos0:\n");
        pos0.dump();
        printf("DistVec:\n");
        distVec0.dump();
        printf("pos2:\n");
        pos2.dump();
        printf("Result: %E\n", result);
        exit(1);

        printf("hcPlane3D::set(Vec3D, Vec3D, Vec3D): Points are collinear!\n");
        printf("Point0:\n");
        pos0.dump();
        Vec3D pos0t = pos0.convCoordCart2Spher();
        pos0t.dump();
        printf("Point1:\n");
        pos1.dump();
        Vec3D pos1t = pos1.convCoordCart2Spher();
        pos1t.dump();
        printf("Point2:\n");
        pos2.dump();
        Vec3D pos2t = pos2.convCoordCart2Spher();
        pos2t.dump();
        return false;
    }

    distVec0.cp(distVec1, &normal);
    normal.normalize();

    hcFloat angle 	= distVec0.getAngle(normal);
    this->pos		= pos0;

    return true;
}

void hcPlane3D::dump() const{

    printf("Dumping hcPlane:\n");
    printf("Position:\n");
    pos.dump();
    printf("Normal:\n");
    normal.dump();
}
