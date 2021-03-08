#ifndef HCCONSTANTS_H
#define HCCONSTANTS_H

#include "engine/math/hcVec.h"

//____naturals____________________________________________________________________________
#define r_sol 6.955E8                    // solar radius
const double PI                     = 3.1415926535897932384626433832795028841971693993751058209749445923;
const double AU                     = 149597870700;
const float solarSiderialRotSpeed   = 14.713 * PI/180;
const double deg2rad                = PI / 180;
const double rad2deg                = 180 / PI;
const double eclipObliquity			= 84381.448 / 3600 * PI / 180; 	// J. Laskar (1986) - ecliptic obliquity of J2000

const uint overlayNumTheta			= 50;
const uint overlayNumPhi			= 100;

extern Vec3D e3x, e3y, e3z;

#endif
