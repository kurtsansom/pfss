#ifndef PLANES_H
#define PLANES_H

#include "engine/math/hcVec.h"

class hcPlaneND{
public:

};

class hcPlane3D{
public:

    Vec3D pos;
    Vec3D normal;

    hcPlane3D(const Vec3D &pos = Vec3D(0.0, 0.0, 0.0),
              const Vec3D &normal = Vec3D(0.0, 0.0, 1.0));  /*!< \brief std constructor     */

    hcPlane3D(const hcPlane3D &other);                      /*!< \brief cpy constructor     */

    hcPlane3D(const Vec3D &pos0, const Vec3D &pos1, const Vec3D &pos2);

    hcPlane3D &operator=(const hcPlane3D &other);

    bool operator==(const hcPlane3D &other);
        /*!< \brief tests identity of this plane with another								*/

    void clear();
    void initNULL();

    bool init(const Vec3D &pos, const Vec3D &normal);
    	/*!< \brief initializes plane by giving position and normal vectors				*/

    bool init(const Vec3D &pos0, const Vec3D &pos1, const Vec3D &pos2);
    	/*!< \brief initializes plane by giving three points on the plane				*/

    void dump() const;
};

#endif // PLANES_H
