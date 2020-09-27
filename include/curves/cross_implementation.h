/**
 * \file cubic_hermite_spline.h
 * \brief class allowing to create a cubic hermite spline of any dimension.
 * \author Justin Carpentier <jcarpent@laas.fr> modified by Jason Chemin <jchemin@laas.fr>
 * \date 05/2019
 */

#ifndef _CLASS_CROSSIMP
#define _CLASS_CROSSIMP

#include "curves/fwd.h"

namespace curves {
    template<typename Point>
Eigen::Vector3d cross(const Point& a, const Point& b){
    Eigen::Vector3d c;
    c << a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0] ;
    return c;
}
}  // namespace curves
#endif  //_CLASS_CROSSIMP
