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
inline
Eigen::Vector3d cross(const Eigen::VectorXd& a, const Eigen::VectorXd& b){
    Eigen::Vector3d c;
    c << a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0] ;
    return c;
}

inline
Eigen::Vector3d cross(const Eigen::Vector3d& a, const Eigen::Vector3d& b){
    return a.cross(b);
}

inline
Eigen::Vector3f cross(const Eigen::Vector3f& a, const Eigen::Vector3f& b){
    return a.cross(b);
}

template<typename N, bool S>
linear_variable<N,S> cross(const linear_variable<N,S>& a, const linear_variable<N,S>& b){
    return a.cross(b);
}
}  // namespace curves
#endif  //_CLASS_CROSSIMP
