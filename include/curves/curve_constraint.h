/**
 * \file curve_constraint.h
 * \brief struct to define constraints on start / end velocities and acceleration
 * on a curve
 * \author Steve T.
 * \version 0.1
 * \date 04/05/2017
 *
 */

#ifndef _CLASS_CURVE_CONSTRAINT
#define _CLASS_CURVE_CONSTRAINT

#include "MathDefs.h"

#include <functional>
#include <vector>

namespace curves {
template <typename Point>
struct curve_constraints {
  typedef Point point_t;
  curve_constraints() {}

  curve_constraints(const curve_constraints& other)
    :init_vel (other.init_vel )
    ,init_acc (other.init_acc )
    ,init_jerk(other.init_jerk)
    ,end_vel  (other.end_vel  )
    ,end_acc  (other.end_acc  )
    ,end_jerk (other.end_jerk )
  {}


   ~curve_constraints(){}
    point_t init_vel;
    point_t init_acc;
    point_t init_jerk;
    point_t end_vel;
    point_t end_acc;
    point_t end_jerk;
};
}  // namespace curves
#endif  //_CLASS_CUBICZEROVELACC

