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

namespace curves
{
  template <typename Point, std::size_t Dim=3>
  struct curve_constraints
  {
    typedef Point   point_t;
    curve_constraints():
    init_vel(point_t::Zero(Dim)),init_acc(init_vel),end_vel(init_vel),end_acc(init_vel){}

    ~curve_constraints(){}

    point_t init_vel;
    point_t init_acc;
    point_t end_vel;
    point_t end_acc;
  };
} // namespace curves
#endif //_CLASS_CUBICZEROVELACC
