/**
 * \file curve_constraint.h
 * \brief struct to define constraints on start / end velocities and
 * acceleration on a curve \author Steve T. \version 0.1 \date 04/05/2017
 *
 */

#ifndef _CLASS_CURVE_CONSTRAINT
#define _CLASS_CURVE_CONSTRAINT

#include <functional>
#include <vector>

#include "MathDefs.h"
#include "serialization/archive.hpp"
#include "serialization/eigen-matrix.hpp"

namespace ndcurves {
template <typename Point>
struct curve_constraints : serialization::Serializable {
  typedef Point point_t;
  curve_constraints(const size_t dim = 3)
      : init_vel(point_t::Zero(dim)),
        init_acc(point_t::Zero(dim)),
        init_jerk(point_t::Zero(dim)),
        end_vel(point_t::Zero(dim)),
        end_acc(point_t::Zero(dim)),
        end_jerk(point_t::Zero(dim)),
        dim_(dim) {}

  curve_constraints(const curve_constraints& other)
      : init_vel(other.init_vel),
        init_acc(other.init_acc),
        init_jerk(other.init_jerk),
        end_vel(other.end_vel),
        end_acc(other.end_acc),
        end_jerk(other.end_jerk),
        dim_(other.dim_) {}

  /// \brief Check if actual curve_constraints and other are equal.
  /// \param other : the other curve_constraints to check.
  /// \return true if the two curve_constraints are equals.
  virtual bool operator==(const curve_constraints& other) const {
    return dim_ == other.dim_ && init_vel == other.init_vel &&
           init_acc == other.init_acc && init_jerk == other.init_jerk &&
           end_vel == other.end_vel && end_acc == other.end_acc &&
           end_jerk == other.end_jerk;
  }

  /// \brief Check if actual curve_constraint and other are different.
  /// \param other : the other curve_constraint to check.
  /// \return true if the two curve_constraint are different.
  virtual bool operator!=(const curve_constraints& other) const {
    return !(*this == other);
  }

  virtual ~curve_constraints() {}
  point_t init_vel;
  point_t init_acc;
  point_t init_jerk;
  point_t end_vel;
  point_t end_acc;
  point_t end_jerk;
  size_t dim_;

  // Serialization of the class
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar& boost::serialization::make_nvp("init_vel", init_vel);
    ar& boost::serialization::make_nvp("init_acc", init_acc);
    ar& boost::serialization::make_nvp("init_jerk", init_jerk);
    ar& boost::serialization::make_nvp("end_vel", end_vel);
    ar& boost::serialization::make_nvp("end_acc", end_acc);
    ar& boost::serialization::make_nvp("end_jerk", end_jerk);
    ar& boost::serialization::make_nvp("dim", dim_);
  }
};
}  // namespace ndcurves
#endif  //_CLASS_CUBICZEROVELACC
