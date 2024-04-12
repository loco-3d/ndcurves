#ifndef _STRUCT_SO3_LINEAR_H
#define _STRUCT_SO3_LINEAR_H

#include <Eigen/Geometry>
#include <boost/math/constants/constants.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/vector.hpp>

#include "MathDefs.h"
#include "constant_curve.h"
#include "curve_abc.h"

namespace ndcurves {

/// \class SO3Linear.
/// \brief Represents a linear interpolation in SO3, using the slerp method
/// provided by Eigen::Quaternion
///
///
template <typename Time = double, typename Numeric = Time, bool Safe = false>
struct SO3Linear : public curve_abc<Time, Numeric, Safe, matrix3_t, point3_t> {
  typedef Numeric Scalar;
  typedef matrix3_t point_t;
  typedef point3_t point_derivate_t;
  typedef Eigen::Quaternion<Scalar> quaternion_t;
  typedef Time time_t;
  typedef curve_abc<Time, Numeric, Safe, point_t, point_derivate_t> curve_abc_t;
  typedef constant_curve<Time, Numeric, Safe, point_derivate_t>
      curve_derivate_t;
  typedef SO3Linear<Time, Numeric, Safe> SO3Linear_t;

 public:
  /* Constructors - destructors */
  /// \brief Empty constructor. Curve obtained this way can not perform other
  /// class functions.
  ///
  SO3Linear()
      : curve_abc_t(),
        dim_(3),
        init_rot_(),
        end_rot_(),
        angular_vel_(),
        T_min_(0),
        T_max_(0) {}

  /// \brief constructor with initial and final rotation and time bounds
  SO3Linear(const quaternion_t& init_rot, const quaternion_t& end_rot,
            const time_t t_min, const time_t t_max)
      : curve_abc_t(),
        dim_(3),
        init_rot_(init_rot),
        end_rot_(end_rot),
        angular_vel_(computeAngularVelocity(init_rot.toRotationMatrix(),
                                            end_rot.toRotationMatrix(), t_min,
                                            t_max)),
        T_min_(t_min),
        T_max_(t_max) {
    safe_check();
  }

  /// \brief constructor with initial and final rotation expressed as rotation
  /// matrix and time bounds
  SO3Linear(const matrix3_t& init_rot, const matrix3_t& end_rot,
            const time_t t_min, const time_t t_max)
      : curve_abc_t(),
        dim_(3),
        init_rot_(quaternion_t(init_rot)),
        end_rot_(quaternion_t(end_rot)),
        angular_vel_(computeAngularVelocity(init_rot, end_rot, t_min, t_max)),
        T_min_(t_min),
        T_max_(t_max) {
    safe_check();
  }

  /// \brief constructor with initial and final rotation, time bounds are set to
  /// [0;1]
  SO3Linear(const quaternion_t& init_rot, const quaternion_t& end_rot)
      : curve_abc_t(),
        dim_(3),
        init_rot_(init_rot),
        end_rot_(end_rot),
        angular_vel_(computeAngularVelocity(
            init_rot.toRotationMatrix(), end_rot.toRotationMatrix(), 0., 1.)),
        T_min_(0.),
        T_max_(1.) {
    safe_check();
  }

  /// \brief constructor with initial and final rotation expressed as rotation
  /// matrix, time bounds are set to [0;1]
  SO3Linear(const matrix3_t& init_rot, const matrix3_t& end_rot)
      : curve_abc_t(),
        dim_(3),
        init_rot_(quaternion_t(init_rot)),
        end_rot_(quaternion_t(end_rot)),
        angular_vel_(computeAngularVelocity(init_rot, end_rot, 0., 1.)),
        T_min_(0.),
        T_max_(1.) {
    safe_check();
  }

  /// \brief Destructor
  virtual ~SO3Linear() {}

  SO3Linear(const SO3Linear& other)
      : dim_(other.dim_),
        init_rot_(other.init_rot_),
        end_rot_(other.end_rot_),
        angular_vel_(other.angular_vel_),
        T_min_(other.T_min_),
        T_max_(other.T_max_) {}

  point3_t computeAngularVelocity(const matrix3_t& init_rot,
                                  const matrix3_t& end_rot, const double t_min,
                                  const double t_max) {
    if (t_min == t_max) {
      return point3_t::Zero();
    } else {
      return log3(init_rot.transpose() * end_rot) / (t_max - t_min);
    }
  }

  quaternion_t computeAsQuaternion(const time_t t) const {
    if (Safe & !(T_min_ <= t && t <= T_max_)) {
      throw std::invalid_argument(
          "can't evaluate bezier curve, time t is out of range");  // TODO
    }
    if (t >= T_max_) return end_rot_;
    if (t <= T_min_) return init_rot_;
    Scalar u = (t - T_min_) / (T_max_ - T_min_);
    return init_rot_.slerp(u, end_rot_);
  }

  ///  \brief Evaluation of the SO3Linear at time t using Eigen slerp.
  ///  \param t : time when to evaluate the spline.
  ///  \return \f$x(t)\f$ point corresponding on spline at time t.
  virtual point_t operator()(const time_t t) const {
    return computeAsQuaternion(t).toRotationMatrix();
  }

  /**
   * @brief isApprox check if other and *this are approximately equals.
   * Only two curves of the same class can be approximately equals, for
   * comparison between different type of curves see isEquivalent
   * @param other the other curve to check
   * @param prec the precision threshold, default
   * Eigen::NumTraits<Numeric>::dummy_precision()
   * @return true is the two curves are approximately equals
   */
  bool isApprox(
      const SO3Linear_t& other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    return ndcurves::isApprox<Numeric>(T_min_, other.min()) &&
           ndcurves::isApprox<Numeric>(T_max_, other.max()) &&
           dim_ == other.dim() &&
           init_rot_.toRotationMatrix().isApprox(
               other.init_rot_.toRotationMatrix(), prec) &&
           end_rot_.toRotationMatrix().isApprox(
               other.end_rot_.toRotationMatrix(), prec);
  }

  virtual bool isApprox(
      const curve_abc_t* other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    const SO3Linear_t* other_cast = dynamic_cast<const SO3Linear_t*>(other);
    if (other_cast)
      return isApprox(*other_cast, prec);
    else
      return false;
  }

  virtual bool operator==(const SO3Linear_t& other) const {
    return isApprox(other);
  }

  virtual bool operator!=(const SO3Linear_t& other) const {
    return !(*this == other);
  }

  ///  \brief Evaluation of the derivative of order N of spline at time t.
  ///  \param t : the time when to evaluate the spline.
  ///  \param order : order of derivative.
  ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$ point corresponding on derivative
  ///  spline at time t.
  virtual point_derivate_t derivate(const time_t t,
                                    const std::size_t order) const {
    if ((t < T_min_ || t > T_max_) && Safe) {
      throw std::invalid_argument(
          "error in SO3_linear : time t to evaluate derivative should be in "
          "range [Tmin, Tmax] of the curve");
    }
    if (order > 1 || t > T_max_ || t < T_min_) {
      return point3_t::Zero(3);
    } else if (order == 1) {
      return angular_vel_;
    } else {
      throw std::invalid_argument("Order must be > 0 ");
    }
  }

  curve_derivate_t compute_derivate(const std::size_t order) const {
    return curve_derivate_t(derivate(T_min_, order), T_min_, T_max_);
  }

  ///  \brief Compute the derived curve at order N.
  ///  \param order : order of derivative.
  ///  \return A pointer to \f$\frac{d^Nx(t)}{dt^N}\f$ derivative order N of the
  ///  curve.
  curve_derivate_t* compute_derivate_ptr(const std::size_t order) const {
    return new curve_derivate_t(compute_derivate(order));
  }

  /*Helpers*/
  /// \brief Get dimension of curve.
  /// \return dimension of curve.
  std::size_t virtual dim() const { return dim_; };
  /// \brief Get the minimum time for which the curve is defined
  /// \return \f$t_{min}\f$ lower bound of time range.
  time_t min() const { return T_min_; }
  /// \brief Get the maximum time for which the curve is defined.
  /// \return \f$t_{max}\f$ upper bound of time range.
  time_t max() const { return T_max_; }
  /// \brief Get the degree of the curve.
  /// \return \f$degree\f$, the degree of the curve.
  virtual std::size_t degree() const { return 1; }
  matrix3_t getInitRotation() const { return init_rot_.toRotationMatrix(); }
  matrix3_t getEndRotation() const { return end_rot_.toRotationMatrix(); }
  matrix3_t getInitRotation() { return init_rot_.toRotationMatrix(); }
  matrix3_t getEndRotation() { return end_rot_.toRotationMatrix(); }

  /*Helpers*/

  /*Attributes*/
  std::size_t dim_;  // const
  quaternion_t init_rot_, end_rot_;
  point3_t angular_vel_;  // const
  time_t T_min_, T_max_;  // const
  /*Attributes*/

  // Serialization of the class
  friend class boost::serialization::access;

  template <class Archive>
  void load(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar >> BOOST_SERIALIZATION_BASE_OBJECT_NVP(curve_abc_t);
    ar >> boost::serialization::make_nvp("dim", dim_);
    matrix3_t init, end;
    ar >> boost::serialization::make_nvp("init_rotation", init);
    ar >> boost::serialization::make_nvp("end_rotation", end);
    init_rot_ = quaternion_t(init);
    end_rot_ = quaternion_t(end);
    ar >> boost::serialization::make_nvp("angular_vel", angular_vel_);
    ar >> boost::serialization::make_nvp("T_min", T_min_);
    ar >> boost::serialization::make_nvp("T_max", T_max_);
  }

  template <class Archive>
  void save(Archive& ar, const unsigned int version) const {
    if (version) {
      // Do something depending on version ?
    }
    ar << BOOST_SERIALIZATION_BASE_OBJECT_NVP(curve_abc_t);
    ar << boost::serialization::make_nvp("dim", dim_);
    matrix3_t init_matrix(getInitRotation());
    matrix3_t end_matrix(getEndRotation());
    ar << boost::serialization::make_nvp("init_rotation", init_matrix);
    ar << boost::serialization::make_nvp("end_rotation", end_matrix);
    ar << boost::serialization::make_nvp("angular_vel", angular_vel_);
    ar << boost::serialization::make_nvp("T_min", T_min_);
    ar << boost::serialization::make_nvp("T_max", T_max_);
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()

  /// \brief Log: SO3 -> so3.
  ///
  /// Pseudo-inverse of log from \f$ SO3 -> { v \in so3, ||v|| \le pi } \f$.
  ///
  /// Code from pinocchio.
  ///
  /// \param[in] R the rotation matrix.    ///
  /// \return The angular velocity vector associated to the rotation matrix.
  ///
  point3_t log3(const matrix3_t& R) {
    Scalar theta;
    static const Scalar PI_value = boost::math::constants::pi<Scalar>();

    point3_t res;
    const Scalar tr = R.trace();
    if (tr > Scalar(3))
      theta = 0;  // acos((3-1)/2)
    else if (tr < Scalar(-1))
      theta = PI_value;  // acos((-1-1)/2)
    else
      theta = acos((tr - Scalar(1)) / Scalar(2));
    if (!std::isfinite(theta)) {
      throw std::runtime_error("theta contains some NaN");
    }

    // From runs of hpp-constraints/tests/logarithm.cc: 1e-6 is too small.
    if (theta < PI_value - 1e-2) {
      const Scalar t =
          ((theta >
            std::pow(std::numeric_limits<Scalar>::epsilon(),
                     Scalar(1) /
                         Scalar(4)))  // precision of taylor serie of degree 3
               ? theta / sin(theta)
               : Scalar(1)) /
          Scalar(2);
      res(0) = t * (R(2, 1) - R(1, 2));
      res(1) = t * (R(0, 2) - R(2, 0));
      res(2) = t * (R(1, 0) - R(0, 1));
    } else {
      // 1e-2: A low value is not required since the computation is
      // using explicit formula. However, the precision of this method
      // is the square root of the precision with the antisymmetric
      // method (Nominal case).
      const Scalar cphi = cos(theta - PI_value);
      const Scalar beta = theta * theta / (Scalar(1) + cphi);
      point3_t tmp((R.diagonal().array() + cphi) * beta);
      res(0) = (R(2, 1) > R(1, 2) ? Scalar(1) : Scalar(-1)) *
               (tmp[0] > Scalar(0) ? sqrt(tmp[0]) : Scalar(0));
      res(1) = (R(0, 2) > R(2, 0) ? Scalar(1) : Scalar(-1)) *
               (tmp[1] > Scalar(0) ? sqrt(tmp[1]) : Scalar(0));
      res(2) = (R(1, 0) > R(0, 1) ? Scalar(1) : Scalar(-1)) *
               (tmp[2] > Scalar(0) ? sqrt(tmp[2]) : Scalar(0));
    }

    return res;
  }

 private:
  void safe_check() {
    if (Safe) {
      if (T_min_ > T_max_) {
        throw std::invalid_argument("Tmin should be inferior to Tmax");
      }
    }
  }

};  // struct SO3Linear

}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(
    SINGLE_ARG(typename Time, typename Numeric, bool Safe),
    SINGLE_ARG(ndcurves::SO3Linear<Time, Numeric, Safe>))

#endif  // _STRUCT_SO3_LINEAR_H
