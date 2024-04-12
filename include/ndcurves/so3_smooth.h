#ifdef CURVES_WITH_PINOCCHIO_SUPPORT

#ifndef NDCURVES_SO3_SMOOTH_H
#define NDCURVES_SO3_SMOOTH_H

#include <Eigen/Geometry>
#include <boost/math/constants/constants.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/vector.hpp>
#include <pinocchio/fwd.hpp>
#include <pinocchio/spatial/explog.hpp>

#include "ndcurves/MathDefs.h"
#include "ndcurves/constant_curve.h"
#include "ndcurves/curve_abc.h"
#include "ndcurves/polynomial.h"

namespace ndcurves {

/// \class SO3Smooth.
/// \brief Represents a smooth interpolation in SO3, using an analytical
///        interpolation between 2 rotation.
///
template <typename Time = double, typename Numeric = Time, bool Safe = false>
struct SO3Smooth : public curve_abc<Time, Numeric, Safe, matrix3_t, point3_t> {
  typedef Numeric Scalar;
  typedef matrix3_t point_t;
  typedef matrix3_t_cst_ref point_t_ref;
  typedef point3_t point_derivate_t;
  typedef Time time_t;
  typedef curve_abc<Time, Numeric, Safe, point_t, point_derivate_t> curve_abc_t;
  typedef constant_curve<Time, Numeric, Safe, point_derivate_t>
      curve_derivate_t;
  typedef SO3Smooth<Time, Numeric, Safe> SO3Smooth_t;
  typedef polynomial<Time, Numeric, Safe, point1_t> min_jerk_t;

 public:
  /* Constructors - destructors */
  /// \brief Empty constructor. Curve obtained this way can not perform other
  /// class functions.
  ///
  SO3Smooth()
      : curve_abc_t(),
        init_rot_(point_t::Identity()),
        end_rot_(point_t::Identity()),
        t_min_(0.0),
        t_max_(1.0),
        min_jerk_(min_jerk_t::MinimumJerk(point1_t(0.0), point1_t(1.0), t_min_,
                                          t_max_)),
        rot_diff_(pinocchio::log3(init_rot_.transpose() * end_rot_)),
        dt_(1e-3) {}

  /// \brief constructor with initial and final rotation and time bounds
  SO3Smooth(const quaternion_t& init_quat, const quaternion_t& end_quat,
            const time_t t_min, const time_t t_max)
      : curve_abc_t(),
        t_min_(t_min),
        t_max_(t_max),
        min_jerk_(min_jerk_t::MinimumJerk(point1_t(0.0), point1_t(1.0), t_min_,
                                          t_max_)),
        dt_(1e-3) {
    quaternion_t tmp_init_quat = init_quat;
    quaternion_t tmp_end_quat = end_quat;
    tmp_init_quat.normalized();
    tmp_end_quat.normalized();
    init_rot_ = tmp_init_quat.toRotationMatrix();
    end_rot_ = tmp_end_quat.toRotationMatrix();
    rot_diff_ = pinocchio::log3(init_rot_.transpose() * end_rot_);
    safe_check();
  }

  /// \brief constructor with initial and final rotation expressed as rotation
  /// matrix and time bounds
  SO3Smooth(const matrix3_t& init_rot, const matrix3_t& end_rot,
            const time_t t_min, const time_t t_max)
      : curve_abc_t(),
        init_rot_(init_rot),
        end_rot_(end_rot),
        t_min_(t_min),
        t_max_(t_max),
        min_jerk_(min_jerk_t::MinimumJerk(point1_t(0.0), point1_t(1.0), t_min_,
                                          t_max_)),
        rot_diff_(pinocchio::log3(init_rot_.transpose() * end_rot_)),
        dt_(1e-3) {
    safe_check();
  }

  /// \brief constructor with initial and final rotation, time bounds are set to
  /// [0;1]
  SO3Smooth(const quaternion_t& init_rot, const quaternion_t& end_rot)
      : curve_abc_t(),
        t_min_(0.),
        t_max_(1.),
        min_jerk_(min_jerk_t::MinimumJerk(point1_t(0.0), point1_t(1.0), t_min_,
                                          t_max_)),
        dt_(1e-3) {
    quaternion_t init_rot_normalized = init_rot.normalized();
    quaternion_t end_rot_normalized = end_rot.normalized();
    init_rot_ = init_rot_normalized.toRotationMatrix();
    end_rot_ = end_rot_normalized.toRotationMatrix();
    rot_diff_ = pinocchio::log3(init_rot_.transpose() * end_rot_);
    safe_check();
  }

  /// \brief constructor with initial and final rotation expressed as rotation
  /// matrix, time bounds are set to [0;1]
  SO3Smooth(const matrix3_t& init_rot, const matrix3_t& end_rot)
      : curve_abc_t(),
        init_rot_(init_rot),
        end_rot_(end_rot),
        t_min_(0.),
        t_max_(1.),
        min_jerk_(min_jerk_t::MinimumJerk(point1_t(0.0), point1_t(1.0), t_min_,
                                          t_max_)),
        dt_(1e-3) {
    rot_diff_ = pinocchio::log3(init_rot_.transpose() * end_rot_);
    safe_check();
  }

  // Generate functions:

  /* Generators. */
  /// \brief Empty constructor. Curve obtained this way can not perform other
  /// class functions.
  ///
  void generate() {
    init_rot_.setIdentity();
    end_rot_.setIdentity();
    t_min_ = 0.0;
    t_max_ = 1.0;
    min_jerk_t::MinimumJerk(min_jerk_, point1_t(0.0), point1_t(1.0), t_min_,
                            t_max_);
    rot_diff_ = pinocchio::log3(init_rot_.transpose() * end_rot_);
    safe_check();
  }

  /// \brief constructor with initial and final rotation and time bounds
  void generate(const matrix3_t& init_rot, const matrix3_t& end_rot,
                const time_t t_min, const time_t t_max) {
    init_rot_ = init_rot;
    end_rot_ = end_rot;
    t_min_ = t_min;
    t_max_ = t_max;
    min_jerk_t::MinimumJerk(min_jerk_, point1_t(0.0), point1_t(1.0), t_min_,
                            t_max_);
    rot_diff_ = pinocchio::log3(init_rot_.transpose() * end_rot_);
    safe_check();
  }

  /// \brief constructor with initial and final rotation expressed as rotation
  /// matrix and time bounds
  void generate(const quaternion_t& init_rot, const quaternion_t& end_rot,
                const time_t t_min, const time_t t_max) {
    quaternion_t init_rot_normalized = init_rot.normalized();
    quaternion_t end_rot_normalized = end_rot.normalized();
    generate(init_rot_normalized.toRotationMatrix(),
             end_rot_normalized.toRotationMatrix(), t_min, t_max);
  }

  void generate(const quaternion_t& init_rot, const quaternion_t& end_rot) {
    quaternion_t init_rot_normalized = init_rot.normalized();
    quaternion_t end_rot_normalized = end_rot.normalized();
    generate(init_rot_normalized.toRotationMatrix(),
             end_rot_normalized.toRotationMatrix(), 0.0, 1.0);
  }

  /// \brief constructor with initial and final rotation expressed as rotation
  /// matrix, time bounds are set to [0;1]
  void generate(const matrix3_t& init_rot, const matrix3_t& end_rot) {
    generate(init_rot, end_rot, 0.0, 1.0);
  }

  /// \brief Destructor
  virtual ~SO3Smooth() {}

  /// \brief Copy constructor
  SO3Smooth(const SO3Smooth& other)
      : curve_abc_t(),
        init_rot_(other.init_rot_),
        end_rot_(other.end_rot_),
        t_min_(other.t_min_),
        t_max_(other.t_max_),
        min_jerk_(other.min_jerk_),
        rot_diff_(other.rot_diff_),
        dt_(1e-3) {}

  ///  \brief Evaluation of the SO3Smooth at time t using \$f R(t) = R1 *
  ///  exp3(min_jerk(t) log3(R1.T * R2)) \$f. \param t : time when to evaluate
  ///  the spline. \return \f$x(t)\f$ point corresponding on spline at time t.
  virtual point_t operator()(const time_t t) const {
    if ((t < t_min_ || t > t_max_) && Safe) {
      throw std::invalid_argument(
          "error in SO3Smooth : time t to evaluate derivative should be in "
          "range [Tmin, Tmax] of the curve");
    }
    return init_rot_ * pinocchio::exp3(min_jerk_(t)[0] * rot_diff_);
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
      const SO3Smooth_t& other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    return ndcurves::isApprox<Numeric>(t_min_, other.min()) &&
           ndcurves::isApprox<Numeric>(t_max_, other.max()) &&
           init_rot_.isApprox(other.init_rot_, prec) &&
           end_rot_.isApprox(other.end_rot_, prec);
  }

  virtual bool isApprox(
      const curve_abc_t* other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    const SO3Smooth_t* other_cast = dynamic_cast<const SO3Smooth_t*>(other);
    if (other_cast)
      return isApprox(*other_cast, prec);
    else
      return false;
  }

  virtual bool operator==(const SO3Smooth_t& other) const {
    return isApprox(other);
  }

  virtual bool operator!=(const SO3Smooth_t& other) const {
    return !(*this == other);
  }

  ///  \brief Evaluation of the derivative of order N of spline at time t.
  ///  \param t : the time when to evaluate the spline.
  ///  \param order : order of derivative.
  ///  \return \f$ \frac{d^N x(t)}{dt^N} \f$ point corresponding on derivative
  ///          spline at time t. Specifically
  ///          \f$ \frac{d R(t)}{dt} = R1 * \frac{d exp^{min_jerk(t) diff}{dw} *
  ///          \frac{d min_jerk(t)}{dt} diff \f$ with \f$ R^3 \f$ -> \f$ SO3 \f$
  ///          : \f$ R = exp(w) \f$ or \frac{d R(t)}{dt} =  = R1 *
  ///          Jexp3(min_jerk(t) diff) * \frac{d min_jerk(t)}{dt} diff \f$. We
  ///          use the pinocchio library in order to compute the exp3, log3, and
  ///          Jexp3 fonctions.
  virtual point_derivate_t derivate(const time_t t,
                                    const std::size_t order) const {
    point_derivate_t out = point_derivate_t::Zero();
    if ((t < t_min_ || t > t_max_) && Safe) {
      throw std::invalid_argument(
          "error in SO3Smooth : time t to evaluate derivative should be in "
          "range [Tmin, Tmax] of the curve");
    }
    if (order > 2) {
      throw std::invalid_argument("Order must be in {1, 2}.");
    } else if (order == 2) {
      double t1 = t;
      double t2 = t + dt_;
      if (t2 > t_max_) {
        t1 = t - dt_;
        t2 = t;
      }
      if (t1 < t_min_) {
        t1 = t_min_;
        t2 = t_max_;
      }
      if (t1 == t2) {
        out.setZero();
      } else {
        out = (derivate(t1, 1) - derivate(t2, 1)) / (t2 - t1);
      }
    } else if (order == 1) {
      matrix3_t jexp;
      pinocchio::Jexp3(min_jerk_(t)[0] * rot_diff_, jexp);
      out.noalias() =
          init_rot_ * jexp * min_jerk_.derivate(t, 1)[0] * rot_diff_;
    } else {
      throw std::invalid_argument("Order must be > 0 ");
    }
    return out;
  }

  curve_derivate_t compute_derivate(const std::size_t /*order*/) const {
    throw std::runtime_error("This function has not been implemented yet.");
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
  std::size_t virtual dim() const { return 3; };
  /// \brief Get the minimum time for which the curve is defined
  /// \return \f$t_{min}\f$ lower bound of time range.
  time_t min() const { return t_min_; }
  /// \brief Get the maximum time for which the curve is defined.
  /// \return \f$t_{max}\f$ upper bound of time range.
  time_t max() const { return t_max_; }
  /// \brief Get the degree of the curve.
  /// \return \f$degree\f$, the degree of the curve.
  virtual std::size_t degree() const { return 5; }
  matrix3_t_cst_ref get_init_rotation() const { return init_rot_; }
  matrix3_t_cst_ref get_end_rotation() const { return end_rot_; }

  /*Attributes*/
  point_t init_rot_; /*! @brief Initial rotation. */
  point_t end_rot_;  /*! @brief End rotation. */
  time_t t_min_;     /*! @brief Start time of the trajectory. */
  time_t t_max_;     /*! @brief End time of the trajectory. */
  min_jerk_t
      min_jerk_; /*! @brief Min jerk trajectory from t_min_ and t_max_. */
  point3_t
      rot_diff_;    /*! @brief The difference between init_rot_ and end_rot_. */
  const double dt_; /*! @brief Finite differentiation delta time for the
                       acceleration computation. */

  // Serialization of the class
  friend class boost::serialization::access;

  template <class Archive>
  void load(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar >> BOOST_SERIALIZATION_BASE_OBJECT_NVP(curve_abc_t);
    ar >> boost::serialization::make_nvp("init_rotation", init_rot_);
    ar >> boost::serialization::make_nvp("end_rotation", end_rot_);
    ar >> boost::serialization::make_nvp("t_min", t_min_);
    ar >> boost::serialization::make_nvp("t_max", t_max_);
    ar >> boost::serialization::make_nvp("min_jerk", min_jerk_);
    ar >> boost::serialization::make_nvp("rot_diff", rot_diff_);
  }

  template <class Archive>
  void save(Archive& ar, const unsigned int version) const {
    if (version) {
      // Do something depending on version ?
    }
    ar << BOOST_SERIALIZATION_BASE_OBJECT_NVP(curve_abc_t);
    ar << boost::serialization::make_nvp("init_rotation", init_rot_);
    ar << boost::serialization::make_nvp("end_rotation", end_rot_);
    ar << boost::serialization::make_nvp("t_min", t_min_);
    ar << boost::serialization::make_nvp("t_max", t_max_);
    ar << boost::serialization::make_nvp("min_jerk", min_jerk_);
    ar << boost::serialization::make_nvp("rot_diff", rot_diff_);
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
      if (t_min_ > t_max_) {
        throw std::invalid_argument("Tmin should be inferior to Tmax");
      }
    }
  }

};  // struct SO3Smooth

}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(
    SINGLE_ARG(typename Time, typename Numeric, bool Safe),
    SINGLE_ARG(ndcurves::SO3Smooth<Time, Numeric, Safe>))

#endif  // NDCURVES_SO3_SMOOTH_H
#endif  // CURVES_WITH_PINOCCHIO_SUPPORT
