/**
 * \file sinusoidal.h
 * \brief class allowing to create a sinusoidal curve.
 * \author Pierre Fernbach
 * \version 0.4
 * \date 29/04/2020
 */

#ifndef _CLASS_SINUSOIDALCURVE
#define _CLASS_SINUSOIDALCURVE

#include <cmath>

#include "curve_abc.h"

namespace ndcurves {
/// \class sinusoidal.
/// \brief Represents a sinusoidal curve, evaluating the following equation:
///  p0 + amplitude * (sin(2pi/T + phi)
///
template <typename Time = double, typename Numeric = Time, bool Safe = false,
          typename Point = Eigen::Matrix<Numeric, Eigen::Dynamic, 1> >
struct sinusoidal : public curve_abc<Time, Numeric, Safe, Point> {
  typedef Point point_t;
  typedef Point point_derivate_t;
  typedef Time time_t;
  typedef Numeric num_t;
  typedef sinusoidal<Time, Numeric, Safe, Point> sinusoidal_t;
  typedef curve_abc<Time, Numeric, Safe, Point> curve_abc_t;  // parent class

  /* Constructors - destructors */
 public:
  /// \brief Empty constructor. Curve obtained this way can not perform other
  /// class functions.
  ///
  sinusoidal() : T_min_(0), T_max_(0), dim_(0) {}

  /// \brief Constructor
  /// \param p0       : Offset of the sinusoidal
  /// \param amplitude: Amplitude
  /// \param T        : The period
  /// \param phi      : the phase
  /// \param T_min    : lower bound of the time interval (default to 0)
  /// \param T_max    : upper bound of the time interval (default to +inf)
  ///
  sinusoidal(const Point& p0, const Point& amplitude, const time_t T,
             const time_t phi, const time_t T_min = 0.,
             const time_t T_max = std::numeric_limits<time_t>::max())
      : p0_(p0),
        amplitude_(amplitude),
        T_(T),
        phi_(std::fmod(phi, 2. * M_PI)),
        T_min_(T_min),
        T_max_(T_max),
        dim_(p0_.size()) {
    if (Safe && T_min_ > T_max_) {
      throw std::invalid_argument(
          "can't create constant curve: min bound is higher than max bound");
    }
    if (T_ <= 0.)
      throw std::invalid_argument("The period must be strictly positive");
    if (static_cast<size_t>(amplitude_.size()) != dim_)
      throw std::invalid_argument(
          "The offset and the amplitude must have the same dimension");
  }

  /// \brief Constructor from stationary points
  /// \param traj_time: duration to go from p_init to p_final (half a period)
  /// \param p_init   : first stationary point, either minimum or maximum
  /// \param p_final  : second stationary point, either minimum or maximum
  /// \param T_min    : lower bound of the time interval (default to 0)
  /// \param T_max    : upper bound of the time interval (default to +inf)
  ///
  sinusoidal(const time_t traj_time, const Point& p_init, const Point& p_final,
             const time_t T_min = 0.,
             const time_t T_max = std::numeric_limits<time_t>::max())
      : T_(2. * traj_time),
        phi_(M_PI / 2.),
        T_min_(T_min),
        T_max_(T_max),
        dim_(p_init.size()) {
    if (Safe && T_min_ > T_max_) {
      throw std::invalid_argument(
          "can't create constant curve: min bound is higher than max bound");
    }
    if (T_ <= 0)
      throw std::invalid_argument("The period must be strictly positive");
    if (p_init.size() != p_final.size())
      throw std::invalid_argument(
          "The two stationary points must have the same dimension");
    p0_ = (p_init + p_final) / 2.;
    amplitude_ = (p_init - p_final) / 2.;
  }

  /// \brief Copy constructor
  /// \param other
  sinusoidal(const sinusoidal_t& other)
      : p0_(other.p0_),
        amplitude_(other.amplitude_),
        T_(other.T_),
        phi_(other.phi_),
        T_min_(other.T_min_),
        T_max_(other.T_max_),
        dim_(other.dim_) {}

  /// \brief Destructor.
  virtual ~sinusoidal() {}
  /* Constructors - destructors */

  /*Operations*/
  ///  \brief Evaluation of the cubic spline at time t.
  ///  \param t : time when to evaluate the spine
  ///  \return \f$x(t)\f$, point corresponding on curve at time t.
  virtual point_t operator()(const time_t t) const {
    if (Safe && (t < T_min_ || t > T_max_)) {
      throw std::invalid_argument(
          "error in sinusoidal curve : time t to evaluate should be in range "
          "[Tmin, Tmax] of the curve");
    }
    return p0_ + amplitude_ * sin(two_pi_f(t) + phi_);
  }

  /// \brief Evaluate the derivative of order N of curve at time t.
  /// \param t : time when to evaluate the spline.
  /// \param order : order of derivative.
  /// \return \f$\frac{d^Nx(t)}{dt^N}\f$, point corresponding on derivative
  /// curve of order N at time t.
  virtual point_derivate_t derivate(const time_t t,
                                    const std::size_t order) const {
    if (Safe && (t < T_min_ || t > T_max_)) {
      throw std::invalid_argument(
          "error in constant curve : time t to derivate should be in range "
          "[Tmin, Tmax] of the curve");
    }
    if (order <= 0)
      throw std::invalid_argument("Order must be strictly positive");
    return amplitude_ * pow(2. * M_PI / T_, static_cast<num_t>(order)) *
           sin(two_pi_f(t) + phi_ + (M_PI * static_cast<num_t>(order) / 2.));
  }

  ///  \brief Compute the derived curve at order N.
  ///  Computes the derivative order N, \f$\frac{d^Nx(t)}{dt^N}\f$ of bezier
  ///  curve of parametric equation x(t). \param order : order of derivative.
  ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$ derivative order N of the curve.
  sinusoidal_t compute_derivate(const std::size_t order) const {
    if (order <= 0)
      throw std::invalid_argument("Order must be strictly positive");
    const point_t amplitude =
        amplitude_ * pow(2. * M_PI / T_, static_cast<num_t>(order));
    const time_t phi = phi_ + (M_PI * static_cast<num_t>(order) / 2.);
    return sinusoidal_t(point_t::Zero(dim_), amplitude, T_, phi, T_min_,
                        T_max_);
  }

  ///  \brief Compute the derived curve at orderN.
  ///  \param order : order of derivative.
  ///  \return A pointer to \f$\frac{d^Nx(t)}{dt^N}\f$ derivative order N of the
  ///  curve.
  virtual sinusoidal_t* compute_derivate_ptr(const std::size_t order) const {
    return new sinusoidal_t(compute_derivate(order));
  }

  /**
   * @brief isApprox check if other and *this are approximately equals given a
   * precision threshold Only two curves of the same class can be approximately
   * equals, for comparison between different type of curves see isEquivalent.
   * @param other the other curve to check
   * @param prec the precision threshold, default
   * Eigen::NumTraits<Numeric>::dummy_precision()
   * @return true is the two curves are approximately equals
   */
  virtual bool isApprox(
      const sinusoidal_t& other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    return ndcurves::isApprox<time_t>(T_min_, other.min()) &&
           ndcurves::isApprox<time_t>(T_max_, other.max()) &&
           dim_ == other.dim() && p0_.isApprox(other.p0_, prec) &&
           amplitude_.isApprox(other.amplitude_, prec) &&
           ndcurves::isApprox<time_t>(T_, other.T_) &&
           ndcurves::isApprox<time_t>(phi_, other.phi_);
  }

  virtual bool isApprox(
      const curve_abc_t* other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    const sinusoidal_t* other_cast = dynamic_cast<const sinusoidal_t*>(other);
    if (other_cast)
      return isApprox(*other_cast, prec);
    else
      return false;
  }

  virtual bool operator==(const sinusoidal_t& other) const {
    return isApprox(other);
  }

  virtual bool operator!=(const sinusoidal_t& other) const {
    return !(*this == other);
  }

  /*Helpers*/
  /// \brief Get dimension of curve.
  /// \return dimension of curve.
  std::size_t virtual dim() const { return dim_; }
  /// \brief Get the minimum time for which the curve is defined
  /// \return \f$t_{min}\f$ lower bound of time range.
  num_t virtual min() const { return T_min_; }
  /// \brief Get the maximum time for which the curve is defined.
  /// \return \f$t_{max}\f$ upper bound of time range.
  num_t virtual max() const { return T_max_; }
  /// \brief Get the degree of the curve.
  /// \return \f$degree\f$, the degree of the curve.
  virtual std::size_t degree() const { return 1; }
  /*Helpers*/

  /*Attributes*/
  Point p0_;  // offset
  Point amplitude_;
  time_t T_;              // period
  time_t phi_;            // phase
  time_t T_min_, T_max_;  // const
  std::size_t dim_;       // const
  /*Attributes*/

  // Serialization of the class
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(curve_abc_t);
    ar& boost::serialization::make_nvp("p0", p0_);
    ar& boost::serialization::make_nvp("amplitude_", amplitude_);
    ar& boost::serialization::make_nvp("T_", T_);
    ar& boost::serialization::make_nvp("phi_", phi_);
    ar& boost::serialization::make_nvp("T_min", T_min_);
    ar& boost::serialization::make_nvp("T_max", T_max_);
    ar& boost::serialization::make_nvp("dim", dim_);
  }

 private:
  inline const num_t two_pi_f(const time_t& t) const {
    return (2 * M_PI / T_) * t;
  }

};  // struct sinusoidal
}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(
    SINGLE_ARG(typename Time, typename Numeric, bool Safe, typename Point),
    SINGLE_ARG(ndcurves::sinusoidal<Time, Numeric, Safe, Point>))

#endif  // _CLASS_SINUSOIDALCURVE
