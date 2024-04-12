// #ifdef CURVES_WITH_PINOCCHIO_SUPPORT

#define EIGEN_RUNTIME_NO_MALLOC
#define BOOST_TEST_MODULE test_polynomial

#include <boost/test/included/unit_test.hpp>
#include <random>

#include "ndcurves/polynomial.h"

using namespace ndcurves;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

double generateRandomNumber(double lower, double upper) {
  // Some random number
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(lower, upper);
  return dist(gen);
}

BOOST_AUTO_TEST_CASE(default_constructor) {
  polynomial1_t traj;
  BOOST_CHECK_EQUAL(traj.dim(), 0);
  BOOST_CHECK_EQUAL(traj.min(), 0.0);
  BOOST_CHECK_EQUAL(traj.max(), 1.0);
  BOOST_CHECK_EQUAL(traj.degree(), 0);
}

BOOST_AUTO_TEST_CASE(min_jerk_constructor) {
  double t_min = generateRandomNumber(0.0, 100.0);
  double t_max = t_min + generateRandomNumber(0.1, 10.0);
  double init = generateRandomNumber(-100.0, 100.0);
  double end = generateRandomNumber(-100.0, 100.0);

  polynomial1_t traj =
      polynomial1_t::MinimumJerk(point1_t(init), point1_t(end), t_min, t_max);

  // Real-time critical
  Eigen::internal::set_is_malloc_allowed(false);
  BOOST_CHECK_EQUAL(traj.dim(), 1);
  BOOST_CHECK_EQUAL(traj.min(), t_min);
  BOOST_CHECK_EQUAL(traj.max(), t_max);
  BOOST_CHECK_EQUAL(traj.degree(), 5);
  BOOST_CHECK_EQUAL(traj(traj.min())[0], init);
  BOOST_CHECK_CLOSE(traj(traj.max())[0], end, 1e-8);
  BOOST_CHECK_EQUAL(traj.derivate(traj.min(), 1)[0], 0.0);
  BOOST_CHECK_LE(std::abs(traj.derivate(traj.max(), 1)[0]), 1e-8);
  BOOST_CHECK_EQUAL(traj.derivate(traj.min(), 2)[0], 0.0);
  BOOST_CHECK_LE(std::abs(traj.derivate(traj.max(), 2)[0]), 1e-8);
  Eigen::internal::set_is_malloc_allowed(true);
}

BOOST_AUTO_TEST_SUITE_END()

// #endif  // CURVES_WITH_PINOCCHIO_SUPPORT
