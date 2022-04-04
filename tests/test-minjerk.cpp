#define BOOST_TEST_MODULE test_minjerk

#include <boost/test/included/unit_test.hpp>

#include "ndcurves/fwd.h"
#include "ndcurves/polynomial.h"

using namespace ndcurves;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(without_timing) {
  point3_t a(1.5, -2, 3.7);
  point3_t b(2, 3, -4.);

  polynomial_t c = polynomial_t::MinimumJerk(a, b);
  BOOST_CHECK_EQUAL(c.min(), 0.);
  BOOST_CHECK_EQUAL(c.max(), 1.);
  BOOST_CHECK_EQUAL(c.dim(), 3);
  BOOST_CHECK_EQUAL(c.degree(), 5);

  BOOST_CHECK(c(0.).isApprox(a));
  BOOST_CHECK(c(1.).isApprox(b));
}

BOOST_AUTO_TEST_CASE(with_timing) {
  point3_t a(1.5, -2, 3.7);
  point3_t b(2, 3, -4.);

  polynomial_t c = polynomial_t::MinimumJerk(a, b, 1., 5.);
  BOOST_CHECK_EQUAL(c.min(), 1.);
  BOOST_CHECK_EQUAL(c.max(), 5.);
  BOOST_CHECK_EQUAL(c.dim(), 3);
  BOOST_CHECK_EQUAL(c.degree(), 5);

  BOOST_CHECK(c(1.).isApprox(a));
  BOOST_CHECK(c(5.).isApprox(b));
}

BOOST_AUTO_TEST_CASE(constructor_error) {
  pointX_t a(3);
  a << 1.5, -2, 3.7;
  pointX_t b(4);
  b << 2, 3, -4., 5.;

  BOOST_CHECK_THROW(polynomial_t::MinimumJerk(a, b);, std::invalid_argument);
  point3_t a1(1.5, -2, 3.7);
  point3_t b1(2, 3, -4.);
  BOOST_CHECK_THROW(polynomial_t::MinimumJerk(a1, b1, 1., 0.5);
                    , std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(evaluate) {
  point3_t a(1.5, -2, 3.7);
  point3_t b(2, 3, -4.);

  polynomial_t c = polynomial_t::MinimumJerk(a, b, 0., 1.5);
  BOOST_CHECK(c.derivate(0., 1).isZero());
  BOOST_CHECK(c.derivate(0., 2).isZero());
  BOOST_CHECK(c.derivate(1.5, 1).isZero());
  BOOST_CHECK(c.derivate(1.5, 2).isZero());
  BOOST_CHECK(c.derivate(1., 6).isZero());
}

BOOST_AUTO_TEST_SUITE_END()
