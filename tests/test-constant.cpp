#define BOOST_TEST_MODULE test_constant

#include <boost/test/included/unit_test.hpp>

#include "ndcurves/constant_curve.h"
#include "ndcurves/fwd.h"
#include "ndcurves/serialization/curves.hpp"

using namespace ndcurves;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(constructors) {
  pointX_t p(5);
  p << -1, 0.5, 2., 1.2, -3.6;

  // default constructor
  constant_t cX(p);
  BOOST_CHECK_EQUAL(cX.value_, p);
  BOOST_CHECK_EQUAL(cX.min(), 0.);
  BOOST_CHECK_EQUAL(cX.max(), std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(cX.dim(), 5);
  BOOST_CHECK_EQUAL(cX.degree(), 0);

  // constructor with time interval:
  constant_t cT(p, 1., 3.5);
  BOOST_CHECK_EQUAL(cT.value_, p);
  BOOST_CHECK_EQUAL(cT.min(), 1.);
  BOOST_CHECK_EQUAL(cT.max(), 3.5);
  BOOST_CHECK_EQUAL(cT.dim(), 5);
  BOOST_CHECK_EQUAL(cT.degree(), 0);

  // tmin > tmax
  BOOST_CHECK_THROW(constant_t(p, 2., 1.), std::invalid_argument);

  // copy constructor
  constant_t cC(cT);
  BOOST_CHECK_EQUAL(cC.value_, p);
  BOOST_CHECK_EQUAL(cC.min(), 1.);
  BOOST_CHECK_EQUAL(cC.max(), 3.5);
  BOOST_CHECK_EQUAL(cC.dim(), 5);
  BOOST_CHECK_EQUAL(cC.degree(), 0);

  // Check with fixed_size point:
  point3_t p3(1, -5, 2.);
  constant_t c3(p3);
  BOOST_CHECK_EQUAL(c3.value_, p3);
  BOOST_CHECK_EQUAL(c3.min(), 0.);
  BOOST_CHECK_EQUAL(c3.max(), std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(c3.dim(), 3);
  BOOST_CHECK_EQUAL(c3.degree(), 0);
}

BOOST_AUTO_TEST_CASE(call) {
  pointX_t p(5);
  p << -1, 0.5, 2., 1.2, -3.6;

  constant_t c(p);

  BOOST_CHECK_EQUAL(c(0.), p);
  BOOST_CHECK_EQUAL(c(1.), p);
  BOOST_CHECK_EQUAL(c(0.5), p);
  BOOST_CHECK_EQUAL(c(7.9), p);

  // check outside of the bounds:
  BOOST_CHECK_THROW(c(-0.1), std::invalid_argument);

  constant_t c2(p, 1., 3.);
  BOOST_CHECK_EQUAL(c2(1.), p);
  BOOST_CHECK_EQUAL(c2(1.5), p);
  BOOST_CHECK_EQUAL(c2(3.), p);

  // check outside of the bounds:
  BOOST_CHECK_THROW(c2(0.2), std::invalid_argument);
  BOOST_CHECK_THROW(c2(7.), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(derivate) {
  pointX_t p(4);
  p << -1, 0.5, 2., 1.2;

  pointX_t p0 = pointX_t::Zero(4);
  constant_t c(p, 1., 3.);
  BOOST_CHECK_EQUAL(c.derivate(1., 1), p0);
  BOOST_CHECK_EQUAL(c.derivate(2., 2), p0);
  BOOST_CHECK_EQUAL(c.derivate(3., 15), p0);

  // check outside of the bounds:
  BOOST_CHECK_THROW(c.derivate(0.2, 1), std::invalid_argument);
  BOOST_CHECK_THROW(c.derivate(7., 1), std::invalid_argument);

  constant_t c_derivate = c.compute_derivate();
  BOOST_CHECK_EQUAL(c_derivate.min(), c.min());
  BOOST_CHECK_EQUAL(c_derivate.max(), c.max());
  BOOST_CHECK_EQUAL(c_derivate.dim(), c.dim());
  BOOST_CHECK_EQUAL(c_derivate.value_, p0);
  BOOST_CHECK_EQUAL(c_derivate(1.5), p0);

  constant_t* c_derivate_ptr = c.compute_derivate_ptr(2);
  BOOST_CHECK_EQUAL(c_derivate_ptr->min(), c.min());
  BOOST_CHECK_EQUAL(c_derivate_ptr->max(), c.max());
  BOOST_CHECK_EQUAL(c_derivate_ptr->dim(), c.dim());
  BOOST_CHECK_EQUAL(c_derivate_ptr->value_, p0);
  BOOST_CHECK_EQUAL(c_derivate_ptr->operator()(1.5), p0);

  // check with different size for the derivative:
  point3_t p03 = point3_t::Zero();
  constant_curve<double, double, true, pointX_t, point3_t> cX3(p, 1., 3.);
  BOOST_CHECK_EQUAL(cX3.derivate(1., 1), p03);
  BOOST_CHECK_EQUAL(cX3.derivate(2., 2), p03);
  BOOST_CHECK_EQUAL(cX3.derivate(3., 15), p03);

  // check outside of the bounds:
  BOOST_CHECK_THROW(cX3.derivate(0.2, 1), std::invalid_argument);
  BOOST_CHECK_THROW(cX3.derivate(7., 1), std::invalid_argument);

  constant_curve<double, double, true, point3_t> cX3_derivate =
      cX3.compute_derivate();
  BOOST_CHECK_EQUAL(cX3_derivate.min(), cX3.min());
  BOOST_CHECK_EQUAL(cX3_derivate.max(), cX3.max());
  BOOST_CHECK_EQUAL(cX3_derivate.dim(), 3);
  BOOST_CHECK_EQUAL(cX3_derivate.value_, p03);
  BOOST_CHECK_EQUAL(cX3_derivate(1.5), p03);

  constant_curve<double, double, true, point3_t>* cX3_derivate_ptr =
      cX3.compute_derivate_ptr(2);
  BOOST_CHECK_EQUAL(cX3_derivate_ptr->min(), cX3.min());
  BOOST_CHECK_EQUAL(cX3_derivate_ptr->max(), cX3.max());
  BOOST_CHECK_EQUAL(cX3_derivate_ptr->dim(), 3);
  BOOST_CHECK_EQUAL(cX3_derivate_ptr->value_, p03);
  BOOST_CHECK_EQUAL(cX3_derivate_ptr->operator()(1.5), p03);
}

BOOST_AUTO_TEST_CASE(comparison) {
  pointX_t p(4);
  p << -1, 0.5, 2., 1.2;

  constant_t c1(p, 1.5, 3.);
  constant_t c2(p, 1.5, 3.);
  constant_t c3(c1);
  curve_abc_t* c_ptr = new constant_t(p, 1.5, 3.);

  pointX_t p3(3);
  p3 << -1, 0.5, 2.;
  constant_t cn1(p3, 1.5, 3.);
  constant_t cn2(p, 1.4, 3.);
  constant_t cn3(p, 1.5, 2.);

  BOOST_CHECK(c1.isApprox(c2));
  BOOST_CHECK(c1.isEquivalent(c_ptr));
  BOOST_CHECK(c1 == c2);
  BOOST_CHECK(c1 == c3);

  BOOST_CHECK(!c1.isApprox(cn1));
  BOOST_CHECK(c1 != cn1);
  BOOST_CHECK(c1 != cn2);
  BOOST_CHECK(c1 != cn3);
}

BOOST_AUTO_TEST_CASE(serialization) {
  std::string fileName("fileTest_constant");
  pointX_t p(4);
  p << -1, 0.5, 2., 1.2;

  constant_t c(p, 1.5, 3.);

  c.saveAsText<constant_t>(fileName + ".txt");
  c.saveAsXML<constant_t>(fileName + ".xml", "constant");
  c.saveAsBinary<constant_t>(fileName);

  constant_t c_txt, c_xml, c_binary;
  c_txt.loadFromText<constant_t>(fileName + ".txt");
  c_xml.loadFromXML<constant_t>(fileName + ".xml", "constant");
  c_binary.loadFromBinary<constant_t>(fileName);

  BOOST_CHECK(c == c_txt);
  BOOST_CHECK(c == c_xml);
  BOOST_CHECK(c == c_binary);

  // } BOOST_AUTO_TEST_CASE(serialization_3) {
  // This test used to be split in two, but this raise
  // a weird "free(): invalid pointer" after test completion
  // without failures on Boost 1.72.0

  std::string fileName3("fileTest_constant3");
  point3_t p3;
  p3 << -1, 0.5, 2.;

  constant3_t c3(p3);

  c3.saveAsText<constant3_t>(fileName3 + ".txt");
  c3.saveAsXML<constant3_t>(fileName3 + ".xml", "constant");
  c3.saveAsBinary<constant3_t>(fileName3);

  constant3_t c3_txt, c3_xml, c3_binary;
  c3_txt.loadFromText<constant3_t>(fileName3 + ".txt");
  c3_xml.loadFromXML<constant3_t>(fileName3 + ".xml", "constant");
  c3_binary.loadFromBinary<constant3_t>(fileName3);

  BOOST_CHECK(c3 == c3_txt);
  BOOST_CHECK(c3 == c3_xml);
  BOOST_CHECK(c3 == c3_binary);
}

BOOST_AUTO_TEST_SUITE_END()
