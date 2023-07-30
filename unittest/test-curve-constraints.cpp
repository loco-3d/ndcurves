#define BOOST_TEST_MODULE test_curve_constraints

#include <boost/test/included/unit_test.hpp>

#include "ndcurves/bezier_curve.h"
#include "ndcurves/fwd.h"

using namespace ndcurves;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(copy_constructor) {
  bezier_t::curve_constraints_t constraints(3);
  constraints.init_vel = point3_t(-1, -1, -1);
  constraints.init_acc = point3_t(-2, -2, -2);
  constraints.init_jerk = point3_t(1, 2, 3);
  constraints.end_vel = point3_t(-10, -10, -10);
  constraints.end_acc = point3_t(-20, -20, -20);
  constraints.end_jerk = point3_t(-1, -2, -3);

  bezier_t::curve_constraints_t constraints2(constraints);
  BOOST_CHECK_EQUAL(constraints.dim_, constraints2.dim_);
  BOOST_CHECK_EQUAL(constraints.init_vel, constraints2.init_vel);
  BOOST_CHECK_EQUAL(constraints.init_acc, constraints2.init_acc);
  BOOST_CHECK_EQUAL(constraints.init_jerk, constraints2.init_jerk);
  BOOST_CHECK_EQUAL(constraints.end_vel, constraints2.end_vel);
  BOOST_CHECK_EQUAL(constraints.end_vel, constraints2.end_vel);
  BOOST_CHECK_EQUAL(constraints.end_jerk, constraints2.end_jerk);
}

BOOST_AUTO_TEST_CASE(operator_equal) {
  bezier_t::curve_constraints_t constraints(3);
  constraints.init_vel = point3_t(-1, -1, -1);
  constraints.init_acc = point3_t(-2, -2, -2);
  constraints.init_jerk = point3_t(1, 2, 3);
  constraints.end_vel = point3_t(-10, -10, -10);
  constraints.end_acc = point3_t(-20, -20, -20);
  constraints.end_jerk = point3_t(-1, -2, -3);

  bezier_t::curve_constraints_t constraints2(constraints);
  BOOST_CHECK(constraints == constraints2);
  constraints2.init_vel = point3_t(1, 1, 1);
  BOOST_CHECK(constraints != constraints2);
  constraints2.init_vel = constraints.init_vel;
  constraints2.init_acc = point3_t(1, 1, 1);
  BOOST_CHECK(constraints != constraints2);
  constraints2.init_acc = constraints.init_acc;
  constraints2.init_jerk = point3_t(1, 1, 1);
  BOOST_CHECK(constraints != constraints2);
  constraints2.init_jerk = constraints.init_jerk;
  constraints2.end_vel = point3_t(1, 1, 1);
  BOOST_CHECK(constraints != constraints2);
  constraints2.end_vel = constraints.end_vel;
  constraints2.end_acc = point3_t(1, 1, 1);
  BOOST_CHECK(constraints != constraints2);
  constraints2.end_acc = constraints.end_acc;
  constraints2.end_jerk = point3_t(1, 1, 1);
  BOOST_CHECK(constraints != constraints2);

  bezier_t::curve_constraints_t constraints3(2);
  constraints.init_vel = pointX_t(2);
  constraints.init_vel << -1, -1;
  constraints.init_acc = pointX_t(2);
  constraints.init_acc << -2, -2;
  constraints.init_jerk = pointX_t(2);
  constraints.init_jerk << 1, 2;
  constraints.end_vel = pointX_t(2);
  constraints.end_vel << -10, -10;
  constraints.end_acc = pointX_t(2);
  constraints.end_acc << -20, -20;
  constraints.end_jerk = pointX_t(2);
  constraints.end_jerk << -1, -2;
  BOOST_CHECK(constraints != constraints3);
}

BOOST_AUTO_TEST_CASE(serialization) {
  bezier_t::curve_constraints_t constraints(3);
  constraints.init_vel = point3_t(-1, -1, -1);
  constraints.init_acc = point3_t(-2, -2, -2);
  constraints.init_jerk = point3_t(1, 2, 3);
  constraints.end_vel = point3_t(-10, -10, -10);
  constraints.end_acc = point3_t(-20, -20, -20);
  constraints.end_jerk = point3_t(-1, -2, -3);

  std::string fileName("curve_constraints");

  constraints.saveAsText<bezier_t::curve_constraints_t>(fileName + ".txt");
  constraints.saveAsXML<bezier_t::curve_constraints_t>(fileName + ".xml",
                                                       "curve_constraints");
  constraints.saveAsBinary<bezier_t::curve_constraints_t>(fileName);

  bezier_t::curve_constraints_t c_txt, c_xml, c_binary;
  c_txt.loadFromText<bezier_t::curve_constraints_t>(fileName + ".txt");
  c_xml.loadFromXML<bezier_t::curve_constraints_t>(fileName + ".xml",
                                                   "curve_constraints");
  c_binary.loadFromBinary<bezier_t::curve_constraints_t>(fileName);

  BOOST_CHECK(constraints == c_txt);
  BOOST_CHECK(constraints == c_xml);
  BOOST_CHECK(constraints == c_binary);
}

BOOST_AUTO_TEST_SUITE_END()
