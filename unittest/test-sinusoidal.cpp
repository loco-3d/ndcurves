#define BOOST_TEST_MODULE test_sinusoidal

#include <boost/test/included/unit_test.hpp>

#include "ndcurves/fwd.h"
#include "ndcurves/serialization/curves.hpp"
#include "ndcurves/sinusoidal.h"

using namespace ndcurves;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(constructors) {
  pointX_t p0(3), amp(3);
  p0 << -1, 0.5, 2.;
  amp << 2, -0.8, -1;
  double T = 1.5;
  double phi = 0.;

  // default constructor
  sinusoidal_t c(p0, amp, T, phi);
  BOOST_CHECK_EQUAL(c.p0_, p0);
  BOOST_CHECK_EQUAL(c.amplitude_, amp);
  BOOST_CHECK_EQUAL(c.T_, T);
  BOOST_CHECK_EQUAL(c.phi_, phi);
  BOOST_CHECK_EQUAL(c.min(), 0.);
  BOOST_CHECK_EQUAL(c.max(), std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(c.dim(), 3);
  BOOST_CHECK_EQUAL(c.degree(), 1);

  // constructor with time interval:
  sinusoidal_t cT(p0, amp, T, phi, 1., 3.5);
  BOOST_CHECK_EQUAL(cT.p0_, p0);
  BOOST_CHECK_EQUAL(cT.amplitude_, amp);
  BOOST_CHECK_EQUAL(cT.T_, T);
  BOOST_CHECK_EQUAL(cT.phi_, phi);
  BOOST_CHECK_EQUAL(cT.min(), 1.);
  BOOST_CHECK_EQUAL(cT.max(), 3.5);
  BOOST_CHECK_EQUAL(cT.dim(), 3);
  BOOST_CHECK_EQUAL(cT.degree(), 1);

  // check if errors are correctly raised
  BOOST_CHECK_THROW(sinusoidal_t(p0, amp, T, phi, 2., 1.),
                    std::invalid_argument);
  BOOST_CHECK_THROW(sinusoidal_t(p0, amp, -1., phi), std::invalid_argument);
  BOOST_CHECK_THROW(sinusoidal_t(p0, amp, 0., phi), std::invalid_argument);
  pointX_t amp4(4);
  amp4 << 2, -0.8, -1, 5.;
  BOOST_CHECK_THROW(sinusoidal_t(p0, amp4, T, phi), std::invalid_argument);

  // copy constructor
  sinusoidal_t cC(cT);
  BOOST_CHECK_EQUAL(cC.p0_, p0);
  BOOST_CHECK_EQUAL(cC.amplitude_, amp);
  BOOST_CHECK_EQUAL(cC.T_, T);
  BOOST_CHECK_EQUAL(cC.phi_, phi);
  BOOST_CHECK_EQUAL(cC.min(), 1.);
  BOOST_CHECK_EQUAL(cC.max(), 3.5);
  BOOST_CHECK_EQUAL(cC.dim(), 3);
  BOOST_CHECK_EQUAL(cC.degree(), 1);

  // Check modulo of phi:
  sinusoidal_t c2(p0, amp, T, 0.59);
  BOOST_CHECK_EQUAL(c2.phi_, 0.59);
  sinusoidal_t c3(p0, amp, T, -1.2);
  BOOST_CHECK_EQUAL(c3.phi_, -1.2);
  sinusoidal_t c4(p0, amp, T, 2. * M_PI + 0.5);
  BOOST_CHECK_EQUAL(c4.phi_, 0.5);

  // constructor from stationary points :
  pointX_t p_init(5), p_end(5);
  p_init << -1, 0.5, 2., 1.2, -3.6;
  p_end << -2, 5.6, -2, 0., 2.8;
  double duration = 0.7;

  sinusoidal_t cS(duration, p_init, p_end);
  BOOST_CHECK_EQUAL(cS.p0_, (p_init + p_end) / 2.);
  BOOST_CHECK_EQUAL(cS.amplitude_, (p_init - p_end) / 2.);
  BOOST_CHECK_EQUAL(cS.T_, duration * 2.);
  BOOST_CHECK_EQUAL(cS.phi_, M_PI / 2.);
  BOOST_CHECK_EQUAL(cS.min(), 0.);
  BOOST_CHECK_EQUAL(cS.max(), std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(cS.dim(), 5);
  BOOST_CHECK_EQUAL(cS.degree(), 1);

  sinusoidal_t cST(duration, p_init, p_end, 1., 5.);
  BOOST_CHECK_EQUAL(cST.p0_, (p_init + p_end) / 2.);
  BOOST_CHECK_EQUAL(cST.amplitude_, (p_init - p_end) / 2.);
  BOOST_CHECK_EQUAL(cST.T_, duration * 2.);
  BOOST_CHECK_EQUAL(cST.phi_, M_PI / 2.);
  BOOST_CHECK_EQUAL(cST.min(), 1.);
  BOOST_CHECK_EQUAL(cST.max(), 5.);
  BOOST_CHECK_EQUAL(cST.dim(), 5);
  BOOST_CHECK_EQUAL(cST.degree(), 1);

  // check if errors are correctly raised
  BOOST_CHECK_THROW(sinusoidal_t(duration, p_init, p_end, 1., 0.2);
                    , std::invalid_argument);
  BOOST_CHECK_THROW(sinusoidal_t(-1.2, p_init, p_end), std::invalid_argument);
  BOOST_CHECK_THROW(sinusoidal_t(0., p_init, p_end), std::invalid_argument);
  pointX_t p_init3(3);
  p_init3 << 2, -0.8, -1.;
  BOOST_CHECK_THROW(sinusoidal_t(duration, p_init3, p_end);
                    , std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(call) {
  pointX_t p0(3), amp(3);
  p0 << -1, 0.5, 2.;
  amp << 2, -0.8, -1;
  double T = 1.5;
  double phi = 0.;

  sinusoidal_t c(p0, amp, T, phi);
  BOOST_CHECK(c(0.).isApprox(p0));
  BOOST_CHECK(c(1.5).isApprox(p0));
  BOOST_CHECK(c(6.).isApprox(p0));
  BOOST_CHECK(c(T / 2.).isApprox(p0));
  BOOST_CHECK(c(T / 4.).isApprox(p0 + amp));
  BOOST_CHECK(c(3. * T / 4.).isApprox(p0 - amp));

  // check outside of the bounds:
  BOOST_CHECK_THROW(c(-0.1), std::invalid_argument);

  sinusoidal_t c2(p0, amp, T, phi, 1., 3.);
  // check outside of the bounds:
  BOOST_CHECK_THROW(c2(0.2), std::invalid_argument);
  BOOST_CHECK_THROW(c2(7.), std::invalid_argument);

  phi = 0.3;
  double two_pi_phi = T * phi / (2. * M_PI);
  sinusoidal_t c3(p0, amp, T, phi);
  BOOST_CHECK(c3(1.5 - two_pi_phi).isApprox(p0));
  BOOST_CHECK(c3(6. - two_pi_phi).isApprox(p0));
  BOOST_CHECK(c3(3. * T / 4. - two_pi_phi).isApprox(p0 - amp));

  // check from stationary points :
  pointX_t p_init(5), p_end(5);
  p_init << -1, 0.5, 2., 1.2, -3.6;
  p_end << -2, 5.6, -2, 0., 2.8;
  double duration = 0.7;

  sinusoidal_t cS(duration, p_init, p_end);
  BOOST_CHECK(cS(0.).isApprox(p_init));
  BOOST_CHECK(cS(duration).isApprox(p_end));
  BOOST_CHECK(cS(duration * 4.).isApprox(p_init));
  BOOST_CHECK(cS(duration * 7.).isApprox(p_end));
  BOOST_CHECK(cS(duration / 2.).isApprox(cS.p0_));
  BOOST_CHECK(cS(3. * duration / 2.).isApprox(cS.p0_));
}

BOOST_AUTO_TEST_CASE(derivate) {
  pointX_t p0(3), amp(3);
  p0 << -1, 0.5, 2.;
  amp << 2, -0.8, -1;
  double T = 1.5;
  double phi = 0.;

  sinusoidal_t c(p0, amp, T, phi, 0., 20.);
  BOOST_CHECK(c.derivate(0., 1).isApprox(amp * 2. * M_PI / T));
  BOOST_CHECK(c.derivate(T, 1).isApprox(amp * 2. * M_PI / T));
  BOOST_CHECK(c.derivate(T / 2., 1).isApprox(-amp * 2. * M_PI / T));
  BOOST_CHECK(c.derivate(T / 4., 1).isZero());
  BOOST_CHECK(c.derivate(3. * T / 4., 1).isZero());

  BOOST_CHECK(
      c.derivate(T / 4., 2).isApprox(-(amp * 2. * M_PI / T) * (2. * M_PI / T)));
  BOOST_CHECK(c.derivate(3. * T / 4., 2)
                  .isApprox((amp * 2. * M_PI / T) * (2. * M_PI / T)));
  BOOST_CHECK(c.derivate(0, 2).isZero());
  BOOST_CHECK(c.derivate(T, 2).isZero());
  BOOST_CHECK(c.derivate(2. * T, 2).isZero());

  phi = 0.3;
  double two_pi_phi = T * phi / (2. * M_PI);
  sinusoidal_t c2(p0, amp, T, phi, 0., 20.);
  BOOST_CHECK(c2.derivate(T - two_pi_phi, 1).isApprox(amp * 2. * M_PI / T));
  BOOST_CHECK(
      c2.derivate(3. * T / 2. - two_pi_phi, 1).isApprox(-amp * 2. * M_PI / T));
  BOOST_CHECK(c2.derivate(3. * T / 4. - two_pi_phi, 1).isZero());

  // check outside of the bounds:
  BOOST_CHECK_THROW(c.derivate(-0.1, 1), std::invalid_argument);
  BOOST_CHECK_THROW(c.derivate(21., 1), std::invalid_argument);
  BOOST_CHECK_THROW(c.derivate(3., 0), std::invalid_argument);

  for (size_t i = 1; i < 10; ++i) {
    sinusoidal_t c_derivate = c2.compute_derivate(i);
    BOOST_CHECK_EQUAL(c_derivate.min(), c2.min());
    BOOST_CHECK_EQUAL(c_derivate.max(), c2.max());
    BOOST_CHECK_EQUAL(c_derivate.dim(), c2.dim());
    BOOST_CHECK(c_derivate(0).isApprox(c2.derivate(0., i)));
    BOOST_CHECK(c_derivate(0.2).isApprox(c2.derivate(0.2, i)));
    BOOST_CHECK(c_derivate(0.9).isApprox(c2.derivate(0.9, i)));
    BOOST_CHECK(c_derivate(T).isApprox(c2.derivate(T, i)));
    BOOST_CHECK(c_derivate(T * 1.5).isApprox(c2.derivate(T * 1.5, i)));
  }

  sinusoidal_t* c_derivate_ptr = c.compute_derivate_ptr(1);
  BOOST_CHECK_EQUAL(c_derivate_ptr->min(), c.min());
  BOOST_CHECK_EQUAL(c_derivate_ptr->max(), c.max());
  BOOST_CHECK_EQUAL(c_derivate_ptr->dim(), c.dim());
  BOOST_CHECK(c_derivate_ptr->operator()(0.).isApprox(amp * 2. * M_PI / T));
}

BOOST_AUTO_TEST_CASE(comparison) {
  pointX_t p0(3), amp(3);
  p0 << -1, 0.5, 2.;
  amp << 2, -0.8, -1;
  double T = 1.5;
  double phi = 0.;

  sinusoidal_t c01(p0, amp, T, phi);
  sinusoidal_t c02(p0, amp, T, phi);
  BOOST_CHECK(c01 == c02);

  sinusoidal_t c1(p0, amp, T, phi, 0., 20.);
  sinusoidal_t c2(p0, amp, T, phi, 0., 20.);
  sinusoidal_t c3(c1);
  curve_abc_t* c_ptr = new sinusoidal_t(p0, amp, T, phi, 0., 20.);

  pointX_t p02(3), amp2(3);
  p02 << -1, 0.4, 2.;
  amp2 << -2, -0.8, -1;
  sinusoidal_t cn0(p0, amp, T, phi, 0.5, 20.);
  sinusoidal_t cn1(p0, amp, T, phi, 0., 10.);
  sinusoidal_t cn2(p0, amp, T, phi);
  sinusoidal_t cn3(p02, amp, T, phi, 0., 20.);
  sinusoidal_t cn4(p0, amp2, T, phi, 0., 20.);
  sinusoidal_t cn5(p0, amp, 1.7, phi, 0., 20.);
  sinusoidal_t cn6(p0, amp, T, 0.3, 0., 20.);

  BOOST_CHECK(c1.isApprox(c2));
  BOOST_CHECK(c1.isEquivalent(c_ptr));
  BOOST_CHECK(c1 == c2);
  BOOST_CHECK(c1 == c3);

  BOOST_CHECK(!c1.isApprox(cn1));
  BOOST_CHECK(c1 != cn0);
  BOOST_CHECK(c1 != cn1);
  BOOST_CHECK(c1 != cn2);
  BOOST_CHECK(c1 != cn3);
  BOOST_CHECK(c1 != cn4);
  BOOST_CHECK(c1 != cn5);
  BOOST_CHECK(c1 != cn6);
}

BOOST_AUTO_TEST_SUITE_END()
