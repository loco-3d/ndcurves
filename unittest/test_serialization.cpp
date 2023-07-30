#define BOOST_TEST_MODULE test_sinusoidal

#include <boost/test/included/unit_test.hpp>

#include "ndcurves/fwd.h"
#include "ndcurves/serialization/curves.hpp"
#include "ndcurves/sinusoidal.h"

using namespace ndcurves;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(serialization) {
  std::string fileName("fileTest_sinusoidal");
  pointX_t p0(3), amp(3);
  p0 << -1, 0.5, 2.;
  amp << 2, -0.8, -1;
  double T = 1.5;
  double phi = 0.;

  sinusoidal_t c(p0, amp, T, phi, 0., 20.);

  c.saveAsText<sinusoidal_t>(fileName + ".txt");
  c.saveAsXML<sinusoidal_t>(fileName + ".xml", "sinusoidal");
  c.saveAsBinary<sinusoidal_t>(fileName);

  sinusoidal_t c_txt, c_xml, c_binary;
  c_txt.loadFromText<sinusoidal_t>(fileName + ".txt");
  c_xml.loadFromXML<sinusoidal_t>(fileName + ".xml", "sinusoidal");
  c_binary.loadFromBinary<sinusoidal_t>(fileName);

  BOOST_CHECK(c == c_txt);
  BOOST_CHECK(c == c_xml);
  BOOST_CHECK(c == c_binary);
}

BOOST_AUTO_TEST_SUITE_END()
