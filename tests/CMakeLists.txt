ADD_UNIT_TEST(
  curves_tests Main.cpp
  )
PKG_CONFIG_USE_DEPENDENCY(curves_tests eigen3)
TARGET_LINK_LIBRARIES(curves_tests ${Boost_LIBRARIES})
