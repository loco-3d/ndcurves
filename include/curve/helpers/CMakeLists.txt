SET(${PROJECT_NAME}_HELPERS_HEADERS
  effector_spline.h
  effector_spline_rotation.h
  )

INSTALL(FILES
  ${${PROJECT_NAME}_HELPERS_HEADERS}
  DESTINATION include/curve/helpers
  )
