SET(${PROJECT_NAME}_HEADERS
  bernstein.h
  curve_conversion.h
  curve_abc.h
  exact_cubic.h
  MathDefs.h
  polynomial.h
  bezier_curve.h
  cubic_spline.h
  curve_constraint.h
  quintic_spline.h
  linear_variable.h
  cubic_hermite_spline.h
  piecewise_curve.h
  )

INSTALL(FILES
  ${${PROJECT_NAME}_HEADERS}
  DESTINATION include/curves
  )

ADD_SUBDIRECTORY(helpers)
ADD_SUBDIRECTORY(proto)
