SET(${PROJECT_NAME}_HEADERS
  bernstein.h
  bezier_polynomial_conversion.h
  curve_abc.h
  exact_cubic.h
  MathDefs.h
  polynomial.h
  spline_deriv_constraint.h
  bezier_curve.h
  cubic_spline.h
  curve_constraint.h
  quintic_spline.h
  linear_variable.h
  cubic_hermite_spline.h
  )

INSTALL(FILES
  ${${PROJECT_NAME}_HEADERS}
  DESTINATION include/curves
  )

ADD_SUBDIRECTORY(helpers)
