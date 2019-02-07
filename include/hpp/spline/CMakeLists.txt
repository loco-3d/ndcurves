SET(${PROJECT_NAME}_HEADERS
  bernstein.h
  bezier_polynom_conversion.h
  curve_abc.h
  exact_cubic.h
  MathDefs.h
  polynom.h
  spline_deriv_constraint.h
  bezier_curve.h
  cubic_spline.h
  curve_constraint.h
  quintic_spline.h
  linear_variable.h
  )

INSTALL(FILES
  ${${PROJECT_NAME}_HEADERS}
  DESTINATION include/hpp/spline
  )

ADD_SUBDIRECTORY(helpers)
