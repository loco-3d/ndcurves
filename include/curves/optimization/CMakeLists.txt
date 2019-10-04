SET(${PROJECT_NAME}_HELPERS_HEADERS
  definitions.h
  details.h
  quadratic_problem.h
  integral_cost.h
  )

INSTALL(FILES
  ${${PROJECT_NAME}_HELPERS_HEADERS}
  DESTINATION include/hpp/spline/optimization
  )
