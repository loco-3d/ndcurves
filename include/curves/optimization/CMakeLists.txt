SET(${PROJECT_NAME}_OPTIMIZATION_HEADERS
  definitions.h
  details.h
  quadratic_problem.h
  integral_cost.h
  )

INSTALL(FILES
  ${${PROJECT_NAME}_OPTIMIZATION_HEADERS}
  DESTINATION include/curves/optimization
  )
