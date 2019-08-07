SET(${PROJECT_NAME}_SERIALIZATION_HEADERS
  aligned-vector.hpp
  archive.hpp
  eigen-matrix.hpp
  fwd.hpp
  )

INSTALL(FILES
  ${${PROJECT_NAME}_SERIALIZATION_HEADERS}
  DESTINATION include/curves/serialization
  )