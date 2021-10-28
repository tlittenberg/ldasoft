# add library
add_library(lisa STATIC LISA.c LISA.h Constants.h)

include_directories ("${PROJECT_SOURCE_DIR}/lisa/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
target_link_libraries(lisa m ${GSL_LIBRARIES} ${HDF5_LIBRARIES})

install(TARGETS lisa DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

