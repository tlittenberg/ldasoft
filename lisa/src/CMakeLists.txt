# add library
add_library(lisa STATIC LISA.c LISA.h Constants.h)

include_directories ("${PROJECT_SOURCE_DIR}/lisa/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
target_link_libraries(lisa hdf5 ${GSL_LIBRARIES})

install(TARGETS lisa DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

