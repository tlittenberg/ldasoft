# add library
add_library(util STATIC util.c util.h)

include_directories ("${PROJECT_SOURCE_DIR}/util/src/")
#include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
#target_link_libraries(lisa m ${GSL_LIBRARIES} ${HDF5_LIBRARIES})

install(TARGETS util DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

