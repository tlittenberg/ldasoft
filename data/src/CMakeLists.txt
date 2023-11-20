# add library
add_library(data STATIC lisa.c lisa.h data.h data.c constants.h)

include_directories ("${PROJECT_SOURCE_DIR}/data/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

install(TARGETS data DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

