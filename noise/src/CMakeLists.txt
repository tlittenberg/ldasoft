# add library
add_library(noise STATIC noise_sampler.c noise_sampler.h noise_io.c noise_io.h noise_model.c noise_model.h)

include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")
include_directories ("${PROJECT_SOURCE_DIR}/data/src/")
include_directories ("${PROJECT_SOURCE_DIR}/noise/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

install(TARGETS noise DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

