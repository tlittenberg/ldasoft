# add library
add_library(glass_utils STATIC 
            glass_lisa.c glass_lisa.h 
            glass_data.h glass_data.c 
            glass_math.c glass_math.h 
            glass_gmm.c glass_gmm.h 
            glass_wavelet.c glass_wavelet.h 
            pix2ang.c pix2ang.h
            astrometry/astrometry.c astrometry/astrometry.h
            glass_galaxy.c glass_galaxy.h 
            glass_constants.h glass_utils.c glass_utils.h) 

include_directories ("${PROJECT_SOURCE_DIR}/utils/src/")
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS} ${BLAS_INCLUDE_DIRS})

install(TARGETS glass_utils DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

