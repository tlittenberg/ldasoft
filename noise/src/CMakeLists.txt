# add library
add_library(noise STATIC Noise.c Noise.h)

include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")
include_directories ("${PROJECT_SOURCE_DIR}/lisa/src/")
include_directories ("${PROJECT_SOURCE_DIR}/gbmcmc/src/")
include_directories ("${PROJECT_SOURCE_DIR}/noise/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
target_link_libraries(noise gbmcmc tools m ${GSL_LIBRARIES} ${HDF5_LIBRARIES})

install(TARGETS noise DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

add_executable(noise_spline_mcmc noise_spline_mcmc.c)

target_link_libraries(noise_spline_mcmc ${GSL_LIBRARIES})
target_link_libraries(noise_spline_mcmc m)
target_link_libraries(noise_spline_mcmc tools)
target_link_libraries(noise_spline_mcmc lisa)
target_link_libraries(noise_spline_mcmc gbmcmc)
target_link_libraries(noise_spline_mcmc noise)
target_link_libraries(noise_spline_mcmc ${HDF5_LIBRARIES})
if(APPLE)
   target_link_libraries(noise_spline_mcmc omp)
endif()
install(TARGETS noise_spline_mcmc DESTINATION bin)

add_executable(noise_mcmc noise_mcmc.c)

target_link_libraries(noise_mcmc ${GSL_LIBRARIES})
target_link_libraries(noise_mcmc m)
target_link_libraries(noise_mcmc tools)
target_link_libraries(noise_mcmc lisa)
target_link_libraries(noise_mcmc gbmcmc)
target_link_libraries(noise_mcmc noise)
target_link_libraries(noise_mcmc ${HDF5_LIBRARIES})
if(APPLE)
   target_link_libraries(noise_mcmc omp)
endif()
install(TARGETS noise_mcmc DESTINATION bin)

